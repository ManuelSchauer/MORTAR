"""
MORTAR - MOlecule fRagmenTAtion fRamework, test data generator (MIT licence, see the MORTAR repository).

Independent reference implementation of the RECAP fragmentation (Lewell et al., J. Chem. Inf. Comput. Sci. 1998, 38,
511-522) as implemented in the MORTAR class RECAPFragmenter: the same eleven bond rules (SMARTS), the same assignment of
bond types (lowest bond type number first), the same selection of the bonds to cleave (no small terminal groups are
cleaved off), and simultaneous cleavage. It uses the RDKit substructure matcher, so that the CDK-based implementation
can be cross-validated against an independent SMARTS engine, see RECAPFragmenterReferenceTest. Each dummy atom carries
the number of the bond type as isotope, e.g. [5*]Oc1ccc(F)cc1.

Usage: python -I RECAP_Reference_Generator.py <input .sdf|.smi> <output .tsv> [--every N] [--max N] [--rdkit]
       [--small-group-carbons N] [--comment TEXT]
  --every N              use only every N-th molecule of the input
  --max N                stop after N molecules
  --rdkit                add a column with the leaf fragments of the RECAP implementation of the RDKit
  --small-group-carbons N  maximum number of carbon atoms of small terminal groups (default 4)
  --comment              add a comment line to the header of the output, e.g. the source and licence of the input data

The reference file of the BIOFACQUIM data set in the MORTAR test resources was generated with RDKit 2025.03.3 by
  python -I RECAP_Reference_Generator.py Tutorial/biofacquim-06-2025.sdf RECAP_Reference_Fragments_BIOFACQUIM.tsv
      --comment "..."
"""
import os
import sys
from collections import deque

from rdkit import Chem
from rdkit import RDLogger
from rdkit import rdBase

RDLogger.DisableLog('rdApp.*')

AMIDE_CARBON = '[C;$(C=O);!$(C(~[#7])~[#7])]'
AMIDE_NITROGEN = '[#7;+0]'
ESTER_CARBON = '[C;$(C=O)]'
ESTER_OXYGEN = '[O;D2;$(O([#6])[#6])]'
SP3_CARBON = '[C;X4]'
AMINE_NITROGEN = '[N;+0;X3;!$(N[#6]=[#7,#8,#16]);!$(N[#7,#8,#15,#16])]'
UREA_CARBON = '[C;$(C(=O)(~[#7])~[#7])]'
ETHER_OXYGEN = '[O;D2;+0;$(O([#6])[#6]);!$(O[#6]=[#7,#8,#16])]'
OLEFIN_CARBON = '[C;!$(C(=*)=*)]'
QUATERNARY_NITROGEN = '[N+;D4;!$([N+]~[!#6])]'
ALIPHATIC_CARBON = '[C]'
LACTAM_NITROGEN = '[N;R;$(N@C=O)]'
SULFONYL_SULFUR = '[S;X4;$(S(=O)=O)]'

SINGLE = '-;!@'
DOUBLE = '=;!@'
# (bond type, name, left, bond, right)
RULES = [
    (1, 'amide', AMIDE_CARBON, SINGLE, AMIDE_NITROGEN),
    (2, 'ester', ESTER_CARBON, SINGLE, ESTER_OXYGEN),
    (3, 'amine', SP3_CARBON, SINGLE, AMINE_NITROGEN),
    (4, 'urea', AMIDE_NITROGEN, SINGLE, UREA_CARBON),
    (5, 'ether', SP3_CARBON, SINGLE, ETHER_OXYGEN),
    (6, 'olefin', OLEFIN_CARBON, DOUBLE, OLEFIN_CARBON),
    (7, 'quaternary nitrogen', QUATERNARY_NITROGEN, SINGLE, '[#6]'),
    (8, 'aromatic nitrogen - aliphatic carbon', '[n]', SINGLE, ALIPHATIC_CARBON),
    (9, 'lactam nitrogen - aliphatic carbon', LACTAM_NITROGEN, SINGLE, ALIPHATIC_CARBON),
    (10, 'aromatic carbon - aromatic carbon', '[c]', SINGLE, '[c]'),
    (11, 'sulphonamide', AMIDE_NITROGEN, SINGLE, SULFONYL_SULFUR),
]
PATTERNS = []
for bond_type, name, left, bond, right in RULES:
    patt = Chem.MolFromSmarts(left + bond + right)
    if patt is None:
        raise ValueError('cannot parse rule %d' % bond_type)
    PATTERNS.append((bond_type, patt))


def find_recap_bonds(mol):
    """returns {bond index: bond type}, the lowest bond type number of all matching rules"""
    found = {}
    for bond_type, patt in PATTERNS:
        for match in mol.GetSubstructMatches(patt, uniquify=False, maxMatches=1000000):
            idx = mol.GetBondBetweenAtoms(match[0], match[1]).GetIdx()
            if idx not in found:
                found[idx] = bond_type
    return found


def is_small_terminal_group(mol, point, cleaved, max_carbons):
    """True if the part containing the cleavage point is a small terminal group: one attachment point, and an amino
    group (NH2) for nitrogen, an acyclic hydrocarbon group with at most max_carbons carbon atoms for an alkyl carbon
    atom, or otherwise (oxygen, carbonyl, sulfonyl) the cleavage point with its double-bonded oxygen atoms and hydrogen
    or such a hydrocarbon group"""
    atom = mol.GetAtomWithIdx(point)
    head = {point}
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetIdx() not in cleaved and bond.GetBondType() == Chem.BondType.DOUBLE and other.GetAtomicNum() == 8:
            head.add(other.GetIdx())
    is_nitrogen = atom.GetAtomicNum() == 7
    is_alkyl_carbon = atom.GetAtomicNum() == 6 and len(head) == 1
    seen = {point}
    queue = deque([point])
    part = []
    attachments = 0
    while queue:
        current = mol.GetAtomWithIdx(queue.popleft())
        if current.GetAtomicNum() == 0:
            attachments += 1
        else:
            part.append(current)
        for bond in current.GetBonds():
            other = bond.GetOtherAtomIdx(current.GetIdx())
            if bond.GetIdx() in cleaved:
                attachments += 1
            elif other not in seen:
                seen.add(other)
                queue.append(other)
    if attachments != 1:
        return False
    if is_nitrogen:
        return len(part) == len(head)
    tail = [a for a in part if is_alkyl_carbon or a.GetIdx() not in head]
    return all(a.GetAtomicNum() == 6 and not a.IsInRing() for a in tail) and len(tail) <= max_carbons


def select_bonds(mol, bonds, max_carbons):
    cleaved = set()
    for idx in sorted(bonds, key=lambda i: (bonds[i], i)):
        cleaved.add(idx)
        bond = mol.GetBondWithIdx(idx)
        if is_small_terminal_group(mol, bond.GetBeginAtomIdx(), cleaved, max_carbons) \
                or is_small_terminal_group(mol, bond.GetEndAtomIdx(), cleaved, max_carbons):
            cleaved.discard(idx)
    return {idx: bonds[idx] for idx in cleaved}


def fragment(mol, bonds):
    if not bonds:
        return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mol, asMols=True))
    indices = sorted(bonds)
    labels = [(bonds[i], bonds[i]) for i in indices]
    frag_mol = Chem.FragmentOnBonds(mol, indices, addDummies=True, dummyLabels=labels)
    return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False))


def rdkit_recap(mol):
    from rdkit.Chem import Recap
    return sorted(Recap.RecapDecompose(mol).GetLeaves().keys())


def molecules(path, every, maximum):
    count = 0
    if path.endswith('.sdf'):
        for i, mol in enumerate(Chem.ForwardSDMolSupplier(path)):
            if i % every or mol is None:
                continue
            name = ''
            for key in ('_Name', 'identifier', 'chembl_id'):
                if mol.HasProp(key) and mol.GetProp(key).strip():
                    name = mol.GetProp(key).strip()
                    break
            yield name or str(i), mol
            count += 1
            if maximum and count >= maximum:
                return
    else:
        with open(path) as handle:
            for i, line in enumerate(handle):
                parts = line.split()
                if i % every or not parts:
                    continue
                mol = Chem.MolFromSmiles(parts[0])
                if mol is None:
                    continue
                yield (parts[1] if len(parts) > 1 else str(i)), mol
                count += 1
                if maximum and count >= maximum:
                    return


def main():
    args = sys.argv[1:]
    every = int(args[args.index('--every') + 1]) if '--every' in args else 1
    maximum = int(args[args.index('--max') + 1]) if '--max' in args else 0
    max_carbons = int(args[args.index('--small-group-carbons') + 1]) if '--small-group-carbons' in args else 4
    with_rdkit = '--rdkit' in args
    with open(args[1], 'w') as out:
        out.write('# RECAP reference fragments (rules of MORTAR RECAPFragmenter, small terminal groups up to %d carbon atoms)'
                  ' generated by RECAP_Reference_Generator.py with RDKit %s from %s\n'
                  % (max_carbons, rdBase.rdkitVersion, os.path.basename(args[0])))
        if '--comment' in args:
            out.write('# %s\n' % args[args.index('--comment') + 1])
        out.write('# name\tinput SMILES\tRECAP fragments\tcleaved bonds (type:begin-end)'
                  + ('\tRDKit RECAP leaves' if with_rdkit else '') + '\n')
        for name, mol in molecules(args[0], every, maximum):
            smiles = Chem.MolToSmiles(mol)
            # re-parse the written SMILES, so that both implementations start from the identical input
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            bonds = select_bonds(mol, find_recap_bonds(mol), max_carbons)
            cut = ' '.join('%d:%d-%d' % (bonds[i], mol.GetBondWithIdx(i).GetBeginAtomIdx(),
                                         mol.GetBondWithIdx(i).GetEndAtomIdx()) for i in sorted(bonds))
            row = [name, smiles, '.'.join(fragment(mol, bonds)), cut]
            if with_rdkit:
                row.append('.'.join(rdkit_recap(mol)))
            out.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    main()
