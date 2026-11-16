"""
MORTAR - MOlecule fRagmenTAtion fRamework, test data generator (MIT licence, see the MORTAR repository).

Independent reference implementation of the ORIGINAL BRICS shredding rules (Degen et al., ChemMedChem 2008, 3,
1503-1507) as specified in SMARTS syntax for the program Recore in the dissertation of J. Degen (Universitaet Hamburg,
2009), appendix C.1.1.3, "BRICS-Shredding-Regelwerk (Teil 1 und 2)". Uses the RDKit substructure matcher, so that the
CDK-based implementation in MORTAR (BRICSFragmenter) can be cross-validated against an independent SMARTS engine, see
BRICSFragmenterReferenceTest. The rules, their names, and their order are the same as in BRICSFragmenter.

Adaptations of the published SMARTS (identical in the Java implementation):
 - '!H' (meaning "not a hydrogen atom" with the explicit hydrogens of the original MOL2 input) is written as '!#1'.
 - the chirality-like primitive '!@' inside the atom brackets of rule C7 ('[C,O,S;!@]') is removed.
 - the right side of rule C7 ("repetition of the expression on the left side") is the left expression.
All bonds are simultaneously cleaved (as in Recore). Each dummy atom carries the link type of the atom it is attached
to as isotope, e.g. [16*]c1ccccc1. Rule O3 has no link type for the P/S side, its dummy atom is unlabelled ([*]).

Usage: python -I BRICS_Reference_Generator.py <input .sdf|.smi> <output .tsv> [--every N] [--max N] [--rdkit]
       [--comment TEXT]
  --every N    use only every N-th molecule of the input
  --max N      stop after N molecules
  --rdkit      add a column with the fragments of the BRICS implementation of the RDKit (BreakBRICSBonds)
  --comment    add a comment line to the header of the output, e.g. the source and licence of the input data

The reference file of the BIOFACQUIM data set in the MORTAR test resources was generated with RDKit 2025.03.3 by
  python -I BRICS_Reference_Generator.py Tutorial/biofacquim-06-2025.sdf BRICS_Reference_Fragments_BIOFACQUIM.tsv
      --comment "..."
"""
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

R13 = '[C;D3;R;$([C]@[N,n,S,s,O,o])]'
R14 = '[c;X3;R;$([c]@[N,n,S,s,O,o])]'
R15 = '[C;D3;R;!$([C]@[N,n,S,s,O,o])]'
R16 = '[c;X3;R;!$([c]@[N,n,S,s,O,o])]'
RING = [(R13, 13), (R14, 14), (R15, 15), (R16, 16)]

C_R8 = '[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]'
C_R6 = '[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]'
C_R7 = '[#6;!D1;!$([#6]-;!@[#6]=;!@[C,O,S]);$([#6][!#1;!D1]!=[!#1;!D1])]'
O1_C = '[C;!D4;!D1;!R;!$(C(=O));$(C([O;D2])[#6;!D1]~[!D1])]'
O1_O = '[O;D2;$(O([C;!D1;!R])[#6;!D1][!D1]);!$(OCO);!$(OC=*);!$(O[P,S])]'
O1R_O = '[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]'
O2_C = '[C;!D4;!D1;!R;$(C(=O));$(C([O;D2])[#6;!D1]!=[!D1])]'
O2_O = '[O;D2;$(O(C(=O))[#6;!D1][#6;!D1]);!$(OCO);!$(O[P,S])]'
O3_PS = '[P,S]'
O3_O = '[O;D2;$(O([P,S])[#6;!D1][!D1]);!$(OCO);!$(OC=*)]'
N1_C = '[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!S;!D1][!D1])]'
N1_N = '[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!S;!D1][!D1])]'
N3_S = '[S;X4;$(S(=O)(=O))]'
N3_N = '[#7;X3;!D1;$([#7](S(=O)(=O))[#6;!D1][!D1]);!$([#7]C(=O))]'
N6_C = '[C;!D4;!D1;!R;$(C(=O));$(C([n])[!S;!D1][!D1])]'
N_R9 = '[n;!$(n@[#6](=O))]'
N2_C = '[C;!D4;!D1;!R;!$(C(=O));$(C([N;X3;!D1])[#6;!D1][!D1])]'
N2_N = '[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([C;!D1;!R])[#6;!D1]!=[!D1])]'
N2R_N = '[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]'
N4_C = '[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([n])!=[#6;!D1]~[!D1])]'
N7_C = '[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([#7;R;D3])!=[#6;!D1]~[!D1])]'
N_R10 = '[#7;R;D3;$([#7]@[#6](=O))]'
N9_C = '[C;!D4;!D1;!R;$(C(=O));$(C([#7;R;D3])[!S;!D1][!D1])]'
S1_C = '[C;!D4;!D1;!R;!$(C(=S));$(C([S;D2])[#6;!D1]~[!D1])]'
S1_S = '[S;D2;$(S([C;!D1;!R])[#6;!D1][!D1]);!$(SCS);!$(SC=*);!$(S[P,S])]'
S1R_S = '[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]'
S2_C = '[C;!D4;!D1;!R;$(C(S=O)[#6;!D1][!D1])]'
S2_S = '[S;$(S(=O)([C;!D1;!R])[#6;!D1][!D1])]'
S2R_S = '[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]'

SINGLE = '-;!@'
DOUBLE = '=;!@'
# (name, left, bond, right, left link type, right link type)
RULES = []
# alkylated ring carbon atoms (named C1a, C1b, C2a, C2b in the dissertation, renamed to avoid duplicate names)
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('C1' + suffix, C_R8, SINGLE, ring, 8, label))
# acyclic bonds between ring carbon atoms, C2a-d, C3a-c, C4a-b, C5
for i, (ring1, label1) in enumerate(RING):
    for j, (ring2, label2) in enumerate(RING[i:]):
        name = 'C5' if i == 3 else 'C%d%s' % (i + 2, 'abcd'[j])
        RULES.append((name, ring1, SINGLE, ring2, label1, label2))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('C6' + suffix, C_R6, SINGLE, ring, 6, label))
RULES.append(('C7', C_R7, DOUBLE, C_R7, 7, 7))
RULES.append(('O1', O1_C, SINGLE, O1_O, 4, 3))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('O1' + suffix, ring, SINGLE, O1R_O, label, 3))
RULES.append(('O2', O2_C, SINGLE, O2_O, 1, 3))
RULES.append(('O3', O3_PS, SINGLE, O3_O, 0, 3))
RULES.append(('N1', N1_C, SINGLE, N1_N, 1, 2))
RULES.append(('N3', N3_S, SINGLE, N3_N, 12, 2))
RULES.append(('N6', N6_C, SINGLE, N_R9, 1, 9))
RULES.append(('N2', N2_C, SINGLE, N2_N, 4, 5))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('N2' + suffix, ring, SINGLE, N2R_N, label, 5 if label in (13, 15) else 2))
RULES.append(('N4', N4_C, SINGLE, N_R9, 8, 9))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('N5' + suffix, ring, SINGLE, N_R9, label, 9))
RULES.append(('N7', N7_C, SINGLE, N_R10, 8, 10))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('N8' + suffix, ring, SINGLE, N_R10, label, 10))
RULES.append(('N9', N9_C, SINGLE, N_R10, 1, 10))
RULES.append(('S1', S1_C, SINGLE, S1_S, 4, 11))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('S1' + suffix, ring, SINGLE, S1R_S, label, 11))
RULES.append(('S2', S2_C, SINGLE, S2_S, 4, 11))
for suffix, (ring, label) in zip('abcd', RING):
    RULES.append(('S2' + suffix, ring, SINGLE, S2R_S, label, 11))

PATTERNS = []
for name, left, bond, right, ll, rl in RULES:
    patt = Chem.MolFromSmarts(left + bond + right)
    if patt is None:
        raise ValueError('cannot parse rule %s' % name)
    PATTERNS.append((name, patt, ll, rl))


def find_brics_bonds(mol):
    """returns {bond index: (rule name, label of begin atom, label of end atom)} and a list of conflicts"""
    found = {}
    conflicts = []
    for name, patt, ll, rl in PATTERNS:
        for match in mol.GetSubstructMatches(patt, uniquify=False, maxMatches=1000000):
            a, b = match[0], match[1]
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond.GetBeginAtomIdx() == a:
                labels = (ll, rl)
            else:
                labels = (rl, ll)
            entry = (name, labels[0], labels[1])
            idx = bond.GetIdx()
            if idx in found:
                if (found[idx][1], found[idx][2]) != labels:
                    conflicts.append((idx, found[idx][0], name))
                continue
            found[idx] = entry
    return found, conflicts


def fragment(mol, bonds):
    if not bonds:
        return [Chem.MolToSmiles(mol)]
    indices = sorted(bonds)
    labels = [(bonds[i][2], bonds[i][1]) for i in indices]  # (label for dummy on end atom, ... on begin atom)
    frag_mol = Chem.FragmentOnBonds(mol, indices, addDummies=True, dummyLabels=labels)
    return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False))


def rdkit_brics(mol):
    from rdkit.Chem import BRICS
    res = BRICS.BreakBRICSBonds(mol)
    return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(res, asMols=True, sanitizeFrags=False))


def molecules(path, every, maximum):
    count = 0
    if path.endswith('.sdf'):
        supplier = Chem.ForwardSDMolSupplier(path)
        for i, mol in enumerate(supplier):
            if i % every:
                continue
            if mol is None:
                continue
            name = ''
            for key in ('_Name', 'identifier', 'chembl_id'):
                if mol.HasProp(key) and mol.GetProp(key).strip():
                    name = mol.GetProp(key).strip()
                    break
            name = name or str(i)
            yield name, mol
            count += 1
            if maximum and count >= maximum:
                return
    else:
        with open(path) as handle:
            for i, line in enumerate(handle):
                if i % every:
                    continue
                parts = line.split()
                if not parts:
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
    with_rdkit = '--rdkit' in args
    n_conflicts = 0
    import os
    from rdkit import rdBase
    with open(args[1], 'w') as out:
        out.write('# BRICS reference fragments (original rule set of Degen et al.) generated by BRICS_Reference_Generator.py'
                  ' with RDKit %s from %s\n' % (rdBase.rdkitVersion, os.path.basename(args[0])))
        if '--comment' in args:
            out.write('# %s\n' % args[args.index('--comment') + 1])
        out.write('# name\tinput SMILES\toriginal BRICS fragments\tcleaved bonds (rule:begin-end)'
                  + ('\tRDKit BRICS fragments' if with_rdkit else '') + '\n')
        for name, mol in molecules(args[0], every, maximum):
            smiles = Chem.MolToSmiles(mol)
            # re-parse the written SMILES, so that both implementations start from the identical input
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            bonds, conflicts = find_brics_bonds(mol)
            n_conflicts += len(conflicts)
            for c in conflicts:
                sys.stderr.write('conflict %s %s\n' % (name, c))
            frags = fragment(mol, bonds)
            cut = ' '.join('%s:%d-%d' % (bonds[i][0], mol.GetBondWithIdx(i).GetBeginAtomIdx(),
                                         mol.GetBondWithIdx(i).GetEndAtomIdx()) for i in sorted(bonds))
            row = [name, smiles, '.'.join(frags), cut]
            if with_rdkit:
                row.append('.'.join(rdkit_brics(mol)))
            out.write('\t'.join(row) + '\n')
    sys.stderr.write('conflicts: %d\n' % n_conflicts)


if __name__ == '__main__':
    main()
