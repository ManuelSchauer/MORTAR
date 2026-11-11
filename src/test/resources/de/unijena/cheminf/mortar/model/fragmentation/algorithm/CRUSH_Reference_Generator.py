"""
MORTAR - MOlecule fRagmenTAtion fRamework, test data generator (MIT licence, see the MORTAR repository).

Reference values for the MORTAR class CRUSHFragmenter, an implementation of CRUSH (Cleavage Rules Using SMIRKS
Heuristics, E. Lopez-Lopez, J. L. Medina-Franco, F. Miljkovic, J. Cheminform. 18, 78 (2026),
https://doi.org/10.1186/s13321-026-01239-w). The 33 rules are read as data from CRUSH_SMIRKS_Rules.txt, a verbatim copy
of the rule file of the authors (https://github.com/EdgL2/CRUSH, MIT licence, see the header of that file). This script
is independent code that uses the RDKit; it does not contain or run code of the authors. It has three parts:

1. original_steps: the single cleavage steps as the authors' published procedure defines them (reaction SMIRKS applied
   with RunReactants, products written as canonical SMILES and parsed again, fragments with fewer than 3 heavy atoms or
   unparsable products rejected, atoms in at least two rings and their neighbours protected, rules with a match on a
   protected atom skipped as a whole). Every product set is classified: V valid, S rejected by size, P rule skipped
   because of protected atoms, R ring bond cleaved by the reaction (overlapping products with duplicated atoms), U
   unparsable product (e.g. cleavage at a neutral aromatic nitrogen atom).
2. mortar_steps: the cleavage steps as CRUSHFragmenter defines them, written independently of part 1: the reactant side
   of each SMIRKS is matched as SMARTS, the bonds between mapped atoms that are missing on the product side are cleaved
   (never ring bonds), and the open valences are saturated with hydrogen atoms. The generator checks itself: every V
   step of part 1 must be a step of part 2 with the same products, and every step of part 2 that is not V in part 1
   must be a U case of part 1. Valid steps of part 1 whose products differ from part 2 only in hydrogen counts are
   category D: RunReactants recomputes the hydrogen count of mapped atoms from the default valences of the RDKit,
   which saturates atoms with radical electrons of the input and changes the valence of hypervalent atoms with hydrogen
   atoms (e.g. O=[PH](O)C gives O=PO), whereas CRUSHFragmenter adds one hydrogen atom per order of a cleaved bond.
3. mortar_route: the exact search for the decomposition with the maximum number of final fragments (memoised per atom
   set, at most a given number of evaluated fragments), with the list of all optimal final fragment multisets.

Usage: python -I CRUSH_Reference_Generator.py steps <input> <output .tsv> [options]
       python -I CRUSH_Reference_Generator.py routes <input> <output .tsv> [options]
  input                  .sdf, .smi (SMILES and name per line), or .tsv (rule cases: name, SMILES, rule, expectation,
                         description; lines starting with '#' are ignored)
  --every N              use only every N-th molecule of the input
  --max N                stop after N molecules
  --intermediates N      (steps) add up to N intermediate fragments of the route search per molecule (default 0)
  --budget N             maximum number of evaluated fragments per molecule in the route search (default 3000)
  --rules PATH           rule file (default: CRUSH_SMIRKS_Rules.txt next to this script)
  --comment TEXT         add a comment line to the header of the output, e.g. source and licence of the input data
  --allow-rdkit-version  run with another RDKit version than REFERENCE_RDKIT_VERSION

Stereochemistry: the default stereo perception of the RDKit is used. It removes the configuration of atoms that are not
stereogenic because of the configuration of other atoms (e.g. C3 of the chiral 1,5-anhydroarabinitol formed from
arabinopyranosides by rules 6 and 7), which the CDK classifies as potential (para) stereocentres and keeps;
CRUSHFragmenterReferenceTest reports such differences separately.

The reference files in the MORTAR test resources were generated with RDKit 2026.09.1 (Python 3.12) by
  python -I CRUSH_Reference_Generator.py steps CRUSH_Rule_Cases_Input.tsv CRUSH_Rule_Cases.tsv --comment "..."
  python -I CRUSH_Reference_Generator.py steps Tutorial/biofacquim-06-2025.sdf CRUSH_Reference_Steps_BIOFACQUIM.tsv
      --intermediates 20 --comment "..."
  python -I CRUSH_Reference_Generator.py routes Tutorial/biofacquim-06-2025.sdf CRUSH_Reference_Fragments_BIOFACQUIM.tsv
      --comment "..."
"""
import ast
import hashlib
import io
import os
import sys
from collections import Counter

from rdkit import Chem
from rdkit import RDLogger
from rdkit import rdBase
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')

REFERENCE_RDKIT_VERSION = '2026.09.1'
RULE_FILE_MD5 = '9584dbe50ad42934d71c0656c77e53b8'
MINIMUM_FRAGMENT_SIZE = 3
DEFAULT_BUDGET = 3000
MAXIMUM_OPTIMAL_SETS = 50
MAXIMUM_OPTIMAL_CUT_SETS = 2000


class Rule:
    """a CRUSH rule: number (order of Table 1 of the publication), name, SMIRKS, compiled reaction, reactant query,
    and the cleaved bonds as pairs of map numbers and of query atom indices"""

    def __init__(self, number, name, smirks):
        self.number = number
        self.name = name
        self.smirks = smirks
        reactant, product = smirks.split('>>')
        self.smarts = reactant
        self.reaction = AllChem.ReactionFromSmarts(smirks)
        self.query = Chem.MolFromSmarts(reactant)
        if self.query is None or self.reaction is None:
            raise ValueError('cannot parse rule %d' % number)
        product_query = Chem.MolFromSmarts(product)
        map_to_index = {a.GetAtomMapNum(): a.GetIdx() for a in self.query.GetAtoms() if a.GetAtomMapNum()}

        def mapped_bonds(mol):
            return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())))
                    for b in mol.GetBonds() if b.GetBeginAtom().GetAtomMapNum() and b.GetEndAtom().GetAtomMapNum()}
        self.cleaved_map_pairs = sorted(mapped_bonds(self.query) - mapped_bonds(product_query))
        self.cleaved_query_pairs = [(map_to_index[a], map_to_index[b]) for a, b in self.cleaved_map_pairs]


def read_rules(path):
    """reads the rules from the verbatim copy of the authors' rule file (header lines starting with '#' are skipped,
    the rest must be identical to the original file; the original has CRLF line endings, Git may store the copy with
    LF line endings, so the line endings are converted back before the check)"""
    with open(path, 'rb') as handle:
        data = handle.read()
    lines = data.replace(b'\r\n', b'\n').split(b'\n')
    body = b'\r\n'.join(lines[next(i for i, l in enumerate(lines) if not l.startswith(b'#')):])
    if hashlib.md5(body).hexdigest() != RULE_FILE_MD5:
        raise ValueError('the rule file is not the verbatim copy of CRUSH_SMIRKS_Rules.txt')
    text = body.decode('utf-8')
    pairs = ast.literal_eval(text[text.index('=') + 1:].strip())
    return [Rule(i + 1, name, smirks) for i, (name, smirks) in enumerate(pairs)]


def protected_atoms(mol):
    """atoms in at least two rings of the ring information of the RDKit (symmetrised SSSR), and their neighbours"""
    counts = [0] * mol.GetNumAtoms()
    for ring in mol.GetRingInfo().AtomRings():
        for index in ring:
            counts[index] += 1
    fused = {i for i, count in enumerate(counts) if count >= 2}
    neighbours = {n.GetIdx() for i in fused for n in mol.GetAtomWithIdx(i).GetNeighbors()}
    return fused | neighbours


def heavy_atom_count(mol):
    return mol.GetNumHeavyAtoms() if mol is not None else 0


def format_bonds(pairs):
    return ' '.join('%d-%d' % pair for pair in sorted(pairs))


# ----------------------------------------------------------------------------------------------------------------------
# part 1: single steps as defined by the published procedure (RunReactants semantics)
# ----------------------------------------------------------------------------------------------------------------------

def original_steps(mol, rules, minimum_size=MINIMUM_FRAGMENT_SIZE):
    """returns a list of (rule number, category, cleaved bonds, product SMILES) of the published procedure for the
    given (sanitised) molecule, each distinct entry once; category P has no bonds and products"""
    protected = protected_atoms(mol)
    entries = []
    seen = set()
    for rule in rules:
        template = rule.reaction.GetReactantTemplate(0)
        matches = mol.GetSubstructMatches(template, useChirality=False)
        if not matches:
            continue
        if any(index in protected for match in matches for index in match):
            entries.append((rule.number, 'P', (), ()))
            continue
        for products in rule.reaction.RunReactants((mol,)):
            map_to_atom = {}
            for product in products:
                for atom in product.GetAtoms():
                    if atom.HasProp('old_mapno') and atom.HasProp('react_atom_idx'):
                        map_to_atom[atom.GetIntProp('old_mapno')] = atom.GetIntProp('react_atom_idx')
            bonds = tuple(sorted(tuple(sorted((map_to_atom[a], map_to_atom[b]))) for a, b in rule.cleaved_map_pairs))
            smiles = [Chem.MolToSmiles(product, True) for product in products]
            parsed = [Chem.MolFromSmiles(s) for s in smiles]
            if any(mol.GetBondBetweenAtoms(a, b).IsInRing() for a, b in bonds):
                category = 'R'
            elif any(p is None for p in parsed):
                category = 'U'
            elif any(heavy_atom_count(p) < minimum_size for p in parsed):
                category = 'S'
            else:
                category = 'V'
            key = (rule.number, category, bonds, tuple(sorted(smiles)))
            if key not in seen:
                seen.add(key)
                entries.append(key)
    return entries


# ----------------------------------------------------------------------------------------------------------------------
# part 2: single steps as defined by CRUSHFragmenter (bond cleavage, hydrogen saturation)
# ----------------------------------------------------------------------------------------------------------------------

def cleave(mol, pairs):
    """cleaves the given bonds (atom index pairs) of the molecule and saturates the open valences with hydrogen atoms;
    stereo elements are kept, a cleaved ligand of a stereocentre or stereo double bond is replaced by the hydrogen atom.
    Returns the new molecule; the indices of the atoms of the molecule are kept."""
    frozen = Chem.RWMol(mol)
    for atom in frozen.GetAtoms():
        hydrogens = atom.GetTotalNumHs()
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(hydrogens)
    frozen = frozen.GetMol()
    if not pairs:
        return frozen
    bond_indices = [mol.GetBondBetweenAtoms(a, b).GetIdx() for a, b in pairs]
    fragmented = Chem.RWMol(Chem.FragmentOnBonds(frozen, bond_indices, addDummies=True))
    removed = []
    for atom in fragmented.GetAtoms():
        if atom.GetIdx() < mol.GetNumAtoms():
            continue
        bond = atom.GetBonds()[0]
        neighbour = bond.GetOtherAtom(atom)
        if bond.GetBondType() == Chem.BondType.SINGLE:
            # an explicit hydrogen atom at the position of the former ligand, merged by RemoveHs
            atom.SetAtomicNum(1)
            atom.SetIsotope(0)
        else:
            neighbour.SetNumExplicitHs(neighbour.GetNumExplicitHs() + int(bond.GetBondTypeAsDouble()))
            removed.append(atom.GetIdx())
    for index in sorted(removed, reverse=True):
        fragmented.RemoveAtom(index)
    result = fragmented.GetMol()
    Chem.SanitizeMol(result)
    parameters = Chem.RemoveHsParameters()
    # hydrogen atoms that define the configuration of a double bond are removed as well, the configuration is then
    # defined with respect to the other substituent
    parameters.removeDefiningBondStereo = True
    result = Chem.RemoveHs(result, parameters)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return result


def fragments_with_atoms(mol, atom_sets):
    """returns the canonical isomeric SMILES of the connected parts of the molecule that contain the given atom sets"""
    mapping = []
    parts = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping)
    result = []
    for atoms in atom_sets:
        for part, indices in zip(parts, mapping):
            if next(iter(atoms)) in indices:
                result.append(Chem.MolToSmiles(part))
                break
    return result


def components_of(mol, atoms, cut_pairs):
    """connected parts of the given atom set after removing the given bonds"""
    cut = {frozenset(pair) for pair in cut_pairs}
    remaining = set(atoms)
    parts = []
    while remaining:
        start = remaining.pop()
        part = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbour in mol.GetAtomWithIdx(current).GetNeighbors():
                index = neighbour.GetIdx()
                if index in remaining and frozenset((current, index)) not in cut:
                    remaining.discard(index)
                    part.add(index)
                    stack.append(index)
        parts.append(frozenset(part))
    return parts


def mortar_step_candidates(mol, rules, minimum_size=MINIMUM_FRAGMENT_SIZE):
    """returns {(rule number, cleaved bonds): parts} of the valid steps of CRUSHFragmenter for the given (sanitised,
    hydrogen-saturated) molecule; every connected component is a fragment of its own"""
    protected = protected_atoms(mol)
    component_of = {}
    for number, component in enumerate(Chem.GetMolFrags(mol)):
        for index in component:
            component_of[index] = (number, frozenset(component))
    steps = {}
    for rule in rules:
        matches = mol.GetSubstructMatches(rule.query, uniquify=False, useChirality=False, maxMatches=10000000)
        skipped_components = {component_of[match[0]][0] for match in matches
                              if any(index in protected for index in match)}
        for match in matches:
            if component_of[match[0]][0] in skipped_components:
                continue
            pairs = tuple(sorted(tuple(sorted((match[a], match[b]))) for a, b in rule.cleaved_query_pairs))
            if any(mol.GetBondBetweenAtoms(a, b).IsInRing() for a, b in pairs):
                continue
            key = (rule.number, pairs)
            if key in steps:
                continue
            parts = components_of(mol, component_of[match[0]][1], pairs)
            steps[key] = parts if all(sum(1 for i in part if mol.GetAtomWithIdx(i).GetAtomicNum() > 1) >= minimum_size
                                      for part in parts) else None
    return {key: parts for key, parts in steps.items() if parts is not None}


def mortar_steps(mol, rules, minimum_size=MINIMUM_FRAGMENT_SIZE):
    """returns a list of (rule number, cleaved bonds, product SMILES) of the valid steps of CRUSHFragmenter"""
    result = []
    for (number, pairs), parts in sorted(mortar_step_candidates(mol, rules, minimum_size).items()):
        cleaved = cleave(mol, pairs)
        result.append((number, pairs, tuple(sorted(fragments_with_atoms(cleaved, parts)))))
    return result


# ----------------------------------------------------------------------------------------------------------------------
# part 3: exact search for the decomposition with the maximum number of final fragments
# ----------------------------------------------------------------------------------------------------------------------

class BudgetExceeded(Exception):
    pass


class RouteSearch:
    """exact search on the atom sets of the given (sanitised) molecule"""

    def __init__(self, mol, rules, budget=DEFAULT_BUDGET, minimum_size=MINIMUM_FRAGMENT_SIZE):
        self.mol = mol
        self.rules = rules
        self.budget = budget
        self.minimum_size = minimum_size
        self.steps_cache = {}
        self.best_cache = {}
        self.cut_sets_cache = {}
        self.evaluations = 0
        self.fragment_smiles = {}

    def fragment(self, atoms):
        """the hydrogen-saturated fragment with the given atoms and the map of its atom indices to the molecule"""
        boundary = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol.GetBonds()
                    if (b.GetBeginAtomIdx() in atoms) != (b.GetEndAtomIdx() in atoms)]
        cleaved = cleave(self.mol, boundary)
        mapping = []
        parts = Chem.GetMolFrags(cleaved, asMols=True, sanitizeFrags=True, fragsMolAtomMapping=mapping)
        start = next(iter(atoms))
        for part, indices in zip(parts, mapping):
            if start in indices:
                return part, list(indices)
        raise ValueError('fragment not found')

    def steps(self, atoms):
        """unique steps of the fragment (lowest rule number per set of cleaved bonds), sorted by rule number, as
        (rule number, cleaved bonds in molecule indices, parts as atom sets of the molecule)"""
        if atoms in self.steps_cache:
            return self.steps_cache[atoms]
        self.evaluations += 1
        if self.evaluations > self.budget:
            raise BudgetExceeded()
        fragment, mapping = self.fragment(atoms)
        self.fragment_smiles[atoms] = Chem.MolToSmiles(fragment)
        unique = {}
        for (number, pairs), parts in mortar_step_candidates(fragment, self.rules, self.minimum_size).items():
            global_pairs = tuple(sorted(tuple(sorted((mapping[a], mapping[b]))) for a, b in pairs))
            if global_pairs not in unique or unique[global_pairs][0] > number:
                unique[global_pairs] = (number, global_pairs,
                                        tuple(frozenset(mapping[i] for i in part) for part in parts))
        result = sorted(unique.values(), key=lambda step: (step[0], step[1]))
        self.steps_cache[atoms] = result
        return result

    def best(self, atoms):
        if atoms in self.best_cache:
            return self.best_cache[atoms]
        steps = self.steps(atoms)
        value = 1 if not steps else max(sum(self.best(part) for part in step[2]) for step in steps)
        self.best_cache[atoms] = value
        return value

    def optimal_cut_sets(self, atoms):
        """all sets of cleaved bonds of optimal decompositions of the fragment, None if there are too many"""
        if atoms in self.cut_sets_cache:
            return self.cut_sets_cache[atoms]
        steps = self.steps(atoms)
        if not steps:
            result = {frozenset()}
        else:
            target = self.best(atoms)
            result = set()
            for number, pairs, parts in steps:
                if sum(self.best(part) for part in parts) != target:
                    continue
                combinations = {frozenset(pairs)}
                for part in parts:
                    part_sets = self.optimal_cut_sets(part)
                    if part_sets is None:
                        result = None
                        break
                    combinations = {c | p for c in combinations for p in part_sets}
                    if len(combinations) > MAXIMUM_OPTIMAL_CUT_SETS:
                        result = None
                        break
                if result is None:
                    break
                result |= combinations
                if len(result) > MAXIMUM_OPTIMAL_CUT_SETS:
                    result = None
                    break
        self.cut_sets_cache[atoms] = result
        return result


def mortar_route(mol, rules, budget=DEFAULT_BUDGET):
    """returns (maximum number of final fragments or None, exhaustive, evaluations, optimal fragment multisets or None
    if there are more than MAXIMUM_OPTIMAL_SETS, search) for the (sanitised) molecule"""
    search = RouteSearch(mol, rules, budget)
    components = [frozenset(c) for c in Chem.GetMolFrags(mol)]
    try:
        maximum = sum(search.best(component) for component in components)
        cut_sets = {frozenset()}
        for component in components:
            component_sets = search.optimal_cut_sets(component)
            if component_sets is None:
                cut_sets = None
                break
            cut_sets = {c | p for c in cut_sets for p in component_sets}
            if len(cut_sets) > MAXIMUM_OPTIMAL_CUT_SETS:
                cut_sets = None
                break
    except BudgetExceeded:
        return None, False, search.evaluations, None, search
    if cut_sets is None:
        return maximum, True, search.evaluations, None, search
    multisets = set()
    for cut_set in cut_sets:
        cleaved = cleave(mol, sorted(cut_set))
        multisets.add(tuple(sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(cleaved, asMols=True,
                                                                                  sanitizeFrags=False))))
    if len(multisets) > MAXIMUM_OPTIMAL_SETS:
        return maximum, True, search.evaluations, None, search
    for multiset in multisets:
        if len(multiset) != maximum:
            raise AssertionError('wrong number of fragments %s' % (multiset,))
    return maximum, True, search.evaluations, sorted(multisets), search


# ----------------------------------------------------------------------------------------------------------------------
# input and output
# ----------------------------------------------------------------------------------------------------------------------

def sdf_records(path, every):
    """yields (record index, mol block) of every N-th record of an SD file, without parsing the other records"""
    index = 0
    lines = []
    with open(path, encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('$$$$'):
                if index % every == 0:
                    yield index, ''.join(lines)
                index += 1
                lines = []
            elif index % every == 0:
                lines.append(line)


def molecules(path, every, maximum):
    """yields (name, SMILES, case information) of the input molecules"""
    count = 0
    if path.endswith('.sdf'):
        for index, block in sdf_records(path, every):
            mol = next(Chem.ForwardSDMolSupplier(io.BytesIO((block + '$$$$\n').encode('utf-8'))), None)
            if mol is None:
                continue
            name = ''
            for key in ('_Name', 'identifier', 'chembl_id'):
                if not name and mol.HasProp(key):
                    name = mol.GetProp(key).strip()
            yield name or str(index), Chem.MolToSmiles(mol), None
            count += 1
            if maximum and count >= maximum:
                return
    else:
        with open(path, encoding='utf-8') as handle:
            for index, line in enumerate(handle):
                if index % every or not line.strip() or line.startswith('#'):
                    continue
                if path.endswith('.tsv'):
                    columns = line.rstrip('\n').split('\t')
                    yield columns[0], columns[1], (int(columns[2]), columns[3], columns[4])
                else:
                    parts = line.split()
                    yield (parts[1] if len(parts) > 1 else str(index)), parts[0], None
                count += 1
                if maximum and count >= maximum:
                    return


def without_hydrogen_counts(smiles_list):
    """canonical SMILES codes (without stereochemistry) of the given SMILES codes in which no atom has hydrogen atoms
    or radical electrons, i.e. the heavy-atom skeletons with charges and bond orders"""
    result = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result.append(smiles)
            continue
        Chem.RemoveStereochemistry(mol)
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(0)
        mol.UpdatePropertyCache(strict=False)
        result.append(Chem.MolToSmiles(mol))
    return sorted(result)


def check_steps(name, mol, original, mortar, counts, problems):
    """self-check of parts 1 and 2, see module documentation; counts the categories and returns the entries of part 1,
    in which V entries whose products differ from part 2 only in hydrogen counts and radical electrons (recomputed by
    the reaction from the default valences) are reclassified as D"""
    mortar_keys = {(number, pairs): products for number, pairs, products in mortar}
    valid_keys = set()
    unparsable_keys = set()
    skipped_rules = set()
    result = []
    for entry in original:
        number, category, pairs, products = entry
        if category == 'V':
            valid_keys.add((number, pairs))
            expected = mortar_keys.get((number, pairs))
            if expected is None:
                problems.append('%s: V step of rule %d (%s) is not a step of part 2' % (name, number, format_bonds(pairs)))
            elif canonical_without_stereo(expected) != canonical_without_stereo(products):
                # the reactions of the RDKit recompute the hydrogen count of mapped atoms from the default valences,
                # which removes radical electrons of the input and changes the valence of hypervalent atoms with
                # hydrogen atoms (e.g. O=[PH](O)C gives O=PO instead of O=[PH2]O); such products differ from part 2
                # only in the hydrogen counts
                if without_hydrogen_counts(expected) == without_hydrogen_counts(products):
                    category = 'D'
                    entry = (number, category, pairs, products)
                else:
                    problems.append('%s: rule %d (%s): products %s (part 1) and %s (part 2) differ'
                                    % (name, number, format_bonds(pairs), products, expected))
            elif sorted(expected) != sorted(products):
                counts['V with different stereo in part 1'] += 1
        counts[category] += 1
        result.append(entry)
        if category == 'U':
            unparsable_keys.add((number, pairs))
        elif category == 'P':
            skipped_rules.add(number)
    for number, pairs, products in mortar:
        counts['part 2'] += 1
        if (number, pairs) in valid_keys:
            continue
        if (number, pairs) in unparsable_keys:
            counts['part 2 step that is U in part 1'] += 1
        elif number in skipped_rules:
            # only possible for molecules with several components (protection checked per component)
            counts['part 2 step of a rule skipped in part 1 (several components)'] += 1
        else:
            problems.append('%s: step of rule %d (%s) of part 2 is neither V nor U in part 1'
                            % (name, number, format_bonds(pairs)))
    return result


def canonical_without_stereo(smiles_list):
    result = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result.append(smiles)
            continue
        Chem.RemoveStereochemistry(mol)
        result.append(Chem.MolToSmiles(mol))
    return sorted(result)


def write_steps(out, name, smiles, mol, rules, case, counts, problems):
    """writes the molecule row and the step rows of one molecule"""
    protected = sorted(protected_atoms(mol))
    if case is None:
        target, expectation, description = 0, '-', ''
    else:
        target, expectation, description = case
    out.write('\t'.join([name, smiles, 'M', str(target), expectation, ' '.join(map(str, protected)), description])
              + '\n')
    original = original_steps(mol, rules)
    mortar = mortar_steps(mol, rules)
    original = check_steps(name, mol, original, mortar, counts, problems)
    if case is not None and target and expectation in ('positive', 'negative'):
        found = any(number == target for number, pairs, products in mortar)
        if found != (expectation == 'positive'):
            problems.append('%s: expectation %s for rule %d not met' % (name, expectation, target))
    for number, pairs, products in mortar:
        out.write('\t'.join([name, '', '2', str(number), 'V', format_bonds(pairs), '.'.join(products)]) + '\n')
    mortar_products = {(number, pairs): products for number, pairs, products in mortar}
    for number, category, pairs, products in original:
        # the products of part 1 are written as '=' if they are identical to those of part 2, and not at all for
        # categories P, S and R (overlapping products with duplicated atoms)
        if category in ('P', 'S', 'R'):
            text = ''
        elif category == 'V' and tuple(sorted(products)) == mortar_products.get((number, pairs)):
            text = '='
        else:
            text = '.'.join(products)
        out.write('\t'.join([name, '', '1', str(number), category, format_bonds(pairs), text]) + '\n')


def intermediate_fragments(mol, rules, budget, maximum):
    """up to the given number of distinct hydrogen-saturated fragments evaluated in the route search (without the
    connected components of the molecule), the largest first"""
    search = RouteSearch(mol, rules, budget)
    components = [frozenset(c) for c in Chem.GetMolFrags(mol)]
    try:
        for component in components:
            search.best(component)
    except BudgetExceeded:
        pass
    smiles = {}
    for atoms, fragment_smiles in search.fragment_smiles.items():
        if atoms not in components:
            smiles[fragment_smiles] = len(atoms)
    ordered = sorted(smiles, key=lambda s: (-smiles[s], s))
    return ordered[:maximum]


def header(out, mode, args, rules, extra):
    out.write('# CRUSH reference %s (CRUSHFragmenter of MORTAR, 33 rules of CRUSH_SMIRKS_Rules.txt, minimum fragment size'
              ' %d, fused ring atoms and their neighbours protected, conservative mode) generated by'
              ' CRUSH_Reference_Generator.py with RDKit %s from %s\n'
              % (mode, MINIMUM_FRAGMENT_SIZE, rdBase.rdkitVersion, os.path.basename(args[1])))
    if '--comment' in args:
        out.write('# %s\n' % args[args.index('--comment') + 1])
    for line in extra:
        out.write('# %s\n' % line)


def main():
    args = sys.argv[1:]
    if len(args) < 3 or args[0] not in ('steps', 'routes'):
        print(__doc__)
        sys.exit(1)
    if rdBase.rdkitVersion != REFERENCE_RDKIT_VERSION and '--allow-rdkit-version' not in args:
        sys.exit('RDKit %s is installed, the reference values are generated with RDKit %s (use --allow-rdkit-version)'
                 % (rdBase.rdkitVersion, REFERENCE_RDKIT_VERSION))
    every = int(args[args.index('--every') + 1]) if '--every' in args else 1
    maximum = int(args[args.index('--max') + 1]) if '--max' in args else 0
    budget = int(args[args.index('--budget') + 1]) if '--budget' in args else DEFAULT_BUDGET
    intermediates = int(args[args.index('--intermediates') + 1]) if '--intermediates' in args else 0
    rule_path = args[args.index('--rules') + 1] if '--rules' in args \
        else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRUSH_SMIRKS_Rules.txt')
    rules = read_rules(rule_path)
    rows = []
    counts = Counter()
    problems = []
    molecule_count = 0
    if args[0] == 'steps':
        body = io.StringIO()
        for name, smiles, case in molecules(args[1], every, maximum):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                if case is not None:
                    problems.append('%s: rule case %s cannot be parsed' % (name, smiles))
                continue
            # the canonical SMILES of the RDKit is the common input of both implementations
            smiles = Chem.MolToSmiles(mol)
            mol = Chem.MolFromSmiles(smiles)
            molecule_count += 1
            write_steps(body, name, smiles, mol, rules, case, counts, problems)
            if intermediates:
                for k, fragment_smiles in enumerate(intermediate_fragments(mol, rules, budget, intermediates)):
                    fragment = Chem.MolFromSmiles(fragment_smiles)
                    molecule_count += 1
                    write_steps(body, '%s/%d' % (name, k + 1), fragment_smiles, fragment, rules, None, counts, problems)
        extra = ['molecules (including intermediate fragments): %d' % molecule_count,
                 'part 1 categories: ' + ', '.join('%s %d' % (c, counts[c]) for c in 'VSPRUD'),
                 'part 1 V steps whose products differ from part 2 in stereochemistry only: %d'
                 % counts['V with different stereo in part 1'],
                 'part 2 steps: %d, of which U in part 1: %d, of rules skipped in part 1 (several components): %d'
                 % (counts['part 2'], counts['part 2 step that is U in part 1'],
                    counts['part 2 step of a rule skipped in part 1 (several components)']),
                 'self-check problems: %d' % len(problems),
                 'record M: name, input SMILES, M, rule (rule cases), expectation, protected atoms, description',
                 'record 2: name, -, 2, rule, V, cleaved bonds (atom indices), products of CRUSHFragmenter (hydrogen)',
                 'record 1: name, -, 1, rule, category (V valid, S size, P protected, R ring bond, U unparsable, D valid'
                 ' with products that differ from part 2 only in hydrogen counts, since the reaction recomputes them from'
                 ' the default valences, e.g. for atoms with radical electrons or hypervalent atoms with hydrogen atoms),'
                 ' cleaved bonds, products of the published procedure (=: identical to part 2; not given for P, S, R)']
        with open(args[2], 'w', encoding='utf-8') as out:
            header(out, 'steps', args, rules, extra)
            out.write('# name\tSMILES\trecord\trule\tcategory\tatoms\tproducts\n')
            out.write(body.getvalue())
    else:
        import time
        body = io.StringIO()
        exhaustive_count = 0
        start = time.time()
        for name, smiles, case in molecules(args[1], every, maximum):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol)
            mol = Chem.MolFromSmiles(smiles)
            molecule_count += 1
            maximum_count, exhaustive, evaluations, multisets, search = mortar_route(mol, rules, budget)
            exhaustive_count += exhaustive
            body.write('\t'.join([name, smiles, str(maximum_count) if exhaustive else '-',
                                  'exhaustive' if exhaustive else 'budget exceeded', str(evaluations),
                                  str(len(multisets)) if multisets is not None else '-',
                                  ' '.join('.'.join(m) for m in multisets) if multisets is not None else '']) + '\n')
        print('%.0f s' % (time.time() - start))
        extra = ['molecules: %d, exhaustive searches (at most %d evaluated fragments): %d'
                 % (molecule_count, budget, exhaustive_count),
                 'optimal fragment sets are listed if there are at most %d' % MAXIMUM_OPTIMAL_SETS]
        with open(args[2], 'w', encoding='utf-8') as out:
            header(out, 'routes', args, rules, extra)
            out.write('# name\tSMILES\tmaximum number of final fragments\tsearch\tevaluated fragments\t'
                      'number of optimal fragment sets\toptimal fragment sets (hydrogen saturation, sets separated by'
                      ' spaces)\n')
            out.write(body.getvalue())
    for problem in problems:
        print(problem)
    print('%d molecules, %s, %d problems' % (molecule_count, dict(counts), len(problems)))
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
