"""
Generates StereoSmilesReference.tsv, the reference data of ChemUtilStereoSmilesTest.

For every hand-picked skeleton below, all combinations of the tetrahedral configurations marked in the template
(every '@' is replaced by '@' or '@@') and of the double bond configurations (every '/' or '\\' marked as '{/}' is
replaced by '/' or '\\') are written out as input SMILES. Each input SMILES is additionally rewritten by RDKit with
random atom orders. All writings are grouped into compounds by their RDKit canonical isomeric SMILES; the InChIKey is
given for information where RDKit can create one (not for structures with pseudo atoms '*').

The canonical SMILES are created with the legacy stereo perception of RDKit (the default of RDKit 2026.09.1). The
new stereo perception (Chem.SetUseLegacyStereoPerception(False)) is not used, because it keeps the configuration of
the non-stereogenic middle carbon atom of the chiral 1,5-anhydropentitols, so that it splits one compound into two
canonical SMILES (e.g. O[C@H]1[C@@H](O)COC[C@@H]1O and O[C@@H]1[C@@H](O)COC[C@@H]1O, both with the InChIKey
QXAMTEJJAZOINB-IMJSIDKUSA-N). For all structures without pseudo atoms, the script checks that the grouping by the
canonical SMILES and the grouping by the InChIKey agree and fails otherwise.

Usage (RDKit 2026.09.1):
    python -I generate_stereo_smiles_reference.py > StereoSmilesReference.tsv
"""

import itertools
import random
import re
import sys

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

# (name, template); '[C@' marks a tetrahedral centre whose configuration is enumerated, '{/}' a double bond marker
SKELETONS = [
    # 1,5-anhydropentitols: C3 is not stereogenic in the chiral compounds and pseudo-asymmetric in the meso compounds
    ("1,5-anhydropentitol", "O[C@H]1COC[C@H](O)[C@H]1O"),
    ("1,5-anhydropentitol, link atom on O3", "O[C@H]1COC[C@H](O)[C@H]1O*"),
    ("1,5-anhydropentitol, link atoms on O2 and O4", "*O[C@H]1COC[C@H](O*)[C@H]1O"),
    ("1,5-anhydropentitol, link atom on C3", "O[C@H]1COC[C@H](O)[C@]1(O)*"),
    # pentitols (ribitol, xylitol, arabinitol)
    ("pentitol", "OC[C@H](O)[C@H](O)[C@H](O)CO"),
    ("pentitol, link atoms on the terminal O", "*OC[C@H](O)[C@H](O)[C@H](O)CO*"),
    # 1,4- and 1,3-disubstituted cyclohexanes
    ("1,4-dimethylcyclohexane", "C[C@H]1CC[C@H](C)CC1"),
    ("cyclohexane-1,4-diol", "O[C@H]1CC[C@H](O)CC1"),
    ("4-methylcyclohexanol", "C[C@H]1CC[C@H](O)CC1"),
    ("cyclohexane-1,4-diyl link atoms", "*[C@H]1CC[C@H](*)CC1"),
    ("4-hydroxycyclohexyl link atom", "*[C@H]1CC[C@H](O)CC1"),
    ("1,3-dimethylcyclohexane", "C[C@H]1CCC[C@H](C)C1"),
    # inositols (9 stereoisomers)
    ("inositol", "O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),
    ("inositol, link atom on O1", "*O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),
    # tartaric acids
    ("tartaric acid", "OC(=O)[C@H](O)[C@H](O)C(=O)O"),
    ("tartaric acid, link atoms on the acid groups", "*C(=O)[C@H](O)[C@H](O)C(=O)*"),
    # 2,3,4-trihydroxyglutaric acids: C3 pseudo-asymmetric in the meso compounds, not stereogenic in the chiral ones
    ("2,3,4-trihydroxyglutaric acid", "OC(=O)[C@H](O)[C@H](O)[C@H](O)C(=O)O"),
    ("2,3,4-trihydroxyglutaric acid, link atoms on the acid groups", "*OC(=O)[C@H](O)[C@H](O)[C@H](O)C(=O)O*"),
    # atoms with two equivalent ligands that are never stereogenic
    ("two link atoms on one carbon atom", "*[C@H](*)O"),
    ("two link atoms on one carbon atom with a stereocentre", "*[C@H](*)[C@H](O)C"),
    ("propane-1,2,3-triol", "OC[C@H](O)CO"),
    # normal stereocentres
    ("alanine", "N[C@H](C)C(=O)O"),
    ("butan-2-ol with link atom", "*[C@H](O)CC"),
    ("glucopyranose", "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),
    ("menthol", "CC(C)[C@H]1CC[C@H](C)C[C@H]1O"),
    # double bonds
    ("but-2-ene", "C{/}C=C{/}C"),
    ("2-methylbut-2-ene", "C{/}C(C)=C{/}C"),
    ("hex-3-ene-2,5-diol", "C[C@H](O){/}C=C{/}[C@H](C)O"),
]

RANDOM_WRITINGS_PER_INPUT = 2
SEED = 42


def enumerate_writings(aTemplate):
    """Returns all SMILES obtained by setting each marked tetrahedral and double bond configuration."""
    tmpParts = re.split(r"(\[C@\]|\[C@H\]|\{/\})", aTemplate)
    tmpOptions = []
    for tmpPart in tmpParts:
        if tmpPart == "[C@]":
            tmpOptions.append(["[C@]", "[C@@]"])
        elif tmpPart == "[C@H]":
            tmpOptions.append(["[C@H]", "[C@@H]"])
        elif tmpPart == "{/}":
            tmpOptions.append(["/", "\\"])
        else:
            tmpOptions.append([tmpPart])
    return ["".join(tmpCombination) for tmpCombination in itertools.product(*tmpOptions)]


def random_writings(aSmiles, aCount, aRandom):
    """Returns RDKit writings of the given SMILES with random atom orders."""
    tmpMol = Chem.MolFromSmiles(aSmiles)
    tmpWritings = []
    for _ in range(aCount):
        tmpOrder = list(range(tmpMol.GetNumAtoms()))
        aRandom.shuffle(tmpOrder)
        tmpRenumbered = Chem.RenumberAtoms(tmpMol, tmpOrder)
        tmpWritings.append(Chem.MolToSmiles(tmpRenumbered, canonical=False))
    return tmpWritings


def canonical_identity(aSmiles):
    """Returns the RDKit canonical isomeric SMILES and the InChIKey (empty for structures with pseudo atoms)."""
    tmpMol = Chem.MolFromSmiles(aSmiles)
    tmpInChIKey = "" if "*" in aSmiles else Chem.MolToInchiKey(tmpMol)
    return Chem.MolToSmiles(tmpMol), tmpInChIKey


def main():
    Chem.SetUseLegacyStereoPerception(True)
    tmpRandom = random.Random(SEED)
    tmpRows = []
    for tmpName, tmpTemplate in SKELETONS:
        for tmpWriting in enumerate_writings(tmpTemplate):
            for tmpSmiles in [tmpWriting] + random_writings(tmpWriting, RANDOM_WRITINGS_PER_INPUT, tmpRandom):
                tmpCanonicalSmiles, tmpInChIKey = canonical_identity(tmpSmiles)
                tmpRows.append((tmpName, tmpSmiles, tmpCanonicalSmiles, tmpInChIKey))
    # canonical SMILES and InChIKey must group the writings in the same way
    tmpSmilesToKey = {}
    tmpKeyToSmiles = {}
    for _, tmpSmiles, tmpCanonicalSmiles, tmpInChIKey in tmpRows:
        if tmpInChIKey == "":
            continue
        if (tmpSmilesToKey.setdefault(tmpCanonicalSmiles, tmpInChIKey) != tmpInChIKey
                or tmpKeyToSmiles.setdefault(tmpInChIKey, tmpCanonicalSmiles) != tmpCanonicalSmiles):
            sys.exit("Canonical SMILES and InChIKey disagree for " + tmpSmiles)
    print("# generated by generate_stereo_smiles_reference.py with RDKit " + rdBase.rdkitVersion)
    print("name\tinput_smiles\trdkit_canonical_smiles\tinchikey")
    tmpSeen = set()
    for tmpName, tmpSmiles, tmpCanonicalSmiles, tmpInChIKey in tmpRows:
        if (tmpName, tmpSmiles) in tmpSeen:
            continue
        tmpSeen.add((tmpName, tmpSmiles))
        print("\t".join((tmpName, tmpSmiles, tmpCanonicalSmiles, tmpInChIKey)))


if __name__ == "__main__":
    main()
