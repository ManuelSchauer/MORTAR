/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2026  Felix Baensch, Jonas Schaub (felix.j.baensch@gmail.com, jonas.schaub@uni-jena.de)
 *
 * Source code is available at <https://github.com/FelixBaensch/MORTAR>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.mortar.model.fragmentation.algorithm;

/*
 * Third-party content: Biosynfoni
 *
 * The SMARTS patterns, the key labels, and the key order of the original pattern set in this file are taken from
 * Biosynfoni, version "full_1103" (https://github.com/lucinamay/biosynfoni, files
 * src/biosynfoni/subkeys/biosmartfonis.py and src/biosynfoni/subkeys/versionfonis.py). The SMILES codes of the
 * representative structures of the NADH and NADPH keys are identical to the PubChem-based NADH and NADPH patterns of
 * Biosynfoni ("co_nadh_pubchem" and "co_nadph_pubchem" in biosmartfonis.py). The pattern variants are modifications
 * of the original patterns; the other representative structures are not taken from Biosynfoni.
 * Biosynfoni is distributed under the following licence:
 *
 * MIT License
 *
 * Copyright (c) 2023 BioSynFoni (biosynfoni)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.model.util.SimpleIDisplayEnumConstantProperty;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.value.ObservableBooleanValue;

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Fragmentation algorithm based on the substructure definitions of the Biosynfoni fingerprint, version "full_1103"
 * (<a href="https://doi.org/10.1186/s13321-025-01081-6">Biosynfoni paper (Journal of Cheminformatics, 2025)</a>,
 * <a href="https://github.com/lucinamay/biosynfoni">original Python implementation</a>; see below for the differences
 * to the reference implementation). Biosynfoni is a biosynthesis-informed, interpretable fingerprint for natural
 * products that consists of 39 SMARTS patterns describing biosynthetically relevant building blocks, e.g. cofactors,
 * amino acids, sugars, shikimate-derived phenyl units, isoprene and acetyl units, functional groups, and carbon rings.
 * <br>
 * <br>Many of these SMARTS patterns are promiscuous, i.e. they match structurally different substructures (e.g.
 * the ethyl pattern "[#6]~[#6]" matches every carbon-carbon bond regardless of bond order and aromaticity). Returning
 * the matched atoms as fragments would therefore scatter the matches of one Biosynfoni key over many different
 * fragments. Therefore, by default, every accepted match of a key is returned as one fixed representative structure
 * that is stored for this key (see {@link BiosynfoniKey#getRepresentativeSmiles(Set)}). All matches of the same key
 * thus end up as one fragment in MORTAR, and the frequency of this fragment in a molecule equals the count of the
 * respective Biosynfoni count fingerprint position. The representative structures are not parts of the input molecule
 * and partly contain pseudo atoms (e.g. "*C"); as a step of a fragmentation pipeline, the fragmenter should therefore
 * return the matched substructures.
 * <br>Alternatively, the actually matched substructures can be returned (see
 * {@link BiosynfoniReturnedFragmentsOption#MATCHED_SUBSTRUCTURES}). This shows which different substructures are
 * covered by a specific SMARTS pattern. Each matched substructure consists of the matched atoms and all bonds between
 * them, free valences are saturated with implicit hydrogen atoms, and aromaticity is re-detected on the fragment. A
 * cut bond of an aromatic system counts as a single bond, so that completely matched aromatic rings of fused systems
 * remain aromatic, e.g. the pyrrole ring of indole, and a matched atom whose aromatic bonds are all cut is not part
 * of the matched aromatic system, so that e.g. the bond between the two rings of biphenyl remains a single bond. If
 * the matched part of an aromatic system has no Kekulé structure, e.g. three atoms of a benzene ring, the fragment is
 * built from the bond orders of the Kekulé form of the input molecule instead, so that its bonds appear in this
 * Kekulé form. If the input molecule has no Kekulé structure, its unset bond orders are treated as single bonds,
 * which leaves open valences in the fragments, e.g. "[CH]1[CH][CH][CH][CH]1" for "c1cccc1". Stereochemistry is not
 * transferred to the matched substructures, since stereo elements can refer to atoms outside the match.
 * <br>Several SMARTS patterns consist of a single atom whose environment is only checked, e.g. the amino acid patterns,
 * which match only the amino nitrogen atom, or the hydroxyl, ether, terminal nitrogen, methyl, and halogen patterns.
 * Their matched substructure would be this atom alone, e.g. "N" or "O", which does not show where it was matched.
 * Therefore, the fragment of such a match additionally contains the environment of the matched atom that the pattern
 * checks: the atoms of the matches of its recursive SMARTS expressions ("$(...)") that start at the matched atom, e.g.
 * the backbone and the side chain of an amino acid, so that different amino acids give different fragments, with atoms
 * that these expressions only require as any atom ("*") as attachment points ("*"), or, if
 * the pattern has no such expression, the direct neighbours of the matched atom: neighbours in a common ring as atoms,
 * e.g. the ring of an epoxide, all other neighbours as attachment points ("*"), since the pattern only requires their
 * presence, e.g. "*C" for a methyl group or "*O*" for an ether oxygen atom. Every match still yields exactly one
 * fragment, and the overlap filters only regard the matched atoms.
 * <br><b>Known limitation of the matched substructures option:</b> MORTAR identifies fragments by their unique SMILES
 * code across the whole data set. Different keys can produce identical fragments, e.g. "C1CO1" for an epoxide and for an
 * acetyl unit matched in the epoxide ring, or, in the original pattern set, alanine for the standard amino acid
 * alanine and for the non-standard amino acid alanine methyl ester, whose ester oxygen atom belongs to the pattern.
 * Such fragments are combined into one fragment, whose frequency in a molecule is the sum over these keys. In the
 * items tab, it is shown once per key that produced it in the respective molecule, with the number of matches of this
 * key; in the fragments tab, it is grouped by all keys that produced it anywhere in the data set; its properties (e.g.
 * the Biosynfoni label in exported files) are those of only one of these keys. The number of matches per key and the
 * representative structures option are not affected.
 * <br>The label and the SMARTS pattern of the key and the method that found its matches (SMARTS pattern or ring
 * perception) are attached to every returned fragment as properties (see
 * {@link #BIOSYNFONI_LABEL_PROPERTY_KEY}, {@link #BIOSYNFONI_SMARTS_PROPERTY_KEY},
 * {@link #BIOSYNFONI_MATCHING_METHOD_PROPERTY_KEY},
 * {@link IMoleculeFragmenter#FRAGMENT_CATEGORY_PROPERTY_KEY}, and
 * {@link IMoleculeFragmenter#FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY}); the latter is used to group the fragments by
 * their key in the GUI. Like all fragment properties, they are also written to exported files, the display label in
 * the language of the GUI.
 * <br>The fragmenter has no fragment saturation setting (see {@link IMoleculeFragmenter}), like e.g. the Sugar Removal
 * Utility fragmenter: the representative structures are fixed, and the free valences of the matched substructures are
 * saturated with implicit hydrogen atoms.
 * <br>
 * <br>As in the reference implementation, two optional overlap filters are available as settings:
 * <ul>
 *   <li><b>Intra-substructure overlap filter</b>: prevents multiple matches of the same SMARTS pattern from sharing
 *   atoms. Matches are accepted greedily, as in the reference implementation. To make the result independent of the
 *   input atom order, in contrast to the reference implementation, matches are processed in the order of the
 *   canonical atom ranks of their atoms (see {@link #getCanonicalRanks(IAtomContainer)}), which are determined as in
 *   the Biosynfoni fingerprinter that is being developed for the CDK.</li>
 *   <li><b>Inter-substructure overlap filter</b>: prevents matches of a SMARTS pattern from sharing atoms with
 *   accepted matches of SMARTS patterns that come earlier in the predefined key order.</li>
 * </ul>
 * Prior to matching, ring perception and Daylight aromaticity detection are applied to a copy of the input molecule
 * (see {@link SmartsPattern#prepare(IAtomContainer)}).
 * <br>
 * <br><b>Pattern sets</b>: Several SMARTS patterns of the reference implementation (Python, RDKit, version
 * "full_1103") match more substructures than their name and their explanation in the reference implementation
 * describe, and some patterns leave room for chemical interpretation. Therefore, the pattern set is a setting (see
 * {@link BiosynfoniPatternSetOption}): the original patterns, identical to the reference implementation, or the
 * pattern variants (default), which follow the key names more closely. The variant pattern set differs from the
 * original one in 15 SMARTS patterns and in the matching of the 8 ring keys; the other 16 keys are identical. The
 * changes are grouped into 13 variants (see {@link BiosynfoniVariant}), and a third option, the custom pattern set,
 * applies only those variants whose setting is turned on, e.g. only the variants that follow the key names without
 * the interpretations. The documentation of each key describes its variant; "variant pattern set" there means that
 * the respective variant is applied.
 * <ul>
 *   <li><b>Patterns that match more than their name describes</b>: the open-chain sugar patterns also match sugar rings
 *   ({@link BiosynfoniKey#S_OPENPYR_C6O6}, {@link BiosynfoniKey#S_OPENFUR_C5O5}); the non-standard amino acid
 *   pattern lists the standard side chains instead of a generic one ({@link BiosynfoniKey#NON_STD_AMINOS}); the
 *   phenethylamine pattern matches the pyrrole ring of indoles ({@link BiosynfoniKey#D_PHENYL_C2N}); the hydroxyl
 *   pattern ignores hemiacetal hydroxyl groups and matches the oxygen atoms of nitro groups and N-oxides depending on
 *   their notation ({@link BiosynfoniKey#O_HYDROXYL}); the phosphate and sulfonate patterns count every
 *   phosphorus-oxygen and sulfur-oxygen bond instead of every group ({@link BiosynfoniKey#PHOSPHATE},
 *   {@link BiosynfoniKey#SULFONATE}); the NADH pattern also matches NADPH ({@link BiosynfoniKey#CO_NADH}); the side
 *   chain condition of the phenyl-C3 and phenyl-C2 patterns depends on the smallest ring of an atom
 *   ({@link BiosynfoniKey#D_PHENYL_C3}); the terminal nitrogen pattern also matches nitriles and imines
 *   ({@link BiosynfoniKey#N_NITRATE}); the ether pattern matches peroxides ({@link BiosynfoniKey#O_ETHER}); and the
 *   ring patterns also match the perimeter of fused and bridged ring systems, e.g. a ten-membered ring in
 *   naphthalene ({@link BiosynfoniKey#getRingSize()}).</li>
 *   <li><b>Interpretations</b>: the rings of the phenyl units from the shikimate pathway have to be aromatic
 *   ({@link BiosynfoniKey#D_PHENYL_C3}), the methyl group of a methylmalonyl unit has to be attached to a branching
 *   point ({@link BiosynfoniKey#D2_METHYLMALONYL}), and standard amino acids that are N-alkylated or esterified count
 *   as standard amino acids ({@link BiosynfoniKey#ALL_STD_AMINOS}).</li>
 * </ul>
 * The variants of the phenyl side chain condition and of the nitro and N-oxide oxygen atoms go back to the
 * Biosynfoni fingerprinter that is being developed for the CDK. Deliberately broad patterns of the reference
 * implementation, e.g. for isoprene, acetyl, and ethyl units, are kept, as is the definition of the sugar rings by
 * their number of oxygen atoms.
 * <br>
 * <br><b>Further differences to the reference implementation</b>, independent of the pattern set:
 * <ul>
 *   <li><b>Intra-substructure overlap filter</b>: the reference implementation selects the matches greedily in the
 *   order in which RDKit finds them, so that its counts depend on the input atom order of the molecule. Here, the
 *   matches are selected in a canonical order (see above), so the counts do not depend on the input atom order but can
 *   differ from the reference implementation where both select different but equally valid sets of non-overlapping
 *   matches.</li>
 *   <li><b>Cheminformatics toolkit</b>: the substructure matching is done with the CDK instead of RDKit. Explicit
 *   hydrogen atoms count as neighbours for the SMARTS primitive D and are not counted by h, so that e.g. the oxygen
 *   atom of ethanol with an explicit hydrogen atom on it is matched as ether instead of hydroxyl oxygen. RDKit removes
 *   most explicit hydrogen atoms when reading molecules. MORTAR converts explicit hydrogen atoms to implicit ones on
 *   import where possible; the remaining explicit hydrogen atoms, e.g. deuterium, tritium, or bridging hydrogen atoms,
 *   are matched as described (RDKit keeps deuterium and tritium explicit as well). If the fragmenter is used directly,
 *   explicit hydrogen atoms in the given molecules are matched as they are.</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class BiosynfoniFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum BiosynfoniKey">
    /**
     * Biosynfoni substructure keys in the order of the fingerprint positions. Each constant holds the label of the
     * key in the reference implementation, its original and, if it has a variant, its variant SMARTS pattern and representative
     * SMILES code, and the variant that switches between them (see {@link BiosynfoniVariant}). The phenyl unit and
     * hydroxyl keys are affected by two variants each, so their SMARTS patterns are composed of the original or
     * variant parts. The representative structure is matched by the SMARTS pattern of the key for every combination of
     * applied variants and is returned as fragment for every accepted match of the key. For every combination, all
     * representative structures are pairwise distinct so that every key yields exactly one fragment.
     */
    public enum BiosynfoniKey {
        //cofactors
        /**
         * Coenzyme A (CoA), matches the CoA backbone including the pantetheine chain, diphosphate linker, and
         * adenosine moiety, also in acyl-CoA thioesters.
         */
        CO_COA("co_coa",
                "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O",
                "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"),
        /**
         * Reduced nicotinamide adenine dinucleotide (NADH); the pattern also matches the oxidised form
         * NAD<sup>+</sup>. In the variant pattern set, the 2'-oxygen atom of the adenosine ribose must not be bonded
         * to a phosphorus atom ("[O;!$(*P)]"). The original pattern also matches NADPH, so that NADPH counts as NADH
         * as well, and with the inter-substructure overlap filter, the following NADPH key can never match, contrary
         * to the key order of the reference implementation that puts the largest cofactors first.
         */
        CO_NADH("co_nadh",
                BiosynfoniVariant.NADH,
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~[O;!$(*P)])~O)~O)~O",
                "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"),
        /**
         * Reduced nicotinamide adenine dinucleotide phosphate (NADPH); the pattern also matches the oxidised form
         * NADP<sup>+</sup>.
         */
        CO_NADPH("co_nadph",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O",
                "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"),
        //amino acids
        /**
         * Standard proteinogenic amino acids (free or peptide-bound), pattern taken from the Daylight SMARTS examples;
         * the pattern matches only the amino nitrogen atom. Represented by glycine. According to the explanation of
         * the reference implementation, "standard" refers to the 20 standard side chains. The original pattern,
         * however, does not match residues with a standard side chain whose amino nitrogen atom carries a further
         * carbon substituent, e.g. N-methylated residues in peptides, or whose carboxy group is esterified, e.g. in
         * depsipeptides. In the variant pattern set, these residues are matched, since N-methylation and
         * esterification are modifications of standard amino acids.
         */
        ALL_STD_AMINOS("allstnd_aminos",
                BiosynfoniVariant.AMINO_ACIDS,
                "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]",
                "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H2][CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2,OX1-,N])]",
                "NCC(=O)O"),
        /**
         * Non-standard amino acids, i.e. alpha-amino acids that are not matched by {@link #ALL_STD_AMINOS}; the
         * pattern matches only the amino nitrogen atom. According to the explanation of the reference implementation,
         * the first part of the pattern describes a generic amino acid ("generic amino acid but not a standard amino
         * acid"). In the original pattern, however, it lists the 20 standard side chains, so that only amino acids
         * with a standard side chain but a modified amino or carboxy group are matched, e.g. alanine methyl ester, and
         * no amino acid with a non-standard side chain, e.g. ornithine or 2-aminobutyric acid. In the variant pattern
         * set, the first part describes a generic alpha-amino acid, including alpha,alpha-disubstituted ones like
         * 2-aminoisobutyric acid. Represented by alanine methyl ester in the original and by 2-aminobutyric acid in
         * the variant pattern set.
         */
        NON_STD_AMINOS("nonstnd_aminos",
                BiosynfoniVariant.AMINO_ACIDS,
                "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]",
                "[$([NX3,NX4+][CX4;H1,H0]([*])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H2][CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2,OX1-,N])])]",
                "COC(=O)C(C)N",
                "CCC(N)C(=O)O"),
        //sugar derivatives
        /**
         * Open-chain hexose, linear chain of six carbon atoms that each carry an oxygen atom. In the variant pattern
         * set, the carbon-carbon bonds must not be ring bonds ("!@"), as the name requires. The original pattern also
         * matches cyclic hexoses, e.g. glucopyranose, so that with the inter-substructure overlap filter,
         * hexopyranoses are counted as open-chain hexoses and the following pyranose key cannot match them.
         */
        S_OPENPYR_C6O6("s_openpyr_C6O6",
                BiosynfoniVariant.OPEN_CHAIN_SUGARS,
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])",
                "OCC(O)C(O)C(O)C(O)C=O"),
        /**
         * Open-chain pentose, linear chain of five carbon atoms that each carry an oxygen atom. In the variant pattern
         * set, the carbon-carbon bonds must not be ring bonds ("!@"), as the name and the explanation of the reference
         * implementation ("linear furanose, in open form") require. The original pattern also matches furanoses and
         * hexopyranoses.
         */
        S_OPENFUR_C5O5("s_openfur_C5O5",
                BiosynfoniVariant.OPEN_CHAIN_SUGARS,
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])",
                "OCC(O)C(O)C(O)C=O"),
        /**
         * Pyranose, six-membered sugar ring containing one oxygen atom, with oxygen atoms at three neighbouring ring
         * carbon atoms. In accordance with the composition "C5O4" in the label, 2-deoxy sugars are not matched.
         */
        S_PYRANOSE_C5O4("s_pyranose_C5O4",
                "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1",
                "C1OCC(O)C(O)C1O"),
        /**
         * Furanose, five-membered sugar ring containing one oxygen atom, with oxygen atoms at two neighbouring ring
         * carbon atoms. In accordance with the composition "C4O3" in the label, 2-deoxyribose is not matched.
         */
        S_FURANOSE_C4O3("s_furanose_C4O3",
                "C~1~[#8]~C~C(~[#8])~C(~[#8])~1",
                "C1OCC(O)C1O"),
        //nitrogen-containing building blocks
        /**
         * Indole scaffold with a two-carbon side chain terminating in a nitrogen atom (tryptamine-like).
         */
        D_INDOLE("d_indoleC2N_12",
                "c1cccc2c1c(~[#6]~[#6]~[#7])cn2",
                "NCCc1c[nH]c2ccccc12"),
        /**
         * Phenethylamine, benzene ring connected to a two-carbon chain ending in a nitrogen atom. In the original
         * pattern, the chain is written without bond symbols, i.e. with single or aromatic bonds, and may form a
         * five-membered ring with the benzene ring, so that the pyrrole ring of every indole (C3a-C3-C2-N1) is matched
         * as phenethylamine unit. In the variant pattern set, the first chain carbon atom must not be part of a
         * five-membered ring, and the chain is written with any-bond expressions. Isoquinolines, whose phenethylamine
         * unit forms a six-membered ring with the benzene ring, are still matched.
         */
        D_PHENYL_C2N("d_phenylC2N_9",
                BiosynfoniVariant.PHENETHYLAMINE,
                "c1ccccc1[#6][#6][#7]",
                "c1ccccc1[#6;!$(*~1~*~*~*~*~1)]~[#6]~[#7]",
                "NCCc1ccccc1"),
        /**
         * Six-membered ring of five carbon atoms and one nitrogen atom, regardless of aromaticity.
         */
        D_C5N("d_c5n_6",
                "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1",
                "C1CCNCC1"),
        /**
         * Five-membered ring of four carbon atoms and one nitrogen atom, regardless of aromaticity.
         */
        D_C4N("d_c4n_5",
                "[#6]~1~[#6]~[#6]~[#6]~[#7]~1",
                "C1CCNC1"),
        //phenyls from the shikimate pathway
        /**
         * Phenyl ring attached to a three-carbon side chain, i.e. a C6-C3 unit from the shikimate pathway. The ring
         * atoms must not be part of further rings ("R1") to prevent fused terpenoid rings from being identified, and
         * the side chain carbon atoms after the first one must not be part of a six-membered ring.
         * <br>In the original pattern set, the ring does not have to be aromatic, and the side chain atoms are
         * excluded via "[#6;!$([r6])]" only if their smallest ring is six-membered, so that, e.g., carbon atoms in a
         * six-membered ring fused to the five-membered lactone ring of aryltetralin lignans are not excluded.
         * <br>The two variants are applied independently, see {@link BiosynfoniVariant#AROMATIC_PHENYL_RINGS} and
         * {@link BiosynfoniVariant#PHENYL_SIDE_CHAIN}.
         * <br>In the variant pattern set, the ring has to be aromatic, so that cyclohexyl groups and monocyclic
         * terpenes, e.g. menthane derivatives, are not matched; non-aromatic shikimate derivatives like shikimic acid
         * are therefore not matched either. The side chain atoms must not be part of any six-membered ring, which is
         * written with any-bond expressions ("~") so that rings with non-aromatic double bonds, e.g. cyclohexene, are
         * recognised as well; without bond symbols, a bond means single or aromatic in SMARTS.
         */
        D_PHENYL_C3("d_phenylC3_9_strict", "CCCc1ccccc1"),
        /**
         * Phenyl ring attached to a two-carbon side chain, i.e. a C6-C2 unit. The second side chain carbon atom must
         * not be part of a six-membered ring. The differences between the pattern sets are described at {@link
         * #D_PHENYL_C3}.
         */
        D_PHENYL_C2("d_phenylC2_8_strict", "CCc1ccccc1"),
        /**
         * Phenyl ring attached to a single carbon atom, i.e. a C6-C1 unit. In the variant pattern set, the ring has to
         * be aromatic, see {@link #D_PHENYL_C3}.
         */
        D_PHENYL_C1("d_phenylC1_7_strict", "Cc1ccccc1"),
        //carbon chains
        /**
         * Isoprene unit, branched five-carbon motif. According to the explanation of the reference implementation, the
         * pattern deliberately ignores double bonds to detect isoprene units within fused rings. Represented by
         * isoprene (2-methylbuta-1,3-diene).
         */
        D_ISOPRENE("d_isoprene_5",
                "[#6]~[#6](~[#6])~[#6]~[#6]",
                "C=CC(C)=C"),
        /**
         * Acetyl unit, two-carbon chain connected to an oxygen atom. Represented by acetaldehyde, i.e. with the
         * carbonyl group of the acetyl unit.
         */
        D2_ACETYL("d2_acetyl_C2O1",
                "[#6]~[#6]~[#8]",
                "CC=O"),
        /**
         * Methylmalonyl unit, three-carbon motif terminating in a methyl carbon atom. The original pattern matches
         * every C-C-CH<sub>3</sub> motif, including chain ends, methyl groups on aromatic rings, and methyl ketones,
         * i.e. typically acetate-derived structures. In the variant pattern set, the methyl group must be attached to
         * an aliphatic carbon atom with at least three heavy-atom neighbours that is not a carbonyl carbon atom, i.e.
         * to a branching point as in polyketide chains built from methylmalonyl units; terpenoid methyl groups are
         * still matched, but with the inter-substructure overlap filter, they are taken by the isoprene key first.
         * Represented by propane in the original and by an isopropyl group attached to a pseudo atom in the variant
         * pattern set, whose pattern requires a further neighbour.
         */
        D2_METHYLMALONYL("d2_methylmalonyl_C3",
                BiosynfoniVariant.METHYLMALONYL,
                "[#6]~[#6][C;D1;h3]",
                "[#6]~[C;D3,D4;!$(*=[#8])]-[C;D1;h3]",
                "CCC",
                "*C(C)C"),
        /**
         * Ethyl unit, two connected carbon atoms.
         */
        D_ETHYL("d_ethyl_2",
                "[#6]~[#6]",
                "CC"),
        /**
         * Terminal methyl carbon. Represented by a methyl group attached to a pseudo atom because the pattern requires
         * exactly one neighbour.
         */
        D_METHYL("d_methyl_1",
                "[C;D1;h3]",
                "*C"),
        //phosphorus- and sulfur-containing groups
        /**
         * Phosphorus atom bonded to an oxygen atom. The original pattern "P~O" matches every phosphorus-oxygen pair,
         * so that one phosphate group is counted four times, although the explanation of the reference implementation
         * states that additional oxygen atoms are not matched. In the variant pattern set, only the phosphorus atom is
         * matched, so that every phosphorus atom bonded to at least one oxygen atom is counted once.
         */
        PHOSPHATE("phosphate_2",
                BiosynfoniVariant.PHOSPHATE_SULFONATE,
                "P~O",
                "[P;$(*~[#8])]",
                "OP"),
        /**
         * Sulfur atom bonded to an oxygen atom, e.g. in sulfoxides, sulfones, sulfonates, and sulfates. The original
         * pattern "S~O" matches every sulfur-oxygen pair, e.g. four times for one sulfate group. In the variant
         * pattern set, only the sulfur atom is matched, so that every sulfur atom bonded to at least one oxygen atom
         * is counted once.
         */
        SULFONATE("sulfonate_2",
                BiosynfoniVariant.PHOSPHATE_SULFONATE,
                "S~O",
                "[S;$(*~[#8])]",
                "OS"),
        //halides
        /**
         * Fluorine atom.
         */
        HAL_F("hal_f", "[#9]", "F"),
        /**
         * Chlorine atom.
         */
        HAL_CL("hal_cl", "[#17]", "Cl"),
        /**
         * Bromine atom.
         */
        HAL_BR("hal_br", "[#35]", "Br"),
        /**
         * Iodine atom.
         */
        HAL_I("hal_i", "[#53]", "I"),
        //functional groups
        /**
         * Terminal aliphatic nitrogen atom with exactly one neighbour (not restricted to nitrate groups despite the
         * name). Represented by an amino group attached to a pseudo atom because the pattern requires exactly one
         * neighbour. In the variant pattern set, the nitrogen atom must be bonded by a single bond, as stated in the
         * explanation of the reference implementation ("N with one single bond"), so that nitriles and imines are not
         * matched.
         */
        N_NITRATE("n_nitrate_1",
                BiosynfoniVariant.TERMINAL_NITROGEN,
                "[N;D1]",
                "[N;D1;!$(*=*);!$(*#*)]",
                "*N"),
        /**
         * Oxygen atom of an epoxide (three-membered ring).
         */
        O_EPOXY("o_epoxy_1",
                "[O;x2;r3]",
                "C1CO1"),
        /**
         * Acyclic ether oxygen atom, excluding esters, phosphates, and sulfates. In the variant pattern set, oxygen
         * atoms bonded to a further oxygen atom, i.e. acyclic peroxides, are excluded as well.
         */
        O_ETHER("o_ether_1",
                BiosynfoniVariant.ETHER,
                "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]",
                "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S);!$(*[#8])]",
                "COC"),
        /**
         * Hydroxyl oxygen atom (-OH) or deprotonated hydroxyl oxygen atom (-O<sup>-</sup>) attached to a carbon or
         * nitrogen atom, excluding carboxylic acids. The original pattern excludes oxygen atoms bonded to a carbon
         * atom that carries a further oxygen atom ("!$(*C~O)"), which also excludes the hydroxyl groups of hemiacetals
         * and hemiketals, e.g. the anomeric hydroxyl group of reducing sugars; its parts "!$(P)" and "!$(S)" have no
         * effect because they test the oxygen atom itself.
         * <br>In the variant pattern set, only oxygen atoms bonded to a carbonyl carbon atom are excluded
         * ("!$(*C=O)"). Additionally, the oxygen atoms of nitro groups, nitrate esters, and N-oxides are excluded:
         * oxygen atoms bonded to a nitrogen atom that carries a further oxygen atom ("!$(*[#7]~[#8])") and negatively
         * charged oxygen atoms bonded to a positively charged nitrogen atom ("!$([#8-][#7+])"). These groups can be
         * written in a pentavalent form ("N(=O)=O", "N=O"), whose doubly bonded oxygen atom is not matched anyway, or
         * in a charge-separated form ("[N+](=O)[O-]", "[N+][O-]"), whose singly bonded "[O-]" would otherwise be
         * matched as deprotonated hydroxyl group, so that the result depended on the notation of the same molecule.
         * The reference implementation matches this "[O-]" because RDKit converts nitro groups and N-oxides into the
         * charge-separated form. Hydroxylamines, oximes, and hydroxamic acids, whose uncharged nitrogen atom carries
         * only one oxygen atom, are still matched. The exclusion of these oxygen atoms is taken from the Biosynfoni
         * fingerprinter that is being developed for the CDK.
         */
        O_HYDROXYL("o_hydroxyl_1", "CO"),
        //carbon rings
        /**
         * Three-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by
         * ring perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1", "C1CC1", 3),
        /**
         * Four-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by ring
         * perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1", "C1CCC1", 4),
        /**
         * Five-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by ring
         * perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1", "C1CCCC1", 5),
        /**
         * Six-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by ring
         * perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCC1", 6),
        /**
         * Seven-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by
         * ring perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCC1", 7),
        /**
         * Eight-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by
         * ring perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCC1", 8),
        /**
         * Nine-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by ring
         * perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCCC1", 9),
        /**
         * Ten-membered carbon ring, regardless of aromaticity. In the variant pattern set, the key is matched by ring
         * perception instead of its SMARTS pattern, see {@link #getRingSize()}.
         */
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCCCC1", 10);

        /**
         * Phenyl ring of the original phenyl unit patterns, six carbon atoms that are not part of further rings.
         */
        private static final String PHENYL_RING_SMARTS = "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1";
        /**
         * Phenyl ring of the phenyl unit pattern variants, six aromatic carbon atoms that are not part of further
         * rings, see {@link BiosynfoniVariant#AROMATIC_PHENYL_RINGS}.
         */
        private static final String AROMATIC_PHENYL_RING_SMARTS = "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1";
        /**
         * Side chain carbon atom of the original phenyl unit patterns that is not part of a ring whose smallest ring
         * is six-membered.
         */
        private static final String PHENYL_SIDE_CHAIN_ATOM_SMARTS = "[#6;!$([r6])]";
        /**
         * Side chain carbon atom of the phenyl unit pattern variants that is not part of any six-membered ring, see
         * {@link BiosynfoniVariant#PHENYL_SIDE_CHAIN}.
         */
        private static final String VARIANT_PHENYL_SIDE_CHAIN_ATOM_SMARTS = "[#6;!$(*~1~*~*~*~*~*~1)]";
        /**
         * Label of the key as used in the reference implementation.
         */
        private final String label;
        /**
         * Variant that switches between the original and the variant SMARTS pattern and representative structure of
         * the key; null if the key is unchanged or composed of several variants.
         */
        private final BiosynfoniVariant variant;
        /**
         * Original SMARTS pattern of the key; null for keys whose SMARTS pattern is composed.
         */
        private final String originalSmarts;
        /**
         * Variant SMARTS pattern of the key; null for keys whose SMARTS pattern is composed.
         */
        private final String variantSmarts;
        /**
         * SMILES code of the original representative structure of the key.
         */
        private final String originalRepresentativeSmiles;
        /**
         * SMILES code of the variant representative structure of the key.
         */
        private final String variantRepresentativeSmiles;
        /**
         * Number of ring atoms of a ring key, 0 for all other keys.
         */
        private final int ringSize;
        /**
         * Constructor for keys whose SMARTS pattern is composed of original or variant parts depending on several
         * variants (see {@link #getSmarts(Set)}).
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aRepresentativeSmiles SMILES code of the representative structure of the key
         */
        private BiosynfoniKey(String aLabel, String aRepresentativeSmiles) {
            this(aLabel, null, null, null, aRepresentativeSmiles, aRepresentativeSmiles, 0);
        }
        //
        /**
         * Constructor for keys without a variant.
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aSmarts SMARTS pattern of the key
         * @param aRepresentativeSmiles SMILES code of the representative structure of the key
         */
        private BiosynfoniKey(String aLabel, String aSmarts, String aRepresentativeSmiles) {
            this(aLabel, null, aSmarts, aSmarts, aRepresentativeSmiles, aRepresentativeSmiles, 0);
        }
        //
        /**
         * Constructor for keys whose SMARTS pattern has a variant while the representative structure fits both patterns.
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aVariant variant that switches between the original and the variant SMARTS pattern
         * @param anOriginalSmarts original SMARTS pattern of the key
         * @param aVariantSmarts variant SMARTS pattern of the key
         * @param aRepresentativeSmiles SMILES code of the representative structure of the key
         */
        private BiosynfoniKey(String aLabel, BiosynfoniVariant aVariant, String anOriginalSmarts, String aVariantSmarts,
                      String aRepresentativeSmiles) {
            this(aLabel, aVariant, anOriginalSmarts, aVariantSmarts, aRepresentativeSmiles, aRepresentativeSmiles, 0);
        }
        //
        /**
         * Constructor for keys whose SMARTS pattern and representative structure have a variant.
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aVariant variant that switches between the original and the variant SMARTS pattern
         * @param anOriginalSmarts original SMARTS pattern of the key
         * @param aVariantSmarts variant SMARTS pattern of the key
         * @param anOriginalRepresentativeSmiles SMILES code of the original representative structure
         * @param aVariantRepresentativeSmiles SMILES code of the variant representative structure
         */
        private BiosynfoniKey(String aLabel, BiosynfoniVariant aVariant, String anOriginalSmarts, String aVariantSmarts,
                      String anOriginalRepresentativeSmiles, String aVariantRepresentativeSmiles) {
            this(aLabel, aVariant, anOriginalSmarts, aVariantSmarts, anOriginalRepresentativeSmiles,
                    aVariantRepresentativeSmiles, 0);
        }
        //
        /**
         * Constructor for ring keys, whose SMARTS pattern and representative structure have no variant, but which are
         * matched by ring perception if {@link BiosynfoniVariant#RING_PERCEPTION} is applied (see
         * {@link #getRingSize()}).
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aSmarts SMARTS pattern of the key
         * @param aRepresentativeSmiles SMILES code of the representative structure of the key
         * @param aRingSize number of ring atoms
         */
        private BiosynfoniKey(String aLabel, String aSmarts, String aRepresentativeSmiles, int aRingSize) {
            this(aLabel, BiosynfoniVariant.RING_PERCEPTION, aSmarts, aSmarts, aRepresentativeSmiles,
                    aRepresentativeSmiles, aRingSize);
        }
        //
        /**
         * Constructor that sets all fields.
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aVariant variant that switches between the original and the variant of the key; null if
         *                  the key is unchanged or its SMARTS pattern is composed
         * @param anOriginalSmarts original SMARTS pattern of the key; null if composed
         * @param aVariantSmarts variant SMARTS pattern of the key; null if composed
         * @param anOriginalRepresentativeSmiles SMILES code of the original representative structure
         * @param aVariantRepresentativeSmiles SMILES code of the variant representative structure
         * @param aRingSize number of ring atoms of a ring key, 0 for all other keys
         */
        private BiosynfoniKey(String aLabel, BiosynfoniVariant aVariant, String anOriginalSmarts, String aVariantSmarts,
                      String anOriginalRepresentativeSmiles, String aVariantRepresentativeSmiles, int aRingSize) {
            this.label = aLabel;
            this.variant = aVariant;
            this.originalSmarts = anOriginalSmarts;
            this.variantSmarts = aVariantSmarts;
            this.originalRepresentativeSmiles = anOriginalRepresentativeSmiles;
            this.variantRepresentativeSmiles = aVariantRepresentativeSmiles;
            this.ringSize = aRingSize;
        }
        //
        /**
         * Returns the label of the key as used in the reference implementation.
         *
         * @return label
         */
        public String getLabel() {
            return this.label;
        }
        //
        /**
         * Returns the SMARTS pattern of the key for the given applied variants. No variants give the original
         * pattern of the reference implementation, all variants the variant pattern.
         *
         * @param aVariantSet the applied variants
         * @return SMARTS pattern
         * @throws NullPointerException if the given set is null
         */
        public String getSmarts(Set<BiosynfoniVariant> aVariantSet) throws NullPointerException {
            Objects.requireNonNull(aVariantSet, "Given set of variants is null.");
            return switch (this) {
                case D_PHENYL_C3 -> BiosynfoniKey.composePhenylSmarts(2, aVariantSet);
                case D_PHENYL_C2 -> BiosynfoniKey.composePhenylSmarts(1, aVariantSet);
                case D_PHENYL_C1 -> BiosynfoniKey.composePhenylSmarts(0, aVariantSet);
                case O_HYDROXYL -> BiosynfoniKey.composeHydroxylSmarts(aVariantSet);
                default -> this.isVariantApplied(aVariantSet) ? this.variantSmarts : this.originalSmarts;
            };
        }
        //
        /**
         * Returns the SMILES code of the representative structure that is returned as fragment for every match of
         * the key with the given applied variants.
         *
         * @param aVariantSet the applied variants
         * @return representative SMILES code
         * @throws NullPointerException if the given set is null
         */
        public String getRepresentativeSmiles(Set<BiosynfoniVariant> aVariantSet) throws NullPointerException {
            Objects.requireNonNull(aVariantSet, "Given set of variants is null.");
            return this.isVariantApplied(aVariantSet) ? this.variantRepresentativeSmiles : this.originalRepresentativeSmiles;
        }
        //
        /**
         * Returns the number of ring atoms of a ring key, or 0 for all other keys.
         * <br>A ring SMARTS pattern like "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1" matches every closed path of the given
         * size, also along the perimeter of fused or bridged ring systems, e.g. the ten carbon atoms of the perimeter
         * of naphthalene, decalin, or steroids, or the six carbon atoms of the perimeter of bicyclo[3.1.0]hexane. If
         * {@link BiosynfoniVariant#RING_PERCEPTION} is applied, ring keys are therefore matched by ring perception
         * instead of their SMARTS pattern: every relevant cycle (see {@link Cycles#relevant(IAtomContainer)}) of the
         * given size that consists only of carbon atoms is a match. The relevant cycles are the union of all minimum
         * cycle bases, so they do not depend on the input atom order, and they include all faces of symmetric ring
         * systems like cubane.
         *
         * @return number of ring atoms, 0 if the key is not a ring key
         */
        public int getRingSize() {
            return this.ringSize;
        }
        //
        /**
         * Returns the language-specific, readable name of the key that is displayed in the GUI.
         *
         * @return display name
         */
        public String getDisplayName() {
            return Message.get("BiosynfoniFragmenter.BiosynfoniKey." + this.name() + ".displayName");
        }
        //
        /**
         * Checks whether the single variant of this key is contained in the given set of applied variants.
         *
         * @param aVariantSet the applied variants
         * @return true if the key has a single variant and it is applied
         */
        private boolean isVariantApplied(Set<BiosynfoniVariant> aVariantSet) {
            return !Objects.isNull(this.variant) && aVariantSet.contains(this.variant);
        }
        //
        /**
         * Composes the SMARTS pattern of a phenyl unit key from the original or variant phenyl ring (see
         * {@link BiosynfoniVariant#AROMATIC_PHENYL_RINGS}) and the original or variant side chain condition (see
         * {@link BiosynfoniVariant#PHENYL_SIDE_CHAIN}).
         *
         * @param aNumberOfConditionedSideChainAtoms number of side chain atoms after the first one, i.e. 2 for the
         *                                           phenyl-C3, 1 for the phenyl-C2, and 0 for the phenyl-C1 unit
         * @param aVariantSet the applied variants
         * @return SMARTS pattern
         */
        private static String composePhenylSmarts(int aNumberOfConditionedSideChainAtoms, Set<BiosynfoniVariant> aVariantSet) {
            StringBuilder tmpSmarts = new StringBuilder(aVariantSet.contains(BiosynfoniVariant.AROMATIC_PHENYL_RINGS)
                    ? BiosynfoniKey.AROMATIC_PHENYL_RING_SMARTS : BiosynfoniKey.PHENYL_RING_SMARTS);
            tmpSmarts.append("~[#6]");
            String tmpSideChainAtom = aVariantSet.contains(BiosynfoniVariant.PHENYL_SIDE_CHAIN)
                    ? BiosynfoniKey.VARIANT_PHENYL_SIDE_CHAIN_ATOM_SMARTS : BiosynfoniKey.PHENYL_SIDE_CHAIN_ATOM_SMARTS;
            for (int i = 0; i < aNumberOfConditionedSideChainAtoms; i++) {
                tmpSmarts.append('~').append(tmpSideChainAtom);
            }
            return tmpSmarts.toString();
        }
        //
        /**
         * Composes the SMARTS pattern of the hydroxyl key from the original or variant exclusion of oxygen atoms bonded
         * to a carbon atom with a further oxygen atom (see {@link BiosynfoniVariant#HYDROXYL_HEMIACETAL}) and the
         * original or variant remaining exclusions (see {@link BiosynfoniVariant#HYDROXYL_NITRO_N_OXIDE}). The
         * original remaining exclusions "!$(P)" and "!$(S)" have no effect but are kept to reproduce the original
         * pattern exactly.
         *
         * @param aVariantSet the applied variants
         * @return SMARTS pattern
         */
        private static String composeHydroxylSmarts(Set<BiosynfoniVariant> aVariantSet) {
            return "[#8;D1;h,!v2;$(*[#6,#7]);"
                    + (aVariantSet.contains(BiosynfoniVariant.HYDROXYL_HEMIACETAL) ? "!$(*C=O)" : "!$(*C~O)")
                    + (aVariantSet.contains(BiosynfoniVariant.HYDROXYL_NITRO_N_OXIDE)
                            ? ";!$(*[#7]~[#8]);!$([#8-][#7+])" : ";!$(P);!$(S)")
                    + "]";
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum BiosynfoniVariant">
    /**
     * Variants of the original Biosynfoni SMARTS patterns. The variant pattern set applies all variants, the
     * original pattern set none, and the custom pattern set those whose setting is turned on (see
     * {@link BiosynfoniPatternSetOption}). Related changes form one variant so that no contradictory combinations
     * arise, e.g. both amino acid keys or both open-chain sugar keys. The documentation of the affected keys describes
     * each variant in detail.
     */
    public enum BiosynfoniVariant {
        /**
         * NADH pattern that does not match NADPH ({@link BiosynfoniKey#CO_NADH}).
         */
        NADH("NADH pattern variant setting"),
        /**
         * Standard amino acids including N-alkylated and esterified residues, and non-standard amino acids with a
         * generic side chain ({@link BiosynfoniKey#ALL_STD_AMINOS}, {@link BiosynfoniKey#NON_STD_AMINOS}).
         */
        AMINO_ACIDS("Amino acid pattern variants setting"),
        /**
         * Open-chain sugar patterns that do not match sugar rings ({@link BiosynfoniKey#S_OPENPYR_C6O6},
         * {@link BiosynfoniKey#S_OPENFUR_C5O5}).
         */
        OPEN_CHAIN_SUGARS("Open-chain sugar pattern variants setting"),
        /**
         * Phenethylamine pattern that does not match the pyrrole ring of indoles ({@link BiosynfoniKey#D_PHENYL_C2N}).
         */
        PHENETHYLAMINE("Phenethylamine pattern variant setting"),
        /**
         * Side chain atoms of the phenyl-C3 and phenyl-C2 units that must not be part of any six-membered ring instead
         * of a six-membered smallest ring ({@link BiosynfoniKey#D_PHENYL_C3}, {@link BiosynfoniKey#D_PHENYL_C2}).
         */
        PHENYL_SIDE_CHAIN("Phenyl side chain variant setting"),
        /**
         * Aromatic rings of the phenyl units ({@link BiosynfoniKey#D_PHENYL_C3}, {@link BiosynfoniKey#D_PHENYL_C2},
         * {@link BiosynfoniKey#D_PHENYL_C1}).
         */
        AROMATIC_PHENYL_RINGS("Aromatic phenyl rings variant setting"),
        /**
         * Methylmalonyl units whose methyl group is attached to a branching point
         * ({@link BiosynfoniKey#D2_METHYLMALONYL}).
         */
        METHYLMALONYL("Methylmalonyl pattern variant setting"),
        /**
         * One match per phosphorus or sulfur atom instead of per bond to an oxygen atom ({@link BiosynfoniKey#PHOSPHATE},
         * {@link BiosynfoniKey#SULFONATE}).
         */
        PHOSPHATE_SULFONATE("Phosphate and sulfonate pattern variants setting"),
        /**
         * Terminal nitrogen atoms with a single bond only ({@link BiosynfoniKey#N_NITRATE}).
         */
        TERMINAL_NITROGEN("Terminal nitrogen pattern variant setting"),
        /**
         * Ether pattern that does not match peroxides ({@link BiosynfoniKey#O_ETHER}).
         */
        ETHER("Ether pattern variant setting"),
        /**
         * Hydroxyl pattern that matches the hydroxyl groups of hemiacetals and hemiketals
         * ({@link BiosynfoniKey#O_HYDROXYL}).
         */
        HYDROXYL_HEMIACETAL("Hydroxyl hemiacetal variant setting"),
        /**
         * Hydroxyl pattern that does not match the oxygen atoms of nitro groups, nitrate esters, and N-oxides
         * ({@link BiosynfoniKey#O_HYDROXYL}).
         */
        HYDROXYL_NITRO_N_OXIDE("Hydroxyl nitro and N-oxide variant setting"),
        /**
         * Ring keys matched by ring perception instead of their SMARTS pattern (see {@link BiosynfoniKey#getRingSize()}).
         */
        RING_PERCEPTION("Ring perception variant setting");
        /**
         * Name of the setting that turns the variant on or off in the custom pattern set.
         */
        private final String settingName;
        /**
         * Constructor.
         *
         * @param aSettingName name of the setting that turns the variant on or off in the custom pattern set
         */
        private BiosynfoniVariant(String aSettingName) {
            this.settingName = aSettingName;
        }
        //
        /**
         * Returns the name of the setting that turns the variant on or off in the custom pattern set.
         *
         * @return setting name
         */
        public String getSettingName() {
            return this.settingName;
        }
        //
        /**
         * Returns the language-specific display name of the setting of the variant.
         *
         * @return display name
         */
        public String getDisplayName() {
            return Message.get("BiosynfoniFragmenter.BiosynfoniVariant." + this.name() + ".displayName");
        }
        //
        /**
         * Returns the language-specific tooltip text of the setting of the variant.
         *
         * @return tooltip text
         */
        public String getTooltipText() {
            return Message.get("BiosynfoniFragmenter.BiosynfoniVariant." + this.name() + ".tooltip");
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum BiosynfoniReturnedFragmentsOption">
    /**
     * Enum for the options of which fragments are returned for the matches of the Biosynfoni SMARTS patterns.
     */
    public enum BiosynfoniReturnedFragmentsOption implements IDisplayEnum {
        /**
         * Every match of a key is returned as the stored representative structure of the key, so all matches of a key
         * are combined into one fragment.
         */
        REPRESENTATIVE_STRUCTURES(
                Message.get("BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.RepresentativeStructures.displayName"),
                Message.get("BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.RepresentativeStructures.tooltip")),
        /**
         * Every match of a key is returned as the substructure formed by the matched atoms, so the different
         * substructures matched by a key become visible. A match of a single-atom pattern additionally contains the
         * environment of the atom that the pattern checks, and identical fragments of different keys are combined, see
         * the class documentation.
         */
        MATCHED_SUBSTRUCTURES(
                Message.get("BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MatchedSubstructures.displayName"),
                Message.get("BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MatchedSubstructures.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        private BiosynfoniReturnedFragmentsOption(String aDisplayName, String aTooltip) {
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        //
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        //
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum BiosynfoniPatternSetOption">
    /**
     * Enum for the options of which SMARTS patterns are used: the original ones of Biosynfoni version "full_1103", the
     * variant ones, or a custom selection of variants (see class documentation and {@link BiosynfoniVariant}).
     */
    public enum BiosynfoniPatternSetOption implements IDisplayEnum {
        /**
         * Pattern variants, which follow the key names more closely than the original patterns and interpret some of
         * them chemically, i.e. all variants are applied.
         */
        VARIANT(
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Variant.displayName"),
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Variant.tooltip")),
        /**
         * Original SMARTS patterns of Biosynfoni version "full_1103", identical to the reference implementation.
         */
        ORIGINAL(
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Original.displayName"),
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Original.tooltip")),
        /**
         * Custom selection of variants: a variant is applied if its setting is turned on (see
         * {@link BiosynfoniVariant}).
         */
        CUSTOM(
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Custom.displayName"),
                Message.get("BiosynfoniFragmenter.BiosynfoniPatternSetOption.Custom.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        private BiosynfoniPatternSetOption(String aDisplayName, String aTooltip) {
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        //
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        //
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm used in this fragmenter.
     */
    public static final String ALGORITHM_NAME = "Biosynfoni";

    /**
     * Property key under which the Biosynfoni label of the matched key is stored in the returned fragments.
     */
    public static final String BIOSYNFONI_LABEL_PROPERTY_KEY = "BiosynfoniFragmenter.Label";

    /**
     * Property key under which the SMARTS pattern of the matched key is stored in the returned fragments. For ring keys
     * that are matched by ring perception (see {@link #BIOSYNFONI_MATCHING_METHOD_PROPERTY_KEY}), it is the ring
     * pattern of the key, which matches the fragment but was not used to find the matches.
     */
    public static final String BIOSYNFONI_SMARTS_PROPERTY_KEY = "BiosynfoniFragmenter.SMARTS";

    /**
     * Property key under which the method that found the matches of the key is stored in the returned fragments,
     * either {@link #MATCHING_METHOD_SMARTS} or {@link #MATCHING_METHOD_RING_PERCEPTION}.
     */
    public static final String BIOSYNFONI_MATCHING_METHOD_PROPERTY_KEY = "BiosynfoniFragmenter.MatchingMethod";

    /**
     * Value of the matching method property for keys that are matched by their SMARTS pattern.
     */
    public static final String MATCHING_METHOD_SMARTS = "SMARTS pattern";

    /**
     * Value of the matching method property for ring keys that are matched by ring perception, i.e. by the relevant
     * carbocycles of their size (see {@link BiosynfoniKey#getRingSize()}).
     */
    public static final String MATCHING_METHOD_RING_PERCEPTION = "ring perception (relevant carbocycles)";

    /**
     * Default value of the intra-substructure overlap filter setting, taken from the default configuration of the
     * Biosynfoni fingerprinter.
     */
    public static final boolean INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT = false;

    /**
     * Default value of the inter-substructure overlap filter setting, taken from the default configuration of the
     * Biosynfoni fingerprinter.
     */
    public static final boolean INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT = false;

    /**
     * Default option for the returned fragments setting.
     */
    public static final BiosynfoniReturnedFragmentsOption RETURNED_FRAGMENTS_OPTION_DEFAULT =
            BiosynfoniReturnedFragmentsOption.REPRESENTATIVE_STRUCTURES;

    /**
     * Default option for the pattern set setting.
     */
    public static final BiosynfoniPatternSetOption PATTERN_SET_OPTION_DEFAULT = BiosynfoniPatternSetOption.VARIANT;

    /**
     * Default value of the settings of the variants in the custom pattern set (see {@link BiosynfoniVariant}), so
     * that the custom pattern set initially corresponds to the variant pattern set.
     */
    public static final boolean VARIANT_SETTING_DEFAULT = true;
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * Applied variants for which the key data below were created; null if they have not been created yet. The key data
     * only depend on the applied variants, so they are created once and reused as long as the variants do not change.
     */
    private Set<BiosynfoniVariant> cachedVariants;

    /**
     * SMARTS codes of the keys for the cached variants, in key order.
     */
    private String[] cachedSmartsCodes;

    /**
     * Compiled SMARTS patterns of the keys for the cached variants, in key order.
     */
    private SmartsPattern[] cachedSmartsPatterns;

    /**
     * Annotated representative structures of the keys for the cached variants, in key order; returned as one clone per
     * match.
     */
    private IAtomContainer[] cachedRepresentatives;

    /**
     * Atom expressions of the keys whose SMARTS pattern consists of a single atom for the cached variants, in key
     * order, used to determine the environment of their matches in the matched substructures (see
     * {@link #getSingleAtomEnvironment(Expr, int, IAtomContainer)}); null for keys whose SMARTS pattern consists of
     * more than one atom.
     */
    private Expr[] cachedSingleAtomExpressions;
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Compiled SMARTS patterns by their SMARTS code, created when a SMARTS pattern is used for the first time. Every
     * fragmenter instance holds its own patterns so that instances used in parallel do not share them. The automatic
     * preparation of the target molecule is disabled because it is prepared once per molecule before all patterns are
     * matched.
     */
    private final Map<String, SmartsPattern> smartsPatternCache;

    /**
     * Compiled patterns of the recursive SMARTS expressions of single-atom SMARTS patterns, by their sub-query, created
     * when an expression is used for the first time; used to determine the environment of single-atom matches in the
     * matched substructures. Every fragmenter instance holds its own patterns, like {@link #smartsPatternCache}.
     */
    private final Map<IAtomContainer, Pattern> environmentPatternCache;

    /**
     * Settings of the variants, which are applied in the custom pattern set if turned on (see
     * {@link BiosynfoniVariant}).
     */
    private final Map<BiosynfoniVariant, SimpleBooleanProperty> variantSettings;

    /**
     * SMILES parser for creating the representative fragment structures.
     */
    private final SmilesParser smilesParser;

    //note: since Java 21, the javadoc build complains about "double comments" when there is a comment
    // for the get() method of the property and the private property itself as well
    private final SimpleIDisplayEnumConstantProperty patternSetSetting;

    private final SimpleBooleanProperty intraSubstructureOverlapFilterSetting;

    private final SimpleBooleanProperty interSubstructureOverlapFilterSetting;

    private final SimpleIDisplayEnumConstantProperty returnedFragmentsSetting;

    /**
     * Aromaticity model that is applied to the matched substructure fragments, corresponding to the Daylight model
     * used for SMARTS matching.
     */
    private final Aromaticity fragmentAromaticityModel;

    /**
     * All settings of this fragmenter, encapsulated in JavaFX properties for binding in GUI.
     */
    private final List<Property<?>> settings;

    /**
     * Map to store pairs of {@literal <setting name, tooltip text>}.
     */
    private final HashMap<String, String> settingNameToTooltipTextMap;

    /**
     * Map to store pairs of {@literal <setting name, display name>}.
     */
    private final HashMap<String, String> settingNameToDisplayNameMap;

    /**
     * Map to store pairs of {@literal <setting name, visibility condition>}; the variant settings are only shown
     * while the custom pattern set is chosen because they are not regarded with the other pattern sets.
     */
    private final HashMap<String, ObservableBooleanValue> settingNameToVisibilityConditionMap;

    /**
     * Map to store pairs of {@literal <setting name, section title>}; the variant settings form a section of their own
     * that starts with the setting of the first variant.
     */
    private final Map<String, String> settingNameToSectionTitleMap;

    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(BiosynfoniFragmenter.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, all settings are initialised with their default values as declared in the respective public
     * constants.
     */
    public BiosynfoniFragmenter() {
        this.smartsPatternCache = new HashMap<>();
        this.environmentPatternCache = new IdentityHashMap<>();
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.fragmentAromaticityModel = new Aromaticity(Aromaticity.Model.Daylight,
                Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
        int tmpNumberOfSettings = 4 + BiosynfoniVariant.values().length;
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        int tmpInitialCapacityForSettingMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);

        this.patternSetSetting = this.createEnumSetting("Pattern set setting",
                BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT, BiosynfoniPatternSetOption.class);
        this.settings.add(this.patternSetSetting);
        this.settingNameToTooltipTextMap.put(this.patternSetSetting.getName(),
                Message.get("BiosynfoniFragmenter.patternSetSetting.tooltip"));
        this.settingNameToDisplayNameMap.put(this.patternSetSetting.getName(),
                Message.get("BiosynfoniFragmenter.patternSetSetting.displayName"));

        this.intraSubstructureOverlapFilterSetting = new SimpleBooleanProperty(this,
                "Intra-substructure overlap filter setting",
                BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        this.settings.add(this.intraSubstructureOverlapFilterSetting);
        this.settingNameToTooltipTextMap.put(this.intraSubstructureOverlapFilterSetting.getName(),
                Message.get("BiosynfoniFragmenter.intraSubstructureOverlapFilterSetting.tooltip"));
        this.settingNameToDisplayNameMap.put(this.intraSubstructureOverlapFilterSetting.getName(),
                Message.get("BiosynfoniFragmenter.intraSubstructureOverlapFilterSetting.displayName"));

        this.interSubstructureOverlapFilterSetting = new SimpleBooleanProperty(this,
                "Inter-substructure overlap filter setting",
                BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        this.settings.add(this.interSubstructureOverlapFilterSetting);
        this.settingNameToTooltipTextMap.put(this.interSubstructureOverlapFilterSetting.getName(),
                Message.get("BiosynfoniFragmenter.interSubstructureOverlapFilterSetting.tooltip"));
        this.settingNameToDisplayNameMap.put(this.interSubstructureOverlapFilterSetting.getName(),
                Message.get("BiosynfoniFragmenter.interSubstructureOverlapFilterSetting.displayName"));

        this.returnedFragmentsSetting = this.createEnumSetting("Returned fragments setting",
                BiosynfoniFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, BiosynfoniReturnedFragmentsOption.class);
        this.settings.add(this.returnedFragmentsSetting);
        this.settingNameToTooltipTextMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("BiosynfoniFragmenter.returnedFragmentsSetting.tooltip"));
        this.settingNameToDisplayNameMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("BiosynfoniFragmenter.returnedFragmentsSetting.displayName"));

        this.variantSettings = new EnumMap<>(BiosynfoniVariant.class);
        this.settingNameToVisibilityConditionMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                BiosynfoniVariant.values().length, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR),
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        ObservableBooleanValue tmpIsCustomPatternSetChosen = this.patternSetSetting.isEqualTo(BiosynfoniPatternSetOption.CUSTOM);
        for (BiosynfoniVariant tmpVariant : BiosynfoniVariant.values()) {
            SimpleBooleanProperty tmpVariantSetting = new SimpleBooleanProperty(this, tmpVariant.getSettingName(),
                    BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT);
            this.variantSettings.put(tmpVariant, tmpVariantSetting);
            this.settings.add(tmpVariantSetting);
            this.settingNameToTooltipTextMap.put(tmpVariantSetting.getName(), tmpVariant.getTooltipText());
            this.settingNameToDisplayNameMap.put(tmpVariantSetting.getName(), tmpVariant.getDisplayName());
            this.settingNameToVisibilityConditionMap.put(tmpVariantSetting.getName(), tmpIsCustomPatternSetChosen);
        }
        this.settingNameToSectionTitleMap = Map.of(BiosynfoniVariant.values()[0].getSettingName(),
                Message.get("BiosynfoniFragmenter.variantSettings.sectionTitle"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the currently set option of the pattern set setting.
     *
     * @return enum constant of the set option
     */
    public BiosynfoniPatternSetOption getPatternSetSetting() {
        return (BiosynfoniPatternSetOption) this.patternSetSetting.get();
    }

    /**
     * Returns the property object of the pattern set setting that can be used to configure this setting.
     *
     * @return property object of the pattern set setting
     */
    public SimpleIDisplayEnumConstantProperty patternSetSettingProperty() {
        return this.patternSetSetting;
    }

    /**
     * Returns the current state of the setting of the given variant, which is only regarded in the custom pattern
     * set.
     *
     * @param aVariant the variant
     * @return true if the variant is applied in the custom pattern set
     * @throws NullPointerException if the given variant is null
     */
    public boolean getVariantSetting(BiosynfoniVariant aVariant) throws NullPointerException {
        Objects.requireNonNull(aVariant, "Given variant is null.");
        return this.variantSettings.get(aVariant).get();
    }

    /**
     * Returns the property object of the setting of the given variant that can be used to configure this setting.
     *
     * @param aVariant the variant
     * @return property object of the setting of the variant
     * @throws NullPointerException if the given variant is null
     */
    public SimpleBooleanProperty variantSettingProperty(BiosynfoniVariant aVariant) throws NullPointerException {
        Objects.requireNonNull(aVariant, "Given variant is null.");
        return this.variantSettings.get(aVariant);
    }

    /**
     * Returns the variants that are applied according to the pattern set setting: none for the original pattern
     * set, all for the variant pattern set, and those whose setting is turned on for the custom pattern set.
     *
     * @return the applied variants, a new set
     */
    public Set<BiosynfoniVariant> getAppliedVariants() {
        BiosynfoniPatternSetOption tmpPatternSet = (BiosynfoniPatternSetOption) this.patternSetSetting.get();
        Set<BiosynfoniVariant> tmpAppliedVariants = EnumSet.noneOf(BiosynfoniVariant.class);
        for (BiosynfoniVariant tmpVariant : BiosynfoniVariant.values()) {
            if (tmpPatternSet == BiosynfoniPatternSetOption.VARIANT
                    || (tmpPatternSet == BiosynfoniPatternSetOption.CUSTOM && this.variantSettings.get(tmpVariant).get())) {
                tmpAppliedVariants.add(tmpVariant);
            }
        }
        return tmpAppliedVariants;
    }

    /**
     * Returns the current state of the intra-substructure overlap filter setting.
     *
     * @return true if multiple matches of the same SMARTS pattern are prevented from sharing atoms
     */
    public boolean getIntraSubstructureOverlapFilterSetting() {
        return this.intraSubstructureOverlapFilterSetting.get();
    }

    /**
     * Returns the property object of the intra-substructure overlap filter setting that can be used to configure this
     * setting.
     *
     * @return property object of the intra-substructure overlap filter setting
     */
    public SimpleBooleanProperty intraSubstructureOverlapFilterSettingProperty() {
        return this.intraSubstructureOverlapFilterSetting;
    }

    /**
     * Returns the current state of the inter-substructure overlap filter setting.
     *
     * @return true if matches are prevented from sharing atoms with accepted matches of earlier SMARTS patterns
     */
    public boolean getInterSubstructureOverlapFilterSetting() {
        return this.interSubstructureOverlapFilterSetting.get();
    }

    /**
     * Returns the property object of the inter-substructure overlap filter setting that can be used to configure this
     * setting.
     *
     * @return property object of the inter-substructure overlap filter setting
     */
    public SimpleBooleanProperty interSubstructureOverlapFilterSettingProperty() {
        return this.interSubstructureOverlapFilterSetting;
    }

    /**
     * Returns the currently set option of the returned fragments setting.
     *
     * @return enum constant of the set option
     */
    public BiosynfoniReturnedFragmentsOption getReturnedFragmentsSetting() {
        return (BiosynfoniReturnedFragmentsOption) this.returnedFragmentsSetting.get();
    }

    /**
     * Returns the property object of the returned fragments setting that can be used to configure this setting.
     *
     * @return property object of the returned fragments setting
     */
    public SimpleIDisplayEnumConstantProperty returnedFragmentsSettingProperty() {
        return this.returnedFragmentsSetting;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set">
    /**
     * Sets the pattern set setting, defining whether the original SMARTS patterns of Biosynfoni version "full_1103",
     * all pattern variants, or a custom selection of the pattern variants (see
     * {@link #setVariantSetting(BiosynfoniVariant, boolean)}) are used.
     *
     * @param anOption a constant from the {@link BiosynfoniPatternSetOption} enum
     * @throws NullPointerException if the given option is null
     */
    public void setPatternSetSetting(BiosynfoniPatternSetOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given option is null.");
        this.patternSetSetting.set(anOption);
    }

    /**
     * Sets the setting of the given variant, which is only regarded in the custom pattern set.
     *
     * @param aVariant the variant
     * @param aBoolean true if the variant should be applied in the custom pattern set
     * @throws NullPointerException if the given variant is null
     */
    public void setVariantSetting(BiosynfoniVariant aVariant, boolean aBoolean) throws NullPointerException {
        Objects.requireNonNull(aVariant, "Given variant is null.");
        this.variantSettings.get(aVariant).set(aBoolean);
    }

    /**
     * Sets the intra-substructure overlap filter setting. If true, multiple matches of the same SMARTS pattern are
     * prevented from sharing atoms; matches are accepted in the order of the canonical ranks of their atoms.
     *
     * @param aBoolean true if the intra-substructure overlap filter should be applied
     */
    public void setIntraSubstructureOverlapFilterSetting(boolean aBoolean) {
        this.intraSubstructureOverlapFilterSetting.set(aBoolean);
    }

    /**
     * Sets the inter-substructure overlap filter setting. If true, matches of a SMARTS pattern are prevented from
     * sharing atoms with accepted matches of SMARTS patterns that come earlier in the key order.
     *
     * @param aBoolean true if the inter-substructure overlap filter should be applied
     */
    public void setInterSubstructureOverlapFilterSetting(boolean aBoolean) {
        this.interSubstructureOverlapFilterSetting.set(aBoolean);
    }

    /**
     * Sets the returned fragments setting, defining whether the matches of a key are returned as the stored
     * representative structure of the key or as the actually matched substructures.
     *
     * @param anOption a constant from the {@link BiosynfoniReturnedFragmentsOption} enum
     * @throws NullPointerException if the given option is null
     */
    public void setReturnedFragmentsSetting(BiosynfoniReturnedFragmentsOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given option is null.");
        this.returnedFragmentsSetting.set(anOption);
    }
    //</editor-fold>
    //
    //<editor-fold desc="IMoleculeFragmenter methods">
    //without the empty line, the code folding does not work properly here...

    @Override
    public List<Property<?>> settingsProperties() {
        return this.settings;
    }

    @Override
    public Map<String, String> getSettingNameToTooltipTextMap() {
        return this.settingNameToTooltipTextMap;
    }

    @Override
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return this.settingNameToDisplayNameMap;
    }

    @Override
    public Map<String, ObservableBooleanValue> getSettingNameToVisibilityConditionMap() {
        return this.settingNameToVisibilityConditionMap;
    }

    @Override
    public Map<String, String> getSettingNameToSectionTitleMap() {
        return this.settingNameToSectionTitleMap;
    }

    @Override
    public String getFragmentationAlgorithmName() {
        return BiosynfoniFragmenter.ALGORITHM_NAME;
    }

    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("BiosynfoniFragmenter.displayName");
    }

    @Override
    public IMoleculeFragmenter copy() {
        BiosynfoniFragmenter tmpCopy = new BiosynfoniFragmenter();
        tmpCopy.setPatternSetSetting((BiosynfoniPatternSetOption) this.patternSetSetting.get());
        for (BiosynfoniVariant tmpVariant : BiosynfoniVariant.values()) {
            tmpCopy.setVariantSetting(tmpVariant, this.variantSettings.get(tmpVariant).get());
        }
        tmpCopy.setIntraSubstructureOverlapFilterSetting(this.intraSubstructureOverlapFilterSetting.get());
        tmpCopy.setInterSubstructureOverlapFilterSetting(this.interSubstructureOverlapFilterSetting.get());
        tmpCopy.setReturnedFragmentsSetting((BiosynfoniReturnedFragmentsOption) this.returnedFragmentsSetting.get());
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {
        this.patternSetSetting.set(BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT);
        for (SimpleBooleanProperty tmpVariantSetting : this.variantSettings.values()) {
            tmpVariantSetting.set(BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT);
        }
        this.intraSubstructureOverlapFilterSetting.set(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        this.interSubstructureOverlapFilterSetting.set(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        this.returnedFragmentsSetting.set(BiosynfoniFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT);
    }

    @Override
    public List<IAtomContainer> fragmentMolecule(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Parameter tests">
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpCanBeFragmented = this.canBeFragmented(aMolecule);
        if (!tmpCanBeFragmented) {
            throw new IllegalArgumentException("Given molecule cannot be fragmented but should be filtered or preprocessed first.");
        }
        //</editor-fold>
        List<IAtomContainer> tmpFragments = new ArrayList<>();
        boolean tmpAreMatchedSubstructuresReturned = this.returnedFragmentsSetting.get().equals(
                BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Set<BiosynfoniVariant> tmpVariants = this.getAppliedVariants();
        try {
            this.updateKeyDataCache(tmpVariants);
            IAtomContainer tmpPreparedMolecule = BiosynfoniFragmenter.prepareMolecule(aMolecule);
            int[] tmpCanonicalRanks = this.intraSubstructureOverlapFilterSetting.get()
                    ? BiosynfoniFragmenter.getCanonicalRanks(tmpPreparedMolecule) : null;
            //ring keys are matched by ring perception if this variant is applied, see BiosynfoniKey.getRingSize()
            Map<Integer, List<int[]>> tmpCarbocyclesBySize = tmpVariants.contains(BiosynfoniVariant.RING_PERCEPTION)
                    ? BiosynfoniFragmenter.getRelevantCarbocyclesBySize(tmpPreparedMolecule) : null;
            //atoms of accepted matches of previous keys, only used by the inter-substructure overlap filter
            boolean[] tmpIsAtomBlockedByPreviousKeys = new boolean[tmpPreparedMolecule.getAtomCount()];
            BiosynfoniKey[] tmpKeys = BiosynfoniKey.values();
            for (int i = 0; i < tmpKeys.length; i++) {
                List<int[]> tmpMatches = this.getFilteredMatches(
                        BiosynfoniFragmenter.getMatches(tmpKeys[i], this.cachedSmartsPatterns[i], tmpPreparedMolecule,
                                tmpCarbocyclesBySize),
                        tmpCanonicalRanks, tmpIsAtomBlockedByPreviousKeys);
                if (tmpAreMatchedSubstructuresReturned) {
                    for (int[] tmpMatch : tmpMatches) {
                        //a single matched atom is returned together with its environment that the pattern checks, see
                        // getSingleAtomEnvironment(); the overlap filters above only regard the matched atoms
                        int[][] tmpFragmentAtoms = Objects.isNull(this.cachedSingleAtomExpressions[i])
                                ? new int[][] {tmpMatch, new int[0]}
                                : this.getSingleAtomEnvironment(this.cachedSingleAtomExpressions[i], tmpMatch[0],
                                        tmpPreparedMolecule);
                        tmpFragments.add(this.createMatchedSubstructureFragment(tmpPreparedMolecule, tmpFragmentAtoms[0],
                                tmpFragmentAtoms[1], tmpKeys[i], this.cachedSmartsCodes[i], tmpVariants));
                    }
                } else {
                    //all matches of one key are represented by the same structure and thereby collapse into one fragment;
                    // every match gets its own clone, so that the returned fragments can be altered independently
                    for (int j = 0; j < tmpMatches.size(); j++) {
                        tmpFragments.add(this.cachedRepresentatives[i].clone());
                    }
                }
            }
        } catch (Exception anException) {
            throw new IllegalArgumentException("An error occurred during fragmentation: " + anException.toString()
                    + " Molecule Name: " + aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException);
        }
        return tmpFragments;
    }

    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        return (Objects.isNull(aMolecule) || aMolecule.isEmpty());
    }

    @Override
    public boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return false;
    }

    @Override
    public boolean canBeFragmented(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered = this.shouldBeFiltered(aMolecule);
        boolean tmpShouldBePreprocessed = this.shouldBePreprocessed(aMolecule);
        return !(tmpShouldBeFiltered || tmpShouldBePreprocessed);
    }

    @Override
    public boolean providesDisplayLabels() {
        return true;
    }

    @Override
    public IAtomContainer applyPreprocessing(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered = this.shouldBeFiltered(aMolecule);
        if (tmpShouldBeFiltered) {
            throw new IllegalArgumentException("The given molecule cannot be preprocessed but should be filtered.");
        }
        return aMolecule.clone();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Returns the compiled SMARTS pattern of the given SMARTS code, from the cache or newly created. The automatic
     * preparation of the target molecule is disabled because it is prepared once per molecule in
     * {@link #prepareMolecule(IAtomContainer)}.
     *
     * @param aSmarts the SMARTS code
     * @return compiled SMARTS pattern
     */
    private SmartsPattern getSmartsPattern(String aSmarts) {
        return this.smartsPatternCache.computeIfAbsent(aSmarts, aSmartsCode -> SmartsPattern.create(aSmartsCode).setPrepare(false));
    }

    /**
     * Returns the atoms of a matched substructure of a key whose SMARTS pattern consists of a single atom: the matched
     * atom and its environment that the pattern checks. A single matched atom alone, e.g. "N" or "O", would not show
     * where it was matched. The environment consists of the atoms of all matches of the recursive SMARTS expressions
     * ("$(...)") of the pattern that start at the matched atom, including the matches of recursive expressions nested
     * in them, e.g. the backbone and the side chain of an amino acid for the amino acid keys or the carbon or nitrogen
     * atom of a hydroxyl group; negated recursive expressions ("!$(...)") are not regarded, since they describe what
     * must not be present. Atoms that a recursive expression only requires as any atom ("*"), e.g. the side chain of
     * the non-standard amino acid variant, become attachment points, i.e. pseudo atoms ("*"), because the pattern does
     * not specify them; e.g. the side chain of phenylglycine would otherwise appear as a methyl group, as in alanine.
     * If the pattern has no recursive expression that adds an atom, e.g. the ether, terminal
     * nitrogen, epoxide, methyl, and halogen patterns, whose conditions refer to the neighbours or rings of the atom
     * via primitives like D, x, or r, the direct neighbours of the matched atom are the environment: neighbours in a
     * common ring with the matched atom as atoms, e.g. the ring of an epoxide, all other neighbours as attachment
     * points, i.e. pseudo atoms ("*"), because the pattern only requires their presence, not their element, e.g.
     * "*C" for a methyl group instead of "CC", which would look like an ethyl unit.
     *
     * @param anExpression atom expression of the single-atom SMARTS pattern
     * @param anAtomIndex index of the matched atom
     * @param aMolecule the molecule the match was found in, prepared for SMARTS matching
     * @return two arrays of atom indices: the matched atom (first) and the atoms of its environment (in ascending
     *         order), and the neighbours that become attachment points (in ascending order, possibly empty)
     */
    private int[][] getSingleAtomEnvironment(Expr anExpression, int anAtomIndex, IAtomContainer aMolecule) {
        Set<Integer> tmpEnvironment = new TreeSet<>();
        Set<Integer> tmpAttachmentPoints = new TreeSet<>();
        this.collectRecursiveEnvironment(anExpression, anAtomIndex, aMolecule, tmpEnvironment, tmpAttachmentPoints);
        tmpEnvironment.remove(anAtomIndex);
        //an atom that is specified by another expression is shown as atom
        tmpAttachmentPoints.removeAll(tmpEnvironment);
        tmpAttachmentPoints.remove(anAtomIndex);
        if (tmpEnvironment.isEmpty() && tmpAttachmentPoints.isEmpty()) {
            IAtom tmpMatchedAtom = aMolecule.getAtom(anAtomIndex);
            for (IBond tmpBond : aMolecule.getConnectedBondsList(tmpMatchedAtom)) {
                int tmpNeighbourIndex = tmpBond.getOther(tmpMatchedAtom).getIndex();
                (tmpBond.isInRing() ? tmpEnvironment : tmpAttachmentPoints).add(tmpNeighbourIndex);
            }
        }
        int[] tmpAtoms = new int[tmpEnvironment.size() + 1];
        tmpAtoms[0] = anAtomIndex;
        int tmpPosition = 1;
        for (int tmpAtomIndex : tmpEnvironment) {
            tmpAtoms[tmpPosition++] = tmpAtomIndex;
        }
        return new int[][] {tmpAtoms, tmpAttachmentPoints.stream().mapToInt(Integer::intValue).toArray()};
    }

    /**
     * Adds the atoms of all matches of the non-negated recursive SMARTS expressions in the given atom expression that
     * start at the given atom to the given sets, including the atoms of the recursive expressions nested in these
     * expressions, see {@link #getSingleAtomEnvironment(Expr, int, IAtomContainer)}: atoms that an expression only
     * requires as any atom ("*") to the attachment points, all other atoms to the environment. Alternatives of an OR
     * expression that do not match do not add atoms.
     *
     * @param anExpression atom expression
     * @param anAtomIndex index of the atom the expression refers to
     * @param aMolecule the molecule, prepared for SMARTS matching
     * @param anEnvironment set to add the indices of the specified atoms to
     * @param anAttachmentPointSet set to add the indices of the atoms required as any atom to
     */
    private void collectRecursiveEnvironment(Expr anExpression, int anAtomIndex, IAtomContainer aMolecule,
                                             Set<Integer> anEnvironment, Set<Integer> anAttachmentPointSet) {
        switch (anExpression.type()) {
            case AND, OR -> {
                this.collectRecursiveEnvironment(anExpression.left(), anAtomIndex, aMolecule, anEnvironment,
                        anAttachmentPointSet);
                this.collectRecursiveEnvironment(anExpression.right(), anAtomIndex, aMolecule, anEnvironment,
                        anAttachmentPointSet);
            }
            case RECURSIVE -> {
                IAtomContainer tmpSubquery = anExpression.subquery();
                Pattern tmpPattern = this.environmentPatternCache.computeIfAbsent(tmpSubquery, Pattern::findSubstructure);
                for (int[] tmpMapping : tmpPattern.matchAll(aMolecule)) {
                    //the first atom of a recursive expression is the atom it refers to
                    if (tmpMapping[0] != anAtomIndex) {
                        continue;
                    }
                    for (int j = 0; j < tmpMapping.length; j++) {
                        Expr tmpAtomExpression = ((QueryAtom) AtomRef.deref(tmpSubquery.getAtom(j))).getExpression();
                        //the first atom is the one the expression refers to, any further atom given as "*" is not
                        // specified by the pattern
                        if (j > 0 && tmpAtomExpression.type() == Expr.Type.TRUE) {
                            anAttachmentPointSet.add(tmpMapping[j]);
                            continue;
                        }
                        anEnvironment.add(tmpMapping[j]);
                        this.collectRecursiveEnvironment(tmpAtomExpression, tmpMapping[j], aMolecule, anEnvironment,
                                anAttachmentPointSet);
                    }
                }
            }
            //negated expressions describe what must not be present, all other types are primitives without atoms
            default -> {
            }
        }
    }

    /**
     * Creates the SMARTS codes, compiled SMARTS patterns, and annotated representative structures of all keys for the
     * given applied variants, unless they have already been created for these variants.
     *
     * @param aVariantSet the applied variants
     * @throws CDKException if a stored representative SMILES code cannot be parsed
     */
    private void updateKeyDataCache(Set<BiosynfoniVariant> aVariantSet) throws CDKException {
        if (aVariantSet.equals(this.cachedVariants)) {
            return;
        }
        BiosynfoniKey[] tmpKeys = BiosynfoniKey.values();
        String[] tmpSmartsCodes = new String[tmpKeys.length];
        SmartsPattern[] tmpSmartsPatterns = new SmartsPattern[tmpKeys.length];
        IAtomContainer[] tmpRepresentatives = new IAtomContainer[tmpKeys.length];
        Expr[] tmpSingleAtomExpressions = new Expr[tmpKeys.length];
        for (int i = 0; i < tmpKeys.length; i++) {
            tmpSmartsCodes[i] = tmpKeys[i].getSmarts(aVariantSet);
            tmpSmartsPatterns[i] = this.getSmartsPattern(tmpSmartsCodes[i]);
            tmpRepresentatives[i] = this.smilesParser.parseSmiles(tmpKeys[i].getRepresentativeSmiles(aVariantSet));
            BiosynfoniFragmenter.annotateFragment(tmpRepresentatives[i], tmpKeys[i], tmpSmartsCodes[i], aVariantSet);
            IAtomContainer tmpQuery = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
            if (!Smarts.parse(tmpQuery, tmpSmartsCodes[i])) {
                throw new CDKException("SMARTS pattern cannot be parsed: " + tmpSmartsCodes[i]);
            }
            if (tmpQuery.getAtomCount() == 1) {
                tmpSingleAtomExpressions[i] = ((QueryAtom) AtomRef.deref(tmpQuery.getAtom(0))).getExpression();
            }
        }
        this.cachedSmartsCodes = tmpSmartsCodes;
        this.cachedSmartsPatterns = tmpSmartsPatterns;
        this.cachedRepresentatives = tmpRepresentatives;
        this.cachedSingleAtomExpressions = tmpSingleAtomExpressions;
        //set last, so that an exception above leaves no partially updated cache behind
        this.cachedVariants = EnumSet.copyOf(aVariantSet);
    }

    /**
     * Creates a setting property for the constants of the given enum class whose setter shows an alert in the GUI and
     * re-throws the exception if an illegal value is given, which properly resets the binding.
     *
     * @param aName name of the setting
     * @param aDefaultOption default option of the setting
     * @param anEnumClass the enum class of the options
     * @return the setting property
     */
    private SimpleIDisplayEnumConstantProperty createEnumSetting(String aName, IDisplayEnum aDefaultOption,
                                                                 Class<? extends IDisplayEnum> anEnumClass) {
        return new SimpleIDisplayEnumConstantProperty(this, aName, aDefaultOption, anEnumClass) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    BiosynfoniFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
    }

    /**
     * Returns a clone of the given molecule that is prepared for matching all SMARTS patterns, i.e. ring atoms and
     * bonds are marked and the Daylight aromaticity model is applied once (see
     * {@link SmartsPattern#prepare(IAtomContainer)}), since the automatic preparation of the patterns is disabled. If
     * the clone has unset bond orders, it is kekulized beforehand in both returned fragments modes, because the
     * aromaticity perception of the SMARTS preparation and the matched substructures need explicit bond orders, and
     * the counts must not depend on the returned fragments setting. If no Kekulé structure can be assigned, the clone
     * is used with its unset bond orders, which are then treated as single bonds in the matched substructures.
     *
     * @param aMolecule the molecule to fragment, it is not altered
     * @return prepared clone of the given molecule
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    private static IAtomContainer prepareMolecule(IAtomContainer aMolecule) throws CloneNotSupportedException {
        //SMARTS preparation alters the aromaticity and ring flags of the atom container, so work on a clone
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        if (BiosynfoniFragmenter.hasUnsetBondOrders(tmpMoleculeClone)) {
            try {
                Kekulization.kekulize(tmpMoleculeClone);
            } catch (CDKException | IllegalArgumentException anException) {
                //IllegalArgumentException e.g. for unset implicit hydrogen counts
                BiosynfoniFragmenter.LOGGER.log(Level.FINE, "Molecule cannot be kekulized, unset bond orders are"
                        + " treated as single bonds in the matched substructures.", anException);
            }
        }
        SmartsPattern.prepare(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Checks whether the given molecule has bonds with unset bond order.
     *
     * @param aMolecule the molecule to check
     * @return true if at least one bond order is null or unset
     */
    private static boolean hasUnsetBondOrders(IAtomContainer aMolecule) {
        for (IBond tmpBond : aMolecule.bonds()) {
            if (Objects.isNull(tmpBond.getOrder()) || tmpBond.getOrder() == IBond.Order.UNSET) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the unique matches of the given key in the prepared molecule, before the overlap filters are applied.
     * Ring keys are matched by the given relevant carbocycles if these are given, i.e. if
     * {@link BiosynfoniVariant#RING_PERCEPTION} is applied (see {@link BiosynfoniKey#getRingSize()}); all other keys
     * are matched by the given SMARTS pattern.
     *
     * @param aKey the Biosynfoni key
     * @param aPattern SMARTS pattern of the key for the applied variants
     * @param aPreparedMolecule molecule prepared for SMARTS matching, see {@link #prepareMolecule(IAtomContainer)}
     * @param aCarbocyclesBySize relevant cycles consisting only of carbon atoms, by ring size; null if ring keys are
     *                           matched by their SMARTS pattern
     * @return matches, each given as atom indices
     */
    private static List<int[]> getMatches(BiosynfoniKey aKey, SmartsPattern aPattern, IAtomContainer aPreparedMolecule,
                                          Map<Integer, List<int[]>> aCarbocyclesBySize) {
        if (aKey.getRingSize() > 0 && !Objects.isNull(aCarbocyclesBySize)) {
            return aCarbocyclesBySize.getOrDefault(aKey.getRingSize(), Collections.emptyList());
        }
        return Arrays.asList(aPattern.matchAll(aPreparedMolecule).uniqueAtoms().toArray());
    }

    /**
     * Returns the relevant cycles of the given molecule that consist only of carbon atoms, grouped by their size (see
     * {@link BiosynfoniKey#getRingSize()}).
     *
     * @param aMolecule the molecule
     * @return atom indices of the relevant carbocycles by ring size
     */
    private static Map<Integer, List<int[]>> getRelevantCarbocyclesBySize(IAtomContainer aMolecule) {
        Map<Integer, List<int[]>> tmpCarbocyclesBySize = new HashMap<>();
        for (int[] tmpPath : Cycles.relevant(aMolecule).paths()) {
            //the path is closed, i.e. its last atom index repeats the first one
            int[] tmpCycle = Arrays.copyOf(tmpPath, tmpPath.length - 1);
            boolean tmpIsCarbocycle = true;
            for (int tmpAtomIndex : tmpCycle) {
                if (!Objects.equals(aMolecule.getAtom(tmpAtomIndex).getAtomicNumber(), IElement.C)) {
                    tmpIsCarbocycle = false;
                    break;
                }
            }
            if (tmpIsCarbocycle) {
                tmpCarbocyclesBySize.computeIfAbsent(tmpCycle.length, aSize -> new ArrayList<>()).add(tmpCycle);
            }
        }
        return tmpCarbocyclesBySize;
    }

    /**
     * Filters the given matches of a key by the overlap filters that are turned on. If the inter-substructure overlap
     * filter is turned on, the atoms of the returned matches are additionally flagged as blocked for the following
     * keys.
     *
     * @param aMatchesList unique matches of the current key, each given as atom indices; the list is not altered
     * @param aCanonicalRanks canonical rank of every atom, indexed by atom index; only used and may be null if the
     *                        intra-substructure overlap filter is turned off
     * @param anIsAtomBlockedArray flags of atoms used by accepted matches of previous keys, indexed by atom index; only
     *                             used by the inter-substructure overlap filter
     * @return filtered matches, each given as atom indices
     */
    private List<int[]> getFilteredMatches(List<int[]> aMatchesList, int[] aCanonicalRanks, boolean[] anIsAtomBlockedArray) {
        List<int[]> tmpMatches = aMatchesList;
        if (this.intraSubstructureOverlapFilterSetting.get()) {
            tmpMatches = BiosynfoniFragmenter.applyIntraSubstructureOverlapFilter(tmpMatches, aCanonicalRanks);
        }
        if (this.interSubstructureOverlapFilterSetting.get()) {
            tmpMatches = BiosynfoniFragmenter.applyInterSubstructureOverlapFilter(tmpMatches, anIsAtomBlockedArray);
            for (int[] tmpMatch : tmpMatches) {
                for (int tmpAtomIndex : tmpMatch) {
                    anIsAtomBlockedArray[tmpAtomIndex] = true;
                }
            }
        }
        return tmpMatches;
    }

    /**
     * Creates a fragment from the given match that consists of the matched atoms and all bonds between them. Free
     * valences resulting from bonds to non-matched atoms are saturated with implicit hydrogen atoms, and aromaticity
     * is re-detected on the fragment. A cut bond of an aromatic system counts as a single bond and a Kekulé structure
     * is assigned to the matched part of the aromatic system, so that completely matched aromatic rings remain
     * aromatic, e.g. the pyrrole ring of indole. A matched atom whose aromatic bonds are all cut is not part of the
     * matched aromatic system, so that e.g. the bond between the two rings of biphenyl remains a single bond. If no
     * Kekulé structure can be assigned, e.g. to three atoms of a benzene ring, the fragment is built from the bond
     * orders of the given molecule instead. Ring atoms and bonds of the fragment are marked, and the fragment is
     * annotated with the key label and SMARTS pattern.
     *
     * @param aMolecule the molecule the match was found in, prepared for SMARTS matching; it should have explicit
     *                  (Kekulé) bond orders, unset bond orders are treated as single bonds
     * @param aMatch atom indices of the match, for single-atom patterns including the environment of the matched atom
     *               (see {@link #getSingleAtomEnvironment(Expr, int, IAtomContainer)})
     * @param anAttachmentPointArray atom indices of neighbours that are added as attachment points, i.e. pseudo atoms
     *                               ("*") bonded to the matched atom; empty except for single-atom patterns
     * @param aKey the Biosynfoni key that produced the match
     * @param aSmarts SMARTS pattern of the key for the applied variants
     * @param aVariantSet the applied variants
     * @return new atom container of the matched substructure
     * @throws CloneNotSupportedException if an atom cannot be cloned
     * @throws CDKException if aromaticity detection on the fragment fails
     */
    private IAtomContainer createMatchedSubstructureFragment(IAtomContainer aMolecule, int[] aMatch,
                                                             int[] anAttachmentPointArray, BiosynfoniKey aKey,
                                                             String aSmarts, Set<BiosynfoniVariant> aVariantSet)
            throws CloneNotSupportedException, CDKException {
        int[] tmpAtoms = aMatch;
        if (anAttachmentPointArray.length > 0) {
            tmpAtoms = Arrays.copyOf(aMatch, aMatch.length + anAttachmentPointArray.length);
            System.arraycopy(anAttachmentPointArray, 0, tmpAtoms, aMatch.length, anAttachmentPointArray.length);
        }
        IAtomContainer tmpFragment = BiosynfoniFragmenter.copyMatchedAtomsAndBonds(aMolecule, tmpAtoms, true);
        try {
            Kekulization.kekulize(tmpFragment);
        } catch (CDKException anException) {
            //the matched part of the aromatic system has no Kekulé structure, the bond orders of the molecule are used
            tmpFragment = BiosynfoniFragmenter.copyMatchedAtomsAndBonds(aMolecule, tmpAtoms, false);
        }
        //the atoms are copied in the given order, so the attachment points are the last atoms of the fragment
        for (int i = aMatch.length; i < tmpAtoms.length; i++) {
            IPseudoAtom tmpAttachmentPoint = tmpFragment.getBuilder().newInstance(IPseudoAtom.class, "*");
            tmpAttachmentPoint.setImplicitHydrogenCount(0);
            tmpFragment.setAtom(i, tmpAttachmentPoint);
        }
        Cycles.markRingAtomsAndBonds(tmpFragment);
        this.fragmentAromaticityModel.apply(tmpFragment);
        BiosynfoniFragmenter.annotateFragment(tmpFragment, aKey, aSmarts, aVariantSet);
        return tmpFragment;
    }

    /**
     * Copies the matched atoms and all bonds between them into a new atom container and saturates the free valences
     * resulting from bonds to non-matched atoms with implicit hydrogen atoms. If the aromaticity is kept, the matched
     * aromatic bonds and their atoms are flagged as aromatic, the bond orders of these bonds are left unset so that a
     * Kekulé structure can be assigned to them afterwards, and every cut aromatic bond of these atoms counts as a
     * single bond. A matched aromatic atom whose aromatic bonds are all cut is not part of the matched aromatic system:
     * it is not flagged and its cut bonds count with the bond orders of the given molecule, so that it is saturated and
     * a single bond between two such atoms, e.g. between the two rings of biphenyl, remains a single bond. Otherwise,
     * all bond orders of the given molecule are used and no aromatic flags are set. Unset bond orders of non-aromatic
     * bonds are treated as single bonds.
     *
     * @param aMolecule the molecule the match was found in, prepared for SMARTS matching
     * @param aMatch atom indices of the match
     * @param isAromaticityKept whether the matched aromatic bonds are flagged and cut aromatic bonds of their atoms
     *                          count as single bonds
     * @return new atom container of the matched atoms and the bonds between them
     * @throws CloneNotSupportedException if an atom cannot be cloned
     */
    private static IAtomContainer copyMatchedAtomsAndBonds(IAtomContainer aMolecule, int[] aMatch, boolean isAromaticityKept)
            throws CloneNotSupportedException {
        IChemObjectBuilder tmpBuilder = aMolecule.getBuilder();
        IAtomContainer tmpFragment = tmpBuilder.newAtomContainer();
        boolean[] tmpIsMatchedAtom = new boolean[aMolecule.getAtomCount()];
        for (int tmpAtomIndex : aMatch) {
            tmpIsMatchedAtom[tmpAtomIndex] = true;
        }
        //only atoms with a matched aromatic bond belong to the matched part of an aromatic system; an aromatic atom whose
        // aromatic bonds are all cut is saturated instead, otherwise the Kekulé structure assignment could turn a
        // single bond between two such atoms into a double bond, e.g. the bond between the two rings of biphenyl
        boolean[] tmpHasMatchedAromaticBond = new boolean[aMolecule.getAtomCount()];
        if (isAromaticityKept) {
            for (IBond tmpBond : aMolecule.bonds()) {
                int tmpBeginIndex = tmpBond.getBegin().getIndex();
                int tmpEndIndex = tmpBond.getEnd().getIndex();
                if (tmpBond.isAromatic() && tmpIsMatchedAtom[tmpBeginIndex] && tmpIsMatchedAtom[tmpEndIndex]) {
                    tmpHasMatchedAromaticBond[tmpBeginIndex] = true;
                    tmpHasMatchedAromaticBond[tmpEndIndex] = true;
                }
            }
        }
        IAtom[] tmpOriginalToFragmentAtom = new IAtom[aMolecule.getAtomCount()];
        for (int tmpAtomIndex : aMatch) {
            IAtom tmpOriginalAtom = aMolecule.getAtom(tmpAtomIndex);
            IAtom tmpFragmentAtom = tmpOriginalAtom.clone();
            int tmpImplicitHydrogenCount = Objects.isNull(tmpOriginalAtom.getImplicitHydrogenCount())
                    ? 0 : tmpOriginalAtom.getImplicitHydrogenCount();
            tmpFragmentAtom.setImplicitHydrogenCount(tmpImplicitHydrogenCount);
            tmpFragmentAtom.setIsAromatic(tmpHasMatchedAromaticBond[tmpAtomIndex]);
            tmpFragmentAtom.setIsInRing(false);
            tmpFragment.addAtom(tmpFragmentAtom);
            //the container may wrap the added atom, so the stored reference has to be taken from it
            tmpOriginalToFragmentAtom[tmpAtomIndex] = tmpFragment.getAtom(tmpFragment.getAtomCount() - 1);
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            IAtom tmpBegin = tmpOriginalToFragmentAtom[tmpBond.getBegin().getIndex()];
            IAtom tmpEnd = tmpOriginalToFragmentAtom[tmpBond.getEnd().getIndex()];
            boolean tmpIsAromaticBondKept = isAromaticityKept && tmpBond.isAromatic();
            IBond.Order tmpBondOrder = Objects.isNull(tmpBond.getOrder()) || tmpBond.getOrder() == IBond.Order.UNSET
                    ? IBond.Order.SINGLE : tmpBond.getOrder();
            if (!Objects.isNull(tmpBegin) && !Objects.isNull(tmpEnd)) {
                IBond tmpFragmentBond = tmpBuilder.newInstance(IBond.class, tmpBegin, tmpEnd,
                        tmpIsAromaticBondKept ? IBond.Order.UNSET : tmpBondOrder);
                tmpFragmentBond.setIsAromatic(tmpIsAromaticBondKept);
                tmpFragment.addBond(tmpFragmentBond);
            } else if (!Objects.isNull(tmpBegin) || !Objects.isNull(tmpEnd)) {
                //bond to a non-matched atom is cut, saturate the free valence with hydrogen atoms
                IAtom tmpMatchedAtom = Objects.isNull(tmpBegin) ? tmpEnd : tmpBegin;
                int tmpMatchedAtomIndex = Objects.isNull(tmpBegin) ? tmpBond.getEnd().getIndex() : tmpBond.getBegin().getIndex();
                boolean tmpCountsAsSingleBond = tmpIsAromaticBondKept && tmpHasMatchedAromaticBond[tmpMatchedAtomIndex];
                tmpMatchedAtom.setImplicitHydrogenCount(tmpMatchedAtom.getImplicitHydrogenCount()
                        + (tmpCountsAsSingleBond ? 1 : tmpBondOrder.numeric()));
            }
        }
        return tmpFragment;
    }

    /**
     * Sets the title and the label, SMARTS, matching method, category, and display label properties of the given
     * fragment according to the given key and applied variants. The display label is the readable display name of
     * the key.
     *
     * @param aFragment the fragment to annotate
     * @param aKey the Biosynfoni key that produced the fragment
     * @param aSmarts SMARTS pattern of the key for the applied variants
     * @param aVariantSet the applied variants
     */
    private static void annotateFragment(IAtomContainer aFragment, BiosynfoniKey aKey, String aSmarts,
                                         Set<BiosynfoniVariant> aVariantSet) {
        aFragment.setTitle(aKey.getLabel());
        aFragment.setProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY, aKey.getLabel());
        aFragment.setProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY, aSmarts);
        aFragment.setProperty(BiosynfoniFragmenter.BIOSYNFONI_MATCHING_METHOD_PROPERTY_KEY,
                aKey.getRingSize() > 0 && aVariantSet.contains(BiosynfoniVariant.RING_PERCEPTION)
                        ? BiosynfoniFragmenter.MATCHING_METHOD_RING_PERCEPTION : BiosynfoniFragmenter.MATCHING_METHOD_SMARTS);
        aFragment.setProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY, aKey.getLabel());
        aFragment.setProperty(IMoleculeFragmenter.FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY, aKey.getDisplayName());
    }

    /**
     * Removes matches that share atoms with previously accepted matches of the same SMARTS pattern. Matches are
     * processed in lexicographic order of their sorted canonical atom ranks so that the result does not depend on the
     * input atom order. The atom indices of the returned matches are not altered.
     *
     * @param aMatchesList matches of one SMARTS pattern
     * @param aCanonicalRanks canonical rank of every atom, indexed by atom index
     * @return atom-disjoint subset of the given matches
     */
    private static List<int[]> applyIntraSubstructureOverlapFilter(List<int[]> aMatchesList, int[] aCanonicalRanks) {
        if (aMatchesList.size() < 2) {
            return aMatchesList;
        }
        int tmpNumberOfMatches = aMatchesList.size();
        int[][] tmpSortedRanksOfMatches = new int[tmpNumberOfMatches][];
        Integer[] tmpOrder = new Integer[tmpNumberOfMatches];
        for (int i = 0; i < tmpNumberOfMatches; i++) {
            int[] tmpMatch = aMatchesList.get(i);
            int[] tmpRanks = new int[tmpMatch.length];
            for (int j = 0; j < tmpMatch.length; j++) {
                tmpRanks[j] = aCanonicalRanks[tmpMatch[j]];
            }
            Arrays.sort(tmpRanks);
            tmpSortedRanksOfMatches[i] = tmpRanks;
            tmpOrder[i] = i;
        }
        Arrays.sort(tmpOrder, (anIndex1, anIndex2) -> Arrays.compare(tmpSortedRanksOfMatches[anIndex1],
                tmpSortedRanksOfMatches[anIndex2]));
        boolean[] tmpIsAtomBlocked = new boolean[aCanonicalRanks.length];
        List<int[]> tmpFilteredMatches = new ArrayList<>(tmpNumberOfMatches);
        for (int tmpIndex : tmpOrder) {
            int[] tmpMatch = aMatchesList.get(tmpIndex);
            if (BiosynfoniFragmenter.hasOverlap(tmpMatch, tmpIsAtomBlocked)) {
                continue;
            }
            tmpFilteredMatches.add(tmpMatch);
            for (int tmpAtomIndex : tmpMatch) {
                tmpIsAtomBlocked[tmpAtomIndex] = true;
            }
        }
        return tmpFilteredMatches;
    }

    /**
     * Removes matches that share atoms with accepted matches of previously processed SMARTS patterns.
     *
     * @param aMatchesList matches of the current SMARTS pattern
     * @param anIsAtomBlockedArray flags of atoms used by accepted matches of previous SMARTS patterns, indexed by atom index
     * @return matches that are atom-disjoint from all accepted matches of previous SMARTS patterns
     */
    private static List<int[]> applyInterSubstructureOverlapFilter(List<int[]> aMatchesList, boolean[] anIsAtomBlockedArray) {
        List<int[]> tmpFilteredMatches = new ArrayList<>(aMatchesList.size());
        for (int[] tmpMatch : aMatchesList) {
            if (!BiosynfoniFragmenter.hasOverlap(tmpMatch, anIsAtomBlockedArray)) {
                tmpFilteredMatches.add(tmpMatch);
            }
        }
        return tmpFilteredMatches;
    }

    /**
     * Checks whether any atom of the given match is flagged as blocked.
     *
     * @param aMatch atom indices of the match
     * @param anIsAtomBlockedArray flags of blocked atoms, indexed by atom index
     * @return true if at least one atom of the match is blocked
     */
    private static boolean hasOverlap(int[] aMatch, boolean[] anIsAtomBlockedArray) {
        for (int tmpAtomIndex : aMatch) {
            if (anIsAtomBlockedArray[tmpAtomIndex]) {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes the canonical rank of every atom, i.e. a numbering of the atoms that does not depend on the input atom
     * order. The atoms are numbered in the order of a depth-first traversal along the bonds that starts at the atom
     * with the smallest canonical label of the CDK canonical labelling ({@link Canon#label(IAtomContainer, int[][])})
     * and visits the neighbours of an atom in the order of their canonical labels; disconnected components are
     * traversed one after the other, in the order of their smallest canonical label. This is the same numbering as in
     * the Biosynfoni fingerprinter that is being developed for the CDK, so that both select the same matches with the
     * intra-substructure overlap filter. It follows the same principle as the atom order of the unique SMILES code and
     * is identical to it for molecules without explicit hydrogen atoms. Because neighbouring atoms get close ranks, the
     * intra-substructure overlap filter selects the matches along chains and rings, similar to the reference
     * implementation, which processes the matches in the order of the (SMILES-based) input atom indices. The matches
     * are thereby packed about as densely, but the selected matches and their number can still differ from the
     * reference implementation (see the class documentation). The canonical labelling does not regard bond orders, so
     * atoms that differ only in the orders of their bonds, e.g. the atoms next to an unset bond in a ring of otherwise
     * identical atoms, can be ranked depending on their input order, and the selected matches can then depend on it as
     * well. If an implicit hydrogen count is unset, for which the default invariants of the canonical labelling throw
     * an exception, the labelling starts from the atom invariants of the unique SMILES code instead, which regard unset
     * implicit hydrogen counts as zero. The method is package-private for testing.
     *
     * @param aMolecule the molecule
     * @return canonical rank of every atom, indexed by atom index
     */
    static int[] getCanonicalRanks(IAtomContainer aMolecule) {
        int[][] tmpGraph = GraphUtil.toAdjList(aMolecule);
        long[] tmpLabels;
        try {
            tmpLabels = Canon.label(aMolecule, tmpGraph);
        } catch (NullPointerException anException) {
            //unset implicit hydrogen counts
            BiosynfoniFragmenter.LOGGER.log(Level.FINE, "Default invariants of the canonical labelling not available,"
                    + " invariants of the unique SMILES code are used instead.", anException);
            tmpLabels = Canon.label(aMolecule, tmpGraph, SmilesGenerator.createComparator(aMolecule, SmiFlavor.Unique));
        }
        long[] tmpFinalLabels = tmpLabels;
        Integer[] tmpAtomsByLabel = new Integer[tmpFinalLabels.length];
        for (int i = 0; i < tmpFinalLabels.length; i++) {
            tmpAtomsByLabel[i] = i;
        }
        Arrays.sort(tmpAtomsByLabel, Comparator.comparingLong((Integer anIndex) -> tmpFinalLabels[anIndex]));
        //neighbours are pushed in descending label order, so that the one with the smallest label is visited first
        Comparator<Integer> tmpDescendingLabel = Comparator.comparingLong((Integer anIndex) -> tmpFinalLabels[anIndex])
                .reversed();
        int[] tmpRanks = new int[tmpFinalLabels.length];
        Arrays.fill(tmpRanks, -1);
        int tmpNextRank = 0;
        Deque<Integer> tmpStack = new ArrayDeque<>();
        for (int tmpStart : tmpAtomsByLabel) {
            if (tmpRanks[tmpStart] >= 0) {
                continue;
            }
            tmpStack.push(tmpStart);
            while (!tmpStack.isEmpty()) {
                int tmpAtom = tmpStack.pop();
                if (tmpRanks[tmpAtom] >= 0) {
                    continue;
                }
                tmpRanks[tmpAtom] = tmpNextRank++;
                Integer[] tmpNeighbours = new Integer[tmpGraph[tmpAtom].length];
                for (int k = 0; k < tmpNeighbours.length; k++) {
                    tmpNeighbours[k] = tmpGraph[tmpAtom][k];
                }
                Arrays.sort(tmpNeighbours, tmpDescendingLabel);
                for (int tmpNeighbour : tmpNeighbours) {
                    if (tmpRanks[tmpNeighbour] < 0) {
                        tmpStack.push(tmpNeighbour);
                    }
                }
            }
        }
        return tmpRanks;
    }
    //</editor-fold>
}
