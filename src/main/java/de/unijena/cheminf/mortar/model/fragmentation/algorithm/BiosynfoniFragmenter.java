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

import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.model.util.SimpleIDisplayEnumConstantProperty;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
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
 * that is stored for this key (see {@link BiosynfoniKey#getRepresentativeSmiles()}). All matches of the same key thus
 * end up as one fragment in MORTAR, and the frequency of this fragment in a molecule equals the count of the
 * respective Biosynfoni count fingerprint position.
 * <br>Alternatively, the actually matched substructures can be returned (see
 * {@link BiosynfoniReturnedFragmentsOption#MATCHED_SUBSTRUCTURES}). This shows which different substructures are
 * covered by a specific SMARTS pattern. Each matched substructure consists of the matched atoms and all bonds between
 * them, free valences are saturated with implicit hydrogen atoms, and aromaticity is re-detected on the fragment.
 * Bonds of aromatic systems that are only partially matched therefore appear in the Kekulé form of the input molecule.
 * <br><b>Known limitation of the matched substructures option:</b> MORTAR identifies fragments by their unique SMILES
 * code across the whole data set. Keys whose SMARTS pattern matches only a single atom (ether, hydroxyl, and epoxide
 * oxygen; terminal nitrogen and both amino acid keys; methyl) can produce identical fragments, e.g. "O" for an ether
 * and a hydroxyl oxygen atom. Such fragments are combined into one fragment that carries the labels of all keys that
 * produced it anywhere in the data set, and its frequency in a molecule is the sum over these keys. The number of
 * matches per key and the representative structures option are not affected.
 * <br>The label and the SMARTS pattern of the key are attached to every returned fragment as properties (see
 * {@link #BIOSYNFONI_LABEL_PROPERTY_KEY}, {@link #BIOSYNFONI_SMARTS_PROPERTY_KEY},
 * {@link IMoleculeFragmenter#FRAGMENT_CATEGORY_PROPERTY_KEY}, and
 * {@link IMoleculeFragmenter#FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY}); the latter is used to group the fragments by
 * their key in the GUI.
 * <br>
 * <br>As in the reference implementation, two optional overlap filters are available as settings:
 * <ul>
 *   <li><b>Intra-substructure overlap filter</b>: prevents multiple matches of the same SMARTS pattern from sharing
 *   atoms. Matches are accepted greedily, as in the reference implementation. To make the result independent of the
 *   input atom order, in contrast to the reference implementation, matches are processed in the order of the
 *   canonical atom ranks of their atoms, i.e. their positions in the unique SMILES code of the molecule.</li>
 *   <li><b>Inter-substructure overlap filter</b>: prevents matches of a SMARTS pattern from sharing atoms with
 *   accepted matches of SMARTS patterns that come earlier in the predefined key order.</li>
 * </ul>
 * Prior to matching, ring perception and Daylight aromaticity detection are applied to a copy of the input molecule
 * (see {@link SmartsPattern#prepare(IAtomContainer)}).
 * <br>
 * <br><b>Differences to the reference implementation</b> (Python, RDKit, version "full_1103"):
 * <ul>
 *   <li><b>Revised SMARTS patterns</b>: 36 of the 39 SMARTS patterns are identical to the reference implementation.
 *   In the patterns of the phenyl-C3 and phenyl-C2 units, the side chain carbon atoms must not be part of any
 *   six-membered ring, while the reference implementation only excludes atoms whose smallest ring is six-membered (see
 *   {@link BiosynfoniKey#D_PHENYL_C3}). This follows the intention of the original patterns not to detect phenyl units
 *   in fused ring systems, e.g. it is no longer relevant whether a six-membered ring is additionally fused to a
 *   five-membered ring, as in aryltetralin lignans. The hydroxyl pattern does not match the oxygen atoms of nitro
 *   groups, nitrate esters, and N-oxides (see {@link BiosynfoniKey#O_HYDROXYL}), which the reference implementation
 *   matches as hydroxyl groups in their charge-separated form. Like this, both notations of these groups give the
 *   same result. The revised patterns are the same as in the BiosynfoniFingerprinter of the CDK.</li>
 *   <li><b>Intra-substructure overlap filter</b>: the reference implementation selects the matches greedily in the
 *   order in which RDKit finds them, so that its counts depend on the input atom order of the molecule. Here, the
 *   matches are selected in a canonical order (see above), so the counts do not depend on the input atom order but can
 *   differ from the reference implementation where both select different but equally valid sets of non-overlapping
 *   matches.</li>
 *   <li><b>Cheminformatics toolkit</b>: the substructure matching is done with the CDK instead of RDKit. RDKit removes
 *   explicit hydrogen atoms, which are counted by the SMARTS primitives D and h. MORTAR suppresses explicit hydrogen
 *   atoms when importing molecules, so this only matters if the fragmenter is used directly with molecules containing
 *   explicit hydrogen atoms.</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class BiosynfoniFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum BiosynfoniKey">
    /**
     * Biosynfoni substructure keys in the order of the fingerprint positions. Each constant holds the label of the
     * key in the reference implementation, its SMARTS pattern, and a representative SMILES code. The representative
     * structure is matched by the SMARTS pattern itself and is returned as fragment for every accepted match of the
     * key. All representative structures are pairwise distinct so that every key yields exactly one fragment.
     */
    public enum BiosynfoniKey {
        //cofactors
        /**
         * Coenzyme A (CoA), matches the CoA backbone including the pantetheine chain, diphosphate linker, and adenosine
         * moiety.
         */
        CO_COA("co_coa",
                "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O",
                "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"),
        /**
         * Reduced nicotinamide adenine dinucleotide (NADH).
         */
        CO_NADH("co_nadh",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O",
                "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"),
        /**
         * Reduced nicotinamide adenine dinucleotide phosphate (NADPH).
         */
        CO_NADPH("co_nadph",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O",
                "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"),
        //amino acids
        /**
         * Standard proteinogenic amino acids (free or peptide-bound); the pattern matches only the amino nitrogen atom.
         * Represented by glycine.
         */
        ALL_STD_AMINOS("allstnd_aminos",
                "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]",
                "NCC(=O)O"),
        /**
         * Amino acids with a standard side chain but a modified amino or carboxy group that are not matched by
         * {@link #ALL_STD_AMINOS}; the pattern matches only the amino nitrogen atom. Represented by alanine methyl
         * ester.
         */
        NON_STD_AMINOS("nonstnd_aminos",
                "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]",
                "COC(=O)C(C)N"),
        //sugar derivatives
        /**
         * Open-chain hexose, linear six-carbon chain with an oxygen atom on every carbon atom.
         */
        S_OPENPYR_C6O6("s_openpyr_C6O6",
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "OCC(O)C(O)C(O)C(O)C=O"),
        /**
         * Open-chain pentose, linear five-carbon chain with an oxygen atom on every carbon atom.
         */
        S_OPENFUR_C5O5("s_openfur_C5O5",
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "OCC(O)C(O)C(O)C=O"),
        /**
         * Pyranose, six-membered sugar ring containing one oxygen atom. Represented by glucopyranose.
         */
        S_PYRANOSE_C5O4("s_pyranose_C5O4",
                "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1",
                "OCC1OC(O)C(O)C(O)C1O"),
        /**
         * Furanose, five-membered sugar ring containing one oxygen atom. Represented by ribofuranose.
         */
        S_FURANOSE_C4O3("s_furanose_C4O3",
                "C~1~[#8]~C~C(~[#8])~C(~[#8])~1",
                "OCC1OC(O)C(O)C1O"),
        //nitrogen-containing building blocks
        /**
         * Indole scaffold with a two-carbon side chain terminating in a nitrogen atom (tryptamine-like).
         */
        D_INDOLE("d_indoleC2N_12",
                "c1cccc2c1c(~[#6]~[#6]~[#7])cn2",
                "NCCc1c[nH]c2ccccc12"),
        /**
         * Phenethylamine, benzene ring connected to a two-carbon chain ending in a nitrogen atom.
         */
        D_PHENYL_C2N("d_phenylC2N_9",
                "c1ccccc1[#6][#6][#7]",
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
         * Phenyl ring attached to a three-carbon side chain. The side chain carbon atoms after the first one must not be
         * part of any six-membered ring. This is a revision of the original pattern, which excludes them via
         * "[#6;!$([r6])]" only if their smallest ring is six-membered, so that, e.g., carbon atoms in a six-membered ring
         * fused to a five-membered ring are not excluded there. Represented by trans-cinnamic acid, the precursor
         * of the phenylpropanoids (C6-C3 units).
         */
        D_PHENYL_C3("d_phenylC3_9_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]~[#6!$(*1*****1)]",
                "OC(=O)/C=C/c1ccccc1"),
        /**
         * Phenyl ring attached to a two-carbon side chain. The second side chain carbon atom must not be part of any
         * six-membered ring, see the revision of the pattern described at {@link #D_PHENYL_C3}. Represented by
         * phenylacetic acid (C6-C2 unit).
         */
        D_PHENYL_C2("d_phenylC2_8_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]",
                "OC(=O)Cc1ccccc1"),
        /**
         * Phenyl ring attached to a single carbon atom. Represented by benzoic acid (C6-C1 unit).
         */
        D_PHENYL_C1("d_phenylC1_7_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]",
                "OC(=O)c1ccccc1"),
        //carbon chains
        /**
         * Isoprene unit, branched five-carbon motif. Represented by isoprene (2-methylbuta-1,3-diene).
         */
        D_ISOPRENE("d_isoprene_5",
                "[#6]~[#6](~[#6])~[#6]~[#6]",
                "C=CC(C)=C"),
        /**
         * Acetyl unit, two-carbon chain connected to an oxygen atom. Represented by acetic acid.
         */
        D2_ACETYL("d2_acetyl_C2O1",
                "[#6]~[#6]~[#8]",
                "CC(=O)O"),
        /**
         * Methylmalonyl unit, three-carbon motif terminating in a methyl carbon. Represented by methylmalonic acid.
         */
        D2_METHYLMALONYL("d2_methylmalonyl_C3",
                "[#6]~[#6][C;D1;h3]",
                "CC(C(=O)O)C(=O)O"),
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
         * Phosphorus atom bonded to an oxygen atom. Represented by phosphoric acid.
         */
        PHOSPHATE("phosphate_2",
                "P~O",
                "OP(=O)(O)O"),
        /**
         * Sulfur atom bonded to an oxygen atom. Represented by methanesulfonic acid.
         */
        SULFONATE("sulfonate_2",
                "S~O",
                "CS(=O)(=O)O"),
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
         * neighbour.
         */
        N_NITRATE("n_nitrate_1",
                "[N;D1]",
                "*N"),
        /**
         * Oxygen atom of an epoxide (three-membered ring).
         */
        O_EPOXY("o_epoxy_1",
                "[O;x2;r3]",
                "C1CO1"),
        /**
         * Acyclic ether oxygen atom, excluding esters, phosphates, and sulfates.
         */
        O_ETHER("o_ether_1",
                "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]",
                "COC"),
        /**
         * Hydroxyl oxygen atom (-OH) or deprotonated hydroxyl oxygen atom (-O<sup>-</sup>) attached to a carbon or
         * nitrogen atom, excluding carboxylic acids, phosphates, and sulfates. This is a revision of the original
         * pattern that additionally excludes the oxygen atoms of nitro groups, nitrate esters, and N-oxides: oxygen
         * atoms bonded to a nitrogen atom that carries a further oxygen atom ("!$(*[#7]~[#8])", the rule of the original
         * pattern for carboxylic acids, "!$(*C~O)", applied to nitrogen) and negatively charged oxygen atoms bonded to a
         * positively charged nitrogen atom ("!$([#8-][#7+])"). These groups can be written in a pentavalent form
         * ("N(=O)=O", "N=O"), whose doubly bonded oxygen atom is not matched anyway, or in a charge-separated form
         * ("[N+](=O)[O-]", "[N+][O-]"), whose singly bonded "[O-]" would otherwise be matched as deprotonated hydroxyl
         * group, so that the result depended on the notation of the same molecule. The reference implementation matches
         * this "[O-]" because RDKit converts nitro groups and N-oxides into the charge-separated form. Hydroxylamines,
         * oximes, and hydroxamic acids, whose uncharged nitrogen atom carries only one oxygen atom, are still matched.
         * Revision taken from the BiosynfoniFingerprinter of the CDK.
         */
        O_HYDROXYL("o_hydroxyl_1",
                "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(*[#7]~[#8]);!$([#8-][#7+]);!$(P);!$(S)]",
                "CO"),
        //carbon rings
        /**
         * Three-membered carbon ring, regardless of aromaticity.
         */
        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1", "C1CC1"),
        /**
         * Four-membered carbon ring, regardless of aromaticity.
         */
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1", "C1CCC1"),
        /**
         * Five-membered carbon ring, regardless of aromaticity.
         */
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1", "C1CCCC1"),
        /**
         * Six-membered carbon ring, regardless of aromaticity.
         */
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCC1"),
        /**
         * Seven-membered carbon ring, regardless of aromaticity.
         */
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCC1"),
        /**
         * Eight-membered carbon ring, regardless of aromaticity.
         */
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCC1"),
        /**
         * Nine-membered carbon ring, regardless of aromaticity.
         */
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCCC1"),
        /**
         * Ten-membered carbon ring, regardless of aromaticity.
         */
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", "C1CCCCCCCCC1");
        /**
         * Label of the key as used in the reference implementation.
         */
        private final String label;
        /**
         * SMARTS pattern of the key.
         */
        private final String smarts;
        /**
         * SMILES code of the representative structure that is returned as fragment for every match of the key.
         */
        private final String representativeSmiles;
        /**
         * Constructor.
         *
         * @param aLabel label of the key as used in the reference implementation
         * @param aSmarts SMARTS pattern of the key
         * @param aRepresentativeSmiles SMILES code of the representative structure of the key
         */
        BiosynfoniKey(String aLabel, String aSmarts, String aRepresentativeSmiles) {
            this.label = aLabel;
            this.smarts = aSmarts;
            this.representativeSmiles = aRepresentativeSmiles;
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
         * Returns the SMARTS pattern of the key.
         *
         * @return SMARTS pattern
         */
        public String getSmarts() {
            return this.smarts;
        }
        //
        /**
         * Returns the SMILES code of the representative structure that is returned as fragment for every match of
         * the key.
         *
         * @return representative SMILES code
         */
        public String getRepresentativeSmiles() {
            return this.representativeSmiles;
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
         * substructures matched by a key become visible. Note that fragments of keys matching only a single atom may
         * coincide and are then combined, see the class documentation.
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
     * Property key under which the SMARTS pattern of the matched key is stored in the returned fragments.
     */
    public static final String BIOSYNFONI_SMARTS_PROPERTY_KEY = "BiosynfoniFragmenter.SMARTS";

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
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Compiled SMARTS patterns in the order of {@link BiosynfoniKey#values()}. Every fragmenter instance holds its
     * own patterns so that instances used in parallel do not share them. The automatic preparation of the target
     * molecule is disabled because it is prepared once per molecule before all patterns are matched.
     */
    private final SmartsPattern[] smartsPatterns;

    /**
     * SMILES parser for creating the representative fragment structures.
     */
    private final SmilesParser smilesParser;

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
        BiosynfoniKey[] tmpKeys = BiosynfoniKey.values();
        this.smartsPatterns = new SmartsPattern[tmpKeys.length];
        for (int i = 0; i < tmpKeys.length; i++) {
            //preparation (ring perception and aromaticity) is done once per molecule in fragmentMolecule()
            this.smartsPatterns[i] = SmartsPattern.create(tmpKeys[i].getSmarts()).setPrepare(false);
        }
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.fragmentAromaticityModel = new Aromaticity(ElectronDonation.daylight(),
                Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
        int tmpNumberOfSettings = 3;
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        int tmpInitialCapacityForSettingMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);

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

        this.returnedFragmentsSetting = new SimpleIDisplayEnumConstantProperty(this, "Returned fragments setting",
                BiosynfoniFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, BiosynfoniReturnedFragmentsOption.class) {
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
        this.settings.add(this.returnedFragmentsSetting);
        this.settingNameToTooltipTextMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("BiosynfoniFragmenter.returnedFragmentsSetting.tooltip"));
        this.settingNameToDisplayNameMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("BiosynfoniFragmenter.returnedFragmentsSetting.displayName"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
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
        tmpCopy.setIntraSubstructureOverlapFilterSetting(this.intraSubstructureOverlapFilterSetting.get());
        tmpCopy.setInterSubstructureOverlapFilterSetting(this.interSubstructureOverlapFilterSetting.get());
        tmpCopy.setReturnedFragmentsSetting((BiosynfoniReturnedFragmentsOption) this.returnedFragmentsSetting.get());
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {
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
        //SMARTS preparation alters the aromaticity and ring flags of the atom container, so work on a clone
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        List<IAtomContainer> tmpFragments = new ArrayList<>();
        boolean tmpAreMatchedSubstructuresReturned = this.returnedFragmentsSetting.get().equals(
                BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        try {
            if (tmpAreMatchedSubstructuresReturned) {
                //the matched substructures are built from the explicit bond orders of the Kekulé form
                for (IBond tmpBond : tmpMoleculeClone.bonds()) {
                    if (Objects.isNull(tmpBond.getOrder()) || tmpBond.getOrder() == IBond.Order.UNSET) {
                        Kekulization.kekulize(tmpMoleculeClone);
                        break;
                    }
                }
            }
            //marks ring atoms and bonds and applies the Daylight aromaticity model once for all patterns, whose
            // automatic preparation is disabled
            SmartsPattern.prepare(tmpMoleculeClone);
            int[] tmpCanonicalRanks = this.intraSubstructureOverlapFilterSetting.get()
                    ? BiosynfoniFragmenter.getCanonicalRanks(tmpMoleculeClone) : null;
            //atoms of accepted matches of previous keys, only used by the inter-substructure overlap filter
            boolean[] tmpIsAtomBlockedByPreviousKeys = new boolean[tmpMoleculeClone.getAtomCount()];
            BiosynfoniKey[] tmpKeys = BiosynfoniKey.values();
            for (int i = 0; i < tmpKeys.length; i++) {
                List<int[]> tmpMatches = new ArrayList<>(Arrays.asList(
                        this.smartsPatterns[i].matchAll(tmpMoleculeClone).uniqueAtoms().toArray()));
                if (this.intraSubstructureOverlapFilterSetting.get()) {
                    tmpMatches = BiosynfoniFragmenter.applyIntraSubstructureOverlapFilter(tmpMatches, tmpCanonicalRanks);
                }
                if (this.interSubstructureOverlapFilterSetting.get()) {
                    tmpMatches = BiosynfoniFragmenter.applyInterSubstructureOverlapFilter(tmpMatches, tmpIsAtomBlockedByPreviousKeys);
                    for (int[] tmpMatch : tmpMatches) {
                        for (int tmpAtomIndex : tmpMatch) {
                            tmpIsAtomBlockedByPreviousKeys[tmpAtomIndex] = true;
                        }
                    }
                }
                if (tmpMatches.isEmpty()) {
                    continue;
                }
                if (tmpAreMatchedSubstructuresReturned) {
                    //TODO: keys matching only a single atom (O_ETHER, O_HYDROXYL, O_EPOXY, N_NITRATE, ALL_STD_AMINOS,
                    // NON_STD_AMINOS, D_METHYL) can produce identical fragments ("O", "N") whose labels and frequencies
                    // are then combined across keys, see class documentation; to be discussed: return the representative
                    // structure for these keys also in this mode, since a single matched atom does not show the
                    // promiscuity of the pattern anyway
                    for (int[] tmpMatch : tmpMatches) {
                        tmpFragments.add(this.createMatchedSubstructureFragment(tmpMoleculeClone, tmpMatch, tmpKeys[i]));
                    }
                } else {
                    //all matches of one key are represented by the same structure and thereby collapse into one fragment
                    IAtomContainer tmpRepresentative = this.createRepresentativeFragment(tmpKeys[i]);
                    for (int j = 0; j < tmpMatches.size(); j++) {
                        tmpFragments.add(tmpRepresentative);
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
     * Creates the representative fragment structure of the given key from its stored SMILES code and annotates it
     * with the key label and SMARTS pattern.
     *
     * @param aKey the Biosynfoni key
     * @return new atom container of the representative structure
     * @throws CDKException if the stored SMILES code cannot be parsed
     */
    private IAtomContainer createRepresentativeFragment(BiosynfoniKey aKey) throws CDKException {
        IAtomContainer tmpRepresentative = this.smilesParser.parseSmiles(aKey.getRepresentativeSmiles());
        BiosynfoniFragmenter.annotateFragment(tmpRepresentative, aKey);
        return tmpRepresentative;
    }

    /**
     * Creates a fragment from the given match that consists of the matched atoms and all bonds between them. Free
     * valences resulting from bonds to non-matched atoms are saturated with implicit hydrogen atoms, and aromaticity
     * is re-detected on the fragment. The fragment is annotated with the key label and SMARTS pattern.
     *
     * @param aMolecule the molecule the match was found in; must have explicit (Kekulé) bond orders
     * @param aMatch atom indices of the match
     * @param aKey the Biosynfoni key that produced the match
     * @return new atom container of the matched substructure
     * @throws CloneNotSupportedException if an atom cannot be cloned
     * @throws CDKException if aromaticity detection on the fragment fails
     */
    private IAtomContainer createMatchedSubstructureFragment(IAtomContainer aMolecule, int[] aMatch, BiosynfoniKey aKey)
            throws CloneNotSupportedException, CDKException {
        IChemObjectBuilder tmpBuilder = aMolecule.getBuilder();
        IAtomContainer tmpFragment = tmpBuilder.newAtomContainer();
        IAtom[] tmpOriginalToFragmentAtom = new IAtom[aMolecule.getAtomCount()];
        for (int tmpAtomIndex : aMatch) {
            IAtom tmpOriginalAtom = aMolecule.getAtom(tmpAtomIndex);
            IAtom tmpFragmentAtom = tmpOriginalAtom.clone();
            int tmpImplicitHydrogenCount = Objects.isNull(tmpOriginalAtom.getImplicitHydrogenCount())
                    ? 0 : tmpOriginalAtom.getImplicitHydrogenCount();
            tmpFragmentAtom.setImplicitHydrogenCount(tmpImplicitHydrogenCount);
            tmpFragmentAtom.setIsAromatic(false);
            tmpFragmentAtom.setIsInRing(false);
            tmpFragment.addAtom(tmpFragmentAtom);
            //the container may wrap the added atom, so the stored reference has to be taken from it
            tmpOriginalToFragmentAtom[tmpAtomIndex] = tmpFragment.getAtom(tmpFragment.getAtomCount() - 1);
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            IAtom tmpBegin = tmpOriginalToFragmentAtom[tmpBond.getBegin().getIndex()];
            IAtom tmpEnd = tmpOriginalToFragmentAtom[tmpBond.getEnd().getIndex()];
            if (!Objects.isNull(tmpBegin) && !Objects.isNull(tmpEnd)) {
                tmpFragment.addBond(tmpBuilder.newInstance(IBond.class, tmpBegin, tmpEnd, tmpBond.getOrder()));
            } else if (!Objects.isNull(tmpBegin) || !Objects.isNull(tmpEnd)) {
                //bond to a non-matched atom is cut, saturate the free valence with hydrogen atoms
                IAtom tmpMatchedAtom = Objects.isNull(tmpBegin) ? tmpEnd : tmpBegin;
                tmpMatchedAtom.setImplicitHydrogenCount(tmpMatchedAtom.getImplicitHydrogenCount() + tmpBond.getOrder().numeric());
            }
        }
        this.fragmentAromaticityModel.apply(tmpFragment);
        BiosynfoniFragmenter.annotateFragment(tmpFragment, aKey);
        return tmpFragment;
    }

    /**
     * Sets the title and the label, SMARTS, category, and display label properties of the given fragment according
     * to the given key. The display label is the readable display name of the key.
     *
     * @param aFragment the fragment to annotate
     * @param aKey the Biosynfoni key that produced the fragment
     */
    private static void annotateFragment(IAtomContainer aFragment, BiosynfoniKey aKey) {
        aFragment.setTitle(aKey.getLabel());
        aFragment.setProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY, aKey.getLabel());
        aFragment.setProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY, aKey.getSmarts());
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
        Arrays.sort(tmpOrder, (a, b) -> Arrays.compare(tmpSortedRanksOfMatches[a], tmpSortedRanksOfMatches[b]));
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
     * Computes the canonical rank of every atom as its position in the unique (canonical) SMILES code of the molecule.
     * This numbering does not depend on the input atom order. Because the SMILES code is written by a depth-first
     * traversal along the bonds, neighbouring atoms get close ranks, so that the intra-substructure overlap filter
     * selects the matches along chains and rings. This packs the matches as densely as the reference implementation,
     * which processes the matches in the order of the (SMILES-based) input atom indices, but without depending on the
     * input atom order. Should the SMILES code not be creatable, e.g. for unset bond orders, the ranks of the CDK
     * canonical labelling are used instead, which are input order independent as well.
     *
     * @param aMolecule the molecule
     * @return canonical rank of every atom, indexed by atom index
     */
    private static int[] getCanonicalRanks(IAtomContainer aMolecule) {
        int[] tmpSmilesOrder = new int[aMolecule.getAtomCount()];
        try {
            new SmilesGenerator(SmiFlavor.Unique).create(aMolecule, tmpSmilesOrder);
            return tmpSmilesOrder;
        } catch (CDKException anException) {
            BiosynfoniFragmenter.LOGGER.log(Level.FINE, "Unique SMILES order not available, canonical labels are used instead.", anException);
        }
        long[] tmpLabels = Canon.label(aMolecule, GraphUtil.toAdjList(aMolecule));
        Integer[] tmpIndices = new Integer[tmpLabels.length];
        for (int i = 0; i < tmpLabels.length; i++) {
            tmpIndices[i] = i;
        }
        Arrays.sort(tmpIndices, Comparator.comparingLong((Integer i) -> tmpLabels[i]).thenComparingInt(i -> i));
        int[] tmpRanks = new int[tmpIndices.length];
        for (int tmpRank = 0; tmpRank < tmpIndices.length; tmpRank++) {
            tmpRanks[tmpIndices[tmpRank]] = tmpRank;
        }
        return tmpRanks;
    }
    //</editor-fold>
}
