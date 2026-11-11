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
import javafx.beans.property.SimpleIntegerProperty;

import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Implementation of the RECAP fragmentation ("Retrosynthetic Combinatorial Analysis Procedure") by Lewell et al., see
 * the references below. RECAP cleaves acyclic bonds that are formed by common chemical reactions, so that the
 * resulting fragments are building blocks for combinatorial chemistry. There are eleven bond types; the attachment
 * points of the fragments are marked with link atoms that carry the number of the bond type:
 * <ol>
 *     <li>amide: carbonyl carbon - nitrogen (no ureas)</li>
 *     <li>ester: carbonyl carbon - ester oxygen</li>
 *     <li>amine: sp3 carbon - amine nitrogen (not acylated, sulfonylated, bonded to an imine or nitrile carbon atom, or
 *     bonded to another heteroatom)</li>
 *     <li>urea: nitrogen - urea carbonyl carbon (both bonds)</li>
 *     <li>ether: sp3 carbon - ether oxygen (not of esters, carbamates or cyanates)</li>
 *     <li>olefin: carbon = carbon (the double bond is cleaved)</li>
 *     <li>quaternary nitrogen: quaternary ammonium nitrogen - carbon</li>
 *     <li>aromatic nitrogen - aliphatic carbon</li>
 *     <li>lactam nitrogen - aliphatic carbon</li>
 *     <li>aromatic carbon - aromatic carbon (biaryls)</li>
 *     <li>sulphonamide: nitrogen - sulfonyl sulfur</li>
 * </ol>
 * Ring bonds are never cleaved. If a bond matches more than one rule, the rule with the lowest number is applied.
 * Protonated amines, e.g. of salts drawn as [NH+] or [NH3+], are not cleaved as amines, since rule 3 requires a neutral
 * nitrogen atom (only quaternary ammonium ions with four carbon substituents are cleaved, by rule 7).
 *
 * <p><b>Rules of the publication and their formulation.</b> The publication (Lewell et al. 1998) shows the eleven
 * bond types as structure drawings (Figure 2) and states that the SMARTS of the program can be defined by the user;
 * the SMARTS are therefore formulated in this class. The following details follow from the text and the fragments
 * shown in the publication:
 * <ul>
 *     <li>All susceptible bonds are cleaved in a single pass, so no intermediate structures occur; only acyclic bonds
 *     are cleaved.</li>
 *     <li>All bonds of the central atom of a group can be cleaved (both bonds of ureas and dialkyl ethers, all bonds
 *     of tertiary amines and quaternary ammonium ions), so fragments that link two or more other fragments, e.g. a
 *     single nitrogen or carbonyl group, are retained.</li>
 *     <li>Ethers and amines are only cleaved at the bond to an sp3 carbon atom: oxygen atoms of aryl ethers remain at
 *     the aromatic ring (e.g. in the cisapride example, Figure 3), and aryl-nitrogen bonds are not cleaved (e.g. the
 *     fragment 1-(2-methoxyphenyl)piperazine in Figure 9), but bonds of aliphatic ring carbon atoms to ether oxygen
 *     atoms are (e.g. sugars in Figures 8 and 9).</li>
 *     <li>Small terminal groups are not cleaved off: according to the text, a terminal fragment that contains only
 *     small functional groups (hydrogen, methyl, ethyl, propyl, butyl) is left uncleaved, so that e.g. the methoxy
 *     groups of cisapride remain at the fragments. The fragments of the publication show that this also applies to
 *     ethyl and methyl esters and acetates, which are not cleaved (Figure 9), and to N-methyl amides (the fragment
 *     indole-CH2-CH([1]NH)-C(=O)-NH-CH3 in Figure 12), but not to alkylated amino groups at amine bonds:
 *     dimethylamino is the most frequent fragment, and diethylamino, isopropylamino and tert-butylamino are frequent
 *     fragments (Figures 9 and 10). Therefore, a terminal part (one attachment point) is not cleaved off if it is an
 *     acyclic hydrocarbon group with at most four carbon atoms (setting, default
 *     {@value #SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT}), if it is an oxygen atom, a carbonyl or a sulfonyl group, or
 *     the nitrogen atom of an amide, urea or sulphonamide bond, with hydrogen or such a hydrocarbon group (e.g.
 *     methoxy, ethoxy, formyl, acetyl, mesyl, the N-methyl and N,N-diethyl groups of amides), or if it is an amino
 *     group (NH2) at any other bond type. Explicit hydrogen atoms that remain after the suppression of hydrogen atoms
 *     (isotopes) are treated like implicit ones, so e.g. a deuterated methoxy group is kept like a methoxy group.
 *     Linking fragments of any size, like the propyl linker of cisapride, are allowed.</li>
 *     <li>The attachment points are labelled with the number of the bond type ("isotopic" labels).</li>
 *     <li>If a bond matches more than one rule, the lowest bond type number is assigned; e.g. bonds of the aromatic
 *     nitrogen atoms of pyrimidinones (uracil) are of type 8, as in Figure 9.</li>
 * </ul>
 *
 * <p><b>Rules in SMIRKS notation.</b> Every rule can be written as a SMIRKS with two mapped atoms, the atoms of the
 * cleaved bond: the left side is the acyclic bond with the environments of its two atoms, the right side replaces it
 * by two bonds to link atoms that carry the number of the bond type. The eleven rules of the list above are (see
 * {@link #getBondTypeSmirks(boolean)}, which also returns the rules with the strict functional group exclusions
 * described below):
 * <pre>{@code
 *  1  [C;$(C=O);!$(C(~[#7])~[#7]):1]-!@[#7;+0:2]>>[*:1]-[1*].[*:2]-[1*]
 *  2  [C;$(C=O):1]-!@[O;D2;$(O([#6])[#6]):2]>>[*:1]-[2*].[*:2]-[2*]
 *  3  [C;X4:1]-!@[N;+0;X3;!$(N[#6]=,#[#7,#8,#16]);!$(N[#7,#8,#15,#16]):2]>>[*:1]-[3*].[*:2]-[3*]
 *  4  [#7;+0:1]-!@[C;$(C(=O)(~[#7])~[#7]):2]>>[*:1]-[4*].[*:2]-[4*]
 *  5  [C;X4:1]-!@[O;D2;+0;$(O([#6])[#6]);!$(O[#6]=,#[#7,#8,#16]):2]>>[*:1]-[5*].[*:2]-[5*]
 *  6  [C;!$(C(=*)=*):1]=!@[C;!$(C(=*)=*):2]>>[*:1]=[6*].[*:2]=[6*]
 *  7  [N+;D4;!$([N+]~[!#6]):1]-!@[#6:2]>>[*:1]-[7*].[*:2]-[7*]
 *  8  [n:1]-!@[C:2]>>[*:1]-[8*].[*:2]-[8*]
 *  9  [N;R;$(N@C=O):1]-!@[C:2]>>[*:1]-[9*].[*:2]-[9*]
 * 10  [c:1]-!@[c:2]>>[*:1]-[10*].[*:2]-[10*]
 * 11  [#7;+0:1]-!@[S;X4;$(S(=O)=O):2]>>[*:1]-[11*].[*:2]-[11*]
 * }</pre>
 * The SMIRKS are a notation of the rules; they are not applied one after the other as reactions. All rules are
 * matched on the intact molecule, a bond that matches several rules is assigned the lowest bond type number, bonds
 * whose cleavage would cut off a small terminal group are skipped, and the remaining bonds are cleaved simultaneously
 * (see below). A sequential application as reactions gives other fragments, since every cleavage changes the
 * environments of the neighbouring atoms: e.g. the nitrogen atom of a cleaved amide is not acylated anymore, so that
 * rule 3 would cleave it off its sp3 carbon atom afterwards.
 *
 * <p><b>Algorithm.</b>
 * <ol>
 *     <li>Explicit hydrogen atoms are suppressed, and ring membership and aromaticity are perceived with the Daylight
 *     model, on which the SMARTS semantics of the rules are based. Since the aromaticity perception needs the bond
 *     orders, a Kekulé structure is assigned first if bond orders are unset (e.g. the aromatic bonds of a molecule
 *     whose kekulization failed on import); molecules that cannot be kekulized are filtered.</li>
 *     <li>Every acyclic, non-aromatic single or double bond is tested against the rules in the order of their numbers;
 *     the atom environments are evaluated lazily and only once per atom.</li>
 *     <li>The cleavable bonds are processed in the order of their bond type number (and bond index); a bond is cleaved
 *     unless its cleavage, together with the bonds cleaved before, cuts off a small terminal group. Since further
 *     cleavages only split fragments, no small terminal group is cut off. All accepted bonds are cleaved
 *     simultaneously, so every heavy atom of the molecule is contained in exactly one fragment.</li>
 *     <li>Each cleaved bond is replaced by two bonds of the same order to new link atoms. A link atom is a pseudo atom
 *     carrying the number of the bond type as mass number (resulting in SMILES like [5*]Oc1ccc(F)cc1, if the mass
 *     number is encoded) and as property {@link #BOND_TYPE_PROPERTY_KEY}. Alternatively, open valences can be
 *     saturated with hydrogen atoms.</li>
 *     <li>Stereo elements are transferred to the fragments; ligands that were cleaved off are replaced by the link atom
 *     (or hydrogen atom) at their position, and stereo elements that lost their stereogenicity are removed.</li>
 * </ol>
 * Pseudo atoms of the input, e.g. attachment points of a preceding fragmentation in a pipeline, count as attachment
 * points in the check for small terminal groups.
 *
 * <p><b>Strict functional groups (setting, default off).</b> The rules of the publication also cleave bonds that are
 * part of larger functional groups. With the strict setting, the following bonds are not cleaved in addition: amides
 * and esters of carbamates, carbonates, haloformates and 1,2-dicarbonyl compounds, anhydrides, imides, acyl ureas,
 * sulphonyl ureas and ureas at imide nitrogen atoms, hydrazides, hydroxamic acids, N-acyl sulphonamides, sulphamides
 * and sulphamates, acyl ammonium ions, bonds of lactam nitrogen atoms in cyclic ureas, carbamates, imides, amidines and
 * hydrazides, bonds of aromatic and lactam nitrogen atoms to carbonyl, imidoyl and thiocarbonyl carbon atoms, olefins
 * with a nitrogen, oxygen or sulfur substituent (e.g. enol ethers, enamines), and carbon-carbon double bonds that are
 * conjugated with another carbon-carbon double or triple bond (polyenes, which would otherwise be cut into two-carbon
 * fragments).
 * Double bonds conjugated with a carbonyl group or an aromatic ring are still cleaved, since they are formed by
 * olefinations like the Horner-Wadsworth-Emmons or Wittig reaction. The exclusions follow the test cases of the RECAP
 * re-implementation for the CDK by Jonas Schaub (see the test resource RECAP_Rule_Cases.tsv).
 *
 * <p><b>RECAP hierarchy (fragment output setting).</b> As an extension of the publication, all fragments of the RECAP
 * hierarchy can be returned instead of the single-pass fragments, see {@link #buildHierarchy(IAtomContainer)}: the
 * children of a fragment are the two fragments that result from the cleavage of one more of its RECAP bonds, unless a
 * small terminal group is cut off. Every set of atoms is represented by one node only, so the hierarchy is a directed
 * acyclic graph whose leaves usually are the single-pass fragments; fragments that occur several times in the molecule
 * are kept as separate nodes. Since the number of nodes grows exponentially with the number of RECAP bonds around a
 * common core, and the total size of the fragments grows with the third power of the number of RECAP bonds in a
 * chain, molecules whose hierarchy has more than {@value #MAXIMUM_HIERARCHY_NODE_COUNT} nodes or more than
 * {@value #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in the fragments of all nodes together are rejected before the fragments
 * are created, which bounds the memory consumption (about 1 KB per fragment atom).
 *
 * <p><b>Differences to the RECAP implementation of the RDKit.</b> The RDKit implementation removes the oxygen atom of
 * cleaved ethers, the nitrogen atom of cleaved amines and the carbonyl group of cleaved ureas, cleaves aryl ethers
 * and aryl amines, cleaves the hydroxy group off carboxylic acids, has no rule for quaternary nitrogen atoms, and only
 * forbids fragments that consist of one to three carbon atoms (so e.g. methoxy, acetyl and butyl groups are cleaved
 * off). Its function RecapDecompose applies the rules hierarchically and returns all unique leaf fragments, which can
 * overlap, whereas the publication cleaves all bonds in a single pass; the multiplicity of fragments is lost. Its
 * hierarchy merges fragments by their SMILES code, whereas the hierarchy of this class merges fragments by their atoms,
 * so that fragments occurring several times are kept.
 *
 * <p><b>Note on bond types in MORTAR.</b> MORTAR identifies fragments by their unique SMILES code without isotope
 * information, so link atoms are written as '*' and fragments that differ only in their bond types are counted as
 * the same fragment.
 *
 * <p>References:
 * <ul>
 *     <li>X. Q. Lewell, D. B. Judd, S. P. Watson, M. M. Hann, "RECAP - Retrosynthetic Combinatorial Analysis
 *     Procedure: A Powerful New Technique for Identifying Privileged Molecular Fragments with Useful Applications in
 *     Combinatorial Chemistry", J. Chem. Inf. Comput. Sci. 1998, 38, 511-522,
 *     <a href="https://doi.org/10.1021/ci970429i">doi:10.1021/ci970429i</a></li>
 *     <li>X. Q. Lewell, D. Judd, S. Watson, M. Hann, "RECAP - Retrosynthetic Combinatorial Analysis Procedure",
 *     presentation, Daylight EMUG 1997, <a href="https://www.daylight.com/meetings/emug97/Lewell/">
 *     www.daylight.com/meetings/emug97/Lewell/</a> (slides "Bond Cleavage Rules", "RECAP Illustration",
 *     "An Example of Cleaving Cisapride")</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class RECAPFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum SaturationOption">
    /**
     * Enumeration of the ways to saturate the open valences that result from the cleavage of a bond.
     */
    public enum SaturationOption implements IDisplayEnum {
        /**
         * Every open valence is saturated with a link atom (pseudo atom) that carries the number of the RECAP bond
         * type.
         */
        LINK_ATOMS(Message.get("RECAPFragmenter.SaturationOption.linkAtoms.displayName"),
                Message.get("RECAPFragmenter.SaturationOption.linkAtoms.tooltip")),
        /**
         * Every open valence is saturated with (implicit) hydrogen atoms, i.e. both atoms of a cleaved double bond
         * (olefin) receive two hydrogen atoms.
         */
        HYDROGEN_SATURATION(Message.get("RECAPFragmenter.SaturationOption.hydrogenSaturation.displayName"),
                Message.get("RECAPFragmenter.SaturationOption.hydrogenSaturation.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor setting the display name and tooltip.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        private SaturationOption(String aDisplayName, String aTooltip) {
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
    //<editor-fold desc="Enum FragmentOutputOption">
    /**
     * Enumeration of the fragments that are returned for a molecule.
     */
    public enum FragmentOutputOption implements IDisplayEnum {
        /**
         * The fragments of the single-pass cleavage of all selected bonds, as in the original publication; every heavy
         * atom of the molecule is contained in exactly one fragment.
         */
        TERMINAL_FRAGMENTS(Message.get("RECAPFragmenter.FragmentOutputOption.terminalFragments.displayName"),
                Message.get("RECAPFragmenter.FragmentOutputOption.terminalFragments.tooltip")),
        /**
         * All fragments of the RECAP hierarchy (see {@link RECAPFragmenter#buildHierarchy(IAtomContainer)}), i.e. the
         * terminal fragments and all intermediate fragments that result from the cleavage of only a part of the bonds;
         * these fragments overlap.
         */
        ALL_HIERARCHY_FRAGMENTS(Message.get("RECAPFragmenter.FragmentOutputOption.allHierarchyFragments.displayName"),
                Message.get("RECAPFragmenter.FragmentOutputOption.allHierarchyFragments.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor setting the display name and tooltip.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        private FragmentOutputOption(String aDisplayName, String aTooltip) {
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
    //<editor-fold desc="Public static nested class RECAPHierarchyNode">
    /**
     * A node of the RECAP hierarchy of a molecule, see {@link RECAPFragmenter#buildHierarchy(IAtomContainer)}: a
     * fragment with references to the fragments it results from by the cleavage of one RECAP bond (parents) and to the
     * fragments that result from it by the cleavage of one more RECAP bond (children). Every set of atoms of the
     * molecule is represented by one node only, so a node can have several parents.
     */
    public static final class RECAPHierarchyNode {
        /**
         * The fragment of this node.
         */
        private final IAtomContainer fragment;
        /**
         * The nodes whose fragments yield this fragment by the cleavage of one RECAP bond.
         */
        private final List<RECAPHierarchyNode> parents;
        /**
         * The nodes whose fragments result from this fragment by the cleavage of one more RECAP bond, in pairs.
         */
        private final List<RECAPHierarchyNode> children;
        /**
         * Constructor.
         *
         * @param aFragment the fragment of the node
         */
        private RECAPHierarchyNode(IAtomContainer aFragment) {
            this.fragment = aFragment;
            this.parents = new ArrayList<>(1);
            this.children = new ArrayList<>(2);
        }
        //
        /**
         * Returns the fragment of this node; for the root, it is the prepared molecule.
         *
         * @return the fragment
         */
        public IAtomContainer getFragment() {
            return this.fragment;
        }
        //
        /**
         * Returns the parent nodes, i.e. the nodes whose fragments yield this fragment by the cleavage of one RECAP
         * bond.
         *
         * @return unmodifiable list of the parent nodes, empty for the root
         */
        public List<RECAPHierarchyNode> getParents() {
            return Collections.unmodifiableList(this.parents);
        }
        //
        /**
         * Returns the child nodes, i.e. the nodes whose fragments result from this fragment by the cleavage of one more
         * RECAP bond. The two fragments of each cleaved bond are adjacent in the list; for the root of a molecule with
         * several disconnected parts, the parts without a cleavable bond follow at the end.
         *
         * @return unmodifiable list of the child nodes, empty for leaves
         */
        public List<RECAPHierarchyNode> getChildren() {
            return Collections.unmodifiableList(this.children);
        }
        //
        /**
         * Returns whether this node is the root of the hierarchy, i.e. has no parents.
         *
         * @return true for the root
         */
        public boolean isRoot() {
            return this.parents.isEmpty();
        }
        //
        /**
         * Returns whether this node is a leaf of the hierarchy, i.e. its fragment cannot be cleaved further.
         *
         * @return true for leaves
         */
        public boolean isLeaf() {
            return this.children.isEmpty();
        }
        //
        /**
         * Returns all nodes below this node, each node once, in breadth-first order.
         *
         * @return list of all descendant nodes
         */
        public List<RECAPHierarchyNode> getAllDescendants() {
            Set<RECAPHierarchyNode> tmpVisitedNodes = Collections.newSetFromMap(new IdentityHashMap<>());
            List<RECAPHierarchyNode> tmpDescendants = new ArrayList<>();
            Deque<RECAPHierarchyNode> tmpQueue = new ArrayDeque<>(this.children);
            while (!tmpQueue.isEmpty()) {
                RECAPHierarchyNode tmpNode = tmpQueue.poll();
                if (tmpVisitedNodes.add(tmpNode)) {
                    tmpDescendants.add(tmpNode);
                    tmpQueue.addAll(tmpNode.children);
                }
            }
            return tmpDescendants;
        }
        //
        /**
         * Returns all leaves below this node, each leaf once, in breadth-first order.
         *
         * @return list of the descendant leaves, empty if this node is a leaf itself
         */
        public List<RECAPHierarchyNode> getLeaves() {
            List<RECAPHierarchyNode> tmpLeaves = new ArrayList<>();
            for (RECAPHierarchyNode tmpNode : this.getAllDescendants()) {
                if (tmpNode.isLeaf()) {
                    tmpLeaves.add(tmpNode);
                }
            }
            return tmpLeaves;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm.
     */
    public static final String ALGORITHM_NAME = "RECAP";
    /**
     * Default option for the saturation of open valences.
     */
    public static final SaturationOption SATURATION_OPTION_DEFAULT = SaturationOption.LINK_ATOMS;
    /**
     * Default setting for returning molecules without cleavable bonds as one unchanged fragment.
     */
    public static final boolean RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT = true;
    /**
     * Default maximum number of carbon atoms of the small terminal groups that are not cleaved off: methyl, ethyl,
     * propyl, and butyl, as in the original publication.
     */
    public static final int SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT = 4;
    /**
     * Default setting for the strict functional group rules: off, i.e. the rules of the original publication are applied.
     */
    public static final boolean STRICT_FUNCTIONAL_GROUPS_OPTION_DEFAULT = false;
    /**
     * Default option for the returned fragments: the terminal fragments of the single-pass cleavage, as in the original
     * publication.
     */
    public static final FragmentOutputOption FRAGMENT_OUTPUT_OPTION_DEFAULT = FragmentOutputOption.TERMINAL_FRAGMENTS;
    /**
     * Maximum number of nodes of a RECAP hierarchy (see {@link #buildHierarchy(IAtomContainer)}). The number of nodes
     * grows exponentially with the number of cleavable bonds around a common core (2^k + k nodes, including the root,
     * for k substituents), so molecules whose hierarchy exceeds this number are rejected.
     */
    public static final int MAXIMUM_HIERARCHY_NODE_COUNT = 10000;
    /**
     * Maximum total number of atoms of the fragments of all nodes of a RECAP hierarchy (see
     * {@link #buildHierarchy(IAtomContainer)}), without link atoms and including the root. The fragments need about
     * 1 KB of memory per atom, and their total size grows with the third power of the number of cleavable bonds in a
     * chain (e.g. about 136,000 atoms for a peptide of 40 amino acids) and exponentially with the number of cleavable
     * bonds around a common core, so molecules whose hierarchy exceeds this number are rejected before its fragments
     * are created.
     */
    public static final int MAXIMUM_HIERARCHY_ATOM_COUNT = 200000;
    /**
     * Key of the atom property of link atoms that holds the number of their RECAP bond type (Integer, 1 to 11).
     */
    public static final String BOND_TYPE_PROPERTY_KEY = "RECAPFragmenter.BondType";
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(RECAPFragmenter.class.getName());
    /**
     * All RECAP bond rules in the order of their bond type numbers.
     */
    private static final List<BondRule> RULES = RECAPFragmenter.createRules(false);
    /**
     * All RECAP bond rules with the strict functional group exclusions, in the order of their bond type numbers.
     */
    private static final List<BondRule> STRICT_RULES = RECAPFragmenter.createRules(true);
    /**
     * Order in which cleavable bonds are processed: by bond type number, then by bond index.
     */
    private static final Comparator<RECAPBond> BOND_PROCESSING_ORDER = Comparator
            .comparingInt(RECAPBond::bondType)
            .thenComparingInt(aRECAPBond -> aRECAPBond.bond().getIndex());
    //</editor-fold>
    //
    //<editor-fold desc="Private records">
    /**
     * A RECAP bond rule: an acyclic bond of the given order between an atom matching the left environment and an atom
     * matching the right environment is cleavable; both atoms are labelled with the bond type.
     *
     * @param bondType number of the RECAP bond type, 1 to 11
     * @param name name of the bond type
     * @param leftEnvironment SMARTS of the left atom, one atom in square brackets
     * @param order order of the bond, single or double
     * @param rightEnvironment SMARTS of the right atom, one atom in square brackets
     */
    private record BondRule(int bondType, String name, String leftEnvironment, IBond.Order order, String rightEnvironment) {
        /**
         * Returns this rule in SMIRKS notation, see {@link RECAPFragmenter#getBondTypeSmirks(boolean)}: the atom maps 1
         * and 2 are inserted into the environments of the left and the right atom.
         *
         * @return SMIRKS of the rule
         */
        String toSmirks() {
            String tmpBond = this.order == IBond.Order.DOUBLE ? "=" : "-";
            return String.format("%s:1]%s!@%s:2]>>[*:1]%s[%d*].[*:2]%s[%d*]",
                    this.leftEnvironment.substring(0, this.leftEnvironment.length() - 1), tmpBond,
                    this.rightEnvironment.substring(0, this.rightEnvironment.length() - 1),
                    tmpBond, this.bondType, tmpBond, this.bondType);
        }
    }
    //
    /**
     * A bond rule whose environments are referenced by their index in the compiled environment patterns.
     *
     * @param rule the rule
     * @param leftEnvironmentIndex index of the left environment pattern
     * @param rightEnvironmentIndex index of the right environment pattern
     */
    private record CompiledRule(BondRule rule, int leftEnvironmentIndex, int rightEnvironmentIndex) {
    }
    //
    /**
     * A bond that matches a RECAP rule.
     *
     * @param bond the bond in the prepared molecule
     * @param bondType number of the RECAP bond type
     */
    record RECAPBond(IBond bond, int bondType) {
    }
    //
    /**
     * A bond after its cleavage, with the new bonds of its former begin and end atom to their link atoms.
     *
     * @param begin former begin atom
     * @param end former end atom
     * @param newBeginBond new bond of the begin atom to its link atom, null if saturated with implicit hydrogen atoms
     * @param newEndBond new bond of the end atom to its link atom, null if saturated with implicit hydrogen atoms
     */
    private record CleavedBond(IAtom begin, IAtom end, IBond newBeginBond, IBond newEndBond) {
        /**
         * Returns the new bond to the link atom of the given atom of the cleaved bond.
         *
         * @param anAtom begin or end atom of the cleaved bond
         * @return new bond of the atom to its link atom, null if saturated with implicit hydrogen atoms
         */
        IBond newBondOf(IAtom anAtom) {
            return anAtom.equals(this.begin) ? this.newBeginBond : this.newEndBond;
        }
        //
        /**
         * Returns the other atom of the cleaved bond.
         *
         * @param anAtom begin or end atom of the cleaved bond
         * @return the other atom
         */
        IAtom other(IAtom anAtom) {
            return anAtom.equals(this.begin) ? this.end : this.begin;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Compiled patterns of all distinct atom environments of the rules, matched with their root atom.
     */
    private final DfPattern[] environmentPatterns;
    /**
     * The rules with references to the compiled environment patterns.
     */
    private final CompiledRule[] compiledRules;
    /**
     * The rules with the strict functional group exclusions, with references to the compiled environment patterns.
     */
    private final CompiledRule[] compiledStrictRules;
    //the settings are documented at their property methods (a Javadoc comment here would duplicate it for javadoc)
    //setting for the saturation of open valences
    private final SimpleIDisplayEnumConstantProperty saturationSetting;
    //setting for the strict functional group rules
    private final SimpleBooleanProperty strictFunctionalGroupsSetting;
    //setting for the returned fragments (terminal fragments or all fragments of the hierarchy)
    private final SimpleIDisplayEnumConstantProperty fragmentOutputSetting;
    //setting for returning molecules without cleavable bonds as one unchanged fragment
    private final SimpleBooleanProperty returnUnfragmentedMoleculesSetting;
    //setting for the maximum number of carbon atoms of small terminal groups that are not cleaved off
    private final SimpleIntegerProperty smallGroupMaximumCarbonAtomsSetting;
    /**
     * All settings of this fragmenter, encapsulated in JavaFX properties for binding in GUI.
     */
    private final List<Property<?>> settings;
    /**
     * Map to store pairs of {@literal <setting name, tooltip text>}.
     */
    private final HashMap<String, String> settingNameTooltipTextMap;
    /**
     * Map to store pairs of {@literal <setting name, display name>}.
     */
    private final HashMap<String, String> settingNameDisplayNameMap;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, all settings are initialised with their default values as declared in the respective public
     * constants. The SMARTS patterns of the rules are compiled for every instance, so that instances can be used in
     * different threads.
     */
    public RECAPFragmenter() {
        //both rule sets share the compiled environment patterns
        Map<String, Integer> tmpEnvironmentIndexMap = new LinkedHashMap<>();
        this.compiledRules = RECAPFragmenter.compileRules(RECAPFragmenter.RULES, tmpEnvironmentIndexMap);
        this.compiledStrictRules = RECAPFragmenter.compileRules(RECAPFragmenter.STRICT_RULES, tmpEnvironmentIndexMap);
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        this.environmentPatterns = new DfPattern[tmpEnvironmentIndexMap.size()];
        for (Map.Entry<String, Integer> tmpEntry : tmpEnvironmentIndexMap.entrySet()) {
            QueryAtomContainer tmpQuery = new QueryAtomContainer(tmpBuilder);
            if (!Smarts.parse(tmpQuery, tmpEntry.getKey())) {
                throw new IllegalStateException(String.format("RECAP environment SMARTS %s cannot be parsed: %s",
                        tmpEntry.getKey(), Smarts.getLastErrorMesg()));
            }
            this.environmentPatterns[tmpEntry.getValue()] = DfPattern.findSubstructure(tmpQuery);
        }
        int tmpNumberOfSettings = 5;
        int tmpInitialCapacityForSettingNameMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.saturationSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Fragment saturation setting",
                RECAPFragmenter.SATURATION_OPTION_DEFAULT,
                RECAPFragmenter.SaturationOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    RECAPFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.saturationSetting.getName(),
                Message.get("RECAPFragmenter.saturationSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.saturationSetting.getName(),
                Message.get("RECAPFragmenter.saturationSetting.displayName"));
        this.smallGroupMaximumCarbonAtomsSetting = new SimpleIntegerProperty(this,
                "Maximum size of small terminal groups setting",
                RECAPFragmenter.SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT) {
            @Override
            public void set(int newValue) {
                if (newValue >= 0) {
                    super.set(newValue);
                } else {
                    IllegalArgumentException tmpException = new IllegalArgumentException(
                            "The maximum size of small terminal groups must not be negative.");
                    RECAPFragmenter.LOGGER.log(Level.WARNING, tmpException.toString(), tmpException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            tmpException.toString(),
                            tmpException);
                    //re-throws the exception to properly reset the binding
                    throw tmpException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.smallGroupMaximumCarbonAtomsSetting.getName(),
                Message.get("RECAPFragmenter.smallGroupMaximumCarbonAtomsSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.smallGroupMaximumCarbonAtomsSetting.getName(),
                Message.get("RECAPFragmenter.smallGroupMaximumCarbonAtomsSetting.displayName"));
        this.returnUnfragmentedMoleculesSetting = new SimpleBooleanProperty(this,
                "Return unfragmented molecules setting",
                RECAPFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("RECAPFragmenter.returnUnfragmentedMoleculesSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("RECAPFragmenter.returnUnfragmentedMoleculesSetting.displayName"));
        this.strictFunctionalGroupsSetting = new SimpleBooleanProperty(this,
                "Strict functional groups setting",
                RECAPFragmenter.STRICT_FUNCTIONAL_GROUPS_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.strictFunctionalGroupsSetting.getName(),
                Message.get("RECAPFragmenter.strictFunctionalGroupsSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.strictFunctionalGroupsSetting.getName(),
                Message.get("RECAPFragmenter.strictFunctionalGroupsSetting.displayName"));
        this.fragmentOutputSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Fragment output setting",
                RECAPFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT,
                RECAPFragmenter.FragmentOutputOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    RECAPFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.fragmentOutputSetting.getName(),
                Message.get("RECAPFragmenter.fragmentOutputSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.fragmentOutputSetting.getName(),
                Message.get("RECAPFragmenter.fragmentOutputSetting.displayName"));
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        this.settings.add(this.saturationSetting);
        this.settings.add(this.smallGroupMaximumCarbonAtomsSetting);
        this.settings.add(this.strictFunctionalGroupsSetting);
        this.settings.add(this.fragmentOutputSetting);
        this.settings.add(this.returnUnfragmentedMoleculesSetting);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the currently set option for the saturation of open valences.
     *
     * @return saturation option
     */
    public SaturationOption getSaturationSetting() {
        return (SaturationOption) this.saturationSetting.get();
    }
    //
    /**
     * Returns the property object of the saturation setting that can be used to configure this setting.
     *
     * @return property object of the saturation setting
     */
    public SimpleIDisplayEnumConstantProperty saturationSettingProperty() {
        return this.saturationSetting;
    }
    //
    /**
     * Returns the maximum number of carbon atoms of the small terminal groups that are not cleaved off (see class
     * documentation).
     *
     * @return maximum size of small terminal groups in carbon atoms
     */
    public int getSmallGroupMaximumCarbonAtomsSetting() {
        return this.smallGroupMaximumCarbonAtomsSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the maximum size of small terminal groups.
     *
     * @return property object of the setting
     */
    public SimpleIntegerProperty smallGroupMaximumCarbonAtomsSettingProperty() {
        return this.smallGroupMaximumCarbonAtomsSetting;
    }
    //
    /**
     * Returns whether molecules without cleavable bonds are returned as one unchanged fragment.
     *
     * @return true if molecules without cleavable bonds are returned unchanged
     */
    public boolean getReturnUnfragmentedMoleculesSetting() {
        return this.returnUnfragmentedMoleculesSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for returning molecules without cleavable bonds.
     *
     * @return property object of the setting
     */
    public SimpleBooleanProperty returnUnfragmentedMoleculesSettingProperty() {
        return this.returnUnfragmentedMoleculesSetting;
    }
    //
    /**
     * Returns whether the strict functional group rules are applied (see class documentation).
     *
     * @return true if the strict functional group rules are applied
     */
    public boolean getStrictFunctionalGroupsSetting() {
        return this.strictFunctionalGroupsSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the strict functional group rules.
     *
     * @return property object of the setting
     */
    public SimpleBooleanProperty strictFunctionalGroupsSettingProperty() {
        return this.strictFunctionalGroupsSetting;
    }
    //
    /**
     * Returns the currently set option for the returned fragments.
     *
     * @return fragment output option
     */
    public FragmentOutputOption getFragmentOutputSetting() {
        return (FragmentOutputOption) this.fragmentOutputSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the returned fragments.
     *
     * @return property object of the setting
     */
    public SimpleIDisplayEnumConstantProperty fragmentOutputSettingProperty() {
        return this.fragmentOutputSetting;
    }
    //
    /**
     * Returns the names of all RECAP bond types in the order of their numbers (index 0 is bond type 1).
     *
     * @return unmodifiable list of bond type names
     */
    public static List<String> getBondTypeNames() {
        List<String> tmpNames = new ArrayList<>(RECAPFragmenter.RULES.size());
        for (BondRule tmpRule : RECAPFragmenter.RULES) {
            tmpNames.add(tmpRule.name());
        }
        return Collections.unmodifiableList(tmpNames);
    }
    //
    /**
     * Returns all RECAP bond rules in SMIRKS notation, in the order of their bond type numbers (index 0 is bond type
     * 1), e.g. {@code [C;X4:1]-!@[N;...:2]>>[*:1]-[3*].[*:2]-[3*]} for amines. The atoms with the maps 1 and 2 are the
     * atoms of the cleaved bond with their environments; the right side replaces the bond by two bonds to link atoms
     * that carry the bond type number. The SMIRKS are created from the rule definitions that the fragmenter matches, so
     * they always describe the rules in use. They are a notation only: the fragmenter matches all rules on the intact
     * molecule, assigns the lowest bond type number to a bond that matches several rules, skips bonds that would cut
     * off a small terminal group, and cleaves all remaining bonds simultaneously (see class documentation), which a
     * sequential application of the SMIRKS as reactions does not reproduce.
     *
     * @param isStrict true for the rules with the strict functional group exclusions (see class documentation)
     * @return unmodifiable list of the eleven rules as SMIRKS
     */
    public static List<String> getBondTypeSmirks(boolean isStrict) {
        List<BondRule> tmpRules = isStrict ? RECAPFragmenter.STRICT_RULES : RECAPFragmenter.RULES;
        List<String> tmpSmirksList = new ArrayList<>(tmpRules.size());
        for (BondRule tmpRule : tmpRules) {
            tmpSmirksList.add(tmpRule.toSmirks());
        }
        return Collections.unmodifiableList(tmpSmirksList);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set">
    /**
     * Sets the option for the saturation of open valences.
     *
     * @param anOption saturation option
     * @throws NullPointerException if the given option is null
     */
    public void setSaturationSetting(SaturationOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given saturation option is null.");
        this.saturationSetting.set(anOption);
    }
    //
    /**
     * Sets the maximum number of carbon atoms of the small terminal groups that are not cleaved off (see class
     * documentation). With 0, only hydrogen-terminated groups (e.g. amino, formyl) are protected.
     *
     * @param aSize maximum size of small terminal groups in carbon atoms, must not be negative
     * @throws IllegalArgumentException if the given size is negative
     */
    public void setSmallGroupMaximumCarbonAtomsSetting(int aSize) throws IllegalArgumentException {
        if (aSize < 0) {
            throw new IllegalArgumentException("The maximum size of small terminal groups must not be negative.");
        }
        this.smallGroupMaximumCarbonAtomsSetting.set(aSize);
    }
    //
    /**
     * Sets whether molecules without cleavable bonds are returned as one unchanged fragment. If false, no fragments
     * are returned for such molecules.
     *
     * @param aBoolean true if molecules without cleavable bonds should be returned unchanged
     */
    public void setReturnUnfragmentedMoleculesSetting(boolean aBoolean) {
        this.returnUnfragmentedMoleculesSetting.set(aBoolean);
    }
    //
    /**
     * Sets whether the strict functional group rules are applied (see class documentation). If false (default), the
     * rules of the original publication are applied.
     *
     * @param aBoolean true if the strict functional group rules should be applied
     */
    public void setStrictFunctionalGroupsSetting(boolean aBoolean) {
        this.strictFunctionalGroupsSetting.set(aBoolean);
    }
    //
    /**
     * Sets the option for the returned fragments.
     *
     * @param anOption fragment output option
     * @throws NullPointerException if the given option is null
     */
    public void setFragmentOutputSetting(FragmentOutputOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given fragment output option is null.");
        this.fragmentOutputSetting.set(anOption);
    }
    //</editor-fold>
    //
    //<editor-fold desc="IMoleculeFragmenter methods">
    @Override
    public List<Property<?>> settingsProperties() {
        return this.settings;
    }
    //
    @Override
    public Map<String, String> getSettingNameToTooltipTextMap() {
        return this.settingNameTooltipTextMap;
    }
    //
    @Override
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return this.settingNameDisplayNameMap;
    }
    //
    @Override
    public String getFragmentationAlgorithmName() {
        return RECAPFragmenter.ALGORITHM_NAME;
    }
    //
    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("RECAPFragmenter.displayName");
    }
    //
    @Override
    public IMoleculeFragmenter copy() {
        RECAPFragmenter tmpCopy = new RECAPFragmenter();
        tmpCopy.setSaturationSetting(this.getSaturationSetting());
        tmpCopy.setSmallGroupMaximumCarbonAtomsSetting(this.getSmallGroupMaximumCarbonAtomsSetting());
        tmpCopy.setReturnUnfragmentedMoleculesSetting(this.getReturnUnfragmentedMoleculesSetting());
        tmpCopy.setStrictFunctionalGroupsSetting(this.getStrictFunctionalGroupsSetting());
        tmpCopy.setFragmentOutputSetting(this.getFragmentOutputSetting());
        return tmpCopy;
    }
    //
    @Override
    public void restoreDefaultSettings() {
        this.saturationSetting.set(RECAPFragmenter.SATURATION_OPTION_DEFAULT);
        this.smallGroupMaximumCarbonAtomsSetting.set(RECAPFragmenter.SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT);
        this.returnUnfragmentedMoleculesSetting.set(RECAPFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
        this.strictFunctionalGroupsSetting.set(RECAPFragmenter.STRICT_FUNCTIONAL_GROUPS_OPTION_DEFAULT);
        this.fragmentOutputSetting.set(RECAPFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT);
    }
    //
    /**
     * {@inheritDoc}
     *
     * <p>The fragments are the disconnected parts of a clone of the given molecule, in which all selected RECAP bonds
     * have been cleaved and the open valences saturated according to the saturation setting. With the fragment output
     * option {@link FragmentOutputOption#ALL_HIERARCHY_FRAGMENTS}, all fragments of the RECAP hierarchy except its root
     * (the molecule) are returned instead, see {@link #buildHierarchy(IAtomContainer)}. If no bond is cleaved, the
     * clone (or its disconnected parts, e.g. counter ions, separately) is returned or, depending on the respective
     * setting, no fragment is returned.</p>
     *
     * @throws IllegalArgumentException also if the molecule has bonds with unset order and cannot be kekulized, or,
     *                                  with the fragment output option
     *                                  {@link FragmentOutputOption#ALL_HIERARCHY_FRAGMENTS}, if the RECAP hierarchy of
     *                                  the molecule has more than {@link #MAXIMUM_HIERARCHY_NODE_COUNT} nodes or more
     *                                  than {@link #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in all fragments together
     */
    @Override
    public List<IAtomContainer> fragmentMolecule(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Parameter tests">
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (!this.canBeFragmented(aMolecule)) {
            throw new IllegalArgumentException("Given molecule cannot be fragmented but should be filtered or preprocessed first.");
        }
        //</editor-fold>
        IAtomContainer tmpMolecule = aMolecule.clone();
        try {
            this.prepareMolecule(tmpMolecule);
            List<RECAPBond> tmpRECAPBonds = this.findRECAPBonds(tmpMolecule);
            List<RECAPBond> tmpBondsToCleave = this.selectBondsToCleave(tmpMolecule, tmpRECAPBonds);
            //if the single pass cleaves no bond, no bond can be cleaved in the hierarchy either, since the cleavage of
            //each bond alone cuts off a small terminal group
            if (!tmpBondsToCleave.isEmpty()) {
                if (this.fragmentOutputSetting.get().equals(FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS)) {
                    List<RECAPHierarchyNode> tmpNodes = this.createHierarchy(tmpMolecule, tmpRECAPBonds);
                    List<IAtomContainer> tmpFragments = new ArrayList<>(tmpNodes.size() - 1);
                    //the root, i.e. the molecule itself, is not a fragment
                    for (int i = 1; i < tmpNodes.size(); i++) {
                        tmpFragments.add(tmpNodes.get(i).getFragment());
                    }
                    return tmpFragments;
                }
                return this.cleaveBonds(tmpMolecule, tmpBondsToCleave);
            }
            if (!this.returnUnfragmentedMoleculesSetting.get()) {
                return new ArrayList<>(0);
            }
            if (ConnectivityChecker.isConnected(tmpMolecule)) {
                List<IAtomContainer> tmpResult = new ArrayList<>(1);
                tmpResult.add(tmpMolecule);
                return tmpResult;
            }
            //disconnected parts, e.g. counter ions, are returned separately, as for molecules with cleaved bonds
            IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpMolecule);
            List<IAtomContainer> tmpResult = new ArrayList<>(tmpComponents.getAtomContainerCount());
            for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
                tmpResult.add(tmpComponent);
            }
            return tmpResult;
        } catch (CDKException | RuntimeException anException) {
            throw new IllegalArgumentException(String.format("An error occurred during fragmentation: %s Molecule Name: %s",
                    anException, aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)), anException);
        }
    }
    //
    /**
     * {@inheritDoc}
     *
     * <p>Molecules are filtered if they are null or empty, or if the order of a bond is unset and cannot be determined
     * by kekulization (see {@link #prepareMolecule(IAtomContainer)}), e.g. the aromatic remnants of opened rings in
     * fragments of a preceding fragmentation, which MORTAR parses without kekulization because their kekulization
     * fails. The given molecule is not altered.</p>
     */
    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        return Objects.isNull(aMolecule) || aMolecule.isEmpty() || this.hasUnresolvableBondOrders(aMolecule);
    }
    //
    @Override
    public boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return false;
    }
    //
    @Override
    public boolean canBeFragmented(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered = this.shouldBeFiltered(aMolecule);
        boolean tmpShouldBePreprocessed = this.shouldBePreprocessed(aMolecule);
        return !(tmpShouldBeFiltered || tmpShouldBePreprocessed);
    }
    //
    @Override
    public IAtomContainer applyPreprocessing(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (this.shouldBeFiltered(aMolecule)) {
            throw new IllegalArgumentException("The given molecule cannot be preprocessed but should be filtered.");
        }
        return aMolecule.clone();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Builds the RECAP hierarchy of the given molecule with the current settings (an extension of the original
     * publication, which cleaves all bonds in a single pass, comparable to the function RecapDecompose of the RDKit).
     * The root of the hierarchy is the molecule; the children of a node are the two fragments that result from the
     * cleavage of one more RECAP bond of its fragment, unless this cleavage cuts off a small terminal group (with
     * respect to the bonds already cleaved for this fragment). Every fragment, i.e. every set of atoms of the molecule,
     * is represented by exactly one node, which has a parent for every bond whose cleavage yields it, so the hierarchy
     * is a directed acyclic graph whose size grows polynomially with the number of RECAP bonds for chains, but
     * exponentially for many substituents around a common core. Fragments that occur several times in the molecule
     * are represented by several nodes, so their multiplicity is kept. If the molecule consists of several disconnected
     * parts, e.g. a salt, the parts without a cleavable bond are children of the root as well. The leaves of the
     * hierarchy usually are the fragments of the single-pass cleavage. The link atoms (or hydrogen atoms) and stereo
     * elements of the fragments are created as in {@link #fragmentMolecule(IAtomContainer)}. To bound the memory
     * consumption, the size of the hierarchy is checked before its fragments are created: it must not have more than
     * {@link #MAXIMUM_HIERARCHY_NODE_COUNT} nodes, and the fragments of all nodes together must not have more than
     * {@link #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms (without link atoms).
     *
     * @param aMolecule the molecule to fragment, it is not altered
     * @return the root node of the hierarchy, its fragment is the prepared clone of the molecule (explicit hydrogen
     *         atoms suppressed, Kekulé structure assigned if bond orders were unset, Daylight aromaticity)
     * @throws NullPointerException if the given molecule is null
     * @throws IllegalArgumentException if the molecule cannot be fragmented (see
     *                                  {@link #canBeFragmented(IAtomContainer)}), if an error occurs during its
     *                                  preparation (e.g. bonds with unset order that cannot be kekulized), or if its
     *                                  hierarchy has more than {@link #MAXIMUM_HIERARCHY_NODE_COUNT} nodes or more than
     *                                  {@link #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in all fragments together
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    public RECAPHierarchyNode buildHierarchy(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (!this.canBeFragmented(aMolecule)) {
            throw new IllegalArgumentException("Given molecule cannot be fragmented but should be filtered or preprocessed first.");
        }
        IAtomContainer tmpMolecule = aMolecule.clone();
        try {
            this.prepareMolecule(tmpMolecule);
        } catch (CDKException anException) {
            throw new IllegalArgumentException(String.format("The molecule could not be prepared: %s", anException),
                    anException);
        }
        return this.createHierarchy(tmpMolecule, this.findRECAPBonds(tmpMolecule)).get(0);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Package-private methods">
    /**
     * Prepares the given molecule (in place) for the matching of the RECAP rules: missing implicit hydrogen counts are
     * set, a Kekulé structure is assigned if the order of a bond is unset (e.g. the aromatic bonds of a molecule
     * parsed without kekulization or whose kekulization failed on import), explicit hydrogen atoms are suppressed, and
     * ring membership and aromaticity (Daylight model) are perceived. The Kekulé structure is needed because the
     * aromaticity perception clears the aromaticity flags of the input and determines the aromaticity from the bond
     * orders, so aromatic rings with unset bond orders would otherwise be regarded as non-aromatic.
     *
     * @param aMolecule the molecule to prepare, it is altered
     * @throws CDKException if implicit hydrogen counts are missing and cannot be determined, or if the order of a bond
     *                      is unset and no Kekulé structure can be assigned
     */
    void prepareMolecule(IAtomContainer aMolecule) throws CDKException {
        boolean tmpIsHydrogenCountMissing = false;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom.getImplicitHydrogenCount())) {
                if (tmpAtom instanceof IPseudoAtom) {
                    tmpAtom.setImplicitHydrogenCount(0);
                } else {
                    tmpIsHydrogenCountMissing = true;
                }
            }
        }
        if (tmpIsHydrogenCountMissing) {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
            CDKHydrogenAdder.getInstance(aMolecule.getBuilder()).addImplicitHydrogens(aMolecule);
        }
        boolean tmpIsBondOrderUnset = false;
        for (IBond tmpBond : aMolecule.bonds()) {
            if (Objects.isNull(tmpBond.getOrder()) || tmpBond.getOrder() == IBond.Order.UNSET) {
                //the kekulization expects UNSET instead of null
                tmpBond.setOrder(IBond.Order.UNSET);
                tmpIsBondOrderUnset = true;
            }
        }
        if (tmpIsBondOrderUnset) {
            //assigns single and double bonds to the bonds of the aromatic atoms, throws a CDKException if not possible
            Kekulization.kekulize(aMolecule);
            for (IBond tmpBond : aMolecule.bonds()) {
                //the kekulization leaves bonds unset that are not flagged as aromatic and receive no double bond
                if (tmpBond.getOrder() == IBond.Order.UNSET) {
                    throw new CDKException(String.format("The order of the bond between the atoms %d and %d is unset "
                            + "and cannot be determined by kekulization.", tmpBond.getBegin().getIndex(),
                            tmpBond.getEnd().getIndex()));
                }
            }
        }
        AtomContainerManipulator.suppressHydrogens(aMolecule);
        //marks ring atoms and bonds, and applies the Daylight aromaticity model after clearing existing aromaticity flags
        SmartsPattern.prepare(aMolecule);
    }
    //
    /**
     * Identifies all bonds of the given prepared molecule that match a RECAP rule (with the strict functional group
     * exclusions if the respective setting is true). If a bond matches several rules, the rule with the lowest bond
     * type number is assigned.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @return the cleavable bonds with their bond types, in the order of the bonds in the molecule
     */
    List<RECAPBond> findRECAPBonds(IAtomContainer aPreparedMolecule) {
        CompiledRule[] tmpRules = this.strictFunctionalGroupsSetting.get() ? this.compiledStrictRules : this.compiledRules;
        //lazily filled results of the environment patterns per atom index, 0: unknown, 1: match, 2: no match
        byte[][] tmpEnvironmentMatches = new byte[this.environmentPatterns.length][];
        int tmpAtomCount = aPreparedMolecule.getAtomCount();
        List<RECAPBond> tmpRECAPBonds = new ArrayList<>();
        for (IBond tmpBond : aPreparedMolecule.bonds()) {
            IBond.Order tmpOrder = tmpBond.getOrder();
            if (tmpBond.isInRing() || tmpBond.isAromatic()
                    || (tmpOrder != IBond.Order.SINGLE && tmpOrder != IBond.Order.DOUBLE)) {
                continue;
            }
            IAtom tmpBegin = tmpBond.getBegin();
            IAtom tmpEnd = tmpBond.getEnd();
            for (CompiledRule tmpCompiledRule : tmpRules) {
                if (tmpCompiledRule.rule().order() != tmpOrder) {
                    continue;
                }
                boolean tmpIsMatch =
                        (this.matches(tmpCompiledRule.leftEnvironmentIndex(), tmpBegin, tmpAtomCount, tmpEnvironmentMatches)
                        && this.matches(tmpCompiledRule.rightEnvironmentIndex(), tmpEnd, tmpAtomCount, tmpEnvironmentMatches))
                        || (this.matches(tmpCompiledRule.leftEnvironmentIndex(), tmpEnd, tmpAtomCount, tmpEnvironmentMatches)
                        && this.matches(tmpCompiledRule.rightEnvironmentIndex(), tmpBegin, tmpAtomCount, tmpEnvironmentMatches));
                if (tmpIsMatch) {
                    tmpRECAPBonds.add(new RECAPBond(tmpBond, tmpCompiledRule.rule().bondType()));
                    break;
                }
            }
        }
        return tmpRECAPBonds;
    }
    //
    /**
     * Selects the bonds to cleave from the given cleavable bonds, so that no small terminal group is cleaved off. The
     * bonds are processed in the order of their bond type (and bond index); a bond is selected unless one of the two
     * parts that result from its cleavage, together with the bonds selected before, is a small terminal group, see
     * {@link #isSmallTerminalGroup(IAtomContainer, int, int, int[][], GraphUtil.EdgeToBondMap, boolean[], int)}.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @param aRECAPBondList cleavable bonds of the molecule
     * @return the bonds to cleave
     */
    List<RECAPBond> selectBondsToCleave(IAtomContainer aPreparedMolecule, List<RECAPBond> aRECAPBondList) {
        if (aRECAPBondList.isEmpty()) {
            return aRECAPBondList;
        }
        List<RECAPBond> tmpCandidates = new ArrayList<>(aRECAPBondList);
        tmpCandidates.sort(RECAPFragmenter.BOND_PROCESSING_ORDER);
        GraphUtil.EdgeToBondMap tmpEdgeToBondMap = GraphUtil.EdgeToBondMap.withSpaceFor(aPreparedMolecule);
        int[][] tmpAdjacencyList = GraphUtil.toAdjList(aPreparedMolecule, tmpEdgeToBondMap);
        boolean[] tmpIsCleaved = new boolean[aPreparedMolecule.getBondCount()];
        int tmpMaximumCarbonAtoms = this.smallGroupMaximumCarbonAtomsSetting.get();
        List<RECAPBond> tmpSelectedBonds = new ArrayList<>(tmpCandidates.size());
        for (RECAPBond tmpCandidate : tmpCandidates) {
            int tmpBondIndex = tmpCandidate.bond().getIndex();
            tmpIsCleaved[tmpBondIndex] = true;
            if (RECAPFragmenter.isSmallTerminalGroup(aPreparedMolecule, tmpCandidate.bond().getBegin().getIndex(),
                    tmpCandidate.bondType(), tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved, tmpMaximumCarbonAtoms)
                    || RECAPFragmenter.isSmallTerminalGroup(aPreparedMolecule, tmpCandidate.bond().getEnd().getIndex(),
                    tmpCandidate.bondType(), tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved, tmpMaximumCarbonAtoms)) {
                tmpIsCleaved[tmpBondIndex] = false;
            } else {
                tmpSelectedBonds.add(tmpCandidate);
            }
        }
        return tmpSelectedBonds;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Checks whether the order of a bond of the given molecule is unset and cannot be determined by kekulization in the
     * preparation of the molecule, see {@link #prepareMolecule(IAtomContainer)}. Molecules without bonds of unset
     * order are not cloned; otherwise, a clone of the molecule is prepared, so that the check and the fragmentation
     * cannot deviate.
     *
     * @param aMolecule the molecule to check, it is not altered
     * @return true if the order of a bond is unset (or null) and the preparation of a clone fails, e.g. because no
     *         Kekulé structure can be assigned, or if the molecule cannot be cloned
     */
    private boolean hasUnresolvableBondOrders(IAtomContainer aMolecule) {
        boolean tmpIsBondOrderUnset = false;
        for (IBond tmpBond : aMolecule.bonds()) {
            if (Objects.isNull(tmpBond.getOrder()) || tmpBond.getOrder() == IBond.Order.UNSET) {
                tmpIsBondOrderUnset = true;
                break;
            }
        }
        if (!tmpIsBondOrderUnset) {
            return false;
        }
        try {
            this.prepareMolecule(aMolecule.clone());
            return false;
        } catch (CDKException | CloneNotSupportedException anException) {
            return true;
        }
    }
    //
    /**
     * Returns whether the given atom matches the environment pattern with the given index; the result is cached.
     *
     * @param anEnvironmentIndex index of the environment pattern
     * @param anAtom atom of the prepared molecule
     * @param anAtomCount number of atoms of the prepared molecule
     * @param aCache cache of the matching results per environment and atom index
     * @return true if the atom matches the environment
     */
    private boolean matches(int anEnvironmentIndex, IAtom anAtom, int anAtomCount, byte[][] aCache) {
        byte[] tmpResults = aCache[anEnvironmentIndex];
        if (tmpResults == null) {
            tmpResults = new byte[anAtomCount];
            aCache[anEnvironmentIndex] = tmpResults;
        }
        int tmpIndex = anAtom.getIndex();
        if (tmpResults[tmpIndex] == 0) {
            tmpResults[tmpIndex] = this.environmentPatterns[anEnvironmentIndex].matchesRoot(anAtom) ? (byte) 1 : (byte) 2;
        }
        return tmpResults[tmpIndex] == 1;
    }
    //
    /**
     * Determines the part of the molecule that contains the given atom of a cleaved bond (the cleavage point) when the
     * marked bonds are cleaved, and returns whether it is a small terminal group that must not be cleaved off. The
     * part must be terminal, i.e. have exactly one attachment point (cleaved bond or pseudo atom of the input), and
     * <ul>
     *     <li>if the cleavage point is a nitrogen atom at a bond of the types 3 (amine), 7, 8 or 9: consist of the
     *     nitrogen atom with hydrogen atoms only (amino group); nitrogen atoms with alkyl groups, e.g. dimethylamino,
     *     isopropylamino or tert-butylamino, are cleaved off as in the publication (Figures 9 and 10);</li>
     *     <li>if the cleavage point is a carbon atom without double-bonded oxygen atom: be an acyclic hydrocarbon group
     *     with at most the given number of carbon atoms (methyl, ethyl, propyl, butyl);</li>
     *     <li>otherwise, i.e. for oxygen atoms, carbonyl or sulfonyl groups, and nitrogen atoms at a bond of the types
     *     1 (amide), 4 (urea) or 11 (sulphonamide): consist of the cleavage point (with its double-bonded oxygen
     *     atoms) and hydrogen atoms or an acyclic hydrocarbon group with at most the given number of carbon atoms in
     *     total (e.g. methoxy, ethoxy, formyl, acetyl, mesyl, NHCH3 or N(CH2CH3)2 of an amide, as the N-methyl amide in
     *     Figure 12 of the publication).</li>
     * </ul>
     * Pseudo atoms of the input count as attachment points. Explicit hydrogen atoms, i.e. hydrogen isotopes that are
     * not suppressed, are ignored like implicit ones.
     *
     * @param aMolecule the prepared molecule
     * @param aCleavagePointIndex index of the atom of the cleaved bond in the part
     * @param aBondType RECAP bond type of the cleaved bond
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @param anIsCleavedArray cleaved bonds by bond index
     * @param aMaximumCarbonAtoms maximum number of carbon atoms of a small terminal group
     * @return true if the part is a small terminal group
     */
    private static boolean isSmallTerminalGroup(IAtomContainer aMolecule, int aCleavagePointIndex, int aBondType,
                                                int[][] anAdjacencyList, GraphUtil.EdgeToBondMap anEdgeToBondMap,
                                                boolean[] anIsCleavedArray, int aMaximumCarbonAtoms) {
        IAtom tmpCleavagePoint = aMolecule.getAtom(aCleavagePointIndex);
        //the cleavage point and the oxygen atoms double-bonded to it (carbonyl, sulfonyl) form the "head" of the group
        Set<Integer> tmpHeadAtoms = new HashSet<>(4);
        tmpHeadAtoms.add(aCleavagePointIndex);
        for (int tmpNeighbour : anAdjacencyList[aCleavagePointIndex]) {
            IBond tmpBond = anEdgeToBondMap.get(aCleavagePointIndex, tmpNeighbour);
            IAtom tmpNeighbourAtom = aMolecule.getAtom(tmpNeighbour);
            if (!anIsCleavedArray[tmpBond.getIndex()] && tmpBond.getOrder() == IBond.Order.DOUBLE
                    && Objects.equals(tmpNeighbourAtom.getAtomicNumber(), 8)) {
                tmpHeadAtoms.add(tmpNeighbour);
            }
        }
        //the nitrogen atom of an amide, urea or sulphonamide bond is a head like an oxygen atom, i.e. small alkyl
        //groups at it are kept (N-methyl amide in Figure 12), whereas at the other bond types, especially amines, only
        //an amino group is small and alkylated amino groups are cleaved off (Figures 9 and 10)
        boolean tmpIsAminoNitrogen = Objects.equals(tmpCleavagePoint.getAtomicNumber(), 7)
                && aBondType != 1 && aBondType != 4 && aBondType != 11;
        boolean tmpIsAlkylCarbon = Objects.equals(tmpCleavagePoint.getAtomicNumber(), 6) && tmpHeadAtoms.size() == 1;
        //a small group has at most the head atoms and the carbon atoms of the alkyl group
        int tmpMaximumPartSize = tmpHeadAtoms.size() + aMaximumCarbonAtoms;
        boolean[] tmpIsVisited = new boolean[anAdjacencyList.length];
        Deque<Integer> tmpQueue = new ArrayDeque<>();
        tmpQueue.add(aCleavagePointIndex);
        tmpIsVisited[aCleavagePointIndex] = true;
        int tmpPartSize = 0;
        int tmpAttachmentPoints = 0;
        //atoms of the group except the head atoms (for an alkyl carbon cleavage point, including it)
        int tmpTailCarbonAtoms = 0;
        boolean tmpIsTailHydrocarbon = true;
        while (!tmpQueue.isEmpty()) {
            int tmpAtomIndex = tmpQueue.poll();
            IAtom tmpAtom = aMolecule.getAtom(tmpAtomIndex);
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAttachmentPoints++;
            } else if (!Objects.equals(tmpAtom.getAtomicNumber(), 1)) {
                //explicit hydrogen atoms (isotopes that are not suppressed) are ignored like implicit ones
                tmpPartSize++;
                if (tmpIsAlkylCarbon || !tmpHeadAtoms.contains(tmpAtomIndex)) {
                    if (Objects.equals(tmpAtom.getAtomicNumber(), 6) && !tmpAtom.isInRing()) {
                        tmpTailCarbonAtoms++;
                    } else {
                        tmpIsTailHydrocarbon = false;
                    }
                }
            }
            for (int tmpNeighbour : anAdjacencyList[tmpAtomIndex]) {
                if (anIsCleavedArray[anEdgeToBondMap.get(tmpAtomIndex, tmpNeighbour).getIndex()]) {
                    tmpAttachmentPoints++;
                } else if (!tmpIsVisited[tmpNeighbour]) {
                    tmpIsVisited[tmpNeighbour] = true;
                    tmpQueue.add(tmpNeighbour);
                }
            }
            //the part is not a small terminal group, no further search needed
            if (tmpAttachmentPoints > 1 || !tmpIsTailHydrocarbon || tmpPartSize > tmpMaximumPartSize) {
                return false;
            }
        }
        if (tmpAttachmentPoints != 1) {
            return false;
        }
        if (tmpIsAminoNitrogen) {
            return tmpPartSize == tmpHeadAtoms.size();
        }
        return tmpTailCarbonAtoms <= aMaximumCarbonAtoms;
    }
    //
    /**
     * Creates the RECAP hierarchy of the given prepared molecule, see {@link #buildHierarchy(IAtomContainer)}. The
     * nodes are first determined as atom sets, with the bonds connecting an atom set to the rest of the molecule
     * regarded as cleaved, and the size of the hierarchy is checked while they are added; then their fragments are
     * created.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}; it becomes the fragment
     *                          of the root and must not be altered afterwards
     * @param aRECAPBondList all cleavable bonds of the molecule
     * @return all nodes of the hierarchy in breadth-first order, the first one is the root
     * @throws IllegalArgumentException if the hierarchy has more than {@link #MAXIMUM_HIERARCHY_NODE_COUNT} nodes or
     *                                  more than {@link #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in all nodes together
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    private List<RECAPHierarchyNode> createHierarchy(IAtomContainer aPreparedMolecule, List<RECAPBond> aRECAPBondList)
            throws IllegalArgumentException, CloneNotSupportedException {
        int tmpAtomCount = aPreparedMolecule.getAtomCount();
        int tmpBondCount = aPreparedMolecule.getBondCount();
        List<RECAPBond> tmpCandidates = new ArrayList<>(aRECAPBondList);
        tmpCandidates.sort(RECAPFragmenter.BOND_PROCESSING_ORDER);
        int[] tmpBondTypes = new int[tmpBondCount];
        for (RECAPBond tmpCandidate : tmpCandidates) {
            tmpBondTypes[tmpCandidate.bond().getIndex()] = tmpCandidate.bondType();
        }
        GraphUtil.EdgeToBondMap tmpEdgeToBondMap = GraphUtil.EdgeToBondMap.withSpaceFor(aPreparedMolecule);
        int[][] tmpAdjacencyList = GraphUtil.toAdjList(aPreparedMolecule, tmpEdgeToBondMap);
        int tmpMaximumCarbonAtoms = this.smallGroupMaximumCarbonAtomsSetting.get();
        List<BitSet> tmpAtomSets = new ArrayList<>();
        Map<BitSet, Integer> tmpNodeIndexMap = new HashMap<>();
        List<List<Integer>> tmpChildIndexLists = new ArrayList<>();
        BitSet tmpRootAtomSet = new BitSet(tmpAtomCount);
        tmpRootAtomSet.set(0, tmpAtomCount);
        RECAPFragmenter.addHierarchyNode(tmpRootAtomSet, tmpAtomSets, tmpNodeIndexMap, tmpChildIndexLists);
        //total number of atoms of the nodes whose atoms have been counted, i.e. of the first tmpCountedNodeCount nodes
        long tmpHierarchyAtomCount = 0;
        int tmpCountedNodeCount = 0;
        ConnectedComponents tmpConnectedComponents = new ConnectedComponents(tmpAdjacencyList);
        int[] tmpComponents = tmpConnectedComponents.components();
        boolean[] tmpIsComponentCleaved = new boolean[tmpConnectedComponents.nComponents() + 1];
        //breadth-first: the list of atom sets grows while it is processed
        for (int tmpNodeIndex = 0; tmpNodeIndex < tmpAtomSets.size(); tmpNodeIndex++) {
            BitSet tmpAtomSet = tmpAtomSets.get(tmpNodeIndex);
            //the bonds that connect the fragment to the rest of the molecule have been cleaved
            boolean[] tmpIsCleaved = new boolean[tmpBondCount];
            for (IBond tmpBond : aPreparedMolecule.bonds()) {
                tmpIsCleaved[tmpBond.getIndex()] = tmpAtomSet.get(tmpBond.getBegin().getIndex())
                        != tmpAtomSet.get(tmpBond.getEnd().getIndex());
            }
            for (RECAPBond tmpCandidate : tmpCandidates) {
                int tmpBeginIndex = tmpCandidate.bond().getBegin().getIndex();
                int tmpEndIndex = tmpCandidate.bond().getEnd().getIndex();
                if (!tmpAtomSet.get(tmpBeginIndex) || !tmpAtomSet.get(tmpEndIndex)) {
                    continue;
                }
                int tmpBondIndex = tmpCandidate.bond().getIndex();
                tmpIsCleaved[tmpBondIndex] = true;
                if (!RECAPFragmenter.isSmallTerminalGroup(aPreparedMolecule, tmpBeginIndex, tmpCandidate.bondType(),
                        tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved, tmpMaximumCarbonAtoms)
                        && !RECAPFragmenter.isSmallTerminalGroup(aPreparedMolecule, tmpEndIndex,
                        tmpCandidate.bondType(), tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved,
                        tmpMaximumCarbonAtoms)) {
                    tmpIsComponentCleaved[tmpComponents[tmpBeginIndex]] = true;
                    for (int tmpStartIndex : new int[] {tmpBeginIndex, tmpEndIndex}) {
                        BitSet tmpPart = RECAPFragmenter.collectPart(tmpStartIndex, tmpAtomCount, tmpAdjacencyList,
                                tmpEdgeToBondMap, tmpIsCleaved);
                        tmpChildIndexLists.get(tmpNodeIndex).add(RECAPFragmenter.addHierarchyNode(tmpPart, tmpAtomSets,
                                tmpNodeIndexMap, tmpChildIndexLists));
                    }
                    tmpHierarchyAtomCount = RECAPFragmenter.checkHierarchySize(tmpAtomSets, tmpCountedNodeCount,
                            tmpHierarchyAtomCount);
                    tmpCountedNodeCount = tmpAtomSets.size();
                }
                tmpIsCleaved[tmpBondIndex] = false;
            }
        }
        //disconnected parts without a cleavable bond, e.g. counter ions, are children of the root (and leaves)
        if (tmpConnectedComponents.nComponents() > 1) {
            for (int tmpComponent = 1; tmpComponent <= tmpConnectedComponents.nComponents(); tmpComponent++) {
                if (tmpIsComponentCleaved[tmpComponent]) {
                    continue;
                }
                BitSet tmpPart = new BitSet(tmpAtomCount);
                for (int i = 0; i < tmpAtomCount; i++) {
                    if (tmpComponents[i] == tmpComponent) {
                        tmpPart.set(i);
                    }
                }
                tmpChildIndexLists.get(0).add(RECAPFragmenter.addHierarchyNode(tmpPart, tmpAtomSets, tmpNodeIndexMap,
                        tmpChildIndexLists));
            }
        }
        RECAPFragmenter.checkHierarchySize(tmpAtomSets, tmpCountedNodeCount, tmpHierarchyAtomCount);
        //the stereogenicity of the atoms of the molecule is determined once for all nodes, by atom index, since the
        //fragment of every node is created from a new clone of the molecule
        BitSet tmpNonStereogenicAtomIndices = aPreparedMolecule.stereoElements().iterator().hasNext()
                ? RECAPFragmenter.findNonStereogenicAtomIndices(aPreparedMolecule) : new BitSet(0);
        List<RECAPHierarchyNode> tmpNodes = new ArrayList<>(tmpAtomSets.size());
        //placeholder for the root, whose fragment is the prepared molecule that must remain unaltered for the clones
        tmpNodes.add(null);
        for (int i = 1; i < tmpAtomSets.size(); i++) {
            tmpNodes.add(new RECAPHierarchyNode(this.createHierarchyFragment(aPreparedMolecule, tmpAtomSets.get(i),
                    tmpBondTypes, tmpNonStereogenicAtomIndices)));
        }
        tmpNodes.set(0, new RECAPHierarchyNode(aPreparedMolecule));
        for (int i = 0; i < tmpNodes.size(); i++) {
            for (int tmpChildIndex : tmpChildIndexLists.get(i)) {
                tmpNodes.get(i).children.add(tmpNodes.get(tmpChildIndex));
                tmpNodes.get(tmpChildIndex).parents.add(tmpNodes.get(i));
            }
        }
        return tmpNodes;
    }
    //
    /**
     * Returns the index of the hierarchy node with the given atom set, adding a new node if there is none.
     *
     * @param anAtomSet atom set of the node
     * @param anAtomSetList atom sets of all nodes by node index
     * @param aNodeIndexMap map of the atom sets to the node indices
     * @param aChildIndexLists child node indices of all nodes by node index
     * @return index of the node
     */
    private static int addHierarchyNode(BitSet anAtomSet, List<BitSet> anAtomSetList, Map<BitSet, Integer> aNodeIndexMap,
                                        List<List<Integer>> aChildIndexLists) {
        Integer tmpIndex = aNodeIndexMap.get(anAtomSet);
        if (tmpIndex != null) {
            return tmpIndex;
        }
        int tmpNewIndex = anAtomSetList.size();
        anAtomSetList.add(anAtomSet);
        aNodeIndexMap.put(anAtomSet, tmpNewIndex);
        aChildIndexLists.add(new ArrayList<>(2));
        return tmpNewIndex;
    }
    //
    /**
     * Adds the atoms of the hierarchy nodes that have not been counted yet to the total number of atoms of the
     * hierarchy and checks the size of the hierarchy against the limits {@link #MAXIMUM_HIERARCHY_NODE_COUNT} and
     * {@link #MAXIMUM_HIERARCHY_ATOM_COUNT}, which bound the memory needed for the fragments of the nodes.
     *
     * @param anAtomSetList atom sets of all nodes by node index
     * @param aCountedNodeCount number of nodes whose atoms have already been counted (the first nodes of the list)
     * @param aCountedAtomCount total number of atoms of these nodes
     * @return total number of atoms of all nodes of the list
     * @throws IllegalArgumentException if the hierarchy has more than {@link #MAXIMUM_HIERARCHY_NODE_COUNT} nodes or
     *                                  more than {@link #MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in all nodes together
     */
    private static long checkHierarchySize(List<BitSet> anAtomSetList, int aCountedNodeCount, long aCountedAtomCount)
            throws IllegalArgumentException {
        if (anAtomSetList.size() > RECAPFragmenter.MAXIMUM_HIERARCHY_NODE_COUNT) {
            throw new IllegalArgumentException(String.format(
                    "The RECAP hierarchy of the molecule has more than %d nodes.",
                    RECAPFragmenter.MAXIMUM_HIERARCHY_NODE_COUNT));
        }
        long tmpAtomCount = aCountedAtomCount;
        for (int i = aCountedNodeCount; i < anAtomSetList.size(); i++) {
            tmpAtomCount += anAtomSetList.get(i).cardinality();
        }
        if (tmpAtomCount > RECAPFragmenter.MAXIMUM_HIERARCHY_ATOM_COUNT) {
            throw new IllegalArgumentException(String.format(
                    "The fragments of the RECAP hierarchy of the molecule have more than %d atoms in total.",
                    RECAPFragmenter.MAXIMUM_HIERARCHY_ATOM_COUNT));
        }
        return tmpAtomCount;
    }
    //
    /**
     * Returns the atoms that are connected to the given start atom without crossing a cleaved bond.
     *
     * @param aStartIndex index of the start atom
     * @param anAtomCount number of atoms of the molecule
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @param anIsCleavedArray cleaved bonds by bond index
     * @return the atom indices of the part
     */
    private static BitSet collectPart(int aStartIndex, int anAtomCount, int[][] anAdjacencyList,
                                      GraphUtil.EdgeToBondMap anEdgeToBondMap, boolean[] anIsCleavedArray) {
        BitSet tmpPart = new BitSet(anAtomCount);
        Deque<Integer> tmpQueue = new ArrayDeque<>();
        tmpQueue.add(aStartIndex);
        tmpPart.set(aStartIndex);
        while (!tmpQueue.isEmpty()) {
            int tmpAtomIndex = tmpQueue.poll();
            for (int tmpNeighbour : anAdjacencyList[tmpAtomIndex]) {
                if (!tmpPart.get(tmpNeighbour)
                        && !anIsCleavedArray[anEdgeToBondMap.get(tmpAtomIndex, tmpNeighbour).getIndex()]) {
                    tmpPart.set(tmpNeighbour);
                    tmpQueue.add(tmpNeighbour);
                }
            }
        }
        return tmpPart;
    }
    //
    /**
     * Creates the fragment of a hierarchy node: a clone of the prepared molecule is cleaved at all bonds that connect
     * the given atom set to the rest of the molecule, and the part with the atom set is returned.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}, it is not altered
     * @param anAtomSet atom indices of the fragment
     * @param aBondTypes RECAP bond types by bond index (0 for bonds that are not cleavable)
     * @param aNonStereogenicAtomIndices indices of the atoms that are not stereogenic in the prepared molecule (empty
     *                                   if it has no stereo elements), or null if they could not be determined
     * @return the fragment
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    private IAtomContainer createHierarchyFragment(IAtomContainer aPreparedMolecule, BitSet anAtomSet, int[] aBondTypes,
                                                   BitSet aNonStereogenicAtomIndices)
            throws CloneNotSupportedException {
        IAtomContainer tmpClone = aPreparedMolecule.clone();
        List<RECAPBond> tmpBoundaryBonds = new ArrayList<>();
        for (IBond tmpBond : tmpClone.bonds()) {
            if (anAtomSet.get(tmpBond.getBegin().getIndex()) != anAtomSet.get(tmpBond.getEnd().getIndex())) {
                tmpBoundaryBonds.add(new RECAPBond(tmpBond, aBondTypes[tmpBond.getIndex()]));
            }
        }
        //the clone has the atom order of the prepared molecule; only the atoms of the fragment are needed
        Set<IAtom> tmpNonStereogenicAtoms = null;
        if (aNonStereogenicAtomIndices != null) {
            tmpNonStereogenicAtoms = new HashSet<>();
            for (int i = aNonStereogenicAtomIndices.nextSetBit(0); i >= 0;
                 i = aNonStereogenicAtomIndices.nextSetBit(i + 1)) {
                if (anAtomSet.get(i)) {
                    tmpNonStereogenicAtoms.add(tmpClone.getAtom(i));
                }
            }
        }
        IAtom tmpMarkerAtom = tmpClone.getAtom(anAtomSet.nextSetBit(0));
        List<IAtomContainer> tmpFragments = this.cleaveBonds(tmpClone, tmpBoundaryBonds, tmpNonStereogenicAtoms,
                tmpMarkerAtom);
        if (tmpFragments.size() != 1) {
            throw new IllegalStateException("The fragment of a RECAP hierarchy node could not be determined.");
        }
        return tmpFragments.getFirst();
    }
    //
    /**
     * Cleaves the given bonds of the prepared molecule simultaneously, saturates the open valences according to the
     * saturation setting, transfers the stereo elements, and returns the disconnected parts as fragments, see
     * {@link #cleaveBonds(IAtomContainer, List, Set, IAtom)}.
     *
     * @param aMolecule the prepared molecule, it is altered
     * @param aRECAPBondList the bonds to cleave
     * @return list of fragments
     */
    private List<IAtomContainer> cleaveBonds(IAtomContainer aMolecule, List<RECAPBond> aRECAPBondList) {
        //the stereogenicity of the atoms is determined before the cleavage
        Set<IAtom> tmpNonStereogenicAtomsOfMolecule = aMolecule.stereoElements().iterator().hasNext()
                ? RECAPFragmenter.findNonStereogenicAtoms(aMolecule) : Collections.emptySet();
        return this.cleaveBonds(aMolecule, aRECAPBondList, tmpNonStereogenicAtomsOfMolecule, null);
    }
    //
    /**
     * Cleaves the given bonds of the prepared molecule simultaneously, saturates the open valences according to the
     * saturation setting, transfers the stereo elements, and returns the disconnected parts as fragments; stereo
     * elements that lost their stereogenicity by the cleavage are removed from the returned fragments. This is the
     * only implementation of the cleavage, used for the single-pass fragments and for the fragments of the hierarchy.
     *
     * @param aMolecule the prepared molecule, it is altered
     * @param aRECAPBondList the bonds to cleave
     * @param aNonStereogenicAtomSetOfMolecule atoms of the molecule that are not stereogenic before the cleavage (at
     *                                         least those of the returned fragments; empty if the molecule has no
     *                                         stereo elements), or null if they could not be determined
     * @param aPartAtom if not null, only the fragment that contains this atom is returned
     * @return list of fragments
     */
    @SuppressWarnings("rawtypes")
    private List<IAtomContainer> cleaveBonds(IAtomContainer aMolecule, List<RECAPBond> aRECAPBondList,
                                             Set<IAtom> aNonStereogenicAtomSetOfMolecule, IAtom aPartAtom) {
        boolean tmpIsHydrogenSaturation = this.saturationSetting.get().equals(SaturationOption.HYDROGEN_SATURATION);
        //the stereo elements are detached before the bonds are removed, because the CDK atom container would otherwise
        //adapt them on bond removal, e.g. by replacing a cleaved ligand of a stereocentre with an implicit hydrogen atom
        List<IStereoElement> tmpOriginalStereoElements = new ArrayList<>();
        for (IStereoElement tmpElement : aMolecule.stereoElements()) {
            tmpOriginalStereoElements.add(tmpElement);
        }
        aMolecule.setStereoElements(new ArrayList<>(0));
        Map<IBond, CleavedBond> tmpCleavedBondMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                aRECAPBondList.size(), BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR),
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        for (RECAPBond tmpRECAPBond : aRECAPBondList) {
            IBond tmpBond = tmpRECAPBond.bond();
            IAtom tmpBegin = tmpBond.getBegin();
            IAtom tmpEnd = tmpBond.getEnd();
            IBond.Order tmpOrder = tmpBond.getOrder();
            aMolecule.removeBond(tmpBond);
            IBond tmpNewBeginBond = this.saturate(aMolecule, tmpBegin, tmpEnd, tmpRECAPBond.bondType(), tmpOrder,
                    tmpIsHydrogenSaturation);
            IBond tmpNewEndBond = this.saturate(aMolecule, tmpEnd, tmpBegin, tmpRECAPBond.bondType(), tmpOrder,
                    tmpIsHydrogenSaturation);
            tmpCleavedBondMap.put(tmpBond, new CleavedBond(tmpBegin, tmpEnd, tmpNewBeginBond, tmpNewEndBond));
        }
        this.transferStereoElements(aMolecule, tmpOriginalStereoElements, tmpCleavedBondMap, tmpIsHydrogenSaturation);
        if (tmpIsHydrogenSaturation) {
            //the new explicit hydrogen atoms become implicit ones, stereo elements are updated accordingly
            AtomContainerManipulator.suppressHydrogens(aMolecule);
        }
        IAtomContainerSet tmpFragmentSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        List<IAtomContainer> tmpFragments = new ArrayList<>(tmpFragmentSet.getAtomContainerCount());
        for (IAtomContainer tmpFragment : tmpFragmentSet.atomContainers()) {
            if (aPartAtom != null && !tmpFragment.contains(aPartAtom)) {
                continue;
            }
            RECAPFragmenter.removeStereoElementsMadeNonStereogenic(tmpFragment, aNonStereogenicAtomSetOfMolecule);
            tmpFragments.add(tmpFragment);
        }
        return tmpFragments;
    }
    //
    /**
     * Saturates the open valence of the given atom that results from the cleavage of its bond to the given former
     * neighbour, either with a link atom or with a hydrogen atom (explicit for single bonds, so that stereo elements
     * can refer to it, two implicit ones for double bonds).
     *
     * @param aMolecule the molecule
     * @param anAtom atom of the cleaved bond
     * @param aFormerNeighbour the other atom of the cleaved bond
     * @param aBondType RECAP bond type of the cleaved bond
     * @param anOrder order of the cleaved bond
     * @param isHydrogenSaturation true if the open valence should be saturated with hydrogen atoms
     * @return the new bond of the atom to its link atom or explicit hydrogen atom, null for implicit hydrogen atoms
     */
    private IBond saturate(IAtomContainer aMolecule, IAtom anAtom, IAtom aFormerNeighbour, int aBondType,
                           IBond.Order anOrder, boolean isHydrogenSaturation) {
        IChemObjectBuilder tmpBuilder = aMolecule.getBuilder();
        IAtom tmpNewAtom;
        if (isHydrogenSaturation) {
            if (anOrder != IBond.Order.SINGLE) {
                anAtom.setImplicitHydrogenCount(anAtom.getImplicitHydrogenCount() + anOrder.numeric());
                return null;
            }
            tmpNewAtom = tmpBuilder.newInstance(IAtom.class, "H");
        } else {
            IPseudoAtom tmpLinkAtom = tmpBuilder.newInstance(IPseudoAtom.class, "*");
            tmpLinkAtom.setMassNumber(aBondType);
            tmpLinkAtom.setProperty(RECAPFragmenter.BOND_TYPE_PROPERTY_KEY, aBondType);
            tmpNewAtom = tmpLinkAtom;
        }
        tmpNewAtom.setImplicitHydrogenCount(0);
        //the new atom takes the position of the former neighbour
        if (aFormerNeighbour.getPoint2d() != null) {
            tmpNewAtom.setPoint2d(new Point2d(aFormerNeighbour.getPoint2d()));
        }
        if (aFormerNeighbour.getPoint3d() != null) {
            tmpNewAtom.setPoint3d(new Point3d(aFormerNeighbour.getPoint3d()));
        }
        aMolecule.addAtom(tmpNewAtom);
        IAtom tmpAddedAtom = aMolecule.getAtom(aMolecule.getAtomCount() - 1);
        aMolecule.addBond(tmpBuilder.newInstance(IBond.class, anAtom, tmpAddedAtom, anOrder));
        return aMolecule.getBond(aMolecule.getBondCount() - 1);
    }
    //
    /**
     * Transfers the stereo elements of the molecule before the cleavage to the molecule after the cleavage of the
     * bonds. Ligands of tetrahedral stereocentres and of stereo double bonds that were cleaved off are replaced by the
     * link atom (or explicit hydrogen atom) at their former position. Stereo elements whose focus bond was cleaved,
     * tetrahedral stereocentres that obtain two hydrogen atoms, and stereo elements that would span more than one
     * fragment are removed.
     *
     * @param aMolecule the molecule after the cleavage of the bonds, its stereo elements are set
     * @param anOriginalStereoElementList the stereo elements of the molecule before the cleavage
     * @param aCleavedBondMap map of the cleaved bonds to their cleavage results
     * @param isHydrogenSaturation true if open valences are saturated with hydrogen atoms
     */
    @SuppressWarnings("rawtypes")
    private void transferStereoElements(IAtomContainer aMolecule, List<IStereoElement> anOriginalStereoElementList,
                                        Map<IBond, CleavedBond> aCleavedBondMap, boolean isHydrogenSaturation) {
        if (anOriginalStereoElementList.isEmpty()) {
            return;
        }
        Map<IAtom, List<CleavedBond>> tmpAtomToCleavedBondsMap = new HashMap<>();
        for (CleavedBond tmpCleavedBond : aCleavedBondMap.values()) {
            tmpAtomToCleavedBondsMap.computeIfAbsent(tmpCleavedBond.begin(), aKey -> new ArrayList<>(2)).add(tmpCleavedBond);
            tmpAtomToCleavedBondsMap.computeIfAbsent(tmpCleavedBond.end(), aKey -> new ArrayList<>(2)).add(tmpCleavedBond);
        }
        int[] tmpComponents = new ConnectedComponents(GraphUtil.toAdjList(aMolecule)).components();
        List<IStereoElement> tmpTransferredElements = new ArrayList<>(anOriginalStereoElementList.size());
        for (IStereoElement tmpElement : anOriginalStereoElementList) {
            IStereoElement tmpTransferredElement;
            if (tmpElement.getConfigClass() == IStereoElement.TH) {
                tmpTransferredElement = this.transferTetrahedralElement(tmpElement, tmpAtomToCleavedBondsMap,
                        isHydrogenSaturation);
            } else if (tmpElement.getConfigClass() == IStereoElement.CT) {
                tmpTransferredElement = this.transferDoubleBondElement(tmpElement, aCleavedBondMap,
                        isHydrogenSaturation);
            } else {
                //other stereo element types are kept only if they are not affected by a cleavage, see below
                tmpTransferredElement = tmpElement;
            }
            if (tmpTransferredElement != null
                    && RECAPFragmenter.isInOneComponent(aMolecule, tmpTransferredElement, aCleavedBondMap.keySet(), tmpComponents)) {
                tmpTransferredElements.add(tmpTransferredElement);
            }
        }
        aMolecule.setStereoElements(tmpTransferredElements);
    }
    //
    /**
     * Transfers a tetrahedral stereo element, see
     * {@link #transferStereoElements(IAtomContainer, List, Map, boolean)}.
     *
     * @param anElement tetrahedral stereo element
     * @param anAtomToCleavedBondsMap map of atoms to the cleaved bonds they were part of
     * @param isHydrogenSaturation true if open valences are saturated with hydrogen atoms
     * @return the transferred element or null if it is removed
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private IStereoElement transferTetrahedralElement(IStereoElement anElement,
                                                      Map<IAtom, List<CleavedBond>> anAtomToCleavedBondsMap,
                                                      boolean isHydrogenSaturation) {
        IAtom tmpFocus = (IAtom) anElement.getFocus();
        List<CleavedBond> tmpCleavedBonds = anAtomToCleavedBondsMap.get(tmpFocus);
        if (tmpCleavedBonds == null) {
            return anElement;
        }
        List<IAtom> tmpCarriers = (List<IAtom>) anElement.getCarriers();
        Map<IAtom, IAtom> tmpAtomMap = new HashMap<>();
        tmpAtomMap.put(tmpFocus, tmpFocus);
        for (IAtom tmpCarrier : tmpCarriers) {
            tmpAtomMap.put(tmpCarrier, tmpCarrier);
        }
        int tmpHydrogenCount = Objects.requireNonNullElse(tmpFocus.getImplicitHydrogenCount(), 0);
        for (CleavedBond tmpCleavedBond : tmpCleavedBonds) {
            IAtom tmpFormerLigand = tmpCleavedBond.other(tmpFocus);
            IBond tmpNewBond = tmpCleavedBond.newBondOf(tmpFocus);
            if (!tmpCarriers.contains(tmpFormerLigand) || tmpNewBond == null) {
                return null;
            }
            tmpAtomMap.put(tmpFormerLigand, tmpNewBond.getOther(tmpFocus));
            if (isHydrogenSaturation) {
                tmpHydrogenCount++;
            }
        }
        if (tmpHydrogenCount > 1) {
            return null;
        }
        return anElement.map(tmpAtomMap, Collections.emptyMap());
    }
    //
    /**
     * Transfers a double bond stereo element, see
     * {@link #transferStereoElements(IAtomContainer, List, Map, boolean)}.
     *
     * @param anElement double bond stereo element
     * @param aCleavedBondMap map of the cleaved bonds to their cleavage results
     * @param isHydrogenSaturation true if open valences are saturated with hydrogen atoms
     * @return the transferred element or null if it is removed
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private IStereoElement transferDoubleBondElement(IStereoElement anElement, Map<IBond, CleavedBond> aCleavedBondMap,
                                                     boolean isHydrogenSaturation) {
        IBond tmpFocus = (IBond) anElement.getFocus();
        if (aCleavedBondMap.containsKey(tmpFocus)) {
            return null;
        }
        List<IBond> tmpCarriers = (List<IBond>) anElement.getCarriers();
        Map<IBond, IBond> tmpBondMap = new HashMap<>();
        tmpBondMap.put(tmpFocus, tmpFocus);
        boolean tmpIsAffected = false;
        for (IBond tmpCarrier : tmpCarriers) {
            CleavedBond tmpCleavedBond = aCleavedBondMap.get(tmpCarrier);
            if (tmpCleavedBond == null) {
                tmpBondMap.put(tmpCarrier, tmpCarrier);
                continue;
            }
            tmpIsAffected = true;
            IAtom tmpDoubleBondAtom = tmpFocus.contains(tmpCleavedBond.begin()) ? tmpCleavedBond.begin() : tmpCleavedBond.end();
            IBond tmpNewBond = tmpCleavedBond.newBondOf(tmpDoubleBondAtom);
            //with hydrogen saturation, a double bond atom that already carries a hydrogen atom is not stereogenic anymore
            if (tmpNewBond == null
                    || (isHydrogenSaturation && Objects.requireNonNullElse(tmpDoubleBondAtom.getImplicitHydrogenCount(), 0) > 0)) {
                return null;
            }
            tmpBondMap.put(tmpCarrier, tmpNewBond);
        }
        if (!tmpIsAffected) {
            return anElement;
        }
        Map<IAtom, IAtom> tmpAtomMap = new HashMap<>();
        tmpAtomMap.put(tmpFocus.getBegin(), tmpFocus.getBegin());
        tmpAtomMap.put(tmpFocus.getEnd(), tmpFocus.getEnd());
        return anElement.map(tmpAtomMap, tmpBondMap);
    }
    //
    /**
     * Checks whether all atoms that the given stereo element refers to are in the same connected component and that
     * none of its bonds was cleaved.
     *
     * @param aMolecule the molecule after the cleavage of the bonds
     * @param anElement stereo element
     * @param aCleavedBondSet the cleaved bonds
     * @param aComponents connected component labels of the atoms of the molecule
     * @return true if the stereo element lies within one fragment
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private static boolean isInOneComponent(IAtomContainer aMolecule, IStereoElement anElement, Set<IBond> aCleavedBondSet,
                                            int[] aComponents) {
        Set<IAtom> tmpAtoms = new HashSet<>();
        List<Object> tmpObjects = new ArrayList<>(anElement.getCarriers());
        tmpObjects.add(anElement.getFocus());
        for (Object tmpObject : tmpObjects) {
            if (tmpObject instanceof IAtom tmpAtom) {
                tmpAtoms.add(tmpAtom);
            } else if (tmpObject instanceof IBond tmpBond) {
                if (aCleavedBondSet.contains(tmpBond)) {
                    return false;
                }
                tmpAtoms.add(tmpBond.getBegin());
                tmpAtoms.add(tmpBond.getEnd());
            }
        }
        int tmpComponent = -1;
        for (IAtom tmpAtom : tmpAtoms) {
            int tmpIndex = aMolecule.indexOf(tmpAtom);
            if (tmpIndex < 0) {
                return false;
            }
            if (tmpComponent == -1) {
                tmpComponent = aComponents[tmpIndex];
            } else if (tmpComponent != aComponents[tmpIndex]) {
                return false;
            }
        }
        return true;
    }
    //
    /**
     * Removes tetrahedral and double bond stereo elements whose atoms were stereogenic in the molecule but are not
     * stereogenic in the given fragment anymore, i.e. whose stereogenicity was lost by the fragmentation, e.g. if two
     * ligands of a former stereocentre are link atoms now. Stereo elements of atoms that are not perceived as
     * stereogenic in the molecule either, e.g. pseudo-asymmetric centres like C3 of tropanes, are kept unchanged. The
     * bond types are not considered, since MORTAR identifies fragments by their SMILES codes without isotope
     * information, in which all link atoms are equal.
     *
     * @param aFragment fragment
     * @param aNonStereogenicAtomSetOfMolecule atoms that are not stereogenic in the molecule before the cleavage, or
     *                                         null if they could not be determined
     */
    @SuppressWarnings("rawtypes")
    private static void removeStereoElementsMadeNonStereogenic(IAtomContainer aFragment,
                                                                Set<IAtom> aNonStereogenicAtomSetOfMolecule) {
        if (!aFragment.stereoElements().iterator().hasNext()) {
            return;
        }
        Set<IAtom> tmpNonStereogenicAtoms = RECAPFragmenter.findNonStereogenicAtoms(aFragment);
        if (tmpNonStereogenicAtoms == null) {
            return;
        }
        List<IStereoElement> tmpRemainingElements = new ArrayList<>();
        for (IStereoElement tmpElement : aFragment.stereoElements()) {
            List<IAtom> tmpStereoAtoms = new ArrayList<>(2);
            if (tmpElement.getConfigClass() == IStereoElement.TH) {
                tmpStereoAtoms.add((IAtom) tmpElement.getFocus());
            } else if (tmpElement.getConfigClass() == IStereoElement.CT) {
                tmpStereoAtoms.add(((IBond) tmpElement.getFocus()).getBegin());
                tmpStereoAtoms.add(((IBond) tmpElement.getFocus()).getEnd());
            }
            boolean tmpIsStereogenicityLost = false;
            for (IAtom tmpAtom : tmpStereoAtoms) {
                if (tmpNonStereogenicAtoms.contains(tmpAtom) && (aNonStereogenicAtomSetOfMolecule == null
                        || !aNonStereogenicAtomSetOfMolecule.contains(tmpAtom))) {
                    tmpIsStereogenicityLost = true;
                    break;
                }
            }
            if (!tmpIsStereogenicityLost) {
                tmpRemainingElements.add(tmpElement);
            }
        }
        aFragment.setStereoElements(tmpRemainingElements);
    }
    //
    /**
     * Returns the atoms of the given molecule that are not stereogenic, see
     * {@link #findNonStereogenicAtomIndices(IAtomContainer)}.
     *
     * @param aMolecule molecule
     * @return set of non-stereogenic atoms, or null if the stereocentres cannot be determined
     */
    private static Set<IAtom> findNonStereogenicAtoms(IAtomContainer aMolecule) {
        BitSet tmpIndices = RECAPFragmenter.findNonStereogenicAtomIndices(aMolecule);
        if (tmpIndices == null) {
            return null;
        }
        Set<IAtom> tmpNonStereogenicAtoms = new HashSet<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpIndices.cardinality(), BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR),
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        for (int i = tmpIndices.nextSetBit(0); i >= 0; i = tmpIndices.nextSetBit(i + 1)) {
            tmpNonStereogenicAtoms.add(aMolecule.getAtom(i));
        }
        return tmpNonStereogenicAtoms;
    }
    //
    /**
     * Returns the indices of the atoms of the given molecule that are not stereogenic according to the CDK class
     * {@link Stereocenters}, which perceives stereocentres by constitutional symmetry. Mass numbers of pseudo atoms,
     * i.e. the bond types, are ignored.
     *
     * @param aMolecule molecule
     * @return indices of the non-stereogenic atoms, or null if the stereocentres cannot be determined
     */
    private static BitSet findNonStereogenicAtomIndices(IAtomContainer aMolecule) {
        //bond types are given as mass numbers and must be ignored for the symmetry perception
        List<IAtom> tmpLinkAtoms = new ArrayList<>();
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom instanceof IPseudoAtom && tmpAtom.getMassNumber() != null) {
                tmpLinkAtoms.add(tmpAtom);
            }
        }
        List<Integer> tmpMassNumbers = new ArrayList<>(tmpLinkAtoms.size());
        for (IAtom tmpLinkAtom : tmpLinkAtoms) {
            tmpMassNumbers.add(tmpLinkAtom.getMassNumber());
            tmpLinkAtom.setMassNumber(null);
        }
        try {
            Stereocenters tmpStereocenters = Stereocenters.of(aMolecule);
            BitSet tmpNonStereogenicAtomIndices = new BitSet(aMolecule.getAtomCount());
            for (int i = 0; i < aMolecule.getAtomCount(); i++) {
                if (!tmpStereocenters.isStereocenter(i)) {
                    tmpNonStereogenicAtomIndices.set(i);
                }
            }
            return tmpNonStereogenicAtomIndices;
        } catch (RuntimeException anException) {
            RECAPFragmenter.LOGGER.log(Level.FINE, "Stereocentres could not be determined", anException);
            return null;
        } finally {
            for (int i = 0; i < tmpLinkAtoms.size(); i++) {
                tmpLinkAtoms.get(i).setMassNumber(tmpMassNumbers.get(i));
            }
        }
    }
    //
    /**
     * Compiles the given rules: the atom environments are added to the given map of environment SMARTS to pattern
     * indices (if not contained yet), and the rules reference their environments by these indices.
     *
     * @param aRuleList rules to compile
     * @param anEnvironmentIndexMap map of environment SMARTS to pattern indices, extended by this method
     * @return the compiled rules in the order of the given list
     */
    private static CompiledRule[] compileRules(List<BondRule> aRuleList, Map<String, Integer> anEnvironmentIndexMap) {
        CompiledRule[] tmpCompiledRules = new CompiledRule[aRuleList.size()];
        for (int i = 0; i < aRuleList.size(); i++) {
            BondRule tmpRule = aRuleList.get(i);
            int tmpLeftIndex = anEnvironmentIndexMap.computeIfAbsent(tmpRule.leftEnvironment(),
                    aKey -> anEnvironmentIndexMap.size());
            int tmpRightIndex = anEnvironmentIndexMap.computeIfAbsent(tmpRule.rightEnvironment(),
                    aKey -> anEnvironmentIndexMap.size());
            tmpCompiledRules[i] = new CompiledRule(tmpRule, tmpLeftIndex, tmpRightIndex);
        }
        return tmpCompiledRules;
    }
    //
    /**
     * Creates the eleven RECAP bond rules in the order of their numbers, with or without the strict functional group
     * exclusions (see class documentation).
     *
     * @param isStrict true for the rules with the strict functional group exclusions
     * @return list of all rules
     */
    private static List<BondRule> createRules(boolean isStrict) {
        IBond.Order tmpSingle = IBond.Order.SINGLE;
        //carbonyl carbon atom of amides, not of ureas (bond type 4)
        String tmpAmideCarbon = "[C;$(C=O);!$(C(~[#7])~[#7])]";
        //nitrogen atom of amides, ureas and sulphonamides (also aromatic nitrogen atoms)
        String tmpAmideNitrogen = "[#7;+0]";
        String tmpUreaNitrogen = tmpAmideNitrogen;
        String tmpSulphonamideNitrogen = tmpAmideNitrogen;
        //carbonyl carbon atom of esters (also carbonates and carbamates)
        String tmpEsterCarbon = "[C;$(C=O)]";
        //ester oxygen atom, bonded to two carbon atoms (no carboxylic acids)
        String tmpEsterOxygen = "[O;D2;$(O([#6])[#6])]";
        //sp3 carbon atom of amines and ethers
        String tmpSp3Carbon = "[C;X4]";
        //amine nitrogen atom: neutral, three-connected, not acylated or bonded to an imine, thiocarbonyl or nitrile carbon
        //atom (amides, ureas, amidines, cyanamides), and not bonded to another heteroatom (e.g. sulfonamides, hydrazines)
        String tmpAmineNitrogen = "[N;+0;X3;!$(N[#6]=,#[#7,#8,#16]);!$(N[#7,#8,#15,#16])]";
        //carbonyl carbon atom of ureas
        String tmpUreaCarbon = "[C;$(C(=O)(~[#7])~[#7])]";
        //ether oxygen atom, bonded to two carbon atoms, not part of an ester, carbamate or cyanate
        String tmpEtherOxygen = "[O;D2;+0;$(O([#6])[#6]);!$(O[#6]=,#[#7,#8,#16])]";
        //olefinic carbon atom, not part of a cumulene
        String tmpOlefinCarbon = "[C;!$(C(=*)=*)]";
        //quaternary ammonium nitrogen atom with four carbon substituents
        String tmpQuaternaryNitrogen = "[N+;D4;!$([N+]~[!#6])]";
        //carbon atom bonded to a quaternary ammonium nitrogen atom
        String tmpQuaternaryNitrogenCarbon = "[#6]";
        //aliphatic (non-aromatic) carbon atom
        String tmpAliphaticCarbon = "[C]";
        //lactam nitrogen atom: ring nitrogen atom bonded to a ring carbonyl carbon atom
        String tmpLactamNitrogen = "[N;R;$(N@C=O)]";
        //sulfonyl sulfur atom
        String tmpSulfonylSulfur = "[S;X4;$(S(=O)=O)]";
        if (isStrict) {
            //no carbamates, carbamoyl halides, or 1,2-dicarbonyl compounds (e.g. alpha-keto amides, oxamides)
            tmpAmideCarbon = "[C;$(C=O);!$(C(~[#7])~[#7]);!$(C(=O)([#7])[!#6;!#1;!#0]);!$(C(=O)C=[#7,#8,#16])]";
            //no hydrazides, hydroxamic acids or N-acyl sulphonamides (non-aromatic single bond to a heteroatom), no
            //imides or acyl ureas (second acyl, imidoyl or thioacyl group), and no nitrogen atoms with a non-aromatic
            //double or triple bond
            tmpAmideNitrogen = "[#7;+0;!$([#7]-[!#6;!#1;!#0]);!$([#7]=,#*);!$([#7](C=[#7,#8,#16])C=[#7,#8,#16])]";
            tmpUreaNitrogen = tmpAmideNitrogen;
            //no acyl ureas, biurets or ureas at imide nitrogen atoms (nitrogen atom with a second acyl, imidoyl or
            //thioacyl group) and no sulphonyl ureas (nitrogen atom with a sulfonyl group), so that neither of their
            //urea bonds is cleaved
            tmpUreaCarbon = "[C;$(C(=O)(~[#7])~[#7]);!$(C(=O)~[#7]C=[#7,#8,#16]);!$(C(=O)~[#7][#16](=O)=O)]";
            //no N-acyl sulphonamides, sulphonyl hydrazides or disulphonimides
            tmpSulphonamideNitrogen = "[#7;+0;!$([#7]C=[#7,#8,#16]);!$([#7](S)[!#6;!#1;!#0]);!$([#7]=,#*)]";
            //no carbonates, carbamates, haloformates, or 1,2-dicarbonyl compounds
            tmpEsterCarbon = "[C;$(C=O);!$(C(=O)([#8])[!#6;!#1;!#0]);!$(C(=O)C=[#7,#8,#16])]";
            //no anhydrides
            tmpEsterOxygen = "[O;D2;$(O([#6])[#6]);!$(O(C=[#7,#8,#16])C=[#7,#8,#16])]";
            //no enols, enol ethers, enamines or other olefins with a nitrogen, oxygen or sulfur substituent, and no
            //double bonds conjugated with another carbon-carbon double or triple bond (polyenes)
            tmpOlefinCarbon = "[C;!$(C(=*)=*);!$(C-[#7,#8,#16]);!$(C-[#6]=,#[#6])]";
            //no acyl ammonium ions, i.e. no bond of a quaternary nitrogen atom that carries an acyl, imidoyl or cyano
            //group is cleaved
            tmpQuaternaryNitrogen = "[N+;D4;!$([N+]~[!#6]);!$([N+][#6]=,#[!#6])]";
            //no carbonyl, imidoyl or thiocarbonyl carbon atoms, e.g. of carbamates or amidines at the nitrogen atom
            tmpAliphaticCarbon = "[C;!$(C=[#7,#8,#16])]";
            //only lactams in the strict sense: the carbonyl carbon atom is bonded to a ring carbon atom (no cyclic ureas
            //or carbamates), and no cyclic imides, amidines or hydrazides
            tmpLactamNitrogen = "[N;R;$(N@[C;R](=O)@[#6]);!$(N-[!#6;!#1;!#0]);!$(N(C=[#7,#8,#16])C=[#7,#8,#16])]";
            //no sulphamides, sulphamates or sulphamoyl halides
            tmpSulfonylSulfur = "[S;X4;$(S(=O)=O);!$(S(=O)(=O)([#7])[!#6;!#1;!#0])]";
        }
        List<BondRule> tmpRules = new ArrayList<>(11);
        tmpRules.add(new BondRule(1, "amide", tmpAmideCarbon, tmpSingle, tmpAmideNitrogen));
        tmpRules.add(new BondRule(2, "ester", tmpEsterCarbon, tmpSingle, tmpEsterOxygen));
        tmpRules.add(new BondRule(3, "amine", tmpSp3Carbon, tmpSingle, tmpAmineNitrogen));
        tmpRules.add(new BondRule(4, "urea", tmpUreaNitrogen, tmpSingle, tmpUreaCarbon));
        tmpRules.add(new BondRule(5, "ether", tmpSp3Carbon, tmpSingle, tmpEtherOxygen));
        tmpRules.add(new BondRule(6, "olefin", tmpOlefinCarbon, IBond.Order.DOUBLE, tmpOlefinCarbon));
        tmpRules.add(new BondRule(7, "quaternary nitrogen", tmpQuaternaryNitrogen, tmpSingle, tmpQuaternaryNitrogenCarbon));
        tmpRules.add(new BondRule(8, "aromatic nitrogen - aliphatic carbon", "[n]", tmpSingle, tmpAliphaticCarbon));
        tmpRules.add(new BondRule(9, "lactam nitrogen - aliphatic carbon", tmpLactamNitrogen, tmpSingle, tmpAliphaticCarbon));
        tmpRules.add(new BondRule(10, "aromatic carbon - aromatic carbon", "[c]", tmpSingle, "[c]"));
        tmpRules.add(new BondRule(11, "sulphonamide", tmpSulphonamideNitrogen, tmpSingle, tmpSulfonylSulfur));
        return Collections.unmodifiableList(tmpRules);
    }
    //</editor-fold>
}
