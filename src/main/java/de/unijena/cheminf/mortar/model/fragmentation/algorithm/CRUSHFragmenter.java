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
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.model.util.SimpleIDisplayEnumConstantProperty;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;

import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
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
 * Implementation of the CRUSH fragmentation ("Cleavage Rules Using SMIRKS Heuristics") by López-López et al., see
 * the references below. CRUSH defines 33 bond cleavage rules as reaction SMIRKS (Table 1 of the publication), from
 * functional groups as in RECAP and BRICS (amides, esters, ureas, carbamates, ethers, amines, sulphonamides) to
 * aliphatic carbon-carbon and carbon-heteroatom bonds, aryl-alkyl, aryl-heteroatom and biaryl bonds, and sulfur-containing
 * groups. A molecule is fragmented recursively, until no rule can be applied anymore, under three constraints: every
 * fragment has at least three heavy atoms, atoms in fused ring systems (and their neighbours) are protected, and a rule
 * with a match on a protected atom is not applied at all (conservative mode). Among all fragmentation routes, the one
 * with the maximum number of fragments is selected.
 *
 * <p><b>Rules.</b> The rules are numbered in the order of Table 1 of the publication (and of the rule file of the
 * authors); the reactant side of each SMIRKS is used verbatim as SMARTS, and the bonds between mapped atoms that are
 * missing on the product side are cleaved (one bond per match; two bonds for rule 19, Carbamate_2). The attachment
 * points of the fragments are marked with link atoms that carry the number of the rule. Unspecified bonds of the SMARTS
 * match single and aromatic bonds and '!@' matches any acyclic bond, so rules 25 and 26 also cleave acyclic double and
 * triple bonds (e.g. the C=N bond of oximes), as the reaction SMIRKS of the authors do.
 *
 * <p><b>Algorithm.</b>
 * <ol>
 *     <li>Explicit hydrogen atoms are suppressed, and ring membership and aromaticity are perceived with the Daylight
 *     model, on which the SMARTS semantics of the rules are based. Every connected component of the molecule (e.g. a
 *     counter ion) is a start fragment of its own.</li>
 *     <li>A fragment is a connected set of atoms of the molecule. Its cleavage steps are determined on a copy that is
 *     saturated with hydrogen atoms at the bonds cleaved before (as the products of the authors' reactions are), with
 *     ring membership and aromaticity perceived again:
 *     <ul>
 *         <li>protected atoms are the atoms in at least two rings (fused ring systems, spiro and bridgehead atoms, ring
 *         set: relevant cycles) and, by default, their neighbours (setting);</li>
 *         <li>for every rule, all matches are determined; in conservative mode (setting, default), the rule yields no
 *         step if any match contains a protected atom, otherwise only the matches with protected atoms are dropped;</li>
 *         <li>a match yields a step if none of its cleaved bonds is a ring bond and all parts after the cleavage have at
 *         least the minimum number of heavy atoms (setting, default {@value #MINIMUM_FRAGMENT_SIZE_DEFAULT});</li>
 *         <li>steps are unique per set of cleaved bonds; if several rules cleave the same bonds, the rule with the lowest
 *         number is assigned.</li>
 *     </ul></li>
 *     <li>The route is a complete decomposition tree of each start fragment: the children of a fragment are the parts
 *     of one of its steps, the leaves are fragments without steps. The route with the maximum number of leaves is
 *     determined by an exact search that is memoised per atom set. Ties between steps are broken by the lowest rule
 *     number and then by the lexicographically smallest sorted list of the unique SMILES codes of the parts (with link
 *     atoms, rule numbers, and stereochemistry), so that the result does not depend on the order of the atoms. At most
 *     {@value #MAXIMUM_SEARCH_FRAGMENT_COUNT} fragments are evaluated per molecule; if this number is exceeded (e.g. for
 *     large glycosides or peptides), the remaining fragments are decomposed greedily (first step in the order above),
 *     see {@link CRUSHRouteNode#isExhaustive()}.</li>
 *     <li>The final fragments (leaves of the route) are returned, so every heavy atom of the molecule is contained in
 *     exactly one fragment; optionally, all fragments of the route are returned (setting). Each cleaved bond is replaced
 *     by two bonds of the same order to new link atoms, pseudo atoms carrying the rule number as mass number (resulting
 *     in SMILES like [26*]OCC, if the mass number is encoded) and as property {@link #RULE_NUMBER_PROPERTY_KEY}.
 *     Alternatively, the open valences can be saturated with hydrogen atoms (setting), as in the output of the authors'
 *     script.</li>
 *     <li>Stereo elements are transferred to the fragments; ligands that were cleaved off are replaced by the link atom
 *     (or hydrogen atom) at their position, and stereo elements that lost their stereogenicity are removed. As in the
 *     RECAP fragmenter, stereogenicity is perceived from the constitution (CDK class {@link Stereocenters}), so the
 *     configuration of a potential (para) stereocentre is kept even if it is not stereogenic because of the
 *     configuration of other atoms (e.g. C3 of the 1,5-anhydroarabinitol unit that rule 7 forms from xylopyranoses,
 *     whose configuration the RDKit removes). A configuration-aware perception would be needed to remove it (in the meso
 *     isomers, this atom is a genuine pseudo-asymmetric centre); like for any molecule with such an atom in MORTAR, the
 *     unique SMILES code of such a fragment then depends on the configuration given for the atom.</li>
 * </ol>
 *
 * <p><b>Differences to the script of the authors.</b> The script applies the reaction SMIRKS with the RDKit and
 * follows, in each recursion step, only one product further; the products of all steps of the best path are its
 * output. The single steps are reproduced exactly, except for the following deliberate differences:
 * <ul>
 *     <li>Ring bonds are never cleaved. The script applies the rules without the '!@' restriction (rules 3, 7, 28,
 *     29 to 33) also to ring bonds, which yields two overlapping products with duplicated atoms (e.g. open-chain
 *     amines from N-alkyl piperidines, and from tetrahydrofurans and sugars); this contradicts the complementary
 *     fragments described in the publication.</li>
 *     <li>Bonds of neutral aromatic nitrogen atoms are cleaved (rules 10, 22, 23), the nitrogen atom obtains a hydrogen
 *     atom. In the script, the product SMILES lacks this hydrogen atom and cannot be parsed, so the cleavage is dropped
 *     silently and only cationic aromatic nitrogen atoms are cleaved.</li>
 *     <li>The ether rule (6) is applied without changing bond orders. Its product template turns the bonds next to the
 *     cleaved bond into single bonds, which makes the products of some heteroaromatic ethers (e.g. thienyl ethers)
 *     unparsable in the script.</li>
 *     <li>Instead of one chain of steps whose selection depends on the order of the rules (shuffled with a fixed seed),
 *     the order of the matches of the RDKit and the pruning of the search, the complete decomposition with the maximum
 *     number of fragments is determined, the criterion stated in the publication.</li>
 *     <li>The final fragments are returned by default; the intermediate fragments of the route, which the script
 *     returns as well, are only returned with the respective setting.</li>
 *     <li>The attachment points are marked with link atoms by default, as described in the publication; the script
 *     saturates them with hydrogen atoms (available as setting).</li>
 *     <li>Counter ions and other disconnected parts are kept as fragments; in the script, they disappear after the first
 *     step.</li>
 *     <li>With hydrogen saturation, each atom of a cleaved bond obtains one hydrogen atom per bond order. The reactions
 *     of the script recompute the hydrogen counts from the default valences of the RDKit, which differs for atoms
 *     with radical electrons in the input (they become saturated) and for hypervalent atoms with hydrogen atoms
 *     (e.g. a phosphinate O=[PH](O)C yields O=PO instead of O=[PH2]O).</li>
 * </ul>
 * Because of the overlapping ring-opened products and the intermediate fragments in the output of the script, the
 * fragment statistics of the publication cannot be reproduced exactly with this class.
 *
 * <p><b>Note on rule numbers in MORTAR.</b> MORTAR identifies fragments by their unique SMILES code without isotope
 * information, so link atoms are written as '*' and fragments that differ only in their rule numbers are counted as the
 * same fragment.
 *
 * <p>References:
 * <ul>
 *     <li>E. López-López, J. L. Medina-Franco, F. Miljković, "CRUSH: Cleavage Rules Using SMIRKS Heuristics",
 *     J. Cheminform. 2026, 18, 78, <a href="https://doi.org/10.1186/s13321-026-01239-w">
 *     doi:10.1186/s13321-026-01239-w</a> (Table 1: rules; Methods: recursive fragmentation and constraints)</li>
 *     <li>Rule file and script of the authors: <a href="https://github.com/EdgL2/CRUSH">github.com/EdgL2/CRUSH</a>
 *     (MIT licence, commit 08e0040c84), also <a href="https://doi.org/10.5281/zenodo.18665641">
 *     doi:10.5281/zenodo.18665641</a></li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class CRUSHFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum SaturationOption">
    /**
     * Enumeration of the ways to saturate the open valences that result from the cleavage of a bond.
     */
    public enum SaturationOption implements IDisplayEnum {
        /**
         * Every open valence is saturated with a link atom (pseudo atom) that carries the number of the CRUSH rule.
         */
        LINK_ATOMS(Message.get("CRUSHFragmenter.SaturationOption.linkAtoms.displayName"),
                Message.get("CRUSHFragmenter.SaturationOption.linkAtoms.tooltip")),
        /**
         * Every open valence is saturated with (implicit) hydrogen atoms, i.e. both atoms of a cleaved double bond
         * receive two hydrogen atoms, as in the output of the script of the authors.
         */
        HYDROGEN_SATURATION(Message.get("CRUSHFragmenter.SaturationOption.hydrogenSaturation.displayName"),
                Message.get("CRUSHFragmenter.SaturationOption.hydrogenSaturation.tooltip"));
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
    //<editor-fold desc="Enum RingProtectionOption">
    /**
     * Enumeration of the atoms that are protected from cleavage.
     */
    public enum RingProtectionOption implements IDisplayEnum {
        /**
         * Atoms in at least two rings (fused ring systems, spiro and bridgehead atoms) and their neighbours are
         * protected, as in the publication.
         */
        FUSED_RING_ATOMS_AND_NEIGHBOURS(Message.get("CRUSHFragmenter.RingProtectionOption.fusedRingAtomsAndNeighbours.displayName"),
                Message.get("CRUSHFragmenter.RingProtectionOption.fusedRingAtomsAndNeighbours.tooltip")),
        /**
         * Only atoms in at least two rings are protected.
         */
        FUSED_RING_ATOMS(Message.get("CRUSHFragmenter.RingProtectionOption.fusedRingAtoms.displayName"),
                Message.get("CRUSHFragmenter.RingProtectionOption.fusedRingAtoms.tooltip")),
        /**
         * No atom is protected.
         */
        NONE(Message.get("CRUSHFragmenter.RingProtectionOption.none.displayName"),
                Message.get("CRUSHFragmenter.RingProtectionOption.none.tooltip"));
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
        private RingProtectionOption(String aDisplayName, String aTooltip) {
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
         * The final fragments of the route, i.e. its leaves; every heavy atom of the molecule is contained in exactly one
         * fragment.
         */
        FINAL_FRAGMENTS(Message.get("CRUSHFragmenter.FragmentOutputOption.finalFragments.displayName"),
                Message.get("CRUSHFragmenter.FragmentOutputOption.finalFragments.tooltip")),
        /**
         * All fragments of the route (see {@link CRUSHFragmenter#buildRoute(IAtomContainer)}) except the molecule
         * itself, i.e. the final fragments and all intermediate fragments, as in the output of the script of the
         * authors; these fragments overlap.
         */
        ALL_ROUTE_FRAGMENTS(Message.get("CRUSHFragmenter.FragmentOutputOption.allRouteFragments.displayName"),
                Message.get("CRUSHFragmenter.FragmentOutputOption.allRouteFragments.tooltip"));
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
    //<editor-fold desc="Public static nested class CRUSHRouteNode">
    /**
     * A node of the fragmentation route of a molecule, see {@link CRUSHFragmenter#buildRoute(IAtomContainer)}: a
     * fragment with the fragments that result from it by one cleavage step (children). The fragment is created when it
     * is requested for the first time.
     */
    public static final class CRUSHRouteNode {
        /**
         * The route this node belongs to.
         */
        private final Route route;
        /**
         * Indices of the atoms of the fragment in the prepared molecule.
         */
        private final BitSet atomSet;
        /**
         * The parent node, null for the root.
         */
        private final CRUSHRouteNode parent;
        /**
         * The nodes of the parts that result from the cleavage step of this fragment.
         */
        private final List<CRUSHRouteNode> children;
        /**
         * Number of the rule of the cleavage step of this fragment, 0 for leaves and for the root of a molecule with
         * several connected components.
         */
        private int ruleNumber;
        /**
         * The fragment, created when requested for the first time.
         */
        private IAtomContainer fragment;
        /**
         * Constructor.
         *
         * @param aRoute the route
         * @param anAtomSet indices of the atoms of the fragment in the prepared molecule
         * @param aParent the parent node, null for the root
         */
        private CRUSHRouteNode(Route aRoute, BitSet anAtomSet, CRUSHRouteNode aParent) {
            this.route = aRoute;
            this.atomSet = anAtomSet;
            this.parent = aParent;
            this.children = new ArrayList<>(2);
            this.ruleNumber = 0;
            this.fragment = null;
        }
        //
        /**
         * Returns the fragment of this node; for the root, it is the prepared molecule. The fragments of the other
         * nodes are created from the prepared molecule when they are requested for the first time, with link atoms or
         * hydrogen atoms at the attachment points (according to the saturation setting at the time the route was built)
         * and with the stereo elements transferred.
         *
         * @return the fragment
         */
        public IAtomContainer getFragment() {
            if (this.fragment == null) {
                if (this.parent == null) {
                    this.fragment = this.route.molecule();
                } else {
                    try {
                        this.fragment = CRUSHFragmenter.createFragment(this.route.molecule(), this.atomSet,
                                this.route.bondLabels(), this.route.isHydrogenSaturation());
                    } catch (CloneNotSupportedException anException) {
                        throw new IllegalStateException("The fragment of a route node could not be created.", anException);
                    }
                }
            }
            return this.fragment;
        }
        //
        /**
         * Returns the parent node, i.e. the node whose cleavage step yields this fragment.
         *
         * @return the parent node, null for the root
         */
        public CRUSHRouteNode getParent() {
            return this.parent;
        }
        //
        /**
         * Returns the child nodes, i.e. the nodes of the parts that result from the cleavage step of this fragment (two
         * parts, three for rule 19). For the root of a molecule with several connected components, the children are the
         * components.
         *
         * @return unmodifiable list of the child nodes, empty for leaves
         */
        public List<CRUSHRouteNode> getChildren() {
            return Collections.unmodifiableList(this.children);
        }
        //
        /**
         * Returns the number of the rule of the cleavage step of this fragment.
         *
         * @return rule number (1 to 33), 0 for leaves and for the root of a molecule with several connected components
         */
        public int getRuleNumber() {
            return this.ruleNumber;
        }
        //
        /**
         * Returns whether this node is the root of the route, i.e. represents the molecule.
         *
         * @return true for the root
         */
        public boolean isRoot() {
            return this.parent == null;
        }
        //
        /**
         * Returns whether this node is a leaf of the route, i.e. a final fragment.
         *
         * @return true for leaves
         */
        public boolean isLeaf() {
            return this.children.isEmpty();
        }
        //
        /**
         * Returns whether the route was determined by an exhaustive search, i.e. it has the maximum number of final
         * fragments. If the search exceeded {@link #MAXIMUM_SEARCH_FRAGMENT_COUNT} evaluated fragments, some fragments
         * were decomposed greedily and the number of final fragments may be lower than the maximum.
         *
         * @return true if the search was exhaustive
         */
        public boolean isExhaustive() {
            return this.route.isExhaustive();
        }
        //
        /**
         * Returns all nodes below this node in breadth-first order.
         *
         * @return list of all descendant nodes
         */
        public List<CRUSHRouteNode> getAllDescendants() {
            List<CRUSHRouteNode> tmpDescendants = new ArrayList<>();
            Deque<CRUSHRouteNode> tmpQueue = new ArrayDeque<>(this.children);
            while (!tmpQueue.isEmpty()) {
                CRUSHRouteNode tmpNode = tmpQueue.poll();
                tmpDescendants.add(tmpNode);
                tmpQueue.addAll(tmpNode.children);
            }
            return tmpDescendants;
        }
        //
        /**
         * Returns all leaves below this node in breadth-first order, i.e. the final fragments of its decomposition.
         *
         * @return list of the descendant leaves, empty if this node is a leaf itself
         */
        public List<CRUSHRouteNode> getLeaves() {
            List<CRUSHRouteNode> tmpLeaves = new ArrayList<>();
            for (CRUSHRouteNode tmpNode : this.getAllDescendants()) {
                if (tmpNode.isLeaf()) {
                    tmpLeaves.add(tmpNode);
                }
            }
            return tmpLeaves;
        }
        //
        /**
         * Returns the indices of the atoms of the fragment in the prepared molecule.
         *
         * @return copy of the atom set
         */
        BitSet getAtomSet() {
            return (BitSet) this.atomSet.clone();
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm.
     */
    public static final String ALGORITHM_NAME = "CRUSH";
    /**
     * Default option for the saturation of open valences.
     */
    public static final SaturationOption SATURATION_OPTION_DEFAULT = SaturationOption.LINK_ATOMS;
    /**
     * Default minimum number of heavy atoms of a fragment, as in the publication.
     */
    public static final int MINIMUM_FRAGMENT_SIZE_DEFAULT = 3;
    /**
     * Default option for the protected atoms: atoms in fused ring systems and their neighbours, as used by the authors.
     */
    public static final RingProtectionOption RING_PROTECTION_OPTION_DEFAULT = RingProtectionOption.FUSED_RING_ATOMS_AND_NEIGHBOURS;
    /**
     * Default setting for the conservative mode: a rule with a match on a protected atom is not applied at all, as in
     * the publication.
     */
    public static final boolean CONSERVATIVE_MODE_OPTION_DEFAULT = true;
    /**
     * Default option for the returned fragments: the final fragments of the route.
     */
    public static final FragmentOutputOption FRAGMENT_OUTPUT_OPTION_DEFAULT = FragmentOutputOption.FINAL_FRAGMENTS;
    /**
     * Default setting for returning molecules without cleavage steps as one unchanged fragment.
     */
    public static final boolean RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT = true;
    /**
     * Maximum number of fragments whose cleavage steps are evaluated in the exact search for the route with the
     * maximum number of fragments, per molecule; the remaining fragments are decomposed greedily. Determined on the
     * BIOFACQUIM data set: 99 % of the molecules need at most this number of evaluations.
     */
    public static final int MAXIMUM_SEARCH_FRAGMENT_COUNT = 3000;
    /**
     * Key of the atom property of link atoms that holds the number of the CRUSH rule of their cleaved bond (Integer, 1
     * to 33).
     */
    public static final String RULE_NUMBER_PROPERTY_KEY = "CRUSHFragmenter.RuleNumber";
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(CRUSHFragmenter.class.getName());
    /**
     * All CRUSH rules in the order of their numbers.
     */
    private static final List<CleavageRule> RULES = CRUSHFragmenter.createRules();
    //</editor-fold>
    //
    //<editor-fold desc="Package-private records">
    /**
     * A CRUSH rule: the reactant side of the reaction SMIRKS of the authors (used as SMARTS) and the bonds it cleaves.
     *
     * @param number number of the rule, 1 to 33, in the order of Table 1 of the publication
     * @param name name of the rule as given by the authors
     * @param smarts reactant side of the reaction SMIRKS
     * @param cleavedBonds the cleaved bonds as pairs of atom map numbers of the SMARTS
     */
    record CleavageRule(int number, String name, String smarts, int[][] cleavedBonds) {
    }
    //
    /**
     * A cleavage step of a fragment.
     *
     * @param ruleNumber number of the rule
     * @param cleavedBonds indices of the cleaved bonds in the prepared molecule, ascending
     * @param parts atom sets of the parts after the cleavage, ordered by their lowest atom index
     */
    record CleavageStep(int ruleNumber, int[] cleavedBonds, BitSet[] parts) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private records and classes">
    /**
     * A rule with its compiled SMARTS pattern.
     *
     * @param rule the rule
     * @param pattern the compiled pattern of the SMARTS
     * @param cleavedQueryAtomPairs the cleaved bonds as pairs of query atom indices
     */
    private record CompiledRule(CleavageRule rule, DfPattern pattern, int[][] cleavedQueryAtomPairs) {
    }
    //
    /**
     * A bond to cleave with the rule number of its link atoms.
     *
     * @param bond the bond
     * @param ruleNumber number of the rule, 0 for link atoms without a rule number
     */
    private record CRUSHBond(IBond bond, int ruleNumber) {
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
    //
    /**
     * The data shared by the nodes of a route.
     *
     * @param molecule the prepared molecule, it must not be altered
     * @param bondLabels rule numbers of the cleaved bonds by bond index, 0 for bonds that are not cleaved
     * @param isHydrogenSaturation true if the fragments are saturated with hydrogen atoms instead of link atoms
     * @param isExhaustive true if the search was exhaustive
     */
    private record Route(IAtomContainer molecule, int[] bondLabels, boolean isHydrogenSaturation, boolean isExhaustive) {
    }
    //
    /**
     * The arguments that determine the tie-breaking key of a part, see
     * {@link #createTieBreakingKey(IAtomContainer, BitSet, int[], boolean)}.
     *
     * @param atomSet atom set of the part
     * @param boundaryBondLabels rule numbers of the bonds that connect the part to the rest of the molecule, in the order
     *                           of their bond indices
     * @param isStereoEncoded whether stereochemistry is encoded
     */
    private record TieBreakingKeyArguments(BitSet atomSet, List<Integer> boundaryBondLabels, boolean isStereoEncoded) {
    }
    //
    /**
     * The search state of a fragment: its unique cleavage steps, its number of final fragments, and the step chosen
     * greedily, if the search budget was exhausted.
     */
    private static final class SearchEntry {
        /**
         * Unique cleavage steps of the fragment.
         */
        private final List<CleavageStep> steps;
        /**
         * Number of final fragments of the best decomposition found, 0 if not determined yet.
         */
        private int fragmentCount;
        /**
         * The step chosen greedily, null if the fragment count was determined exactly or the fragment has no step.
         */
        private CleavageStep greedyStep;
        /**
         * Constructor.
         *
         * @param aStepList unique cleavage steps of the fragment
         */
        private SearchEntry(List<CleavageStep> aStepList) {
            this.steps = aStepList;
            this.fragmentCount = 0;
            this.greedyStep = null;
        }
    }
    //
    /**
     * The search for the route with the maximum number of final fragments of a prepared molecule.
     */
    private final class RouteSearch {
        /**
         * The prepared molecule.
         */
        private final IAtomContainer molecule;
        /**
         * Adjacency list of the molecule.
         */
        private final int[][] adjacencyList;
        /**
         * Map of the edges of the adjacency list to the bonds.
         */
        private final GraphUtil.EdgeToBondMap edgeToBondMap;
        /**
         * Atoms of the molecule in at least two rings.
         */
        private final BitSet fusedRingAtoms;
        /**
         * Search states of all evaluated fragments by atom set.
         */
        private final Map<BitSet, SearchEntry> entries;
        /**
         * Rule numbers of the cleaved bonds of the route by bond index, filled when the route is assembled.
         */
        private final int[] bondLabels;
        /**
         * Tie-breaking keys of the parts of cleavage steps that were already created, since the same part occurs in
         * several steps.
         */
        private final Map<TieBreakingKeyArguments, String> tieBreakingKeys;
        /**
         * Number of evaluated fragments.
         */
        private int evaluationCount;
        /**
         * Whether all fragments were evaluated exactly.
         */
        private boolean isExhaustive;
        /**
         * Constructor.
         *
         * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
         */
        private RouteSearch(IAtomContainer aPreparedMolecule) {
            this.molecule = aPreparedMolecule;
            this.edgeToBondMap = GraphUtil.EdgeToBondMap.withSpaceFor(aPreparedMolecule);
            this.adjacencyList = GraphUtil.toAdjList(aPreparedMolecule, this.edgeToBondMap);
            this.fusedRingAtoms = CRUSHFragmenter.findFusedRingAtoms(aPreparedMolecule);
            this.entries = new HashMap<>();
            this.bondLabels = new int[aPreparedMolecule.getBondCount()];
            this.tieBreakingKeys = new HashMap<>();
            this.evaluationCount = 0;
            this.isExhaustive = true;
        }
        //
        /**
         * Returns the search state of the given fragment; its cleavage steps are determined if it is new.
         *
         * @param anAtomSet atom set of the fragment
         * @return search state
         */
        private SearchEntry entryOf(BitSet anAtomSet) {
            SearchEntry tmpEntry = this.entries.get(anAtomSet);
            if (tmpEntry == null) {
                this.evaluationCount++;
                tmpEntry = new SearchEntry(CRUSHFragmenter.this.findUniqueSteps(this.molecule, anAtomSet,
                        this.adjacencyList, this.edgeToBondMap, this.fusedRingAtoms));
                this.entries.put(anAtomSet, tmpEntry);
            }
            return tmpEntry;
        }
        //
        /**
         * Returns the maximum number of final fragments of the given fragment; exact as long as the search budget is not
         * exhausted, afterwards by greedy decomposition.
         *
         * @param anAtomSet atom set of the fragment
         * @return number of final fragments
         * @throws CloneNotSupportedException if a fragment for the tie-breaking cannot be created
         */
        private int countFragments(BitSet anAtomSet) throws CloneNotSupportedException {
            SearchEntry tmpEntry = this.entries.get(anAtomSet);
            if (tmpEntry != null && tmpEntry.fragmentCount > 0) {
                return tmpEntry.fragmentCount;
            }
            boolean tmpIsBudgetExhausted = tmpEntry == null
                    && this.evaluationCount >= CRUSHFragmenter.MAXIMUM_SEARCH_FRAGMENT_COUNT;
            tmpEntry = this.entryOf(anAtomSet);
            if (tmpEntry.steps.isEmpty()) {
                tmpEntry.fragmentCount = 1;
                return 1;
            }
            if (tmpIsBudgetExhausted) {
                this.isExhaustive = false;
                CleavageStep tmpStep = this.selectStep(tmpEntry.steps, null);
                int tmpCount = 0;
                for (BitSet tmpPart : tmpStep.parts()) {
                    tmpCount += this.countFragments(tmpPart);
                }
                tmpEntry.greedyStep = tmpStep;
                tmpEntry.fragmentCount = tmpCount;
                return tmpCount;
            }
            int tmpMaximum = 0;
            for (CleavageStep tmpStep : tmpEntry.steps) {
                int tmpCount = 0;
                for (BitSet tmpPart : tmpStep.parts()) {
                    tmpCount += this.countFragments(tmpPart);
                }
                tmpMaximum = Math.max(tmpMaximum, tmpCount);
            }
            tmpEntry.fragmentCount = tmpMaximum;
            return tmpMaximum;
        }
        //
        /**
         * Creates the route node of the given fragment and, recursively, the nodes of its decomposition: for fragments
         * that were evaluated exactly, the step with the maximum number of final fragments is selected (ties are broken
         * as described in the class documentation); for fragments that were decomposed greedily, the greedy step. The
         * rule numbers of the cleaved bonds are set in the bond labels.
         *
         * @param aNode the node of the fragment, its children are added
         * @throws CloneNotSupportedException if a fragment for the tie-breaking cannot be created
         */
        private void assemble(CRUSHRouteNode aNode) throws CloneNotSupportedException {
            SearchEntry tmpEntry = this.entries.get(aNode.atomSet);
            if (tmpEntry == null || tmpEntry.steps.isEmpty()) {
                return;
            }
            CleavageStep tmpStep = tmpEntry.greedyStep;
            if (tmpStep == null) {
                List<CleavageStep> tmpCandidates = new ArrayList<>(tmpEntry.steps.size());
                for (CleavageStep tmpCandidate : tmpEntry.steps) {
                    int tmpCount = 0;
                    for (BitSet tmpPart : tmpCandidate.parts()) {
                        tmpCount += this.entries.get(tmpPart).fragmentCount;
                    }
                    if (tmpCount == tmpEntry.fragmentCount) {
                        tmpCandidates.add(tmpCandidate);
                    }
                }
                tmpStep = this.selectStep(tmpCandidates, this.bondLabels);
            }
            for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                this.bondLabels[tmpBondIndex] = tmpStep.ruleNumber();
            }
            aNode.ruleNumber = tmpStep.ruleNumber();
            for (BitSet tmpPart : tmpStep.parts()) {
                CRUSHRouteNode tmpChild = new CRUSHRouteNode(aNode.route, tmpPart, aNode);
                aNode.children.add(tmpChild);
            }
            for (CRUSHRouteNode tmpChild : aNode.children) {
                this.assemble(tmpChild);
            }
        }
        //
        /**
         * Selects the first of the given steps in the tie-breaking order: lowest rule number, then the lexicographically
         * smallest sorted list of the unique SMILES codes of the parts with link atoms. If bond labels are given, the
         * link atoms carry the rule numbers of the bonds cleaved before and of the step, and stereochemistry is encoded;
         * otherwise (greedy decomposition), the link atoms of the bonds cleaved before carry no rule number.
         *
         * @param aStepList steps, sorted by rule number
         * @param aBondLabelArray rule numbers of the bonds cleaved before by bond index, or null
         * @return the selected step
         * @throws CloneNotSupportedException if a fragment cannot be created
         */
        private CleavageStep selectStep(List<CleavageStep> aStepList, int[] aBondLabelArray)
                throws CloneNotSupportedException {
            int tmpLowestRuleNumber = aStepList.getFirst().ruleNumber();
            List<CleavageStep> tmpTies = new ArrayList<>(2);
            for (CleavageStep tmpStep : aStepList) {
                if (tmpStep.ruleNumber() == tmpLowestRuleNumber) {
                    tmpTies.add(tmpStep);
                }
            }
            if (tmpTies.size() == 1) {
                return tmpTies.getFirst();
            }
            CleavageStep tmpSelectedStep = null;
            List<String> tmpSelectedKey = null;
            for (CleavageStep tmpStep : tmpTies) {
                int[] tmpLabels = aBondLabelArray == null ? new int[this.molecule.getBondCount()] : aBondLabelArray;
                for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                    tmpLabels[tmpBondIndex] = tmpStep.ruleNumber();
                }
                List<String> tmpKey = new ArrayList<>(tmpStep.parts().length);
                for (BitSet tmpPart : tmpStep.parts()) {
                    tmpKey.add(this.tieBreakingKeyOf(tmpPart, tmpLabels, aBondLabelArray != null));
                }
                for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                    tmpLabels[tmpBondIndex] = 0;
                }
                Collections.sort(tmpKey);
                if (tmpSelectedKey == null || CRUSHFragmenter.compare(tmpKey, tmpSelectedKey) < 0) {
                    tmpSelectedKey = tmpKey;
                    tmpSelectedStep = tmpStep;
                }
            }
            return tmpSelectedStep;
        }
        //
        /**
         * Returns the tie-breaking key of the given part, see
         * {@link #createTieBreakingKey(IAtomContainer, BitSet, int[], boolean)}; keys are created once per part and
         * rule numbers of its boundary bonds. Without stereochemistry, the key is created from a fragment that consists
         * of the atoms of the part only, see {@link #createTieBreakingKeyWithoutStereo(IAtomContainer, BitSet, int[])}.
         *
         * @param anAtomSet atom set of the part
         * @param aBondLabelArray rule numbers by bond index for the link atoms
         * @param isStereoEncoded whether stereochemistry is encoded
         * @return key
         * @throws CloneNotSupportedException if an atom or the molecule cannot be cloned
         */
        private String tieBreakingKeyOf(BitSet anAtomSet, int[] aBondLabelArray, boolean isStereoEncoded)
                throws CloneNotSupportedException {
            List<Integer> tmpBoundaryBondIndices = new ArrayList<>();
            for (int i = anAtomSet.nextSetBit(0); i >= 0; i = anAtomSet.nextSetBit(i + 1)) {
                for (int tmpNeighbour : this.adjacencyList[i]) {
                    if (!anAtomSet.get(tmpNeighbour)) {
                        tmpBoundaryBondIndices.add(this.edgeToBondMap.get(i, tmpNeighbour).getIndex());
                    }
                }
            }
            Collections.sort(tmpBoundaryBondIndices);
            List<Integer> tmpBoundaryBondLabels = new ArrayList<>(tmpBoundaryBondIndices.size());
            for (int tmpBondIndex : tmpBoundaryBondIndices) {
                tmpBoundaryBondLabels.add(aBondLabelArray[tmpBondIndex]);
            }
            TieBreakingKeyArguments tmpArguments = new TieBreakingKeyArguments(anAtomSet, tmpBoundaryBondLabels,
                    isStereoEncoded);
            String tmpKey = this.tieBreakingKeys.get(tmpArguments);
            if (tmpKey == null) {
                tmpKey = isStereoEncoded
                        ? CRUSHFragmenter.createTieBreakingKey(this.molecule, anAtomSet, aBondLabelArray, true)
                        : CRUSHFragmenter.createTieBreakingKeyWithoutStereo(this.molecule, anAtomSet, aBondLabelArray);
                this.tieBreakingKeys.put(tmpArguments, tmpKey);
            }
            return tmpKey;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * The rules with their compiled SMARTS patterns, in the order of their numbers.
     */
    private final CompiledRule[] compiledRules;
    //the settings are documented at their property methods (a Javadoc comment here would duplicate it for javadoc)
    //setting for the saturation of open valences
    private final SimpleIDisplayEnumConstantProperty saturationSetting;
    //setting for the minimum number of heavy atoms of a fragment
    private final SimpleIntegerProperty minimumFragmentSizeSetting;
    //setting for the protected atoms
    private final SimpleIDisplayEnumConstantProperty ringProtectionSetting;
    //setting for the conservative mode
    private final SimpleBooleanProperty conservativeModeSetting;
    //setting for the returned fragments (final fragments or all fragments of the route)
    private final SimpleIDisplayEnumConstantProperty fragmentOutputSetting;
    //setting for returning molecules without cleavage steps as one unchanged fragment
    private final SimpleBooleanProperty returnUnfragmentedMoleculesSetting;
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
    public CRUSHFragmenter() {
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        this.compiledRules = new CompiledRule[CRUSHFragmenter.RULES.size()];
        for (int i = 0; i < this.compiledRules.length; i++) {
            CleavageRule tmpRule = CRUSHFragmenter.RULES.get(i);
            QueryAtomContainer tmpQuery = new QueryAtomContainer(tmpBuilder);
            if (!Smarts.parse(tmpQuery, tmpRule.smarts())) {
                throw new IllegalStateException(String.format("CRUSH rule SMARTS %s cannot be parsed: %s",
                        tmpRule.smarts(), Smarts.getLastErrorMesg()));
            }
            Map<Integer, Integer> tmpMapNumberToQueryIndex = new HashMap<>();
            for (int j = 0; j < tmpQuery.getAtomCount(); j++) {
                tmpMapNumberToQueryIndex.put(tmpQuery.getAtom(j).getMapIdx(), j);
            }
            int[][] tmpCleavedQueryAtomPairs = new int[tmpRule.cleavedBonds().length][];
            for (int j = 0; j < tmpCleavedQueryAtomPairs.length; j++) {
                tmpCleavedQueryAtomPairs[j] = new int[] {
                        tmpMapNumberToQueryIndex.get(tmpRule.cleavedBonds()[j][0]),
                        tmpMapNumberToQueryIndex.get(tmpRule.cleavedBonds()[j][1])};
            }
            this.compiledRules[i] = new CompiledRule(tmpRule, DfPattern.findSubstructure(tmpQuery), tmpCleavedQueryAtomPairs);
        }
        int tmpNumberOfSettings = 6;
        int tmpInitialCapacityForSettingNameMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.saturationSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Fragment saturation setting",
                CRUSHFragmenter.SATURATION_OPTION_DEFAULT,
                CRUSHFragmenter.SaturationOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    CRUSHFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
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
                Message.get("CRUSHFragmenter.saturationSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.saturationSetting.getName(),
                Message.get("CRUSHFragmenter.saturationSetting.displayName"));
        this.minimumFragmentSizeSetting = new SimpleIntegerProperty(this,
                "Minimum fragment size setting",
                CRUSHFragmenter.MINIMUM_FRAGMENT_SIZE_DEFAULT) {
            @Override
            public void set(int newValue) {
                if (newValue >= 1) {
                    super.set(newValue);
                } else {
                    IllegalArgumentException tmpException = new IllegalArgumentException(
                            "The minimum fragment size must be at least 1.");
                    CRUSHFragmenter.LOGGER.log(Level.WARNING, tmpException.toString(), tmpException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            tmpException.toString(),
                            tmpException);
                    //re-throws the exception to properly reset the binding
                    throw tmpException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.minimumFragmentSizeSetting.getName(),
                Message.get("CRUSHFragmenter.minimumFragmentSizeSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.minimumFragmentSizeSetting.getName(),
                Message.get("CRUSHFragmenter.minimumFragmentSizeSetting.displayName"));
        this.ringProtectionSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Ring protection setting",
                CRUSHFragmenter.RING_PROTECTION_OPTION_DEFAULT,
                CRUSHFragmenter.RingProtectionOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    CRUSHFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.ringProtectionSetting.getName(),
                Message.get("CRUSHFragmenter.ringProtectionSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.ringProtectionSetting.getName(),
                Message.get("CRUSHFragmenter.ringProtectionSetting.displayName"));
        this.conservativeModeSetting = new SimpleBooleanProperty(this,
                "Conservative mode setting",
                CRUSHFragmenter.CONSERVATIVE_MODE_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.conservativeModeSetting.getName(),
                Message.get("CRUSHFragmenter.conservativeModeSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.conservativeModeSetting.getName(),
                Message.get("CRUSHFragmenter.conservativeModeSetting.displayName"));
        this.fragmentOutputSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Fragment output setting",
                CRUSHFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT,
                CRUSHFragmenter.FragmentOutputOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    CRUSHFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
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
                Message.get("CRUSHFragmenter.fragmentOutputSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.fragmentOutputSetting.getName(),
                Message.get("CRUSHFragmenter.fragmentOutputSetting.displayName"));
        this.returnUnfragmentedMoleculesSetting = new SimpleBooleanProperty(this,
                "Return unfragmented molecules setting",
                CRUSHFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("CRUSHFragmenter.returnUnfragmentedMoleculesSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("CRUSHFragmenter.returnUnfragmentedMoleculesSetting.displayName"));
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        this.settings.add(this.saturationSetting);
        this.settings.add(this.minimumFragmentSizeSetting);
        this.settings.add(this.ringProtectionSetting);
        this.settings.add(this.conservativeModeSetting);
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
     * Returns the minimum number of heavy atoms of a fragment.
     *
     * @return minimum fragment size
     */
    public int getMinimumFragmentSizeSetting() {
        return this.minimumFragmentSizeSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the minimum number of heavy atoms of a fragment.
     *
     * @return property object of the setting
     */
    public SimpleIntegerProperty minimumFragmentSizeSettingProperty() {
        return this.minimumFragmentSizeSetting;
    }
    //
    /**
     * Returns the currently set option for the protected atoms.
     *
     * @return ring protection option
     */
    public RingProtectionOption getRingProtectionSetting() {
        return (RingProtectionOption) this.ringProtectionSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the protected atoms.
     *
     * @return property object of the setting
     */
    public SimpleIDisplayEnumConstantProperty ringProtectionSettingProperty() {
        return this.ringProtectionSetting;
    }
    //
    /**
     * Returns whether the conservative mode is applied, i.e. a rule with a match on a protected atom is not applied at
     * all.
     *
     * @return true for the conservative mode
     */
    public boolean getConservativeModeSetting() {
        return this.conservativeModeSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the conservative mode.
     *
     * @return property object of the setting
     */
    public SimpleBooleanProperty conservativeModeSettingProperty() {
        return this.conservativeModeSetting;
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
     * Returns whether molecules without cleavage steps are returned as one unchanged fragment.
     *
     * @return true if molecules without cleavage steps are returned unchanged
     */
    public boolean getReturnUnfragmentedMoleculesSetting() {
        return this.returnUnfragmentedMoleculesSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for returning molecules without cleavage steps.
     *
     * @return property object of the setting
     */
    public SimpleBooleanProperty returnUnfragmentedMoleculesSettingProperty() {
        return this.returnUnfragmentedMoleculesSetting;
    }
    //
    /**
     * Returns the names of all CRUSH rules in the order of their numbers (index 0 is rule 1), as given by the authors.
     *
     * @return unmodifiable list of rule names
     */
    public static List<String> getRuleNames() {
        List<String> tmpNames = new ArrayList<>(CRUSHFragmenter.RULES.size());
        for (CleavageRule tmpRule : CRUSHFragmenter.RULES) {
            tmpNames.add(tmpRule.name());
        }
        return Collections.unmodifiableList(tmpNames);
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
     * Sets the minimum number of heavy atoms of a fragment; steps that yield a smaller part are not applied. With the
     * default value 3, rule 19 (Carbamate_2) never applies, since it cuts off a carbonyl group with two heavy atoms.
     *
     * @param aSize minimum fragment size, at least 1
     * @throws IllegalArgumentException if the given size is smaller than 1
     */
    public void setMinimumFragmentSizeSetting(int aSize) throws IllegalArgumentException {
        if (aSize < 1) {
            throw new IllegalArgumentException("The minimum fragment size must be at least 1.");
        }
        this.minimumFragmentSizeSetting.set(aSize);
    }
    //
    /**
     * Sets the option for the protected atoms.
     *
     * @param anOption ring protection option
     * @throws NullPointerException if the given option is null
     */
    public void setRingProtectionSetting(RingProtectionOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given ring protection option is null.");
        this.ringProtectionSetting.set(anOption);
    }
    //
    /**
     * Sets whether the conservative mode is applied. If true (default), a rule with a match on a protected atom is not
     * applied to the fragment at all; if false, only the matches on protected atoms are ignored.
     *
     * @param aBoolean true for the conservative mode
     */
    public void setConservativeModeSetting(boolean aBoolean) {
        this.conservativeModeSetting.set(aBoolean);
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
    //
    /**
     * Sets whether molecules without cleavage steps are returned as one unchanged fragment. If false, no fragments are
     * returned for such molecules.
     *
     * @param aBoolean true if molecules without cleavage steps should be returned unchanged
     */
    public void setReturnUnfragmentedMoleculesSetting(boolean aBoolean) {
        this.returnUnfragmentedMoleculesSetting.set(aBoolean);
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
        return CRUSHFragmenter.ALGORITHM_NAME;
    }
    //
    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("CRUSHFragmenter.displayName");
    }
    //
    @Override
    public IMoleculeFragmenter copy() {
        CRUSHFragmenter tmpCopy = new CRUSHFragmenter();
        tmpCopy.setSaturationSetting(this.getSaturationSetting());
        tmpCopy.setMinimumFragmentSizeSetting(this.getMinimumFragmentSizeSetting());
        tmpCopy.setRingProtectionSetting(this.getRingProtectionSetting());
        tmpCopy.setConservativeModeSetting(this.getConservativeModeSetting());
        tmpCopy.setFragmentOutputSetting(this.getFragmentOutputSetting());
        tmpCopy.setReturnUnfragmentedMoleculesSetting(this.getReturnUnfragmentedMoleculesSetting());
        return tmpCopy;
    }
    //
    @Override
    public void restoreDefaultSettings() {
        this.saturationSetting.set(CRUSHFragmenter.SATURATION_OPTION_DEFAULT);
        this.minimumFragmentSizeSetting.set(CRUSHFragmenter.MINIMUM_FRAGMENT_SIZE_DEFAULT);
        this.ringProtectionSetting.set(CRUSHFragmenter.RING_PROTECTION_OPTION_DEFAULT);
        this.conservativeModeSetting.set(CRUSHFragmenter.CONSERVATIVE_MODE_OPTION_DEFAULT);
        this.fragmentOutputSetting.set(CRUSHFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT);
        this.returnUnfragmentedMoleculesSetting.set(CRUSHFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
    }
    //
    /**
     * {@inheritDoc}
     *
     * <p>The fragments are the final fragments of the route with the maximum number of fragments (see
     * {@link #buildRoute(IAtomContainer)}), i.e. the disconnected parts of a clone of the given molecule in which all
     * bonds cleaved on the route have been cleaved and the open valences saturated according to the saturation setting.
     * With the fragment output option {@link FragmentOutputOption#ALL_ROUTE_FRAGMENTS}, all fragments of the route
     * except its root (the molecule) are returned instead. If no bond is cleaved, the clone (or its disconnected parts,
     * e.g. counter ions, separately) is returned or, depending on the respective setting, no fragment is returned.</p>
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
            CRUSHRouteNode tmpRoot = this.createRoute(tmpMolecule);
            int[] tmpBondLabels = tmpRoot.route.bondLabels();
            boolean tmpIsCleaved = false;
            for (int tmpLabel : tmpBondLabels) {
                if (tmpLabel != 0) {
                    tmpIsCleaved = true;
                    break;
                }
            }
            if (!tmpIsCleaved) {
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
            }
            if (this.fragmentOutputSetting.get().equals(FragmentOutputOption.ALL_ROUTE_FRAGMENTS)) {
                List<CRUSHRouteNode> tmpNodes = tmpRoot.getAllDescendants();
                List<IAtomContainer> tmpFragments = new ArrayList<>(tmpNodes.size());
                for (CRUSHRouteNode tmpNode : tmpNodes) {
                    tmpFragments.add(tmpNode.getFragment());
                }
                return tmpFragments;
            }
            //all bonds of the route are cleaved at once, the molecule is not needed anymore
            List<CRUSHBond> tmpBondsToCleave = new ArrayList<>();
            for (IBond tmpBond : tmpMolecule.bonds()) {
                int tmpLabel = tmpBondLabels[tmpBond.getIndex()];
                if (tmpLabel != 0) {
                    tmpBondsToCleave.add(new CRUSHBond(tmpBond, tmpLabel));
                }
            }
            return CRUSHFragmenter.cleaveBonds(tmpMolecule, tmpBondsToCleave,
                    this.saturationSetting.get().equals(SaturationOption.HYDROGEN_SATURATION));
        } catch (CDKException | RuntimeException anException) {
            throw new IllegalArgumentException(String.format("An error occurred during fragmentation: %s Molecule Name: %s",
                    anException, aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)), anException);
        }
    }
    //
    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        return Objects.isNull(aMolecule) || aMolecule.isEmpty();
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
     * Determines the fragmentation route of the given molecule with the current settings: the decomposition tree with
     * the maximum number of final fragments, see the class documentation. The root of the route is the molecule; its
     * children are the parts of the first cleavage step, or, for molecules with several connected components (e.g.
     * salts), the components, whose children are the parts of their first cleavage steps. The leaves are the final
     * fragments. The fragments of the nodes are created on request, with link atoms or hydrogen atoms according to the
     * saturation setting.
     *
     * @param aMolecule the molecule to fragment, it is not altered
     * @return the root node of the route, its fragment is the prepared clone of the molecule (explicit hydrogen atoms
     *         suppressed, Daylight aromaticity)
     * @throws NullPointerException if the given molecule is null
     * @throws IllegalArgumentException if the molecule cannot be fragmented (see
     *                                  {@link #canBeFragmented(IAtomContainer)}) or if an error occurs during its
     *                                  preparation
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    public CRUSHRouteNode buildRoute(IAtomContainer aMolecule)
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
        return this.createRoute(tmpMolecule);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Package-private methods">
    /**
     * Returns all CRUSH rules in the order of their numbers.
     *
     * @return unmodifiable list of the rules
     */
    static List<CleavageRule> getRules() {
        return CRUSHFragmenter.RULES;
    }
    //
    /**
     * Prepares the given molecule (in place) for the matching of the CRUSH rules: missing implicit hydrogen counts are
     * set, explicit hydrogen atoms are suppressed, a Kekulé structure is assigned if the order of a bond is unset (e.g.
     * the aromatic bonds of a molecule whose kekulization failed on import), and ring membership and aromaticity
     * (Daylight model) are perceived. The Kekulé structure is needed because the aromaticity perception determines the
     * aromaticity from the bond orders, so aromatic rings with unset bond orders would otherwise be regarded as
     * non-aromatic.
     *
     * @param aMolecule the molecule to prepare, it is altered
     * @throws CDKException if implicit hydrogen counts are missing and cannot be determined, or if the order of a bond is
     *                      unset and no Kekulé structure can be assigned
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
        AtomContainerManipulator.suppressHydrogens(aMolecule);
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
        //marks ring atoms and bonds, and applies the Daylight aromaticity model after clearing existing aromaticity flags
        SmartsPattern.prepare(aMolecule);
    }
    //
    /**
     * Returns the protected atoms of the given prepared molecule according to the ring protection setting: atoms in at
     * least two rings (ring set: relevant cycles) and, by default, their neighbours.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @return indices of the protected atoms
     */
    BitSet findProtectedAtoms(IAtomContainer aPreparedMolecule) {
        BitSet tmpAllAtoms = new BitSet(aPreparedMolecule.getAtomCount());
        tmpAllAtoms.set(0, aPreparedMolecule.getAtomCount());
        return this.findProtectedAtoms(tmpAllAtoms, GraphUtil.toAdjList(aPreparedMolecule),
                CRUSHFragmenter.findFusedRingAtoms(aPreparedMolecule));
    }
    //
    /**
     * Determines the cleavage steps of the fragment with the given atoms of the prepared molecule, for every rule
     * separately (see class documentation). A step of one rule is contained once per set of cleaved bonds; if several
     * rules cleave the same bonds, the step is contained once for each of them.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @param anAtomSet atom indices of a connected fragment, e.g. a connected component of the molecule
     * @return the steps sorted by rule number and cleaved bonds
     */
    List<CleavageStep> findCleavageSteps(IAtomContainer aPreparedMolecule, BitSet anAtomSet) {
        GraphUtil.EdgeToBondMap tmpEdgeToBondMap = GraphUtil.EdgeToBondMap.withSpaceFor(aPreparedMolecule);
        int[][] tmpAdjacencyList = GraphUtil.toAdjList(aPreparedMolecule, tmpEdgeToBondMap);
        return this.findCleavageSteps(aPreparedMolecule, anAtomSet, tmpAdjacencyList, tmpEdgeToBondMap,
                CRUSHFragmenter.findFusedRingAtoms(aPreparedMolecule));
    }
    //
    /**
     * Creates the fragment with the given atoms: a clone of the prepared molecule is cleaved at all bonds that connect
     * the atom set to the rest of the molecule, and the part with the atom set is returned. The open valences are
     * saturated with link atoms carrying the given rule numbers or with hydrogen atoms, and the stereo elements are
     * transferred.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}, it is not altered
     * @param anAtomSet atom indices of the fragment
     * @param aBondLabelArray rule numbers by bond index for the link atoms of the cleaved bonds (0: link atoms without
     *                        rule number)
     * @param isHydrogenSaturation true if the open valences should be saturated with hydrogen atoms
     * @return the fragment
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    static IAtomContainer createFragment(IAtomContainer aPreparedMolecule, BitSet anAtomSet, int[] aBondLabelArray,
                                         boolean isHydrogenSaturation) throws CloneNotSupportedException {
        IAtomContainer tmpClone = aPreparedMolecule.clone();
        List<CRUSHBond> tmpBoundaryBonds = new ArrayList<>();
        for (IBond tmpBond : tmpClone.bonds()) {
            if (anAtomSet.get(tmpBond.getBegin().getIndex()) != anAtomSet.get(tmpBond.getEnd().getIndex())) {
                tmpBoundaryBonds.add(new CRUSHBond(tmpBond, aBondLabelArray[tmpBond.getIndex()]));
            }
        }
        IAtom tmpMarkerAtom = tmpClone.getAtom(anAtomSet.nextSetBit(0));
        for (IAtomContainer tmpFragment : CRUSHFragmenter.cleaveBonds(tmpClone, tmpBoundaryBonds, isHydrogenSaturation)) {
            if (tmpFragment.contains(tmpMarkerAtom)) {
                return tmpFragment;
            }
        }
        throw new IllegalStateException("The fragment with the given atoms could not be determined.");
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Determines the route of the given prepared molecule, see {@link #buildRoute(IAtomContainer)}.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}; it becomes the fragment
     *                          of the root and must not be altered afterwards
     * @return the root node of the route
     * @throws CloneNotSupportedException if a fragment for the tie-breaking cannot be created
     */
    private CRUSHRouteNode createRoute(IAtomContainer aPreparedMolecule) throws CloneNotSupportedException {
        RouteSearch tmpSearch = new RouteSearch(aPreparedMolecule);
        int tmpAtomCount = aPreparedMolecule.getAtomCount();
        ConnectedComponents tmpConnectedComponents = new ConnectedComponents(tmpSearch.adjacencyList);
        int[] tmpComponentLabels = tmpConnectedComponents.components();
        List<BitSet> tmpComponents = new ArrayList<>(tmpConnectedComponents.nComponents());
        for (int tmpComponent = 1; tmpComponent <= tmpConnectedComponents.nComponents(); tmpComponent++) {
            BitSet tmpAtomSet = new BitSet(tmpAtomCount);
            for (int i = 0; i < tmpAtomCount; i++) {
                if (tmpComponentLabels[i] == tmpComponent) {
                    tmpAtomSet.set(i);
                }
            }
            tmpComponents.add(tmpAtomSet);
        }
        for (BitSet tmpComponent : tmpComponents) {
            tmpSearch.countFragments(tmpComponent);
        }
        Route tmpRoute = new Route(aPreparedMolecule, tmpSearch.bondLabels,
                this.saturationSetting.get().equals(SaturationOption.HYDROGEN_SATURATION), tmpSearch.isExhaustive);
        BitSet tmpAllAtoms = new BitSet(tmpAtomCount);
        tmpAllAtoms.set(0, tmpAtomCount);
        CRUSHRouteNode tmpRoot = new CRUSHRouteNode(tmpRoute, tmpAllAtoms, null);
        if (tmpComponents.size() == 1) {
            tmpSearch.assemble(tmpRoot);
        } else {
            for (BitSet tmpComponent : tmpComponents) {
                CRUSHRouteNode tmpComponentNode = new CRUSHRouteNode(tmpRoute, tmpComponent, tmpRoot);
                tmpRoot.children.add(tmpComponentNode);
                tmpSearch.assemble(tmpComponentNode);
            }
        }
        return tmpRoot;
    }
    //
    /**
     * Determines the cleavage steps of a fragment that are unique per set of cleaved bonds, with the lowest rule number.
     *
     * @param aPreparedMolecule the prepared molecule
     * @param anAtomSet atom set of the fragment
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @param aFusedRingAtomSet atoms of the molecule in at least two rings
     * @return the unique steps sorted by rule number and cleaved bonds
     */
    private List<CleavageStep> findUniqueSteps(IAtomContainer aPreparedMolecule, BitSet anAtomSet, int[][] anAdjacencyList,
                                               GraphUtil.EdgeToBondMap anEdgeToBondMap, BitSet aFusedRingAtomSet) {
        List<CleavageStep> tmpSteps = this.findCleavageSteps(aPreparedMolecule, anAtomSet, anAdjacencyList,
                anEdgeToBondMap, aFusedRingAtomSet);
        Set<List<Integer>> tmpCleavedBondSets = new HashSet<>();
        List<CleavageStep> tmpUniqueSteps = new ArrayList<>(tmpSteps.size());
        for (CleavageStep tmpStep : tmpSteps) {
            List<Integer> tmpBondList = new ArrayList<>(tmpStep.cleavedBonds().length);
            for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                tmpBondList.add(tmpBondIndex);
            }
            if (tmpCleavedBondSets.add(tmpBondList)) {
                tmpUniqueSteps.add(tmpStep);
            }
        }
        return tmpUniqueSteps;
    }
    //
    /**
     * Determines the cleavage steps of the fragment with the given atoms of the prepared molecule, see
     * {@link #findCleavageSteps(IAtomContainer, BitSet)}. The rules are matched on a copy of the fragment that is
     * saturated with hydrogen atoms at its boundary bonds, with ring membership and aromaticity perceived again; the
     * protected atoms are determined from the rings of the molecule, since ring bonds are never cleaved.
     *
     * @param aPreparedMolecule the prepared molecule
     * @param anAtomSet atom set of the fragment
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @param aFusedRingAtomSet atoms of the molecule in at least two rings
     * @return the steps sorted by rule number and cleaved bonds
     */
    private List<CleavageStep> findCleavageSteps(IAtomContainer aPreparedMolecule, BitSet anAtomSet,
                                                 int[][] anAdjacencyList, GraphUtil.EdgeToBondMap anEdgeToBondMap,
                                                 BitSet aFusedRingAtomSet) {
        int[] tmpCopyToMoleculeIndex = new int[anAtomSet.cardinality()];
        IAtomContainer tmpCopy = CRUSHFragmenter.createSaturatedCopy(aPreparedMolecule, anAtomSet, tmpCopyToMoleculeIndex);
        BitSet tmpProtectedAtoms = this.findProtectedAtoms(anAtomSet, anAdjacencyList, aFusedRingAtomSet);
        boolean tmpIsConservative = this.conservativeModeSetting.get();
        int tmpMinimumSize = this.minimumFragmentSizeSetting.get();
        List<CleavageStep> tmpSteps = new ArrayList<>();
        for (CompiledRule tmpCompiledRule : this.compiledRules) {
            List<int[]> tmpMatches = new ArrayList<>();
            boolean tmpIsAnyMatchProtected = false;
            for (int[] tmpMatch : tmpCompiledRule.pattern().matchAll(tmpCopy)) {
                boolean tmpIsProtected = false;
                for (int tmpCopyIndex : tmpMatch) {
                    if (tmpProtectedAtoms.get(tmpCopyToMoleculeIndex[tmpCopyIndex])) {
                        tmpIsProtected = true;
                        break;
                    }
                }
                if (tmpIsProtected) {
                    tmpIsAnyMatchProtected = true;
                    if (tmpIsConservative) {
                        break;
                    }
                } else {
                    //the mappings reuse their array
                    tmpMatches.add(tmpMatch.clone());
                }
            }
            if (tmpIsConservative && tmpIsAnyMatchProtected) {
                continue;
            }
            Map<List<Integer>, CleavageStep> tmpRuleSteps = new LinkedHashMap<>();
            for (int[] tmpMatch : tmpMatches) {
                int[][] tmpPairs = tmpCompiledRule.cleavedQueryAtomPairs();
                int[] tmpCleavedBonds = new int[tmpPairs.length];
                boolean tmpIsRingBond = false;
                for (int i = 0; i < tmpPairs.length; i++) {
                    int tmpBegin = tmpCopyToMoleculeIndex[tmpMatch[tmpPairs[i][0]]];
                    int tmpEnd = tmpCopyToMoleculeIndex[tmpMatch[tmpPairs[i][1]]];
                    IBond tmpBond = anEdgeToBondMap.get(tmpBegin, tmpEnd);
                    if (tmpBond.isInRing()) {
                        tmpIsRingBond = true;
                        break;
                    }
                    tmpCleavedBonds[i] = tmpBond.getIndex();
                }
                if (tmpIsRingBond) {
                    continue;
                }
                Arrays.sort(tmpCleavedBonds);
                List<Integer> tmpKey = new ArrayList<>(tmpCleavedBonds.length);
                for (int tmpBondIndex : tmpCleavedBonds) {
                    tmpKey.add(tmpBondIndex);
                }
                if (tmpRuleSteps.containsKey(tmpKey)) {
                    continue;
                }
                BitSet[] tmpParts = CRUSHFragmenter.splitFragment(aPreparedMolecule, anAtomSet, tmpCleavedBonds,
                        anAdjacencyList, anEdgeToBondMap);
                boolean tmpIsTooSmall = false;
                for (BitSet tmpPart : tmpParts) {
                    if (CRUSHFragmenter.countHeavyAtoms(aPreparedMolecule, tmpPart) < tmpMinimumSize) {
                        tmpIsTooSmall = true;
                        break;
                    }
                }
                //steps rejected because of their size are kept as null, so that they are not evaluated again
                tmpRuleSteps.put(tmpKey, tmpIsTooSmall ? null
                        : new CleavageStep(tmpCompiledRule.rule().number(), tmpCleavedBonds, tmpParts));
            }
            List<CleavageStep> tmpValidSteps = new ArrayList<>(tmpRuleSteps.size());
            for (CleavageStep tmpStep : tmpRuleSteps.values()) {
                if (tmpStep != null) {
                    tmpValidSteps.add(tmpStep);
                }
            }
            tmpValidSteps.sort((aStep1, aStep2) -> Arrays.compare(aStep1.cleavedBonds(), aStep2.cleavedBonds()));
            tmpSteps.addAll(tmpValidSteps);
        }
        return tmpSteps;
    }
    //
    /**
     * Returns the protected atoms of the given fragment according to the ring protection setting.
     *
     * @param anAtomSet atom set of the fragment
     * @param anAdjacencyList adjacency list of the molecule
     * @param aFusedRingAtomSet atoms of the molecule in at least two rings
     * @return indices of the protected atoms of the fragment
     */
    private BitSet findProtectedAtoms(BitSet anAtomSet, int[][] anAdjacencyList, BitSet aFusedRingAtomSet) {
        BitSet tmpProtectedAtoms = new BitSet(anAdjacencyList.length);
        RingProtectionOption tmpOption = (RingProtectionOption) this.ringProtectionSetting.get();
        if (tmpOption == RingProtectionOption.NONE) {
            return tmpProtectedAtoms;
        }
        for (int i = aFusedRingAtomSet.nextSetBit(0); i >= 0; i = aFusedRingAtomSet.nextSetBit(i + 1)) {
            if (!anAtomSet.get(i)) {
                continue;
            }
            tmpProtectedAtoms.set(i);
            if (tmpOption == RingProtectionOption.FUSED_RING_ATOMS_AND_NEIGHBOURS) {
                for (int tmpNeighbour : anAdjacencyList[i]) {
                    if (anAtomSet.get(tmpNeighbour)) {
                        tmpProtectedAtoms.set(tmpNeighbour);
                    }
                }
            }
        }
        return tmpProtectedAtoms;
    }
    //
    /**
     * Returns the atoms of the given molecule that are contained in at least two rings of its relevant cycles (fused
     * ring systems, spiro and bridgehead atoms), corresponding to the symmetrised smallest set of smallest rings of the
     * RDKit used by the authors. If the relevant cycles cannot be determined, the minimum cycle basis is used.
     *
     * @param aMolecule molecule
     * @return indices of the atoms in at least two rings
     */
    private static BitSet findFusedRingAtoms(IAtomContainer aMolecule) {
        int tmpAtomCount = aMolecule.getAtomCount();
        BitSet tmpFusedRingAtoms = new BitSet(tmpAtomCount);
        Cycles tmpCycles;
        try {
            tmpCycles = Cycles.or(Cycles.relevant(), Cycles.mcb()).find(aMolecule);
        } catch (Intractable anException) {
            CRUSHFragmenter.LOGGER.log(Level.FINE, "Relevant cycles could not be determined", anException);
            tmpCycles = Cycles.mcb(aMolecule);
        }
        int[] tmpRingCounts = new int[tmpAtomCount];
        for (int[] tmpPath : tmpCycles.paths()) {
            //the paths are closed, i.e. the first atom is repeated at the end
            for (int i = 0; i < tmpPath.length - 1; i++) {
                tmpRingCounts[tmpPath[i]]++;
            }
        }
        for (int i = 0; i < tmpAtomCount; i++) {
            if (tmpRingCounts[i] >= 2) {
                tmpFusedRingAtoms.set(i);
            }
        }
        return tmpFusedRingAtoms;
    }
    //
    /**
     * Creates a copy of the fragment with the given atoms of the prepared molecule whose open valences at the boundary
     * bonds are saturated with implicit hydrogen atoms (one per bond order), with ring membership and aromaticity
     * perceived again. Stereo elements are not copied.
     *
     * @param aPreparedMolecule the prepared molecule
     * @param anAtomSet atom set of the fragment
     * @param aCopyToMoleculeIndexArray array that is filled with the atom indices in the molecule of the atoms of the
     *                                  copy
     * @return the copy
     */
    private static IAtomContainer createSaturatedCopy(IAtomContainer aPreparedMolecule, BitSet anAtomSet,
                                                      int[] aCopyToMoleculeIndexArray) {
        IChemObjectBuilder tmpBuilder = aPreparedMolecule.getBuilder();
        IAtomContainer tmpCopy = tmpBuilder.newAtomContainer();
        Map<IAtom, IAtom> tmpAtomMap = new IdentityHashMap<>(anAtomSet.cardinality() * 2);
        int tmpCopyIndex = 0;
        for (int i = anAtomSet.nextSetBit(0); i >= 0; i = anAtomSet.nextSetBit(i + 1)) {
            IAtom tmpAtom = aPreparedMolecule.getAtom(i);
            IAtom tmpAtomCopy;
            if (tmpAtom instanceof IPseudoAtom tmpPseudoAtom) {
                tmpAtomCopy = tmpBuilder.newInstance(IPseudoAtom.class, tmpPseudoAtom.getLabel());
            } else {
                tmpAtomCopy = tmpBuilder.newAtom();
                tmpAtomCopy.setAtomicNumber(tmpAtom.getAtomicNumber());
                tmpAtomCopy.setSymbol(tmpAtom.getSymbol());
                tmpAtomCopy.setFormalCharge(tmpAtom.getFormalCharge());
            }
            tmpAtomCopy.setMassNumber(tmpAtom.getMassNumber());
            int tmpHydrogenCount = Objects.requireNonNullElse(tmpAtom.getImplicitHydrogenCount(), 0);
            for (IBond tmpBond : tmpAtom.bonds()) {
                if (!anAtomSet.get(tmpBond.getOther(tmpAtom).getIndex())) {
                    tmpHydrogenCount += tmpBond.getOrder().numeric();
                }
            }
            tmpAtomCopy.setImplicitHydrogenCount(tmpHydrogenCount);
            tmpCopy.addAtom(tmpAtomCopy);
            tmpAtomMap.put(tmpAtom, tmpCopy.getAtom(tmpCopyIndex));
            aCopyToMoleculeIndexArray[tmpCopyIndex] = i;
            tmpCopyIndex++;
        }
        for (IBond tmpBond : aPreparedMolecule.bonds()) {
            IAtom tmpBegin = tmpAtomMap.get(tmpBond.getBegin());
            IAtom tmpEnd = tmpAtomMap.get(tmpBond.getEnd());
            if (tmpBegin != null && tmpEnd != null) {
                tmpCopy.addBond(tmpBuilder.newInstance(IBond.class, tmpBegin, tmpEnd, tmpBond.getOrder()));
            }
        }
        SmartsPattern.prepare(tmpCopy);
        return tmpCopy;
    }
    //
    /**
     * Returns the parts of the given fragment after the cleavage of the given bonds.
     *
     * @param aPreparedMolecule the prepared molecule
     * @param anAtomSet atom set of the fragment
     * @param aCleavedBondArray indices of the cleaved bonds
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @return atom sets of the parts, ordered by their lowest atom index
     */
    private static BitSet[] splitFragment(IAtomContainer aPreparedMolecule, BitSet anAtomSet, int[] aCleavedBondArray,
                                          int[][] anAdjacencyList, GraphUtil.EdgeToBondMap anEdgeToBondMap) {
        List<BitSet> tmpParts = new ArrayList<>(aCleavedBondArray.length + 1);
        BitSet tmpAssigned = new BitSet(anAdjacencyList.length);
        for (int tmpBondIndex : aCleavedBondArray) {
            IBond tmpBond = aPreparedMolecule.getBond(tmpBondIndex);
            for (int tmpStart : new int[] {tmpBond.getBegin().getIndex(), tmpBond.getEnd().getIndex()}) {
                if (tmpAssigned.get(tmpStart)) {
                    continue;
                }
                BitSet tmpPart = new BitSet(anAdjacencyList.length);
                Deque<Integer> tmpQueue = new ArrayDeque<>();
                tmpQueue.add(tmpStart);
                tmpPart.set(tmpStart);
                while (!tmpQueue.isEmpty()) {
                    int tmpAtomIndex = tmpQueue.poll();
                    for (int tmpNeighbour : anAdjacencyList[tmpAtomIndex]) {
                        if (!anAtomSet.get(tmpNeighbour) || tmpPart.get(tmpNeighbour)) {
                            continue;
                        }
                        int tmpEdgeBondIndex = anEdgeToBondMap.get(tmpAtomIndex, tmpNeighbour).getIndex();
                        boolean tmpIsCleaved = false;
                        for (int tmpCleavedBondIndex : aCleavedBondArray) {
                            if (tmpCleavedBondIndex == tmpEdgeBondIndex) {
                                tmpIsCleaved = true;
                                break;
                            }
                        }
                        if (!tmpIsCleaved) {
                            tmpPart.set(tmpNeighbour);
                            tmpQueue.add(tmpNeighbour);
                        }
                    }
                }
                tmpAssigned.or(tmpPart);
                tmpParts.add(tmpPart);
            }
        }
        tmpParts.sort((aPart1, aPart2) -> Integer.compare(aPart1.nextSetBit(0), aPart2.nextSetBit(0)));
        return tmpParts.toArray(new BitSet[0]);
    }
    //
    /**
     * Returns the number of heavy atoms (no hydrogen and pseudo atoms) of the given atom set.
     *
     * @param aMolecule the molecule
     * @param anAtomSet atom set
     * @return number of heavy atoms
     */
    private static int countHeavyAtoms(IAtomContainer aMolecule, BitSet anAtomSet) {
        int tmpCount = 0;
        for (int i = anAtomSet.nextSetBit(0); i >= 0; i = anAtomSet.nextSetBit(i + 1)) {
            IAtom tmpAtom = aMolecule.getAtom(i);
            if (!(tmpAtom instanceof IPseudoAtom) && tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() > 1) {
                tmpCount++;
            }
        }
        return tmpCount;
    }
    //
    /**
     * Creates the key of a part for the tie-breaking of cleavage steps: the unique SMILES code of the fragment with
     * link atoms carrying the given rule numbers as mass numbers and, if stereochemistry is encoded, additionally the
     * unique SMILES code with stereochemistry and without rule numbers.
     *
     * @param aPreparedMolecule the prepared molecule
     * @param anAtomSet atom set of the part
     * @param aBondLabelArray rule numbers by bond index for the link atoms
     * @param isStereoEncoded whether stereochemistry is encoded
     * @return key
     * @throws CloneNotSupportedException if the fragment cannot be created
     */
    private static String createTieBreakingKey(IAtomContainer aPreparedMolecule, BitSet anAtomSet, int[] aBondLabelArray,
                                               boolean isStereoEncoded) throws CloneNotSupportedException {
        IAtomContainer tmpFragment = CRUSHFragmenter.createFragment(aPreparedMolecule, anAtomSet, aBondLabelArray, false);
        String tmpSmiles = CRUSHFragmenter.createUniqueSmilesWithRuleNumbers(tmpFragment);
        if (!isStereoEncoded || !tmpFragment.stereoElements().iterator().hasNext()) {
            return tmpSmiles;
        }
        //the canonical SMILES code with stereochemistry cannot be created with mass numbers of pseudo atoms
        for (IAtom tmpAtom : tmpFragment.atoms()) {
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAtom.setMassNumber(null);
            }
        }
        return tmpSmiles + " " + Objects.requireNonNullElse(ChemUtil.createUniqueSmiles(tmpFragment, true), "");
    }
    //
    /**
     * Creates the key of a part for the tie-breaking of cleavage steps without stereochemistry, i.e. the result of
     * {@link #createTieBreakingKey(IAtomContainer, BitSet, int[], boolean)} without stereochemistry, from a fragment
     * that is built from the atoms of the part only, so that its cost does not depend on the size of the molecule. The
     * fragment has the same atoms, bonds, and link atoms in the same order as the fragment of
     * {@link #createFragment(IAtomContainer, BitSet, int[], boolean)} (the atoms of the part in the order of their
     * indices, followed by the link atoms in the order of the indices of the cleaved bonds), so that also its unique
     * SMILES code is the same; stereo elements are not needed.
     *
     * @param aPreparedMolecule the prepared molecule, it is not altered
     * @param anAtomSet atom set of the part
     * @param aBondLabelArray rule numbers by bond index for the link atoms
     * @return key
     * @throws CloneNotSupportedException if an atom cannot be cloned
     */
    private static String createTieBreakingKeyWithoutStereo(IAtomContainer aPreparedMolecule, BitSet anAtomSet,
                                                            int[] aBondLabelArray) throws CloneNotSupportedException {
        IChemObjectBuilder tmpBuilder = aPreparedMolecule.getBuilder();
        IAtomContainer tmpFragment = tmpBuilder.newAtomContainer();
        Map<IAtom, IAtom> tmpAtomMap = new IdentityHashMap<>(anAtomSet.cardinality() * 2);
        for (int i = anAtomSet.nextSetBit(0); i >= 0; i = anAtomSet.nextSetBit(i + 1)) {
            IAtom tmpAtom = aPreparedMolecule.getAtom(i);
            tmpFragment.addAtom(tmpAtom.clone());
            tmpAtomMap.put(tmpAtom, tmpFragment.getAtom(tmpFragment.getAtomCount() - 1));
        }
        List<IBond> tmpBoundaryBonds = new ArrayList<>();
        for (IBond tmpBond : aPreparedMolecule.bonds()) {
            IAtom tmpBegin = tmpAtomMap.get(tmpBond.getBegin());
            IAtom tmpEnd = tmpAtomMap.get(tmpBond.getEnd());
            if (tmpBegin != null && tmpEnd != null) {
                IBond tmpBondCopy = tmpBuilder.newInstance(IBond.class, tmpBegin, tmpEnd, tmpBond.getOrder());
                tmpBondCopy.setIsAromatic(tmpBond.isAromatic());
                tmpBondCopy.setIsInRing(tmpBond.isInRing());
                tmpFragment.addBond(tmpBondCopy);
            } else if (tmpBegin != null || tmpEnd != null) {
                tmpBoundaryBonds.add(tmpBond);
            }
        }
        for (IBond tmpBond : tmpBoundaryBonds) {
            IAtom tmpPartAtom = tmpAtomMap.containsKey(tmpBond.getBegin()) ? tmpBond.getBegin() : tmpBond.getEnd();
            CRUSHFragmenter.saturate(tmpFragment, tmpAtomMap.get(tmpPartAtom), tmpBond.getOther(tmpPartAtom),
                    aBondLabelArray[tmpBond.getIndex()], tmpBond.getOrder(), false);
        }
        for (ISingleElectron tmpSingleElectron : aPreparedMolecule.singleElectrons()) {
            IAtom tmpAtom = tmpAtomMap.get(tmpSingleElectron.getAtom());
            if (tmpAtom != null) {
                tmpFragment.addSingleElectron(tmpBuilder.newInstance(ISingleElectron.class, tmpAtom));
            }
        }
        for (ILonePair tmpLonePair : aPreparedMolecule.lonePairs()) {
            IAtom tmpAtom = tmpAtomMap.get(tmpLonePair.getAtom());
            if (tmpAtom != null) {
                tmpFragment.addLonePair(tmpBuilder.newInstance(ILonePair.class, tmpAtom));
            }
        }
        return CRUSHFragmenter.createUniqueSmilesWithRuleNumbers(tmpFragment);
    }
    //
    /**
     * Returns the unique SMILES code of the given fragment with aromatic symbols and the rule numbers of the link atoms
     * as mass numbers.
     *
     * @param aFragment fragment
     * @return unique SMILES code, an empty string if it cannot be created
     */
    private static String createUniqueSmilesWithRuleNumbers(IAtomContainer aFragment) {
        try {
            return SmilesGenerator.create(aFragment, SmiFlavor.Unique | SmiFlavor.AtomicMass | SmiFlavor.UseAromaticSymbols,
                    new int[aFragment.getAtomCount()]);
        } catch (CDKException anException) {
            CRUSHFragmenter.LOGGER.log(Level.FINE, "Unique SMILES code could not be created", anException);
            return Objects.requireNonNullElse(ChemUtil.createUniqueSmiles(aFragment, false), "");
        }
    }
    //
    /**
     * Compares two sorted lists of strings lexicographically.
     *
     * @param aList1 first list
     * @param aList2 second list
     * @return negative, zero, or positive, if the first list is smaller, equal, or larger
     */
    private static int compare(List<String> aList1, List<String> aList2) {
        for (int i = 0; i < Math.min(aList1.size(), aList2.size()); i++) {
            int tmpComparison = aList1.get(i).compareTo(aList2.get(i));
            if (tmpComparison != 0) {
                return tmpComparison;
            }
        }
        return Integer.compare(aList1.size(), aList2.size());
    }
    //
    /**
     * Cleaves the given bonds of the molecule simultaneously, saturates the open valences with link atoms or hydrogen
     * atoms, transfers the stereo elements, and returns the disconnected parts as fragments.
     *
     * @param aMolecule the molecule, it is altered
     * @param aCRUSHBondList the bonds to cleave with the rule numbers of their link atoms
     * @param isHydrogenSaturation true if the open valences should be saturated with hydrogen atoms
     * @return list of fragments
     */
    @SuppressWarnings("rawtypes")
    private static List<IAtomContainer> cleaveBonds(IAtomContainer aMolecule, List<CRUSHBond> aCRUSHBondList,
                                                    boolean isHydrogenSaturation) {
        //the stereo elements are detached before the bonds are removed, because the CDK atom container would otherwise
        //adapt them on bond removal, e.g. by replacing a cleaved ligand of a stereocentre with an implicit hydrogen atom
        List<IStereoElement> tmpOriginalStereoElements = new ArrayList<>();
        for (IStereoElement tmpElement : aMolecule.stereoElements()) {
            tmpOriginalStereoElements.add(tmpElement);
        }
        Set<IAtom> tmpNonStereogenicAtomsOfMolecule = tmpOriginalStereoElements.isEmpty()
                ? Collections.emptySet() : CRUSHFragmenter.findNonStereogenicAtoms(aMolecule);
        aMolecule.setStereoElements(new ArrayList<>(0));
        Map<IBond, CleavedBond> tmpCleavedBondMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                aCRUSHBondList.size(), BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR),
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        for (CRUSHBond tmpCRUSHBond : aCRUSHBondList) {
            IBond tmpBond = tmpCRUSHBond.bond();
            IAtom tmpBegin = tmpBond.getBegin();
            IAtom tmpEnd = tmpBond.getEnd();
            IBond.Order tmpOrder = tmpBond.getOrder();
            aMolecule.removeBond(tmpBond);
            IBond tmpNewBeginBond = CRUSHFragmenter.saturate(aMolecule, tmpBegin, tmpEnd, tmpCRUSHBond.ruleNumber(),
                    tmpOrder, isHydrogenSaturation);
            IBond tmpNewEndBond = CRUSHFragmenter.saturate(aMolecule, tmpEnd, tmpBegin, tmpCRUSHBond.ruleNumber(),
                    tmpOrder, isHydrogenSaturation);
            tmpCleavedBondMap.put(tmpBond, new CleavedBond(tmpBegin, tmpEnd, tmpNewBeginBond, tmpNewEndBond));
        }
        CRUSHFragmenter.transferStereoElements(aMolecule, tmpOriginalStereoElements, tmpCleavedBondMap, isHydrogenSaturation);
        if (isHydrogenSaturation) {
            //the new explicit hydrogen atoms become implicit ones, stereo elements are updated accordingly
            AtomContainerManipulator.suppressHydrogens(aMolecule);
        }
        IAtomContainerSet tmpFragmentSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        List<IAtomContainer> tmpFragments = new ArrayList<>(tmpFragmentSet.getAtomContainerCount());
        for (IAtomContainer tmpFragment : tmpFragmentSet.atomContainers()) {
            CRUSHFragmenter.removeStereoElementsMadeNonStereogenic(tmpFragment, tmpNonStereogenicAtomsOfMolecule);
            tmpFragments.add(tmpFragment);
        }
        return tmpFragments;
    }
    //
    /**
     * Saturates the open valence of the given atom that results from the cleavage of its bond to the given former
     * neighbour, either with a link atom or with a hydrogen atom (explicit for single bonds, so that stereo elements
     * can refer to it, implicit ones for bonds of higher order).
     *
     * @param aMolecule the molecule
     * @param anAtom atom of the cleaved bond
     * @param aFormerNeighbour the other atom of the cleaved bond
     * @param aRuleNumber rule number of the cleaved bond, 0 for a link atom without rule number
     * @param anOrder order of the cleaved bond
     * @param isHydrogenSaturation true if the open valence should be saturated with hydrogen atoms
     * @return the new bond of the atom to its link atom or explicit hydrogen atom, null for implicit hydrogen atoms
     */
    private static IBond saturate(IAtomContainer aMolecule, IAtom anAtom, IAtom aFormerNeighbour, int aRuleNumber,
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
            if (aRuleNumber > 0) {
                tmpLinkAtom.setMassNumber(aRuleNumber);
                tmpLinkAtom.setProperty(CRUSHFragmenter.RULE_NUMBER_PROPERTY_KEY, aRuleNumber);
            }
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
    private static void transferStereoElements(IAtomContainer aMolecule, List<IStereoElement> anOriginalStereoElementList,
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
                tmpTransferredElement = CRUSHFragmenter.transferTetrahedralElement(tmpElement, tmpAtomToCleavedBondsMap,
                        isHydrogenSaturation);
            } else if (tmpElement.getConfigClass() == IStereoElement.CT) {
                tmpTransferredElement = CRUSHFragmenter.transferDoubleBondElement(tmpElement, aCleavedBondMap,
                        isHydrogenSaturation);
            } else {
                //other stereo element types are kept only if they are not affected by a cleavage, see below
                tmpTransferredElement = tmpElement;
            }
            if (tmpTransferredElement != null
                    && CRUSHFragmenter.isInOneComponent(aMolecule, tmpTransferredElement, aCleavedBondMap.keySet(), tmpComponents)) {
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
    private static IStereoElement transferTetrahedralElement(IStereoElement anElement,
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
    private static IStereoElement transferDoubleBondElement(IStereoElement anElement, Map<IBond, CleavedBond> aCleavedBondMap,
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
     * stereogenic in the molecule either, e.g. pseudo-asymmetric centres, are kept unchanged. The rule numbers are not
     * considered, since MORTAR identifies fragments by their SMILES codes without isotope information, in which all
     * link atoms are equal.
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
        Set<IAtom> tmpNonStereogenicAtoms = CRUSHFragmenter.findNonStereogenicAtoms(aFragment);
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
     * Returns the atoms of the given molecule that are not stereogenic according to the CDK class
     * {@link Stereocenters}, which perceives stereocentres by constitutional symmetry. Mass numbers of pseudo atoms,
     * i.e. the rule numbers, are ignored.
     *
     * @param aMolecule molecule
     * @return set of non-stereogenic atoms, or null if the stereocentres cannot be determined
     */
    private static Set<IAtom> findNonStereogenicAtoms(IAtomContainer aMolecule) {
        //rule numbers are given as mass numbers and must be ignored for the symmetry perception
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
            Set<IAtom> tmpNonStereogenicAtoms = new HashSet<>();
            for (int i = 0; i < aMolecule.getAtomCount(); i++) {
                if (!tmpStereocenters.isStereocenter(i)) {
                    tmpNonStereogenicAtoms.add(aMolecule.getAtom(i));
                }
            }
            return tmpNonStereogenicAtoms;
        } catch (RuntimeException anException) {
            CRUSHFragmenter.LOGGER.log(Level.FINE, "Stereocentres could not be determined", anException);
            return null;
        } finally {
            for (int i = 0; i < tmpLinkAtoms.size(); i++) {
                tmpLinkAtoms.get(i).setMassNumber(tmpMassNumbers.get(i));
            }
        }
    }
    //
    /**
     * Creates the 33 CRUSH rules in the order of their numbers. The SMARTS are the reactant sides of the reaction SMIRKS
     * of the authors (file CRUSH_SMIRKS_Rules.txt, identical to Table 1 of the publication), the cleaved bonds are the
     * bonds between mapped atoms that are missing on the product sides; the names are those of the authors.
     *
     * @return list of all rules
     */
    private static List<CleavageRule> createRules() {
        int[][] tmpBond12 = {{1, 2}};
        int[][] tmpBond24 = {{2, 4}};
        List<CleavageRule> tmpRules = new ArrayList<>(33);
        tmpRules.add(new CleavageRule(1, "Acyclic carbonyl-N cleavage", "[#6:1]!@[C!R:2](=[O:3])!@[N:4]", tmpBond24));
        tmpRules.add(new CleavageRule(2, "Ester", "[#6:1]!@[C!R:2](=[O:3])!@[O!R:4][#6:5]", tmpBond24));
        tmpRules.add(new CleavageRule(3, "Tertiary/quaternary aliphatic amine C–N cleavage",
                "[NX3,N+X4:1]([CX4:2])([CX4:3])[CX4:4]", new int[][] {{1, 4}}));
        tmpRules.add(new CleavageRule(4, "Ring tertiary amine substituent disconnect", "[N;D3;R:1]!@[C;!$(C=*):2]",
                tmpBond12));
        tmpRules.add(new CleavageRule(5, "Urea", "[N:1]!@[C!R:2](=[O:3])!@[N:4]", tmpBond24));
        tmpRules.add(new CleavageRule(6, "Ether", "[#6,#1:1][#6:2]!@[O!R:3]!@[#6:4][#6,#1:5]", new int[][] {{3, 4}}));
        tmpRules.add(new CleavageRule(7, "Hetero-substituted alkyl C–X cleavage",
                "[C;$(C(-;@[C,N,O,S])-;@[N,O,S]):1]-[*:2]", tmpBond12));
        tmpRules.add(new CleavageRule(8, "Olefin", "[C:1]!@=[C:2]", tmpBond12));
        tmpRules.add(new CleavageRule(9, "Quaternary nitrogen", "[N+X4&H0:1]!@[C:2]", tmpBond12));
        tmpRules.add(new CleavageRule(10, "Aromatic nitrogen – aliphatic carbon", "[n,n+:1]!@[C:2]", tmpBond12));
        tmpRules.add(new CleavageRule(11, "Lactam nitrogen – aliphatic carbon",
                "[NR:1](@[CR;!$(C=[O,N]):2])(!@[C:3])@[CR:4]=[O:5]", new int[][] {{1, 3}}));
        tmpRules.add(new CleavageRule(12, "Aromatic carbon – aromatic carbon", "[c:1]!@[c:2]", tmpBond12));
        tmpRules.add(new CleavageRule(13, "Sulphonamide", "[S!R:1](=[O:2])(=[O:3])!@[N:4]!@[#6;!$(C=[O,N]):5]",
                new int[][] {{1, 4}}));
        tmpRules.add(new CleavageRule(14, "Aromatic carbon – aliphatic carbon", "[c:1]!@[CX4;!$(C~[!#1&!#6]):2]",
                tmpBond12));
        tmpRules.add(new CleavageRule(15, "Aromatic carbon – aliphatic amine", "[c:1]!@[N:2]", tmpBond12));
        tmpRules.add(new CleavageRule(16, "Alkyne", "[C:1]!@#[C:2]", tmpBond12));
        tmpRules.add(new CleavageRule(17, "Thioether", "[#6;!R:1]!@[S;X2:2]!@[#6;!R:3]", tmpBond12));
        tmpRules.add(new CleavageRule(18, "Carbamate_1", "[O:1]!@[C!R:2](=[O:3])!@[N:4]", tmpBond24));
        tmpRules.add(new CleavageRule(19, "Carbamate_2", "[O;!R:1]!@[C;!R:2](=[O:3])!@[N;!R;!$([N]-[C](=O)-[N]):4]",
                new int[][] {{1, 2}, {2, 4}}));
        tmpRules.add(new CleavageRule(20, "Acyl-sulphamide", "[#6:1][S!R:2](=[O:3])(=[O:4])!@[N:5]!@[C!R:6](=[O:7])",
                new int[][] {{5, 6}}));
        tmpRules.add(new CleavageRule(21, "Sulfacylation-O", "[S!R:1](=[O:2])(=[O:3])!@[O:4]!@[c:5]",
                new int[][] {{1, 4}}));
        tmpRules.add(new CleavageRule(22, "Aryl-C–N cleavage", "[c:1]!@[n:2]", tmpBond12));
        tmpRules.add(new CleavageRule(23, "Inter-Aryl-N–N cleavage", "[n:1]!@[n:2]", tmpBond12));
        tmpRules.add(new CleavageRule(24, "Ar–alkyl cleavage", "[c:1]-[C;!R:2]", tmpBond12));
        tmpRules.add(new CleavageRule(25, "C–C break (aliphatic)", "[C;!R;!$(C=[O,N]):1]!@[C;!R;!$(C=[O,N]):2]",
                tmpBond12));
        tmpRules.add(new CleavageRule(26, "C–hetero (aliphatic)", "[C;!R:1]!@[!#6;!R:2]", tmpBond12));
        tmpRules.add(new CleavageRule(27, "Sulfonate ester cleavage", "[S!R:1](=[O:2])(=[O:3])!@[O:4]!@[C:5]",
                new int[][] {{4, 5}}));
        tmpRules.add(new CleavageRule(28, "Disulfide bond cleavage", "[S:1]-[S:2]", tmpBond12));
        tmpRules.add(new CleavageRule(29, "Ar-vinyl-core cleavage", "[c:1]-[C:2]=[C:3]", tmpBond12));
        tmpRules.add(new CleavageRule(30, "Ar-O cleavage", "[c:1]-[O:2]", tmpBond12));
        tmpRules.add(new CleavageRule(31, "Ar–N exocyclic cleavage", "[c:1]-[N:2]", tmpBond12));
        tmpRules.add(new CleavageRule(32, "Ar-S exocyclic cleavage", "[c:1]-[S:2]", tmpBond12));
        tmpRules.add(new CleavageRule(33, "Ar–sp–sp C cleavage", "[c:1]-[C:2]#[C:3]", tmpBond12));
        return Collections.unmodifiableList(tmpRules);
    }
    //</editor-fold>
}
