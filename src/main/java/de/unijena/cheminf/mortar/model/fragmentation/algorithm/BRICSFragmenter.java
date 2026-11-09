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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Implementation of the BRICS fragmentation ("Breaking of Retrosynthetically Interesting Chemical Substructures")
 * by Degen et al., see the reference below. BRICS cleaves acyclic bonds whose chemical environment corresponds to a
 * common, synthetically accessible chemical motif (e.g. amides, esters, ethers, amines, sulfonamides, biaryls, olefins)
 * and marks every resulting attachment point with a typed link atom. There are 16 link types (L1 to L16) that describe
 * the chemical environment of the attachment point:
 * <ul>
 *     <li>L1: C (carbonyl carbon of amides and esters)</li>
 *     <li>L2: N (amide nitrogen or nitrogen bonded to an aromatic ring, N.pl3)</li>
 *     <li>L3: O (ether and ester oxygen)</li>
 *     <li>L4: C (alkyl carbon of amines, ethers and thioethers)</li>
 *     <li>L5: N (amine nitrogen, N.sp3)</li>
 *     <li>L6: C (carbonyl carbon of acylated rings)</li>
 *     <li>L7: C (olefinic carbon, the double bond is cleaved)</li>
 *     <li>L8: C (alkyl carbon of ring substituents)</li>
 *     <li>L9: n (aromatic nitrogen of pyrrole-like rings)</li>
 *     <li>L10: N (lactam nitrogen)</li>
 *     <li>L11: S (sulfur of thioethers, sulfoxides and sulfones)</li>
 *     <li>L12: S (sulfonyl sulfur of sulfonamides)</li>
 *     <li>L13: C (aliphatic ring carbon next to N, O, or S in the ring, "activated")</li>
 *     <li>L14: c (aromatic ring carbon next to N, O, or S in the ring, "activated")</li>
 *     <li>L15: C (aliphatic ring carbon not next to N, O, or S in the ring)</li>
 *     <li>L16: c (aromatic ring carbon not next to N, O, or S in the ring)</li>
 * </ul>
 *
 * <p><b>Rule set.</b> The cleavage rules are the complete original BRICS rule set, i.e. the 55 bond rules in SMARTS
 * syntax that were used with the shredding program Recore to create the BRICS fragment spaces. They are published in
 * the dissertation of the first author, appendix C.1.1.3 ("BRICS-Shredding-Regelwerk", see references). Each rule
 * specifies the environments of both atoms of a cleavable acyclic bond and the link types they receive. The rule names
 * are those of the dissertation, except for the four rules for ring-substituting alkyl carbons (L8 - L13 to L16),
 * which are named C1a to C1d here because the dissertation names two of them "C2a" and "C2b" a second time. The
 * environments are formulated so that the generation of small terminal fragments is avoided: bonds to hydrogen atoms,
 * halogens, hydroxy, nitro, carboxylate, methoxy, methyl, ethyl and isopropyl groups are not cleaved. They also
 * prevent the creation of acetals, thioacetals, imides, and Michael acceptors, and ring bonds are never cleaved.
 * <br>Three adaptations of the published SMARTS were necessary, all of them restore the intended meaning:
 * <ul>
 *     <li>Recore worked on structures with explicit hydrogen atoms, but the degree primitive 'D' of the rules
 *     counts heavy atom neighbours only (otherwise, rules like [C;!D4;!D1] for alkyl carbons would never match).
 *     Explicit hydrogen atoms are therefore suppressed before matching, and '!H' (a non-hydrogen atom in the
 *     original rule C7) is written as '!#1', because 'H' denotes the hydrogen count in Daylight SMARTS.</li>
 *     <li>The atom expression [C,O,S;!@] in rule C7 contains a misplaced bond primitive, '!@' is removed.</li>
 *     <li>The right atom of rule C7 is given in the dissertation as "repetition of the expression on the left side".</li>
 * </ul>
 * Two rules of the shredding rule set have no entry in the link compatibility matrix of the BRICS fragment space:
 * rule N6 (acylated aromatic nitrogen atoms, L1 - L9) and rule O3 (esters of phosphorus and sulfur acids). Both are
 * applied here, since they define cleavable bonds. Rule O3 assigns no link type to the phosphorus or sulfur side,
 * which is therefore marked by an untyped link atom (link type {@link #UNTYPED_LINK}). The filters that were
 * additionally applied for the compilation of the fragment spaces (fragments with at most 16 heavy atoms, no rings
 * with more than eight atoms, no reactive or toxic groups, no molecules without cleavable bonds) are not applied,
 * because they curate a fragment library and do not define the decomposition of a molecule.
 *
 * <p><b>Algorithm.</b> As in the original shredding procedure, all cleavable bonds are first identified on the input
 * molecule and then cleaved simultaneously:
 * <ol>
 *     <li>Explicit hydrogen atoms are suppressed, and ring membership and aromaticity are perceived with the Daylight
 *     model, which is the model that the SMARTS semantics of the rules are based on.</li>
 *     <li>Every acyclic, non-aromatic single or double bond is tested against the rules in their published order, and
 *     the first matching rule determines the link types. The atom environments are evaluated lazily and only once per
 *     atom, i.e. only for atoms of candidate bonds. The environments of the rules are disjoint, so a bond is never
 *     matched by two rules with different link types (none were found in the test data sets).</li>
 *     <li>Each cleaved bond is replaced by two bonds of the same order to new link atoms, one on each side. A link atom
 *     is a pseudo atom carrying the link type of the atom it is bonded to, as mass number (resulting in the common
 *     SMILES notation of BRICS fragments, e.g. [16*]c1ccccc1, if the mass number is encoded) and as property
 *     {@link #LINK_TYPE_PROPERTY_KEY}. Alternatively, open valences can be saturated with hydrogen atoms.</li>
 *     <li>Stereo elements are transferred to the fragments. Tetrahedral and double bond stereo elements whose ligands
 *     were cleaved off refer to the link atom (or hydrogen atom) at the former position afterwards. Stereo elements
 *     whose focus was cleaved and stereo elements that lost their stereogenicity by the fragmentation (e.g. with two
 *     link atoms on a former stereocentre) are removed.</li>
 *     <li>The disconnected parts are returned as fragments; every heavy atom of the molecule is contained in exactly
 *     one fragment, so the frequencies of identical fragments in a molecule are retained. Disconnected parts of the
 *     input molecule, e.g. counter ions, are returned as separate fragments.</li>
 * </ol>
 *
 * <p><b>Differences to the BRICS implementation of the RDKit.</b> The RDKit implementation, which is widely used,
 * deviates from the original rule set in several aspects. It uses less specific atom environments (so e.g. methyl,
 * ethyl, methoxy, trifluoromethyl and acetyl groups, N-methyl amides and benzylic carbon atoms are cleaved off), it
 * merges the link type L2 (amide and aniline nitrogen) into L5 (amine nitrogen), and it does not cleave the bonds of
 * the link type combinations L13 - L13 and L15 - L15 (e.g. bicyclohexyl). Its function BRICSDecompose applies the
 * rules recursively to fragments and returns the set of unique fragment SMILES, so the multiplicity of fragments in
 * a molecule is lost. Fragments of this implementation are therefore usually larger than the RDKit BRICS fragments.
 *
 * <p><b>Note on link types in MORTAR.</b> MORTAR identifies fragments by their unique SMILES code without isotope
 * information, so link atoms are written as '*' and fragments that differ only in their link types are counted as
 * the same fragment.
 *
 * <p>References:
 * <ul>
 *     <li>J. Degen, C. Wegscheid-Gerlach, A. Zaliani, M. Rarey, "On the Art of Compiling and Using 'Drug-Like'
 *     Chemical Fragment Spaces", ChemMedChem 2008, 3, 1503-1507, <a href="https://doi.org/10.1002/cmdc.200800178">
 *     doi:10.1002/cmdc.200800178</a></li>
 *     <li>J. Degen, "Entwicklung eines strukturbasierten de-novo-Design-Verfahrens und chemischer Fragmenträume für
 *     den rationalen Wirkstoffentwurf", dissertation, Universität Hamburg, 2009, chapter 6.1 and appendix C.1.1,
 *     <a href="https://ediss.sub.uni-hamburg.de/handle/ediss/2543">ediss.sub.uni-hamburg.de/handle/ediss/2543</a></li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class BRICSFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum SaturationOption">
    /**
     * Enumeration of the ways to saturate the open valences that result from the cleavage of a bond.
     */
    public enum SaturationOption implements IDisplayEnum {
        /**
         * Every open valence is saturated with a link atom (pseudo atom) that carries the BRICS link type of the
         * attachment point.
         */
        LINK_ATOMS(Message.get("BRICSFragmenter.SaturationOption.linkAtoms.displayName"),
                Message.get("BRICSFragmenter.SaturationOption.linkAtoms.tooltip")),
        /**
         * Every open valence is saturated with (implicit) hydrogen atoms, i.e. both atoms of a cleaved double bond
         * (rule C7) receive two hydrogen atoms.
         */
        HYDROGEN_SATURATION(Message.get("BRICSFragmenter.SaturationOption.hydrogenSaturation.displayName"),
                Message.get("BRICSFragmenter.SaturationOption.hydrogenSaturation.tooltip"));
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
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm.
     */
    public static final String ALGORITHM_NAME = "BRICS";
    /**
     * Default option for the saturation of open valences.
     */
    public static final SaturationOption SATURATION_OPTION_DEFAULT = SaturationOption.LINK_ATOMS;
    /**
     * Default setting for returning molecules without cleavable bonds as one unchanged fragment.
     */
    public static final boolean RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT = true;
    /**
     * Key of the atom property of link atoms that holds their BRICS link type (Integer, 1 to 16, or
     * {@link #UNTYPED_LINK}).
     */
    public static final String LINK_TYPE_PROPERTY_KEY = "BRICSFragmenter.LinkType";
    /**
     * Link type of the link atoms on the phosphorus or sulfur side of a bond cleaved by rule O3, which has no link
     * type in the original rule set.
     */
    public static final int UNTYPED_LINK = 0;
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(BRICSFragmenter.class.getName());
    /**
     * Environment of link type L13, aliphatic ring carbon next to N, O, or S in the ring.
     */
    private static final String L13 = "[C;D3;R;$([C]@[N,n,S,s,O,o])]";
    /**
     * Environment of link type L14, aromatic ring carbon next to N, O, or S in the ring.
     */
    private static final String L14 = "[c;X3;R;$([c]@[N,n,S,s,O,o])]";
    /**
     * Environment of link type L15, aliphatic ring carbon not next to N, O, or S in the ring.
     */
    private static final String L15 = "[C;D3;R;!$([C]@[N,n,S,s,O,o])]";
    /**
     * Environment of link type L16, aromatic ring carbon not next to N, O, or S in the ring.
     */
    private static final String L16 = "[c;X3;R;!$([c]@[N,n,S,s,O,o])]";
    /**
     * Environments of the ring carbon link types L13 to L16, in this order.
     */
    private static final String[] RING_CARBON_ENVIRONMENTS = {L13, L14, L15, L16};
    /**
     * Link types L13 to L16, in the order of {@link #RING_CARBON_ENVIRONMENTS}.
     */
    private static final int[] RING_CARBON_LINK_TYPES = {13, 14, 15, 16};
    /**
     * Letters of the rule names for the ring carbon link types L13 to L16.
     */
    private static final String[] RING_CARBON_RULE_SUFFIXES = {"a", "b", "c", "d"};
    /**
     * All BRICS bond rules in the order of their publication.
     */
    private static final List<BondRule> RULES = BRICSFragmenter.createRules();
    //</editor-fold>
    //
    //<editor-fold desc="Private records">
    /**
     * A BRICS bond rule: an acyclic bond of the given order between an atom matching the left environment and an atom
     * matching the right environment is cleaved, and the atoms receive the given link types.
     *
     * @param name name of the rule in the original publication
     * @param leftEnvironment SMARTS of the left atom
     * @param order order of the bond, single or double
     * @param rightEnvironment SMARTS of the right atom
     * @param leftLinkType link type of the left atom
     * @param rightLinkType link type of the right atom
     */
    private record BondRule(String name, String leftEnvironment, IBond.Order order, String rightEnvironment,
                            int leftLinkType, int rightLinkType) {
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
     * A bond that is cleaved by BRICS.
     *
     * @param bond the bond in the prepared molecule
     * @param ruleName name of the rule that matched the bond
     * @param beginLinkType link type of the begin atom of the bond
     * @param endLinkType link type of the end atom of the bond
     */
    record BRICSBond(IBond bond, String ruleName, int beginLinkType, int endLinkType) {
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
     * Setting for the saturation of open valences.
     */
    private final SimpleIDisplayEnumConstantProperty saturationSetting;
    /**
     * Setting for returning molecules without cleavable bonds as one unchanged fragment.
     */
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
    public BRICSFragmenter() {
        Map<String, Integer> tmpEnvironmentIndexMap = new LinkedHashMap<>();
        this.compiledRules = new CompiledRule[BRICSFragmenter.RULES.size()];
        for (int i = 0; i < BRICSFragmenter.RULES.size(); i++) {
            BondRule tmpRule = BRICSFragmenter.RULES.get(i);
            int tmpLeftIndex = tmpEnvironmentIndexMap.computeIfAbsent(tmpRule.leftEnvironment(), aKey -> tmpEnvironmentIndexMap.size());
            int tmpRightIndex = tmpEnvironmentIndexMap.computeIfAbsent(tmpRule.rightEnvironment(), aKey -> tmpEnvironmentIndexMap.size());
            this.compiledRules[i] = new CompiledRule(tmpRule, tmpLeftIndex, tmpRightIndex);
        }
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        this.environmentPatterns = new DfPattern[tmpEnvironmentIndexMap.size()];
        for (Map.Entry<String, Integer> tmpEntry : tmpEnvironmentIndexMap.entrySet()) {
            QueryAtomContainer tmpQuery = new QueryAtomContainer(tmpBuilder);
            if (!Smarts.parse(tmpQuery, tmpEntry.getKey())) {
                throw new IllegalStateException(String.format("BRICS environment SMARTS %s cannot be parsed: %s",
                        tmpEntry.getKey(), Smarts.getLastErrorMesg()));
            }
            this.environmentPatterns[tmpEntry.getValue()] = DfPattern.findSubstructure(tmpQuery);
        }
        int tmpNumberOfSettings = 2;
        int tmpInitialCapacityForSettingNameMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.saturationSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Fragment saturation setting",
                BRICSFragmenter.SATURATION_OPTION_DEFAULT,
                BRICSFragmenter.SaturationOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    BRICSFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
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
                Message.get("BRICSFragmenter.saturationSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.saturationSetting.getName(),
                Message.get("BRICSFragmenter.saturationSetting.displayName"));
        this.returnUnfragmentedMoleculesSetting = new SimpleBooleanProperty(this,
                "Return unfragmented molecules setting",
                BRICSFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("BRICSFragmenter.returnUnfragmentedMoleculesSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("BRICSFragmenter.returnUnfragmentedMoleculesSetting.displayName"));
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        this.settings.add(this.saturationSetting);
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
     * Returns the names of all BRICS bond rules in the order in which they are applied.
     *
     * @return unmodifiable list of rule names
     */
    public static List<String> getRuleNames() {
        List<String> tmpNames = new ArrayList<>(BRICSFragmenter.RULES.size());
        for (BondRule tmpRule : BRICSFragmenter.RULES) {
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
     * Sets whether molecules without cleavable bonds are returned as one unchanged fragment. If false, no fragments
     * are returned for such molecules.
     *
     * @param aBoolean true if molecules without cleavable bonds should be returned unchanged
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
        return BRICSFragmenter.ALGORITHM_NAME;
    }
    //
    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("BRICSFragmenter.displayName");
    }
    //
    @Override
    public IMoleculeFragmenter copy() {
        BRICSFragmenter tmpCopy = new BRICSFragmenter();
        tmpCopy.setSaturationSetting(this.getSaturationSetting());
        tmpCopy.setReturnUnfragmentedMoleculesSetting(this.getReturnUnfragmentedMoleculesSetting());
        return tmpCopy;
    }
    //
    @Override
    public void restoreDefaultSettings() {
        this.saturationSetting.set(BRICSFragmenter.SATURATION_OPTION_DEFAULT);
        this.returnUnfragmentedMoleculesSetting.set(BRICSFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
    }
    //
    /**
     * {@inheritDoc}
     *
     * <p>The fragments are the disconnected parts of a clone of the given molecule, in which all BRICS bonds have been
     * cleaved and the open valences saturated according to the saturation setting. If the molecule contains no
     * cleavable bond, the clone (or its disconnected parts, e.g. counter ions, separately) is returned or, depending on
     * the respective setting, no fragment is returned.</p>
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
            List<BRICSBond> tmpBRICSBonds = this.findBRICSBonds(tmpMolecule);
            if (!tmpBRICSBonds.isEmpty()) {
                return this.cleaveBonds(tmpMolecule, tmpBRICSBonds);
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
    //<editor-fold desc="Package-private methods">
    /**
     * Prepares the given molecule (in place) for the matching of the BRICS rules: missing implicit hydrogen counts are
     * set, explicit hydrogen atoms are suppressed, and ring membership and aromaticity (Daylight model) are perceived.
     *
     * @param aMolecule the molecule to prepare, it is altered
     * @throws CDKException if implicit hydrogen counts are missing and cannot be determined
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
        //marks ring atoms and bonds, and applies the Daylight aromaticity model after clearing existing aromaticity flags
        SmartsPattern.prepare(aMolecule);
    }
    //
    /**
     * Identifies all bonds of the given prepared molecule that are cleaved by BRICS.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @return the bonds to cleave with their link types, in the order of the bonds in the molecule
     */
    List<BRICSBond> findBRICSBonds(IAtomContainer aPreparedMolecule) {
        //lazily filled results of the environment patterns per atom index, 0: unknown, 1: match, 2: no match
        byte[][] tmpEnvironmentMatches = new byte[this.environmentPatterns.length][];
        int tmpAtomCount = aPreparedMolecule.getAtomCount();
        List<BRICSBond> tmpBRICSBonds = new ArrayList<>();
        for (IBond tmpBond : aPreparedMolecule.bonds()) {
            IBond.Order tmpOrder = tmpBond.getOrder();
            if (tmpBond.isInRing() || tmpBond.isAromatic()
                    || (tmpOrder != IBond.Order.SINGLE && tmpOrder != IBond.Order.DOUBLE)) {
                continue;
            }
            IAtom tmpBegin = tmpBond.getBegin();
            IAtom tmpEnd = tmpBond.getEnd();
            for (CompiledRule tmpCompiledRule : this.compiledRules) {
                BondRule tmpRule = tmpCompiledRule.rule();
                if (tmpRule.order() != tmpOrder) {
                    continue;
                }
                if (this.matches(tmpCompiledRule.leftEnvironmentIndex(), tmpBegin, tmpAtomCount, tmpEnvironmentMatches)
                        && this.matches(tmpCompiledRule.rightEnvironmentIndex(), tmpEnd, tmpAtomCount, tmpEnvironmentMatches)) {
                    tmpBRICSBonds.add(new BRICSBond(tmpBond, tmpRule.name(), tmpRule.leftLinkType(), tmpRule.rightLinkType()));
                    break;
                }
                if (this.matches(tmpCompiledRule.leftEnvironmentIndex(), tmpEnd, tmpAtomCount, tmpEnvironmentMatches)
                        && this.matches(tmpCompiledRule.rightEnvironmentIndex(), tmpBegin, tmpAtomCount, tmpEnvironmentMatches)) {
                    tmpBRICSBonds.add(new BRICSBond(tmpBond, tmpRule.name(), tmpRule.rightLinkType(), tmpRule.leftLinkType()));
                    break;
                }
            }
        }
        return tmpBRICSBonds;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
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
     * Cleaves the given bonds of the prepared molecule simultaneously, saturates the open valences according to the
     * saturation setting, transfers the stereo elements, and returns the disconnected parts as fragments.
     *
     * @param aMolecule the prepared molecule, it is altered
     * @param aBRICSBondList the bonds to cleave
     * @return list of fragments
     */
    @SuppressWarnings("rawtypes")
    private List<IAtomContainer> cleaveBonds(IAtomContainer aMolecule, List<BRICSBond> aBRICSBondList) {
        boolean tmpIsHydrogenSaturation = this.saturationSetting.get().equals(SaturationOption.HYDROGEN_SATURATION);
        //the stereo elements are detached before the bonds are removed, because the CDK atom container would otherwise
        //adapt them on bond removal, e.g. by replacing a cleaved ligand of a stereocentre with an implicit hydrogen atom
        List<IStereoElement> tmpOriginalStereoElements = new ArrayList<>();
        for (IStereoElement tmpElement : aMolecule.stereoElements()) {
            tmpOriginalStereoElements.add(tmpElement);
        }
        Set<IAtom> tmpNonStereogenicAtomsOfMolecule = tmpOriginalStereoElements.isEmpty()
                ? Collections.emptySet() : BRICSFragmenter.findNonStereogenicAtoms(aMolecule);
        aMolecule.setStereoElements(new ArrayList<>(0));
        Map<IBond, CleavedBond> tmpCleavedBondMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                aBRICSBondList.size(), BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR),
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        for (BRICSBond tmpBRICSBond : aBRICSBondList) {
            IBond tmpBond = tmpBRICSBond.bond();
            IAtom tmpBegin = tmpBond.getBegin();
            IAtom tmpEnd = tmpBond.getEnd();
            IBond.Order tmpOrder = tmpBond.getOrder();
            aMolecule.removeBond(tmpBond);
            IBond tmpNewBeginBond = this.saturate(aMolecule, tmpBegin, tmpEnd, tmpBRICSBond.beginLinkType(), tmpOrder,
                    tmpIsHydrogenSaturation);
            IBond tmpNewEndBond = this.saturate(aMolecule, tmpEnd, tmpBegin, tmpBRICSBond.endLinkType(), tmpOrder,
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
            BRICSFragmenter.removeStereoElementsMadeNonStereogenic(tmpFragment, tmpNonStereogenicAtomsOfMolecule);
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
     * @param aLinkType link type of the atom
     * @param anOrder order of the cleaved bond
     * @param isHydrogenSaturation true if the open valence should be saturated with hydrogen atoms
     * @return the new bond of the atom to its link atom or explicit hydrogen atom, null for implicit hydrogen atoms
     */
    private IBond saturate(IAtomContainer aMolecule, IAtom anAtom, IAtom aFormerNeighbour, int aLinkType,
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
            if (aLinkType != BRICSFragmenter.UNTYPED_LINK) {
                tmpLinkAtom.setMassNumber(aLinkType);
            }
            tmpLinkAtom.setProperty(BRICSFragmenter.LINK_TYPE_PROPERTY_KEY, aLinkType);
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
                    && BRICSFragmenter.isInOneComponent(aMolecule, tmpTransferredElement, aCleavedBondMap.keySet(), tmpComponents)) {
                tmpTransferredElements.add(tmpTransferredElement);
            }
        }
        aMolecule.setStereoElements(tmpTransferredElements);
    }
    //
    /**
     * Transfers a tetrahedral stereo element, see {@link #transferStereoElements(IAtomContainer, Map, boolean)}.
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
     * Transfers a double bond stereo element, see {@link #transferStereoElements(IAtomContainer, Map, boolean)}.
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
     * link types are not considered, since MORTAR identifies fragments by their SMILES codes without isotope
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
        Set<IAtom> tmpNonStereogenicAtoms = BRICSFragmenter.findNonStereogenicAtoms(aFragment);
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
     * i.e. the link types, are ignored.
     *
     * @param aMolecule molecule
     * @return set of non-stereogenic atoms, or null if the stereocentres cannot be determined
     */
    private static Set<IAtom> findNonStereogenicAtoms(IAtomContainer aMolecule) {
        //link types are given as mass numbers and must be ignored for the symmetry perception
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
            BRICSFragmenter.LOGGER.log(Level.FINE, "Stereocentres could not be determined", anException);
            return null;
        } finally {
            for (int i = 0; i < tmpLinkAtoms.size(); i++) {
                tmpLinkAtoms.get(i).setMassNumber(tmpMassNumbers.get(i));
            }
        }
    }
    //
    /**
     * Creates the BRICS bond rules in the order of their publication in the dissertation of J. Degen, appendix
     * C.1.1.3. The SMARTS are the published ones with the adaptations described in the class documentation.
     *
     * @return list of all rules
     */
    private static List<BondRule> createRules() {
        IBond.Order tmpSingle = IBond.Order.SINGLE;
        //alkylated ring carbon atoms (named C1a, C1b, C2a, C2b in the dissertation)
        String tmpAlkylL8 = "[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]";
        //acylated ring carbon atoms
        String tmpAcylL6 = "[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]";
        //olefins (without Michael acceptors), '!H' -> '!#1' and '!@' removed from [C,O,S;!@]
        String tmpOlefinL7 = "[#6;!D1;!$([#6]-;!@[#6]=;!@[C,O,S]);$([#6][!#1;!D1]!=[!#1;!D1])]";
        //ethers (no lactones, acetals, or Michael acceptors)
        String tmpEtherCarbonL4 = "[C;!D4;!D1;!R;!$(C(=O));$(C([O;D2])[#6;!D1]~[!D1])]";
        String tmpEtherOxygenL3 = "[O;D2;$(O([C;!D1;!R])[#6;!D1][!D1]);!$(OCO);!$(OC=*);!$(O[P,S])]";
        String tmpRingEtherOxygenL3 = "[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]";
        //esters, urethanes (no anhydrides or Michael acceptors)
        String tmpEsterCarbonL1 = "[C;!D4;!D1;!R;$(C(=O));$(C([O;D2])[#6;!D1]!=[!D1])]";
        String tmpEsterOxygenL3 = "[O;D2;$(O(C(=O))[#6;!D1][#6;!D1]);!$(OCO);!$(O[P,S])]";
        String tmpPhosphorusSulfurEsterOxygenL3 = "[O;D2;$(O([P,S])[#6;!D1][!D1]);!$(OCO);!$(OC=*)]";
        //amides, sulfonamides
        String tmpAmideCarbonL1 = "[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!S;!D1][!D1])]";
        String tmpAmideNitrogenL2 = "[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!S;!D1][!D1])]";
        String tmpSulfonamideSulfurL12 = "[S;X4;$(S(=O)(=O))]";
        String tmpSulfonamideNitrogenL2 = "[#7;X3;!D1;$([#7](S(=O)(=O))[#6;!D1][!D1]);!$([#7]C(=O))]";
        String tmpAcylatedAromaticNitrogenCarbonL1 = "[C;!D4;!D1;!R;$(C(=O));$(C([n])[!S;!D1][!D1])]";
        //primary and secondary amines
        String tmpAmineCarbonL4 = "[C;!D4;!D1;!R;!$(C(=O));$(C([N;X3;!D1])[#6;!D1][!D1])]";
        String tmpAmineNitrogenL5 = "[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([C;!D1;!R])[#6;!D1]!=[!D1])]";
        String tmpRingAmineNitrogen = "[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]";
        //pyrrole-like aromatic nitrogen atoms
        String tmpPyrroleNitrogenL9 = "[n;!$(n@[#6](=O))]";
        String tmpPyrroleAlkylL8 = "[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([n])!=[#6;!D1]~[!D1])]";
        //lactam nitrogen atoms
        String tmpLactamNitrogenL10 = "[#7;R;D3;$([#7]@[#6](=O))]";
        String tmpLactamAlkylL8 = "[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([#7;R;D3])!=[#6;!D1]~[!D1])]";
        String tmpLactamAcylL1 = "[C;!D4;!D1;!R;$(C(=O));$(C([#7;R;D3])[!S;!D1][!D1])]";
        //thioethers (no thiolactones or thioacetals)
        String tmpThioetherCarbonL4 = "[C;!D4;!D1;!R;!$(C(=S));$(C([S;D2])[#6;!D1]~[!D1])]";
        String tmpThioetherSulfurL11 = "[S;D2;$(S([C;!D1;!R])[#6;!D1][!D1]);!$(SCS);!$(SC=*);!$(S[P,S])]";
        String tmpRingThioetherSulfurL11 = "[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]";
        //sulfoxides, sulfones (no sulfonamides)
        String tmpSulfoxideCarbonL4 = "[C;!D4;!D1;!R;$(C(S=O)[#6;!D1][!D1])]";
        String tmpSulfoxideSulfurL11 = "[S;$(S(=O)([C;!D1;!R])[#6;!D1][!D1])]";
        String tmpRingSulfoxideSulfurL11 = "[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]";
        //
        List<BondRule> tmpRules = new ArrayList<>(55);
        //alkylated ring carbon atoms, L8 - L13..L16
        for (int i = 0; i < 4; i++) {
            tmpRules.add(new BondRule("C1" + BRICSFragmenter.RING_CARBON_RULE_SUFFIXES[i], tmpAlkylL8, tmpSingle,
                    BRICSFragmenter.RING_CARBON_ENVIRONMENTS[i], 8, BRICSFragmenter.RING_CARBON_LINK_TYPES[i]));
        }
        //acyclic bonds between ring carbon atoms, C2a-d: L13 - L13..L16, C3a-c: L14 - L14..L16, C4a-b: L15 - L15..L16,
        //C5: L16 - L16
        for (int i = 0; i < 4; i++) {
            for (int j = i; j < 4; j++) {
                String tmpName = (i == 3) ? "C5" : "C" + (i + 2) + BRICSFragmenter.RING_CARBON_RULE_SUFFIXES[j - i];
                tmpRules.add(new BondRule(tmpName, BRICSFragmenter.RING_CARBON_ENVIRONMENTS[i], tmpSingle,
                        BRICSFragmenter.RING_CARBON_ENVIRONMENTS[j], BRICSFragmenter.RING_CARBON_LINK_TYPES[i],
                        BRICSFragmenter.RING_CARBON_LINK_TYPES[j]));
            }
        }
        //acylated ring carbon atoms, L6 - L13..L16
        BRICSFragmenter.addRingCarbonRules(tmpRules, "C6", tmpAcylL6, 6, null);
        //olefins, L7 = L7
        tmpRules.add(new BondRule("C7", tmpOlefinL7, IBond.Order.DOUBLE, tmpOlefinL7, 7, 7));
        //ethers, L4 - L3 and L13..L16 - L3
        tmpRules.add(new BondRule("O1", tmpEtherCarbonL4, tmpSingle, tmpEtherOxygenL3, 4, 3));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "O1", tmpRingEtherOxygenL3, -1, new int[] {3, 3, 3, 3});
        //esters, L1 - L3
        tmpRules.add(new BondRule("O2", tmpEsterCarbonL1, tmpSingle, tmpEsterOxygenL3, 1, 3));
        //esters of phosphorus and sulfur acids, untyped - L3
        tmpRules.add(new BondRule("O3", "[P,S]", tmpSingle, tmpPhosphorusSulfurEsterOxygenL3, BRICSFragmenter.UNTYPED_LINK, 3));
        //amides, L1 - L2
        tmpRules.add(new BondRule("N1", tmpAmideCarbonL1, tmpSingle, tmpAmideNitrogenL2, 1, 2));
        //sulfonamides, L12 - L2
        tmpRules.add(new BondRule("N3", tmpSulfonamideSulfurL12, tmpSingle, tmpSulfonamideNitrogenL2, 12, 2));
        //acylated aromatic nitrogen atoms, L1 - L9
        tmpRules.add(new BondRule("N6", tmpAcylatedAromaticNitrogenCarbonL1, tmpSingle, tmpPyrroleNitrogenL9, 1, 9));
        //amines, L4 - L5, L13 - L5, L14 - L2, L15 - L5, L16 - L2
        tmpRules.add(new BondRule("N2", tmpAmineCarbonL4, tmpSingle, tmpAmineNitrogenL5, 4, 5));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "N2", tmpRingAmineNitrogen, -1, new int[] {5, 2, 5, 2});
        //pyrrole-like aromatic nitrogen atoms, L8 - L9 and L13..L16 - L9
        tmpRules.add(new BondRule("N4", tmpPyrroleAlkylL8, tmpSingle, tmpPyrroleNitrogenL9, 8, 9));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "N5", tmpPyrroleNitrogenL9, -1, new int[] {9, 9, 9, 9});
        //lactam nitrogen atoms, L8 - L10, L13..L16 - L10, and L1 - L10
        tmpRules.add(new BondRule("N7", tmpLactamAlkylL8, tmpSingle, tmpLactamNitrogenL10, 8, 10));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "N8", tmpLactamNitrogenL10, -1, new int[] {10, 10, 10, 10});
        tmpRules.add(new BondRule("N9", tmpLactamAcylL1, tmpSingle, tmpLactamNitrogenL10, 1, 10));
        //thioethers, L4 - L11 and L13..L16 - L11
        tmpRules.add(new BondRule("S1", tmpThioetherCarbonL4, tmpSingle, tmpThioetherSulfurL11, 4, 11));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "S1", tmpRingThioetherSulfurL11, -1, new int[] {11, 11, 11, 11});
        //sulfoxides and sulfones, L4 - L11 and L13..L16 - L11
        tmpRules.add(new BondRule("S2", tmpSulfoxideCarbonL4, tmpSingle, tmpSulfoxideSulfurL11, 4, 11));
        BRICSFragmenter.addRingCarbonRules(tmpRules, "S2", tmpRingSulfoxideSulfurL11, -1, new int[] {11, 11, 11, 11});
        return Collections.unmodifiableList(tmpRules);
    }
    //
    /**
     * Adds the four rules for bonds between an atom of the given environment and the ring carbon atoms of the link
     * types L13 to L16 (rule name suffixes a to d).
     *
     * @param aRuleList list to add the rules to
     * @param aNamePrefix prefix of the rule names
     * @param anEnvironment SMARTS of the atom bonded to the ring carbon atom
     * @param aLinkType link type of the atom bonded to the ring carbon atom if the ring carbon atom is the right atom
     *                  of the rule (the environment is the left one), or -1 if the ring carbon atom is the left atom
     * @param aLinkTypePerRingCarbon link types of the atom bonded to the ring carbon atom for L13 to L16 if the ring
     *                               carbon atom is the left atom of the rule, otherwise null
     */
    private static void addRingCarbonRules(List<BondRule> aRuleList, String aNamePrefix, String anEnvironment,
                                           int aLinkType, int[] aLinkTypePerRingCarbon) {
        for (int i = 0; i < 4; i++) {
            String tmpName = aNamePrefix + BRICSFragmenter.RING_CARBON_RULE_SUFFIXES[i];
            if (aLinkTypePerRingCarbon == null) {
                aRuleList.add(new BondRule(tmpName, anEnvironment, IBond.Order.SINGLE,
                        BRICSFragmenter.RING_CARBON_ENVIRONMENTS[i], aLinkType, BRICSFragmenter.RING_CARBON_LINK_TYPES[i]));
            } else {
                aRuleList.add(new BondRule(tmpName, BRICSFragmenter.RING_CARBON_ENVIRONMENTS[i], IBond.Order.SINGLE,
                        anEnvironment, BRICSFragmenter.RING_CARBON_LINK_TYPES[i], aLinkTypePerRingCarbon[i]));
            }
        }
    }
    //</editor-fold>
}
