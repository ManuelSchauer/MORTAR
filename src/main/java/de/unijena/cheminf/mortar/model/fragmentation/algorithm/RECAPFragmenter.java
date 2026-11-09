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
import java.util.Collections;
import java.util.Deque;
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
 * Implementation of the RECAP fragmentation ("Retrosynthetic Combinatorial Analysis Procedure") by Lewell et al., see
 * the references below. RECAP cleaves acyclic bonds that are formed by common chemical reactions, so that the
 * resulting fragments are building blocks for combinatorial chemistry. There are eleven bond types; the attachment
 * points of the fragments are marked with link atoms that carry the number of the bond type:
 * <ol>
 *     <li>amide: carbonyl carbon - nitrogen (no ureas)</li>
 *     <li>ester: carbonyl carbon - ester oxygen</li>
 *     <li>amine: sp3 carbon - amine nitrogen (not acylated, sulfonylated, or bonded to another heteroatom)</li>
 *     <li>urea: nitrogen - urea carbonyl carbon (both bonds)</li>
 *     <li>ether: sp3 carbon - ether oxygen</li>
 *     <li>olefin: carbon = carbon (the double bond is cleaved)</li>
 *     <li>quaternary nitrogen: quaternary ammonium nitrogen - carbon</li>
 *     <li>aromatic nitrogen - aliphatic carbon</li>
 *     <li>lactam nitrogen - aliphatic carbon</li>
 *     <li>aromatic carbon - aromatic carbon (biaryls)</li>
 *     <li>sulphonamide: nitrogen - sulfonyl sulfur</li>
 * </ol>
 * Ring bonds are never cleaved. If a bond matches more than one rule, the rule with the lowest number is applied.
 *
 * <p><b>Sources and reconstruction of the rules.</b> The publication (Lewell et al. 1998) is not freely available.
 * The rules are therefore implemented according to the presentation of the authors (Daylight EMUG 1997), which
 * shows the eleven bond types with the cleaved bonds, states that small fragments are avoided, and gives the
 * fragmentation of cisapride as example (see references). The SMARTS of the rules are formulated in this class; the
 * following details are derived from the presentation:
 * <ul>
 *     <li>All bonds of the central atom of a group can be cleaved (both bonds of ureas and dialkyl ethers, all bonds
 *     of tertiary amines and quaternary ammonium ions), so fragments that link two or more other fragments, e.g. a
 *     single nitrogen or carbonyl group, are retained.</li>
 *     <li>Ethers and amines are only cleaved at the bond to an sp3 carbon atom: in cisapride, the oxygen atom of the
 *     aryl ether remains at the fluorophenyl ring, and the methoxy groups are not cleaved off.</li>
 *     <li>Small terminal fragments are avoided: a bond is not cleaved if one of the two resulting parts is a terminal
 *     fragment (only one attachment point) with fewer heavy atoms than the minimum size of terminal fragments
 *     (setting, default {@value #MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT}). In cisapride, methyl, methoxy and amino
 *     groups are not cleaved off, and in the ranitidine example of the presentation, the dimethylamino group is
 *     kept. The default protects alkyl groups up to butyl. Linking fragments of any size, like the propyl linker of
 *     cisapride, are allowed.</li>
 *     <li>The attachment points are labelled with the number of the bond type, as in the cisapride example.</li>
 * </ul>
 * The details of the publication that the presentation does not show may differ from this reconstruction.
 *
 * <p><b>Algorithm.</b>
 * <ol>
 *     <li>Explicit hydrogen atoms are suppressed, and ring membership and aromaticity are perceived with the Daylight
 *     model, on which the SMARTS semantics of the rules are based.</li>
 *     <li>Every acyclic, non-aromatic single or double bond is tested against the rules in the order of their numbers;
 *     the atom environments are evaluated lazily and only once per atom.</li>
 *     <li>The cleavable bonds are processed in the order of their bond type number (and bond index); a bond is cleaved
 *     unless its cleavage, together with the bonds cleaved before, creates a terminal fragment below the minimum size.
 *     Since further cleavages only split fragments, no small terminal fragment is created. All accepted bonds are
 *     cleaved simultaneously, so every heavy atom of the molecule is contained in exactly one fragment.</li>
 *     <li>Each cleaved bond is replaced by two bonds of the same order to new link atoms. A link atom is a pseudo atom
 *     carrying the number of the bond type as mass number (resulting in SMILES like [5*]Oc1ccc(F)cc1, if the mass
 *     number is encoded) and as property {@link #BOND_TYPE_PROPERTY_KEY}. Alternatively, open valences can be
 *     saturated with hydrogen atoms.</li>
 *     <li>Stereo elements are transferred to the fragments; ligands that were cleaved off are replaced by the link atom
 *     (or hydrogen atom) at their position, and stereo elements that lost their stereogenicity are removed.</li>
 * </ol>
 * Pseudo atoms of the input, e.g. attachment points of a preceding fragmentation in a pipeline, count as attachment
 * points and not as heavy atoms in the size check.
 *
 * <p><b>Differences to the RECAP implementation of the RDKit.</b> The RDKit implementation removes the oxygen atom of
 * cleaved ethers, the nitrogen atom of cleaved amines and the carbonyl group of cleaved ureas, cleaves aryl ethers
 * and aryl amines, cleaves the hydroxy group off carboxylic acids, has no rule for quaternary nitrogen atoms, and only
 * forbids fragments that consist of one to three carbon atoms. Its function RecapDecompose applies the rules
 * hierarchically and returns all unique leaf fragments, which can overlap; the multiplicity of fragments is lost.
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
     * Default minimum number of heavy atoms of a terminal fragment (fragment with one attachment point); smaller
     * terminal fragments are not cleaved off. The default protects alkyl groups up to butyl.
     */
    public static final int MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT = 5;
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
    private static final List<BondRule> RULES = RECAPFragmenter.createRules();
    //</editor-fold>
    //
    //<editor-fold desc="Private records">
    /**
     * A RECAP bond rule: an acyclic bond of the given order between an atom matching the left environment and an atom
     * matching the right environment is cleavable; both atoms are labelled with the bond type.
     *
     * @param bondType number of the RECAP bond type, 1 to 11
     * @param name name of the bond type
     * @param leftEnvironment SMARTS of the left atom
     * @param order order of the bond, single or double
     * @param rightEnvironment SMARTS of the right atom
     */
    private record BondRule(int bondType, String name, String leftEnvironment, IBond.Order order, String rightEnvironment) {
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
     * Setting for the saturation of open valences.
     */
    private final SimpleIDisplayEnumConstantProperty saturationSetting;
    /**
     * Setting for returning molecules without cleavable bonds as one unchanged fragment.
     */
    private final SimpleBooleanProperty returnUnfragmentedMoleculesSetting;
    /**
     * Setting for the minimum number of heavy atoms of terminal fragments.
     */
    private final SimpleIntegerProperty minimumTerminalFragmentSizeSetting;
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
        Map<String, Integer> tmpEnvironmentIndexMap = new LinkedHashMap<>();
        this.compiledRules = new CompiledRule[RECAPFragmenter.RULES.size()];
        for (int i = 0; i < RECAPFragmenter.RULES.size(); i++) {
            BondRule tmpRule = RECAPFragmenter.RULES.get(i);
            int tmpLeftIndex = tmpEnvironmentIndexMap.computeIfAbsent(tmpRule.leftEnvironment(), aKey -> tmpEnvironmentIndexMap.size());
            int tmpRightIndex = tmpEnvironmentIndexMap.computeIfAbsent(tmpRule.rightEnvironment(), aKey -> tmpEnvironmentIndexMap.size());
            this.compiledRules[i] = new CompiledRule(tmpRule, tmpLeftIndex, tmpRightIndex);
        }
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
        int tmpNumberOfSettings = 3;
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
        this.minimumTerminalFragmentSizeSetting = new SimpleIntegerProperty(this,
                "Minimum size of terminal fragments setting",
                RECAPFragmenter.MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT) {
            @Override
            public void set(int newValue) {
                if (newValue > 0) {
                    super.set(newValue);
                } else {
                    IllegalArgumentException tmpException = new IllegalArgumentException(
                            "The minimum size of terminal fragments must be positive.");
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
        this.settingNameTooltipTextMap.put(this.minimumTerminalFragmentSizeSetting.getName(),
                Message.get("RECAPFragmenter.minimumTerminalFragmentSizeSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.minimumTerminalFragmentSizeSetting.getName(),
                Message.get("RECAPFragmenter.minimumTerminalFragmentSizeSetting.displayName"));
        this.returnUnfragmentedMoleculesSetting = new SimpleBooleanProperty(this,
                "Return unfragmented molecules setting",
                RECAPFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
        this.settingNameTooltipTextMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("RECAPFragmenter.returnUnfragmentedMoleculesSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.returnUnfragmentedMoleculesSetting.getName(),
                Message.get("RECAPFragmenter.returnUnfragmentedMoleculesSetting.displayName"));
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        this.settings.add(this.saturationSetting);
        this.settings.add(this.minimumTerminalFragmentSizeSetting);
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
     * Returns the minimum number of heavy atoms of terminal fragments; smaller terminal fragments are not cleaved off.
     *
     * @return minimum size of terminal fragments
     */
    public int getMinimumTerminalFragmentSizeSetting() {
        return this.minimumTerminalFragmentSizeSetting.get();
    }
    //
    /**
     * Returns the property object of the setting for the minimum size of terminal fragments.
     *
     * @return property object of the setting
     */
    public SimpleIntegerProperty minimumTerminalFragmentSizeSettingProperty() {
        return this.minimumTerminalFragmentSizeSetting;
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
     * Sets the minimum number of heavy atoms of terminal fragments; smaller terminal fragments are not cleaved off.
     *
     * @param aSize minimum size, must be positive
     * @throws IllegalArgumentException if the given size is not positive
     */
    public void setMinimumTerminalFragmentSizeSetting(int aSize) throws IllegalArgumentException {
        if (aSize <= 0) {
            throw new IllegalArgumentException("The minimum size of terminal fragments must be positive.");
        }
        this.minimumTerminalFragmentSizeSetting.set(aSize);
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
        tmpCopy.setMinimumTerminalFragmentSizeSetting(this.getMinimumTerminalFragmentSizeSetting());
        tmpCopy.setReturnUnfragmentedMoleculesSetting(this.getReturnUnfragmentedMoleculesSetting());
        return tmpCopy;
    }
    //
    @Override
    public void restoreDefaultSettings() {
        this.saturationSetting.set(RECAPFragmenter.SATURATION_OPTION_DEFAULT);
        this.minimumTerminalFragmentSizeSetting.set(RECAPFragmenter.MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT);
        this.returnUnfragmentedMoleculesSetting.set(RECAPFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT);
    }
    //
    /**
     * {@inheritDoc}
     *
     * <p>The fragments are the disconnected parts of a clone of the given molecule, in which all selected RECAP bonds
     * have been cleaved and the open valences saturated according to the saturation setting. If no bond is cleaved,
     * the clone (or its disconnected parts, e.g. counter ions, separately) is returned or, depending on the respective
     * setting, no fragment is returned.</p>
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
            List<RECAPBond> tmpBondsToCleave = this.selectBondsToCleave(tmpMolecule, this.findRECAPBonds(tmpMolecule));
            if (!tmpBondsToCleave.isEmpty()) {
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
     * Prepares the given molecule (in place) for the matching of the RECAP rules: missing implicit hydrogen counts are
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
     * Identifies all bonds of the given prepared molecule that match a RECAP rule. If a bond matches several rules,
     * the rule with the lowest bond type number is assigned.
     *
     * @param aPreparedMolecule molecule prepared with {@link #prepareMolecule(IAtomContainer)}
     * @return the cleavable bonds with their bond types, in the order of the bonds in the molecule
     */
    List<RECAPBond> findRECAPBonds(IAtomContainer aPreparedMolecule) {
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
            for (CompiledRule tmpCompiledRule : this.compiledRules) {
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
     * Selects the bonds to cleave from the given cleavable bonds, so that no small terminal fragment is created. The
     * bonds are processed in the order of their bond type (and bond index); a bond is selected unless one of the two
     * parts that result from its cleavage, together with the bonds selected before, is a terminal fragment (one
     * attachment point) with fewer heavy atoms than the minimum size of terminal fragments.
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
        tmpCandidates.sort((aFirst, aSecond) -> aFirst.bondType() != aSecond.bondType()
                ? Integer.compare(aFirst.bondType(), aSecond.bondType())
                : Integer.compare(aFirst.bond().getIndex(), aSecond.bond().getIndex()));
        GraphUtil.EdgeToBondMap tmpEdgeToBondMap = GraphUtil.EdgeToBondMap.withSpaceFor(aPreparedMolecule);
        int[][] tmpAdjacencyList = GraphUtil.toAdjList(aPreparedMolecule, tmpEdgeToBondMap);
        boolean[] tmpIsCleaved = new boolean[aPreparedMolecule.getBondCount()];
        int tmpMinimumSize = this.minimumTerminalFragmentSizeSetting.get();
        List<RECAPBond> tmpSelectedBonds = new ArrayList<>(tmpCandidates.size());
        for (RECAPBond tmpCandidate : tmpCandidates) {
            int tmpBondIndex = tmpCandidate.bond().getIndex();
            tmpIsCleaved[tmpBondIndex] = true;
            if (RECAPFragmenter.isSmallTerminalPart(aPreparedMolecule, tmpCandidate.bond().getBegin().getIndex(),
                    tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved, tmpMinimumSize)
                    || RECAPFragmenter.isSmallTerminalPart(aPreparedMolecule, tmpCandidate.bond().getEnd().getIndex(),
                    tmpAdjacencyList, tmpEdgeToBondMap, tmpIsCleaved, tmpMinimumSize)) {
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
     * Determines the part of the molecule that contains the given atom when the marked bonds are cleaved, and returns
     * whether it is a small terminal fragment: it has exactly one attachment point (cleaved bond or pseudo atom of the
     * input) and fewer heavy atoms (pseudo atoms not counted) than the given minimum size.
     *
     * @param aMolecule the prepared molecule
     * @param aStartAtomIndex index of an atom of the part
     * @param anAdjacencyList adjacency list of the molecule
     * @param anEdgeToBondMap map of the edges of the adjacency list to the bonds
     * @param anIsCleavedArray cleaved bonds by bond index
     * @param aMinimumSize minimum number of heavy atoms of a terminal fragment
     * @return true if the part is a terminal fragment below the minimum size
     */
    private static boolean isSmallTerminalPart(IAtomContainer aMolecule, int aStartAtomIndex, int[][] anAdjacencyList,
                                               GraphUtil.EdgeToBondMap anEdgeToBondMap, boolean[] anIsCleavedArray,
                                               int aMinimumSize) {
        boolean[] tmpIsVisited = new boolean[anAdjacencyList.length];
        Deque<Integer> tmpQueue = new ArrayDeque<>();
        tmpQueue.add(aStartAtomIndex);
        tmpIsVisited[aStartAtomIndex] = true;
        int tmpHeavyAtoms = 0;
        int tmpAttachmentPoints = 0;
        while (!tmpQueue.isEmpty()) {
            int tmpAtomIndex = tmpQueue.poll();
            IAtom tmpAtom = aMolecule.getAtom(tmpAtomIndex);
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAttachmentPoints++;
            } else if (tmpAtom.getAtomicNumber() == null || tmpAtom.getAtomicNumber() != 1) {
                tmpHeavyAtoms++;
            }
            for (int tmpNeighbour : anAdjacencyList[tmpAtomIndex]) {
                if (anIsCleavedArray[anEdgeToBondMap.get(tmpAtomIndex, tmpNeighbour).getIndex()]) {
                    tmpAttachmentPoints++;
                } else if (!tmpIsVisited[tmpNeighbour]) {
                    tmpIsVisited[tmpNeighbour] = true;
                    tmpQueue.add(tmpNeighbour);
                }
            }
            //the part is neither terminal nor small anymore, no further search needed
            if (tmpAttachmentPoints > 1 || tmpHeavyAtoms >= aMinimumSize) {
                return false;
            }
        }
        return tmpAttachmentPoints == 1 && tmpHeavyAtoms < aMinimumSize;
    }
    //
    /**
     * Cleaves the given bonds of the prepared molecule simultaneously, saturates the open valences according to the
     * saturation setting, transfers the stereo elements, and returns the disconnected parts as fragments.
     *
     * @param aMolecule the prepared molecule, it is altered
     * @param aRECAPBondList the bonds to cleave
     * @return list of fragments
     */
    @SuppressWarnings("rawtypes")
    private List<IAtomContainer> cleaveBonds(IAtomContainer aMolecule, List<RECAPBond> aRECAPBondList) {
        boolean tmpIsHydrogenSaturation = this.saturationSetting.get().equals(SaturationOption.HYDROGEN_SATURATION);
        //the stereo elements are detached before the bonds are removed, because the CDK atom container would otherwise
        //adapt them on bond removal, e.g. by replacing a cleaved ligand of a stereocentre with an implicit hydrogen atom
        List<IStereoElement> tmpOriginalStereoElements = new ArrayList<>();
        for (IStereoElement tmpElement : aMolecule.stereoElements()) {
            tmpOriginalStereoElements.add(tmpElement);
        }
        Set<IAtom> tmpNonStereogenicAtomsOfMolecule = tmpOriginalStereoElements.isEmpty()
                ? Collections.emptySet() : RECAPFragmenter.findNonStereogenicAtoms(aMolecule);
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
            RECAPFragmenter.removeStereoElementsMadeNonStereogenic(tmpFragment, tmpNonStereogenicAtomsOfMolecule);
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
     * Returns the atoms of the given molecule that are not stereogenic according to the CDK class
     * {@link Stereocenters}, which perceives stereocentres by constitutional symmetry. Mass numbers of pseudo atoms,
     * i.e. the bond types, are ignored.
     *
     * @param aMolecule molecule
     * @return set of non-stereogenic atoms, or null if the stereocentres cannot be determined
     */
    private static Set<IAtom> findNonStereogenicAtoms(IAtomContainer aMolecule) {
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
            Set<IAtom> tmpNonStereogenicAtoms = new HashSet<>();
            for (int i = 0; i < aMolecule.getAtomCount(); i++) {
                if (!tmpStereocenters.isStereocenter(i)) {
                    tmpNonStereogenicAtoms.add(aMolecule.getAtom(i));
                }
            }
            return tmpNonStereogenicAtoms;
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
     * Creates the eleven RECAP bond rules in the order of their numbers.
     *
     * @return list of all rules
     */
    private static List<BondRule> createRules() {
        IBond.Order tmpSingle = IBond.Order.SINGLE;
        //carbonyl carbon atom of amides, not of ureas (bond type 4)
        String tmpAmideCarbon = "[C;$(C=O);!$(C(~[#7])~[#7])]";
        //nitrogen atom of amides and ureas (also aromatic nitrogen atoms)
        String tmpAmideNitrogen = "[#7;+0]";
        //carbonyl carbon atom of esters (also carbonates and carbamates)
        String tmpEsterCarbon = "[C;$(C=O)]";
        //ester oxygen atom, bonded to two carbon atoms (no carboxylic acids)
        String tmpEsterOxygen = "[O;D2;$(O([#6])[#6])]";
        //sp3 carbon atom of amines and ethers
        String tmpSp3Carbon = "[C;X4]";
        //amine nitrogen atom: neutral, three-connected, not acylated or bonded to an imine or thiocarbonyl carbon atom
        //(amides, ureas, amidines), and not bonded to another heteroatom (e.g. sulfonamides, hydrazines)
        String tmpAmineNitrogen = "[N;+0;X3;!$(N[#6]=[#7,#8,#16]);!$(N[#7,#8,#15,#16])]";
        //carbonyl carbon atom of ureas
        String tmpUreaCarbon = "[C;$(C(=O)(~[#7])~[#7])]";
        //ether oxygen atom, bonded to two carbon atoms, not part of an ester
        String tmpEtherOxygen = "[O;D2;+0;$(O([#6])[#6]);!$(O[#6]=[#7,#8,#16])]";
        //olefinic carbon atom, not part of a cumulene
        String tmpOlefinCarbon = "[C;!$(C(=*)=*)]";
        //quaternary ammonium nitrogen atom with four carbon substituents
        String tmpQuaternaryNitrogen = "[N+;D4;!$([N+]~[!#6])]";
        //aliphatic (non-aromatic) carbon atom
        String tmpAliphaticCarbon = "[C]";
        //lactam nitrogen atom: ring nitrogen atom bonded to a ring carbonyl carbon atom
        String tmpLactamNitrogen = "[N;R;$(N@C=O)]";
        //sulfonyl sulfur atom
        String tmpSulfonylSulfur = "[S;X4;$(S(=O)=O)]";
        List<BondRule> tmpRules = new ArrayList<>(11);
        tmpRules.add(new BondRule(1, "amide", tmpAmideCarbon, tmpSingle, tmpAmideNitrogen));
        tmpRules.add(new BondRule(2, "ester", tmpEsterCarbon, tmpSingle, tmpEsterOxygen));
        tmpRules.add(new BondRule(3, "amine", tmpSp3Carbon, tmpSingle, tmpAmineNitrogen));
        tmpRules.add(new BondRule(4, "urea", tmpAmideNitrogen, tmpSingle, tmpUreaCarbon));
        tmpRules.add(new BondRule(5, "ether", tmpSp3Carbon, tmpSingle, tmpEtherOxygen));
        tmpRules.add(new BondRule(6, "olefin", tmpOlefinCarbon, IBond.Order.DOUBLE, tmpOlefinCarbon));
        tmpRules.add(new BondRule(7, "quaternary nitrogen", tmpQuaternaryNitrogen, tmpSingle, "[#6]"));
        tmpRules.add(new BondRule(8, "aromatic nitrogen - aliphatic carbon", "[n]", tmpSingle, tmpAliphaticCarbon));
        tmpRules.add(new BondRule(9, "lactam nitrogen - aliphatic carbon", tmpLactamNitrogen, tmpSingle, tmpAliphaticCarbon));
        tmpRules.add(new BondRule(10, "aromatic carbon - aromatic carbon", "[c]", tmpSingle, "[c]"));
        tmpRules.add(new BondRule(11, "sulphonamide", tmpAmideNitrogen, tmpSingle, tmpSulfonylSulfur));
        return Collections.unmodifiableList(tmpRules);
    }
    //</editor-fold>
}
