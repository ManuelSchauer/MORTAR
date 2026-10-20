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
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Fragmenter for pure hydrocarbon (alkyl) structures, intended as downstream step after other fragmentations, e.g.
 * of the alkane fragments of the {@link ErtlFunctionalGroupsFinderFragmenter}.
 * <p>
 * The atoms of a molecule are classified into substructures and every bond is either kept or cut according to the
 * following rules; the fragments are the connected parts of the kept bonds. The examples use the default settings
 * (isolation of tertiary/quaternary carbons active, separation from rings inactive, maximum chain length 6, hydrogen
 * saturation) unless stated otherwise.
 * </p>
 * <ol>
 *     <li><b>Rings</b>: ring bonds are never cut, fused and spiro ring systems stay one fragment; rings connected by a
 *     non-conjugated single bond are separated. Example: {@code C1CCC(CC1)C1CCCCC1 -> C1CCCCC1, C1CCCCC1}.</li>
 *     <li><b>Conjugated pi systems</b>: bonds of conjugated systems are never cut. A single bond is conjugated if both
 *     its atoms are part of a double bond or are aromatic (double-single-double), or both are part of a triple bond
 *     (triple-single-triple); a double and a triple bond separated by a single bond are not conjugated. Aromatic atoms
 *     are always part of a conjugated system, also if they are remnants of aromatic rings without aromatic bonds or
 *     bond orders (e.g. after the removal of hetero atoms by a previous fragmentation); bonds between aromatic atoms
 *     are never cut and aromatic atoms are neither tertiary/quaternary carbons nor chain atoms. Rings and conjugated
 *     systems form the "core" substructures. Examples: {@code C=Cc1ccccc1 -> C=Cc1ccccc1};
 *     {@code C=CC#C -> C=C, C#C}; Ertl remnant {@code cCCC -> c, CCC}.</li>
 *     <li><b>Multiple bonds</b>: double and triple bonds are never cut. A multiple bond attached to a core atom
 *     (e.g. an exocyclic double bond) belongs to that core, all other multiple bonds form isolated fragments.
 *     Examples: {@code C=C1CCCCC1 -> C=C1CCCCC1}; {@code CC=CC -> C=C, C, C}.</li>
 *     <li><b>Tertiary and quaternary carbons</b>: non-cyclic, non-aromatic carbons with only single bonds to three or
 *     four other atoms.
 *     <ul>
 *         <li>With the isolation setting, such a carbon is cut from all its neighbors, unless it is bonded to a ring
 *         and not to be separated from it. The neighbors of an isolated carbon are treated like any other atoms, e.g.
 *         saturated neighbors are ordinary chain atoms. Example: {@code CC(C)(C)CCCC(C)C -> *C(*)(*)*, *C(*)*, C, C,
 *         C, C, C, CCC}; with separation from rings: {@code c1ccc(cc1)C(C)CC -> c1ccccc1, *C(*)*, C, CC}.</li>
 *         <li>A non-isolated tertiary/quaternary carbon (isolation setting inactive, or kept at a ring) forms a unit
 *         with its saturated neighbors (like tert-butyl) and keeps its bonds to other non-isolated
 *         tertiary/quaternary carbons, to conjugated systems, and to rings unless separated from them. Bonds of the
 *         saturated neighbors to further atoms are cut. Examples: isolation inactive
 *         {@code CC(C)(C)CCCC(C)C -> CC(C)(C)C, CC(C)C, C}; kept at ring {@code CC(C)c1ccccc1 -> CC(C)c1ccccc1};
 *         bonded to ring and conjugated system {@code CC(C1CCCCC1)C=CC=C -> C=CC=CC(C)C1CCCCC1}.</li>
 *     </ul></li>
 *     <li><b>Pseudo atoms</b>: every cut bond of a tertiary/quaternary carbon is replaced by a pseudo atom on that carbon,
 *     so that the carbon keeps its original number of bonds; the other side of the cut is saturated with hydrogen.
 *     Examples: isolation inactive {@code C=CC(C)(C)C -> C=C, *C(C)(C)C}; separation from rings and isolation
 *     inactive {@code CC(C)c1ccccc1 -> *C(C)C, c1ccccc1}; kept at ring next to an isolated carbon
 *     {@code CC(C)C(C)C1CCCCC1 -> *C(*)*, *C(C)C1CCCCC1, C, C}.</li>
 *     <li><b>Linear chains</b>: the remaining saturated carbons form chains that are cut from all other substructures and,
 *     if activated, dissected into chains of the set maximum length, walked from the chain end with the lowest atom
 *     index. Examples: {@code CC1CCCCC1 -> C1CCCCC1, C}; maximum length 3:
 *     {@code CCCCCCCC -> CCC, CCC, CC}.</li>
 * </ol>
 * <p>
 * Afterwards, cut bonds are saturated with implicit hydrogen atoms (if the respective setting is active) and the
 * stereo chemistry of the input molecule is transferred to the fragments where the stereo element is still complete
 * and still a stereocenter (if the respective setting is active).
 * </p>
 * <p>
 * Every atom of a returned fragment (except pseudo atoms) carries the index of the corresponding atom in the
 * fragmented (hydrogen-suppressed) molecule as property {@link #INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY}.
 * </p>
 *
 * @author Maximilian Rottmann (maximilian.rottmann@studmail.w-hs.de), Manuel Schauer
 * @version 2.0.0.0
 */
public class AlkylStructureFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the fragmenter.
     */
    public static final String ALGORITHM_NAME = "Alkyl Fragmenter";
    /**
     * Default value of the setting for keeping non-fragmentable molecules in the fragmentation pipeline. Active by
     * default, since the fragmenter is mainly applied as downstream step in a pipeline (e.g. after the
     * {@link ErtlFunctionalGroupsFinderFragmenter}), where fragments of the previous steps would be lost otherwise.
     */
    public static final boolean KEEP_NON_FRAGMENTABLE_MOLECULES_SETTING_DEFAULT = true;
    /**
     * Default value of the maximum length of linear chain fragments.
     */
    public static final int MAX_CHAIN_LENGTH_SETTING_DEFAULT = 6;
    /**
     * Default value of the setting for dissecting linear chains according to the maximum chain length.
     */
    public static final boolean FRAGMENT_SIDE_CHAINS_SETTING_DEFAULT = true;
    /**
     * Default value of the setting for isolating tertiary and quaternary carbons.
     */
    public static final boolean ISOLATE_TERT_QUAT_CARBONS_SETTING_DEFAULT = true;
    /**
     * Default value of the setting for separating tertiary and quaternary carbons from rings.
     */
    public static final boolean SEPARATE_TERT_QUAT_CARBON_FROM_RING_SETTING_DEFAULT = false;
    /**
     * Default value of the setting for transferring the stereo chemistry of the input molecule to the fragments.
     */
    public static final boolean PRESERVE_STEREO_CHEMISTRY_SETTING_DEFAULT = true;
    /**
     * Key of the atom property holding the index of the corresponding atom in the fragmented molecule.
     */
    public static final String INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY = "ASF.ATOM_INDEX";
    /**
     * Key of the bond property holding the index of the corresponding bond in the fragmented molecule.
     */
    public static final String INTERNAL_ASF_BOND_INDEX_PROPERTY_KEY = "ASF.BOND_INDEX";
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    //note: settings properties are not commented here to avoid "double comments" with their get methods in javadoc
    private final SimpleIDisplayEnumConstantProperty fragmentSaturationSetting;

    private final SimpleBooleanProperty keepNonFragmentableMoleculesSetting;

    private final SimpleBooleanProperty fragmentSideChainsSetting;

    private final SimpleIntegerProperty maxChainLengthSetting;

    private final SimpleBooleanProperty isolateTertQuatCarbonSetting;

    private final SimpleBooleanProperty separateTertQuatCarbonFromRingSetting;

    private final SimpleBooleanProperty preserveStereoChemistrySetting;

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
     * Builder for the fragment atoms and bonds.
     */
    private final IChemObjectBuilder builder;

    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(AlkylStructureFragmenter.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, all settings are initialized with their default values as declared in the respective public constants.
     */
    public AlkylStructureFragmenter() {
        int tmpNumberOfSettings = 7;
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        int tmpInitialCapacity = CollectionUtil.calculateInitialHashCollectionCapacity(tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToTooltipTextMap = new HashMap<>(tmpInitialCapacity, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameToDisplayNameMap = new HashMap<>(tmpInitialCapacity, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.fragmentSaturationSetting = new SimpleIDisplayEnumConstantProperty(this, "Fragment saturation setting",
                IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT, IMoleculeFragmenter.FragmentSaturationOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    AlkylStructureFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.addSetting(this.fragmentSaturationSetting, "fragmentSaturationSetting");
        this.keepNonFragmentableMoleculesSetting = new SimpleBooleanProperty(this, "Keep non-fragmentable molecules in Pipeline",
                AlkylStructureFragmenter.KEEP_NON_FRAGMENTABLE_MOLECULES_SETTING_DEFAULT);
        this.addSetting(this.keepNonFragmentableMoleculesSetting, "keepNonFragmentableMoleculesSetting");
        this.fragmentSideChainsSetting = new SimpleBooleanProperty(this, "Fragment linear chains",
                AlkylStructureFragmenter.FRAGMENT_SIDE_CHAINS_SETTING_DEFAULT);
        this.addSetting(this.fragmentSideChainsSetting, "fragmentSideChainsSetting");
        this.maxChainLengthSetting = new SimpleIntegerProperty(this, "Limit length of returned chains",
                AlkylStructureFragmenter.MAX_CHAIN_LENGTH_SETTING_DEFAULT) {
            @Override
            public void set(int newValue) throws IllegalArgumentException {
                if (newValue < 1) {
                    IllegalArgumentException tmpException = new IllegalArgumentException(
                            "Maximum length of alkyl chain fragments setting only accepts positive, non-zero values, got: " + newValue);
                    AlkylStructureFragmenter.LOGGER.log(Level.WARNING, tmpException.toString(), tmpException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            tmpException.toString(),
                            tmpException);
                    //throws the exception to properly reset the binding
                    throw tmpException;
                }
                super.set(newValue);
            }
        };
        this.addSetting(this.maxChainLengthSetting, "maxChainLengthSetting");
        this.isolateTertQuatCarbonSetting = new SimpleBooleanProperty(this, "Isolate non-cyclic tertiary and quaternary carbons",
                AlkylStructureFragmenter.ISOLATE_TERT_QUAT_CARBONS_SETTING_DEFAULT);
        this.addSetting(this.isolateTertQuatCarbonSetting, "isolateTertQuatCarbonsSetting");
        this.separateTertQuatCarbonFromRingSetting = new SimpleBooleanProperty(this, "Separate tertiary and quaternary carbon atoms from rings",
                AlkylStructureFragmenter.SEPARATE_TERT_QUAT_CARBON_FROM_RING_SETTING_DEFAULT);
        this.addSetting(this.separateTertQuatCarbonFromRingSetting, "separateTertQuatCarbonFromRingSetting");
        this.preserveStereoChemistrySetting = new SimpleBooleanProperty(this, "Preserve stereo chemistry in fragmentation",
                AlkylStructureFragmenter.PRESERVE_STEREO_CHEMISTRY_SETTING_DEFAULT);
        this.addSetting(this.preserveStereoChemistrySetting, "preserveStereoChemistrySetting");
        this.builder = SilentChemObjectBuilder.getInstance();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
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
        return AlkylStructureFragmenter.ALGORITHM_NAME;
    }

    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("AlkylStructureFragmenter.displayName");
    }

    /**
     * Returns the currently set option for saturating free valences of the returned fragments.
     *
     * @return the set option
     */
    public IMoleculeFragmenter.FragmentSaturationOption getFragmentSaturationSetting() {
        return (IMoleculeFragmenter.FragmentSaturationOption) this.fragmentSaturationSetting.get();
    }

    /**
     * Returns the property of the setting for saturating free valences of the returned fragments.
     *
     * @return setting property for fragment saturation
     */
    public SimpleIDisplayEnumConstantProperty fragmentSaturationSettingProperty() {
        return this.fragmentSaturationSetting;
    }

    /**
     * Returns whether non-fragmentable molecules are kept in the fragmentation pipeline.
     *
     * @return value of the setting
     */
    public boolean getKeepNonFragmentableMoleculesSetting() {
        return this.keepNonFragmentableMoleculesSetting.get();
    }

    /**
     * Returns the property of the setting for keeping non-fragmentable molecules in the fragmentation pipeline.
     *
     * @return setting property
     */
    public SimpleBooleanProperty getKeepNonFragmentableMoleculesSettingProperty() {
        return this.keepNonFragmentableMoleculesSetting;
    }

    /**
     * Returns whether linear chains are dissected according to the maximum chain length.
     *
     * @return value of the setting
     */
    public boolean getFragmentSideChainsSetting() {
        return this.fragmentSideChainsSetting.get();
    }

    /**
     * Returns the property of the setting for dissecting linear chains.
     *
     * @return setting property
     */
    public SimpleBooleanProperty getFragmentSideChainsSettingProperty() {
        return this.fragmentSideChainsSetting;
    }

    /**
     * Returns the maximum length of linear chain fragments.
     *
     * @return value of the setting
     */
    public int getMaxChainLengthSetting() {
        return this.maxChainLengthSetting.get();
    }

    /**
     * Returns the property of the setting for the maximum length of linear chain fragments.
     *
     * @return setting property
     */
    public SimpleIntegerProperty getMaxChainLengthSettingProperty() {
        return this.maxChainLengthSetting;
    }

    /**
     * Returns whether non-cyclic tertiary and quaternary carbons are isolated.
     *
     * @return value of the setting
     */
    public boolean getIsolateTertQuatCarbonSetting() {
        return this.isolateTertQuatCarbonSetting.get();
    }

    /**
     * Returns the property of the setting for isolating non-cyclic tertiary and quaternary carbons.
     *
     * @return setting property
     */
    public SimpleBooleanProperty getIsolateTertQuatCarbonsSettingProperty() {
        return this.isolateTertQuatCarbonSetting;
    }

    /**
     * Returns whether tertiary and quaternary carbons are separated from rings.
     *
     * @return value of the setting
     */
    public boolean getSeparateTertQuatCarbonFromRingSetting() {
        return this.separateTertQuatCarbonFromRingSetting.get();
    }

    /**
     * Returns the property of the setting for separating tertiary and quaternary carbons from rings.
     *
     * @return setting property
     */
    public SimpleBooleanProperty getSeparateTertQuatCarbonFromRingSettingProperty() {
        return this.separateTertQuatCarbonFromRingSetting;
    }

    /**
     * Returns whether the stereo chemistry of the input molecule is transferred to the fragments.
     *
     * @return value of the setting
     */
    public boolean getPreserveStereoChemistrySetting() {
        return this.preserveStereoChemistrySetting.get();
    }

    /**
     * Returns the property of the setting for transferring the stereo chemistry to the fragments.
     *
     * @return setting property
     */
    public SimpleBooleanProperty getPreserveStereoChemistrySettingProperty() {
        return this.preserveStereoChemistrySetting;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set">
    /**
     * Sets the option for saturating free valences of the returned fragments.
     *
     * @param anOption the saturation option to use
     * @throws NullPointerException if the given option is null
     */
    public void setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given saturation option is null.");
        this.fragmentSaturationSetting.set(anOption);
    }

    /**
     * Sets whether non-fragmentable molecules are kept in the fragmentation pipeline.
     *
     * @param aBoolean new value
     */
    public void setKeepNonFragmentableMoleculesSetting(boolean aBoolean) {
        this.keepNonFragmentableMoleculesSetting.set(aBoolean);
    }

    /**
     * Sets whether linear chains are dissected according to the maximum chain length.
     *
     * @param aBoolean new value
     */
    public void setFragmentSideChainsSetting(boolean aBoolean) {
        this.fragmentSideChainsSetting.set(aBoolean);
    }

    /**
     * Sets the maximum length of linear chain fragments.
     *
     * @param aValue new value, has to be at least 1
     * @throws IllegalArgumentException if the given value is smaller than 1
     */
    public void setMaxChainLengthSetting(int aValue) throws IllegalArgumentException {
        if (aValue < 1) {
            throw new IllegalArgumentException("Given chain length cannot be smaller than 1.");
        }
        this.maxChainLengthSetting.set(aValue);
    }

    /**
     * Sets whether non-cyclic tertiary and quaternary carbons are isolated.
     *
     * @param aBoolean new value
     */
    public void setIsolateTertQuatCarbonsSetting(boolean aBoolean) {
        this.isolateTertQuatCarbonSetting.set(aBoolean);
    }

    /**
     * Sets whether tertiary and quaternary carbons are separated from rings.
     *
     * @param aBoolean new value
     */
    public void setSeparateTertQuatCarbonFromRingSetting(boolean aBoolean) {
        this.separateTertQuatCarbonFromRingSetting.set(aBoolean);
    }

    /**
     * Sets whether the stereo chemistry of the input molecule is transferred to the fragments.
     *
     * @param aBoolean new value
     */
    public void setPreserveStereoChemistrySetting(boolean aBoolean) {
        this.preserveStereoChemistrySetting.set(aBoolean);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    @Override
    public IMoleculeFragmenter copy() {
        AlkylStructureFragmenter tmpCopy = new AlkylStructureFragmenter();
        tmpCopy.setFragmentSaturationSetting(this.getFragmentSaturationSetting());
        tmpCopy.setKeepNonFragmentableMoleculesSetting(this.getKeepNonFragmentableMoleculesSetting());
        tmpCopy.setFragmentSideChainsSetting(this.getFragmentSideChainsSetting());
        tmpCopy.setMaxChainLengthSetting(this.getMaxChainLengthSetting());
        tmpCopy.setIsolateTertQuatCarbonsSetting(this.getIsolateTertQuatCarbonSetting());
        tmpCopy.setSeparateTertQuatCarbonFromRingSetting(this.getSeparateTertQuatCarbonFromRingSetting());
        tmpCopy.setPreserveStereoChemistrySetting(this.getPreserveStereoChemistrySetting());
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {
        this.fragmentSaturationSetting.set(IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT);
        this.keepNonFragmentableMoleculesSetting.set(AlkylStructureFragmenter.KEEP_NON_FRAGMENTABLE_MOLECULES_SETTING_DEFAULT);
        this.fragmentSideChainsSetting.set(AlkylStructureFragmenter.FRAGMENT_SIDE_CHAINS_SETTING_DEFAULT);
        this.maxChainLengthSetting.set(AlkylStructureFragmenter.MAX_CHAIN_LENGTH_SETTING_DEFAULT);
        this.isolateTertQuatCarbonSetting.set(AlkylStructureFragmenter.ISOLATE_TERT_QUAT_CARBONS_SETTING_DEFAULT);
        this.separateTertQuatCarbonFromRingSetting.set(AlkylStructureFragmenter.SEPARATE_TERT_QUAT_CARBON_FROM_RING_SETTING_DEFAULT);
        this.preserveStereoChemistrySetting.set(AlkylStructureFragmenter.PRESERVE_STEREO_CHEMISTRY_SETTING_DEFAULT);
    }

    /**
     * {@inheritDoc}
     * <p>
     * A molecule is fragmentable if it contains at least one carbon atom, only carbon and hydrogen atoms (no hetero or
     * pseudo atoms) and only single, double, triple or quadruple bonds, or aromatic bonds without assigned bond order.
     * Non-fragmentable molecules are not filtered if they are to be kept in the pipeline (see
     * {@link #setKeepNonFragmentableMoleculesSetting(boolean)}); they are then returned unchanged as single fragment.
     * Null or empty molecules and molecules containing null atoms are always filtered.
     * </p>
     */
    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        if (Objects.isNull(aMolecule) || aMolecule.isEmpty()) {
            return true;
        }
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom)) {
                return true;
            }
        }
        if (this.isFragmentable(aMolecule)) {
            return false;
        }
        return !this.keepNonFragmentableMoleculesSetting.get();
    }

    /**
     * {@inheritDoc}
     * Always returns false, no preprocessing is needed.
     */
    @Override
    public boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return false;
    }

    @Override
    public boolean canBeFragmented(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return !this.shouldBeFiltered(aMolecule) && this.isFragmentable(aMolecule);
    }

    /**
     * {@inheritDoc}
     * No preprocessing is needed, a clone of the given molecule is returned.
     */
    @Override
    public IAtomContainer applyPreprocessing(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (this.shouldBeFiltered(aMolecule)) {
            throw new IllegalArgumentException("The given molecule cannot be preprocessed but should be filtered.");
        }
        return aMolecule.clone();
    }

    /**
     * {@inheritDoc}
     * <p>
     * See the class documentation for the fragmentation rules. A non-fragmentable molecule (see
     * {@link #shouldBeFiltered(IAtomContainer)}) is returned unchanged as single fragment if non-fragmentable molecules
     * are kept in the pipeline.
     * </p>
     */
    @Override
    public List<IAtomContainer> fragmentMolecule(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (!this.isFragmentable(aMolecule)) {
            if (this.keepNonFragmentableMoleculesSetting.get() && !this.shouldBeFiltered(aMolecule)) {
                List<IAtomContainer> tmpUnchangedMoleculeList = new ArrayList<>(1);
                tmpUnchangedMoleculeList.add(aMolecule);
                return tmpUnchangedMoleculeList;
            }
            throw new IllegalArgumentException("The given molecule cannot be fragmented but should be filtered.");
        }
        //explicit hydrogens are made implicit, otherwise they would count as neighbors of tertiary/quaternary carbons
        IAtomContainer tmpMolecule = AtomContainerManipulator.suppressHydrogens(aMolecule.clone());
        this.ensureImplicitHydrogenCounts(tmpMolecule);
        MoleculeClassification tmpClassification = this.classify(tmpMolecule);
        boolean[] tmpIsBondKept = this.determineKeptBonds(tmpMolecule, tmpClassification);
        if (this.fragmentSideChainsSetting.get()) {
            this.dissectLinearChains(tmpMolecule, tmpClassification, tmpIsBondKept, this.maxChainLengthSetting.get());
        }
        return this.buildFragments(tmpMolecule, tmpClassification, tmpIsBondKept);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Protected methods">
    /**
     * Classifies the atoms and bonds of the given (hydrogen-suppressed) molecule into the substructures relevant for the
     * fragmentation rules (see class documentation).
     *
     * @param aMolecule molecule to classify
     * @return the classification
     */
    protected MoleculeClassification classify(IAtomContainer aMolecule) {
        int tmpAtomCount = aMolecule.getAtomCount();
        int tmpBondCount = aMolecule.getBondCount();
        MoleculeClassification tmpClassification = new MoleculeClassification(tmpAtomCount, tmpBondCount);
        Cycles.markRingAtomsAndBonds(aMolecule);
        for (int i = 0; i < tmpBondCount; i++) {
            IBond tmpBond = aMolecule.getBond(i);
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            tmpClassification.isRingBond[i] = tmpBond.isInRing();
            if (this.isAromaticBondWithoutOrder(tmpBond)) {
                tmpClassification.hasAromaticBondWithoutOrder[tmpBegin] = true;
                tmpClassification.hasAromaticBondWithoutOrder[tmpEnd] = true;
            } else if (tmpBond.getOrder() == IBond.Order.DOUBLE) {
                tmpClassification.hasDoubleBond[tmpBegin] = true;
                tmpClassification.hasDoubleBond[tmpEnd] = true;
            } else if (tmpBond.getOrder() == IBond.Order.TRIPLE || tmpBond.getOrder() == IBond.Order.QUADRUPLE) {
                tmpClassification.hasTripleBond[tmpBegin] = true;
                tmpClassification.hasTripleBond[tmpEnd] = true;
            }
        }
        for (int i = 0; i < tmpAtomCount; i++) {
            IAtom tmpAtom = aMolecule.getAtom(i);
            tmpClassification.isRingAtom[i] = tmpAtom.isInRing();
            tmpClassification.isAromaticAtom[i] = tmpAtom.isAromatic();
            tmpClassification.isConjugatedAtom[i] = tmpAtom.isAromatic();
        }
        //conjugated bonds: aromatic bonds and single bonds between two pi atoms (D-S-D incl. aromatic, T-S-T)
        for (int i = 0; i < tmpBondCount; i++) {
            IBond tmpBond = aMolecule.getBond(i);
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            boolean tmpIsAromatic = this.isAromaticBondWithoutOrder(tmpBond)
                    || (tmpClassification.isAromaticAtom[tmpBegin] && tmpClassification.isAromaticAtom[tmpEnd]);
            boolean tmpIsConjugatedSingleBond = tmpBond.getOrder() == IBond.Order.SINGLE
                    && ((this.isDoubleBondOrAromaticAtom(tmpClassification, tmpBegin) && this.isDoubleBondOrAromaticAtom(tmpClassification, tmpEnd))
                    || (tmpClassification.hasTripleBond[tmpBegin] && tmpClassification.hasTripleBond[tmpEnd]));
            if (tmpIsAromatic || tmpIsConjugatedSingleBond) {
                tmpClassification.isConjugatedBond[i] = true;
                tmpClassification.isConjugatedAtom[tmpBegin] = true;
                tmpClassification.isConjugatedAtom[tmpEnd] = true;
            }
        }
        //multiple bonds of conjugated atoms are conjugated as well
        for (int i = 0; i < tmpBondCount; i++) {
            IBond tmpBond = aMolecule.getBond(i);
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            if (this.isMultipleBond(tmpBond) && (tmpClassification.isConjugatedAtom[tmpBegin] || tmpClassification.isConjugatedAtom[tmpEnd])) {
                tmpClassification.isConjugatedBond[i] = true;
                tmpClassification.isConjugatedAtom[tmpBegin] = true;
                tmpClassification.isConjugatedAtom[tmpEnd] = true;
            }
        }
        //tertiary and quaternary carbons
        for (int i = 0; i < tmpAtomCount; i++) {
            IAtom tmpAtom = aMolecule.getAtom(i);
            int tmpBondCountOfAtom = tmpAtom.getBondCount();
            boolean tmpHasOnlySingleBonds = true;
            for (IBond tmpBond : tmpAtom.bonds()) {
                tmpHasOnlySingleBonds &= tmpBond.getOrder() == IBond.Order.SINGLE && !this.isAromaticBondWithoutOrder(tmpBond);
            }
            if (!tmpClassification.isRingAtom[i] && !tmpClassification.isAromaticAtom[i] && tmpHasOnlySingleBonds
                    && (tmpBondCountOfAtom == 3 || tmpBondCountOfAtom == 4)) {
                tmpClassification.isTertiaryCarbon[i] = tmpBondCountOfAtom == 3;
                tmpClassification.isQuaternaryCarbon[i] = tmpBondCountOfAtom == 4;
            }
        }
        boolean tmpIsolate = this.isolateTertQuatCarbonSetting.get();
        boolean tmpSeparateFromRing = this.separateTertQuatCarbonFromRingSetting.get();
        for (int i = 0; i < tmpAtomCount; i++) {
            if (!tmpClassification.isTertQuatCarbon(i)) {
                continue;
            }
            for (IAtom tmpNeighbor : aMolecule.getAtom(i).neighbors()) {
                tmpClassification.isTertQuatCarbonBondedToRing[i] |= tmpClassification.isRingAtom[aMolecule.indexOf(tmpNeighbor)];
            }
            boolean tmpIsKeptAtRing = tmpClassification.isTertQuatCarbonBondedToRing[i] && !tmpSeparateFromRing;
            tmpClassification.isIsolatedTertQuatCarbon[i] = tmpIsolate && !tmpIsKeptAtRing;
        }
        //chain atoms: saturated, non-core, non-tertiary/quaternary atoms
        for (int i = 0; i < tmpAtomCount; i++) {
            tmpClassification.isChainAtom[i] = !tmpClassification.isCoreAtom(i) && !tmpClassification.isTertQuatCarbon(i)
                    && !tmpClassification.hasDoubleBond[i] && !tmpClassification.hasTripleBond[i]
                    && !tmpClassification.hasAromaticBondWithoutOrder[i];
        }
        //members of non-isolated tertiary/quaternary units: chain atoms bonded to a non-isolated tertiary/quaternary carbon
        for (int i = 0; i < tmpAtomCount; i++) {
            if (!tmpClassification.isTertQuatCarbon(i) || tmpClassification.isIsolatedTertQuatCarbon[i]) {
                continue;
            }
            for (IAtom tmpNeighbor : aMolecule.getAtom(i).neighbors()) {
                int tmpNeighborIndex = aMolecule.indexOf(tmpNeighbor);
                if (tmpClassification.isChainAtom[tmpNeighborIndex]) {
                    tmpClassification.isTertQuatUnitMember[tmpNeighborIndex] = true;
                }
            }
        }
        return tmpClassification;
    }

    /**
     * Decides for every bond of the given molecule whether it is kept in the fragments (see class documentation).
     *
     * @param aMolecule the molecule to fragment
     * @param aClassification classification of the molecule's atoms and bonds
     * @return array indexed by bond index, true if the bond is kept
     */
    protected boolean[] determineKeptBonds(IAtomContainer aMolecule, MoleculeClassification aClassification) {
        boolean tmpSeparateFromRing = this.separateTertQuatCarbonFromRingSetting.get();
        boolean[] tmpIsBondKept = new boolean[aMolecule.getBondCount()];
        for (int i = 0; i < aMolecule.getBondCount(); i++) {
            IBond tmpBond = aMolecule.getBond(i);
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            if (aClassification.isRingBond[i] || aClassification.isConjugatedBond[i] || this.isMultipleBond(tmpBond)
                    || this.isAromaticBondWithoutOrder(tmpBond)) {
                //rules 1-3: ring, conjugated and multiple bonds are never cut
                tmpIsBondKept[i] = true;
            } else if (aClassification.isTertQuatCarbon(tmpBegin) || aClassification.isTertQuatCarbon(tmpEnd)) {
                //rule 4: bonds of tertiary/quaternary carbons
                tmpIsBondKept[i] = this.isTertQuatBondKept(aClassification, tmpBegin, tmpEnd, tmpSeparateFromRing)
                        && this.isTertQuatBondKept(aClassification, tmpEnd, tmpBegin, tmpSeparateFromRing);
            } else {
                //rule 6: bonds between chain atoms that are not part of a tertiary/quaternary unit
                tmpIsBondKept[i] = aClassification.isChainAtom[tmpBegin] && aClassification.isChainAtom[tmpEnd]
                        && !aClassification.isTertQuatUnitMember[tmpBegin] && !aClassification.isTertQuatUnitMember[tmpEnd];
            }
        }
        return tmpIsBondKept;
    }

    /**
     * Dissects the linear chains (connected chain atoms, see class documentation) into chains of the given maximum
     * length by marking the respective bonds as cut. Every chain is walked from its terminal atom with the lowest index.
     *
     * @param aMolecule the molecule to fragment
     * @param aClassification classification of the molecule's atoms and bonds
     * @param anIsBondKeptArray kept bonds, changed by this method
     * @param aMaxChainLength maximum number of atoms per chain fragment
     */
    protected void dissectLinearChains(IAtomContainer aMolecule, MoleculeClassification aClassification,
                                       boolean[] anIsBondKeptArray, int aMaxChainLength) {
        int tmpAtomCount = aMolecule.getAtomCount();
        int[] tmpPositionInSegment = new int[tmpAtomCount];
        for (int tmpStart = 0; tmpStart < tmpAtomCount; tmpStart++) {
            //chains are started at the terminal chain atom (at most one kept chain bond) with the lowest index
            if (!aClassification.isChainAtom[tmpStart] || aClassification.isTertQuatUnitMember[tmpStart]
                    || tmpPositionInSegment[tmpStart] != 0 || this.getKeptChainBonds(aMolecule, anIsBondKeptArray, tmpStart).size() > 1) {
                continue;
            }
            tmpPositionInSegment[tmpStart] = 1;
            int tmpCurrent = tmpStart;
            boolean tmpHasNext = true;
            while (tmpHasNext) {
                tmpHasNext = false;
                for (IBond tmpBond : this.getKeptChainBonds(aMolecule, anIsBondKeptArray, tmpCurrent)) {
                    int tmpNext = aMolecule.indexOf(tmpBond.getOther(aMolecule.getAtom(tmpCurrent)));
                    if (tmpPositionInSegment[tmpNext] != 0) {
                        continue;
                    }
                    if (tmpPositionInSegment[tmpCurrent] < aMaxChainLength) {
                        tmpPositionInSegment[tmpNext] = tmpPositionInSegment[tmpCurrent] + 1;
                    } else {
                        anIsBondKeptArray[aMolecule.indexOf(tmpBond)] = false;
                        tmpPositionInSegment[tmpNext] = 1;
                    }
                    tmpCurrent = tmpNext;
                    tmpHasNext = true;
                    break;
                }
            }
        }
    }

    /**
     * Builds the fragments from the kept bonds: the connected parts of the kept bonds are copied into new atom
     * containers, cut bonds of tertiary/quaternary carbons are replaced by pseudo atoms, free valences are saturated
     * with implicit hydrogens (if set) and the stereo chemistry is transferred (if set).
     *
     * @param aMolecule the molecule to fragment
     * @param aClassification classification of the molecule's atoms and bonds
     * @param anIsBondKeptArray kept bonds
     * @return list of fragments, ordered by their atom with the lowest index
     */
    protected List<IAtomContainer> buildFragments(IAtomContainer aMolecule, MoleculeClassification aClassification,
                                                  boolean[] anIsBondKeptArray) {
        int tmpAtomCount = aMolecule.getAtomCount();
        //connected parts of the kept bonds (union-find)
        int[] tmpParent = new int[tmpAtomCount];
        for (int i = 0; i < tmpAtomCount; i++) {
            tmpParent[i] = i;
        }
        for (int i = 0; i < aMolecule.getBondCount(); i++) {
            if (anIsBondKeptArray[i]) {
                IBond tmpBond = aMolecule.getBond(i);
                int tmpBeginRoot = this.findRoot(tmpParent, aMolecule.indexOf(tmpBond.getBegin()));
                int tmpEndRoot = this.findRoot(tmpParent, aMolecule.indexOf(tmpBond.getEnd()));
                tmpParent[Math.max(tmpBeginRoot, tmpEndRoot)] = Math.min(tmpBeginRoot, tmpEndRoot);
            }
        }
        boolean tmpSaturate = this.fragmentSaturationSetting.get().equals(IMoleculeFragmenter.FragmentSaturationOption.HYDROGEN_SATURATION);
        Map<Integer, IAtomContainer> tmpFragmentByRoot = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(tmpAtomCount));
        List<IAtomContainer> tmpFragments = new ArrayList<>();
        IAtom[] tmpAtomCopies = new IAtom[tmpAtomCount];
        for (int i = 0; i < tmpAtomCount; i++) {
            IAtomContainer tmpFragment = tmpFragmentByRoot.computeIfAbsent(this.findRoot(tmpParent, i), aRoot -> {
                IAtomContainer tmpNewFragment = this.builder.newAtomContainer();
                tmpFragments.add(tmpNewFragment);
                return tmpNewFragment;
            });
            tmpFragment.addAtom(this.copyAtom(aMolecule.getAtom(i), i));
            tmpAtomCopies[i] = tmpFragment.getAtom(tmpFragment.getAtomCount() - 1);
        }
        IBond[] tmpBondCopies = new IBond[aMolecule.getBondCount()];
        int[] tmpCutBondOrderSum = new int[tmpAtomCount];
        for (int i = 0; i < aMolecule.getBondCount(); i++) {
            IBond tmpBond = aMolecule.getBond(i);
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            if (anIsBondKeptArray[i]) {
                IAtomContainer tmpFragment = tmpFragmentByRoot.get(this.findRoot(tmpParent, tmpBegin));
                IBond tmpBondCopy = this.builder.newInstance(IBond.class, tmpAtomCopies[tmpBegin], tmpAtomCopies[tmpEnd],
                        tmpBond.getOrder() == null ? IBond.Order.UNSET : tmpBond.getOrder());
                tmpBondCopy.setIsAromatic(tmpBond.isAromatic());
                tmpBondCopy.setProperty(AlkylStructureFragmenter.INTERNAL_ASF_BOND_INDEX_PROPERTY_KEY, i);
                tmpFragment.addBond(tmpBondCopy);
                tmpBondCopies[i] = tmpFragment.getBond(tmpFragment.getBondCount() - 1);
                continue;
            }
            //rule 5: cut bonds of tertiary/quaternary carbons are replaced by pseudo atoms, all others by hydrogens
            for (int tmpAtomIndex : new int[]{tmpBegin, tmpEnd}) {
                if (aClassification.isTertQuatCarbon(tmpAtomIndex)) {
                    IAtomContainer tmpFragment = tmpFragmentByRoot.get(this.findRoot(tmpParent, tmpAtomIndex));
                    IPseudoAtom tmpPseudoAtom = this.builder.newInstance(IPseudoAtom.class);
                    tmpPseudoAtom.setImplicitHydrogenCount(0);
                    tmpFragment.addAtom(tmpPseudoAtom);
                    tmpFragment.addBond(this.builder.newInstance(IBond.class, tmpAtomCopies[tmpAtomIndex],
                            tmpFragment.getAtom(tmpFragment.getAtomCount() - 1), IBond.Order.SINGLE));
                } else {
                    tmpCutBondOrderSum[tmpAtomIndex] += tmpBond.getOrder() == null || tmpBond.getOrder() == IBond.Order.UNSET
                            ? 1 : tmpBond.getOrder().numeric();
                }
            }
        }
        if (tmpSaturate) {
            for (int i = 0; i < tmpAtomCount; i++) {
                tmpAtomCopies[i].setImplicitHydrogenCount(tmpAtomCopies[i].getImplicitHydrogenCount() + tmpCutBondOrderSum[i]);
            }
        }
        if (this.preserveStereoChemistrySetting.get()) {
            this.transferStereoChemistry(aMolecule, anIsBondKeptArray, tmpAtomCopies, tmpBondCopies);
            for (IAtomContainer tmpFragment : tmpFragments) {
                this.removeNonStereocenterElements(tmpFragment);
            }
        }
        return tmpFragments;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Adds the given setting to the settings list and its tooltip and display name to the respective maps.
     *
     * @param aSetting the setting property
     * @param aMessageKeyPart part of the message keys of the setting's tooltip and display name
     */
    private void addSetting(Property<?> aSetting, String aMessageKeyPart) {
        this.settings.add(aSetting);
        this.settingNameToTooltipTextMap.put(aSetting.getName(),
                Message.get("AlkylStructureFragmenter." + aMessageKeyPart + ".tooltip"));
        this.settingNameToDisplayNameMap.put(aSetting.getName(),
                Message.get("AlkylStructureFragmenter." + aMessageKeyPart + ".displayName"));
    }

    /**
     * Checks whether the given molecule can be fragmented (see {@link #shouldBeFiltered(IAtomContainer)}),
     * independent of the setting for keeping non-fragmentable molecules.
     *
     * @param aMolecule molecule to check
     * @return true if the molecule can be fragmented
     */
    private boolean isFragmentable(IAtomContainer aMolecule) {
        if (Objects.isNull(aMolecule) || aMolecule.isEmpty()) {
            return false;
        }
        boolean tmpContainsCarbon = false;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom) || tmpAtom instanceof IPseudoAtom || Objects.isNull(tmpAtom.getAtomicNumber())) {
                return false;
            }
            int tmpAtomicNumber = tmpAtom.getAtomicNumber();
            if (tmpAtomicNumber == IElement.C) {
                tmpContainsCarbon = true;
            } else if (tmpAtomicNumber != IElement.H) {
                return false;
            }
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            if (this.isAromaticBondWithoutOrder(tmpBond)) {
                continue;
            }
            IBond.Order tmpOrder = tmpBond.getOrder();
            if (tmpOrder != IBond.Order.SINGLE && tmpOrder != IBond.Order.DOUBLE && tmpOrder != IBond.Order.TRIPLE
                    && tmpOrder != IBond.Order.QUADRUPLE) {
                return false;
            }
        }
        return tmpContainsCarbon;
    }

    /**
     * Sets missing implicit hydrogen counts of the given molecule's atoms.
     *
     * @param aMolecule molecule to check
     */
    private void ensureImplicitHydrogenCounts(IAtomContainer aMolecule) {
        boolean tmpIsCountMissing = false;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpIsCountMissing |= Objects.isNull(tmpAtom.getImplicitHydrogenCount());
        }
        if (!tmpIsCountMissing) {
            return;
        }
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
            CDKHydrogenAdder.getInstance(this.builder).addImplicitHydrogens(aMolecule);
        } catch (CDKException anException) {
            AlkylStructureFragmenter.LOGGER.log(Level.WARNING, String.format("Implicit hydrogen counts could not be determined for molecule %s: %s",
                    aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException));
        }
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom.getImplicitHydrogenCount())) {
                tmpAtom.setImplicitHydrogenCount(0);
            }
        }
    }

    /**
     * Checks whether the bond between the given tertiary/quaternary carbon and the given other atom is kept (rule 4).
     * Returns true if the first atom is no tertiary/quaternary carbon, so that both directions of a bond can be checked.
     *
     * @param aClassification classification of the molecule's atoms and bonds
     * @param aTertQuatIndex index of the (possibly) tertiary/quaternary carbon
     * @param anOtherIndex index of the other bond atom
     * @param isSeparatedFromRing value of the setting for separating tertiary/quaternary carbons from rings
     * @return true if the bond is kept from the point of view of the first atom
     */
    private boolean isTertQuatBondKept(MoleculeClassification aClassification, int aTertQuatIndex, int anOtherIndex,
                                       boolean isSeparatedFromRing) {
        if (!aClassification.isTertQuatCarbon(aTertQuatIndex)) {
            return true;
        }
        if (aClassification.isIsolatedTertQuatCarbon[aTertQuatIndex]) {
            return false;
        }
        if (aClassification.isRingAtom[anOtherIndex]) {
            return !isSeparatedFromRing;
        }
        if (aClassification.isConjugatedAtom[anOtherIndex]) {
            return true;
        }
        if (aClassification.isTertQuatCarbon(anOtherIndex)) {
            return !aClassification.isIsolatedTertQuatCarbon[anOtherIndex];
        }
        return aClassification.isChainAtom[anOtherIndex];
    }

    /**
     * Returns the kept bonds of the given atom to other chain atoms that are not part of a tertiary/quaternary unit.
     *
     * @param aMolecule the molecule to fragment
     * @param anIsBondKeptArray kept bonds
     * @param anAtomIndex index of the chain atom
     * @return list of the kept chain bonds of the atom
     */
    private List<IBond> getKeptChainBonds(IAtomContainer aMolecule, boolean[] anIsBondKeptArray, int anAtomIndex) {
        List<IBond> tmpChainBonds = new ArrayList<>(2);
        for (IBond tmpBond : aMolecule.getAtom(anAtomIndex).bonds()) {
            if (anIsBondKeptArray[aMolecule.indexOf(tmpBond)]) {
                tmpChainBonds.add(tmpBond);
            }
        }
        return tmpChainBonds;
    }

    /**
     * Creates a copy of the given atom for a fragment, carrying the given index as property.
     *
     * @param anAtom atom to copy
     * @param anIndex index of the atom in the fragmented molecule
     * @return the copy
     */
    private IAtom copyAtom(IAtom anAtom, int anIndex) {
        IAtom tmpCopy = this.builder.newInstance(IAtom.class, anAtom.getSymbol());
        tmpCopy.setAtomicNumber(anAtom.getAtomicNumber());
        tmpCopy.setImplicitHydrogenCount(anAtom.getImplicitHydrogenCount());
        tmpCopy.setFormalCharge(anAtom.getFormalCharge());
        tmpCopy.setMassNumber(anAtom.getMassNumber());
        tmpCopy.setIsAromatic(anAtom.isAromatic());
        if (Objects.nonNull(anAtom.getPoint2d())) {
            tmpCopy.setPoint2d(new Point2d(anAtom.getPoint2d()));
        }
        if (Objects.nonNull(anAtom.getPoint3d())) {
            tmpCopy.setPoint3d(new Point3d(anAtom.getPoint3d()));
        }
        tmpCopy.setProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY, anIndex);
        return tmpCopy;
    }

    /**
     * Transfers the tetrahedral and double bond stereo elements of the given molecule to the fragments if the stereo
     * element is complete in a fragment, i.e. the focus and all carriers are present and bonded as in the molecule.
     *
     * @param aMolecule the fragmented molecule
     * @param anIsBondKeptArray kept bonds
     * @param anAtomCopyArray fragment atoms, indexed by the index of the original atom
     * @param aBondCopyArray fragment bonds, indexed by the index of the original bond (null if cut)
     */
    private void transferStereoChemistry(IAtomContainer aMolecule, boolean[] anIsBondKeptArray, IAtom[] anAtomCopyArray,
                                         IBond[] aBondCopyArray) {
        for (IStereoElement tmpElement : aMolecule.stereoElements()) {
            try {
                if (tmpElement instanceof ITetrahedralChirality tmpChirality) {
                    IAtom tmpFocus = tmpChirality.getChiralAtom();
                    IAtom[] tmpLigands = tmpChirality.getLigands();
                    IAtom[] tmpLigandCopies = new IAtom[tmpLigands.length];
                    boolean tmpIsComplete = true;
                    for (int i = 0; i < tmpLigands.length; i++) {
                        IBond tmpBond = aMolecule.getBond(tmpFocus, tmpLigands[i]);
                        //implicit hydrogen is represented by the focus itself
                        tmpIsComplete &= tmpLigands[i].equals(tmpFocus) || (tmpBond != null && anIsBondKeptArray[aMolecule.indexOf(tmpBond)]);
                        tmpLigandCopies[i] = anAtomCopyArray[aMolecule.indexOf(tmpLigands[i])];
                    }
                    if (tmpIsComplete) {
                        IAtom tmpFocusCopy = anAtomCopyArray[aMolecule.indexOf(tmpFocus)];
                        tmpFocusCopy.getContainer().addStereoElement(
                                new TetrahedralChirality(tmpFocusCopy, tmpLigandCopies, tmpChirality.getStereo()));
                    }
                } else if (tmpElement instanceof IDoubleBondStereochemistry tmpDoubleBondStereo) {
                    IBond tmpFocusBond = tmpDoubleBondStereo.getStereoBond();
                    IBond[] tmpCarrierBonds = tmpDoubleBondStereo.getBonds();
                    int tmpFocusBondIndex = aMolecule.indexOf(tmpFocusBond);
                    int tmpFirstCarrierIndex = aMolecule.indexOf(tmpCarrierBonds[0]);
                    int tmpSecondCarrierIndex = aMolecule.indexOf(tmpCarrierBonds[1]);
                    if (anIsBondKeptArray[tmpFocusBondIndex] && anIsBondKeptArray[tmpFirstCarrierIndex] && anIsBondKeptArray[tmpSecondCarrierIndex]) {
                        IBond tmpFocusBondCopy = aBondCopyArray[tmpFocusBondIndex];
                        tmpFocusBondCopy.getContainer().addStereoElement(new DoubleBondStereochemistry(tmpFocusBondCopy,
                                new IBond[]{aBondCopyArray[tmpFirstCarrierIndex], aBondCopyArray[tmpSecondCarrierIndex]},
                                tmpDoubleBondStereo.getStereo()));
                    }
                }
                //other stereo element types are not transferred
            } catch (Exception anException) {
                AlkylStructureFragmenter.LOGGER.log(Level.WARNING, String.format("Stereo element could not be transferred for molecule %s: %s",
                        aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException), anException);
            }
        }
    }

    /**
     * Removes all stereo elements of the given fragment whose focus is no longer a stereocenter, because fragmentation
     * made formerly different substituents identical (e.g. a sec-butyl stereocenter becomes an isopropyl group).
     * Otherwise, identical fragments would have different stereo SMILES codes.
     *
     * @param aFragment the fragment
     */
    private void removeNonStereocenterElements(IAtomContainer aFragment) {
        if (!aFragment.stereoElements().iterator().hasNext()) {
            return;
        }
        Stereocenters tmpStereocenters;
        try {
            tmpStereocenters = Stereocenters.of(aFragment);
        } catch (Exception anException) {
            AlkylStructureFragmenter.LOGGER.log(Level.WARNING, "Stereocenters of a fragment could not be determined, stereo elements are kept.", anException);
            return;
        }
        List<IStereoElement> tmpKeptElements = new ArrayList<>();
        for (IStereoElement tmpElement : aFragment.stereoElements()) {
            boolean tmpIsStereocenter = true;
            if (tmpElement instanceof ITetrahedralChirality tmpChirality) {
                tmpIsStereocenter = tmpStereocenters.isStereocenter(aFragment.indexOf(tmpChirality.getChiralAtom()));
            } else if (tmpElement instanceof IDoubleBondStereochemistry tmpDoubleBondStereo) {
                tmpIsStereocenter = tmpStereocenters.isStereocenter(aFragment.indexOf(tmpDoubleBondStereo.getStereoBond().getBegin()))
                        && tmpStereocenters.isStereocenter(aFragment.indexOf(tmpDoubleBondStereo.getStereoBond().getEnd()));
            }
            if (tmpIsStereocenter) {
                tmpKeptElements.add(tmpElement);
            }
        }
        aFragment.setStereoElements(tmpKeptElements);
    }

    /**
     * Checks whether the given bond is aromatic but has no assigned bond order, e.g. in remnants of aromatic rings that
     * cannot be kekulized after a previous fragmentation.
     *
     * @param aBond the bond to check
     * @return true if the bond is aromatic and its order is unset
     */
    private boolean isAromaticBondWithoutOrder(IBond aBond) {
        return aBond.isAromatic() && (Objects.isNull(aBond.getOrder()) || aBond.getOrder() == IBond.Order.UNSET);
    }

    /**
     * Checks whether the given bond is a double, triple or quadruple bond.
     *
     * @param aBond the bond to check
     * @return true if the bond is a multiple bond
     */
    private boolean isMultipleBond(IBond aBond) {
        return aBond.getOrder() == IBond.Order.DOUBLE || aBond.getOrder() == IBond.Order.TRIPLE
                || aBond.getOrder() == IBond.Order.QUADRUPLE;
    }

    /**
     * Checks whether the given atom is part of a double bond (and no triple bond), aromatic or part of an aromatic bond
     * without bond order, i.e. whether a single bond to it can be part of a conjugated system of double bonds.
     *
     * @param aClassification classification of the molecule's atoms and bonds
     * @param anAtomIndex index of the atom
     * @return true if the atom is a double bond or aromatic atom
     */
    private boolean isDoubleBondOrAromaticAtom(MoleculeClassification aClassification, int anAtomIndex) {
        return (aClassification.hasDoubleBond[anAtomIndex] && !aClassification.hasTripleBond[anAtomIndex])
                || aClassification.isAromaticAtom[anAtomIndex] || aClassification.hasAromaticBondWithoutOrder[anAtomIndex];
    }

    /**
     * Returns the root of the given element in the given union-find array (with path halving).
     *
     * @param aParentArray union-find array
     * @param anElement element
     * @return root of the element
     */
    private int findRoot(int[] aParentArray, int anElement) {
        int tmpElement = anElement;
        while (aParentArray[tmpElement] != tmpElement) {
            aParentArray[tmpElement] = aParentArray[aParentArray[tmpElement]];
            tmpElement = aParentArray[tmpElement];
        }
        return tmpElement;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Inner class MoleculeClassification">
    /**
     * Classification of the atoms and bonds of a molecule into the substructures relevant for the fragmentation rules,
     * indexed by the atom and bond indices of the classified molecule.
     */
    protected static final class MoleculeClassification {
        /**
         * Atoms in rings.
         */
        final boolean[] isRingAtom;
        /**
         * Bonds in rings.
         */
        final boolean[] isRingBond;
        /**
         * Atoms with the aromaticity flag.
         */
        final boolean[] isAromaticAtom;
        /**
         * Atoms of conjugated pi systems, including all aromatic atoms.
         */
        final boolean[] isConjugatedAtom;
        /**
         * Bonds of conjugated pi systems.
         */
        final boolean[] isConjugatedBond;
        /**
         * Atoms of double bonds.
         */
        final boolean[] hasDoubleBond;
        /**
         * Atoms of triple (or quadruple) bonds.
         */
        final boolean[] hasTripleBond;
        /**
         * Atoms of aromatic bonds without bond order.
         */
        final boolean[] hasAromaticBondWithoutOrder;
        /**
         * Tertiary carbons (non-cyclic, non-aromatic, three single bonds).
         */
        final boolean[] isTertiaryCarbon;
        /**
         * Quaternary carbons (non-cyclic, non-aromatic, four single bonds).
         */
        final boolean[] isQuaternaryCarbon;
        /**
         * Tertiary/quaternary carbons bonded to a ring atom.
         */
        final boolean[] isTertQuatCarbonBondedToRing;
        /**
         * Tertiary/quaternary carbons cut from all neighbors (isolation setting, not kept at a ring).
         */
        final boolean[] isIsolatedTertQuatCarbon;
        /**
         * Saturated carbons that are neither core atoms nor tertiary/quaternary carbons.
         */
        final boolean[] isChainAtom;
        /**
         * Chain atoms bonded to a non-isolated tertiary/quaternary carbon, they form a unit with it.
         */
        final boolean[] isTertQuatUnitMember;

        /**
         * Constructor, all arrays are initialized with false.
         *
         * @param anAtomCount number of atoms
         * @param aBondCount number of bonds
         */
        MoleculeClassification(int anAtomCount, int aBondCount) {
            this.isRingAtom = new boolean[anAtomCount];
            this.isRingBond = new boolean[aBondCount];
            this.isAromaticAtom = new boolean[anAtomCount];
            this.isConjugatedAtom = new boolean[anAtomCount];
            this.isConjugatedBond = new boolean[aBondCount];
            this.hasDoubleBond = new boolean[anAtomCount];
            this.hasTripleBond = new boolean[anAtomCount];
            this.hasAromaticBondWithoutOrder = new boolean[anAtomCount];
            this.isTertiaryCarbon = new boolean[anAtomCount];
            this.isQuaternaryCarbon = new boolean[anAtomCount];
            this.isTertQuatCarbonBondedToRing = new boolean[anAtomCount];
            this.isIsolatedTertQuatCarbon = new boolean[anAtomCount];
            this.isChainAtom = new boolean[anAtomCount];
            this.isTertQuatUnitMember = new boolean[anAtomCount];
        }

        /**
         * Returns whether the atom with the given index is a tertiary or quaternary carbon.
         *
         * @param anAtomIndex atom index
         * @return true if tertiary or quaternary carbon
         */
        boolean isTertQuatCarbon(int anAtomIndex) {
            return this.isTertiaryCarbon[anAtomIndex] || this.isQuaternaryCarbon[anAtomIndex];
        }

        /**
         * Returns whether the atom with the given index belongs to a core substructure (ring or conjugated system).
         *
         * @param anAtomIndex atom index
         * @return true if ring or conjugated atom
         */
        boolean isCoreAtom(int anAtomIndex) {
            return this.isRingAtom[anAtomIndex] || this.isConjugatedAtom[anAtomIndex];
        }
    }
    //</editor-fold>
}
