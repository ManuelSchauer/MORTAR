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
import de.unijena.cheminf.mortar.model.util.SimpleSmartsPatternListProperty;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;

import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Fragmenter that extracts substructures defined by user-supplied SMARTS patterns. Every match of a pattern becomes a
 * fragment that consists of the matched atoms and all bonds between them. Optionally, the remainder of the molecule,
 * i.e. the connected parts of the atoms that are not part of any returned match, is returned as well, so that the
 * molecule is partitioned like by the Ertl fragmenter.
 * <br>
 * <br><b>Patterns</b>: The patterns are given as text with one entry per line, a SMARTS pattern optionally followed by
 * whitespace and a name (see {@link SimpleSmartsPatternListProperty}). Every match fragment is annotated with the name
 * of its pattern, or with the SMARTS pattern itself if the pattern is unnamed, as category (see
 * {@link IMoleculeFragmenter#FRAGMENT_CATEGORY_PROPERTY_KEY}); remainder fragments get the category
 * {@link #FRAGMENT_CATEGORY_REMAINDER_VALUE}. The patterns are matched with the CDK SMARTS implementation
 * ({@link SmartsPattern}). Prior to matching, ring atoms and bonds are marked and the Daylight aromaticity model is
 * applied to a copy of the input molecule (see {@link SmartsPattern#prepare(IAtomContainer)}), which is what the
 * SMARTS semantics assume; existing aromaticity flags of the input are not used. Recursive SMARTS and component-level
 * grouping are supported. A pattern whose components are separated by "." matches atoms that do not have to be
 * connected; such a match is returned as one fragment, i.e. it can consist of several unconnected parts, because the
 * pattern defines these parts as one unit.
 * <br>
 * <br><b>Overlapping matches</b>: If the overlap filter is turned on (default), the matches are accepted greedily and
 * no atom is part of more than one returned match. The order of the patterns in the list is their priority, i.e. a
 * match is only accepted if none of its atoms is part of an accepted match of the same or an earlier pattern. The
 * matches of one pattern are processed in the lexicographic order of the sorted canonical ranks of their atoms, so
 * that the result does not depend on the input atom order. If the filter is turned off, all unique matches (unique
 * sets of matched atoms) are returned, and they may overlap.
 * <br>
 * <br><b>Fragment construction</b>: Every fragment is built by copying the respective atoms and the bonds between
 * them into a new atom container, without removing bonds from a molecule, so the input molecule is not altered. Free
 * valences resulting from bonds to atoms outside the fragment are saturated with implicit hydrogen atoms; a cut bond
 * of an aromatic system counts as a single bond. A Kekulé structure is assigned to the aromatic part of this saturated
 * fragment, and aromaticity is detected again (Daylight model), so that completely contained aromatic rings stay
 * aromatic, e.g. the pyrrole ring of indole, while the rest of a partially contained aromatic system gets a Kekulé
 * structure. If the aromatic part of a fragment has no Kekulé structure on its own, e.g. three atoms of a benzene
 * ring, the bond orders of the (Kekulé form of the) input molecule are used instead; since the Kekulé structure of the
 * input may depend on its atom order, such rare fragments of fused aromatic systems may then depend on it as well. If
 * no saturation is set, the added hydrogen atoms are removed again afterwards, so that both saturation options agree
 * in bond orders and aromaticity. Stereo elements of the input molecule are only transferred to a fragment if all of
 * their atoms and bonds are part of the fragment, i.e. stereocentres and stereo double bonds at the border of a
 * fragment lose their configuration.
 * <br>
 * <br><b>Input restrictions</b>: Molecules with unset bond orders that cannot be resolved by assigning a Kekulé
 * structure are filtered, because SMARTS bond primitives cannot be evaluated on them. If the pattern list is empty,
 * no molecule is filtered: every molecule is returned as remainder, split into its unconnected parts, or no fragments
 * are returned if only the matched substructures are requested.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SmartsFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Enum SmartsFragmenterReturnedFragmentsOption">
    /**
     * Enum for defining which fragments are returned, only the matched substructures or additionally the remainder of
     * the molecule.
     */
    public enum SmartsFragmenterReturnedFragmentsOption implements IDisplayEnum {
        /**
         * Return only the substructures matched by the SMARTS patterns.
         */
        MATCHED_SUBSTRUCTURES(
                Message.get("SmartsFragmenter.ReturnedFragmentsOption.MatchedSubstructures.displayName"),
                Message.get("SmartsFragmenter.ReturnedFragmentsOption.MatchedSubstructures.tooltip")),
        /**
         * Return the matched substructures and the connected parts of the atoms in no returned match.
         */
        MATCHED_SUBSTRUCTURES_AND_REMAINDER(
                Message.get("SmartsFragmenter.ReturnedFragmentsOption.MatchedSubstructuresAndRemainder.displayName"),
                Message.get("SmartsFragmenter.ReturnedFragmentsOption.MatchedSubstructuresAndRemainder.tooltip"));
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
        private SmartsFragmenterReturnedFragmentsOption(String aDisplayName, String aTooltip) {
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
    //<editor-fold desc="Record CompiledPattern">
    /**
     * A compiled SMARTS pattern together with the category assigned to its match fragments.
     *
     * @param pattern the compiled SMARTS pattern, automatic preparation of the target molecule is disabled
     * @param category name of the pattern or, if it is unnamed, the SMARTS pattern
     */
    private record CompiledPattern(SmartsPattern pattern, String category) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm used in this fragmenter.
     */
    public static final String ALGORITHM_NAME = "Custom SMARTS Fragmenter";

    /**
     * Default pattern list, empty.
     */
    public static final String PATTERNS_SETTING_DEFAULT = "";

    /**
     * Default returned fragments option, the matched substructures and the remainder.
     */
    public static final SmartsFragmenterReturnedFragmentsOption RETURNED_FRAGMENTS_OPTION_DEFAULT =
            SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES_AND_REMAINDER;

    /**
     * Default value of the filter overlapping matches setting.
     */
    public static final boolean FILTER_OVERLAPPING_MATCHES_SETTING_DEFAULT = true;

    /**
     * Remainder fragments, i.e. the connected parts of the atoms in no returned match, are assigned this value for
     * the property with key {@link IMoleculeFragmenter#FRAGMENT_CATEGORY_PROPERTY_KEY}.
     */
    public static final String FRAGMENT_CATEGORY_REMAINDER_VALUE = "SmartsFragmenter.Remainder";
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(SmartsFragmenter.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    //note: since Java 21, the javadoc build complains about "double comments" when there is a comment
    // for the get() method of the property and the private property itself as well
    private final SimpleSmartsPatternListProperty patternsSetting;

    private final SimpleIDisplayEnumConstantProperty returnedFragmentsSetting;

    private final SimpleBooleanProperty filterOverlappingMatchesSetting;

    private final SimpleIDisplayEnumConstantProperty fragmentSaturationSetting;

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
    //<editor-fold desc="Private variables">
    /**
     * Compiled SMARTS patterns of the pattern list text stored in {@link #compiledPatternsText}; compiled lazily when
     * the patterns are needed for the first time after the pattern list has changed.
     */
    private List<CompiledPattern> compiledPatterns;

    /**
     * Pattern list text that the compiled patterns were created from, null if no patterns were compiled yet.
     */
    private String compiledPatternsText;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, all settings are initialised with their default values as declared in the respective public constants.
     */
    public SmartsFragmenter() {
        int tmpNumberOfSettings = 4;
        int tmpInitialCapacityForSettingMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        //these names are for internal use, the language-specific display names for the GUI are stored in the map
        this.patternsSetting = new SimpleSmartsPatternListProperty(this, "SMARTS patterns setting",
                SmartsFragmenter.PATTERNS_SETTING_DEFAULT) {
            @Override
            public void set(String newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    SmartsFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                            Message.get("Fragmenter.IllegalSettingValue.Header"),
                            anException.toString(),
                            anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.patternsSetting.getName(),
                Message.get("SmartsFragmenter.patternsSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.patternsSetting.getName(),
                Message.get("SmartsFragmenter.patternsSetting.displayName"));
        this.returnedFragmentsSetting = this.createEnumSetting("Returned fragments setting",
                SmartsFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, SmartsFragmenterReturnedFragmentsOption.class);
        this.settingNameTooltipTextMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("SmartsFragmenter.returnedFragmentsSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.returnedFragmentsSetting.getName(),
                Message.get("SmartsFragmenter.returnedFragmentsSetting.displayName"));
        this.filterOverlappingMatchesSetting = new SimpleBooleanProperty(this, "Filter overlapping matches setting",
                SmartsFragmenter.FILTER_OVERLAPPING_MATCHES_SETTING_DEFAULT);
        this.settingNameTooltipTextMap.put(this.filterOverlappingMatchesSetting.getName(),
                Message.get("SmartsFragmenter.filterOverlappingMatchesSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.filterOverlappingMatchesSetting.getName(),
                Message.get("SmartsFragmenter.filterOverlappingMatchesSetting.displayName"));
        this.fragmentSaturationSetting = this.createEnumSetting("Fragment saturation setting",
                IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT, IMoleculeFragmenter.FragmentSaturationOption.class);
        this.settingNameTooltipTextMap.put(this.fragmentSaturationSetting.getName(),
                Message.get("SmartsFragmenter.fragmentSaturationSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.fragmentSaturationSetting.getName(),
                Message.get("SmartsFragmenter.fragmentSaturationSetting.displayName"));
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        this.settings.add(this.patternsSetting);
        this.settings.add(this.returnedFragmentsSetting);
        this.settings.add(this.filterOverlappingMatchesSetting);
        this.settings.add(this.fragmentSaturationSetting);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the pattern list text, one SMARTS pattern per line, optionally followed by whitespace and a name (see
     * {@link SimpleSmartsPatternListProperty}).
     *
     * @return pattern list text
     */
    public String getPatternsSetting() {
        return this.patternsSetting.get();
    }

    /**
     * Returns the property object of the patterns setting that can be used to configure this setting.
     *
     * @return property object of the patterns setting
     */
    public SimpleSmartsPatternListProperty patternsSettingProperty() {
        return this.patternsSetting;
    }

    /**
     * Returns the currently set option for the returned fragments setting.
     *
     * @return enum constant of the set option
     */
    public SmartsFragmenterReturnedFragmentsOption getReturnedFragmentsSetting() {
        return (SmartsFragmenterReturnedFragmentsOption) this.returnedFragmentsSetting.get();
    }

    /**
     * Returns the property object of the returned fragments setting that can be used to configure this setting.
     *
     * @return property object of the returned fragments setting
     */
    public SimpleIDisplayEnumConstantProperty returnedFragmentsSettingProperty() {
        return this.returnedFragmentsSetting;
    }

    /**
     * Returns the current state of the filter overlapping matches setting.
     *
     * @return true if overlapping matches are filtered
     */
    public boolean getFilterOverlappingMatchesSetting() {
        return this.filterOverlappingMatchesSetting.get();
    }

    /**
     * Returns the property object of the filter overlapping matches setting that can be used to configure this setting.
     *
     * @return property object of the filter overlapping matches setting
     */
    public SimpleBooleanProperty filterOverlappingMatchesSettingProperty() {
        return this.filterOverlappingMatchesSetting;
    }

    /**
     * Returns the currently set option for saturating free valences on returned fragment molecules.
     *
     * @return the set option
     */
    public IMoleculeFragmenter.FragmentSaturationOption getFragmentSaturationSetting() {
        return (IMoleculeFragmenter.FragmentSaturationOption) this.fragmentSaturationSetting.get();
    }

    /**
     * Returns the property representing the setting for fragment saturation.
     *
     * @return setting property for fragment saturation
     */
    public SimpleIDisplayEnumConstantProperty fragmentSaturationSettingProperty() {
        return this.fragmentSaturationSetting;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set">
    /**
     * Sets the pattern list text, one SMARTS pattern per line, optionally followed by whitespace and a name; blank
     * lines and lines starting with "#" are ignored (see {@link SimpleSmartsPatternListProperty}). The order of the
     * patterns is their priority if overlapping matches are filtered.
     *
     * @param aPatternListText the pattern list text
     * @throws NullPointerException if the given text is null
     * @throws IllegalArgumentException if a line of the given text holds an invalid SMARTS pattern
     */
    public void setPatternsSetting(String aPatternListText) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aPatternListText, "Given pattern list text is null.");
        List<SimpleSmartsPatternListProperty.ParsedLine> tmpInvalidLines = SimpleSmartsPatternListProperty.validate(aPatternListText);
        if (!tmpInvalidLines.isEmpty()) {
            //checked here to not trigger the GUI alert of the property
            throw new IllegalArgumentException("The given SMARTS pattern list contains invalid lines:"
                    + SimpleSmartsPatternListProperty.LINE_SEPARATOR
                    + SimpleSmartsPatternListProperty.describeInvalidLines(tmpInvalidLines));
        }
        this.patternsSetting.set(aPatternListText);
    }

    /**
     * Sets the returned fragments setting, defining whether only the matched substructures or additionally the
     * remainder of the molecule are returned.
     *
     * @param anOption a constant from the {@link SmartsFragmenterReturnedFragmentsOption} enum
     * @throws NullPointerException if the given option is null
     */
    public void setReturnedFragmentsSetting(SmartsFragmenterReturnedFragmentsOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given option is null.");
        this.returnedFragmentsSetting.set(anOption);
    }

    /**
     * Sets the filter overlapping matches setting. If true, no atom is part of more than one returned match, and the
     * order of the patterns is their priority.
     *
     * @param aBoolean true if overlapping matches should be filtered
     */
    public void setFilterOverlappingMatchesSetting(boolean aBoolean) {
        this.filterOverlappingMatchesSetting.set(aBoolean);
    }

    /**
     * Sets the option for saturating free valences on returned fragment molecules.
     *
     * @param anOption the saturation option to use
     * @throws NullPointerException if the given option is null
     */
    public void setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given saturation option is null.");
        this.fragmentSaturationSetting.set(anOption);
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
        return this.settingNameTooltipTextMap;
    }

    @Override
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return this.settingNameDisplayNameMap;
    }

    @Override
    public String getFragmentationAlgorithmName() {
        return SmartsFragmenter.ALGORITHM_NAME;
    }

    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("SmartsFragmenter.displayName");
    }

    @Override
    public IMoleculeFragmenter copy() {
        //the copy compiles its own patterns when it needs them
        SmartsFragmenter tmpCopy = new SmartsFragmenter();
        tmpCopy.setPatternsSetting(this.patternsSetting.get());
        tmpCopy.setReturnedFragmentsSetting((SmartsFragmenterReturnedFragmentsOption) this.returnedFragmentsSetting.get());
        tmpCopy.setFilterOverlappingMatchesSetting(this.filterOverlappingMatchesSetting.get());
        tmpCopy.setFragmentSaturationSetting((IMoleculeFragmenter.FragmentSaturationOption) this.fragmentSaturationSetting.get());
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {
        this.patternsSetting.set(SmartsFragmenter.PATTERNS_SETTING_DEFAULT);
        this.returnedFragmentsSetting.set(SmartsFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT);
        this.filterOverlappingMatchesSetting.set(SmartsFragmenter.FILTER_OVERLAPPING_MATCHES_SETTING_DEFAULT);
        this.fragmentSaturationSetting.set(IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT);
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
        boolean tmpIsSaturated = this.fragmentSaturationSetting.get().equals(FragmentSaturationOption.HYDROGEN_SATURATION);
        boolean tmpIsOverlapFiltered = this.filterOverlappingMatchesSetting.get();
        List<IAtomContainer> tmpFragments = new ArrayList<>();
        try {
            List<CompiledPattern> tmpPatterns = this.getCompiledPatterns();
            IAtomContainer tmpPreparedMolecule = SmartsFragmenter.prepareMolecule(aMolecule);
            int[] tmpCanonicalRanks = SmartsFragmenter.getCanonicalRanks(tmpPreparedMolecule);
            //atoms that are part of a returned match
            boolean[] tmpIsAtomMatched = new boolean[tmpPreparedMolecule.getAtomCount()];
            for (CompiledPattern tmpPattern : tmpPatterns) {
                int[][] tmpMatches = tmpPattern.pattern().matchAll(tmpPreparedMolecule).uniqueAtoms().toArray();
                for (int[] tmpMatch : SmartsFragmenter.sortByCanonicalRanks(tmpMatches, tmpCanonicalRanks)) {
                    if (tmpIsOverlapFiltered && SmartsFragmenter.isAnyAtomFlagged(tmpMatch, tmpIsAtomMatched)) {
                        continue;
                    }
                    for (int tmpAtomIndex : tmpMatch) {
                        tmpIsAtomMatched[tmpAtomIndex] = true;
                    }
                    IAtomContainer tmpFragment = SmartsFragmenter.createFragment(tmpPreparedMolecule, tmpMatch, tmpIsSaturated);
                    tmpFragment.setProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY, tmpPattern.category());
                    tmpFragments.add(tmpFragment);
                }
            }
            if (this.returnedFragmentsSetting.get().equals(SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES_AND_REMAINDER)) {
                for (int[] tmpComponent : SmartsFragmenter.getRemainderComponents(tmpPreparedMolecule, tmpIsAtomMatched, tmpCanonicalRanks)) {
                    IAtomContainer tmpFragment = SmartsFragmenter.createFragment(tmpPreparedMolecule, tmpComponent, tmpIsSaturated);
                    tmpFragment.setProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY,
                            SmartsFragmenter.FRAGMENT_CATEGORY_REMAINDER_VALUE);
                    tmpFragments.add(tmpFragment);
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
        if (Objects.isNull(aMolecule) || aMolecule.isEmpty()) {
            return true;
        }
        if (!SmartsFragmenter.hasUnsetBondOrders(aMolecule)) {
            return false;
        }
        //unset bond orders of aromatic systems are resolved by kekulization, see prepareMolecule()
        try {
            IAtomContainer tmpMoleculeClone = aMolecule.clone();
            SmartsFragmenter.setUnsetImplicitHydrogenCountsToZero(tmpMoleculeClone);
            Kekulization.kekulize(tmpMoleculeClone);
            return SmartsFragmenter.hasUnsetBondOrders(tmpMoleculeClone);
        } catch (CDKException | CloneNotSupportedException | RuntimeException anException) {
            return true;
        }
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
                    SmartsFragmenter.LOGGER.log(Level.WARNING, anException.toString(), anException);
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
     * Returns the compiled SMARTS patterns of the current pattern list, in their order. The patterns are compiled
     * again if the pattern list has changed since the last call. Automatic preparation of the target molecule is
     * disabled because the molecule is prepared once for all patterns.
     *
     * @return compiled patterns
     */
    private List<CompiledPattern> getCompiledPatterns() {
        String tmpPatternListText = this.patternsSetting.get();
        if (Objects.isNull(this.compiledPatterns) || !tmpPatternListText.equals(this.compiledPatternsText)) {
            List<SimpleSmartsPatternListProperty.SmartsPatternEntry> tmpEntries =
                    SimpleSmartsPatternListProperty.parseEntries(tmpPatternListText);
            List<CompiledPattern> tmpCompiledPatterns = new ArrayList<>(tmpEntries.size());
            for (SimpleSmartsPatternListProperty.SmartsPatternEntry tmpEntry : tmpEntries) {
                SmartsPattern tmpPattern = SmartsPattern.create(tmpEntry.smarts(), SilentChemObjectBuilder.getInstance())
                        .setPrepare(false);
                tmpCompiledPatterns.add(new CompiledPattern(tmpPattern, tmpEntry.getLabel()));
            }
            this.compiledPatterns = tmpCompiledPatterns;
            this.compiledPatternsText = tmpPatternListText;
        }
        return this.compiledPatterns;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static methods">
    /**
     * Creates a fragment that consists of the given atoms and all bonds between them, see the class documentation.
     *
     * @param aMolecule the prepared molecule, see {@link #prepareMolecule(IAtomContainer)}
     * @param anAtomIndices indices of the atoms of the fragment
     * @param isSaturated whether free valences resulting from cut bonds are saturated with implicit hydrogen atoms
     * @return new atom container of the fragment
     * @throws CloneNotSupportedException if an atom cannot be cloned
     */
    private static IAtomContainer createFragment(IAtomContainer aMolecule, int[] anAtomIndices, boolean isSaturated)
            throws CloneNotSupportedException {
        int[] tmpAddedHydrogenCounts = new int[anAtomIndices.length];
        IAtomContainer tmpFragment = SmartsFragmenter.copyAtomsAndBonds(aMolecule, anAtomIndices, true, tmpAddedHydrogenCounts);
        try {
            Kekulization.kekulize(tmpFragment);
        } catch (CDKException anException) {
            //the aromatic part of the fragment has no Kekulé structure on its own, the bond orders of the molecule are used
            tmpFragment = SmartsFragmenter.copyAtomsAndBonds(aMolecule, anAtomIndices, false, tmpAddedHydrogenCounts);
        }
        //ring perception and Daylight aromaticity detection, as for the molecule prior to SMARTS matching
        SmartsPattern.prepare(tmpFragment);
        if (!isSaturated) {
            //bond orders and aromaticity were determined for the saturated fragment, so that both options agree in them
            for (int i = 0; i < anAtomIndices.length; i++) {
                IAtom tmpAtom = tmpFragment.getAtom(i);
                tmpAtom.setImplicitHydrogenCount(tmpAtom.getImplicitHydrogenCount() - tmpAddedHydrogenCounts[i]);
            }
        }
        SmartsFragmenter.addContainedStereoElements(aMolecule, anAtomIndices, tmpFragment);
        return tmpFragment;
    }

    /**
     * Returns a clone of the given molecule that is prepared for matching all SMARTS patterns: unset implicit hydrogen
     * counts are set to zero, unset bond orders of aromatic systems are resolved by kekulization, and ring atoms and
     * bonds are marked and the Daylight aromaticity model is applied (see {@link SmartsPattern#prepare(IAtomContainer)}).
     *
     * @param aMolecule the molecule to fragment, it is not altered
     * @return prepared clone of the given molecule
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     * @throws CDKException if unset bond orders cannot be resolved; such molecules should have been filtered
     */
    private static IAtomContainer prepareMolecule(IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        SmartsFragmenter.setUnsetImplicitHydrogenCountsToZero(tmpMoleculeClone);
        if (SmartsFragmenter.hasUnsetBondOrders(tmpMoleculeClone)) {
            Kekulization.kekulize(tmpMoleculeClone);
        }
        SmartsPattern.prepare(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Sets the implicit hydrogen counts of all atoms of the given molecule that have none to zero, as the importer does.
     *
     * @param aMolecule the molecule to alter
     */
    private static void setUnsetImplicitHydrogenCountsToZero(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom.getImplicitHydrogenCount())) {
                tmpAtom.setImplicitHydrogenCount(0);
            }
        }
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
     * Computes the canonical rank of every atom as its position in the unique SMILES code of the molecule. This
     * numbering does not depend on the input atom order, and neighbouring atoms get close ranks because the SMILES code
     * is written by a depth-first traversal. Should the SMILES code not be creatable, the ranks of the CDK canonical
     * labelling are used instead, which do not depend on the input atom order either.
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
            SmartsFragmenter.LOGGER.log(Level.FINE, "Unique SMILES order not available, canonical labels are used instead.", anException);
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

    /**
     * Sorts the given matches in the lexicographic order of the sorted canonical ranks of their atoms, so that the
     * matches are processed independently of the input atom order. The atom indices of the matches are not altered.
     *
     * @param aMatches matches, each given as atom indices
     * @param aCanonicalRanks canonical rank of every atom, indexed by atom index
     * @return the matches in canonical order
     */
    private static List<int[]> sortByCanonicalRanks(int[][] aMatches, int[] aCanonicalRanks) {
        int[][] tmpSortedRanks = new int[aMatches.length][];
        Integer[] tmpOrder = new Integer[aMatches.length];
        for (int i = 0; i < aMatches.length; i++) {
            int[] tmpRanks = new int[aMatches[i].length];
            for (int j = 0; j < aMatches[i].length; j++) {
                tmpRanks[j] = aCanonicalRanks[aMatches[i][j]];
            }
            Arrays.sort(tmpRanks);
            tmpSortedRanks[i] = tmpRanks;
            tmpOrder[i] = i;
        }
        Arrays.sort(tmpOrder, (a, b) -> Arrays.compare(tmpSortedRanks[a], tmpSortedRanks[b]));
        List<int[]> tmpSortedMatches = new ArrayList<>(aMatches.length);
        for (int tmpIndex : tmpOrder) {
            tmpSortedMatches.add(aMatches[tmpIndex]);
        }
        return tmpSortedMatches;
    }

    /**
     * Checks whether any atom of the given match is flagged.
     *
     * @param aMatch atom indices of the match
     * @param anAtomFlagArray flags, indexed by atom index
     * @return true if at least one atom of the match is flagged
     */
    private static boolean isAnyAtomFlagged(int[] aMatch, boolean[] anAtomFlagArray) {
        for (int tmpAtomIndex : aMatch) {
            if (anAtomFlagArray[tmpAtomIndex]) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the connected components of the atoms that are not part of any returned match, considering only the
     * bonds between these atoms. The components are ordered by the smallest canonical rank of their atoms.
     *
     * @param aMolecule the prepared molecule
     * @param anIsAtomMatchedArray flags of the atoms that are part of a returned match, indexed by atom index
     * @param aCanonicalRanks canonical rank of every atom, indexed by atom index
     * @return atom indices of every remainder component
     */
    private static List<int[]> getRemainderComponents(IAtomContainer aMolecule, boolean[] anIsAtomMatchedArray,
                                                      int[] aCanonicalRanks) {
        int[][] tmpAdjacencyList = GraphUtil.toAdjList(aMolecule);
        boolean[] tmpIsVisited = new boolean[tmpAdjacencyList.length];
        //start atoms in canonical order, so that the components are ordered by their smallest canonical rank
        Integer[] tmpStartAtoms = new Integer[tmpAdjacencyList.length];
        for (int i = 0; i < tmpStartAtoms.length; i++) {
            tmpStartAtoms[i] = i;
        }
        Arrays.sort(tmpStartAtoms, Comparator.comparingInt((Integer i) -> aCanonicalRanks[i]));
        List<int[]> tmpComponents = new ArrayList<>();
        for (int tmpStartAtom : tmpStartAtoms) {
            if (anIsAtomMatchedArray[tmpStartAtom] || tmpIsVisited[tmpStartAtom]) {
                continue;
            }
            List<Integer> tmpComponent = new ArrayList<>();
            Deque<Integer> tmpQueue = new ArrayDeque<>();
            tmpQueue.add(tmpStartAtom);
            tmpIsVisited[tmpStartAtom] = true;
            while (!tmpQueue.isEmpty()) {
                int tmpAtom = tmpQueue.poll();
                tmpComponent.add(tmpAtom);
                for (int tmpNeighbour : tmpAdjacencyList[tmpAtom]) {
                    if (!anIsAtomMatchedArray[tmpNeighbour] && !tmpIsVisited[tmpNeighbour]) {
                        tmpIsVisited[tmpNeighbour] = true;
                        tmpQueue.add(tmpNeighbour);
                    }
                }
            }
            tmpComponents.add(tmpComponent.stream().mapToInt(Integer::intValue).toArray());
        }
        return tmpComponents;
    }

    /**
     * Copies the given atoms and all bonds between them into a new atom container and saturates the free valences
     * resulting from bonds to other atoms with implicit hydrogen atoms. The number of hydrogen atoms added to every
     * atom is stored in the given array, so that they can be removed again. If the aromaticity is kept, the aromatic
     * flags of the atoms and bonds are copied, the orders of the aromatic bonds are left unset so that a Kekulé
     * structure can be assigned to them afterwards, and every cut aromatic bond counts as a single bond. Otherwise, the
     * bond orders of the given molecule are used and no aromatic flags are set.
     *
     * @param aMolecule the prepared molecule
     * @param anAtomIndices indices of the atoms to copy; the i-th atom of the new container is a copy of the atom with
     *                      the i-th index
     * @param isAromaticityKept whether the aromatic flags are copied and cut aromatic bonds count as single bonds
     * @param anAddedHydrogenCountArray array of the same length as the atom indices that is filled with the number of
     *                                  implicit hydrogen atoms added to each atom
     * @return new atom container of the given atoms and the bonds between them
     * @throws CloneNotSupportedException if an atom cannot be cloned
     */
    private static IAtomContainer copyAtomsAndBonds(IAtomContainer aMolecule, int[] anAtomIndices, boolean isAromaticityKept,
                                                    int[] anAddedHydrogenCountArray) throws CloneNotSupportedException {
        IChemObjectBuilder tmpBuilder = aMolecule.getBuilder();
        IAtomContainer tmpFragment = tmpBuilder.newAtomContainer();
        int[] tmpFragmentIndexOfAtom = new int[aMolecule.getAtomCount()];
        Arrays.fill(tmpFragmentIndexOfAtom, -1);
        for (int i = 0; i < anAtomIndices.length; i++) {
            IAtom tmpOriginalAtom = aMolecule.getAtom(anAtomIndices[i]);
            IAtom tmpFragmentAtom = tmpOriginalAtom.clone();
            tmpFragmentAtom.setImplicitHydrogenCount(tmpOriginalAtom.getImplicitHydrogenCount());
            tmpFragmentAtom.setIsAromatic(isAromaticityKept && tmpOriginalAtom.isAromatic());
            tmpFragmentAtom.setIsInRing(false);
            tmpFragment.addAtom(tmpFragmentAtom);
            tmpFragmentIndexOfAtom[anAtomIndices[i]] = i;
            anAddedHydrogenCountArray[i] = 0;
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            int tmpBegin = tmpFragmentIndexOfAtom[aMolecule.indexOf(tmpBond.getBegin())];
            int tmpEnd = tmpFragmentIndexOfAtom[aMolecule.indexOf(tmpBond.getEnd())];
            boolean tmpIsAromaticBondKept = isAromaticityKept && tmpBond.isAromatic();
            if (tmpBegin >= 0 && tmpEnd >= 0) {
                IBond tmpFragmentBond = tmpBuilder.newInstance(IBond.class, tmpFragment.getAtom(tmpBegin),
                        tmpFragment.getAtom(tmpEnd), tmpIsAromaticBondKept ? IBond.Order.UNSET : tmpBond.getOrder());
                tmpFragmentBond.setIsAromatic(tmpIsAromaticBondKept);
                tmpFragment.addBond(tmpFragmentBond);
            } else if (tmpBegin >= 0 || tmpEnd >= 0) {
                //bond to an atom outside the fragment is cut, saturate the free valence with hydrogen atoms
                int tmpFragmentAtomIndex = tmpBegin >= 0 ? tmpBegin : tmpEnd;
                int tmpHydrogenCount = tmpIsAromaticBondKept ? 1 : tmpBond.getOrder().numeric();
                IAtom tmpAtom = tmpFragment.getAtom(tmpFragmentAtomIndex);
                tmpAtom.setImplicitHydrogenCount(tmpAtom.getImplicitHydrogenCount() + tmpHydrogenCount);
                anAddedHydrogenCountArray[tmpFragmentAtomIndex] += tmpHydrogenCount;
            }
        }
        return tmpFragment;
    }

    /**
     * Adds those stereo elements of the given molecule to the fragment whose atoms and bonds are all part of the
     * fragment.
     *
     * @param aMolecule the prepared molecule
     * @param anAtomIndices indices of the atoms of the fragment; the i-th atom of the fragment is a copy of the atom
     *                      with the i-th index
     * @param aFragment the fragment
     */
    private static void addContainedStereoElements(IAtomContainer aMolecule, int[] anAtomIndices, IAtomContainer aFragment) {
        if (!aMolecule.stereoElements().iterator().hasNext()) {
            return;
        }
        Map<IAtom, IAtom> tmpAtomMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(anAtomIndices.length));
        for (int i = 0; i < anAtomIndices.length; i++) {
            tmpAtomMap.put(aMolecule.getAtom(anAtomIndices[i]), aFragment.getAtom(i));
        }
        Map<IBond, IBond> tmpBondMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(aFragment.getBondCount()));
        for (IBond tmpBond : aMolecule.bonds()) {
            IAtom tmpBegin = tmpAtomMap.get(tmpBond.getBegin());
            IAtom tmpEnd = tmpAtomMap.get(tmpBond.getEnd());
            if (!Objects.isNull(tmpBegin) && !Objects.isNull(tmpEnd)) {
                tmpBondMap.put(tmpBond, aFragment.getBond(tmpBegin, tmpEnd));
            }
        }
        for (IStereoElement tmpStereoElement : aMolecule.stereoElements()) {
            boolean tmpIsContained = SmartsFragmenter.isMapped(tmpStereoElement.getFocus(), tmpAtomMap, tmpBondMap);
            for (Object tmpCarrier : tmpStereoElement.getCarriers()) {
                tmpIsContained = tmpIsContained && SmartsFragmenter.isMapped(tmpCarrier, tmpAtomMap, tmpBondMap);
            }
            if (tmpIsContained) {
                aFragment.addStereoElement(tmpStereoElement.map(tmpAtomMap, tmpBondMap));
            }
        }
    }

    /**
     * Checks whether the given focus or carrier of a stereo element is an atom or bond that is part of the fragment.
     *
     * @param anObject focus or carrier of a stereo element
     * @param anAtomMap map of the molecule atoms to the fragment atoms
     * @param aBondMap map of the molecule bonds to the fragment bonds
     * @return true if the object is a mapped atom or bond
     */
    private static boolean isMapped(Object anObject, Map<IAtom, IAtom> anAtomMap, Map<IBond, IBond> aBondMap) {
        if (anObject instanceof IAtom tmpAtom) {
            return anAtomMap.containsKey(tmpAtom);
        }
        if (anObject instanceof IBond tmpBond) {
            return aBondMap.containsKey(tmpBond);
        }
        return false;
    }
    //</editor-fold>
}
