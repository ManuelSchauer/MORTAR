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

package de.unijena.cheminf.mortar.model.filter;

import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.FileUtil;
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.model.util.SimpleIDisplayEnumConstantProperty;
import de.unijena.cheminf.mortar.preference.PreferenceContainer;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Container for the settings of the molecule filters that are applied at import (disconnected structures, elements,
 * heavy atom count, duplicates). Externally, the settings can be accessed via JavaFX properties in the same way as the
 * global settings ({@link de.unijena.cheminf.mortar.model.settings.SettingsContainer}) and the fragmenter settings,
 * so that the generic settings view can render them from {@link #settingsProperties()},
 * {@link #getSettingNameToDisplayNameMap()}, and {@link #getSettingNameToTooltipTextMap()}. For persistence, they are
 * translated into {@link de.unijena.cheminf.mortar.preference.IPreference} objects and written to an own file in the
 * settings directory, separate from the global settings file.
 * <br>All default values are neutral, i.e. with them, an import behaves exactly as without molecule filters.
 * <br>Thread safety and GUI independence: in contrast to the global settings container, this class never raises GUI
 * alerts, neither in the setters nor at persistence or reload; invalid values are logged and rejected with an
 * IllegalArgumentException (re-thrown so that a bidirectional binding is reset), and the persistence methods report
 * failures by their return value. Therefore, all methods can be called from non-GUI threads and in headless tests.
 * The setters of the two heavy atom count bounds, {@link #restoreDefaultSettings()}, and {@link #createSnapshot()}
 * synchronise on an internal lock, so a snapshot taken by an import thread always contains a consistent range. Apart
 * from that, the JavaFX properties are not thread-safe and should only be changed from one thread at a time.
 * <br>The invalid case "minimum heavy atom count larger than maximum heavy atom count" cannot occur: if the minimum is
 * set above the current maximum, the maximum is raised to the same value, and if the maximum is set below the current
 * minimum, the minimum is lowered to the same value. Thereby, every non-negative value is accepted in any order of
 * setting the two bounds (e.g. when restoring recent values in the settings view or reloading persisted settings), and
 * setting a valid pair of bounds in any order always yields exactly that pair.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class MoleculeFilterSettings {
    //<editor-fold desc="Enum DisconnectedStructuresOption">
    /**
     * Enum for the options how molecules consisting of more than one unconnected component (salts, hydrates, solvates,
     * mixtures) are treated at import.
     */
    public enum DisconnectedStructuresOption implements IDisplayEnum {
        /**
         * Molecules with several components are imported unchanged (behaviour of MORTAR without molecule filters).
         */
        KEEP_UNCHANGED(
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.KeepUnchanged.displayName"),
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.KeepUnchanged.tooltip")),
        /**
         * Molecules with several components are reduced to their largest component, determined by heavy atom count;
         * see {@link MoleculeFilterUtil#extractLargestComponent(org.openscience.cdk.interfaces.IAtomContainer)} for
         * the tie-break.
         */
        KEEP_LARGEST_COMPONENT(
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.KeepLargestComponent.displayName"),
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.KeepLargestComponent.tooltip")),
        /**
         * Molecules with several components are removed from the data set.
         */
        REMOVE_MOLECULE(
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.RemoveMolecule.displayName"),
                Message.get("MoleculeFilterSettings.DisconnectedStructuresOption.RemoveMolecule.tooltip"));
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
        private DisconnectedStructuresOption(String aDisplayName, String aTooltip) {
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
     * Number of different settings wrapped by the container.
     */
    public static final int NUMBER_OF_SETTINGS = 6;

    /**
     * Name of the file (without extension) in the settings directory that persists the molecule filter settings.
     */
    public static final String SETTINGS_FILE_NAME = "MORTAR_Molecule_Filter_Settings";

    /**
     * Default value of the disconnected structures setting, molecules are kept unchanged.
     */
    public static final DisconnectedStructuresOption DISCONNECTED_STRUCTURES_SETTING_DEFAULT =
            DisconnectedStructuresOption.KEEP_UNCHANGED;

    /**
     * Default value of whether the element filter is active.
     */
    public static final boolean ELEMENT_FILTER_SETTING_DEFAULT = false;

    /**
     * Default value of whether the heavy atom count filter is active.
     */
    public static final boolean HEAVY_ATOM_COUNT_FILTER_SETTING_DEFAULT = false;

    /**
     * Default value of the inclusive minimum heavy atom count; only used if the heavy atom count filter is activated.
     * With 2, molecules consisting of a single heavy atom (lone ions, water, methane) are removed.
     */
    public static final int MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT = 2;

    /**
     * Default value of the inclusive maximum heavy atom count; only used if the heavy atom count filter is activated.
     */
    public static final int MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT = 100;

    /**
     * Default value of whether the duplicate filter is active.
     */
    public static final boolean DUPLICATE_FILTER_SETTING_DEFAULT = false;
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(MoleculeFilterSettings.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Lock object for the setters of the heavy atom count bounds, restoring the defaults, and taking snapshots.
     */
    private final Object settingsLock;

    //note: no Javadoc comments on the private properties because the javadoc build complains about "double comments"
    // for the properties and their get() methods
    private final SimpleIDisplayEnumConstantProperty disconnectedStructuresSetting;

    private final SimpleBooleanProperty elementFilterSetting;

    private final SimpleBooleanProperty heavyAtomCountFilterSetting;

    private final SimpleIntegerProperty minimumHeavyAtomCountSetting;

    private final SimpleIntegerProperty maximumHeavyAtomCountSetting;

    private final SimpleBooleanProperty duplicateFilterSetting;

    /**
     * List of all settings for display in the settings view, in the order in which the filters are applied.
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
     * Constructor that sets all settings to their default values. Persisted settings are not loaded here, call
     * {@link #reloadSettings()} for that.
     */
    public MoleculeFilterSettings() {
        this.settingsLock = new Object();
        int tmpInitialCapacityForSettingNameMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                MoleculeFilterSettings.NUMBER_OF_SETTINGS,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingNameMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        //note: these names are for persistence and de-persistence! The map values are for display
        this.disconnectedStructuresSetting = new SimpleIDisplayEnumConstantProperty(this,
                "Disconnected structures setting",
                MoleculeFilterSettings.DISCONNECTED_STRUCTURES_SETTING_DEFAULT,
                MoleculeFilterSettings.DisconnectedStructuresOption.class) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    //call to super.set() for parameter checks
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    MoleculeFilterSettings.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    //no GUI alert here to keep this class usable from non-GUI threads;
                    // re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.disconnectedStructuresSetting.getName(),
                Message.get("MoleculeFilterSettings.disconnectedStructuresSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.disconnectedStructuresSetting.getName(),
                Message.get("MoleculeFilterSettings.disconnectedStructuresSetting.displayName"));
        this.elementFilterSetting = new SimpleBooleanProperty(this,
                "Element filter setting",
                MoleculeFilterSettings.ELEMENT_FILTER_SETTING_DEFAULT);
        this.settingNameTooltipTextMap.put(this.elementFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.elementFilterSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.elementFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.elementFilterSetting.displayName"));
        this.heavyAtomCountFilterSetting = new SimpleBooleanProperty(this,
                "Heavy atom count filter setting",
                MoleculeFilterSettings.HEAVY_ATOM_COUNT_FILTER_SETTING_DEFAULT);
        this.settingNameTooltipTextMap.put(this.heavyAtomCountFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.heavyAtomCountFilterSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.heavyAtomCountFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.heavyAtomCountFilterSetting.displayName"));
        this.minimumHeavyAtomCountSetting = new SimpleIntegerProperty(this,
                "Minimum heavy atom count setting",
                MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT) {
            @Override
            public void set(int newValue) throws IllegalArgumentException {
                MoleculeFilterSettings.this.checkHeavyAtomCountBound(newValue, this.getName());
                synchronized (MoleculeFilterSettings.this.settingsLock) {
                    SimpleIntegerProperty tmpMaximumSetting = MoleculeFilterSettings.this.maximumHeavyAtomCountSetting;
                    //raises the maximum first so that the range is never empty, also not for concurrent readers
                    if (Objects.nonNull(tmpMaximumSetting) && newValue > tmpMaximumSetting.get()) {
                        tmpMaximumSetting.set(newValue);
                    }
                    super.set(newValue);
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.minimumHeavyAtomCountSetting.getName(),
                Message.get("MoleculeFilterSettings.minimumHeavyAtomCountSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.minimumHeavyAtomCountSetting.getName(),
                Message.get("MoleculeFilterSettings.minimumHeavyAtomCountSetting.displayName"));
        this.maximumHeavyAtomCountSetting = new SimpleIntegerProperty(this,
                "Maximum heavy atom count setting",
                MoleculeFilterSettings.MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT) {
            @Override
            public void set(int newValue) throws IllegalArgumentException {
                MoleculeFilterSettings.this.checkHeavyAtomCountBound(newValue, this.getName());
                synchronized (MoleculeFilterSettings.this.settingsLock) {
                    SimpleIntegerProperty tmpMinimumSetting = MoleculeFilterSettings.this.minimumHeavyAtomCountSetting;
                    //lowers the minimum first so that the range is never empty, also not for concurrent readers
                    if (Objects.nonNull(tmpMinimumSetting) && newValue < tmpMinimumSetting.get()) {
                        tmpMinimumSetting.set(newValue);
                    }
                    super.set(newValue);
                }
            }
        };
        this.settingNameTooltipTextMap.put(this.maximumHeavyAtomCountSetting.getName(),
                Message.get("MoleculeFilterSettings.maximumHeavyAtomCountSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.maximumHeavyAtomCountSetting.getName(),
                Message.get("MoleculeFilterSettings.maximumHeavyAtomCountSetting.displayName"));
        this.duplicateFilterSetting = new SimpleBooleanProperty(this,
                "Duplicate filter setting",
                MoleculeFilterSettings.DUPLICATE_FILTER_SETTING_DEFAULT);
        this.settingNameTooltipTextMap.put(this.duplicateFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.duplicateFilterSetting.tooltip"));
        this.settingNameDisplayNameMap.put(this.duplicateFilterSetting.getName(),
                Message.get("MoleculeFilterSettings.duplicateFilterSetting.displayName"));
        //order of the filters' application
        this.settings = new ArrayList<>(MoleculeFilterSettings.NUMBER_OF_SETTINGS);
        this.settings.add(this.disconnectedStructuresSetting);
        this.settings.add(this.elementFilterSetting);
        this.settings.add(this.heavyAtomCountFilterSetting);
        this.settings.add(this.minimumHeavyAtomCountSetting);
        this.settings.add(this.maximumHeavyAtomCountSetting);
        this.settings.add(this.duplicateFilterSetting);
        try {
            PreferenceUtil.checkPropertiesForPreferenceRestrictions(this.settings);
        } catch (UnsupportedOperationException anException) {
            //can only happen if a developer introduces an invalid setting name; logged, no GUI alert (see class doc)
            MoleculeFilterSettings.LOGGER.log(Level.SEVERE, anException.toString(), anException);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties">
    /**
     * Returns a list of all settings as JavaFX properties, in the order in which the filters are applied, for display
     * in the settings view.
     *
     * @return list of settings as properties
     */
    public List<Property<?>> settingsProperties() {
        return this.settings;
    }

    /**
     * Returns a map containing descriptive texts (values) for the settings with the given names (keys) to be used as
     * tooltips in the GUI.
     *
     * @return map with tooltip texts
     */
    public Map<String, String> getSettingNameToTooltipTextMap() {
        return this.settingNameTooltipTextMap;
    }

    /**
     * Returns a map containing language-specific names for display (values) for the settings with the given names
     * (keys) to be used in the GUI.
     *
     * @return map with display names
     */
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return this.settingNameDisplayNameMap;
    }

    /**
     * Returns the language-specific display name of the molecule filter settings as a whole, e.g. as title of a
     * settings tab.
     *
     * @return display name
     */
    public String getDisplayName() {
        return Message.get("MoleculeFilterSettings.displayName");
    }

    /**
     * Returns the currently set option for the treatment of disconnected structures.
     *
     * @return disconnected structures setting value
     */
    public MoleculeFilterSettings.DisconnectedStructuresOption getDisconnectedStructuresSetting() {
        return (MoleculeFilterSettings.DisconnectedStructuresOption) this.disconnectedStructuresSetting.get();
    }

    /**
     * Returns the property wrapping the disconnected structures setting.
     *
     * @return disconnected structures setting property
     */
    public SimpleIDisplayEnumConstantProperty disconnectedStructuresSettingProperty() {
        return this.disconnectedStructuresSetting;
    }

    /**
     * Returns whether the element filter is active.
     *
     * @return element filter setting value
     */
    public boolean getElementFilterSetting() {
        return this.elementFilterSetting.get();
    }

    /**
     * Returns the property wrapping the element filter setting.
     *
     * @return element filter setting property
     */
    public SimpleBooleanProperty elementFilterSettingProperty() {
        return this.elementFilterSetting;
    }

    /**
     * Returns whether the heavy atom count filter is active.
     *
     * @return heavy atom count filter setting value
     */
    public boolean getHeavyAtomCountFilterSetting() {
        return this.heavyAtomCountFilterSetting.get();
    }

    /**
     * Returns the property wrapping the heavy atom count filter setting.
     *
     * @return heavy atom count filter setting property
     */
    public SimpleBooleanProperty heavyAtomCountFilterSettingProperty() {
        return this.heavyAtomCountFilterSetting;
    }

    /**
     * Returns the inclusive minimum heavy atom count.
     *
     * @return minimum heavy atom count setting value
     */
    public int getMinimumHeavyAtomCountSetting() {
        return this.minimumHeavyAtomCountSetting.get();
    }

    /**
     * Returns the property wrapping the minimum heavy atom count setting.
     *
     * @return minimum heavy atom count setting property
     */
    public SimpleIntegerProperty minimumHeavyAtomCountSettingProperty() {
        return this.minimumHeavyAtomCountSetting;
    }

    /**
     * Returns the inclusive maximum heavy atom count.
     *
     * @return maximum heavy atom count setting value
     */
    public int getMaximumHeavyAtomCountSetting() {
        return this.maximumHeavyAtomCountSetting.get();
    }

    /**
     * Returns the property wrapping the maximum heavy atom count setting.
     *
     * @return maximum heavy atom count setting property
     */
    public SimpleIntegerProperty maximumHeavyAtomCountSettingProperty() {
        return this.maximumHeavyAtomCountSetting;
    }

    /**
     * Returns whether the duplicate filter is active.
     *
     * @return duplicate filter setting value
     */
    public boolean getDuplicateFilterSetting() {
        return this.duplicateFilterSetting.get();
    }

    /**
     * Returns the property wrapping the duplicate filter setting.
     *
     * @return duplicate filter setting property
     */
    public SimpleBooleanProperty duplicateFilterSettingProperty() {
        return this.duplicateFilterSetting;
    }

    /**
     * Sets the option for the treatment of disconnected structures.
     *
     * @param anOption a constant of the disconnected structures option enum
     * @throws NullPointerException if the given option is null
     */
    public void setDisconnectedStructuresSetting(MoleculeFilterSettings.DisconnectedStructuresOption anOption) throws NullPointerException {
        Objects.requireNonNull(anOption, "Given disconnected structures option is null.");
        this.disconnectedStructuresSetting.set(anOption);
    }

    /**
     * Sets whether the element filter is active.
     *
     * @param isActive whether only molecules consisting of the allowed elements should be imported
     */
    public void setElementFilterSetting(boolean isActive) {
        this.elementFilterSetting.set(isActive);
    }

    /**
     * Sets whether the heavy atom count filter is active.
     *
     * @param isActive whether only molecules with a heavy atom count in the set range should be imported
     */
    public void setHeavyAtomCountFilterSetting(boolean isActive) {
        this.heavyAtomCountFilterSetting.set(isActive);
    }

    /**
     * Sets the inclusive minimum heavy atom count. If the given value is larger than the current maximum, the maximum
     * is raised to the same value.
     *
     * @param aCount the minimum heavy atom count
     * @throws IllegalArgumentException if the given count is negative
     */
    public void setMinimumHeavyAtomCountSetting(int aCount) throws IllegalArgumentException {
        this.minimumHeavyAtomCountSetting.set(aCount);
    }

    /**
     * Sets the inclusive maximum heavy atom count. If the given value is smaller than the current minimum, the minimum
     * is lowered to the same value.
     *
     * @param aCount the maximum heavy atom count
     * @throws IllegalArgumentException if the given count is negative
     */
    public void setMaximumHeavyAtomCountSetting(int aCount) throws IllegalArgumentException {
        this.maximumHeavyAtomCountSetting.set(aCount);
    }

    /**
     * Sets whether the duplicate filter is active.
     *
     * @param isActive whether molecules with a unique SMILES string that occurred before should be removed
     */
    public void setDuplicateFilterSetting(boolean isActive) {
        this.duplicateFilterSetting.set(isActive);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Restores all settings to their default values according to the respective public constants in this class. With
     * the defaults, all filters are inactive.
     */
    public void restoreDefaultSettings() {
        synchronized (this.settingsLock) {
            this.disconnectedStructuresSetting.set(MoleculeFilterSettings.DISCONNECTED_STRUCTURES_SETTING_DEFAULT);
            this.elementFilterSetting.set(MoleculeFilterSettings.ELEMENT_FILTER_SETTING_DEFAULT);
            this.heavyAtomCountFilterSetting.set(MoleculeFilterSettings.HEAVY_ATOM_COUNT_FILTER_SETTING_DEFAULT);
            //order does not matter, see the class doc
            this.minimumHeavyAtomCountSetting.set(MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT);
            this.maximumHeavyAtomCountSetting.set(MoleculeFilterSettings.MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT);
            this.duplicateFilterSetting.set(MoleculeFilterSettings.DUPLICATE_FILTER_SETTING_DEFAULT);
        }
    }

    /**
     * Takes an immutable snapshot of the current settings. To be called once at the start of an import run, so that
     * settings changes during the run do not affect it.
     *
     * @return snapshot of the current settings
     */
    public MoleculeFilterSettingsSnapshot createSnapshot() {
        synchronized (this.settingsLock) {
            return new MoleculeFilterSettingsSnapshot(
                    this.getDisconnectedStructuresSetting(),
                    this.elementFilterSetting.get(),
                    this.heavyAtomCountFilterSetting.get(),
                    this.minimumHeavyAtomCountSetting.get(),
                    this.maximumHeavyAtomCountSetting.get(),
                    this.duplicateFilterSetting.get());
        }
    }

    /**
     * Persists the current settings in the file {@link #SETTINGS_FILE_NAME} in the MORTAR settings directory
     * ({@link FileUtil#getSettingsDirPath()}) for reloading at the next application session. Failures are logged and
     * reported by the return value, no GUI alert is raised.
     *
     * @return true if the settings were written successfully
     */
    public boolean preserveSettings() {
        String tmpSettingsDirectoryPath;
        try {
            tmpSettingsDirectoryPath = FileUtil.getSettingsDirPath();
        } catch (SecurityException anException) {
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, String.format("Molecule filter settings persistence went wrong, " +
                    "settings directory not available: %s", anException.toString()), anException);
            return false;
        }
        return this.preserveSettings(tmpSettingsDirectoryPath);
    }

    /**
     * Persists the current settings in the file {@link #SETTINGS_FILE_NAME} in the given directory, which is created
     * if it does not exist. Intended for the application settings directory (see {@link #preserveSettings()}) and for
     * tests that must not touch it. Failures are logged and reported by the return value, no GUI alert is raised.
     *
     * @param aDirectoryPath path of the directory to write the settings file to
     * @return true if the settings were written successfully
     * @throws NullPointerException if the given path is null
     */
    public boolean preserveSettings(String aDirectoryPath) throws NullPointerException {
        Objects.requireNonNull(aDirectoryPath, "Given directory path is null.");
        File tmpDirectory = new File(aDirectoryPath);
        boolean tmpMKDirsSuccessful = true;
        if (!tmpDirectory.exists()) {
            tmpMKDirsSuccessful = tmpDirectory.mkdirs();
        }
        if (!tmpMKDirsSuccessful || !tmpDirectory.isDirectory() || !tmpDirectory.canWrite()) {
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, "Molecule filter settings persistence went wrong, " +
                    "cannot write to directory {0}.", aDirectoryPath);
            return false;
        }
        String tmpFilePath = MoleculeFilterSettings.getSettingsFilePath(aDirectoryPath);
        try {
            PreferenceContainer tmpPreferenceContainer = PreferenceUtil.translateJavaFxPropertiesToPreferences(
                    this.settings, tmpFilePath);
            tmpPreferenceContainer.writeRepresentation();
            return true;
        } catch (NullPointerException | IllegalArgumentException | IOException | SecurityException anException) {
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, String.format("Molecule filter settings persistence " +
                    "went wrong, exception: %s", anException.toString()), anException);
            return false;
        }
    }

    /**
     * Reloads the settings from the file {@link #SETTINGS_FILE_NAME} in the MORTAR settings directory
     * ({@link FileUtil#getSettingsDirPath()}). Settings that cannot be found or read remain in their current value
     * (the default, if called right after construction). Failures are logged and reported by the return value, no GUI
     * alert is raised.
     *
     * @return true if the settings file was found and read; false if it does not exist or could not be read
     */
    public boolean reloadSettings() {
        String tmpSettingsDirectoryPath;
        try {
            tmpSettingsDirectoryPath = FileUtil.getSettingsDirPath();
        } catch (SecurityException anException) {
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, String.format("Unable to reload molecule filter settings, " +
                    "settings directory not available: %s", anException.toString()), anException);
            return false;
        }
        return this.reloadSettings(tmpSettingsDirectoryPath);
    }

    /**
     * Reloads the settings from the file {@link #SETTINGS_FILE_NAME} in the given directory. Settings that cannot be
     * found or read remain in their current value. Failures are logged and reported by the return value, no GUI alert
     * is raised.
     *
     * @param aDirectoryPath path of the directory to read the settings file from
     * @return true if the settings file was found and read; false if it does not exist or could not be read
     * @throws NullPointerException if the given path is null
     */
    public boolean reloadSettings(String aDirectoryPath) throws NullPointerException {
        Objects.requireNonNull(aDirectoryPath, "Given directory path is null.");
        File tmpSettingsFile = new File(MoleculeFilterSettings.getSettingsFilePath(aDirectoryPath));
        if (!tmpSettingsFile.exists() || !tmpSettingsFile.isFile() || !tmpSettingsFile.canRead()) {
            MoleculeFilterSettings.LOGGER.log(Level.INFO, "No persisted molecule filter settings could be found, " +
                    "they remain in their current values (defaults at start-up).");
            return false;
        }
        PreferenceContainer tmpPreferenceContainer;
        try {
            tmpPreferenceContainer = new PreferenceContainer(tmpSettingsFile);
        } catch (IllegalArgumentException | IOException | SecurityException anException) {
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, String.format("Unable to reload molecule filter settings: %s",
                    anException.toString()), anException);
            return false;
        }
        //invalid single values are logged and skipped by the utility method
        PreferenceUtil.updatePropertiesFromPreferences(this.settings, tmpPreferenceContainer);
        return true;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Returns the path of the settings file in the given directory.
     *
     * @param aDirectoryPath path of the directory, with or without trailing file separator
     * @return path of the settings file
     * @throws NullPointerException if the given path is null
     */
    public static String getSettingsFilePath(String aDirectoryPath) throws NullPointerException {
        Objects.requireNonNull(aDirectoryPath, "Given directory path is null.");
        return new File(aDirectoryPath, MoleculeFilterSettings.SETTINGS_FILE_NAME
                + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION).getPath();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Checks whether the given value is a legal heavy atom count bound, i.e. not negative, and throws an
     * IllegalArgumentException after logging it if not.
     *
     * @param aCount the value to check
     * @param aSettingName name of the setting to report
     * @throws IllegalArgumentException if the given value is negative
     */
    private void checkHeavyAtomCountBound(int aCount, String aSettingName) throws IllegalArgumentException {
        if (aCount < 0) {
            IllegalArgumentException tmpException = new IllegalArgumentException(String.format(
                    "An illegal value was given for the %s: %d; it must not be negative.", aSettingName, aCount));
            MoleculeFilterSettings.LOGGER.log(Level.WARNING, tmpException.toString(), tmpException);
            //no GUI alert here to keep this class usable from non-GUI threads (the settings view only accepts
            // non-negative integers anyway); re-throws the exception to properly reset the binding
            throw tmpException;
        }
    }
    //</editor-fold>
}
