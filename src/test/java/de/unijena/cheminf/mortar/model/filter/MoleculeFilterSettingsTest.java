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

import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.preference.BooleanPreference;
import de.unijena.cheminf.mortar.preference.PreferenceContainer;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;
import de.unijena.cheminf.mortar.preference.SingleIntegerPreference;
import de.unijena.cheminf.mortar.preference.SingleTermPreference;

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Tests the {@link MoleculeFilterSettings} container and the {@link MoleculeFilterSettingsSnapshot}: neutral defaults,
 * the settings view contract (properties, display names, tooltips, preference restrictions), the handling of the heavy
 * atom count range, snapshots, and the persistence round trip into a temporary directory (the real MORTAR settings
 * directory is not touched). All tests run headless.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class MoleculeFilterSettingsTest {
    //<editor-fold desc="Static initializer">
    /**
     * Sets the default locale to British English, required for the message texts of the settings.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Tests">
    /**
     * Tests that a new settings container holds the neutral defaults (all filters inactive) and that its snapshot
     * equals the default snapshot.
     */
    @Test
    public void testDefaultsAreNeutral() {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        Assertions.assertEquals(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED,
                tmpSettings.getDisconnectedStructuresSetting());
        Assertions.assertFalse(tmpSettings.getElementFilterSetting());
        Assertions.assertFalse(tmpSettings.getHeavyAtomCountFilterSetting());
        Assertions.assertFalse(tmpSettings.getDuplicateFilterSetting());
        Assertions.assertEquals(MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                tmpSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertEquals(MoleculeFilterSettings.MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                tmpSettings.getMaximumHeavyAtomCountSetting());
        MoleculeFilterSettingsSnapshot tmpSnapshot = tmpSettings.createSnapshot();
        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(), tmpSnapshot);
        Assertions.assertFalse(tmpSnapshot.isAnyFilterActive());
        Assertions.assertFalse(tmpSnapshot.isAnyStructuralFilterActive());
    }

    /**
     * Tests the contract of the generic settings view: six properties of renderable types with unique, persistable
     * names, a display name and a tooltip for each, and display names and tooltips for all enum constants. Missing
     * message keys would show up as "Key '...' not found.".
     */
    @Test
    public void testSettingsViewContract() {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        List<Property<?>> tmpProperties = tmpSettings.settingsProperties();
        Assertions.assertEquals(MoleculeFilterSettings.NUMBER_OF_SETTINGS, tmpProperties.size());
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpProperties));
        Set<String> tmpNames = new HashSet<>();
        for (Property<?> tmpProperty : tmpProperties) {
            Assertions.assertTrue(tmpNames.add(tmpProperty.getName()));
            Assertions.assertNotNull(tmpProperty.getValue());
            this.assertMessageText(tmpSettings.getSettingNameToDisplayNameMap().get(tmpProperty.getName()));
            this.assertMessageText(tmpSettings.getSettingNameToTooltipTextMap().get(tmpProperty.getName()));
        }
        this.assertMessageText(tmpSettings.getDisplayName());
        //order of application
        Assertions.assertSame(tmpSettings.disconnectedStructuresSettingProperty(), tmpProperties.get(0));
        Assertions.assertSame(tmpSettings.elementFilterSettingProperty(), tmpProperties.get(1));
        Assertions.assertSame(tmpSettings.heavyAtomCountFilterSettingProperty(), tmpProperties.get(2));
        Assertions.assertSame(tmpSettings.minimumHeavyAtomCountSettingProperty(), tmpProperties.get(3));
        Assertions.assertSame(tmpSettings.maximumHeavyAtomCountSettingProperty(), tmpProperties.get(4));
        Assertions.assertSame(tmpSettings.duplicateFilterSettingProperty(), tmpProperties.get(5));
        for (IDisplayEnum tmpConstant : MoleculeFilterSettings.DisconnectedStructuresOption.values()) {
            this.assertMessageText(tmpConstant.getDisplayName());
            this.assertMessageText(tmpConstant.getTooltipText());
        }
        for (IDisplayEnum tmpConstant : MoleculeFilterReason.values()) {
            this.assertMessageText(tmpConstant.getDisplayName());
            this.assertMessageText(tmpConstant.getTooltipText());
        }
        Assertions.assertTrue(tmpSettings.getSettingNameToTooltipTextMap()
                .get(tmpSettings.elementFilterSettingProperty().getName())
                .contains("H, B, C, N, O, F, Si, P, S, Cl, Se, Br, and I"));
    }

    /**
     * Tests the setters and that invalid values are rejected with an IllegalArgumentException without any GUI alert
     * (which would fail with a different exception in this headless test), also when called from a non-GUI thread.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testSettersAndInvalidValuesHeadless() throws Exception {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        tmpSettings.setDisconnectedStructuresSetting(MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE);
        tmpSettings.setElementFilterSetting(true);
        tmpSettings.setHeavyAtomCountFilterSetting(true);
        tmpSettings.setDuplicateFilterSetting(true);
        Assertions.assertEquals(MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE,
                tmpSettings.getDisconnectedStructuresSetting());
        Assertions.assertTrue(tmpSettings.getElementFilterSetting());
        Assertions.assertTrue(tmpSettings.getHeavyAtomCountFilterSetting());
        Assertions.assertTrue(tmpSettings.getDuplicateFilterSetting());
        Assertions.assertTrue(tmpSettings.createSnapshot().isAnyFilterActive());
        Assertions.assertThrows(NullPointerException.class, () -> tmpSettings.setDisconnectedStructuresSetting(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpSettings.disconnectedStructuresSettingProperty().set(null));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpSettings.setMinimumHeavyAtomCountSetting(-1));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpSettings.maximumHeavyAtomCountSettingProperty().setValue(-5));
        Assertions.assertEquals(MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                tmpSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertEquals(MoleculeFilterSettings.MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                tmpSettings.getMaximumHeavyAtomCountSetting());
        AtomicReference<Throwable> tmpThrown = new AtomicReference<>();
        Thread tmpThread = new Thread(() -> {
            try {
                tmpSettings.setMaximumHeavyAtomCountSetting(-1);
            } catch (Throwable aThrowable) {
                tmpThrown.set(aThrowable);
            }
            tmpSettings.setMaximumHeavyAtomCountSetting(50);
            tmpSettings.restoreDefaultSettings();
        });
        tmpThread.start();
        tmpThread.join();
        Assertions.assertInstanceOf(IllegalArgumentException.class, tmpThrown.get());
        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(), tmpSettings.createSnapshot());
    }

    /**
     * Tests that the heavy atom count range can never become empty: setting the minimum above the maximum raises the
     * maximum and vice versa, and setting any valid pair of bounds in any order from any start state gives exactly
     * this pair (as needed for restoring recent values in the settings view and for reloading persisted settings).
     */
    @Test
    public void testHeavyAtomCountRangeNeverInvalid() {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        tmpSettings.setMinimumHeavyAtomCountSetting(10);
        tmpSettings.setMaximumHeavyAtomCountSetting(20);
        tmpSettings.setMinimumHeavyAtomCountSetting(30);
        Assertions.assertEquals(30, tmpSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertEquals(30, tmpSettings.getMaximumHeavyAtomCountSetting());
        tmpSettings.setMaximumHeavyAtomCountSetting(5);
        Assertions.assertEquals(5, tmpSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertEquals(5, tmpSettings.getMaximumHeavyAtomCountSetting());
        //the settings view sets 0 if a text field is cleared
        tmpSettings.minimumHeavyAtomCountSettingProperty().setValue(null);
        Assertions.assertEquals(0, tmpSettings.getMinimumHeavyAtomCountSetting());
        int[][] tmpPairs = {{0, 0}, {0, 1}, {2, 100}, {50, 60}, {150, 300}, {7, 7}, {1000, 100000}};
        for (int[] tmpStart : tmpPairs) {
            for (int[] tmpTarget : tmpPairs) {
                for (boolean tmpIsMinimumFirst : new boolean[] {true, false}) {
                    MoleculeFilterSettings tmpTestSettings = new MoleculeFilterSettings();
                    tmpTestSettings.setMaximumHeavyAtomCountSetting(Integer.MAX_VALUE);
                    tmpTestSettings.setMinimumHeavyAtomCountSetting(tmpStart[0]);
                    tmpTestSettings.setMaximumHeavyAtomCountSetting(tmpStart[1]);
                    if (tmpIsMinimumFirst) {
                        tmpTestSettings.setMinimumHeavyAtomCountSetting(tmpTarget[0]);
                        tmpTestSettings.setMaximumHeavyAtomCountSetting(tmpTarget[1]);
                    } else {
                        tmpTestSettings.setMaximumHeavyAtomCountSetting(tmpTarget[1]);
                        tmpTestSettings.setMinimumHeavyAtomCountSetting(tmpTarget[0]);
                    }
                    Assertions.assertEquals(tmpTarget[0], tmpTestSettings.getMinimumHeavyAtomCountSetting());
                    Assertions.assertEquals(tmpTarget[1], tmpTestSettings.getMaximumHeavyAtomCountSetting());
                    Assertions.assertDoesNotThrow(tmpTestSettings::createSnapshot);
                }
            }
        }
    }

    /**
     * Tests restoring the defaults and that a snapshot is not affected by later settings changes.
     */
    @Test
    public void testRestoreDefaultsAndSnapshotIsolation() {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        tmpSettings.setDisconnectedStructuresSetting(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT);
        tmpSettings.setElementFilterSetting(true);
        tmpSettings.setHeavyAtomCountFilterSetting(true);
        tmpSettings.setMinimumHeavyAtomCountSetting(500);
        tmpSettings.setDuplicateFilterSetting(true);
        MoleculeFilterSettingsSnapshot tmpSnapshot = tmpSettings.createSnapshot();
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(tmpSettings);
        tmpSettings.restoreDefaultSettings();
        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(), tmpSettings.createSnapshot());
        MoleculeFilterSettingsSnapshot tmpExpected = new MoleculeFilterSettingsSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, true, 500, 500, true);
        Assertions.assertEquals(tmpExpected, tmpSnapshot);
        Assertions.assertEquals(tmpExpected, tmpRun.getSettingsSnapshot());
    }

    /**
     * Tests the argument checks of the snapshot record.
     */
    @Test
    public void testSnapshotValidation() {
        Assertions.assertThrows(NullPointerException.class, () -> new MoleculeFilterSettingsSnapshot(null,
                false, false, 1, 2, false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> new MoleculeFilterSettingsSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, -1, 2, false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> new MoleculeFilterSettingsSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 3, 2, false));
        MoleculeFilterSettingsSnapshot tmpOnlyDuplicates = new MoleculeFilterSettingsSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 3, 3, true);
        Assertions.assertTrue(tmpOnlyDuplicates.isAnyFilterActive());
        Assertions.assertFalse(tmpOnlyDuplicates.isAnyStructuralFilterActive());
    }

    /**
     * Tests the persistence round trip (write, change, reload) into a temporary directory instead of the MORTAR
     * settings directory, including a bound pair that only reloads correctly because of the range handling, a missing
     * file, and a file with invalid values.
     *
     * @param aTempDir temporary directory provided by JUnit
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testPersistenceRoundTrip(@TempDir Path aTempDir) throws Exception {
        String tmpDirectoryPath = aTempDir.toString();
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        //nothing persisted yet, settings remain unchanged
        Assertions.assertFalse(tmpSettings.reloadSettings(tmpDirectoryPath));
        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(), tmpSettings.createSnapshot());
        tmpSettings.setDisconnectedStructuresSetting(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT);
        tmpSettings.setElementFilterSetting(true);
        tmpSettings.setHeavyAtomCountFilterSetting(true);
        //both bounds above the default maximum: reloading into a default container must raise the maximum first
        tmpSettings.setMaximumHeavyAtomCountSetting(300);
        tmpSettings.setMinimumHeavyAtomCountSetting(150);
        tmpSettings.setDuplicateFilterSetting(true);
        MoleculeFilterSettingsSnapshot tmpPersistedSnapshot = tmpSettings.createSnapshot();
        Assertions.assertTrue(tmpSettings.preserveSettings(tmpDirectoryPath));
        File tmpFile = new File(MoleculeFilterSettings.getSettingsFilePath(tmpDirectoryPath));
        Assertions.assertTrue(tmpFile.isFile());
        Assertions.assertEquals(MoleculeFilterSettings.SETTINGS_FILE_NAME + ".txt", tmpFile.getName());
        //change and reload
        tmpSettings.restoreDefaultSettings();
        tmpSettings.setMaximumHeavyAtomCountSetting(3);
        Assertions.assertTrue(tmpSettings.reloadSettings(tmpDirectoryPath));
        Assertions.assertEquals(tmpPersistedSnapshot, tmpSettings.createSnapshot());
        //fresh container
        MoleculeFilterSettings tmpReloadedSettings = new MoleculeFilterSettings();
        Assertions.assertTrue(tmpReloadedSettings.reloadSettings(tmpDirectoryPath));
        Assertions.assertEquals(tmpPersistedSnapshot, tmpReloadedSettings.createSnapshot());
        //lower pair into a container with a higher range
        tmpReloadedSettings.setMinimumHeavyAtomCountSetting(0);
        tmpReloadedSettings.setMaximumHeavyAtomCountSetting(1);
        Assertions.assertTrue(tmpReloadedSettings.preserveSettings(tmpDirectoryPath));
        Assertions.assertTrue(tmpSettings.reloadSettings(tmpDirectoryPath));
        Assertions.assertEquals(0, tmpSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertEquals(1, tmpSettings.getMaximumHeavyAtomCountSetting());
        //a file with invalid values: they are skipped, the valid ones are applied, nothing is thrown
        PreferenceContainer tmpInvalidContainer = new PreferenceContainer(MoleculeFilterSettings.getSettingsFilePath(tmpDirectoryPath));
        tmpInvalidContainer.add(new SingleTermPreference(tmpSettings.disconnectedStructuresSettingProperty().getName(), "NO_SUCH_OPTION"));
        tmpInvalidContainer.add(new SingleIntegerPreference(tmpSettings.minimumHeavyAtomCountSettingProperty().getName(), -7));
        tmpInvalidContainer.add(new BooleanPreference(tmpSettings.elementFilterSettingProperty().getName(), false));
        tmpInvalidContainer.writeRepresentation();
        MoleculeFilterSettings tmpThirdSettings = new MoleculeFilterSettings();
        tmpThirdSettings.setElementFilterSetting(true);
        Assertions.assertTrue(tmpThirdSettings.reloadSettings(tmpDirectoryPath));
        Assertions.assertEquals(MoleculeFilterSettings.DISCONNECTED_STRUCTURES_SETTING_DEFAULT,
                tmpThirdSettings.getDisconnectedStructuresSetting());
        Assertions.assertEquals(MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                tmpThirdSettings.getMinimumHeavyAtomCountSetting());
        Assertions.assertFalse(tmpThirdSettings.getElementFilterSetting());
        //a path that cannot be a directory
        Assertions.assertFalse(tmpSettings.preserveSettings(tmpFile.getPath()));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Asserts that the given message text exists, i.e. is neither null, blank, nor the placeholder that is returned
     * for missing message keys.
     *
     * @param aText the text to check
     */
    private void assertMessageText(String aText) {
        Assertions.assertNotNull(aText);
        Assertions.assertFalse(aText.isBlank());
        Assertions.assertFalse(aText.contains("not found"), aText);
    }
    //</editor-fold>
}
