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

package de.unijena.cheminf.mortar.model.clustering;

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Tests for {@link ClusteringService}, {@link ClusteringTask#runPipeline} and {@link ClusteringRunResult}.
 *
 * @author Manuel Schauer
 */
class ClusteringServiceTest {
    /**
     * Name of the test fragmentation.
     */
    private static final String FRAGMENTATION_NAME = "TestFragmentation";
    //
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    ClusteringServiceTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests that all settings have valid preference names (needed for persistence), a display name and a tooltip, and
     * that the defaults are MORTAR's standard aromaticity model and the cosine kernel.
     */
    @Test
    void testSettingsAndDefaults() {
        ClusteringService tmpService = new ClusteringService();
        List<Property<?>> tmpAllSettings = new ArrayList<>(tmpService.settingsProperties());
        tmpAllSettings.addAll(tmpService.descriptorSettingsProperties());
        PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpAllSettings);
        for (Property<?> tmpSetting : tmpAllSettings) {
            Assertions.assertTrue(PreferenceUtil.isValidName(tmpSetting.getName()), tmpSetting.getName());
            Assertions.assertNotNull(tmpService.getSettingNameToDisplayNameMap().get(tmpSetting.getName()), tmpSetting.getName());
            Assertions.assertNotNull(tmpService.getSettingNameToTooltipTextMap().get(tmpSetting.getName()), tmpSetting.getName());
        }
        //7 preset and fingerprint settings plus one setting per descriptor that is no fingerprint
        long tmpNumberOfDescriptors = java.util.Arrays.stream(Descriptor.getAllDescriptors()).filter(aDescriptor -> !aDescriptor.isFingerprint()).count();
        Assertions.assertEquals(tmpNumberOfDescriptors + 7, tmpService.descriptorSettingsProperties().size());
        ClusteringParameters tmpParameters = tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false);
        Assertions.assertEquals(IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, tmpParameters.electronDonationModelOption());
        Assertions.assertEquals(IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET, tmpParameters.cycleFinderOption());
        Assertions.assertEquals(Art2aClusterer.Art2aAlgorithm.COSINE, tmpParameters.algorithm());
        Assertions.assertEquals(Art2aClusterer.CONVERGENCE_THRESHOLD_COSINE_DEFAULT, tmpParameters.convergenceThreshold());
        for (Descriptor tmpDescriptor : tmpParameters.descriptors()) {
            Assertions.assertTrue(ClusteringService.isDescriptorSelectedByDefault(tmpDescriptor));
            Assertions.assertFalse(tmpDescriptor.isFingerprint());
        }
    }
    //
    /**
     * Tests that changing the algorithm sets the convergence threshold default of the kernel and that restoring the
     * defaults works.
     */
    @Test
    void testAlgorithmChangeSetsConvergenceDefault() {
        ClusteringService tmpService = new ClusteringService();
        SimpleDoubleProperty tmpConvergence = (SimpleDoubleProperty) ClusteringServiceTest.getSetting(tmpService, "Convergence threshold setting");
        Property tmpAlgorithm = ClusteringServiceTest.getSetting(tmpService, "Clustering algorithm setting");
        tmpAlgorithm.setValue(ClusteringService.AlgorithmOption.ART_2A_EUCLID);
        Assertions.assertEquals(Art2aClusterer.CONVERGENCE_THRESHOLD_EUCLID_DEFAULT, tmpConvergence.get(), 1.0e-6);
        ((SimpleIntegerProperty) ClusteringServiceTest.getSetting(tmpService, "Random seed setting")).set(7);
        tmpService.restoreDefaultSettings();
        ClusteringParameters tmpParameters = tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false);
        Assertions.assertEquals(Art2aClusterer.Art2aAlgorithm.COSINE, tmpParameters.algorithm());
        Assertions.assertEquals(Art2aClusterer.CONVERGENCE_THRESHOLD_COSINE_DEFAULT, tmpParameters.convergenceThreshold());
        Assertions.assertEquals(Art2aClusterer.RANDOM_SEED_DEFAULT, tmpParameters.randomSeed());
    }
    //
    /**
     * Tests that the float defaults of the library are shown as exact decimal values in the settings (0.99, not
     * 0.9900000095367432) and still result in the library defaults when passed to the kernel.
     */
    @Test
    void testDefaultsAreExactDecimals() {
        ClusteringService tmpService = new ClusteringService();
        Assertions.assertEquals(0.99, ((SimpleDoubleProperty) ClusteringServiceTest.getSetting(tmpService, "Convergence threshold setting")).get());
        Assertions.assertEquals(0.01, ((SimpleDoubleProperty) ClusteringServiceTest.getSetting(tmpService, "Learning parameter setting")).get());
        Assertions.assertEquals(1.0, ((SimpleDoubleProperty) ClusteringServiceTest.getSetting(tmpService, "Offset for contrast enhancement setting")).get());
        Property tmpAlgorithm = ClusteringServiceTest.getSetting(tmpService, "Clustering algorithm setting");
        tmpAlgorithm.setValue(ClusteringService.AlgorithmOption.ART_2A_EUCLID);
        Assertions.assertEquals(0.1, ((SimpleDoubleProperty) ClusteringServiceTest.getSetting(tmpService, "Convergence threshold setting")).get());
        ClusteringParameters tmpParameters = tmpService.createClusteringParameters("Test", false);
        Assertions.assertEquals(Art2aClusterer.CONVERGENCE_THRESHOLD_EUCLID_DEFAULT, tmpParameters.convergenceThreshold());
        Assertions.assertEquals(Art2aClusterer.LEARNING_PARAMETER_DEFAULT, tmpParameters.learningParameter());
    }
    //
    /**
     * Tests that the correlation filter is switched off by default, that the setting is passed to the parameters, and
     * that the pipeline only removes correlated features if it is switched on.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testCorrelationFilterSwitch() throws Exception {
        ClusteringService tmpService = new ClusteringService();
        Assertions.assertFalse(tmpService.correlationFilterProperty().get());
        Assertions.assertFalse(tmpService.createClusteringParameters("Test", false).isCorrelationFilterEnabled());
        Assertions.assertTrue(tmpService.filterSettingsProperties().contains(tmpService.correlationFilterProperty()));
        //two perfectly correlated descriptors (heavy atom count and carbon count of hydrocarbons)
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        String[] tmpSmiles = {"C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"};
        for (int i = 0; i < tmpSmiles.length; i++) {
            //the fragment counts alternate and are not correlated with the descriptors
            tmpMolecules.add(ClusteringServiceTest.createMolecule(tmpSmiles[i], "m" + i, Map.of("C", 1 + i % 2)));
        }
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.ATOM_COUNT_HEAVY, true);
        tmpService.setDescriptorSelected(Descriptor.ATOM_COUNT_C, true);
        ((SimpleIntegerProperty) ClusteringServiceTest.getSetting(tmpService, "Minimum fragment molecule frequency setting")).set(1);
        ClusteringRunResult tmpWithoutFilter = ClusteringTask.runPipeline(tmpMolecules,
                tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false),
                tmpService.getPriorityTable(), () -> false, (aStep, aMessage) -> { });
        Assertions.assertTrue(tmpWithoutFilter.removedFeatures().stream()
                .noneMatch(aFeature -> aFeature.reason() == RemovedFeature.RemovalReason.CORRELATION));
        tmpService.correlationFilterProperty().set(true);
        ClusteringRunResult tmpWithFilter = ClusteringTask.runPipeline(tmpMolecules,
                tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false),
                tmpService.getPriorityTable(), () -> false, (aStep, aMessage) -> { });
        Assertions.assertTrue(tmpWithFilter.removedFeatures().stream()
                .anyMatch(aFeature -> aFeature.reason() == RemovedFeature.RemovalReason.CORRELATION));
        Assertions.assertTrue(tmpWithFilter.processedMatrix().getColumnCount() < tmpWithoutFilter.processedMatrix().getColumnCount());
    }
    //
    /**
     * Tests the descriptor selection.
     */
    @Test
    void testDescriptorSelection() {
        ClusteringService tmpService = new ClusteringService();
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.ATOM_COUNT_C, true);
        Assertions.assertArrayEquals(new Descriptor[] {Descriptor.ATOM_COUNT_C}, tmpService.getSelectedDescriptors());
        tmpService.restoreDefaultDescriptorSelection();
        Assertions.assertTrue(tmpService.getSelectedDescriptors().length > 1);
    }
    //
    /**
     * Tests that a preset sets the check boxes, that a deviating manual change switches the selection to custom and that
     * restoring the defaults selects the default preset again.
     */
    @Test
    void testDescriptorPresets() {
        ClusteringService tmpService = new ClusteringService();
        tmpService.selectDescriptorPreset(ClusteringService.DescriptorPresetOption.BASIC.name());
        Assertions.assertEquals(ClusteringService.DescriptorPresetOption.BASIC_DESCRIPTORS,
                java.util.Set.of(tmpService.getSelectedDescriptors()));
        //selecting a descriptor of the preset again keeps the preset
        tmpService.setDescriptorSelected(Descriptor.TPSA, true);
        Assertions.assertEquals(ClusteringService.DescriptorPresetOption.BASIC.name(), tmpService.descriptorPresetKeyProperty().get());
        //a deviation switches to custom and keeps the check boxes
        tmpService.setDescriptorSelected(Descriptor.CHI_PATH, true);
        Assertions.assertEquals(ClusteringService.CUSTOM_PRESET_KEY, tmpService.descriptorPresetKeyProperty().get());
        Assertions.assertEquals(ClusteringService.DescriptorPresetOption.BASIC_DESCRIPTORS.size() + 1, tmpService.getSelectedDescriptors().length);
        Assertions.assertEquals("Basic (custom)", tmpService.getProposedUserPresetName());
        tmpService.descriptorPresetKeyProperty().set(ClusteringService.DescriptorPresetOption.ALL_2D.name());
        for (Descriptor tmpDescriptor : tmpService.getSelectedDescriptors()) {
            Assertions.assertFalse(tmpDescriptor.requires3DCoordinates());
        }
        tmpService.restoreDefaultDescriptorSelection();
        Assertions.assertEquals(ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.name(), tmpService.descriptorPresetKeyProperty().get());
        for (Descriptor tmpDescriptor : tmpService.getSelectedDescriptors()) {
            Assertions.assertTrue(tmpDescriptor.isFast() && tmpDescriptor.isSafe() && !tmpDescriptor.requires3DCoordinates());
        }
    }
    //
    /**
     * Tests saving a modified selection as own preset: it is listed first, survives a reload in a new service instance
     * (including the selection of the preset) and can be deleted. Invalid names are rejected.
     *
     * @param aTempDirectory temporary settings directory provided by JUnit
     * @throws Exception if the presets cannot be written
     */
    @Test
    void testUserDescriptorPresets(@org.junit.jupiter.api.io.TempDir java.nio.file.Path aTempDirectory) throws Exception {
        ClusteringService tmpService = new ClusteringService(aTempDirectory.toString());
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.TPSA, true);
        tmpService.setDescriptorSelected(Descriptor.X_LOG_P, true);
        String tmpKey = tmpService.saveCurrentDescriptorSelectionAsPreset("My polarity set");
        Assertions.assertEquals(tmpKey, tmpService.descriptorPresetKeyProperty().get());
        DescriptorPreset tmpFirstPreset = tmpService.getDescriptorPresets().getFirst();
        Assertions.assertTrue(tmpFirstPreset.isUserDefined());
        Assertions.assertEquals(java.util.Set.of(Descriptor.TPSA, Descriptor.X_LOG_P), tmpFirstPreset.descriptors());
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpService.saveCurrentDescriptorSelectionAsPreset("Fast"));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpService.saveCurrentDescriptorSelectionAsPreset("  "));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpService.saveCurrentDescriptorSelectionAsPreset("a=b"));
        //a new service reloads the preset and its selection
        tmpService.persistClusteringSettings();
        ClusteringService tmpReloadedService = new ClusteringService(aTempDirectory.toString());
        tmpReloadedService.reloadClusteringSettings();
        Assertions.assertEquals(tmpKey, tmpReloadedService.descriptorPresetKeyProperty().get());
        Assertions.assertArrayEquals(new Descriptor[] {Descriptor.TPSA, Descriptor.X_LOG_P}, tmpReloadedService.getSelectedDescriptors());
        tmpReloadedService.deleteUserDescriptorPreset(tmpKey);
        Assertions.assertEquals(ClusteringService.CUSTOM_PRESET_KEY, tmpReloadedService.descriptorPresetKeyProperty().get());
        Assertions.assertFalse(tmpReloadedService.getDescriptorPresets().getFirst().isUserDefined());
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpReloadedService.deleteUserDescriptorPreset(ClusteringService.DescriptorPresetOption.FAST.name()));
    }
    //
    /**
     * Tests that the fingerprints are selected via their own settings and that type and radius are mapped to the
     * correct circular fingerprint (radius 2 = ECFP4).
     */
    @Test
    void testFingerprintSettings() {
        ClusteringService tmpService = new ClusteringService();
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        ClusteringServiceTest.getDescriptorSetting(tmpService, "Circular fingerprint setting").setValue(true);
        ClusteringServiceTest.getDescriptorSetting(tmpService, "Circular fingerprint type setting")
                .setValue(ClusteringService.CircularFingerprintTypeOption.FCFP);
        ClusteringServiceTest.getDescriptorSetting(tmpService, "MACCS fingerprint setting").setValue(true);
        Assertions.assertArrayEquals(new Descriptor[] {Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4, Descriptor.MACCS_FINGERPRINTER},
                tmpService.getSelectedDescriptors());
        for (ClusteringService.CircularFingerprintTypeOption tmpType : ClusteringService.CircularFingerprintTypeOption.values()) {
            for (ClusteringService.CircularFingerprintRadiusOption tmpRadius : ClusteringService.CircularFingerprintRadiusOption.values()) {
                Descriptor tmpDescriptor = ClusteringService.getCircularFingerprintDescriptor(tmpType, tmpRadius);
                Assertions.assertTrue(tmpDescriptor.isFingerprint());
                Assertions.assertTrue(tmpDescriptor.name().endsWith(tmpType.name() + "_" + 2 * tmpRadius.getRadius()));
            }
        }
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpService.setDescriptorSelected(Descriptor.MACCS_FINGERPRINTER, true));
        Assertions.assertEquals("Chi Path (slow)", ClusteringService.getDescriptorDisplayName(Descriptor.CHI_PATH));
        Assertions.assertEquals("TPSA", ClusteringService.getDescriptorDisplayName(Descriptor.TPSA));
    }
    //
    /**
     * Returns the descriptor tab setting with the given name.
     *
     * @param aService clustering service
     * @param aName setting name
     * @return setting
     */
    private static Property getDescriptorSetting(ClusteringService aService, String aName) {
        return aService.descriptorSettingsProperties().stream().filter(aSetting -> aSetting.getName().equals(aName)).findFirst().orElseThrow();
    }
    //
    /**
     * Tests the complete pipeline on two groups of molecules: aromatic molecules sharing the benzene fragment and
     * aliphatic alcohols sharing the hydroxy fragment. The groups must end up in different clusters, the summaries must
     * list the representative first and the matrices of all stages must be available.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testRunPipelineSeparatesTwoGroups() throws Exception {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        String[] tmpAromatics = {"c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1Cl", "c1ccccc1Br"};
        String[] tmpAlcohols = {"CCO", "CCCO", "CCCCO", "CC(C)O", "CCC(C)O"};
        for (int i = 0; i < tmpAromatics.length; i++) {
            tmpMolecules.add(ClusteringServiceTest.createMolecule(tmpAromatics[i], "aromatic" + i, Map.of("c1ccccc1", 1, "C", 1 + i % 2)));
            tmpMolecules.add(ClusteringServiceTest.createMolecule(tmpAlcohols[i], "alcohol" + i, Map.of("O", 1, "CC", 1 + i % 2)));
        }
        ClusteringService tmpService = new ClusteringService();
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.AROMATIC_ATOMS_COUNT, true);
        tmpService.setDescriptorSelected(Descriptor.ATOM_COUNT_O, true);
        //vigilance 0.5 with the default fragment frequency filter of 2 molecules
        ClusteringParameters tmpParameters = tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false);
        List<Integer> tmpSteps = new ArrayList<>();
        ClusteringRunResult tmpRunResult = ClusteringTask.runPipeline(tmpMolecules, tmpParameters, tmpService.getPriorityTable(),
                () -> false, (aStep, aMessage) -> tmpSteps.add(aStep));
        Assertions.assertEquals(List.of(0, 1, 2, 3, 4, 5), tmpSteps);
        //4 fragments and 2 descriptors
        Assertions.assertEquals(6, tmpRunResult.rawMatrix().getColumnCount());
        //the descriptors follow the fragments in enum order; toluene has 6 aromatic atoms
        FeatureMatrix tmpRawMatrix = tmpRunResult.rawMatrix();
        int tmpAromaticAtomsColumn = -1;
        for (int j = 0; j < tmpRawMatrix.getColumnCount(); j++) {
            if (tmpRawMatrix.getFeatureGroupName(j).equals(Descriptor.AROMATIC_ATOMS_COUNT.name())) {
                tmpAromaticAtomsColumn = j;
            }
        }
        Assertions.assertEquals(5, tmpAromaticAtomsColumn);
        Assertions.assertEquals(6.0f, tmpRawMatrix.getValue(0, tmpAromaticAtomsColumn));
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        for (int i = 0; i < tmpMolecules.size(); i += 2) {
            //aromatic and alcohol molecules are never in the same cluster
            for (int j = 1; j < tmpMolecules.size(); j += 2) {
                Assertions.assertNotEquals(tmpLabels[i], tmpLabels[j]);
            }
        }
        List<ClusteringRunResult.ClusterSummary> tmpSummaries = tmpRunResult.getClusterSummaries();
        int tmpTotal = 0;
        for (ClusteringRunResult.ClusterSummary tmpSummary : tmpSummaries) {
            tmpTotal += tmpSummary.getSize();
            if (tmpSummary.clusterIndex() >= 0) {
                Assertions.assertSame(tmpSummary.representative(), tmpSummary.members().getFirst());
            }
        }
        Assertions.assertEquals(tmpMolecules.size(), tmpTotal);
    }
    //
    /**
     * Tests that the pipeline stops if it is cancelled.
     */
    @Test
    void testRunPipelineCancellation() {
        List<MoleculeDataModel> tmpMolecules = List.of(ClusteringServiceTest.createMolecule("CCO", "m1", Map.of("O", 1)));
        ClusteringService tmpService = new ClusteringService();
        Assertions.assertThrows(java.util.concurrent.CancellationException.class, () -> ClusteringTask.runPipeline(
                tmpMolecules, tmpService.createClusteringParameters(ClusteringServiceTest.FRAGMENTATION_NAME, false),
                tmpService.getPriorityTable(), () -> true, (aStep, aMessage) -> { }));
    }
    //
    /**
     * Returns the general setting with the given name.
     *
     * @param aService clustering service
     * @param aName setting name
     * @return setting
     */
    private static Property<?> getSetting(ClusteringService aService, String aName) {
        return aService.settingsProperties().stream().filter(aSetting -> aSetting.getName().equals(aName)).findFirst().orElseThrow();
    }
    //
    /**
     * Creates a molecule data model with the given fragment frequencies for the test fragmentation.
     *
     * @param aSmiles SMILES of the molecule
     * @param aName name of the molecule
     * @param aFrequencies fragment frequencies
     * @return molecule data model
     */
    private static MoleculeDataModel createMolecule(String aSmiles, String aName, Map<String, Integer> aFrequencies) {
        MoleculeDataModel tmpMolecule = new MoleculeDataModel(aSmiles, aName, new HashMap<>());
        tmpMolecule.getFragmentFrequencies().put(ClusteringServiceTest.FRAGMENTATION_NAME, aFrequencies);
        return tmpMolecule;
    }
}
