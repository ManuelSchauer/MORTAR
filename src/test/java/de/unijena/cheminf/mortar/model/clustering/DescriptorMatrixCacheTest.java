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
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

/**
 * Tests for {@link DescriptorMatrixCache} and for the clustering parameters with given descriptors.
 *
 * @author Manuel Schauer
 */
class DescriptorMatrixCacheTest {
    /**
     * Electron donation model of the tests.
     */
    private static final IMoleculeFragmenter.ElectronDonationModelOption ELECTRON_DONATION =
            IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT;
    /**
     * Cycle finder of the tests.
     */
    private static final IMoleculeFragmenter.CycleFinderOption CYCLE_FINDER = IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET;
    //
    /**
     * Constructor that sets the default locale to British English, which is needed for the message resources of the
     * fragmenter enums.
     */
    DescriptorMatrixCacheTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates test molecules.
     *
     * @return molecules
     */
    private static List<MoleculeDataModel> createMolecules() {
        return List.of(new MoleculeDataModel("c1ccccc1O", "phenol", new HashMap<>()),
                new MoleculeDataModel("CCO", "ethanol", new HashMap<>()),
                new MoleculeDataModel("CC(=O)O", "acetic acid", new HashMap<>()));
    }
    //
    /**
     * Tests that the cache returns the same values as a direct calculation, in the requested descriptor and molecule
     * order, and that the descriptors are cached afterwards.
     *
     * @throws Exception if the calculation fails
     */
    @Test
    void testCachedValuesEqualDirectCalculation() throws Exception {
        List<MoleculeDataModel> tmpMolecules = DescriptorMatrixCacheTest.createMolecules();
        Descriptor[] tmpDescriptors = {Descriptor.TPSA, Descriptor.ATOM_COUNT_HEAVY};
        FeatureMatrix tmpDirect = FeatureMatrixBuilder.buildDescriptorMatrix(tmpMolecules, tmpDescriptors,
                ChemUtil.createAromaticityModel(DescriptorMatrixCacheTest.ELECTRON_DONATION, DescriptorMatrixCacheTest.CYCLE_FINDER), false);
        DescriptorMatrixCache tmpCache = new DescriptorMatrixCache();
        FeatureMatrix tmpCached = tmpCache.getDescriptorMatrix(tmpMolecules, tmpDescriptors, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false);
        Assertions.assertArrayEquals(tmpDirect.getDataCopy(), tmpCached.getDataCopy());
        Assertions.assertTrue(tmpCache.isCached(tmpMolecules, Descriptor.TPSA, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
        //another order of descriptors and a subset of the molecules in another order come from the cache
        FeatureMatrix tmpSubset = tmpCache.getDescriptorMatrix(List.of(tmpMolecules.get(2), tmpMolecules.get(0)),
                new Descriptor[] {Descriptor.ATOM_COUNT_HEAVY, Descriptor.TPSA}, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false);
        Assertions.assertEquals("acetic acid", tmpSubset.getRowName(0));
        Assertions.assertEquals(tmpDirect.getValue(2, 1), tmpSubset.getValue(0, 0));
        Assertions.assertEquals(tmpDirect.getValue(0, 0), tmpSubset.getValue(1, 1));
        Assertions.assertEquals("TPSA", tmpSubset.getFeatureName(1));
        //the cache still holds all molecules
        Assertions.assertTrue(tmpCache.isCached(tmpMolecules, Descriptor.TPSA, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
    }
    //
    /**
     * Tests that missing descriptors are added, that another aromaticity model or a missing descriptor for a subset of
     * the molecules resets the cache and that clear empties it.
     *
     * @throws Exception if the calculation fails
     */
    @Test
    void testResetAndClear() throws Exception {
        List<MoleculeDataModel> tmpMolecules = DescriptorMatrixCacheTest.createMolecules();
        DescriptorMatrixCache tmpCache = new DescriptorMatrixCache();
        tmpCache.getDescriptorMatrix(tmpMolecules, new Descriptor[] {Descriptor.TPSA}, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false);
        FeatureMatrix tmpExtended = tmpCache.getDescriptorMatrix(tmpMolecules, new Descriptor[] {Descriptor.TPSA, Descriptor.ATOM_COUNT_HEAVY},
                DescriptorMatrixCacheTest.ELECTRON_DONATION, DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false);
        Assertions.assertEquals(2, tmpExtended.getColumnCount());
        Assertions.assertTrue(tmpCache.isCached(tmpMolecules, Descriptor.TPSA, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
        //another aromaticity model is not served from the cache
        IMoleculeFragmenter.ElectronDonationModelOption tmpOtherModel = Arrays.stream(IMoleculeFragmenter.ElectronDonationModelOption.values())
                .filter(anOption -> anOption != DescriptorMatrixCacheTest.ELECTRON_DONATION).findFirst().orElseThrow();
        Assertions.assertFalse(tmpCache.isCached(tmpMolecules, Descriptor.TPSA, tmpOtherModel, DescriptorMatrixCacheTest.CYCLE_FINDER));
        //a new descriptor for a subset resets the cache to the subset
        List<MoleculeDataModel> tmpSubset = tmpMolecules.subList(0, 2);
        tmpCache.getDescriptorMatrix(tmpSubset, new Descriptor[] {Descriptor.X_LOG_P}, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false);
        Assertions.assertTrue(tmpCache.isCached(tmpSubset, Descriptor.X_LOG_P, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
        Assertions.assertFalse(tmpCache.isCached(tmpMolecules, Descriptor.TPSA, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
        tmpCache.clear();
        Assertions.assertFalse(tmpCache.isCached(tmpSubset, Descriptor.X_LOG_P, DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER));
        //no descriptors: an empty matrix with one row per molecule
        Assertions.assertEquals(3, tmpCache.getDescriptorMatrix(tmpMolecules, new Descriptor[0], DescriptorMatrixCacheTest.ELECTRON_DONATION,
                DescriptorMatrixCacheTest.CYCLE_FINDER, 1024, false).getRowCount());
    }
    //
    /**
     * Tests that a descriptor calculation of the service fills the cache used by its clustering tasks and that the
     * clustering parameters with given descriptors contain exactly these descriptors and fingerprints.
     *
     * @throws Exception if the calculation fails
     */
    @Test
    void testServiceCalculationAndParametersWithGivenDescriptors() throws Exception {
        ClusteringService tmpService = new ClusteringService();
        List<MoleculeDataModel> tmpMolecules = DescriptorMatrixCacheTest.createMolecules();
        tmpService.calculateDescriptorMatrix(tmpMolecules, new Descriptor[] {Descriptor.TPSA}, false);
        Assertions.assertTrue(tmpService.getDescriptorMatrixCache().isCached(tmpMolecules, Descriptor.TPSA,
                ClusteringService.ELECTRON_DONATION_MODEL_OPTION_DEFAULT, ClusteringService.CYCLE_FINDER_OPTION_DEFAULT));
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.MOLECULAR_WEIGHT, true);
        ClusteringParameters tmpParameters = tmpService.createClusteringParameters("Frag",
                new Descriptor[] {Descriptor.X_LOG_P, Descriptor.TPSA, Descriptor.MACCS_FINGERPRINTER}, false);
        //exactly the given descriptors and fingerprints, e.g. the features calculated in a descriptor tab
        Assertions.assertEquals(List.of(Descriptor.X_LOG_P, Descriptor.TPSA, Descriptor.MACCS_FINGERPRINTER), Arrays.asList(tmpParameters.descriptors()));
        tmpService.clearDescriptorMatrixCache();
        Assertions.assertFalse(tmpService.getDescriptorMatrixCache().isCached(tmpMolecules, Descriptor.TPSA,
                ClusteringService.ELECTRON_DONATION_MODEL_OPTION_DEFAULT, ClusteringService.CYCLE_FINDER_OPTION_DEFAULT));
    }
}
