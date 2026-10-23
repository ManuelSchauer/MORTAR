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
import de.unijena.cheminf.mortar.model.io.ClusteringDataCsvWriter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

/**
 * Tests for {@link ClusterProfiler}.
 *
 * @author Manuel Schauer
 */
class ClusterProfilerTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    ClusterProfilerTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the profile values against hand-calculated values: four molecules, two clusters, one fragment column, one
     * fingerprint bit column and one descriptor column with a missing value.
     */
    @Test
    void testProfileValues() {
        //fragment "F1": rows 0, 1, 2; bit "MACCS [1]": row 3; descriptor "D": 1, 3, NaN, 8
        float[][] tmpData = {{1, 0, 1}, {2, 0, 3}, {1, 0, Float.NaN}, {0, 1, 8}};
        FeatureMatrix tmpMatrix = new FeatureMatrix(tmpData, new String[] {"m0", "m1", "m2", "m3"},
                new String[] {"F1", "MACCS [1]", "D"},
                new FeatureType[] {FeatureType.COUNT, FeatureType.BIT, FeatureType.CONTINUOUS},
                new String[] {"Frag", Descriptor.MACCS_FINGERPRINTER.name(), Descriptor.TPSA.name()});
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        for (int i = 0; i < 4; i++) {
            tmpMolecules.add(new MoleculeDataModel("C", "m" + i, new HashMap<>()));
        }
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", true, 1, 0, new Descriptor[0],
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        //cluster 0: rows 0, 1, 2; cluster 1: row 3
        ClusteringRunResult tmpRunResult = new ClusteringRunResult(tmpParameters, tmpMolecules, tmpMatrix, tmpMatrix, tmpMatrix,
                List.of(), new ClusteringResult(new int[] {0, 0, 0, 1}, 2, new int[] {0, 3}, false, true, 2, 0.5f, new float[2][3]));
        ClusterProfiler.ClusterProfile tmpProfile = new ClusterProfiler(tmpRunResult).getProfile(0);
        Assertions.assertEquals(3, tmpProfile.size());
        ClusterProfiler.FeatureProfile tmpFragment = tmpProfile.features().get(0);
        Assertions.assertEquals(ClusterProfiler.FeatureCategory.FRAGMENT, tmpFragment.category());
        //present in 3 of 3 cluster molecules and 3 of 4 molecules overall
        Assertions.assertEquals(1.0, tmpFragment.clusterValue(), 1.0e-9);
        Assertions.assertEquals(0.75, tmpFragment.overallValue(), 1.0e-9);
        Assertions.assertEquals(0.25, tmpFragment.score(), 1.0e-9);
        ClusterProfiler.FeatureProfile tmpBit = tmpProfile.features().get(1);
        Assertions.assertEquals(ClusterProfiler.FeatureCategory.FINGERPRINT_BIT, tmpBit.category());
        Assertions.assertEquals(-0.25, tmpBit.score(), 1.0e-9);
        ClusterProfiler.FeatureProfile tmpDescriptor = tmpProfile.features().get(2);
        Assertions.assertEquals(ClusterProfiler.FeatureCategory.DESCRIPTOR, tmpDescriptor.category());
        //NaN is ignored: cluster mean (1 + 3) / 2 = 2, overall mean (1 + 3 + 8) / 3 = 4, SD sqrt((9 + 1 + 16) / 3)
        Assertions.assertEquals(2.0, tmpDescriptor.clusterValue(), 1.0e-9);
        Assertions.assertEquals(4.0, tmpDescriptor.overallValue(), 1.0e-9);
        Assertions.assertEquals(-2.0 / Math.sqrt(26.0 / 3.0), tmpDescriptor.score(), 1.0e-9);
        Assertions.assertEquals(List.of(tmpFragment), tmpProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, 5, true));
        Assertions.assertTrue(tmpProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, 5, false).isEmpty());
        Assertions.assertThrows(IllegalArgumentException.class, () -> new ClusterProfiler(tmpRunResult).getProfile(2));
    }
    //
    /**
     * Tests on the aromatic/alcohol data set that the benzene fragment and the aromatic atom count characterise the
     * cluster of the aromatic molecules and the hydroxy fragment the cluster of the alcohols; and that the profiles can
     * be written as CSV.
     *
     * @throws Exception if the clustering or writing fails
     */
    @Test
    void testProfilesOfAromaticAndAlcoholClusters() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        ClusterProfiler tmpProfiler = new ClusterProfiler(tmpRunResult);
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        //row 0 is aromatic, row 1 an alcohol
        ClusterProfiler.ClusterProfile tmpAromaticProfile = tmpProfiler.getProfile(tmpLabels[0]);
        ClusterProfiler.ClusterProfile tmpAlcoholProfile = tmpProfiler.getProfile(tmpLabels[1]);
        //all aromatic molecules contain "c1ccccc1" and "C" (100 % vs. 50 % overall), all alcohols "O" and "CC"
        List<ClusterProfiler.FeatureProfile> tmpAromaticFragments = tmpAromaticProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, 2, true);
        Assertions.assertTrue(tmpAromaticFragments.stream().anyMatch(aFeature -> aFeature.featureName().equals("c1ccccc1")));
        Assertions.assertEquals(0.5, tmpAromaticFragments.getFirst().score(), 1.0e-9);
        List<ClusterProfiler.FeatureProfile> tmpAlcoholFragments = tmpAlcoholProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, 2, true);
        Assertions.assertTrue(tmpAlcoholFragments.stream().anyMatch(aFeature -> aFeature.featureName().equals("O")));
        //the benzene fragment is depleted in the alcohol cluster
        Assertions.assertTrue(tmpAlcoholProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, 2, false).stream()
                .anyMatch(aFeature -> aFeature.featureName().equals("c1ccccc1")));
        ClusterProfiler.FeatureProfile tmpAromaticAtoms = tmpAromaticProfile.getMostDistinctFeatures(
                ClusterProfiler.FeatureCategory.DESCRIPTOR, 5).stream()
                .filter(aFeature -> aFeature.featureName().equals("naAromAtom")).findFirst().orElseThrow();
        Assertions.assertTrue(tmpAromaticAtoms.score() > 0.0);
        StringWriter tmpWriter = new StringWriter();
        ClusteringDataCsvWriter.writeClusterProfiles(tmpWriter, List.of(tmpAromaticProfile, tmpAlcoholProfile), 3, ",");
        String[] tmpLines = tmpWriter.toString().split(System.lineSeparator());
        Assertions.assertEquals("Cluster,ClusterSize,Category,Feature,ClusterValue,OverallValue,Score", tmpLines[0]);
        Assertions.assertTrue(tmpLines.length > 2);
    }
}
