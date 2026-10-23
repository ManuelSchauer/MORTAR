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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests for {@link ClusterDescriptorStatistics} (medians and range filter) and for the restriction of the training/test
 * partitioning and the diverse subset to marked clusters.
 *
 * @author Manuel Schauer
 */
class ClusterDescriptorStatisticsTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the message resources of the
     * fragmenter enums.
     */
    ClusterDescriptorStatisticsTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates a clustering result with six molecules: cluster 0 = rows 0-2, cluster 1 = rows 3-4, row 5 not clusterable.
     * Raw columns: a fragment bit, "LogP" (1, 3, NaN | 5, 7 | 100) and "TPSA" (10, 20, 30 | 40, 50 | 60); "TPSA" is not
     * part of the processed matrix, i.e. removed by the filter.
     *
     * @return clustering result
     */
    private static ClusteringRunResult createRunResult() {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            tmpMolecules.add(new MoleculeDataModel("C", "m" + i, new HashMap<>()));
        }
        String[] tmpRowNames = {"m0", "m1", "m2", "m3", "m4", "m5"};
        FeatureMatrix tmpRaw = new FeatureMatrix(new float[][] {
                {1.0f, 1.0f, 10.0f}, {1.0f, 3.0f, 20.0f}, {0.0f, Float.NaN, 30.0f},
                {0.0f, 5.0f, 40.0f}, {1.0f, 7.0f, 50.0f}, {0.0f, 100.0f, 60.0f}},
                tmpRowNames, new String[] {"C", "LogP", "TPSA"},
                new FeatureType[] {FeatureType.BIT, FeatureType.CONTINUOUS, FeatureType.CONTINUOUS}, new String[] {"Frag", "X_LOG_P", "TPSA"});
        FeatureMatrix tmpProcessed = tmpRaw.selectColumns(new int[] {0, 1});
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", false, 1, 0, new Descriptor[0],
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        return new ClusteringRunResult(tmpParameters, tmpMolecules, tmpRaw, tmpRaw, tmpProcessed, List.of(),
                new ClusteringResult(new int[] {0, 0, 0, 1, 1, -1}, 2, new int[] {0, 3}, false, true, 2, 0.5f, new float[2][2]));
    }
    //
    /**
     * Tests the medians per cluster and overall from the raw values: NaN values are ignored, an even number of values
     * gives the mean of the two middle values, the not clusterable molecules have no median but count overall.
     */
    @Test
    void testMedians() {
        ClusterDescriptorStatistics tmpStatistics = new ClusterDescriptorStatistics(ClusterDescriptorStatisticsTest.createRunResult());
        Assertions.assertEquals(List.of("LogP", "TPSA"), tmpStatistics.getDescriptorColumnNames());
        Assertions.assertEquals(2.0, tmpStatistics.getMedian(0, "LogP"));
        Assertions.assertEquals(6.0, tmpStatistics.getMedian(1, "LogP"));
        Assertions.assertEquals(5.0, tmpStatistics.getOverallMedian("LogP"));
        Assertions.assertEquals(20.0, tmpStatistics.getMedian(0, "TPSA"));
        Assertions.assertEquals(45.0, tmpStatistics.getMedian(1, "TPSA"));
        Assertions.assertEquals(35.0, tmpStatistics.getOverallMedian("TPSA"));
        Assertions.assertTrue(Double.isNaN(tmpStatistics.getMedian(-1, "LogP")));
        Assertions.assertTrue(Double.isNaN(tmpStatistics.getValue(2, "LogP")));
        Assertions.assertEquals(100.0, tmpStatistics.getValue(5, "LogP"));
        Assertions.assertTrue(tmpStatistics.isUsedForClustering("LogP"));
        Assertions.assertFalse(tmpStatistics.isUsedForClustering("TPSA"));
        Assertions.assertFalse(tmpStatistics.hasDescriptorColumn("C"));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpStatistics.getMedian(0, "C"));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpStatistics.getMedian(2, "LogP"));
        Assertions.assertTrue(Double.isNaN(ClusterDescriptorStatistics.getMedianOfSorted(new double[0], 0)));
        Assertions.assertEquals(2.0, ClusterDescriptorStatistics.getMedianOfSorted(new double[] {1.0, 2.0, 9.0}, 3));
    }
    //
    /**
     * Tests the range filter: bounds are included, a missing bound is open, all conditions must hold and without
     * conditions no cluster matches.
     */
    @Test
    void testRangeFilter() {
        ClusterDescriptorStatistics tmpStatistics = new ClusterDescriptorStatistics(ClusterDescriptorStatisticsTest.createRunResult());
        Assertions.assertEquals(Set.of(0), tmpStatistics.getMatchingClusters(List.of(
                new ClusterDescriptorStatistics.RangeCondition("LogP", 0.0, 2.0))));
        Assertions.assertEquals(Set.of(0, 1), tmpStatistics.getMatchingClusters(List.of(
                new ClusterDescriptorStatistics.RangeCondition("LogP", Double.NaN, Double.NaN))));
        Assertions.assertEquals(Set.of(1), tmpStatistics.getMatchingClusters(List.of(
                new ClusterDescriptorStatistics.RangeCondition("LogP", 4.0, Double.NaN),
                new ClusterDescriptorStatistics.RangeCondition("TPSA", Double.NaN, 45.0))));
        Assertions.assertEquals(Set.of(), tmpStatistics.getMatchingClusters(List.of(
                new ClusterDescriptorStatistics.RangeCondition("LogP", Double.NaN, 4.0),
                new ClusterDescriptorStatistics.RangeCondition("TPSA", 30.0, Double.NaN))));
        Assertions.assertTrue(tmpStatistics.getMatchingClusters(List.of()).isEmpty());
        Assertions.assertFalse(new ClusterDescriptorStatistics.RangeCondition("LogP", Double.NaN, Double.NaN).isSatisfiedBy(Double.NaN));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpStatistics.getMatchingClusters(List.of(
                new ClusterDescriptorStatistics.RangeCondition("Unknown", 0.0, 1.0))));
    }
    //
    /**
     * Tests the row indices of marked clusters and the training/test partitioning restricted to marked clusters: the
     * molecules of other clusters and the not clusterable molecules are in neither set.
     */
    @Test
    void testRestrictionToMarkedClusters() {
        ClusteringRunResult tmpRunResult = ClusterDescriptorStatisticsTest.createRunResult();
        Assertions.assertArrayEquals(new int[] {3, 4}, tmpRunResult.getRowIndicesOfClusters(Set.of(1)));
        Assertions.assertArrayEquals(new int[] {0, 1, 2, 3, 4, 5}, tmpRunResult.getRowIndicesOfClusters(Set.of()));
        Assertions.assertArrayEquals(new int[] {0, 1, 2, 3, 4, 5}, tmpRunResult.getRowIndicesOfClusters(null));
        TrainTestSplitter.TrainTestSplit tmpSplit = TrainTestSplitter.split(tmpRunResult, 0.5,
                TrainTestSplitter.SelectionMode.RANDOM, 1L, Set.of(0));
        int[] tmpAll = java.util.stream.IntStream.concat(Arrays.stream(tmpSplit.trainingIndices()), Arrays.stream(tmpSplit.testIndices()))
                .sorted().toArray();
        Assertions.assertArrayEquals(new int[] {0, 1, 2}, tmpAll);
        Assertions.assertTrue(tmpSplit.trainingIndices().length >= 1 && tmpSplit.testIndices().length >= 1);
        //without marks all molecules are partitioned, the not clusterable one into the training set
        TrainTestSplitter.TrainTestSplit tmpFullSplit = TrainTestSplitter.split(tmpRunResult, 0.5,
                TrainTestSplitter.SelectionMode.RANDOM, 1L, Set.of());
        Assertions.assertEquals(6, tmpFullSplit.trainingIndices().length + tmpFullSplit.testIndices().length);
        Assertions.assertTrue(tmpFullSplit.isTraining(5));
    }
    //
    /**
     * Tests the diverse subset restricted to a marked cluster: only molecules of this cluster are selected and the not
     * clusterable molecules are not added.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testDiverseSubsetOfMarkedClusters() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        int tmpCluster = tmpRunResult.clusteringResult().clusterLabels()[0];
        List<DiverseSubsetSelector.SelectedMolecule> tmpSelection = DiverseSubsetSelector.select(tmpRunResult, 2, true,
                DiverseSubsetSelector.SelectionMode.DIVERSE, Set.of(tmpCluster));
        Assertions.assertFalse(tmpSelection.isEmpty());
        Assertions.assertTrue(tmpSelection.stream().allMatch(aSelected -> aSelected.clusterIndex() == tmpCluster));
        Assertions.assertTrue(tmpSelection.size() <= 2);
    }
}
