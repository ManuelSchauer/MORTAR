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
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests for {@link ClusterQuality} and {@link ClusterStability}.
 *
 * @author Manuel Schauer
 */
class ClusterQualityTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the message resources of the
     * fragmenter enums.
     */
    ClusterQualityTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates a clustering result with two well separated groups of three molecules and one group that lies between
     * them, with hand-made processed features (two features) and labels; the cosine kernel compares the kernel vectors.
     *
     * @return clustering result
     */
    private static ClusteringRunResult createRunResult() {
        float[][] tmpData = {{1.0f, 0.0f}, {0.95f, 0.05f}, {0.9f, 0.1f}, {0.0f, 1.0f}, {0.05f, 0.95f}, {0.1f, 0.9f}};
        String[] tmpNames = {"a1", "a2", "a3", "b1", "b2", "b3"};
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        for (String tmpName : tmpNames) {
            tmpMolecules.add(new MoleculeDataModel("C", tmpName, new HashMap<>()));
        }
        FeatureMatrix tmpMatrix = new FeatureMatrix(tmpData, tmpNames, new String[] {"x", "y"},
                new FeatureType[] {FeatureType.CONTINUOUS, FeatureType.CONTINUOUS}, new String[] {"X", "Y"});
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", false, 1, 0, new Descriptor[0],
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        return new ClusteringRunResult(tmpParameters, tmpMolecules, tmpMatrix, tmpMatrix, tmpMatrix, List.of(),
                new ClusteringResult(new int[] {0, 0, 0, 1, 1, 1}, 2, new int[] {1, 4}, false, true, 2, 0.5f,
                        new float[][] {{0.95f, 0.05f}, {0.05f, 0.95f}}));
    }
    //
    /**
     * Tests the silhouette formula and the number of neighbours depending on the cluster size.
     */
    @Test
    void testSilhouetteFormula() {
        //number of neighbours: at most 10, at most the other molecules of the cluster
        Assertions.assertEquals(0, ClusterQuality.getNumberOfNeighbours(1));
        Assertions.assertEquals(2, ClusterQuality.getNumberOfNeighbours(3));
        Assertions.assertEquals(10, ClusterQuality.getNumberOfNeighbours(11));
        Assertions.assertEquals(10, ClusterQuality.getNumberOfNeighbours(500));
        Assertions.assertEquals(0.75, ClusterQuality.getSilhouette(0.1, 0.4), 1.0e-12);
        Assertions.assertEquals(-0.75, ClusterQuality.getSilhouette(0.4, 0.1), 1.0e-12);
        Assertions.assertEquals(0.0, ClusterQuality.getSilhouette(0.2, 0.2), 1.0e-12);
        Assertions.assertEquals(0.0, ClusterQuality.getSilhouette(0.0, 0.0), 1.0e-12);
    }
    //
    /**
     * Tests two well separated clusters: high silhouette, all nearest neighbours in the own cluster, symmetric mean
     * distances that are larger between than within the clusters, and a silhouette computed by hand for one molecule.
     */
    @Test
    void testSeparatedClusters() {
        ClusteringRunResult tmpRunResult = ClusterQualityTest.createRunResult();
        ClusterQuality tmpQuality = new ClusterQuality(tmpRunResult);
        Assertions.assertFalse(tmpQuality.isSampled());
        Assertions.assertEquals(6, tmpQuality.getNumberOfEvaluatedMolecules());
        Assertions.assertTrue(tmpQuality.getOverallSilhouette() > 0.8);
        //clusters of 3 molecules use k = 2 neighbours, the other two molecules of the cluster: 100 % instead of 2 of 10
        Assertions.assertEquals(1.0, tmpQuality.getOverallNeighbourPurity(), 1.0e-12);
        for (ClusterQuality.ClusterQualityEntry tmpEntry : tmpQuality.getEntries()) {
            Assertions.assertEquals(3, tmpEntry.size());
            Assertions.assertEquals(2, tmpEntry.numberOfNeighbours());
            Assertions.assertEquals(1.0, tmpEntry.neighbourPurity(), 1.0e-12);
            Assertions.assertEquals(0, tmpEntry.numberOfNegativeSilhouettes());
            Assertions.assertEquals(1 - tmpEntry.clusterIndex(), tmpEntry.nearestClusterIndex());
            Assertions.assertTrue(tmpEntry.meanDistanceWithinCluster() < tmpEntry.meanDistanceToNearestCluster());
        }
        Assertions.assertEquals(tmpQuality.getMeanDistance(0, 1), tmpQuality.getMeanDistance(1, 0), 1.0e-12);
        //silhouette of molecule a1 by hand: cosine kernel vectors are the normalised rows
        SimilaritySearch tmpSearch = new SimilaritySearch(tmpRunResult);
        double tmpOwn = 0.0;
        double tmpOther = 0.0;
        for (int j = 1; j < 3; j++) {
            tmpOwn += 1.0 - tmpSearch.getSimilarityBetweenRows(0, j);
        }
        for (int j = 3; j < 6; j++) {
            tmpOther += 1.0 - tmpSearch.getSimilarityBetweenRows(0, j);
        }
        Assertions.assertEquals(ClusterQuality.getSilhouette(tmpOwn / 2.0, tmpOther / 3.0), tmpQuality.getMoleculeSilhouette(0), 1.0e-9);
        List<ClusterQuality.ClusterDistance> tmpDistances = tmpQuality.getDistancesToOtherClusters(0);
        Assertions.assertEquals(1, tmpDistances.size());
        Assertions.assertEquals(1, tmpDistances.getFirst().clusterIndex());
        Assertions.assertEquals(tmpQuality.getMeanDistance(0, 1), tmpDistances.getFirst().meanDistance(), 1.0e-12);
        Assertions.assertEquals(tmpQuality.getCentreDistance(0, 1), tmpDistances.getFirst().centreDistance(), 1.0e-12);
    }
    //
    /**
     * Tests the quality of a clustering result of the pipeline: one entry per cluster with consistent sizes and values
     * in their ranges.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testQualityOfPipelineResult() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        ClusterQuality tmpQuality = new ClusterQuality(tmpRunResult);
        Assertions.assertEquals(tmpRunResult.clusteringResult().numberOfClusters(), tmpQuality.getEntries().size());
        for (ClusterQuality.ClusterQualityEntry tmpEntry : tmpQuality.getEntries()) {
            Assertions.assertEquals(tmpRunResult.getRowIndicesOfClusters(Set.of(tmpEntry.clusterIndex())).length, tmpEntry.size());
            if (!Double.isNaN(tmpEntry.silhouette())) {
                Assertions.assertTrue(tmpEntry.silhouette() >= -1.0 && tmpEntry.silhouette() <= 1.0);
            }
            Assertions.assertTrue(tmpEntry.neighbourPurity() >= 0.0 && tmpEntry.neighbourPurity() <= 1.0);
        }
    }
    //
    /**
     * Tests the Jaccard indices and the adjusted Rand index: identical and relabelled assignments give 1, a split
     * cluster gives the expected Jaccard index, and the stability check of a pipeline result gives values in range.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testStability() throws Exception {
        int[] tmpLabels = {0, 0, 0, 1, 1, 1};
        int[] tmpRelabelled = {1, 1, 1, 0, 0, 0};
        Assertions.assertEquals(1.0, ClusterStability.getAdjustedRandIndex(tmpLabels, tmpRelabelled), 1.0e-12);
        Assertions.assertArrayEquals(new double[] {1.0, 1.0}, ClusterStability.getBestJaccardIndices(tmpLabels, 2, tmpRelabelled, 2), 1.0e-12);
        //cluster 0 split into {0,1} and {2}: best Jaccard 2/3
        int[] tmpSplit = {0, 0, 2, 1, 1, 1};
        Assertions.assertEquals(2.0 / 3.0, ClusterStability.getBestJaccardIndices(tmpLabels, 2, tmpSplit, 3)[0], 1.0e-12);
        Assertions.assertTrue(ClusterStability.getAdjustedRandIndex(tmpLabels, tmpSplit) < 1.0);
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        ClusterStability.StabilityResult tmpResult = ClusterStability.check(tmpRunResult, 3, () -> false, (aDone, aTotal) -> { });
        Assertions.assertEquals(3, tmpResult.numberOfRuns());
        Assertions.assertEquals(tmpRunResult.clusteringResult().numberOfClusters(), tmpResult.clusterJaccardIndices().length);
        for (double tmpJaccard : tmpResult.clusterJaccardIndices()) {
            Assertions.assertTrue(tmpJaccard >= 0.0 && tmpJaccard <= 1.0);
        }
        Assertions.assertTrue(tmpResult.meanAdjustedRandIndex() <= 1.0 + 1.0e-12);
        Assertions.assertThrows(java.util.concurrent.CancellationException.class,
                () -> ClusterStability.check(tmpRunResult, 3, () -> true, (aDone, aTotal) -> { }));
    }
}
