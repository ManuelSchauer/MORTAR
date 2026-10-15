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

/**
 * Tests for {@link TrainTestSplitter}.
 *
 * @author Manuel Schauer
 */
class TrainTestSplitterTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    TrainTestSplitterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates a clustering result with the given cluster labels (9 molecules, random features).
     *
     * @param aLabels cluster label of every molecule, -1 for not clusterable
     * @param aNumberOfClusters number of clusters
     * @return clustering run result
     */
    private static ClusteringRunResult createRunResult(int[] aLabels, int aNumberOfClusters) {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        String[] tmpRowNames = new String[aLabels.length];
        float[][] tmpData = new float[aLabels.length][2];
        for (int i = 0; i < aLabels.length; i++) {
            tmpMolecules.add(new MoleculeDataModel("C", "m" + i, new HashMap<>()));
            tmpRowNames[i] = "m" + i;
            tmpData[i][0] = i;
            tmpData[i][1] = aLabels.length - i;
        }
        FeatureMatrix tmpMatrix = new FeatureMatrix(tmpData, tmpRowNames, new String[] {"a", "b"},
                new FeatureType[] {FeatureType.CONTINUOUS, FeatureType.CONTINUOUS}, new String[] {"A", "B"});
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", false, 1, 0, new Descriptor[0],
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        int[] tmpRepresentatives = new int[aNumberOfClusters];
        for (int k = 0; k < aNumberOfClusters; k++) {
            for (int i = 0; i < aLabels.length; i++) {
                if (aLabels[i] == k) {
                    tmpRepresentatives[k] = i;
                    break;
                }
            }
        }
        return new ClusteringRunResult(tmpParameters, tmpMolecules, tmpMatrix, tmpMatrix, tmpMatrix, List.of(),
                new ClusteringResult(aLabels, aNumberOfClusters, tmpRepresentatives, false, true, 2, 0.5f, new float[aNumberOfClusters][2]));
    }
    //
    /**
     * Tests the parsing of the training/test ratio.
     */
    @Test
    void testParseTrainingFraction() {
        Assertions.assertEquals(0.4, TrainTestSplitter.parseTrainingFraction("40:60"), 1.0e-9);
        Assertions.assertEquals(0.5, TrainTestSplitter.parseTrainingFraction(" 1 : 1 "), 1.0e-9);
        Assertions.assertEquals(0.8, TrainTestSplitter.parseTrainingFraction("80"), 1.0e-9);
        Assertions.assertEquals(0.8, TrainTestSplitter.parseTrainingFraction("80%"), 1.0e-9);
        Assertions.assertEquals(0.7, TrainTestSplitter.parseTrainingFraction("0.7"), 1.0e-9);
        Assertions.assertEquals(0.7, TrainTestSplitter.parseTrainingFraction("0,7"), 1.0e-9);
        for (String tmpInvalid : new String[] {"", "abc", "0:0", "100", "0", "100:0", "1:2:3", "-10:20"}) {
            Assertions.assertThrows(IllegalArgumentException.class, () -> TrainTestSplitter.parseTrainingFraction(tmpInvalid), tmpInvalid);
        }
    }
    //
    /**
     * Tests the rules: a single-molecule cluster and the not clusterable molecule go into the training set, every
     * larger cluster gets round(n * fraction) training molecules but at least one in each set.
     */
    @Test
    void testSplitRules() {
        //cluster 0: row 0; cluster 1: rows 1-5; cluster 2: rows 6-7; row 8 not clusterable
        int[] tmpLabels = {0, 1, 1, 1, 1, 1, 2, 2, -1};
        ClusteringRunResult tmpRunResult = TrainTestSplitterTest.createRunResult(tmpLabels, 3);
        TrainTestSplitter.TrainTestSplit tmpSplit = TrainTestSplitter.split(tmpRunResult, 0.4, TrainTestSplitter.SelectionMode.RANDOM, 1L);
        Assertions.assertTrue(tmpSplit.isTraining(0));
        Assertions.assertTrue(tmpSplit.isTraining(8));
        Assertions.assertEquals(2, TrainTestSplitterTest.countTraining(tmpSplit, tmpLabels, 1));
        Assertions.assertEquals(1, TrainTestSplitterTest.countTraining(tmpSplit, tmpLabels, 2));
        Assertions.assertEquals(5, tmpSplit.trainingIndices().length);
        Assertions.assertEquals(4, tmpSplit.testIndices().length);
        Assertions.assertEquals(5.0 / 9.0, tmpSplit.getActualTrainingFraction(), 1.0e-9);
        //a high fraction still leaves one test molecule per cluster, a low one still one training molecule
        TrainTestSplitter.TrainTestSplit tmpHighSplit = TrainTestSplitter.split(tmpRunResult, 0.95, TrainTestSplitter.SelectionMode.RANDOM, 1L);
        Assertions.assertEquals(4, TrainTestSplitterTest.countTraining(tmpHighSplit, tmpLabels, 1));
        Assertions.assertEquals(1, TrainTestSplitterTest.countTraining(tmpHighSplit, tmpLabels, 2));
        TrainTestSplitter.TrainTestSplit tmpLowSplit = TrainTestSplitter.split(tmpRunResult, 0.05, TrainTestSplitter.SelectionMode.RANDOM, 1L);
        Assertions.assertEquals(1, TrainTestSplitterTest.countTraining(tmpLowSplit, tmpLabels, 1));
        //every molecule is in exactly one set and the random selection is reproducible
        Assertions.assertEquals(tmpLabels.length, tmpLowSplit.trainingIndices().length + tmpLowSplit.testIndices().length);
        Assertions.assertArrayEquals(tmpSplit.trainingIndices(),
                TrainTestSplitter.split(tmpRunResult, 0.4, TrainTestSplitter.SelectionMode.RANDOM, 1L).trainingIndices());
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> TrainTestSplitter.split(tmpRunResult, 1.0, TrainTestSplitter.SelectionMode.RANDOM, 1L));
    }
    //
    /**
     * Tests that with the selection closest to the cluster centre, every training molecule of a cluster is at least as
     * similar to the cluster vector as every test molecule of that cluster.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testClosestToClusterCenterSelection() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        TrainTestSplitter.TrainTestSplit tmpSplit = TrainTestSplitter.split(tmpRunResult, 0.5,
                TrainTestSplitter.SelectionMode.CLOSEST_TO_CLUSTER_CENTER, 1L);
        SimilaritySearch tmpSearch = new SimilaritySearch(tmpRunResult);
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        for (int tmpTraining : tmpSplit.trainingIndices()) {
            for (int tmpTest : tmpSplit.testIndices()) {
                if (tmpLabels[tmpTraining] == tmpLabels[tmpTest]) {
                    Assertions.assertTrue(tmpSearch.getSimilarityToOwnClusterVector(tmpTraining)
                            >= tmpSearch.getSimilarityToOwnClusterVector(tmpTest));
                }
            }
        }
        Assertions.assertTrue(tmpSplit.testIndices().length > 0);
    }
    //
    /**
     * Returns the number of training molecules of the given cluster.
     *
     * @param aSplit split
     * @param aLabels cluster labels
     * @param aCluster cluster index
     * @return number of training molecules of the cluster
     */
    private static int countTraining(TrainTestSplitter.TrainTestSplit aSplit, int[] aLabels, int aCluster) {
        int tmpCount = 0;
        for (int tmpIndex : aSplit.trainingIndices()) {
            if (aLabels[tmpIndex] == aCluster) {
                tmpCount++;
            }
        }
        return tmpCount;
    }
}
