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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.CancellationException;

/**
 * Tests for {@link VigilanceScan} and {@link DiverseSubsetSelector}.
 *
 * @author Manuel Schauer
 */
class VigilanceScanAndSubsetTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    VigilanceScanAndSubsetTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the generation of the vigilance values including the end value and invalid ranges.
     */
    @Test
    void testGetVigilanceValues() {
        Assertions.assertArrayEquals(new float[] {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f},
                VigilanceScan.getVigilanceValues(0.1, 0.9, 0.1));
        Assertions.assertArrayEquals(new float[] {0.5f}, VigilanceScan.getVigilanceValues(0.5, 0.5, 0.1));
        Assertions.assertArrayEquals(new float[] {0.9f, 0.95f}, VigilanceScan.getVigilanceValues(0.9, 0.99, 0.05));
        Assertions.assertThrows(IllegalArgumentException.class, () -> VigilanceScan.getVigilanceValues(0.0, 0.9, 0.1));
        Assertions.assertThrows(IllegalArgumentException.class, () -> VigilanceScan.getVigilanceValues(0.5, 0.4, 0.1));
        Assertions.assertThrows(IllegalArgumentException.class, () -> VigilanceScan.getVigilanceValues(0.1, 0.9, 0.0));
        Assertions.assertThrows(IllegalArgumentException.class, () -> VigilanceScan.getVigilanceValues(0.1, 0.9, 0.00001));
    }
    //
    /**
     * Tests the scan on the aromatic/alcohol data: one entry per value, the results use the scanned vigilance and the
     * features of the base run, consistent counts, more clusters for a very high vigilance, and that the result of the
     * base vigilance is reproduced.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testScan() throws Exception {
        ClusteringRunResult tmpBase = SimilaritySearchTest.createRunResult();
        float[] tmpValues = {0.1f, tmpBase.clusteringResult().vigilance(), 0.99f};
        List<Integer> tmpProgress = new ArrayList<>();
        List<VigilanceScan.ScanEntry> tmpEntries = VigilanceScan.scan(tmpBase, tmpValues, () -> false,
                (aDone, aTotal) -> tmpProgress.add(aDone));
        Assertions.assertEquals(List.of(0, 1, 2, 3), tmpProgress);
        Assertions.assertEquals(3, tmpEntries.size());
        for (int v = 0; v < tmpValues.length; v++) {
            VigilanceScan.ScanEntry tmpEntry = tmpEntries.get(v);
            Assertions.assertEquals(tmpValues[v], tmpEntry.vigilance());
            Assertions.assertEquals(tmpValues[v], tmpEntry.runResult().parameters().vigilance());
            Assertions.assertSame(tmpBase.processedMatrix(), tmpEntry.runResult().processedMatrix());
            Assertions.assertTrue(tmpEntry.numberOfSingletons() <= tmpEntry.numberOfClusters());
            Assertions.assertTrue(tmpEntry.largestClusterSize() <= tmpBase.molecules().size());
            Assertions.assertTrue(tmpEntry.meanSimilarityToClusterVector() > 0.0 && tmpEntry.meanSimilarityToClusterVector() <= 1.0 + 1.0e-6);
        }
        Assertions.assertTrue(tmpEntries.get(2).numberOfClusters() >= tmpEntries.get(0).numberOfClusters());
        Assertions.assertArrayEquals(tmpBase.clusteringResult().clusterLabels(), tmpEntries.get(1).runResult().clusteringResult().clusterLabels());
        Assertions.assertThrows(CancellationException.class, () -> VigilanceScan.scan(tmpBase, tmpValues, () -> true, (aDone, aTotal) -> { }));
    }
    //
    /**
     * Tests the diverse subset: with one molecule per cluster every cluster contributes its centre (most typical
     * molecule) in both modes, the clusters are ordered by their index, and with many molecules per cluster all
     * clustered molecules are selected.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testDiverseSubset() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        SimilaritySearch tmpSearch = new SimilaritySearch(tmpRunResult);
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        for (DiverseSubsetSelector.SelectionMode tmpMode : DiverseSubsetSelector.SelectionMode.values()) {
            List<DiverseSubsetSelector.SelectedMolecule> tmpOnePerCluster = DiverseSubsetSelector.select(tmpRunResult, 1, false, tmpMode);
            Assertions.assertEquals(tmpRunResult.clusteringResult().numberOfClusters(), tmpOnePerCluster.size());
            for (int k = 0; k < tmpOnePerCluster.size(); k++) {
                DiverseSubsetSelector.SelectedMolecule tmpSelected = tmpOnePerCluster.get(k);
                //the clusters appear in the order of their index
                Assertions.assertEquals(k, tmpSelected.clusterIndex());
                Assertions.assertEquals(1, tmpSelected.rankInCluster());
                //no member of the cluster is more typical than the selected centre
                for (int i = 0; i < tmpLabels.length; i++) {
                    if (tmpLabels[i] == k) {
                        Assertions.assertTrue(tmpSearch.getSimilarityToOwnClusterVector(i) <= tmpSelected.similarityToClusterVector() + 1.0e-9);
                    }
                }
            }
            Assertions.assertEquals(tmpRunResult.molecules().size(), DiverseSubsetSelector.select(tmpRunResult, 1000, true, tmpMode).size());
        }
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> DiverseSubsetSelector.select(tmpRunResult, 0, false, DiverseSubsetSelector.SelectionMode.DIVERSE));
    }
    //
    /**
     * Tests the selection from one cluster: typical takes the members in the order of typicality, diverse takes the centre
     * and then the member least similar to it (MaxMin); both select different molecules only.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testSelectFromCluster() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        SimilaritySearch tmpSearch = new SimilaritySearch(tmpRunResult);
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        //the cluster of the first aromatic molecule
        int tmpCluster = tmpLabels[0];
        List<Integer> tmpMembers = new ArrayList<>();
        for (int i = 0; i < tmpLabels.length; i++) {
            if (tmpLabels[i] == tmpCluster) {
                tmpMembers.add(i);
            }
        }
        Assertions.assertTrue(tmpMembers.size() >= 3, "test data must give a cluster with at least 3 members");
        List<DiverseSubsetSelector.SelectedMolecule> tmpTypical = DiverseSubsetSelector.selectFromCluster(tmpRunResult, tmpCluster, 3,
                DiverseSubsetSelector.SelectionMode.TYPICAL);
        List<DiverseSubsetSelector.SelectedMolecule> tmpDiverse = DiverseSubsetSelector.selectFromCluster(tmpRunResult, tmpCluster, 3,
                DiverseSubsetSelector.SelectionMode.DIVERSE);
        Assertions.assertEquals(3, tmpTypical.size());
        Assertions.assertEquals(3, tmpDiverse.size());
        //typical: decreasing similarity to the cluster vector
        Assertions.assertTrue(tmpTypical.get(0).similarityToClusterVector() >= tmpTypical.get(1).similarityToClusterVector());
        Assertions.assertTrue(tmpTypical.get(1).similarityToClusterVector() >= tmpTypical.get(2).similarityToClusterVector());
        //both start with the centre
        Assertions.assertEquals(tmpTypical.getFirst().rowIndex(), tmpDiverse.getFirst().rowIndex());
        //diverse: the second molecule is the member least similar to the centre
        int tmpCentre = tmpDiverse.getFirst().rowIndex();
        double tmpLowestSimilarity = Double.POSITIVE_INFINITY;
        for (int tmpMember : tmpMembers) {
            if (tmpMember != tmpCentre) {
                tmpLowestSimilarity = Math.min(tmpLowestSimilarity, tmpSearch.getSimilarityBetweenRows(tmpMember, tmpCentre));
            }
        }
        Assertions.assertEquals(tmpLowestSimilarity, tmpSearch.getSimilarityBetweenRows(tmpDiverse.get(1).rowIndex(), tmpCentre), 1.0e-9);
        Assertions.assertEquals(3, tmpDiverse.stream().map(DiverseSubsetSelector.SelectedMolecule::rowIndex).distinct().count());
        Assertions.assertEquals(List.of(1, 2, 3), tmpDiverse.stream().map(DiverseSubsetSelector.SelectedMolecule::rankInCluster).toList());
        Assertions.assertEquals(tmpMembers.size(), DiverseSubsetSelector.selectFromCluster(tmpRunResult, tmpCluster, 1000,
                DiverseSubsetSelector.SelectionMode.DIVERSE).size());
        Assertions.assertThrows(IllegalArgumentException.class, () -> DiverseSubsetSelector.selectFromCluster(tmpRunResult,
                tmpRunResult.clusteringResult().numberOfClusters(), 1, DiverseSubsetSelector.SelectionMode.TYPICAL));
    }
}
