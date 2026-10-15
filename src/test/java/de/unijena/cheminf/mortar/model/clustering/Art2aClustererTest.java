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

import java.util.Arrays;

/**
 * Tests for {@link Art2aClusterer}.
 *
 * @author Manuel Schauer
 */
class Art2aClustererTest {
    /**
     * Creates a bit matrix with two groups of molecules that have disjoint sets of fragments.
     *
     * @return test matrix with 20 rows and 10 columns
     */
    private static FeatureMatrix createTwoGroupMatrix() {
        int tmpRowCount = 20;
        int tmpColumnCount = 10;
        float[][] tmpData = new float[tmpRowCount][tmpColumnCount];
        String[] tmpRowNames = new String[tmpRowCount];
        for (int i = 0; i < tmpRowCount; i++) {
            tmpRowNames[i] = "m" + i;
            //group 1 uses fragments 0-4, group 2 uses fragments 5-9; one fragment varies per molecule
            int tmpOffset = i < tmpRowCount / 2 ? 0 : 5;
            for (int j = 0; j < 4; j++) {
                tmpData[i][tmpOffset + j] = 1.0f;
            }
            tmpData[i][tmpOffset + 4] = i % 2;
        }
        String[] tmpNames = new String[tmpColumnCount];
        FeatureType[] tmpTypes = new FeatureType[tmpColumnCount];
        String[] tmpGroups = new String[tmpColumnCount];
        for (int j = 0; j < tmpColumnCount; j++) {
            tmpNames[j] = "f" + j;
            tmpTypes[j] = FeatureType.BIT;
            tmpGroups[j] = "F";
        }
        return new FeatureMatrix(tmpData, tmpRowNames, tmpNames, tmpTypes, tmpGroups);
    }
    //
    /**
     * Tests that both kernels separate two disjoint groups and that the labels are consistent with the groups.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testTwoGroupsAreSeparated() throws Exception {
        for (Art2aClusterer.Art2aAlgorithm tmpAlgorithm : Art2aClusterer.Art2aAlgorithm.values()) {
            float tmpConvergence = tmpAlgorithm == Art2aClusterer.Art2aAlgorithm.COSINE
                    ? Art2aClusterer.CONVERGENCE_THRESHOLD_COSINE_DEFAULT
                    : Art2aClusterer.CONVERGENCE_THRESHOLD_EUCLID_DEFAULT;
            ClusteringResult tmpResult = Art2aClusterer.cluster(Art2aClustererTest.createTwoGroupMatrix(), tmpAlgorithm,
                    0.5f, 20, Art2aClusterer.MAXIMUM_NUMBER_OF_EPOCHS_DEFAULT, tmpConvergence,
                    Art2aClusterer.LEARNING_PARAMETER_DEFAULT, Art2aClusterer.OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT,
                    Art2aClusterer.RANDOM_SEED_DEFAULT, false);
            int[] tmpLabels = tmpResult.clusterLabels();
            Assertions.assertEquals(2, tmpResult.numberOfClusters(), tmpAlgorithm.name());
            Assertions.assertFalse(tmpResult.isOverflow());
            Assertions.assertEquals(0, tmpResult.getNumberOfUnclusteredRows());
            //all molecules of a group share one label and the groups have different labels
            Assertions.assertEquals(1, Arrays.stream(tmpLabels, 0, 10).distinct().count());
            Assertions.assertEquals(1, Arrays.stream(tmpLabels, 10, 20).distinct().count());
            Assertions.assertNotEquals(tmpLabels[0], tmpLabels[10]);
            Assertions.assertEquals(10, tmpResult.getClusterSize(tmpLabels[0]));
            Assertions.assertEquals(2, tmpResult.representativeIndices().length);
        }
    }
    //
    /**
     * Tests that the clustering is reproducible for a fixed seed.
     *
     * @throws Exception if the clustering fails
     */
    @Test
    void testReproducibleWithFixedSeed() throws Exception {
        ClusteringResult tmpFirst = Art2aClusterer.cluster(Art2aClustererTest.createTwoGroupMatrix(),
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.9f, 20, 10, 0.99f, 0.01f, 1.0f, 5L, false);
        ClusteringResult tmpSecond = Art2aClusterer.cluster(Art2aClustererTest.createTwoGroupMatrix(),
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.9f, 20, 10, 0.99f, 0.01f, 1.0f, 5L, true);
        Assertions.assertArrayEquals(tmpFirst.clusterLabels(), tmpSecond.clusterLabels());
    }
}
