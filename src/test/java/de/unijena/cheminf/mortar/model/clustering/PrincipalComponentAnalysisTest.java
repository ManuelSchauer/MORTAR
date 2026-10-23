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

import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.concurrent.CancellationException;

/**
 * Tests for {@link PrincipalComponentAnalysis}.
 *
 * @author Manuel Schauer
 */
class PrincipalComponentAnalysisTest {
    /**
     * Constructor that sets the default locale to British English, which is needed for the message resources of the
     * fragmenter enums.
     */
    PrincipalComponentAnalysisTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the PCA of data with known variances along the axes (x: 9, y: 4, z: 1 in expectation, rotated randomly is
     * not needed): the eigenvalues equal the variances of the scores, are descending, the scores are uncorrelated, the
     * explained variance ratios sum up to 1 for all components, and a data row projected as additional vector gives its
     * score.
     */
    @Test
    void testPropertiesOfComponents() {
        Random tmpRandom = new Random(42L);
        int tmpRowCount = 500;
        double[][] tmpData = new double[tmpRowCount][3];
        for (int i = 0; i < tmpRowCount; i++) {
            double tmpA = 3.0 * tmpRandom.nextGaussian();
            double tmpB = 2.0 * tmpRandom.nextGaussian();
            double tmpC = tmpRandom.nextGaussian();
            //mixed into the three features, so that the components are not the axes
            tmpData[i][0] = tmpA + 0.5 * tmpB + 5.0;
            tmpData[i][1] = -0.5 * tmpA + tmpB + tmpC;
            tmpData[i][2] = 0.3 * tmpA - tmpC - 2.0;
        }
        PrincipalComponentAnalysis.PcaResult tmpResult = PrincipalComponentAnalysis.calculate(tmpData, new double[][] {tmpData[7]},
                3, () -> false);
        Assertions.assertEquals(3, tmpResult.getNumberOfComponents());
        double[] tmpEigenvalues = tmpResult.eigenvalues();
        Assertions.assertTrue(tmpEigenvalues[0] >= tmpEigenvalues[1] && tmpEigenvalues[1] >= tmpEigenvalues[2]);
        double tmpRatioSum = 0.0;
        for (double tmpRatio : tmpResult.explainedVarianceRatios()) {
            tmpRatioSum += tmpRatio;
        }
        Assertions.assertEquals(1.0, tmpRatioSum, 1.0e-8);
        double[][] tmpScores = tmpResult.scores();
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                double tmpCovariance = 0.0;
                for (double[] tmpScore : tmpScores) {
                    tmpCovariance += tmpScore[a] * tmpScore[b];
                }
                tmpCovariance /= tmpRowCount - 1;
                Assertions.assertEquals(a == b ? tmpEigenvalues[a] : 0.0, tmpCovariance, 1.0e-6 * tmpEigenvalues[0]);
            }
        }
        for (int c = 0; c < 3; c++) {
            Assertions.assertEquals(tmpScores[7][c], tmpResult.clusterVectorScores()[0][c], 1.0e-9);
        }
        //reproducible
        PrincipalComponentAnalysis.PcaResult tmpSecond = PrincipalComponentAnalysis.calculate(tmpData, new double[0][], 3, () -> false);
        Assertions.assertArrayEquals(tmpScores[0], tmpSecond.scores()[0], 1.0e-9);
    }
    //
    /**
     * Tests that collinear data has one component with the whole variance, that the number of components is limited by
     * the data and that invalid input and cancellation are handled.
     */
    @Test
    void testCollinearDataAndLimits() {
        double[][] tmpData = new double[10][];
        for (int i = 0; i < tmpData.length; i++) {
            tmpData[i] = new double[] {i, 2.0 * i, -i};
        }
        PrincipalComponentAnalysis.PcaResult tmpResult = PrincipalComponentAnalysis.calculate(tmpData, new double[0][], 2, () -> false);
        Assertions.assertEquals(1.0, tmpResult.explainedVarianceRatios()[0], 1.0e-9);
        Assertions.assertEquals(0.0, tmpResult.explainedVarianceRatios()[1], 1.0e-9);
        //two rows allow only one component
        Assertions.assertEquals(1, PrincipalComponentAnalysis.calculate(new double[][] {{1.0, 2.0}, {3.0, 1.0}}, new double[0][], 3,
                () -> false).getNumberOfComponents());
        Assertions.assertThrows(IllegalArgumentException.class, () -> PrincipalComponentAnalysis.calculate(new double[][] {{1.0}},
                new double[0][], 2, () -> false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> PrincipalComponentAnalysis.calculate(tmpData, new double[0][], 0,
                () -> false));
        Assertions.assertThrows(CancellationException.class, () -> PrincipalComponentAnalysis.calculate(tmpData, new double[0][], 2,
                () -> true));
    }
    //
    /**
     * Tests the PCA of a clustering result: one score per molecule and one per cluster vector.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testPcaOfClusteringResult() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        PrincipalComponentAnalysis.PcaResult tmpResult = PrincipalComponentAnalysis.calculate(tmpRunResult,
                PrincipalComponentAnalysis.NUMBER_OF_COMPONENTS_DEFAULT, () -> false);
        Assertions.assertEquals(tmpRunResult.molecules().size(), tmpResult.scores().length);
        Assertions.assertEquals(tmpRunResult.clusteringResult().numberOfClusters(), tmpResult.clusterVectorScores().length);
        Assertions.assertTrue(tmpResult.explainedVarianceRatios()[0] > 0.0);
    }
    //
    /**
     * Tests the loadings: the component vectors are orthonormal, the scores are the centred data projected onto them,
     * and the feature with by far the largest variance has the largest weight in the first component.
     */
    @Test
    void testLoadings() {
        Random tmpRandom = new Random(7L);
        int tmpRowCount = 300;
        double[][] tmpData = new double[tmpRowCount][4];
        for (int i = 0; i < tmpRowCount; i++) {
            tmpData[i][0] = 0.1 * tmpRandom.nextGaussian();
            tmpData[i][1] = 5.0 * tmpRandom.nextGaussian();
            tmpData[i][2] = 0.5 * tmpRandom.nextGaussian();
            tmpData[i][3] = 0.2 * tmpRandom.nextGaussian();
        }
        PrincipalComponentAnalysis.PcaResult tmpResult = PrincipalComponentAnalysis.calculate(tmpData, new double[0][], 3, () -> false);
        double[][] tmpLoadings = tmpResult.loadings();
        Assertions.assertEquals(3, tmpLoadings.length);
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                double tmpDot = 0.0;
                for (int j = 0; j < 4; j++) {
                    tmpDot += tmpLoadings[a][j] * tmpLoadings[b][j];
                }
                Assertions.assertEquals(a == b ? 1.0 : 0.0, tmpDot, 1.0e-9);
            }
        }
        //scores = centred data x loadings
        double[] tmpMeans = new double[4];
        for (double[] tmpRow : tmpData) {
            for (int j = 0; j < 4; j++) {
                tmpMeans[j] += tmpRow[j] / tmpRowCount;
            }
        }
        for (int c = 0; c < 3; c++) {
            double tmpScore = 0.0;
            for (int j = 0; j < 4; j++) {
                tmpScore += (tmpData[5][j] - tmpMeans[j]) * tmpLoadings[c][j];
            }
            Assertions.assertEquals(tmpScore, tmpResult.scores()[5][c], 1.0e-9);
        }
        List<PrincipalComponentAnalysis.FeatureLoading> tmpTop = tmpResult.getTopLoadings(0, 2);
        Assertions.assertEquals(2, tmpTop.size());
        Assertions.assertEquals(1, tmpTop.getFirst().featureIndex());
        Assertions.assertTrue(Math.abs(tmpTop.getFirst().loading()) > 0.99);
        Assertions.assertTrue(Math.abs(tmpTop.getFirst().loading()) >= Math.abs(tmpTop.get(1).loading()));
    }
}
