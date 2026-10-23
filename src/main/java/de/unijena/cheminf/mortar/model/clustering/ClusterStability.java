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

import java.util.Objects;
import java.util.concurrent.CancellationException;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;

/**
 * Stability of the clusters of a clustering result: ART-2A depends on the order in which the molecules are presented,
 * which is shuffled with the random seed. The processed features of the result are therefore clustered again several
 * times with other seeds (all other parameters are kept) and every run is compared with the result:
 * <ul>
 *     <li>Per cluster: the Jaccard index (size of the intersection divided by the size of the union) of its molecules
 *     with the best matching cluster of the run, averaged over the runs. Near 1 the cluster is found again in every run, a low value means it falls apart
 *     differently from run to run and is rather an artefact of the order.</li>
 *     <li>Overall: the adjusted Rand index (ARI) of the cluster assignments, averaged over the runs; 1 means identical
 *     assignments, about 0 an agreement by chance. Only molecules clustered in both assignments are compared.</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class ClusterStability {
    //<editor-fold desc="Record StabilityResult" defaultstate="collapsed">
    /**
     * Result of the stability check.
     *
     * @param clusterJaccardIndices mean Jaccard index of every cluster of the result with its best match, by cluster
     *                              index
     * @param meanAdjustedRandIndex mean adjusted Rand index of the runs
     * @param numberOfRuns number of additional runs
     */
    public record StabilityResult(double[] clusterJaccardIndices, double meanAdjustedRandIndex, int numberOfRuns) {
        /**
         * Returns the mean Jaccard index of the given cluster.
         *
         * @param aClusterIndex cluster index
         * @return mean Jaccard index in [0, 1]
         */
        public double getClusterJaccardIndex(int aClusterIndex) {
            return this.clusterJaccardIndices[aClusterIndex];
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Default number of additional runs.
     */
    public static final int NUMBER_OF_RUNS_DEFAULT = 10;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private ClusterStability() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Clusters the processed features of the given result again with the given number of other seeds (seed of the
     * result + 1, + 2, ...) and compares every run with the result.
     *
     * @param aRunResult clustering result
     * @param aNumberOfRuns number of additional runs, at least 1
     * @param anIsCancelledSupplier returns true if the check should be stopped
     * @param aProgressConsumer receives the number of finished runs and the number of runs
     * @return stability of the clusters
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of runs is smaller than 1
     * @throws CancellationException if the check was cancelled
     * @throws Exception if a clustering fails
     */
    public static StabilityResult check(ClusteringRunResult aRunResult, int aNumberOfRuns, BooleanSupplier anIsCancelledSupplier,
                                        BiConsumer<Integer, Integer> aProgressConsumer)
            throws NullPointerException, IllegalArgumentException, CancellationException, Exception {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        Objects.requireNonNull(aProgressConsumer, "Given progress consumer is null.");
        if (aNumberOfRuns < 1) {
            throw new IllegalArgumentException("The number of runs must be at least 1.");
        }
        ClusteringParameters tmpParameters = aRunResult.parameters();
        int[] tmpBaseLabels = aRunResult.clusteringResult().clusterLabels();
        int tmpBaseClusterCount = aRunResult.clusteringResult().numberOfClusters();
        double[] tmpJaccardSums = new double[tmpBaseClusterCount];
        double tmpRandSum = 0.0;
        for (int r = 0; r < aNumberOfRuns; r++) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The stability check was cancelled.");
            }
            aProgressConsumer.accept(r, aNumberOfRuns);
            ClusteringResult tmpRun = Art2aClusterer.cluster(aRunResult.processedMatrix(), tmpParameters.algorithm(),
                    tmpParameters.vigilance(), tmpParameters.maximumNumberOfClusters(), tmpParameters.maximumNumberOfEpochs(),
                    tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(), tmpParameters.offsetForContrastEnhancement(),
                    tmpParameters.randomSeed() + r + 1, tmpParameters.isParallelCalculation());
            double[] tmpJaccards = ClusterStability.getBestJaccardIndices(tmpBaseLabels, tmpBaseClusterCount,
                    tmpRun.clusterLabels(), tmpRun.numberOfClusters());
            for (int k = 0; k < tmpBaseClusterCount; k++) {
                tmpJaccardSums[k] += tmpJaccards[k];
            }
            tmpRandSum += ClusterStability.getAdjustedRandIndex(tmpBaseLabels, tmpRun.clusterLabels());
        }
        aProgressConsumer.accept(aNumberOfRuns, aNumberOfRuns);
        for (int k = 0; k < tmpBaseClusterCount; k++) {
            tmpJaccardSums[k] /= aNumberOfRuns;
        }
        return new StabilityResult(tmpJaccardSums, tmpRandSum / aNumberOfRuns, aNumberOfRuns);
    }
    //
    /**
     * Returns for every cluster of the first assignment the Jaccard index with the best matching cluster of the second
     * assignment. Molecules with label -1 (not clusterable) belong to no cluster.
     *
     * @param aFirstLabels cluster labels of the first assignment
     * @param aFirstClusterCount number of clusters of the first assignment
     * @param aSecondLabels cluster labels of the second assignment, same length
     * @param aSecondClusterCount number of clusters of the second assignment
     * @return best Jaccard index of every cluster of the first assignment, 0 for an empty cluster
     * @throws IllegalArgumentException if the label arrays have different lengths
     */
    public static double[] getBestJaccardIndices(int[] aFirstLabels, int aFirstClusterCount, int[] aSecondLabels, int aSecondClusterCount)
            throws IllegalArgumentException {
        if (aFirstLabels.length != aSecondLabels.length) {
            throw new IllegalArgumentException("The label arrays have different lengths.");
        }
        int[][] tmpIntersections = new int[aFirstClusterCount][aSecondClusterCount];
        int[] tmpFirstSizes = new int[aFirstClusterCount];
        int[] tmpSecondSizes = new int[aSecondClusterCount];
        for (int i = 0; i < aFirstLabels.length; i++) {
            if (aFirstLabels[i] >= 0) {
                tmpFirstSizes[aFirstLabels[i]]++;
            }
            if (aSecondLabels[i] >= 0) {
                tmpSecondSizes[aSecondLabels[i]]++;
            }
            if (aFirstLabels[i] >= 0 && aSecondLabels[i] >= 0) {
                tmpIntersections[aFirstLabels[i]][aSecondLabels[i]]++;
            }
        }
        double[] tmpBest = new double[aFirstClusterCount];
        for (int k = 0; k < aFirstClusterCount; k++) {
            for (int m = 0; m < aSecondClusterCount; m++) {
                int tmpUnion = tmpFirstSizes[k] + tmpSecondSizes[m] - tmpIntersections[k][m];
                if (tmpUnion > 0) {
                    tmpBest[k] = Math.max(tmpBest[k], (double) tmpIntersections[k][m] / tmpUnion);
                }
            }
        }
        return tmpBest;
    }
    //
    /**
     * Returns the adjusted Rand index of two cluster assignments, compared on the molecules clustered in both (label not
     * -1). The labels need not be numbered alike.
     *
     * @param aFirstLabels cluster labels of the first assignment
     * @param aSecondLabels cluster labels of the second assignment, same length
     * @return adjusted Rand index, 1 for identical assignments; 1 if both assignments put all molecules into one cluster
     * @throws IllegalArgumentException if the label arrays have different lengths
     */
    public static double getAdjustedRandIndex(int[] aFirstLabels, int[] aSecondLabels) throws IllegalArgumentException {
        if (aFirstLabels.length != aSecondLabels.length) {
            throw new IllegalArgumentException("The label arrays have different lengths.");
        }
        int tmpFirstMax = -1;
        int tmpSecondMax = -1;
        for (int i = 0; i < aFirstLabels.length; i++) {
            tmpFirstMax = Math.max(tmpFirstMax, aFirstLabels[i]);
            tmpSecondMax = Math.max(tmpSecondMax, aSecondLabels[i]);
        }
        long[][] tmpContingency = new long[tmpFirstMax + 1][tmpSecondMax + 1];
        long[] tmpRowSums = new long[tmpFirstMax + 1];
        long[] tmpColumnSums = new long[tmpSecondMax + 1];
        long tmpCount = 0;
        for (int i = 0; i < aFirstLabels.length; i++) {
            if (aFirstLabels[i] < 0 || aSecondLabels[i] < 0) {
                continue;
            }
            tmpContingency[aFirstLabels[i]][aSecondLabels[i]]++;
            tmpRowSums[aFirstLabels[i]]++;
            tmpColumnSums[aSecondLabels[i]]++;
            tmpCount++;
        }
        double tmpIndex = 0.0;
        for (long[] tmpRow : tmpContingency) {
            for (long tmpValue : tmpRow) {
                tmpIndex += ClusterStability.pairs(tmpValue);
            }
        }
        double tmpRowPairs = 0.0;
        for (long tmpValue : tmpRowSums) {
            tmpRowPairs += ClusterStability.pairs(tmpValue);
        }
        double tmpColumnPairs = 0.0;
        for (long tmpValue : tmpColumnSums) {
            tmpColumnPairs += ClusterStability.pairs(tmpValue);
        }
        double tmpTotalPairs = ClusterStability.pairs(tmpCount);
        if (tmpTotalPairs == 0.0) {
            return 1.0;
        }
        double tmpExpected = tmpRowPairs * tmpColumnPairs / tmpTotalPairs;
        double tmpMaximum = (tmpRowPairs + tmpColumnPairs) / 2.0;
        return tmpMaximum == tmpExpected ? 1.0 : (tmpIndex - tmpExpected) / (tmpMaximum - tmpExpected);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns the number of pairs of n elements, n * (n - 1) / 2.
     *
     * @param aCount number of elements
     * @return number of pairs
     */
    private static double pairs(long aCount) {
        return aCount * (aCount - 1) / 2.0;
    }
    //</editor-fold>
}
