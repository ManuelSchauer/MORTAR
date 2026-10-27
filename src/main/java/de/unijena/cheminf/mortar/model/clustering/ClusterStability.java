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

import java.util.Arrays;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;

/**
 * Stability of the clusters of a clustering result. The processed features of the result are clustered again several
 * times (all parameters are kept) and every run is compared with the result, in one of two modes:
 * <ul>
 *     <li>{@link StabilityMode#SEEDS}: all molecules with other random seeds. ART-2A depends on the order in which the
 *     molecules are presented, which is shuffled with the seed; the check shows whether a cluster is an artefact of
 *     the order.</li>
 *     <li>{@link StabilityMode#SUBSAMPLES}: a random subsample of {@link #SUBSAMPLE_FRACTION} of the molecules, drawn
 *     without replacement (bootstrap-like subsampling, as for the cluster-wise stability of Hennig), with another seed.
 *     The check shows whether a cluster depends on single molecules of the data set; it is the stricter check, because
 *     the data change, too.</li>
 * </ul>
 * The comparison of a run with the result uses the molecules of the run that are clustered in both:
 * <ul>
 *     <li>Per cluster: the Jaccard index (size of the intersection divided by the size of the union) of its molecules
 *     with the best matching cluster of the run, averaged over the runs in which the cluster has molecules. Near 1 the
 *     cluster is found again in every run, a low value means it falls apart differently from run to run.</li>
 *     <li>Overall: the adjusted Rand index (ARI) of the cluster assignments, averaged over the runs; 1 means identical
 *     assignments, about 0 an agreement by chance.</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class ClusterStability {
    //<editor-fold desc="Enum StabilityMode" defaultstate="collapsed">
    /**
     * What is changed in the additional runs of the stability check.
     */
    public enum StabilityMode {
        /**
         * All molecules, other random seeds (other order of the molecules).
         */
        SEEDS,
        /**
         * Random subsamples of {@link #SUBSAMPLE_FRACTION} of the molecules with other random seeds.
         */
        SUBSAMPLES
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record StabilityResult" defaultstate="collapsed">
    /**
     * Result of the stability check.
     *
     * @param mode what was changed in the additional runs
     * @param clusterJaccardIndices mean Jaccard index of every cluster of the result with its best match, by cluster
     *                              index; NaN for a cluster that had no molecule in any subsample
     * @param meanAdjustedRandIndex mean adjusted Rand index of the runs
     * @param numberOfRuns number of additional runs
     */
    public record StabilityResult(StabilityMode mode, double[] clusterJaccardIndices, double meanAdjustedRandIndex, int numberOfRuns) {
        /**
         * Returns the mean Jaccard index of the given cluster.
         *
         * @param aClusterIndex cluster index
         * @return mean Jaccard index in [0, 1], NaN if the cluster was in no run
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
    /**
     * Maximum number of additional runs that can be chosen.
     */
    public static final int MAXIMUM_NUMBER_OF_RUNS = 100;
    /**
     * Share of the molecules in a subsample of the mode {@link StabilityMode#SUBSAMPLES}.
     */
    public static final double SUBSAMPLE_FRACTION = 0.8;
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
     * Clusters all molecules of the given result again with the given number of other seeds, see
     * {@link #check(ClusteringRunResult, StabilityMode, int, BooleanSupplier, BiConsumer)} with
     * {@link StabilityMode#SEEDS}.
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
        return ClusterStability.check(aRunResult, StabilityMode.SEEDS, aNumberOfRuns, anIsCancelledSupplier, aProgressConsumer);
    }
    //
    /**
     * Clusters the processed features of the given result again the given number of times, with other seeds (seed of
     * the result + 1, + 2, ...) and, in the mode {@link StabilityMode#SUBSAMPLES}, on a random subsample of the
     * molecules drawn with the same seed, and compares every run with the result on the molecules of the run.
     *
     * @param aRunResult clustering result
     * @param aMode what is changed in the additional runs
     * @param aNumberOfRuns number of additional runs, at least 1
     * @param anIsCancelledSupplier returns true if the check should be stopped
     * @param aProgressConsumer receives the number of finished runs and the number of runs
     * @return stability of the clusters
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of runs is smaller than 1
     * @throws CancellationException if the check was cancelled
     * @throws Exception if a clustering fails
     */
    public static StabilityResult check(ClusteringRunResult aRunResult, StabilityMode aMode, int aNumberOfRuns,
                                        BooleanSupplier anIsCancelledSupplier, BiConsumer<Integer, Integer> aProgressConsumer)
            throws NullPointerException, IllegalArgumentException, CancellationException, Exception {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aMode, "Given stability mode is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        Objects.requireNonNull(aProgressConsumer, "Given progress consumer is null.");
        if (aNumberOfRuns < 1) {
            throw new IllegalArgumentException("The number of runs must be at least 1.");
        }
        ClusteringParameters tmpParameters = aRunResult.parameters();
        int[] tmpBaseLabels = aRunResult.clusteringResult().clusterLabels();
        int tmpBaseClusterCount = aRunResult.clusteringResult().numberOfClusters();
        int tmpRowCount = tmpBaseLabels.length;
        int tmpSubsampleSize = (int) Math.round(tmpRowCount * ClusterStability.SUBSAMPLE_FRACTION);
        double[] tmpJaccardSums = new double[tmpBaseClusterCount];
        //runs in which a cluster had molecules; in a subsample a small cluster may be missing
        int[] tmpJaccardCounts = new int[tmpBaseClusterCount];
        double tmpRandSum = 0.0;
        for (int r = 0; r < aNumberOfRuns; r++) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The stability check was cancelled.");
            }
            aProgressConsumer.accept(r, aNumberOfRuns);
            long tmpSeed = tmpParameters.randomSeed() + r + 1;
            //rows of the run: all molecules or a random subsample, ascending
            int[] tmpRows = aMode == StabilityMode.SEEDS ? null
                    : ClusterStability.drawSubsample(tmpRowCount, tmpSubsampleSize, new Random(tmpSeed));
            ClusteringResult tmpRun = Art2aClusterer.cluster(tmpRows == null ? aRunResult.processedMatrix()
                            : aRunResult.processedMatrix().selectRows(tmpRows), tmpParameters.algorithm(),
                    tmpParameters.vigilance(), tmpParameters.maximumNumberOfClusters(), tmpParameters.maximumNumberOfEpochs(),
                    tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(), tmpParameters.offsetForContrastEnhancement(),
                    tmpSeed, tmpParameters.isParallelCalculation());
            //labels of the result for the molecules of the run
            int[] tmpBaseRunLabels = tmpRows == null ? tmpBaseLabels : Arrays.stream(tmpRows).map(i -> tmpBaseLabels[i]).toArray();
            boolean[] tmpIsPresent = new boolean[tmpBaseClusterCount];
            for (int tmpLabel : tmpBaseRunLabels) {
                if (tmpLabel >= 0) {
                    tmpIsPresent[tmpLabel] = true;
                }
            }
            double[] tmpJaccards = ClusterStability.getBestJaccardIndices(tmpBaseRunLabels, tmpBaseClusterCount,
                    tmpRun.clusterLabels(), tmpRun.numberOfClusters());
            for (int k = 0; k < tmpBaseClusterCount; k++) {
                if (tmpIsPresent[k]) {
                    tmpJaccardSums[k] += tmpJaccards[k];
                    tmpJaccardCounts[k]++;
                }
            }
            tmpRandSum += ClusterStability.getAdjustedRandIndex(tmpBaseRunLabels, tmpRun.clusterLabels());
        }
        aProgressConsumer.accept(aNumberOfRuns, aNumberOfRuns);
        for (int k = 0; k < tmpBaseClusterCount; k++) {
            tmpJaccardSums[k] = tmpJaccardCounts[k] == 0 ? Double.NaN : tmpJaccardSums[k] / tmpJaccardCounts[k];
        }
        return new StabilityResult(aMode, tmpJaccardSums, tmpRandSum / aNumberOfRuns, aNumberOfRuns);
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
     * Draws a random subsample of the given size from the row indices 0 to the row count - 1, without replacement.
     *
     * @param aRowCount number of rows
     * @param aSubsampleSize number of rows to draw
     * @param aRandom random generator
     * @return drawn row indices, ascending
     */
    private static int[] drawSubsample(int aRowCount, int aSubsampleSize, Random aRandom) {
        int[] tmpRows = new int[aRowCount];
        for (int i = 0; i < aRowCount; i++) {
            tmpRows[i] = i;
        }
        //partial Fisher-Yates shuffle: the first rows become the subsample
        for (int i = 0; i < aSubsampleSize; i++) {
            int tmpOther = i + aRandom.nextInt(aRowCount - i);
            int tmpRow = tmpRows[i];
            tmpRows[i] = tmpRows[tmpOther];
            tmpRows[tmpOther] = tmpRow;
        }
        int[] tmpSubsample = Arrays.copyOf(tmpRows, aSubsampleSize);
        Arrays.sort(tmpSubsample);
        return tmpSubsample;
    }
    //
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
