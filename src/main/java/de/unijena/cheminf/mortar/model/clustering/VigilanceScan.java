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

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.concurrent.CancellationException;
import java.util.concurrent.CompletionException;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;
import java.util.stream.IntStream;

/**
 * Vigilance scan: clusters the processed feature matrix of an existing clustering run again for a series of vigilance
 * values and reports for each value the number of clusters, the number of single-molecule clusters, the largest cluster
 * and the mean similarity of the molecules to their cluster vector. The scan helps to choose the vigilance, which has no
 * general default. The features (fragment fingerprint, descriptors, transformation, filter) do not depend on the
 * vigilance and are taken from the given run, so all values of the scan are comparable and only the clustering itself is
 * repeated; all other parameters (kernel, epochs, learning parameter, seed, ...) are kept.
 * <p>
 * Instead of a range, a target number of clusters can be given ({@link #searchForNumberOfClusters}): the vigilance is
 * searched with two decimals by bisection from both ends of (0,1), see the method description.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class VigilanceScan {
    //<editor-fold desc="Record ScanEntry" defaultstate="collapsed">
    /**
     * Result of the clustering with one vigilance value.
     *
     * @param vigilance vigilance value
     * @param numberOfClusters number of clusters
     * @param numberOfSingletons number of clusters with one molecule
     * @param largestClusterSize size of the largest cluster
     * @param numberOfUnclusteredMolecules number of molecules that could not be clustered
     * @param meanSimilarityToClusterVector mean similarity of the clustered molecules to their cluster vector (compactness
     *                                      of the clusters, 1 = every molecule equals its cluster vector)
     * @param overallSilhouette mean silhouette of the molecules (see {@link ClusterQuality}), NaN with fewer than two
     *                          clusters
     * @param overallNeighbourPurity mean neighbour purity of the molecules (see {@link ClusterQuality})
     * @param isOverflow true if the maximum number of clusters was reached
     * @param isConverged true if the clustering converged
     * @param runResult complete clustering result with this vigilance, e.g. to show it as result tab
     */
    public record ScanEntry(float vigilance, int numberOfClusters, int numberOfSingletons, int largestClusterSize,
                            int numberOfUnclusteredMolecules, double meanSimilarityToClusterVector, double overallSilhouette,
                            double overallNeighbourPurity, boolean isOverflow, boolean isConverged, ClusteringRunResult runResult) {
    }
    //
    /**
     * Result of the search for a target number of clusters.
     *
     * @param entries entries of all vigilance values clustered during the search, in ascending order of the vigilance
     * @param bestEntry entry whose number of clusters is closest to the target
     * @param targetNumberOfClusters requested number of clusters
     */
    public record TargetSearchResult(List<ScanEntry> entries, ScanEntry bestEntry, int targetNumberOfClusters) {
        /**
         * Returns whether the best entry has exactly the requested number of clusters.
         *
         * @return true if the target was hit exactly
         */
        public boolean isExact() {
            return this.bestEntry.numberOfClusters() == this.targetNumberOfClusters;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Number of steps of the vigilance grid of the target search: the vigilance is searched with two decimals, i.e. in
     * steps of 1 / 100.
     */
    public static final int SEARCH_GRID_STEPS = 100;
    /**
     * Upper bound of the number of clusterings of a target search: the two ends plus the bisection of 97 grid points.
     */
    public static final int MAXIMUM_NUMBER_OF_SEARCH_STEPS = 9;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private VigilanceScan() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns the vigilance values from the given start to the given end (both included) with the given step size,
     * rounded to 4 decimals to avoid floating point artefacts like 0.30000001.
     *
     * @param aStart first vigilance value in (0,1)
     * @param anEnd last vigilance value in (0,1), not smaller than the start
     * @param aStep step size, greater than 0
     * @return vigilance values in ascending order
     * @throws IllegalArgumentException if a parameter is invalid or more than 1000 values would result
     */
    public static float[] getVigilanceValues(double aStart, double anEnd, double aStep) throws IllegalArgumentException {
        if (!(aStart > 0.0 && aStart < 1.0) || !(anEnd > 0.0 && anEnd < 1.0) || anEnd < aStart) {
            throw new IllegalArgumentException("Start and end must be in (0,1) and the end must not be smaller than the start.");
        }
        if (!(aStep > 0.0) || (anEnd - aStart) / aStep > 1000) {
            throw new IllegalArgumentException("The step size must be greater than 0 and result in at most 1000 values.");
        }
        List<Float> tmpValues = new ArrayList<>();
        //a small tolerance includes the end despite rounding errors
        for (int i = 0; aStart + i * aStep <= anEnd + 1.0e-9; i++) {
            tmpValues.add((float) (Math.round((aStart + i * aStep) * 10000.0) / 10000.0));
        }
        float[] tmpArray = new float[tmpValues.size()];
        for (int i = 0; i < tmpArray.length; i++) {
            tmpArray[i] = tmpValues.get(i);
        }
        return tmpArray;
    }
    //
    /**
     * Clusters the processed matrix of the given run with every given vigilance value; the other parameters of the run
     * are kept. The clusterings are independent of each other and run in parallel if the parallel calculation of the run
     * is switched on; the result does not depend on it. The silhouette and neighbour purity of all clusterings are
     * calculated at the end in one pass over the distances of the molecules (see {@link ClusterQuality#getOverallValues}).
     *
     * @param aRunResult clustering run whose features and parameters are used
     * @param aVigilanceValues vigilance values, each in (0,1)
     * @param anIsCancelledSupplier returns true if the scan should be stopped
     * @param aProgressConsumer receives the number of finished values and the total number
     * @return one entry per vigilance value in the given order
     * @throws NullPointerException if a parameter is null
     * @throws CancellationException if the scan was cancelled
     * @throws Exception if a clustering fails
     */
    public static List<ScanEntry> scan(ClusteringRunResult aRunResult, float[] aVigilanceValues, BooleanSupplier anIsCancelledSupplier,
                                       BiConsumer<Integer, Integer> aProgressConsumer)
            throws NullPointerException, CancellationException, Exception {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aVigilanceValues, "Given vigilance values are null.");
        ScanEntry[] tmpEntries = new ScanEntry[aVigilanceValues.length];
        AtomicInteger tmpFinishedCount = new AtomicInteger();
        aProgressConsumer.accept(0, aVigilanceValues.length);
        IntStream tmpIndices = IntStream.range(0, aVigilanceValues.length);
        try {
            (aRunResult.parameters().isParallelCalculation() ? tmpIndices.parallel() : tmpIndices).forEach(v -> {
                if (anIsCancelledSupplier.getAsBoolean()) {
                    throw new CancellationException("Vigilance scan was cancelled.");
                }
                try {
                    tmpEntries[v] = VigilanceScan.clusterWithVigilance(aRunResult, aVigilanceValues[v]);
                } catch (Exception anException) {
                    throw new CompletionException(anException);
                }
                aProgressConsumer.accept(tmpFinishedCount.incrementAndGet(), aVigilanceValues.length);
            });
        } catch (CompletionException anException) {
            //the exception of a clustering, wrapped to leave the stream
            throw anException.getCause() instanceof Exception tmpCause ? tmpCause : anException;
        }
        //the last clustering has already reported the end of the progress
        return VigilanceScan.addQuality(aRunResult, List.of(tmpEntries), anIsCancelledSupplier);
    }
    //
    /**
     * Formats a vigilance with as many decimals as it has, at most four (the precision of the range scan): 0.43, 0.5 or
     * 0.125; values of the target search have at most two decimals.
     *
     * @param aVigilance vigilance
     * @return formatted vigilance with a decimal point
     */
    public static String formatVigilance(float aVigilance) {
        //the float is converted via its decimal representation, e.g. 0.3 instead of 0.30000001192092896
        BigDecimal tmpValue = BigDecimal.valueOf(ClusteringService.toDecimalDouble(aVigilance)).setScale(4, RoundingMode.HALF_UP)
                .stripTrailingZeros();
        return tmpValue.scale() < 1 ? tmpValue.setScale(1, RoundingMode.UNNECESSARY).toPlainString() : tmpValue.toPlainString();
    }
    //
    /**
     * Searches the vigilance (two decimals, 0.01 to 0.99) whose clustering has the given number of clusters. The number
     * of clusters rises with the vigilance, so the search narrows the vigilance from both sides (bisection):
     * <ol>
     *     <li>Both ends 0.01 and 0.99 are clustered. If the target is not above the number of clusters at 0.01 or not
     *     below the one at 0.99, the search ends there.</li>
     *     <li>Otherwise, the middle of the interval between the lower bound (fewer clusters than the target) and the upper
     *     bound (more clusters) is clustered and replaces the bound on its side, until the target is hit exactly or no
     *     vigilance with two decimals lies between the bounds.</li>
     * </ol>
     * Since ART-2A is not strictly monotone, the best entry is chosen among all clustered values: the smallest deviation
     * from the target, then the fewest singletons, then the highest mean similarity to the cluster vector. At most
     * {@link #MAXIMUM_NUMBER_OF_SEARCH_STEPS} clusterings are needed.
     *
     * @param aRunResult clustering run whose features and parameters are used
     * @param aTargetNumberOfClusters requested number of clusters, at least 1
     * @param anIsCancelledSupplier returns true if the search should be stopped
     * @param aProgressConsumer receives the number of finished clusterings and the maximum number
     * @return clustered entries and the best one
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the target is smaller than 1
     * @throws CancellationException if the search was cancelled
     * @throws Exception if a clustering fails
     */
    public static TargetSearchResult searchForNumberOfClusters(ClusteringRunResult aRunResult, int aTargetNumberOfClusters,
                                                               BooleanSupplier anIsCancelledSupplier, BiConsumer<Integer, Integer> aProgressConsumer)
            throws NullPointerException, IllegalArgumentException, CancellationException, Exception {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        Objects.requireNonNull(aProgressConsumer, "Given progress consumer is null.");
        if (aTargetNumberOfClusters < 1) {
            throw new IllegalArgumentException("The target number of clusters must be at least 1.");
        }
        //clustered grid points (vigilance = index / SEARCH_GRID_STEPS) in ascending order
        TreeMap<Integer, ScanEntry> tmpEntryOfGridPointMap = new TreeMap<>();
        int tmpLower = 1;
        int tmpUpper = VigilanceScan.SEARCH_GRID_STEPS - 1;
        //both ends first, at the same time if the calculation is parallel; if the target lies outside, the search ends at
        // the corresponding end
        if (aRunResult.parameters().isParallelCalculation()) {
            VigilanceScan.evaluateGridPointsInParallel(aRunResult, new int[] {tmpLower, tmpUpper}, tmpEntryOfGridPointMap,
                    anIsCancelledSupplier, aProgressConsumer);
        }
        if (VigilanceScan.evaluateGridPoint(aRunResult, tmpLower, tmpEntryOfGridPointMap, anIsCancelledSupplier, aProgressConsumer)
                .numberOfClusters() < aTargetNumberOfClusters
                && VigilanceScan.evaluateGridPoint(aRunResult, tmpUpper, tmpEntryOfGridPointMap, anIsCancelledSupplier, aProgressConsumer)
                .numberOfClusters() > aTargetNumberOfClusters) {
            //bisection: the lower bound has fewer, the upper bound more clusters than the target
            while (tmpUpper - tmpLower > 1) {
                int tmpMiddle = (tmpLower + tmpUpper) / 2;
                int tmpNumberOfClusters = VigilanceScan.evaluateGridPoint(aRunResult, tmpMiddle, tmpEntryOfGridPointMap,
                        anIsCancelledSupplier, aProgressConsumer).numberOfClusters();
                if (tmpNumberOfClusters == aTargetNumberOfClusters) {
                    break;
                } else if (tmpNumberOfClusters < aTargetNumberOfClusters) {
                    tmpLower = tmpMiddle;
                } else {
                    tmpUpper = tmpMiddle;
                }
            }
        }
        //silhouette and neighbour purity of all clusterings of the search in one pass
        List<ScanEntry> tmpEntries = VigilanceScan.addQuality(aRunResult, new ArrayList<>(tmpEntryOfGridPointMap.values()),
                anIsCancelledSupplier);
        aProgressConsumer.accept(VigilanceScan.MAXIMUM_NUMBER_OF_SEARCH_STEPS, VigilanceScan.MAXIMUM_NUMBER_OF_SEARCH_STEPS);
        ScanEntry tmpBest = tmpEntries.stream().min(Comparator
                .comparingInt((ScanEntry anEntry) -> Math.abs(anEntry.numberOfClusters() - aTargetNumberOfClusters))
                .thenComparingInt(ScanEntry::numberOfSingletons)
                .thenComparing(Comparator.comparingDouble(ScanEntry::meanSimilarityToClusterVector).reversed())).orElseThrow();
        return new TargetSearchResult(List.copyOf(tmpEntries), tmpBest, aTargetNumberOfClusters);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Clusters the given grid point of the target search (vigilance = grid point / {@link #SEARCH_GRID_STEPS}) unless
     * it was already clustered, and stores its entry.
     *
     * @param aRunResult clustering run whose features and parameters are used
     * @param aGridPoint grid point between 1 and SEARCH_GRID_STEPS - 1
     * @param anEntryOfGridPointMap entries of the grid points clustered so far (IS CHANGED)
     * @param anIsCancelledSupplier returns true if the search should be stopped
     * @param aProgressConsumer receives the number of finished clusterings and the maximum number
     * @return entry of the grid point
     * @throws CancellationException if the search was cancelled
     * @throws Exception if the clustering fails
     */
    private static ScanEntry evaluateGridPoint(ClusteringRunResult aRunResult, int aGridPoint, Map<Integer, ScanEntry> anEntryOfGridPointMap,
                                               BooleanSupplier anIsCancelledSupplier, BiConsumer<Integer, Integer> aProgressConsumer)
            throws CancellationException, Exception {
        ScanEntry tmpEntry = anEntryOfGridPointMap.get(aGridPoint);
        if (tmpEntry != null) {
            return tmpEntry;
        }
        if (anIsCancelledSupplier.getAsBoolean()) {
            throw new CancellationException("Vigilance search was cancelled.");
        }
        aProgressConsumer.accept(anEntryOfGridPointMap.size(), VigilanceScan.MAXIMUM_NUMBER_OF_SEARCH_STEPS);
        tmpEntry = VigilanceScan.clusterWithVigilance(aRunResult, (float) aGridPoint / VigilanceScan.SEARCH_GRID_STEPS);
        anEntryOfGridPointMap.put(aGridPoint, tmpEntry);
        return tmpEntry;
    }
    //
    /**
     * Clusters the given grid points of the target search at the same time and stores their entries, see
     * {@link #evaluateGridPoint}.
     *
     * @param aRunResult clustering run whose features and parameters are used
     * @param aGridPoints grid points between 1 and SEARCH_GRID_STEPS - 1
     * @param anEntryOfGridPointMap entries of the grid points clustered so far (IS CHANGED)
     * @param anIsCancelledSupplier returns true if the search should be stopped
     * @param aProgressConsumer receives the number of finished clusterings and the maximum number
     * @throws CancellationException if the search was cancelled
     * @throws Exception if a clustering fails
     */
    private static void evaluateGridPointsInParallel(ClusteringRunResult aRunResult, int[] aGridPoints, Map<Integer, ScanEntry> anEntryOfGridPointMap,
                                                     BooleanSupplier anIsCancelledSupplier, BiConsumer<Integer, Integer> aProgressConsumer)
            throws CancellationException, Exception {
        if (anIsCancelledSupplier.getAsBoolean()) {
            throw new CancellationException("Vigilance search was cancelled.");
        }
        aProgressConsumer.accept(anEntryOfGridPointMap.size(), VigilanceScan.MAXIMUM_NUMBER_OF_SEARCH_STEPS);
        ScanEntry[] tmpEntries = new ScanEntry[aGridPoints.length];
        try {
            IntStream.range(0, aGridPoints.length).parallel().forEach(g -> {
                try {
                    tmpEntries[g] = VigilanceScan.clusterWithVigilance(aRunResult, (float) aGridPoints[g] / VigilanceScan.SEARCH_GRID_STEPS);
                } catch (Exception anException) {
                    throw new CompletionException(anException);
                }
            });
        } catch (CompletionException anException) {
            throw anException.getCause() instanceof Exception tmpCause ? tmpCause : anException;
        }
        //the map is only changed in this thread
        for (int g = 0; g < aGridPoints.length; g++) {
            anEntryOfGridPointMap.put(aGridPoints[g], tmpEntries[g]);
        }
    }
    //
    /**
     * Clusters the processed matrix of the given run with the given vigilance; the other parameters of the run are kept.
     *
     * @param aRunResult clustering run whose features and parameters are used
     * @param aVigilance vigilance in (0,1)
     * @return scan entry of the clustering
     * @throws Exception if the clustering fails
     */
    private static ScanEntry clusterWithVigilance(ClusteringRunResult aRunResult, float aVigilance) throws Exception {
        ClusteringParameters tmpParameters = aRunResult.parameters();
        ClusteringResult tmpResult = Art2aClusterer.cluster(aRunResult.processedMatrix(), tmpParameters.algorithm(),
                aVigilance, tmpParameters.maximumNumberOfClusters(), tmpParameters.maximumNumberOfEpochs(),
                tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(),
                tmpParameters.offsetForContrastEnhancement(), tmpParameters.randomSeed(), tmpParameters.isParallelCalculation());
        ClusteringParameters tmpScanParameters = VigilanceScan.withVigilance(tmpParameters, aVigilance);
        ClusteringRunResult tmpRunResult = new ClusteringRunResult(tmpScanParameters, aRunResult.molecules(), aRunResult.rawMatrix(),
                aRunResult.transformedMatrix(), aRunResult.processedMatrix(), aRunResult.removedFeatures(), tmpResult);
        return VigilanceScan.createEntry(tmpRunResult);
    }
    //
    /**
     * Creates the scan entry of the given clustering result.
     *
     * @param aRunResult clustering result
     * @return scan entry
     */
    private static ScanEntry createEntry(ClusteringRunResult aRunResult) {
        ClusteringResult tmpResult = aRunResult.clusteringResult();
        int[] tmpSizes = new int[tmpResult.numberOfClusters()];
        for (int tmpLabel : tmpResult.clusterLabels()) {
            if (tmpLabel >= 0) {
                tmpSizes[tmpLabel]++;
            }
        }
        int tmpSingletons = 0;
        int tmpLargest = 0;
        for (int tmpSize : tmpSizes) {
            if (tmpSize == 1) {
                tmpSingletons++;
            }
            tmpLargest = Math.max(tmpLargest, tmpSize);
        }
        SimilaritySearch tmpSearch = new SimilaritySearch(aRunResult);
        double tmpSimilaritySum = 0.0;
        int tmpClustered = 0;
        for (int i = 0; i < tmpResult.clusterLabels().length; i++) {
            if (tmpResult.clusterLabels()[i] >= 0) {
                tmpSimilaritySum += tmpSearch.getSimilarityToOwnClusterVector(i);
                tmpClustered++;
            }
        }
        //silhouette and neighbour purity are added for all entries at once, see addQuality
        return new ScanEntry(tmpResult.vigilance(), tmpResult.numberOfClusters(), tmpSingletons, tmpLargest,
                tmpResult.getNumberOfUnclusteredRows(), tmpClustered == 0 ? Double.NaN : tmpSimilaritySum / tmpClustered,
                Double.NaN, Double.NaN, tmpResult.isOverflow(), tmpResult.isConverged(), aRunResult);
    }
    //
    /**
     * Returns copies of the given entries with their overall silhouette and neighbour purity, calculated for all
     * entries in one pass over the distances of the molecules of the given run, whose features they share.
     *
     * @param aRunResult clustering run whose features the entries share
     * @param anEntryList entries without quality
     * @param anIsCancelledSupplier returns true if the calculation should be stopped
     * @return entries with quality, in the same order
     * @throws CancellationException if the calculation was cancelled
     */
    private static List<ScanEntry> addQuality(ClusteringRunResult aRunResult, List<ScanEntry> anEntryList, BooleanSupplier anIsCancelledSupplier)
            throws CancellationException {
        List<int[]> tmpLabelArrays = new ArrayList<>(anEntryList.size());
        for (ScanEntry tmpEntry : anEntryList) {
            tmpLabelArrays.add(tmpEntry.runResult().clusteringResult().clusterLabels());
        }
        double[][] tmpValues = ClusterQuality.getOverallValues(aRunResult, tmpLabelArrays, anIsCancelledSupplier);
        List<ScanEntry> tmpEntries = new ArrayList<>(anEntryList.size());
        for (int i = 0; i < anEntryList.size(); i++) {
            ScanEntry tmpEntry = anEntryList.get(i);
            tmpEntries.add(new ScanEntry(tmpEntry.vigilance(), tmpEntry.numberOfClusters(), tmpEntry.numberOfSingletons(),
                    tmpEntry.largestClusterSize(), tmpEntry.numberOfUnclusteredMolecules(), tmpEntry.meanSimilarityToClusterVector(),
                    tmpValues[i][0], tmpValues[i][1], tmpEntry.isOverflow(), tmpEntry.isConverged(), tmpEntry.runResult()));
        }
        return tmpEntries;
    }
    //
    /**
     * Returns a copy of the given parameters with another vigilance.
     *
     * @param aParameters parameters
     * @param aVigilance new vigilance
     * @return parameters with the new vigilance
     */
    private static ClusteringParameters withVigilance(ClusteringParameters aParameters, float aVigilance) {
        return new ClusteringParameters(aParameters.fragmentationName(), aParameters.isCountFingerprint(),
                aParameters.minimumFragmentMoleculeFrequency(), aParameters.maximumNumberOfFragmentFeatures(),
                aParameters.descriptors(), aParameters.electronDonationModelOption(), aParameters.cycleFinderOption(),
                aParameters.countTransformation(), aParameters.continuousTransformation(), aParameters.maximumNaNFraction(),
                aParameters.isCorrelationFilterEnabled(), aParameters.correlationThreshold(), aParameters.correlationMethod(),
                aParameters.algorithm(), aVigilance, aParameters.maximumNumberOfClusters(), aParameters.maximumNumberOfEpochs(),
                aParameters.convergenceThreshold(), aParameters.learningParameter(), aParameters.offsetForContrastEnhancement(),
                aParameters.randomSeed(), aParameters.circularFingerprintSize(), aParameters.isParallelCalculation(),
                aParameters.isFragmentFingerprintUsed());
    }
    //</editor-fold>
}
