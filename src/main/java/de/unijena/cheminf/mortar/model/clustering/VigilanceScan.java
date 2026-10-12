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

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.CancellationException;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;

/**
 * Vigilance scan: clusters the processed feature matrix of an existing clustering run again for a series of vigilance
 * values and reports for each value the number of clusters, the number of single-molecule clusters, the largest cluster
 * and the mean similarity of the molecules to their cluster vector. The scan helps to choose the vigilance, which has no
 * general default. The features (fragment fingerprint, descriptors, transformation, filter) do not depend on the
 * vigilance and are taken from the given run, so all values of the scan are comparable and only the clustering itself is
 * repeated; all other parameters (kernel, epochs, learning parameter, seed, ...) are kept.
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
     * @param isOverflow true if the maximum number of clusters was reached
     * @param isConverged true if the clustering converged
     * @param runResult complete clustering result with this vigilance, e.g. to show it as result tab
     */
    public record ScanEntry(float vigilance, int numberOfClusters, int numberOfSingletons, int largestClusterSize,
                            int numberOfUnclusteredMolecules, double meanSimilarityToClusterVector, boolean isOverflow,
                            boolean isConverged, ClusteringRunResult runResult) {
    }
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
     * are kept.
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
        ClusteringParameters tmpParameters = aRunResult.parameters();
        List<ScanEntry> tmpEntries = new ArrayList<>(aVigilanceValues.length);
        for (int v = 0; v < aVigilanceValues.length; v++) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("Vigilance scan was cancelled.");
            }
            aProgressConsumer.accept(v, aVigilanceValues.length);
            float tmpVigilance = aVigilanceValues[v];
            ClusteringResult tmpResult = Art2aClusterer.cluster(aRunResult.processedMatrix(), tmpParameters.algorithm(),
                    tmpVigilance, tmpParameters.maximumNumberOfClusters(), tmpParameters.maximumNumberOfEpochs(),
                    tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(),
                    tmpParameters.offsetForContrastEnhancement(), tmpParameters.randomSeed(), tmpParameters.isParallelCalculation());
            ClusteringParameters tmpScanParameters = VigilanceScan.withVigilance(tmpParameters, tmpVigilance);
            ClusteringRunResult tmpRunResult = new ClusteringRunResult(tmpScanParameters, aRunResult.molecules(), aRunResult.rawMatrix(),
                    aRunResult.transformedMatrix(), aRunResult.processedMatrix(), aRunResult.removedFeatures(), tmpResult);
            tmpEntries.add(VigilanceScan.createEntry(tmpRunResult));
        }
        aProgressConsumer.accept(aVigilanceValues.length, aVigilanceValues.length);
        return tmpEntries;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
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
        return new ScanEntry(tmpResult.vigilance(), tmpResult.numberOfClusters(), tmpSingletons, tmpLargest,
                tmpResult.getNumberOfUnclusteredRows(), tmpClustered == 0 ? Double.NaN : tmpSimilaritySum / tmpClustered,
                tmpResult.isOverflow(), tmpResult.isConverged(), aRunResult);
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
                aParameters.randomSeed(), aParameters.circularFingerprintSize(), aParameters.isParallelCalculation());
    }
    //</editor-fold>
}
