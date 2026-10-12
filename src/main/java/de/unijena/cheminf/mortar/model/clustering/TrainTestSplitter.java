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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Random;

/**
 * Creates a training/test partitioning of the clustered molecules based on a clustering result, so that both sets cover
 * the chemical space of the data in the same way (cluster-stratified split). The rules follow the training/test split of
 * the ART2a-Clustering-for-Java library (Art2aResult.getTrainingAndTestIndices) and are available for both kernels:
 * <ul>
 *     <li>A cluster with one molecule goes completely into the training set.</li>
 *     <li>A cluster with n &ge; 2 molecules contributes round(n * training fraction) molecules to the training set,
 *     but at least one and at most n - 1, so that every such cluster is represented in both sets.</li>
 *     <li>Molecules that could not be clustered go into the training set, like single-molecule clusters.</li>
 * </ul>
 * Because of these rules, the actual training fraction can deviate from the requested one, especially with many small
 * clusters; it is part of the result.
 * <p>
 * Which molecules of a cluster go into the training set is decided by the selection mode: the molecules closest to the
 * cluster vector (as in the library; the test set then contains the less typical molecules) or a random selection
 * with a fixed seed.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class TrainTestSplitter {
    //<editor-fold desc="Enum SelectionMode" defaultstate="collapsed">
    /**
     * Enumeration of the ways to select the training molecules of a cluster.
     */
    public enum SelectionMode {
        /**
         * The molecules closest to the cluster vector go into the training set (as in the ART2a library).
         */
        CLOSEST_TO_CLUSTER_CENTER,
        /**
         * The training molecules of every cluster are selected randomly with a fixed seed.
         */
        RANDOM;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record TrainTestSplit" defaultstate="collapsed">
    /**
     * Result of a training/test partitioning.
     *
     * @param trainingIndices row indices of the training molecules in the clustered data, ascending
     * @param testIndices row indices of the test molecules in the clustered data, ascending
     * @param requestedTrainingFraction requested fraction of training molecules
     * @param selectionMode way the training molecules of a cluster were selected
     */
    public record TrainTestSplit(int[] trainingIndices, int[] testIndices, double requestedTrainingFraction, SelectionMode selectionMode) {
        /**
         * Returns the actual fraction of training molecules.
         *
         * @return number of training molecules divided by the number of all molecules
         */
        public double getActualTrainingFraction() {
            int tmpTotal = this.trainingIndices.length + this.testIndices.length;
            return tmpTotal == 0 ? 0.0 : (double) this.trainingIndices.length / tmpTotal;
        }
        //
        /**
         * Returns whether the molecule with the given row index is in the training set.
         *
         * @param aRowIndex row index
         * @return true for training, false for test
         */
        public boolean isTraining(int aRowIndex) {
            return Arrays.binarySearch(this.trainingIndices, aRowIndex) >= 0;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private TrainTestSplitter() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Creates the cluster-based training/test partitioning, see the class description for the rules.
     *
     * @param aRunResult clustering result
     * @param aTrainingFraction requested fraction of training molecules in (0,1)
     * @param aSelectionMode way to select the training molecules of a cluster
     * @param aRandomSeed seed for {@link SelectionMode#RANDOM}
     * @return training/test partitioning
     * @throws NullPointerException if the result or mode is null
     * @throws IllegalArgumentException if the fraction is not in (0,1)
     */
    public static TrainTestSplit split(ClusteringRunResult aRunResult, double aTrainingFraction, SelectionMode aSelectionMode,
                                       long aRandomSeed) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aSelectionMode, "Given selection mode is null.");
        if (!(aTrainingFraction > 0.0 && aTrainingFraction < 1.0)) {
            throw new IllegalArgumentException("The training fraction must be in (0,1).");
        }
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        int tmpNumberOfClusters = aRunResult.clusteringResult().numberOfClusters();
        List<List<Integer>> tmpMembers = new ArrayList<>(tmpNumberOfClusters);
        for (int k = 0; k < tmpNumberOfClusters; k++) {
            tmpMembers.add(new ArrayList<>());
        }
        List<Integer> tmpTraining = new ArrayList<>();
        List<Integer> tmpTest = new ArrayList<>();
        for (int i = 0; i < tmpLabels.length; i++) {
            if (tmpLabels[i] < 0) {
                //not clusterable molecules are treated like single-molecule clusters
                tmpTraining.add(i);
            } else {
                tmpMembers.get(tmpLabels[i]).add(i);
            }
        }
        SimilaritySearch tmpSearch = aSelectionMode == SelectionMode.CLOSEST_TO_CLUSTER_CENTER ? new SimilaritySearch(aRunResult) : null;
        Random tmpRandom = new Random(aRandomSeed);
        for (List<Integer> tmpCluster : tmpMembers) {
            int tmpSize = tmpCluster.size();
            if (tmpSize == 0) {
                continue;
            }
            if (tmpSize == 1) {
                tmpTraining.add(tmpCluster.getFirst());
                continue;
            }
            if (tmpSearch != null) {
                //most typical molecules first, ties by row index
                tmpCluster.sort(Comparator.comparingDouble((Integer anIndex) -> -tmpSearch.getSimilarityToOwnClusterVector(anIndex))
                        .thenComparingInt(anIndex -> anIndex));
            } else {
                Collections.shuffle(tmpCluster, tmpRandom);
            }
            int tmpNumberOfTraining = (int) Math.max(1, Math.min(tmpSize - 1, Math.round(tmpSize * aTrainingFraction)));
            tmpTraining.addAll(tmpCluster.subList(0, tmpNumberOfTraining));
            tmpTest.addAll(tmpCluster.subList(tmpNumberOfTraining, tmpSize));
        }
        int[] tmpTrainingIndices = tmpTraining.stream().mapToInt(Integer::intValue).sorted().toArray();
        int[] tmpTestIndices = tmpTest.stream().mapToInt(Integer::intValue).sorted().toArray();
        return new TrainTestSplit(tmpTrainingIndices, tmpTestIndices, aTrainingFraction, aSelectionMode);
    }
    //
    /**
     * Parses a training/test ratio: "40:60" (training:test, any positive numbers), "40" (percent of training) or "0.4"
     * (fraction of training).
     *
     * @param aRatio ratio text
     * @return training fraction
     * @throws IllegalArgumentException if the text cannot be parsed or the fraction is not in (0,1)
     */
    public static double parseTrainingFraction(String aRatio) throws IllegalArgumentException {
        if (aRatio == null || aRatio.isBlank()) {
            throw new IllegalArgumentException("No training/test ratio given.");
        }
        double tmpFraction;
        try {
            String tmpRatio = aRatio.trim().replace(',', '.').replace("%", "");
            if (tmpRatio.contains(":") || tmpRatio.contains("/")) {
                String[] tmpParts = tmpRatio.split("[:/]");
                if (tmpParts.length != 2) {
                    throw new IllegalArgumentException("A ratio must have the form training:test, e.g. 80:20.");
                }
                double tmpTraining = Double.parseDouble(tmpParts[0].trim());
                double tmpTest = Double.parseDouble(tmpParts[1].trim());
                if (tmpTraining < 0 || tmpTest < 0 || tmpTraining + tmpTest == 0) {
                    throw new IllegalArgumentException("The parts of the ratio must not be negative.");
                }
                tmpFraction = tmpTraining / (tmpTraining + tmpTest);
            } else {
                double tmpValue = Double.parseDouble(tmpRatio);
                tmpFraction = tmpValue >= 1.0 ? tmpValue / 100.0 : tmpValue;
            }
        } catch (NumberFormatException anException) {
            throw new IllegalArgumentException("The training/test ratio could not be read: " + aRatio, anException);
        }
        if (!(tmpFraction > 0.0 && tmpFraction < 1.0)) {
            throw new IllegalArgumentException("The training fraction must be between 0 and 1 (exclusive): " + aRatio);
        }
        return tmpFraction;
    }
    //</editor-fold>
}
