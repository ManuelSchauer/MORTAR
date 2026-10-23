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
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;

/**
 * Creates cluster profiles: for every cluster, the features that distinguish it from the whole data set, so that a
 * cluster can be interpreted ("92 % of the molecules contain a glycoside, compared to 12 % overall"; "TPSA is 1.8
 * standard deviations above the mean"). The profiles are computed on the raw feature matrix, so that the values are in
 * the original units.
 * <ul>
 *     <li>Fragments and fingerprint bits (presence features): share of the molecules that contain the feature (value
 *     &gt; 0) in the cluster and in the whole data set; the score is the difference of the shares (-1 to 1).</li>
 *     <li>Descriptors (continuous features): mean in the cluster and in the whole data set; the score is the difference
 *     of the means divided by the standard deviation in the whole data set. NaN values are ignored.</li>
 * </ul>
 * Molecules that could not be clustered belong to no cluster but are part of the whole data set.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusterProfiler {
    //<editor-fold desc="Enum FeatureCategory" defaultstate="collapsed">
    /**
     * Enumeration of the feature categories of a profile.
     */
    public enum FeatureCategory {
        /**
         * Fragment of the fragmentation (bit or count), profiled by its presence.
         */
        FRAGMENT,
        /**
         * Descriptor component, profiled by its mean.
         */
        DESCRIPTOR,
        /**
         * Bit of a CDK fingerprint, profiled by its presence.
         */
        FINGERPRINT_BIT;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Records" defaultstate="collapsed">
    /**
     * Profile value of one feature in one cluster.
     *
     * @param columnIndex column index of the feature in the raw matrix
     * @param featureName name of the feature, e.g. the unique SMILES of a fragment or the descriptor name
     * @param category category of the feature
     * @param clusterValue share of molecules containing the feature (presence features) or mean (descriptors) in the
     *                     cluster; NaN if no value is available
     * @param overallValue the same for the whole data set
     * @param score difference of the shares (presence features) or of the means in standard deviations (descriptors);
     *              positive if the feature is enriched or higher in the cluster
     */
    public record FeatureProfile(int columnIndex, String featureName, FeatureCategory category, double clusterValue,
                                 double overallValue, double score) {
    }
    //
    /**
     * Profile of one cluster.
     *
     * @param clusterIndex index of the cluster
     * @param size number of molecules in the cluster
     * @param features profile values of all features in the order of the raw matrix columns
     */
    public record ClusterProfile(int clusterIndex, int size, List<FeatureProfile> features) {
        /**
         * Returns the features of the given category with the highest (positive) or lowest (negative) scores, i.e. the
         * most enriched or most depleted ones. Features with a score of 0 are not returned.
         *
         * @param aCategory feature category
         * @param aNumber maximum number of features
         * @param anIsPositive true for the highest scores (enriched), false for the lowest (depleted)
         * @return features sorted by score, the most extreme first
         */
        public List<FeatureProfile> getTopFeatures(FeatureCategory aCategory, int aNumber, boolean anIsPositive) {
            return this.features.stream()
                    .filter(aFeature -> aFeature.category() == aCategory)
                    .filter(aFeature -> anIsPositive ? aFeature.score() > 0.0 : aFeature.score() < 0.0)
                    .sorted(Comparator.comparingDouble((FeatureProfile aFeature) -> anIsPositive ? -aFeature.score() : aFeature.score())
                            .thenComparingInt(FeatureProfile::columnIndex))
                    .limit(aNumber)
                    .toList();
        }
        //
        /**
         * Returns the features of the given category with the largest absolute scores, e.g. the descriptors that deviate
         * most from the whole data set in either direction.
         *
         * @param aCategory feature category
         * @param aNumber maximum number of features
         * @return features sorted by absolute score, the most extreme first
         */
        public List<FeatureProfile> getMostDistinctFeatures(FeatureCategory aCategory, int aNumber) {
            return this.features.stream()
                    .filter(aFeature -> aFeature.category() == aCategory && aFeature.score() != 0.0)
                    .sorted(Comparator.comparingDouble((FeatureProfile aFeature) -> -Math.abs(aFeature.score()))
                            .thenComparingInt(FeatureProfile::columnIndex))
                    .limit(aNumber)
                    .toList();
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Clustering result that is profiled.
     */
    private final ClusteringRunResult runResult;
    /**
     * Category of every raw matrix column.
     */
    private final FeatureCategory[] categories;
    /**
     * Value of every column in the whole data set: share of molecules containing the feature or mean.
     */
    private final double[] overallValues;
    /**
     * Standard deviation of every descriptor column in the whole data set (0 for presence features).
     */
    private final double[] overallStandardDeviations;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, computes the values of the whole data set.
     *
     * @param aRunResult clustering result to profile
     * @throws NullPointerException if the result is null
     */
    public ClusterProfiler(ClusteringRunResult aRunResult) throws NullPointerException {
        this.runResult = Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        FeatureMatrix tmpRaw = aRunResult.rawMatrix();
        int tmpColumnCount = tmpRaw.getColumnCount();
        this.categories = new FeatureCategory[tmpColumnCount];
        this.overallValues = new double[tmpColumnCount];
        this.overallStandardDeviations = new double[tmpColumnCount];
        int[] tmpAllRows = new int[tmpRaw.getRowCount()];
        for (int i = 0; i < tmpAllRows.length; i++) {
            tmpAllRows[i] = i;
        }
        for (int j = 0; j < tmpColumnCount; j++) {
            if (FeatureMatrixBuilder.isMolecularDescriptorColumn(tmpRaw, j)) {
                this.categories[j] = FeatureCategory.DESCRIPTOR;
            } else if (tmpRaw.getFeatureGroupName(j).equals(aRunResult.parameters().fragmentationName())) {
                this.categories[j] = FeatureCategory.FRAGMENT;
            } else {
                this.categories[j] = FeatureCategory.FINGERPRINT_BIT;
            }
            this.overallValues[j] = this.getValue(j, tmpAllRows);
            if (this.categories[j] == FeatureCategory.DESCRIPTOR) {
                this.overallStandardDeviations[j] = this.getStandardDeviation(j, this.overallValues[j]);
            }
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the profile of the given cluster.
     *
     * @param aClusterIndex index of the cluster
     * @return cluster profile
     * @throws IllegalArgumentException if the cluster index is out of range
     */
    public ClusterProfile getProfile(int aClusterIndex) throws IllegalArgumentException {
        if (aClusterIndex < 0 || aClusterIndex >= this.runResult.clusteringResult().numberOfClusters()) {
            throw new IllegalArgumentException("Cluster index out of range: " + aClusterIndex);
        }
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        int[] tmpMembers = IntStream.range(0, tmpLabels.length).filter(i -> tmpLabels[i] == aClusterIndex).toArray();
        FeatureMatrix tmpRaw = this.runResult.rawMatrix();
        List<FeatureProfile> tmpFeatures = new ArrayList<>(tmpRaw.getColumnCount());
        for (int j = 0; j < tmpRaw.getColumnCount(); j++) {
            double tmpClusterValue = this.getValue(j, tmpMembers);
            double tmpScore;
            if (Double.isNaN(tmpClusterValue) || Double.isNaN(this.overallValues[j])) {
                tmpScore = 0.0;
            } else if (this.categories[j] == FeatureCategory.DESCRIPTOR) {
                tmpScore = this.overallStandardDeviations[j] > 0.0
                        ? (tmpClusterValue - this.overallValues[j]) / this.overallStandardDeviations[j] : 0.0;
            } else {
                tmpScore = tmpClusterValue - this.overallValues[j];
            }
            tmpFeatures.add(new FeatureProfile(j, tmpRaw.getFeatureName(j), this.categories[j], tmpClusterValue,
                    this.overallValues[j], tmpScore));
        }
        return new ClusterProfile(aClusterIndex, tmpMembers.length, tmpFeatures);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns the value of the given column over the given rows: the share of rows with a value &gt; 0 for presence
     * features, the mean of the finite values for descriptors.
     *
     * @param aColumnIndex raw matrix column
     * @param aRowIndices rows
     * @return share or mean; NaN if there is no (finite) value
     */
    private double getValue(int aColumnIndex, int[] aRowIndices) {
        FeatureMatrix tmpRaw = this.runResult.rawMatrix();
        double tmpSum = 0.0;
        int tmpCount = 0;
        for (int tmpRow : aRowIndices) {
            float tmpValue = tmpRaw.getValue(tmpRow, aColumnIndex);
            if (!Float.isFinite(tmpValue)) {
                continue;
            }
            tmpSum += this.categories[aColumnIndex] == FeatureCategory.DESCRIPTOR ? tmpValue : (tmpValue > 0.0f ? 1.0 : 0.0);
            tmpCount++;
        }
        return tmpCount == 0 ? Double.NaN : tmpSum / tmpCount;
    }
    //
    /**
     * Returns the population standard deviation of the finite values of the given column in the whole data set.
     *
     * @param aColumnIndex raw matrix column
     * @param aMean mean of the column
     * @return standard deviation; 0 if there is no finite value
     */
    private double getStandardDeviation(int aColumnIndex, double aMean) {
        if (Double.isNaN(aMean)) {
            return 0.0;
        }
        FeatureMatrix tmpRaw = this.runResult.rawMatrix();
        double tmpSumOfSquares = 0.0;
        int tmpCount = 0;
        for (int i = 0; i < tmpRaw.getRowCount(); i++) {
            float tmpValue = tmpRaw.getValue(i, aColumnIndex);
            if (Float.isFinite(tmpValue)) {
                tmpSumOfSquares += (tmpValue - aMean) * (tmpValue - aMean);
                tmpCount++;
            }
        }
        return tmpCount == 0 ? 0.0 : Math.sqrt(tmpSumOfSquares / tmpCount);
    }
    //</editor-fold>
}
