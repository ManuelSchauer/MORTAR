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
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Statistics of the descriptor values of the clusters of a clustering result: per descriptor column the median in every
 * cluster and the median of the whole data set. The values are taken from the raw feature matrix,
 * i.e. untransformed and in the original units of the descriptors; NaN values (descriptor could not be calculated) are
 * ignored. Descriptor columns are the CONTINUOUS columns of the raw matrix, including the ones removed by the filter.
 * The statistics of a column are calculated on first use and cached.
 * <p>
 * The range filter marks the clusters whose medians lie within given bounds for all given descriptors.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusterDescriptorStatistics {
    //<editor-fold desc="Records" defaultstate="collapsed">
    /**
     * Condition of the range filter: the median of the cluster must lie within the bounds (both included).
     *
     * @param columnName name of the descriptor column in the raw matrix
     * @param lowerBound lower bound, NaN for no lower bound
     * @param upperBound upper bound, NaN for no upper bound
     */
    public record RangeCondition(String columnName, double lowerBound, double upperBound) {
        /**
         * Compact constructor that checks the name for null.
         *
         * @throws NullPointerException if the column name is null
         */
        public RangeCondition {
            Objects.requireNonNull(columnName, "Given column name is null.");
        }
        //
        /**
         * Returns whether the given value lies within the bounds; NaN never does.
         *
         * @param aValue value, e.g. the median of a cluster
         * @return true if the value satisfies the condition
         */
        public boolean isSatisfiedBy(double aValue) {
            return !Double.isNaN(aValue)
                    && (Double.isNaN(this.lowerBound) || aValue >= this.lowerBound)
                    && (Double.isNaN(this.upperBound) || aValue <= this.upperBound);
        }
    }
    //
    /**
     * Statistics of one descriptor column.
     *
     * @param clusterMedians median per cluster (index = cluster index), NaN if the cluster has no finite value
     * @param overallMedian median of all molecules, NaN if there is no finite value
     */
    private record ColumnStatistics(double[] clusterMedians, double overallMedian) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Clustering result.
     */
    private final ClusteringRunResult runResult;
    /**
     * Raw matrix column index of every descriptor column name, in column order.
     */
    private final Map<String, Integer> columnIndexOfNameMap;
    /**
     * Cache of the statistics of the columns that were already used, by raw matrix column index.
     */
    private final Map<Integer, ColumnStatistics> statisticsOfColumnMap;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aRunResult clustering result
     * @throws NullPointerException if the result is null
     */
    public ClusterDescriptorStatistics(ClusteringRunResult aRunResult) throws NullPointerException {
        this.runResult = Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        FeatureMatrix tmpRaw = aRunResult.rawMatrix();
        this.columnIndexOfNameMap = new HashMap<>();
        for (int j = 0; j < tmpRaw.getColumnCount(); j++) {
            if (tmpRaw.getFeatureType(j) == FeatureType.CONTINUOUS) {
                this.columnIndexOfNameMap.putIfAbsent(tmpRaw.getFeatureName(j), j);
            }
        }
        this.statisticsOfColumnMap = new HashMap<>();
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the names of all descriptor columns in the order of the raw matrix.
     *
     * @return unmodifiable list of descriptor column names
     */
    public List<String> getDescriptorColumnNames() {
        List<String> tmpNames = new ArrayList<>(this.columnIndexOfNameMap.keySet());
        tmpNames.sort((aName1, aName2) -> Integer.compare(this.columnIndexOfNameMap.get(aName1), this.columnIndexOfNameMap.get(aName2)));
        return Collections.unmodifiableList(tmpNames);
    }
    //
    /**
     * Returns whether the clustering result contains a descriptor column with the given name.
     *
     * @param aColumnName descriptor column name
     * @return true if the column exists
     */
    public boolean hasDescriptorColumn(String aColumnName) {
        return this.columnIndexOfNameMap.containsKey(aColumnName);
    }
    //
    /**
     * Returns whether the descriptor column with the given name was used for the clustering, i.e. not removed by the
     * filter.
     *
     * @param aColumnName descriptor column name
     * @return true if the column is part of the processed matrix
     */
    public boolean isUsedForClustering(String aColumnName) {
        FeatureMatrix tmpProcessed = this.runResult.processedMatrix();
        for (int j = 0; j < tmpProcessed.getColumnCount(); j++) {
            if (tmpProcessed.getFeatureType(j) == FeatureType.CONTINUOUS && tmpProcessed.getFeatureName(j).equals(aColumnName)) {
                return true;
            }
        }
        return false;
    }
    //
    /**
     * Returns the median of the given descriptor in the given cluster.
     *
     * @param aClusterIndex cluster index; -1 (molecules that could not be clustered) returns NaN
     * @param aColumnName descriptor column name
     * @return median, NaN if the cluster has no finite value
     * @throws IllegalArgumentException if the column does not exist or the cluster index is out of range
     */
    public double getMedian(int aClusterIndex, String aColumnName) throws IllegalArgumentException {
        return aClusterIndex < 0 ? Double.NaN : this.getStatistics(aColumnName).clusterMedians()[this.checkClusterIndex(aClusterIndex)];
    }
    //
    /**
     * Returns the median of the given descriptor over all molecules of the data set.
     *
     * @param aColumnName descriptor column name
     * @return median, NaN if there is no finite value
     * @throws IllegalArgumentException if the column does not exist
     */
    public double getOverallMedian(String aColumnName) throws IllegalArgumentException {
        return this.getStatistics(aColumnName).overallMedian();
    }
    //
    /**
     * Returns the value of the given descriptor for the molecule in the given row, untransformed.
     *
     * @param aRowIndex row index of the molecule in the clustered data
     * @param aColumnName descriptor column name
     * @return raw value, NaN if it could not be calculated
     * @throws IllegalArgumentException if the column does not exist
     */
    public double getValue(int aRowIndex, String aColumnName) throws IllegalArgumentException {
        return this.runResult.rawMatrix().getValue(aRowIndex, this.getColumnIndex(aColumnName));
    }
    //
    /**
     * Range filter: returns the clusters whose median satisfies all given conditions. Clusters without a finite value for
     * a descriptor of the conditions do not match. Without conditions, no cluster is returned.
     *
     * @param aConditionList conditions, all must be satisfied
     * @return indices of the matching clusters in ascending order
     * @throws NullPointerException if the list is null
     * @throws IllegalArgumentException if a condition refers to a column that does not exist
     */
    public Set<Integer> getMatchingClusters(List<RangeCondition> aConditionList) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aConditionList, "Given condition list is null.");
        Set<Integer> tmpMatchingClusters = new LinkedHashSet<>();
        if (aConditionList.isEmpty()) {
            return tmpMatchingClusters;
        }
        for (int k = 0; k < this.runResult.clusteringResult().numberOfClusters(); k++) {
            boolean tmpIsMatching = true;
            for (RangeCondition tmpCondition : aConditionList) {
                if (!tmpCondition.isSatisfiedBy(this.getMedian(k, tmpCondition.columnName()))) {
                    tmpIsMatching = false;
                    break;
                }
            }
            if (tmpIsMatching) {
                tmpMatchingClusters.add(k);
            }
        }
        return tmpMatchingClusters;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns the median of the given values, which must be sorted ascending; the mean of the two middle values for an
     * even number of values.
     *
     * @param aSortedValues sorted values
     * @param aCount number of values to use from the start of the array
     * @return median, NaN if the count is 0
     */
    public static double getMedianOfSorted(double[] aSortedValues, int aCount) {
        if (aCount == 0) {
            return Double.NaN;
        }
        int tmpMiddle = aCount / 2;
        return aCount % 2 == 1 ? aSortedValues[tmpMiddle] : (aSortedValues[tmpMiddle - 1] + aSortedValues[tmpMiddle]) / 2.0;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns the raw matrix column index of the given descriptor column.
     *
     * @param aColumnName descriptor column name
     * @return column index
     * @throws IllegalArgumentException if the column does not exist
     */
    private int getColumnIndex(String aColumnName) throws IllegalArgumentException {
        Integer tmpIndex = this.columnIndexOfNameMap.get(aColumnName);
        if (tmpIndex == null) {
            throw new IllegalArgumentException("No descriptor column with name " + aColumnName);
        }
        return tmpIndex;
    }
    //
    /**
     * Checks the given cluster index.
     *
     * @param aClusterIndex cluster index
     * @return the cluster index
     * @throws IllegalArgumentException if the index is out of range
     */
    private int checkClusterIndex(int aClusterIndex) throws IllegalArgumentException {
        if (aClusterIndex >= this.runResult.clusteringResult().numberOfClusters()) {
            throw new IllegalArgumentException("Cluster index out of range: " + aClusterIndex);
        }
        return aClusterIndex;
    }
    //
    /**
     * Returns the statistics of the given descriptor column, calculating them on first use.
     *
     * @param aColumnName descriptor column name
     * @return column statistics
     * @throws IllegalArgumentException if the column does not exist
     */
    private ColumnStatistics getStatistics(String aColumnName) throws IllegalArgumentException {
        int tmpColumnIndex = this.getColumnIndex(aColumnName);
        return this.statisticsOfColumnMap.computeIfAbsent(tmpColumnIndex, this::calculateStatistics);
    }
    //
    /**
     * Calculates the statistics of the given raw matrix column.
     *
     * @param aColumnIndex raw matrix column index
     * @return column statistics
     */
    private ColumnStatistics calculateStatistics(int aColumnIndex) {
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        int tmpNumberOfClusters = this.runResult.clusteringResult().numberOfClusters();
        float[] tmpColumn = this.runResult.rawMatrix().getColumn(aColumnIndex);
        //finite values grouped by cluster
        double[][] tmpClusterValues = new double[tmpNumberOfClusters][];
        int[] tmpClusterCounts = new int[tmpNumberOfClusters];
        for (int tmpLabel : tmpLabels) {
            if (tmpLabel >= 0) {
                tmpClusterCounts[tmpLabel]++;
            }
        }
        for (int k = 0; k < tmpNumberOfClusters; k++) {
            tmpClusterValues[k] = new double[tmpClusterCounts[k]];
        }
        Arrays.fill(tmpClusterCounts, 0);
        double[] tmpAllValues = new double[tmpColumn.length];
        int tmpAllCount = 0;
        for (int i = 0; i < tmpColumn.length; i++) {
            if (!Float.isFinite(tmpColumn[i])) {
                continue;
            }
            tmpAllValues[tmpAllCount++] = tmpColumn[i];
            if (tmpLabels[i] >= 0) {
                tmpClusterValues[tmpLabels[i]][tmpClusterCounts[tmpLabels[i]]++] = tmpColumn[i];
            }
        }
        double[] tmpMedians = new double[tmpNumberOfClusters];
        for (int k = 0; k < tmpNumberOfClusters; k++) {
            Arrays.sort(tmpClusterValues[k], 0, tmpClusterCounts[k]);
            tmpMedians[k] = ClusterDescriptorStatistics.getMedianOfSorted(tmpClusterValues[k], tmpClusterCounts[k]);
        }
        Arrays.sort(tmpAllValues, 0, tmpAllCount);
        return new ColumnStatistics(tmpMedians, ClusterDescriptorStatistics.getMedianOfSorted(tmpAllValues, tmpAllCount));
    }
    //</editor-fold>
}
