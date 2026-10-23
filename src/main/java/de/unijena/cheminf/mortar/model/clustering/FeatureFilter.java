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
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

/**
 * Filters the columns of a {@link FeatureMatrix} before clustering. The steps are performed in this order:
 * <ol>
 *     <li>Columns whose fraction of NaN values exceeds the given maximum are removed; remaining NaN values are replaced
 *     by the median of their column (the clustering kernels do not accept NaN).</li>
 *     <li>Constant columns are removed (they would still count in the contrast enhancement threshold of the kernels).</li>
 *     <li>Correlation filter with priority lookup: COUNT and CONTINUOUS columns are ordered by descending priority
 *     (see {@link DescriptorPriorityTable}) and, for equal priority, by their original column order. A column is kept
 *     if its absolute correlation with every already kept column is not above the threshold; otherwise it is removed
 *     and the kept partner is recorded. BIT columns do not take part in the correlation filter.</li>
 * </ol>
 * The filter should be applied after the {@link FeatureTransformer}, because Pearson's r changes under non-linear
 * transformations.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class FeatureFilter {
    //<editor-fold desc="Enum CorrelationMethod" defaultstate="collapsed">
    /**
     * Enumeration of the available correlation coefficients.
     */
    public enum CorrelationMethod {
        /**
         * Pearson product-moment correlation, measures linear dependence.
         */
        PEARSON,
        /**
         * Spearman rank correlation, measures monotone dependence and is robust against skewed values and outliers.
         */
        SPEARMAN;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record FilterResult" defaultstate="collapsed">
    /**
     * Result of the filter.
     *
     * @param filteredMatrix feature matrix that only contains the kept columns, NaN values imputed
     * @param removedFeatures report of the removed features in the order of removal
     */
    public record FilterResult(FeatureMatrix filteredMatrix, List<RemovedFeature> removedFeatures) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Columns whose range (maximum minus minimum) is not larger than this value are regarded as constant.
     */
    private static final double CONSTANT_RANGE_TOLERANCE = 1.0e-9;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private FeatureFilter() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Filters the given feature matrix, see the class description for the steps.
     *
     * @param aFeatureMatrix matrix to filter (IS NOT CHANGED)
     * @param aPriorityTable priority lookup table for the correlation filter
     * @param aMaximumNaNFraction maximum allowed fraction of NaN values per column in [0,1]
     * @param aCorrelationThreshold columns with an absolute correlation above this value are regarded as redundant; a
     *                              value of 1 or larger disables the correlation filter
     * @param aCorrelationMethod correlation coefficient to use
     * @return filter result with the filtered matrix and the report of removed features
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the maximum NaN fraction is not in [0,1] or the threshold is negative
     */
    public static FilterResult filter(
            FeatureMatrix aFeatureMatrix,
            DescriptorPriorityTable aPriorityTable,
            double aMaximumNaNFraction,
            double aCorrelationThreshold,
            CorrelationMethod aCorrelationMethod)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(aPriorityTable, "Given priority table is null.");
        Objects.requireNonNull(aCorrelationMethod, "Given correlation method is null.");
        if (aMaximumNaNFraction < 0.0 || aMaximumNaNFraction > 1.0) {
            throw new IllegalArgumentException("Maximum NaN fraction must be in [0,1].");
        }
        if (aCorrelationThreshold < 0.0) {
            throw new IllegalArgumentException("Correlation threshold must not be negative.");
        }
        int tmpRowCount = aFeatureMatrix.getRowCount();
        int tmpColumnCount = aFeatureMatrix.getColumnCount();
        List<RemovedFeature> tmpRemovedFeatures = new ArrayList<>();
        boolean[] tmpIsRemoved = new boolean[tmpColumnCount];
        float[][] tmpImputedColumns = new float[tmpColumnCount][];
        //step 1 and 2: missing values and constant columns
        for (int j = 0; j < tmpColumnCount; j++) {
            float[] tmpColumn = aFeatureMatrix.getColumn(j);
            int tmpNaNCount = 0;
            for (float tmpValue : tmpColumn) {
                if (!Float.isFinite(tmpValue)) {
                    tmpNaNCount++;
                }
            }
            if (tmpRowCount == 0 || tmpNaNCount > aMaximumNaNFraction * tmpRowCount || tmpNaNCount == tmpRowCount) {
                tmpIsRemoved[j] = true;
                tmpRemovedFeatures.add(new RemovedFeature(aFeatureMatrix.getFeatureName(j),
                        RemovedFeature.RemovalReason.MISSING_VALUES, null, Double.NaN));
                continue;
            }
            if (tmpNaNCount > 0) {
                FeatureFilter.imputeWithMedian(tmpColumn);
            }
            float tmpMinimum = Float.POSITIVE_INFINITY;
            float tmpMaximum = Float.NEGATIVE_INFINITY;
            for (float tmpValue : tmpColumn) {
                tmpMinimum = Math.min(tmpMinimum, tmpValue);
                tmpMaximum = Math.max(tmpMaximum, tmpValue);
            }
            if (tmpMaximum - tmpMinimum <= FeatureFilter.CONSTANT_RANGE_TOLERANCE) {
                tmpIsRemoved[j] = true;
                tmpRemovedFeatures.add(new RemovedFeature(aFeatureMatrix.getFeatureName(j),
                        RemovedFeature.RemovalReason.CONSTANT, null, Double.NaN));
                continue;
            }
            tmpImputedColumns[j] = tmpColumn;
        }
        //step 3: greedy correlation filter by priority
        if (aCorrelationThreshold < 1.0) {
            List<Integer> tmpCandidates = new ArrayList<>();
            for (int j = 0; j < tmpColumnCount; j++) {
                if (!tmpIsRemoved[j] && aFeatureMatrix.getFeatureType(j) != FeatureType.BIT) {
                    tmpCandidates.add(j);
                }
            }
            tmpCandidates.sort(Comparator.comparingInt((Integer anIndex) -> -aPriorityTable.getPriority(aFeatureMatrix.getFeatureGroupName(anIndex)))
                    .thenComparingInt(anIndex -> anIndex));
            List<Integer> tmpKeptIndices = new ArrayList<>();
            List<double[]> tmpKeptStandardizedColumns = new ArrayList<>();
            for (int tmpCandidate : tmpCandidates) {
                double[] tmpStandardized = FeatureFilter.getStandardizedColumn(tmpImputedColumns[tmpCandidate], aCorrelationMethod);
                int tmpPartner = -1;
                double tmpPartnerCorrelation = 0.0;
                for (int k = 0; k < tmpKeptIndices.size(); k++) {
                    double tmpCorrelation = FeatureFilter.getDotProduct(tmpStandardized, tmpKeptStandardizedColumns.get(k));
                    if (Math.abs(tmpCorrelation) > aCorrelationThreshold) {
                        tmpPartner = tmpKeptIndices.get(k);
                        tmpPartnerCorrelation = tmpCorrelation;
                        break;
                    }
                }
                if (tmpPartner >= 0) {
                    tmpIsRemoved[tmpCandidate] = true;
                    tmpRemovedFeatures.add(new RemovedFeature(aFeatureMatrix.getFeatureName(tmpCandidate),
                            RemovedFeature.RemovalReason.CORRELATION, aFeatureMatrix.getFeatureName(tmpPartner), tmpPartnerCorrelation));
                } else {
                    tmpKeptIndices.add(tmpCandidate);
                    tmpKeptStandardizedColumns.add(tmpStandardized);
                }
            }
        }
        //assemble the filtered matrix in the original column order
        int[] tmpKeptColumnIndices = new int[tmpColumnCount];
        int tmpKeptCount = 0;
        for (int j = 0; j < tmpColumnCount; j++) {
            if (!tmpIsRemoved[j]) {
                tmpKeptColumnIndices[tmpKeptCount++] = j;
            }
        }
        tmpKeptColumnIndices = Arrays.copyOf(tmpKeptColumnIndices, tmpKeptCount);
        float[][] tmpData = new float[tmpRowCount][tmpKeptCount];
        for (int k = 0; k < tmpKeptCount; k++) {
            float[] tmpColumn = tmpImputedColumns[tmpKeptColumnIndices[k]];
            for (int i = 0; i < tmpRowCount; i++) {
                tmpData[i][k] = tmpColumn[i];
            }
        }
        FeatureMatrix tmpFilteredMatrix = aFeatureMatrix.selectColumns(tmpKeptColumnIndices).withData(tmpData);
        return new FilterResult(tmpFilteredMatrix, tmpRemovedFeatures);
    }
    //
    /**
     * Calculates the correlation coefficient of the two given columns. The columns must not contain NaN values.
     *
     * @param aFirstColumn first column
     * @param aSecondColumn second column of the same length
     * @param aCorrelationMethod correlation coefficient to use
     * @return correlation coefficient in [-1,1]; 0 if one of the columns is constant
     */
    public static double getCorrelation(float[] aFirstColumn, float[] aSecondColumn, CorrelationMethod aCorrelationMethod) {
        return FeatureFilter.getDotProduct(
                FeatureFilter.getStandardizedColumn(aFirstColumn, aCorrelationMethod),
                FeatureFilter.getStandardizedColumn(aSecondColumn, aCorrelationMethod));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns the centered column scaled to unit Euclidean length, so that the dot product of two standardized columns
     * is their Pearson correlation. For Spearman, the values are replaced by their ranks first.
     *
     * @param aColumn column values without NaN
     * @param aCorrelationMethod correlation coefficient to use
     * @return standardized column; all zeros if the column is constant
     */
    private static double[] getStandardizedColumn(float[] aColumn, CorrelationMethod aCorrelationMethod) {
        double[] tmpValues = aCorrelationMethod == CorrelationMethod.SPEARMAN
                ? FeatureFilter.getRanks(aColumn)
                : FeatureFilter.toDoubleArray(aColumn);
        double tmpMean = 0.0;
        for (double tmpValue : tmpValues) {
            tmpMean += tmpValue;
        }
        tmpMean /= tmpValues.length;
        double tmpSumOfSquares = 0.0;
        for (int i = 0; i < tmpValues.length; i++) {
            tmpValues[i] -= tmpMean;
            tmpSumOfSquares += tmpValues[i] * tmpValues[i];
        }
        double tmpLength = Math.sqrt(tmpSumOfSquares);
        if (tmpLength > 0.0) {
            for (int i = 0; i < tmpValues.length; i++) {
                tmpValues[i] /= tmpLength;
            }
        }
        return tmpValues;
    }
    //
    /**
     * Returns the ranks of the given values (1-based); tied values get their average rank.
     *
     * @param aColumn column values
     * @return ranks
     */
    private static double[] getRanks(float[] aColumn) {
        Integer[] tmpOrder = new Integer[aColumn.length];
        for (int i = 0; i < aColumn.length; i++) {
            tmpOrder[i] = i;
        }
        Arrays.sort(tmpOrder, Comparator.comparingDouble(anIndex -> aColumn[anIndex]));
        double[] tmpRanks = new double[aColumn.length];
        int tmpStart = 0;
        while (tmpStart < aColumn.length) {
            int tmpEnd = tmpStart;
            while (tmpEnd + 1 < aColumn.length && aColumn[tmpOrder[tmpEnd + 1]] == aColumn[tmpOrder[tmpStart]]) {
                tmpEnd++;
            }
            double tmpAverageRank = (tmpStart + tmpEnd) / 2.0 + 1.0;
            for (int k = tmpStart; k <= tmpEnd; k++) {
                tmpRanks[tmpOrder[k]] = tmpAverageRank;
            }
            tmpStart = tmpEnd + 1;
        }
        return tmpRanks;
    }
    //
    /**
     * Replaces all non-finite values of the given column by the median of its finite values.
     *
     * @param aColumn column values (IS CHANGED), at least one finite value
     */
    private static void imputeWithMedian(float[] aColumn) {
        float[] tmpFinite = new float[aColumn.length];
        int tmpCount = 0;
        for (float tmpValue : aColumn) {
            if (Float.isFinite(tmpValue)) {
                tmpFinite[tmpCount++] = tmpValue;
            }
        }
        tmpFinite = Arrays.copyOf(tmpFinite, tmpCount);
        Arrays.sort(tmpFinite);
        float tmpMedian = FeatureTransformer.getPercentile(tmpFinite, 0.5);
        for (int i = 0; i < aColumn.length; i++) {
            if (!Float.isFinite(aColumn[i])) {
                aColumn[i] = tmpMedian;
            }
        }
    }
    //
    /**
     * Returns the dot product of the two given vectors.
     *
     * @param aFirstVector first vector
     * @param aSecondVector second vector of the same length
     * @return dot product
     */
    private static double getDotProduct(double[] aFirstVector, double[] aSecondVector) {
        double tmpSum = 0.0;
        for (int i = 0; i < aFirstVector.length; i++) {
            tmpSum += aFirstVector[i] * aSecondVector[i];
        }
        return tmpSum;
    }
    //
    /**
     * Converts the given float array into a double array.
     *
     * @param anArray float array
     * @return double array
     */
    private static double[] toDoubleArray(float[] anArray) {
        double[] tmpArray = new double[anArray.length];
        for (int i = 0; i < anArray.length; i++) {
            tmpArray[i] = anArray[i];
        }
        return tmpArray;
    }
    //</editor-fold>
}
