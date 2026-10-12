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

/**
 * Applies {@link FeatureTransformation}s to the columns of a {@link FeatureMatrix} depending on their
 * {@link FeatureType}. BIT features are never transformed. NaN values stay NaN.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class FeatureTransformer {
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Lower percentile used by {@link FeatureTransformation#WINSORIZE}.
     */
    public static final double WINSORIZE_LOWER_PERCENTILE = 0.01;
    /**
     * Upper percentile used by {@link FeatureTransformation#WINSORIZE}.
     */
    public static final double WINSORIZE_UPPER_PERCENTILE = 0.99;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private FeatureTransformer() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns a new feature matrix in which the COUNT columns are transformed with the given count transformation and
     * the CONTINUOUS columns with the given continuous transformation. BIT columns are copied unchanged.
     *
     * @param aFeatureMatrix matrix to transform (IS NOT CHANGED)
     * @param aCountTransformation transformation for COUNT columns
     * @param aContinuousTransformation transformation for CONTINUOUS columns
     * @return new transformed feature matrix
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if LOG1P or SQRT is requested for a column that contains negative values
     */
    public static FeatureMatrix transform(
            FeatureMatrix aFeatureMatrix,
            FeatureTransformation aCountTransformation,
            FeatureTransformation aContinuousTransformation)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(aCountTransformation, "Given count transformation is null.");
        Objects.requireNonNull(aContinuousTransformation, "Given continuous transformation is null.");
        float[][] tmpData = aFeatureMatrix.getDataCopy();
        for (int j = 0; j < aFeatureMatrix.getColumnCount(); j++) {
            FeatureTransformation tmpTransformation = switch (aFeatureMatrix.getFeatureType(j)) {
                case BIT -> FeatureTransformation.NONE;
                case COUNT -> aCountTransformation;
                case CONTINUOUS -> aContinuousTransformation;
            };
            float[] tmpColumn = aFeatureMatrix.getColumn(j);
            FeatureTransformer.transformColumn(tmpColumn, tmpTransformation, aFeatureMatrix.getFeatureName(j));
            for (int i = 0; i < tmpData.length; i++) {
                tmpData[i][j] = tmpColumn[i];
            }
        }
        return aFeatureMatrix.withData(tmpData);
    }
    //
    /**
     * Transforms the given column in place.
     *
     * @param aColumn column values (IS CHANGED)
     * @param aTransformation transformation to apply
     * @param aFeatureName name of the feature, used in exception messages
     * @throws IllegalArgumentException if LOG1P or SQRT is requested and the column contains negative values
     */
    public static void transformColumn(float[] aColumn, FeatureTransformation aTransformation, String aFeatureName)
            throws IllegalArgumentException {
        switch (aTransformation) {
            case NONE:
                break;
            case LOG1P:
            case SQRT:
                for (int i = 0; i < aColumn.length; i++) {
                    if (aColumn[i] < 0.0f) {
                        throw new IllegalArgumentException(aTransformation + " is not applicable to negative values of feature " + aFeatureName);
                    }
                    //NaN stays NaN in both functions
                    aColumn[i] = aTransformation == FeatureTransformation.LOG1P
                            ? (float) Math.log1p(aColumn[i])
                            : (float) Math.sqrt(aColumn[i]);
                }
                break;
            case BINARIZE:
                for (int i = 0; i < aColumn.length; i++) {
                    if (!Float.isNaN(aColumn[i])) {
                        aColumn[i] = aColumn[i] > 0.0f ? 1.0f : 0.0f;
                    }
                }
                break;
            case WINSORIZE:
                float[] tmpFiniteValues = FeatureTransformer.getSortedFiniteValues(aColumn);
                if (tmpFiniteValues.length == 0) {
                    break;
                }
                float tmpLowerBound = FeatureTransformer.getPercentile(tmpFiniteValues, FeatureTransformer.WINSORIZE_LOWER_PERCENTILE);
                float tmpUpperBound = FeatureTransformer.getPercentile(tmpFiniteValues, FeatureTransformer.WINSORIZE_UPPER_PERCENTILE);
                for (int i = 0; i < aColumn.length; i++) {
                    if (!Float.isNaN(aColumn[i])) {
                        aColumn[i] = Math.min(Math.max(aColumn[i], tmpLowerBound), tmpUpperBound);
                    }
                }
                break;
            default:
                throw new IllegalArgumentException("Undefined feature transformation.");
        }
    }
    //
    /**
     * Returns the percentile of the given sorted values using linear interpolation between the closest ranks.
     *
     * @param aSortedValues sorted values, at least one
     * @param aPercentile percentile in [0,1]
     * @return percentile value
     */
    public static float getPercentile(float[] aSortedValues, double aPercentile) {
        double tmpPosition = aPercentile * (aSortedValues.length - 1);
        int tmpLowerIndex = (int) Math.floor(tmpPosition);
        int tmpUpperIndex = (int) Math.ceil(tmpPosition);
        double tmpFraction = tmpPosition - tmpLowerIndex;
        return (float) (aSortedValues[tmpLowerIndex] + tmpFraction * (aSortedValues[tmpUpperIndex] - aSortedValues[tmpLowerIndex]));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns the finite values of the given column in ascending order.
     *
     * @param aColumn column values
     * @return sorted finite values
     */
    private static float[] getSortedFiniteValues(float[] aColumn) {
        float[] tmpValues = new float[aColumn.length];
        int tmpCount = 0;
        for (float tmpValue : aColumn) {
            if (Float.isFinite(tmpValue)) {
                tmpValues[tmpCount++] = tmpValue;
            }
        }
        float[] tmpFiniteValues = Arrays.copyOf(tmpValues, tmpCount);
        Arrays.sort(tmpFiniteValues);
        return tmpFiniteValues;
    }
    //</editor-fold>
}
