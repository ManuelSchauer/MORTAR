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
 * Data matrix of molecules (rows) and features (columns) together with the metadata of every column: its name, its
 * value type and the name of the group it belongs to (the descriptor enum constant name for descriptor components,
 * the fragmentation name for fragment fingerprint bits/counts). The group name is used to look up the priority of a
 * feature in the {@link DescriptorPriorityTable}.
 * <p>
 * The matrix may contain {@link Float#NaN} values for failed descriptor calculations. Instances are not changed by
 * any method of this class; methods that select or transform data return new instances.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class FeatureMatrix {
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Data matrix, data[i][j] is the value of feature j for molecule i.
     */
    private final float[][] data;
    /**
     * Identifiers of the rows, i.e. the molecule names.
     */
    private final String[] rowNames;
    /**
     * Names of the features (columns), e.g. "ALOGP [1]" or the unique SMILES of a fragment.
     */
    private final String[] featureNames;
    /**
     * Value types of the features (columns).
     */
    private final FeatureType[] featureTypes;
    /**
     * Names of the groups the features (columns) belong to, e.g. the descriptor enum constant name.
     */
    private final String[] featureGroupNames;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor. The given arrays are not copied.
     *
     * @param aData data matrix with one row per molecule, all rows must have the same length
     * @param aRowNameArray names of the rows (molecules), length must equal the number of rows
     * @param aFeatureNameArray names of the features, length must equal the number of columns
     * @param aFeatureTypeArray types of the features, length must equal the number of columns
     * @param aFeatureGroupNameArray group names of the features, length must equal the number of columns
     * @throws NullPointerException if any parameter is null
     * @throws IllegalArgumentException if the dimensions of the given arrays do not match
     */
    public FeatureMatrix(
            float[][] aData,
            String[] aRowNameArray,
            String[] aFeatureNameArray,
            FeatureType[] aFeatureTypeArray,
            String[] aFeatureGroupNameArray)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aData, "Given data matrix is null.");
        Objects.requireNonNull(aRowNameArray, "Given row name array is null.");
        Objects.requireNonNull(aFeatureNameArray, "Given feature name array is null.");
        Objects.requireNonNull(aFeatureTypeArray, "Given feature type array is null.");
        Objects.requireNonNull(aFeatureGroupNameArray, "Given feature group name array is null.");
        if (aData.length != aRowNameArray.length) {
            throw new IllegalArgumentException("Number of rows and number of row names differ.");
        }
        int tmpColumnCount = aFeatureNameArray.length;
        if (aFeatureTypeArray.length != tmpColumnCount || aFeatureGroupNameArray.length != tmpColumnCount) {
            throw new IllegalArgumentException("Lengths of feature name, type and group name arrays differ.");
        }
        for (float[] tmpRow : aData) {
            if (tmpRow == null || tmpRow.length != tmpColumnCount) {
                throw new IllegalArgumentException("A row of the data matrix is null or has a wrong length.");
            }
        }
        this.data = aData;
        this.rowNames = aRowNameArray;
        this.featureNames = aFeatureNameArray;
        this.featureTypes = aFeatureTypeArray;
        this.featureGroupNames = aFeatureGroupNameArray;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the number of rows (molecules).
     *
     * @return number of rows
     */
    public int getRowCount() {
        return this.data.length;
    }
    //
    /**
     * Returns the number of columns (features).
     *
     * @return number of columns
     */
    public int getColumnCount() {
        return this.featureNames.length;
    }
    //
    /**
     * Returns the value of the given feature for the given molecule.
     *
     * @param aRowIndex index of the molecule
     * @param aColumnIndex index of the feature
     * @return value, may be NaN
     */
    public float getValue(int aRowIndex, int aColumnIndex) {
        return this.data[aRowIndex][aColumnIndex];
    }
    //
    /**
     * Returns a copy of the given column.
     *
     * @param aColumnIndex index of the feature
     * @return copy of the column values
     */
    public float[] getColumn(int aColumnIndex) {
        float[] tmpColumn = new float[this.data.length];
        for (int i = 0; i < this.data.length; i++) {
            tmpColumn[i] = this.data[i][aColumnIndex];
        }
        return tmpColumn;
    }
    //
    /**
     * Returns a deep copy of the data matrix that can be passed to the clustering kernels.
     *
     * @return copy of the data matrix
     */
    public float[][] getDataCopy() {
        float[][] tmpCopy = new float[this.data.length][];
        for (int i = 0; i < this.data.length; i++) {
            tmpCopy[i] = Arrays.copyOf(this.data[i], this.data[i].length);
        }
        return tmpCopy;
    }
    //
    /**
     * Returns the name of the given row (molecule).
     *
     * @param aRowIndex index of the molecule
     * @return row name
     */
    public String getRowName(int aRowIndex) {
        return this.rowNames[aRowIndex];
    }
    //
    /**
     * Returns the name of the given feature.
     *
     * @param aColumnIndex index of the feature
     * @return feature name
     */
    public String getFeatureName(int aColumnIndex) {
        return this.featureNames[aColumnIndex];
    }
    //
    /**
     * Returns the type of the given feature.
     *
     * @param aColumnIndex index of the feature
     * @return feature type
     */
    public FeatureType getFeatureType(int aColumnIndex) {
        return this.featureTypes[aColumnIndex];
    }
    //
    /**
     * Returns the group name of the given feature, e.g. the descriptor enum constant name.
     *
     * @param aColumnIndex index of the feature
     * @return feature group name
     */
    public String getFeatureGroupName(int aColumnIndex) {
        return this.featureGroupNames[aColumnIndex];
    }
    //
    /**
     * Returns a new feature matrix that only contains the given columns in the given order.
     *
     * @param aColumnIndexArray indices of the columns to keep
     * @return new feature matrix with copied data
     * @throws NullPointerException if the given array is null
     */
    public FeatureMatrix selectColumns(int[] aColumnIndexArray) throws NullPointerException {
        Objects.requireNonNull(aColumnIndexArray, "Given column index array is null.");
        int tmpNewColumnCount = aColumnIndexArray.length;
        float[][] tmpData = new float[this.data.length][tmpNewColumnCount];
        String[] tmpNames = new String[tmpNewColumnCount];
        FeatureType[] tmpTypes = new FeatureType[tmpNewColumnCount];
        String[] tmpGroupNames = new String[tmpNewColumnCount];
        for (int j = 0; j < tmpNewColumnCount; j++) {
            int tmpOldIndex = aColumnIndexArray[j];
            tmpNames[j] = this.featureNames[tmpOldIndex];
            tmpTypes[j] = this.featureTypes[tmpOldIndex];
            tmpGroupNames[j] = this.featureGroupNames[tmpOldIndex];
            for (int i = 0; i < this.data.length; i++) {
                tmpData[i][j] = this.data[i][tmpOldIndex];
            }
        }
        return new FeatureMatrix(tmpData, this.rowNames.clone(), tmpNames, tmpTypes, tmpGroupNames);
    }
    //
    /**
     * Returns a new feature matrix with the same metadata but the given data, e.g. after a transformation.
     *
     * @param aData new data matrix with the same dimensions as this matrix
     * @return new feature matrix
     * @throws NullPointerException if the given data is null
     * @throws IllegalArgumentException if the dimensions do not match
     */
    public FeatureMatrix withData(float[][] aData) throws NullPointerException, IllegalArgumentException {
        return new FeatureMatrix(aData, this.rowNames.clone(), this.featureNames.clone(),
                this.featureTypes.clone(), this.featureGroupNames.clone());
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Concatenates the given feature matrices horizontally, i.e. the columns of all matrices are put side by side.
     * All matrices must have the same rows in the same order.
     *
     * @param aFeatureMatrices feature matrices to concatenate, at least one
     * @return new feature matrix containing the columns of all given matrices
     * @throws NullPointerException if the given array or one of its elements is null
     * @throws IllegalArgumentException if no matrix is given or the row names differ
     */
    public static FeatureMatrix concatenate(FeatureMatrix... aFeatureMatrices) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aFeatureMatrices, "Given feature matrix array is null.");
        if (aFeatureMatrices.length == 0) {
            throw new IllegalArgumentException("No feature matrix given.");
        }
        FeatureMatrix tmpFirst = Objects.requireNonNull(aFeatureMatrices[0], "A given feature matrix is null.");
        int tmpTotalColumnCount = 0;
        for (FeatureMatrix tmpMatrix : aFeatureMatrices) {
            Objects.requireNonNull(tmpMatrix, "A given feature matrix is null.");
            if (!Arrays.equals(tmpMatrix.rowNames, tmpFirst.rowNames)) {
                throw new IllegalArgumentException("Feature matrices have different rows.");
            }
            tmpTotalColumnCount += tmpMatrix.getColumnCount();
        }
        float[][] tmpData = new float[tmpFirst.getRowCount()][tmpTotalColumnCount];
        String[] tmpNames = new String[tmpTotalColumnCount];
        FeatureType[] tmpTypes = new FeatureType[tmpTotalColumnCount];
        String[] tmpGroupNames = new String[tmpTotalColumnCount];
        int tmpOffset = 0;
        for (FeatureMatrix tmpMatrix : aFeatureMatrices) {
            int tmpColumnCount = tmpMatrix.getColumnCount();
            System.arraycopy(tmpMatrix.featureNames, 0, tmpNames, tmpOffset, tmpColumnCount);
            System.arraycopy(tmpMatrix.featureTypes, 0, tmpTypes, tmpOffset, tmpColumnCount);
            System.arraycopy(tmpMatrix.featureGroupNames, 0, tmpGroupNames, tmpOffset, tmpColumnCount);
            for (int i = 0; i < tmpData.length; i++) {
                System.arraycopy(tmpMatrix.data[i], 0, tmpData[i], tmpOffset, tmpColumnCount);
            }
            tmpOffset += tmpColumnCount;
        }
        return new FeatureMatrix(tmpData, tmpFirst.rowNames.clone(), tmpNames, tmpTypes, tmpGroupNames);
    }
    //</editor-fold>
}
