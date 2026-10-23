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
 * Similarity search in the feature space of a clustering run: finds the molecules most similar to a query molecule and
 * the cluster the query belongs to.
 * <p>
 * The ART-2A library offers no nearest-neighbour search and does not assign new vectors publicly. This class therefore
 * reproduces the internal preprocessing of the kernel that was used (per-column min-max scaling with the minima and
 * maxima of the clustered data, for the cosine kernel also L2 normalisation, contrast enhancement with
 * θ = 1/sqrt(n + offset) and renormalisation), so that "similar" means exactly what the clustering saw:
 * <ul>
 *     <li>cosine kernel: similarity = cosine of the preprocessed vectors (the value compared with the vigilance)</li>
 *     <li>Euclid kernel: similarity = 1 - squared distance / n (a new cluster is opened by the kernel if this is below
 *     the vigilance)</li>
 * </ul>
 * The cluster of a query is found like the ART-2A winner search: the cluster vector (prototype) with the highest
 * similarity; it is reported together with whether the similarity passes the vigilance. Molecules of the clustered data
 * keep their own cluster.
 * <p>
 * A query that is not part of the clustered data is given as raw feature vector (see
 * {@link FeatureMatrixBuilder#buildRawFeatureVector}); it is transformed like the data (same transformations, winsorizing
 * bounds of the data, NaN replaced by the median of the data) before the preprocessing.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SimilaritySearch {
    //<editor-fold desc="Records" defaultstate="collapsed">
    /**
     * A molecule found by the search.
     *
     * @param rowIndex row index of the molecule in the clustered data
     * @param similarity similarity to the query in [0,1]
     */
    public record Neighbor(int rowIndex, double similarity) {
    }
    //
    /**
     * Result of a search.
     *
     * @param neighbors most similar molecules, most similar first
     * @param clusterIndex cluster of the query; -1 if the query has no features in the kernel space (zero vector)
     * @param clusterSimilarity similarity of the query to the cluster vector of its cluster
     * @param isWithinVigilance true if the cluster similarity passes the vigilance, i.e. the kernel would have put the
     *                          query into this cluster; false if it would have opened a new cluster
     */
    public record SearchResult(List<Neighbor> neighbors, int clusterIndex, double clusterSimilarity, boolean isWithinVigilance) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Clustering run whose feature space is searched.
     */
    private final ClusteringRunResult runResult;
    /**
     * True for the cosine kernel, false for the Euclid kernel.
     */
    private final boolean isCosineKernel;
    /**
     * Threshold of the contrast enhancement, 1/sqrt(n + offset).
     */
    private final float contrastEnhancementThreshold;
    /**
     * Minimum of every processed column, used for the min-max scaling.
     */
    private final float[] columnMinima;
    /**
     * Maximum of every processed column, used for the min-max scaling.
     */
    private final float[] columnMaxima;
    /**
     * Index of the raw (and transformed) column of every processed column.
     */
    private final int[] rawColumnIndices;
    /**
     * Median of every processed column in the transformed data, used to replace missing values of a query.
     */
    private final float[] transformedColumnMedians;
    /**
     * Lower winsorizing bound of every processed column, computed from the raw data; NaN if not winsorized.
     */
    private final float[] winsorizingLowerBounds;
    /**
     * Upper winsorizing bound of every processed column, computed from the raw data; NaN if not winsorized.
     */
    private final float[] winsorizingUpperBounds;
    /**
     * Preprocessed vector of every clustered molecule in the kernel space.
     */
    private final float[][] kernelVectors;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, preprocesses the clustered data.
     *
     * @param aRunResult clustering run whose feature space is searched
     * @throws NullPointerException if the run result is null
     */
    public SimilaritySearch(ClusteringRunResult aRunResult) throws NullPointerException {
        this.runResult = Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        FeatureMatrix tmpProcessed = aRunResult.processedMatrix();
        FeatureMatrix tmpRaw = aRunResult.rawMatrix();
        int tmpColumnCount = tmpProcessed.getColumnCount();
        this.isCosineKernel = aRunResult.parameters().algorithm() == Art2aClusterer.Art2aAlgorithm.COSINE;
        this.contrastEnhancementThreshold = (float) (1.0 / Math.sqrt(tmpColumnCount + aRunResult.parameters().offsetForContrastEnhancement()));
        this.columnMinima = new float[tmpColumnCount];
        this.columnMaxima = new float[tmpColumnCount];
        this.rawColumnIndices = new int[tmpColumnCount];
        this.transformedColumnMedians = new float[tmpColumnCount];
        this.winsorizingLowerBounds = new float[tmpColumnCount];
        this.winsorizingUpperBounds = new float[tmpColumnCount];
        //the processed columns are an order-preserving subset of the raw columns
        int tmpRawIndex = 0;
        for (int j = 0; j < tmpColumnCount; j++) {
            while (!tmpRaw.getFeatureName(tmpRawIndex).equals(tmpProcessed.getFeatureName(j))) {
                tmpRawIndex++;
            }
            this.rawColumnIndices[j] = tmpRawIndex;
            float[] tmpColumn = tmpProcessed.getColumn(j);
            float tmpMinimum = Float.POSITIVE_INFINITY;
            float tmpMaximum = Float.NEGATIVE_INFINITY;
            for (float tmpValue : tmpColumn) {
                tmpMinimum = Math.min(tmpMinimum, tmpValue);
                tmpMaximum = Math.max(tmpMaximum, tmpValue);
            }
            this.columnMinima[j] = tmpMinimum;
            this.columnMaxima[j] = tmpMaximum;
            float[] tmpSortedTransformed = SimilaritySearch.getSortedFiniteValues(aRunResult.transformedMatrix().getColumn(tmpRawIndex));
            this.transformedColumnMedians[j] = tmpSortedTransformed.length == 0
                    ? 0.0f : FeatureTransformer.getPercentile(tmpSortedTransformed, 0.5);
            this.winsorizingLowerBounds[j] = Float.NaN;
            this.winsorizingUpperBounds[j] = Float.NaN;
            if (this.getTransformation(tmpRawIndex) == FeatureTransformation.WINSORIZE) {
                float[] tmpSortedRaw = SimilaritySearch.getSortedFiniteValues(tmpRaw.getColumn(tmpRawIndex));
                if (tmpSortedRaw.length > 0) {
                    this.winsorizingLowerBounds[j] = FeatureTransformer.getPercentile(tmpSortedRaw, FeatureTransformer.WINSORIZE_LOWER_PERCENTILE);
                    this.winsorizingUpperBounds[j] = FeatureTransformer.getPercentile(tmpSortedRaw, FeatureTransformer.WINSORIZE_UPPER_PERCENTILE);
                }
            }
            tmpRawIndex++;
        }
        this.kernelVectors = new float[tmpProcessed.getRowCount()][];
        for (int i = 0; i < this.kernelVectors.length; i++) {
            float[] tmpRow = new float[tmpColumnCount];
            for (int j = 0; j < tmpColumnCount; j++) {
                tmpRow[j] = tmpProcessed.getValue(i, j);
            }
            this.kernelVectors[i] = this.preprocess(tmpRow);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Searches the molecules most similar to the clustered molecule with the given row index; the molecule itself is not
     * part of the result and its cluster is the one it was assigned to.
     *
     * @param aRowIndex row index of the query molecule in the clustered data
     * @param aNumberOfNeighbors maximum number of molecules to return, at least 1
     * @return search result
     * @throws IllegalArgumentException if the row index or number of neighbours is invalid
     */
    public SearchResult searchByRowIndex(int aRowIndex, int aNumberOfNeighbors) throws IllegalArgumentException {
        if (aRowIndex < 0 || aRowIndex >= this.kernelVectors.length) {
            throw new IllegalArgumentException("Row index out of range: " + aRowIndex);
        }
        float[] tmpQuery = this.kernelVectors[aRowIndex];
        int tmpCluster = this.runResult.clusteringResult().clusterLabels()[aRowIndex];
        double tmpClusterSimilarity = tmpCluster < 0 ? 0.0
                : this.getSimilarity(tmpQuery, this.runResult.clusteringResult().clusterVectors()[tmpCluster]);
        return new SearchResult(this.getNeighbors(tmpQuery, aNumberOfNeighbors, aRowIndex), tmpCluster, tmpClusterSimilarity,
                tmpCluster >= 0 && tmpClusterSimilarity >= this.runResult.clusteringResult().vigilance());
    }
    //
    /**
     * Searches the molecules most similar to a query given as raw feature vector (columns of the raw matrix of the run)
     * and assigns it to the cluster with the most similar cluster vector.
     *
     * @param aRawFeatureVector raw feature vector of the query, may contain NaN
     * @param aNumberOfNeighbors maximum number of molecules to return, at least 1
     * @return search result
     * @throws NullPointerException if the vector is null
     * @throws IllegalArgumentException if the vector length does not match the raw matrix or the number is invalid
     */
    public SearchResult searchByRawFeatureVector(float[] aRawFeatureVector, int aNumberOfNeighbors)
            throws NullPointerException, IllegalArgumentException {
        float[] tmpQuery = this.getKernelVector(aRawFeatureVector);
        int tmpBestCluster = -1;
        double tmpBestSimilarity = 0.0;
        if (!SimilaritySearch.isZeroVector(tmpQuery)) {
            float[][] tmpClusterVectors = this.runResult.clusteringResult().clusterVectors();
            tmpBestSimilarity = Double.NEGATIVE_INFINITY;
            for (int k = 0; k < tmpClusterVectors.length; k++) {
                double tmpSimilarity = this.getSimilarity(tmpQuery, tmpClusterVectors[k]);
                if (tmpSimilarity > tmpBestSimilarity) {
                    tmpBestSimilarity = tmpSimilarity;
                    tmpBestCluster = k;
                }
            }
        }
        return new SearchResult(this.getNeighbors(tmpQuery, aNumberOfNeighbors, -1), tmpBestCluster,
                Math.max(tmpBestSimilarity, 0.0),
                tmpBestCluster >= 0 && tmpBestSimilarity >= this.runResult.clusteringResult().vigilance());
    }
    //
    /**
     * Returns the similarity of two clustered molecules in the feature space of the clustering.
     *
     * @param aFirstRowIndex row index of the first molecule
     * @param aSecondRowIndex row index of the second molecule
     * @return similarity, 1 for identical feature vectors
     * @throws IllegalArgumentException if a row index is out of range
     */
    public double getSimilarityBetweenRows(int aFirstRowIndex, int aSecondRowIndex) throws IllegalArgumentException {
        if (aFirstRowIndex < 0 || aFirstRowIndex >= this.kernelVectors.length
                || aSecondRowIndex < 0 || aSecondRowIndex >= this.kernelVectors.length) {
            throw new IllegalArgumentException("Row index out of range.");
        }
        return this.getSimilarity(this.kernelVectors[aFirstRowIndex], this.kernelVectors[aSecondRowIndex]);
    }
    //
    /**
     * Returns the similarity of the clustered molecule with the given row index to the cluster vector (prototype) of its
     * cluster, i.e. how typical the molecule is for its cluster.
     *
     * @param aRowIndex row index of the molecule in the clustered data
     * @return similarity to the own cluster vector; 0 for molecules that could not be clustered
     * @throws IllegalArgumentException if the row index is out of range
     */
    public double getSimilarityToOwnClusterVector(int aRowIndex) throws IllegalArgumentException {
        if (aRowIndex < 0 || aRowIndex >= this.kernelVectors.length) {
            throw new IllegalArgumentException("Row index out of range: " + aRowIndex);
        }
        int tmpCluster = this.runResult.clusteringResult().clusterLabels()[aRowIndex];
        return tmpCluster < 0 ? 0.0
                : this.getSimilarity(this.kernelVectors[aRowIndex], this.runResult.clusteringResult().clusterVectors()[tmpCluster]);
    }
    //
    /**
     * Returns the vector of the clustered molecule with the given row index in the kernel space, i.e. the vector the
     * ART-2A clustering compared with the cluster vectors.
     *
     * @param aRowIndex row index of the molecule in the clustered data
     * @return copy of the kernel vector
     * @throws IllegalArgumentException if the row index is out of range
     */
    public float[] getKernelVectorOfRow(int aRowIndex) throws IllegalArgumentException {
        if (aRowIndex < 0 || aRowIndex >= this.kernelVectors.length) {
            throw new IllegalArgumentException("Row index out of range: " + aRowIndex);
        }
        return this.kernelVectors[aRowIndex].clone();
    }
    //
    /**
     * Transforms and preprocesses the given raw feature vector into the kernel space, like the clustered data.
     *
     * @param aRawFeatureVector raw feature vector (columns of the raw matrix of the run), may contain NaN
     * @return vector in the kernel space
     * @throws NullPointerException if the vector is null
     * @throws IllegalArgumentException if the vector length does not match the raw matrix
     */
    public float[] getKernelVector(float[] aRawFeatureVector) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aRawFeatureVector, "Given raw feature vector is null.");
        if (aRawFeatureVector.length != this.runResult.rawMatrix().getColumnCount()) {
            throw new IllegalArgumentException("The raw feature vector does not match the columns of the clustered data.");
        }
        float[] tmpProcessed = new float[this.rawColumnIndices.length];
        for (int j = 0; j < tmpProcessed.length; j++) {
            int tmpRawIndex = this.rawColumnIndices[j];
            float tmpValue = aRawFeatureVector[tmpRawIndex];
            FeatureTransformation tmpTransformation = this.getTransformation(tmpRawIndex);
            if (tmpTransformation == FeatureTransformation.WINSORIZE) {
                if (!Float.isNaN(tmpValue) && !Float.isNaN(this.winsorizingLowerBounds[j])) {
                    tmpValue = Math.min(Math.max(tmpValue, this.winsorizingLowerBounds[j]), this.winsorizingUpperBounds[j]);
                }
            } else if ((tmpTransformation == FeatureTransformation.LOG1P || tmpTransformation == FeatureTransformation.SQRT)
                    && tmpValue < 0.0f) {
                //not defined for negative values, treated as missing value
                tmpValue = Float.NaN;
            } else {
                float[] tmpSingleValue = new float[] {tmpValue};
                FeatureTransformer.transformColumn(tmpSingleValue, tmpTransformation, this.runResult.rawMatrix().getFeatureName(tmpRawIndex));
                tmpValue = tmpSingleValue[0];
            }
            //missing values are replaced by the median of the data, as done by the filter
            tmpProcessed[j] = Float.isFinite(tmpValue) ? tmpValue : this.transformedColumnMedians[j];
        }
        return this.preprocess(tmpProcessed);
    }
    //
    /**
     * Returns the similarity of two vectors in the kernel space: the cosine (dot product of unit vectors) for the
     * cosine kernel and 1 - squared distance / n for the Euclid kernel.
     *
     * @param aFirstVector first vector in the kernel space
     * @param aSecondVector second vector in the kernel space
     * @return similarity, 1 for identical vectors
     */
    public double getSimilarity(float[] aFirstVector, float[] aSecondVector) {
        double tmpSum = 0.0;
        if (this.isCosineKernel) {
            for (int j = 0; j < aFirstVector.length; j++) {
                tmpSum += (double) aFirstVector[j] * aSecondVector[j];
            }
            return tmpSum;
        }
        for (int j = 0; j < aFirstVector.length; j++) {
            double tmpDifference = aFirstVector[j] - aSecondVector[j];
            tmpSum += tmpDifference * tmpDifference;
        }
        return 1.0 - tmpSum / aFirstVector.length;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns the transformation that was applied to the given raw column.
     *
     * @param aRawColumnIndex raw column index
     * @return transformation
     */
    private FeatureTransformation getTransformation(int aRawColumnIndex) {
        return switch (this.runResult.rawMatrix().getFeatureType(aRawColumnIndex)) {
            case BIT -> FeatureTransformation.NONE;
            case COUNT -> this.runResult.parameters().countTransformation();
            case CONTINUOUS -> this.runResult.parameters().continuousTransformation();
        };
    }
    //
    /**
     * Applies the internal preprocessing of the kernel to a vector of processed (transformed and filtered) values:
     * min-max scaling with the data's minima and maxima (values outside the range of the data are clamped), then for
     * the cosine kernel L2 normalisation, contrast enhancement and renormalisation, for the Euclid kernel contrast
     * enhancement only. A vector that is zero after scaling stays zero.
     *
     * @param aProcessedVector processed values (IS CHANGED)
     * @return vector in the kernel space (the given array)
     */
    private float[] preprocess(float[] aProcessedVector) {
        for (int j = 0; j < aProcessedVector.length; j++) {
            float tmpRange = this.columnMaxima[j] - this.columnMinima[j];
            float tmpScaled = tmpRange > 0.0f ? (aProcessedVector[j] - this.columnMinima[j]) / tmpRange : 0.0f;
            aProcessedVector[j] = Math.min(1.0f, Math.max(0.0f, tmpScaled));
        }
        if (SimilaritySearch.isZeroVector(aProcessedVector)) {
            return aProcessedVector;
        }
        if (this.isCosineKernel) {
            SimilaritySearch.normalize(aProcessedVector);
        }
        boolean tmpIsChanged = false;
        for (int j = 0; j < aProcessedVector.length; j++) {
            if (aProcessedVector[j] > 0.0f && aProcessedVector[j] <= this.contrastEnhancementThreshold) {
                aProcessedVector[j] = 0.0f;
                tmpIsChanged = true;
            }
        }
        if (this.isCosineKernel && tmpIsChanged && !SimilaritySearch.isZeroVector(aProcessedVector)) {
            SimilaritySearch.normalize(aProcessedVector);
        }
        return aProcessedVector;
    }
    //
    /**
     * Returns the clustered molecules most similar to the given kernel vector, most similar first (ties by row index).
     *
     * @param aQueryVector query in the kernel space
     * @param aNumberOfNeighbors maximum number of molecules to return, at least 1
     * @param anExcludedRowIndex row index that is not returned (the query itself), -1 for none
     * @return list of neighbours
     * @throws IllegalArgumentException if the number of neighbours is smaller than 1
     */
    private List<Neighbor> getNeighbors(float[] aQueryVector, int aNumberOfNeighbors, int anExcludedRowIndex)
            throws IllegalArgumentException {
        if (aNumberOfNeighbors < 1) {
            throw new IllegalArgumentException("The number of similar molecules must be at least 1.");
        }
        List<Neighbor> tmpCandidates = new ArrayList<>(this.kernelVectors.length);
        for (int i = 0; i < this.kernelVectors.length; i++) {
            if (i != anExcludedRowIndex) {
                tmpCandidates.add(new Neighbor(i, this.getSimilarity(aQueryVector, this.kernelVectors[i])));
            }
        }
        tmpCandidates.sort(Comparator.comparingDouble(Neighbor::similarity).reversed().thenComparingInt(Neighbor::rowIndex));
        return List.copyOf(tmpCandidates.subList(0, Math.min(aNumberOfNeighbors, tmpCandidates.size())));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns whether all components of the given vector are zero.
     *
     * @param aVector vector
     * @return true for a zero vector
     */
    private static boolean isZeroVector(float[] aVector) {
        for (float tmpValue : aVector) {
            if (tmpValue != 0.0f) {
                return false;
            }
        }
        return true;
    }
    //
    /**
     * Scales the given non-zero vector to unit length.
     *
     * @param aVector vector (IS CHANGED)
     */
    private static void normalize(float[] aVector) {
        double tmpSum = 0.0;
        for (float tmpValue : aVector) {
            tmpSum += (double) tmpValue * tmpValue;
        }
        float tmpLength = (float) Math.sqrt(tmpSum);
        for (int j = 0; j < aVector.length; j++) {
            aVector[j] /= tmpLength;
        }
    }
    //
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
        float[] tmpSorted = Arrays.copyOf(tmpValues, tmpCount);
        Arrays.sort(tmpSorted);
        return tmpSorted;
    }
    //</editor-fold>
}
