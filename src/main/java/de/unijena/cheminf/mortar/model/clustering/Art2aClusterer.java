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

import de.unijena.cheminf.clustering.art2a.Art2aEuclidKernel;
import de.unijena.cheminf.clustering.art2a.Art2aEuclidResult;
import de.unijena.cheminf.clustering.art2a.Art2aKernel;
import de.unijena.cheminf.clustering.art2a.Art2aResult;

import java.util.Arrays;
import java.util.Objects;

/**
 * Runs an ART-2A clustering of a {@link FeatureMatrix} with the kernels of the ART2a-Clustering-for-Java library and
 * converts the library result into a {@link ClusteringResult}.
 * <p>
 * Note: Both kernels min-max scale every column to [0,1] internally and then apply a contrast enhancement that sets
 * small components to zero, see the library documentation. The feature matrix must not contain NaN values and must
 * have at least two columns, i.e. it should be the output of {@link FeatureFilter}.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class Art2aClusterer {
    //<editor-fold desc="Enum Art2aAlgorithm" defaultstate="collapsed">
    /**
     * Enumeration of the available ART-2A kernels.
     */
    public enum Art2aAlgorithm {
        /**
         * Classic ART-2A: cosine similarity of L2-normalised, contrast-enhanced vectors. For binary vectors the cosine
         * equals the Ochiai coefficient. Vigilance is the minimum cosine between a vector and its cluster vector.
         */
        COSINE,
        /**
         * ART-2A-Euclid: squared Euclidean distance of contrast-enhanced vectors. Not recommended for sparse
         * fingerprint data, because the vigilance has nearly no effect there.
         */
        EUCLID;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Default maximum number of epochs (library default).
     */
    public static final int MAXIMUM_NUMBER_OF_EPOCHS_DEFAULT = 10;
    /**
     * Default convergence threshold of the cosine kernel, i.e. the minimum cosine between a cluster vector before and
     * after an epoch (library default).
     */
    public static final float CONVERGENCE_THRESHOLD_COSINE_DEFAULT = 0.99f;
    /**
     * Default convergence threshold of the Euclid kernel, i.e. the maximum Euclidean shift of a cluster vector during
     * an epoch (library default).
     */
    public static final float CONVERGENCE_THRESHOLD_EUCLID_DEFAULT = 0.1f;
    /**
     * Default learning parameter (library default).
     */
    public static final float LEARNING_PARAMETER_DEFAULT = 0.01f;
    /**
     * Default offset for the contrast enhancement threshold 1/sqrt(n + offset) (library default).
     */
    public static final float OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT = 1.0f;
    /**
     * Default seed of the random number generator that shuffles the order of the data vectors in every epoch.
     */
    public static final long RANDOM_SEED_DEFAULT = 1L;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private Art2aClusterer() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Clusters the rows of the given feature matrix.
     *
     * @param aFeatureMatrix feature matrix without NaN values and with at least two columns (IS NOT CHANGED)
     * @param anAlgorithm ART-2A kernel to use
     * @param aVigilance vigilance parameter in (0,1); higher values result in more and smaller clusters
     * @param aMaximumNumberOfClusters maximum number of clusters, at least 2; if it is reached, the remaining vectors
     *                                 are assigned to the most similar cluster regardless of the vigilance (overflow)
     * @param aMaximumNumberOfEpochs maximum number of epochs, at least 1
     * @param aConvergenceThreshold convergence threshold, meaning depends on the kernel (see default constants)
     * @param aLearningParameter learning parameter in (0,1)
     * @param anOffsetForContrastEnhancement offset for the contrast enhancement threshold, greater than 0
     * @param aRandomSeed seed of the random number generator, greater than 0
     * @param anIsParallelCalculation true: the winner cluster search is parallelized
     * @return clustering result
     * @throws NullPointerException if the feature matrix or algorithm is null
     * @throws IllegalArgumentException if a parameter is invalid or the matrix contains non-finite values
     * @throws Exception if the clustering fails
     */
    public static ClusteringResult cluster(
            FeatureMatrix aFeatureMatrix,
            Art2aAlgorithm anAlgorithm,
            float aVigilance,
            int aMaximumNumberOfClusters,
            int aMaximumNumberOfEpochs,
            float aConvergenceThreshold,
            float aLearningParameter,
            float anOffsetForContrastEnhancement,
            long aRandomSeed,
            boolean anIsParallelCalculation)
            throws NullPointerException, IllegalArgumentException, Exception {
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(anAlgorithm, "Given algorithm is null.");
        float[][] tmpData = aFeatureMatrix.getDataCopy();
        int tmpRowCount = tmpData.length;
        int[] tmpLabels = new int[tmpRowCount];
        //-1 marks vectors that could not be clustered (zero-length vectors after scaling)
        Arrays.fill(tmpLabels, -1);
        int tmpNumberOfClusters;
        int[] tmpRepresentatives;
        boolean tmpIsOverflow;
        boolean tmpIsConverged;
        int tmpNumberOfEpochs;
        float[][] tmpClusterVectors;
        switch (anAlgorithm) {
            case COSINE -> {
                Art2aKernel tmpKernel = new Art2aKernel(tmpData, aMaximumNumberOfClusters, aMaximumNumberOfEpochs,
                        aConvergenceThreshold, aLearningParameter, anOffsetForContrastEnhancement, aRandomSeed, true);
                Art2aResult tmpResult = tmpKernel.getClusterResult(aVigilance, anIsParallelCalculation);
                tmpNumberOfClusters = tmpResult.getNumberOfDetectedClusters();
                for (int k = 0; k < tmpNumberOfClusters; k++) {
                    for (int tmpIndex : tmpResult.getDataVectorIndicesOfCluster(k)) {
                        tmpLabels[tmpIndex] = k;
                    }
                }
                tmpRepresentatives = tmpResult.getRepresentativeIndicesOfClusters();
                tmpIsOverflow = tmpResult.isClusterOverflow();
                tmpIsConverged = tmpResult.isConverged();
                tmpNumberOfEpochs = tmpResult.getNumberOfEpochs();
                //copies, because the library returns its internal arrays
                tmpClusterVectors = new float[tmpNumberOfClusters][];
                for (int k = 0; k < tmpNumberOfClusters; k++) {
                    tmpClusterVectors[k] = tmpResult.getClusterVector(k).clone();
                }
            }
            case EUCLID -> {
                Art2aEuclidKernel tmpKernel = new Art2aEuclidKernel(tmpData, aMaximumNumberOfClusters, aMaximumNumberOfEpochs,
                        aConvergenceThreshold, aLearningParameter, anOffsetForContrastEnhancement, aRandomSeed, true);
                Art2aEuclidResult tmpResult = tmpKernel.getClusterResult(aVigilance, anIsParallelCalculation);
                tmpNumberOfClusters = tmpResult.getNumberOfDetectedClusters();
                for (int k = 0; k < tmpNumberOfClusters; k++) {
                    for (int tmpIndex : tmpResult.getDataVectorIndicesOfCluster(k)) {
                        tmpLabels[tmpIndex] = k;
                    }
                }
                tmpRepresentatives = tmpResult.getRepresentativeIndicesOfClusters();
                tmpIsOverflow = tmpResult.isClusterOverflow();
                tmpIsConverged = tmpResult.isConverged();
                tmpNumberOfEpochs = tmpResult.getNumberOfEpochs();
                tmpClusterVectors = new float[tmpNumberOfClusters][];
                for (int k = 0; k < tmpNumberOfClusters; k++) {
                    tmpClusterVectors[k] = tmpResult.getClusterVector(k).clone();
                }
            }
            default -> throw new IllegalArgumentException("Undefined ART-2A algorithm.");
        }
        return new ClusteringResult(tmpLabels, tmpNumberOfClusters, tmpRepresentatives, tmpIsOverflow,
                tmpIsConverged, tmpNumberOfEpochs, aVigilance, tmpClusterVectors);
    }
    //</editor-fold>
}
