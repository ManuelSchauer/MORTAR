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
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BooleanSupplier;
import java.util.stream.IntStream;

/**
 * Principal component analysis (PCA) of a clustering result for its visualisation. The PCA is calculated on the vectors
 * the ART-2A clustering compared, i.e. the molecules in the kernel space (see {@link SimilaritySearch}), so that the
 * projection shows the space in which the clusters were formed; the cluster vectors are projected into the same space.
 * <p>
 * Only the first few components are needed, therefore no covariance matrix (number of features squared) is set up: the
 * components are calculated by subspace iteration on the centred data with a Rayleigh-Ritz step (eigen decomposition of
 * the small projected matrix by the Jacobi method), which needs only products of the data matrix with a few vectors.
 * The start vectors come from a random number generator with a fixed seed and the sign of every component is fixed
 * (its largest loading is positive), so the result is reproducible.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class PrincipalComponentAnalysis {
    //<editor-fold desc="Records" defaultstate="collapsed">
    /**
     * Result of a PCA.
     *
     * @param scores coordinates of every row (molecule) on the components, rows x components
     * @param clusterVectorScores coordinates of every cluster vector on the components, clusters x components; empty if
     *                            no cluster vectors were given
     * @param eigenvalues variance along every component, descending
     * @param explainedVarianceRatios share of the total variance explained by every component, in [0,1]
     * @param loadings weight of every feature (column of the data) in every component, components x features; the
     *                 orthonormal component vectors, with the sign fixed so that the largest weight is positive
     */
    public record PcaResult(double[][] scores, double[][] clusterVectorScores, double[] eigenvalues, double[] explainedVarianceRatios,
                            double[][] loadings) {
        /**
         * Returns the number of components.
         *
         * @return number of components
         */
        public int getNumberOfComponents() {
            return this.eigenvalues.length;
        }
        //
        /**
         * Returns the features with the largest absolute weight in the given component, i.e. the features the component
         * mainly stands for; a positive loading means that molecules with a higher value of the feature lie further in
         * the positive direction of the component.
         *
         * @param aComponent component index
         * @param aNumber maximum number of features
         * @return feature loadings, the largest absolute weight first
         */
        public List<FeatureLoading> getTopLoadings(int aComponent, int aNumber) {
            double[] tmpLoadings = this.loadings[aComponent];
            return IntStream.range(0, tmpLoadings.length).boxed()
                    .sorted(Comparator.comparingDouble((Integer aFeature) -> -Math.abs(tmpLoadings[aFeature])).thenComparingInt(aFeature -> aFeature))
                    .limit(aNumber)
                    .map(aFeature -> new FeatureLoading(aFeature, tmpLoadings[aFeature]))
                    .toList();
        }
    }
    //
    /**
     * Weight of one feature in a component.
     *
     * @param featureIndex index of the feature (column of the data)
     * @param loading weight of the feature in the component, in [-1, 1]
     */
    public record FeatureLoading(int featureIndex, double loading) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Default number of components calculated for the visualisation.
     */
    public static final int NUMBER_OF_COMPONENTS_DEFAULT = 3;
    //</editor-fold>
    //
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Maximum number of iterations of the subspace iteration.
     */
    private static final int MAXIMUM_NUMBER_OF_ITERATIONS = 500;
    /**
     * Relative change of the eigenvalues below which the iteration is converged.
     */
    private static final double CONVERGENCE_TOLERANCE = 1.0e-10;
    /**
     * Seed of the start vectors.
     */
    private static final long RANDOM_SEED = 1L;
    /**
     * Maximum number of sweeps of the Jacobi eigen decomposition.
     */
    private static final int MAXIMUM_NUMBER_OF_JACOBI_SWEEPS = 100;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private PrincipalComponentAnalysis() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Calculates the PCA of the given clustering result on the kernel vectors of the clustered molecules and projects its
     * cluster vectors into the same space.
     *
     * @param aRunResult clustering result
     * @param aNumberOfComponents number of components, at least 1; fewer are calculated if the data allows fewer
     * @param anIsCancelledSupplier returns true if the calculation should be stopped
     * @return PCA result
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of components is smaller than 1 or there are fewer than two molecules
     * @throws CancellationException if the calculation was cancelled
     */
    public static PcaResult calculate(ClusteringRunResult aRunResult, int aNumberOfComponents, BooleanSupplier anIsCancelledSupplier)
            throws NullPointerException, IllegalArgumentException, CancellationException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        SimilaritySearch tmpSearch = new SimilaritySearch(aRunResult);
        int tmpRowCount = aRunResult.molecules().size();
        double[][] tmpData = new double[tmpRowCount][];
        for (int i = 0; i < tmpRowCount; i++) {
            tmpData[i] = PrincipalComponentAnalysis.toDoubles(tmpSearch.getKernelVectorOfRow(i));
        }
        float[][] tmpClusterVectors = aRunResult.clusteringResult().clusterVectors();
        double[][] tmpClusterData = new double[tmpClusterVectors == null ? 0 : tmpClusterVectors.length][];
        for (int k = 0; k < tmpClusterData.length; k++) {
            tmpClusterData[k] = PrincipalComponentAnalysis.toDoubles(tmpClusterVectors[k]);
        }
        return PrincipalComponentAnalysis.calculate(tmpData, tmpClusterData, aNumberOfComponents, anIsCancelledSupplier);
    }
    //
    /**
     * Calculates the PCA of the given data and projects the given additional vectors (e.g. cluster vectors) into the same
     * space. The data is centred by its column means; the additional vectors are centred by the same means.
     *
     * @param aData data, rows x features, at least two rows of equal length
     * @param anAdditionalVectors vectors to project, each with the number of features; may be empty
     * @param aNumberOfComponents number of components, at least 1; at most the number of features and the number of rows
     *                            minus one are calculated
     * @param anIsCancelledSupplier returns true if the calculation should be stopped
     * @return PCA result
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of components is smaller than 1, there are fewer than two rows or
     *                                  the lengths of the vectors differ
     * @throws CancellationException if the calculation was cancelled
     */
    public static PcaResult calculate(double[][] aData, double[][] anAdditionalVectors, int aNumberOfComponents,
                                      BooleanSupplier anIsCancelledSupplier)
            throws NullPointerException, IllegalArgumentException, CancellationException {
        Objects.requireNonNull(aData, "Given data is null.");
        Objects.requireNonNull(anAdditionalVectors, "Given additional vectors are null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        if (aNumberOfComponents < 1) {
            throw new IllegalArgumentException("The number of components must be at least 1.");
        }
        if (aData.length < 2) {
            throw new IllegalArgumentException("At least two rows are needed for a PCA.");
        }
        int tmpRowCount = aData.length;
        int tmpColumnCount = aData[0].length;
        for (double[] tmpRow : aData) {
            if (tmpRow.length != tmpColumnCount) {
                throw new IllegalArgumentException("The rows of the data have different lengths.");
            }
        }
        for (double[] tmpVector : anAdditionalVectors) {
            if (tmpVector.length != tmpColumnCount) {
                throw new IllegalArgumentException("An additional vector does not match the number of features.");
            }
        }
        int tmpComponentCount = Math.max(1, Math.min(aNumberOfComponents, Math.min(tmpColumnCount, tmpRowCount - 1)));
        //centred data and total variance (sum of the column variances)
        double[] tmpMeans = new double[tmpColumnCount];
        for (double[] tmpRow : aData) {
            for (int j = 0; j < tmpColumnCount; j++) {
                tmpMeans[j] += tmpRow[j];
            }
        }
        for (int j = 0; j < tmpColumnCount; j++) {
            tmpMeans[j] /= tmpRowCount;
        }
        double[][] tmpCentred = new double[tmpRowCount][tmpColumnCount];
        double tmpTotalVariance = 0.0;
        for (int i = 0; i < tmpRowCount; i++) {
            for (int j = 0; j < tmpColumnCount; j++) {
                tmpCentred[i][j] = aData[i][j] - tmpMeans[j];
                tmpTotalVariance += tmpCentred[i][j] * tmpCentred[i][j];
            }
        }
        tmpTotalVariance /= tmpRowCount - 1;
        //subspace iteration with Rayleigh-Ritz step
        double[][] tmpBasis = new double[tmpComponentCount][tmpColumnCount];
        Random tmpRandom = new Random(PrincipalComponentAnalysis.RANDOM_SEED);
        for (double[] tmpVector : tmpBasis) {
            for (int j = 0; j < tmpColumnCount; j++) {
                tmpVector[j] = tmpRandom.nextGaussian();
            }
        }
        PrincipalComponentAnalysis.orthonormalize(tmpBasis);
        double[] tmpEigenvalues = new double[tmpComponentCount];
        double[] tmpPreviousEigenvalues = new double[tmpComponentCount];
        for (int tmpIteration = 0; tmpIteration < PrincipalComponentAnalysis.MAXIMUM_NUMBER_OF_ITERATIONS; tmpIteration++) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The PCA was cancelled.");
            }
            //W = X^T X V
            double[][] tmpProjected = PrincipalComponentAnalysis.multiplyCovariance(tmpCentred, tmpBasis);
            PrincipalComponentAnalysis.orthonormalize(tmpProjected);
            //Rayleigh-Ritz: eigen decomposition of the small matrix V^T X^T X V
            double[][] tmpSmall = PrincipalComponentAnalysis.projectCovariance(tmpCentred, tmpProjected);
            double[][] tmpRotation = new double[tmpComponentCount][tmpComponentCount];
            double[] tmpSmallEigenvalues = PrincipalComponentAnalysis.jacobiEigen(tmpSmall, tmpRotation);
            Integer[] tmpOrder = new Integer[tmpComponentCount];
            for (int c = 0; c < tmpComponentCount; c++) {
                tmpOrder[c] = c;
            }
            Arrays.sort(tmpOrder, (aFirst, aSecond) -> Double.compare(tmpSmallEigenvalues[aSecond], tmpSmallEigenvalues[aFirst]));
            for (int c = 0; c < tmpComponentCount; c++) {
                int tmpSource = tmpOrder[c];
                Arrays.fill(tmpBasis[c], 0.0);
                for (int r = 0; r < tmpComponentCount; r++) {
                    double tmpFactor = tmpRotation[r][tmpSource];
                    for (int j = 0; j < tmpColumnCount; j++) {
                        tmpBasis[c][j] += tmpFactor * tmpProjected[r][j];
                    }
                }
                tmpEigenvalues[c] = Math.max(0.0, tmpSmallEigenvalues[tmpSource]) / (tmpRowCount - 1);
            }
            boolean tmpIsConverged = tmpIteration > 0;
            for (int c = 0; c < tmpComponentCount && tmpIsConverged; c++) {
                double tmpScale = Math.max(Math.abs(tmpEigenvalues[c]), Double.MIN_NORMAL);
                tmpIsConverged = Math.abs(tmpEigenvalues[c] - tmpPreviousEigenvalues[c]) / tmpScale < PrincipalComponentAnalysis.CONVERGENCE_TOLERANCE;
            }
            System.arraycopy(tmpEigenvalues, 0, tmpPreviousEigenvalues, 0, tmpComponentCount);
            if (tmpIsConverged) {
                break;
            }
        }
        //fixed sign: the largest loading of every component is positive
        for (double[] tmpComponent : tmpBasis) {
            int tmpLargest = 0;
            for (int j = 1; j < tmpColumnCount; j++) {
                if (Math.abs(tmpComponent[j]) > Math.abs(tmpComponent[tmpLargest])) {
                    tmpLargest = j;
                }
            }
            if (tmpComponent[tmpLargest] < 0.0) {
                for (int j = 0; j < tmpColumnCount; j++) {
                    tmpComponent[j] = -tmpComponent[j];
                }
            }
        }
        double[] tmpRatios = new double[tmpComponentCount];
        for (int c = 0; c < tmpComponentCount; c++) {
            tmpRatios[c] = tmpTotalVariance > 0.0 ? tmpEigenvalues[c] / tmpTotalVariance : 0.0;
        }
        double[][] tmpScores = new double[tmpRowCount][];
        for (int i = 0; i < tmpRowCount; i++) {
            tmpScores[i] = PrincipalComponentAnalysis.project(tmpCentred[i], null, tmpBasis);
        }
        double[][] tmpAdditionalScores = new double[anAdditionalVectors.length][];
        for (int k = 0; k < anAdditionalVectors.length; k++) {
            tmpAdditionalScores[k] = PrincipalComponentAnalysis.project(anAdditionalVectors[k], tmpMeans, tmpBasis);
        }
        return new PcaResult(tmpScores, tmpAdditionalScores, tmpEigenvalues, tmpRatios, tmpBasis);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Converts the given float vector to doubles.
     *
     * @param aVector vector
     * @return double vector
     */
    private static double[] toDoubles(float[] aVector) {
        double[] tmpVector = new double[aVector.length];
        for (int j = 0; j < aVector.length; j++) {
            tmpVector[j] = aVector[j];
        }
        return tmpVector;
    }
    //
    /**
     * Returns the coordinates of the given vector on the given orthonormal components.
     *
     * @param aVector vector
     * @param aMeans means to subtract first, null if the vector is already centred
     * @param aBasis components, each of the length of the vector
     * @return coordinates
     */
    private static double[] project(double[] aVector, double[] aMeans, double[][] aBasis) {
        double[] tmpScores = new double[aBasis.length];
        for (int c = 0; c < aBasis.length; c++) {
            double tmpSum = 0.0;
            for (int j = 0; j < aVector.length; j++) {
                tmpSum += (aMeans == null ? aVector[j] : aVector[j] - aMeans[j]) * aBasis[c][j];
            }
            tmpScores[c] = tmpSum;
        }
        return tmpScores;
    }
    //
    /**
     * Returns X^T X v for every given vector v, via the products X v, without setting up X^T X.
     *
     * @param aData centred data X, rows x features
     * @param aVectors vectors v, each with the number of features
     * @return products, one per vector
     */
    private static double[][] multiplyCovariance(double[][] aData, double[][] aVectors) {
        int tmpColumnCount = aData[0].length;
        double[][] tmpResult = new double[aVectors.length][tmpColumnCount];
        for (double[] tmpRow : aData) {
            for (int c = 0; c < aVectors.length; c++) {
                double tmpDot = 0.0;
                for (int j = 0; j < tmpColumnCount; j++) {
                    tmpDot += tmpRow[j] * aVectors[c][j];
                }
                if (tmpDot != 0.0) {
                    for (int j = 0; j < tmpColumnCount; j++) {
                        tmpResult[c][j] += tmpDot * tmpRow[j];
                    }
                }
            }
        }
        return tmpResult;
    }
    //
    /**
     * Returns the small matrix V^T X^T X V of the given orthonormal vectors V, i.e. the scalar products of the projections
     * X V.
     *
     * @param aData centred data X, rows x features
     * @param aVectors orthonormal vectors V
     * @return symmetric matrix, vectors x vectors
     */
    private static double[][] projectCovariance(double[][] aData, double[][] aVectors) {
        int tmpCount = aVectors.length;
        double[][] tmpResult = new double[tmpCount][tmpCount];
        for (double[] tmpRow : aData) {
            double[] tmpProjection = PrincipalComponentAnalysis.project(tmpRow, null, aVectors);
            for (int a = 0; a < tmpCount; a++) {
                for (int b = 0; b < tmpCount; b++) {
                    tmpResult[a][b] += tmpProjection[a] * tmpProjection[b];
                }
            }
        }
        return tmpResult;
    }
    //
    /**
     * Orthonormalizes the given vectors in place by the modified Gram-Schmidt method; a vector that becomes (nearly) zero
     * is replaced by a unit vector orthogonal to the previous ones where possible.
     *
     * @param aVectors vectors (ARE CHANGED)
     */
    private static void orthonormalize(double[][] aVectors) {
        int tmpLength = aVectors[0].length;
        for (int c = 0; c < aVectors.length; c++) {
            for (int b = 0; b < c; b++) {
                double tmpDot = PrincipalComponentAnalysis.dot(aVectors[c], aVectors[b]);
                for (int j = 0; j < tmpLength; j++) {
                    aVectors[c][j] -= tmpDot * aVectors[b][j];
                }
            }
            double tmpNorm = Math.sqrt(PrincipalComponentAnalysis.dot(aVectors[c], aVectors[c]));
            if (tmpNorm < 1.0e-12) {
                //degenerate direction (e.g. no variance left): a unit vector, made orthogonal to the previous ones
                Arrays.fill(aVectors[c], 0.0);
                aVectors[c][c % tmpLength] = 1.0;
                for (int b = 0; b < c; b++) {
                    double tmpDot = PrincipalComponentAnalysis.dot(aVectors[c], aVectors[b]);
                    for (int j = 0; j < tmpLength; j++) {
                        aVectors[c][j] -= tmpDot * aVectors[b][j];
                    }
                }
                tmpNorm = Math.sqrt(PrincipalComponentAnalysis.dot(aVectors[c], aVectors[c]));
                if (tmpNorm < 1.0e-12) {
                    continue;
                }
            }
            for (int j = 0; j < tmpLength; j++) {
                aVectors[c][j] /= tmpNorm;
            }
        }
    }
    //
    /**
     * Returns the scalar product of the given vectors.
     *
     * @param aFirst first vector
     * @param aSecond second vector of the same length
     * @return scalar product
     */
    private static double dot(double[] aFirst, double[] aSecond) {
        double tmpSum = 0.0;
        for (int j = 0; j < aFirst.length; j++) {
            tmpSum += aFirst[j] * aSecond[j];
        }
        return tmpSum;
    }
    //
    /**
     * Eigen decomposition of a small symmetric matrix by the cyclic Jacobi method.
     *
     * @param aMatrix symmetric matrix (IS CHANGED)
     * @param anEigenvectors receives the eigenvectors as columns, must have the size of the matrix (IS CHANGED)
     * @return eigenvalues in the order of the eigenvector columns
     */
    private static double[] jacobiEigen(double[][] aMatrix, double[][] anEigenvectors) {
        int tmpSize = aMatrix.length;
        for (int i = 0; i < tmpSize; i++) {
            Arrays.fill(anEigenvectors[i], 0.0);
            anEigenvectors[i][i] = 1.0;
        }
        for (int tmpSweep = 0; tmpSweep < PrincipalComponentAnalysis.MAXIMUM_NUMBER_OF_JACOBI_SWEEPS; tmpSweep++) {
            double tmpOffDiagonal = 0.0;
            for (int p = 0; p < tmpSize; p++) {
                for (int q = p + 1; q < tmpSize; q++) {
                    tmpOffDiagonal += aMatrix[p][q] * aMatrix[p][q];
                }
            }
            if (tmpOffDiagonal < 1.0e-30) {
                break;
            }
            for (int p = 0; p < tmpSize; p++) {
                for (int q = p + 1; q < tmpSize; q++) {
                    if (Math.abs(aMatrix[p][q]) < 1.0e-300) {
                        continue;
                    }
                    double tmpTheta = (aMatrix[q][q] - aMatrix[p][p]) / (2.0 * aMatrix[p][q]);
                    double tmpT = Math.signum(tmpTheta) / (Math.abs(tmpTheta) + Math.sqrt(tmpTheta * tmpTheta + 1.0));
                    if (tmpTheta == 0.0) {
                        tmpT = 1.0;
                    }
                    double tmpCos = 1.0 / Math.sqrt(tmpT * tmpT + 1.0);
                    double tmpSin = tmpT * tmpCos;
                    for (int k = 0; k < tmpSize; k++) {
                        double tmpKp = aMatrix[k][p];
                        double tmpKq = aMatrix[k][q];
                        aMatrix[k][p] = tmpCos * tmpKp - tmpSin * tmpKq;
                        aMatrix[k][q] = tmpSin * tmpKp + tmpCos * tmpKq;
                    }
                    for (int k = 0; k < tmpSize; k++) {
                        double tmpPk = aMatrix[p][k];
                        double tmpQk = aMatrix[q][k];
                        aMatrix[p][k] = tmpCos * tmpPk - tmpSin * tmpQk;
                        aMatrix[q][k] = tmpSin * tmpPk + tmpCos * tmpQk;
                    }
                    for (int k = 0; k < tmpSize; k++) {
                        double tmpKp = anEigenvectors[k][p];
                        double tmpKq = anEigenvectors[k][q];
                        anEigenvectors[k][p] = tmpCos * tmpKp - tmpSin * tmpKq;
                        anEigenvectors[k][q] = tmpSin * tmpKp + tmpCos * tmpKq;
                    }
                }
            }
        }
        double[] tmpEigenvalues = new double[tmpSize];
        for (int i = 0; i < tmpSize; i++) {
            tmpEigenvalues[i] = aMatrix[i][i];
        }
        return tmpEigenvalues;
    }
    //</editor-fold>
}
