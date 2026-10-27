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

/**
 * Result of an ART-2A clustering of the rows of a {@link FeatureMatrix}, independent of the kernel used.
 *
 * @param clusterLabels cluster index of every row; -1 for rows that could not be clustered (zero-length vectors)
 * @param numberOfClusters number of detected clusters
 * @param representativeIndices row index of the representative of every cluster (the member most similar to the
 *                              cluster vector)
 * @param isOverflow true if the maximum number of clusters was reached; then the remaining vectors were assigned to the
 *                   most similar cluster regardless of the vigilance, and the user should be warned
 * @param isConverged true if the clustering converged before the maximum number of epochs was reached
 * @param numberOfEpochs number of epochs that were performed
 * @param vigilance vigilance parameter that was used
 * @param clusterVectors cluster vector (prototype) of every cluster in the kernel's internal feature space (min-max
 *                       scaled, for the cosine kernel also unit length and contrast-enhanced); used to assign new
 *                       molecules to a cluster like the ART-2A winner search
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record ClusteringResult(
        int[] clusterLabels,
        int numberOfClusters,
        int[] representativeIndices,
        boolean isOverflow,
        boolean isConverged,
        int numberOfEpochs,
        float vigilance,
        float[][] clusterVectors) {
    /**
     * Returns the number of rows in the given cluster.
     *
     * @param aClusterIndex cluster index
     * @return cluster size
     */
    public int getClusterSize(int aClusterIndex) {
        int tmpSize = 0;
        for (int tmpLabel : this.clusterLabels) {
            if (tmpLabel == aClusterIndex) {
                tmpSize++;
            }
        }
        return tmpSize;
    }
    //
    /**
     * Returns the number of rows that could not be clustered.
     *
     * @return number of rows with label -1
     */
    public int getNumberOfUnclusteredRows() {
        return this.getClusterSize(-1);
    }
}
