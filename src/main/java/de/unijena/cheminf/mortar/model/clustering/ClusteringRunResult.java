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

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

/**
 * Complete result of one clustering run: the parameters, the clustered molecules, the feature matrix of every stage
 * of the pipeline (for the export of raw, transformed and processed data), the report of the removed features and the
 * clustering result itself. Row i of every matrix and label i of the clustering result belong to molecule i.
 *
 * @param parameters parameters of the run
 * @param molecules clustered molecules
 * @param rawMatrix feature matrix as calculated (fragment fingerprint block followed by descriptor block), may contain NaN
 * @param transformedMatrix feature matrix after the non-linear transformations, may contain NaN
 * @param processedMatrix feature matrix after filtering, i.e. the input of the clustering kernel
 * @param removedFeatures report of the features removed by the filter
 * @param clusteringResult result of the ART-2A clustering
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record ClusteringRunResult(
        ClusteringParameters parameters,
        List<MoleculeDataModel> molecules,
        FeatureMatrix rawMatrix,
        FeatureMatrix transformedMatrix,
        FeatureMatrix processedMatrix,
        List<RemovedFeature> removedFeatures,
        ClusteringResult clusteringResult) {
    //<editor-fold desc="Record ClusterSummary" defaultstate="collapsed">
    /**
     * Summary of one cluster for the display in the GUI.
     *
     * @param clusterIndex index of the cluster, -1 for the molecules that could not be clustered
     * @param members molecules of the cluster, the representative first
     * @param representative representative of the cluster, null for the not clusterable molecules
     */
    public record ClusterSummary(int clusterIndex, List<MoleculeDataModel> members, MoleculeDataModel representative) {
        /**
         * Returns the number of molecules in the cluster.
         *
         * @return cluster size
         */
        public int getSize() {
            return this.members.size();
        }
    }
    //</editor-fold>
    //
    /**
     * Compact constructor that checks the parameters for null.
     *
     * @throws NullPointerException if a parameter is null
     */
    public ClusteringRunResult {
        Objects.requireNonNull(parameters, "Given parameters are null.");
        Objects.requireNonNull(molecules, "Given molecule list is null.");
        Objects.requireNonNull(rawMatrix, "Given raw matrix is null.");
        Objects.requireNonNull(transformedMatrix, "Given transformed matrix is null.");
        Objects.requireNonNull(processedMatrix, "Given processed matrix is null.");
        Objects.requireNonNull(removedFeatures, "Given removed feature list is null.");
        Objects.requireNonNull(clusteringResult, "Given clustering result is null.");
        molecules = List.copyOf(molecules);
        removedFeatures = List.copyOf(removedFeatures);
    }
    //
    /**
     * Returns the cluster summaries ranked from the largest to the smallest cluster (ties in the order of the cluster
     * indices); the summary of the molecules that could not be clustered, if any, comes last.
     *
     * @return ranked list of cluster summaries
     */
    public List<ClusterSummary> getClusterSummariesRankedBySize() {
        List<ClusterSummary> tmpRankedSummaries = new ArrayList<>(this.getClusterSummaries());
        tmpRankedSummaries.sort(Comparator.comparing((ClusterSummary aCluster) -> aCluster.clusterIndex() < 0)
                .thenComparing(Comparator.comparingInt(ClusterSummary::getSize).reversed())
                .thenComparingInt(ClusterSummary::clusterIndex));
        return tmpRankedSummaries;
    }
    //
    /**
     * Returns one summary per cluster in the order of the cluster indices; if there are molecules that could not be
     * clustered, a summary with index -1 and without representative is appended.
     *
     * @return list of cluster summaries
     */
    public List<ClusterSummary> getClusterSummaries() {
        int[] tmpLabels = this.clusteringResult.clusterLabels();
        int tmpNumberOfClusters = this.clusteringResult.numberOfClusters();
        List<List<MoleculeDataModel>> tmpMembers = new ArrayList<>(tmpNumberOfClusters + 1);
        for (int k = 0; k <= tmpNumberOfClusters; k++) {
            tmpMembers.add(new ArrayList<>());
        }
        int[] tmpRepresentatives = this.clusteringResult.representativeIndices();
        //representatives first
        for (int k = 0; k < tmpNumberOfClusters; k++) {
            tmpMembers.get(k).add(this.molecules.get(tmpRepresentatives[k]));
        }
        for (int i = 0; i < tmpLabels.length; i++) {
            int tmpLabel = tmpLabels[i];
            if (tmpLabel < 0) {
                tmpMembers.get(tmpNumberOfClusters).add(this.molecules.get(i));
            } else if (i != tmpRepresentatives[tmpLabel]) {
                tmpMembers.get(tmpLabel).add(this.molecules.get(i));
            }
        }
        List<ClusterSummary> tmpSummaries = new ArrayList<>(tmpNumberOfClusters + 1);
        for (int k = 0; k < tmpNumberOfClusters; k++) {
            tmpSummaries.add(new ClusterSummary(k, tmpMembers.get(k), this.molecules.get(tmpRepresentatives[k])));
        }
        if (!tmpMembers.get(tmpNumberOfClusters).isEmpty()) {
            tmpSummaries.add(new ClusterSummary(-1, tmpMembers.get(tmpNumberOfClusters), null));
        }
        return tmpSummaries;
    }
}
