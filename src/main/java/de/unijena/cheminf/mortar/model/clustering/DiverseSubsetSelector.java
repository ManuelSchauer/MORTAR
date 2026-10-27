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
import java.util.Set;

/**
 * Selects molecules of a clustering result: a subset with the given number of molecules from every cluster, so that the
 * subset covers the whole chemical space of the data with few molecules (e.g. for screening or a reduced data set), or
 * the given number of molecules from one cluster. Two selection modes are available within a cluster:
 * <ul>
 *     <li>{@link SelectionMode#TYPICAL}: the molecules most similar to the cluster vector; they represent the cluster
 *     best but are similar to each other.</li>
 *     <li>{@link SelectionMode#DIVERSE}: the most typical molecule (cluster centre) first, then MaxMin selection: every
 *     further molecule is the one whose highest similarity to the already selected molecules is lowest, i.e. the
 *     molecules spread to the edges of the cluster in different directions. Note that outliers at the edge of a
 *     cluster are preferred.</li>
 * </ul>
 * All similarities are calculated in the feature space of the clustering (see {@link SimilaritySearch}). With one
 * molecule per cluster both modes select the cluster centre. Clusters smaller than the requested number contribute all
 * their molecules. Molecules that could not be clustered can be added.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class DiverseSubsetSelector {
    //<editor-fold desc="Enum SelectionMode" defaultstate="collapsed">
    /**
     * Enumeration of the ways to select molecules within a cluster.
     */
    public enum SelectionMode {
        /**
         * The molecules closest to the cluster centre (most typical).
         */
        TYPICAL,
        /**
         * The cluster centre, then MaxMin selection towards the edges of the cluster (most diverse).
         */
        DIVERSE;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record SelectedMolecule" defaultstate="collapsed">
    /**
     * A selected molecule.
     *
     * @param rowIndex row index of the molecule in the clustered data
     * @param clusterIndex cluster of the molecule, -1 if it could not be clustered
     * @param rankInCluster position of the molecule in the selection of its cluster, 1 = first selected (cluster centre)
     * @param similarityToClusterVector similarity to the cluster vector; 0 for molecules that could not be clustered
     */
    public record SelectedMolecule(int rowIndex, int clusterIndex, int rankInCluster, double similarityToClusterVector) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Default selection mode.
     */
    public static final SelectionMode SELECTION_MODE_DEFAULT = SelectionMode.DIVERSE;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private DiverseSubsetSelector() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Selects the given number of molecules from every cluster. The molecules are returned cluster by cluster in the
     * order of the cluster indices (0, 1, 2, ...), within a cluster in the order of selection; the not clusterable
     * molecules come last.
     *
     * @param aRunResult clustering result
     * @param aNumberPerCluster number of molecules per cluster, at least 1
     * @param anIsIncludingUnclusteredMolecules true to add the molecules that could not be clustered
     * @param aSelectionMode way to select the molecules within a cluster
     * @return selected molecules
     * @throws NullPointerException if the result or mode is null
     * @throws IllegalArgumentException if the number per cluster is smaller than 1
     */
    public static List<SelectedMolecule> select(ClusteringRunResult aRunResult, int aNumberPerCluster,
                                                boolean anIsIncludingUnclusteredMolecules, SelectionMode aSelectionMode)
            throws NullPointerException, IllegalArgumentException {
        return DiverseSubsetSelector.select(aRunResult, aNumberPerCluster, anIsIncludingUnclusteredMolecules, aSelectionMode, null);
    }
    //
    /**
     * Selects the given number of molecules from every given cluster, e.g. from the clusters marked by the user, in the
     * order of the cluster indices; see {@link #select(ClusteringRunResult, int, boolean, SelectionMode)}. If clusters
     * are given, the molecules that could not be clustered are not added.
     *
     * @param aRunResult clustering result
     * @param aNumberPerCluster number of molecules per cluster, at least 1
     * @param anIsIncludingUnclusteredMolecules true to add the molecules that could not be clustered (only without given
     *                                          clusters)
     * @param aSelectionMode way to select the molecules within a cluster
     * @param aClusterIndices indices of the clusters to select from; null or empty for all clusters
     * @return selected molecules
     * @throws NullPointerException if the result or mode is null
     * @throws IllegalArgumentException if the number per cluster is smaller than 1
     */
    public static List<SelectedMolecule> select(ClusteringRunResult aRunResult, int aNumberPerCluster,
                                                boolean anIsIncludingUnclusteredMolecules, SelectionMode aSelectionMode,
                                                Set<Integer> aClusterIndices)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aSelectionMode, "Given selection mode is null.");
        if (aNumberPerCluster < 1) {
            throw new IllegalArgumentException("The number of molecules per cluster must be at least 1.");
        }
        SimilaritySearch tmpSearch = new SimilaritySearch(aRunResult);
        List<SelectedMolecule> tmpSelection = new ArrayList<>();
        boolean tmpIsRestricted = aClusterIndices != null && !aClusterIndices.isEmpty();
        for (int k = 0; k < aRunResult.clusteringResult().numberOfClusters(); k++) {
            if (!tmpIsRestricted || aClusterIndices.contains(k)) {
                tmpSelection.addAll(DiverseSubsetSelector.selectFromCluster(aRunResult, tmpSearch, k, aNumberPerCluster, aSelectionMode));
            }
        }
        if (anIsIncludingUnclusteredMolecules && !tmpIsRestricted) {
            int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
            for (int i = 0; i < tmpLabels.length; i++) {
                if (tmpLabels[i] < 0) {
                    tmpSelection.add(new SelectedMolecule(i, -1, 1, 0.0));
                }
            }
        }
        return tmpSelection;
    }
    //
    /**
     * Selects the given number of molecules from one cluster, in the order of selection.
     *
     * @param aRunResult clustering result
     * @param aClusterIndex index of the cluster
     * @param aNumber number of molecules, at least 1; all members are returned if the cluster is smaller
     * @param aSelectionMode way to select the molecules
     * @return selected molecules
     * @throws NullPointerException if the result or mode is null
     * @throws IllegalArgumentException if the cluster index is out of range or the number is smaller than 1
     */
    public static List<SelectedMolecule> selectFromCluster(ClusteringRunResult aRunResult, int aClusterIndex, int aNumber,
                                                           SelectionMode aSelectionMode)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aSelectionMode, "Given selection mode is null.");
        if (aClusterIndex < 0 || aClusterIndex >= aRunResult.clusteringResult().numberOfClusters()) {
            throw new IllegalArgumentException("Cluster index out of range: " + aClusterIndex);
        }
        if (aNumber < 1) {
            throw new IllegalArgumentException("The number of molecules must be at least 1.");
        }
        return DiverseSubsetSelector.selectFromCluster(aRunResult, new SimilaritySearch(aRunResult), aClusterIndex, aNumber, aSelectionMode);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Selects the given number of molecules from one cluster with an existing similarity search.
     *
     * @param aRunResult clustering result
     * @param aSearch similarity search of the clustering result
     * @param aClusterIndex index of the cluster
     * @param aNumber number of molecules
     * @param aSelectionMode way to select the molecules
     * @return selected molecules in the order of selection
     */
    private static List<SelectedMolecule> selectFromCluster(ClusteringRunResult aRunResult, SimilaritySearch aSearch,
                                                            int aClusterIndex, int aNumber, SelectionMode aSelectionMode) {
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        List<Integer> tmpMembers = new ArrayList<>();
        for (int i = 0; i < tmpLabels.length; i++) {
            if (tmpLabels[i] == aClusterIndex) {
                tmpMembers.add(i);
            }
        }
        //most typical first, ties by row index
        tmpMembers.sort(Comparator.comparingDouble((Integer anIndex) -> -aSearch.getSimilarityToOwnClusterVector(anIndex))
                .thenComparingInt(anIndex -> anIndex));
        int tmpNumber = Math.min(aNumber, tmpMembers.size());
        List<Integer> tmpSelectedRows;
        if (aSelectionMode == SelectionMode.TYPICAL || tmpNumber <= 1) {
            tmpSelectedRows = tmpMembers.subList(0, tmpNumber);
        } else {
            tmpSelectedRows = DiverseSubsetSelector.selectMaxMin(aSearch, tmpMembers, tmpNumber);
        }
        List<SelectedMolecule> tmpSelection = new ArrayList<>(tmpNumber);
        for (int r = 0; r < tmpSelectedRows.size(); r++) {
            int tmpRow = tmpSelectedRows.get(r);
            tmpSelection.add(new SelectedMolecule(tmpRow, aClusterIndex, r + 1, aSearch.getSimilarityToOwnClusterVector(tmpRow)));
        }
        return tmpSelection;
    }
    //
    /**
     * MaxMin selection: starts with the first member (the cluster centre) and adds the member whose highest similarity to
     * the already selected members is lowest, until the given number is reached. Ties are broken by the order of the
     * members, i.e. by typicality.
     *
     * @param aSearch similarity search of the clustering result
     * @param aMembersByTypicality members of the cluster, most typical first
     * @param aNumber number of members to select, at least 1 and at most the number of members
     * @return selected row indices in the order of selection
     */
    private static List<Integer> selectMaxMin(SimilaritySearch aSearch, List<Integer> aMembersByTypicality, int aNumber) {
        int tmpSize = aMembersByTypicality.size();
        boolean[] tmpIsSelected = new boolean[tmpSize];
        //highest similarity of every member to the selected members
        double[] tmpMaximumSimilarity = new double[tmpSize];
        Arrays.fill(tmpMaximumSimilarity, Double.NEGATIVE_INFINITY);
        List<Integer> tmpSelected = new ArrayList<>(aNumber);
        int tmpNext = 0;
        while (tmpSelected.size() < aNumber) {
            tmpIsSelected[tmpNext] = true;
            int tmpNextRow = aMembersByTypicality.get(tmpNext);
            tmpSelected.add(tmpNextRow);
            int tmpBest = -1;
            for (int m = 0; m < tmpSize; m++) {
                if (tmpIsSelected[m]) {
                    continue;
                }
                tmpMaximumSimilarity[m] = Math.max(tmpMaximumSimilarity[m],
                        aSearch.getSimilarityBetweenRows(aMembersByTypicality.get(m), tmpNextRow));
                if (tmpBest < 0 || tmpMaximumSimilarity[m] < tmpMaximumSimilarity[tmpBest]) {
                    tmpBest = m;
                }
            }
            tmpNext = tmpBest;
        }
        return tmpSelected;
    }
    //</editor-fold>
}
