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
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BooleanSupplier;
import java.util.stream.IntStream;

/**
 * Quality of the clusters of a clustering result and distances between them, measured on the molecules, not only on the
 * cluster centres. All values are calculated in the space the ART-2A clustering worked in: the kernel vectors of the
 * molecules are compared with the similarity of the kernel (see {@link SimilaritySearch#getSimilarity}); the distance is
 * 1 - similarity. All pairwise distances are calculated once (in parallel, without storing the distance matrix) and give:
 * <ul>
 *     <li>Silhouette (classic definition): for a molecule, a is the mean distance to the other molecules of its cluster
 *     and b the mean distance to the molecules of the nearest other cluster (the one with the smallest mean distance);
 *     s = (b - a) / max(a, b) in [-1, 1]. Near 1 the molecule is clearly assigned, near 0 it lies between two clusters,
 *     below 0 it is on average closer to another cluster. Molecules of single-molecule clusters get 0, as usual. Unlike
 *     a comparison of the centres, the silhouette takes the shape and extent of the clusters into account.</li>
 *     <li>Neighbour purity: share of the k nearest neighbours of a molecule (among all clustered molecules) that belong
 *     to its own cluster; a purely local measure, independent of the positions of the centres. The number of neighbours
 *     depends on the size of the cluster: k = min({@link #NUMBER_OF_NEIGHBOURS}, cluster size - 1). A cluster with s
 *     molecules can have at most s - 1 own neighbours, so with a fixed k small clusters could never reach 100 % (a
 *     perfectly separated cluster of 3 molecules would get 2 of 10 = 20 %); with k = s - 1 it gets 100 % if the other
 *     molecules of the cluster are its nearest neighbours. Clusters with at least {@link #NUMBER_OF_NEIGHBOURS} + 1
 *     molecules use the full k; single-molecule clusters have no neighbours and no purity (NaN).</li>
 *     <li>Contacts: for every pair of clusters, how many of the nearest neighbours of the molecules of one cluster lie in
 *     the other, i.e. which clusters touch or mix.</li>
 *     <li>Mean distance between clusters (average linkage): mean distance of all pairs of molecules of two clusters.</li>
 * </ul>
 * In addition, the mean similarity of the molecules to their own cluster vector (the cohesion as seen by ART-2A) and the
 * distances between the cluster vectors are given.
 * <p>
 * Random baseline: silhouette and neighbour purity are also calculated for {@link #NUMBER_OF_RANDOM_ASSIGNMENTS} random
 * assignments, in which the cluster labels are shuffled among the clustered molecules, so that the clusters keep their
 * sizes. Because the molecules and their distances stay the same, the baseline needs no further distance calculation.
 * In the chemical space, the silhouette is often low in absolute terms; compared with the baseline, it shows how much
 * better than chance the clusters are. The random neighbour purity of a cluster grows with its size (a cluster with half
 * of the molecules gets about 50 % by chance), so it is given per cluster, too. The same single pass over the distances
 * also evaluates other assignments of the same molecules, e.g. the clusterings of a vigilance scan, see
 * {@link #getOverallValues(ClusteringRunResult, List, BooleanSupplier)}. Molecules that could not be clustered are not considered. For more
 * than {@link #MAXIMUM_NUMBER_OF_QUERY_MOLECULES} clustered molecules, the per-molecule values are calculated for a
 * stratified random sample of the molecules (the same share of every cluster, at least one molecule), compared with all
 * molecules; the result tells whether it was sampled. With only one cluster, the values that need another cluster are
 * NaN.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusterQuality {
    //<editor-fold desc="Records" defaultstate="collapsed">
    /**
     * Quality of one cluster.
     *
     * @param clusterIndex index of the cluster
     * @param size number of molecules
     * @param meanSimilarityToOwnCentre mean similarity of the molecules to the cluster vector (cohesion seen by ART-2A)
     * @param meanDistanceWithinCluster mean distance between the molecules of the cluster, NaN for a single molecule
     * @param silhouette mean silhouette of the molecules, NaN if there is no other cluster
     * @param numberOfNegativeSilhouettes number of molecules with negative silhouette, i.e. on average closer to another
     *                                    cluster (counted among the evaluated molecules)
     * @param numberOfNeighbours number k of nearest neighbours used for the molecules of this cluster, min(
     *                           {@link #NUMBER_OF_NEIGHBOURS}, size - 1); 0 for a single molecule
     * @param neighbourPurity mean share of the k nearest neighbours of the molecules that belong to the cluster, NaN for
     *                        a single molecule
     * @param nearestClusterIndex cluster with the smallest mean distance to this cluster, -1 if there is none
     * @param meanDistanceToNearestCluster mean distance to the molecules of the nearest cluster, NaN if there is none
     * @param mostContactsClusterIndex other cluster in which most of the nearest neighbours of the molecules lie that are
     *                                 not in this cluster, -1 if all neighbours are in this cluster
     * @param contactShareWithMostContacts share of all nearest neighbours of the molecules that lie in that cluster
     * @param randomNeighbourPurity mean neighbour purity of a cluster of the same size in the random assignments, NaN for
     *                              a single molecule
     */
    public record ClusterQualityEntry(int clusterIndex, int size, double meanSimilarityToOwnCentre, double meanDistanceWithinCluster,
                                      double silhouette, int numberOfNegativeSilhouettes, int numberOfNeighbours, double neighbourPurity,
                                      int nearestClusterIndex, double meanDistanceToNearestCluster, int mostContactsClusterIndex,
                                      double contactShareWithMostContacts, double randomNeighbourPurity) {
    }
    //
    /**
     * Distance of a cluster to another cluster.
     *
     * @param clusterIndex index of the other cluster
     * @param meanDistance mean distance of all pairs of molecules of both clusters (average linkage)
     * @param centreDistance distance between the cluster vectors
     * @param neighbourContactShare share of the nearest neighbours of the molecules of the first cluster that lie in
     *                              the other cluster
     */
    public record ClusterDistance(int clusterIndex, double meanDistance, double centreDistance, double neighbourContactShare) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Maximum number of nearest neighbours of the neighbour purity and the contacts; the molecules of a cluster with s
     * molecules use k = min(NUMBER_OF_NEIGHBOURS, s - 1) neighbours, see the class description.
     */
    public static final int NUMBER_OF_NEIGHBOURS = 10;
    /**
     * Maximum number of molecules whose per-molecule values are calculated; above, a stratified sample is used.
     */
    public static final int MAXIMUM_NUMBER_OF_QUERY_MOLECULES = 5000;
    /**
     * Number of random assignments of the baseline.
     */
    public static final int NUMBER_OF_RANDOM_ASSIGNMENTS = 10;
    //</editor-fold>
    //
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Seed of the sample of the molecules.
     */
    private static final long SAMPLE_SEED = 1L;
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Quality of every cluster, by cluster index.
     */
    private final List<ClusterQualityEntry> entries;
    /**
     * Mean distances between the clusters (average linkage), clusters x clusters; NaN if not available.
     */
    private final double[][] meanDistances;
    /**
     * Distances between the cluster vectors, clusters x clusters.
     */
    private final double[][] centreDistances;
    /**
     * Share of the nearest neighbours of the molecules of a cluster (row) that lie in a cluster (column).
     */
    private final double[][] contactShares;
    /**
     * Silhouette of every molecule; NaN if not evaluated, not clustered or if there is only one cluster.
     */
    private final double[] moleculeSilhouettes;
    /**
     * Mean silhouette of all evaluated molecules.
     */
    private final double overallSilhouette;
    /**
     * Mean neighbour purity of all evaluated molecules.
     */
    private final double overallNeighbourPurity;
    /**
     * Mean of the overall silhouettes of the random assignments.
     */
    private final double randomOverallSilhouette;
    /**
     * Standard deviation of the overall silhouettes of the random assignments.
     */
    private final double randomOverallSilhouetteDeviation;
    /**
     * Mean of the overall neighbour purities of the random assignments.
     */
    private final double randomOverallNeighbourPurity;
    /**
     * Standard deviation of the overall neighbour purities of the random assignments.
     */
    private final double randomOverallNeighbourPurityDeviation;
    /**
     * Number of molecules whose values were calculated.
     */
    private final int numberOfEvaluatedMolecules;
    /**
     * True if the per-molecule values were calculated for a sample of the molecules.
     */
    private final boolean isSampled;
    //</editor-fold>
    //
    //<editor-fold desc="Constructors" defaultstate="collapsed">
    /**
     * Constructor, calculates the quality of all clusters of the given result.
     *
     * @param aRunResult clustering result
     * @throws NullPointerException if the result is null
     */
    public ClusterQuality(ClusteringRunResult aRunResult) throws NullPointerException {
        this(aRunResult, () -> false);
    }
    //
    /**
     * Constructor, calculates the quality of all clusters of the given result; the calculation can be cancelled.
     *
     * @param aRunResult clustering result
     * @param anIsCancelledSupplier returns true if the calculation should be stopped
     * @throws NullPointerException if a parameter is null
     * @throws CancellationException if the calculation was cancelled
     */
    public ClusterQuality(ClusteringRunResult aRunResult, BooleanSupplier anIsCancelledSupplier)
            throws NullPointerException, CancellationException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        SimilaritySearch tmpSearch = new SimilaritySearch(aRunResult);
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        float[][] tmpCentres = aRunResult.clusteringResult().clusterVectors();
        int tmpClusterCount = aRunResult.clusteringResult().numberOfClusters();
        //clustered molecules and their kernel vectors
        int[] tmpClusteredRows = IntStream.range(0, tmpLabels.length).filter(i -> tmpLabels[i] >= 0).toArray();
        float[][] tmpVectors = new float[tmpLabels.length][];
        for (int tmpRow : tmpClusteredRows) {
            tmpVectors[tmpRow] = tmpSearch.getKernelVectorOfRow(tmpRow);
        }
        int[] tmpSizes = new int[tmpClusterCount];
        for (int tmpRow : tmpClusteredRows) {
            tmpSizes[tmpLabels[tmpRow]]++;
        }
        //distances between the cluster vectors
        this.centreDistances = new double[tmpClusterCount][tmpClusterCount];
        for (int k = 0; k < tmpClusterCount; k++) {
            for (int l = k + 1; l < tmpClusterCount; l++) {
                double tmpDistance = 1.0 - tmpSearch.getSimilarity(tmpCentres[k], tmpCentres[l]);
                this.centreDistances[k][l] = tmpDistance;
                this.centreDistances[l][k] = tmpDistance;
            }
        }
        //molecules whose values are calculated, a stratified sample for large data sets
        int[] tmpQueryRows = ClusterQuality.selectQueryRows(tmpClusteredRows, tmpLabels, tmpSizes);
        this.isSampled = tmpQueryRows.length < tmpClusteredRows.length;
        this.numberOfEvaluatedMolecules = tmpQueryRows.length;
        //random assignments of the baseline: the labels are shuffled among the clustered molecules, the sizes stay
        int[][] tmpRandomLabels = new int[ClusterQuality.NUMBER_OF_RANDOM_ASSIGNMENTS][];
        int[][] tmpRandomSizes = new int[ClusterQuality.NUMBER_OF_RANDOM_ASSIGNMENTS][];
        for (int r = 0; r < tmpRandomLabels.length; r++) {
            tmpRandomLabels[r] = ClusterQuality.shuffleLabels(tmpLabels, tmpClusteredRows, new Random(ClusterQuality.SAMPLE_SEED + r + 1));
            tmpRandomSizes[r] = tmpSizes;
        }
        //per query molecule: mean distance to every cluster, the clusters of its nearest neighbours and the values of the
        // random assignments
        MoleculeValues[] tmpValues = IntStream.range(0, tmpQueryRows.length).parallel().mapToObj(q -> {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The calculation of the cluster quality was cancelled.");
            }
            return ClusterQuality.calculateMoleculeValues(tmpQueryRows[q], tmpClusteredRows, tmpLabels, tmpSizes, tmpVectors,
                    tmpSearch, tmpClusterCount, tmpRandomLabels, tmpRandomSizes);
        }).toArray(MoleculeValues[]::new);
        //aggregation per cluster
        double[] tmpOwnSimilaritySums = new double[tmpClusterCount];
        for (int tmpRow : tmpClusteredRows) {
            tmpOwnSimilaritySums[tmpLabels[tmpRow]] += tmpSearch.getSimilarity(tmpVectors[tmpRow], tmpCentres[tmpLabels[tmpRow]]);
        }
        int[] tmpEvaluatedCounts = new int[tmpClusterCount];
        double[] tmpSilhouetteSums = new double[tmpClusterCount];
        int[] tmpNegativeCounts = new int[tmpClusterCount];
        double[] tmpPuritySums = new double[tmpClusterCount];
        double[][] tmpLinkageSums = new double[tmpClusterCount][tmpClusterCount];
        double[][] tmpContactCounts = new double[tmpClusterCount][tmpClusterCount];
        double[] tmpNeighbourTotals = new double[tmpClusterCount];
        this.moleculeSilhouettes = new double[tmpLabels.length];
        Arrays.fill(this.moleculeSilhouettes, Double.NaN);
        double tmpOverallSilhouetteSum = 0.0;
        double tmpOverallPuritySum = 0.0;
        int tmpPurityCount = 0;
        int tmpSilhouetteCount = 0;
        for (int q = 0; q < tmpQueryRows.length; q++) {
            int tmpRow = tmpQueryRows[q];
            int tmpOwn = tmpLabels[tmpRow];
            MoleculeValues tmpMolecule = tmpValues[q];
            tmpEvaluatedCounts[tmpOwn]++;
            for (int k = 0; k < tmpClusterCount; k++) {
                tmpLinkageSums[tmpOwn][k] += tmpMolecule.meanDistances()[k];
            }
            for (int tmpNeighbourLabel : tmpMolecule.neighbourLabels()) {
                tmpContactCounts[tmpOwn][tmpNeighbourLabel]++;
                tmpNeighbourTotals[tmpOwn]++;
            }
            double tmpPurity = tmpMolecule.neighbourLabels().length == 0 ? Double.NaN
                    : (double) Arrays.stream(tmpMolecule.neighbourLabels()).filter(aLabel -> aLabel == tmpOwn).count()
                    / tmpMolecule.neighbourLabels().length;
            if (!Double.isNaN(tmpPurity)) {
                tmpPuritySums[tmpOwn] += tmpPurity;
                tmpOverallPuritySum += tmpPurity;
                tmpPurityCount++;
            }
            if (tmpClusterCount > 1) {
                double tmpSilhouette = ClusterQuality.getMoleculeSilhouette(tmpMolecule.meanDistances(), tmpOwn, tmpSizes);
                this.moleculeSilhouettes[tmpRow] = tmpSilhouette;
                tmpSilhouetteSums[tmpOwn] += tmpSilhouette;
                if (tmpSilhouette < 0.0) {
                    tmpNegativeCounts[tmpOwn]++;
                }
                tmpOverallSilhouetteSum += tmpSilhouette;
                tmpSilhouetteCount++;
            }
        }
        this.overallSilhouette = tmpSilhouetteCount == 0 ? Double.NaN : tmpOverallSilhouetteSum / tmpSilhouetteCount;
        //random baseline: overall values of every random assignment, purity per cluster over all random assignments
        double[][] tmpRandomOverallValues = ClusterQuality.aggregateOverallValues(tmpValues, tmpQueryRows, tmpRandomLabels, tmpClusterCount > 1);
        double[] tmpRandomPuritySums = new double[tmpClusterCount];
        int[] tmpRandomPurityCounts = new int[tmpClusterCount];
        for (int q = 0; q < tmpQueryRows.length; q++) {
            for (int r = 0; r < tmpRandomLabels.length; r++) {
                double tmpPurity = tmpValues[q].purities()[r];
                if (!Double.isNaN(tmpPurity)) {
                    tmpRandomPuritySums[tmpRandomLabels[r][tmpQueryRows[q]]] += tmpPurity;
                    tmpRandomPurityCounts[tmpRandomLabels[r][tmpQueryRows[q]]]++;
                }
            }
        }
        double[] tmpRandomSilhouettes = new double[tmpRandomOverallValues.length];
        double[] tmpRandomPurities = new double[tmpRandomOverallValues.length];
        for (int r = 0; r < tmpRandomOverallValues.length; r++) {
            tmpRandomSilhouettes[r] = tmpRandomOverallValues[r][0];
            tmpRandomPurities[r] = tmpRandomOverallValues[r][1];
        }
        this.randomOverallSilhouette = ClusterQuality.getMean(tmpRandomSilhouettes);
        this.randomOverallSilhouetteDeviation = ClusterQuality.getStandardDeviation(tmpRandomSilhouettes);
        this.randomOverallNeighbourPurity = ClusterQuality.getMean(tmpRandomPurities);
        this.randomOverallNeighbourPurityDeviation = ClusterQuality.getStandardDeviation(tmpRandomPurities);
        //molecules of single-molecule clusters have no neighbours and are left out
        this.overallNeighbourPurity = tmpPurityCount == 0 ? Double.NaN : tmpOverallPuritySum / tmpPurityCount;
        //mean distances between clusters, averaged over both directions (equal without sampling)
        this.meanDistances = new double[tmpClusterCount][tmpClusterCount];
        for (int k = 0; k < tmpClusterCount; k++) {
            for (int l = 0; l < tmpClusterCount; l++) {
                double tmpForward = tmpEvaluatedCounts[k] == 0 ? Double.NaN : tmpLinkageSums[k][l] / tmpEvaluatedCounts[k];
                double tmpBackward = tmpEvaluatedCounts[l] == 0 ? Double.NaN : tmpLinkageSums[l][k] / tmpEvaluatedCounts[l];
                this.meanDistances[k][l] = k == l ? tmpForward
                        : Double.isNaN(tmpForward) ? tmpBackward : Double.isNaN(tmpBackward) ? tmpForward : (tmpForward + tmpBackward) / 2.0;
            }
        }
        this.contactShares = new double[tmpClusterCount][tmpClusterCount];
        for (int k = 0; k < tmpClusterCount; k++) {
            for (int l = 0; l < tmpClusterCount; l++) {
                this.contactShares[k][l] = tmpNeighbourTotals[k] == 0.0 ? Double.NaN : tmpContactCounts[k][l] / tmpNeighbourTotals[k];
            }
        }
        this.entries = new ArrayList<>(tmpClusterCount);
        for (int k = 0; k < tmpClusterCount; k++) {
            int tmpNearest = -1;
            int tmpMostContacts = -1;
            for (int l = 0; l < tmpClusterCount; l++) {
                if (l == k) {
                    continue;
                }
                if (!Double.isNaN(this.meanDistances[k][l]) && (tmpNearest < 0 || this.meanDistances[k][l] < this.meanDistances[k][tmpNearest])) {
                    tmpNearest = l;
                }
                if (tmpContactCounts[k][l] > 0.0 && (tmpMostContacts < 0 || tmpContactCounts[k][l] > tmpContactCounts[k][tmpMostContacts])) {
                    tmpMostContacts = l;
                }
            }
            int tmpEvaluated = tmpEvaluatedCounts[k];
            boolean tmpHasOther = tmpClusterCount > 1 && tmpEvaluated > 0;
            this.entries.add(new ClusterQualityEntry(k, tmpSizes[k],
                    tmpSizes[k] == 0 ? Double.NaN : tmpOwnSimilaritySums[k] / tmpSizes[k],
                    tmpSizes[k] < 2 ? Double.NaN : this.meanDistances[k][k],
                    tmpHasOther ? tmpSilhouetteSums[k] / tmpEvaluated : Double.NaN,
                    tmpNegativeCounts[k],
                    ClusterQuality.getNumberOfNeighbours(tmpSizes[k]),
                    tmpEvaluated == 0 || ClusterQuality.getNumberOfNeighbours(tmpSizes[k]) == 0 ? Double.NaN : tmpPuritySums[k] / tmpEvaluated,
                    tmpNearest,
                    tmpNearest < 0 ? Double.NaN : this.meanDistances[k][tmpNearest],
                    tmpMostContacts,
                    tmpMostContacts < 0 ? 0.0 : this.contactShares[k][tmpMostContacts],
                    tmpRandomPurityCounts[k] == 0 ? Double.NaN : tmpRandomPuritySums[k] / tmpRandomPurityCounts[k]));
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the quality of every cluster in the order of the cluster indices.
     *
     * @return unmodifiable list of quality entries
     */
    public List<ClusterQualityEntry> getEntries() {
        return Collections.unmodifiableList(this.entries);
    }
    //
    /**
     * Returns the quality of the given cluster.
     *
     * @param aClusterIndex cluster index
     * @return quality entry
     * @throws IndexOutOfBoundsException if the index is out of range
     */
    public ClusterQualityEntry getEntry(int aClusterIndex) throws IndexOutOfBoundsException {
        return this.entries.get(aClusterIndex);
    }
    //
    /**
     * Returns the mean silhouette of all evaluated molecules, a measure of the quality of the whole clustering.
     *
     * @return overall silhouette in [-1, 1], NaN if there are fewer than two clusters
     */
    public double getOverallSilhouette() {
        return this.overallSilhouette;
    }
    //
    /**
     * Returns the mean neighbour purity of all evaluated molecules that have neighbours (molecules of single-molecule
     * clusters are left out); every molecule uses k = min({@link #NUMBER_OF_NEIGHBOURS}, size of its cluster - 1).
     *
     * @return overall neighbour purity in [0, 1], NaN if no molecule has neighbours
     */
    public double getOverallNeighbourPurity() {
        return this.overallNeighbourPurity;
    }
    //
    /**
     * Returns the mean of the overall silhouettes of the random assignments with the same cluster sizes.
     *
     * @return random overall silhouette, NaN if there are fewer than two clusters
     */
    public double getRandomOverallSilhouette() {
        return this.randomOverallSilhouette;
    }
    //
    /**
     * Returns the standard deviation of the overall silhouettes of the random assignments.
     *
     * @return standard deviation of the random overall silhouette
     */
    public double getRandomOverallSilhouetteDeviation() {
        return this.randomOverallSilhouetteDeviation;
    }
    //
    /**
     * Returns the mean of the overall neighbour purities of the random assignments with the same cluster sizes.
     *
     * @return random overall neighbour purity in [0, 1]
     */
    public double getRandomOverallNeighbourPurity() {
        return this.randomOverallNeighbourPurity;
    }
    //
    /**
     * Returns the standard deviation of the overall neighbour purities of the random assignments.
     *
     * @return standard deviation of the random overall neighbour purity
     */
    public double getRandomOverallNeighbourPurityDeviation() {
        return this.randomOverallNeighbourPurityDeviation;
    }
    //
    /**
     * Returns the silhouette of the molecule with the given row index.
     *
     * @param aRowIndex row index of the molecule in the clustered data
     * @return silhouette, NaN if the molecule was not evaluated, could not be clustered or there is only one cluster
     */
    public double getMoleculeSilhouette(int aRowIndex) {
        return this.moleculeSilhouettes[aRowIndex];
    }
    //
    /**
     * Returns the number of molecules whose values were calculated.
     *
     * @return number of evaluated molecules
     */
    public int getNumberOfEvaluatedMolecules() {
        return this.numberOfEvaluatedMolecules;
    }
    //
    /**
     * Returns whether the per-molecule values were calculated for a sample of the molecules.
     *
     * @return true if sampled
     */
    public boolean isSampled() {
        return this.isSampled;
    }
    //
    /**
     * Returns the mean distance between the molecules of the given clusters (average linkage).
     *
     * @param aFirstClusterIndex first cluster
     * @param aSecondClusterIndex second cluster
     * @return mean distance, NaN if not available
     */
    public double getMeanDistance(int aFirstClusterIndex, int aSecondClusterIndex) {
        return this.meanDistances[aFirstClusterIndex][aSecondClusterIndex];
    }
    //
    /**
     * Returns the distance between the cluster vectors of the given clusters.
     *
     * @param aFirstClusterIndex first cluster
     * @param aSecondClusterIndex second cluster
     * @return distance, 1 - similarity of the cluster vectors; 0 for the same cluster
     */
    public double getCentreDistance(int aFirstClusterIndex, int aSecondClusterIndex) {
        return this.centreDistances[aFirstClusterIndex][aSecondClusterIndex];
    }
    //
    /**
     * Returns the share of the nearest neighbours of the molecules of the first cluster that lie in the second cluster.
     *
     * @param aFromClusterIndex cluster whose molecules' neighbours are counted
     * @param aToClusterIndex cluster the neighbours lie in
     * @return contact share in [0, 1], NaN if not available
     */
    public double getContactShare(int aFromClusterIndex, int aToClusterIndex) {
        return this.contactShares[aFromClusterIndex][aToClusterIndex];
    }
    //
    /**
     * Returns the distances of the given cluster to all other clusters, the nearest (smallest mean distance) first.
     *
     * @param aClusterIndex cluster index
     * @return distances to the other clusters, ascending by mean distance
     */
    public List<ClusterDistance> getDistancesToOtherClusters(int aClusterIndex) {
        List<ClusterDistance> tmpDistances = new ArrayList<>(this.meanDistances.length);
        for (int l = 0; l < this.meanDistances.length; l++) {
            if (l != aClusterIndex) {
                tmpDistances.add(new ClusterDistance(l, this.meanDistances[aClusterIndex][l], this.centreDistances[aClusterIndex][l],
                        this.contactShares[aClusterIndex][l]));
            }
        }
        tmpDistances.sort(Comparator.comparingDouble((ClusterDistance aDistance) -> Double.isNaN(aDistance.meanDistance())
                        ? Double.POSITIVE_INFINITY : aDistance.meanDistance()).thenComparingInt(ClusterDistance::clusterIndex));
        return tmpDistances;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns the silhouette of a molecule from its mean distance to its own cluster and to the nearest other cluster.
     *
     * @param anOwnDistance mean distance a to the other molecules of the own cluster
     * @param anOtherDistance mean distance b to the molecules of the nearest other cluster
     * @return (b - a) / max(a, b), 0 if both distances are 0
     */
    public static double getSilhouette(double anOwnDistance, double anOtherDistance) {
        double tmpMaximum = Math.max(anOwnDistance, anOtherDistance);
        return tmpMaximum <= 0.0 ? 0.0 : (anOtherDistance - anOwnDistance) / tmpMaximum;
    }
    //
    /**
     * Returns the number of nearest neighbours used for the molecules of a cluster of the given size: at most
     * {@link #NUMBER_OF_NEIGHBOURS}, but not more than the other molecules of the cluster, so that a perfectly separated
     * small cluster reaches a neighbour purity of 100 %.
     *
     * @param aClusterSize number of molecules of the cluster
     * @return k = min(NUMBER_OF_NEIGHBOURS, cluster size - 1), at least 0
     */
    public static int getNumberOfNeighbours(int aClusterSize) {
        return Math.max(0, Math.min(ClusterQuality.NUMBER_OF_NEIGHBOURS, aClusterSize - 1));
    }
    //
    /**
     * Returns the overall silhouette and neighbour purity of other assignments of the molecules of the given result to
     * clusters, e.g. of the clusterings of a vigilance scan, which share the processed features of the result. All
     * assignments are evaluated in one pass over the distances, in the space of the result, on the molecules clustered in
     * the result (the same sample as the quality of the result for large data sets); a molecule that is not clustered
     * in an assignment is left out for it.
     *
     * @param aRunResult clustering result whose features and molecules the assignments share
     * @param aLabelArrayList cluster labels of every assignment, one per molecule of the result, -1 for not clustered
     * @param anIsCancelledSupplier returns true if the calculation should be stopped
     * @return per assignment the overall silhouette (NaN with fewer than two clusters) and the overall neighbour purity
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if a label array does not have one label per molecule
     * @throws CancellationException if the calculation was cancelled
     */
    public static double[][] getOverallValues(ClusteringRunResult aRunResult, List<int[]> aLabelArrayList, BooleanSupplier anIsCancelledSupplier)
            throws NullPointerException, IllegalArgumentException, CancellationException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aLabelArrayList, "Given label array list is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        int tmpClusterCount = aRunResult.clusteringResult().numberOfClusters();
        int[] tmpClusteredRows = IntStream.range(0, tmpLabels.length).filter(i -> tmpLabels[i] >= 0).toArray();
        int[] tmpSizes = new int[tmpClusterCount];
        for (int tmpRow : tmpClusteredRows) {
            tmpSizes[tmpLabels[tmpRow]]++;
        }
        //sizes of the clusters of every assignment, counted on the molecules of the comparison
        int[][] tmpAssignmentLabels = new int[aLabelArrayList.size()][];
        int[][] tmpAssignmentSizes = new int[aLabelArrayList.size()][];
        boolean[] tmpHasSeveralClusters = new boolean[aLabelArrayList.size()];
        for (int a = 0; a < tmpAssignmentLabels.length; a++) {
            int[] tmpAssignment = aLabelArrayList.get(a);
            if (tmpAssignment.length != tmpLabels.length) {
                throw new IllegalArgumentException("A label array does not have one label per molecule.");
            }
            int tmpMaximum = Arrays.stream(tmpAssignment).max().orElse(-1);
            int[] tmpAssignmentSize = new int[tmpMaximum + 1];
            for (int tmpRow : tmpClusteredRows) {
                if (tmpAssignment[tmpRow] >= 0) {
                    tmpAssignmentSize[tmpAssignment[tmpRow]]++;
                }
            }
            tmpAssignmentLabels[a] = tmpAssignment;
            tmpAssignmentSizes[a] = tmpAssignmentSize;
            tmpHasSeveralClusters[a] = Arrays.stream(tmpAssignmentSize).filter(aSize -> aSize > 0).count() > 1;
        }
        SimilaritySearch tmpSearch = new SimilaritySearch(aRunResult);
        float[][] tmpVectors = new float[tmpLabels.length][];
        for (int tmpRow : tmpClusteredRows) {
            tmpVectors[tmpRow] = tmpSearch.getKernelVectorOfRow(tmpRow);
        }
        int[] tmpQueryRows = ClusterQuality.selectQueryRows(tmpClusteredRows, tmpLabels, tmpSizes);
        MoleculeValues[] tmpValues = IntStream.range(0, tmpQueryRows.length).parallel().mapToObj(q -> {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The calculation of the cluster quality was cancelled.");
            }
            return ClusterQuality.calculateMoleculeValues(tmpQueryRows[q], tmpClusteredRows, tmpLabels, tmpSizes, tmpVectors,
                    tmpSearch, tmpClusterCount, tmpAssignmentLabels, tmpAssignmentSizes);
        }).toArray(MoleculeValues[]::new);
        double[][] tmpOverallValues = new double[tmpAssignmentLabels.length][];
        for (int a = 0; a < tmpAssignmentLabels.length; a++) {
            tmpOverallValues[a] = ClusterQuality.aggregateOverallValues(tmpValues, tmpQueryRows, new int[][] {tmpAssignmentLabels[a]},
                    tmpHasSeveralClusters[a], a)[0];
        }
        return tmpOverallValues;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private records and static methods" defaultstate="collapsed">
    /**
     * Values of one evaluated molecule.
     *
     * @param meanDistances mean distance to the molecules of every cluster (for the own cluster without the molecule
     *                      itself), NaN for a cluster without other molecules
     * @param neighbourLabels clusters of the nearest neighbours
     * @param silhouettes silhouette of the molecule in every additional assignment, NaN if it is not clustered there
     * @param purities neighbour purity of the molecule in every additional assignment, NaN if it is not clustered there
     *                 or has no neighbours
     */
    private record MoleculeValues(double[] meanDistances, int[] neighbourLabels, double[] silhouettes, double[] purities) {
    }
    //
    /**
     * Selects the molecules whose values are calculated: all clustered molecules, or for more than
     * {@link #MAXIMUM_NUMBER_OF_QUERY_MOLECULES} a stratified random sample (the same share of every cluster, at least one
     * molecule per cluster).
     *
     * @param aClusteredRows row indices of the clustered molecules
     * @param aLabels cluster labels
     * @param aSizes cluster sizes
     * @return row indices of the molecules to evaluate, ascending
     */
    private static int[] selectQueryRows(int[] aClusteredRows, int[] aLabels, int[] aSizes) {
        if (aClusteredRows.length <= ClusterQuality.MAXIMUM_NUMBER_OF_QUERY_MOLECULES) {
            return aClusteredRows;
        }
        double tmpFraction = (double) ClusterQuality.MAXIMUM_NUMBER_OF_QUERY_MOLECULES / aClusteredRows.length;
        List<List<Integer>> tmpMembers = new ArrayList<>(aSizes.length);
        for (int k = 0; k < aSizes.length; k++) {
            tmpMembers.add(new ArrayList<>(aSizes[k]));
        }
        for (int tmpRow : aClusteredRows) {
            tmpMembers.get(aLabels[tmpRow]).add(tmpRow);
        }
        Random tmpRandom = new Random(ClusterQuality.SAMPLE_SEED);
        List<Integer> tmpSelected = new ArrayList<>(ClusterQuality.MAXIMUM_NUMBER_OF_QUERY_MOLECULES + aSizes.length);
        for (List<Integer> tmpCluster : tmpMembers) {
            if (tmpCluster.isEmpty()) {
                continue;
            }
            Collections.shuffle(tmpCluster, tmpRandom);
            int tmpCount = Math.max(1, (int) Math.round(tmpCluster.size() * tmpFraction));
            tmpSelected.addAll(tmpCluster.subList(0, Math.min(tmpCount, tmpCluster.size())));
        }
        return tmpSelected.stream().mapToInt(Integer::intValue).sorted().toArray();
    }
    //
    /**
     * Calculates the distances of one molecule to all clustered molecules: the mean distance to every cluster and the
     * clusters of its nearest neighbours, and from the same distances its silhouette and neighbour purity in every given
     * additional assignment (e.g. a random one). The nearest neighbours do not depend on the assignment; every assignment
     * uses k = min({@link #NUMBER_OF_NEIGHBOURS}, size of the molecule's cluster in it - 1) of them.
     *
     * @param aRow row index of the molecule
     * @param aClusteredRows row indices of all clustered molecules
     * @param aLabels cluster labels
     * @param aSizes cluster sizes
     * @param aVectors kernel vectors by row index (null for not clustered molecules)
     * @param aSearch similarity search providing the similarity of the kernel
     * @param aClusterCount number of clusters
     * @param anAssignmentLabels cluster labels of every additional assignment, -1 for not clustered
     * @param anAssignmentSizes cluster sizes of every additional assignment
     * @return values of the molecule
     */
    private static MoleculeValues calculateMoleculeValues(int aRow, int[] aClusteredRows, int[] aLabels, int[] aSizes, float[][] aVectors,
                                                          SimilaritySearch aSearch, int aClusterCount, int[][] anAssignmentLabels,
                                                          int[][] anAssignmentSizes) {
        double[] tmpDistanceSums = new double[aClusterCount];
        double[][] tmpAssignmentDistanceSums = new double[anAssignmentLabels.length][];
        for (int a = 0; a < anAssignmentLabels.length; a++) {
            tmpAssignmentDistanceSums[a] = new double[anAssignmentSizes[a].length];
        }
        //nearest neighbours of all assignments, sorted ascending by distance
        double[] tmpNeighbourDistances = new double[ClusterQuality.NUMBER_OF_NEIGHBOURS];
        int[] tmpNeighbourRows = new int[ClusterQuality.NUMBER_OF_NEIGHBOURS];
        Arrays.fill(tmpNeighbourDistances, Double.POSITIVE_INFINITY);
        int tmpFound = 0;
        float[] tmpVector = aVectors[aRow];
        for (int tmpOther : aClusteredRows) {
            if (tmpOther == aRow) {
                continue;
            }
            double tmpDistance = 1.0 - aSearch.getSimilarity(tmpVector, aVectors[tmpOther]);
            tmpDistanceSums[aLabels[tmpOther]] += tmpDistance;
            for (int a = 0; a < anAssignmentLabels.length; a++) {
                int tmpLabel = anAssignmentLabels[a][tmpOther];
                if (tmpLabel >= 0) {
                    tmpAssignmentDistanceSums[a][tmpLabel] += tmpDistance;
                }
            }
            if (tmpDistance < tmpNeighbourDistances[ClusterQuality.NUMBER_OF_NEIGHBOURS - 1]) {
                int tmpPosition = Math.min(tmpFound, ClusterQuality.NUMBER_OF_NEIGHBOURS - 1);
                while (tmpPosition > 0 && tmpNeighbourDistances[tmpPosition - 1] > tmpDistance) {
                    tmpNeighbourDistances[tmpPosition] = tmpNeighbourDistances[tmpPosition - 1];
                    tmpNeighbourRows[tmpPosition] = tmpNeighbourRows[tmpPosition - 1];
                    tmpPosition--;
                }
                tmpNeighbourDistances[tmpPosition] = tmpDistance;
                tmpNeighbourRows[tmpPosition] = tmpOther;
                tmpFound = Math.min(tmpFound + 1, ClusterQuality.NUMBER_OF_NEIGHBOURS);
            }
        }
        //k depends on the size of the own cluster, so that a perfectly separated small cluster reaches 100 %
        int tmpOwn = aLabels[aRow];
        int tmpNeighbourCount = Math.min(tmpFound, ClusterQuality.getNumberOfNeighbours(aSizes[tmpOwn]));
        int[] tmpNeighbourLabels = new int[tmpNeighbourCount];
        for (int n = 0; n < tmpNeighbourCount; n++) {
            tmpNeighbourLabels[n] = aLabels[tmpNeighbourRows[n]];
        }
        double[] tmpMeanDistances = ClusterQuality.getMeanDistances(tmpDistanceSums, tmpOwn, aSizes);
        //the additional assignments
        double[] tmpSilhouettes = new double[anAssignmentLabels.length];
        double[] tmpPurities = new double[anAssignmentLabels.length];
        for (int a = 0; a < anAssignmentLabels.length; a++) {
            int tmpAssignmentOwn = anAssignmentLabels[a][aRow];
            if (tmpAssignmentOwn < 0) {
                tmpSilhouettes[a] = Double.NaN;
                tmpPurities[a] = Double.NaN;
                continue;
            }
            int[] tmpAssignmentSizes = anAssignmentSizes[a];
            tmpSilhouettes[a] = ClusterQuality.getMoleculeSilhouette(
                    ClusterQuality.getMeanDistances(tmpAssignmentDistanceSums[a], tmpAssignmentOwn, tmpAssignmentSizes),
                    tmpAssignmentOwn, tmpAssignmentSizes);
            //the first k neighbours that are clustered in the assignment
            int tmpK = ClusterQuality.getNumberOfNeighbours(tmpAssignmentSizes[tmpAssignmentOwn]);
            int tmpUsed = 0;
            int tmpSame = 0;
            for (int n = 0; n < tmpFound && tmpUsed < tmpK; n++) {
                int tmpLabel = anAssignmentLabels[a][tmpNeighbourRows[n]];
                if (tmpLabel >= 0) {
                    tmpUsed++;
                    if (tmpLabel == tmpAssignmentOwn) {
                        tmpSame++;
                    }
                }
            }
            tmpPurities[a] = tmpUsed == 0 ? Double.NaN : (double) tmpSame / tmpUsed;
        }
        return new MoleculeValues(tmpMeanDistances, tmpNeighbourLabels, tmpSilhouettes, tmpPurities);
    }
    //
    /**
     * Returns the mean distance of a molecule to every cluster from the sums of its distances.
     *
     * @param aDistanceSums sum of the distances to the molecules of every cluster
     * @param anOwnCluster own cluster of the molecule, whose molecule count excludes the molecule itself
     * @param aSizes cluster sizes
     * @return mean distance to every cluster, NaN for a cluster without other molecules
     */
    private static double[] getMeanDistances(double[] aDistanceSums, int anOwnCluster, int[] aSizes) {
        double[] tmpMeanDistances = new double[aDistanceSums.length];
        for (int k = 0; k < aDistanceSums.length; k++) {
            int tmpOthers = k == anOwnCluster ? aSizes[k] - 1 : aSizes[k];
            tmpMeanDistances[k] = tmpOthers <= 0 ? Double.NaN : aDistanceSums[k] / tmpOthers;
        }
        return tmpMeanDistances;
    }
    //
    /**
     * Returns the overall silhouette and neighbour purity of every additional assignment from the values of the
     * evaluated molecules: the means over the molecules clustered in the assignment.
     *
     * @param aValues values of the evaluated molecules
     * @param aQueryRows row indices of the evaluated molecules
     * @param anAssignmentLabels cluster labels of the assignments
     * @param anIsSilhouetteDefined false if the assignments have fewer than two clusters (silhouette NaN)
     * @return per assignment the overall silhouette and the overall neighbour purity
     */
    private static double[][] aggregateOverallValues(MoleculeValues[] aValues, int[] aQueryRows, int[][] anAssignmentLabels,
                                                     boolean anIsSilhouetteDefined) {
        return ClusterQuality.aggregateOverallValues(aValues, aQueryRows, anAssignmentLabels, anIsSilhouetteDefined, 0);
    }
    //
    /**
     * Returns the overall silhouette and neighbour purity of the given assignments, whose values are stored in the
     * molecule values from the given offset on.
     *
     * @param aValues values of the evaluated molecules
     * @param aQueryRows row indices of the evaluated molecules
     * @param anAssignmentLabels cluster labels of the assignments
     * @param anIsSilhouetteDefined false if the assignments have fewer than two clusters (silhouette NaN)
     * @param anOffset index of the first assignment in the arrays of the molecule values
     * @return per assignment the overall silhouette and the overall neighbour purity
     */
    private static double[][] aggregateOverallValues(MoleculeValues[] aValues, int[] aQueryRows, int[][] anAssignmentLabels,
                                                     boolean anIsSilhouetteDefined, int anOffset) {
        double[][] tmpOverallValues = new double[anAssignmentLabels.length][2];
        for (int a = 0; a < anAssignmentLabels.length; a++) {
            double tmpSilhouetteSum = 0.0;
            int tmpSilhouetteCount = 0;
            double tmpPuritySum = 0.0;
            int tmpPurityCount = 0;
            for (int q = 0; q < aQueryRows.length; q++) {
                if (anAssignmentLabels[a][aQueryRows[q]] < 0) {
                    continue;
                }
                double tmpSilhouette = aValues[q].silhouettes()[anOffset + a];
                if (!Double.isNaN(tmpSilhouette)) {
                    tmpSilhouetteSum += tmpSilhouette;
                    tmpSilhouetteCount++;
                }
                double tmpPurity = aValues[q].purities()[anOffset + a];
                if (!Double.isNaN(tmpPurity)) {
                    tmpPuritySum += tmpPurity;
                    tmpPurityCount++;
                }
            }
            tmpOverallValues[a][0] = !anIsSilhouetteDefined || tmpSilhouetteCount == 0 ? Double.NaN : tmpSilhouetteSum / tmpSilhouetteCount;
            tmpOverallValues[a][1] = tmpPurityCount == 0 ? Double.NaN : tmpPuritySum / tmpPurityCount;
        }
        return tmpOverallValues;
    }
    //
    /**
     * Returns a copy of the given labels in which the labels of the clustered molecules are shuffled, so that every
     * cluster keeps its size.
     *
     * @param aLabels cluster labels
     * @param aClusteredRows row indices of the clustered molecules
     * @param aRandom random generator
     * @return shuffled labels
     */
    private static int[] shuffleLabels(int[] aLabels, int[] aClusteredRows, Random aRandom) {
        int[] tmpShuffled = aLabels.clone();
        //Fisher-Yates shuffle of the labels of the clustered rows
        for (int k = aClusteredRows.length - 1; k > 0; k--) {
            int tmpOther = aRandom.nextInt(k + 1);
            int tmpLabel = tmpShuffled[aClusteredRows[k]];
            tmpShuffled[aClusteredRows[k]] = tmpShuffled[aClusteredRows[tmpOther]];
            tmpShuffled[aClusteredRows[tmpOther]] = tmpLabel;
        }
        return tmpShuffled;
    }
    //
    /**
     * Returns the mean of the given values, ignoring NaN.
     *
     * @param aValues values
     * @return mean, NaN if there is no value
     */
    private static double getMean(double[] aValues) {
        return Arrays.stream(aValues).filter(aValue -> !Double.isNaN(aValue)).average().orElse(Double.NaN);
    }
    //
    /**
     * Returns the standard deviation (population) of the given values, ignoring NaN.
     *
     * @param aValues values
     * @return standard deviation, NaN if there is no value
     */
    private static double getStandardDeviation(double[] aValues) {
        double tmpMean = ClusterQuality.getMean(aValues);
        if (Double.isNaN(tmpMean)) {
            return Double.NaN;
        }
        return Math.sqrt(Arrays.stream(aValues).filter(aValue -> !Double.isNaN(aValue))
                .map(aValue -> (aValue - tmpMean) * (aValue - tmpMean)).average().orElse(0.0));
    }
    //
    /**
     * Returns the silhouette of a molecule from its mean distances to all clusters; 0 for a molecule of a
     * single-molecule cluster or if there is no other non-empty cluster.
     *
     * @param aMeanDistances mean distance to every cluster
     * @param anOwnCluster own cluster
     * @param aSizes cluster sizes
     * @return silhouette
     */
    private static double getMoleculeSilhouette(double[] aMeanDistances, int anOwnCluster, int[] aSizes) {
        if (aSizes[anOwnCluster] < 2) {
            return 0.0;
        }
        double tmpNearestOther = Double.POSITIVE_INFINITY;
        for (int k = 0; k < aMeanDistances.length; k++) {
            if (k != anOwnCluster && !Double.isNaN(aMeanDistances[k])) {
                tmpNearestOther = Math.min(tmpNearestOther, aMeanDistances[k]);
            }
        }
        return Double.isInfinite(tmpNearestOther) ? 0.0 : ClusterQuality.getSilhouette(aMeanDistances[anOwnCluster], tmpNearestOther);
    }
    //</editor-fold>
}
