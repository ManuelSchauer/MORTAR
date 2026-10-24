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
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorCalculator;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.scaffold.ScaffoldGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.CancellationException;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * Chemical check of the clusters of a clustering result that is independent of the features the molecules were
 * clustered with: it answers whether the molecules of a cluster are also chemically similar when measured otherwise.
 * <ul>
 *     <li>Mean Tanimoto similarity of the molecules of a cluster on a bit fingerprint chosen by the user (e.g. ECFP4),
 *     over all pairs of molecules of the cluster (a random sample of {@link #MAXIMUM_NUMBER_OF_MOLECULES_PER_CLUSTER}
 *     molecules for larger clusters).</li>
 *     <li>Share of the molecules of a cluster that have its most frequent Bemis-Murcko scaffold (CDK scaffold generator,
 *     mode SCAFFOLD: rings and linkers with exocyclic double bonds); molecules without ring have the empty, "acyclic"
 *     scaffold.</li>
 * </ul>
 * Both values are compared with a random baseline: the cluster labels are shuffled among the clustered molecules, so
 * that the clusters keep their sizes, and the values are calculated again; this is repeated with different seeds.
 * The baseline matters especially for the scaffold share, which is high for small clusters by chance alone (a cluster
 * of two molecules has a share of at least 50 %).
 * <p>
 * The check is only independent if the fingerprint was not used for the clustering, see
 * {@link #getIndependence(Descriptor, ClusteringParameters)}. Clusters with fewer than two molecules with valid
 * values have no value (NaN); the overall values are means over the clusters weighted by their number of molecules.
 * Molecules whose fingerprint or scaffold cannot be calculated are left out. The molecules that could not be clustered
 * are not part of the check.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class ClusterChemistryCheck {
    //<editor-fold desc="Enum Independence" defaultstate="collapsed">
    /**
     * Independence of a fingerprint from the features of a clustering run.
     */
    public enum Independence {
        /**
         * The fingerprint was not used for the clustering and is not related to a fingerprint used for it.
         */
        INDEPENDENT,
        /**
         * The fingerprint was not used for the clustering, but a related fingerprint was: both are circular fingerprints
         * (ECFP or FCFP of another radius), which share many substructures.
         */
        RELATED,
        /**
         * The fingerprint was used for the clustering, so the check only repeats what the clustering optimised.
         */
        USED
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record ChemistryCheckResult" defaultstate="collapsed">
    /**
     * Result of the chemical check, the arrays by cluster index.
     *
     * @param fingerprint fingerprint the Tanimoto similarities were calculated with
     * @param meanTanimotos mean Tanimoto similarity of the molecules of every cluster, NaN if it has fewer than two
     *                      molecules with a fingerprint
     * @param randomMeanTanimotos mean Tanimoto similarity of every cluster in the random assignments with the same cluster
     *                            sizes, averaged over the random runs
     * @param scaffoldShares share of the molecules of every cluster with its most frequent scaffold, NaN if it has fewer
     *                       than two molecules with a scaffold
     * @param randomScaffoldShares scaffold share of every cluster in the random assignments, averaged over the runs
     * @param mostFrequentScaffolds unique SMILES of the most frequent scaffold of every cluster, "" for acyclic molecules,
     *                              null if the cluster has no molecule with a scaffold
     * @param overallMeanTanimoto mean Tanimoto similarity over all clusters, weighted by their number of molecules
     * @param overallRandomMeanTanimoto mean of the overall mean Tanimoto similarity of the random runs
     * @param overallRandomMeanTanimotoDeviation standard deviation of the overall mean Tanimoto similarity of the random runs
     * @param overallScaffoldShare scaffold share over all clusters, weighted by their number of molecules
     * @param overallRandomScaffoldShare mean of the overall scaffold share of the random runs
     * @param overallRandomScaffoldShareDeviation standard deviation of the overall scaffold share of the random runs
     * @param numberOfRandomRuns number of random assignments
     * @param numberOfFailedMolecules number of clustered molecules whose fingerprint or scaffold could not be calculated
     */
    public record ChemistryCheckResult(Descriptor fingerprint, double[] meanTanimotos, double[] randomMeanTanimotos,
                                       double[] scaffoldShares, double[] randomScaffoldShares, String[] mostFrequentScaffolds,
                                       double overallMeanTanimoto, double overallRandomMeanTanimoto,
                                       double overallRandomMeanTanimotoDeviation, double overallScaffoldShare,
                                       double overallRandomScaffoldShare, double overallRandomScaffoldShareDeviation,
                                       int numberOfRandomRuns, int numberOfFailedMolecules) {
        /**
         * Returns the mean Tanimoto similarity of the given cluster.
         *
         * @param aClusterIndex cluster index
         * @return mean Tanimoto similarity in [0, 1] or NaN
         */
        public double getMeanTanimoto(int aClusterIndex) {
            return this.meanTanimotos[aClusterIndex];
        }
        //
        /**
         * Returns the mean Tanimoto similarity of the given cluster in the random assignments.
         *
         * @param aClusterIndex cluster index
         * @return mean Tanimoto similarity in [0, 1] or NaN
         */
        public double getRandomMeanTanimoto(int aClusterIndex) {
            return this.randomMeanTanimotos[aClusterIndex];
        }
        //
        /**
         * Returns the share of the molecules of the given cluster with its most frequent scaffold.
         *
         * @param aClusterIndex cluster index
         * @return share in [0, 1] or NaN
         */
        public double getScaffoldShare(int aClusterIndex) {
            return this.scaffoldShares[aClusterIndex];
        }
        //
        /**
         * Returns the scaffold share of the given cluster in the random assignments.
         *
         * @param aClusterIndex cluster index
         * @return share in [0, 1] or NaN
         */
        public double getRandomScaffoldShare(int aClusterIndex) {
            return this.randomScaffoldShares[aClusterIndex];
        }
        //
        /**
         * Returns the unique SMILES of the most frequent scaffold of the given cluster.
         *
         * @param aClusterIndex cluster index
         * @return SMILES, "" for acyclic molecules, null if the cluster has no molecule with a scaffold
         */
        public String getMostFrequentScaffold(int aClusterIndex) {
            return this.mostFrequentScaffolds[aClusterIndex];
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record ClusterValues" defaultstate="collapsed">
    /**
     * Values of one assignment of the molecules to the clusters.
     *
     * @param meanTanimotos mean Tanimoto similarity by cluster index
     * @param scaffoldShares share of the most frequent scaffold by cluster index
     * @param mostFrequentScaffoldIds id of the most frequent scaffold by cluster index, -1 if there is none
     * @param overallMeanTanimoto weighted mean of the Tanimoto similarities
     * @param overallScaffoldShare weighted mean of the scaffold shares
     */
    private record ClusterValues(double[] meanTanimotos, double[] scaffoldShares, int[] mostFrequentScaffoldIds,
                                 double overallMeanTanimoto, double overallScaffoldShare) {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Default number of random assignments of the baseline.
     */
    public static final int NUMBER_OF_RANDOM_RUNS_DEFAULT = 10;
    /**
     * Maximum number of molecules of a cluster whose pairwise Tanimoto similarities are calculated; a random sample is
     * taken from larger clusters, which limits the number of pairs to about half a million per cluster.
     */
    public static final int MAXIMUM_NUMBER_OF_MOLECULES_PER_CLUSTER = 1000;
    /**
     * Fingerprints proposed for the check in this order, the first one that was not used for the clustering is the
     * default (see {@link #getDefaultFingerprint(ClusteringParameters)}).
     */
    public static final List<Descriptor> PREFERRED_FINGERPRINTS = List.of(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4,
            Descriptor.MACCS_FINGERPRINTER, Descriptor.PUBCHEM_FINGERPRINTER);
    //</editor-fold>
    //
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusterChemistryCheck.class.getName());
    /**
     * Number of molecules whose fingerprints are calculated at once, so that the float matrix of the descriptor engine
     * stays small; the fingerprints are kept as bits.
     */
    private static final int FINGERPRINT_BATCH_SIZE = 1000;
    /**
     * Number of molecules whose scaffolds are generated by one scaffold generator instance (instances are not
     * thread-safe, the batches run in parallel).
     */
    private static final int SCAFFOLD_BATCH_SIZE = 500;
    /**
     * Prefix of the enum constant names of the circular fingerprints.
     */
    private static final String CIRCULAR_FINGERPRINT_PREFIX = "CIRCULAR_FINGERPRINTER_";
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private ClusterChemistryCheck() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns the bit fingerprints that can be chosen for the check (all fingerprints without 3D coordinates).
     *
     * @return fingerprints
     */
    public static Descriptor[] getBitFingerprints() {
        return Descriptor.getAllFingerprints();
    }
    //
    /**
     * Returns whether the given fingerprint is independent of the features of the given clustering run: USED if it is
     * one of its descriptors, RELATED if both it and one of its descriptors are circular fingerprints, else INDEPENDENT.
     *
     * @param aFingerprint fingerprint for the check
     * @param aParameters parameters of the clustering run
     * @return independence of the fingerprint
     * @throws NullPointerException if a parameter is null
     */
    public static Independence getIndependence(Descriptor aFingerprint, ClusteringParameters aParameters) throws NullPointerException {
        Objects.requireNonNull(aFingerprint, "Given fingerprint is null.");
        Objects.requireNonNull(aParameters, "Given clustering parameters are null.");
        boolean tmpIsRelated = false;
        for (Descriptor tmpDescriptor : aParameters.descriptors()) {
            if (tmpDescriptor == aFingerprint) {
                return Independence.USED;
            }
            if (ClusterChemistryCheck.isCircularFingerprint(tmpDescriptor) && ClusterChemistryCheck.isCircularFingerprint(aFingerprint)) {
                tmpIsRelated = true;
            }
        }
        return tmpIsRelated ? Independence.RELATED : Independence.INDEPENDENT;
    }
    //
    /**
     * Returns the fingerprint proposed for the check of the given clustering run: the first of
     * {@link #PREFERRED_FINGERPRINTS} that is independent, otherwise the first independent bit fingerprint, otherwise
     * the first preferred one.
     *
     * @param aParameters parameters of the clustering run
     * @return proposed fingerprint
     * @throws NullPointerException if the parameters are null
     */
    public static Descriptor getDefaultFingerprint(ClusteringParameters aParameters) throws NullPointerException {
        for (Descriptor tmpFingerprint : ClusterChemistryCheck.PREFERRED_FINGERPRINTS) {
            if (ClusterChemistryCheck.getIndependence(tmpFingerprint, aParameters) == Independence.INDEPENDENT) {
                return tmpFingerprint;
            }
        }
        for (Descriptor tmpFingerprint : ClusterChemistryCheck.getBitFingerprints()) {
            if (ClusterChemistryCheck.getIndependence(tmpFingerprint, aParameters) == Independence.INDEPENDENT) {
                return tmpFingerprint;
            }
        }
        return ClusterChemistryCheck.PREFERRED_FINGERPRINTS.getFirst();
    }
    //
    /**
     * Checks the clusters of the given result chemically with the given bit fingerprint and the Bemis-Murcko scaffolds
     * and compares them with the given number of random assignments with the same cluster sizes (seeds of the result + 1,
     * + 2, ...). The fingerprints are calculated with the aromaticity model and the circular fingerprint size of the run,
     * as bit fingerprints.
     *
     * @param aRunResult clustering result
     * @param aFingerprint bit fingerprint for the Tanimoto similarities
     * @param aNumberOfRandomRuns number of random assignments, at least 1
     * @param anIsCancelledSupplier returns true if the check should be stopped
     * @param aProgressConsumer receives the number of finished steps and the number of steps
     * @return result of the check
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the descriptor is no fingerprint or the number of runs is smaller than 1
     * @throws CancellationException if the check was cancelled
     * @throws InterruptedException if the fingerprint calculation is interrupted
     */
    public static ChemistryCheckResult check(ClusteringRunResult aRunResult, Descriptor aFingerprint, int aNumberOfRandomRuns,
                                             BooleanSupplier anIsCancelledSupplier, BiConsumer<Integer, Integer> aProgressConsumer)
            throws NullPointerException, IllegalArgumentException, CancellationException, InterruptedException {
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aFingerprint, "Given fingerprint is null.");
        Objects.requireNonNull(anIsCancelledSupplier, "Given cancellation supplier is null.");
        Objects.requireNonNull(aProgressConsumer, "Given progress consumer is null.");
        if (!aFingerprint.isFingerprint()) {
            throw new IllegalArgumentException("The given descriptor is no fingerprint.");
        }
        if (aNumberOfRandomRuns < 1) {
            throw new IllegalArgumentException("The number of random runs must be at least 1.");
        }
        ClusteringParameters tmpParameters = aRunResult.parameters();
        List<MoleculeDataModel> tmpMolecules = aRunResult.molecules();
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        int tmpClusterCount = aRunResult.clusteringResult().numberOfClusters();
        //steps: fingerprints, scaffolds, the result and every random run
        int tmpNumberOfSteps = 3 + aNumberOfRandomRuns;
        aProgressConsumer.accept(0, tmpNumberOfSteps);
        //only the clustered molecules are needed
        boolean[] tmpIsClustered = new boolean[tmpMolecules.size()];
        for (int i = 0; i < tmpLabels.length; i++) {
            tmpIsClustered[i] = tmpLabels[i] >= 0;
        }
        long[][] tmpBits = ClusterChemistryCheck.calculateFingerprintBits(tmpMolecules, tmpIsClustered, aFingerprint,
                FeatureMatrixBuilder.getDescriptorCalculator(tmpParameters.circularFingerprintSize()),
                ChemUtil.createAromaticityModel(tmpParameters.electronDonationModelOption(), tmpParameters.cycleFinderOption()),
                tmpParameters.isParallelCalculation(), anIsCancelledSupplier);
        aProgressConsumer.accept(1, tmpNumberOfSteps);
        List<String> tmpScaffolds = new ArrayList<>();
        int[] tmpScaffoldIds = ClusterChemistryCheck.calculateScaffoldIds(tmpMolecules, tmpIsClustered, tmpScaffolds,
                tmpParameters.isParallelCalculation(), anIsCancelledSupplier);
        aProgressConsumer.accept(2, tmpNumberOfSteps);
        int tmpFailedCount = 0;
        for (int i = 0; i < tmpLabels.length; i++) {
            if (tmpIsClustered[i] && (tmpBits[i] == null || tmpScaffoldIds[i] < 0)) {
                tmpFailedCount++;
            }
        }
        ClusterValues tmpValues = ClusterChemistryCheck.evaluate(tmpLabels, tmpClusterCount, tmpBits, tmpScaffoldIds,
                tmpScaffolds.size(), new Random(tmpParameters.randomSeed()), tmpParameters.isParallelCalculation());
        aProgressConsumer.accept(3, tmpNumberOfSteps);
        //random baseline: the labels are shuffled among the clustered molecules, so the clusters keep their sizes
        int[] tmpClusteredRows = IntStream.range(0, tmpLabels.length).filter(i -> tmpIsClustered[i]).toArray();
        double[] tmpRandomTanimotoSums = new double[tmpClusterCount];
        double[] tmpRandomScaffoldSums = new double[tmpClusterCount];
        int[] tmpRandomTanimotoCounts = new int[tmpClusterCount];
        int[] tmpRandomScaffoldCounts = new int[tmpClusterCount];
        double[] tmpOverallRandomTanimotos = new double[aNumberOfRandomRuns];
        double[] tmpOverallRandomScaffoldShares = new double[aNumberOfRandomRuns];
        for (int r = 0; r < aNumberOfRandomRuns; r++) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The chemical check was cancelled.");
            }
            Random tmpRandom = new Random(tmpParameters.randomSeed() + r + 1);
            int[] tmpShuffledLabels = tmpLabels.clone();
            //Fisher-Yates shuffle of the labels of the clustered rows
            for (int k = tmpClusteredRows.length - 1; k > 0; k--) {
                int tmpOther = tmpRandom.nextInt(k + 1);
                int tmpLabel = tmpShuffledLabels[tmpClusteredRows[k]];
                tmpShuffledLabels[tmpClusteredRows[k]] = tmpShuffledLabels[tmpClusteredRows[tmpOther]];
                tmpShuffledLabels[tmpClusteredRows[tmpOther]] = tmpLabel;
            }
            ClusterValues tmpRandomValues = ClusterChemistryCheck.evaluate(tmpShuffledLabels, tmpClusterCount, tmpBits,
                    tmpScaffoldIds, tmpScaffolds.size(), tmpRandom, tmpParameters.isParallelCalculation());
            for (int c = 0; c < tmpClusterCount; c++) {
                if (!Double.isNaN(tmpRandomValues.meanTanimotos()[c])) {
                    tmpRandomTanimotoSums[c] += tmpRandomValues.meanTanimotos()[c];
                    tmpRandomTanimotoCounts[c]++;
                }
                if (!Double.isNaN(tmpRandomValues.scaffoldShares()[c])) {
                    tmpRandomScaffoldSums[c] += tmpRandomValues.scaffoldShares()[c];
                    tmpRandomScaffoldCounts[c]++;
                }
            }
            tmpOverallRandomTanimotos[r] = tmpRandomValues.overallMeanTanimoto();
            tmpOverallRandomScaffoldShares[r] = tmpRandomValues.overallScaffoldShare();
            aProgressConsumer.accept(4 + r, tmpNumberOfSteps);
        }
        double[] tmpRandomTanimotos = new double[tmpClusterCount];
        double[] tmpRandomScaffoldShares = new double[tmpClusterCount];
        String[] tmpMostFrequentScaffolds = new String[tmpClusterCount];
        for (int c = 0; c < tmpClusterCount; c++) {
            tmpRandomTanimotos[c] = tmpRandomTanimotoCounts[c] == 0 ? Double.NaN : tmpRandomTanimotoSums[c] / tmpRandomTanimotoCounts[c];
            tmpRandomScaffoldShares[c] = tmpRandomScaffoldCounts[c] == 0 ? Double.NaN : tmpRandomScaffoldSums[c] / tmpRandomScaffoldCounts[c];
            int tmpScaffoldId = tmpValues.mostFrequentScaffoldIds()[c];
            tmpMostFrequentScaffolds[c] = tmpScaffoldId < 0 ? null : tmpScaffolds.get(tmpScaffoldId);
        }
        return new ChemistryCheckResult(aFingerprint, tmpValues.meanTanimotos(), tmpRandomTanimotos, tmpValues.scaffoldShares(),
                tmpRandomScaffoldShares, tmpMostFrequentScaffolds, tmpValues.overallMeanTanimoto(),
                ClusterChemistryCheck.getMean(tmpOverallRandomTanimotos), ClusterChemistryCheck.getStandardDeviation(tmpOverallRandomTanimotos),
                tmpValues.overallScaffoldShare(), ClusterChemistryCheck.getMean(tmpOverallRandomScaffoldShares),
                ClusterChemistryCheck.getStandardDeviation(tmpOverallRandomScaffoldShares), aNumberOfRandomRuns, tmpFailedCount);
    }
    //
    /**
     * Returns the Tanimoto similarity of two bit fingerprints given as bits packed into long values: the number of bits
     * set in both divided by the number of bits set in at least one.
     *
     * @param aFirstBits first fingerprint
     * @param aSecondBits second fingerprint, same length
     * @return Tanimoto similarity in [0, 1]; 1 if no bit is set in both fingerprints
     */
    public static double getTanimoto(long[] aFirstBits, long[] aSecondBits) {
        int tmpIntersection = 0;
        int tmpUnion = 0;
        for (int w = 0; w < aFirstBits.length; w++) {
            tmpIntersection += Long.bitCount(aFirstBits[w] & aSecondBits[w]);
            tmpUnion += Long.bitCount(aFirstBits[w] | aSecondBits[w]);
        }
        return tmpUnion == 0 ? 1.0 : (double) tmpIntersection / tmpUnion;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns whether the given descriptor is a circular fingerprint (ECFP or FCFP).
     *
     * @param aDescriptor descriptor
     * @return true for a circular fingerprint
     */
    private static boolean isCircularFingerprint(Descriptor aDescriptor) {
        return aDescriptor.name().startsWith(ClusterChemistryCheck.CIRCULAR_FINGERPRINT_PREFIX);
    }
    //
    /**
     * Calculates the given fingerprint of the clustered molecules in batches and packs its bits into long values.
     *
     * @param aMoleculeList molecules
     * @param anIsClusteredArray true for the molecules whose fingerprint is needed
     * @param aFingerprint fingerprint
     * @param aCalculator descriptor calculator
     * @param anAromaticityModel aromaticity model applied before the calculation
     * @param anIsParallelCalculation true: molecules are processed in parallel
     * @param anIsCancelledSupplier returns true if the check should be stopped
     * @return bits of every molecule, null for molecules without fingerprint (not needed or failed)
     * @throws InterruptedException if the calculation is interrupted
     * @throws CancellationException if the check was cancelled
     */
    private static long[][] calculateFingerprintBits(List<MoleculeDataModel> aMoleculeList, boolean[] anIsClusteredArray,
                                                     Descriptor aFingerprint, DescriptorCalculator aCalculator,
                                                     Aromaticity anAromaticityModel, boolean anIsParallelCalculation,
                                                     BooleanSupplier anIsCancelledSupplier) throws InterruptedException {
        long[][] tmpBits = new long[aMoleculeList.size()][];
        int[] tmpRows = IntStream.range(0, aMoleculeList.size()).filter(i -> anIsClusteredArray[i]).toArray();
        for (int tmpStart = 0; tmpStart < tmpRows.length; tmpStart += ClusterChemistryCheck.FINGERPRINT_BATCH_SIZE) {
            if (anIsCancelledSupplier.getAsBoolean()) {
                throw new CancellationException("The chemical check was cancelled.");
            }
            int tmpEnd = Math.min(tmpRows.length, tmpStart + ClusterChemistryCheck.FINGERPRINT_BATCH_SIZE);
            List<MoleculeDataModel> tmpBatch = new ArrayList<>(tmpEnd - tmpStart);
            for (int k = tmpStart; k < tmpEnd; k++) {
                tmpBatch.add(aMoleculeList.get(tmpRows[k]));
            }
            FeatureMatrix tmpMatrix = FeatureMatrixBuilder.buildDescriptorMatrix(tmpBatch, new Descriptor[] {aFingerprint},
                    anAromaticityModel, aCalculator, anIsParallelCalculation);
            int tmpNumberOfBits = tmpMatrix.getColumnCount();
            for (int k = 0; k < tmpBatch.size(); k++) {
                long[] tmpWords = new long[(tmpNumberOfBits + 63) / 64];
                boolean tmpIsValid = true;
                for (int j = 0; j < tmpNumberOfBits && tmpIsValid; j++) {
                    float tmpValue = tmpMatrix.getValue(k, j);
                    if (Float.isNaN(tmpValue)) {
                        tmpIsValid = false;
                    } else if (tmpValue > 0.0f) {
                        tmpWords[j >> 6] |= 1L << (j & 63);
                    }
                }
                tmpBits[tmpRows[tmpStart + k]] = tmpIsValid ? tmpWords : null;
            }
        }
        return tmpBits;
    }
    //
    /**
     * Generates the Bemis-Murcko scaffold of every clustered molecule and numbers the distinct scaffolds.
     *
     * @param aMoleculeList molecules
     * @param anIsClusteredArray true for the molecules whose scaffold is needed
     * @param aScaffoldList receives the unique SMILES of the distinct scaffolds, index = scaffold id ("" = acyclic)
     * @param anIsParallelCalculation true: batches of molecules are processed in parallel
     * @param anIsCancelledSupplier returns true if the check should be stopped
     * @return scaffold id of every molecule, -1 for molecules without scaffold (not needed or failed)
     * @throws CancellationException if the check was cancelled
     */
    private static int[] calculateScaffoldIds(List<MoleculeDataModel> aMoleculeList, boolean[] anIsClusteredArray,
                                              List<String> aScaffoldList, boolean anIsParallelCalculation,
                                              BooleanSupplier anIsCancelledSupplier) {
        String[] tmpScaffolds = new String[aMoleculeList.size()];
        int tmpNumberOfBatches = (aMoleculeList.size() + ClusterChemistryCheck.SCAFFOLD_BATCH_SIZE - 1) / ClusterChemistryCheck.SCAFFOLD_BATCH_SIZE;
        IntStream tmpBatches = IntStream.range(0, tmpNumberOfBatches);
        (anIsParallelCalculation ? tmpBatches.parallel() : tmpBatches).forEach(b -> {
            //a scaffold generator instance per batch, the instances are not thread-safe
            ScaffoldGenerator tmpGenerator = new ScaffoldGenerator();
            tmpGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.SCAFFOLD);
            int tmpEnd = Math.min(aMoleculeList.size(), (b + 1) * ClusterChemistryCheck.SCAFFOLD_BATCH_SIZE);
            for (int i = b * ClusterChemistryCheck.SCAFFOLD_BATCH_SIZE; i < tmpEnd; i++) {
                if (!anIsClusteredArray[i] || anIsCancelledSupplier.getAsBoolean()) {
                    continue;
                }
                try {
                    IAtomContainer tmpScaffold = tmpGenerator.getScaffold(aMoleculeList.get(i).getAtomContainer().clone(), true);
                    tmpScaffolds[i] = tmpScaffold.getAtomCount() == 0 ? "" : ChemUtil.createUniqueSmiles(tmpScaffold, false);
                } catch (Exception anException) {
                    ClusterChemistryCheck.LOGGER.log(Level.WARNING, String.format("No scaffold could be generated for molecule %s.",
                            aMoleculeList.get(i).getName()), anException);
                }
            }
        });
        if (anIsCancelledSupplier.getAsBoolean()) {
            throw new CancellationException("The chemical check was cancelled.");
        }
        Map<String, Integer> tmpIdOfScaffoldMap = new HashMap<>();
        int[] tmpIds = new int[aMoleculeList.size()];
        for (int i = 0; i < tmpIds.length; i++) {
            if (tmpScaffolds[i] == null) {
                tmpIds[i] = -1;
                continue;
            }
            Integer tmpId = tmpIdOfScaffoldMap.get(tmpScaffolds[i]);
            if (tmpId == null) {
                tmpId = aScaffoldList.size();
                tmpIdOfScaffoldMap.put(tmpScaffolds[i], tmpId);
                aScaffoldList.add(tmpScaffolds[i]);
            }
            tmpIds[i] = tmpId;
        }
        return tmpIds;
    }
    //
    /**
     * Calculates the mean Tanimoto similarity and the share of the most frequent scaffold of every cluster of the given
     * assignment and their means over the clusters, weighted by the number of molecules with a value.
     *
     * @param aLabels cluster label of every molecule, -1 for not clustered
     * @param aClusterCount number of clusters
     * @param aBits fingerprint bits of every molecule, null if missing
     * @param aScaffoldIds scaffold id of every molecule, -1 if missing
     * @param aNumberOfScaffolds number of distinct scaffolds
     * @param aRandom random generator for the sample of large clusters
     * @param anIsParallelCalculation true: clusters are processed in parallel
     * @return values of the assignment
     */
    private static ClusterValues evaluate(int[] aLabels, int aClusterCount, long[][] aBits, int[] aScaffoldIds,
                                          int aNumberOfScaffolds, Random aRandom, boolean anIsParallelCalculation) {
        List<List<Integer>> tmpFingerprintMembers = new ArrayList<>(aClusterCount);
        List<List<Integer>> tmpScaffoldMembers = new ArrayList<>(aClusterCount);
        for (int c = 0; c < aClusterCount; c++) {
            tmpFingerprintMembers.add(new ArrayList<>());
            tmpScaffoldMembers.add(new ArrayList<>());
        }
        for (int i = 0; i < aLabels.length; i++) {
            if (aLabels[i] < 0) {
                continue;
            }
            if (aBits[i] != null) {
                tmpFingerprintMembers.get(aLabels[i]).add(i);
            }
            if (aScaffoldIds[i] >= 0) {
                tmpScaffoldMembers.get(aLabels[i]).add(i);
            }
        }
        //the samples of large clusters are drawn before the parallel part, so that the result does not depend on it
        int[][] tmpSamples = new int[aClusterCount][];
        for (int c = 0; c < aClusterCount; c++) {
            List<Integer> tmpMembers = new ArrayList<>(tmpFingerprintMembers.get(c));
            if (tmpMembers.size() > ClusterChemistryCheck.MAXIMUM_NUMBER_OF_MOLECULES_PER_CLUSTER) {
                Collections.shuffle(tmpMembers, aRandom);
                tmpMembers = tmpMembers.subList(0, ClusterChemistryCheck.MAXIMUM_NUMBER_OF_MOLECULES_PER_CLUSTER);
            }
            tmpSamples[c] = tmpMembers.stream().mapToInt(Integer::intValue).toArray();
        }
        double[] tmpTanimotos = new double[aClusterCount];
        IntStream tmpClusters = IntStream.range(0, aClusterCount);
        (anIsParallelCalculation ? tmpClusters.parallel() : tmpClusters).forEach(c -> {
            int[] tmpSample = tmpSamples[c];
            if (tmpSample.length < 2) {
                tmpTanimotos[c] = Double.NaN;
                return;
            }
            double tmpSum = 0.0;
            for (int a = 0; a < tmpSample.length; a++) {
                for (int b = a + 1; b < tmpSample.length; b++) {
                    tmpSum += ClusterChemistryCheck.getTanimoto(aBits[tmpSample[a]], aBits[tmpSample[b]]);
                }
            }
            tmpTanimotos[c] = tmpSum / ((double) tmpSample.length * (tmpSample.length - 1) / 2.0);
        });
        double[] tmpShares = new double[aClusterCount];
        int[] tmpMostFrequentIds = new int[aClusterCount];
        int[] tmpCounts = new int[aNumberOfScaffolds];
        for (int c = 0; c < aClusterCount; c++) {
            List<Integer> tmpMembers = tmpScaffoldMembers.get(c);
            int tmpBestId = -1;
            int tmpBestCount = 0;
            for (int tmpRow : tmpMembers) {
                int tmpId = aScaffoldIds[tmpRow];
                tmpCounts[tmpId]++;
                //ties are broken by the smaller id, so the result is reproducible
                if (tmpCounts[tmpId] > tmpBestCount || (tmpCounts[tmpId] == tmpBestCount && tmpId < tmpBestId)) {
                    tmpBestCount = tmpCounts[tmpId];
                    tmpBestId = tmpId;
                }
            }
            for (int tmpRow : tmpMembers) {
                tmpCounts[aScaffoldIds[tmpRow]] = 0;
            }
            tmpMostFrequentIds[c] = tmpBestId;
            tmpShares[c] = tmpMembers.size() < 2 ? Double.NaN : (double) tmpBestCount / tmpMembers.size();
        }
        double tmpTanimotoSum = 0.0;
        double tmpTanimotoWeight = 0.0;
        double tmpShareSum = 0.0;
        double tmpShareWeight = 0.0;
        for (int c = 0; c < aClusterCount; c++) {
            if (!Double.isNaN(tmpTanimotos[c])) {
                tmpTanimotoSum += tmpTanimotos[c] * tmpFingerprintMembers.get(c).size();
                tmpTanimotoWeight += tmpFingerprintMembers.get(c).size();
            }
            if (!Double.isNaN(tmpShares[c])) {
                tmpShareSum += tmpShares[c] * tmpScaffoldMembers.get(c).size();
                tmpShareWeight += tmpScaffoldMembers.get(c).size();
            }
        }
        return new ClusterValues(tmpTanimotos, tmpShares, tmpMostFrequentIds,
                tmpTanimotoWeight == 0.0 ? Double.NaN : tmpTanimotoSum / tmpTanimotoWeight,
                tmpShareWeight == 0.0 ? Double.NaN : tmpShareSum / tmpShareWeight);
    }
    //
    /**
     * Returns the mean of the given values, ignoring NaN.
     *
     * @param aValues values
     * @return mean or NaN if there is no value
     */
    private static double getMean(double[] aValues) {
        return Arrays.stream(aValues).filter(aValue -> !Double.isNaN(aValue)).average().orElse(Double.NaN);
    }
    //
    /**
     * Returns the standard deviation (population) of the given values, ignoring NaN.
     *
     * @param aValues values
     * @return standard deviation, 0 for a single value, NaN if there is no value
     */
    private static double getStandardDeviation(double[] aValues) {
        double tmpMean = ClusterChemistryCheck.getMean(aValues);
        if (Double.isNaN(tmpMean)) {
            return Double.NaN;
        }
        return Math.sqrt(Arrays.stream(aValues).filter(aValue -> !Double.isNaN(aValue))
                .map(aValue -> (aValue - tmpMean) * (aValue - tmpMean)).average().orElse(0.0));
    }
    //</editor-fold>
}
