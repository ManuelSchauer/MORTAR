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

package de.unijena.cheminf.mortar.model.io;

import de.unijena.cheminf.mortar.model.clustering.ClusterChemistryCheck;
import de.unijena.cheminf.mortar.model.clustering.ClusterProfiler;
import de.unijena.cheminf.mortar.model.clustering.ClusterQuality;
import de.unijena.cheminf.mortar.model.clustering.ClusterStability;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.DiverseSubsetSelector;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.FeatureTransformation;
import de.unijena.cheminf.mortar.model.clustering.RemovedFeature;
import de.unijena.cheminf.mortar.model.clustering.TrainTestSplitter;
import de.unijena.cheminf.mortar.model.clustering.VigilanceScan;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/**
 * Writes the data of the clustering pipeline as CSV: a feature matrix of any stage (raw, transformed, filtered =
 * clustering input), optionally with the cluster label of every molecule, the description of the feature columns
 * (group, type and transformation) and the report of the removed features.
 * Values are written with "." as decimal separator, NaN as "NaN".
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class ClusteringDataCsvWriter {
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Header of the molecule name column.
     */
    public static final String NAME_COLUMN_HEADER = "Name";
    /**
     * Header of the cluster label column.
     */
    public static final String CLUSTER_COLUMN_HEADER = "Cluster";
    /**
     * Header of the SMILES column.
     */
    public static final String SMILES_COLUMN_HEADER = "SMILES";
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private ClusteringDataCsvWriter() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Writes the given feature matrix with one line per molecule. The first column contains the molecule name, the
     * second column the cluster label if labels are given, followed by one column per feature.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aFeatureMatrix feature matrix to write
     * @param aClusterLabelArray cluster label of every row, may be null
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if the writer, matrix or separator is null
     * @throws IllegalArgumentException if the length of the label array does not match the number of rows
     * @throws IOException if writing fails
     */
    public static void writeFeatureMatrix(
            Writer aWriter,
            FeatureMatrix aFeatureMatrix,
            int[] aClusterLabelArray,
            String aSeparator)
            throws NullPointerException, IllegalArgumentException, IOException {
        ClusteringDataCsvWriter.writeFeatureMatrix(aWriter, aFeatureMatrix, aClusterLabelArray, null, aSeparator);
    }
    //
    /**
     * Writes the given rows of the given feature matrix with one line per molecule, see
     * {@link #writeFeatureMatrix(Writer, FeatureMatrix, int[], String)}, e.g. only the molecules of the selected clusters.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aFeatureMatrix feature matrix to write
     * @param aClusterLabelArray cluster label of every row of the matrix, may be null
     * @param aRowIndices indices of the rows to write in the given order; null for all rows
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if the writer, matrix or separator is null
     * @throws IllegalArgumentException if the length of the label array does not match the number of rows
     * @throws IOException if writing fails
     */
    public static void writeFeatureMatrix(
            Writer aWriter,
            FeatureMatrix aFeatureMatrix,
            int[] aClusterLabelArray,
            int[] aRowIndices,
            String aSeparator)
            throws NullPointerException, IllegalArgumentException, IOException {
        ClusteringDataCsvWriter.writeFeatureMatrix(aWriter, aFeatureMatrix, aClusterLabelArray, aRowIndices, null, aSeparator);
    }
    //
    /**
     * Writes the given rows of the given feature matrix with one line per molecule, see
     * {@link #writeFeatureMatrix(Writer, FeatureMatrix, int[], String)}; if molecules are given, their SMILES code is
     * written in the second column, after the name.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aFeatureMatrix feature matrix to write
     * @param aClusterLabelArray cluster label of every row of the matrix, may be null
     * @param aRowIndices indices of the rows to write in the given order; null for all rows
     * @param aMoleculeList molecules in the order of the matrix rows, whose SMILES codes are written; null for none
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if the writer, matrix or separator is null
     * @throws IllegalArgumentException if the length of the label array or of the molecule list does not match the
     *                                  number of rows
     * @throws IOException if writing fails
     */
    public static void writeFeatureMatrix(
            Writer aWriter,
            FeatureMatrix aFeatureMatrix,
            int[] aClusterLabelArray,
            int[] aRowIndices,
            List<MoleculeDataModel> aMoleculeList,
            String aSeparator)
            throws NullPointerException, IllegalArgumentException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        if (aClusterLabelArray != null && aClusterLabelArray.length != aFeatureMatrix.getRowCount()) {
            throw new IllegalArgumentException("Number of cluster labels does not match the number of rows.");
        }
        if (aMoleculeList != null && aMoleculeList.size() != aFeatureMatrix.getRowCount()) {
            throw new IllegalArgumentException("Number of molecules does not match the number of rows.");
        }
        StringBuilder tmpLine = new StringBuilder(ClusteringDataCsvWriter.NAME_COLUMN_HEADER);
        if (aMoleculeList != null) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.SMILES_COLUMN_HEADER);
        }
        if (aClusterLabelArray != null) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER);
        }
        for (int j = 0; j < aFeatureMatrix.getColumnCount(); j++) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getFeatureName(j), aSeparator));
        }
        aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        int tmpNumberOfRows = aRowIndices == null ? aFeatureMatrix.getRowCount() : aRowIndices.length;
        for (int r = 0; r < tmpNumberOfRows; r++) {
            int i = aRowIndices == null ? r : aRowIndices[r];
            tmpLine.setLength(0);
            tmpLine.append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getRowName(i), aSeparator));
            if (aMoleculeList != null) {
                tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.escape(aMoleculeList.get(i).getUniqueSmiles(), aSeparator));
            }
            if (aClusterLabelArray != null) {
                tmpLine.append(aSeparator).append(aClusterLabelArray[i]);
            }
            for (int j = 0; j < aFeatureMatrix.getColumnCount(); j++) {
                tmpLine.append(aSeparator).append(aFeatureMatrix.getValue(i, j));
            }
            aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the quality of the clusters with one line per cluster (see {@link ClusterQuality}): size, mean similarity to
     * the own cluster vector, mean distance within the cluster, silhouette, number of molecules with negative
     * silhouette, number k of nearest neighbours (min(10, size - 1)), neighbour purity with its random baseline, nearest
     * cluster with the mean
     * distance to it, cluster with the most neighbour
     * contacts and their share, and, if given, the stability (mean Jaccard index) of the stability check and the values
     * of the chemical check (mean Tanimoto similarity and share of the most frequent scaffold, each with its random
     * baseline, and the most frequent scaffold; the fingerprint is named in the Tanimoto headers).
     *
     * @param aWriter writer to write to (is not closed)
     * @param aQuality cluster quality
     * @param aStabilityResult result of the stability check, null if it was not run (the column is left out)
     * @param aChemistryCheckResult result of the chemical check, null if it was not run (the columns are left out)
     * @param aClusterIndices clusters to write; null or empty for all
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if the writer, quality or separator is null
     * @throws IOException if writing fails
     */
    public static void writeClusterQuality(Writer aWriter, ClusterQuality aQuality, ClusterStability.StabilityResult aStabilityResult,
                                           ClusterChemistryCheck.ChemistryCheckResult aChemistryCheckResult,
                                           Set<Integer> aClusterIndices, String aSeparator) throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aQuality, "Given cluster quality is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        List<String> tmpHeader = new ArrayList<>(List.of(ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER, "Size", "MeanSimilarityToOwnCentre",
                "MeanDistanceWithinCluster", "Silhouette", "NegativeSilhouettes", "NumberOfNeighbours", "NeighbourPurity", "RandomNeighbourPurity", "NearestCluster",
                "MeanDistanceToNearestCluster", "MostContactsCluster", "ContactShareWithMostContacts"));
        if (aStabilityResult != null) {
            //the mode tells what was changed in the runs (seeds or subsamples of the molecules)
            tmpHeader.add("StabilityJaccard_" + aStabilityResult.mode().name());
        }
        if (aChemistryCheckResult != null) {
            String tmpFingerprintName = aChemistryCheckResult.fingerprint().name();
            tmpHeader.addAll(List.of("MeanTanimoto_" + tmpFingerprintName, "RandomMeanTanimoto_" + tmpFingerprintName,
                    "MostFrequentScaffoldShare", "RandomMostFrequentScaffoldShare", "MostFrequentScaffold"));
        }
        aWriter.write(String.join(aSeparator, tmpHeader) + System.lineSeparator());
        for (ClusterQuality.ClusterQualityEntry tmpEntry : aQuality.getEntries()) {
            if (aClusterIndices != null && !aClusterIndices.isEmpty() && !aClusterIndices.contains(tmpEntry.clusterIndex())) {
                continue;
            }
            List<String> tmpValues = new ArrayList<>(List.of(Integer.toString(tmpEntry.clusterIndex()), Integer.toString(tmpEntry.size()),
                    Double.toString(tmpEntry.meanSimilarityToOwnCentre()), Double.toString(tmpEntry.meanDistanceWithinCluster()),
                    Double.toString(tmpEntry.silhouette()), Integer.toString(tmpEntry.numberOfNegativeSilhouettes()),
                    Integer.toString(tmpEntry.numberOfNeighbours()), Double.toString(tmpEntry.neighbourPurity()),
                    Double.toString(tmpEntry.randomNeighbourPurity()), Integer.toString(tmpEntry.nearestClusterIndex()),
                    Double.toString(tmpEntry.meanDistanceToNearestCluster()), Integer.toString(tmpEntry.mostContactsClusterIndex()),
                    Double.toString(tmpEntry.contactShareWithMostContacts())));
            if (aStabilityResult != null) {
                tmpValues.add(Double.toString(aStabilityResult.getClusterJaccardIndex(tmpEntry.clusterIndex())));
            }
            if (aChemistryCheckResult != null) {
                int tmpIndex = tmpEntry.clusterIndex();
                String tmpScaffold = aChemistryCheckResult.getMostFrequentScaffold(tmpIndex);
                tmpValues.addAll(List.of(Double.toString(aChemistryCheckResult.getMeanTanimoto(tmpIndex)),
                        Double.toString(aChemistryCheckResult.getRandomMeanTanimoto(tmpIndex)),
                        Double.toString(aChemistryCheckResult.getScaffoldShare(tmpIndex)),
                        Double.toString(aChemistryCheckResult.getRandomScaffoldShare(tmpIndex)),
                        //"" stands for acyclic molecules, written as such to tell it from a missing scaffold
                        tmpScaffold == null ? "" : ClusteringDataCsvWriter.escape(tmpScaffold.isEmpty() ? "acyclic" : tmpScaffold, aSeparator)));
            }
            aWriter.write(String.join(aSeparator, tmpValues) + System.lineSeparator());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes a table of descriptor values with one line per molecule: name, SMILES and one column per descriptor
     * component of the given matrix.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aMoleculeList molecules in the order of the matrix rows
     * @param aDescriptorMatrix descriptor values
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of molecules does not match the number of rows
     * @throws IOException if writing fails
     */
    public static void writeDescriptorTable(Writer aWriter, List<MoleculeDataModel> aMoleculeList, FeatureMatrix aDescriptorMatrix,
                                            String aSeparator) throws NullPointerException, IllegalArgumentException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aMoleculeList, "Given molecule list is null.");
        Objects.requireNonNull(aDescriptorMatrix, "Given descriptor matrix is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        if (aMoleculeList.size() != aDescriptorMatrix.getRowCount()) {
            throw new IllegalArgumentException("Number of molecules does not match the number of rows.");
        }
        StringBuilder tmpLine = new StringBuilder(ClusteringDataCsvWriter.NAME_COLUMN_HEADER).append(aSeparator).append(ClusteringDataCsvWriter.SMILES_COLUMN_HEADER);
        for (int j = 0; j < aDescriptorMatrix.getColumnCount(); j++) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.escape(aDescriptorMatrix.getFeatureName(j), aSeparator));
        }
        aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        for (int i = 0; i < aMoleculeList.size(); i++) {
            MoleculeDataModel tmpMolecule = aMoleculeList.get(i);
            tmpLine.setLength(0);
            tmpLine.append(ClusteringDataCsvWriter.escape(tmpMolecule.getName(), aSeparator)).append(aSeparator)
                    .append(ClusteringDataCsvWriter.escape(tmpMolecule.getUniqueSmiles(), aSeparator));
            for (int j = 0; j < aDescriptorMatrix.getColumnCount(); j++) {
                tmpLine.append(aSeparator).append(aDescriptorMatrix.getValue(i, j));
            }
            aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the description of the columns of the given feature matrix, one line per column in column order: the
     * feature name (the column header in the matrix files, e.g. the CDK component name "MW" or a fragment SMILES), the
     * group (descriptor enum constant name, e.g. "MOLECULAR_WEIGHT", or fragmentation name) and the type (BIT, COUNT or
     * CONTINUOUS), which determines the transformation of the column. If transformations are given, the column
     * "Transformation" gives the transformation applied to the column according to its type (BIT columns are never
     * transformed); if a processed matrix is given, the column "UsedForClustering" gives whether the column is part of
     * it, i.e. was not removed by the filter.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aFeatureMatrix feature matrix whose columns are described, e.g. the raw matrix
     * @param aCountTransformation transformation of the COUNT columns; null (together with the continuous
     *                             transformation) to omit the transformation column
     * @param aContinuousTransformation transformation of the CONTINUOUS columns; null (together with the count
     *                                  transformation) to omit the transformation column
     * @param aProcessedMatrix matrix the clustering was performed on; null to omit the column "UsedForClustering"
     * @param aSeparator column separator, e.g. "," or ";"
     * @throws NullPointerException if the writer, the feature matrix, the separator or only one of the transformations
     *                              is null
     * @throws IOException if writing fails
     */
    public static void writeFeatureDescription(Writer aWriter, FeatureMatrix aFeatureMatrix, FeatureTransformation aCountTransformation,
                                               FeatureTransformation aContinuousTransformation, FeatureMatrix aProcessedMatrix,
                                               String aSeparator) throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        boolean tmpIsTransformationWritten = aCountTransformation != null || aContinuousTransformation != null;
        if (tmpIsTransformationWritten) {
            Objects.requireNonNull(aCountTransformation, "Given count transformation is null.");
            Objects.requireNonNull(aContinuousTransformation, "Given continuous transformation is null.");
        }
        //columns of the processed matrix by group and name, since a name is only unique within its group
        Set<String> tmpProcessedKeys = null;
        if (aProcessedMatrix != null) {
            tmpProcessedKeys = new HashSet<>(aProcessedMatrix.getColumnCount() * 2);
            for (int j = 0; j < aProcessedMatrix.getColumnCount(); j++) {
                tmpProcessedKeys.add(aProcessedMatrix.getFeatureGroupName(j) + "\n" + aProcessedMatrix.getFeatureName(j));
            }
        }
        List<String> tmpHeader = new ArrayList<>(List.of("Feature", "Group", "Type"));
        if (tmpIsTransformationWritten) {
            tmpHeader.add("Transformation");
        }
        if (tmpProcessedKeys != null) {
            tmpHeader.add("UsedForClustering");
        }
        aWriter.write(String.join(aSeparator, tmpHeader) + System.lineSeparator());
        StringBuilder tmpLine = new StringBuilder();
        for (int j = 0; j < aFeatureMatrix.getColumnCount(); j++) {
            tmpLine.setLength(0);
            tmpLine.append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getFeatureName(j), aSeparator)).append(aSeparator)
                    .append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getFeatureGroupName(j), aSeparator)).append(aSeparator)
                    .append(aFeatureMatrix.getFeatureType(j).name());
            if (tmpIsTransformationWritten) {
                FeatureTransformation tmpTransformation = switch (aFeatureMatrix.getFeatureType(j)) {
                    case BIT -> FeatureTransformation.NONE;
                    case COUNT -> aCountTransformation;
                    case CONTINUOUS -> aContinuousTransformation;
                };
                tmpLine.append(aSeparator).append(tmpTransformation.name());
            }
            if (tmpProcessedKeys != null) {
                tmpLine.append(aSeparator).append(tmpProcessedKeys.contains(
                        aFeatureMatrix.getFeatureGroupName(j) + "\n" + aFeatureMatrix.getFeatureName(j)));
            }
            aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the report of the removed features with the columns feature, reason, kept partner feature and correlation.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aRemovedFeatureList removed features
     * @param aSeparator column separator
     * @throws NullPointerException if a parameter is null
     * @throws IOException if writing fails
     */
    public static void writeRemovedFeatures(Writer aWriter, List<RemovedFeature> aRemovedFeatureList, String aSeparator)
            throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aRemovedFeatureList, "Given removed feature list is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        aWriter.write(String.join(aSeparator, "Feature", "Reason", "KeptPartner", "Correlation") + System.lineSeparator());
        for (RemovedFeature tmpFeature : aRemovedFeatureList) {
            aWriter.write(String.join(aSeparator,
                    ClusteringDataCsvWriter.escape(tmpFeature.featureName(), aSeparator),
                    tmpFeature.reason().name(),
                    tmpFeature.partnerFeatureName() == null ? "" : ClusteringDataCsvWriter.escape(tmpFeature.partnerFeatureName(), aSeparator),
                    Double.isNaN(tmpFeature.correlation()) ? "" : Double.toString(tmpFeature.correlation()))
                    + System.lineSeparator());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the cluster profiles with the columns cluster, cluster size, category, feature, cluster value, overall value
     * and score. To keep the file small, per cluster only the given number of most enriched and most depleted fragments
     * and fingerprint bits is written, but all descriptors.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aProfileList profiles of the clusters
     * @param aNumberOfPresenceFeatures number of most enriched and of most depleted fragments and fingerprint bits per cluster
     * @param aSeparator column separator
     * @throws NullPointerException if a parameter is null
     * @throws IOException if writing fails
     */
    public static void writeClusterProfiles(Writer aWriter, List<ClusterProfiler.ClusterProfile> aProfileList,
                                            int aNumberOfPresenceFeatures, String aSeparator) throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aProfileList, "Given profile list is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        aWriter.write(String.join(aSeparator, ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER, "ClusterSize", "Category", "Feature",
                "ClusterValue", "OverallValue", "Score") + System.lineSeparator());
        for (ClusterProfiler.ClusterProfile tmpProfile : aProfileList) {
            List<ClusterProfiler.FeatureProfile> tmpFeatures = new ArrayList<>();
            for (ClusterProfiler.FeatureCategory tmpCategory : List.of(ClusterProfiler.FeatureCategory.FRAGMENT,
                    ClusterProfiler.FeatureCategory.FINGERPRINT_BIT)) {
                tmpFeatures.addAll(tmpProfile.getTopFeatures(tmpCategory, aNumberOfPresenceFeatures, true));
                tmpFeatures.addAll(tmpProfile.getTopFeatures(tmpCategory, aNumberOfPresenceFeatures, false));
            }
            tmpFeatures.addAll(tmpProfile.features().stream()
                    .filter(aFeature -> aFeature.category() == ClusterProfiler.FeatureCategory.DESCRIPTOR).toList());
            for (ClusterProfiler.FeatureProfile tmpFeature : tmpFeatures) {
                aWriter.write(String.join(aSeparator,
                        Integer.toString(tmpProfile.clusterIndex()),
                        Integer.toString(tmpProfile.size()),
                        tmpFeature.category().name(),
                        ClusteringDataCsvWriter.escape(tmpFeature.featureName(), aSeparator),
                        Double.toString(tmpFeature.clusterValue()),
                        Double.toString(tmpFeature.overallValue()),
                        Double.toString(tmpFeature.score())) + System.lineSeparator());
            }
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the result of a vigilance scan with one line per vigilance value.
     *
     * @param aWriter writer to write to (is not closed)
     * @param anEntryList entries of the scan
     * @param aSeparator column separator
     * @throws NullPointerException if a parameter is null
     * @throws IOException if writing fails
     */
    public static void writeVigilanceScan(Writer aWriter, List<VigilanceScan.ScanEntry> anEntryList, String aSeparator)
            throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(anEntryList, "Given entry list is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        aWriter.write(String.join(aSeparator, "Vigilance", "NumberOfClusters", "NumberOfSingletons", "LargestClusterSize",
                "NumberOfUnclusteredMolecules", "MeanSimilarityToClusterVector", "Silhouette", "NeighbourPurity", "IsOverflow",
                "IsConverged") + System.lineSeparator());
        for (VigilanceScan.ScanEntry tmpEntry : anEntryList) {
            aWriter.write(String.join(aSeparator, VigilanceScan.formatVigilance(tmpEntry.vigilance()), Integer.toString(tmpEntry.numberOfClusters()),
                    Integer.toString(tmpEntry.numberOfSingletons()), Integer.toString(tmpEntry.largestClusterSize()),
                    Integer.toString(tmpEntry.numberOfUnclusteredMolecules()), Double.toString(tmpEntry.meanSimilarityToClusterVector()),
                    Double.toString(tmpEntry.overallSilhouette()), Double.toString(tmpEntry.overallNeighbourPurity()),
                    Boolean.toString(tmpEntry.isOverflow()), Boolean.toString(tmpEntry.isConverged())) + System.lineSeparator());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes a diverse subset with the columns name, SMILES, cluster, rank in the cluster and similarity to the cluster
     * vector.
     *
     * @param aWriter writer to write to (is not closed)
     * @param aRunResult clustering result the subset was selected from
     * @param aSelection selected molecules
     * @param aSeparator column separator
     * @throws NullPointerException if a parameter is null
     * @throws IOException if writing fails
     */
    public static void writeDiverseSubset(Writer aWriter, ClusteringRunResult aRunResult,
                                          List<DiverseSubsetSelector.SelectedMolecule> aSelection, String aSeparator)
            throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aSelection, "Given selection is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        aWriter.write(String.join(aSeparator, ClusteringDataCsvWriter.NAME_COLUMN_HEADER, ClusteringDataCsvWriter.SMILES_COLUMN_HEADER,
                ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER, "RankInCluster", "SimilarityToClusterVector") + System.lineSeparator());
        for (DiverseSubsetSelector.SelectedMolecule tmpSelected : aSelection) {
            MoleculeDataModel tmpMolecule = aRunResult.molecules().get(tmpSelected.rowIndex());
            aWriter.write(String.join(aSeparator,
                    ClusteringDataCsvWriter.escape(tmpMolecule.getName(), aSeparator),
                    ClusteringDataCsvWriter.escape(tmpMolecule.getUniqueSmiles(), aSeparator),
                    Integer.toString(tmpSelected.clusterIndex()), Integer.toString(tmpSelected.rankInCluster()),
                    Double.toString(tmpSelected.similarityToClusterVector())) + System.lineSeparator());
        }
        aWriter.flush();
    }
    //
    /**
     * Writes the molecules of the given rows with the columns name, SMILES, cluster and set (training or test).
     *
     * @param aWriter writer to write to (is not closed)
     * @param aRunResult clustering result
     * @param aSplit training/test partitioning of the clustering result
     * @param aRowIndices rows to write, e.g. all, the training or the test indices
     * @param aSeparator column separator
     * @throws NullPointerException if a parameter is null
     * @throws IOException if writing fails
     */
    public static void writeTrainTestSplit(Writer aWriter, ClusteringRunResult aRunResult, TrainTestSplitter.TrainTestSplit aSplit,
                                           int[] aRowIndices, String aSeparator) throws NullPointerException, IOException {
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        Objects.requireNonNull(aSplit, "Given split is null.");
        Objects.requireNonNull(aRowIndices, "Given row indices are null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        aWriter.write(String.join(aSeparator, ClusteringDataCsvWriter.NAME_COLUMN_HEADER, ClusteringDataCsvWriter.SMILES_COLUMN_HEADER,
                ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER, "Set") + System.lineSeparator());
        for (int tmpRow : aRowIndices) {
            MoleculeDataModel tmpMolecule = aRunResult.molecules().get(tmpRow);
            aWriter.write(String.join(aSeparator,
                    ClusteringDataCsvWriter.escape(tmpMolecule.getName(), aSeparator),
                    ClusteringDataCsvWriter.escape(tmpMolecule.getUniqueSmiles(), aSeparator),
                    Integer.toString(aRunResult.clusteringResult().clusterLabels()[tmpRow]),
                    aSplit.isTraining(tmpRow) ? "training" : "test") + System.lineSeparator());
        }
        aWriter.flush();
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Encloses the given text in double quotes if it contains the separator, a double quote or a line break; contained
     * double quotes are doubled (RFC 4180). SMILES may contain characters like "," that need escaping. Package-private,
     * also used by the {@link Exporter}.
     *
     * @param aText text to escape, may be null
     * @param aSeparator column separator
     * @return escaped text, empty string for null
     */
    static String escape(String aText, String aSeparator) {
        if (aText == null) {
            return "";
        }
        if (aText.contains(aSeparator) || aText.contains("\"") || aText.contains("\n") || aText.contains("\r")) {
            return "\"" + aText.replace("\"", "\"\"") + "\"";
        }
        return aText;
    }
    //</editor-fold>
}
