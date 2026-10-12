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

import de.unijena.cheminf.mortar.model.clustering.ClusterProfiler;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.DiverseSubsetSelector;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.RemovedFeature;
import de.unijena.cheminf.mortar.model.clustering.TrainTestSplitter;
import de.unijena.cheminf.mortar.model.clustering.VigilanceScan;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Writes the data of the clustering pipeline as CSV: a feature matrix of any stage (raw, transformed, filtered =
 * clustering input), optionally with the cluster label of every molecule, and the report of the removed features.
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
        Objects.requireNonNull(aWriter, "Given writer is null.");
        Objects.requireNonNull(aFeatureMatrix, "Given feature matrix is null.");
        Objects.requireNonNull(aSeparator, "Given separator is null.");
        if (aClusterLabelArray != null && aClusterLabelArray.length != aFeatureMatrix.getRowCount()) {
            throw new IllegalArgumentException("Number of cluster labels does not match the number of rows.");
        }
        StringBuilder tmpLine = new StringBuilder(ClusteringDataCsvWriter.NAME_COLUMN_HEADER);
        if (aClusterLabelArray != null) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.CLUSTER_COLUMN_HEADER);
        }
        for (int j = 0; j < aFeatureMatrix.getColumnCount(); j++) {
            tmpLine.append(aSeparator).append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getFeatureName(j), aSeparator));
        }
        aWriter.write(tmpLine.append(System.lineSeparator()).toString());
        for (int i = 0; i < aFeatureMatrix.getRowCount(); i++) {
            tmpLine.setLength(0);
            tmpLine.append(ClusteringDataCsvWriter.escape(aFeatureMatrix.getRowName(i), aSeparator));
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
                "NumberOfUnclusteredMolecules", "MeanSimilarityToClusterVector", "IsOverflow", "IsConverged") + System.lineSeparator());
        for (VigilanceScan.ScanEntry tmpEntry : anEntryList) {
            aWriter.write(String.join(aSeparator, Float.toString(tmpEntry.vigilance()), Integer.toString(tmpEntry.numberOfClusters()),
                    Integer.toString(tmpEntry.numberOfSingletons()), Integer.toString(tmpEntry.largestClusterSize()),
                    Integer.toString(tmpEntry.numberOfUnclusteredMolecules()), Double.toString(tmpEntry.meanSimilarityToClusterVector()),
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
        aWriter.write(String.join(aSeparator, ClusteringDataCsvWriter.NAME_COLUMN_HEADER, "SMILES",
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
        aWriter.write(String.join(aSeparator, ClusteringDataCsvWriter.NAME_COLUMN_HEADER, "SMILES",
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
     * double quotes are doubled (RFC 4180). SMILES may contain characters like "," that need escaping.
     *
     * @param aText text to escape, may be null
     * @param aSeparator column separator
     * @return escaped text, empty string for null
     */
    private static String escape(String aText, String aSeparator) {
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
