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

import de.unijena.cheminf.mortar.model.clustering.Art2aClusterer;
import de.unijena.cheminf.mortar.model.clustering.ClusteringParameters;
import de.unijena.cheminf.mortar.model.clustering.ClusteringResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.FeatureFilter;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.FeatureTransformation;
import de.unijena.cheminf.mortar.model.clustering.FeatureType;
import de.unijena.cheminf.mortar.model.clustering.RemovedFeature;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

/**
 * Tests for {@link ClusteringDataCsvWriter}.
 *
 * @author Manuel Schauer
 */
class ClusteringDataCsvWriterTest {
    /**
     * Tests the written feature matrix including cluster labels, NaN values and escaping of a SMILES with a separator.
     *
     * @throws Exception if writing fails
     */
    @Test
    void testWriteFeatureMatrix() throws Exception {
        FeatureMatrix tmpMatrix = new FeatureMatrix(new float[][] {{1.0f, Float.NaN}, {0.5f, 2.0f}},
                new String[] {"m1", "m2"}, new String[] {"C,C", "ALOGP"},
                new FeatureType[] {FeatureType.BIT, FeatureType.CONTINUOUS}, new String[] {"F", "A_LOG_P"});
        StringWriter tmpWriter = new StringWriter();
        ClusteringDataCsvWriter.writeFeatureMatrix(tmpWriter, tmpMatrix, new int[] {0, -1}, ",");
        String[] tmpLines = tmpWriter.toString().split(System.lineSeparator());
        Assertions.assertEquals("Name,Cluster,\"C,C\",ALOGP", tmpLines[0]);
        Assertions.assertEquals("m1,0,1.0,NaN", tmpLines[1]);
        Assertions.assertEquals("m2,-1,0.5,2.0", tmpLines[2]);
    }
    //
    /**
     * Tests the written report of removed features.
     *
     * @throws Exception if writing fails
     */
    @Test
    void testWriteRemovedFeatures() throws Exception {
        StringWriter tmpWriter = new StringWriter();
        ClusteringDataCsvWriter.writeRemovedFeatures(tmpWriter, List.of(
                new RemovedFeature("A", RemovedFeature.RemovalReason.CORRELATION, "B", 0.97),
                new RemovedFeature("K", RemovedFeature.RemovalReason.CONSTANT, null, Double.NaN)), ";");
        String[] tmpLines = tmpWriter.toString().split(System.lineSeparator());
        Assertions.assertEquals("Feature;Reason;KeptPartner;Correlation", tmpLines[0]);
        Assertions.assertEquals("A;CORRELATION;B;0.97", tmpLines[1]);
        Assertions.assertEquals("K;CONSTANT;;", tmpLines[2]);
    }
    //
    /**
     * Tests that the Exporter writes all six clustering files into the given directory.
     *
     * @param aTempDirectory temporary directory provided by JUnit
     * @throws Exception if writing fails
     */
    @Test
    void testExporterWritesClusteringFiles(@TempDir Path aTempDirectory) throws Exception {
        Locale.setDefault(Locale.of("en", "GB"));
        FeatureMatrix tmpMatrix = new FeatureMatrix(new float[][] {{1.0f, 2.0f}, {0.0f, 3.0f}},
                new String[] {"m1", "m2"}, new String[] {"C", "ALOGP"},
                new FeatureType[] {FeatureType.BIT, FeatureType.CONTINUOUS}, new String[] {"F", "A_LOG_P"});
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag 1", false, 1, 0, new Descriptor[] {Descriptor.A_LOG_P},
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, true, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        ClusteringRunResult tmpRunResult = new ClusteringRunResult(tmpParameters,
                List.of(new MoleculeDataModel("C", "m1", new HashMap<>()), new MoleculeDataModel("CC", "m2", new HashMap<>())),
                tmpMatrix, tmpMatrix, tmpMatrix, List.of(),
                new ClusteringResult(new int[] {0, -1}, 1, new int[] {0}, false, true, 2, 0.5f, new float[1][2]));
        List<String> tmpFailed = new Exporter(new SettingsContainer()).exportClusteringCsvFiles(aTempDirectory.toFile(), tmpRunResult, ',');
        Assertions.assertTrue(tmpFailed.isEmpty());
        File[] tmpFiles = aTempDirectory.toFile().listFiles();
        Assertions.assertNotNull(tmpFiles);
        Assertions.assertEquals(7, tmpFiles.length);
        List<String> tmpClusterLines = Files.readAllLines(aTempDirectory.resolve("Clustering_Frag_1_clusters.csv"));
        Assertions.assertEquals("m1,C,0,true", tmpClusterLines.get(1));
        Assertions.assertEquals("m2,CC,-1,false", tmpClusterLines.get(2));
        Assertions.assertTrue(Files.readAllLines(aTempDirectory.resolve("Clustering_Frag_1_parameters.txt")).contains("descriptors=A_LOG_P"));
    }
    //
    /**
     * Tests that the cluster summaries are ranked by size (largest first, ties by cluster index), that every summary
     * starts with its representative and that the not clusterable molecules come last.
     */
    @Test
    void testClusterSummariesRankedBySize() {
        Locale.setDefault(Locale.of("en", "GB"));
        List<MoleculeDataModel> tmpMolecules = new java.util.ArrayList<>();
        for (int i = 0; i < 6; i++) {
            tmpMolecules.add(new MoleculeDataModel("C", "m" + i, new HashMap<>()));
        }
        FeatureMatrix tmpMatrix = new FeatureMatrix(new float[6][2], new String[] {"m0", "m1", "m2", "m3", "m4", "m5"},
                new String[] {"a", "b"}, new FeatureType[] {FeatureType.BIT, FeatureType.BIT}, new String[] {"F", "F"});
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", false, 1, 0, new Descriptor[0],
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        //cluster 0: m0; cluster 1: m1, m2, m3 (representative m2); cluster 2: m4; m5 not clusterable
        ClusteringRunResult tmpRunResult = new ClusteringRunResult(tmpParameters, tmpMolecules, tmpMatrix, tmpMatrix, tmpMatrix,
                List.of(), new ClusteringResult(new int[] {0, 1, 1, 1, 2, -1}, 3, new int[] {0, 2, 4}, false, true, 2, 0.5f, new float[3][2]));
        List<ClusteringRunResult.ClusterSummary> tmpRanked = tmpRunResult.getClusterSummariesRankedBySize();
        Assertions.assertEquals(List.of(1, 0, 2, -1), tmpRanked.stream().map(ClusteringRunResult.ClusterSummary::clusterIndex).toList());
        Assertions.assertEquals(3, tmpRanked.getFirst().getSize());
        Assertions.assertSame(tmpMolecules.get(2), tmpRanked.getFirst().members().getFirst());
        Assertions.assertSame(tmpMolecules.get(2), tmpRanked.getFirst().representative());
        Assertions.assertNull(tmpRanked.getLast().representative());
    }
}
