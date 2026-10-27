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

import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.io.ClusteringDataCsvWriter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.io.StringWriter;
import java.util.Locale;
import java.util.concurrent.CancellationException;

/**
 * Test class for {@link ClusterChemistryCheck}: Tanimoto similarity, independence of the fingerprint, the check on the
 * aromatic/alcohol data set against the random baseline and the export of its values.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
class ClusterChemistryCheckTest {
    /**
     * Tests the Tanimoto similarity of packed bit fingerprints against hand-calculated values.
     */
    @Test
    void testTanimoto() {
        //bits 0, 1, 2 and bits 1, 2, 3: intersection 2, union 4
        Assertions.assertEquals(0.5, ClusterChemistryCheck.getTanimoto(new long[] {0b0111L}, new long[] {0b1110L}), 1.0e-12);
        Assertions.assertEquals(1.0, ClusterChemistryCheck.getTanimoto(new long[] {0b0111L, 1L}, new long[] {0b0111L, 1L}), 1.0e-12);
        Assertions.assertEquals(0.0, ClusterChemistryCheck.getTanimoto(new long[] {0b0001L}, new long[] {0b0010L}), 1.0e-12);
        //bits in the second word count, too
        Assertions.assertEquals(1.0 / 3.0, ClusterChemistryCheck.getTanimoto(new long[] {1L, 1L}, new long[] {1L, 2L}), 1.0e-12);
        //two empty fingerprints are identical
        Assertions.assertEquals(1.0, ClusterChemistryCheck.getTanimoto(new long[] {0L}, new long[] {0L}), 1.0e-12);
    }
    //
    /**
     * Tests that a fingerprint used for the clustering is recognised as used, another circular fingerprint as related
     * and MACCS as independent, and that an independent fingerprint is proposed.
     */
    @Test
    void testIndependence() {
        ClusteringParameters tmpParameters = new ClusteringParameters("Frag", true, 1, 0,
                new Descriptor[] {Descriptor.TPSA, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4},
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET,
                FeatureTransformation.NONE, FeatureTransformation.NONE, 0.1, false, 0.95, FeatureFilter.CorrelationMethod.PEARSON,
                Art2aClusterer.Art2aAlgorithm.COSINE, 0.5f, 10, 10, 0.99f, 0.01f, 1.0f, 1L, 1024, false);
        Assertions.assertEquals(ClusterChemistryCheck.Independence.USED,
                ClusterChemistryCheck.getIndependence(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, tmpParameters));
        Assertions.assertEquals(ClusterChemistryCheck.Independence.RELATED,
                ClusterChemistryCheck.getIndependence(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2, tmpParameters));
        Assertions.assertEquals(ClusterChemistryCheck.Independence.INDEPENDENT,
                ClusterChemistryCheck.getIndependence(Descriptor.MACCS_FINGERPRINTER, tmpParameters));
        //ECFP4 was used, so the next preferred fingerprint is proposed
        Assertions.assertEquals(Descriptor.MACCS_FINGERPRINTER, ClusterChemistryCheck.getDefaultFingerprint(tmpParameters));
        for (Descriptor tmpFingerprint : ClusterChemistryCheck.getBitFingerprints()) {
            Assertions.assertTrue(tmpFingerprint.isFingerprint(), tmpFingerprint.name());
        }
    }
    //
    /**
     * Tests the check on the aromatic/alcohol data set: all aromatic molecules share the benzene scaffold and all
     * alcohols are acyclic, so the scaffold share and the Tanimoto similarity within the clusters are above the random
     * baseline; the check is reproducible and its values are exported with the cluster quality.
     *
     * @throws Exception if the clustering or the check fails
     */
    @Test
    void testCheckOfAromaticAndAlcoholClusters() throws Exception {
        Locale.setDefault(Locale.of("en", "GB"));
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        int[] tmpLabels = tmpRunResult.clusteringResult().clusterLabels();
        //row 0 is aromatic, row 1 an alcohol; the two groups are separated by the clustering
        int tmpAromaticCluster = tmpLabels[0];
        int tmpAlcoholCluster = tmpLabels[1];
        Assertions.assertNotEquals(tmpAromaticCluster, tmpAlcoholCluster);
        Assertions.assertEquals(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, ClusterChemistryCheck.getDefaultFingerprint(tmpRunResult.parameters()));
        ClusterChemistryCheck.ChemistryCheckResult tmpResult = ClusterChemistryCheck.check(tmpRunResult,
                Descriptor.MACCS_FINGERPRINTER, 5, () -> false, (aDone, aTotal) -> { });
        Assertions.assertEquals(0, tmpResult.numberOfFailedMolecules());
        Assertions.assertEquals("c1ccccc1", tmpResult.getMostFrequentScaffold(tmpAromaticCluster));
        Assertions.assertEquals("", tmpResult.getMostFrequentScaffold(tmpAlcoholCluster));
        Assertions.assertEquals(1.0, tmpResult.getScaffoldShare(tmpAromaticCluster), 1.0e-12);
        Assertions.assertEquals(1.0, tmpResult.getScaffoldShare(tmpAlcoholCluster), 1.0e-12);
        Assertions.assertTrue(tmpResult.getRandomScaffoldShare(tmpAromaticCluster) < 1.0);
        Assertions.assertTrue(tmpResult.overallScaffoldShare() > tmpResult.overallRandomScaffoldShare());
        Assertions.assertTrue(tmpResult.getMeanTanimoto(tmpAromaticCluster) > tmpResult.getRandomMeanTanimoto(tmpAromaticCluster));
        Assertions.assertTrue(tmpResult.overallMeanTanimoto() > tmpResult.overallRandomMeanTanimoto());
        Assertions.assertTrue(tmpResult.overallRandomMeanTanimotoDeviation() >= 0.0);
        Assertions.assertEquals(5, tmpResult.numberOfRandomRuns());
        //same seeds, same result
        ClusterChemistryCheck.ChemistryCheckResult tmpAgain = ClusterChemistryCheck.check(tmpRunResult,
                Descriptor.MACCS_FINGERPRINTER, 5, () -> false, (aDone, aTotal) -> { });
        Assertions.assertArrayEquals(tmpResult.meanTanimotos(), tmpAgain.meanTanimotos());
        Assertions.assertArrayEquals(tmpResult.randomScaffoldShares(), tmpAgain.randomScaffoldShares());
        //export: the chemistry columns follow the quality columns, the acyclic scaffold is named
        StringWriter tmpWriter = new StringWriter();
        ClusteringDataCsvWriter.writeClusterQuality(tmpWriter, new ClusterQuality(tmpRunResult), null, tmpResult, null, ",");
        String[] tmpLines = tmpWriter.toString().split(System.lineSeparator());
        Assertions.assertTrue(tmpLines[0].endsWith(",MeanTanimoto_MACCS_FINGERPRINTER,RandomMeanTanimoto_MACCS_FINGERPRINTER,"
                + "MostFrequentScaffoldShare,RandomMostFrequentScaffoldShare,MostFrequentScaffold"), tmpLines[0]);
        Assertions.assertTrue(tmpWriter.toString().contains(",acyclic" + System.lineSeparator()));
        Assertions.assertTrue(tmpWriter.toString().contains(",c1ccccc1" + System.lineSeparator()));
        //invalid input
        Assertions.assertThrows(IllegalArgumentException.class, () -> ClusterChemistryCheck.check(tmpRunResult, Descriptor.TPSA, 5,
                () -> false, (aDone, aTotal) -> { }));
        Assertions.assertThrows(IllegalArgumentException.class, () -> ClusterChemistryCheck.check(tmpRunResult,
                Descriptor.MACCS_FINGERPRINTER, 0, () -> false, (aDone, aTotal) -> { }));
        Assertions.assertThrows(CancellationException.class, () -> ClusterChemistryCheck.check(tmpRunResult,
                Descriptor.MACCS_FINGERPRINTER, 5, () -> true, (aDone, aTotal) -> { }));
    }
}
