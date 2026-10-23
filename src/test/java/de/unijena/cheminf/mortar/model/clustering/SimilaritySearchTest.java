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
import de.unijena.cheminf.mortar.model.fragmentation.FragmentationService;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Tests for {@link SimilaritySearch}, {@link FeatureMatrixBuilder#buildRawFeatureVector} and the single-molecule
 * fragmentation of the {@link FragmentationService} used by the similarity search.
 *
 * @author Manuel Schauer
 */
class SimilaritySearchTest {
    /**
     * Name of the test fragmentation.
     */
    private static final String FRAGMENTATION_NAME = "TestFragmentation";
    //
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    SimilaritySearchTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Clusters five aromatic molecules (even rows) and five alcohols (odd rows) by fragments and two descriptors.
     *
     * @return clustering run result
     * @throws Exception if the pipeline fails
     */
    static ClusteringRunResult createRunResult() throws Exception {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        String[] tmpAromatics = {"c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1Cl", "c1ccccc1Br"};
        String[] tmpAlcohols = {"CCO", "CCCO", "CCCCO", "CC(C)O", "CCC(C)O"};
        for (int i = 0; i < tmpAromatics.length; i++) {
            tmpMolecules.add(SimilaritySearchTest.createMolecule(tmpAromatics[i], "aromatic" + i, Map.of("c1ccccc1", 1, "C", 1 + i % 2)));
            tmpMolecules.add(SimilaritySearchTest.createMolecule(tmpAlcohols[i], "alcohol" + i, Map.of("O", 1, "CC", 1 + i % 2)));
        }
        ClusteringService tmpService = new ClusteringService();
        tmpService.applyDescriptorPreset(ClusteringService.DescriptorPresetOption.NONE);
        tmpService.setDescriptorSelected(Descriptor.AROMATIC_ATOMS_COUNT, true);
        tmpService.setDescriptorSelected(Descriptor.ATOM_COUNT_O, true);
        return ClusteringTask.runPipeline(tmpMolecules,
                tmpService.createClusteringParameters(SimilaritySearchTest.FRAGMENTATION_NAME, false),
                tmpService.getPriorityTable(), () -> false, (aStep, aMessage) -> { });
    }
    //
    /**
     * Creates a molecule data model with the given fragment frequencies for the test fragmentation.
     *
     * @param aSmiles SMILES of the molecule
     * @param aName name of the molecule
     * @param aFrequencies fragment frequencies
     * @return molecule data model
     */
    private static MoleculeDataModel createMolecule(String aSmiles, String aName, Map<String, Integer> aFrequencies) {
        MoleculeDataModel tmpMolecule = new MoleculeDataModel(aSmiles, aName, new HashMap<>());
        tmpMolecule.getFragmentFrequencies().put(SimilaritySearchTest.FRAGMENTATION_NAME, aFrequencies);
        return tmpMolecule;
    }
    //
    /**
     * Tests the search for a clustered molecule: the molecule itself is excluded, the most similar molecules are of the
     * same kind, the similarities are sorted and the cluster is the molecule's own cluster.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testSearchByRowIndex() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        SimilaritySearch tmpSearch = new SimilaritySearch(tmpRunResult);
        SimilaritySearch.SearchResult tmpResult = tmpSearch.searchByRowIndex(0, 4);
        Assertions.assertEquals(4, tmpResult.neighbors().size());
        double tmpPreviousSimilarity = Double.POSITIVE_INFINITY;
        for (SimilaritySearch.Neighbor tmpNeighbor : tmpResult.neighbors()) {
            Assertions.assertNotEquals(0, tmpNeighbor.rowIndex());
            //all four other aromatic molecules (even rows) are more similar than any alcohol
            Assertions.assertEquals(0, tmpNeighbor.rowIndex() % 2);
            Assertions.assertTrue(tmpNeighbor.similarity() <= tmpPreviousSimilarity);
            tmpPreviousSimilarity = tmpNeighbor.similarity();
        }
        Assertions.assertEquals(tmpRunResult.clusteringResult().clusterLabels()[0], tmpResult.clusterIndex());
        Assertions.assertTrue(tmpResult.isWithinVigilance());
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpSearch.searchByRowIndex(0, 0));
    }
    //
    /**
     * Tests that the raw feature vector of a clustered molecule equals its row of the raw matrix and that searching with
     * it finds the molecule itself with similarity 1 and its cluster.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testSearchByRawFeatureVectorOfClusteredMolecule() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        float[] tmpRawVector = FeatureMatrixBuilder.buildRawFeatureVector(tmpRunResult.molecules().get(3), tmpRunResult);
        for (int j = 0; j < tmpRawVector.length; j++) {
            Assertions.assertEquals(tmpRunResult.rawMatrix().getValue(3, j), tmpRawVector[j], 1.0e-6f);
        }
        SimilaritySearch.SearchResult tmpResult = new SimilaritySearch(tmpRunResult).searchByRawFeatureVector(tmpRawVector, 3);
        Assertions.assertEquals(3, tmpResult.neighbors().getFirst().rowIndex());
        Assertions.assertEquals(1.0, tmpResult.neighbors().getFirst().similarity(), 1.0e-5);
        Assertions.assertEquals(tmpRunResult.clusteringResult().clusterLabels()[3], tmpResult.clusterIndex());
    }
    //
    /**
     * Tests that a new aromatic molecule is assigned to the cluster of the aromatic molecules and finds them first.
     *
     * @throws Exception if the pipeline fails
     */
    @Test
    void testSearchForNewMolecule() throws Exception {
        ClusteringRunResult tmpRunResult = SimilaritySearchTest.createRunResult();
        MoleculeDataModel tmpQuery = SimilaritySearchTest.createMolecule("c1ccccc1CCCC", "butylbenzene", Map.of("c1ccccc1", 1, "C", 1));
        SimilaritySearch.SearchResult tmpResult = new SimilaritySearch(tmpRunResult).searchByRawFeatureVector(
                FeatureMatrixBuilder.buildRawFeatureVector(tmpQuery, tmpRunResult), 5);
        Assertions.assertEquals(tmpRunResult.clusteringResult().clusterLabels()[0], tmpResult.clusterIndex());
        for (SimilaritySearch.Neighbor tmpNeighbor : tmpResult.neighbors()) {
            Assertions.assertEquals(0, tmpNeighbor.rowIndex() % 2);
        }
    }
    //
    /**
     * Tests that a molecule can be fragmented like a single-algorithm fragmentation of the session, and that this is not
     * possible for unknown fragmentations.
     *
     * @throws Exception if the fragmentation fails
     */
    @Test
    void testFragmentSingleMolecule() throws Exception {
        FragmentationService tmpService = new FragmentationService();
        List<MoleculeDataModel> tmpMolecules = List.of(new MoleculeDataModel("OCC1=CC=CC=C1", "benzyl alcohol", new HashMap<>()));
        tmpService.startSingleFragmentation(tmpMolecules, 1, false);
        String tmpFragmentationName = tmpService.getCurrentFragmentationName();
        Assertions.assertTrue(tmpService.canFragmentSingleMolecule(tmpFragmentationName));
        Assertions.assertFalse(tmpService.canFragmentSingleMolecule("unknown"));
        int tmpNumberOfFragmentsBefore = tmpService.getFragments().size();
        MoleculeDataModel tmpQuery = new MoleculeDataModel("OCCC1=CC=CC=C1", "query", new HashMap<>());
        tmpService.fragmentSingleMolecule(tmpQuery, tmpFragmentationName);
        Assertions.assertNotNull(tmpQuery.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName));
        Assertions.assertFalse(tmpQuery.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName).isEmpty());
        //the fragments of the data set are not changed by the query
        Assertions.assertEquals(tmpNumberOfFragmentsBefore, tmpService.getFragments().size());
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpService.fragmentSingleMolecule(tmpQuery, "unknown"));
    }
}
