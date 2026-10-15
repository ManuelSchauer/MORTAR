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

package de.unijena.cheminf.mortar.model.fragmentation;

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.ErtlFunctionalGroupsFinderFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.ScaffoldGeneratorFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.SugarRemovalUtilityFragmenter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Tests the fragmentation of single molecules like an existing fragmentation
 * ({@link FragmentationService#fragmentSingleMolecule}), which is used by the similarity search of the clustering for
 * query molecules that are not part of the data set. For pipelines, every molecule of a data set is fragmented alone
 * and compared with its fragments from the fragmentation of the whole data set.
 *
 * @author Manuel Schauer
 */
class FragmentationServiceSingleMoleculeTest {
    /**
     * Small data set of natural-product-like molecules (glycosides, aromatics, aliphatics) and drugs.
     */
    private static final String[] SMILES_CODES = {
            "OCC1OC(OC2=CC=C(C=C2)C=O)C(O)C(O)C1O",
            "OCC1OC(OCC2=CC=CC=C2)C(O)C(O)C1O",
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "OC1=CC=C(C=C1)CCN",
            "OCC1OC(OC2=CC(O)=C3C(=O)C=C(OC3=C2)C2=CC=C(O)C=C2)C(O)C(O)C1O",
            "CCCCCCCCCC(=O)O",
            "O=C(O)C1=CC=CC=C1O",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "OCC1OC(O)C(O)C(O)C1O"
    };
    //
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle.
     */
    FragmentationServiceSingleMoleculeTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates the test molecules.
     *
     * @return list of molecule data models
     */
    private static List<MoleculeDataModel> createMolecules() {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>();
        for (int i = 0; i < FragmentationServiceSingleMoleculeTest.SMILES_CODES.length; i++) {
            tmpMolecules.add(new MoleculeDataModel(FragmentationServiceSingleMoleculeTest.SMILES_CODES[i], "m" + i, new HashMap<>()));
        }
        return tmpMolecules;
    }
    //
    /**
     * Fragments the data set with the given pipeline, then every molecule alone like this fragmentation, and returns the
     * molecules whose fragment frequencies differ, with both results.
     *
     * @param aPipeline pipeline fragmenters
     * @return descriptions of the differing molecules, empty if all are equal
     * @throws Exception if a fragmentation fails
     */
    private static List<String> getDifferencesBetweenDataSetAndSingleMoleculeFragmentation(IMoleculeFragmenter[] aPipeline)
            throws Exception {
        FragmentationService tmpService = new FragmentationService();
        tmpService.setPipelineFragmenter(aPipeline);
        tmpService.setPipeliningFragmentationName("Test pipeline");
        List<MoleculeDataModel> tmpMolecules = FragmentationServiceSingleMoleculeTest.createMolecules();
        tmpService.startPipelineFragmentation(tmpMolecules, 2, false, false);
        String tmpFragmentationName = tmpService.getCurrentFragmentationName();
        Assertions.assertTrue(tmpService.canFragmentSingleMolecule(tmpFragmentationName));
        int tmpNumberOfFragments = tmpService.getFragments().size();
        List<String> tmpDifferences = new ArrayList<>();
        int tmpNumberOfComparedFragments = 0;
        for (MoleculeDataModel tmpMolecule : tmpMolecules) {
            MoleculeDataModel tmpSingle = new MoleculeDataModel(tmpMolecule.getUniqueSmiles(), tmpMolecule.getName(), new HashMap<>());
            tmpService.fragmentSingleMolecule(tmpSingle, tmpFragmentationName);
            Map<String, Integer> tmpExpected = tmpMolecule.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName);
            Map<String, Integer> tmpActual = tmpSingle.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName);
            Assertions.assertNotNull(tmpExpected, "molecule was not fragmented in the data set: " + tmpMolecule.getUniqueSmiles());
            tmpNumberOfComparedFragments += tmpExpected.size();
            if (!java.util.Objects.equals(tmpExpected, tmpActual)) {
                tmpDifferences.add(tmpMolecule.getUniqueSmiles() + ": data set " + tmpExpected + ", alone " + tmpActual);
            }
        }
        //the comparison must not be trivial
        Assertions.assertTrue(tmpNumberOfComparedFragments > tmpMolecules.size(), "too few fragments compared");
        //the fragmentation of single molecules does not change the fragments of the service
        Assertions.assertEquals(tmpNumberOfFragments, tmpService.getFragments().size());
        Assertions.assertEquals(tmpFragmentationName, tmpService.getCurrentFragmentationName());
        return tmpDifferences;
    }
    //
    /**
     * Sugar removal followed by functional group extraction: every molecule fragmented alone gets the same fragments
     * and frequencies as in the fragmentation of the whole data set.
     *
     * @throws Exception if a fragmentation fails
     */
    @Test
    void testSugarAndFunctionalGroupPipelineEqualsDataSetFragmentation() throws Exception {
        List<String> tmpDifferences = FragmentationServiceSingleMoleculeTest.getDifferencesBetweenDataSetAndSingleMoleculeFragmentation(
                new IMoleculeFragmenter[] {new SugarRemovalUtilityFragmenter(), new ErtlFunctionalGroupsFinderFragmenter()});
        Assertions.assertTrue(tmpDifferences.isEmpty(), String.join("\n", tmpDifferences));
    }
    //
    /**
     * Sugar removal followed by scaffold generation: every molecule fragmented alone gets the same fragments and
     * frequencies as in the fragmentation of the whole data set.
     *
     * @throws Exception if a fragmentation fails
     */
    @Test
    void testSugarAndScaffoldPipelineEqualsDataSetFragmentation() throws Exception {
        List<String> tmpDifferences = FragmentationServiceSingleMoleculeTest.getDifferencesBetweenDataSetAndSingleMoleculeFragmentation(
                new IMoleculeFragmenter[] {new SugarRemovalUtilityFragmenter(), new ScaffoldGeneratorFragmenter()});
        Assertions.assertTrue(tmpDifferences.isEmpty(), String.join("\n", tmpDifferences));
    }
}
