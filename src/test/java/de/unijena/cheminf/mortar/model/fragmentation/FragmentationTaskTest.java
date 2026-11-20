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

import de.unijena.cheminf.mortar.model.data.FragmentDataModel;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.BiosynfoniFragmenter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Class to test the correct working of {@link FragmentationTask}, currently the collection of the fragment display
 * labels.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
class FragmentationTaskTest {
    /**
     * Name of the fragmentation used in the tests.
     */
    private static final String FRAGMENTATION_NAME = "Test fragmentation";
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the display names are imported from the message.properties file.
     */
    public FragmentationTaskTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests that the display labels set by the fragmenter are collected if requested, as done for single
     * fragmentations, so that every Biosynfoni fragment is labelled with the display name of its key.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsCollectedTest() throws Exception {
        Map<String, FragmentDataModel> tmpFragments = FragmentationTaskTest.fragment(new BiosynfoniFragmenter(), true,
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4", "CCO");
        Assertions.assertFalse(tmpFragments.isEmpty());
        for (FragmentDataModel tmpFragment : tmpFragments.values()) {
            String tmpLabel = (String) tmpFragment.getProperties().get(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY);
            String tmpExpectedDisplayName = null;
            for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
                if (tmpKey.getLabel().equals(tmpLabel)) {
                    tmpExpectedDisplayName = tmpKey.getDisplayName();
                }
            }
            Assertions.assertNotNull(tmpExpectedDisplayName, tmpFragment.getUniqueSmiles());
            Assertions.assertEquals(tmpExpectedDisplayName, tmpFragment.getDisplayLabelText(), tmpFragment.getUniqueSmiles());
        }
    }
    //
    /**
     * Tests that no display labels are collected if not requested, as done for pipeline fragmentations.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsNotCollectedTest() throws Exception {
        Map<String, FragmentDataModel> tmpFragments = FragmentationTaskTest.fragment(new BiosynfoniFragmenter(), false,
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4", "CCO");
        Assertions.assertFalse(tmpFragments.isEmpty());
        for (FragmentDataModel tmpFragment : tmpFragments.values()) {
            Assertions.assertNull(tmpFragment.getDisplayLabelText(), tmpFragment.getUniqueSmiles());
        }
    }
    //
    /**
     * Tests that the display labels of different fragmenter outputs resulting in the same fragment are combined in
     * alphabetical order. With the matched substructures option of the Biosynfoni fragmenter, the ether oxygen atom of
     * dimethyl ether and the hydroxyl oxygen atom of ethanol both result in the fragment "O" (known limitation of
     * this option).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsCombinedTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Map<String, FragmentDataModel> tmpFragments = FragmentationTaskTest.fragment(tmpFragmenter, true, "COC", "CCO");
        Assertions.assertTrue(tmpFragments.containsKey("O"), tmpFragments.keySet().toString());
        Assertions.assertEquals(
                BiosynfoniFragmenter.BiosynfoniKey.O_ETHER.getDisplayName() + ", "
                        + BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getDisplayName(),
                tmpFragments.get("O").getDisplayLabelText());
    }
    //
    /**
     * Fragments the given molecules with a single fragmentation task.
     *
     * @param aFragmenter fragmenter to use
     * @param areDisplayLabelsCollected whether the display labels of the fragments should be collected
     * @param aSmilesArray SMILES codes of the molecules to fragment
     * @return map of unique SMILES codes to the resulting fragments
     * @throws Exception if anything goes wrong
     */
    private static Map<String, FragmentDataModel> fragment(BiosynfoniFragmenter aFragmenter, boolean areDisplayLabelsCollected,
                                                           String... aSmilesArray) throws Exception {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(aSmilesArray.length);
        for (String tmpSmiles : aSmilesArray) {
            tmpMolecules.add(new MoleculeDataModel(tmpSmiles, tmpSmiles, new HashMap<>()));
        }
        Map<String, FragmentDataModel> tmpFragments = new ConcurrentHashMap<>();
        FragmentationTaskResult tmpResult = new FragmentationTask(tmpMolecules, aFragmenter, tmpFragments,
                FragmentationTaskTest.FRAGMENTATION_NAME, false, areDisplayLabelsCollected).call();
        Assertions.assertEquals(0, tmpResult.exceptionsCount());
        Assertions.assertEquals(0, tmpResult.unexpectedExceptionsCount());
        return tmpFragments;
    }
}
