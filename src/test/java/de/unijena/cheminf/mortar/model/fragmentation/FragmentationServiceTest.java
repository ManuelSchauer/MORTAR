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
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.SugarRemovalUtilityFragmenter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Class to test the correct working of {@link FragmentationService}, currently the collection of the fragment display
 * labels in single fragmentations.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
class FragmentationServiceTest {
    /**
     * SMILES code of the COCONUT natural product CNP0151033, a glycoside, used for test fragmentations.
     */
    private static final String TEST_MOLECULE_SMILES =
            "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4";
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenters because the display names are imported from the message.properties file.
     */
    public FragmentationServiceTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests that only the Biosynfoni fragmenter declares that it provides display labels.
     */
    @Test
    void providesDisplayLabelsTest() {
        FragmentationService tmpService = new FragmentationService();
        for (IMoleculeFragmenter tmpFragmenter : tmpService.getFragmenters()) {
            Assertions.assertEquals(tmpFragmenter instanceof BiosynfoniFragmenter, tmpFragmenter.providesDisplayLabels(),
                    tmpFragmenter.getFragmentationAlgorithmName());
        }
    }
    //
    /**
     * Tests that a single fragmentation with the Biosynfoni fragmenter collects the display labels of the fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void biosynfoniDisplayLabelsCollectedTest() throws Exception {
        Map<String, FragmentDataModel> tmpFragments = FragmentationServiceTest.fragment(
                new BiosynfoniFragmenter().getFragmentationAlgorithmDisplayName(), new HashMap<>());
        Assertions.assertFalse(tmpFragments.isEmpty());
        for (FragmentDataModel tmpFragment : tmpFragments.values()) {
            Assertions.assertNotNull(tmpFragment.getDisplayLabelText(), tmpFragment.getUniqueSmiles());
        }
    }
    //
    /**
     * Tests that a single fragmentation with a fragmenter that does not provide display labels does not take over a
     * display label from the input molecule, e.g. after re-importing exported Biosynfoni fragments. The Sugar Removal
     * Utility fragmenter returns the aglycone as a copy of the input molecule including its properties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsOfInputMoleculesIgnoredTest() throws Exception {
        Map<Object, Object> tmpProperties = new HashMap<>();
        tmpProperties.put(IMoleculeFragmenter.FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY,
                BiosynfoniFragmenter.BiosynfoniKey.O_ETHER.getDisplayName());
        Map<String, FragmentDataModel> tmpFragments = FragmentationServiceTest.fragment(
                new SugarRemovalUtilityFragmenter().getFragmentationAlgorithmDisplayName(), tmpProperties);
        Assertions.assertFalse(tmpFragments.isEmpty());
        //precondition: the property is passed on to the aglycone fragment, so the label would be collected otherwise
        Assertions.assertTrue(tmpFragments.values().stream().anyMatch(aFragment ->
                aFragment.getProperties().containsKey(IMoleculeFragmenter.FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY)));
        for (FragmentDataModel tmpFragment : tmpFragments.values()) {
            Assertions.assertNull(tmpFragment.getDisplayLabelText(), tmpFragment.getUniqueSmiles());
        }
    }
    //
    /**
     * Fragments the test molecule in a single fragmentation with the fragmenter of the given display name.
     *
     * @param aFragmenterDisplayName display name of the fragmenter to select
     * @param aPropertyMap properties of the test molecule
     * @return map of unique SMILES codes to the resulting fragments
     * @throws Exception if anything goes wrong
     */
    private static Map<String, FragmentDataModel> fragment(String aFragmenterDisplayName, Map<Object, Object> aPropertyMap)
            throws Exception {
        FragmentationService tmpService = new FragmentationService();
        tmpService.setSelectedFragmenter(aFragmenterDisplayName);
        Assertions.assertEquals(aFragmenterDisplayName, tmpService.getSelectedFragmenter().getFragmentationAlgorithmDisplayName());
        MoleculeDataModel tmpMolecule = new MoleculeDataModel(FragmentationServiceTest.TEST_MOLECULE_SMILES, "CNP0151033", aPropertyMap);
        tmpService.startSingleFragmentation(List.of(tmpMolecule), 1, false);
        return tmpService.getFragments();
    }
}
