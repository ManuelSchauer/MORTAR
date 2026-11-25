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
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

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
    /**
     * SMILES code of L-alanine.
     */
    private static final String ALANINE = "C[C@H](N)C(=O)O";
    /**
     * SMILES code of L-alanine methyl ester.
     */
    private static final String ALANINE_METHYL_ESTER = "COC(=O)[C@@H](C)N";
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
            for (MoleculeDataModel tmpParentMolecule : tmpFragment.getParentMolecules()) {
                Assertions.assertEquals(tmpExpectedDisplayName, tmpFragment.getDisplayLabelText(tmpParentMolecule),
                        tmpFragment.getUniqueSmiles());
            }
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
            for (MoleculeDataModel tmpParentMolecule : tmpFragment.getParentMolecules()) {
                Assertions.assertNull(tmpFragment.getDisplayLabelText(tmpParentMolecule), tmpFragment.getUniqueSmiles());
            }
        }
    }
    //
    /**
     * Tests that the constructor without the parameter for the display labels collects no display labels, as before
     * the parameter was added.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsNotCollectedByDefaultTest() throws Exception {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(2);
        for (String tmpSmiles : new String[] {"O=C(O)C=CC1=CC=C(O)C=C1", "CCO"}) {
            tmpMolecules.add(new MoleculeDataModel(tmpSmiles, tmpSmiles, new HashMap<>()));
        }
        Map<String, FragmentDataModel> tmpFragments = new ConcurrentHashMap<>();
        FragmentationTaskResult tmpResult = new FragmentationTask(tmpMolecules, new BiosynfoniFragmenter(), tmpFragments,
                FragmentationTaskTest.FRAGMENTATION_NAME, false).call();
        Assertions.assertEquals(0, tmpResult.exceptionsCount());
        Assertions.assertEquals(0, tmpResult.unexpectedExceptionsCount());
        Assertions.assertFalse(tmpFragments.isEmpty());
        for (FragmentDataModel tmpFragment : tmpFragments.values()) {
            Assertions.assertNull(tmpFragment.getDisplayLabelText(), tmpFragment.getUniqueSmiles());
            for (MoleculeDataModel tmpParentMolecule : tmpFragment.getParentMolecules()) {
                Assertions.assertNull(tmpFragment.getDisplayLabelText(tmpParentMolecule), tmpFragment.getUniqueSmiles());
            }
        }
    }
    //
    /**
     * Tests that the display labels of different fragmenter outputs resulting in the same fragment are combined in
     * alphabetical order, in the whole data set and per molecule. With the matched substructures option of the
     * Biosynfoni fragmenter and the original pattern set, the standard amino acid alanine and the non-standard amino
     * acid alanine methyl ester, whose ester oxygen atom belongs to the pattern, both result in the fragment alanine
     * (known limitation of this option). In each molecule, the fragment carries only the labels of the keys that
     * matched in this molecule.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelsCombinedTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
        Map<String, FragmentDataModel> tmpFragments = FragmentationTaskTest.fragment(tmpFragmenter, true,
                FragmentationTaskTest.ALANINE, FragmentationTaskTest.ALANINE_METHYL_ESTER,
                FragmentationTaskTest.ALANINE + "." + FragmentationTaskTest.ALANINE_METHYL_ESTER);
        String tmpBackbone = ChemUtil.createUniqueSmiles(new SmilesParser(SilentChemObjectBuilder.getInstance())
                .parseSmiles("CC(N)C(=O)O"), false);
        Assertions.assertTrue(tmpFragments.containsKey(tmpBackbone), tmpFragments.keySet().toString());
        String tmpStandard = BiosynfoniFragmenter.BiosynfoniKey.ALL_STD_AMINOS.getDisplayName();
        String tmpNonStandard = BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getDisplayName();
        FragmentDataModel tmpFragment = tmpFragments.get(tmpBackbone);
        Assertions.assertEquals(tmpNonStandard + ", " + tmpStandard, tmpFragment.getDisplayLabelText());
        Map<String, String> tmpExpectedLabelOfMolecule = Map.of(FragmentationTaskTest.ALANINE, tmpStandard,
                FragmentationTaskTest.ALANINE_METHYL_ESTER, tmpNonStandard,
                FragmentationTaskTest.ALANINE + "." + FragmentationTaskTest.ALANINE_METHYL_ESTER, tmpNonStandard + ", " + tmpStandard);
        Assertions.assertEquals(tmpExpectedLabelOfMolecule.size(), tmpFragment.getParentMolecules().size());
        for (MoleculeDataModel tmpParentMolecule : tmpFragment.getParentMolecules()) {
            Assertions.assertEquals(tmpExpectedLabelOfMolecule.get(tmpParentMolecule.getName()),
                    tmpFragment.getDisplayLabelText(tmpParentMolecule), tmpParentMolecule.getName());
            //the frequencies per label are only stored where the fragment has more than one label in a molecule
            if (tmpParentMolecule.getName().contains(".")) {
                Assertions.assertNotNull(tmpFragment.getDisplayLabelFrequencies(tmpParentMolecule));
            } else {
                Assertions.assertNull(tmpFragment.getDisplayLabelFrequencies(tmpParentMolecule), tmpParentMolecule.getName());
            }
        }
        //a fragment produced by only one key has the same label in all molecules
        FragmentDataModel tmpTerminalNitrogen = tmpFragments.get("*N");
        Assertions.assertNotNull(tmpTerminalNitrogen, tmpFragments.keySet().toString());
        for (MoleculeDataModel tmpParentMolecule : tmpTerminalNitrogen.getParentMolecules()) {
            Assertions.assertEquals(tmpTerminalNitrogen.getDisplayLabelText(), tmpTerminalNitrogen.getDisplayLabelText(tmpParentMolecule));
            Assertions.assertNull(tmpTerminalNitrogen.getDisplayLabelFrequencies(tmpParentMolecule));
        }
    }
    //
    /**
     * Tests that the frequencies of the display labels of a fragment in a molecule are counted per label, in the
     * order of the fragmenter output, i.e. the key order: in a mixture of two alanine molecules and one alanine methyl
     * ester molecule, the fragment alanine is produced twice by the standard amino acid key and once by the
     * non-standard amino acid key of the original pattern set, see {@link #displayLabelsCombinedTest()}.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void displayLabelFrequenciesTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
        Map<String, FragmentDataModel> tmpFragments = FragmentationTaskTest.fragment(tmpFragmenter, true,
                FragmentationTaskTest.ALANINE + "." + FragmentationTaskTest.ALANINE + "." + FragmentationTaskTest.ALANINE_METHYL_ESTER);
        String tmpBackbone = ChemUtil.createUniqueSmiles(new SmilesParser(SilentChemObjectBuilder.getInstance())
                .parseSmiles("CC(N)C(=O)O"), false);
        FragmentDataModel tmpFragment = tmpFragments.get(tmpBackbone);
        Assertions.assertNotNull(tmpFragment, tmpFragments.keySet().toString());
        MoleculeDataModel tmpMolecule = tmpFragment.getParentMolecules().iterator().next();
        Assertions.assertEquals(3, tmpMolecule.getFragmentFrequencyOfSpecificFragmentation(
                FragmentationTaskTest.FRAGMENTATION_NAME).get(tmpBackbone));
        Assertions.assertEquals(List.of(
                        Map.entry(BiosynfoniFragmenter.BiosynfoniKey.ALL_STD_AMINOS.getDisplayName(), 2),
                        Map.entry(BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getDisplayName(), 1)),
                List.copyOf(tmpFragment.getDisplayLabelFrequencies(tmpMolecule).entrySet()));
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
