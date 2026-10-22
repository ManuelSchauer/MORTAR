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

package de.unijena.cheminf.mortar.model.fragmentation.algorithm;

import de.unijena.cheminf.mortar.model.util.ChemUtil;

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Class to test the correct working of
 * {@link de.unijena.cheminf.mortar.model.fragmentation.algorithm.BiosynfoniFragmenter}.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
class BiosynfoniFragmenterTest {
    /**
     * SMILES code of the COCONUT natural product CNP0151033 used for test fragmentations.
     */
    private static final String TEST_MOLECULE_SMILES =
            "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4";
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public BiosynfoniFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests instantiation and basic settings retrieval, including verification that all default values match
     * their declared constants and that copies and restoring the defaults work.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void basicTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        Assertions.assertDoesNotThrow(tmpFragmenter::getFragmentationAlgorithmName);
        Assertions.assertDoesNotThrow(tmpFragmenter::getFragmentationAlgorithmDisplayName);
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
        }
        Assertions.assertEquals(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setInterSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        BiosynfoniFragmenter tmpCopy = (BiosynfoniFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(tmpFragmenter.getReturnedFragmentsSetting(), tmpCopy.getReturnedFragmentsSetting());
        Assertions.assertEquals(tmpFragmenter.getIntraSubstructureOverlapFilterSetting(), tmpCopy.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(tmpFragmenter.getInterSubstructureOverlapFilterSetting(), tmpCopy.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getInterSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, tmpFragmenter.getReturnedFragmentsSetting());
    }
    //
    /**
     * Tests that the stored representative structures are pairwise distinct, so that every key results in its own
     * fragment, and that every representative structure is matched by the SMARTS pattern of its key.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void representativeStructuresTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        Set<String> tmpUniqueSmilesSet = new HashSet<>();
        Set<String> tmpDisplayNameSet = new HashSet<>();
        for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
            IAtomContainer tmpRepresentative = tmpSmiPar.parseSmiles(tmpKey.getRepresentativeSmiles());
            Assertions.assertTrue(tmpUniqueSmilesSet.add(tmpSmiGen.create(tmpRepresentative)),
                    "Representative structure is not unique: " + tmpKey.getLabel());
            Assertions.assertTrue(SmartsPattern.create(tmpKey.getSmarts()).matches(tmpRepresentative),
                    "Representative structure is not matched by its own SMARTS: " + tmpKey.getLabel());
            Assertions.assertTrue(tmpDisplayNameSet.add(tmpKey.getDisplayName()), "Display name is not unique: " + tmpKey.getLabel());
            Assertions.assertFalse(tmpKey.getDisplayName().startsWith("Key '"), "Display name is missing: " + tmpKey.getLabel());
        }
    }
    //
    /**
     * Does a test fragmentation on the COCONUT natural product CNP0151033 and checks that all matches of one key are
     * returned as the same representative structure annotated with the label and SMARTS pattern of the key, and that
     * the input molecule is not altered.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void fragmentationTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles(BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES);
        String tmpInputSmilesBefore = tmpSmiGen.create(tmpMolecule);
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        Assertions.assertFalse(tmpFragmenter.shouldBeFiltered(tmpMolecule));
        Assertions.assertFalse(tmpFragmenter.shouldBePreprocessed(tmpMolecule));
        Assertions.assertTrue(tmpFragmenter.canBeFragmented(tmpMolecule));
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertFalse(tmpFragments.isEmpty());
        Assertions.assertEquals(tmpInputSmilesBefore, tmpSmiGen.create(tmpMolecule));
        Map<String, String> tmpLabelToSmilesMap = new HashMap<>();
        for (IAtomContainer tmpFragment : tmpFragments) {
            String tmpLabel = tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY);
            Assertions.assertNotNull(tmpLabel);
            Assertions.assertNotNull(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY));
            Assertions.assertEquals(tmpLabel, tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY));
            String tmpSmiles = tmpSmiGen.create(tmpFragment);
            String tmpPreviousSmiles = tmpLabelToSmilesMap.putIfAbsent(tmpLabel, tmpSmiles);
            if (tmpPreviousSmiles != null) {
                Assertions.assertEquals(tmpPreviousSmiles, tmpSmiles);
            }
        }
        //one distinct fragment per matched key
        Assertions.assertEquals(tmpLabelToSmilesMap.size(), new HashSet<>(tmpLabelToSmilesMap.values()).size());
        //a pyranose, a hydroxyl group, and a phenyl unit are present in the test molecule
        Assertions.assertTrue(tmpLabelToSmilesMap.containsKey(BiosynfoniFragmenter.BiosynfoniKey.S_PYRANOSE_C5O4.getLabel()));
        Assertions.assertTrue(tmpLabelToSmilesMap.containsKey(BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getLabel()));
        Assertions.assertTrue(tmpLabelToSmilesMap.containsKey(BiosynfoniFragmenter.BiosynfoniKey.D_PHENYL_C3.getLabel()));
    }
    //
    /**
     * Tests that the overlap filters reduce the number of returned fragments, and that the result of the
     * intra-substructure overlap filter does not depend on the input atom order.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void overlapFilterTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles(BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES);
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        int tmpUnfilteredCount = tmpFragmenter.fragmentMolecule(tmpMolecule).size();
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        Map<String, Integer> tmpIntraFilteredCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpMolecule));
        int tmpIntraFilteredCount = tmpIntraFilteredCounts.values().stream().mapToInt(Integer::intValue).sum();
        Assertions.assertTrue(tmpIntraFilteredCount < tmpUnfilteredCount);
        //same molecule with a different atom order (canonical SMILES order instead of input order)
        IAtomContainer tmpReorderedMolecule = tmpSmiPar.parseSmiles(new SmilesGenerator(SmiFlavor.Unique).create(tmpMolecule));
        Assertions.assertEquals(tmpIntraFilteredCounts, BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpReorderedMolecule)));
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(false);
        tmpFragmenter.setInterSubstructureOverlapFilterSetting(true);
        int tmpInterFilteredCount = tmpFragmenter.fragmentMolecule(tmpMolecule).size();
        Assertions.assertTrue(tmpInterFilteredCount < tmpUnfilteredCount);
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        int tmpBothFilteredCount = tmpFragmenter.fragmentMolecule(tmpMolecule).size();
        Assertions.assertTrue(tmpBothFilteredCount <= tmpIntraFilteredCount);
    }
    //
    /**
     * Tests the matched substructures option: for all overlap filter combinations, the number of fragments per key
     * equals the one of the representative structures option, every fragment has as many atoms as the SMARTS pattern
     * of its key matches, carries the label of its key, and can be converted to a valid SMILES code that MORTAR can
     * parse again. Additionally, the promiscuous ethyl pattern must produce different substructures.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void matchedSubstructuresTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles(BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES);
        String tmpInputSmilesBefore = tmpSmiGen.create(tmpMolecule);
        BiosynfoniFragmenter tmpRepresentativeFragmenter = new BiosynfoniFragmenter();
        BiosynfoniFragmenter tmpMatchedFragmenter = new BiosynfoniFragmenter();
        tmpMatchedFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        for (boolean tmpIsIntraFilterApplied : new boolean[] {false, true}) {
            for (boolean tmpIsInterFilterApplied : new boolean[] {false, true}) {
                tmpRepresentativeFragmenter.setIntraSubstructureOverlapFilterSetting(tmpIsIntraFilterApplied);
                tmpRepresentativeFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
                tmpMatchedFragmenter.setIntraSubstructureOverlapFilterSetting(tmpIsIntraFilterApplied);
                tmpMatchedFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
                List<IAtomContainer> tmpFragments = tmpMatchedFragmenter.fragmentMolecule(tmpMolecule);
                //SMARTS matching in this test prepares the target, so a clone is used to keep the input molecule untouched
                IAtomContainer tmpMoleculeClone = tmpMolecule.clone();
                Assertions.assertEquals(
                        BiosynfoniFragmenterTest.countPerLabel(tmpRepresentativeFragmenter.fragmentMolecule(tmpMolecule)),
                        BiosynfoniFragmenterTest.countPerLabel(tmpFragments));
                for (IAtomContainer tmpFragment : tmpFragments) {
                    String tmpLabel = tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY);
                    Assertions.assertEquals(BiosynfoniFragmenter.BiosynfoniKey.valueOf(BiosynfoniFragmenterTest.keyNameForLabel(tmpLabel)).getDisplayName(),
                            tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_DISPLAY_LABEL_PROPERTY_KEY));
                    SmartsPattern tmpPattern = SmartsPattern.create(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY));
                    int[] tmpAnyMatch = tmpPattern.match(tmpMoleculeClone);
                    Assertions.assertEquals(tmpAnyMatch.length, tmpFragment.getAtomCount(), tmpLabel);
                    String tmpSmiles = ChemUtil.createUniqueSmiles(tmpFragment, false);
                    Assertions.assertNotNull(tmpSmiles, tmpLabel);
                    Assertions.assertDoesNotThrow(() -> ChemUtil.parseSmilesToAtomContainer(tmpSmiles), tmpSmiles);
                }
            }
        }
        tmpMatchedFragmenter.setIntraSubstructureOverlapFilterSetting(false);
        tmpMatchedFragmenter.setInterSubstructureOverlapFilterSetting(false);
        Set<String> tmpEthylSmilesSet = new HashSet<>();
        for (IAtomContainer tmpFragment : tmpMatchedFragmenter.fragmentMolecule(tmpMolecule)) {
            if (BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel().equals(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY))) {
                tmpEthylSmilesSet.add(ChemUtil.createUniqueSmiles(tmpFragment, false));
            }
        }
        Assertions.assertTrue(tmpEthylSmilesSet.contains("CC"), tmpEthylSmilesSet.toString());
        Assertions.assertTrue(tmpEthylSmilesSet.contains("C=C"), tmpEthylSmilesSet.toString());
        Assertions.assertEquals(tmpInputSmilesBefore, tmpSmiGen.create(tmpMolecule));
    }
    //
    /**
     * Returns the name of the key enum constant with the given label.
     *
     * @param aLabel Biosynfoni label
     * @return enum constant name
     */
    private static String keyNameForLabel(String aLabel) {
        for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
            if (tmpKey.getLabel().equals(aLabel)) {
                return tmpKey.name();
            }
        }
        throw new IllegalArgumentException("Unknown label: " + aLabel);
    }    //
    /**
     * Counts the returned fragments per Biosynfoni label.
     *
     * @param aFragmentList fragments returned by the fragmenter
     * @return map of label to count
     */
    private static Map<String, Integer> countPerLabel(List<IAtomContainer> aFragmentList) {
        Map<String, Integer> tmpCounts = new HashMap<>();
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpCounts.merge(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY), 1, Integer::sum);
        }
        return tmpCounts;
    }
}
