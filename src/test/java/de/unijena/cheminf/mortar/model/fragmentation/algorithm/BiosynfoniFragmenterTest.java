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

import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.preference.PreferenceContainer;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;

import javafx.beans.property.Property;
import javafx.beans.value.ObservableBooleanValue;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
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
    /**
     * Seed of the random atom and bond permutations used to test the independence of the input atom order.
     */
    private static final long PERMUTATION_SEED = 42L;
    /**
     * Number of random atom and bond permutations per molecule and configuration used to test the independence of the
     * input atom order.
     */
    private static final int NUMBER_OF_PERMUTATIONS = 5;
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
        Assertions.assertEquals(BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT, tmpFragmenter.getPatternSetSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setInterSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            Assertions.assertEquals(BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT, tmpFragmenter.getVariantSetting(tmpVariant));
        }
        tmpFragmenter.setVariantSetting(BiosynfoniFragmenter.BiosynfoniVariant.ETHER, !BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT);
        BiosynfoniFragmenter tmpCopy = (BiosynfoniFragmenter) tmpFragmenter.copy();
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            Assertions.assertEquals(tmpFragmenter.getVariantSetting(tmpVariant), tmpCopy.getVariantSetting(tmpVariant));
        }
        Assertions.assertEquals(tmpFragmenter.getPatternSetSetting(), tmpCopy.getPatternSetSetting());
        Assertions.assertEquals(tmpFragmenter.getReturnedFragmentsSetting(), tmpCopy.getReturnedFragmentsSetting());
        Assertions.assertEquals(tmpFragmenter.getIntraSubstructureOverlapFilterSetting(), tmpCopy.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(tmpFragmenter.getInterSubstructureOverlapFilterSetting(), tmpCopy.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT, tmpFragmenter.getPatternSetSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT,
                tmpFragmenter.getVariantSetting(BiosynfoniFragmenter.BiosynfoniVariant.ETHER));
        Assertions.assertEquals(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getInterSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, tmpFragmenter.getReturnedFragmentsSetting());
    }
    //
    /**
     * Tests that the stored representative structures are pairwise distinct, so that every key results in its own
     * fragment, and that every representative structure is matched by the SMARTS pattern of its key. This is checked
     * for no variants, all variants, every single variant, and all variants but one, which covers all combinations
     * of the two variants of the phenyl unit and hydroxyl keys.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void representativeStructuresTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        for (Set<BiosynfoniFragmenter.BiosynfoniVariant> tmpVariants : BiosynfoniFragmenterTest.getTestVariantSets()) {
            Set<String> tmpUniqueSmilesSet = new HashSet<>();
            Set<String> tmpDisplayNameSet = new HashSet<>();
            for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
                String tmpMessage = tmpKey.getLabel() + ", " + tmpVariants;
                IAtomContainer tmpRepresentative = tmpSmiPar.parseSmiles(tmpKey.getRepresentativeSmiles(tmpVariants));
                Assertions.assertTrue(tmpUniqueSmilesSet.add(tmpSmiGen.create(tmpRepresentative)),
                        "Representative structure is not unique: " + tmpMessage);
                Assertions.assertTrue(SmartsPattern.create(tmpKey.getSmarts(tmpVariants)).matches(tmpRepresentative),
                        "Representative structure is not matched by its own SMARTS: " + tmpMessage);
                Assertions.assertTrue(tmpDisplayNameSet.add(tmpKey.getDisplayName()), "Display name is not unique: " + tmpMessage);
                Assertions.assertFalse(tmpKey.getDisplayName().startsWith("Key '"), "Display name is missing: " + tmpMessage);
            }
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
     * Tests that the intra-substructure overlap filter selects the matches along chains and rings, so that they are
     * packed as densely as in the reference implementation. The expected counts are taken from the reference
     * implementation, where they do not depend on the input atom order. Each molecule is given in two different atom
     * orders to check that the result does not depend on it.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void intraSubstructureOverlapFilterDensePackingTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        //p-benzoquinone, three disjoint ethyl units in the six-membered ring
        for (String tmpSmiles : new String[] {"O=C1C=CC(=O)C=C1", "C1=CC(=O)C=CC1=O"}) {
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)));
            Assertions.assertEquals(3, tmpCounts.get(BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel()), tmpSmiles);
            Assertions.assertEquals(2, tmpCounts.get(BiosynfoniFragmenter.BiosynfoniKey.D2_ACETYL.getLabel()), tmpSmiles);
        }
        //thiophene 1,1-dioxide, two disjoint ethyl units in the chain of four carbon atoms
        for (String tmpSmiles : new String[] {"O=S1(=O)C=CC=C1", "C1=CS(=O)(=O)C=C1"}) {
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)));
            Assertions.assertEquals(2, tmpCounts.get(BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel()), tmpSmiles);
        }
    }
    //
    /**
     * Tests that nitro groups, nitrate esters, and N-oxides give the same fragments in their pentavalent and
     * charge-separated notation for all overlap filter combinations, and that their oxygen atoms are not matched as
     * hydroxyl group, while hydroxylamines (also deprotonated), oximes, hydroxamic acids, and deprotonated hydroxyl
     * groups still are.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void nitroGroupAndNOxideNotationTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String tmpHydroxylLabel = BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getLabel();
        //pentavalent and charge-separated notation of nitromethane, nitrobenzene, methyl nitrate, 4-nitrobenzyl alcohol,
        // trimethylamine N-oxide, and pyridine N-oxide
        String[][] tmpSmilesPairs = {
                {"CN(=O)=O", "C[N+](=O)[O-]"},
                {"c1ccccc1N(=O)=O", "c1ccccc1[N+](=O)[O-]"},
                {"CON(=O)=O", "CO[N+](=O)[O-]"},
                {"OCc1ccc(cc1)N(=O)=O", "OCc1ccc(cc1)[N+](=O)[O-]"},
                {"CN(C)(C)=O", "C[N+](C)(C)[O-]"},
                {"O=n1ccccc1", "[O-][n+]1ccccc1"}};
        //only the alcohol group of 4-nitrobenzyl alcohol is a hydroxyl group
        int[] tmpExpectedHydroxylCounts = {0, 0, 0, 1, 0, 0};
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        for (boolean tmpIsIntraFilterApplied : new boolean[] {false, true}) {
            for (boolean tmpIsInterFilterApplied : new boolean[] {false, true}) {
                tmpFragmenter.setIntraSubstructureOverlapFilterSetting(tmpIsIntraFilterApplied);
                tmpFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
                for (int i = 0; i < tmpSmilesPairs.length; i++) {
                    Map<String, Integer> tmpPentavalentCounts = BiosynfoniFragmenterTest.countPerLabel(
                            tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmilesPairs[i][0])));
                    Map<String, Integer> tmpChargeSeparatedCounts = BiosynfoniFragmenterTest.countPerLabel(
                            tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmilesPairs[i][1])));
                    Assertions.assertEquals(tmpPentavalentCounts, tmpChargeSeparatedCounts, tmpSmilesPairs[i][1]);
                    Assertions.assertEquals(tmpExpectedHydroxylCounts[i], tmpChargeSeparatedCounts.getOrDefault(tmpHydroxylLabel, 0),
                            tmpSmilesPairs[i][1]);
                }
            }
        }
        //N-methylhydroxylamine, deprotonated N,N-dimethylhydroxylamine, acetaldoxime, acetohydroxamic acid, and phenolate
        tmpFragmenter.restoreDefaultSettings();
        for (String tmpSmiles : new String[] {"CNO", "CN(C)[O-]", "CC=NO", "CC(=O)NO", "[O-]c1ccccc1"}) {
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)));
            Assertions.assertEquals(1, tmpCounts.getOrDefault(tmpHydroxylLabel, 0), tmpSmiles);
        }
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
     * Tests that the side chain carbon atoms of the phenyl-C3 and phenyl-C2 units must not be part of any six-membered
     * ring, regardless of the bond orders in this ring and of further fused rings like the lactone ring of
     * aryltetralin lignans, while open side chains are still matched.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void phenylSideChainRingTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        String tmpPhenylC3Label = BiosynfoniFragmenter.BiosynfoniKey.D_PHENYL_C3.getLabel();
        String tmpPhenylC2Label = BiosynfoniFragmenter.BiosynfoniKey.D_PHENYL_C2.getLabel();
        //side chain in a six-membered ring: phenylcyclohexane, 4-phenylcyclohexene (non-aromatic double bond in the
        // ring), aryltetralin lignan skeletons with and without lactone ring, and aryldihydronaphthalene lignan
        // skeletons (double bond in the ring) with and without lactone ring
        for (String tmpSmiles : new String[] {"c1ccccc1C1CCCCC1", "c1ccccc1C1CC=CCC1",
                "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "OCC1Cc2cc(OC)c(O)cc2C(c2ccc(O)c(OC)c2)C1CO",
                "COc1cc(cc(OC)c1OC)C1=C2C(COC2=O)Cc2cc3OCOc3cc12", "OCC1=Cc2cc(OC)c(O)cc2C(c2ccc(O)c(OC)c2)C1CO"}) {
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)));
            Assertions.assertEquals(0, tmpCounts.getOrDefault(tmpPhenylC3Label, 0), tmpSmiles);
            Assertions.assertEquals(0, tmpCounts.getOrDefault(tmpPhenylC2Label, 0), tmpSmiles);
        }
        //open side chains: propylbenzene, cinnamyl alcohol, and cannabidiol, whose only phenyl-C3 unit is in the pentyl
        // chain, while its aryl group is attached to a cyclohexene ring
        for (String tmpSmiles : new String[] {"CCCc1ccccc1", "OC/C=C/c1ccccc1", "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"}) {
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)));
            Assertions.assertEquals(1, tmpCounts.getOrDefault(tmpPhenylC3Label, 0), tmpSmiles);
        }
    }
    //
    /**
     * Tests that completely matched aromatic rings of fused ring systems remain aromatic in the matched substructures
     * mode: the pyrrole ring of indole, the pyridine ring of quinoline, and both rings of naphthalene.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void matchedSubstructuresOfFusedAromaticRingsTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        //molecule, key of the ring, expected fragment of the ring
        String[][] tmpCases = {
                {"c1ccc2[nH]ccc2c1", BiosynfoniFragmenter.BiosynfoniKey.D_C4N.getLabel(), "c1cc[nH]c1"},
                {"c1ccc2ncccc2c1", BiosynfoniFragmenter.BiosynfoniKey.D_C5N.getLabel(), "c1ccncc1"},
                {"c1ccc2ccccc2c1", BiosynfoniFragmenter.BiosynfoniKey.R_C6.getLabel(), "c1ccccc1"}};
        for (String[] tmpCase : tmpCases) {
            String tmpExpectedSmiles = tmpSmiGen.create(tmpSmiPar.parseSmiles(tmpCase[2]));
            int tmpRingCount = 0;
            for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpCase[0]))) {
                if (tmpCase[1].equals(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY))) {
                    Assertions.assertEquals(tmpExpectedSmiles, tmpSmiGen.create(tmpFragment), tmpCase[0]);
                    tmpRingCount++;
                }
            }
            Assertions.assertTrue(tmpRingCount > 0, tmpCase[0]);
        }
    }
    //
    /**
     * Tests that a molecule without Kekulé structure, which MORTAR reads with unset bond orders if necessary, is also
     * fragmented in the matched substructures mode, with the same number of fragments per key as in the
     * representative structures mode.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void matchedSubstructuresOfNonKekulizableMoleculeTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpSmiPar.kekulise(false);
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("OCCc1cccc1");
        BiosynfoniFragmenter tmpRepresentativeFragmenter = new BiosynfoniFragmenter();
        BiosynfoniFragmenter tmpMatchedFragmenter = new BiosynfoniFragmenter();
        tmpMatchedFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        List<IAtomContainer> tmpFragments = tmpMatchedFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertFalse(tmpFragments.isEmpty());
        Assertions.assertEquals(BiosynfoniFragmenterTest.countPerLabel(tmpRepresentativeFragmenter.fragmentMolecule(tmpMolecule)),
                BiosynfoniFragmenterTest.countPerLabel(tmpFragments));
        for (IAtomContainer tmpFragment : tmpFragments) {
            Assertions.assertNotNull(ChemUtil.createUniqueSmiles(tmpFragment, false),
                    (String) tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY));
        }
    }
    //
    /**
     * Tests every variant of the variant pattern set against the original pattern set: for each case, the number of
     * matches of the given key in the given molecule is checked in both pattern sets, without overlap filters. The
     * cases cover the patterns that match more than their name describes (e.g. open-chain sugar patterns matching sugar
     * rings, the phenethylamine pattern matching indoles, ring patterns matching the perimeter of fused ring systems)
     * and the interpretations (e.g. aromatic phenyl rings, methylmalonyl units at branching points), as well as
     * unchanged behaviour.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void patternSetTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //key, SMILES code, number of matches with the original patterns, number of matches with the pattern variants
        String[][] tmpCases = {
                //NADH must not match NADPH
                {"CO_NADH", "NC(=O)C1=CN(C=CC1)C1OC(COP(O)(=O)OP(O)(=O)OCC2OC(C(O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", "1", "1"},
                {"CO_NADH", "NC(=O)C1=CN(C=CC1)C1OC(COP(O)(=O)OP(O)(=O)OCC2OC(C(OP(O)(O)=O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", "1", "0"},
                //standard amino acids, also N-methylated and esterified
                {"ALL_STD_AMINOS", "NCC(=O)NC(C)C(=O)O", "2", "2"},
                {"ALL_STD_AMINOS", "COC(=O)C(C)N", "0", "1"},
                {"ALL_STD_AMINOS", "CC(=O)N(C)C(C)C(N)=O", "0", "1"},
                //non-standard amino acids: ornithine, 2-aminobutyric acid, 2-aminoisobutyric acid, alanine methyl ester
                {"NON_STD_AMINOS", "NCCCC(N)C(=O)O", "0", "1"},
                {"NON_STD_AMINOS", "CCC(N)C(=O)O", "0", "1"},
                {"NON_STD_AMINOS", "CC(C)(N)C(=O)O", "0", "1"},
                {"NON_STD_AMINOS", "COC(=O)C(C)N", "1", "0"},
                {"NON_STD_AMINOS", "CC(N)C(=O)O", "0", "0"},
                //open-chain sugars: sorbitol and ribitol, but not glucopyranose and ribofuranose
                {"S_OPENPYR_C6O6", "OCC(O)C(O)C(O)C(O)CO", "1", "1"},
                {"S_OPENPYR_C6O6", "OCC1OC(O)C(O)C(O)C1O", "1", "0"},
                {"S_OPENFUR_C5O5", "OCC(O)C(O)C(O)CO", "1", "1"},
                {"S_OPENFUR_C5O5", "OCC1OC(O)C(O)C1O", "1", "0"},
                //phenethylamine units: phenethylamine and isoquinoline, but not indole and tryptamine
                {"D_PHENYL_C2N", "NCCc1ccccc1", "1", "1"},
                {"D_PHENYL_C2N", "c1ccc2cnccc2c1", "1", "1"},
                {"D_PHENYL_C2N", "c1ccc2[nH]ccc2c1", "1", "0"},
                {"D_PHENYL_C2N", "NCCc1c[nH]c2ccccc12", "1", "0"},
                //phenyl units: aromatic rings only, side chain in no six-membered ring
                {"D_PHENYL_C3", "CCCc1ccccc1", "1", "1"},
                {"D_PHENYL_C3", "CCCC1CCCCC1", "1", "0"},
                {"D_PHENYL_C3", "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "2", "0"},
                {"D_PHENYL_C2", "CC(C)C1CCC(C)CC1O", "2", "0"},
                {"D_PHENYL_C1", "Cc1ccccc1", "1", "1"},
                {"D_PHENYL_C1", "CC(C)C1CCC(C)CC1O", "2", "0"},
                //methylmalonyl units at branching points, not at chain ends, aromatic rings, or methyl ketones
                {"D2_METHYLMALONYL", "CCC(C)C(=O)O", "3", "2"},
                {"D2_METHYLMALONYL", "CCCCCCCCCCCCCCCC(=O)O", "1", "0"},
                {"D2_METHYLMALONYL", "Cc1ccccc1", "2", "0"},
                {"D2_METHYLMALONYL", "CC(C)=O", "1", "0"},
                //one match per phosphorus and sulfur atom
                {"PHOSPHATE", "OP(O)(=O)OP(O)(O)=O", "8", "2"},
                {"SULFONATE", "COS(=O)(=O)O", "4", "1"},
                //terminal nitrogen atoms with a single bond only
                {"N_NITRATE", "CN", "1", "1"},
                {"N_NITRATE", "CC#N", "1", "0"},
                //no peroxides as ethers
                {"O_ETHER", "CCOCC", "1", "1"},
                {"O_ETHER", "COOC", "2", "0"},
                //hemiacetal hydroxyl groups and nitro groups
                {"O_HYDROXYL", "OCC1OC(O)C(O)C(O)C1O", "4", "5"},
                {"O_HYDROXYL", "CC(=O)O", "0", "0"},
                {"O_HYDROXYL", "C[N+](=O)[O-]", "1", "0"},
                //rings in the chemical sense only, not the perimeter of fused or bridged ring systems
                {"R_C6", "c1ccc2ccccc2c1", "2", "2"},
                {"R_C10", "c1ccc2ccccc2c1", "1", "0"},
                {"R_C10", "C1CCC2CCCCC2C1", "1", "0"},
                {"R_C10", "C1CCCCCCCCC1", "1", "1"},
                {"R_C9", "c1ccc2CCCc2c1", "1", "0"},
                {"R_C6", "C1CC2CC2C1", "1", "0"},
                {"R_C5", "C1CC2CC2C1", "1", "1"},
                {"R_C3", "C1CC2CC2C1", "1", "1"},
                //unchanged keys
                {"S_PYRANOSE_C5O4", "OCC1OC(O)C(O)C(O)C1O", "1", "1"},
                {"D_ISOPRENE", "C=CC(C)=C", "1", "1"}};
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        for (String[] tmpCase : tmpCases) {
            String tmpLabel = BiosynfoniFragmenter.BiosynfoniKey.valueOf(tmpCase[0]).getLabel();
            tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
            Map<String, Integer> tmpOriginalCounts = BiosynfoniFragmenterTest.countPerLabel(
                    tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpCase[1])));
            tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT);
            Map<String, Integer> tmpVariantCounts = BiosynfoniFragmenterTest.countPerLabel(
                    tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpCase[1])));
            Assertions.assertEquals(Integer.parseInt(tmpCase[2]), tmpOriginalCounts.getOrDefault(tmpLabel, 0),
                    tmpCase[0] + ", original: " + tmpCase[1]);
            Assertions.assertEquals(Integer.parseInt(tmpCase[3]), tmpVariantCounts.getOrDefault(tmpLabel, 0),
                    tmpCase[0] + ", variant: " + tmpCase[1]);
        }
    }
    //
    /**
     * Tests that the pattern set setting changes the SMARTS pattern and the representative structure that are attached
     * to the returned fragments, for the two keys whose representative structures differ between the pattern sets, and
     * that the original pattern set uses the original patterns of the reference implementation for the phenyl-C3 and
     * hydroxyl keys, which were variant first.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void patternSetFragmentPropertiesTest() throws Exception {
        Assertions.assertEquals("[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",
                BiosynfoniFragmenter.BiosynfoniKey.D_PHENYL_C3.getSmarts(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniVariant.class)));
        Assertions.assertEquals("[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]",
                BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getSmarts(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniVariant.class)));
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        //isovaleric acid for the methylmalonyl key, ornithyl-alanine methyl ester for the non-standard amino acid key,
        // which is matched in both pattern sets, by the alanine methyl ester residue or by the ornithine residue
        BiosynfoniFragmenter.BiosynfoniKey[] tmpKeys = {BiosynfoniFragmenter.BiosynfoniKey.D2_METHYLMALONYL,
                BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS};
        String[] tmpSmilesArray = {"CC(C)CC(=O)O", "COC(=O)C(C)NC(=O)C(N)CCCN"};
        for (BiosynfoniFragmenter.BiosynfoniPatternSetOption tmpPatternSet : BiosynfoniFragmenter.BiosynfoniPatternSetOption.values()) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            Set<BiosynfoniFragmenter.BiosynfoniVariant> tmpVariants = tmpFragmenter.getAppliedVariants();
            for (int i = 0; i < tmpKeys.length; i++) {
                String tmpExpectedSmiles = tmpSmiGen.create(tmpSmiPar.parseSmiles(tmpKeys[i].getRepresentativeSmiles(tmpVariants)));
                int tmpCount = 0;
                for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmilesArray[i]))) {
                    if (tmpKeys[i].getLabel().equals(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY))) {
                        Assertions.assertEquals(tmpKeys[i].getSmarts(tmpVariants),
                                tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY), tmpPatternSet.name());
                        Assertions.assertEquals(tmpExpectedSmiles, tmpSmiGen.create(tmpFragment), tmpPatternSet.name());
                        tmpCount++;
                    }
                }
                Assertions.assertTrue(tmpCount > 0, tmpKeys[i].name() + ", " + tmpPatternSet.name());
            }
        }
        Assertions.assertNotEquals(
                BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getRepresentativeSmiles(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniVariant.class)),
                BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getRepresentativeSmiles(EnumSet.allOf(BiosynfoniFragmenter.BiosynfoniVariant.class)));
    }
    //
    /**
     * Tests the matching method property: ring keys matched by ring perception are marked as such, all other matches
     * as matched by their SMARTS pattern. Also checks that switching the pattern set on the same instance updates the
     * cached key data and that the cached representative structures are returned as new instances on every call.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void matchingMethodPropertyTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //naphthalene, with a methyl group for a key that is not a ring key
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("Cc1ccc2ccccc2c1");
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        for (BiosynfoniFragmenter.BiosynfoniPatternSetOption tmpPatternSet : new BiosynfoniFragmenter.BiosynfoniPatternSetOption[] {
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT, BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL,
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT}) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            boolean tmpIsRingPerceptionApplied = tmpPatternSet == BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT;
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(tmpFragments);
            Assertions.assertEquals(2, tmpCounts.get("r_c6"), tmpPatternSet.name());
            Assertions.assertEquals(tmpIsRingPerceptionApplied ? null : 1, tmpCounts.get("r_c10"), tmpPatternSet.name());
            for (IAtomContainer tmpFragment : tmpFragments) {
                String tmpLabel = tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY);
                boolean tmpIsRingKey = tmpLabel.startsWith("r_c");
                Assertions.assertEquals(tmpIsRingKey && tmpIsRingPerceptionApplied
                                ? BiosynfoniFragmenter.MATCHING_METHOD_RING_PERCEPTION : BiosynfoniFragmenter.MATCHING_METHOD_SMARTS,
                        tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_MATCHING_METHOD_PROPERTY_KEY),
                        tmpLabel + ", " + tmpPatternSet.name());
            }
            List<IAtomContainer> tmpFragmentsOfSecondCall = tmpFragmenter.fragmentMolecule(tmpMolecule);
            for (int i = 0; i < tmpFragments.size(); i++) {
                Assertions.assertNotSame(tmpFragments.get(i), tmpFragmentsOfSecondCall.get(i));
            }
        }
    }
    //
    /**
     * Tests the custom pattern set: with all variants turned on, it gives the same fragments as the variant pattern
     * set, and with all variants turned off the same as the original one; every variant changes exactly the SMARTS
     * patterns of the keys it is meant for; and the two variants of the phenyl unit and hydroxyl keys as well as the
     * ring perception work independently of each other.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void customPatternSetTest() throws Exception {
        //keys whose SMARTS pattern each variant changes; the ring perception does not change any SMARTS pattern
        Map<BiosynfoniFragmenter.BiosynfoniVariant, Set<String>> tmpExpectedKeys = new HashMap<>();
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.NADH, Set.of("CO_NADH"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.AMINO_ACIDS, Set.of("ALL_STD_AMINOS", "NON_STD_AMINOS"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.OPEN_CHAIN_SUGARS, Set.of("S_OPENPYR_C6O6", "S_OPENFUR_C5O5"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.PHENETHYLAMINE, Set.of("D_PHENYL_C2N"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.PHENYL_SIDE_CHAIN, Set.of("D_PHENYL_C3", "D_PHENYL_C2"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.AROMATIC_PHENYL_RINGS, Set.of("D_PHENYL_C3", "D_PHENYL_C2", "D_PHENYL_C1"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.METHYLMALONYL, Set.of("D2_METHYLMALONYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.PHOSPHATE_SULFONATE, Set.of("PHOSPHATE", "SULFONATE"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.TERMINAL_NITROGEN, Set.of("N_NITRATE"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.ETHER, Set.of("O_ETHER"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.HYDROXYL_HEMIACETAL, Set.of("O_HYDROXYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.HYDROXYL_NITRO_N_OXIDE, Set.of("O_HYDROXYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniVariant.RING_PERCEPTION, Set.of());
        Set<BiosynfoniFragmenter.BiosynfoniVariant> tmpNoVariants = EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniVariant.class);
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            Set<String> tmpChangedKeys = new HashSet<>();
            for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
                if (!tmpKey.getSmarts(tmpNoVariants).equals(tmpKey.getSmarts(EnumSet.of(tmpVariant)))) {
                    tmpChangedKeys.add(tmpKey.name());
                }
            }
            Assertions.assertEquals(tmpExpectedKeys.get(tmpVariant), tmpChangedKeys, tmpVariant.name());
        }
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpPresetFragmenter = new BiosynfoniFragmenter();
        BiosynfoniFragmenter tmpCustomFragmenter = new BiosynfoniFragmenter();
        tmpCustomFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.CUSTOM);
        String[] tmpSmilesArray = {BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES, "c1ccc2ccccc2c1",
                "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "COC(=O)C(C)NC(=O)C(N)CCCN"};
        for (boolean tmpIsApplied : new boolean[] {true, false}) {
            tmpPresetFragmenter.setPatternSetSetting(tmpIsApplied
                    ? BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT : BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
            for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
                tmpCustomFragmenter.setVariantSetting(tmpVariant, tmpIsApplied);
            }
            for (String tmpSmiles : tmpSmilesArray) {
                Assertions.assertEquals(
                        BiosynfoniFragmenterTest.countPerLabel(tmpPresetFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles))),
                        BiosynfoniFragmenterTest.countPerLabel(tmpCustomFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles))),
                        tmpIsApplied + ": " + tmpSmiles);
            }
        }
        //single variants: variant, key, SMILES code, number of matches
        String[][] tmpCases = {
                {"AROMATIC_PHENYL_RINGS", "D_PHENYL_C3", "CCCC1CCCCC1", "0"},
                {"AROMATIC_PHENYL_RINGS", "D_PHENYL_C3", "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "2"},
                {"PHENYL_SIDE_CHAIN", "D_PHENYL_C3", "CCCC1CCCCC1", "1"},
                {"PHENYL_SIDE_CHAIN", "D_PHENYL_C3", "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "0"},
                {"HYDROXYL_HEMIACETAL", "O_HYDROXYL", "OCC1OC(O)C(O)C(O)C1O", "5"},
                {"HYDROXYL_HEMIACETAL", "O_HYDROXYL", "C[N+](=O)[O-]", "1"},
                {"HYDROXYL_NITRO_N_OXIDE", "O_HYDROXYL", "OCC1OC(O)C(O)C(O)C1O", "4"},
                {"HYDROXYL_NITRO_N_OXIDE", "O_HYDROXYL", "C[N+](=O)[O-]", "0"},
                {"RING_PERCEPTION", "R_C10", "c1ccc2ccccc2c1", "0"},
                {"ETHER", "R_C10", "c1ccc2ccccc2c1", "1"}};
        for (String[] tmpCase : tmpCases) {
            for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
                tmpCustomFragmenter.setVariantSetting(tmpVariant, tmpVariant.name().equals(tmpCase[0]));
            }
            Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(
                    tmpCustomFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpCase[2])));
            Assertions.assertEquals(Integer.parseInt(tmpCase[3]),
                    tmpCounts.getOrDefault(BiosynfoniFragmenter.BiosynfoniKey.valueOf(tmpCase[1]).getLabel(), 0),
                    tmpCase[0] + ", " + tmpCase[1] + ": " + tmpCase[2]);
        }
    }
    //
    /**
     * Tests that the variant settings are only shown in the GUI while the custom pattern set is chosen, that all other
     * settings are always shown, and that a copy of the fragmenter has its own visibility conditions.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void variantSettingsVisibilityTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        Map<String, ObservableBooleanValue> tmpConditions = tmpFragmenter.getSettingNameToVisibilityConditionMap();
        Set<String> tmpVariantSettingNames = new HashSet<>();
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            tmpVariantSettingNames.add(tmpVariant.getSettingName());
        }
        Assertions.assertEquals(tmpVariantSettingNames, tmpConditions.keySet());
        for (BiosynfoniFragmenter.BiosynfoniPatternSetOption tmpPatternSet : BiosynfoniFragmenter.BiosynfoniPatternSetOption.values()) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            for (ObservableBooleanValue tmpCondition : tmpConditions.values()) {
                Assertions.assertEquals(tmpPatternSet == BiosynfoniFragmenter.BiosynfoniPatternSetOption.CUSTOM,
                        tmpCondition.get(), tmpPatternSet.name());
            }
        }
        tmpFragmenter.restoreDefaultSettings();
        for (ObservableBooleanValue tmpCondition : tmpConditions.values()) {
            Assertions.assertFalse(tmpCondition.get());
        }
        BiosynfoniFragmenter tmpCopy = (BiosynfoniFragmenter) tmpFragmenter.copy();
        tmpCopy.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.CUSTOM);
        for (String tmpSettingName : tmpVariantSettingNames) {
            Assertions.assertTrue(tmpCopy.getSettingNameToVisibilityConditionMap().get(tmpSettingName).get());
            Assertions.assertFalse(tmpConditions.get(tmpSettingName).get());
        }
    }
    //
    /**
     * Tests that the variant settings form a section of their own in the GUI: the only section title is attached to
     * the first variant setting, which directly follows the other settings and is hidden unless the custom pattern set
     * is chosen, so that the header bar of the section is hidden together with the variant settings.
     */
    @Test
    void variantSettingsSectionTitleTest() {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        Map<String, String> tmpSectionTitles = tmpFragmenter.getSettingNameToSectionTitleMap();
        String tmpFirstVariantSettingName = BiosynfoniFragmenter.BiosynfoniVariant.values()[0].getSettingName();
        Assertions.assertEquals(Set.of(tmpFirstVariantSettingName), tmpSectionTitles.keySet());
        Assertions.assertFalse(tmpSectionTitles.get(tmpFirstVariantSettingName).isBlank());
        List<String> tmpSettingNames = tmpFragmenter.settingsProperties().stream().map(Property::getName).toList();
        int tmpIndex = tmpSettingNames.indexOf(tmpFirstVariantSettingName);
        Assertions.assertEquals(tmpSettingNames.size() - BiosynfoniFragmenter.BiosynfoniVariant.values().length, tmpIndex);
        Assertions.assertTrue(tmpFragmenter.getSettingNameToVisibilityConditionMap().containsKey(tmpFirstVariantSettingName));
        Assertions.assertEquals(tmpSectionTitles, ((BiosynfoniFragmenter) tmpFragmenter.copy()).getSettingNameToSectionTitleMap());
    }
    //
    /**
     * Tests that every match of a key is returned as a fragment instance of its own in the representative structures
     * mode, so that altering one returned fragment does not alter the others. Butane contains three ethyl units.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void representativeFragmentInstancesTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        List<IAtomContainer> tmpEthylFragments = new ArrayList<>(3);
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("CCCC"))) {
            if (BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel().equals(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY))) {
                tmpEthylFragments.add(tmpFragment);
            }
        }
        Assertions.assertEquals(3, tmpEthylFragments.size());
        for (int i = 0; i < tmpEthylFragments.size(); i++) {
            for (int j = i + 1; j < tmpEthylFragments.size(); j++) {
                Assertions.assertNotSame(tmpEthylFragments.get(i), tmpEthylFragments.get(j));
                Assertions.assertNotSame(tmpEthylFragments.get(i).getAtom(0), tmpEthylFragments.get(j).getAtom(0));
            }
        }
        tmpEthylFragments.get(0).setProperty("Test property", "altered");
        Assertions.assertNull(tmpEthylFragments.get(1).getProperty("Test property"));
    }
    //
    /**
     * Tests that the original pattern set gives exactly the counts of the Python reference implementation of Biosynfoni
     * (version "full_1103") for ten diverse natural products, without overlap filters and with the inter-substructure
     * overlap filter only, the two configurations in which the counts are identical to the reference implementation
     * (see the class documentation of {@link BiosynfoniFragmenter}). This fails if an original SMARTS pattern or the key
     * order is changed by accident. The expected counts must not be adapted to the results of the fragmenter but only
     * be replaced by new results of the reference implementation.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void originalPatternSetReferenceCountsTest() throws Exception {
        //labels of the reference implementation in the order of the fingerprint positions
        String[] tmpReferenceLabels = {"co_coa", "co_nadh", "co_nadph", "allstnd_aminos", "nonstnd_aminos",
                "s_openpyr_C6O6", "s_openfur_C5O5", "s_pyranose_C5O4", "s_furanose_C4O3", "d_indoleC2N_12", "d_phenylC2N_9",
                "d_c5n_6", "d_c4n_5", "d_phenylC3_9_strict", "d_phenylC2_8_strict", "d_phenylC1_7_strict", "d_isoprene_5",
                "d2_acetyl_C2O1", "d2_methylmalonyl_C3", "d_ethyl_2", "d_methyl_1", "phosphate_2", "sulfonate_2", "hal_f",
                "hal_cl", "hal_br", "hal_i", "n_nitrate_1", "o_epoxy_1", "o_ether_1", "o_hydroxyl_1", "r_c3", "r_c4", "r_c5",
                "r_c6", "r_c7", "r_c8", "r_c9", "r_c10"};
        //sucrose, glutathione, chlorogenic acid, cholesterol, ATP, yohimbine, emodin, erythromycin A, griseofulvin, and
        // quercetin 3-sulfate
        String[] tmpSmilesArray = {
                "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
                "NC(CCC(=O)NC(CS)C(=O)NCC(=O)O)C(=O)O",
                "OC(=O)C1(O)CC(O)C(O)C(C1)OC(=O)C=Cc1ccc(O)c(O)c1",
                "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
                "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
                "COC(=O)C1C(O)CCC2CN3CCc4c([nH]c5ccccc45)C3CC12",
                "Cc1cc(O)c2C(=O)c3c(O)cc(O)cc3C(=O)c2c1",
                "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O",
                "COC1=CC(=O)CC(C)C12Oc1c(Cl)c(OC)cc(OC)c1C2=O",
                "O=c1c(OS(=O)(=O)O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"};
        //count vectors of the Python reference implementation (https://github.com/lucinamay/biosynfoni, commit 0ba4620,
        // RDKit 2025.03.6), in the order of the fingerprint positions, computed as
        // biosynfoni.concerto_fp.Biosynfoni(Chem.MolFromSmiles(smiles), version_name="full_1103", intersub_overlap=True,
        // intrasub_overlap=True).fingerprint without overlap filters and with intersub_overlap=False for the
        // inter-substructure overlap filter only
        int[][] tmpNoFilterCounts = {
                {0, 0, 0, 0, 0, 2, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 23, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 7, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 5, 17, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 2, 11, 30, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 1, 2},
                {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 6, 0, 12, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 19, 4, 0, 21, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 10, 2, 17, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 2},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 33, 17, 31, 13, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 15, 2, 15, 4, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6, 16, 0, 16, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2, 0, 0, 0, 0}};
        int[][] tmpInterFilterCounts = {
                {0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 6, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0}};
        BiosynfoniFragmenter.BiosynfoniKey[] tmpKeys = BiosynfoniFragmenter.BiosynfoniKey.values();
        Assertions.assertEquals(tmpReferenceLabels.length, tmpKeys.length);
        for (int i = 0; i < tmpKeys.length; i++) {
            Assertions.assertEquals(tmpReferenceLabels[i], tmpKeys[i].getLabel());
        }
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
        for (boolean tmpIsInterFilterApplied : new boolean[] {false, true}) {
            tmpFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
            int[][] tmpExpectedCounts = tmpIsInterFilterApplied ? tmpInterFilterCounts : tmpNoFilterCounts;
            for (int i = 0; i < tmpSmilesArray.length; i++) {
                Map<String, Integer> tmpCounts = BiosynfoniFragmenterTest.countPerLabel(
                        tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmilesArray[i])));
                int[] tmpCountVector = new int[tmpKeys.length];
                for (int j = 0; j < tmpKeys.length; j++) {
                    tmpCountVector[j] = tmpCounts.getOrDefault(tmpKeys[j].getLabel(), 0);
                }
                Assertions.assertArrayEquals(tmpExpectedCounts[i], tmpCountVector,
                        (tmpIsInterFilterApplied ? "inter-substructure overlap filter: " : "no overlap filter: ") + tmpSmilesArray[i]);
            }
        }
    }
    //
    /**
     * Tests that the counts with the intra-substructure overlap filter do not depend on the input atom and bond order,
     * for both preset pattern sets, with and without the inter-substructure overlap filter: the counts of several seeded
     * random permutations of the atoms and bonds must equal those of the parsed molecule. The molecules are symmetric
     * (mesitylene, cubane), fused (cholesterol), disconnected with two identical components (two butyric acid
     * molecules), and with many overlapping ethyl, acetyl, and isoprene matches (farnesol and the COCONUT natural
     * product CNP0151033).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void intraSubstructureOverlapFilterAtomOrderTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String[] tmpSmilesArray = {"Cc1cc(C)cc(C)c1", "C12C3C4C1C5C2C3C45", "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
                "CCCC(=O)O.CCCC(=O)O", "CC(C)=CCCC(C)=CCCC(C)=CCO", BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES};
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        for (BiosynfoniFragmenter.BiosynfoniPatternSetOption tmpPatternSet : new BiosynfoniFragmenter.BiosynfoniPatternSetOption[] {
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT, BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL}) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            for (boolean tmpIsInterFilterApplied : new boolean[] {false, true}) {
                tmpFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
                for (String tmpSmiles : tmpSmilesArray) {
                    BiosynfoniFragmenterTest.assertCountsIndependentOfAtomOrder(tmpFragmenter, tmpSmiPar.parseSmiles(tmpSmiles),
                            tmpPatternSet.name() + ", inter-substructure overlap filter " + tmpIsInterFilterApplied + ": " + tmpSmiles);
                }
            }
        }
    }
    //
    /**
     * Tests the fallback of the canonical atom ranks for the intra-substructure overlap filter, which is used if no
     * unique SMILES code can be created, here because of an unset order of a non-aromatic bond: the filter still
     * reduces the number of fragments, and the counts do not depend on the input atom and bond order, for both preset
     * pattern sets, with and without the inter-substructure overlap filter.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void canonicalRanksFallbackTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //farnesol with an unset order of the bond between its first two carbon atoms
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("CC(C)=CCCC(C)=CCCC(C)=CCO");
        tmpMolecule.getBond(tmpMolecule.getAtom(0), tmpMolecule.getAtom(1)).setOrder(IBond.Order.UNSET);
        //precondition: no unique SMILES code, so the fallback is used
        Assertions.assertThrows(CDKException.class, () -> new SmilesGenerator(SmiFlavor.Unique).create(tmpMolecule));
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        int tmpUnfilteredCount = tmpFragmenter.fragmentMolecule(tmpMolecule).size();
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(tmpMolecule).size() < tmpUnfilteredCount);
        for (BiosynfoniFragmenter.BiosynfoniPatternSetOption tmpPatternSet : new BiosynfoniFragmenter.BiosynfoniPatternSetOption[] {
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.VARIANT, BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL}) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            for (boolean tmpIsInterFilterApplied : new boolean[] {false, true}) {
                tmpFragmenter.setInterSubstructureOverlapFilterSetting(tmpIsInterFilterApplied);
                BiosynfoniFragmenterTest.assertCountsIndependentOfAtomOrder(tmpFragmenter, tmpMolecule,
                        tmpPatternSet.name() + ", inter-substructure overlap filter " + tmpIsInterFilterApplied);
            }
        }
    }
    //
    /**
     * Tests that a molecule with unset implicit hydrogen counts, for which neither a unique SMILES code nor the default
     * invariants of the CDK canonical labelling are available, can be fragmented with the intra-substructure overlap
     * filter by means of the fallback of the canonical atom ranks, and that its counts do not depend on the input atom
     * and bond order. MORTAR sets unset implicit hydrogen counts on import, so this only concerns the direct use of the
     * fragmenter. Molecules with rings cannot be fragmented with unset implicit hydrogen counts at all, regardless of
     * the overlap filters, because the aromaticity perception of the SMARTS matching requires them, so an acyclic
     * molecule is used.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void unsetImplicitHydrogenCountsTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //1-butanol with unset implicit hydrogen counts
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("CCCCO");
        for (IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpAtom.setImplicitHydrogenCount(null);
        }
        //precondition: no unique SMILES code, so the fallback is used
        Assertions.assertThrows(NullPointerException.class, () -> new SmilesGenerator(SmiFlavor.Unique).create(tmpMolecule));
        String tmpEthylLabel = BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel();
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        //the chain of four carbon atoms contains three overlapping ethyl units, of which one or two are disjoint
        Assertions.assertEquals(3, BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpMolecule)).get(tmpEthylLabel));
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(true);
        int tmpFilteredEthylCount = BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpMolecule)).get(tmpEthylLabel);
        Assertions.assertTrue(tmpFilteredEthylCount == 1 || tmpFilteredEthylCount == 2, String.valueOf(tmpFilteredEthylCount));
        BiosynfoniFragmenterTest.assertCountsIndependentOfAtomOrder(tmpFragmenter, tmpMolecule, "CCCCO");
    }
    //
    /**
     * Tests that all 17 settings, each set to a value other than its default, are restored by a new fragmenter after
     * they have been written to and read from a preference container file as MORTAR persists them (see
     * FragmentationService.persistFragmenterSettings() and reloadFragmenterSettings()). The file is written to a
     * temporary directory.
     *
     * @param aTempDirectory temporary directory for the preference container file, provided by JUnit
     * @throws Exception if anything goes wrong
     */
    @Test
    void settingsPersistenceTest(@TempDir Path aTempDirectory) throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.CUSTOM);
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setInterSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            tmpFragmenter.setVariantSetting(tmpVariant, !BiosynfoniFragmenter.VARIANT_SETTING_DEFAULT);
        }
        List<Property<?>> tmpSettings = tmpFragmenter.settingsProperties();
        List<Property<?>> tmpDefaultSettings = new BiosynfoniFragmenter().settingsProperties();
        Assertions.assertEquals(17, tmpSettings.size());
        for (int i = 0; i < tmpSettings.size(); i++) {
            Assertions.assertNotEquals(tmpDefaultSettings.get(i).getValue(), tmpSettings.get(i).getValue(), tmpSettings.get(i).getName());
        }
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpSettings));
        String tmpFilePath = aTempDirectory.resolve(BiosynfoniFragmenter.class.getSimpleName()
                + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION).toString();
        PreferenceContainer tmpContainer = PreferenceUtil.translateJavaFxPropertiesToPreferences(tmpSettings, tmpFilePath);
        Assertions.assertEquals(tmpSettings.size(), tmpContainer.getSize());
        tmpContainer.writeRepresentation();
        BiosynfoniFragmenter tmpReloadedFragmenter = new BiosynfoniFragmenter();
        PreferenceUtil.updatePropertiesFromPreferences(tmpReloadedFragmenter.settingsProperties(), new PreferenceContainer(new File(tmpFilePath)));
        List<Property<?>> tmpReloadedSettings = tmpReloadedFragmenter.settingsProperties();
        Assertions.assertEquals(tmpSettings.size(), tmpReloadedSettings.size());
        for (int i = 0; i < tmpSettings.size(); i++) {
            Assertions.assertEquals(tmpSettings.get(i).getName(), tmpReloadedSettings.get(i).getName());
            Assertions.assertEquals(tmpSettings.get(i).getValue(), tmpReloadedSettings.get(i).getValue(), tmpSettings.get(i).getName());
        }
    }
    //
    /**
     * Tests the behaviour for unusual input molecules: an empty molecule is filtered and cannot be fragmented; single
     * atoms that no key matches give no fragments in all modes; the counts of a disconnected molecule are the sums of
     * the counts of its components; pseudo atoms are not matched as carbon atoms; and explicit hydrogen atoms count as
     * neighbours, so that the oxygen atom of ethanol with an explicit hydrogen atom on it is matched as ether instead of
     * hydroxyl oxygen (see the class documentation of {@link BiosynfoniFragmenter}).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void unusualInputTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        IAtomContainer tmpEmptyMolecule = SilentChemObjectBuilder.getInstance().newAtomContainer();
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpEmptyMolecule));
        Assertions.assertFalse(tmpFragmenter.canBeFragmented(tmpEmptyMolecule));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.fragmentMolecule(tmpEmptyMolecule));
        //methane, water, and an iron atom
        for (boolean tmpIsIntraFilterApplied : new boolean[] {false, true}) {
            tmpFragmenter.setIntraSubstructureOverlapFilterSetting(tmpIsIntraFilterApplied);
            for (BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption tmpOption : BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.values()) {
                tmpFragmenter.setReturnedFragmentsSetting(tmpOption);
                for (String tmpSmiles : new String[] {"C", "O", "[Fe]"}) {
                    Assertions.assertTrue(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles)).isEmpty(),
                            tmpOption.name() + ", intra-substructure overlap filter " + tmpIsIntraFilterApplied + ": " + tmpSmiles);
                }
            }
        }
        tmpFragmenter.restoreDefaultSettings();
        //acetic acid and ethylamine as one disconnected molecule
        Map<String, Integer> tmpExpectedCounts = BiosynfoniFragmenterTest.countPerLabel(
                tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("CC(=O)O")));
        BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("CCN")))
                .forEach((aLabel, aCount) -> tmpExpectedCounts.merge(aLabel, aCount, Integer::sum));
        Assertions.assertEquals(tmpExpectedCounts, BiosynfoniFragmenterTest.countPerLabel(
                tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("CC(=O)O.CCN"))));
        //only the bond between the two carbon atoms is matched, as ethyl unit
        Assertions.assertEquals(Map.of(BiosynfoniFragmenter.BiosynfoniKey.D_ETHYL.getLabel(), 1),
                BiosynfoniFragmenterTest.countPerLabel(tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("*CC*"))));
        String tmpEtherLabel = BiosynfoniFragmenter.BiosynfoniKey.O_ETHER.getLabel();
        String tmpHydroxylLabel = BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getLabel();
        //ethanol with an implicit and with an explicit hydrogen atom on the oxygen atom
        Map<String, Integer> tmpImplicitHydrogenCounts = BiosynfoniFragmenterTest.countPerLabel(
                tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("OCC")));
        Assertions.assertEquals(0, tmpImplicitHydrogenCounts.getOrDefault(tmpEtherLabel, 0));
        Assertions.assertEquals(1, tmpImplicitHydrogenCounts.getOrDefault(tmpHydroxylLabel, 0));
        Map<String, Integer> tmpExplicitHydrogenCounts = BiosynfoniFragmenterTest.countPerLabel(
                tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles("[H]OCC")));
        Assertions.assertEquals(1, tmpExplicitHydrogenCounts.getOrDefault(tmpEtherLabel, 0));
        Assertions.assertEquals(0, tmpExplicitHydrogenCounts.getOrDefault(tmpHydroxylLabel, 0));
    }
    //
    /**
     * Returns sets of variants for testing: no variants, all variants, every single variant, and all variants but
     * one.
     *
     * @return list of variant sets
     */
    private static List<Set<BiosynfoniFragmenter.BiosynfoniVariant>> getTestVariantSets() {
        List<Set<BiosynfoniFragmenter.BiosynfoniVariant>> tmpVariantSets = new ArrayList<>();
        tmpVariantSets.add(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniVariant.class));
        tmpVariantSets.add(EnumSet.allOf(BiosynfoniFragmenter.BiosynfoniVariant.class));
        for (BiosynfoniFragmenter.BiosynfoniVariant tmpVariant : BiosynfoniFragmenter.BiosynfoniVariant.values()) {
            tmpVariantSets.add(EnumSet.of(tmpVariant));
            tmpVariantSets.add(EnumSet.complementOf(EnumSet.of(tmpVariant)));
        }
        return tmpVariantSets;
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
    }
    //
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
    //
    /**
     * Asserts that the given fragmenter returns the same number of fragments per key for the given molecule and for
     * several seeded random permutations of its atoms and bonds.
     *
     * @param aFragmenter the fragmenter with the settings to test
     * @param aMolecule the molecule, it is not altered
     * @param aMessage message for failed assertions
     * @throws Exception if anything goes wrong
     */
    private static void assertCountsIndependentOfAtomOrder(BiosynfoniFragmenter aFragmenter, IAtomContainer aMolecule,
                                                           String aMessage) throws Exception {
        Map<String, Integer> tmpExpectedCounts = BiosynfoniFragmenterTest.countPerLabel(aFragmenter.fragmentMolecule(aMolecule));
        Random tmpRandom = new Random(BiosynfoniFragmenterTest.PERMUTATION_SEED);
        for (int i = 0; i < BiosynfoniFragmenterTest.NUMBER_OF_PERMUTATIONS; i++) {
            IAtomContainer tmpPermutedMolecule = BiosynfoniFragmenterTest.permuteAtomsAndBonds(aMolecule, tmpRandom);
            Assertions.assertEquals(tmpExpectedCounts,
                    BiosynfoniFragmenterTest.countPerLabel(aFragmenter.fragmentMolecule(tmpPermutedMolecule)),
                    aMessage + ", permutation " + i);
        }
    }
    //
    /**
     * Returns a clone of the given molecule whose atoms and bonds are randomly reordered.
     *
     * @param aMolecule the molecule, it is not altered
     * @param aRandom random number generator for the permutations
     * @return reordered clone of the molecule
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    private static IAtomContainer permuteAtomsAndBonds(IAtomContainer aMolecule, Random aRandom) throws CloneNotSupportedException {
        IAtomContainer tmpClone = aMolecule.clone();
        List<IAtom> tmpAtoms = new ArrayList<>(tmpClone.getAtomCount());
        for (IAtom tmpAtom : tmpClone.atoms()) {
            tmpAtoms.add(tmpAtom);
        }
        List<IBond> tmpBonds = new ArrayList<>(tmpClone.getBondCount());
        for (IBond tmpBond : tmpClone.bonds()) {
            tmpBonds.add(tmpBond);
        }
        Collections.shuffle(tmpAtoms, aRandom);
        Collections.shuffle(tmpBonds, aRandom);
        tmpClone.setAtoms(tmpAtoms.toArray(new IAtom[0]));
        tmpClone.setBonds(tmpBonds.toArray(new IBond[0]));
        return tmpClone;
    }
}
