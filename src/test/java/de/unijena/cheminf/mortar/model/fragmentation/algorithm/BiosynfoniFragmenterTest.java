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
import javafx.beans.value.ObservableBooleanValue;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.EnumSet;
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
        Assertions.assertEquals(BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT, tmpFragmenter.getPatternSetSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT,
                tmpFragmenter.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
        tmpFragmenter.setIntraSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTRA_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setInterSubstructureOverlapFilterSetting(!BiosynfoniFragmenter.INTER_SUBSTRUCTURE_OVERLAP_FILTER_SETTING_DEFAULT);
        tmpFragmenter.setReturnedFragmentsSetting(BiosynfoniFragmenter.BiosynfoniReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
            Assertions.assertEquals(BiosynfoniFragmenter.REVISION_SETTING_DEFAULT, tmpFragmenter.getRevisionSetting(tmpRevision));
        }
        tmpFragmenter.setRevisionSetting(BiosynfoniFragmenter.BiosynfoniRevision.ETHER, !BiosynfoniFragmenter.REVISION_SETTING_DEFAULT);
        BiosynfoniFragmenter tmpCopy = (BiosynfoniFragmenter) tmpFragmenter.copy();
        for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
            Assertions.assertEquals(tmpFragmenter.getRevisionSetting(tmpRevision), tmpCopy.getRevisionSetting(tmpRevision));
        }
        Assertions.assertEquals(tmpFragmenter.getPatternSetSetting(), tmpCopy.getPatternSetSetting());
        Assertions.assertEquals(tmpFragmenter.getReturnedFragmentsSetting(), tmpCopy.getReturnedFragmentsSetting());
        Assertions.assertEquals(tmpFragmenter.getIntraSubstructureOverlapFilterSetting(), tmpCopy.getIntraSubstructureOverlapFilterSetting());
        Assertions.assertEquals(tmpFragmenter.getInterSubstructureOverlapFilterSetting(), tmpCopy.getInterSubstructureOverlapFilterSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(BiosynfoniFragmenter.PATTERN_SET_OPTION_DEFAULT, tmpFragmenter.getPatternSetSetting());
        Assertions.assertEquals(BiosynfoniFragmenter.REVISION_SETTING_DEFAULT,
                tmpFragmenter.getRevisionSetting(BiosynfoniFragmenter.BiosynfoniRevision.ETHER));
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
     * for no revisions, all revisions, every single revision, and all revisions but one, which covers all combinations
     * of the two revisions of the phenyl unit and hydroxyl keys.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void representativeStructuresTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        for (Set<BiosynfoniFragmenter.BiosynfoniRevision> tmpRevisions : BiosynfoniFragmenterTest.getTestRevisionSets()) {
            Set<String> tmpUniqueSmilesSet = new HashSet<>();
            Set<String> tmpDisplayNameSet = new HashSet<>();
            for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
                String tmpMessage = tmpKey.getLabel() + ", " + tmpRevisions;
                IAtomContainer tmpRepresentative = tmpSmiPar.parseSmiles(tmpKey.getRepresentativeSmiles(tmpRevisions));
                Assertions.assertTrue(tmpUniqueSmilesSet.add(tmpSmiGen.create(tmpRepresentative)),
                        "Representative structure is not unique: " + tmpMessage);
                Assertions.assertTrue(SmartsPattern.create(tmpKey.getSmarts(tmpRevisions)).matches(tmpRepresentative),
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
     * Tests every revision of the revised pattern set against the original pattern set: for each case, the number of
     * matches of the given key in the given molecule is checked in both pattern sets, without overlap filters. The
     * cases cover the corrected errors (e.g. open-chain sugar patterns matching sugar rings, the phenethylamine pattern
     * matching indoles, ring patterns matching the perimeter of fused ring systems) and the chemical decisions (e.g.
     * aromatic phenyl rings, methylmalonyl units at branching points), as well as unchanged behaviour.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void patternSetTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //key, SMILES code, number of matches with the original patterns, number of matches with the revised patterns
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
            tmpFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.REVISED);
            Map<String, Integer> tmpRevisedCounts = BiosynfoniFragmenterTest.countPerLabel(
                    tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpCase[1])));
            Assertions.assertEquals(Integer.parseInt(tmpCase[2]), tmpOriginalCounts.getOrDefault(tmpLabel, 0),
                    tmpCase[0] + ", original: " + tmpCase[1]);
            Assertions.assertEquals(Integer.parseInt(tmpCase[3]), tmpRevisedCounts.getOrDefault(tmpLabel, 0),
                    tmpCase[0] + ", revised: " + tmpCase[1]);
        }
    }
    //
    /**
     * Tests that the pattern set setting changes the SMARTS pattern and the representative structure that are attached
     * to the returned fragments, for the two keys whose representative structures differ between the pattern sets, and
     * that the original pattern set uses the original patterns of the reference implementation for the phenyl-C3 and
     * hydroxyl keys, which were revised first.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void patternSetFragmentPropertiesTest() throws Exception {
        Assertions.assertEquals("[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",
                BiosynfoniFragmenter.BiosynfoniKey.D_PHENYL_C3.getSmarts(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniRevision.class)));
        Assertions.assertEquals("[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]",
                BiosynfoniFragmenter.BiosynfoniKey.O_HYDROXYL.getSmarts(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniRevision.class)));
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
            Set<BiosynfoniFragmenter.BiosynfoniRevision> tmpRevisions = tmpFragmenter.getAppliedRevisions();
            for (int i = 0; i < tmpKeys.length; i++) {
                String tmpExpectedSmiles = tmpSmiGen.create(tmpSmiPar.parseSmiles(tmpKeys[i].getRepresentativeSmiles(tmpRevisions)));
                int tmpCount = 0;
                for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmilesArray[i]))) {
                    if (tmpKeys[i].getLabel().equals(tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_LABEL_PROPERTY_KEY))) {
                        Assertions.assertEquals(tmpKeys[i].getSmarts(tmpRevisions),
                                tmpFragment.getProperty(BiosynfoniFragmenter.BIOSYNFONI_SMARTS_PROPERTY_KEY), tmpPatternSet.name());
                        Assertions.assertEquals(tmpExpectedSmiles, tmpSmiGen.create(tmpFragment), tmpPatternSet.name());
                        tmpCount++;
                    }
                }
                Assertions.assertTrue(tmpCount > 0, tmpKeys[i].name() + ", " + tmpPatternSet.name());
            }
        }
        Assertions.assertNotEquals(
                BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getRepresentativeSmiles(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniRevision.class)),
                BiosynfoniFragmenter.BiosynfoniKey.NON_STD_AMINOS.getRepresentativeSmiles(EnumSet.allOf(BiosynfoniFragmenter.BiosynfoniRevision.class)));
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
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.REVISED, BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL,
                BiosynfoniFragmenter.BiosynfoniPatternSetOption.REVISED}) {
            tmpFragmenter.setPatternSetSetting(tmpPatternSet);
            boolean tmpIsRingPerceptionApplied = tmpPatternSet == BiosynfoniFragmenter.BiosynfoniPatternSetOption.REVISED;
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
     * Tests the custom pattern set: with all revisions turned on, it gives the same fragments as the revised pattern
     * set, and with all revisions turned off the same as the original one; every revision changes exactly the SMARTS
     * patterns of the keys it is meant for; and the two revisions of the phenyl unit and hydroxyl keys as well as the
     * ring perception work independently of each other.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void customPatternSetTest() throws Exception {
        //keys whose SMARTS pattern each revision changes; the ring perception does not change any SMARTS pattern
        Map<BiosynfoniFragmenter.BiosynfoniRevision, Set<String>> tmpExpectedKeys = new HashMap<>();
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.NADH, Set.of("CO_NADH"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.AMINO_ACIDS, Set.of("ALL_STD_AMINOS", "NON_STD_AMINOS"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.OPEN_CHAIN_SUGARS, Set.of("S_OPENPYR_C6O6", "S_OPENFUR_C5O5"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.PHENETHYLAMINE, Set.of("D_PHENYL_C2N"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.PHENYL_SIDE_CHAIN, Set.of("D_PHENYL_C3", "D_PHENYL_C2"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.AROMATIC_PHENYL_RINGS, Set.of("D_PHENYL_C3", "D_PHENYL_C2", "D_PHENYL_C1"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.METHYLMALONYL, Set.of("D2_METHYLMALONYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.PHOSPHATE_SULFONATE, Set.of("PHOSPHATE", "SULFONATE"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.TERMINAL_NITROGEN, Set.of("N_NITRATE"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.ETHER, Set.of("O_ETHER"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.HYDROXYL_HEMIACETAL, Set.of("O_HYDROXYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.HYDROXYL_NITRO_N_OXIDE, Set.of("O_HYDROXYL"));
        tmpExpectedKeys.put(BiosynfoniFragmenter.BiosynfoniRevision.RING_PERCEPTION, Set.of());
        Set<BiosynfoniFragmenter.BiosynfoniRevision> tmpNoRevisions = EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniRevision.class);
        for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
            Set<String> tmpChangedKeys = new HashSet<>();
            for (BiosynfoniFragmenter.BiosynfoniKey tmpKey : BiosynfoniFragmenter.BiosynfoniKey.values()) {
                if (!tmpKey.getSmarts(tmpNoRevisions).equals(tmpKey.getSmarts(EnumSet.of(tmpRevision)))) {
                    tmpChangedKeys.add(tmpKey.name());
                }
            }
            Assertions.assertEquals(tmpExpectedKeys.get(tmpRevision), tmpChangedKeys, tmpRevision.name());
        }
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        BiosynfoniFragmenter tmpPresetFragmenter = new BiosynfoniFragmenter();
        BiosynfoniFragmenter tmpCustomFragmenter = new BiosynfoniFragmenter();
        tmpCustomFragmenter.setPatternSetSetting(BiosynfoniFragmenter.BiosynfoniPatternSetOption.CUSTOM);
        String[] tmpSmilesArray = {BiosynfoniFragmenterTest.TEST_MOLECULE_SMILES, "c1ccc2ccccc2c1",
                "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "COC(=O)C(C)NC(=O)C(N)CCCN"};
        for (boolean tmpIsApplied : new boolean[] {true, false}) {
            tmpPresetFragmenter.setPatternSetSetting(tmpIsApplied
                    ? BiosynfoniFragmenter.BiosynfoniPatternSetOption.REVISED : BiosynfoniFragmenter.BiosynfoniPatternSetOption.ORIGINAL);
            for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
                tmpCustomFragmenter.setRevisionSetting(tmpRevision, tmpIsApplied);
            }
            for (String tmpSmiles : tmpSmilesArray) {
                Assertions.assertEquals(
                        BiosynfoniFragmenterTest.countPerLabel(tmpPresetFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles))),
                        BiosynfoniFragmenterTest.countPerLabel(tmpCustomFragmenter.fragmentMolecule(tmpSmiPar.parseSmiles(tmpSmiles))),
                        tmpIsApplied + ": " + tmpSmiles);
            }
        }
        //single revisions: revision, key, SMILES code, number of matches
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
            for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
                tmpCustomFragmenter.setRevisionSetting(tmpRevision, tmpRevision.name().equals(tmpCase[0]));
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
     * Tests that the revision settings are only shown in the GUI while the custom pattern set is chosen, that all other
     * settings are always shown, and that a copy of the fragmenter has its own visibility conditions.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void revisionSettingsVisibilityTest() throws Exception {
        BiosynfoniFragmenter tmpFragmenter = new BiosynfoniFragmenter();
        Map<String, ObservableBooleanValue> tmpConditions = tmpFragmenter.getSettingNameToVisibilityConditionMap();
        Set<String> tmpRevisionSettingNames = new HashSet<>();
        for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
            tmpRevisionSettingNames.add(tmpRevision.getSettingName());
        }
        Assertions.assertEquals(tmpRevisionSettingNames, tmpConditions.keySet());
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
        for (String tmpSettingName : tmpRevisionSettingNames) {
            Assertions.assertTrue(tmpCopy.getSettingNameToVisibilityConditionMap().get(tmpSettingName).get());
            Assertions.assertFalse(tmpConditions.get(tmpSettingName).get());
        }
    }
    //
    /**
     * Returns sets of revisions for testing: no revisions, all revisions, every single revision, and all revisions but
     * one.
     *
     * @return list of revision sets
     */
    private static List<Set<BiosynfoniFragmenter.BiosynfoniRevision>> getTestRevisionSets() {
        List<Set<BiosynfoniFragmenter.BiosynfoniRevision>> tmpRevisionSets = new ArrayList<>();
        tmpRevisionSets.add(EnumSet.noneOf(BiosynfoniFragmenter.BiosynfoniRevision.class));
        tmpRevisionSets.add(EnumSet.allOf(BiosynfoniFragmenter.BiosynfoniRevision.class));
        for (BiosynfoniFragmenter.BiosynfoniRevision tmpRevision : BiosynfoniFragmenter.BiosynfoniRevision.values()) {
            tmpRevisionSets.add(EnumSet.of(tmpRevision));
            tmpRevisionSets.add(EnumSet.complementOf(EnumSet.of(tmpRevision)));
        }
        return tmpRevisionSets;
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
}
