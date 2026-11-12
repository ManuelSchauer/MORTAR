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

import de.unijena.cheminf.mortar.model.fragmentation.FragmentationService;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;
import de.unijena.cheminf.mortar.preference.SingleTermPreference;

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

/**
 * Class to test the correct working of {@link SmartsFragmenter}.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SmartsFragmenterTest {
    /**
     * SMILES parser that assigns a Kekulé structure and keeps the aromaticity flags of aromatic SMILES.
     */
    private final SmilesParser smilesParser;
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of
     * the fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public SmartsFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    }
    //
    /**
     * Tests instantiation, default settings, names, and that the settings can be persisted.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void basicTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        Assertions.assertEquals(SmartsFragmenter.ALGORITHM_NAME, tmpFragmenter.getFragmentationAlgorithmName());
        Assertions.assertEquals("Custom SMARTS", tmpFragmenter.getFragmentationAlgorithmDisplayName());
        Assertions.assertTrue(PreferenceUtil.isValidName(tmpFragmenter.getFragmentationAlgorithmName()));
        Assertions.assertTrue(SingleTermPreference.isValidContent(tmpFragmenter.getFragmentationAlgorithmName()));
        Assertions.assertTrue(PreferenceUtil.isValidName(tmpFragmenter.getFragmentationAlgorithmDisplayName()));
        Assertions.assertTrue(SingleTermPreference.isValidContent(tmpFragmenter.getFragmentationAlgorithmDisplayName()));
        Assertions.assertEquals("", tmpFragmenter.getPatternsSetting());
        Assertions.assertEquals(SmartsFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, tmpFragmenter.getReturnedFragmentsSetting());
        Assertions.assertEquals(SmartsFragmenter.FILTER_OVERLAPPING_MATCHES_SETTING_DEFAULT, tmpFragmenter.getFilterOverlappingMatchesSetting());
        Assertions.assertEquals(IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT, tmpFragmenter.getFragmentSaturationSetting());
        Assertions.assertEquals(4, tmpFragmenter.settingsProperties().size());
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
        }
        //the empty default and patterns with special characters must be persistable
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpFragmenter.settingsProperties()));
        tmpFragmenter.setPatternsSetting("[#6;R1]~[$([OX2H]),$([O-])] ring carbon with OH or O-\n*!@[N+](=O)[O-] nitro");
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpFragmenter.settingsProperties()));
    }
    //
    /**
     * Tests that the fragmenter is registered in the fragmentation service, passes its checks of the settings, and can
     * be instantiated by its algorithm name, as done for pipeline fragmenters.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void fragmentationServiceRegistrationTest() throws Exception {
        //the constructor checks all fragmenters and their settings for persistability
        FragmentationService tmpService = new FragmentationService();
        long tmpCount = 0;
        for (IMoleculeFragmenter tmpFragmenter : tmpService.getFragmenters()) {
            if (tmpFragmenter instanceof SmartsFragmenter) {
                tmpCount++;
            }
        }
        Assertions.assertEquals(1, tmpCount);
        Assertions.assertInstanceOf(SmartsFragmenter.class,
                tmpService.createNewFragmenterObjectByAlgorithmName(SmartsFragmenter.ALGORITHM_NAME));
    }
    //
    /**
     * Tests the extraction of named and unnamed matches and of the remainder, and their categories.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void matchingAndCategoriesTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[CX3](=O)[OX2H1] carboxylic acid\n[NX3;H2]");
        //4-aminobutanoic acid (GABA)
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("NCCCC(=O)O");
        String tmpSmilesBefore = this.createSmiles(tmpMolecule);
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertEquals(this.normalisedSmiles("OC=O", "N", "CCC"), this.fragmentSmiles(tmpFragments));
        Map<String, String> tmpSmilesToCategory = new HashMap<>();
        for (IAtomContainer tmpFragment : tmpFragments) {
            tmpSmilesToCategory.put(this.createSmiles(tmpFragment),
                    tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY));
        }
        Assertions.assertEquals("carboxylic acid", tmpSmilesToCategory.get(this.normalisedSmiles("OC=O").getFirst()));
        //unnamed pattern, category is the SMARTS pattern
        Assertions.assertEquals("[NX3;H2]", tmpSmilesToCategory.get(this.normalisedSmiles("N").getFirst()));
        Assertions.assertEquals(SmartsFragmenter.FRAGMENT_CATEGORY_REMAINDER_VALUE,
                tmpSmilesToCategory.get(this.normalisedSmiles("CCC").getFirst()));
        //the input molecule is not altered
        Assertions.assertEquals(tmpSmilesBefore, this.createSmiles(tmpMolecule));
        //only matches
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Assertions.assertEquals(this.normalisedSmiles("OC=O", "N"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        //no match at all
        IAtomContainer tmpPropane = this.smilesParser.parseSmiles("CCC");
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(tmpPropane).isEmpty());
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES_AND_REMAINDER);
        Assertions.assertEquals(this.normalisedSmiles("CCC"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpPropane)));
    }
    //
    /**
     * Tests the overlap filter: with the filter, the order of the patterns is their priority and matches of the same
     * pattern do not overlap; without the filter, all unique matches are returned.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void overlapFilterTest() throws Exception {
        //lactic acid
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("CC(O)C(=O)O");
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[CX3](=O)[OX2H1] acid\n[OX2H1] hydroxyl");
        Assertions.assertEquals(this.normalisedSmiles("OC=O", "O", "CC"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        //reversed priority, the acid cannot match anymore because its OH group is taken
        tmpFragmenter.setPatternsSetting("[OX2H1] hydroxyl\n[CX3](=O)[OX2H1] acid");
        Assertions.assertEquals(this.normalisedSmiles("O", "O", "CCC=O"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        //without the filter, both patterns match independently, the remainder are the atoms of no match
        tmpFragmenter.setFilterOverlappingMatchesSetting(false);
        Assertions.assertEquals(this.normalisedSmiles("O", "O", "OC=O", "CC"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        //matches of the same pattern
        IAtomContainer tmpButane = this.smilesParser.parseSmiles("CCCC");
        tmpFragmenter.setPatternsSetting("CC");
        Assertions.assertEquals(this.normalisedSmiles("CC", "CC", "CC"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpButane)));
        tmpFragmenter.setFilterOverlappingMatchesSetting(true);
        Assertions.assertEquals(this.normalisedSmiles("CC", "CC"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpButane)));
    }
    //
    /**
     * Tests that the result does not depend on the input atom order, for molecules and patterns where a greedy
     * selection of overlapping matches in input order would give different results, also for stereochemistry and
     * aromatic systems.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void atomOrderIndependenceTest() throws Exception {
        String[] tmpSmilesCodes = new String[] {
                //isopentane: "CC" can be taken as C1-C2 + C3-C4 or as C2-C3 + two single carbons
                "CC(C)CC",
                "CCC(C)C(C)CCC(C)(C)C",
                //naphthalene and morphine
                "c1ccc2ccccc2c1",
                "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
                //a sugar with stereocentres
                "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                //E/Z double bonds
                "C/C=C/C=C\\CC(=O)OC/C=C/C"
        };
        String tmpPatterns = String.join("\n", "c1ccccc1 benzene ring", "[OX2H] hydroxyl", "C=C alkene", "CC", "[#6]~[#8]");
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting(tmpPatterns);
        Random tmpRandom = new Random(42);
        for (boolean tmpIsFiltered : new boolean[] {true, false}) {
            tmpFragmenter.setFilterOverlappingMatchesSetting(tmpIsFiltered);
            for (String tmpSmiles : tmpSmilesCodes) {
                IAtomContainer tmpMolecule = this.smilesParser.parseSmiles(tmpSmiles);
                List<String> tmpExpected = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule));
                for (int i = 0; i < 20; i++) {
                    IAtomContainer tmpPermutedMolecule = SmartsFragmenterTest.permuteAtoms(tmpMolecule, tmpRandom);
                    Assertions.assertEquals(tmpExpected, this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpPermutedMolecule)),
                            "Permutation " + i + " of " + tmpSmiles + ", overlap filter " + tmpIsFiltered);
                }
            }
        }
        //the case that motivates the canonical order
        tmpFragmenter.setPatternsSetting("CC");
        tmpFragmenter.setFilterOverlappingMatchesSetting(true);
        List<String> tmpIsopentaneResult = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CC(C)CC")));
        Assertions.assertTrue(tmpIsopentaneResult.equals(this.normalisedSmiles("CC", "CC", "C"))
                || tmpIsopentaneResult.equals(this.normalisedSmiles("CC", "C", "C")), tmpIsopentaneResult.toString());
    }
    //
    /**
     * Tests the fragment saturation setting: without saturation, the atoms at cut bonds keep their open valences.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void saturationTest() throws Exception {
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("CC(=O)O");
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)O");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        List<IAtomContainer> tmpSaturatedFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertEquals(1, tmpSaturatedFragments.size());
        Assertions.assertEquals(2, this.getTotalHydrogenCount(tmpSaturatedFragments.getFirst()));
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        List<IAtomContainer> tmpUnsaturatedFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertEquals(1, this.getTotalHydrogenCount(tmpUnsaturatedFragments.getFirst()));
        //aromatic ring of biphenyl: bond orders and aromaticity are the same with and without saturation
        IAtomContainer tmpBiphenyl = this.smilesParser.parseSmiles("c1ccccc1-c1ccccc1");
        tmpFragmenter.setPatternsSetting("c1ccccc1");
        List<IAtomContainer> tmpPhenylFragments = tmpFragmenter.fragmentMolecule(tmpBiphenyl);
        Assertions.assertEquals(2, tmpPhenylFragments.size());
        for (IAtomContainer tmpFragment : tmpPhenylFragments) {
            Assertions.assertEquals(5, this.getTotalHydrogenCount(tmpFragment));
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Assertions.assertTrue(tmpAtom.isAromatic());
            }
            Assertions.assertNotNull(this.createSmiles(tmpFragment));
        }
        Assertions.assertEquals(this.createSmiles(tmpPhenylFragments.get(0)), this.createSmiles(tmpPhenylFragments.get(1)));
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.HYDROGEN_SATURATION);
        Assertions.assertEquals(this.normalisedSmiles("c1ccccc1", "c1ccccc1"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpBiphenyl)));
    }
    //
    /**
     * Tests patterns with unconnected parts, component-level grouping, and recursive SMARTS.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void disconnectedAndRecursiveSmartsTest() throws Exception {
        //ethylene glycol
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("OCCO");
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[OX2H].[OX2H] two hydroxyl groups");
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertEquals(this.normalisedSmiles("O.O", "CC"), this.fragmentSmiles(tmpFragments));
        for (IAtomContainer tmpFragment : tmpFragments) {
            if ("two hydroxyl groups".equals(tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY))) {
                //one fragment of two unconnected parts
                Assertions.assertEquals(2, tmpFragment.getAtomCount());
                Assertions.assertFalse(ConnectivityChecker.isConnected(tmpFragment));
            }
        }
        //component-level grouping: both oxygen atoms have to be in different components
        tmpFragmenter.setPatternsSetting("([OX2H]).([OX2H])");
        Assertions.assertEquals(this.normalisedSmiles("OCCO"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        Assertions.assertEquals(this.normalisedSmiles("O.O", "C", "C"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CO.OC"))));
        //recursive SMARTS: carbonyl carbon atoms of carboxylic acids only
        tmpFragmenter.setPatternsSetting("[C;$(C(=O)[OH])] acid carbon");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Assertions.assertEquals(this.normalisedSmiles("C"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CCC(=O)O"))));
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CCC(=O)C")).isEmpty());
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CCC(=O)OC")).isEmpty());
    }
    //
    /**
     * Tests that stereo elements are kept if they are completely contained in a fragment, with the correct
     * configuration, and dropped otherwise.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void stereoTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        //L-alanine, the alpha carbon atom and all its neighbours are matched
        tmpFragmenter.setPatternsSetting("[NX3][CX4]([CH3])[CX3]=O");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("N[C@@H](C)C(=O)O"));
        Assertions.assertEquals(1, tmpFragments.size());
        Assertions.assertTrue(tmpFragments.getFirst().stereoElements().iterator().hasNext());
        Assertions.assertEquals(this.normalisedSmiles("N[C@@H](C)C=O"), this.fragmentSmiles(tmpFragments));
        Assertions.assertNotEquals(this.normalisedSmiles("N[C@H](C)C=O"), this.fragmentSmiles(tmpFragments));
        //E double bond in the remainder whose neighbouring bond to the carboxyl group is cut
        tmpFragmenter.setPatternsSetting("C(=O)[OH]");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES_AND_REMAINDER);
        List<String> tmpSmiles = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C/CC(=O)O")));
        Assertions.assertEquals(this.normalisedSmiles("C/C=C/C", "OC=O"), tmpSmiles);
        Assertions.assertNotEquals(this.normalisedSmiles("C/C=C\\C", "OC=O"), tmpSmiles);
        tmpSmiles = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C\\CC(=O)O")));
        Assertions.assertEquals(this.normalisedSmiles("C/C=C\\C", "OC=O"), tmpSmiles);
        //stereocentre at the border of a fragment: one of its neighbours is cut, so the stereo element is dropped
        tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(CC)C(=O)O"));
        for (IAtomContainer tmpFragment : tmpFragments) {
            for (IStereoElement tmpStereoElement : tmpFragment.stereoElements()) {
                Assertions.fail("Unexpected stereo element " + tmpStereoElement);
            }
        }
        Assertions.assertEquals(this.normalisedSmiles("CC(O)CC", "OC=O"), this.fragmentSmiles(tmpFragments));
    }
    //
    /**
     * Tests that aromatic rings are matched identically in molecules given in Kekulé form, with aromaticity flags, or
     * with unset bond orders, and that partially matched aromatic systems get a Kekulé structure.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void aromaticityTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("c1ccccc1 benzene ring\n[OX2H] hydroxyl");
        List<String> tmpExpected = this.normalisedSmiles("c1ccccc1", "O", "CC");
        //Kekulé form without aromaticity flags
        SmilesParser tmpKekuleParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpKekule = tmpKekuleParser.parseSmiles("OC1=CC=C(CC)C=C1");
        for (IAtom tmpAtom : tmpKekule.atoms()) {
            Assertions.assertFalse(tmpAtom.isAromatic());
        }
        Assertions.assertEquals(tmpExpected, this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpKekule)));
        //aromatic SMILES, Kekulé structure and aromaticity flags
        Assertions.assertEquals(tmpExpected, this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("Oc1ccc(CC)cc1"))));
        //aromatic SMILES with unset bond orders
        SmilesParser tmpNonKekulisingParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpNonKekulisingParser.kekulise(false);
        IAtomContainer tmpUnsetBondOrders = tmpNonKekulisingParser.parseSmiles("Oc1ccc(CC)cc1");
        Assertions.assertTrue(tmpFragmenter.canBeFragmented(tmpUnsetBondOrders));
        Assertions.assertEquals(tmpExpected, this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpUnsetBondOrders)));
        //the input with unset bond orders is not altered
        for (IBond tmpBond : tmpUnsetBondOrders.bonds()) {
            if (tmpBond.isAromatic()) {
                Assertions.assertEquals(IBond.Order.UNSET, tmpBond.getOrder());
            }
        }
        //fused ring system: one benzene ring is matched, the rest of the other ring gets a Kekulé structure
        Assertions.assertEquals(this.normalisedSmiles("c1ccccc1", "C=CC=C"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("c1ccc2ccccc2c1"))));
        //pyrrole ring of indole stays aromatic
        tmpFragmenter.setPatternsSetting("c1cc[nH]c1");
        Assertions.assertEquals(this.normalisedSmiles("c1cc[nH]c1", "C=CC=C"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("c1ccc2[nH]ccc2c1"))));
        //aliphatic SMARTS atoms do not match aromatic atoms
        tmpFragmenter.setPatternsSetting("C=C");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(tmpKekule).isEmpty());
    }
    //
    /**
     * Tests the filtering of molecules: empty molecules and molecules with unset bond orders that cannot be resolved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void filterTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(null));
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
        Assertions.assertFalse(tmpFragmenter.shouldBeFiltered(this.smilesParser.parseSmiles("CCO")));
        Assertions.assertFalse(tmpFragmenter.shouldBePreprocessed(this.smilesParser.parseSmiles("CCO")));
        SmilesParser tmpNonKekulisingParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpNonKekulisingParser.kekulise(false);
        //pyrrole written without the hydrogen atom at the nitrogen atom has no Kekulé structure
        IAtomContainer tmpNotKekulizable = tmpNonKekulisingParser.parseSmiles("c1ccnc1");
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpNotKekulizable));
        Assertions.assertFalse(tmpFragmenter.canBeFragmented(tmpNotKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.fragmentMolecule(tmpNotKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.applyPreprocessing(tmpNotKekulizable));
        //unset order of a non-aromatic bond
        IAtomContainer tmpUnsetBond = this.smilesParser.parseSmiles("CCO");
        tmpUnsetBond.getBond(0).setOrder(IBond.Order.UNSET);
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpUnsetBond));
    }
    //
    /**
     * Tests the behaviour for an empty pattern list: no exception, the molecule is returned as remainder or no
     * fragments are returned.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void emptyPatternListTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("# only a comment\n\n");
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("C[C@@H](N)C(=O)O");
        Assertions.assertTrue(tmpFragmenter.canBeFragmented(tmpMolecule));
        Assertions.assertEquals(this.normalisedSmiles("C[C@@H](N)C(=O)O"), this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule)));
        //unconnected parts are split
        Assertions.assertEquals(this.normalisedSmiles("CC(=O)[O-]", "[Na+]"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("CC(=O)[O-].[Na+]"))));
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(tmpMolecule).isEmpty());
    }
    //
    /**
     * Tests that invalid SMARTS patterns are rejected without changing the setting.
     */
    @Test
    public void invalidSmartsTest() {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)O acid");
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setPatternsSetting("C(=O)O acid\nC(C broken"));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setPatternsSetting("[C"));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setPatternsSetting(null));
        Assertions.assertEquals("C(=O)O acid", tmpFragmenter.getPatternsSetting());
    }
    //
    /**
     * Tests that a copy has the same settings, produces the same fragments, and is independent of the original, and
     * that changed patterns are compiled again.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void copyTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[OX2H] hydroxyl\nC=O carbonyl");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        tmpFragmenter.setFilterOverlappingMatchesSetting(false);
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles("OCC(=O)CO");
        List<String> tmpOriginalResult = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule));
        SmartsFragmenter tmpCopy = (SmartsFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(tmpFragmenter.getPatternsSetting(), tmpCopy.getPatternsSetting());
        Assertions.assertEquals(tmpFragmenter.getReturnedFragmentsSetting(), tmpCopy.getReturnedFragmentsSetting());
        Assertions.assertEquals(tmpFragmenter.getFilterOverlappingMatchesSetting(), tmpCopy.getFilterOverlappingMatchesSetting());
        Assertions.assertEquals(tmpFragmenter.getFragmentSaturationSetting(), tmpCopy.getFragmentSaturationSetting());
        Assertions.assertEquals(tmpOriginalResult, this.fragmentSmiles(tmpCopy.fragmentMolecule(tmpMolecule)));
        //changing the original does not change the copy, and the original compiles its changed patterns
        tmpFragmenter.setPatternsSetting("C=O carbonyl");
        Assertions.assertEquals(1, tmpFragmenter.fragmentMolecule(tmpMolecule).size());
        Assertions.assertEquals("[OX2H] hydroxyl\nC=O carbonyl", tmpCopy.getPatternsSetting());
        Assertions.assertEquals(tmpOriginalResult, this.fragmentSmiles(tmpCopy.fragmentMolecule(tmpMolecule)));
    }
    //
    /**
     * Tests that the default settings are restored.
     */
    @Test
    public void restoreDefaultSettingsTest() {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[OX2H] hydroxyl");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        tmpFragmenter.setFilterOverlappingMatchesSetting(false);
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(SmartsFragmenter.PATTERNS_SETTING_DEFAULT, tmpFragmenter.getPatternsSetting());
        Assertions.assertEquals(SmartsFragmenter.RETURNED_FRAGMENTS_OPTION_DEFAULT, tmpFragmenter.getReturnedFragmentsSetting());
        Assertions.assertEquals(SmartsFragmenter.FILTER_OVERLAPPING_MATCHES_SETTING_DEFAULT, tmpFragmenter.getFilterOverlappingMatchesSetting());
        Assertions.assertEquals(IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT, tmpFragmenter.getFragmentSaturationSetting());
    }
    //
    /**
     * Tests that the fragments can be written as SMILES codes, parsed again, and passed to the fragmenter again, as in
     * a pipeline.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void fragmentRoundTripTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("c1ccccc1\n[OX2H]\nC(=O)O");
        //CNP0151033, see ErtlFunctionalGroupsFinderFragmenterTest
        IAtomContainer tmpMolecule = this.smilesParser.parseSmiles(
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4");
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
        Assertions.assertFalse(tmpFragments.isEmpty());
        for (IAtomContainer tmpFragment : tmpFragments) {
            String tmpSmiles = ChemUtil.createUniqueSmiles(tmpFragment, true);
            Assertions.assertNotNull(tmpSmiles);
            IAtomContainer tmpReparsedFragment = ChemUtil.parseSmilesToAtomContainer(tmpSmiles);
            Assertions.assertTrue(tmpFragmenter.canBeFragmented(tmpReparsedFragment), tmpSmiles);
            Assertions.assertDoesNotThrow(() -> tmpFragmenter.fragmentMolecule(tmpReparsedFragment));
        }
    }
    //
    //<editor-fold desc="Private methods">
    /**
     * Creates the unique SMILES code with stereochemistry that MORTAR uses to identify fragments.
     *
     * @param aMolecule the molecule
     * @return unique SMILES code
     */
    private String createSmiles(IAtomContainer aMolecule) {
        return ChemUtil.createUniqueSmiles(aMolecule, true);
    }
    //
    /**
     * Returns the sorted unique SMILES codes of the given fragments.
     *
     * @param aFragmentList fragments
     * @return sorted SMILES codes
     */
    private List<String> fragmentSmiles(List<IAtomContainer> aFragmentList) {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(this.createSmiles(tmpFragment));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Parses the given SMILES codes, applies the Daylight aromaticity model as the fragmenter does, and returns their
     * sorted unique SMILES codes.
     *
     * @param aSmilesCodes SMILES codes
     * @return sorted unique SMILES codes
     * @throws Exception if a SMILES code cannot be parsed
     */
    private List<String> normalisedSmiles(String... aSmilesCodes) throws Exception {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesCodes.length);
        for (String tmpSmiles : aSmilesCodes) {
            IAtomContainer tmpMolecule = this.smilesParser.parseSmiles(tmpSmiles);
            SmartsPattern.prepare(tmpMolecule);
            tmpSmilesList.add(this.createSmiles(tmpMolecule));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the number of implicit hydrogen atoms of the given molecule.
     *
     * @param aMolecule the molecule
     * @return sum of the implicit hydrogen counts
     */
    private int getTotalHydrogenCount(IAtomContainer aMolecule) {
        int tmpCount = 0;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpCount += tmpAtom.getImplicitHydrogenCount();
        }
        return tmpCount;
    }
    //
    /**
     * Returns a copy of the given molecule with randomly permuted atoms and bonds, including its stereo elements.
     *
     * @param aMolecule the molecule
     * @param aRandom random number generator
     * @return permuted copy
     * @throws Exception if an atom cannot be cloned
     */
    private static IAtomContainer permuteAtoms(IAtomContainer aMolecule, Random aRandom) throws Exception {
        List<Integer> tmpAtomOrder = new ArrayList<>();
        for (int i = 0; i < aMolecule.getAtomCount(); i++) {
            tmpAtomOrder.add(i);
        }
        Collections.shuffle(tmpAtomOrder, aRandom);
        IAtomContainer tmpPermuted = aMolecule.getBuilder().newAtomContainer();
        Map<IAtom, IAtom> tmpAtomMap = new HashMap<>();
        for (int tmpIndex : tmpAtomOrder) {
            tmpPermuted.addAtom(aMolecule.getAtom(tmpIndex).clone());
            tmpAtomMap.put(aMolecule.getAtom(tmpIndex), tmpPermuted.getAtom(tmpPermuted.getAtomCount() - 1));
        }
        List<IBond> tmpBonds = new ArrayList<>();
        aMolecule.bonds().forEach(tmpBonds::add);
        Collections.shuffle(tmpBonds, aRandom);
        Map<IBond, IBond> tmpBondMap = new HashMap<>();
        for (IBond tmpBond : tmpBonds) {
            IBond tmpNewBond = aMolecule.getBuilder().newInstance(IBond.class,
                    tmpAtomMap.get(tmpBond.getBegin()), tmpAtomMap.get(tmpBond.getEnd()), tmpBond.getOrder());
            tmpNewBond.setIsAromatic(tmpBond.isAromatic());
            tmpPermuted.addBond(tmpNewBond);
            tmpBondMap.put(tmpBond, tmpPermuted.getBond(tmpPermuted.getBondCount() - 1));
        }
        for (IStereoElement tmpStereoElement : aMolecule.stereoElements()) {
            tmpPermuted.addStereoElement(tmpStereoElement.map(tmpAtomMap, tmpBondMap));
        }
        //sanity check
        Assertions.assertEquals(new SmilesGenerator(SmiFlavor.Absolute).create(aMolecule),
                new SmilesGenerator(SmiFlavor.Absolute).create(tmpPermuted));
        return tmpPermuted;
    }
    //</editor-fold>
}
