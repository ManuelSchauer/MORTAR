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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

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
        //stereocentre at the border of a fragment that lost one neighbour, see borderStereocentreTest()
        Assertions.assertEquals(this.normalisedSmiles("C[C@](O)(CC)[H]", "OC=O"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(CC)C(=O)O"))));
    }
    //
    /**
     * Tests that a stereocentre at the border of a fragment keeps its configuration if it lost one neighbour and had no
     * hydrogen atom, with the new hydrogen atom at the position of the cut neighbour, and that both enantiomers give
     * different fragments. The expected fragment is the parent SMILES code with the cut neighbour replaced by an
     * explicit hydrogen atom at the same position, i.e. the geometry is compared, not only the presence of a descriptor.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void borderStereocentreTest() throws Exception {
        //cross-check of the expected SMILES codes: an explicit hydrogen atom at the position of the third ligand
        Assertions.assertEquals(this.normalisedSmiles("C[C@H](O)CC"), this.normalisedSmiles("C[C@](O)(CC)[H]"));
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid");
        //2-hydroxy-2-methylbutanoic acid, the remainder is butan-2-ol
        List<String> tmpRemainderOfS = this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(CC)C(=O)O")));
        List<String> tmpRemainderOfR = this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@@](O)(CC)C(=O)O")));
        Assertions.assertEquals(this.normalisedSmiles("C[C@](O)(CC)[H]"), tmpRemainderOfS);
        Assertions.assertEquals(this.normalisedSmiles("C[C@@](O)(CC)[H]"), tmpRemainderOfR);
        Assertions.assertNotEquals(tmpRemainderOfS, tmpRemainderOfR);
        Assertions.assertNotEquals(this.normalisedSmiles("CC(O)CC"), tmpRemainderOfS);
        //the cut neighbour written at another position of the parent SMILES code
        Assertions.assertEquals(this.normalisedSmiles("OC(=O)[C@@](O)(C)CC".replace("OC(=O)", "[H]")),
                this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("OC(=O)[C@@](O)(C)CC"))));
        //without saturation, the centre has an open valence and loses its configuration
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(CC)C(=O)O"))) {
            Assertions.assertFalse(tmpFragment.stereoElements().iterator().hasNext());
        }
        //a completely contained stereocentre is kept without saturation
        tmpFragmenter.setPatternsSetting("O[CX4](C)CC");
        tmpFragmenter.setReturnedFragmentsSetting(SmartsFragmenter.SmartsFragmenterReturnedFragmentsOption.MATCHED_SUBSTRUCTURES);
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@H](O)CCC(=O)O"));
        Assertions.assertEquals(1, tmpFragments.size());
        Assertions.assertTrue(tmpFragments.getFirst().stereoElements().iterator().hasNext());
        this.assertRoundTrip(tmpFragments);
    }
    //
    /**
     * Tests that stereo elements that are no longer stereogenic in the fragment are dropped, so that both enantiomers
     * of the parent give the same fragment SMILES code without descriptor, while pseudo-asymmetric centres are kept.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void lostStereogenicityTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid");
        String[][] tmpEnantiomersAndRemainder = new String[][] {
                //citramalic acid: border centre whose new hydrogen atom comes with two methyl groups
                {"C[C@](O)(CC(=O)O)C(=O)O", "C[C@@](O)(CC(=O)O)C(=O)O", "CC(C)O"},
                //completely contained centre that gets two methyl groups
                {"C[C@](O)(CC)CC(=O)O", "C[C@@](O)(CC)CC(=O)O", "CCC(C)(C)O"},
                //completely contained centre that gets two ethyl groups (constitutionally equivalent, not terminal)
                {"CC[C@](C)(O)CCC(=O)O", "CC[C@@](C)(O)CCC(=O)O", "CCC(C)(O)CC"},
                //border centre that gets a hydrogen atom and two ethyl groups
                {"CC[C@](O)(CCC(=O)O)C(=O)O", "CC[C@@](O)(CCC(=O)O)C(=O)O", "CCC(O)CC"},
                //double bond whose end gets two methyl groups
                {"C/C(CC(=O)O)=C/C", "C/C(CC(=O)O)=C\\C", "CC(C)=CC"}
        };
        for (String[] tmpCase : tmpEnantiomersAndRemainder) {
            List<IAtomContainer> tmpFragmentsOfFirst = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles(tmpCase[0]));
            List<IAtomContainer> tmpFragmentsOfSecond = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles(tmpCase[1]));
            Assertions.assertEquals(this.remainderSmiles(tmpFragmentsOfFirst), this.remainderSmiles(tmpFragmentsOfSecond), tmpCase[0]);
            Assertions.assertEquals(this.normalisedSmiles(tmpCase[2]), this.remainderSmiles(tmpFragmentsOfFirst), tmpCase[0]);
            for (IAtomContainer tmpFragment : tmpFragmentsOfFirst) {
                Assertions.assertFalse(tmpFragment.stereoElements().iterator().hasNext(), tmpCase[0]);
            }
        }
        //pentane-2,3,4-triol (no match, whole molecule): the pseudo-asymmetric middle carbon atom is kept, the two
        // diastereomers that differ only in its configuration give different fragments
        List<String> tmpFirst = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@H](O)[C@H](O)[C@H](O)C")));
        List<String> tmpSecond = this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@H](O)[C@@H](O)[C@H](O)C")));
        Assertions.assertNotEquals(tmpFirst, tmpSecond);
        Assertions.assertEquals(this.normalisedSmiles("C[C@H](O)[C@H](O)[C@H](O)C"), tmpFirst);
        Assertions.assertEquals(this.normalisedSmiles("C[C@H](O)[C@@H](O)[C@H](O)C"), tmpSecond);
        //the same with a cut: both outer carbon atoms are border centres, the middle one stays pseudo-asymmetric
        tmpFirst = this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(C(=O)O)[C@H](O)[C@](O)(C(=O)O)C")));
        tmpSecond = this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C[C@](O)(C(=O)O)[C@@H](O)[C@](O)(C(=O)O)C")));
        Assertions.assertNotEquals(tmpFirst, tmpSecond);
        Assertions.assertEquals(this.normalisedSmiles("C[C@](O)([H])[C@H](O)[C@](O)([H])C"), tmpFirst);
        Assertions.assertEquals(this.normalisedSmiles("C[C@](O)([H])[C@@H](O)[C@](O)([H])C"), tmpSecond);
    }
    //
    /**
     * Tests stereocentres with two cut neighbours and stereocentres that already had a hydrogen atom: they end up with
     * at least two hydrogen atoms, so their elements are dropped without exception.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void droppedBorderStereocentreTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid\n[NX3;H2] amine");
        //isovaline: two neighbours of the centre are cut
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("N[C@@](CC)(C)C(=O)O"));
        Assertions.assertEquals(this.normalisedSmiles("CCCC", "N", "OC=O"), this.fragmentSmiles(tmpFragments));
        //alanine: the centre had a hydrogen atom
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid");
        tmpFragments.addAll(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("N[C@@H](C)C(=O)O")));
        Assertions.assertEquals(this.normalisedSmiles("CCN", "OC=O"),
                this.fragmentSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("N[C@@H](C)C(=O)O"))));
        for (IAtomContainer tmpFragment : tmpFragments) {
            Assertions.assertFalse(tmpFragment.stereoElements().iterator().hasNext());
        }
    }
    //
    /**
     * Tests stereo double bonds whose referenced substituent was cut: the configuration is kept via the other
     * substituent on this end, or the element is dropped if this end carries only hydrogen atoms afterwards.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void borderDoubleBondTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid");
        //2-methylbut-2-enoic acids, the parsed element refers to the carboxyl group, which is cut
        for (String tmpSmiles : new String[] {"C/C=C(/C(=O)O)C", "C/C=C(\\C(=O)O)C"}) {
            IAtomContainer tmpMolecule = this.smilesParser.parseSmiles(tmpSmiles);
            boolean tmpIsCarboxylReferenced = false;
            for (IStereoElement tmpElement : tmpMolecule.stereoElements()) {
                for (IBond tmpCarrier : ((IDoubleBondStereochemistry) tmpElement).getBonds()) {
                    tmpIsCarboxylReferenced = tmpIsCarboxylReferenced || tmpCarrier.contains(tmpMolecule.getAtom(3));
                }
            }
            Assertions.assertTrue(tmpIsCarboxylReferenced, "precondition: the carboxyl carbon atom is referenced");
            List<String> tmpRemainder = this.remainderSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule));
            //the parent SMILES code with the carboxyl group replaced by a hydrogen atom at the same position
            Assertions.assertEquals(this.normalisedSmiles(tmpSmiles.replace("C(=O)O", "[H]")), tmpRemainder, tmpSmiles);
        }
        //cross-check: the carboxyl group trans to the methyl group of the other end means cis-but-2-ene
        Assertions.assertEquals(this.normalisedSmiles("C/C=C\\C"),
                this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C(/C(=O)O)C"))));
        Assertions.assertEquals(this.normalisedSmiles("C/C=C/C"),
                this.remainderSmiles(tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C(\\C(=O)O)C"))));
        //crotonic acid: the end becomes =CH2, the element is dropped
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C/C(=O)O"));
        Assertions.assertEquals(this.normalisedSmiles("C=CC", "OC=O"), this.fragmentSmiles(tmpFragments));
        for (IAtomContainer tmpFragment : tmpFragments) {
            Assertions.assertFalse(tmpFragment.stereoElements().iterator().hasNext());
        }
        //without saturation, the element at the border is dropped
        tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles("C/C=C(/C(=O)O)C"))) {
            Assertions.assertFalse(tmpFragment.stereoElements().iterator().hasNext());
        }
    }
    //
    /**
     * Tests a stereocentre at the border of a matched fragment, with the overlap filter turned on and off.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void borderStereoInMatchedFragmentTest() throws Exception {
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("[CH3][CX4]([OH])[CH2][CH3] sec-butanol\n[CX4]C(=O)[OH] alpha carbon and acid");
        for (boolean tmpIsFiltered : new boolean[] {true, false}) {
            tmpFragmenter.setFilterOverlappingMatchesSetting(tmpIsFiltered);
            for (String tmpSmiles : new String[] {"C[C@](O)(CC)C(=O)O", "C[C@@](O)(CC)C(=O)O"}) {
                List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(this.smilesParser.parseSmiles(tmpSmiles));
                List<String> tmpButanol = new ArrayList<>();
                for (IAtomContainer tmpFragment : tmpFragments) {
                    if ("sec-butanol".equals(tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY))) {
                        tmpButanol.add(this.createSmiles(tmpFragment));
                    }
                }
                Assertions.assertEquals(this.normalisedSmiles(tmpSmiles.replace("C(=O)O", "[H]")), tmpButanol);
                //with the filter, the second pattern cannot match because its alpha carbon atom is taken, the rest of
                // the acid is the remainder; without the filter, the alpha carbon atom of the acid loses three neighbours
                Assertions.assertEquals(this.normalisedSmiles(tmpSmiles.replace("C(=O)O", "[H]"), tmpIsFiltered ? "OC=O" : "CC(=O)O"),
                        this.fragmentSmiles(tmpFragments), tmpSmiles + ", overlap filter " + tmpIsFiltered);
                this.assertRoundTrip(tmpFragments);
            }
        }
    }
    //
    /**
     * Tests that the configurations at the borders of fragments do not depend on the input atom order and that the
     * fragment SMILES codes can be parsed again to the same SMILES codes, as in a pipeline.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void borderStereoAtomOrderIndependenceTest() throws Exception {
        String[] tmpSmilesCodes = new String[] {
                "C[C@](O)(CC)C(=O)O",
                "C[C@](O)(CC(=O)O)C(=O)O",
                "N[C@@](CC)(C)C(=O)O",
                "C/C=C(/C(=O)O)C",
                "CC/C(C(=O)O)=C(/C)CC(=O)O",
                "C[C@](O)(C(=O)O)[C@H](O)[C@](O)(C(=O)O)C",
                //morphine and a glucuronic acid
                "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
                "O[C@@H]1O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]1O"
        };
        SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
        tmpFragmenter.setPatternsSetting("C(=O)[OH] carboxylic acid\n[NX3;H2] amine\nc1ccccc1 benzene ring\n[CX4]O[CX4] ether");
        Random tmpRandom = new Random(7);
        for (boolean tmpIsFiltered : new boolean[] {true, false}) {
            tmpFragmenter.setFilterOverlappingMatchesSetting(tmpIsFiltered);
            for (String tmpSmiles : tmpSmilesCodes) {
                IAtomContainer tmpMolecule = this.smilesParser.parseSmiles(tmpSmiles);
                List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                this.assertRoundTrip(tmpFragments);
                List<String> tmpExpected = this.fragmentSmiles(tmpFragments);
                for (int i = 0; i < 20; i++) {
                    IAtomContainer tmpPermutedMolecule = SmartsFragmenterTest.permuteAtoms(tmpMolecule, tmpRandom);
                    Assertions.assertEquals(tmpExpected, this.fragmentSmiles(tmpFragmenter.fragmentMolecule(tmpPermutedMolecule)),
                            "Permutation " + i + " of " + tmpSmiles + ", overlap filter " + tmpIsFiltered);
                }
            }
        }
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
     * Returns the sorted unique SMILES codes of the remainder fragments among the given fragments.
     *
     * @param aFragmentList fragments
     * @return sorted SMILES codes of the remainder fragments
     */
    private List<String> remainderSmiles(List<IAtomContainer> aFragmentList) {
        List<IAtomContainer> tmpRemainderFragments = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            if (SmartsFragmenter.FRAGMENT_CATEGORY_REMAINDER_VALUE.equals(
                    tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY))) {
                tmpRemainderFragments.add(tmpFragment);
            }
        }
        return this.fragmentSmiles(tmpRemainderFragments);
    }
    //
    /**
     * Asserts that the unique SMILES code of every given fragment is parsed again to an atom container with the same
     * unique SMILES code, as the fragment data model does it in MORTAR, e.g. for the next step of a pipeline.
     *
     * @param aFragmentList fragments
     * @throws Exception if a SMILES code cannot be parsed
     */
    private void assertRoundTrip(List<IAtomContainer> aFragmentList) throws Exception {
        for (IAtomContainer tmpFragment : aFragmentList) {
            String tmpSmiles = this.createSmiles(tmpFragment);
            Assertions.assertNotNull(tmpSmiles);
            IAtomContainer tmpParsedFragment;
            try {
                tmpParsedFragment = ChemUtil.parseSmilesToAtomContainer(tmpSmiles, true, true);
            } catch (CDKException anException) {
                tmpParsedFragment = ChemUtil.parseSmilesToAtomContainer(tmpSmiles, false, false);
            }
            Assertions.assertEquals(tmpSmiles, this.createSmiles(tmpParsedFragment), tmpSmiles);
        }
    }
    //
    /**
     * Parses the given SMILES codes, makes explicit hydrogen atoms implicit, applies the Daylight aromaticity model as
     * the fragmenter does, and returns their sorted unique SMILES codes.
     *
     * @param aSmilesCodes SMILES codes
     * @return sorted unique SMILES codes
     * @throws Exception if a SMILES code cannot be parsed
     */
    private List<String> normalisedSmiles(String... aSmilesCodes) throws Exception {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesCodes.length);
        for (String tmpSmiles : aSmilesCodes) {
            //explicit hydrogen atoms become implicit ones, the CDK adapts the stereo elements accordingly
            IAtomContainer tmpMolecule = AtomContainerManipulator.suppressHydrogens(this.smilesParser.parseSmiles(tmpSmiles));
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
