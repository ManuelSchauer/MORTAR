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

package de.unijena.cheminf.mortar.model.util;

import de.unijena.cheminf.mortar.model.data.FragmentDataModel;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.fragmentation.FragmentationTask;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.CDKCircularFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.SugarRemovalUtilityFragmenter;
import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.graph.AtomContainerAtomPermutor;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import io.github.dan2097.jnainchi.InchiFlag;
import io.github.dan2097.jnainchi.InchiStatus;

/**
 * Tests that the unique SMILES codes created by {@link ChemUtil#createUniqueSmiles(IAtomContainer, boolean, boolean)}
 * do not depend on the double bond placement of conjugated rings without aromaticity flags (Kekulé structures, as
 * imported from SD files or cut out of such molecules) or on the atom order, in both the aromatic and the Kekulé mode,
 * while aromaticity flags that are given are still encoded with aromatic symbols.
 *
 * @author Manuel Schauer
 */
class ChemUtilKekuleSmilesTest {
    /**
     * Sets the default locale to British English, required for the settings container of the importer.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Seed of the random atom permutations.
     */
    private static final long SEED = 42L;
    //
    /**
     * Different Kekulé notations (different double bond placement and atom order) of one compound, parsed without
     * aromaticity flags, must give one SMILES code in all four modes (with and without stereochemistry and
     * aromaticity), also for random atom permutations. The expected codes are those of the Kekulé mode, which the
     * CDK normalises.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testKekuleNotationsOfOneCompound() throws Exception {
        Map<String, String[]> tmpCompounds = new LinkedHashMap<>();
        tmpCompounds.put("C=1C=CC(=CC1)C", new String[] {"CC1=CC=CC=C1", "CC1C=CC=CC=1", "C1=CC=C(C)C=C1"}); //toluene
        tmpCompounds.put("N=1C=CC=CC1", new String[] {"N1=CC=CC=C1", "N=1C=CC=CC1", "C1=CC=NC=C1"}); //pyridine
        tmpCompounds.put("C=1C=CC=2C=CC=CC2C1", new String[] {"C1=CC=C2C=CC=CC2=C1", "C1=CC2=CC=CC=C2C=C1",
                "C=1C=CC=2C=CC=CC2C1"}); //naphthalene, all three Kekulé structures
        tmpCompounds.put("O=C(NC1=CC=C(O)C=C1)C", new String[] {"CC(=O)NC1=CC=C(O)C=C1", "CC(=O)NC=1C=CC(O)=CC1"}); //paracetamol
        tmpCompounds.put("O=C1NC(=O)N(C=C1C)CC=2C=CC=CC2", new String[] {"CC1=CN(CC2=CC=CC=C2)C(=O)NC1=O",
                "CC1=CN(CC2C=CC=CC=2)C(=O)NC1=O"}); //1-benzylthymine
        tmpCompounds.put("C=1C=CC(=C(C1)C)C", new String[] {"CC1=C(C)C=CC=C1", "CC1C(C)=CC=CC=1"}); //o-xylene
        tmpCompounds.put("C=1C=CC=2NC=CC2C1", new String[] {"C1=CC=C2NC=CC2=C1", "C1=CC2=C(N1)C=CC=C2"}); //indole
        tmpCompounds.put("C=1C=CNC1", new String[] {"N1C=CC=C1", "C1=CNC=C1"}); //pyrrole, one Kekulé structure
        tmpCompounds.put("C=1C=C[N+](=CC1)C", new String[] {"C[N+]1=CC=CC=C1", "C[N+]=1C=CC=CC1"}); //N-methylpyridinium
        tmpCompounds.put("C=1C=C[CH-]C1", new String[] {"[CH-]1C=CC=C1", "C1=C[CH-]C=C1"}); //cyclopentadienyl anion
        tmpCompounds.put("C1=CC=C[CH+]C=C1", new String[] {"[CH+]1C=CC=CC=C1", "C1=C[CH+]C=CC=C1"}); //tropylium
        tmpCompounds.put("C=1C=CC2=CC=CC2=CC1", new String[] {"C1=CC=C2C=CC=C2C=C1", "C1=CC2=CC=CC=CC2=C1"}); //azulene
        tmpCompounds.put("C=1C=C2C=CC=C2C1", new String[] {"C1=CC2=CC=CC2=C1", "C=1C=C2C=CC=C2C1"}); //pentalene
        tmpCompounds.put("C=1C=CC=CC=CC1", new String[] {"C1=CC=CC=CC=C1", "C=1C=CC=CC=CC1"}); //cyclooctatetraene
        tmpCompounds.put("C=1C=CC2=CC=CC=C(C1)C2", new String[] {"C12=CC=CC=C(C2)C=CC=C1",
                "C12C=CC=CC(C2)=CC=CC=1"}); //1,6-methano[10]annulene, bridged
        tmpCompounds.put("C=1C=CC=CC=CC=CC=CC=CC=CC1", new String[] {"C1=CC=CC=CC=CC=CC=CC=CC=C1",
                "C=1C=CC=CC=CC=CC=CC=CC=CC1"}); //[18]annulene without double bond configurations
        tmpCompounds.put("N=1C=CC=CC1.C=1C=CC(=CC1)C", new String[] {"CC1=CC=CC=C1.C1=CC=NC=C1",
                "CC1C=CC=CC=1.N=1C=CC=CC1"}); //disconnected
        tmpCompounds.put("*C=1C=CC=CC1", new String[] {"*C1=CC=CC=C1", "*C1C=CC=CC=1"}); //pseudo atom
        tmpCompounds.put("[H]C=1C=CC=CC1", new String[] {"[2H]C1=CC=CC=C1", "[H]C1C=CC=CC=1"}); //explicit hydrogen atom
        tmpCompounds.put("[CH2]C=1C=CC=CC1", new String[] {"[CH2]C1=CC=CC=C1", "[CH2]C1C=CC=CC=1"}); //radical
        Random tmpRandom = new Random(ChemUtilKekuleSmilesTest.SEED);
        for (Map.Entry<String, String[]> tmpEntry : tmpCompounds.entrySet()) {
            for (String tmpNotation : tmpEntry.getValue()) {
                IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpNotation, false, false);
                for (boolean tmpIsStereo : new boolean[] {true, false}) {
                    for (boolean tmpIsAromatic : new boolean[] {true, false}) {
                        Assertions.assertEquals(tmpEntry.getKey(), ChemUtil.createUniqueSmiles(tmpMolecule, tmpIsStereo, tmpIsAromatic),
                                tmpNotation + " stereo " + tmpIsStereo + " aromatic " + tmpIsAromatic);
                        for (int i = 0; i < 5; i++) {
                            IAtomContainer tmpPermutation = ChemUtilKekuleSmilesTest.permuteAtoms(tmpMolecule, tmpRandom);
                            Assertions.assertEquals(tmpEntry.getKey(), ChemUtil.createUniqueSmiles(tmpPermutation, tmpIsStereo, tmpIsAromatic),
                                    tmpNotation + " permuted, stereo " + tmpIsStereo + " aromatic " + tmpIsAromatic);
                        }
                    }
                }
            }
        }
    }
    //
    /**
     * Stereochemistry is encoded as before: two Kekulé notations of one stereoisomer give one stereo SMILES code,
     * the other stereoisomer gives a different one.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testStereoChemistryIsKept() throws Exception {
        String tmpS1 = this.createSmiles("C[C@H](N)C1=CC=CC=C1", true, true);
        Assertions.assertEquals("C[C@@H](C1=CC=CC=C1)N", tmpS1);
        Assertions.assertEquals(tmpS1, this.createSmiles("C[C@H](N)C1C=CC=CC=1", true, true));
        Assertions.assertEquals(tmpS1, this.createSmiles("C[C@H](N)c1ccccc1", true, false));
        Assertions.assertNotEquals(tmpS1, this.createSmiles("C[C@@H](N)C1=CC=CC=C1", true, true));
        String tmpE = this.createSmiles("C/C=C/C1=CC=CC=C1", true, true);
        Assertions.assertEquals(tmpE, this.createSmiles("C/C=C/C1C=CC=CC=1", true, true));
        Assertions.assertNotEquals(tmpE, this.createSmiles("C/C=C\\C1=CC=CC=C1", true, true));
        //pseudo atom in a stereo SMILES code (canonised with InChI numbers, pseudo atoms are given as Rf)
        Assertions.assertEquals(this.createSmiles("*[C@H](N)C1=CC=CC=C1", true, true),
                this.createSmiles("*[C@H](N)C1C=CC=CC=1", true, true));
    }
    //
    /**
     * Atoms and bonds with aromaticity flags are still written with aromatic symbols, and the distinction between
     * aromatic and aliphatic fragments of the Ertl fragmenter is kept.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testAromaticityFlagsAreKept() throws Exception {
        for (boolean tmpIsKekulized : new boolean[] {true, false}) {
            Assertions.assertEquals("c1ccccc1", ChemUtil.createUniqueSmiles(
                    ChemUtil.parseSmilesToAtomContainer("c1ccccc1", tmpIsKekulized, false), true, true));
            Assertions.assertEquals("c1ccc(cc1)C", ChemUtil.createUniqueSmiles(
                    ChemUtil.parseSmilesToAtomContainer("Cc1ccccc1", tmpIsKekulized, false), true, true));
            Assertions.assertEquals("C=1C=CC(=CC1)C", ChemUtil.createUniqueSmiles(
                    ChemUtil.parseSmilesToAtomContainer("Cc1ccccc1", tmpIsKekulized, false), true, false));
        }
        String[][] tmpPairs = {{"[*]O[*]", "*O*", "*o*"}, {"[*]N([*])[*]", "*N(*)*", "*n(*)*"}, {"C=C", "C=C", "cc"},
                {"O=C", "O=C", "O=[cH2]"}};
        String[] tmpAromaticNotations = {"[*]o[*]", "[*]n([*])[*]", "cc", "O=[cH2]"};
        for (int i = 0; i < tmpPairs.length; i++) {
            String tmpAliphatic = this.createSmiles(tmpPairs[i][0], true, true);
            String tmpAromatic = this.createSmiles(tmpAromaticNotations[i], true, true);
            Assertions.assertEquals(tmpPairs[i][1], tmpAliphatic);
            Assertions.assertEquals(tmpPairs[i][2], tmpAromatic);
            Assertions.assertNotEquals(tmpAliphatic, tmpAromatic);
        }
        //aromatic atoms outside of complete rings, as in the environments of the Ertl fragmenter
        Assertions.assertEquals("c-c1ccccc1", this.createSmiles("c-c1ccccc1", true, true));
        Assertions.assertEquals("cc(-c1ccccc1)c2ccccc2", this.createSmiles("cc(-c1ccccc1)c1ccccc1", true, true));
    }
    //
    /**
     * Molecules with aromaticity flags in some rings only: the rings without flags get the canonical double bond
     * placement, the flagged rings keep their aromatic symbols. Note: the two rings must not be constitutionally
     * equivalent here (e.g. biphenyl with one flagged ring), since the canonical atom order does not regard the
     * aromaticity flags and the order of two such rings would then depend on the atom order of the input.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testMixedAromaticityFlags() throws Exception {
        String tmpExpected = "OC=1C=CC(=CC1)c2ccccc2";
        for (boolean tmpIsKekulized : new boolean[] {true, false}) {
            for (String tmpNotation : new String[] {"c1ccccc1C1=CC=C(O)C=C1", "c1ccccc1C1C=CC(O)=CC=1"}) {
                IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpNotation, tmpIsKekulized, false);
                Assertions.assertEquals(tmpExpected, ChemUtil.createUniqueSmiles(tmpMolecule, true, true));
                Random tmpRandom = new Random(ChemUtilKekuleSmilesTest.SEED);
                for (int i = 0; i < 5; i++) {
                    Assertions.assertEquals(tmpExpected, ChemUtil.createUniqueSmiles(
                            ChemUtilKekuleSmilesTest.permuteAtoms(tmpMolecule, tmpRandom), true, true));
                }
            }
        }
        //2-naphthol with aromaticity flags on the atoms and bonds of the hydroxylated ring only, from different Kekulé
        // structures and atom orders, i.e. the fused bond is double in some of them and single in others
        IAtomContainer tmpNaphthol = ChemUtil.parseSmilesToAtomContainer("OC1=CC=C2C=CC=CC2=C1", false, false);
        Set<String> tmpResults = new TreeSet<>();
        Random tmpRandom = new Random(ChemUtilKekuleSmilesTest.SEED);
        for (int i = 0; i < 10; i++) {
            IAtomContainer tmpAlternative = ChemUtilKekuleSmilesTest.createAlternativeKekuleStructure(tmpNaphthol, tmpRandom);
            ChemUtilKekuleSmilesTest.flagHydroxylatedSixMemberedRing(tmpAlternative);
            tmpResults.add(ChemUtil.createUniqueSmiles(tmpAlternative, true, true));
        }
        Assertions.assertEquals(Set.of("Oc1ccc2C=CC=Cc2c1"), tmpResults);
    }
    //
    /**
     * Molecules whose SMILES code does not depend on the double bond placement get the same SMILES code as before,
     * e.g. rings with only one Kekulé structure, conjugated chains, and molecules with aromaticity flags only.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testUnambiguousSmilesCodesAreUnchanged() throws Exception {
        String[][] tmpCases = {
                {"O=C1C=CC(=O)C=C1", "O=C1C=CC(=O)C=C1"}, //p-benzoquinone, one Kekulé structure
                {"C=C1C=CC(=O)C=C1", "O=C1C=CC(=C)C=C1"}, //quinone methide
                {"C=CC=CC=C", "C=CC=CC=C"}, //hexatriene
                {"C1=CCCCC1", "C1=CCCCC1"}, //cyclohexene
                {"C1=CNC=C1", "C=1C=CNC1"}, //pyrrole
                {"O1C=CC=C1", "O1C=CC=C1"}, //furan
                {"c1ccccc1", "c1ccccc1"},
                {"Cc1ccccc1", "c1ccc(cc1)C"},
                {"c1ccc2ccccc2c1", "c1ccc2ccccc2c1"},
                {"c1ccncc1", "n1ccccc1"},
                {"C[C@@H](N)c1ccccc1", "C[C@H](c1ccccc1)N"},
                {"c1ccccc1C1=CCCC=C1", "c1ccc(cc1)C=2C=CCCC2"} //cyclohexadiene, not conjugated around the ring
        };
        for (String[] tmpCase : tmpCases) {
            Assertions.assertEquals(tmpCase[1], this.createSmiles(tmpCase[0], true, true), tmpCase[0]);
        }
    }
    //
    /**
     * Edge cases: aromatic bonds without bond order, molecules that cannot be kekulized, the existing fallbacks for
     * aromatic nitrogen atoms, and a large ring system with many Kekulé structures (fullerene C60).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testEdgeCases() throws Exception {
        //aromatic flags but no bond orders, as parsed without kekulization
        Assertions.assertEquals("OC=1C=CC(=CC1)c2ccccc2", this.createSmiles("c1ccccc1C1=CC=C(O)C=C1", true, true));
        Assertions.assertEquals("OC=1C=CC(=CC1)c2ccccc2", this.createSmiles("c1ccccc1C1C=CC(O)=CC=1", true, true));
        //cannot be kekulized, aromatic symbols are written as given; the Kekulé mode needs the aromatic nitrogen fix
        Assertions.assertEquals("n1cccc1", this.createSmiles("n1cccc1", true, true));
        Assertions.assertEquals("C=1C=CNC1", this.createSmiles("n1cccc1", true, false));
        //aromatic nitrogen fix with stereo (CHEBI:20794), see ChemUtilTest
        Assertions.assertNotNull(this.createSmiles("C[n+]1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c2nc(N)nc(=O)c21", true, true));
        //a bond without bond order and without aromaticity flag cannot be written, null is returned
        IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer("CC1=CC=CC=C1", false, false);
        tmpMolecule.getBond(0).setOrder(IBond.Order.UNSET);
        Assertions.assertNull(ChemUtil.createUniqueSmiles(tmpMolecule, true, true));
        Assertions.assertNull(ChemUtil.createUniqueSmiles(tmpMolecule, true, false));
        //fullerene C60 (Kekulé SMILES code of PubChem CID 123591) and other Kekulé structures of it in other atom orders
        IAtomContainer tmpFullerene = ChemUtil.parseSmilesToAtomContainer("C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3"
                + "=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4"
                + "=C3C3=C5C(=C81)C%10=C23", false, false);
        String tmpFullereneSmiles = ChemUtil.createUniqueSmiles(tmpFullerene, false, true);
        Assertions.assertNotNull(tmpFullereneSmiles);
        Random tmpRandom = new Random(ChemUtilKekuleSmilesTest.SEED);
        for (int i = 0; i < 3; i++) {
            IAtomContainer tmpAlternative = ChemUtilKekuleSmilesTest.permuteAtoms(tmpFullerene, tmpRandom);
            for (IAtom tmpAtom : tmpAlternative.atoms()) {
                tmpAtom.setIsAromatic(true);
            }
            for (IBond tmpBond : tmpAlternative.bonds()) {
                tmpBond.setIsAromatic(true);
                tmpBond.setOrder(IBond.Order.UNSET);
            }
            Kekulization.kekulize(tmpAlternative);
            ChemUtilKekuleSmilesTest.clearAromaticityFlags(tmpAlternative);
            Assertions.assertEquals(tmpFullereneSmiles, ChemUtil.createUniqueSmiles(tmpAlternative, false, true));
        }
    }
    //
    /**
     * Data set check on BIOFACQUIM (Tutorial/biofacquim-06-2025.sdf), imported like in the GUI with the default
     * settings (implicit hydrogen atoms, no kekulization enforced): for every molecule and with and without
     * stereochemistry, the unique SMILES code of the molecule as imported equals the one of three random atom
     * permutations and of an alternative Kekulé structure (aromaticity perceived, re-kekulized in a permuted atom
     * order, aromaticity flags removed again). Then, the CDK circular fragmenter and the Sugar Removal Utility
     * fragmenter (default settings) must not list one compound under two SMILES codes. Compound identity: the
     * non-standard InChI with fixed hydrogen layer (and stereo layers if stereochemistry is regarded), created by the
     * CDK InChI module from the parsed fragment SMILES code.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testBiofacquimDataSet() throws Exception {
        File tmpFile = Paths.get("Tutorial", "biofacquim-06-2025.sdf").toFile();
        for (boolean tmpIsStereo : new boolean[] {true, false}) {
            SettingsContainer tmpSettings = new SettingsContainer();
            tmpSettings.setKeepAtomContainerInDataModelSetting(true);
            List<MoleculeDataModel> tmpMolecules = new Importer(tmpSettings).importMoleculeFile(tmpFile, tmpIsStereo, true, false);
            Assertions.assertEquals(583, tmpMolecules.size());
            Random tmpRandom = new Random(ChemUtilKekuleSmilesTest.SEED);
            List<String> tmpFailures = new ArrayList<>();
            int tmpAlternativesCount = 0;
            for (MoleculeDataModel tmpMoleculeDataModel : tmpMolecules) {
                IAtomContainer tmpMolecule = tmpMoleculeDataModel.getAtomContainer();
                String tmpSmiles = ChemUtil.createUniqueSmiles(tmpMolecule, tmpIsStereo);
                Assertions.assertEquals(tmpMoleculeDataModel.getUniqueSmiles(), tmpSmiles);
                for (int i = 0; i < 3; i++) {
                    String tmpPermutedSmiles = ChemUtil.createUniqueSmiles(ChemUtilKekuleSmilesTest.permuteAtoms(tmpMolecule, tmpRandom), tmpIsStereo);
                    if (!tmpSmiles.equals(tmpPermutedSmiles)) {
                        tmpFailures.add(tmpMoleculeDataModel.getName() + " permutation: " + tmpSmiles + " vs " + tmpPermutedSmiles);
                    }
                }
                IAtomContainer tmpAlternative = ChemUtilKekuleSmilesTest.createAlternativeKekuleStructure(tmpMolecule, tmpRandom);
                if (tmpAlternative != null) {
                    tmpAlternativesCount++;
                    String tmpAlternativeSmiles = ChemUtil.createUniqueSmiles(tmpAlternative, tmpIsStereo);
                    if (!tmpSmiles.equals(tmpAlternativeSmiles)) {
                        tmpFailures.add(tmpMoleculeDataModel.getName() + " Kekulé structure: " + tmpSmiles + " vs " + tmpAlternativeSmiles);
                    }
                }
            }
            Assertions.assertTrue(tmpFailures.isEmpty(), tmpFailures.size() + " differences, e.g. " + tmpFailures.subList(0, Math.min(5, tmpFailures.size())));
            Assertions.assertTrue(tmpAlternativesCount > 500, "too few alternative Kekulé structures: " + tmpAlternativesCount);
            //fragments
            for (IMoleculeFragmenter tmpFragmenter : new IMoleculeFragmenter[] {new CDKCircularFragmenter(), new SugarRemovalUtilityFragmenter()}) {
                tmpFragmenter.restoreDefaultSettings();
                Map<String, FragmentDataModel> tmpFragmentTable = new ConcurrentHashMap<>();
                String tmpFragmentationName = tmpFragmenter.getFragmentationAlgorithmName() + tmpIsStereo;
                new FragmentationTask(tmpMolecules, tmpFragmenter, tmpFragmentTable, tmpFragmentationName, tmpIsStereo).call();
                Assertions.assertFalse(tmpFragmentTable.isEmpty());
                Map<String, List<String>> tmpSmilesByIdentity = new HashMap<>();
                int tmpWithoutInChICount = 0;
                for (String tmpFragmentSmiles : tmpFragmentTable.keySet()) {
                    String tmpIdentity = ChemUtilKekuleSmilesTest.createInChI(tmpFragmentSmiles, tmpIsStereo);
                    if (tmpIdentity == null) {
                        tmpWithoutInChICount++;
                        tmpIdentity = "no InChI: " + tmpFragmentSmiles;
                    }
                    tmpSmilesByIdentity.computeIfAbsent(tmpIdentity, aKey -> new ArrayList<>(2)).add(tmpFragmentSmiles);
                }
                Assertions.assertTrue(tmpWithoutInChICount * 100 < tmpFragmentTable.size(), tmpWithoutInChICount + " fragments without InChI");
                List<List<String>> tmpSplitCompounds = new ArrayList<>();
                for (List<String> tmpSpellings : tmpSmilesByIdentity.values()) {
                    if (tmpSpellings.size() > 1) {
                        Collections.sort(tmpSpellings);
                        tmpSplitCompounds.add(tmpSpellings);
                    }
                }
                Assertions.assertTrue(tmpSplitCompounds.isEmpty(), tmpFragmenter.getFragmentationAlgorithmName() + ", stereo "
                        + tmpIsStereo + ": " + tmpSplitCompounds.size() + " compounds listed under more than one SMILES code, e.g. "
                        + tmpSplitCompounds.subList(0, Math.min(5, tmpSplitCompounds.size())));
            }
        }
    }
    //
    /**
     * Parses the given SMILES code without kekulization and atom typing and creates the unique SMILES code.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether stereochemistry should be encoded
     * @param isAromatic whether aromaticity should be encoded
     * @return unique SMILES code
     * @throws Exception if the SMILES code cannot be parsed
     */
    private String createSmiles(String aSmiles, boolean isStereo, boolean isAromatic) throws Exception {
        return ChemUtil.createUniqueSmiles(ChemUtil.parseSmilesToAtomContainer(aSmiles, false, false), isStereo, isAromatic);
    }
    //
    /**
     * Returns a copy of the given molecule with randomly permuted atom order.
     *
     * @param aMolecule molecule
     * @param aRandom random number generator
     * @return permuted copy
     */
    private static IAtomContainer permuteAtoms(IAtomContainer aMolecule, Random aRandom) {
        List<Integer> tmpIndices = new ArrayList<>(aMolecule.getAtomCount());
        for (int i = 0; i < aMolecule.getAtomCount(); i++) {
            tmpIndices.add(i);
        }
        Collections.shuffle(tmpIndices, aRandom);
        int[] tmpPermutation = new int[tmpIndices.size()];
        for (int i = 0; i < tmpPermutation.length; i++) {
            tmpPermutation[i] = tmpIndices.get(i);
        }
        return new AtomContainerAtomPermutor(aMolecule).containerFromPermutation(tmpPermutation);
    }
    //
    /**
     * Returns a copy of the given molecule without aromaticity flags whose aromatic rings (Daylight model) are
     * kekulized anew in a random atom order, i.e. usually with another double bond placement, or null if the copy
     * cannot be kekulized.
     *
     * @param aMolecule molecule without aromaticity flags
     * @param aRandom random number generator
     * @return copy with another Kekulé structure or null
     */
    private static IAtomContainer createAlternativeKekuleStructure(IAtomContainer aMolecule, Random aRandom) {
        try {
            IAtomContainer tmpCopy = ChemUtilKekuleSmilesTest.permuteAtoms(aMolecule, aRandom);
            new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.relevant())).apply(tmpCopy);
            for (IBond tmpBond : tmpCopy.bonds()) {
                if (tmpBond.isAromatic()) {
                    tmpBond.setOrder(IBond.Order.UNSET);
                }
            }
            Kekulization.kekulize(tmpCopy);
            ChemUtilKekuleSmilesTest.clearAromaticityFlags(tmpCopy);
            return tmpCopy;
        } catch (Exception anException) {
            return null;
        }
    }
    //
    /**
     * Removes all aromaticity flags of the given molecule.
     *
     * @param aMolecule molecule
     */
    private static void clearAromaticityFlags(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpAtom.setIsAromatic(false);
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            tmpBond.setIsAromatic(false);
        }
    }
    //
    /**
     * Sets the aromaticity flags of the atoms and bonds of the six-membered ring of the given molecule that contains
     * an atom bonded to an oxygen atom (one ring of the smallest set of smallest rings).
     *
     * @param aMolecule molecule
     */
    private static void flagHydroxylatedSixMemberedRing(IAtomContainer aMolecule) {
        for (IAtomContainer tmpRing : Cycles.sssr(aMolecule).toRingSet().atomContainers()) {
            boolean tmpIsHydroxylated = false;
            for (IAtom tmpAtom : tmpRing.atoms()) {
                for (IAtom tmpNeighbour : aMolecule.getConnectedAtomsList(tmpAtom)) {
                    if (tmpNeighbour.getAtomicNumber() == 8) {
                        tmpIsHydroxylated = true;
                    }
                }
            }
            if (tmpRing.getAtomCount() == 6 && tmpIsHydroxylated) {
                for (IAtom tmpAtom : tmpRing.atoms()) {
                    tmpAtom.setIsAromatic(true);
                }
                for (IBond tmpBond : tmpRing.bonds()) {
                    tmpBond.setIsAromatic(true);
                }
                return;
            }
        }
    }
    //
    /**
     * Creates the non-standard InChI with fixed hydrogen layer of the molecule parsed from the given SMILES code, with
     * or without stereo layers, or returns null if no InChI can be created.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether the stereo layers are created
     * @return InChI or null
     */
    private static String createInChI(String aSmiles, boolean isStereo) {
        try {
            IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(aSmiles, true, false);
            InChIGenerator tmpGenerator = isStereo
                    ? InChIGeneratorFactory.getInstance().getInChIGenerator(tmpMolecule, InchiFlag.FixedH)
                    : InChIGeneratorFactory.getInstance().getInChIGenerator(tmpMolecule, InchiFlag.FixedH, InchiFlag.SNon);
            if (tmpGenerator.getStatus() == InchiStatus.ERROR) {
                return null;
            }
            return tmpGenerator.getInchi();
        } catch (Exception anException) {
            return null;
        }
    }
}
