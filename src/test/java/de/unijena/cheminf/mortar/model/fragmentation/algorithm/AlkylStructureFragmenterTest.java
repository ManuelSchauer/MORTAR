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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.io.FileInputStream;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Class to test the correct working of
 * {@link AlkylStructureFragmenter}.
 *
 * @author Maximilian Rottmann
 * @version 1.0.0.0
 */
public class AlkylStructureFragmenterTest {
    /**
     * Static setter for Locale.
     * If Locale is set in Test constructor, a MissingResourceException is thrown when user locale diverts from set Locale.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    /**
     * Constructor of AlkylStructureFragmenter test class.
     */
    private AlkylStructureFragmenterTest() {}
    //
    //<editor-fold desc="Example Usage">
    /**
     * Method to showcase the usage of the AlkylStructureFragmenter class.
     * The molecules used in this 'test' are all molecules used in the tests below, complemented by various substructures
     * the class will detect and return intact.
     */
    @Test
    public void exampleUsageTest() throws Exception {
        //get file reader with test structures file
        URL tmpURL = this.getClass().getResource("ASF_Test_Substructure_Library.sdf");
        Assertions.assertNotNull(tmpURL);
        File tmpResourceFile = Paths.get(tmpURL.toURI()).toFile();
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(new FileInputStream(tmpResourceFile),
                SilentChemObjectBuilder.getInstance());
        //read structures into AtomContainerSet
        IAtomContainerSet tmpStructuresSet = new AtomContainerSet();
        while (tmpSDFReader.hasNext()) {
            tmpStructuresSet.addAtomContainer(tmpSDFReader.next());
        }
        //get file reader with test structures expected fragments file
        URL tmpExpectedFragmentsURL = this.getClass().getResource("ASF_Test_Substructure_Library_Fragments_Strings");
        Assertions.assertNotNull(tmpExpectedFragmentsURL);
        List<String> tmpExpectedSMILESStringList = Files.readAllLines(
                Path.of(tmpExpectedFragmentsURL.toURI()),
                StandardCharsets.UTF_8);
        //constructs an AlkylStructureFragmenter with default settings
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        //generate fragments and add SMILES to fragment list
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        for (IAtomContainer tmpAC : tmpStructuresSet.atomContainers()) {
            this.preprocessTestMolecule(tmpASF, tmpAC,
                    false, false, true);
            int tmpPreFragmentationCount = this.countAtoms(tmpAC);
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpAC);
            if (this.checkChemicalFormula(tmpPreFragmentationCount, tmpACList)) {
                tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
            }
        }
        //check if match
        boolean tmpTestSuccessful = this.compareListsIgnoringOrder(tmpExpectedSMILESStringList, tmpFragmentsSMILESList);
        if (!tmpTestSuccessful) {
            Assertions.fail();
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Unit Tests">
    /**
     * Method testing the fragmentation of a conceptual molecule containing several key structural features detected and
     * extracted by the AlkylStructureFragmenter - an isolated double bond, a quaternary C, an aromatic ring system with
     * a connected conjugated system, an isolated cyclohexane ring and a tertiary C - in default settings.
     *
     * @throws Exception if SMILES cannot be parsed correctly or fragmentation fails
     */
    @Test
    public void extractFragmentsTest() throws Exception {
        ArrayList<String> tmpExpectedFragmentList = new ArrayList<>(12);
        tmpExpectedFragmentList.add("C=CC1=CC=C2C=CC=CC2=C1"); //naphthalene with conjugated double bond
        tmpExpectedFragmentList.add("*C(*)*"); //isolated tertiary carbon
        tmpExpectedFragmentList.add("C"); //methyl groups of the tertiary carbon
        tmpExpectedFragmentList.add("C");
        tmpExpectedFragmentList.add("*C(*)(*)*"); //isolated quaternary carbon
        tmpExpectedFragmentList.add("C"); //methyl groups of the quaternary carbon
        tmpExpectedFragmentList.add("C");
        tmpExpectedFragmentList.add("C"); //methylene groups next to the quaternary carbon
        tmpExpectedFragmentList.add("C");
        tmpExpectedFragmentList.add("C=C"); //isolated double bond
        tmpExpectedFragmentList.add("CCC"); //propyl chain
        tmpExpectedFragmentList.add("C1CCCCC1"); //cyclohexane
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestAC = tmpParser.parseSmiles("C12=CC(=CC=C1C(=CC(=C2)C3CCCCC3)CC(C)(C)C/C=C/CCC)/C=C\\C(C)C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        ArrayList<String> tmpActualFragmentList = new ArrayList<>(12);
        for (IAtomContainer tmpAC : tmpASF.fragmentMolecule(tmpTestAC)) {
            tmpActualFragmentList.add(ChemUtil.createUniqueSmiles(tmpAC, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(tmpActualFragmentList, tmpExpectedFragmentList));
    }
    //<editor-fold desc="Classification Tests">
    /**
     * Test for correct classification of tertiary and quaternary carbon atoms and of the neighbors forming a unit with
     * them if they are not isolated. Atom types are deliberately not perceived, the classification has to work on the
     * bonds of the molecule alone.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkTertQuatAndNeighbors() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        //isobutane: atom 1 is tertiary, atoms 0, 2, 3 form a unit with it
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon, 1);
        this.assertAtoms(tmpClassification.isQuaternaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2, 3);
        //neopentyl chain: atom 1 is quaternary, atom 5 is not part of its unit
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)(C)CC"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isQuaternaryCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2, 3, 4);
        //branching atom with a double bond (isobutene) is no tertiary carbon
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)=C"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember);
        //isolated tertiary carbons (default) do not form a unit with their neighbors
        tmpClassification = new AlkylStructureFragmenter().classify(tmpParser.parseSmiles("CC(C)C"));
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember);
    }

    /**
     * Test for correct classification of ring atoms and bonds.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkRings() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //methylcyclopropane: methyl atom 0 and its bond are not part of the ring
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1CC1");
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = new AlkylStructureFragmenter().classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isRingAtom, 1, 2, 3);
        this.assertBonds(tmpMolecule, tmpClassification.isRingBond, new int[][]{{1, 2}, {2, 3}, {3, 1}});
    }

    /**
     * Test for correct classification of double and triple bond atoms.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkMultiBonds() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = new AlkylStructureFragmenter().classify(tmpParser.parseSmiles("CC=CC#C"));
        this.assertAtoms(tmpClassification.hasDoubleBond, 1, 2);
        this.assertAtoms(tmpClassification.hasTripleBond, 3, 4);
        this.assertAtoms(tmpClassification.isChainAtom, 0);
    }

    /**
     * Test for correct classification of tertiary/quaternary carbon atoms directly bonded to a ring, which are not
     * isolated unless they are to be separated from rings.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkConnectedTertQuatRing() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        //isopropylcyclopropane: tertiary atom 1 is bonded to ring atom 3 and kept at the ring by default
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C1CC1"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatCarbonBondedToRing, 1);
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C1CC1"));
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon, 1);
        //ring atoms themselves are never tertiary carbons
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC1CC1"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatCarbonBondedToRing);
    }

    /**
     * Test for the detection of conjugated pi bond systems.
     * Hexatriene has to be classified completely. Double and triple bonds separated by a single bond (D-S-T) are not
     * regarded as conjugated, also not if the double bond is part of an aromatic ring (phenylacetylene).
     *
     * @throws InvalidSmilesException if SMILES-parser is not able to parse SMILES due to incorrect syntax
     */
    @Test
    public void testMarkConjugatedPiSystems() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C=CC=CC=C");
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isConjugatedAtom, 0, 1, 2, 3, 4, 5);
        this.assertBonds(tmpMolecule, tmpClassification.isConjugatedBond, new int[][]{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
        //D-S-T is not conjugated
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("C=CC#C"));
        this.assertAtoms(tmpClassification.isConjugatedAtom);
        //phenylacetylene: only the ring is conjugated
        tmpMolecule = tmpParser.parseSmiles("C#CC1=CC=CC=C1");
        tmpClassification = tmpASF.classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isConjugatedAtom, 2, 3, 4, 5, 6, 7);
        this.assertBonds(tmpMolecule, tmpClassification.isConjugatedBond, new int[][]{{2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 2}});
    }
    //</editor-fold>
    /**
     * Test that the returned fragments consist of new atom and bond instances: changing a fragment must not change the
     * input molecule (formerly a test of the deep copy methods with butene, which used to make problems).
     *
     * @throws Exception if SMILES cannot be parsed or fragmentation fails
     */
    @Test
    public void deepCopyButeneTest() throws Exception {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpButeneContainer = tmpParser.parseSmiles("C=CCC");
        String tmpButeneSMILES = ChemUtil.createUniqueSmiles(tmpButeneContainer, false);
        List<IAtomContainer> tmpFragments = tmpASF.fragmentMolecule(tmpButeneContainer);
        Assertions.assertTrue(this.compareListsIgnoringOrder(this.generateSMILESFromACList(tmpFragments, false),
                new ArrayList<>(List.of("C=C", "CC"))));
        for (IAtomContainer tmpFragment : tmpFragments) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                for (IAtom tmpInputAtom : tmpButeneContainer.atoms()) {
                    Assertions.assertNotSame(tmpInputAtom, tmpAtom);
                }
                tmpAtom.setImplicitHydrogenCount(0);
                tmpAtom.setAtomicNumber(14);
            }
            for (IBond tmpBond : tmpFragment.bonds()) {
                tmpBond.setOrder(IBond.Order.TRIPLE);
            }
        }
        Assertions.assertEquals(tmpButeneSMILES, ChemUtil.createUniqueSmiles(tmpButeneContainer, false));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Basic Fragmentation Tests">
    /**
     * Method testing correct fragmentation with a basic example molecule.
     * This test focuses on tertiary and quaternary carbon systems, tested with a conceptual molecule.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void basicTest01() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CC(C)(C)CCCC(C)C
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC(C)(C)CCCC(C)C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
             tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("*C(*)(*)*");
        tmpExpectedSMILESList.add("*C(*)*");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate setting value
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)(C)C");
        tmpExpectedSMILESList.add("CC(C)C");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a basic example molecule.
     * This test focuses on extraction of isolated double and triple bonds,
     * tested with the corresponding butenes and butynes.
     * Additionally, the non-fragmentation of a cyclo-hexane with an internal "isolated" double bond is also tested.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void basicTest02() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CC=CC
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC=CC");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C=C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: C=CCC
        tmpTestStructureAC = tmpParser.parseSmiles("C=CCC");
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpTestStructureAC));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpTestStructureAC));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpTestStructureAC));
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("C=C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: CC#CC
        tmpTestStructureAC = tmpParser.parseSmiles("CC#CC");
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpTestStructureAC));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpTestStructureAC));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpTestStructureAC));
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C#C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: CC#CC
        tmpTestStructureAC = tmpParser.parseSmiles("C1=CCCCC1");
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpTestStructureAC));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpTestStructureAC));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpTestStructureAC));
        tmpACList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList.clear();
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("C1=CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a basic example molecule.
     * This test focuses on safe retention of rings, tested with Cyclohexane and Benzene.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void basicTest03() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CCCCC1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CCCCC1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: c1ccccc1
        tmpTestStructureAC = tmpParser.parseSmiles("c1ccccc1");
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpTestStructureAC));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpTestStructureAC));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpTestStructureAC));
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccccc1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a basic example molecule.
     * This test focuses on safe retention of more complex ring systems, tested with Pentalene and Biphenyl.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void basicTest04() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CC2CCCC2C1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CC2CCCC2C1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C1CC2CCCC2C1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: C=1C=CC(=CC1)C2=CC=CC=C2
        tmpTestStructureAC = tmpParser.parseSmiles("C=1C=CC(=CC1)C2=CC=CC=C2");
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpTestStructureAC));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpTestStructureAC));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpTestStructureAC));
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("C=1C=CC(=CC1)C2=CC=CC=C2");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a basic example molecule.
     * This test focuses on separation of ring systems, tested with Bicyclohexyl.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void basicTest05() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CC(CCC1)C1CCCCC1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CC(CCC1)C1CCCCC1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C1CCCCC1");
        tmpExpectedSMILESList.add("C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate setting value since two connected rings also have at least tertiary configuration at their joint atoms
        //though there should be no difference in resulting fragments
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing the correct handling of spiro configurated rings with an example molecule.
     * The molecule used in this test is Spiro[5.5]undecane.
     *
     * @throws InvalidSmilesException if SMILES can not be parsed
     * @throws CloneNotSupportedException if cloning fails
     */
    @Test
    public void basicTest06() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CCC2(CCCCC2)CC1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CCC2(CCCCC2)CC1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C1CCC2(CC1)CCCCC2");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing the correct handling of tertiary and quaternary atoms attached to rings with an example molecule.
     * The molecule used in this test is a concept molecule comprised of a cyclohexane bonded to a quaternary carbon system.
     *
     * @throws InvalidSmilesException if SMILES can not be parsed
     * @throws CloneNotSupportedException if cloning fails
     */
    @Test
    public void basicTest07() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CCCCC1C(C)(C)C
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CCCCC1C(C)(C)C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        //for default settings:
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("CC(C)(C)C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //combinations of relevant settings
        //tmpASF.setIsolateTertQuatCarbonsSetting(true); -> default value
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("*C(*)(*)*");
        tmpExpectedSMILESList.add("C1CCCCC1");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                    new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        //tmpASF.setSeparateTertQuatCarbonFromRingSetting(true); -> still correct value from test above
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("*C(C)(C)C");
        tmpExpectedSMILESList.add("C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        //tmpASF.setIsolateTertQuatCarbonsSetting(false); -> still correct value from test above
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)(C)C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method to test correct extraction of allene structures.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the atomcontainer in fragmentation is not supported
     */
    @Test
    public void basicTest08() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC=C=CC");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C(=C)=C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct handling of conjugated aromatic rings, conjoined with a single bond,
     * making the whole molecule conjugated.
     * Tested with molecule Biphenyl.
     *
     * @throws InvalidSmilesException if SMILES can not be parsed
     * @throws CloneNotSupportedException if cloning of ASF is not supported
     */
    @Test
    public void basicTest09() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("c1ccccc1c2ccccc2");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("c1ccc(cc1)c2ccccc2");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing for correct behavior in dissection and separation of linear carbon chains of varying sizes (including 'no restrictions applied').
     */
    @Test
    public void defaultLinearChainDissectionTest() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpCarbonChainAC = tmpParser.parseSmiles("CCCCCCCCCCCCCC");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpCarbonChainAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpCarbonChainAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpCarbonChainAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate setting value -> maximum length of linear chains = 7
        tmpASF.setMaxChainLengthSetting(7);
        tmpFragmentsSMILESList.clear();
        tmpExpectedSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpCarbonChainAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpCarbonChainAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("CCCCCCC");
        tmpExpectedSMILESList.add("CCCCCCC");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setFragmentSideChainsSetting(false);
        tmpExpectedSMILESList.clear();
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpCarbonChainAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpCarbonChainAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("CCCCCCCCCCCCCC");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Specific Fragmentation Tests">
    //all molecules are tested with all relevant settings
    //example: only changing maximum size of linear chains in a linear C6 molecule
    /**
     * Method testing correct fragmentation with a specific example molecule.
     * This test focuses on extraction of ring systems, tested with a derivative of the molecule Hapalindole B after
     * partial fragmentation with the (in MORTAR included) ErtlFunctionalGroupsFinder.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void specificTest01() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CC1CCC2C(C1)C(=C)C3=C(C=CC=C3)C2(C)C (Hapalindole B Derivative)
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC1CCC2C(C1)C(=C)C3=C(C=CC=C3)C2(C)C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C=C1C=2C=CC=CC2CC3CCCCC13");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a specific example molecule.
     * This test focuses on extraction of quaternary carbon systems, tested with the molecule 4,4,6,6-Tetramethylnonane.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void specificTest02() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CCCC(C)(C)CC(C)(C)CCC
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CCCC(C)(C)CC(C)(C)CCC");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("*C(*)(*)*");
        tmpExpectedSMILESList.add("*C(*)(*)*");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate setting values
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("CC(C)(C)CC(C)(C)C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Method testing correct fragmentation with a specific example molecule.
     * This test focuses on extraction of a singular ring and separation of linear, conjugated Pi systems from side chains,
     * tested with a derivative of the molecule Dehydropinguisanin after fragmentation with the (in MORTAR included) ErtlFunctionalGroupsFinder.
     *
     * @throws InvalidSmilesException if SMILES is not correctly parsed of otherwise faulty
     * @throws CloneNotSupportedException if something goes wrong during cloning step in fragmentation
     */
    @Test
    public void specificTest03() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C=CC=CCC1(C)C(C)CCC1C (Dehydropinguisanin)
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C=CC=CCC1(C)C(C)CCC1C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C=CC=C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C1CCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Tests molecule "CC(C)(C)CC1CCC(=C)C2CC21" which showed difficulties in fragmentation in past versions.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    public void specificTest04() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CC(C)(C)CC1CCC(=C)C2CC21
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC(C)(C)CC1CCC(=C)C2CC21");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("*C(*)(*)*");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C=C1CCCC2CC12");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)(C)C");
        tmpExpectedSMILESList.add("C=C1CCCC2CC12");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Tests correct functionality of tetrahedral stereo chemistry preservation with concept molecule containing a central
     * chiral atom, surrounded by four ring structures. Additionally, a second concept molecule is tested, with three
     * ligand groups and an implicit hydrogen connected to the chiral atom.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    public void specificTest05() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C1CCC(C1)[C@](c1cccc2ccccc12)(C1CCCCC1)c1ccccc1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C1CCC(C1)[C@](c1cccc2ccccc12)(C1CCCCC1)c1ccccc1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(false);
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpASF.setPreserveStereoChemistrySetting(true);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, true));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("c1ccc(cc1)[C@@](C2CCCCC2)(C3CCCC3)c4cccc5ccccc54");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setPreserveStereoChemistrySetting(false);
        tmpFragmentsSMILESList.clear();
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccc(cc1)C(c2cccc3ccccc32)(C4CCCCC4)C5CCCC5");
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //test structure: c1ccccc1[C@@H](C1CCC2CCCC12)C1CCCCC1
        tmpTestStructureAC = tmpParser.parseSmiles("c1ccccc1[C@@H](C1CCC2CCCC12)C1CCCCC1");
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(false);
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpASF.setPreserveStereoChemistrySetting(true);
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList.clear();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, true));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccc(cc1)[C@H](C2CCCCC2)C3CCC4CCCC43");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setPreserveStereoChemistrySetting(false);
        tmpFragmentsSMILESList.clear();
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccc(cc1)C(C2CCCCC2)C3CCC4CCCC43");
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Tests correct functionality of double bond stereo chemistry preservation with concept molecule with consecutive double bonds.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest06() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: C=C\C=C/C=C
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C=C\\C=C/C=C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(false);
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpASF.setPreserveStereoChemistrySetting(true);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, true));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("C=C/C=C\\C=C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setPreserveStereoChemistrySetting(false);
        tmpFragmentsSMILESList.clear();
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("C=CC=CC=C");
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }/**
     * Tests correct functionality of specific ring, tertiary and linear chain configuration.
     * It is tested with a derivative of Cholesterin (CNP0174859.1), moreover, the alkyl base structure
     * (SMILES: CCCC[C@H](C)[C@H]1CC[C@]2(C)[C@H](CC[C@@H]2[C@H]1C)[C@H](C)CCCC(C)C ).
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest07() throws InvalidSmilesException, CloneNotSupportedException {
        //test structure: CCCC[C@H](C)[C@H]1CC[C@]2(C)[C@H](CC[C@@H]2[C@H]1C)[C@H](C)CCCC(C)C
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CCCC[C@H](C)[C@H]1CC[C@]2(C)[C@H](CC[C@@H]2[C@H]1C)[C@H](C)CCCC(C)C");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, true));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("*C(*)*");
        //the two former stereocenters next to the ring bear two methyl groups after fragmentation and are therefore no
        // stereocenters anymore, only the three ring stereocenters are kept
        tmpExpectedSMILESList.add("CC(C)[C@H]1CCC2[C@H](CC[C@@H]2C(C)C)C1");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }

    /**
     * Tests a tertiary/quaternary carbon kept at a ring (default settings) that is bonded to a second, isolated
     * tertiary/quaternary carbon. The isolated carbon is returned saturated with pseudo atoms, the cut bond of the
     * carbon kept at the ring is marked with a pseudo atom as well, so that its branching is preserved in the ring
     * fragment instead of being saturated with hydrogen.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest08() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        //test structure: tertiary carbon at a ring bonded to an isopropyl group
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC(C)C(C)C1CCCCC1");
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
        List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
        List<String> tmpExpectedSMILESList = new ArrayList<>(List.of("*C(*)*", "*C(C)C1CCCCC1", "C", "C"));
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)), tmpFragmentsSMILESList::toString);
        //test structure: quaternary carbon at a ring bonded to a tert-butyl group
        tmpTestStructureAC = tmpParser.parseSmiles("CC(C)(C)C(C)(C)C1CCCCC1");
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
        tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
        tmpExpectedSMILESList = new ArrayList<>(List.of("*C(*)(*)*", "*C(C)(C)C1CCCCC1", "C", "C", "C"));
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)), tmpFragmentsSMILESList::toString);
    }

    /**
     * Tests non-isolated tertiary/quaternary carbon systems (isolation setting inactive) that are directly bonded to an
     * isolated double or triple bond. The multi bond is extracted separately, the cut bond of the tertiary/quaternary
     * carbon is marked with a pseudo atom so that the carbon keeps its original number of bonds.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest09() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        Map<String, List<String>> tmpExpectedFragmentsMap = new HashMap<>();
        tmpExpectedFragmentsMap.put("C=CC(C)(C)C", List.of("C=C", "*C(C)(C)C"));
        tmpExpectedFragmentsMap.put("C#CC(C)(C)C", List.of("C#C", "*C(C)(C)C"));
        tmpExpectedFragmentsMap.put("C=CC(C)C", List.of("C=C", "*C(C)C"));
        for (Map.Entry<String, List<String>> tmpEntry : tmpExpectedFragmentsMap.entrySet()) {
            IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles(tmpEntry.getKey());
            this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
            int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
            Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
            List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
            Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                    new ArrayList<>(tmpEntry.getValue())), () -> tmpEntry.getKey() + " -> " + tmpFragmentsSMILESList);
        }
    }

    /**
     * Checks all fragmentation examples given in the class documentation of {@link AlkylStructureFragmenter}, so that
     * the documented rules and the implementation cannot diverge. Expected and actual fragments are compared as unique
     * SMILES with aromatic symbols.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    void documentationExamplesTest() throws Exception {
        //columns: input SMILES, isolate tert/quat, separate from ring, max chain length, expected fragments
        Object[][] tmpExamples = {
                {"C1CCC(CC1)C1CCCCC1", true, false, 6, "C1CCCCC1 C1CCCCC1"},
                {"C=Cc1ccccc1", true, false, 6, "C=Cc1ccccc1"},
                {"C=CC#C", true, false, 6, "C=C C#C"},
                {"cCCC", true, false, 6, "c CCC"},
                {"C=C1CCCCC1", true, false, 6, "C=C1CCCCC1"},
                {"CC=CC", true, false, 6, "C=C C C"},
                {"CC(C)(C)CCCC(C)C", true, false, 6, "*C(*)(*)* *C(*)* C C C C C CCC"},
                {"c1ccc(cc1)C(C)CC", true, true, 6, "c1ccccc1 *C(*)* C CC"},
                {"CC(C)(C)CCCC(C)C", false, false, 6, "CC(C)(C)C CC(C)C C"},
                {"CC(C)c1ccccc1", true, false, 6, "CC(C)c1ccccc1"},
                {"CC(C1CCCCC1)C=CC=C", true, false, 6, "C=CC=CC(C)C1CCCCC1"},
                {"C=CC(C)(C)C", false, false, 6, "C=C *C(C)(C)C"},
                {"CC(C)c1ccccc1", false, true, 6, "*C(C)C c1ccccc1"},
                {"CC(C)C(C)C1CCCCC1", true, false, 6, "*C(*)* *C(C)C1CCCCC1 C C"},
                {"CC1CCCCC1", true, false, 6, "C1CCCCC1 C"},
                {"CCCCCCCC", true, false, 3, "CCC CCC CC"}};
        for (Object[] tmpExample : tmpExamples) {
            AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
            tmpASF.setIsolateTertQuatCarbonsSetting((boolean) tmpExample[1]);
            tmpASF.setSeparateTertQuatCarbonFromRingSetting((boolean) tmpExample[2]);
            tmpASF.setMaxChainLengthSetting((int) tmpExample[3]);
            IAtomContainer tmpMolecule = AlkylTestMoleculeGenerator.parseAsInPipeline((String) tmpExample[0]);
            Assertions.assertNotNull(tmpMolecule, () -> "Invalid example SMILES " + tmpExample[0]);
            this.preprocessTestMolecule(tmpASF, tmpMolecule, false, false, true);
            List<String> tmpActual = this.generateSMILESFromACList(tmpASF.fragmentMolecule(tmpMolecule), false);
            List<String> tmpExpected = new ArrayList<>();
            for (String tmpExpectedSmiles : ((String) tmpExample[4]).split(" ")) {
                tmpExpected.add(ChemUtil.createUniqueSmiles(AlkylTestMoleculeGenerator.parseAsInPipeline(tmpExpectedSmiles), false, true));
            }
            Collections.sort(tmpActual);
            Collections.sort(tmpExpected);
            Assertions.assertEquals(tmpExpected, tmpActual, tmpExample[0] + " " + Arrays.toString(tmpExample));
        }
    }

    /**
     * Tests tertiary/quaternary carbons that are bonded to a ring and separated from it (isolation and separation
     * settings active). The carbon is isolated like any other tertiary/quaternary carbon and its neighbors are treated as
     * ordinary chain atoms, so an ethyl group is returned as ethane. (The previous implementation additionally cut the
     * first chain atom of every neighbor and returned the ethyl group as two methanes.)
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest10() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        Map<String, List<String>> tmpExpectedFragmentsMap = new HashMap<>();
        tmpExpectedFragmentsMap.put("c1ccc(cc1)C(C)CC", List.of("c1ccccc1", "*C(*)*", "C", "CC"));
        tmpExpectedFragmentsMap.put("c1ccc(cc1)C(CC)CC", List.of("c1ccccc1", "*C(*)*", "CC", "CC"));
        tmpExpectedFragmentsMap.put("CC(C)(CCC)C1CCCCC1", List.of("C1CCCCC1", "*C(*)(*)*", "C", "C", "CCC"));
        for (Map.Entry<String, List<String>> tmpEntry : tmpExpectedFragmentsMap.entrySet()) {
            IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles(tmpEntry.getKey());
            this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
            int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
            Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
            List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
            Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                    new ArrayList<>(tmpEntry.getValue())), () -> tmpEntry.getKey() + " -> " + tmpFragmentsSMILESList);
        }
    }

    /**
     * Tests tertiary/quaternary carbons kept at a ring (default settings) that are bonded to a conjugated system. Like all
     * non-isolated tertiary/quaternary carbons (e.g. with the isolation setting deactivated), they keep their bond to
     * the conjugated system, so ring, tertiary/quaternary carbon and conjugated system form one fragment. (The previous
     * implementation cut this bond only in this setting combination, but kept it with isolation deactivated.)
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void specificTest11() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        Map<String, List<String>> tmpExpectedFragmentsMap = new HashMap<>();
        tmpExpectedFragmentsMap.put("CC(C1CCCCC1)C(=C)C1=CC=CC1", List.of("C=C(C1=CC=CC1)C(C)C2CCCCC2"));
        tmpExpectedFragmentsMap.put("CC(C1CCCCC1)C=CC=C", List.of("C=CC=CC(C)C1CCCCC1"));
        for (boolean tmpIsolate : new boolean[]{true, false}) {
            AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
            tmpASF.setIsolateTertQuatCarbonsSetting(tmpIsolate);
            for (Map.Entry<String, List<String>> tmpEntry : tmpExpectedFragmentsMap.entrySet()) {
                IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles(tmpEntry.getKey());
                this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
                int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
                List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
                Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
                List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
                Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                        new ArrayList<>(tmpEntry.getValue())), () -> tmpEntry.getKey() + " isolate=" + tmpIsolate + " -> " + tmpFragmentsSMILESList);
            }
        }
    }

    /**
     * Tests hydrogen atoms that remain explicit after hydrogen suppression (isotopes, molecular hydrogen). They are no
     * neighbors in the sense of the tertiary/quaternary carbon detection, stay bonded to their atom and do not form
     * fragments of their own; saturation does not create molecular hydrogen.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void remainingExplicitHydrogenTest() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        Map<String, List<String>> tmpExpectedFragmentsMap = new HashMap<>();
        //secondary carbon with deuterium, no tertiary carbon
        tmpExpectedFragmentsMap.put("[2H]C(C)C", List.of("[2H]C(C)C"));
        //tertiary carbon with deuterium, the deuterium stays on the carbon
        tmpExpectedFragmentsMap.put("[2H]C(C)(C)C", List.of("[2H]C(*)(*)*", "C", "C", "C"));
        //deuterium in a chain, the carbon is secondary and therefore cut from the ring
        tmpExpectedFragmentsMap.put("[2H]C(CC)C1CCCCC1", List.of("[2H]CCC", "C1CCCCC1"));
        //molecular hydrogen is no fragment of its own after saturation
        tmpExpectedFragmentsMap.put("CC.[H][H]", List.of("CC", "[H][H]"));
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (Map.Entry<String, List<String>> tmpEntry : tmpExpectedFragmentsMap.entrySet()) {
            IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles(tmpEntry.getKey());
            int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
            Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
            List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
            List<String> tmpExpectedSMILESList = new ArrayList<>();
            for (String tmpExpectedSmiles : tmpEntry.getValue()) {
                tmpExpectedSMILESList.add(ChemUtil.createUniqueSmiles(tmpParser.parseSmiles(tmpExpectedSmiles), false, true));
            }
            Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList), tmpExpectedSMILESList),
                    () -> tmpEntry.getKey() + " -> " + tmpFragmentsSMILESList);
        }
    }

    /**
     * Tests that the dissection of linear chains does not depend on the atom order without hydrogen saturation. The
     * chain is walked from its end that is bonded to the rest of the molecule, so the unsaturated (radical) ends of
     * the chain segments are always at the same positions.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void chainDissectionAtomOrderTest() throws InvalidSmilesException, CloneNotSupportedException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        tmpASF.setMaxChainLengthSetting(3);
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        for (String tmpExpectedSmiles : List.of("[CH2]C[CH2]", "[CH2]C", "[CH]1CCCCC1")) {
            tmpExpectedSMILESList.add(ChemUtil.createUniqueSmiles(tmpParser.parseSmiles(tmpExpectedSmiles), false, true));
        }
        //chain atoms before and after the ring atoms, attachment point at the start, in the middle and at the end
        for (String tmpSmiles : List.of("CCCCCC1CCCCC1", "C1CCCCC1CCCCC", "C(C1CCCCC1)CCCC", "C1CCCC(C1)CCCCC")) {
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpParser.parseSmiles(tmpSmiles));
            List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
            Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList), new ArrayList<>(tmpExpectedSMILESList)),
                    () -> tmpSmiles + " -> " + tmpFragmentsSMILESList);
        }
    }

    /**
     * Tests the fragmentation of aromatic remnants as produced by a previous fragmentation in a pipeline, e.g. alkane
     * fragments of the ErtlFunctionalGroupsFinderFragmenter whose aromatic ring was opened by the removal of a hetero
     * atom. These cannot be kekulized, so their aromatic bonds have no bond order. They are treated as conjugated pi
     * systems and kept intact, alkyl substituents are fragmented as usual. Cut bonds of the aromatic atoms are
     * saturated with hydrogen, so identical aromatic remnants result in identical fragments.
     *
     * @throws CDKException if a SMILES cannot be parsed
     * @throws CloneNotSupportedException if cloning of the original molecule is not supported
     */
    @Test
    void aromaticRemnantTest() throws CDKException, CloneNotSupportedException {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        Map<String, List<String>> tmpExpectedFragmentsMap = new HashMap<>();
        tmpExpectedFragmentsMap.put("ccc(c1ccc(cc1)C)CC", List.of("cccc1ccccc1", "C", "CC"));
        tmpExpectedFragmentsMap.put("CC(C)(C)ccc", List.of("*C(*)(*)*", "ccc", "C", "C", "C"));
        tmpExpectedFragmentsMap.put("cc(C)cc(CCCCCCCC)c", List.of("ccccc", "C", "CCCCCC", "CC"));
        tmpExpectedFragmentsMap.put("cc(CC)ccc", List.of("ccccc", "CC"));
        //single aromatic remnant atoms (no aromatic bond) stay conjugated with neighboring aromatic atoms and are
        // neither part of alkyl chains nor tertiary/quaternary carbons
        tmpExpectedFragmentsMap.put("c-c1ccccc1", List.of("c-c1ccccc1"));
        tmpExpectedFragmentsMap.put("cc-c", List.of("cc-c"));
        tmpExpectedFragmentsMap.put("cCCC", List.of("c", "CCC"));
        tmpExpectedFragmentsMap.put("c(C)(C)C", List.of("c", "C", "C", "C"));
        for (Map.Entry<String, List<String>> tmpEntry : tmpExpectedFragmentsMap.entrySet()) {
            //parsed as in the MORTAR pipeline, kekulization of these remnants fails
            Assertions.assertThrows(CDKException.class, () -> ChemUtil.parseSmilesToAtomContainer(tmpEntry.getKey(), true, true));
            IAtomContainer tmpTestStructureAC = ChemUtil.parseSmilesToAtomContainer(tmpEntry.getKey(), false, false);
            this.preprocessTestMolecule(tmpASF, tmpTestStructureAC, false, false, true);
            int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
            List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
            Assertions.assertTrue(this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList));
            List<String> tmpFragmentsSMILESList = this.generateSMILESFromACList(tmpACList, false);
            Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                    new ArrayList<>(tmpEntry.getValue())), () -> tmpEntry.getKey() + " -> " + tmpFragmentsSMILESList);
        }
        //aromatic remnants are only accepted with aromatic flags, bonds without order and flag are still filtered
        IAtomContainer tmpUnsetBondAC = ChemUtil.parseSmilesToAtomContainer("CCCC", false, false);
        tmpUnsetBondAC.getBond(1).setOrder(IBond.Order.UNSET);
        tmpASF.setKeepNonFragmentableMoleculesSetting(false);
        Assertions.assertTrue(tmpASF.shouldBeFiltered(tmpUnsetBondAC));
    }

    /**
     * Method to test a default alkyl structure fragmentation on a concept molecule covering a broad range of resulting fragments.
     *
     * @throws Exception if fragmentation does not result in expected fragments
     */
    @Test
    public void defaultFragmentationTest() throws Exception {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //test structure: CC(CC1C2CC2C(=C)CC1c1ccc(cc1C\C=C/c1ccccc1)C(C)C)(CCCCCCCCC)CC#CC
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC(CC1C2CC2C(=C)CC1c1ccc(cc1C\\C=C/c1ccccc1)C(C)C)(CCCCCCCCC)CC#CC");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>(12);
        tmpExpectedSMILESList.add("c1ccc(cc1)C(C)C"); //ring structure with connected tertiary carbon system
        tmpExpectedSMILESList.add("C=Cc1ccccc1"); //benzene with connected double bond forming conjugated pi bond system
        tmpExpectedSMILESList.add("C"); //methane between the substructures above
        tmpExpectedSMILESList.add("C=C1CCCC2CC12"); //ring system with connected non-cyclic double bond
        tmpExpectedSMILESList.add("C"); //methane connecting ring system above and quaternary carbon below
        tmpExpectedSMILESList.add("*C(*)(*)*"); //quaternary system in pseudoatom representation
        tmpExpectedSMILESList.add("C"); //methane residual of quaternary system above
        tmpExpectedSMILESList.add("CCCCCC"); //hexane residual of fragmented nonane chain
        tmpExpectedSMILESList.add("CCC"); //propane residual of fragmented nonane chain
        tmpExpectedSMILESList.add("C"); //methane connecting triple bond and quaternary system
        tmpExpectedSMILESList.add("C#C"); //triple bond
        tmpExpectedSMILESList.add("C"); //methane residual of triple bond chain
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate settings
        tmpASF.setFragmentSideChainsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccc(cc1)C(C)C"); //ring structure with connected tertiary carbon system
        tmpExpectedSMILESList.add("C=Cc1ccccc1"); //benzene with connected double bond forming conjugated pi bond system
        tmpExpectedSMILESList.add("C"); //methane between the substructures above
        tmpExpectedSMILESList.add("C=C1CCCC2CC12"); //ring system with connected non-cyclic double bond
        tmpExpectedSMILESList.add("C"); //methane connecting ring system above and quaternary carbon below
        tmpExpectedSMILESList.add("*C(*)(*)*"); //quaternary system in pseudoatom representation
        tmpExpectedSMILESList.add("C"); //methane residual of quaternary system above
        tmpExpectedSMILESList.add("CCCCCCCCC"); //nonane chain
        tmpExpectedSMILESList.add("C"); //methane connecting triple bond and quaternary system
        tmpExpectedSMILESList.add("C#C"); //triple bond
        tmpExpectedSMILESList.add("C"); //methane residual of triple bond chain
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setFragmentSideChainsSetting(true);
        tmpASF.setMaxChainLengthSetting(3);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("c1ccc(cc1)C(C)C"); //ring structure with connected tertiary carbon system
        tmpExpectedSMILESList.add("C=Cc1ccccc1"); //benzene with connected double bond forming conjugated pi bond system
        tmpExpectedSMILESList.add("C"); //methane between the substructures above
        tmpExpectedSMILESList.add("C=C1CCCC2CC12"); //ring system with connected non-cyclic double bond
        tmpExpectedSMILESList.add("C"); //methane connecting ring system above and quaternary carbon below
        tmpExpectedSMILESList.add("*C(*)(*)*"); //quaternary system in pseudoatom representation
        tmpExpectedSMILESList.add("C"); //methane residual of quaternary system above
        tmpExpectedSMILESList.add("CCC"); //propane residual of fragmented nonane chain
        tmpExpectedSMILESList.add("CCC"); //propane residual of fragmented nonane chain
        tmpExpectedSMILESList.add("CCC"); //propane residual of fragmented nonane chain
        tmpExpectedSMILESList.add("C"); //methane connecting triple bond and quaternary system
        tmpExpectedSMILESList.add("C#C"); //triple bond
        tmpExpectedSMILESList.add("C"); //methane residual of triple bond chain
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setMaxChainLengthSetting(6); //max chain length default value
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("*C(*)(*)*");
        tmpExpectedSMILESList.add("*C(*)*");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C=C1CCCC2CC12");
        tmpExpectedSMILESList.add("c1ccccc1");
        tmpExpectedSMILESList.add("C=Cc1ccccc1");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C#C");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("CCCCCC");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(false);
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)(C)C");
        tmpExpectedSMILESList.add("C=C1CCCC2CC12");
        tmpExpectedSMILESList.add("c1ccc(cc1)C(C)C");
        tmpExpectedSMILESList.add("C=Cc1ccccc1");
        tmpExpectedSMILESList.add("C#C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CC");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)(C)C");
        tmpExpectedSMILESList.add("C=C1CCCC2CC12");
        tmpExpectedSMILESList.add("c1ccccc1");
        tmpExpectedSMILESList.add("C=Cc1ccccc1");
        tmpExpectedSMILESList.add("*C(C)C");
        tmpExpectedSMILESList.add("C#C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    /**
     * Example fragmentation using natural compounds CNP0321029.1; CNP0509407.0 and CNP0242949.0. The natural compounds
     * used in this test can be found in the natural product database COCONUT <a href="https://coconut.naturalproducts.net/"<a/>.
     *
     * @throws Exception if molecule is not correctly fragmented
     */
    @Test
    public void naturalCompoundFragmentationTest() throws Exception {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //test structure: CNP0321029.1
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("C[C@@H]1C[C@@H](C)[C@@H](C)C1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>(5);
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C1CCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        //test structure: CNP0509407.0
        tmpTestStructureAC = tmpParser.parseSmiles("CCCCCCCCCC(CC)(CCCCCCC)CCCCCCCC");
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        tmpExpectedSMILESList.clear();
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCC");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("*C(*)(*)*");
        System.out.println("Result: " + tmpFragmentsSMILESList);
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //alternate settings
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpExpectedSMILESList.clear();
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("CC");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CC(C)(C)C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpASF.setFragmentSideChainsSetting(false);
        tmpExpectedSMILESList.clear();
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("CCCCCCCC");
        tmpExpectedSMILESList.add("CCCCCCC");
        tmpExpectedSMILESList.add("CCCCCC");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("CC(C)(C)C");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        //
        //test structure: CNP0242949.0
        tmpTestStructureAC = tmpParser.parseSmiles("CC1CCCCCCCCC1");
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        tmpExpectedSMILESList.clear();
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C1CCCCCCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Null Atom Input">
    /**
     * Test for the handling of molecules containing null atoms (formerly a test of the inner class MolecularArrays,
     * which skipped null elements). Since the CDK does not allow the addition of null atoms to existing atom container
     * instances, reflection is used in this test to manipulate the inner fields of the atom container. Such molecules
     * are always filtered, independent of the setting for keeping non-fragmentable molecules.
     *
     * @throws InvalidSmilesException if parser cannot parse the given SMILES
     * @throws NoSuchFieldException if either the atoms or bonds field can not be found
     * @throws IllegalAccessException if reflection can not be done
     */
    @Test
    public void fillMolecularArraysTest() throws InvalidSmilesException, NoSuchFieldException, IllegalAccessException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCCCCCCCCCCC");
        //reflection is used to manually insert/inflict a null value
        Field atomsField = AtomContainer.class.getDeclaredField("atoms");
        Field bondsField = AtomContainer.class.getDeclaredField("bonds");
        atomsField.setAccessible(true);
        bondsField.setAccessible(true);
        //"getter" for BaseAtomRef and BondRef
        Object atoms = atomsField.get(tmpMolecule);
        Object[] newAtoms = (Object[]) Array.newInstance(atoms.getClass().getComponentType(), Array.getLength(atoms));
        Object bonds = bondsField.get(tmpMolecule);
        Object[] newBonds = (Object[]) Array.newInstance(bonds.getClass().getComponentType(), Array.getLength(bonds));
        //fill newAtoms with actual atoms
        for (int i = 0; i < 12; i++) {
               newAtoms[i] = tmpMolecule.getAtom(i);
               if (i < 11) {
                   newBonds[i] = tmpMolecule.getBond(i);
               }
        }
        //set position 5 null
        newAtoms[5] = null;
        newBonds[5] = null;
        //insert the fields back into tmpMolecule
        atomsField.set(tmpMolecule, newAtoms);
        bondsField.set(tmpMolecule, newBonds);
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (boolean tmpKeepNonFragmentable : new boolean[]{false, true}) {
            tmpASF.setKeepNonFragmentableMoleculesSetting(tmpKeepNonFragmentable);
            Assertions.assertTrue(tmpASF.shouldBeFiltered(tmpMolecule));
            Assertions.assertFalse(tmpASF.canBeFragmented(tmpMolecule));
        }
    }

    /**
     * Test to ensure correct detection and handling of allowed/not-allowed molecules.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed by SmilesParser
     */
    @Test
    public void filterTest() throws InvalidSmilesException{
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> tmpTestStructuresACList = new ArrayList<>(4);
        tmpTestStructuresACList.add(tmpParser.parseSmiles("CCC"));
        tmpTestStructuresACList.add(tmpParser.parseSmiles("CN"));
        tmpTestStructuresACList.add(tmpParser.parseSmiles("C*"));
        tmpTestStructuresACList.add(tmpParser.parseSmiles("[H]C([H])([H])[H]"));
        boolean[] tmpFilterArray = new boolean[tmpTestStructuresACList.size()];
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setKeepNonFragmentableMoleculesSetting(true);
        for (int i = 0; i < tmpFilterArray.length; i++) {
            IAtomContainer tmpAC = tmpTestStructuresACList.get(i);
            tmpFilterArray[i] = tmpASF.shouldBeFiltered(tmpAC);
        }
        boolean[] tmpExpectedArray = new boolean[tmpTestStructuresACList.size()];
        //add filter cases with expected value
        tmpExpectedArray[0] = false;
        //since keepNon-FragmentableInPipeline == true, values for "CN" and "C*" are false to allow passing through filtering
        // to allow them to be kept in pipeline
        tmpExpectedArray[1] = false; //"CN"
        tmpExpectedArray[2] = false; //"C*"
        tmpExpectedArray[3] = false; //explicit hydrogen
        for (int i = 0; i < tmpFilterArray.length; i++) {
            if (tmpFilterArray[i] != tmpExpectedArray[i]) {
                Assertions.fail();
            }
        }
        //test for excluding non-fragmentables from pipeline
        tmpASF.setKeepNonFragmentableMoleculesSetting(false);
        for (int i = 0; i < tmpFilterArray.length; i++) {
            IAtomContainer tmpAC = tmpTestStructuresACList.get(i);
            tmpFilterArray[i] = tmpASF.shouldBeFiltered(tmpAC);
        }
        tmpExpectedArray[0] = false;
        tmpExpectedArray[1] = true;
        tmpExpectedArray[2] = true;
        tmpExpectedArray[3] = false;
        for (int i = 0; i < tmpFilterArray.length; i++) {
            if (tmpFilterArray[i] != tmpExpectedArray[i]) {
                Assertions.fail();
            }
        }

    }
    @Test
    public void disconnectedStructureInputTest() throws InvalidSmilesException, CloneNotSupportedException{
        //test structure: CC(C)C.C1CCCCC1C1CCCCC1
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpTestStructureAC = tmpParser.parseSmiles("CC(C)C.C1CCCCC1C1CCCCC1");
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        this.preprocessTestMolecule(tmpASF, tmpTestStructureAC,
                false, false, true);
        int tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        List<IAtomContainer> tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        List<String> tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        List<String> tmpExpectedSMILESList = new ArrayList<>();
        tmpExpectedSMILESList.add("*C(*)*");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C");
        tmpExpectedSMILESList.add("C1CCCCC1");
        tmpExpectedSMILESList.add("C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpFragmentsSMILESList.clear();
        tmpPreFragmenationCount = this.countAtoms(tmpTestStructureAC);
        tmpACList.clear();
        tmpACList = tmpASF.fragmentMolecule(tmpTestStructureAC);
        tmpFragmentsSMILESList = new ArrayList<>();
        if (this.checkChemicalFormula(tmpPreFragmenationCount, tmpACList)) {
            tmpFragmentsSMILESList.addAll(this.generateSMILESFromACList(tmpACList, false));
        }
        tmpExpectedSMILESList.clear();
        tmpExpectedSMILESList.add("CC(C)C");
        tmpExpectedSMILESList.add("C1CCCCC1");
        tmpExpectedSMILESList.add("C1CCCCC1");
        Assertions.assertTrue(this.compareListsIgnoringOrder(new ArrayList<>(tmpFragmentsSMILESList),
                new ArrayList<>(tmpExpectedSMILESList)));
    }
    //</editor-fold
    //
    //<editor-fold desc="Private Utility Methods">
    /**
     * Compares two provided lists for equality while ignoring the lists' orders.
     * The second given list will be empty after the comparison if both lists contain equal objects.
     *
     * @param aList1 First given list to compare
     * @param aList2 Second given list to compare
     * @return boolean whether given lists are equal
     */
    private boolean compareListsIgnoringOrder(List<?> aList1, List<?> aList2) {
        if (aList1 == null || aList2 == null) {
            return false;
        }
        if (aList1.size() != aList2.size()) {
            return false;
        }
        for (Object o : aList1) {
            aList2.remove(o);
        }
        return aList2.isEmpty();
    }
    /**
     * Utility method to generate SMILES notation strings for a given List containing IAtomContainer instances.
     *
     * @param anACList given List of IAtomContainers
     * @return List with SMILES Strings of originally given List of IAtomContainers
     */
    private List<String> generateSMILESFromACList(List<IAtomContainer> anACList, boolean isStereoChemEncoded) {
        List<String> tmpReturnSmilesList = new ArrayList<>(anACList.size());
        for (IAtomContainer tmpAC: anACList) {
            tmpReturnSmilesList.add(ChemUtil.createUniqueSmiles(tmpAC, isStereoChemEncoded, true));
        }
        return tmpReturnSmilesList;
    }
    /**
     * Private utility method responsible for calling preprocessing methods of a given AlkylStructureFragmenter instance.
     *
     * @param anAlkylStructureFragmenterInstance the given ASF instance of which to call preprocessing
     * @param aMolecule to preprocess
     * @param aShouldBeFilteredStatement for shouldBeFiltered()
     * @param aShouldBePreprocessedStatement for shouldBePreprocessed()
     * @param aCanBeFragmentedStatement for canBeFragmented()
     */
    private void preprocessTestMolecule(AlkylStructureFragmenter anAlkylStructureFragmenterInstance,
                                        IAtomContainer aMolecule,
                                        boolean aShouldBeFilteredStatement,
                                        boolean aShouldBePreprocessedStatement,
                                        boolean aCanBeFragmentedStatement) {
        Assertions.assertEquals(aShouldBeFilteredStatement,
                anAlkylStructureFragmenterInstance.shouldBeFiltered(aMolecule));
        Assertions.assertEquals(aShouldBePreprocessedStatement,
                anAlkylStructureFragmenterInstance.shouldBePreprocessed(aMolecule));
        Assertions.assertEquals(aCanBeFragmentedStatement,
                anAlkylStructureFragmenterInstance.canBeFragmented(aMolecule));
    }
    /**
     * Asserts that exactly the atoms with the given indices are classified as true in the given classification array.
     *
     * @param aClassificationArray classification array indexed by atom index
     * @param anExpectedAtomIndices indices of the atoms expected to be classified as true
     */
    private void assertAtoms(boolean[] aClassificationArray, int... anExpectedAtomIndices) {
        Set<Integer> tmpExpectedIndices = new HashSet<>();
        for (int tmpIndex : anExpectedAtomIndices) {
            tmpExpectedIndices.add(tmpIndex);
        }
        for (int i = 0; i < aClassificationArray.length; i++) {
            Assertions.assertEquals(tmpExpectedIndices.contains(i), aClassificationArray[i], String.format("Classification of atom %d", i));
        }
    }
    /**
     * Asserts that exactly the bonds between the given atom index pairs are classified as true in the given
     * classification array.
     *
     * @param aMolecule the classified molecule
     * @param aClassificationArray classification array indexed by bond index
     * @param anExpectedBonds atom index pairs of the bonds expected to be classified as true
     */
    private void assertBonds(IAtomContainer aMolecule, boolean[] aClassificationArray, int[][] anExpectedBonds) {
        Set<Integer> tmpExpectedBondIndices = new HashSet<>();
        for (int[] tmpAtomPair : anExpectedBonds) {
            IBond tmpBond = aMolecule.getBond(aMolecule.getAtom(tmpAtomPair[0]), aMolecule.getAtom(tmpAtomPair[1]));
            Assertions.assertNotNull(tmpBond, "No bond between atoms " + tmpAtomPair[0] + " and " + tmpAtomPair[1]);
            tmpExpectedBondIndices.add(aMolecule.indexOf(tmpBond));
        }
        for (int i = 0; i < aClassificationArray.length; i++) {
            IBond tmpBond = aMolecule.getBond(i);
            Assertions.assertEquals(tmpExpectedBondIndices.contains(i), aClassificationArray[i],
                    String.format("Classification of bond %d-%d", aMolecule.indexOf(tmpBond.getBegin()), aMolecule.indexOf(tmpBond.getEnd())));
        }
    }
    /**
     * Counts the atoms of the given atom container that are no pseudo atoms.
     *
     * @param anAtomContainer atom container to count the atoms of
     * @return number of non-pseudo atoms
     */
    private int countAtoms(IAtomContainer anAtomContainer) {
        int tmpCount = 0;
        for (IAtom tmpAtom : anAtomContainer.atoms()) {
            if (!(tmpAtom instanceof IPseudoAtom) && tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() != 0) {
                tmpCount++;
            }
        }
        return tmpCount;
    }

    /**
     * Utility method checking the persistence of chemical formula before and after fragmentation by comparing their respective atom counts.
     *
     * @param aPreFragAtomCount from the molecule before fragmentation
     * @param anAtomContainerList with fragments of a fragmented molecule
     * @return boolean for persistence (true if persistent)
     */
    private boolean checkChemicalFormula(int aPreFragAtomCount, List<IAtomContainer> anAtomContainerList) {
        Objects.requireNonNull(anAtomContainerList);
        int tmpPostFragmentationAtomCount = 0;
        for (IAtomContainer tmpAC : anAtomContainerList) {
            tmpPostFragmentationAtomCount += this.countAtoms(tmpAC);
        }
        return aPreFragAtomCount == tmpPostFragmentationAtomCount;
    }
    //</editor-fold>
}
