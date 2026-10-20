/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2025  Felix Baensch, Jonas Schaub (felix.j.baensch@gmail.com, jonas.schaub@uni-jena.de)
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
public class AlkylStructureFragmenterTest extends AlkylStructureFragmenter {
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
     * Method for unit testing the internal methods AlkylStructureFragmenter.markSubstructures() and getFragmentationResults().
     * The molecule used for testing is a conceptual molecule containing several key structural features detected and
     * extracted by the AlkylStructureFragmenter - an isolated double bond, a quaternary C, an aromatic ring system with
     * a connected conjugated system, an isolated cyclohexane ring and a tertiary C - in default settings.
     *
     * @throws CDKException if SMILES cannot be parsed correctly
     */
    @Test
    public void extractFragmentsTest() throws CDKException {
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
        MolecularArrays tmpTestMolecularArrays = new MolecularArrays(tmpTestAC);
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.markSubstructures(tmpTestMolecularArrays, tmpTestAC);
        IAtomContainerSet tmpExtractedFrag;
        ArrayList<String> tmpActualFragmentList = new ArrayList<>(12);
        tmpExtractedFrag = tmpASF.getFragmentationResults(tmpTestMolecularArrays);
        for (IAtomContainer tmpAC : tmpExtractedFrag.atomContainers()) {
            ChemUtil.saturateWithHydrogen(tmpAC);
            tmpActualFragmentList.add(ChemUtil.createUniqueSmiles(tmpAC, false));
        }
        Assertions.assertTrue(this.compareListsIgnoringOrder(tmpActualFragmentList, tmpExpectedFragmentList));
    }
    //<editor-fold desc="Mark Tests">
    /**
     * Test for correct marking of tertiary and quaternary carbon atoms and their neighbor atoms and bonds.
     * Atom types are deliberately not perceived, the marking has to work on the bonds of the molecule alone.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkTertQuatAndNeighbors() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //isobutane: atom 1 is tertiary, atoms 0, 2, 3 are its neighbors
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC(C)C");
        this.markTertQuatAndNeighbors(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY, 1);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_QUATERNARY_CARBON_PROPERTY_KEY);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, 0, 2, 3);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, new int[][]{{0, 1}, {1, 2}, {1, 3}});
        //neopentyl chain: atom 1 is quaternary, atom 5 is not a neighbor
        tmpMolecule = tmpParser.parseSmiles("CC(C)(C)CC");
        this.markTertQuatAndNeighbors(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_QUATERNARY_CARBON_PROPERTY_KEY, 1);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, 0, 2, 3, 4);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, new int[][]{{0, 1}, {1, 2}, {1, 3}, {1, 4}});
        //branching atom with a double bond (isobutene) is no tertiary carbon
        tmpMolecule = tmpParser.parseSmiles("CC(C)=C");
        this.markTertQuatAndNeighbors(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY);
    }

    /**
     * Test for correct function of routine marking ring substructures.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkRings() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //methylcyclopropane: methyl atom 0 and its bond are not part of the ring
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1CC1");
        this.markRings(new MolecularArrays(tmpMolecule), tmpMolecule);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_RING_MARKER_KEY, 1, 2, 3);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_RING_MARKER_KEY, new int[][]{{1, 2}, {2, 3}, {3, 1}});
    }

    /**
     * Test for correct function of routine marking multi-bonds, double and triple bonds have to be marked with their
     * respective markers.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkMultiBonds() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC=CC#C");
        this.markMultiBonds(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_DOUBLE_BOND_MARKER_KEY, 1, 2);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_DOUBLE_BOND_MARKER_KEY, new int[][]{{1, 2}});
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TRIPLE_BOND_MARKER_KEY, 3, 4);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TRIPLE_BOND_MARKER_KEY, new int[][]{{3, 4}});
    }

    /**
     * Test for correct marking of tertiary/quaternary carbon atoms directly bonded to a ring and the connecting bond.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void testMarkConnectedTertQuatRing() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //isopropylcyclopropane: tertiary atom 1 is bonded to ring atom 3
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC(C)C1CC1");
        MolecularArrays tmpMolecularArrays = new MolecularArrays(tmpMolecule);
        this.markRings(tmpMolecularArrays, tmpMolecule);
        this.markTertQuatAndNeighbors(tmpMolecularArrays);
        this.markConnectedTertQuatRing(tmpMolecularArrays);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY, 1);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONNECTED_TERTIARY_QUATERNARY_RING_MARKER_KEY, 1, 3);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONNECTED_TERTIARY_QUATERNARY_RING_MARKER_KEY, new int[][]{{1, 3}});
        //ring atoms themselves are never tertiary carbons, no connection marked for methylcyclopropane
        tmpMolecule = tmpParser.parseSmiles("CC1CC1");
        tmpMolecularArrays = new MolecularArrays(tmpMolecule);
        this.markRings(tmpMolecularArrays, tmpMolecule);
        this.markTertQuatAndNeighbors(tmpMolecularArrays);
        this.markConnectedTertQuatRing(tmpMolecularArrays);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY);
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONNECTED_TERTIARY_QUATERNARY_RING_MARKER_KEY);
    }

    /**
     * Method to test the internal algorithm for detecting and marking conjugated pi bond systems.
     * Hexatriene has to be marked completely. Double and triple bonds separated by a single bond (D-S-T) are not
     * regarded as conjugated, also not if the double bond is part of an aromatic ring (phenylacetylene).
     *
     * @throws InvalidSmilesException if SMILES-parser is not able to parse SMILES due to incorrect syntax
     */
    @Test
    public void testMarkConjugatedPiSystems() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C=CC=CC=C");
        this.markConjugatedPiSystems(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY, 0, 1, 2, 3, 4, 5);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY,
                new int[][]{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
        //D-S-T is not conjugated
        tmpMolecule = tmpParser.parseSmiles("C=CC#C");
        this.markConjugatedPiSystems(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY);
        //phenylacetylene: only the ring is conjugated, regardless of the order in which the ring bonds are visited
        tmpMolecule = tmpParser.parseSmiles("C#CC1=CC=CC=C1");
        this.markConjugatedPiSystems(new MolecularArrays(tmpMolecule));
        this.assertAtomMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY, 2, 3, 4, 5, 6, 7);
        this.assertBondMarker(tmpMolecule, AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY,
                new int[][]{{2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 2}});
    }
    //</editor-fold>
    /**
     * Test for correct deepCopy methods by copying a butene molecule which used to make problems in earlier versions.
     * A seemingly unused MolecularArrays instance is instanced here. This is needed for correct filling of internal atom
     * and bond arrays, since they are created and filled in the MolecularArrays constructor.
     *
     * @throws InvalidSmilesException if SMILES parser fails to parse SMILES from String with incorrect syntax
     */
    @Test
    public void deepCopyButeneTest() throws InvalidSmilesException {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpButeneContainer = tmpParser.parseSmiles("C=CCC");
        //two steps below needed for correct internal index handling
        MolecularArrays tmpMolecularArraysInstance = new MolecularArrays(tmpButeneContainer);
        IAtomContainer tmpCopyAC = tmpButeneContainer.getBuilder().newAtomContainer();
        for (IAtom tmpAtom: tmpButeneContainer.atoms()) {
            tmpCopyAC.addAtom(tmpASF.deepCopyAtom(tmpAtom));
        }
        for (IBond tmpBond: tmpButeneContainer.bonds()) {
            tmpCopyAC.addBond(deepCopyBond(tmpBond, tmpCopyAC));
        }
        //Comparison of original and copied AtomContainer
        String tmpButeneSMILES = ChemUtil.createUniqueSmiles(tmpButeneContainer, false);
        String tmpCopySMILES = ChemUtil.createUniqueSmiles(tmpCopyAC, false);
        Assertions.assertEquals(tmpButeneSMILES, tmpCopySMILES);
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
    //<editor-fold desc="Inner Class 'Molecular Arrays'">
    /**
     * Test method for correct functionality of AlkylStructureFragmenter inner class MolecularArrays.
     * It tests the correct filtering of null elements during the array filling step.
     * It also tests correct setting of atom and bond properties.
     * Since the CDK does not allow the addition of null atoms to existing atom container instances, reflection is used
     * in this test to access the inner fields of atoms and bonds to manipulate their values.
     *
     * @throws InvalidSmilesException if parser cannot parse the given SMILES
     * @throws NoSuchFieldException if either the atoms or bonds field can not be found
     * @throws IllegalAccessException if reflection can not be done
     */
    @Test
    public void fillMolecularArraysTest() throws InvalidSmilesException, NoSuchFieldException, IllegalAccessException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCCCCCCCCCCC");
        MolecularArrays refMolecularArrays = new MolecularArrays(tmpMolecule);
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
        //test for expected discrepancy
        MolecularArrays testMolecularArrays = new MolecularArrays(tmpMolecule);
        Assertions.assertNotEquals(testMolecularArrays.getAtomArray().length, refMolecularArrays.getAtomArray().length);
        Assertions.assertNotEquals(testMolecularArrays.getBondArray().length, refMolecularArrays.getBondArray().length);
        if (testMolecularArrays.getAtomArray().length != 11 && testMolecularArrays.getBondArray().length != 10) {
            Assertions.fail("Array lengths were not in expected range.");
        }
        for (int i = 0; i < testMolecularArrays.getAtomArray().length; i++) {
            Map<Object, Object> tmpExpectedAtomPropertyMap = new HashMap<>();
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY, i);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_RING_MARKER_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_TERTIARY_CARBON_PROPERTY_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_QUATERNARY_CARBON_PROPERTY_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_DOUBLE_BOND_MARKER_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_TRIPLE_BOND_MARKER_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, false);
            tmpExpectedAtomPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_CONNECTED_TERTIARY_QUATERNARY_RING_MARKER_KEY, false);
            IAtom tmpAtom = testMolecularArrays.getAtomArray()[i];
            Map<Object, Object> tmpTestAtomPropertyMap = tmpAtom.getProperties();
            for (Map.Entry<Object, Object> entry : tmpTestAtomPropertyMap.entrySet()) {
                if (!tmpExpectedAtomPropertyMap.containsKey(entry.getKey()) || !tmpExpectedAtomPropertyMap.containsValue(entry.getValue())) {
                    Assertions.fail("Properties not set correctly!");
                }
            }
        }
        for (int i = 0; i < testMolecularArrays.getBondArray().length; i++) {
            Map<Object, Object> tmpExpectedBondPropertyMap = new HashMap<>();
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_BOND_INDEX_PROPERTY_KEY, i);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_RING_MARKER_KEY, false);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_CONJ_PI_MARKER_KEY, false);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_DOUBLE_BOND_MARKER_KEY, false);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_TRIPLE_BOND_MARKER_KEY, false);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_NEIGHBOR_MARKER_KEY, false);
            tmpExpectedBondPropertyMap.put(AlkylStructureFragmenter.INTERNAL_ASF_CONNECTED_TERTIARY_QUATERNARY_RING_MARKER_KEY, false);
            IBond tmpBond = testMolecularArrays.getBondArray()[i];
            Map<Object, Object> tmpTestAtomPropertyMap = tmpBond.getProperties();
            for (Map.Entry<Object, Object> entry : tmpTestAtomPropertyMap.entrySet()) {
                if (!tmpExpectedBondPropertyMap.containsKey(entry.getKey()) || !tmpExpectedBondPropertyMap.containsValue(entry.getValue())) {
                    Assertions.fail("Properties not set correctly!");
                }
            }
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
        setKeepNonFragmentableMoleculesSetting(true);
        for (int i = 0; i < tmpFilterArray.length; i++) {
            IAtomContainer tmpAC = tmpTestStructuresACList.get(i);
            tmpFilterArray[i] = shouldBeFiltered(tmpAC);
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
        setKeepNonFragmentableMoleculesSetting(false);
        for (int i = 0; i < tmpFilterArray.length; i++) {
            IAtomContainer tmpAC = tmpTestStructuresACList.get(i);
            tmpFilterArray[i] = shouldBeFiltered(tmpAC);
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
     * Asserts that exactly the atoms with the given indices carry the given boolean marker property.
     *
     * @param aMolecule the marked molecule
     * @param aMarkerKey property key of the marker
     * @param anExpectedMarkedAtomIndices indices of the atoms expected to be marked
     */
    private void assertAtomMarker(IAtomContainer aMolecule, String aMarkerKey, int... anExpectedMarkedAtomIndices) {
        Set<Integer> tmpExpectedIndices = new HashSet<>();
        for (int tmpIndex : anExpectedMarkedAtomIndices) {
            tmpExpectedIndices.add(tmpIndex);
        }
        for (int i = 0; i < aMolecule.getAtomCount(); i++) {
            Assertions.assertEquals(tmpExpectedIndices.contains(i), aMolecule.getAtom(i).getProperty(aMarkerKey),
                    String.format("Marker %s of atom %d", aMarkerKey, i));
        }
    }
    /**
     * Asserts that exactly the bonds between the given atom index pairs carry the given boolean marker property.
     *
     * @param aMolecule the marked molecule
     * @param aMarkerKey property key of the marker
     * @param anExpectedMarkedBonds atom index pairs of the bonds expected to be marked
     */
    private void assertBondMarker(IAtomContainer aMolecule, String aMarkerKey, int[][] anExpectedMarkedBonds) {
        Set<IBond> tmpExpectedBonds = new HashSet<>();
        for (int[] tmpAtomPair : anExpectedMarkedBonds) {
            IBond tmpBond = aMolecule.getBond(aMolecule.getAtom(tmpAtomPair[0]), aMolecule.getAtom(tmpAtomPair[1]));
            Assertions.assertNotNull(tmpBond, "No bond between atoms " + tmpAtomPair[0] + " and " + tmpAtomPair[1]);
            tmpExpectedBonds.add(tmpBond);
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            Assertions.assertEquals(tmpExpectedBonds.contains(tmpBond), tmpBond.getProperty(aMarkerKey),
                    String.format("Marker %s of bond %d-%d", aMarkerKey, aMolecule.indexOf(tmpBond.getBegin()),
                            aMolecule.indexOf(tmpBond.getEnd())));
        }
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
