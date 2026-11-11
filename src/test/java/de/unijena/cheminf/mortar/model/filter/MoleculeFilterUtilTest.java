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

package de.unijena.cheminf.mortar.model.filter;

import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupType;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests the structural filter functions of {@link MoleculeFilterUtil}: atomic numbers and heavy atom count, the element
 * filter, connectivity, the selection of the largest component including its tie-break, and the preservation of
 * container properties, coordinates, charges, isotopes, flags, Sgroups, and stereo elements in the reduced molecule.
 * It also pins what CDK's ConnectivityChecker.partitionIntoMolecules() preserves and what it drops in CDK 2.13.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class MoleculeFilterUtilTest {
    //<editor-fold desc="Static initializer">
    /**
     * Sets the default locale to British English.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * SD record of sodium 2-hydroxy(5-13C)pent-3-enoate with 2D coordinates, a wedge bond defining the stereocentre,
     * an (E)-configured double bond (from the coordinates), a carbon-13 label, charges, a title, and two data fields.
     * Generated with CDK's structure diagram generator and MDL V2000 writer from
     * [O-]C(=O)[C@@H](O)/C=C/[13CH3].[Na+].
     */
    private static final String SALT_SDF_RECORD = """
            test salt title
              CDK     11112619092D

              9  7  0  0  1  0  0  0  0  0999 V2000
                2.6986    3.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
                3.9976    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                3.9976    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
                5.2966    3.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
                5.2966    4.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
                6.5957    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                7.8947    3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                9.1937    2.2500    0.0000 C   1  0  0  0  0  0  0  0  0  0  0  0
                0.7500    1.8750    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0
              1  2  1  0  0  0  0
              2  3  2  0  0  0  0
              2  4  1  0  0  0  0
              4  5  1  1  0  0  0
              4  6  1  0  0  0  0
              6  7  2  0  0  0  0
              7  8  1  0  0  0  0
            M  CHG  1   1  -1
            M  CHG  1   9   1
            M  ISO  1   8  13
            M  END
            > <DB_ID>
            TEST0001

            > <Comment>
            sodium salt with stereocentre, E double bond, and 13C label

            $$$$
            """;

    /**
     * The anion of {@link #SALT_SDF_RECORD} alone, i.e. the same MOL block without the sodium atom.
     */
    private static final String ANION_SDF_RECORD = """
            test anion title
              CDK     11112619092D

              8  7  0  0  1  0  0  0  0  0999 V2000
                2.6986    3.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
                3.9976    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                3.9976    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
                5.2966    3.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
                5.2966    4.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
                6.5957    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                7.8947    3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                9.1937    2.2500    0.0000 C   1  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0  0  0  0
              2  3  2  0  0  0  0
              2  4  1  0  0  0  0
              4  5  1  1  0  0  0
              4  6  1  0  0  0  0
              6  7  2  0  0  0  0
              7  8  1  0  0  0  0
            M  CHG  1   1  -1
            M  ISO  1   8  13
            M  END
            $$$$
            """;

    /**
     * SD record of a 1:1 mixture of toluene and phenol (seven heavy atoms each) as Kekule structures without
     * aromaticity flags, written by hand without coordinates; the double bonds of phenol are shifted compared to
     * toluene.
     */
    private static final String TOLUENE_PHENOL_SDF_RECORD = """
            toluene and phenol
              MORTAR test

             14 14  0  0  0  0  0  0  0  0999 V2000
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  2  0  0  0  0
              2  3  1  0  0  0  0
              3  4  2  0  0  0  0
              4  5  1  0  0  0  0
              5  6  2  0  0  0  0
              6  1  1  0  0  0  0
              1  7  1  0  0  0  0
              8  9  1  0  0  0  0
              9 10  2  0  0  0  0
             10 11  1  0  0  0  0
             11 12  2  0  0  0  0
             12 13  1  0  0  0  0
             13  8  2  0  0  0  0
              8 14  1  0  0  0  0
            M  END
            $$$$
            """;

    /**
     * SD record of toluene with the methyl carbon atom first, a permuted order of the ring atoms, and double bonds at
     * other positions than in a Kekule SMILES string written ring first (same as "toluene_permuted" in the importer's
     * test file MoleculeFiltersDuplicatesTestFile.sdf).
     */
    private static final String TOLUENE_PERMUTED_SDF_RECORD = """
            toluene_permuted
              MORTAR test

              7  7  0  0  0  0  0  0  0  0999 V2000
                2.5714    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -1.9286    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.3214    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.3214   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -1.1786    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                1.0714   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -1.1786   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
              1  6  1  0
              6  4  1  0
              4  7  2  0
              7  2  1  0
              2  5  2  0
              5  3  1  0
              3  6  2  0
            M  END
            $$$$
            """;
    //</editor-fold>
    //
    //<editor-fold desc="Tests">
    /**
     * Tests atomic numbers and the heavy atom count: hydrogen atoms (explicit, implicit, deuterium) and pseudo atoms do
     * not count. Uses unprocessed SMILES parser output for the explicit hydrogen atoms because the importer's
     * preprocessing converts non-isotopic explicit hydrogen atoms into implicit ones.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testAtomicNumbersAndHeavyAtomCount() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        Assertions.assertEquals(3, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("CCO")));
        Assertions.assertEquals(3, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("[H]C([H])([H])C([H])([H])O[H]")));
        Assertions.assertEquals(2, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("[2H]C([2H])([2H])O")));
        Assertions.assertEquals(3, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("*CCO")));
        Assertions.assertEquals(0, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("[H+]")));
        Assertions.assertEquals(0, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("[H][H]")));
        Assertions.assertEquals(1, MoleculeFilterUtil.countHeavyAtoms(tmpSmilesParser.parseSmiles("[Na+]")));
        Assertions.assertEquals(0, MoleculeFilterUtil.countHeavyAtoms(SilentChemObjectBuilder.getInstance().newAtomContainer()));
        //the same after the importer's preprocessing (deuterium stays explicit)
        List<IAtomContainer> tmpMolecules = new PreprocessingTestImporter().parseAndPreprocess(
                "[H]C([H])([H])C([H])([H])O[H]", "[2H]C([2H])([2H])O");
        Assertions.assertEquals(3, MoleculeFilterUtil.countHeavyAtoms(tmpMolecules.get(0)));
        Assertions.assertEquals(5, tmpMolecules.get(1).getAtomCount());
        Assertions.assertEquals(2, MoleculeFilterUtil.countHeavyAtoms(tmpMolecules.get(1)));
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        IPseudoAtom tmpPseudoAtom = tmpBuilder.newInstance(IPseudoAtom.class, "R1");
        Assertions.assertEquals(0, MoleculeFilterUtil.getAtomicNumber(tmpPseudoAtom));
        Assertions.assertFalse(MoleculeFilterUtil.isHeavyAtom(tmpPseudoAtom));
        //an unset atomic number means an unknown element (CDK derives the symbol from the atomic number)
        IAtom tmpUnknownAtom = tmpBuilder.newInstance(IAtom.class, "C");
        tmpUnknownAtom.setAtomicNumber(null);
        Assertions.assertNull(tmpUnknownAtom.getSymbol());
        Assertions.assertEquals(0, MoleculeFilterUtil.getAtomicNumber(tmpUnknownAtom));
        Assertions.assertFalse(MoleculeFilterUtil.isHeavyAtom(tmpUnknownAtom));
        Assertions.assertFalse(MoleculeFilterUtil.isAllowedElement(tmpUnknownAtom));
        Assertions.assertEquals(6, MoleculeFilterUtil.getAtomicNumber(tmpBuilder.newInstance(IAtom.class, "C")));
        Assertions.assertThrows(NullPointerException.class, () -> MoleculeFilterUtil.countHeavyAtoms(null));
    }

    /**
     * Tests the element filter: exactly the 13 allowed elements pass (checked for all elements 1 to 118), isotopes of
     * them pass, pseudo atoms, metals, and metalloids other than B, Si, and Se fail.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testElementFilter() throws Exception {
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        List<String> tmpAllowedSymbols = new ArrayList<>();
        for (int tmpAtomicNumber = 1; tmpAtomicNumber <= 118; tmpAtomicNumber++) {
            String tmpSymbol = Elements.ofNumber(tmpAtomicNumber).symbol();
            if (MoleculeFilterUtil.isAllowedElement(tmpBuilder.newInstance(IAtom.class, tmpSymbol))) {
                tmpAllowedSymbols.add(tmpSymbol);
            }
        }
        Assertions.assertEquals(List.of("H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"), tmpAllowedSymbols);
        Assertions.assertFalse(MoleculeFilterUtil.isAllowedElement(tmpBuilder.newInstance(IPseudoAtom.class, "R")));
        SmilesParser tmpSmilesParser = new SmilesParser(tmpBuilder);
        for (String tmpAllowedSmiles : new String[] {"[2H]C([2H])([2H])O", "[3H]O", "[13CH3]C(=O)O", "[15NH3]", "C[18OH]",
                "[H]C([H])([H])O[H]", "OB(O)c1ccccc1", "C[Si](C)(C)C", "C[Se]C", "FC(Cl)(Br)I", "OP(=O)(O)O",
                "CS(=O)(=O)O", ""}) {
            Assertions.assertTrue(MoleculeFilterUtil.consistsOfAllowedElementsOnly(tmpSmilesParser.parseSmiles(tmpAllowedSmiles)),
                    tmpAllowedSmiles);
        }
        for (String tmpDisallowedSmiles : new String[] {"*CCO", "[Na+]", "CC(=O)[O-].[Na+]", "Cl[Pt](Cl)(N)N",
                "[Fe+2]", "C[Hg]C", "CC[Sn](CC)(CC)CC", "C[As](C)C", "C[Te]C", "[Ge](C)(C)(C)C", "[He]", "[Al+3]",
                "[Li]C", "[K+].[Cl-]"}) {
            Assertions.assertFalse(MoleculeFilterUtil.consistsOfAllowedElementsOnly(tmpSmilesParser.parseSmiles(tmpDisallowedSmiles)),
                    tmpDisallowedSmiles);
        }
    }

    /**
     * Tests the connectivity check; every unconnected atom is an own component, molecules with fewer than two atoms
     * are connected.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testConnectivity() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        Assertions.assertTrue(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("CCO")));
        Assertions.assertTrue(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("[Na+]")));
        Assertions.assertTrue(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("[H][H]")));
        Assertions.assertTrue(MoleculeFilterUtil.isConnected(SilentChemObjectBuilder.getInstance().newAtomContainer()));
        Assertions.assertFalse(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("CCO.O")));
        Assertions.assertFalse(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("[H+].[Cl-]")));
        Assertions.assertFalse(MoleculeFilterUtil.isConnected(tmpSmilesParser.parseSmiles("*.C")));
    }

    /**
     * Tests the selection of the largest component: salts, hydrates, connected and trivial molecules (returned as
     * they are), components without heavy atoms, the documented limitation that a larger counter-ion wins, and the
     * tie-break that does not depend on the order of the components in the input.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testLargestComponentSelection() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        this.assertLargestComponent(tmpImporter, "CC(=O)[O-].[Na+]", "CC(=O)[O-]");
        this.assertLargestComponent(tmpImporter, "[Na+].CC(=O)[O-]", "CC(=O)[O-]");
        this.assertLargestComponent(tmpImporter, "O.O.CCO.O", "CCO");
        this.assertLargestComponent(tmpImporter, "[Na+].[Na+].[O-]C(=O)C(=O)[O-]", "[O-]C(=O)C(=O)[O-]");
        this.assertLargestComponent(tmpImporter, "C[NH3+].[Cl-]", "C[NH3+]");
        //pseudo atoms are no heavy atoms
        this.assertLargestComponent(tmpImporter, "*.*.*.*.C", "C");
        //limitation: the size criterion picks the larger counter-ion (tosylate, 11 heavy atoms) over the cation of
        // interest (tetramethylammonium, 5 heavy atoms)
        this.assertLargestComponent(tmpImporter, "C[N+](C)(C)C.Cc1ccc(cc1)S(=O)(=O)[O-]", "Cc1ccc(cc1)S(=O)(=O)[O-]");
        //ties are broken independently of the input order: 1:1 mixture, ion pair, enantiomers
        this.assertLargestComponent(tmpImporter, "c1ccccc1.C1CCCCC1", "C1CCCCC1");
        this.assertLargestComponent(tmpImporter, "C1CCCCC1.c1ccccc1", "C1CCCCC1");
        this.assertLargestComponent(tmpImporter, "[Na+].[Cl-]", "[Cl-]");
        this.assertLargestComponent(tmpImporter, "[Cl-].[Na+]", "[Cl-]");
        String tmpEnantiomerOne = this.uniqueStereoSmiles(MoleculeFilterUtil.extractLargestComponent(
                tmpImporter.parseAndPreprocess("C[C@H](N)O.C[C@@H](N)O").getFirst()));
        String tmpEnantiomerTwo = this.uniqueStereoSmiles(MoleculeFilterUtil.extractLargestComponent(
                tmpImporter.parseAndPreprocess("C[C@@H](N)O.C[C@H](N)O").getFirst()));
        Assertions.assertEquals(tmpEnantiomerOne, tmpEnantiomerTwo);
        Assertions.assertTrue(tmpEnantiomerOne.contains("@"));
        //identical components: the first one in atom order
        IAtomContainer tmpTwoEthanol = tmpImporter.parseAndPreprocess("CCO.CCO").getFirst();
        IAtomContainer tmpFirstEthanol = MoleculeFilterUtil.extractLargestComponent(tmpTwoEthanol);
        Assertions.assertEquals(3, tmpFirstEthanol.getAtomCount());
        Assertions.assertSame(AtomRef.deref(tmpTwoEthanol.getAtom(0)), AtomRef.deref(tmpFirstEthanol.getAtom(0)));
        //connected, single-atom, and empty molecules are returned as they are
        for (String tmpSmiles : new String[] {"CCO", "[Na+]", ""}) {
            IAtomContainer tmpMolecule = tmpImporter.parseAndPreprocess(tmpSmiles).getFirst();
            Assertions.assertSame(tmpMolecule, MoleculeFilterUtil.extractLargestComponent(tmpMolecule));
        }
    }

    /**
     * Pins the behaviour of ConnectivityChecker.partitionIntoMolecules() in CDK 2.13 that the largest component
     * extraction builds on: for a disconnected molecule, the components are new containers without any container
     * properties (no title, no MORTAR name, no SDF data fields), but with the stereo elements of the component,
     * re-wrapped for the new container, and with the same atom objects as the input; for a connected molecule, the
     * input container itself is returned. If this test fails after a CDK update, check whether
     * {@link MoleculeFilterUtil#extractLargestComponent(IAtomContainer)} still copies everything that is dropped.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testPartitionIntoMoleculesBehaviourOfCdk() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        IAtomContainer tmpSalt = tmpImporter.readAndPreprocessSdfRecord(MoleculeFilterUtilTest.SALT_SDF_RECORD, "salt");
        Assertions.assertEquals(2, this.countStereoElements(tmpSalt));
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpSalt);
        Assertions.assertEquals(2, tmpComponents.getAtomContainerCount());
        IAtomContainer tmpAnion = tmpComponents.getAtomContainer(0);
        Assertions.assertEquals(8, tmpAnion.getAtomCount());
        //dropped: all container properties
        Assertions.assertTrue(tmpAnion.getProperties().isEmpty());
        Assertions.assertNull(tmpAnion.getTitle());
        //kept: stereo elements, re-wrapped for the component
        Assertions.assertEquals(2, this.countStereoElements(tmpAnion));
        for (IStereoElement<?, ?> tmpElement : tmpAnion.stereoElements()) {
            IChemObject tmpFocus = tmpElement.getFocus();
            if (tmpFocus instanceof IAtom tmpAtom) {
                Assertions.assertSame(tmpAnion, tmpAtom.getContainer());
            } else {
                Assertions.assertSame(tmpAnion, ((IBond) tmpFocus).getContainer());
            }
        }
        //shared: the atom objects
        Assertions.assertSame(AtomRef.deref(tmpSalt.getAtom(0)), AtomRef.deref(tmpAnion.getAtom(0)));
        //connected input: the same container
        IAtomContainer tmpEthanol = tmpImporter.parseAndPreprocess("CCO").getFirst();
        Assertions.assertSame(tmpEthanol, ConnectivityChecker.partitionIntoMolecules(tmpEthanol).getAtomContainer(0));
    }

    /**
     * Tests that the largest component of an SD record keeps the title, the MORTAR name, the SDF data fields, charges,
     * the isotope label, 2D coordinates, implicit hydrogen counts, and both stereo elements, and that its stereo
     * SMILES equals the one of the same anion imported alone.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testReducedMoleculeKeepsDataOfSdfRecord() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        IAtomContainer tmpSalt = tmpImporter.readAndPreprocessSdfRecord(MoleculeFilterUtilTest.SALT_SDF_RECORD, "salt name");
        IAtomContainer tmpAnionAlone = tmpImporter.readAndPreprocessSdfRecord(MoleculeFilterUtilTest.ANION_SDF_RECORD, "anion name");
        String tmpSaltStereoSmiles = this.uniqueStereoSmiles(tmpSalt);
        IAtomContainer tmpReduced = MoleculeFilterUtil.extractLargestComponent(tmpSalt);
        Assertions.assertNotSame(tmpSalt, tmpReduced);
        Assertions.assertEquals(8, tmpReduced.getAtomCount());
        Assertions.assertEquals(7, tmpReduced.getBondCount());
        //container properties
        Assertions.assertEquals("test salt title", tmpReduced.getTitle());
        Assertions.assertEquals("salt name", tmpReduced.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY));
        Assertions.assertEquals("TEST0001", tmpReduced.getProperty("DB_ID"));
        Assertions.assertEquals("sodium salt with stereocentre, E double bond, and 13C label", tmpReduced.getProperty("Comment"));
        Assertions.assertEquals(tmpSalt.getProperties(), tmpReduced.getProperties());
        //atoms: charges, isotope, coordinates, implicit hydrogen counts, in the same order as the anion alone
        for (int i = 0; i < tmpReduced.getAtomCount(); i++) {
            IAtom tmpAtom = tmpReduced.getAtom(i);
            IAtom tmpExpectedAtom = tmpAnionAlone.getAtom(i);
            Assertions.assertEquals(tmpExpectedAtom.getSymbol(), tmpAtom.getSymbol());
            Assertions.assertEquals(tmpExpectedAtom.getFormalCharge(), tmpAtom.getFormalCharge());
            Assertions.assertEquals(tmpExpectedAtom.getMassNumber(), tmpAtom.getMassNumber());
            Assertions.assertEquals(tmpExpectedAtom.getImplicitHydrogenCount(), tmpAtom.getImplicitHydrogenCount());
            Assertions.assertEquals(tmpExpectedAtom.getPoint2d(), tmpAtom.getPoint2d());
        }
        Assertions.assertEquals(-1, tmpReduced.getAtom(0).getFormalCharge());
        Assertions.assertEquals(13, tmpReduced.getAtom(7).getMassNumber());
        Assertions.assertEquals(new Point2d(9.1937, 2.25), tmpReduced.getAtom(7).getPoint2d());
        //stereo elements: one tetrahedral centre and one double bond, belonging to the reduced container
        int tmpTetrahedralCount = 0;
        int tmpDoubleBondCount = 0;
        for (IStereoElement<?, ?> tmpElement : tmpReduced.stereoElements()) {
            if (tmpElement instanceof ITetrahedralChirality tmpTetrahedral) {
                tmpTetrahedralCount++;
                Assertions.assertSame(tmpReduced, tmpTetrahedral.getFocus().getContainer());
            } else if (tmpElement instanceof IDoubleBondStereochemistry tmpDoubleBond) {
                tmpDoubleBondCount++;
                Assertions.assertSame(tmpReduced, tmpDoubleBond.getFocus().getContainer());
            }
        }
        Assertions.assertEquals(1, tmpTetrahedralCount);
        Assertions.assertEquals(1, tmpDoubleBondCount);
        //stereo SMILES equal to the anion imported alone and to the anion part of the salt's SMILES
        String tmpReducedStereoSmiles = this.uniqueStereoSmiles(tmpReduced);
        Assertions.assertEquals("C/C=C/[C@@H](C(=O)[O-])O", tmpReducedStereoSmiles);
        Assertions.assertEquals(this.uniqueStereoSmiles(tmpAnionAlone), tmpReducedStereoSmiles);
        Assertions.assertEquals("C/C=C/[C@@H](C(=O)[O-])O.[Na+]", tmpSaltStereoSmiles);
    }

    /**
     * Tests the reduction for a molecule read like in the importer's MOL file import, i.e. with the MDL V2000 reader in
     * relaxed mode into an {@link org.openscience.cdk.AtomContainer} of the cdk-data module instead of the silent one.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testReducedMoleculeFromMolFileReader() throws Exception {
        String tmpMolBlock = MoleculeFilterUtilTest.SALT_SDF_RECORD.substring(0,
                MoleculeFilterUtilTest.SALT_SDF_RECORD.indexOf("M  END") + "M  END".length()) + "\n";
        IAtomContainer tmpSalt;
        try (MDLV2000Reader tmpReader = new MDLV2000Reader(new StringReader(tmpMolBlock), IChemObjectReader.Mode.RELAXED)) {
            tmpSalt = tmpReader.read(new org.openscience.cdk.AtomContainer());
        }
        tmpSalt.setProperty(Importer.MOLECULE_NAME_PROPERTY_KEY, "mol file salt");
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        tmpImporter.preprocess(tmpSalt);
        IAtomContainer tmpReduced = MoleculeFilterUtil.extractLargestComponent(tmpSalt);
        Assertions.assertNotSame(tmpSalt, tmpReduced);
        Assertions.assertEquals(8, tmpReduced.getAtomCount());
        Assertions.assertEquals("test salt title", tmpReduced.getTitle());
        Assertions.assertEquals("mol file salt", tmpReduced.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY));
        Assertions.assertEquals(2, this.countStereoElements(tmpReduced));
        Assertions.assertEquals("C/C=C/[C@@H](C(=O)[O-])O", this.uniqueStereoSmiles(tmpReduced));
    }

    /**
     * Tests that the largest component keeps aromaticity flags of atoms and bonds, atom properties, 3D coordinates,
     * the container ID and flags, and only the Sgroups that lie completely in the component.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testReducedMoleculeKeepsFlagsCoordinatesAndSgroups() throws Exception {
        List<IAtomContainer> tmpMolecules = new PreprocessingTestImporter().parseAndPreprocess(
                "c1ccccc1C(=O)[O-].[Na+]", "c1ccccc1C(=O)[O-]");
        IAtomContainer tmpSalt = tmpMolecules.get(0);
        for (int i = 0; i < tmpSalt.getAtomCount(); i++) {
            tmpSalt.getAtom(i).setPoint3d(new Point3d(i, 2.0 * i, 3.0 * i));
        }
        tmpSalt.getAtom(0).setProperty("TEST_ATOM_PROPERTY", "ring atom");
        tmpSalt.setID("ID42");
        tmpSalt.set(IChemObject.VISITED);
        Sgroup tmpCarboxylateSgroup = new Sgroup();
        tmpCarboxylateSgroup.setType(SgroupType.CtabAbbreviation);
        tmpCarboxylateSgroup.setSubscript("COO");
        for (int i = 6; i <= 8; i++) {
            tmpCarboxylateSgroup.addAtom(tmpSalt.getAtom(i));
        }
        Sgroup tmpSodiumSgroup = new Sgroup();
        tmpSodiumSgroup.setType(SgroupType.CtabAbbreviation);
        tmpSodiumSgroup.setSubscript("Na");
        tmpSodiumSgroup.addAtom(tmpSalt.getAtom(9));
        List<Sgroup> tmpSgroups = new ArrayList<>(2);
        tmpSgroups.add(tmpCarboxylateSgroup);
        tmpSgroups.add(tmpSodiumSgroup);
        tmpSalt.setProperty(CDKConstants.CTAB_SGROUPS, tmpSgroups);
        IAtomContainer tmpReduced = MoleculeFilterUtil.extractLargestComponent(tmpSalt);
        Assertions.assertEquals(9, tmpReduced.getAtomCount());
        for (int i = 0; i < tmpReduced.getAtomCount(); i++) {
            IAtom tmpAtom = tmpReduced.getAtom(i);
            Assertions.assertEquals(new Point3d(i, 2.0 * i, 3.0 * i), tmpAtom.getPoint3d());
            Assertions.assertEquals(i < 6, tmpAtom.isAromatic());
        }
        int tmpAromaticBondCount = 0;
        for (IBond tmpBond : tmpReduced.bonds()) {
            if (tmpBond.isAromatic()) {
                tmpAromaticBondCount++;
            }
        }
        Assertions.assertEquals(6, tmpAromaticBondCount);
        Assertions.assertEquals("ring atom", tmpReduced.getAtom(0).getProperty("TEST_ATOM_PROPERTY"));
        Assertions.assertEquals("ID42", tmpReduced.getID());
        Assertions.assertTrue(tmpReduced.is(IChemObject.VISITED));
        List<Sgroup> tmpReducedSgroups = tmpReduced.getProperty(CDKConstants.CTAB_SGROUPS);
        Assertions.assertEquals(1, tmpReducedSgroups.size());
        Assertions.assertEquals("COO", tmpReducedSgroups.getFirst().getSubscript());
        for (IAtom tmpSgroupAtom : tmpReducedSgroups.getFirst().getAtoms()) {
            Assertions.assertTrue(tmpReduced.contains(tmpSgroupAtom));
        }
        Assertions.assertEquals(this.uniqueStereoSmiles(tmpMolecules.get(1)), this.uniqueStereoSmiles(tmpReduced));
    }

    /**
     * Tests that the tie-break of the largest component does not depend on the notation of the input: a 1:1 mixture of
     * toluene and phenol (seven heavy atoms each) as aromatic SMILES, as Kekule SMILES with different double bond
     * positions, in both orders of the components, and as SD record without aromaticity flags always gives toluene,
     * whose canonical Kekule SMILES is the lexicographically smaller one. Also tests the fallback for components whose
     * SMILES string cannot be generated because no Kekule structure could be assigned (cyclopentadienyl without
     * hydrogen): the other component wins in both orders, and of two such components, the first in atom order.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testLargestComponentTieBreakIsIndependentOfNotation() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        String tmpTolueneKekuleSmiles = this.uniqueKekuleStereoSmiles(tmpImporter.parseAndPreprocess("Cc1ccccc1").getFirst());
        List<IAtomContainer> tmpMixtures = new ArrayList<>(tmpImporter.parseAndPreprocess("Cc1ccccc1.Oc1ccccc1",
                "Oc1ccccc1.Cc1ccccc1", "CC1=CC=CC=C1.OC1=CC=CC=C1", "OC1C=CC=CC=1.CC1C=CC=CC=1"));
        tmpMixtures.add(tmpImporter.readAndPreprocessSdfRecord(MoleculeFilterUtilTest.TOLUENE_PHENOL_SDF_RECORD, "mixture"));
        for (IAtomContainer tmpMixture : tmpMixtures) {
            IAtomContainer tmpLargestComponent = MoleculeFilterUtil.extractLargestComponent(tmpMixture);
            Assertions.assertEquals(7, tmpLargestComponent.getAtomCount());
            Assertions.assertEquals(tmpTolueneKekuleSmiles, this.uniqueKekuleStereoSmiles(tmpLargestComponent));
        }
        //fallback: no SMILES string for cyclopentadienyl, which cannot be kekulized
        IAtomContainer tmpCyclopentadienyl = tmpImporter.parseAndPreprocess("c1cccc1").getFirst();
        Assertions.assertNull(this.uniqueKekuleStereoSmiles(tmpCyclopentadienyl));
        for (String tmpSmiles : new String[] {"c1cccc1.CCCCC", "CCCCC.c1cccc1"}) {
            IAtomContainer tmpLargestComponent = MoleculeFilterUtil.extractLargestComponent(
                    tmpImporter.parseAndPreprocess(tmpSmiles).getFirst());
            Assertions.assertEquals("CCCCC", this.uniqueKekuleStereoSmiles(tmpLargestComponent), tmpSmiles);
        }
        IAtomContainer tmpTwoRings = tmpImporter.parseAndPreprocess("c1cccc1.c1cccc1").getFirst();
        IAtomContainer tmpFirstRing = MoleculeFilterUtil.extractLargestComponent(tmpTwoRings);
        Assertions.assertEquals(5, tmpFirstRing.getAtomCount());
        Assertions.assertSame(AtomRef.deref(tmpTwoRings.getAtom(0)), AtomRef.deref(tmpFirstRing.getAtom(0)));
    }

    /**
     * Tests the duplicate key: one compound in aromatic and different Kekule notations (also from an SD record with
     * another atom order and double bond positions) gives one key, while the stored unique SMILES strings with aromatic
     * symbols differ; different compounds give different keys; stereoisomers differ only with stereochemistry; a stored
     * SMILES string without aromaticity encoding is reused as key; the stored SMILES string is the fallback if no
     * Kekule SMILES string can be generated.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDuplicateKey() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        String[][] tmpNotationGroups = {
                {"Cc1ccccc1", "CC1=CC=CC=C1", "CC1C=CC=CC=1", "c1cc(C)ccc1"},
                {"c1ccncc1", "C1=CC=NC=C1", "N1=CC=CC=C1"},
                {"c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1", "C1=CC2=CC=CC=C2C=C1"},
                {"CC(=O)Nc1ccc(O)cc1", "CC(=O)NC1=CC=C(O)C=C1", "OC1=CC=C(NC(C)=O)C=C1"}};
        Set<String> tmpAllKeys = new HashSet<>();
        for (String[] tmpGroup : tmpNotationGroups) {
            Set<String> tmpKeys = new HashSet<>();
            Set<String> tmpStoredSmiles = new HashSet<>();
            for (IAtomContainer tmpMolecule : tmpImporter.parseAndPreprocess(tmpGroup)) {
                String tmpStored = this.uniqueStereoSmiles(tmpMolecule);
                tmpStoredSmiles.add(tmpStored);
                tmpKeys.add(MoleculeFilterUtil.createDuplicateKey(tmpMolecule, true, tmpStored, true));
            }
            Assertions.assertEquals(1, tmpKeys.size(), tmpKeys.toString());
            Assertions.assertTrue(tmpStoredSmiles.size() > 1, tmpStoredSmiles.toString());
            tmpAllKeys.addAll(tmpKeys);
        }
        Assertions.assertEquals(tmpNotationGroups.length, tmpAllKeys.size());
        //SD record of toluene
        IAtomContainer tmpToluene = tmpImporter.parseAndPreprocess("Cc1ccccc1").getFirst();
        String tmpTolueneKey = MoleculeFilterUtil.createDuplicateKey(tmpToluene, true, this.uniqueStereoSmiles(tmpToluene), true);
        IAtomContainer tmpSdfToluene = tmpImporter.readAndPreprocessSdfRecord(MoleculeFilterUtilTest.TOLUENE_PERMUTED_SDF_RECORD, "toluene");
        Assertions.assertEquals(tmpTolueneKey, MoleculeFilterUtil.createDuplicateKey(tmpSdfToluene, true,
                this.uniqueStereoSmiles(tmpSdfToluene), true));
        //stereoisomers
        List<IAtomContainer> tmpAlanines = tmpImporter.parseAndPreprocess("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O");
        Assertions.assertNotEquals(MoleculeFilterUtil.createDuplicateKey(tmpAlanines.get(0), true, "x", true),
                MoleculeFilterUtil.createDuplicateKey(tmpAlanines.get(1), true, "x", true));
        Assertions.assertEquals(MoleculeFilterUtil.createDuplicateKey(tmpAlanines.get(0), false, "x", true),
                MoleculeFilterUtil.createDuplicateKey(tmpAlanines.get(1), false, "x", true));
        //reuse of a stored SMILES string without aromaticity encoding
        String tmpStoredKekuleSmiles = this.uniqueKekuleStereoSmiles(tmpToluene);
        Assertions.assertEquals(tmpTolueneKey, tmpStoredKekuleSmiles);
        Assertions.assertSame(tmpStoredKekuleSmiles, MoleculeFilterUtil.createDuplicateKey(tmpToluene, true, tmpStoredKekuleSmiles, false));
        //fallback
        IAtomContainer tmpCyclopentadienyl = tmpImporter.parseAndPreprocess("c1cccc1").getFirst();
        String tmpStoredAromaticSmiles = this.uniqueStereoSmiles(tmpCyclopentadienyl);
        Assertions.assertNotNull(tmpStoredAromaticSmiles);
        Assertions.assertSame(tmpStoredAromaticSmiles, MoleculeFilterUtil.createDuplicateKey(tmpCyclopentadienyl, true,
                tmpStoredAromaticSmiles, true));
        //argument checks
        Assertions.assertThrows(NullPointerException.class, () -> MoleculeFilterUtil.createDuplicateKey(null, true, "C", true));
        Assertions.assertThrows(NullPointerException.class, () -> MoleculeFilterUtil.createDuplicateKey(tmpToluene, true, null, true));
        Assertions.assertThrows(IllegalArgumentException.class, () -> MoleculeFilterUtil.createDuplicateKey(tmpToluene, true, " ", true));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Asserts that the largest component of the given SMILES, after the importer's preprocessing, has the same unique
     * stereo SMILES as the given expected component, preprocessed the same way.
     *
     * @param anImporter preprocessing helper
     * @param aSmiles SMILES of the molecule
     * @param anExpectedComponentSmiles SMILES of the expected largest component
     * @throws Exception if anything goes wrong
     */
    private void assertLargestComponent(PreprocessingTestImporter anImporter, String aSmiles, String anExpectedComponentSmiles)
            throws Exception {
        List<IAtomContainer> tmpMolecules = anImporter.parseAndPreprocess(aSmiles, anExpectedComponentSmiles);
        IAtomContainer tmpLargestComponent = MoleculeFilterUtil.extractLargestComponent(tmpMolecules.get(0));
        Assertions.assertTrue(MoleculeFilterUtil.isConnected(tmpLargestComponent));
        Assertions.assertEquals(this.uniqueStereoSmiles(tmpMolecules.get(1)), this.uniqueStereoSmiles(tmpLargestComponent), aSmiles);
    }

    /**
     * Returns MORTAR's unique SMILES with stereochemistry and aromaticity encoded (the importer's default).
     *
     * @param aMolecule the molecule
     * @return unique SMILES
     */
    private String uniqueStereoSmiles(IAtomContainer aMolecule) {
        return ChemUtil.createUniqueSmiles(aMolecule, true, true);
    }

    /**
     * Returns MORTAR's unique SMILES with stereochemistry and without aromaticity encoding, i.e. in the normalised Kekule
     * notation of the duplicate key.
     *
     * @param aMolecule the molecule
     * @return unique SMILES or null if it cannot be generated
     */
    private String uniqueKekuleStereoSmiles(IAtomContainer aMolecule) {
        return ChemUtil.createUniqueSmiles(aMolecule, true, false);
    }

    /**
     * Counts the stereo elements of the given molecule.
     *
     * @param aMolecule the molecule
     * @return number of stereo elements
     */
    private int countStereoElements(IAtomContainer aMolecule) {
        int tmpCount = 0;
        for (IStereoElement<?, ?> tmpElement : aMolecule.stereoElements()) {
            tmpCount++;
        }
        return tmpCount;
    }
    //</editor-fold>
}
