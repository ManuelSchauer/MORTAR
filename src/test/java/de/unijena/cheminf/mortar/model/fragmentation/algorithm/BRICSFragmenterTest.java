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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests the {@link BRICSFragmenter}. The expected fragments of the examples were derived by hand from the original BRICS
 * rules (dissertation of J. Degen, appendix C.1.1.3), the rule that cleaves a bond is given in the comments. Fragments
 * are compared as canonical SMILES with the link types encoded as mass numbers of the link atoms, e.g. [16*]c1ccccc1.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class BRICSFragmenterTest {
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public BRICSFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the settings, the rule set, copying, and the restoration of the default settings.
     */
    @Test
    public void basicSettingsTest() {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        Assertions.assertEquals(BRICSFragmenter.ALGORITHM_NAME, tmpFragmenter.getFragmentationAlgorithmName());
        Assertions.assertFalse(tmpFragmenter.getFragmentationAlgorithmDisplayName().isBlank());
        Assertions.assertEquals(2, tmpFragmenter.settingsProperties().size());
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
        }
        Assertions.assertEquals(BRICSFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(BRICSFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT,
                tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        tmpFragmenter.setSaturationSetting(BRICSFragmenter.SaturationOption.HYDROGEN_SATURATION);
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        BRICSFragmenter tmpCopy = (BRICSFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(BRICSFragmenter.SaturationOption.HYDROGEN_SATURATION, tmpCopy.getSaturationSetting());
        Assertions.assertFalse(tmpCopy.getReturnUnfragmentedMoleculesSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(BRICSFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setSaturationSetting(null));
        //all 55 rules of the original rule set with unique names
        List<String> tmpRuleNames = BRICSFragmenter.getRuleNames();
        Assertions.assertEquals(55, tmpRuleNames.size());
        Assertions.assertEquals(55, new HashSet<>(tmpRuleNames).size());
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(null));
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
    }
    //
    /**
     * Tests the fragmentation of sorafenib, the example molecule of the original publication. The urea is cleaved at
     * both amide bonds (N1, L1 - L2) and the diaryl ether at both ether bonds (O1d, L16 - L3). The aryl-nitrogen bonds
     * of the urea, the N-methyl amide, and the trifluoromethyl group are not cleaved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void sorafenibTest() throws Exception {
        this.assertFragments("CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",
                "[1*]C([1*])=O", "[16*]c1ccnc(C(=O)NC)c1", "[2*]Nc1ccc(Cl)c(C(F)(F)F)c1", "[2*]Nc1ccc([16*])cc1",
                "[3*]O[3*]");
    }
    //
    /**
     * Tests the fragmentation of imatinib: amide (N1), diarylamine (N2b, L14 - L2 and N2d, L16 - L2), biaryl (C3c,
     * L14 - L16), and the benzylic amine (N2, L4 - L5). The benzyl carbon atom stays at the ring, because an alkyl
     * carbon atom of type L8 needs a further carbon neighbour.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void imatinibTest() throws Exception {
        this.assertFragments("Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
                "[1*]C(=O)c1ccc(C[4*])cc1", "[14*]c1ccnc([14*])n1", "[16*]c1cccnc1", "[2*]N[2*]",
                "[2*]Nc1ccc(C)c([16*])c1", "[5*]N1CCN(C)CC1");
    }
    //
    /**
     * Tests that bonds to small terminal groups are not cleaved: hydrogen, halogens, hydroxy, nitro, carboxylate,
     * methoxy, methyl, ethyl, propyl, isopropyl, and acetyl groups, methyl and ethyl esters, and N-methyl and N-ethyl
     * amides. Ring bonds are never cleaved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void smallTerminalGroupsAreNotCleavedTest() throws Exception {
        String[] tmpUnfragmentedMolecules = {
                "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1", "CC(C)c1ccccc1", "Clc1ccccc1", "FC(F)(F)c1ccccc1",
                "Oc1ccccc1", "COc1ccccc1", "CCOc1ccccc1", "O=[N+]([O-])c1ccccc1", "OC(=O)c1ccccc1", "CC(=O)c1ccccc1",
                "CCC(=O)c1ccccc1", "COC(=O)c1ccccc1", "CCOC(=O)c1ccccc1", "CC(=O)Nc1ccccc1", "CNC(=O)c1ccccc1",
                "CCNC(=O)c1ccccc1", "Nc1ccccc1", "NS(=O)(=O)c1ccccc1", "CSc1ccccc1", "C1CCC2CCCCC2C1", "c1ccc2ccccc2c1",
                "CCCCCCCC", "OCCCCO"};
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        for (String tmpSmiles : tmpUnfragmentedMolecules) {
            IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(tmpSmiles);
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            Assertions.assertEquals(1, tmpFragments.size(), tmpSmiles);
            Assertions.assertEquals(BRICSFragmenterTest.canonical(tmpSmiles, false),
                    BRICSFragmenterTest.toCanonicalSmiles(tmpFragments.getFirst(), false), tmpSmiles);
        }
    }
    //
    /**
     * Tests examples for the different rule groups of the original rule set.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleGroupsTest() throws Exception {
        //C1d: alkyl chain at a ring, L8 - L16
        this.assertFragments("CCCCc1ccccc1", "[16*]c1ccccc1", "[8*]CCCC");
        //C1d: bridging methylene group, L8 - L16
        this.assertFragments("c1ccccc1Cc1ccccc1", "[16*]c1ccccc1", "[16*]c1ccccc1", "[8*]C[8*]");
        //C3c: biaryl, L14 - L16, C5: L16 - L16
        this.assertFragments("c1ccccc1-c1ccccn1", "[16*]c1ccccc1", "[14*]c1ccccn1");
        this.assertFragments("c1ccccc1-c1ccccc1", "[16*]c1ccccc1", "[16*]c1ccccc1");
        //C4a: two aliphatic rings, L15 - L15 (not cleaved by the RDKit implementation)
        this.assertFragments("C1CCC(CC1)C1CCCCC1", "[15*]C1CCCCC1", "[15*]C1CCCCC1");
        //C2a: two activated aliphatic rings, L13 - L13
        this.assertFragments("C1CCOC(C1)C1CCCCO1", "[13*]C1CCCCO1", "[13*]C1CCCCO1");
        //C6d: diaryl ketone, L6 - L16
        this.assertFragments("O=C(c1ccccc1)c1ccccc1", "[16*]c1ccccc1", "[16*]c1ccccc1", "[6*]C([6*])=O");
        //C7: olefin, L7 = L7
        this.assertFragments("CCCCC=CCCCC", "[7*]=CCCCC", "[7*]=CCCCC");
        //no cleavage of Michael acceptors, neither the olefin nor the ester
        this.assertFragments("CCCCC=CC(=O)OCCC", "CCCCC=CC(=O)OCCC");
        //O1: dialkyl ether, L4 - L3
        this.assertFragments("CCCCOCCCC", "[3*]O[3*]", "[4*]CCCC", "[4*]CCCC");
        //O1d: diaryl ether, L16 - L3
        this.assertFragments("c1ccccc1Oc1ccccc1", "[16*]c1ccccc1", "[16*]c1ccccc1", "[3*]O[3*]");
        //O1: aryl alkyl ether, only the alkyl side is cleaved
        this.assertFragments("CCCCOc1ccccc1", "[3*]Oc1ccccc1", "[4*]CCCC");
        //O2: ester, L1 - L3
        this.assertFragments("CCCOC(=O)c1ccccc1", "[1*]C(=O)c1ccccc1", "[3*]OCCC");
        //O3: phosphate ester, untyped link on the phosphorus side
        this.assertFragments("CCCCOP(=O)(O)O", "[3*]OCCCC", "*P(=O)(O)O");
        //N1: amide, L1 - L2
        this.assertFragments("CCCNC(=O)c1ccccc1", "[1*]C(=O)c1ccccc1", "[2*]NCCC");
        //N1: tertiary amide with a ring nitrogen atom
        this.assertFragments("O=C(c1ccccc1)N1CCCCC1", "[1*]C(=O)c1ccccc1", "[2*]N1CCCCC1");
        //N3: sulfonamide, L12 - L2
        this.assertFragments("CCCNS(=O)(=O)c1ccccc1", "[12*]S(=O)(=O)c1ccccc1", "[2*]NCCC");
        //N6: acylated aromatic nitrogen atom, L1 - L9
        this.assertFragments("O=C(c1ccccc1)n1cccc1", "[1*]C(=O)c1ccccc1", "[9*]n1cccc1");
        //N2: dialkyl amine, L4 - L5
        this.assertFragments("CCCCNCCCC", "[4*]CCCC", "[4*]CCCC", "[5*]N[5*]");
        //N2c: aliphatic ring carbon atom and ring amine, L15 - L5
        this.assertFragments("C1CCC(CC1)N1CCCCC1", "[15*]C1CCCCC1", "[5*]N1CCCCC1");
        //N2d: aryl piperidine, L16 - L2
        this.assertFragments("c1ccccc1N1CCCCC1", "[16*]c1ccccc1", "[2*]N1CCCCC1");
        //N4: alkylated pyrrole, L8 - L9, N5d: arylated pyrrole, L16 - L9
        this.assertFragments("CCCCn1cccc1", "[8*]CCCC", "[9*]n1cccc1");
        this.assertFragments("c1ccccc1n1cccc1", "[16*]c1ccccc1", "[9*]n1cccc1");
        //no pyrrole-like cleavage of pyridones
        this.assertFragments("O=c1ccccn1CCCC", "[10*]n1ccccc1=O", "[8*]CCCC");
        //N7: alkylated lactam, L8 - L10, N8d: arylated lactam, L16 - L10, N9: acylated lactam, L1 - L10
        this.assertFragments("CCCCN1CCCC1=O", "[8*]CCCC", "[10*]N1CCCC1=O");
        this.assertFragments("O=C1CCCN1c1ccccc1", "[10*]N1CCCC1=O", "[16*]c1ccccc1");
        this.assertFragments("CCCC(=O)N1CCCC1=O", "[1*]C(=O)CCC", "[10*]N1CCCC1=O");
        //S1: dialkyl thioether, L4 - L11
        this.assertFragments("CCCCSCCCC", "[11*]S[11*]", "[4*]CCCC", "[4*]CCCC");
        //S2: dialkyl sulfoxide, L4 - L11
        this.assertFragments("CCCCS(=O)CCCC", "[11*]S([11*])=O", "[4*]CCCC", "[4*]CCCC");
        //S2d: diaryl sulfone, L16 - L11
        this.assertFragments("c1ccccc1S(=O)(=O)c1ccccc1", "[11*]S([11*])(=O)=O", "[16*]c1ccccc1", "[16*]c1ccccc1");
        //S2: aryl alkyl sulfone, only the alkyl side is cleaved (L4 - L11)
        this.assertFragments("CCCCS(=O)(=O)c1ccccc1", "[11*]S(=O)(=O)c1ccccc1", "[4*]CCCC");
    }
    //
    /**
     * Tests the saturation of the open valences with hydrogen atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hydrogenSaturationTest() throws Exception {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        tmpFragmenter.setSaturationSetting(BRICSFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCCOCCCC", false, "O", "CCCC", "CCCC");
        //cleaved double bonds are saturated with two hydrogen atoms
        this.assertFragments(tmpFragmenter, "CCCCC=CCCCC", false, "CCCCC", "CCCCC");
        this.assertFragments(tmpFragmenter, "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1", false,
                "C=O", "CNC(=O)c1ccccn1", "Nc1ccc(Cl)c(C(F)(F)F)c1", "Nc1ccccc1", "O");
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("CCCCOCCCC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Assertions.assertFalse(tmpAtom instanceof IPseudoAtom);
                Assertions.assertNotEquals(1, tmpAtom.getAtomicNumber());
            }
        }
    }
    //
    /**
     * Tests the link atoms: pseudo atoms with the link type as mass number and atom property.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void linkAtomTest() throws Exception {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("CCCCOc1ccccc1"));
        Set<Integer> tmpLinkTypes = new HashSet<>();
        for (IAtomContainer tmpFragment : tmpFragments) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom) {
                    Integer tmpLinkType = tmpAtom.getProperty(BRICSFragmenter.LINK_TYPE_PROPERTY_KEY);
                    Assertions.assertEquals(tmpLinkType, tmpAtom.getMassNumber());
                    Assertions.assertEquals(0, tmpAtom.getImplicitHydrogenCount());
                    Assertions.assertEquals(1, tmpFragment.getConnectedBondsCount(tmpAtom));
                    tmpLinkTypes.add(tmpLinkType);
                }
            }
        }
        Assertions.assertEquals(Set.of(3, 4), tmpLinkTypes);
        //the untyped link of rule O3 has no mass number
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("CCCCOP(=O)(O)O"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom && tmpFragment.getConnectedAtomsList(tmpAtom).getFirst().getAtomicNumber() == 15) {
                    Assertions.assertNull(tmpAtom.getMassNumber());
                    Assertions.assertEquals(BRICSFragmenter.UNTYPED_LINK, (int) tmpAtom.getProperty(BRICSFragmenter.LINK_TYPE_PROPERTY_KEY));
                }
            }
        }
    }
    //
    /**
     * Tests the setting for molecules without cleavable bonds and the handling of multi-component molecules,
     * explicit hydrogen atoms, and pseudo atoms of preceding fragmentations in a pipeline.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void inputHandlingTest() throws Exception {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        Assertions.assertEquals(1, tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("Cc1ccccc1")).size());
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("Cc1ccccc1")).isEmpty());
        Assertions.assertEquals(3, tmpFragmenter.fragmentMolecule(BRICSFragmenterTest.parse("CCCCOCCCC")).size());
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(true);
        //counter ions are returned unchanged
        this.assertFragments("CCCCNCCCC.[Na+].[Cl-]", "[4*]CCCC", "[4*]CCCC", "[5*]N[5*]", "[Na+]", "[Cl-]");
        //explicit hydrogen atoms do not change the result
        IAtomContainer tmpMolecule = BRICSFragmenterTest.parse("CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1");
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        Assertions.assertEquals(
                BRICSFragmenterTest.canonical(List.of("[1*]C([1*])=O", "[16*]c1ccnc(C(=O)NC)c1",
                        "[2*]Nc1ccc(Cl)c(C(F)(F)F)c1", "[2*]Nc1ccc([16*])cc1", "[3*]O[3*]"), false),
                BRICSFragmenterTest.toCanonicalSmiles(new BRICSFragmenter().fragmentMolecule(tmpMolecule), false));
        //pseudo atoms of a preceding fragmentation count as heavy atom neighbours
        this.assertFragments("*CCCOCCCC", "[3*]O[3*]", "[4*]CCC*", "[4*]CCCC");
        //the input molecule is not altered
        IAtomContainer tmpInput = BRICSFragmenterTest.parse("CCCCOCCCC");
        new BRICSFragmenter().fragmentMolecule(tmpInput);
        Assertions.assertEquals(9, tmpInput.getAtomCount());
        Assertions.assertEquals(8, tmpInput.getBondCount());
    }
    //
    /**
     * Tests the transfer of tetrahedral stereo elements to the fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void tetrahedralStereoTest() throws Exception {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        //the ether oxygen atom of the stereocentre is replaced by the link atom at the same position
        this.assertFragments(tmpFragmenter, "CCCC[C@@H](C)OCCCC", true, "CCCC[C@@H](C)[4*]", "[3*]O[3*]", "[4*]CCCC");
        this.assertFragments(tmpFragmenter, "CCCC[C@H](C)OCCCC", true, "CCCC[C@H](C)[4*]", "[3*]O[3*]", "[4*]CCCC");
        //stereo elements of uncleaved parts are kept
        this.assertFragments(tmpFragmenter, "C[C@H](O)CCCCOCCCC", true, "C[C@H](O)CCCC[4*]", "[3*]O[3*]", "[4*]CCCC");
        //two link atoms at a former stereocentre, the stereo element is removed
        this.assertFragments(tmpFragmenter, "CCCCN[C@@H](CCCC)OCCCC", true, "[4*]C([4*])CCCC", "[3*]O[3*]",
                "[4*]CCCC", "[4*]CCCC", "[5*]N[5*]");
        //the pseudo-asymmetric centre C3 of a tropane keeps its configuration
        this.assertFragments(tmpFragmenter, "CN1[C@@H]2CC[C@H]1C[C@@H](C2)OC(=O)CCCC", true,
                "CN1[C@@H]2CC[C@H]1C[C@@H](C2)O[3*]", "[1*]C(=O)CCCC");
        this.assertFragments(tmpFragmenter, "CN1[C@@H]2CC[C@H]1C[C@H](C2)OC(=O)CCCC", true,
                "CN1[C@@H]2CC[C@H]1C[C@H](C2)O[3*]", "[1*]C(=O)CCCC");
        //quaternary carbon atoms are not cleaved (!D4)
        this.assertFragments(tmpFragmenter, "CCCC[C@@](CC)(C)OCCCC", true, "CCCC[C@@](CC)(C)O[3*]", "[4*]CCCC");
        //with hydrogen saturation, the former stereocentre has two hydrogen atoms
        tmpFragmenter.setSaturationSetting(BRICSFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCC[C@@H](C)OCCCC", true, "CCCCCC", "O", "CCCC");
        this.assertFragments(tmpFragmenter, "C[C@H](O)CCCCOCCCC", true, "C[C@H](O)CCCC", "O", "CCCC");
    }
    //
    /**
     * Tests the transfer of double bond stereo elements to the fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void doubleBondStereoTest() throws Exception {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
        //the amine nitrogen atom is a ligand of the stereo double bond (N2, L4 - L5); the double bond itself is not
        //cleaved (C7), because the methyl substituent is too small
        this.assertFragments(tmpFragmenter, "C/C=C(\\NCCCC)CCC", true, "C/C=C(\\[4*])CCC", "[5*]N[5*]", "[4*]CCCC");
        this.assertFragments(tmpFragmenter, "C/C=C(/NCCCC)CCC", true, "C/C=C(/[4*])CCC", "[5*]N[5*]", "[4*]CCCC");
        //cleaved stereo double bond (C7)
        this.assertFragments(tmpFragmenter, "CCCC/C=C/CCCC", true, "[7*]=CCCCC", "[7*]=CCCCC");
        //with hydrogen saturation, the configuration refers to the remaining substituent
        tmpFragmenter.setSaturationSetting(BRICSFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "C/C=C(\\NCCCC)CCC", true, "C/C=C/CCC", "N", "CCCC");
        this.assertFragments(tmpFragmenter, "C/C=C(/NCCCC)CCC", true, "C/C=C\\CCC", "N", "CCCC");
    }
    //
    /**
     * Fragments the given molecule with a fragmenter with default settings and compares the fragments to the given
     * expected fragments, without stereochemistry.
     *
     * @param aSmiles SMILES of the molecule
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(String aSmiles, String... anExpectedFragments) throws Exception {
        this.assertFragments(new BRICSFragmenter(), aSmiles, false, anExpectedFragments);
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param isStereo whether the stereochemistry is compared
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(BRICSFragmenter aFragmenter, String aSmiles, boolean isStereo,
                                 String... anExpectedFragments) throws Exception {
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(BRICSFragmenterTest.parse(aSmiles));
        //link types are compared without stereochemistry, since the stereo SMILES codes do not encode them
        Assertions.assertEquals(BRICSFragmenterTest.canonical(List.of(anExpectedFragments), false),
                BRICSFragmenterTest.toCanonicalSmiles(tmpFragments, false), aSmiles);
        if (isStereo) {
            Assertions.assertEquals(BRICSFragmenterTest.toStereoSmiles(List.of(anExpectedFragments), false),
                    BRICSFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false), aSmiles);
        }
    }
    //
    /**
     * Parses the given SMILES code with kekulization, as MORTAR does.
     *
     * @param aSmiles SMILES code
     * @return molecule
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static IAtomContainer parse(String aSmiles) throws CDKException {
        return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(aSmiles);
    }
    //
    /**
     * Returns the sorted canonical SMILES codes of the given fragments, see {@link #canonical(String, boolean)}.
     *
     * @param aFragmentList fragments
     * @param isStereo whether stereochemistry is encoded
     * @return sorted list of canonical SMILES codes
     * @throws CDKException if a SMILES code cannot be created
     */
    static List<String> toCanonicalSmiles(List<IAtomContainer> aFragmentList, boolean isStereo) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(BRICSFragmenterTest.toCanonicalSmiles(tmpFragment, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the canonical SMILES code of the given fragment, see {@link #canonical(String, boolean)}.
     *
     * @param aFragment fragment
     * @param isStereo whether stereochemistry is encoded
     * @return canonical SMILES code
     * @throws CDKException if a SMILES code cannot be created
     */
    static String toCanonicalSmiles(IAtomContainer aFragment, boolean isStereo) throws CDKException {
        if (isStereo) {
            return BRICSFragmenterTest.toStereoSmiles(aFragment, false);
        }
        String tmpSmiles = SmilesGenerator.create(aFragment, SmiFlavor.Unique | SmiFlavor.AtomicMass,
                new int[aFragment.getAtomCount()]);
        return BRICSFragmenterTest.canonical(tmpSmiles, false);
    }
    //
    /**
     * Returns the unique SMILES code with stereochemistry of the given fragment as MORTAR creates it, i.e. with all link
     * atoms written as '*' (the link types are not encoded). The CDK creates canonical SMILES codes with
     * stereochemistry based on InChI numbers, which do not support pseudo atoms with mass numbers.
     *
     * @param aFragment fragment, it is not altered
     * @param isNonStereogenicElementRemoved whether stereo elements of atoms that are not stereogenic according to the
     *                                       CDK are removed first
     * @return unique SMILES code with stereochemistry
     * @throws CDKException if the fragment cannot be cloned
     */
    @SuppressWarnings("rawtypes")
    static String toStereoSmiles(IAtomContainer aFragment, boolean isNonStereogenicElementRemoved) throws CDKException {
        IAtomContainer tmpFragment;
        try {
            tmpFragment = aFragment.clone();
        } catch (CloneNotSupportedException anException) {
            throw new CDKException(anException.toString(), anException);
        }
        for (IAtom tmpAtom : tmpFragment.atoms()) {
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAtom.setMassNumber(null);
            }
        }
        //the aromaticity is perceived on the fragment for a representation independent of the parent molecule
        SmartsPattern.prepare(tmpFragment);
        if (isNonStereogenicElementRemoved && tmpFragment.stereoElements().iterator().hasNext()) {
            Stereocenters tmpStereocenters = Stereocenters.of(tmpFragment);
            List<IStereoElement> tmpElements = new ArrayList<>();
            for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                if (tmpElement.getConfigClass() == IStereoElement.TH
                        && !tmpStereocenters.isStereocenter(tmpFragment.indexOf((IAtom) tmpElement.getFocus()))) {
                    continue;
                }
                if (tmpElement.getConfigClass() == IStereoElement.CT
                        && (!tmpStereocenters.isStereocenter(tmpFragment.indexOf(((IBond) tmpElement.getFocus()).getBegin()))
                        || !tmpStereocenters.isStereocenter(tmpFragment.indexOf(((IBond) tmpElement.getFocus()).getEnd())))) {
                    continue;
                }
                tmpElements.add(tmpElement);
            }
            tmpFragment.setStereoElements(tmpElements);
        }
        return ChemUtil.createUniqueSmiles(tmpFragment, true);
    }
    //
    /**
     * Returns the sorted unique SMILES codes with stereochemistry of the given SMILES codes, see
     * {@link #toStereoSmiles(IAtomContainer, boolean)}.
     *
     * @param aSmilesList SMILES codes
     * @param isNonStereogenicElementRemoved whether stereo elements of non-stereogenic atoms are removed first
     * @return sorted list of unique SMILES codes with stereochemistry
     * @throws CDKException if a SMILES code cannot be parsed
     */
    static List<String> toStereoSmiles(List<String> aSmilesList, boolean isNonStereogenicElementRemoved) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesList.size());
        for (String tmpSmiles : aSmilesList) {
            IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(tmpSmiles);
            SmartsPattern.prepare(tmpMolecule);
            tmpSmilesList.add(BRICSFragmenterTest.toStereoSmiles(tmpMolecule, isNonStereogenicElementRemoved));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the sorted unique SMILES codes with stereochemistry of the given fragments, see
     * {@link #toStereoSmiles(IAtomContainer, boolean)}.
     *
     * @param aFragmentList fragments
     * @param isNonStereogenicElementRemoved whether stereo elements of non-stereogenic atoms are removed first
     * @return sorted list of unique SMILES codes with stereochemistry
     * @throws CDKException if a fragment cannot be cloned
     */
    static List<String> toStereoSmilesOfFragments(List<IAtomContainer> aFragmentList, boolean isNonStereogenicElementRemoved)
            throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(BRICSFragmenterTest.toStereoSmiles(tmpFragment, isNonStereogenicElementRemoved));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the sorted canonical SMILES codes of the given SMILES codes, see {@link #canonical(String, boolean)}.
     *
     * @param aSmilesList SMILES codes
     * @param isStereo whether stereochemistry is encoded
     * @return sorted list of canonical SMILES codes
     * @throws CDKException if a SMILES code cannot be parsed
     */
    static List<String> canonical(List<String> aSmilesList, boolean isStereo) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesList.size());
        for (String tmpSmiles : aSmilesList) {
            tmpSmilesList.add(BRICSFragmenterTest.canonical(tmpSmiles, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns a canonical SMILES code of the given SMILES code that is independent of the aromaticity notation of the
     * input: the SMILES is parsed and kekulized, aromaticity is perceived with the Daylight model (as by the fragmenter), and a unique SMILES
     * with mass numbers (i.e. link types) and aromatic symbols is created.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether stereochemistry is encoded
     * @return canonical SMILES code
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static String canonical(String aSmiles, boolean isStereo) throws CDKException {
        IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(aSmiles);
        SmartsPattern.prepare(tmpMolecule);
        if (isStereo) {
            return BRICSFragmenterTest.toStereoSmiles(tmpMolecule, false);
        }
        return SmilesGenerator.create(tmpMolecule, SmiFlavor.Unique | SmiFlavor.AtomicMass | SmiFlavor.UseAromaticSymbols,
                new int[tmpMolecule.getAtomCount()]);
    }
}
