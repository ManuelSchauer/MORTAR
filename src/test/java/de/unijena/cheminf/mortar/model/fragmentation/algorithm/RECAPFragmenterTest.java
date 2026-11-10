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
import org.junit.jupiter.api.Assumptions;
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

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests the {@link RECAPFragmenter}. The expected fragments were derived by hand from the RECAP rules; the bond type
 * that cleaves a bond is given in the comments. Fragments are compared as canonical SMILES with the bond types encoded
 * as mass numbers of the link atoms, e.g. [5*]Oc1ccc(F)cc1.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class RECAPFragmenterTest {
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public RECAPFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the settings, the rule set, copying, and the restoration of the default settings.
     */
    @Test
    public void basicSettingsTest() {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        Assertions.assertEquals(RECAPFragmenter.ALGORITHM_NAME, tmpFragmenter.getFragmentationAlgorithmName());
        Assertions.assertFalse(tmpFragmenter.getFragmentationAlgorithmDisplayName().isBlank());
        Assertions.assertEquals(5, tmpFragmenter.settingsProperties().size());
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
        }
        Assertions.assertEquals(RECAPFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(RECAPFragmenter.SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT,
                tmpFragmenter.getSmallGroupMaximumCarbonAtomsSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertFalse(tmpFragmenter.getStrictFunctionalGroupsSetting());
        Assertions.assertEquals(RECAPFragmenter.FragmentOutputOption.TERMINAL_FRAGMENTS, tmpFragmenter.getFragmentOutputSetting());
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(2);
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        tmpFragmenter.setStrictFunctionalGroupsSetting(true);
        tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS);
        RECAPFragmenter tmpCopy = (RECAPFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION, tmpCopy.getSaturationSetting());
        Assertions.assertEquals(2, tmpCopy.getSmallGroupMaximumCarbonAtomsSetting());
        Assertions.assertFalse(tmpCopy.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertTrue(tmpCopy.getStrictFunctionalGroupsSetting());
        Assertions.assertEquals(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS, tmpCopy.getFragmentOutputSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(RECAPFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(RECAPFragmenter.SMALL_GROUP_MAXIMUM_CARBON_ATOMS_DEFAULT,
                tmpFragmenter.getSmallGroupMaximumCarbonAtomsSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertEquals(RECAPFragmenter.STRICT_FUNCTIONAL_GROUPS_OPTION_DEFAULT,
                tmpFragmenter.getStrictFunctionalGroupsSetting());
        Assertions.assertEquals(RECAPFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT, tmpFragmenter.getFragmentOutputSetting());
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setSaturationSetting(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setFragmentOutputSetting(null));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(-1));
        List<String> tmpBondTypeNames = RECAPFragmenter.getBondTypeNames();
        Assertions.assertEquals(11, tmpBondTypeNames.size());
        Assertions.assertEquals(11, new HashSet<>(tmpBondTypeNames).size());
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(null));
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
    }
    //
    /**
     * Tests the fragmentation of cisapride, the example of the authors' presentation, which shows the fragments
     * [1]C(=O)-aryl, [1]NH-piperidine-[3]N, [3]CH2-CH2-[5]CH2, and [5]O-fluorophenyl: the amide (1), the amine (3), and the
     * ether (5) at the propyl linker are cleaved, but not the methoxy and amino groups, and the oxygen atom of the aryl
     * ether stays at the aromatic ring.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void cisaprideTest() throws Exception {
        this.assertFragments("COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC",
                "[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]NC1CCN([3*])CC1OC", "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1");
        //the stereocentres of the piperidine ring are kept
        this.assertFragments(new RECAPFragmenter(), "COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC", true,
                "[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]N[C@@H]1CCN([3*])C[C@@H]1OC", "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1");
    }
    //
    /**
     * Tests examples for all eleven bond types.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bondTypesTest() throws Exception {
        //1 amide
        this.assertFragments("CCCCCC(=O)NCCCCC", "[1*]C(=O)CCCCC", "[1*]NCCCCC");
        //2 ester, the oxygen atom stays at the alcohol part
        this.assertFragments("CCCCCC(=O)OCCCCC", "[2*]C(=O)CCCCC", "[2*]OCCCCC");
        //3 amine, all three bonds of a tertiary amine
        this.assertFragments("CCCCCN(CCCCC)CCCCC", "[3*]N([3*])[3*]", "[3*]CCCCC", "[3*]CCCCC", "[3*]CCCCC");
        //4 urea, both bonds
        this.assertFragments("CCCCCNC(=O)NCCCCC", "[4*]C([4*])=O", "[4*]NCCCCC", "[4*]NCCCCC");
        //5 ether, both bonds of a dialkyl ether
        this.assertFragments("CCCCCOCCCCC", "[5*]O[5*]", "[5*]CCCCC", "[5*]CCCCC");
        //5 ether, aryl alkyl ether only at the sp3 carbon atom
        this.assertFragments("CCCCCOc1ccccc1", "[5*]CCCCC", "[5*]Oc1ccccc1");
        //6 olefin
        this.assertFragments("CCCCCC=CCCCCC", "[6*]=CCCCCC", "[6*]=CCCCCC");
        //7 quaternary nitrogen, all four bonds
        this.assertFragments("CCCCC[N+](CCCCC)(CCCCC)CCCCC", "[7*][N+]([7*])([7*])[7*]", "[7*]CCCCC", "[7*]CCCCC",
                "[7*]CCCCC", "[7*]CCCCC");
        //8 aromatic nitrogen - aliphatic carbon
        this.assertFragments("CCCCCn1ccnc1", "[8*]CCCCC", "[8*]n1ccnc1");
        //9 lactam nitrogen - aliphatic carbon
        this.assertFragments("CCCCCN1CCCC1=O", "[9*]CCCCC", "[9*]N1CCCC1=O");
        //10 aromatic carbon - aromatic carbon
        this.assertFragments("c1ccccc1-c1ccccn1", "[10*]c1ccccc1", "[10*]c1ccccn1");
        //11 sulphonamide
        this.assertFragments("CCCCCNS(=O)(=O)c1ccccc1", "[11*]NCCCCC", "[11*]S(=O)(=O)c1ccccc1");
        //carbamate: the amide bond (1) is cleaved, the tert-butoxy group is a small terminal group (ester, 2)
        this.assertFragments("CC(C)(C)OC(=O)NCCCCC", "[1*]C(=O)OC(C)(C)C", "[1*]NCCCCC");
    }
    //
    /**
     * Tests molecules whose fragments are shown in the original publication (Figures 9 and 10): isopropylamino
     * (propranolol), dimethylamino (the most frequent fragment) and benzhydryl (diphenhydramine), the piperidine fragment
     * of pethidine analogues with the ethyl ester, 1-(2-methoxyphenyl)piperazine, tert-butylamino (salbutamol), and the
     * type 8 bond of uracil.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void publicationFragmentsTest() throws Exception {
        this.assertFragments("CC(C)NCC(O)COc1cccc2ccccc12", "[3*]NC(C)C", "[3*]CC(O)C[5*]", "[5*]Oc1cccc2ccccc12");
        this.assertFragments("CN(C)CCOC(c1ccccc1)c1ccccc1", "[3*]N(C)C", "[3*]CC[5*]", "[5*]O[5*]",
                "[5*]C(c1ccccc1)c1ccccc1");
        this.assertFragments("CCOC(=O)C1(c2ccccc2)CCN(CCc2ccccc2)CC1", "[3*]N1CCC(C(=O)OCC)(c2ccccc2)CC1",
                "[3*]CCc1ccccc1");
        this.assertFragments("COc1ccccc1N1CCN(CCCCc2ccccc2)CC1", "[3*]N1CCN(c2ccccc2OC)CC1", "[3*]CCCCc1ccccc1");
        this.assertFragments("CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "[3*]NC(C)(C)C", "[3*]CC(O)c1ccc(O)c(CO)c1");
        this.assertFragments("O=c1ccn(CCCCC)c(=O)[nH]1", "[8*]n1ccc(=O)[nH]c1=O", "[8*]CCCCC");
    }
    //
    /**
     * Tests molecules whose fragments are shown in Figures 8, 9 and 12 of the original publication:
     * <ul>
     *     <li>Figure 8: oxotremorine, a 2-connection fragment (the butynyl linker between the amine (3) and the lactam
     *     nitrogen atom (9)), and a monobactam with a carbamate and an N-sulfonic acid, a 3-connection fragment
     *     (amide 1, ester 2 and sulphonamide 11 at the azetidinone);</li>
     *     <li>Figure 9: chlorprothixene (exocyclic double bond, 6), retinol (vitamin A, polyene cleaved into
     *     two-carbon fragments), the quaternary ring nitrogen atom of a piperidinium ion (7), imipramine and verapamil
     *     (dimethylamino and methylamino fragments, 3), nicardipine (the methyl ester is kept, the aminoethyl ester is
     *     cleaved), adenosine (aromatic nitrogen atom, 8), cefotaxime (the acetate is kept), arbutin (glycoside, 5),
     *     and prostaglandin F2-alpha (olefins, 6);</li>
     *     <li>Figure 12: N-benzoyl tryptophan derivatives, whose N-methyl amide, ethyl ester and N-formyl group are not
     *     cleaved off (small terminal groups), while the benzamide (1) is cleaved.</li>
     * </ul>
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void publicationFiguresTest() throws Exception {
        //Figure 8
        this.assertFragments("O=C1CCCN1CC#CCN1CCCC1", "[3*]CC#CC[9*]", "[3*]N1CCCC1", "[9*]N1CCCC1=O");
        this.assertFragments("NC(=O)OCC1C(NC(=O)c2ccccc2)C(=O)N1S(=O)(=O)O", "[1*]C(=O)c1ccccc1",
                "[1*]NC1C(=O)N([11*])C1CO[2*]", "[11*]S(=O)(=O)O", "[2*]C(N)=O");
        //Figure 9
        this.assertFragments("CN(C)CCC=C1c2ccccc2Sc2ccc(Cl)cc12", "[3*]CCC=[6*]", "[3*]N(C)C",
                "[6*]=C1c2ccccc2Sc2ccc(Cl)cc21");
        this.assertFragments("CC1=C(C=CC(C)=CC=CC(C)=CCO)C(C)(C)CCC1", "[6*]=CC1=C(C)CCCC1(C)C", "[6*]=CC(=[6*])C",
                "[6*]=CC(=[6*])C", "[6*]=CC=[6*]", "[6*]=CCO");
        this.assertFragments("C[N+]1(CCCCC1)CC(O)c1ccccc1", "[7*]CC(O)c1ccccc1", "[7*][N+]1(C)CCCCC1");
        this.assertFragments("CN(C)CCCN1c2ccccc2CCc2ccccc21", "[3*]CCC[3*]", "[3*]N(C)C",
                "[3*]N1c2ccccc2CCc2ccccc21");
        this.assertFragments("COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
                "[3*]CCCC(C#N)(c1ccc(OC)c(OC)c1)C(C)C", "[3*]CCc1ccc(OC)c(OC)c1", "[3*]N([3*])C");
        this.assertFragments("COC(=O)C1=C(C)NC(C)=C(C(=O)OCCN(C)Cc2ccccc2)C1c1cccc([N+](=O)[O-])c1",
                "[2*]C(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1cccc([N+](=O)[O-])c1", "[2*]OCC[3*]", "[3*]Cc1ccccc1",
                "[3*]N([3*])C");
        this.assertFragments("Nc1ncnc2n(cnc12)C1OC(CO)C(O)C1O", "[8*]C1OC(CO)C(O)C1O", "[8*]n1cnc2c(N)ncnc21");
        this.assertFragments("CC(=O)OCC1=C(N2C(SC1)C(NC(=O)C(=NOC)c1csc(N)n1)C2=O)C(=O)O",
                "[1*]C(=O)C(=NOC)c1csc(N)n1", "[1*]NC1C(=O)N2C(C(=O)O)=C(COC(C)=O)CSC12");
        this.assertFragments("OCC1OC(Oc2ccc(O)cc2)C(O)C(O)C1O", "[5*]C1OC(CO)C(O)C(O)C1O", "[5*]Oc1ccc(O)cc1");
        this.assertFragments("CCCCCC(O)C=CC1C(O)CC(O)C1CC=CCCCC(=O)O", "[6*]=CC(O)CCCCC",
                "[6*]=CCC1C(O)CC(O)C1C=[6*]", "[6*]=CCCCC(=O)O");
        //Figure 12
        this.assertFragments("O=C(c1ccccc1)NC(Cc1c[nH]c2ccccc12)C(=O)NC", "[1*]C(=O)c1ccccc1",
                "[1*]NC(Cc1c[nH]c2ccccc12)C(=O)NC");
        this.assertFragments("O=C(c1ccccc1)NC(Cc1c[nH]c2ccccc12)C(=O)OCC", "[1*]C(=O)c1ccccc1",
                "[1*]NC(Cc1c[nH]c2ccccc12)C(=O)OCC");
        this.assertFragments("O=C(c1ccccc1)NC(Cc1cn(C=O)c2ccccc12)C(=O)N", "[1*]C(=O)c1ccccc1",
                "[1*]NC(Cc1cn(C=O)c2ccccc12)C(N)=O");
    }
    //
    /**
     * Tests the small terminal groups at the nitrogen atoms of amide, urea and sulphonamide bonds: alkyl groups with up
     * to the maximum number of carbon atoms in total are kept, as the N-methyl amide in Figure 12 of the publication,
     * whereas alkylated amino groups are cleaved off at amine bonds (Figures 9 and 10). The hierarchy applies the same
     * rule.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void acylatedNitrogenSmallGroupTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        //N-methyl amide (1), Figure 12
        String tmpMethylAmide = "O=C(c1ccccc1)NC(Cc1c[nH]c2ccccc12)C(=O)NC";
        this.assertFragments(tmpFragmenter, tmpMethylAmide, false, "[1*]C(=O)c1ccccc1",
                "[1*]NC(Cc1c[nH]c2ccccc12)C(=O)NC");
        //N,N-diethyl amide (DEET), not cleaved at all
        this.assertFragments(tmpFragmenter, "CCN(CC)C(=O)c1cccc(C)c1", false, "CCN(CC)C(=O)c1cccc(C)c1");
        //N-methyl sulphonamide (11) of sumatriptan, whose dimethylamino group is cleaved off (3)
        this.assertFragments(tmpFragmenter, "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1", false,
                "[3*]CCc1c[nH]c2ccc(CS(=O)(=O)NC)cc12", "[3*]N(C)C");
        //N-methyl urea (4): only the bond to the pentylamino group is cleaved
        this.assertFragments(tmpFragmenter, "CCCCCNC(=O)NC", false, "[4*]C(=O)NC", "[4*]NCCCCC");
        //dipropylamino has six carbon atoms and is cleaved off
        this.assertFragments(tmpFragmenter, "CCCCCC(=O)N(CCC)CCC", false, "[1*]C(=O)CCCCC", "[1*]N(CCC)CCC");
        //methylamino at an amine bond (3) is cleaved off
        this.assertFragments(tmpFragmenter, "CNC(C)C(O)c1ccccc1", false, "[3*]C(C)C(O)c1ccccc1", "[3*]NC");
        //the leaves of the hierarchy are the single-pass fragments
        this.assertLeavesAreSinglePassFragments(tmpFragmenter, RECAPFragmenterTest.parse(tmpMethylAmide));
        Assertions.assertTrue(
                tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse("CCN(CC)C(=O)c1cccc(C)c1")).isLeaf());
        //the setting applies to the nitrogen atoms as well, the amino group of a primary amide is always kept
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(0);
        this.assertFragments(tmpFragmenter, "CCN(CC)C(=O)c1cccc(C)c1", false, "[1*]C(=O)c1cccc(C)c1", "[1*]N(CC)CC");
        this.assertFragments(tmpFragmenter, "NC(=O)c1ccccc1", false, "NC(=O)c1ccccc1");
    }
    //
    /**
     * Tests that explicit hydrogen atoms that are not suppressed (isotopes) are ignored in the check for small terminal
     * groups, so deuterated groups are treated like their protium forms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hydrogenIsotopeTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        //deuterated methoxy and amino groups are small terminal groups
        this.assertFragments(tmpFragmenter, "[2H]C([2H])([2H])Oc1ccccc1", false, "[2H]C([2H])([2H])Oc1ccccc1");
        this.assertFragments(tmpFragmenter, "[2H]N([2H])CC1CCCCC1", false, "[2H]N([2H])CC1CCCCC1");
        //the deuterated methyl groups are kept at the amine nitrogen atom, which is cleaved off like dimethylamino
        String tmpDeuteratedAmine = "[2H]C([2H])([2H])N(CCc1ccccc1)C([2H])([2H])[2H]";
        this.assertFragments(tmpFragmenter, tmpDeuteratedAmine, false, "[3*]CCc1ccccc1",
                "[3*]N(C([2H])([2H])[2H])C([2H])([2H])[2H]");
        this.assertFragments(tmpFragmenter, "CN(CCc1ccccc1)C", false, "[3*]CCc1ccccc1", "[3*]N(C)C");
        //the hierarchy applies the same rule
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(tmpDeuteratedAmine);
        Assertions.assertEquals(2, tmpFragmenter.buildHierarchy(tmpMolecule).getAllDescendants().size());
        this.assertLeavesAreSinglePassFragments(tmpFragmenter, tmpMolecule);
        Assertions.assertTrue(
                tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse("[2H]C([2H])([2H])Oc1ccccc1")).isLeaf());
    }
    //
    /**
     * Tests molecules with bonds of unset order, e.g. aromatic bonds of molecules parsed without kekulization: a Kekulé
     * structure is assigned, so that their fragments equal those of the kekulized molecules, and molecules that cannot
     * be kekulized are not fragmented.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void unsetBondOrderTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpParser.kekulise(false);
        for (String tmpSmiles : new String[] {"CCCCCn1ccnc1", "c1ccccc1-c1ccccc1", "CCCCCOc1ccccc1",
                "COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC"}) {
            IAtomContainer tmpMolecule = tmpParser.parseSmiles(tmpSmiles);
            Assertions.assertTrue(this.hasUnsetBondOrder(tmpMolecule), tmpSmiles);
            List<IAtomContainer> tmpExpectedFragments = tmpFragmenter.fragmentMolecule(
                    RECAPFragmenterTest.parse(tmpSmiles));
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            Assertions.assertTrue(tmpExpectedFragments.size() > 1, tmpSmiles);
            Assertions.assertEquals(RECAPFragmenterTest.toCanonicalSmiles(tmpExpectedFragments, false),
                    RECAPFragmenterTest.toCanonicalSmiles(tmpFragments, false), tmpSmiles);
            Assertions.assertEquals(RECAPFragmenterTest.toStereoSmilesOfFragments(tmpExpectedFragments, false),
                    RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false), tmpSmiles);
            this.assertLeavesAreSinglePassFragments(tmpFragmenter, tmpMolecule);
            //the input molecule is not altered
            Assertions.assertTrue(this.hasUnsetBondOrder(tmpMolecule), tmpSmiles);
        }
        this.assertFragments(tmpFragmenter, "CCCCCn1ccnc1", false, "[8*]CCCCC", "[8*]n1ccnc1");
        //bonds without order (null) are treated like bonds with unset order
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCCCCn1ccnc1");
        for (IBond tmpBond : tmpMolecule.bonds()) {
            if (tmpBond.getOrder() == IBond.Order.UNSET) {
                tmpBond.setOrder(null);
            }
        }
        Assertions.assertEquals(RECAPFragmenterTest.canonical(List.of("[8*]CCCCC", "[8*]n1ccnc1"), false),
                RECAPFragmenterTest.toCanonicalSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule), false));
        //no Kekulé structure can be assigned
        IAtomContainer tmpNonKekulizable = tmpParser.parseSmiles("c1cccc1");
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.fragmentMolecule(tmpNonKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.buildHierarchy(tmpNonKekulizable));
        //an unset order of a bond that is not aromatic cannot be determined
        IAtomContainer tmpUnsetEther = RECAPFragmenterTest.parse("CCCCCOCCCCC");
        tmpUnsetEther.getBond(0).setOrder(IBond.Order.UNSET);
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.fragmentMolecule(tmpUnsetEther));
    }
    //
    /**
     * Tests that ring bonds and the bonds of small terminal groups are not cleaved, and that no bond outside the
     * eleven types is cleaved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void notCleavedTest() throws Exception {
        String[] tmpUnfragmentedMolecules = {
                //small terminal groups: methyl, ethyl, butyl, methoxy, ethoxy, ethyl and methyl esters, acetate,
                //acetyl and butyryl amides, mesyl, amino
                "CCOc1ccccc1", "COc1ccccc1", "CCCCN(CCCC)CCCC", "CCOC(=O)c1ccccc1", "COC(=O)C1CCCCC1",
                "CC(=O)OCC1CCCCC1", "CC(=O)NCCCCC", "CCCC(=O)NCCCCC", "CS(=O)(=O)NCCCCC", "NCc1ccccc1",
                "NS(=O)(=O)c1ccccc1", "NC(=O)c1ccccc1",
                //no RECAP bond types: alkyl-aryl bonds, alcohols, carboxylic acids
                "Cc1ccccc1", "OCC1CCCCC1", "OC(=O)c1ccccc1",
                //no RECAP bond types: aryl ethers and amines, ketones, thioethers, alkyl-aryl bonds, ring bonds
                "c1ccccc1Oc1ccccc1", "c1ccccc1N1CCCCC1", "O=C(c1ccccc1)c1ccccc1", "CCCCCSCCCCC", "CCCCCc1ccccc1",
                "O=C1CCCCN1", "C1CCOC1", "c1ccc2ccccc2c1"};
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        for (String tmpSmiles : tmpUnfragmentedMolecules) {
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpSmiles));
            Assertions.assertEquals(1, tmpFragments.size(), tmpSmiles);
            Assertions.assertEquals(RECAPFragmenterTest.canonical(tmpSmiles, false),
                    RECAPFragmenterTest.toCanonicalSmiles(tmpFragments.getFirst(), false), tmpSmiles);
        }
    }
    //
    /**
     * Tests the small terminal groups and the setting for their maximum size.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void smallTerminalGroupTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        //the propyl groups are small, but the pentyl group can be cleaved off the dipropylamino part, since alkylated
        //amino groups are not small
        this.assertFragments(tmpFragmenter, "CCCCCN(CCC)CCC", false, "[3*]CCCCC", "[3*]N(CCC)CCC");
        //trimethylammonium is cleaved off, the methyl groups are not
        this.assertFragments(tmpFragmenter, "CCCCCC[N+](C)(C)C", false, "[7*]CCCCCC", "[7*][N+](C)(C)C");
        //rings are no small groups
        this.assertFragments(tmpFragmenter, "C1CC1N(C1CC1)C1CC1", false, "[3*]N([3*])[3*]", "[3*]C1CC1", "[3*]C1CC1",
                "[3*]C1CC1");
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(3);
        this.assertFragments(tmpFragmenter, "CCCCN(CCCC)CCCC", false, "[3*]N([3*])[3*]", "[3*]CCCC", "[3*]CCCC", "[3*]CCCC");
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(0);
        this.assertFragments(tmpFragmenter, "COc1ccccc1", false, "[5*]C", "[5*]Oc1ccccc1");
        //the amino group is protected independent of the setting
        this.assertFragments(tmpFragmenter, "NCC1CCCCC1", false, "NCC1CCCCC1");
        //a linking fragment is allowed independent of its size, here the single oxygen atom
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(7);
        this.assertFragments(tmpFragmenter, "CCCCCCCCOCCCCCCCC", false, "[5*]O[5*]", "[5*]CCCCCCCC", "[5*]CCCCCCCC");
        tmpFragmenter.setSmallGroupMaximumCarbonAtomsSetting(8);
        this.assertFragments(tmpFragmenter, "CCCCCCCCOCCCCCCCC", false, "CCCCCCCCOCCCCCCCC");
    }
    //
    /**
     * Tests the saturation of the open valences with hydrogen atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hydrogenSaturationTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCCCOCCCCC", false, "O", "CCCCC", "CCCCC");
        //cleaved double bonds are saturated with two hydrogen atoms
        this.assertFragments(tmpFragmenter, "CCCCCC=CCCCCC", false, "CCCCCC", "CCCCCC");
        this.assertFragments(tmpFragmenter, "COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC", false,
                "O=Cc1cc(Cl)c(N)cc1OC", "NC1CCNCC1OC", "CCC", "Oc1ccc(F)cc1");
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("CCCCCOCCCCC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Assertions.assertFalse(tmpAtom instanceof IPseudoAtom);
                Assertions.assertNotEquals(1, tmpAtom.getAtomicNumber());
            }
        }
    }
    //
    /**
     * Tests the link atoms: pseudo atoms with the bond type as mass number and atom property.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void linkAtomTest() throws Exception {
        Set<Integer> tmpBondTypes = new HashSet<>();
        for (IAtomContainer tmpFragment : new RECAPFragmenter().fragmentMolecule(
                RECAPFragmenterTest.parse("COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom) {
                    Integer tmpBondType = tmpAtom.getProperty(RECAPFragmenter.BOND_TYPE_PROPERTY_KEY);
                    Assertions.assertEquals(tmpBondType, tmpAtom.getMassNumber());
                    Assertions.assertEquals(0, tmpAtom.getImplicitHydrogenCount());
                    Assertions.assertEquals(1, tmpFragment.getConnectedBondsCount(tmpAtom));
                    tmpBondTypes.add(tmpBondType);
                }
            }
        }
        Assertions.assertEquals(Set.of(1, 3, 5), tmpBondTypes);
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
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        Assertions.assertEquals(1, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("Cc1ccccc1")).size());
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("Cc1ccccc1")).isEmpty());
        Assertions.assertEquals(3, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("CCCCCOCCCCC")).size());
        //counter ions are returned separately, also for molecules without cleaved bonds
        this.assertFragments("CCCCCOCCCCC.[Na+].[Cl-]", "[5*]O[5*]", "[5*]CCCCC", "[5*]CCCCC", "[Na+]", "[Cl-]");
        this.assertFragments("Cc1ccccc1.[Na+]", "Cc1ccccc1", "[Na+]");
        //explicit hydrogen atoms do not change the result
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse("COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC");
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        Assertions.assertEquals(
                RECAPFragmenterTest.canonical(List.of("[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]NC1CCN([3*])CC1OC",
                        "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1"), false),
                RECAPFragmenterTest.toCanonicalSmiles(new RECAPFragmenter().fragmentMolecule(tmpMolecule), false));
        //pseudo atoms of a preceding fragmentation count as attachment points, so *CCC is a linking fragment
        this.assertFragments("*CCCOCCCCC", "[5*]O[5*]", "[5*]CCC*", "[5*]CCCCC");
        //the input molecule is not altered
        IAtomContainer tmpInput = RECAPFragmenterTest.parse("CCCCCOCCCCC");
        new RECAPFragmenter().fragmentMolecule(tmpInput);
        Assertions.assertEquals(11, tmpInput.getAtomCount());
        Assertions.assertEquals(10, tmpInput.getBondCount());
    }
    //
    /**
     * Tests the transfer of tetrahedral and double bond stereo elements to the fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void stereoTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        //the ether oxygen atom of the stereocentre is replaced by the link atom at the same position
        this.assertFragments(tmpFragmenter, "CCCCC[C@@H](C)OCCCCC", true, "CCCCC[C@@H](C)[5*]", "[5*]O[5*]", "[5*]CCCCC");
        this.assertFragments(tmpFragmenter, "CCCCC[C@H](C)OCCCCC", true, "CCCCC[C@H](C)[5*]", "[5*]O[5*]", "[5*]CCCCC");
        //two link atoms at a former stereocentre, the stereo element is removed
        this.assertFragments(tmpFragmenter, "CCCCCN[C@@H](CCCCC)OCCCCC", true, "[3*]C([5*])CCCCC", "[5*]O[5*]",
                "[3*]N[3*]", "[5*]CCCCC", "[3*]CCCCC");
        //the lactam nitrogen atom is a ligand of the stereo double bond (9); the double bond itself is not cleaved (6),
        //because the methyl group is too small
        this.assertFragments(tmpFragmenter, "O=C1CCCN1/C(CCCCC)=C/C", true, "[9*]/C(CCCCC)=C/C", "[9*]N1CCCC1=O");
        //cleaved stereo double bond
        this.assertFragments(tmpFragmenter, "CCCCC/C=C/CCCCC", true, "[6*]=CCCCCC", "[6*]=CCCCCC");
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        //with hydrogen saturation, the former stereocentre has two hydrogen atoms
        this.assertFragments(tmpFragmenter, "CCCCC[C@@H](C)OCCCCC", true, "CCCCCCC", "O", "CCCCC");
        //with hydrogen saturation, the configuration refers to the remaining substituent
        this.assertFragments(tmpFragmenter, "O=C1CCCN1/C(CCCCC)=C/C", true, "CCCCC/C=C\\C", "O=C1CCCN1");
    }
    //
    /**
     * Tests that cyanates and cyanamides are not cleaved as ethers or amines: the nitrile carbon atom is excluded like
     * an imine or carbonyl carbon atom.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void nitrileTest() throws Exception {
        this.assertFragments("CCCCCCCCOC#N", "CCCCCCCCOC#N");
        this.assertFragments("N#CN(CCCCCC)CCCCCC", "N#CN(CCCCCC)CCCCCC");
    }
    //
    /**
     * Tests the strict functional group setting: bonds of extended functional groups are cleaved with the default
     * rules, but not with the strict ones, while plain amides and esters are cleaved with both.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void strictFunctionalGroupsTest() throws Exception {
        RECAPFragmenter tmpStrictFragmenter = new RECAPFragmenter();
        tmpStrictFragmenter.setStrictFunctionalGroupsSetting(true);
        String[][] tmpCases = {
                //smiles, default fragments separated by '.'
                {"CCCCCNC(=O)OCc1ccccc1", "[1*]C([2*])=O.[1*]NCCCCC.[2*]OCc1ccccc1"}, //carbamate (1, 2)
                {"O=C(OC(=O)c1ccccc1)c1ccccc1", "[2*]C(=O)c1ccccc1.[2*]C(=O)c1ccccc1.[2*]O[2*]"}, //anhydride (2)
                {"O=C(NC(=O)c1ccccc1)c1ccccc1", "[1*]C(=O)c1ccccc1.[1*]C(=O)c1ccccc1.[1*]N[1*]"}, //imide (1)
                {"O=C(NC(=O)NCc1ccccc1)c1ccccc1", "[1*]C(=O)c1ccccc1.[1*]N[4*].[4*]C([4*])=O.[4*]NCc1ccccc1"}, //acyl urea
                {"O=S(=O)(NCc1ccccc1)NCc1ccccc1", "[11*]NCc1ccccc1.[11*]NCc1ccccc1.[11*]S([11*])(=O)=O"}, //sulphamide
                {"O=C1N(CCc2ccccc2)CCN1C", "[9*]CCc1ccccc1.[9*]N1C(=O)N(C)CC1"}, //cyclic urea (9)
                //sulphonyl ureas (4, 11): acetohexamide, and tolbutamide, whose butylamino group is a small group
                {"CC(=O)c1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
                        "[11*]S(=O)(=O)c1ccc(C(C)=O)cc1.[4*]C([4*])=O.[4*]NC1CCCCC1.[4*]N[11*]"},
                {"CCCCNC(=O)NS(=O)(=O)c1ccc(C)cc1", "[11*]S(=O)(=O)c1ccc(C)cc1.[4*]C(=O)NCCCC.[4*]N[11*]"},
        };
        for (String[] tmpCase : tmpCases) {
            this.assertFragments(tmpCase[0], tmpCase[1].split("\\."));
            this.assertFragments(tmpStrictFragmenter, tmpCase[0], false, tmpCase[0]);
        }
        //piperine: the polyene is not cleaved with the strict rules, but the amide is
        String tmpPiperine = "O=C(/C=C/C=C/c1ccc2OCOc2c1)N1CCCCC1";
        this.assertFragments(tmpPiperine, "[1*]C(=O)C=[6*]", "[1*]N1CCCCC1", "[6*]=CC=[6*]", "[6*]=Cc1ccc2OCOc2c1");
        this.assertFragments(tmpStrictFragmenter, tmpPiperine, false, "[1*]C(=O)C=CC=Cc1ccc2OCOc2c1", "[1*]N1CCCCC1");
        //enol ether: the olefin with an oxygen substituent is not cleaved, the ether is
        String tmpEnolEther = "CCCCCC/C=C/OCC1CCCCC1";
        this.assertFragments(tmpEnolEther, "[5*]CC1CCCCC1", "[5*]OC=[6*]", "[6*]=CCCCCCC");
        this.assertFragments(tmpStrictFragmenter, tmpEnolEther, false, "[5*]CC1CCCCC1", "[5*]OC=CCCCCCC");
        //plain amides and esters
        this.assertFragments(tmpStrictFragmenter, "O=C(Cc1ccccc1)NCc1ccccc1", false, "[1*]C(=O)Cc1ccccc1", "[1*]NCc1ccccc1");
        this.assertFragments(tmpStrictFragmenter, "O=C(OCc1ccccc1)c1ccccc1", false, "[2*]C(=O)c1ccccc1", "[2*]OCc1ccccc1");
    }
    //
    /**
     * Tests the RECAP hierarchy of 1,2-diphenoxyethane: the two ether bonds can be cleaved in either order, so the
     * hierarchy has six nodes, and the phenoxy fragments and the central fragment have two parents each. Its leaves are
     * the single-pass fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hierarchyTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse("c1ccccc1OCCOc1ccccc1");
        String tmpInputSmiles = RECAPFragmenterTest.toCanonicalSmiles(tmpMolecule, true);
        RECAPFragmenter.RECAPHierarchyNode tmpRoot = tmpFragmenter.buildHierarchy(tmpMolecule);
        Assertions.assertEquals(tmpInputSmiles, RECAPFragmenterTest.toCanonicalSmiles(tmpMolecule, true));
        Assertions.assertTrue(tmpRoot.isRoot());
        Assertions.assertFalse(tmpRoot.isLeaf());
        Assertions.assertEquals(4, tmpRoot.getChildren().size());
        Assertions.assertEquals(5, tmpRoot.getAllDescendants().size());
        Assertions.assertEquals(RECAPFragmenterTest.canonical(List.of("[5*]Oc1ccccc1", "[5*]CCOc1ccccc1",
                        "[5*]Oc1ccccc1", "[5*]CCOc1ccccc1"), false),
                RECAPFragmenterTest.toCanonicalSmiles(RECAPFragmenterTest.fragmentsOf(tmpRoot.getChildren()), false));
        List<RECAPFragmenter.RECAPHierarchyNode> tmpLeaves = tmpRoot.getLeaves();
        Assertions.assertEquals(RECAPFragmenterTest.toCanonicalSmiles(tmpFragmenter.fragmentMolecule(tmpMolecule), false),
                RECAPFragmenterTest.toCanonicalSmiles(RECAPFragmenterTest.fragmentsOf(tmpLeaves), false));
        //the phenoxy fragments result from the molecule and from the larger fragment, the central fragment from both
        //larger fragments
        List<IAtomContainer> tmpFragmentsWithTwoParents = new ArrayList<>();
        for (RECAPFragmenter.RECAPHierarchyNode tmpNode : tmpRoot.getAllDescendants()) {
            Assertions.assertTrue(tmpNode.getParents().size() <= 2);
            if (tmpNode.getParents().size() == 2) {
                tmpFragmentsWithTwoParents.add(tmpNode.getFragment());
            }
        }
        Assertions.assertEquals(RECAPFragmenterTest.canonical(List.of("[5*]Oc1ccccc1", "[5*]Oc1ccccc1", "[5*]CC[5*]"), false),
                RECAPFragmenterTest.toCanonicalSmiles(tmpFragmentsWithTwoParents, false));
        //all hierarchy fragments, i.e. all nodes except the root
        tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS);
        this.assertFragments(tmpFragmenter, "c1ccccc1OCCOc1ccccc1", false, "[5*]Oc1ccccc1", "[5*]Oc1ccccc1",
                "[5*]CCOc1ccccc1", "[5*]CCOc1ccccc1", "[5*]CC[5*]");
        //counter ions are children of the root and returned as fragments
        tmpRoot = tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse("c1ccccc1OCCOc1ccccc1.[Na+].[Cl-]"));
        Assertions.assertEquals(6, tmpRoot.getChildren().size());
        this.assertFragments(tmpFragmenter, "c1ccccc1OCCOc1ccccc1.[Na+].[Cl-]", false, "[5*]Oc1ccccc1", "[5*]Oc1ccccc1",
                "[5*]CCOc1ccccc1", "[5*]CCOc1ccccc1", "[5*]CC[5*]", "[Na+]", "[Cl-]");
        //molecules without cleavable bonds
        tmpRoot = tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse("c1ccccc1"));
        Assertions.assertTrue(tmpRoot.isLeaf());
        Assertions.assertTrue(tmpRoot.getLeaves().isEmpty());
        this.assertFragments(tmpFragmenter, "c1ccccc1", false, "c1ccccc1");
    }
    //
    /**
     * Tests the size of RECAP hierarchies: a chain of 13 para-phenylene rings (12 biaryl bonds) has one node per
     * contiguous section (91 nodes), 1,3,5-triphenylbenzene has 2^3 + 3 = 11 nodes, and a molecule with 14 benzamide
     * groups at a common core exceeds the maximum number of nodes. The leaves of the hierarchy of cisapride are its
     * single-pass fragments, including their stereochemistry.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hierarchySizeTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        String tmpChain = "c1ccc(cc1)-".repeat(12) + "c1ccccc1";
        RECAPFragmenter.RECAPHierarchyNode tmpRoot = tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse(tmpChain));
        Assertions.assertEquals(91, 1 + tmpRoot.getAllDescendants().size());
        Assertions.assertEquals(13, tmpRoot.getLeaves().size());
        tmpRoot = tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse("c1ccc(cc1)-c1cc(cc(c1)-c1ccccc1)-c1ccccc1"));
        Assertions.assertEquals(11, 1 + tmpRoot.getAllDescendants().size());
        String tmpStar = "C" + "C(NC(=O)c1ccccc1)".repeat(14) + "C";
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse(tmpStar)));
        tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS);
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpStar)));
        tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.TERMINAL_FRAGMENTS);
        Assertions.assertEquals(15, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpStar)).size());
        IAtomContainer tmpCisapride = RECAPFragmenterTest.parse(
                "COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC");
        tmpRoot = tmpFragmenter.buildHierarchy(tmpCisapride);
        Assertions.assertEquals(
                RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragmenter.fragmentMolecule(tmpCisapride), false),
                RECAPFragmenterTest.toStereoSmilesOfFragments(RECAPFragmenterTest.fragmentsOf(tmpRoot.getLeaves()), false));
    }
    //
    /**
     * Tests the limit of the total number of atoms of a RECAP hierarchy: a molecule with 13 benzamide groups at a
     * common core has 2^13 + 13 = 8205 nodes, fewer than the maximum number of nodes, but more than
     * {@value RECAPFragmenter#MAXIMUM_HIERARCHY_ATOM_COUNT} atoms in all fragments together (about 655,000), so it is
     * rejected with and without stereochemistry, while its single-pass fragments are returned. The hierarchy of a
     * molecule with ten benzamide groups and ten stereocentres (1034 nodes) is created, and its leaves are the
     * single-pass fragments, including their stereochemistry.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hierarchyAtomCountTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        for (String tmpStar : new String[] {"C" + "[C@H](NC(=O)c1ccccc1)".repeat(13) + "C",
                "C" + "C(NC(=O)c1ccccc1)".repeat(13) + "C"}) {
            tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.TERMINAL_FRAGMENTS);
            Assertions.assertThrows(IllegalArgumentException.class,
                    () -> tmpFragmenter.buildHierarchy(RECAPFragmenterTest.parse(tmpStar)), tmpStar);
            Assertions.assertEquals(14, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpStar)).size());
            tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS);
            Assertions.assertThrows(IllegalArgumentException.class,
                    () -> tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpStar)), tmpStar);
        }
        tmpFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.TERMINAL_FRAGMENTS);
        IAtomContainer tmpStar = RECAPFragmenterTest.parse("C" + "[C@H](NC(=O)c1ccccc1)".repeat(10) + "C");
        RECAPFragmenter.RECAPHierarchyNode tmpRoot = tmpFragmenter.buildHierarchy(tmpStar);
        Assertions.assertEquals(1034, 1 + tmpRoot.getAllDescendants().size());
        Assertions.assertEquals(
                RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragmenter.fragmentMolecule(tmpStar), false),
                RECAPFragmenterTest.toStereoSmilesOfFragments(
                        RECAPFragmenterTest.fragmentsOf(tmpRoot.getLeaves()), false));
    }
    //
    /**
     * Tests the rules in SMIRKS notation: there are eleven SMIRKS for the default and for the strict rules, each with
     * the two mapped atoms of the cleaved bond on the left side and the two link atoms with the bond type number on the
     * right side; the rules for amines, ethers and biaryls have no strict functional group exclusions, all others have.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bondTypeSmirksTest() throws Exception {
        List<String> tmpSmirksList = RECAPFragmenter.getBondTypeSmirks(false);
        List<String> tmpStrictSmirksList = RECAPFragmenter.getBondTypeSmirks(true);
        Assertions.assertEquals(11, tmpSmirksList.size());
        Assertions.assertEquals(11, tmpStrictSmirksList.size());
        Assertions.assertEquals(RECAPFragmenter.getBondTypeNames().size(), tmpSmirksList.size());
        Assertions.assertEquals("[C;$(C=O);!$(C(~[#7])~[#7]):1]-!@[#7;+0:2]>>[*:1]-[1*].[*:2]-[1*]", tmpSmirksList.get(0));
        Assertions.assertEquals("[C;!$(C(=*)=*):1]=!@[C;!$(C(=*)=*):2]>>[*:1]=[6*].[*:2]=[6*]", tmpSmirksList.get(5));
        Assertions.assertEquals("[c:1]-!@[c:2]>>[*:1]-[10*].[*:2]-[10*]", tmpSmirksList.get(9));
        for (int i = 0; i < tmpSmirksList.size(); i++) {
            int tmpBondType = i + 1;
            for (String tmpSmirks : new String[] {tmpSmirksList.get(i), tmpStrictSmirksList.get(i)}) {
                String tmpBond = tmpBondType == 6 ? "=" : "-";
                Assertions.assertTrue(tmpSmirks.endsWith(String.format(">>[*:1]%s[%d*].[*:2]%s[%d*]", tmpBond,
                        tmpBondType, tmpBond, tmpBondType)), tmpSmirks);
                String tmpLeftSide = tmpSmirks.substring(0, tmpSmirks.indexOf(">>"));
                Assertions.assertTrue(tmpLeftSide.contains(":1]" + tmpBond + "!@["), tmpSmirks);
                Assertions.assertTrue(tmpLeftSide.endsWith(":2]"), tmpSmirks);
                //the left side is a valid SMARTS pattern of two atoms
                Assertions.assertDoesNotThrow(() -> SmartsPattern.create(tmpLeftSide), tmpSmirks);
            }
            boolean tmpHasStrictExclusions = tmpBondType != 3 && tmpBondType != 5 && tmpBondType != 10;
            Assertions.assertEquals(tmpHasStrictExclusions, !tmpSmirksList.get(i).equals(tmpStrictSmirksList.get(i)),
                    tmpSmirksList.get(i));
        }
        //the lists cannot be altered
        Assertions.assertThrows(UnsupportedOperationException.class, () -> tmpSmirksList.add(""));
    }
    //
    /**
     * Tests that the SMIRKS of the default rules listed in the class documentation of the RECAPFragmenter are the ones
     * created from the rules in use, so that the documentation cannot deviate from the implementation. Skipped if the
     * source file is not available, i.e. if the test is not run from the project directory.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bondTypeSmirksDocumentationTest() throws Exception {
        Path tmpSourceFile = Path.of("src", "main", "java", "de", "unijena", "cheminf", "mortar", "model",
                "fragmentation", "algorithm", "RECAPFragmenter.java");
        Assumptions.assumeTrue(Files.isRegularFile(tmpSourceFile), "Source file of the RECAPFragmenter not found.");
        List<String> tmpSourceLines = Files.readAllLines(tmpSourceFile, StandardCharsets.UTF_8);
        List<String> tmpSmirksList = RECAPFragmenter.getBondTypeSmirks(false);
        for (int i = 0; i < tmpSmirksList.size(); i++) {
            String tmpDocumentationLine = String.format(" * %2d  %s", i + 1, tmpSmirksList.get(i));
            Assertions.assertTrue(tmpSourceLines.contains(tmpDocumentationLine),
                    "Not in the class documentation: " + tmpDocumentationLine);
        }
        //no further SMIRKS are listed
        long tmpDocumentedSmirks = tmpSourceLines.stream()
                .filter(aLine -> aLine.startsWith(" * ") && aLine.contains(">>[*:1]")).count();
        Assertions.assertEquals(tmpSmirksList.size(), tmpDocumentedSmirks);
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
        this.assertFragments(new RECAPFragmenter(), aSmiles, false, anExpectedFragments);
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments, once with bond types
     * and without stereochemistry, and optionally with stereochemistry as MORTAR encodes it.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param isStereo whether the stereochemistry is compared
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(RECAPFragmenter aFragmenter, String aSmiles, boolean isStereo,
                                 String... anExpectedFragments) throws Exception {
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(aSmiles));
        Assertions.assertEquals(RECAPFragmenterTest.canonical(List.of(anExpectedFragments), false),
                RECAPFragmenterTest.toCanonicalSmiles(tmpFragments, false), aSmiles);
        if (isStereo) {
            Assertions.assertEquals(RECAPFragmenterTest.toStereoSmiles(List.of(anExpectedFragments), false),
                    RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false), aSmiles);
        }
    }
    //
    /**
     * Asserts that the leaves of the RECAP hierarchy of the given molecule are its single-pass fragments, compared with
     * bond types and without stereochemistry.
     *
     * @param aFragmenter fragmenter
     * @param aMolecule molecule
     * @throws Exception if anything goes wrong
     */
    private void assertLeavesAreSinglePassFragments(RECAPFragmenter aFragmenter, IAtomContainer aMolecule)
            throws Exception {
        Assertions.assertEquals(RECAPFragmenterTest.toCanonicalSmiles(aFragmenter.fragmentMolecule(aMolecule), false),
                RECAPFragmenterTest.toCanonicalSmiles(RECAPFragmenterTest.fragmentsOf(
                        aFragmenter.buildHierarchy(aMolecule).getLeaves()), false));
    }
    //
    /**
     * Returns whether the order of a bond of the given molecule is unset.
     *
     * @param aMolecule molecule
     * @return true if a bond has the order UNSET or no order
     */
    private boolean hasUnsetBondOrder(IAtomContainer aMolecule) {
        for (IBond tmpBond : aMolecule.bonds()) {
            if (tmpBond.getOrder() == null || tmpBond.getOrder() == IBond.Order.UNSET) {
                return true;
            }
        }
        return false;
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
     * Returns the fragments of the given hierarchy nodes.
     *
     * @param aNodeList hierarchy nodes
     * @return list of their fragments
     */
    static List<IAtomContainer> fragmentsOf(List<RECAPFragmenter.RECAPHierarchyNode> aNodeList) {
        List<IAtomContainer> tmpFragments = new ArrayList<>(aNodeList.size());
        for (RECAPFragmenter.RECAPHierarchyNode tmpNode : aNodeList) {
            tmpFragments.add(tmpNode.getFragment());
        }
        return tmpFragments;
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
            tmpSmilesList.add(RECAPFragmenterTest.toCanonicalSmiles(tmpFragment, isStereo));
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
            return RECAPFragmenterTest.toStereoSmiles(aFragment, false);
        }
        String tmpSmiles = SmilesGenerator.create(aFragment, SmiFlavor.Unique | SmiFlavor.AtomicMass,
                new int[aFragment.getAtomCount()]);
        return RECAPFragmenterTest.canonical(tmpSmiles, false);
    }
    //
    /**
     * Returns the unique SMILES code with stereochemistry of the given fragment as MORTAR creates it, i.e. with all link
     * atoms written as '*' (the bond types are not encoded). The CDK creates canonical SMILES codes with
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
            tmpSmilesList.add(RECAPFragmenterTest.toStereoSmiles(RECAPFragmenterTest.parse(tmpSmiles),
                    isNonStereogenicElementRemoved));
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
            tmpSmilesList.add(RECAPFragmenterTest.toStereoSmiles(tmpFragment, isNonStereogenicElementRemoved));
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
            tmpSmilesList.add(RECAPFragmenterTest.canonical(tmpSmiles, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns a canonical SMILES code of the given SMILES code that is independent of the aromaticity notation of the
     * input: the SMILES is parsed and kekulized, aromaticity is perceived with the Daylight model (as by the
     * fragmenter), and a unique SMILES with mass numbers (i.e. bond types) and aromatic symbols is created.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether stereochemistry is encoded (as MORTAR encodes it, without bond types)
     * @return canonical SMILES code
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static String canonical(String aSmiles, boolean isStereo) throws CDKException {
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(aSmiles);
        SmartsPattern.prepare(tmpMolecule);
        if (isStereo) {
            return RECAPFragmenterTest.toStereoSmiles(tmpMolecule, false);
        }
        return SmilesGenerator.create(tmpMolecule, SmiFlavor.Unique | SmiFlavor.AtomicMass | SmiFlavor.UseAromaticSymbols,
                new int[tmpMolecule.getAtomCount()]);
    }
}
