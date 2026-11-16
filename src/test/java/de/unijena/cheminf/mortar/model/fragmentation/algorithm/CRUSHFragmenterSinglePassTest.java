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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.List;
import java.util.Locale;
import java.util.Random;

/**
 * Tests the single-pass variant of the {@link CRUSHFragmenter} (default variant). The expected fragments were derived
 * by hand from the definition of the variant (rule groups and their priorities, candidate bonds, processing order, and
 * the size limits for terminal fragments and linkers, see the class documentation of the fragmenter); the rule and the
 * order of the processing are given in the comments. Fragments are compared as canonical SMILES with the rule numbers
 * encoded as mass numbers of the link atoms, e.g. [6*]OCCCCC.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class CRUSHFragmenterSinglePassTest {
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public CRUSHFragmenterSinglePassTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the rule groups and their priorities: functional group rules (priority 1), aromatic ring bond rules
     * (priority 2), and aliphatic skeleton rules (priority 3), as specified for the single-pass variant.
     */
    @Test
    public void ruleGroupTableTest() {
        int[] tmpFunctionalGroupRules = {1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 17, 18, 19, 20, 21, 27, 28};
        int[] tmpAromaticRingBondRules = {10, 12, 14, 15, 22, 23, 24, 29, 30, 31, 32, 33};
        int[] tmpAliphaticSkeletonRules = {8, 16, 25, 26};
        int tmpCount = 0;
        for (int tmpRule : tmpFunctionalGroupRules) {
            Assertions.assertEquals(1, CRUSHFragmenter.getRuleGroupPriority(tmpRule), "rule " + tmpRule);
            tmpCount++;
        }
        for (int tmpRule : tmpAromaticRingBondRules) {
            Assertions.assertEquals(2, CRUSHFragmenter.getRuleGroupPriority(tmpRule), "rule " + tmpRule);
            tmpCount++;
        }
        for (int tmpRule : tmpAliphaticSkeletonRules) {
            Assertions.assertEquals(3, CRUSHFragmenter.getRuleGroupPriority(tmpRule), "rule " + tmpRule);
            tmpCount++;
        }
        Assertions.assertEquals(CRUSHFragmenter.getRules().size(), tmpCount);
        Assertions.assertThrows(IllegalArgumentException.class, () -> CRUSHFragmenter.getRuleGroupPriority(0));
        Assertions.assertThrows(IllegalArgumentException.class, () -> CRUSHFragmenter.getRuleGroupPriority(34));
    }
    //
    /**
     * Tests the functional group rules: the ester bond is cleaved by rule 2; the O-C bond of the alkyl group, matched by
     * the ether rule (6), is not cleaved, since the oxygen atom alone would remain. Without the functional group rules,
     * the molecule is not fragmented.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void functionalGroupRulesTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCCCCC(=O)OCCCCC", "[2*]C(=O)CCCCCCC", "[2*]OCCCCC");
        //thioether (rule 17), amide (rule 1), urea (rule 5), and sulphonamide (rule 13)
        this.assertFragments(tmpFragmenter, "CCCCCCSCCCCC", "[17*]CCCCCC", "[17*]SCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCCC(=O)NCCCCC", "[1*]C(=O)CCCCC", "[1*]NCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCCNC(=O)NCCCCC", "[5*]C(=O)NCCCCC", "[5*]NCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCCS(=O)(=O)NCCCCC", "[13*]S(=O)(=O)CCCCC", "[13*]NCCCCC");
        tmpFragmenter.setFunctionalGroupRulesSetting(false);
        this.assertFragments(tmpFragmenter, "CCCCCCCC(=O)OCCCCC", "CCCCCCCC(=O)OCCCCC");
    }
    //
    /**
     * Tests the aromatic ring bond rules on 4-pentylbiphenyl: the biaryl bond (rule 12) is processed first, then the
     * aryl-alkyl bond (rules 14 and 24, rule 14 is assigned); the phenylene linker has six heavy atoms and two attachment
     * points, which is valid. Without the aromatic ring bond rules, the molecule is not fragmented.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void aromaticRingBondRulesTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCCc1ccc(cc1)-c1ccccc1", "[12*]c1ccccc1", "[12*]c1ccc([14*])cc1",
                "[14*]CCCCC");
        //aryl ether: the O-C(aryl) bond is matched by rules 6 and 30, the functional group rule 6 is assigned; it is
        //processed before the other ether bond, since its smaller part is larger
        this.assertFragments(tmpFragmenter, "CCCCCOc1ccccc1", "[6*]OCCCCC", "[6*]c1ccccc1");
        tmpFragmenter.setAromaticRingBondRulesSetting(false);
        this.assertFragments(tmpFragmenter, "CCCCCc1ccc(cc1)-c1ccccc1", "CCCCCc1ccc(cc1)-c1ccccc1");
        //the aryl ether is still cleaved by the functional group rule 6
        this.assertFragments(tmpFragmenter, "CCCCCOc1ccccc1", "[6*]OCCCCC", "[6*]c1ccccc1");
    }
    //
    /**
     * Tests the aliphatic skeleton rules (off by default): alkanes, olefins, and alkynes are only cleaved if they are
     * enabled; the bonds are processed with the most balanced cleavage first (central bond of dodecane). Double and
     * triple bonds are cleaved by rules 8 and 16, which come before rule 25 (that also matches them).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void aliphaticSkeletonRulesTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCCCCCCCCC", "CCCCCCCCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCCC=CCCCCC", "CCCCCC=CCCCCC");
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        this.assertFragments(tmpFragmenter, "CCCCCCCCCCCC", "[25*]CCCCCC", "[25*]CCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCCC=CCCCCC", "[8*]=CCCCCC", "[8*]=CCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCC#CCCCCC", "[16*]#CCCCC", "[16*]#CCCCCC");
        //with smaller limits, the central bond and then the bonds that leave propyl groups are cleaved
        tmpFragmenter.setMinimumFragmentSizeSetting(3);
        tmpFragmenter.setMinimumLinkerSizeSetting(3);
        this.assertFragments(tmpFragmenter, "CCCCCCCCCCCC", "[25*]CCC", "[25*]CCC[25*]", "[25*]CCC[25*]", "[25*]CCC");
    }
    //
    /**
     * Tests the priority of the rule groups. The C-O bond of a sulfonate ester is matched by rule 26 (aliphatic
     * skeleton) and rule 27 (functional group); rule 27 is assigned, although its number is higher. In
     * 1-phenyl-4-(hexylthio)butane, the thioether bonds (rule 17, functional group) are processed before the aryl-alkyl
     * bond (rule 14, aromatic ring bond): the thioether bond with the larger smaller part is cleaved, and the butyl
     * linker then stays at the phenyl group; processed by rule number, the aryl-alkyl bond would be cleaved instead.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void rulePriorityTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        String tmpSulfonate = "CCCS(=O)(=O)OCCCCC";
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSulfonate);
        tmpFragmenter.prepareMolecule(tmpMolecule);
        int[] tmpCandidates = tmpFragmenter.findSinglePassCandidates(tmpMolecule);
        //O6-C7 is the ester bond (rules 26 and 27), C7-C8 and the other bonds of the alkyl groups are matched by rule 25
        //only, the S3-O6 bond by no rule
        Assertions.assertEquals(27, tmpCandidates[tmpMolecule.getBond(tmpMolecule.getAtom(6), tmpMolecule.getAtom(7))
                .getIndex()]);
        Assertions.assertEquals(25, tmpCandidates[tmpMolecule.getBond(tmpMolecule.getAtom(7), tmpMolecule.getAtom(8))
                .getIndex()]);
        Assertions.assertEquals(0, tmpCandidates[tmpMolecule.getBond(tmpMolecule.getAtom(3), tmpMolecule.getAtom(6))
                .getIndex()]);
        this.assertFragments(tmpFragmenter, tmpSulfonate, "[27*]OS(=O)(=O)CCC", "[27*]CCCCC");
        //the original variant assigns the lowest rule number
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(CRUSHFragmenterTest.createOriginalFragmenter(), tmpSulfonate,
                26, false).contains(CRUSHFragmenterTest.canonical(List.of("[26*]OS(=O)(=O)CCC", "[26*]CCCCC"), false)));
        this.assertFragments(new CRUSHFragmenter(), "c1ccccc1CCCCSCCCCCC", "[17*]CCCCc1ccccc1", "[17*]SCCCCCC");
    }
    //
    /**
     * Tests the minimum fragment size: in propyl butanoate, the ester bond would leave a propoxy group with four heavy
     * atoms, so the molecule is not fragmented with the default minimum fragment size of 5, but with 4 and 3.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void minimumFragmentSizeTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCC(=O)OCCC", "CCCC(=O)OCCC");
        tmpFragmenter.setMinimumFragmentSizeSetting(4);
        this.assertFragments(tmpFragmenter, "CCCC(=O)OCCC", "[2*]C(=O)CCC", "[2*]OCCC");
        //with the minimum size 3, the propyl group is still not cleaved off, since the oxygen atom alone would remain
        tmpFragmenter.setMinimumFragmentSizeSetting(3);
        this.assertFragments(tmpFragmenter, "CCCC(=O)OCCC", "[2*]C(=O)CCC", "[2*]OCCC");
        tmpFragmenter.setMinimumFragmentSizeSetting(6);
        this.assertFragments(tmpFragmenter, "CCCCCCCC(=O)OCCCCC", "[2*]C(=O)CCCCCCC", "[2*]OCCCCC");
        tmpFragmenter.setMinimumFragmentSizeSetting(7);
        this.assertFragments(tmpFragmenter, "CCCCCCCC(=O)OCCCCC", "CCCCCCCC(=O)OCCCCC");
    }
    //
    /**
     * Tests the minimum linker size: in the propanediol diester, the O-C-C-C-O linker has five heavy atoms; with the
     * default minimum linker size of 6, it stays attached to one of the acyl groups (the second ester bond is not
     * cleaved), with a minimum linker size of 5 (or less, i.e. without an extra limit) it is cut out. The minimum
     * fragment size also applies to linkers.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void minimumLinkerSizeTest() throws Exception {
        String tmpDiester = "CCCCCC(=O)OCCCOC(=O)CCCCC";
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, tmpDiester, "[2*]C(=O)CCCCC", "[2*]OCCCOC(=O)CCCCC");
        tmpFragmenter.setMinimumLinkerSizeSetting(5);
        this.assertFragments(tmpFragmenter, tmpDiester, "[2*]C(=O)CCCCC", "[2*]C(=O)CCCCC", "[2*]OCCCO[2*]");
        tmpFragmenter.setMinimumLinkerSizeSetting(1);
        this.assertFragments(tmpFragmenter, tmpDiester, "[2*]C(=O)CCCCC", "[2*]C(=O)CCCCC", "[2*]OCCCO[2*]");
        //the ethanediol linker has four heavy atoms, less than the minimum fragment size
        this.assertFragments(tmpFragmenter, "CCCCCC(=O)OCCOC(=O)CCCCC", "[2*]C(=O)CCCCC", "[2*]OCCOC(=O)CCCCC");
        tmpFragmenter.setMinimumLinkerSizeSetting(7);
        this.assertFragments(tmpFragmenter, "CCCCC(=O)Nc1ccc(cc1)-c1ccccc1", "[1*]C(=O)CCCC", "[1*]Nc1ccc([12*])cc1",
                "[12*]c1ccccc1");
        tmpFragmenter.setMinimumLinkerSizeSetting(8);
        //the amide bond is processed first, then the biaryl bond would leave the aminophenylene linker with seven heavy
        //atoms
        this.assertFragments(tmpFragmenter, "CCCCC(=O)Nc1ccc(cc1)-c1ccccc1", "[1*]C(=O)CCCC",
                "[1*]Nc1ccc(cc1)-c1ccccc1");
    }
    //
    /**
     * Tests that a molecule is not fragmented if all rule groups are disabled, and the setting for returning
     * unfragmented molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleGroupSwitchesTest() throws Exception {
        String tmpSmiles = "CCCCC(=O)Nc1ccc(cc1)-c1ccccc1";
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        Assertions.assertEquals(3, tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(tmpSmiles)).size());
        tmpFragmenter.setFunctionalGroupRulesSetting(false);
        tmpFragmenter.setAromaticRingBondRulesSetting(false);
        tmpFragmenter.setAliphaticSkeletonRulesSetting(false);
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
        tmpFragmenter.prepareMolecule(tmpMolecule);
        for (int tmpRuleNumber : tmpFragmenter.findSinglePassCandidates(tmpMolecule)) {
            Assertions.assertEquals(0, tmpRuleNumber);
        }
        this.assertFragments(tmpFragmenter, tmpSmiles, tmpSmiles);
        this.assertFragments(tmpFragmenter, "CCCCCCCOCCCCC.[Na+].[Cl-]", "CCCCCCCOCCCCC", "[Na+]", "[Cl-]");
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(tmpSmiles)).isEmpty());
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("CCCCCCCOCCCCC.[Na+]")).isEmpty());
        //with a cleaved bond, the fragments are returned regardless of the setting
        CRUSHFragmenter tmpDefaultRulesFragmenter = new CRUSHFragmenter();
        tmpDefaultRulesFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertEquals(2, tmpDefaultRulesFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("CCCCCCCOCCCCC"))
                .size());
    }
    //
    /**
     * Tests the protected atoms: matches that contain a protected atom are dropped, but the other matches of the same
     * rule remain (the conservative mode is not used). In the benzylic ether, the matches of the ether rule (6) at the
     * naphthylmethyl group contain a protected atom, the pentyl ether bonds are cleaved; without protection, the
     * benzylic ether bond is cleaved instead, since its cleavage is more balanced. Propranolol is only fragmented without
     * the protection of the neighbours of the fused ring atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ringProtectionTest() throws Exception {
        String tmpEther = "CCCCCOCCOCc1cccc2ccccc12";
        String tmpPropranolol = "CC(C)NCC(O)COc1cccc2ccccc12";
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, tmpEther, "[6*]OCCCCC", "[6*]CCOCc1cccc2ccccc12");
        this.assertFragments(tmpFragmenter, tmpPropranolol, tmpPropranolol);
        //the conservative mode (setting of the original variant) has no effect
        tmpFragmenter.setConservativeModeSetting(false);
        this.assertFragments(tmpFragmenter, tmpEther, "[6*]OCCCCC", "[6*]CCOCc1cccc2ccccc12");
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.FUSED_RING_ATOMS);
        this.assertFragments(tmpFragmenter, tmpEther, "[6*]OCCOCCCCC", "[6*]Cc1cccc2ccccc12");
        //the match of the ether rule with the fused ring atom as atom 5 is dropped, the other one remains
        this.assertFragments(tmpFragmenter, tmpPropranolol, "[6*]OCC(O)CNC(C)C", "[6*]c1cccc2ccccc12");
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.NONE);
        this.assertFragments(tmpFragmenter, tmpEther, "[6*]OCCOCCCCC", "[6*]Cc1cccc2ccccc12");
        this.assertFragments(tmpFragmenter, tmpPropranolol, "[6*]OCC(O)CNC(C)C", "[6*]c1cccc2ccccc12");
        //ring bonds are never cleaved
        IAtomContainer tmpGlycoside = CRUSHFragmenterTest.parse("CCCCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O");
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        tmpFragmenter.setMinimumFragmentSizeSetting(1);
        tmpFragmenter.setMinimumLinkerSizeSetting(1);
        tmpFragmenter.prepareMolecule(tmpGlycoside);
        int[] tmpCandidates = tmpFragmenter.findSinglePassCandidates(tmpGlycoside);
        for (IBond tmpBond : tmpGlycoside.bonds()) {
            Assertions.assertTrue(tmpCandidates[tmpBond.getIndex()] == 0 || !tmpBond.isInRing());
        }
    }
    //
    /**
     * Tests rule 19 (Carbamate_2), whose two bonds are independent candidates: in hexyl N-ethylcarbamate, the C-N bond
     * (rule 18) would leave an ethylamino group, the O-C bond of rule 19 is cleaved. In pentyl N-pentylcarbamate, the
     * C-N bond is cleaved by rule 18, after which the cleavage of the O-C bond would leave the carbonyl group alone.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void carbamateTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCCCOC(=O)NCC", "[19*]OCCCCCC", "[19*]C(=O)NCC");
        this.assertFragments(tmpFragmenter, "CCCCCOC(=O)NCCCCC", "[18*]C(=O)OCCCCC", "[18*]NCCCCC");
    }
    //
    /**
     * Tests the handling of salts and of the pseudo atoms of the input, which count as attachment points but not as
     * heavy atoms: with a pseudo atom, the butoxy part of butyl pentyl ether becomes a linker with five heavy atoms and
     * is not cleaved off.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void saltAndPseudoAtomTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCCCCOCCCCC.[Na+].[Cl-]", "[6*]CCCCCCC", "[6*]OCCCCC", "[Na+]", "[Cl-]");
        this.assertFragments(tmpFragmenter, "Cc1ccccc1.[Na+]", "Cc1ccccc1", "[Na+]");
        this.assertFragments(tmpFragmenter, "CCCCOCCCCC", "[6*]OCCCC", "[6*]CCCCC");
        this.assertFragments(tmpFragmenter, "*CCCCOCCCCC", "*CCCCOCCCCC");
        this.assertFragments(tmpFragmenter, "*CCCCCOCCCCC", "*CCCCCO[6*]", "[6*]CCCCC");
        tmpFragmenter.setMinimumLinkerSizeSetting(5);
        this.assertFragments(tmpFragmenter, "*CCCCOCCCCC", "*CCCCO[6*]", "[6*]CCCCC");
        //explicit hydrogen atoms do not change the result, the input molecule is not altered
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse("CCCCC(=O)Nc1ccc(cc1)-c1ccccc1");
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        int tmpAtomCount = tmpMolecule.getAtomCount();
        this.assertFragments(new CRUSHFragmenter(), tmpMolecule, "[1*]C(=O)CCCC", "[1*]Nc1ccc([12*])cc1",
                "[12*]c1ccccc1");
        Assertions.assertEquals(tmpAtomCount, tmpMolecule.getAtomCount());
    }
    //
    /**
     * Tests the transfer of tetrahedral and double bond stereo elements: a cleaved ligand is replaced by the link atom
     * (or hydrogen atom) at its position, configurations away from the cleaved bond are kept, and stereo elements of
     * atoms that are not stereogenic anymore are removed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void stereoTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertStereoFragments(tmpFragmenter, "CCCC[C@](C)(CC)OCCCCC", "CCCC[C@](C)(CC)*", "*OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCC[C@@](C)(CC)OCCCCC", "CCCC[C@@](C)(CC)*", "*OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCCC[C@H](O)CCOCCCCC", "CCCCC[C@H](O)CC*", "*OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCCC[C@@H](CCCCCC)OCCCCC", "CCCCC[C@@H](CCCCCC)*", "*OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCCC/C(=C/CCCCC)OCCCCC", "CCCCC/C(=C/CCCCC)*", "*OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCCC/C(=C\\CCCCC)OCCCCC", "CCCCC/C(=C\\CCCCC)*", "*OCCCCC");
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertStereoFragments(tmpFragmenter, "CCCC[C@](C)(CC)OCCCCC", "CCCC[C@H](C)CC", "OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCC[C@@](C)(CC)OCCCCC", "CCCC[C@@H](C)CC", "OCCCCC");
        //two hydrogen atoms: not stereogenic anymore
        this.assertStereoFragments(tmpFragmenter, "CCCCC[C@@H](CCCCCC)OCCCCC", "CCCCCCCCCCCC", "OCCCCC");
        this.assertStereoFragments(tmpFragmenter, "CCCCC/C(=C/CCCCC)OCCCCC", "CCCCC/C=C/CCCCC", "OCCCCC");
    }
    //
    /**
     * Tests the hydrogen saturation of the open valences, also for cleaved double bonds, whose atoms receive two
     * hydrogen atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hydrogenSaturationTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCCCCCC(=O)OCCCCC", "O=CCCCCCCC", "OCCCCC");
        this.assertFragments(tmpFragmenter, "CCCCC(=O)Nc1ccc(cc1)-c1ccccc1", "O=CCCCC", "Nc1ccccc1", "c1ccccc1");
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        this.assertFragments(tmpFragmenter, "CCCCCC=CCCCCC", "CCCCCC", "CCCCCC");
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("CCCCCC=CCCCCC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Assertions.assertFalse(tmpAtom instanceof IPseudoAtom);
            }
        }
    }
    //
    /**
     * Tests molecules with bonds of unset order: if they can be kekulized, the fragments are the same as for the
     * kekulized molecule; otherwise, the molecule is filtered.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void unsetBondOrderTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        for (String tmpSmiles : new String[] {"CCCCC(=O)Nc1ccc(cc1)-c1ccccc1", "CCCCCOc1ccccc1"}) {
            IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpSmiles, false, false);
            boolean tmpIsUnset = false;
            for (IBond tmpBond : tmpMolecule.bonds()) {
                tmpIsUnset |= tmpBond.getOrder() == null || tmpBond.getOrder() == IBond.Order.UNSET;
            }
            Assertions.assertTrue(tmpIsUnset, tmpSmiles);
            Assertions.assertFalse(tmpFragmenter.shouldBeFiltered(tmpMolecule), tmpSmiles);
            List<IAtomContainer> tmpExpected = tmpFragmenter.fragmentMolecule(ChemUtil.parseSmilesToAtomContainer(
                    tmpSmiles, true, true));
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            Assertions.assertTrue(tmpExpected.size() > 1, tmpSmiles);
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpExpected, false),
                    CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false), tmpSmiles);
            Assertions.assertEquals(CRUSHFragmenterReferenceTest.createUniqueSmilesList(tmpExpected),
                    CRUSHFragmenterReferenceTest.createUniqueSmilesList(tmpFragments), tmpSmiles);
        }
        IAtomContainer tmpNonKekulizable = ChemUtil.parseSmilesToAtomContainer("c1ccc(cc1)ccc-c2ccccc2", false, false);
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpNonKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.fragmentMolecule(tmpNonKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.buildRoute(tmpNonKekulizable));
    }
    //
    /**
     * Tests the route of the single-pass variant: the children of the root are the final fragments (including counter
     * ions), the root of a connected molecule without a cleaved bond has no children, the route is exhaustive, and the
     * fragment output setting of the original variant is ignored.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void routeTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        for (String tmpSmiles : new String[] {"CCCCC(=O)Nc1ccc(cc1)-c1ccccc1", "CCCCCCCOCCCCC.[Na+].[Cl-]", "c1ccccc1",
                "Cc1ccccc1.[Na+]"}) {
            IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
            CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(tmpMolecule);
            Assertions.assertTrue(tmpRoot.isRoot());
            Assertions.assertTrue(tmpRoot.isExhaustive());
            Assertions.assertEquals(0, tmpRoot.getRuleNumber());
            List<CRUSHFragmenter.CRUSHRouteNode> tmpLeaves = tmpRoot.isLeaf() ? List.of(tmpRoot) : tmpRoot.getChildren();
            Assertions.assertEquals(tmpRoot.getLeaves().size(), tmpRoot.getAllDescendants().size());
            for (CRUSHFragmenter.CRUSHRouteNode tmpNode : tmpRoot.getChildren()) {
                Assertions.assertTrue(tmpNode.isLeaf());
                Assertions.assertSame(tmpRoot, tmpNode.getParent());
            }
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                    CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(tmpLeaves), false), tmpSmiles);
            CRUSHFragmenter tmpRouteFragmenter = new CRUSHFragmenter();
            tmpRouteFragmenter.setFragmentOutputSetting(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS);
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                    CRUSHFragmenterTest.toCanonicalSmiles(tmpRouteFragmenter.fragmentMolecule(tmpMolecule), false),
                    tmpSmiles);
        }
        Assertions.assertTrue(tmpFragmenter.buildRoute(CRUSHFragmenterTest.parse("c1ccccc1")).isLeaf());
        Assertions.assertEquals(4, tmpFragmenter.buildRoute(CRUSHFragmenterTest.parse("CCCCCCCOCCCCC.[Na+].[Cl-]"))
                .getChildren().size());
    }
    //
    /**
     * Tests the switch between the variants and that a copy of the fragmenter with non-default settings gives the same
     * fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void variantAndCopyTest() throws Exception {
        String tmpSmiles = "CCCCCCCCC";
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        this.assertFragments(tmpFragmenter, tmpSmiles, tmpSmiles);
        tmpFragmenter.setVariantSetting(CRUSHFragmenter.VariantOption.ORIGINAL);
        tmpFragmenter.setMinimumFragmentSizeSetting(3);
        //route with the maximum number of fragments, see CRUSHFragmenterTest.routeTest()
        this.assertFragments(tmpFragmenter, tmpSmiles, "[25*]CCC", "[25*]CCC[25*]", "[25*]CCC");
        //the single pass prefers balanced cleavages: after the cleavage of the central bonds (equivalent by symmetry),
        //a further cleavage would leave a part with less than three heavy atoms
        tmpFragmenter.setVariantSetting(CRUSHFragmenter.VariantOption.SINGLE_PASS);
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        tmpFragmenter.setMinimumLinkerSizeSetting(3);
        this.assertFragments(tmpFragmenter, tmpSmiles, "[25*]CCCC", "[25*]CCCCC");
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        CRUSHFragmenter tmpCopy = (CRUSHFragmenter) tmpFragmenter.copy();
        for (String tmpMolecule : new String[] {tmpSmiles, "CCCCCCOC(=O)NCC", "CCCCCC(=O)OCCCOC(=O)CCCCC"}) {
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragmenter.fragmentMolecule(
                    CRUSHFragmenterTest.parse(tmpMolecule)), false), CRUSHFragmenterTest.toCanonicalSmiles(
                    tmpCopy.fragmentMolecule(CRUSHFragmenterTest.parse(tmpMolecule)), false), tmpMolecule);
        }
        this.assertFragments(tmpCopy, tmpSmiles, "CCCC", "CCCCC");
    }
    //
    /**
     * Tests that the fragments do not depend on the order of the atoms and bonds of the input, including molecules
     * whose candidate bonds are tied in the first three sort keys: symmetric diesters and meso diesters, in which the
     * two constitutionally equivalent ester bonds give enantiomeric fragments (only one of them can be cleaved, since
     * the ethanediol linker is too small).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void atomOrderInvarianceTest() throws Exception {
        String[] tmpMolecules = {"CC[C@@H](C)C(=O)OCCOC(=O)[C@@H](C)CC", "CC[C@@H](C)C(=O)OCCOC(=O)[C@H](C)CC",
                "CCCC[C@@H](C)OCC(=O)OCCOC(=O)CO[C@@H](C)CCCC", "CCCCCC(=O)OCCCOC(=O)CCCCC",
                "CC(C)NCC(O)COc1cccc2ccccc12", "Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O",
                "COC1=CC=CC2=C1C(=O)C1=C(O)C3=C(C[C@@](O)(C(=O)CO)C[C@@H]3O[C@H]3C[C@H](N)C[C@H](C)O3)C(O)=C1C2=O",
                "NC1=NC(N)=C2N=C(CNC3=CC=C(C(=O)N[C@@H](CCCNC(=O)C4=CC=CC=C4C(=O)O)C(=O)O)C=C3)C=NC2=N1",
                "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(COC(N)=O)CS[C@H]12)C1=CC=CO1.[Na+]"};
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        CRUSHFragmenter tmpAliphaticFragmenter = new CRUSHFragmenter();
        tmpAliphaticFragmenter.setAliphaticSkeletonRulesSetting(true);
        tmpAliphaticFragmenter.setMinimumFragmentSizeSetting(3);
        Random tmpRandom = new Random(42);
        for (CRUSHFragmenter tmpCurrentFragmenter : new CRUSHFragmenter[] {tmpFragmenter, tmpAliphaticFragmenter}) {
            for (String tmpSmiles : tmpMolecules) {
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
                List<IAtomContainer> tmpFragments = tmpCurrentFragmenter.fragmentMolecule(tmpMolecule);
                List<String> tmpExpected = CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false);
                List<String> tmpExpectedKeys = CRUSHFragmenterReferenceTest.createUniqueSmilesList(tmpFragments);
                for (int i = 0; i < 10; i++) {
                    List<IAtomContainer> tmpShuffledFragments = tmpCurrentFragmenter.fragmentMolecule(
                            CRUSHFragmenterTest.shuffle(tmpMolecule, tmpRandom));
                    Assertions.assertEquals(tmpExpected, CRUSHFragmenterTest.toCanonicalSmiles(tmpShuffledFragments,
                            false), tmpSmiles);
                    Assertions.assertEquals(tmpExpectedKeys, CRUSHFragmenterReferenceTest.createUniqueSmilesList(
                            tmpShuffledFragments), tmpSmiles);
                }
            }
        }
        //the meso diester: one acyl group is cleaved off, the other one stays at the linker, with the opposite
        //configuration
        this.assertStereoFragments(tmpFragmenter, "CC[C@@H](C)C(=O)OCCOC(=O)[C@@H](C)CC", "CC[C@@H](C)C(=O)*",
                "CC[C@H](C)C(=O)OCCO*");
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments, without
     * stereochemistry and with the rule numbers of the link atoms.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(CRUSHFragmenter aFragmenter, String aSmiles, String... anExpectedFragments)
            throws Exception {
        this.assertFragments(aFragmenter, CRUSHFragmenterTest.parse(aSmiles), anExpectedFragments);
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments, without
     * stereochemistry and with the rule numbers of the link atoms.
     *
     * @param aFragmenter fragmenter
     * @param aMolecule the molecule
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(CRUSHFragmenter aFragmenter, IAtomContainer aMolecule, String... anExpectedFragments)
            throws Exception {
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(aMolecule);
        Assertions.assertEquals(CRUSHFragmenterTest.canonical(List.of(anExpectedFragments), false),
                CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false));
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments with stereochemistry,
     * as MORTAR encodes it (link atoms as '*').
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertStereoFragments(CRUSHFragmenter aFragmenter, String aSmiles, String... anExpectedFragments)
            throws Exception {
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(aSmiles));
        Assertions.assertEquals(CRUSHFragmenterTest.toStereoSmiles(List.of(anExpectedFragments), false),
                CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false), aSmiles);
    }
}
