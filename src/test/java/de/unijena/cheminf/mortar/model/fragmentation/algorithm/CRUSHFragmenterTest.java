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
import de.unijena.cheminf.mortar.preference.PreferenceUtil;

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;

/**
 * Tests the settings and the original variant of the {@link CRUSHFragmenter}, with the minimum fragment size 3 of the
 * publication (see {@link #createOriginalFragmenter()}); the single-pass variant is tested in
 * {@link CRUSHFragmenterSinglePassTest}. The expected fragments were derived by hand from the CRUSH rules; the rule that
 * cleaves a bond is given in the comments. Fragments are compared as canonical SMILES with the rule numbers encoded as
 * mass numbers of the link atoms, e.g. [26*]OCCC.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class CRUSHFragmenterTest {
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public CRUSHFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Returns a fragmenter with the original variant and the minimum fragment size of the publication (3), on which the
     * tests of the original variant and the reference files of the steps, routes, and rule cases are based; all other
     * settings have their default values.
     *
     * @return fragmenter
     */
    static CRUSHFragmenter createOriginalFragmenter() {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        tmpFragmenter.setVariantSetting(CRUSHFragmenter.VariantOption.ORIGINAL);
        tmpFragmenter.setMinimumFragmentSizeSetting(CRUSHFragmenter.PUBLISHED_MINIMUM_FRAGMENT_SIZE);
        return tmpFragmenter;
    }
    //
    /**
     * Tests the settings, their order and persistence names, copying, the restoration of the default settings, the
     * preference restrictions, and the parameter checks.
     */
    @Test
    public void basicSettingsTest() {
        CRUSHFragmenter tmpFragmenter = new CRUSHFragmenter();
        Assertions.assertEquals(CRUSHFragmenter.ALGORITHM_NAME, tmpFragmenter.getFragmentationAlgorithmName());
        Assertions.assertEquals("CRUSH", tmpFragmenter.getFragmentationAlgorithmDisplayName());
        //the names of the settings are used for their persistence and must not change
        List<String> tmpNames = new ArrayList<>();
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            tmpNames.add(tmpSetting.getName());
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
            Assertions.assertFalse(tmpFragmenter.getSettingNameToTooltipTextMap().get(tmpSetting.getName()).startsWith("!"));
            //Message returns "Key ... not found." for missing keys
            Assertions.assertFalse(tmpFragmenter.getSettingNameToTooltipTextMap().get(tmpSetting.getName()).startsWith(
                    "Key"));
            Assertions.assertFalse(tmpFragmenter.getSettingNameToDisplayNameMap().get(tmpSetting.getName()).startsWith(
                    "Key"));
        }
        Assertions.assertEquals(List.of("Algorithm variant setting", "Fragment saturation setting",
                "Minimum fragment size setting", "Minimum linker size setting", "Functional group rules setting",
                "Aromatic ring bond rules setting", "Aliphatic skeleton rules setting", "Ring protection setting",
                "Conservative mode setting", "Fragment output setting", "Return unfragmented molecules setting"), tmpNames);
        Assertions.assertEquals(11, tmpFragmenter.getSettingNameToDisplayNameMap().size());
        Assertions.assertEquals(11, tmpFragmenter.getSettingNameToTooltipTextMap().size());
        for (CRUSHFragmenter.VariantOption tmpOption : CRUSHFragmenter.VariantOption.values()) {
            Assertions.assertFalse(tmpOption.getDisplayName().startsWith("Key"));
            Assertions.assertFalse(tmpOption.getTooltipText().startsWith("Key"));
        }
        //the published procedure is reproduced with the original variant and the minimum fragment size 3, as stated in
        //the tooltips
        Assertions.assertTrue(CRUSHFragmenter.VariantOption.ORIGINAL.getTooltipText().contains(
                "minimum fragment size of 3"));
        Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().get("Minimum fragment size setting")
                .contains("variant original and a minimum fragment size of 3"));
        Assertions.assertEquals(CRUSHFragmenter.VariantOption.SINGLE_PASS, tmpFragmenter.getVariantSetting());
        Assertions.assertEquals(CRUSHFragmenter.SaturationOption.LINK_ATOMS, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(5, tmpFragmenter.getMinimumFragmentSizeSetting());
        Assertions.assertEquals(6, tmpFragmenter.getMinimumLinkerSizeSetting());
        Assertions.assertTrue(tmpFragmenter.getFunctionalGroupRulesSetting());
        Assertions.assertTrue(tmpFragmenter.getAromaticRingBondRulesSetting());
        Assertions.assertFalse(tmpFragmenter.getAliphaticSkeletonRulesSetting());
        Assertions.assertEquals(CRUSHFragmenter.RingProtectionOption.FUSED_RING_ATOMS_AND_NEIGHBOURS,
                tmpFragmenter.getRingProtectionSetting());
        Assertions.assertTrue(tmpFragmenter.getConservativeModeSetting());
        Assertions.assertEquals(CRUSHFragmenter.FragmentOutputOption.FINAL_FRAGMENTS, tmpFragmenter.getFragmentOutputSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(
                tmpFragmenter.settingsProperties()));
        tmpFragmenter.setVariantSetting(CRUSHFragmenter.VariantOption.ORIGINAL);
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        tmpFragmenter.setMinimumFragmentSizeSetting(2);
        tmpFragmenter.setMinimumLinkerSizeSetting(9);
        tmpFragmenter.setFunctionalGroupRulesSetting(false);
        tmpFragmenter.setAromaticRingBondRulesSetting(false);
        tmpFragmenter.setAliphaticSkeletonRulesSetting(true);
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.NONE);
        tmpFragmenter.setConservativeModeSetting(false);
        tmpFragmenter.setFragmentOutputSetting(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS);
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(
                tmpFragmenter.settingsProperties()));
        CRUSHFragmenter tmpCopy = (CRUSHFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(CRUSHFragmenter.VariantOption.ORIGINAL, tmpCopy.getVariantSetting());
        Assertions.assertEquals(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION, tmpCopy.getSaturationSetting());
        Assertions.assertEquals(2, tmpCopy.getMinimumFragmentSizeSetting());
        Assertions.assertEquals(9, tmpCopy.getMinimumLinkerSizeSetting());
        Assertions.assertFalse(tmpCopy.getFunctionalGroupRulesSetting());
        Assertions.assertFalse(tmpCopy.getAromaticRingBondRulesSetting());
        Assertions.assertTrue(tmpCopy.getAliphaticSkeletonRulesSetting());
        Assertions.assertEquals(CRUSHFragmenter.RingProtectionOption.NONE, tmpCopy.getRingProtectionSetting());
        Assertions.assertFalse(tmpCopy.getConservativeModeSetting());
        Assertions.assertEquals(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS, tmpCopy.getFragmentOutputSetting());
        Assertions.assertFalse(tmpCopy.getReturnUnfragmentedMoleculesSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(CRUSHFragmenter.VARIANT_OPTION_DEFAULT, tmpFragmenter.getVariantSetting());
        Assertions.assertEquals(CRUSHFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(CRUSHFragmenter.MINIMUM_FRAGMENT_SIZE_DEFAULT, tmpFragmenter.getMinimumFragmentSizeSetting());
        Assertions.assertEquals(CRUSHFragmenter.MINIMUM_LINKER_SIZE_DEFAULT, tmpFragmenter.getMinimumLinkerSizeSetting());
        Assertions.assertEquals(CRUSHFragmenter.FUNCTIONAL_GROUP_RULES_OPTION_DEFAULT,
                tmpFragmenter.getFunctionalGroupRulesSetting());
        Assertions.assertEquals(CRUSHFragmenter.AROMATIC_RING_BOND_RULES_OPTION_DEFAULT,
                tmpFragmenter.getAromaticRingBondRulesSetting());
        Assertions.assertEquals(CRUSHFragmenter.ALIPHATIC_SKELETON_RULES_OPTION_DEFAULT,
                tmpFragmenter.getAliphaticSkeletonRulesSetting());
        Assertions.assertEquals(CRUSHFragmenter.RING_PROTECTION_OPTION_DEFAULT, tmpFragmenter.getRingProtectionSetting());
        Assertions.assertEquals(CRUSHFragmenter.CONSERVATIVE_MODE_OPTION_DEFAULT, tmpFragmenter.getConservativeModeSetting());
        Assertions.assertEquals(CRUSHFragmenter.FRAGMENT_OUTPUT_OPTION_DEFAULT, tmpFragmenter.getFragmentOutputSetting());
        Assertions.assertEquals(CRUSHFragmenter.RETURN_UNFRAGMENTED_MOLECULES_OPTION_DEFAULT,
                tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        //the copy is independent of the original
        Assertions.assertEquals(2, tmpCopy.getMinimumFragmentSizeSetting());
        Assertions.assertEquals(CRUSHFragmenter.VariantOption.ORIGINAL, tmpCopy.getVariantSetting());
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setVariantSetting(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setSaturationSetting(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setRingProtectionSetting(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setFragmentOutputSetting(null));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setMinimumFragmentSizeSetting(0));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setMinimumLinkerSizeSetting(0));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setMinimumLinkerSizeSetting(-1));
        Assertions.assertEquals(CRUSHFragmenter.MINIMUM_LINKER_SIZE_DEFAULT, tmpFragmenter.getMinimumLinkerSizeSetting());
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(null));
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.fragmentMolecule(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.buildRoute(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.shouldBePreprocessed(null));
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.canBeFragmented(null));
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.fragmentMolecule(SilentChemObjectBuilder.getInstance().newAtomContainer()));
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> tmpFragmenter.applyPreprocessing(SilentChemObjectBuilder.getInstance().newAtomContainer()));
    }
    //
    /**
     * Tests the rule table: 33 rules in the order of Table 1 of the publication, with the SMARTS (reactant sides) and
     * cleaved bonds of the rule file of the authors (verbatim copy CRUSH_SMIRKS_Rules.txt in the test resources), where
     * the cleaved bonds are derived from the SMIRKS as the bonds between mapped atoms that are missing on the product
     * side, and the cleaved bonds as listed in the implementation plan.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleTableTest() throws Exception {
        List<CRUSHFragmenter.CleavageRule> tmpRules = CRUSHFragmenter.getRules();
        Assertions.assertEquals(33, tmpRules.size());
        Assertions.assertEquals(33, new HashSet<>(CRUSHFragmenter.getRuleNames()).size());
        List<String[]> tmpFileRules = new ArrayList<>();
        try (InputStream tmpStream = CRUSHFragmenterTest.class.getResourceAsStream("CRUSH_SMIRKS_Rules.txt")) {
            Objects.requireNonNull(tmpStream, "Rule file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            java.util.regex.Pattern tmpRulePattern = java.util.regex.Pattern.compile("\\('([^']+)',\\s*'([^']+)'\\)");
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.startsWith("#")) {
                    continue;
                }
                Matcher tmpMatcher = tmpRulePattern.matcher(tmpLine);
                if (tmpMatcher.find()) {
                    tmpFileRules.add(new String[] {tmpMatcher.group(1), tmpMatcher.group(2)});
                }
            }
        }
        Assertions.assertEquals(33, tmpFileRules.size());
        //cleaved bonds (map numbers) as listed in the implementation plan (section 3.1)
        int[][][] tmpExpectedCleavedBonds = {{{2, 4}}, {{2, 4}}, {{1, 4}}, {{1, 2}}, {{2, 4}}, {{3, 4}}, {{1, 2}},
                {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 3}}, {{1, 2}}, {{1, 4}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}},
                {{2, 4}}, {{1, 2}, {2, 4}}, {{5, 6}}, {{1, 4}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}},
                {{4, 5}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}};
        for (int i = 0; i < 33; i++) {
            CRUSHFragmenter.CleavageRule tmpRule = tmpRules.get(i);
            String[] tmpSmirks = tmpFileRules.get(i)[1].split(">>");
            Assertions.assertEquals(i + 1, tmpRule.number());
            Assertions.assertEquals(tmpFileRules.get(i)[0], tmpRule.name());
            Assertions.assertEquals(tmpSmirks[0], tmpRule.smarts());
            Set<List<Integer>> tmpDerivedBonds = CRUSHFragmenterTest.mappedBonds(tmpSmirks[0]);
            tmpDerivedBonds.removeAll(CRUSHFragmenterTest.mappedBonds(tmpSmirks[1]));
            Set<List<Integer>> tmpRuleBonds = new HashSet<>();
            for (int[] tmpPair : tmpRule.cleavedBonds()) {
                tmpRuleBonds.add(List.of(tmpPair[0], tmpPair[1]));
            }
            Assertions.assertEquals(tmpDerivedBonds, tmpRuleBonds, tmpRule.name());
            Set<List<Integer>> tmpPlanBonds = new HashSet<>();
            for (int[] tmpPair : tmpExpectedCleavedBonds[i]) {
                tmpPlanBonds.add(List.of(tmpPair[0], tmpPair[1]));
            }
            Assertions.assertEquals(tmpPlanBonds, tmpRuleBonds, tmpRule.name());
        }
    }
    //
    /**
     * Tests every rule on a minimal example: the expected step of the rule with its products, once with link atoms
     * (carrying the rule number) and once saturated with hydrogen atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleExamplesTest() throws Exception {
        //rule number, molecule, products with link atoms, products with hydrogen atoms
        String[][] tmpCases = {
                {"1", "CCCC(=O)NCCC", "[1*]C(=O)CCC.[1*]NCCC", "CCCC=O.NCCC"},
                {"2", "CCCC(=O)OCCC", "[2*]C(=O)CCC.[2*]OCCC", "CCCC=O.OCCC"},
                {"3", "CCCN(CCC)CCC", "[3*]N(CCC)CCC.[3*]CCC", "CCCNCCC.CCC"},
                {"4", "CCCN1CCOCC1", "[4*]N1CCOCC1.[4*]CCC", "C1COCCN1.CCC"},
                {"5", "CCCNC(=O)NCCC", "[5*]C(=O)NCCC.[5*]NCCC", "O=CNCCC.NCCC"},
                {"6", "CCCOCCC", "[6*]OCCC.[6*]CCC", "OCCC.CCC"},
                {"7", "CCCCC1CCCO1", "[7*]C1CCCO1.[7*]CCCC", "C1CCOC1.CCCC"},
                {"8", "CCCC=CCCC", "[8*]=CCCC.[8*]=CCCC", "CCCC.CCCC"},
                {"9", "CCCC[N+](C)(C)CCCC", "[9*][N+](C)(C)CCCC.[9*]CCCC", "C[NH+](C)CCCC.CCCC"},
                {"10", "CCCCn1ccnc1", "[10*]n1ccnc1.[10*]CCCC", "c1c[nH]cn1.CCCC"},
                {"11", "CCCN1CCCC1=O", "[11*]N1CCCC1=O.[11*]CCC", "O=C1CCCN1.CCC"},
                {"12", "c1ccc(-c2ccccc2)cc1", "[12*]c1ccccc1.[12*]c1ccccc1", "c1ccccc1.c1ccccc1"},
                {"13", "CCCS(=O)(=O)NCCC", "[13*]S(=O)(=O)CCC.[13*]NCCC", "CCC[SH](=O)=O.NCCC"},
                {"14", "CCCc1ccccc1", "[14*]c1ccccc1.[14*]CCC", "c1ccccc1.CCC"},
                {"15", "CCCNc1ccccc1", "[15*]c1ccccc1.[15*]NCCC", "c1ccccc1.NCCC"},
                {"16", "CCCC#CCCC", "[16*]#CCCC.[16*]#CCCC", "CCCC.CCCC"},
                {"17", "CCCSCCC", "[17*]CCC.[17*]SCCC", "CCC.SCCC"},
                {"18", "CCCOC(=O)NCCC", "[18*]C(=O)OCCC.[18*]NCCC", "O=COCCC.NCCC"},
                {"20", "CCCS(=O)(=O)NC(=O)CCC", "[20*]NS(=O)(=O)CCC.[20*]C(=O)CCC", "NS(=O)(=O)CCC.O=CCCC"},
                {"21", "CCCS(=O)(=O)Oc1ccccc1", "[21*]S(=O)(=O)CCC.[21*]Oc1ccccc1", "CCC[SH](=O)=O.Oc1ccccc1"},
                {"22", "c1ccc(-n2cccc2)cc1", "[22*]c1ccccc1.[22*]n1cccc1", "c1ccccc1.c1cc[nH]c1"},
                {"23", "c1ccn(-n2cccc2)c1", "[23*]n1cccc1.[23*]n1cccc1", "c1cc[nH]c1.c1cc[nH]c1"},
                {"24", "CCCc1ccccc1", "[24*]c1ccccc1.[24*]CCC", "c1ccccc1.CCC"},
                {"25", "CCCCCC", "[25*]CCC.[25*]CCC", "CCC.CCC"},
                {"26", "CCCOCCC", "[26*]OCCC.[26*]CCC", "OCCC.CCC"},
                {"27", "CCCS(=O)(=O)OCCC", "[27*]OS(=O)(=O)CCC.[27*]CCC", "OS(=O)(=O)CCC.CCC"},
                {"28", "CCCSSCCC", "[28*]SCCC.[28*]SCCC", "SCCC.SCCC"},
                {"29", "CCC=Cc1ccccc1", "[29*]c1ccccc1.[29*]C=CCC", "c1ccccc1.C=CCC"},
                {"30", "CCCOc1ccccc1", "[30*]c1ccccc1.[30*]OCCC", "c1ccccc1.OCCC"},
                {"31", "CCCNc1ccccc1", "[31*]c1ccccc1.[31*]NCCC", "c1ccccc1.NCCC"},
                {"32", "CCCSc1ccccc1", "[32*]c1ccccc1.[32*]SCCC", "c1ccccc1.SCCC"},
                {"33", "CCCC#Cc1ccccc1", "[33*]c1ccccc1.[33*]C#CCCC", "c1ccccc1.C#CCCC"}};
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        Set<Integer> tmpTestedRules = new HashSet<>();
        for (String[] tmpCase : tmpCases) {
            int tmpRule = Integer.parseInt(tmpCase[0]);
            List<List<String>> tmpProductSets = CRUSHFragmenterTest.stepProducts(tmpFragmenter, tmpCase[1], tmpRule, false);
            Assertions.assertTrue(tmpProductSets.contains(CRUSHFragmenterTest.canonical(List.of(tmpCase[2].split("\\.")),
                    false)), tmpCase[1] + " " + tmpProductSets);
            Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpFragmenter, tmpCase[1], tmpRule, true).contains(
                    CRUSHFragmenterTest.canonical(List.of(tmpCase[3].split("\\.")), false)), tmpCase[1]);
            tmpTestedRules.add(tmpRule);
        }
        //rule 19 (Carbamate_2) cleaves two bonds and cuts off a carbonyl group with two heavy atoms, so it never
        //applies with the minimum fragment size of 3 of the publication, but with 2
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCOC(=O)NCCC", 19, false).isEmpty());
        tmpFragmenter.setMinimumFragmentSizeSetting(2);
        Assertions.assertEquals(List.of(CRUSHFragmenterTest.canonical(List.of("[19*]OCCC", "[19*]C([19*])=O",
                "[19*]NCCC"), false)), CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCOC(=O)NCCC", 19, false));
        Assertions.assertEquals(List.of(CRUSHFragmenterTest.canonical(List.of("OCCC", "C=O", "NCCC"), false)),
                CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCOC(=O)NCCC", 19, true));
        tmpTestedRules.add(19);
        Assertions.assertEquals(33, tmpTestedRules.size());
    }
    //
    /**
     * Tests the cleavage of bonds of higher order: double bonds (rule 8 and rule 25), triple bonds (rule 16), and the
     * C=N bond of an oxime ether (rule 26); the link atoms are bonded with the order of the cleaved bond, and the atoms
     * receive one hydrogen atom per bond order with hydrogen saturation.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bondOrderTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        this.assertFragments(tmpFragmenter, "CCCC=CCCC", "[8*]=CCCC", "[8*]=CCCC");
        this.assertFragments(tmpFragmenter, "CCCC#CCCC", "[16*]#CCCC", "[16*]#CCCC");
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCC=NOCCC", 26, false).contains(
                CRUSHFragmenterTest.canonical(List.of("[26*]=CCCC", "[26*]=NOCCC"), false)));
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCC=CCCC", "CCCC", "CCCC");
        this.assertFragments(tmpFragmenter, "CCCC#CCCC", "CCCC", "CCCC");
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCC=NOCCC", 26, true).contains(
                CRUSHFragmenterTest.canonical(List.of("CCCC", "NOCCC"), false)));
        for (IAtomContainer tmpFragment : CRUSHFragmenterTest.createOriginalFragmenter().fragmentMolecule(
                CRUSHFragmenterTest.parse("CCCC=CCCC"))) {
            for (IBond tmpBond : tmpFragment.bonds()) {
                if (tmpBond.getBegin() instanceof IPseudoAtom || tmpBond.getEnd() instanceof IPseudoAtom) {
                    Assertions.assertEquals(IBond.Order.DOUBLE, tmpBond.getOrder());
                }
            }
        }
    }
    //
    /**
     * Tests the route search on molecules with a unique optimal decomposition (from the BIOFACQUIM reference file) and
     * the tie-breaking by the lowest rule number.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void routeTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        //two steps with two fragments each: rule 1 (amide bond, also matched by rule 26) wins against rule 26 (N-propyl)
        this.assertFragments(tmpFragmenter, "CCCC(=O)NCCC", "[1*]C(=O)CCC", "[1*]NCCC");
        //only the central bond of hexane yields two fragments with three heavy atoms
        this.assertFragments(tmpFragmenter, "CCCCCC", "[25*]CCC", "[25*]CCC");
        this.assertFragments(tmpFragmenter, "CCCCCCCCC", "[25*]CCC", "[25*]CCC[25*]", "[25*]CCC");
        //CNP0080333.1: the aryl-alkyl bond is matched by rules 14 and 24, rule 14 is assigned
        this.assertFragments(tmpFragmenter, "Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O", "[14*]c1ccc(C)c(O)c1",
                "[14*]C(C)C[25*]", "[25*]CC(O)[25*]", "[25*]C(C)(C)O");
        CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(CRUSHFragmenterTest.parse(
                "Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O"));
        Assertions.assertTrue(tmpRoot.isExhaustive());
        Assertions.assertTrue(tmpRoot.isRoot());
        Assertions.assertNull(tmpRoot.getParent());
        Assertions.assertEquals(4, tmpRoot.getLeaves().size());
        Assertions.assertEquals(6, tmpRoot.getAllDescendants().size());
        for (CRUSHFragmenter.CRUSHRouteNode tmpNode : tmpRoot.getAllDescendants()) {
            Assertions.assertFalse(tmpNode.isRoot());
            Assertions.assertTrue(tmpNode.getParent().getChildren().contains(tmpNode));
            Assertions.assertEquals(tmpNode.isLeaf(), tmpNode.getRuleNumber() == 0);
        }
        //CNP0250368.0 with hydrogen saturation, as in the reference file
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "Cc1ccc(C(=O)CCC(=O)c2ccc(C)cc2O)c(O)c1", "CC=O", "CC=O", "Cc1cccc(O)c1",
                "Cc1cccc(O)c1");
        this.assertFragments(tmpFragmenter, "CCCCCCCCCCCCCCCc1cccc(O)c1C(=O)O", "CCC", "CCC", "CCC", "CCC", "CCC",
                "O=CO", "Oc1ccccc1");
    }
    //
    /**
     * Tests the output options: the final fragments contain every heavy atom of the molecule exactly once and two link
     * atoms per cleaved bond; all route fragments are the fragments of all nodes of the route except the root and, for
     * molecules with several connected components, except the components that are cleaved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void outputOptionTest() throws Exception {
        String[] tmpMolecules = {"Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O", "CC(C)NCC(O)COc1cccc2ccccc12",
                "CCCCCCCC(=O)OCCCC.[K+].[Br-]", "CCCOC(=O)NCCC"};
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        CRUSHFragmenter tmpRouteFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        tmpRouteFragmenter.setFragmentOutputSetting(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS);
        for (String tmpSmiles : tmpMolecules) {
            IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(tmpMolecule);
            int tmpCleavedBonds = 0;
            for (CRUSHFragmenter.CRUSHRouteNode tmpNode : tmpRoot.getAllDescendants()) {
                if (!tmpNode.isLeaf()) {
                    tmpCleavedBonds += tmpNode.getChildren().size() - 1;
                }
            }
            if (!tmpRoot.isLeaf() && tmpRoot.getRuleNumber() != 0) {
                tmpCleavedBonds += tmpRoot.getChildren().size() - 1;
            }
            int tmpHeavyAtoms = 0;
            int tmpLinkAtoms = 0;
            for (IAtomContainer tmpFragment : tmpFragments) {
                for (IAtom tmpAtom : tmpFragment.atoms()) {
                    if (tmpAtom instanceof IPseudoAtom) {
                        tmpLinkAtoms++;
                        Assertions.assertEquals((Integer) tmpAtom.getProperty(CRUSHFragmenter.RULE_NUMBER_PROPERTY_KEY),
                                tmpAtom.getMassNumber());
                        Assertions.assertEquals(0, tmpAtom.getImplicitHydrogenCount());
                        Assertions.assertEquals(1, tmpFragment.getConnectedBondsCount(tmpAtom));
                    } else {
                        tmpHeavyAtoms++;
                    }
                }
            }
            Assertions.assertEquals(tmpMolecule.getAtomCount(), tmpHeavyAtoms, tmpSmiles);
            Assertions.assertEquals(2 * tmpCleavedBonds, tmpLinkAtoms, tmpSmiles);
            //the leaves of the route are the final fragments
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                    CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(tmpRoot.getLeaves()), false));
            Assertions.assertEquals(CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false),
                    CRUSHFragmenterTest.toStereoSmilesOfFragments(CRUSHFragmenterTest.fragmentsOf(tmpRoot.getLeaves()), false));
            //all route fragments are all nodes except the root and the cleaved components, including the intermediate
            //fragments
            Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(
                    CRUSHFragmenterTest.allRouteFragmentNodes(tmpRoot)), false),
                    CRUSHFragmenterTest.toCanonicalSmiles(tmpRouteFragmenter.fragmentMolecule(tmpMolecule), false));
        }
        //intermediate fragments of the route of nonane (the parts of the first step)
        this.assertFragments(tmpRouteFragmenter, "CCCCCCCCC", "[25*]CCC", "[25*]CCCCCC", "[25*]CCC", "[25*]CCC[25*]");
        //salts: the components are children of the root; like the molecule itself, a cleaved component is not returned,
        //but the components without a cleavage step are
        CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(CRUSHFragmenterTest.parse("CCCCCC.[Na+].[Cl-]"));
        Assertions.assertEquals(3, tmpRoot.getChildren().size());
        Assertions.assertEquals(0, tmpRoot.getRuleNumber());
        this.assertFragments(tmpRouteFragmenter, "CCCCCC.[Na+].[Cl-]", "[25*]CCC", "[25*]CCC", "[Na+]", "[Cl-]");
        //a salt gives the same route fragments as the free base, apart from the counter ion
        for (String tmpSmiles : new String[] {"CCCN1CCCCC1", "CC(C)NCC(O)COc1cccc2ccccc12", "CCCCCCCC(=O)OCCCC"}) {
            List<String> tmpExpectedFragments = new ArrayList<>(CRUSHFragmenterTest.toCanonicalSmiles(
                    tmpRouteFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(tmpSmiles)), false));
            tmpExpectedFragments.add(CRUSHFragmenterTest.canonical("Cl", false));
            Collections.sort(tmpExpectedFragments);
            Assertions.assertEquals(tmpExpectedFragments, CRUSHFragmenterTest.toCanonicalSmiles(
                    tmpRouteFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(tmpSmiles + ".Cl")), false), tmpSmiles);
        }
    }
    //
    /**
     * Tests the protected atoms (atoms in at least two rings and their neighbours) and the conservative mode: in
     * 1-naphthylmethyl propyl ether, the benzylic carbon atom is bonded to a protected atom, so the ether rule (6) is
     * skipped as a whole in conservative mode, while only its matches on protected atoms are dropped otherwise.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void protectionTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        //1-butyl naphthalene: fused ring atoms 8 and 13, their neighbours 4, 7, 9, 12
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse("CCCCc1cccc2ccccc12");
        tmpFragmenter.prepareMolecule(tmpMolecule);
        Assertions.assertEquals(CRUSHFragmenterTest.bitSet(4, 7, 8, 9, 12, 13), tmpFragmenter.findProtectedAtoms(tmpMolecule));
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.FUSED_RING_ATOMS);
        Assertions.assertEquals(CRUSHFragmenterTest.bitSet(8, 13), tmpFragmenter.findProtectedAtoms(tmpMolecule));
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.NONE);
        Assertions.assertTrue(tmpFragmenter.findProtectedAtoms(tmpMolecule).isEmpty());
        //spiro and bridgehead atoms are in two rings
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.FUSED_RING_ATOMS);
        tmpMolecule = CRUSHFragmenterTest.parse("C1CCC2(CC1)CCCC2");
        tmpFragmenter.prepareMolecule(tmpMolecule);
        Assertions.assertEquals(CRUSHFragmenterTest.bitSet(3), tmpFragmenter.findProtectedAtoms(tmpMolecule));
        //norbornane: the bridgehead atoms and the one-carbon bridge are in both five-membered rings (the six-membered
        //ring is not a relevant cycle, as in the symmetrised SSSR of the RDKit)
        tmpMolecule = CRUSHFragmenterTest.parse("C1CC2CCC1C2");
        tmpFragmenter.prepareMolecule(tmpMolecule);
        Assertions.assertEquals(CRUSHFragmenterTest.bitSet(2, 5, 6), tmpFragmenter.findProtectedAtoms(tmpMolecule));
        //bicyclo[2.2.2]octane: every atom is in two of the three relevant six-membered rings
        tmpMolecule = CRUSHFragmenterTest.parse("C1CC2CCC1CC2");
        tmpFragmenter.prepareMolecule(tmpMolecule);
        Assertions.assertEquals(8, tmpFragmenter.findProtectedAtoms(tmpMolecule).cardinality());
        //conservative mode: every match of the ether rule (6) on the benzylic ether contains the protected naphthyl
        //carbon atom (as atom 1 or 5), so the rule is skipped for the whole molecule, also for the propyl ether
        String tmpEther = "CCCOCCOCc1cccc2ccccc12";
        CRUSHFragmenter tmpConservativeFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpConservativeFragmenter, tmpEther, 6, false).isEmpty());
        Assertions.assertFalse(CRUSHFragmenterTest.stepProducts(tmpConservativeFragmenter, tmpEther, 26, false).isEmpty());
        //without conservative mode, only the matches on protected atoms are dropped: both bonds of the propyl ether
        //oxygen atom are cleaved by rule 6
        CRUSHFragmenter tmpNonConservativeFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        tmpNonConservativeFragmenter.setConservativeModeSetting(false);
        List<List<String>> tmpProducts = CRUSHFragmenterTest.stepProducts(tmpNonConservativeFragmenter, tmpEther, 6, false);
        Assertions.assertEquals(2, tmpProducts.size());
        Assertions.assertTrue(tmpProducts.contains(CRUSHFragmenterTest.canonical(List.of("[6*]CCC",
                "[6*]OCCOCc1cccc2ccccc12"), false)));
        Assertions.assertTrue(tmpProducts.contains(CRUSHFragmenterTest.canonical(List.of("[6*]OCCC",
                "[6*]CCOCc1cccc2ccccc12"), false)));
        //without protection, also both bonds of the benzylic ether oxygen atom are cleaved by rule 6
        tmpNonConservativeFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.NONE);
        Assertions.assertEquals(4, CRUSHFragmenterTest.stepProducts(tmpNonConservativeFragmenter, tmpEther, 6, false).size());
        //propranolol: the ether and aryl ether rules are skipped because of the naphthyl ether, the same bond is cleaved
        //by rule 26
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpConservativeFragmenter, "CC(C)NCC(O)COc1cccc2ccccc12",
                30, false).isEmpty());
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpConservativeFragmenter, "CC(C)NCC(O)COc1cccc2ccccc12",
                26, false).contains(CRUSHFragmenterTest.canonical(List.of("[26*]Oc1cccc2ccccc12", "[26*]CC(O)CNC(C)C"),
                false)));
    }
    //
    /**
     * Tests that ring bonds are never cleaved, although rules 3, 7, 28, and 29 to 33 also match ring bonds: the rings
     * of tetrahydrofurans, piperidines, lipoic acid, and sugars stay intact.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ringBondTest() throws Exception {
        String[] tmpMolecules = {"CCCCC1CCCO1", "CCCCN1CCCCC1", "OC(=O)CCCCC1CCSS1",
                "CCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "c1ccc2c(c1)NCC2", "c1ccc2c(c1)OCO2", "C1CCSC1"};
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        tmpFragmenter.setRingProtectionSetting(CRUSHFragmenter.RingProtectionOption.NONE);
        tmpFragmenter.setMinimumFragmentSizeSetting(1);
        for (String tmpSmiles : tmpMolecules) {
            IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
            tmpFragmenter.prepareMolecule(tmpMolecule);
            BitSet tmpAllAtoms = new BitSet();
            tmpAllAtoms.set(0, tmpMolecule.getAtomCount());
            for (CRUSHFragmenter.CleavageStep tmpStep : tmpFragmenter.findCleavageSteps(tmpMolecule, tmpAllAtoms)) {
                for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                    Assertions.assertFalse(tmpMolecule.getBond(tmpBondIndex).isInRing(), tmpSmiles);
                }
            }
            int tmpRingBonds = 0;
            for (IBond tmpBond : tmpMolecule.bonds()) {
                tmpRingBonds += tmpBond.isInRing() ? 1 : 0;
            }
            int tmpFragmentRingBonds = 0;
            for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse(tmpSmiles))) {
                Cycles.markRingAtomsAndBonds(tmpFragment);
                for (IBond tmpBond : tmpFragment.bonds()) {
                    tmpFragmentRingBonds += tmpBond.isInRing() ? 1 : 0;
                }
            }
            Assertions.assertEquals(tmpRingBonds, tmpFragmentRingBonds, tmpSmiles);
        }
        //the bonds of ring carbon atoms next to the ring heteroatom to substituents are cleaved by rule 7
        this.assertFragments(CRUSHFragmenterTest.createOriginalFragmenter(), "CCCCC1CCCO1", "[7*]C1CCCO1", "[7*]CCCC");
        //lipoic acid: the only decomposition into three fragments (the carboxy group cannot be cleaved off)
        this.assertFragments(CRUSHFragmenterTest.createOriginalFragmenter(), "OC(=O)CCCCC1CCSS1", "[7*]C1CCSS1",
                "[7*]CCC[25*]", "[25*]CC(=O)O");
    }
    //
    /**
     * Tests the cleavage of bonds of aromatic nitrogen atoms: neutral aromatic nitrogen atoms obtain a hydrogen atom
     * (in the script of the authors, these cleavages are dropped because the products cannot be parsed), cationic ones
     * keep their charge.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void aromaticNitrogenTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        this.assertFragments(tmpFragmenter, "CCCCn1ccnc1", "[10*]n1ccnc1", "[10*]CCCC");
        this.assertFragments(tmpFragmenter, "CCCC[n+]1ccccc1", "[10*][n+]1ccccc1", "[10*]CCCC");
        this.assertFragments(tmpFragmenter, "c1ccc(-n2cccc2)cc1", "[22*]c1ccccc1", "[22*]n1cccc1");
        tmpFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCCn1ccnc1", "c1c[nH]cn1", "CCCC");
        this.assertFragments(tmpFragmenter, "CCCC[n+]1ccccc1", "c1cc[nH+]cc1", "CCCC");
        this.assertFragments(tmpFragmenter, "c1ccc(-n2cccc2)cc1", "c1ccccc1", "c1cc[nH]c1");
        //quaternary ammonium: both butyl groups are cleaved off (rules 3 and 9 cleave the same bonds), the nitrogen
        //atom obtains a hydrogen atom per cleaved bond; charged hypervalent sulfur
        this.assertFragments(tmpFragmenter, "CCCC[N+](C)(C)CCCC", "CCCC", "CCCC", "C[NH2+]C");
        Assertions.assertTrue(CRUSHFragmenterTest.stepProducts(tmpFragmenter, "CCCCS(=O)(=O)[O-]", 26, true).contains(
                CRUSHFragmenterTest.canonical(List.of("CCCC", "O=[SH](=O)[O-]"), false)));
    }
    //
    /**
     * Tests the transfer of tetrahedral and double bond stereo elements to the fragments: configurations away from the
     * cleaved bond are kept, a cleaved ligand is replaced by the link atom (or hydrogen atom) at its position, and stereo
     * elements of atoms that are not stereogenic anymore are removed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void stereoTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        //the ether oxygen atom of the quaternary stereocentre is replaced by the link atom at the same position
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCCC[C@](C)(CC)OCCC", 6, false)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCCC[C@](C)(CC)*", "*OCCC"), false)));
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCCC[C@@](C)(CC)OCCC", 6, false)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCCC[C@@](C)(CC)*", "*OCCC"), false)));
        //with hydrogen saturation, the configuration refers to the hydrogen atom at the former position of the oxygen atom
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCCC[C@](C)(CC)OCCC", 6, true)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCCC[C@H](C)CC", "OCCC"), false)));
        //a configuration away from the cleaved bond is kept
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCC[C@H](O)CCOCCC", 6, true)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCC[C@H](O)CCO", "CCC"), false)));
        //a stereocentre with a hydrogen atom loses its configuration with hydrogen saturation, but not with a link atom
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCC[C@@H](CCCC)OCCC", 6, true)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCCCCCCC", "OCCC"), false)));
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCC[C@@H](CCCC)OCCC", 6, false)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCC[C@@H](CCCC)*", "*OCCC"), false)));
        //stereo double bond next to the cleaved bond: the configuration refers to the remaining substituent
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCC/C(=C/CC)OCCC", 26, true)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CC/C=C/CCC", "OCCC"), false)));
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCC/C(=C/CC)OCCC", 26, false)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("CCC/C(=C/CC)*", "*OCCC"), false)));
        //cleaved stereo double bond (step of rule 8), and the route of E-dec-5-ene with three fragments, in which the
        //double bond is kept with its configuration
        Assertions.assertTrue(CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "CCCC/C=C/CCCC", 8, false)
                .contains(CRUSHFragmenterTest.toStereoSmiles(List.of("*=CCCCC", "*=CCCCC"), false)));
        this.assertStereoFragments(tmpFragmenter, "CCCC/C=C/CCCC", "*CCC", "*CCCC", "*C/C=C/*");
        this.assertStereoFragments(tmpFragmenter, "CCCC/C=C\\CCCC", "*CCC", "*CCCC", "*C/C=C\\*");
        //two pseudo atoms at a former stereocentre (one of the input, one link atom): the configuration is removed if the
        //bond of the stereocentre is cleaved, and kept if the other bond of the ether oxygen atom is cleaved
        List<List<String>> tmpStereoProducts = CRUSHFragmenterTest.stereoStepProducts(tmpFragmenter, "*[C@@H](CCC)OCCC",
                6, false);
        Assertions.assertEquals(2, tmpStereoProducts.size());
        Assertions.assertTrue(tmpStereoProducts.contains(CRUSHFragmenterTest.toStereoSmiles(List.of("*C(*)CCC", "*OCCC"),
                false)));
        Assertions.assertTrue(tmpStereoProducts.contains(CRUSHFragmenterTest.toStereoSmiles(List.of("*[C@@H](CCC)O*",
                "*CCC"), false)));
        //the stereocentres of the example molecule CNP0080333.1 are kept
        this.assertStereoFragments(tmpFragmenter, "Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O", "*c1ccc(C)c(O)c1",
                "*[C@H](C)C*", "*CC(O)*", "*C(C)(C)O");
    }
    //
    /**
     * Tests the handling of salts, explicit hydrogen atoms, pseudo atoms of the input, and molecules without cleavage
     * steps, and that the input molecule is not altered.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void inputHandlingTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        //counter ions are kept as fragments
        this.assertFragments(tmpFragmenter, "CCCCCC.[Na+].[Cl-]", "[25*]CCC", "[25*]CCC", "[Na+]", "[Cl-]");
        this.assertFragments(tmpFragmenter, "Cc1ccccc1.[Na+]", "Cc1ccccc1", "[Na+]");
        //pseudo atoms of the input do not count as heavy atoms
        this.assertFragments(tmpFragmenter, "*C1CCCO1", "*C1CCCO1");
        this.assertFragments(tmpFragmenter, "*CCCOCCC", "*CCCO[6*]", "[6*]CCC");
        //molecules without cleavage steps
        this.assertFragments(tmpFragmenter, "c1ccccc1", "c1ccccc1");
        this.assertFragments(tmpFragmenter, "CCOCC", "CCOCC");
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("c1ccccc1")).isEmpty());
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("Cc1ccccc1.[Na+]")).isEmpty());
        Assertions.assertEquals(2, tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("CCCCCC")).size());
        tmpFragmenter.setFragmentOutputSetting(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(CRUSHFragmenterTest.parse("c1ccccc1")).isEmpty());
        CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(CRUSHFragmenterTest.parse("c1ccccc1"));
        Assertions.assertTrue(tmpRoot.isLeaf());
        Assertions.assertTrue(tmpRoot.getLeaves().isEmpty());
        Assertions.assertTrue(tmpRoot.isExhaustive());
        //explicit hydrogen atoms do not change the result
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse("Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O");
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        this.assertFragments(CRUSHFragmenterTest.createOriginalFragmenter(), tmpMolecule, "[14*]c1ccc(C)c(O)c1",
                "[14*]C(C)C[25*]", "[25*]CC(O)[25*]", "[25*]C(C)(C)O");
        //the input molecule is not altered
        IAtomContainer tmpInput = CRUSHFragmenterTest.parse("CCCCCC");
        CRUSHFragmenterTest.createOriginalFragmenter().fragmentMolecule(tmpInput);
        CRUSHFragmenterTest.createOriginalFragmenter().buildRoute(tmpInput);
        Assertions.assertEquals(6, tmpInput.getAtomCount());
        Assertions.assertEquals(5, tmpInput.getBondCount());
        //molecules given with aromatic bonds without bond orders are kekulized
        IAtomContainer tmpAromatic = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CCCOc1ccccc1");
        for (IBond tmpBond : tmpAromatic.bonds()) {
            if (tmpBond.isAromatic()) {
                tmpBond.setOrder(IBond.Order.UNSET);
            }
        }
        Assertions.assertEquals(2, CRUSHFragmenterTest.createOriginalFragmenter().fragmentMolecule(tmpAromatic).size());
        //a molecule whose aromatic ring cannot be kekulized is not fragmented as a non-aromatic molecule but rejected
        SmilesParser tmpNonKekulizingParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpNonKekulizingParser.kekulise(false);
        IAtomContainer tmpNonKekulizable = tmpNonKekulizingParser.parseSmiles("c1cccc1CCCCCOCCCCCC");
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> CRUSHFragmenterTest.createOriginalFragmenter().fragmentMolecule(tmpNonKekulizable));
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> CRUSHFragmenterTest.createOriginalFragmenter().buildRoute(tmpNonKekulizable));
    }
    //
    /**
     * Tests the filter for molecules with bonds of unset order that cannot be resolved by kekulization, e.g. the
     * aromatic remnants of opened rings in fragments of the Ertl fragmenter, which MORTAR parses without kekulization
     * because their kekulization fails: they are filtered and not altered, whereas molecules with unset bond orders that
     * can be kekulized are not filtered and give the same fragments as the kekulized molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void unresolvableBondOrderTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        //unset bond orders that can be resolved, also in a molecule with an amide and a biaryl bond
        for (String tmpSmiles : new String[] {"CCCc1ccccc1OCCCC", "CCCCC(=O)Nc1ccc(cc1)-c1ccccc1"}) {
            IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpSmiles, false, false);
            Assertions.assertTrue(this.hasUnsetBondOrder(tmpMolecule), tmpSmiles);
            Assertions.assertFalse(tmpFragmenter.shouldBeFiltered(tmpMolecule), tmpSmiles);
            Assertions.assertTrue(tmpFragmenter.canBeFragmented(tmpMolecule), tmpSmiles);
            List<String> tmpExpectedKeys = this.createUniqueSmilesList(tmpFragmenter.fragmentMolecule(
                    ChemUtil.parseSmilesToAtomContainer(tmpSmiles, true, true)));
            List<String> tmpKeys = this.createUniqueSmilesList(tmpFragmenter.fragmentMolecule(tmpMolecule));
            Assertions.assertTrue(tmpExpectedKeys.size() > 1, tmpSmiles);
            Assertions.assertFalse(tmpKeys.contains(null), tmpSmiles);
            Assertions.assertEquals(tmpExpectedKeys, tmpKeys, tmpSmiles);
        }
        //aromatic remnants of opened rings cannot be kekulized
        for (String tmpSmiles : new String[] {"ccc", "c1ccc(cc1)ccc-c2ccccc2"}) {
            IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpSmiles, false, false);
            List<String> tmpBondOrdersAndAromaticity = this.getBondOrdersAndAromaticity(tmpMolecule);
            Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpMolecule), tmpSmiles);
            Assertions.assertFalse(tmpFragmenter.canBeFragmented(tmpMolecule), tmpSmiles);
            Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.fragmentMolecule(tmpMolecule),
                    tmpSmiles);
            //the input molecule is not altered
            Assertions.assertTrue(this.hasUnsetBondOrder(tmpMolecule), tmpSmiles);
            Assertions.assertEquals(tmpBondOrdersAndAromaticity, this.getBondOrdersAndAromaticity(tmpMolecule), tmpSmiles);
        }
        //a bond of unset or no order between atoms that are not aromatic cannot be resolved
        IAtomContainer tmpUnsetBondMolecule = SilentChemObjectBuilder.getInstance().newAtomContainer();
        for (int i = 0; i < 2; i++) {
            IAtom tmpAtom = SilentChemObjectBuilder.getInstance().newInstance(IAtom.class, "C");
            tmpAtom.setImplicitHydrogenCount(3);
            tmpUnsetBondMolecule.addAtom(tmpAtom);
        }
        tmpUnsetBondMolecule.addBond(0, 1, IBond.Order.UNSET);
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpUnsetBondMolecule));
        Assertions.assertFalse(tmpFragmenter.canBeFragmented(tmpUnsetBondMolecule));
        tmpUnsetBondMolecule.getBond(0).setOrder(null);
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(tmpUnsetBondMolecule));
        Assertions.assertNull(tmpUnsetBondMolecule.getBond(0).getOrder());
        //molecules without unset bond orders are not filtered
        Assertions.assertFalse(tmpFragmenter.shouldBeFiltered(
                ChemUtil.parseSmilesToAtomContainer("CCCc1ccccc1OCCCC", true, true)));
    }
    //
    /**
     * Tests that the fragments do not depend on the order of the atoms and bonds of the input: the molecules are
     * fragmented with their atoms and bonds in random order, and the fragments must be identical, including rule numbers
     * and stereochemistry.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void atomOrderInvarianceTest() throws Exception {
        String[] tmpMolecules = {"Cc1ccc([C@H](C)CCC(O)C(C)(C)O)cc1O", "CC(C)NCC(O)COc1cccc2ccccc12",
                "CCCCCCCCCCCCCCCc1cccc(O)c1C(=O)O", "OCCOCCOCCOCCO",
                "COC1=CC=CC2=C1C(=O)C1=C(O)C3=C(C[C@@](O)(C(=O)CO)C[C@@H]3O[C@H]3C[C@H](N)C[C@H](C)O3)C(O)=C1C2=O",
                "NC1=NC(N)=C2N=C(CNC3=CC=C(C(=O)N[C@@H](CCCNC(=O)C4=CC=CC=C4C(=O)O)C(=O)O)C=C3)C=NC2=N1",
                "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(COC(N)=O)CS[C@H]12)C1=CC=CO1.[Na+]"};
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        Random tmpRandom = new Random(42);
        for (String tmpSmiles : tmpMolecules) {
            IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
            List<String> tmpExpected = CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false);
            List<String> tmpExpectedStereo = CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false);
            for (int i = 0; i < 5; i++) {
                List<IAtomContainer> tmpShuffledFragments = tmpFragmenter.fragmentMolecule(
                        CRUSHFragmenterTest.shuffle(tmpMolecule, tmpRandom));
                Assertions.assertEquals(tmpExpected, CRUSHFragmenterTest.toCanonicalSmiles(tmpShuffledFragments, false),
                        tmpSmiles);
                //the unique SMILES codes of the CDK are not canonical for some fragments with link atoms bonded by
                //different bond orders to the same atom (e.g. *=C(*)C(*)=O), so these are compared by isomorphism
                if (!tmpExpectedStereo.equals(CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpShuffledFragments, false))) {
                    //in both directions, since the stereo elements of the query must be matched
                    Assertions.assertTrue(CRUSHFragmenterTest.isIsomorphicWithStereo(tmpFragments, tmpShuffledFragments)
                            && CRUSHFragmenterTest.isIsomorphicWithStereo(tmpShuffledFragments, tmpFragments), tmpSmiles);
                }
            }
        }
    }
    //
    /**
     * Tests the greedy decomposition after the search budget is exhausted, with the largest molecule of the BIOFACQUIM
     * data set (CNP0085395.1, 189 heavy atoms): the route is not exhaustive, but a complete decomposition, i.e. every
     * heavy atom is in exactly one final fragment, and no final fragment has a cleavage step.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void searchBudgetTest() throws Exception {
        String tmpSmiles = "CCCCCCCCCCCC(=O)O[C@H]1[C@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@H](O[C@H]3[C@H](O[C@@H](CCCCC)"
                + "CCCCCCCCCC(=O)O[C@H]4[C@H](O)[C@@H](O)[C@H](O[C@@H]5[C@@H](O[C@@H]6O[C@@H](C)[C@H](OC(=O)[C@@H](C)CC)"
                + "[C@@H](O)[C@H]6OC(=O)/C=C/c6ccccc6)[C@H](C)O[C@@H](O[C@@H]6[C@@H](O)[C@H]7OC(=O)CCCCCCCCC[C@H](CCCCC)O"
                + "[C@@H]8O[C@H](C)[C@H](O)[C@H](O)[C@H]8O[C@@H]7O[C@H]6C)[C@@H]5OC(=O)CCCCCCCCCCC)O[C@@H]4CO)O[C@H](C)"
                + "[C@H](O)[C@@H]3O)O[C@H]2C)O[C@@H](C)[C@H](O[C@@H]2O[C@@H](C)[C@H](OC(=O)[C@@H](C)CC)[C@@H](OC(=O)/C=C/"
                + "c3ccccc3)[C@H]2O)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O";
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
        CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(tmpMolecule);
        Assertions.assertFalse(tmpRoot.isExhaustive());
        List<CRUSHFragmenter.CRUSHRouteNode> tmpLeaves = tmpRoot.getLeaves();
        BitSet tmpCoveredAtoms = new BitSet();
        int tmpAtomCount = 0;
        IAtomContainer tmpPreparedMolecule = tmpRoot.getFragment();
        for (CRUSHFragmenter.CRUSHRouteNode tmpLeaf : tmpLeaves) {
            BitSet tmpAtomSet = tmpLeaf.getAtomSet();
            Assertions.assertFalse(tmpCoveredAtoms.intersects(tmpAtomSet));
            tmpCoveredAtoms.or(tmpAtomSet);
            tmpAtomCount += tmpAtomSet.cardinality();
            Assertions.assertTrue(tmpFragmenter.findCleavageSteps(tmpPreparedMolecule, tmpAtomSet).isEmpty());
        }
        Assertions.assertEquals(tmpPreparedMolecule.getAtomCount(), tmpAtomCount);
        Assertions.assertTrue(tmpLeaves.size() > 10);
        Assertions.assertEquals(tmpLeaves.size(), tmpFragmenter.fragmentMolecule(tmpMolecule).size());
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
    private void assertFragments(CRUSHFragmenter aFragmenter, String aSmiles, String... anExpectedFragments) throws Exception {
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
     * Returns the aromaticity flags of the atoms and the orders and aromaticity flags of the bonds of the given
     * molecule, in the order of the atoms and bonds.
     *
     * @param aMolecule molecule
     * @return list with one entry per atom and per bond
     */
    private List<String> getBondOrdersAndAromaticity(IAtomContainer aMolecule) {
        List<String> tmpList = new ArrayList<>(aMolecule.getAtomCount() + aMolecule.getBondCount());
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpList.add(String.valueOf(tmpAtom.isAromatic()));
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            tmpList.add(tmpBond.getOrder() + " " + tmpBond.isAromatic());
        }
        return tmpList;
    }
    //
    /**
     * Returns the sorted unique SMILES codes with stereochemistry of the given fragments as MORTAR creates them
     * ({@link ChemUtil#createUniqueSmiles(IAtomContainer, boolean)}), with null for fragments whose SMILES code cannot be
     * created.
     *
     * @param aFragmentList fragments
     * @return sorted list of unique SMILES codes, null values first
     */
    private List<String> createUniqueSmilesList(List<IAtomContainer> aFragmentList) {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(ChemUtil.createUniqueSmiles(tmpFragment, true));
        }
        tmpSmilesList.sort(Comparator.nullsFirst(Comparator.naturalOrder()));
        return tmpSmilesList;
    }
    //
    /**
     * Returns the products of all steps of the given rule for the given molecule (one connected component), as sorted
     * lists of canonical SMILES codes with rule numbers.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param aRuleNumber rule number
     * @param isHydrogenSaturation true for products saturated with hydrogen atoms, false for link atoms
     * @return list of the products of the steps
     * @throws Exception if anything goes wrong
     */
    static List<List<String>> stepProducts(CRUSHFragmenter aFragmenter, String aSmiles, int aRuleNumber,
                                           boolean isHydrogenSaturation) throws Exception {
        List<List<String>> tmpResult = new ArrayList<>();
        for (List<IAtomContainer> tmpProducts : CRUSHFragmenterTest.stepFragments(aFragmenter, aSmiles, aRuleNumber,
                isHydrogenSaturation)) {
            tmpResult.add(CRUSHFragmenterTest.toCanonicalSmiles(tmpProducts, false));
        }
        return tmpResult;
    }
    //
    /**
     * Returns the products of all steps of the given rule for the given molecule (one connected component), as sorted
     * lists of unique SMILES codes with stereochemistry (link atoms as '*').
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param aRuleNumber rule number
     * @param isHydrogenSaturation true for products saturated with hydrogen atoms, false for link atoms
     * @return list of the products of the steps
     * @throws Exception if anything goes wrong
     */
    static List<List<String>> stereoStepProducts(CRUSHFragmenter aFragmenter, String aSmiles, int aRuleNumber,
                                                 boolean isHydrogenSaturation) throws Exception {
        List<List<String>> tmpResult = new ArrayList<>();
        for (List<IAtomContainer> tmpProducts : CRUSHFragmenterTest.stepFragments(aFragmenter, aSmiles, aRuleNumber,
                isHydrogenSaturation)) {
            tmpResult.add(CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpProducts, false));
        }
        return tmpResult;
    }
    //
    /**
     * Returns the products of all steps of the given rule for the given molecule (one connected component).
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param aRuleNumber rule number
     * @param isHydrogenSaturation true for products saturated with hydrogen atoms, false for link atoms
     * @return list of the products of the steps
     * @throws Exception if anything goes wrong
     */
    private static List<List<IAtomContainer>> stepFragments(CRUSHFragmenter aFragmenter, String aSmiles, int aRuleNumber,
                                                            boolean isHydrogenSaturation) throws Exception {
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(aSmiles);
        aFragmenter.prepareMolecule(tmpMolecule);
        BitSet tmpAllAtoms = new BitSet();
        tmpAllAtoms.set(0, tmpMolecule.getAtomCount());
        List<List<IAtomContainer>> tmpResult = new ArrayList<>();
        for (CRUSHFragmenter.CleavageStep tmpStep : aFragmenter.findCleavageSteps(tmpMolecule, tmpAllAtoms)) {
            if (tmpStep.ruleNumber() != aRuleNumber) {
                continue;
            }
            int[] tmpLabels = new int[tmpMolecule.getBondCount()];
            for (int tmpBondIndex : tmpStep.cleavedBonds()) {
                tmpLabels[tmpBondIndex] = aRuleNumber;
            }
            List<IAtomContainer> tmpProducts = new ArrayList<>();
            for (BitSet tmpPart : tmpStep.parts()) {
                tmpProducts.add(CRUSHFragmenter.createFragment(tmpMolecule, tmpPart, tmpLabels, isHydrogenSaturation));
            }
            tmpResult.add(tmpProducts);
        }
        return tmpResult;
    }
    //
    /**
     * Returns the bonds between mapped atoms of the given SMARTS as pairs of map numbers.
     *
     * @param aSmarts SMARTS
     * @return set of sorted pairs of map numbers
     */
    private static Set<List<Integer>> mappedBonds(String aSmarts) {
        QueryAtomContainer tmpQuery = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
        Assertions.assertTrue(Smarts.parse(tmpQuery, aSmarts), aSmarts);
        Set<List<Integer>> tmpBonds = new HashSet<>();
        for (IBond tmpBond : tmpQuery.bonds()) {
            int tmpBegin = tmpBond.getBegin().getMapIdx();
            int tmpEnd = tmpBond.getEnd().getMapIdx();
            if (tmpBegin > 0 && tmpEnd > 0) {
                tmpBonds.add(List.of(Math.min(tmpBegin, tmpEnd), Math.max(tmpBegin, tmpEnd)));
            }
        }
        return tmpBonds;
    }
    //
    /**
     * Returns a bit set with the given indices.
     *
     * @param anIndexArray indices
     * @return bit set
     */
    private static BitSet bitSet(int... anIndexArray) {
        BitSet tmpBitSet = new BitSet();
        for (int tmpIndex : anIndexArray) {
            tmpBitSet.set(tmpIndex);
        }
        return tmpBitSet;
    }
    //
    /**
     * Returns a copy of the given molecule with its atoms and bonds in random order; the stereo elements are mapped to
     * the copy.
     *
     * @param aMolecule molecule
     * @param aRandom random number generator
     * @return copy with shuffled atoms and bonds
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    static IAtomContainer shuffle(IAtomContainer aMolecule, Random aRandom) throws CloneNotSupportedException {
        IAtomContainer tmpClone = aMolecule.clone();
        List<IAtom> tmpAtoms = new ArrayList<>();
        for (IAtom tmpAtom : tmpClone.atoms()) {
            tmpAtoms.add(tmpAtom);
        }
        List<IBond> tmpBonds = new ArrayList<>();
        for (IBond tmpBond : tmpClone.bonds()) {
            tmpBonds.add(tmpBond);
        }
        Collections.shuffle(tmpAtoms, aRandom);
        Collections.shuffle(tmpBonds, aRandom);
        IAtomContainer tmpShuffled = SilentChemObjectBuilder.getInstance().newAtomContainer();
        Map<IAtom, IAtom> tmpAtomMap = new HashMap<>();
        for (IAtom tmpAtom : tmpAtoms) {
            IAtom tmpNewAtom = tmpAtom.clone();
            tmpShuffled.addAtom(tmpNewAtom);
            tmpAtomMap.put(tmpAtom, tmpShuffled.getAtom(tmpShuffled.getAtomCount() - 1));
        }
        Map<IBond, IBond> tmpBondMap = new HashMap<>();
        for (IBond tmpBond : tmpBonds) {
            //the direction of the bond is reversed randomly
            boolean tmpIsReversed = aRandom.nextBoolean();
            IAtom tmpBegin = tmpAtomMap.get(tmpIsReversed ? tmpBond.getEnd() : tmpBond.getBegin());
            IAtom tmpEnd = tmpAtomMap.get(tmpIsReversed ? tmpBond.getBegin() : tmpBond.getEnd());
            IBond tmpNewBond = SilentChemObjectBuilder.getInstance().newInstance(IBond.class, tmpBegin, tmpEnd,
                    tmpBond.getOrder());
            tmpNewBond.setIsAromatic(tmpBond.isAromatic());
            tmpShuffled.addBond(tmpNewBond);
            tmpBondMap.put(tmpBond, tmpShuffled.getBond(tmpShuffled.getBondCount() - 1));
        }
        List<IStereoElement> tmpElements = new ArrayList<>();
        for (IStereoElement tmpElement : tmpClone.stereoElements()) {
            tmpElements.add(tmpElement.map(tmpAtomMap, tmpBondMap));
        }
        tmpShuffled.setStereoElements(tmpElements);
        return tmpShuffled;
    }
    //
    /**
     * Returns whether the given lists of fragments are isomorphic as multisets, including stereochemistry: elements,
     * formal charges, implicit hydrogen counts, mass numbers (rule numbers of link atoms), bond orders (or aromaticity),
     * and the configurations of the stereo elements of the first fragments must match.
     *
     * @param aFragmentList1 first list of fragments
     * @param aFragmentList2 second list of fragments
     * @return true if there is a one-to-one assignment of isomorphic fragments
     * @throws CloneNotSupportedException if a fragment cannot be cloned
     */
    static boolean isIsomorphicWithStereo(List<IAtomContainer> aFragmentList1, List<IAtomContainer> aFragmentList2)
            throws CloneNotSupportedException {
        if (aFragmentList1.size() != aFragmentList2.size()) {
            return false;
        }
        AtomMatcher tmpAtomMatcher = new AtomMatcher() {
            @Override
            public boolean matches(IAtom anAtom1, IAtom anAtom2) {
                return Objects.equals(anAtom1.getAtomicNumber(), anAtom2.getAtomicNumber())
                        && Objects.equals(anAtom1.getFormalCharge(), anAtom2.getFormalCharge())
                        && Objects.equals(anAtom1.getImplicitHydrogenCount(), anAtom2.getImplicitHydrogenCount())
                        && Objects.equals(anAtom1.getMassNumber(), anAtom2.getMassNumber());
            }
        };
        BondMatcher tmpBondMatcher = new BondMatcher() {
            @Override
            public boolean matches(IBond aBond1, IBond aBond2) {
                return aBond1.isAromatic() == aBond2.isAromatic()
                        && (aBond1.isAromatic() || aBond1.getOrder() == aBond2.getOrder());
            }
        };
        List<IAtomContainer> tmpRemaining = new ArrayList<>(aFragmentList2.size());
        for (IAtomContainer tmpFragment : aFragmentList2) {
            IAtomContainer tmpClone = tmpFragment.clone();
            SmartsPattern.prepare(tmpClone);
            tmpRemaining.add(tmpClone);
        }
        for (IAtomContainer tmpFragment : aFragmentList1) {
            IAtomContainer tmpQuery = tmpFragment.clone();
            SmartsPattern.prepare(tmpQuery);
            Pattern tmpPattern = VentoFoggia.findIdentical(tmpQuery, tmpAtomMatcher, tmpBondMatcher);
            IAtomContainer tmpMatch = null;
            for (IAtomContainer tmpCandidate : tmpRemaining) {
                if (tmpCandidate.getAtomCount() == tmpQuery.getAtomCount()
                        && tmpPattern.matchAll(tmpCandidate).stereochemistry().atLeast(1)) {
                    tmpMatch = tmpCandidate;
                    break;
                }
            }
            if (tmpMatch == null) {
                return false;
            }
            tmpRemaining.remove(tmpMatch);
        }
        return true;
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
     * Returns the fragments of the given route nodes.
     *
     * @param aNodeList route nodes
     * @return list of their fragments
     */
    static List<IAtomContainer> fragmentsOf(List<CRUSHFragmenter.CRUSHRouteNode> aNodeList) {
        List<IAtomContainer> tmpFragments = new ArrayList<>(aNodeList.size());
        for (CRUSHFragmenter.CRUSHRouteNode tmpNode : aNodeList) {
            tmpFragments.add(tmpNode.getFragment());
        }
        return tmpFragments;
    }
    //
    /**
     * Returns the nodes of the given route whose fragments are returned with the fragment output option
     * {@link CRUSHFragmenter.FragmentOutputOption#ALL_ROUTE_FRAGMENTS}: all nodes except the root and, for a molecule
     * with several connected components (root with rule number 0), except the components that are cleaved.
     *
     * @param aRoot root of the route
     * @return list of the nodes in breadth-first order
     */
    static List<CRUSHFragmenter.CRUSHRouteNode> allRouteFragmentNodes(CRUSHFragmenter.CRUSHRouteNode aRoot) {
        List<CRUSHFragmenter.CRUSHRouteNode> tmpNodes = new ArrayList<>();
        for (CRUSHFragmenter.CRUSHRouteNode tmpNode : aRoot.getAllDescendants()) {
            if (aRoot.getRuleNumber() == 0 && tmpNode.getParent() == aRoot && !tmpNode.isLeaf()) {
                continue;
            }
            tmpNodes.add(tmpNode);
        }
        return tmpNodes;
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
            tmpSmilesList.add(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragment, isStereo));
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
            return CRUSHFragmenterTest.toStereoSmiles(aFragment, false);
        }
        String tmpSmiles = SmilesGenerator.create(aFragment, SmiFlavor.Unique | SmiFlavor.AtomicMass,
                new int[aFragment.getAtomCount()]);
        return CRUSHFragmenterTest.canonical(tmpSmiles, false);
    }
    //
    /**
     * Returns the unique SMILES code with stereochemistry of the given fragment as MORTAR creates it, i.e. with all link
     * atoms written as '*' (the rule numbers are not encoded). The CDK creates canonical SMILES codes with
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
     * Returns the unique SMILES code with stereochemistry of the given fragment (see
     * {@link #toStereoSmiles(IAtomContainer, boolean)}), without the configurations of atoms that are not stereogenic
     * or only potential (para) stereocentres according to the CDK, i.e. atoms whose stereogenicity depends on the
     * configuration of other atoms.
     *
     * @param aFragment fragment, it is not altered
     * @return unique SMILES code with stereochemistry of the true stereocentres
     * @throws CDKException if the fragment cannot be cloned
     */
    @SuppressWarnings("rawtypes")
    static String toStereoSmilesWithoutParaStereocentres(IAtomContainer aFragment) throws CDKException {
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
        SmartsPattern.prepare(tmpFragment);
        if (tmpFragment.stereoElements().iterator().hasNext()) {
            Stereocenters tmpStereocenters = Stereocenters.of(tmpFragment);
            List<IStereoElement> tmpElements = new ArrayList<>();
            for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                List<IAtom> tmpAtoms = new ArrayList<>(2);
                if (tmpElement.getConfigClass() == IStereoElement.TH) {
                    tmpAtoms.add((IAtom) tmpElement.getFocus());
                } else if (tmpElement.getConfigClass() == IStereoElement.CT) {
                    tmpAtoms.add(((IBond) tmpElement.getFocus()).getBegin());
                    tmpAtoms.add(((IBond) tmpElement.getFocus()).getEnd());
                }
                boolean tmpIsTrueStereocentre = true;
                for (IAtom tmpAtom : tmpAtoms) {
                    if (tmpStereocenters.stereocenterType(tmpFragment.indexOf(tmpAtom)) != Stereocenters.Stereocenter.True) {
                        tmpIsTrueStereocentre = false;
                    }
                }
                if (tmpIsTrueStereocentre) {
                    tmpElements.add(tmpElement);
                }
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
            tmpSmilesList.add(CRUSHFragmenterTest.toStereoSmiles(CRUSHFragmenterTest.parse(tmpSmiles),
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
            tmpSmilesList.add(CRUSHFragmenterTest.toStereoSmiles(tmpFragment, isNonStereogenicElementRemoved));
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
            tmpSmilesList.add(CRUSHFragmenterTest.canonical(tmpSmiles, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns a canonical SMILES code of the given SMILES code that is independent of the aromaticity notation of the
     * input: the SMILES is parsed and kekulized, aromaticity is perceived with the Daylight model (as by the
     * fragmenter), and a unique SMILES with mass numbers (i.e. rule numbers) and aromatic symbols is created.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether stereochemistry is encoded (as MORTAR encodes it, without rule numbers)
     * @return canonical SMILES code
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static String canonical(String aSmiles, boolean isStereo) throws CDKException {
        IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(aSmiles);
        SmartsPattern.prepare(tmpMolecule);
        if (isStereo) {
            return CRUSHFragmenterTest.toStereoSmiles(tmpMolecule, false);
        }
        return SmilesGenerator.create(tmpMolecule, SmiFlavor.Unique | SmiFlavor.AtomicMass | SmiFlavor.UseAromaticSymbols,
                new int[tmpMolecule.getAtomCount()]);
    }
}
