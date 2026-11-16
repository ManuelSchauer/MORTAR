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
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.smarts.SmartsPattern;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

/**
 * Cross-validates the {@link CRUSHFragmenter} against reference values that are computed with the RDKit from the rule
 * file of the authors by the script CRUSH_Reference_Generator.py in the test resources:
 * <ul>
 *     <li>steps (file CRUSH_Reference_Steps_BIOFACQUIM.tsv): for every molecule of the BIOFACQUIM data set and up to 20
 *     intermediate fragments of its route search, the cleavage steps of every rule must be identical to those of the
 *     independent RDKit implementation of the same definition (cleaved bonds and hydrogen-saturated products, including
 *     stereochemistry), and every valid step of the published procedure (RunReactants semantics) must be among them
 *     with the same products;</li>
 *     <li>routes (file CRUSH_Reference_Fragments_BIOFACQUIM.tsv): for every molecule whose search is exhaustive in both
 *     implementations, the number of final fragments must be the maximum found by the exact search of the RDKit
 *     implementation, and the final fragments (hydrogen saturation) must be one of the optimal fragment sets.</li>
 *     <li>single pass (file CRUSH_Reference_SinglePass_BIOFACQUIM.tsv): for every molecule, the fragments of the
 *     single-pass variant with default settings (with link atoms carrying the rule numbers) must be one of the results
 *     that the toolkit-independent part of its processing order allows, as determined by an independent RDKit
 *     implementation.</li>
 * </ul>
 * The steps, routes, and rule cases refer to the original variant with the minimum fragment size 3 of the
 * publication. For the comparison with stereochemistry, stereo elements of atoms that are not stereogenic according to
 * the CDK are removed on both sides, because the RDKit keeps some of them.
 *
 * <p>The BIOFACQUIM data set (natural products from Mexico, CC BY 4.0, see Tutorial/BIOFACQUIM-description-licence.txt)
 * is part of the test resources. A further steps reference file, e.g. for a sample of ChEMBL, can be given with the
 * environment variable {@value #REFERENCE_FILE_ENVIRONMENT_VARIABLE}, a further routes reference file with
 * {@value #ROUTE_REFERENCE_FILE_ENVIRONMENT_VARIABLE}.</p>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class CRUSHFragmenterReferenceTest {
    /**
     * Name of the steps reference file of the BIOFACQUIM data set in the test resources.
     */
    private static final String BIOFACQUIM_STEPS_FILE = "CRUSH_Reference_Steps_BIOFACQUIM.tsv";
    /**
     * Name of the routes reference file of the BIOFACQUIM data set in the test resources.
     */
    private static final String BIOFACQUIM_ROUTES_FILE = "CRUSH_Reference_Fragments_BIOFACQUIM.tsv";
    /**
     * Name of the single-pass reference file of the BIOFACQUIM data set in the test resources.
     */
    private static final String BIOFACQUIM_SINGLE_PASS_FILE = "CRUSH_Reference_SinglePass_BIOFACQUIM.tsv";
    /**
     * Environment variable with the path of an additional steps reference file.
     */
    private static final String REFERENCE_FILE_ENVIRONMENT_VARIABLE = "MORTAR_CRUSH_REFERENCE_FILE";
    /**
     * Environment variable with the path of an additional routes reference file.
     */
    private static final String ROUTE_REFERENCE_FILE_ENVIRONMENT_VARIABLE = "MORTAR_CRUSH_ROUTE_REFERENCE_FILE";
    /**
     * Maximum number of reported mismatches.
     */
    private static final int MAXIMUM_REPORTED_MISMATCHES = 25;
    //
    /**
     * Result of the comparison of a steps reference file.
     *
     * @param molecules number of compared molecules (including intermediate fragments)
     * @param skipped number of molecules whose SMILES could not be parsed or prepared by the CDK
     * @param mismatches number of molecules with different steps (without stereochemistry)
     * @param stereoMismatches number of molecules with steps whose products differ in stereochemistry only
     * @param originalMismatches number of molecules with a valid step of the published procedure that is not a step of
     *                           the fragmenter with the same products (without stereochemistry)
     * @param originalStereoDifferences number of valid steps of the published procedure whose products differ from those
     *                                  of the fragmenter in stereochemistry only (reported, not asserted)
     * @param protectedAtomMismatches number of molecules with different protected atoms
     * @param expectationMismatches number of rule cases whose expectation is not met
     * @param paraStereocentreDifferences number of molecules with steps whose products differ only in the configuration
     *                                    of potential (para) stereocentres (reported, not asserted)
     * @param isomorphismMatches number of molecules with steps whose products have different unique SMILES codes in the
     *                           CDK but are isomorphic (non-canonical unique SMILES codes of radicals)
     * @param stepMilliseconds time spent in the determination of the steps
     * @param report description of the first mismatches
     */
    record StepComparisonResult(int molecules, int skipped, int mismatches, int stereoMismatches, int originalMismatches,
                                int originalStereoDifferences, int protectedAtomMismatches, int expectationMismatches,
                                int paraStereocentreDifferences, int isomorphismMatches, long stepMilliseconds,
                                String report) {
        /**
         * Returns a summary of the numbers.
         *
         * @return summary
         */
        String summary() {
            return String.format("%d molecules, %d skipped, %d mismatches, %d stereo mismatches, %d mismatches against "
                            + "the published procedure, %d stereo differences against the published procedure, %d "
                            + "protected atom mismatches, %d expectation mismatches, %d molecules with para "
                            + "stereocentre differences, %d molecules with products matched by isomorphism, %d ms",
                    this.molecules, this.skipped, this.mismatches, this.stereoMismatches, this.originalMismatches,
                    this.originalStereoDifferences, this.protectedAtomMismatches, this.expectationMismatches,
                    this.paraStereocentreDifferences, this.isomorphismMatches, this.stepMilliseconds);
        }
    }
    //
    /**
     * Result of the comparison of a routes reference file.
     *
     * @param molecules number of compared molecules
     * @param bothExhaustive number of molecules whose search is exhaustive in both implementations
     * @param referenceExhaustive number of molecules whose reference search is exhaustive
     * @param fragmenterExhaustive number of molecules whose search is exhaustive in the fragmenter
     * @param mismatches number of molecules with a different number of final fragments or final fragments that are not
     *                   an optimal set of the reference (without stereochemistry)
     * @param stereoMismatches number of molecules whose final fragments differ from the optimal sets only in
     *                         stereochemistry
     * @param uncheckedSets number of molecules with more optimal sets than listed, for which only the number of final
     *                      fragments is compared
     * @param isomorphismMatches number of molecules whose final fragments have different unique SMILES codes in the CDK
     *                           than the optimal set they are isomorphic to (non-canonical unique SMILES codes, e.g. of
     *                           Kekulé structures of sulfoximine and phosphazene rings)
     * @param milliseconds time spent in the fragmentation of all molecules
     * @param slowestMilliseconds time of the slowest molecule
     * @param slowestMolecule name of the slowest molecule
     * @param report description of the first mismatches
     */
    record RouteComparisonResult(int molecules, int bothExhaustive, int referenceExhaustive, int fragmenterExhaustive,
                                 int mismatches, int stereoMismatches, int uncheckedSets, int isomorphismMatches,
                                 long milliseconds, long slowestMilliseconds, String slowestMolecule, String report) {
    }
    //
    /**
     * Result of the comparison of a single-pass reference file.
     *
     * @param molecules number of compared molecules
     * @param severalResults number of molecules for which the reference lists more than one possible result, i.e. whose
     *                       result depends on the canonical rank of the bonds (fourth sort key of the processing order)
     * @param uncheckedMolecules number of molecules with more possible results than listed, which are not compared
     * @param mismatches number of molecules whose fragments are not one of the possible results (without
     *                   stereochemistry, with the rule numbers of the link atoms)
     * @param stereoMismatches number of molecules whose fragments differ from the possible results only in
     *                         stereochemistry
     * @param paraStereocentreDifferences number of molecules whose fragments differ from a possible result only in the
     *                                    configuration of potential (para) stereocentres (reported, not asserted)
     * @param isomorphismMatches number of molecules whose fragments have different unique SMILES codes in the CDK than
     *                           the possible result they are isomorphic to
     * @param milliseconds time spent in the fragmentation of all molecules
     * @param report description of the first mismatches
     */
    record SinglePassComparisonResult(int molecules, int severalResults, int uncheckedMolecules, int mismatches,
                                      int stereoMismatches, int paraStereocentreDifferences, int isomorphismMatches,
                                      long milliseconds, String report) {
    }
    //
    /**
     * A molecule of a steps reference file with its reference rows.
     */
    private static final class StepBlock {
        /**
         * Name of the molecule.
         */
        private final String name;
        /**
         * SMILES code of the molecule.
         */
        private final String smiles;
        /**
         * Tested rule of a rule case, 0 otherwise.
         */
        private final int targetRule;
        /**
         * Expectation of a rule case ("positive", "negative") or another text.
         */
        private final String expectation;
        /**
         * Protected atom indices, separated by spaces.
         */
        private final String protectedAtoms;
        /**
         * Rows of the steps of the CRUSHFragmenter definition (part 2).
         */
        private final List<String[]> mortarRows;
        /**
         * Rows of the steps of the published procedure (part 1).
         */
        private final List<String[]> originalRows;
        /**
         * Constructor.
         *
         * @param aColumnArray columns of the molecule row
         */
        private StepBlock(String[] aColumnArray) {
            this.name = aColumnArray[0];
            this.smiles = aColumnArray[1];
            this.targetRule = Integer.parseInt(aColumnArray[3]);
            this.expectation = aColumnArray[4];
            this.protectedAtoms = aColumnArray[5];
            this.mortarRows = new ArrayList<>();
            this.originalRows = new ArrayList<>();
        }
    }
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public CRUSHFragmenterReferenceTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Compares the steps of all molecules of the BIOFACQUIM data set and of up to 20 intermediate fragments per molecule
     * with the reference steps, which must match exactly, including stereochemistry and protected atoms. Two kinds of
     * differences are counted separately, with the numbers of this data set:
     * <ul>
     *     <li>16 molecules (intermediate fragments of CNP0383187.1 and CNP0436828.1 and CNP0436828.1 itself) with steps of
     *     rules 6 and 7 at arabinopyranosyl and 1,5-anhydrohexitol units whose products, 1,5-anhydroarabinitol units
     *     (saturated with hydrogen atoms), differ only in the configuration of C3: it is not stereogenic in these
     *     products because of the configuration of C2 and C4, so the RDKit removes it, while the CDK classifies it as a
     *     potential (para) stereocentre and keeps it;</li>
     *     <li>5 molecules (intermediate fragments of CNP0561292.0, CNP0560805.0 and CNP0520177.1, which contain atoms with
     *     radical electrons) with products that are identical but have different unique SMILES codes in the CDK
     *     (e.g. C=C[CH2] and [CH2]C=C); they are compared by isomorphism.</li>
     * </ul>
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimStepReferenceTest() throws Exception {
        StepComparisonResult tmpResult;
        try (InputStream tmpStream = CRUSHFragmenterReferenceTest.class.getResourceAsStream(
                CRUSHFragmenterReferenceTest.BIOFACQUIM_STEPS_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            tmpResult = CRUSHFragmenterReferenceTest.compareSteps(new BufferedReader(new InputStreamReader(tmpStream,
                    StandardCharsets.UTF_8)), CRUSHFragmenterTest.createOriginalFragmenter(),
                    CRUSHFragmenterReferenceTest.MAXIMUM_REPORTED_MISMATCHES);
        }
        System.out.println(tmpResult.summary());
        Assertions.assertTrue(tmpResult.molecules() > 3500);
        Assertions.assertEquals(0, tmpResult.skipped(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.mismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.stereoMismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.originalMismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.protectedAtomMismatches(), tmpResult.report());
        Assertions.assertEquals(16, tmpResult.paraStereocentreDifferences(), tmpResult.summary());
        Assertions.assertEquals(5, tmpResult.isomorphismMatches(), tmpResult.summary());
    }
    //
    /**
     * Compares the routes of all molecules of the BIOFACQUIM data set with the reference routes: for all molecules whose
     * search is exhaustive in both implementations (at least 95 %), the final fragments must be an optimal fragment
     * set of the reference, including stereochemistry.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimRouteReferenceTest() throws Exception {
        RouteComparisonResult tmpResult;
        try (InputStream tmpStream = CRUSHFragmenterReferenceTest.class.getResourceAsStream(
                CRUSHFragmenterReferenceTest.BIOFACQUIM_ROUTES_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            tmpResult = CRUSHFragmenterReferenceTest.compareRoutes(new BufferedReader(new InputStreamReader(tmpStream,
                    StandardCharsets.UTF_8)), CRUSHFragmenterTest.createOriginalFragmenter());
        }
        System.out.println(tmpResult);
        Assertions.assertTrue(tmpResult.molecules() > 550);
        Assertions.assertTrue(tmpResult.bothExhaustive() >= 0.95 * tmpResult.molecules(), tmpResult.toString());
        Assertions.assertEquals(0, tmpResult.mismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.stereoMismatches(), tmpResult.report());
    }
    //
    /**
     * Checks the consistency of the routes and fragments of the original variant for all molecules of the BIOFACQUIM
     * data set, with link atoms: every heavy atom of the molecule is contained in exactly one final fragment, there are
     * two link atoms per cleaved bond carrying the rule number of the step, the final fragments are the leaves of the
     * route, and all route fragments are the fragments of all nodes of the route except the root and the cleaved
     * components of molecules with several connected components. For every tenth molecule whose search is exhaustive,
     * the fragments must not depend on the order of the atoms and bonds of the input.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimRouteConsistencyTest() throws Exception {
        CRUSHFragmenter tmpFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        CRUSHFragmenter tmpRouteFragmenter = CRUSHFragmenterTest.createOriginalFragmenter();
        tmpRouteFragmenter.setFragmentOutputSetting(CRUSHFragmenter.FragmentOutputOption.ALL_ROUTE_FRAGMENTS);
        java.util.Random tmpRandom = new java.util.Random(42);
        int tmpMolecules = 0;
        int tmpShuffledMolecules = 0;
        try (InputStream tmpStream = CRUSHFragmenterReferenceTest.class.getResourceAsStream(
                CRUSHFragmenterReferenceTest.BIOFACQUIM_ROUTES_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t", -1);
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpColumns[1]);
                List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                CRUSHFragmenter.CRUSHRouteNode tmpRoot = tmpFragmenter.buildRoute(tmpMolecule);
                int tmpCleavedBonds = 0;
                List<CRUSHFragmenter.CRUSHRouteNode> tmpNodes = new ArrayList<>(tmpRoot.getAllDescendants());
                tmpNodes.add(tmpRoot);
                for (CRUSHFragmenter.CRUSHRouteNode tmpNode : tmpNodes) {
                    if (tmpNode.getRuleNumber() != 0) {
                        tmpCleavedBonds += tmpNode.getChildren().size() - 1;
                    }
                }
                int tmpHeavyAtoms = 0;
                int tmpLinkAtoms = 0;
                for (IAtomContainer tmpFragment : tmpFragments) {
                    for (IAtom tmpAtom : tmpFragment.atoms()) {
                        Integer tmpRuleNumber = tmpAtom.getProperty(CRUSHFragmenter.RULE_NUMBER_PROPERTY_KEY);
                        if (tmpRuleNumber != null) {
                            tmpLinkAtoms++;
                            Assertions.assertTrue(tmpRuleNumber >= 1 && tmpRuleNumber <= 33);
                            Assertions.assertEquals(tmpRuleNumber, tmpAtom.getMassNumber());
                        } else if (tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() > 1) {
                            tmpHeavyAtoms++;
                        }
                    }
                }
                int tmpExpectedHeavyAtoms = 0;
                for (IAtom tmpAtom : tmpMolecule.atoms()) {
                    if (tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() > 1) {
                        tmpExpectedHeavyAtoms++;
                    }
                }
                Assertions.assertEquals(tmpExpectedHeavyAtoms, tmpHeavyAtoms, tmpColumns[0]);
                Assertions.assertEquals(2 * tmpCleavedBonds, tmpLinkAtoms, tmpColumns[0]);
                if (tmpRoot.isLeaf()) {
                    Assertions.assertEquals(1, tmpFragments.size(), tmpColumns[0]);
                } else {
                    Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                            CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(tmpRoot.getLeaves()),
                                    false), tmpColumns[0]);
                    Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(
                                    CRUSHFragmenterTest.allRouteFragmentNodes(tmpRoot)), false),
                            CRUSHFragmenterTest.toCanonicalSmiles(tmpRouteFragmenter.fragmentMolecule(tmpMolecule), false),
                            tmpColumns[0]);
                }
                if (tmpMolecules % 10 == 0 && tmpRoot.isExhaustive()) {
                    tmpShuffledMolecules++;
                    List<IAtomContainer> tmpShuffledFragments = tmpFragmenter.fragmentMolecule(
                            CRUSHFragmenterTest.shuffle(tmpMolecule, tmpRandom));
                    Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                            CRUSHFragmenterTest.toCanonicalSmiles(tmpShuffledFragments, false), tmpColumns[0]);
                    if (!CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false).equals(
                            CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpShuffledFragments, false))) {
                        Assertions.assertTrue(CRUSHFragmenterTest.isIsomorphicWithStereo(tmpFragments, tmpShuffledFragments)
                                && CRUSHFragmenterTest.isIsomorphicWithStereo(tmpShuffledFragments, tmpFragments),
                                tmpColumns[0]);
                    }
                }
                tmpMolecules++;
            }
        }
        System.out.printf("%d molecules, %d with shuffled atoms%n", tmpMolecules, tmpShuffledMolecules);
        Assertions.assertTrue(tmpMolecules > 550);
    }
    //
    /**
     * Compares the fragments of the single-pass variant with default settings for all molecules of the BIOFACQUIM data
     * set with the single-pass reference: the fragments (with link atoms carrying the rule numbers, and with
     * stereochemistry) must be one of the results that the toolkit-independent sort keys of the processing order allow.
     * For 4 molecules, the reference lists two possible results, i.e. their result depends on the canonical rank of the
     * bonds. Toolkit differences are handled as in {@link #biofacquimStepReferenceTest()}.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimSinglePassReferenceTest() throws Exception {
        SinglePassComparisonResult tmpResult;
        try (InputStream tmpStream = CRUSHFragmenterReferenceTest.class.getResourceAsStream(
                CRUSHFragmenterReferenceTest.BIOFACQUIM_SINGLE_PASS_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            tmpResult = CRUSHFragmenterReferenceTest.compareSinglePassResults(new BufferedReader(new InputStreamReader(
                    tmpStream, StandardCharsets.UTF_8)), new CRUSHFragmenter());
        }
        System.out.println(tmpResult);
        Assertions.assertEquals(583, tmpResult.molecules());
        Assertions.assertEquals(0, tmpResult.uncheckedMolecules(), tmpResult.toString());
        Assertions.assertEquals(0, tmpResult.mismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.stereoMismatches(), tmpResult.report());
        Assertions.assertEquals(4, tmpResult.severalResults(), tmpResult.toString());
    }
    //
    /**
     * Checks the fragments of the single-pass variant with default settings for all molecules of the BIOFACQUIM data
     * set: every heavy atom of the molecule is contained in exactly one fragment; there are two link atoms per cleaved
     * bond, carrying the rule number; every fragment with a link atom has at least the minimum fragment size of heavy
     * atoms, and at least the minimum linker size if it has two or more attachment points; the fragments are the leaves
     * of the route; and the fragments (unique SMILES codes with stereochemistry as MORTAR creates them, and with the rule
     * numbers) do not depend on the order of the atoms and bonds of the input. The same is checked with hydrogen
     * saturation (without link atoms) and with the aliphatic skeleton rules. With default settings, 259 of the 583
     * molecules are fragmented.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimSinglePassConsistencyTest() throws Exception {
        CRUSHFragmenter tmpAliphaticFragmenter = new CRUSHFragmenter();
        tmpAliphaticFragmenter.setAliphaticSkeletonRulesSetting(true);
        CRUSHFragmenter tmpHydrogenFragmenter = new CRUSHFragmenter();
        tmpHydrogenFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        java.util.Random tmpRandom = new java.util.Random(42);
        int tmpMolecules = 0;
        int tmpFragmentedMolecules = 0;
        int tmpIsomorphismMatches = 0;
        try (InputStream tmpStream = CRUSHFragmenterReferenceTest.class.getResourceAsStream(
                CRUSHFragmenterReferenceTest.BIOFACQUIM_SINGLE_PASS_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t", -1);
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpColumns[1]);
                tmpMolecules++;
                for (CRUSHFragmenter tmpFragmenter : new CRUSHFragmenter[] {new CRUSHFragmenter(), tmpAliphaticFragmenter,
                        tmpHydrogenFragmenter}) {
                    boolean tmpIsDefault = tmpFragmenter.getSaturationSetting() == CRUSHFragmenter.SaturationOption.LINK_ATOMS
                            && !tmpFragmenter.getAliphaticSkeletonRulesSetting();
                    int tmpCleavedBonds = CRUSHFragmenterReferenceTest.checkSinglePassFragments(tmpFragmenter, tmpMolecule,
                            tmpColumns[0]);
                    if (tmpIsDefault && tmpCleavedBonds > 0) {
                        tmpFragmentedMolecules++;
                    }
                    List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                    List<IAtomContainer> tmpShuffledFragments = tmpFragmenter.fragmentMolecule(
                            CRUSHFragmenterTest.shuffle(tmpMolecule, tmpRandom));
                    Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                            CRUSHFragmenterTest.toCanonicalSmiles(tmpShuffledFragments, false), tmpColumns[0]);
                    if (!CRUSHFragmenterReferenceTest.createUniqueSmilesList(tmpFragments).equals(
                            CRUSHFragmenterReferenceTest.createUniqueSmilesList(tmpShuffledFragments))) {
                        //the unique SMILES codes of the CDK are not canonical for some fragments, see
                        //CRUSHFragmenterTest.atomOrderInvarianceTest(), these are compared by isomorphism
                        tmpIsomorphismMatches++;
                        Assertions.assertTrue(CRUSHFragmenterTest.isIsomorphicWithStereo(tmpFragments, tmpShuffledFragments)
                                && CRUSHFragmenterTest.isIsomorphicWithStereo(tmpShuffledFragments, tmpFragments),
                                tmpColumns[0]);
                    }
                }
            }
        }
        System.out.printf("%d molecules, %d fragmented with default settings, %d fragment lists of shuffled molecules"
                + " matched by isomorphism%n", tmpMolecules, tmpFragmentedMolecules, tmpIsomorphismMatches);
        Assertions.assertEquals(583, tmpMolecules);
        Assertions.assertEquals(259, tmpFragmentedMolecules);
    }
    //
    /**
     * Compares the steps with the steps reference file given by the environment variable
     * {@value #REFERENCE_FILE_ENVIRONMENT_VARIABLE}; skipped if it is not set. Large and diverse data sets contain a
     * few molecules for which the toolkits differ (e.g. in the aromaticity perception), so up to one mismatch per
     * thousand molecules is tolerated; the mismatches are printed. On a sample of every 100th molecule of ChEMBL 37
     * (28,978 molecules), all differences are caused by the toolkits:
     * <ul>
     *     <li>3 molecules with different steps: aromaticity of tellurium-containing rings (CHEMBL336433, CHEMBL516219)
     *     and of a cyclic imide in a chlorin (CHEMBL437550), which only one of the toolkits perceives as aromatic;</li>
     *     <li>15 molecules with different stereochemistry: 9 with products that are identical including stereochemistry
     *     (isomorphism) but whose unique SMILES codes of the CDK differ, i.e. symmetric cations drawn with the charge at
     *     the other atom (imidazolium, cyanine and rhodamine dyes) and naphthoquinone dioximes; 3 sulfoxides with a ring
     *     closure at the sulfur atom, whose SMILES the CDK and the RDKit interpret differently (the CDK reads the RDKit
     *     SMILES of the unchanged molecule as the other stereoisomer); 2 squaric acid bis-imines, for which the CDK SMILES
     *     parser creates no stereo element at the aromatic four-membered ring; and 1 macrocycle (CHEMBL5415612) whose
     *     reference product has the wrong configuration (the double bond is Z in the molecule and must stay Z);</li>
     *     <li>14 molecules with para stereocentre differences and 3 with products matched by isomorphism, see
     *     {@link #biofacquimStepReferenceTest()}.</li>
     * </ul>
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void additionalReferenceTest() throws Exception {
        String tmpPath = System.getenv(CRUSHFragmenterReferenceTest.REFERENCE_FILE_ENVIRONMENT_VARIABLE);
        Assumptions.assumeTrue(tmpPath != null && !tmpPath.isBlank(), "No additional reference file given.");
        StepComparisonResult tmpResult;
        try (BufferedReader tmpReader = Files.newBufferedReader(Path.of(tmpPath), StandardCharsets.UTF_8)) {
            tmpResult = CRUSHFragmenterReferenceTest.compareSteps(tmpReader, CRUSHFragmenterTest.createOriginalFragmenter(),
                    CRUSHFragmenterReferenceTest.MAXIMUM_REPORTED_MISMATCHES);
        }
        System.out.println(tmpResult.summary());
        System.out.println(tmpResult.report());
        int tmpDifferentMolecules = tmpResult.skipped() + tmpResult.mismatches() + tmpResult.stereoMismatches()
                + tmpResult.protectedAtomMismatches();
        Assertions.assertTrue(tmpDifferentMolecules <= tmpResult.molecules() / 1000, tmpResult.report());
    }
    //
    /**
     * Compares the routes with the routes reference file given by the environment variable
     * {@value #ROUTE_REFERENCE_FILE_ENVIRONMENT_VARIABLE}; skipped if it is not set. Up to one mismatch per thousand
     * molecules is tolerated; the mismatches are printed. On a sample of every 2000th molecule of ChEMBL 37 (1449
     * molecules, 1435 searches exhaustive in both implementations), 2 molecules have final fragments that are identical
     * to an optimal set but have different unique SMILES codes in the CDK (a sulfoximine and a phosphazene ring with
     * another Kekulé structure; compared by isomorphism), and 1 molecule differs in stereochemistry (CHEMBL19398: the
     * CDK SMILES parser creates no stereo element for the exocyclic C=N bonds of the aromatic four-membered ring, the
     * RDKit does).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void additionalRouteReferenceTest() throws Exception {
        String tmpPath = System.getenv(CRUSHFragmenterReferenceTest.ROUTE_REFERENCE_FILE_ENVIRONMENT_VARIABLE);
        Assumptions.assumeTrue(tmpPath != null && !tmpPath.isBlank(), "No additional route reference file given.");
        RouteComparisonResult tmpResult;
        try (BufferedReader tmpReader = Files.newBufferedReader(Path.of(tmpPath), StandardCharsets.UTF_8)) {
            tmpResult = CRUSHFragmenterReferenceTest.compareRoutes(tmpReader, CRUSHFragmenterTest.createOriginalFragmenter());
        }
        System.out.println(tmpResult);
        Assertions.assertTrue(tmpResult.mismatches() + tmpResult.stereoMismatches() <= tmpResult.molecules() / 1000,
                tmpResult.report());
    }
    //
    /**
     * Compares the steps of all molecules of the given steps reference file (format of CRUSH_Reference_Generator.py,
     * mode "steps") with the steps of the given fragmenter.
     *
     * @param aReader reader of the reference file
     * @param aFragmenter fragmenter with the settings of the reference (original variant, minimum fragment size 3)
     * @param aMaximumReportedCount maximum number of reported mismatches
     * @return comparison result
     * @throws IOException if the file cannot be read
     */
    static StepComparisonResult compareSteps(BufferedReader aReader, CRUSHFragmenter aFragmenter, int aMaximumReportedCount)
            throws IOException {
        int[] tmpCounts = new int[10];
        long[] tmpNanoseconds = new long[1];
        StringBuilder tmpReport = new StringBuilder();
        StepBlock tmpBlock = null;
        String tmpLine;
        while ((tmpLine = aReader.readLine()) != null) {
            if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                continue;
            }
            String[] tmpColumns = tmpLine.split("\t", -1);
            if (tmpColumns[2].equals("M")) {
                if (tmpBlock != null) {
                    CRUSHFragmenterReferenceTest.compareStepBlock(tmpBlock, aFragmenter, tmpCounts, tmpNanoseconds,
                            tmpReport, aMaximumReportedCount);
                }
                tmpBlock = new StepBlock(tmpColumns);
            } else if (tmpColumns[2].equals("2")) {
                Objects.requireNonNull(tmpBlock).mortarRows.add(tmpColumns);
            } else {
                Objects.requireNonNull(tmpBlock).originalRows.add(tmpColumns);
            }
        }
        if (tmpBlock != null) {
            CRUSHFragmenterReferenceTest.compareStepBlock(tmpBlock, aFragmenter, tmpCounts, tmpNanoseconds, tmpReport,
                    aMaximumReportedCount);
        }
        return new StepComparisonResult(tmpCounts[0], tmpCounts[1], tmpCounts[2], tmpCounts[3], tmpCounts[4],
                tmpCounts[5], tmpCounts[6], tmpCounts[7], tmpCounts[8], tmpCounts[9], tmpNanoseconds[0] / 1_000_000L,
                tmpReport.toString());
    }
    //
    /**
     * Compares the steps of one molecule of a steps reference file with the steps of the given fragmenter. For every
     * step (rule and cleaved bonds), the products are compared by their canonical SMILES codes; if these differ, the
     * products are compared by graph isomorphism, since the unique SMILES codes of the CDK are not canonical for some
     * radicals (e.g. C=C[CH2] and [CH2]C=C). Products are compared with stereochemistry for the steps whose products
     * are identical without stereochemistry; a difference that disappears when the configurations of potential (para)
     * stereocentres are removed on both sides is counted separately (para stereocentres, e.g. C3 of 1,5-anhydropentitols:
     * the RDKit removes their configuration if it is not stereogenic because of the configuration of other atoms, the
     * CDK keeps it).
     *
     * @param aBlock the molecule with its reference rows
     * @param aFragmenter the fragmenter
     * @param aCountArray counts: molecules, skipped, mismatches, stereo mismatches, mismatches against the published
     *                    procedure, stereo differences against the published procedure, protected atom mismatches,
     *                    expectation mismatches, molecules with para stereocentre differences, molecules with products
     *                    matched by isomorphism
     * @param aNanosecondArray time spent in the determination of the steps
     * @param aReport report
     * @param aMaximumReportedCount maximum number of reported mismatches
     */
    private static void compareStepBlock(StepBlock aBlock, CRUSHFragmenter aFragmenter, int[] aCountArray,
                                         long[] aNanosecondArray, StringBuilder aReport, int aMaximumReportedCount) {
        aCountArray[0]++;
        Map<String, String> tmpActualByStep = new HashMap<>();
        Map<String, String> tmpActualStereoByStep = new HashMap<>();
        Map<String, List<IAtomContainer>> tmpActualFragmentsByStep = new HashMap<>();
        Map<String, String> tmpExpectedRawByStep = new HashMap<>();
        Map<String, String> tmpExpectedByStep = new HashMap<>();
        Map<String, String> tmpExpectedStereoByStep = new HashMap<>();
        String tmpActualProtectedAtoms;
        try {
            IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(aBlock.smiles);
            int tmpParsedAtomCount = tmpMolecule.getAtomCount();
            aFragmenter.prepareMolecule(tmpMolecule);
            if (tmpParsedAtomCount != tmpMolecule.getAtomCount()) {
                throw new CDKException("Atoms were removed by the preparation, the atom indices cannot be compared.");
            }
            tmpActualProtectedAtoms = CRUSHFragmenterReferenceTest.join(aFragmenter.findProtectedAtoms(tmpMolecule));
            int[] tmpComponentLabels = new ConnectedComponents(GraphUtil.toAdjList(tmpMolecule)).components();
            int tmpComponentCount = Arrays.stream(tmpComponentLabels).max().orElse(0);
            int[] tmpNoLabels = new int[tmpMolecule.getBondCount()];
            for (int tmpComponent = 1; tmpComponent <= tmpComponentCount; tmpComponent++) {
                BitSet tmpAtomSet = new BitSet(tmpMolecule.getAtomCount());
                for (int i = 0; i < tmpComponentLabels.length; i++) {
                    if (tmpComponentLabels[i] == tmpComponent) {
                        tmpAtomSet.set(i);
                    }
                }
                long tmpStart = System.nanoTime();
                List<CRUSHFragmenter.CleavageStep> tmpSteps = aFragmenter.findCleavageSteps(tmpMolecule, tmpAtomSet);
                aNanosecondArray[0] += System.nanoTime() - tmpStart;
                for (CRUSHFragmenter.CleavageStep tmpStep : tmpSteps) {
                    List<IAtomContainer> tmpProducts = new ArrayList<>(tmpStep.parts().length);
                    for (BitSet tmpPart : tmpStep.parts()) {
                        tmpProducts.add(CRUSHFragmenter.createFragment(tmpMolecule, tmpPart, tmpNoLabels, true));
                    }
                    String tmpStepKey = tmpStep.ruleNumber() + "\t" + CRUSHFragmenterReferenceTest.describeBonds(tmpMolecule,
                            tmpStep.cleavedBonds());
                    tmpActualByStep.put(tmpStepKey, String.join(".", CRUSHFragmenterTest.toCanonicalSmiles(tmpProducts, false)));
                    tmpActualStereoByStep.put(tmpStepKey, String.join(".",
                            CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpProducts, true)));
                    tmpActualFragmentsByStep.put(tmpStepKey, tmpProducts);
                }
            }
            for (String[] tmpRow : aBlock.mortarRows) {
                List<String> tmpProducts = Arrays.asList(tmpRow[6].split("\\."));
                String tmpStepKey = tmpRow[3] + "\t" + tmpRow[5];
                tmpExpectedRawByStep.put(tmpStepKey, tmpRow[6]);
                tmpExpectedByStep.put(tmpStepKey, String.join(".", CRUSHFragmenterTest.canonical(tmpProducts, false)));
                tmpExpectedStereoByStep.put(tmpStepKey, String.join(".", CRUSHFragmenterTest.toStereoSmiles(tmpProducts,
                        true)));
            }
        } catch (CDKException | CloneNotSupportedException | RuntimeException anException) {
            aCountArray[1]++;
            CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                    String.format("%s %s: skipped, %s%n", aBlock.name, aBlock.smiles, anException));
            return;
        }
        if (!tmpActualProtectedAtoms.equals(aBlock.protectedAtoms)) {
            aCountArray[6]++;
            CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                    String.format("%s %s: protected atoms expected [%s], actual [%s]%n", aBlock.name, aBlock.smiles,
                            aBlock.protectedAtoms, tmpActualProtectedAtoms));
        }
        List<String> tmpOnlyReference = new ArrayList<>();
        List<String> tmpOnlyFragmenter = new ArrayList<>();
        List<String> tmpStereoDifferences = new ArrayList<>();
        boolean tmpIsParaDifference = false;
        boolean tmpIsIsomorphismMatch = false;
        List<String> tmpStepKeys = new ArrayList<>(tmpExpectedByStep.keySet());
        for (String tmpStepKey : tmpActualByStep.keySet()) {
            if (!tmpExpectedByStep.containsKey(tmpStepKey)) {
                tmpStepKeys.add(tmpStepKey);
            }
        }
        Collections.sort(tmpStepKeys);
        for (String tmpStepKey : tmpStepKeys) {
            String tmpActual = tmpActualByStep.get(tmpStepKey);
            String tmpExpected = tmpExpectedByStep.get(tmpStepKey);
            if (tmpActual == null || tmpExpected == null) {
                if (tmpExpected != null) {
                    tmpOnlyReference.add(tmpStepKey + "\t" + tmpExpected);
                } else {
                    tmpOnlyFragmenter.add(tmpStepKey + "\t" + tmpActual);
                }
                continue;
            }
            if (!tmpActual.equals(tmpExpected)) {
                if (CRUSHFragmenterReferenceTest.isIsomorphic(tmpExpectedRawByStep.get(tmpStepKey),
                        tmpActualFragmentsByStep.get(tmpStepKey))) {
                    tmpIsIsomorphismMatch = true;
                } else {
                    tmpOnlyReference.add(tmpStepKey + "\t" + tmpExpected);
                    tmpOnlyFragmenter.add(tmpStepKey + "\t" + tmpActual);
                }
                continue;
            }
            String tmpActualStereo = tmpActualStereoByStep.get(tmpStepKey);
            String tmpExpectedStereo = tmpExpectedStereoByStep.get(tmpStepKey);
            if (!tmpActualStereo.equals(tmpExpectedStereo)) {
                if (CRUSHFragmenterReferenceTest.isEqualWithoutParaStereocentres(tmpExpectedRawByStep.get(tmpStepKey),
                        tmpActualFragmentsByStep.get(tmpStepKey))) {
                    tmpIsParaDifference = true;
                    if (aCountArray[8] < 5) {
                        System.out.printf("para stereocentre difference: %s %s, step %s%n  reference:  %s%n  fragmenter: %s%n",
                                aBlock.name, aBlock.smiles, tmpStepKey, tmpExpectedStereo, tmpActualStereo);
                    }
                } else {
                    tmpStereoDifferences.add(String.format("%s: reference %s, fragmenter %s", tmpStepKey,
                            tmpExpectedStereo, tmpActualStereo));
                }
            }
        }
        if (tmpIsParaDifference) {
            aCountArray[8]++;
        }
        if (tmpIsIsomorphismMatch) {
            aCountArray[9]++;
        }
        if (!tmpOnlyReference.isEmpty() || !tmpOnlyFragmenter.isEmpty()) {
            aCountArray[2]++;
            CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                    String.format("%s %s%n  only reference: %s%n  only fragmenter: %s%n", aBlock.name, aBlock.smiles,
                            tmpOnlyReference, tmpOnlyFragmenter));
        } else if (!tmpStereoDifferences.isEmpty()) {
            aCountArray[3]++;
            CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                    String.format("%s %s (stereo)%n  %s%n", aBlock.name, aBlock.smiles,
                            String.join(String.format("%n  "), tmpStereoDifferences)));
        }
        boolean tmpIsOriginalMismatch = false;
        for (String[] tmpRow : aBlock.originalRows) {
            String tmpStepKey = tmpRow[3] + "\t" + tmpRow[5];
            if (tmpRow[4].equals("D")) {
                //valid step of the published procedure whose products differ only in hydrogen counts, since the
                //reactions of the RDKit recompute them from the default valences (atoms with radical electrons,
                //hypervalent atoms with hydrogen atoms); the step must exist, its products are compared with part 2
                if (!tmpActualByStep.containsKey(tmpStepKey)) {
                    tmpIsOriginalMismatch = true;
                    CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                            String.format("%s %s: valid step of the published procedure %s (radical) is missing%n",
                                    aBlock.name, aBlock.smiles, tmpStepKey));
                }
                continue;
            }
            if (!tmpRow[4].equals("V")) {
                continue;
            }
            String tmpProductText = tmpRow[6].equals("=") ? tmpExpectedRawByStep.get(tmpStepKey) : tmpRow[6];
            String tmpActualProducts = tmpActualByStep.get(tmpStepKey);
            try {
                List<String> tmpProducts = Arrays.asList(Objects.requireNonNull(tmpProductText).split("\\."));
                if (tmpActualProducts == null
                        || (!tmpActualProducts.equals(String.join(".", CRUSHFragmenterTest.canonical(tmpProducts, false)))
                        && !CRUSHFragmenterReferenceTest.isIsomorphic(tmpProductText,
                        tmpActualFragmentsByStep.get(tmpStepKey)))) {
                    tmpIsOriginalMismatch = true;
                    CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                            String.format("%s %s: valid step of the published procedure %s -> %s, fragmenter: %s%n",
                                    aBlock.name, aBlock.smiles, tmpStepKey, tmpProductText, tmpActualProducts));
                } else if (!tmpActualStereoByStep.get(tmpStepKey).equals(String.join(".",
                        CRUSHFragmenterTest.toStereoSmiles(tmpProducts, true)))) {
                    aCountArray[5]++;
                }
            } catch (CDKException | RuntimeException anException) {
                tmpIsOriginalMismatch = true;
                CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                        String.format("%s %s: products of the published procedure %s not comparable, %s%n", aBlock.name,
                                aBlock.smiles, tmpProductText, anException));
            }
        }
        if (tmpIsOriginalMismatch) {
            aCountArray[4]++;
        }
        if (aBlock.targetRule > 0 && (aBlock.expectation.equals("positive") || aBlock.expectation.equals("negative"))) {
            boolean tmpIsFound = false;
            for (String tmpStepKey : tmpActualByStep.keySet()) {
                if (tmpStepKey.startsWith(aBlock.targetRule + "\t")) {
                    tmpIsFound = true;
                    break;
                }
            }
            if (tmpIsFound != aBlock.expectation.equals("positive")) {
                aCountArray[7]++;
                CRUSHFragmenterReferenceTest.report(aReport, aCountArray, aMaximumReportedCount,
                        String.format("%s %s: expectation %s for rule %d not met%n", aBlock.name, aBlock.smiles,
                                aBlock.expectation, aBlock.targetRule));
            }
        }
    }
    //
    /**
     * Returns whether the molecules of the given SMILES codes (separated by '.') are isomorphic to the given fragments,
     * as multisets: elements, formal charges, implicit hydrogen counts, mass numbers, and bond orders (or aromaticity,
     * Daylight model) must match.
     *
     * @param aSmiles SMILES codes separated by '.'
     * @param aFragmentList fragments
     * @return true if there is a one-to-one assignment of isomorphic molecules
     */
    static boolean isIsomorphic(String aSmiles, List<IAtomContainer> aFragmentList) {
        try {
            List<IAtomContainer> tmpRemaining = new ArrayList<>(aFragmentList.size());
            for (IAtomContainer tmpFragment : aFragmentList) {
                IAtomContainer tmpClone = tmpFragment.clone();
                SmartsPattern.prepare(tmpClone);
                tmpRemaining.add(tmpClone);
            }
            String[] tmpSmilesArray = aSmiles.split("\\.");
            if (tmpSmilesArray.length != tmpRemaining.size()) {
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
            for (String tmpSmiles : tmpSmilesArray) {
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpSmiles);
                SmartsPattern.prepare(tmpMolecule);
                Pattern tmpPattern = VentoFoggia.findIdentical(tmpMolecule, tmpAtomMatcher, tmpBondMatcher);
                IAtomContainer tmpMatch = null;
                for (IAtomContainer tmpCandidate : tmpRemaining) {
                    if (tmpCandidate.getAtomCount() == tmpMolecule.getAtomCount() && tmpPattern.matches(tmpCandidate)) {
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
        } catch (CDKException | CloneNotSupportedException anException) {
            return false;
        }
    }
    //
    /**
     * Returns whether the molecules of the given SMILES codes (separated by '.') and the given fragments have the same
     * unique SMILES codes with stereochemistry if the configurations of potential (para) stereocentres according to the
     * CDK, and of atoms that are not stereogenic, are removed on both sides.
     *
     * @param aSmiles SMILES codes separated by '.'
     * @param aFragmentList fragments
     * @return true if the SMILES codes are equal without the configurations of para stereocentres
     */
    static boolean isEqualWithoutParaStereocentres(String aSmiles, List<IAtomContainer> aFragmentList) {
        try {
            List<String> tmpExpected = new ArrayList<>();
            for (String tmpSmiles : aSmiles.split("\\.")) {
                tmpExpected.add(CRUSHFragmenterTest.toStereoSmilesWithoutParaStereocentres(
                        CRUSHFragmenterTest.parse(tmpSmiles)));
            }
            List<String> tmpActual = new ArrayList<>();
            for (IAtomContainer tmpFragment : aFragmentList) {
                tmpActual.add(CRUSHFragmenterTest.toStereoSmilesWithoutParaStereocentres(tmpFragment));
            }
            Collections.sort(tmpExpected);
            Collections.sort(tmpActual);
            return tmpExpected.equals(tmpActual);
        } catch (CDKException anException) {
            return false;
        }
    }
    //
    /**
     * Compares the final fragments of all molecules of the given routes reference file (format of
     * CRUSH_Reference_Generator.py, mode "routes") with the final fragments of the given fragmenter. If the canonical
     * SMILES codes differ, the fragments are compared by graph isomorphism, as in the step comparison, since the unique
     * SMILES codes of the CDK are not canonical for some Kekulé structures (e.g. sulfoximine and phosphazene rings);
     * such molecules are counted separately.
     *
     * @param aReader reader of the reference file
     * @param aFragmenter fragmenter with the settings of the reference (original variant, minimum fragment size 3);
     *                    its saturation setting is changed to hydrogen saturation
     * @return comparison result
     * @throws IOException if the file cannot be read
     */
    static RouteComparisonResult compareRoutes(BufferedReader aReader, CRUSHFragmenter aFragmenter) throws IOException {
        aFragmenter.setSaturationSetting(CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION);
        int tmpMolecules = 0;
        int tmpBothExhaustive = 0;
        int tmpReferenceExhaustive = 0;
        int tmpFragmenterExhaustive = 0;
        int tmpMismatches = 0;
        int tmpStereoMismatches = 0;
        int tmpUncheckedSets = 0;
        int tmpIsomorphismMatches = 0;
        long tmpNanoseconds = 0;
        long tmpSlowestNanoseconds = 0;
        String tmpSlowestMolecule = "";
        StringBuilder tmpReport = new StringBuilder();
        String tmpLine;
        while ((tmpLine = aReader.readLine()) != null) {
            if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                continue;
            }
            String[] tmpColumns = tmpLine.split("\t", -1);
            tmpMolecules++;
            boolean tmpIsReferenceExhaustive = tmpColumns[3].equals("exhaustive");
            if (tmpIsReferenceExhaustive) {
                tmpReferenceExhaustive++;
            }
            try {
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpColumns[1]);
                long tmpStart = System.nanoTime();
                List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(tmpMolecule);
                long tmpTime = System.nanoTime() - tmpStart;
                CRUSHFragmenter.CRUSHRouteNode tmpRoot = aFragmenter.buildRoute(tmpMolecule);
                tmpNanoseconds += tmpTime;
                if (tmpTime > tmpSlowestNanoseconds) {
                    tmpSlowestNanoseconds = tmpTime;
                    tmpSlowestMolecule = tmpColumns[0];
                }
                if (tmpRoot.isExhaustive()) {
                    tmpFragmenterExhaustive++;
                }
                if (!tmpIsReferenceExhaustive || !tmpRoot.isExhaustive()) {
                    continue;
                }
                tmpBothExhaustive++;
                int tmpExpectedCount = Integer.parseInt(tmpColumns[2]);
                List<String> tmpActual = CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false);
                List<String> tmpActualStereo = CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, true);
                if (tmpActual.size() != tmpExpectedCount) {
                    tmpMismatches++;
                    CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                            String.format("%s %s: %d final fragments, expected %d%n  actual: %s%n", tmpColumns[0],
                                    tmpColumns[1], tmpActual.size(), tmpExpectedCount, tmpActual));
                    continue;
                }
                if (tmpColumns[5].equals("-")) {
                    tmpUncheckedSets++;
                    continue;
                }
                boolean tmpIsFound = false;
                boolean tmpIsStereoFound = false;
                boolean tmpIsIsomorphismMatch = false;
                for (String tmpSet : tmpColumns[6].split(" ")) {
                    List<String> tmpSmilesList = Arrays.asList(tmpSet.split("\\."));
                    if (CRUSHFragmenterTest.canonical(tmpSmilesList, false).equals(tmpActual)) {
                        tmpIsFound = true;
                        if (CRUSHFragmenterTest.toStereoSmiles(tmpSmilesList, true).equals(tmpActualStereo)) {
                            tmpIsStereoFound = true;
                            break;
                        }
                    } else if (CRUSHFragmenterReferenceTest.isIsomorphic(tmpSet, tmpFragments)) {
                        //identical molecules whose unique SMILES codes of the CDK differ (as in the step comparison);
                        //their stereochemistry is compared by isomorphism as well
                        tmpIsFound = true;
                        tmpIsIsomorphismMatch = true;
                        List<IAtomContainer> tmpReferenceFragments = new ArrayList<>(tmpSmilesList.size());
                        for (String tmpSmiles : tmpSmilesList) {
                            tmpReferenceFragments.add(CRUSHFragmenterTest.parse(tmpSmiles));
                        }
                        if (CRUSHFragmenterTest.isIsomorphicWithStereo(tmpReferenceFragments, tmpFragments)
                                && CRUSHFragmenterTest.isIsomorphicWithStereo(tmpFragments, tmpReferenceFragments)) {
                            tmpIsStereoFound = true;
                            break;
                        }
                    }
                }
                if (tmpIsFound && tmpIsIsomorphismMatch) {
                    tmpIsomorphismMatches++;
                }
                if (!tmpIsFound) {
                    tmpMismatches++;
                    CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                            String.format("%s %s: final fragments not among the %s optimal sets%n  actual: %s%n  "
                                    + "reference: %s%n", tmpColumns[0], tmpColumns[1], tmpColumns[5], tmpActual,
                                    tmpColumns[6]));
                } else if (!tmpIsStereoFound) {
                    tmpStereoMismatches++;
                    CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                            String.format("%s %s (stereo): final fragments not among the %s optimal sets%n  actual: "
                                    + "%s%n  reference: %s%n", tmpColumns[0], tmpColumns[1], tmpColumns[5],
                                    tmpActualStereo, tmpColumns[6]));
                }
            } catch (CDKException | CloneNotSupportedException | IllegalArgumentException anException) {
                tmpMismatches++;
                CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                        String.format("%s %s: error, %s%n", tmpColumns[0], tmpColumns[1], anException));
            }
        }
        return new RouteComparisonResult(tmpMolecules, tmpBothExhaustive, tmpReferenceExhaustive, tmpFragmenterExhaustive,
                tmpMismatches, tmpStereoMismatches, tmpUncheckedSets, tmpIsomorphismMatches, tmpNanoseconds / 1_000_000L,
                tmpSlowestNanoseconds / 1_000_000L, tmpSlowestMolecule, tmpReport.toString());
    }
    //
    /**
     * Compares the fragments of all molecules of the given single-pass reference file (format of
     * CRUSH_Reference_Generator.py, mode "singlepass") with the fragments of the given fragmenter: they must be one of
     * the listed possible results, compared by canonical SMILES codes with the rule numbers of the link atoms and by
     * unique SMILES codes with stereochemistry (link atoms without rule numbers, stereo elements of atoms that are not
     * stereogenic according to the CDK removed on both sides). As in the step comparison, fragments whose canonical
     * SMILES codes differ are compared by graph isomorphism (non-canonical unique SMILES codes of the CDK), and
     * differences only in the configuration of potential (para) stereocentres are counted separately.
     *
     * @param aReader reader of the reference file
     * @param aFragmenter fragmenter with the settings of the reference (single-pass variant, link atoms)
     * @return comparison result
     * @throws IOException if the file cannot be read
     */
    static SinglePassComparisonResult compareSinglePassResults(BufferedReader aReader, CRUSHFragmenter aFragmenter)
            throws IOException {
        int tmpMolecules = 0;
        int tmpSeveralResults = 0;
        int tmpUncheckedMolecules = 0;
        int tmpMismatches = 0;
        int tmpStereoMismatches = 0;
        int tmpParaStereocentreDifferences = 0;
        int tmpIsomorphismMatches = 0;
        long tmpNanoseconds = 0;
        StringBuilder tmpReport = new StringBuilder();
        String tmpLine;
        while ((tmpLine = aReader.readLine()) != null) {
            if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                continue;
            }
            String[] tmpColumns = tmpLine.split("\t", -1);
            tmpMolecules++;
            if (tmpColumns[3].equals("-")) {
                tmpUncheckedMolecules++;
                continue;
            }
            if (Integer.parseInt(tmpColumns[3]) > 1) {
                tmpSeveralResults++;
            }
            try {
                IAtomContainer tmpMolecule = CRUSHFragmenterTest.parse(tmpColumns[1]);
                long tmpStart = System.nanoTime();
                List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(tmpMolecule);
                tmpNanoseconds += System.nanoTime() - tmpStart;
                List<String> tmpActual = CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false);
                List<String> tmpActualStereo = CRUSHFragmenterTest.toStereoSmilesOfFragments(tmpFragments, true);
                boolean tmpIsFound = false;
                boolean tmpIsStereoFound = false;
                boolean tmpIsParaDifference = false;
                boolean tmpIsIsomorphismMatch = false;
                for (String tmpSet : tmpColumns[4].split(" ")) {
                    List<String> tmpSmilesList = Arrays.asList(tmpSet.split("\\."));
                    if (CRUSHFragmenterTest.canonical(tmpSmilesList, false).equals(tmpActual)) {
                        tmpIsFound = true;
                        if (CRUSHFragmenterTest.toStereoSmiles(tmpSmilesList, true).equals(tmpActualStereo)) {
                            tmpIsStereoFound = true;
                            break;
                        }
                        if (CRUSHFragmenterReferenceTest.isEqualWithoutParaStereocentres(tmpSet, tmpFragments)) {
                            tmpIsParaDifference = true;
                        }
                    } else if (CRUSHFragmenterReferenceTest.isIsomorphic(tmpSet, tmpFragments)) {
                        tmpIsFound = true;
                        tmpIsIsomorphismMatch = true;
                        List<IAtomContainer> tmpReferenceFragments = new ArrayList<>(tmpSmilesList.size());
                        for (String tmpSmiles : tmpSmilesList) {
                            tmpReferenceFragments.add(CRUSHFragmenterTest.parse(tmpSmiles));
                        }
                        if (CRUSHFragmenterTest.isIsomorphicWithStereo(tmpReferenceFragments, tmpFragments)
                                && CRUSHFragmenterTest.isIsomorphicWithStereo(tmpFragments, tmpReferenceFragments)) {
                            tmpIsStereoFound = true;
                            break;
                        }
                    }
                }
                if (tmpIsIsomorphismMatch) {
                    tmpIsomorphismMatches++;
                }
                if (!tmpIsFound) {
                    tmpMismatches++;
                    CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                            String.format("%s %s: fragments not among the %s possible results%n  actual: %s%n  "
                                    + "reference: %s%n", tmpColumns[0], tmpColumns[1], tmpColumns[3], tmpActual,
                                    tmpColumns[4]));
                } else if (!tmpIsStereoFound && tmpIsParaDifference) {
                    tmpParaStereocentreDifferences++;
                } else if (!tmpIsStereoFound) {
                    tmpStereoMismatches++;
                    CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                            String.format("%s %s (stereo): fragments not among the %s possible results%n  actual: "
                                    + "%s%n  reference: %s%n", tmpColumns[0], tmpColumns[1], tmpColumns[3],
                                    tmpActualStereo, tmpColumns[4]));
                }
            } catch (CDKException | CloneNotSupportedException | IllegalArgumentException anException) {
                tmpMismatches++;
                CRUSHFragmenterReferenceTest.report(tmpReport, tmpMismatches + tmpStereoMismatches,
                        String.format("%s %s: error, %s%n", tmpColumns[0], tmpColumns[1], anException));
            }
        }
        return new SinglePassComparisonResult(tmpMolecules, tmpSeveralResults, tmpUncheckedMolecules, tmpMismatches,
                tmpStereoMismatches, tmpParaStereocentreDifferences, tmpIsomorphismMatches, tmpNanoseconds / 1_000_000L,
                tmpReport.toString());
    }
    //
    /**
     * Checks the fragments of the given molecule of the given fragmenter (single-pass variant): the fragments are the
     * leaves of the route, every heavy atom of the molecule is contained in exactly one fragment, there are two link
     * atoms per cleaved bond carrying a rule number (none with hydrogen saturation), and every fragment with a cleaved
     * bond has at least the minimum fragment size of heavy atoms, and at least the minimum linker size if it has two or
     * more attachment points (cleaved bonds and pseudo atoms of the input). Connected components without a cleaved bond
     * are not checked for their size.
     *
     * @param aFragmenter fragmenter with the single-pass variant
     * @param aMolecule molecule, it is not altered
     * @param aName name of the molecule for the messages
     * @return number of cleaved bonds
     * @throws Exception if anything goes wrong
     */
    static int checkSinglePassFragments(CRUSHFragmenter aFragmenter, IAtomContainer aMolecule, String aName)
            throws Exception {
        Assertions.assertEquals(CRUSHFragmenter.VariantOption.SINGLE_PASS, aFragmenter.getVariantSetting());
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(aMolecule);
        CRUSHFragmenter.CRUSHRouteNode tmpRoot = aFragmenter.buildRoute(aMolecule);
        Assertions.assertTrue(tmpRoot.isExhaustive(), aName);
        Assertions.assertEquals(0, tmpRoot.getRuleNumber(), aName);
        IAtomContainer tmpPreparedMolecule = tmpRoot.getFragment();
        List<CRUSHFragmenter.CRUSHRouteNode> tmpLeaves = tmpRoot.isLeaf() ? List.of(tmpRoot) : tmpRoot.getChildren();
        for (CRUSHFragmenter.CRUSHRouteNode tmpLeaf : tmpLeaves) {
            Assertions.assertTrue(tmpLeaf.isLeaf(), aName);
        }
        Assertions.assertEquals(CRUSHFragmenterTest.toCanonicalSmiles(tmpFragments, false),
                CRUSHFragmenterTest.toCanonicalSmiles(CRUSHFragmenterTest.fragmentsOf(tmpLeaves), false), aName);
        int[] tmpLeafOfAtom = new int[tmpPreparedMolecule.getAtomCount()];
        Arrays.fill(tmpLeafOfAtom, -1);
        for (int i = 0; i < tmpLeaves.size(); i++) {
            BitSet tmpAtomSet = tmpLeaves.get(i).getAtomSet();
            for (int j = tmpAtomSet.nextSetBit(0); j >= 0; j = tmpAtomSet.nextSetBit(j + 1)) {
                Assertions.assertEquals(-1, tmpLeafOfAtom[j], aName);
                tmpLeafOfAtom[j] = i;
            }
        }
        int[] tmpAttachmentPoints = new int[tmpLeaves.size()];
        int[] tmpCleavedBondsOfLeaf = new int[tmpLeaves.size()];
        int[] tmpHeavyAtoms = new int[tmpLeaves.size()];
        int tmpCleavedBonds = 0;
        for (IBond tmpBond : tmpPreparedMolecule.bonds()) {
            int tmpBeginLeaf = tmpLeafOfAtom[tmpBond.getBegin().getIndex()];
            int tmpEndLeaf = tmpLeafOfAtom[tmpBond.getEnd().getIndex()];
            Assertions.assertTrue(tmpBeginLeaf >= 0 && tmpEndLeaf >= 0, aName);
            if (tmpBeginLeaf != tmpEndLeaf) {
                tmpCleavedBonds++;
                Assertions.assertFalse(tmpBond.isInRing(), aName);
                for (int tmpLeaf : new int[] {tmpBeginLeaf, tmpEndLeaf}) {
                    tmpAttachmentPoints[tmpLeaf]++;
                    tmpCleavedBondsOfLeaf[tmpLeaf]++;
                }
            }
        }
        int tmpMoleculeHeavyAtoms = 0;
        for (IAtom tmpAtom : tmpPreparedMolecule.atoms()) {
            int tmpLeaf = tmpLeafOfAtom[tmpAtom.getIndex()];
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAttachmentPoints[tmpLeaf]++;
            } else if (tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() > 1) {
                tmpHeavyAtoms[tmpLeaf]++;
                tmpMoleculeHeavyAtoms++;
            }
        }
        int tmpMinimumSize = aFragmenter.getMinimumFragmentSizeSetting();
        int tmpMinimumLinkerSize = Math.max(tmpMinimumSize, aFragmenter.getMinimumLinkerSizeSetting());
        for (int i = 0; i < tmpLeaves.size(); i++) {
            if (tmpCleavedBondsOfLeaf[i] == 0) {
                continue;
            }
            Assertions.assertTrue(tmpHeavyAtoms[i] >= (tmpAttachmentPoints[i] >= 2 ? tmpMinimumLinkerSize : tmpMinimumSize),
                    String.format("%s: fragment with %d heavy atoms and %d attachment points", aName, tmpHeavyAtoms[i],
                            tmpAttachmentPoints[i]));
        }
        int tmpFragmentHeavyAtoms = 0;
        int tmpLinkAtoms = 0;
        for (IAtomContainer tmpFragment : tmpFragments) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Integer tmpRuleNumber = tmpAtom.getProperty(CRUSHFragmenter.RULE_NUMBER_PROPERTY_KEY);
                if (tmpRuleNumber != null) {
                    tmpLinkAtoms++;
                    Assertions.assertTrue(tmpRuleNumber >= 1 && tmpRuleNumber <= 33, aName);
                    Assertions.assertEquals(tmpRuleNumber, tmpAtom.getMassNumber(), aName);
                } else if (!(tmpAtom instanceof IPseudoAtom) && tmpAtom.getAtomicNumber() != null
                        && tmpAtom.getAtomicNumber() > 1) {
                    tmpFragmentHeavyAtoms++;
                }
            }
        }
        Assertions.assertEquals(tmpMoleculeHeavyAtoms, tmpFragmentHeavyAtoms, aName);
        boolean tmpIsHydrogenSaturation = aFragmenter.getSaturationSetting()
                == CRUSHFragmenter.SaturationOption.HYDROGEN_SATURATION;
        Assertions.assertEquals(tmpIsHydrogenSaturation ? 0 : 2 * tmpCleavedBonds, tmpLinkAtoms, aName);
        return tmpCleavedBonds;
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
    static List<String> createUniqueSmilesList(List<IAtomContainer> aFragmentList) {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(ChemUtil.createUniqueSmiles(tmpFragment, true));
        }
        tmpSmilesList.sort(Comparator.nullsFirst(Comparator.naturalOrder()));
        return tmpSmilesList;
    }
    //
    /**
     * Describes the given bonds as sorted pairs of atom indices, separated by spaces, as the reference files do.
     *
     * @param aMolecule molecule
     * @param aBondIndexArray bond indices
     * @return description
     */
    private static String describeBonds(IAtomContainer aMolecule, int[] aBondIndexArray) {
        List<int[]> tmpPairs = new ArrayList<>(aBondIndexArray.length);
        for (int tmpBondIndex : aBondIndexArray) {
            IBond tmpBond = aMolecule.getBond(tmpBondIndex);
            int tmpBegin = tmpBond.getBegin().getIndex();
            int tmpEnd = tmpBond.getEnd().getIndex();
            tmpPairs.add(new int[] {Math.min(tmpBegin, tmpEnd), Math.max(tmpBegin, tmpEnd)});
        }
        tmpPairs.sort(Arrays::compare);
        List<String> tmpDescriptions = new ArrayList<>(tmpPairs.size());
        for (int[] tmpPair : tmpPairs) {
            tmpDescriptions.add(tmpPair[0] + "-" + tmpPair[1]);
        }
        return String.join(" ", tmpDescriptions);
    }
    //
    /**
     * Joins the indices of the given set with spaces.
     *
     * @param aBitSet set of indices
     * @return text
     */
    private static String join(BitSet aBitSet) {
        List<String> tmpIndices = new ArrayList<>(aBitSet.cardinality());
        for (int i = aBitSet.nextSetBit(0); i >= 0; i = aBitSet.nextSetBit(i + 1)) {
            tmpIndices.add(String.valueOf(i));
        }
        return String.join(" ", tmpIndices);
    }
    //
    /**
     * Returns the elements of the first list that are not contained in the second list (as multisets).
     *
     * @param aList1 first list
     * @param aList2 second list
     * @return difference
     */
    private static List<String> difference(List<String> aList1, List<String> aList2) {
        List<String> tmpDifference = new ArrayList<>(aList1);
        for (String tmpElement : aList2) {
            tmpDifference.remove(tmpElement);
        }
        return tmpDifference;
    }
    //
    /**
     * Appends the given text to the report if the maximum number of reported mismatches is not reached.
     *
     * @param aReport report
     * @param aCountArray counts of the step comparison; the number of reported cases is the sum of the skipped and
     *                    mismatching molecules
     * @param aMaximumReportedCount maximum number of reported mismatches
     * @param aText text to append
     */
    private static void report(StringBuilder aReport, int[] aCountArray, int aMaximumReportedCount, String aText) {
        int tmpCount = aCountArray[1] + aCountArray[2] + aCountArray[3] + aCountArray[4] + aCountArray[6] + aCountArray[7];
        if (tmpCount <= aMaximumReportedCount) {
            aReport.append(aText);
        }
    }
    //
    /**
     * Appends the given text to the report if the maximum number of reported mismatches is not reached.
     *
     * @param aReport report
     * @param aCount number of the mismatch
     * @param aText text to append
     */
    private static void report(StringBuilder aReport, int aCount, String aText) {
        if (aCount <= CRUSHFragmenterReferenceTest.MAXIMUM_REPORTED_MISMATCHES) {
            aReport.append(aText);
        }
    }
}
