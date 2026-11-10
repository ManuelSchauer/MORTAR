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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

/**
 * Cross-validates the {@link RECAPFragmenter} against an independent reference implementation of the same RECAP rules
 * and bond selection that uses the substructure matcher of the RDKit (script RECAP_Reference_Generator.py in the test
 * resources). The comparison tests the SMARTS semantics of the CDK, the assignment of the bond types, the avoidance of
 * small terminal fragments, and the simultaneous cleavage. The fragments are compared once with bond types and without
 * stereochemistry, and once with stereochemistry as MORTAR identifies fragments, i.e. with all link atoms written as
 * '*'; for the latter, stereo elements of atoms that are not stereogenic according to the CDK are removed on both
 * sides, because the RDKit keeps them. In addition, the data set is used to check that the leaves of the RECAP
 * hierarchy are the single-pass fragments and that the fragments do not depend on the atom order of the input.
 *
 * <p>The reference fragments of the BIOFACQUIM data set (natural products from Mexico, CC BY 4.0, see
 * Tutorial/BIOFACQUIM-description-licence.txt) are part of the test resources. A further reference file, e.g. for a
 * sample of ChEMBL, can be given with the environment variable {@value #REFERENCE_FILE_ENVIRONMENT_VARIABLE}.</p>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class RECAPFragmenterReferenceTest {
    /**
     * Name of the reference file of the BIOFACQUIM data set in the test resources.
     */
    private static final String BIOFACQUIM_REFERENCE_FILE = "RECAP_Reference_Fragments_BIOFACQUIM.tsv";
    /**
     * Name of the reference file of the BIOFACQUIM data set with strict functional group rules in the test resources.
     */
    private static final String BIOFACQUIM_STRICT_REFERENCE_FILE = "RECAP_Reference_Fragments_BIOFACQUIM_strict.tsv";
    /**
     * Environment variable with the path of an additional reference file.
     */
    private static final String REFERENCE_FILE_ENVIRONMENT_VARIABLE = "MORTAR_RECAP_REFERENCE_FILE";
    /**
     * Environment variable with the path of an additional reference file with strict functional group rules.
     */
    private static final String STRICT_REFERENCE_FILE_ENVIRONMENT_VARIABLE = "MORTAR_RECAP_STRICT_REFERENCE_FILE";
    /**
     * Maximum number of reported mismatches.
     */
    private static final int MAXIMUM_REPORTED_MISMATCHES = 25;
    //
    /**
     * Result of the comparison of a reference file.
     *
     * @param molecules number of compared molecules
     * @param skipped number of molecules whose input or reference SMILES could not be parsed by the CDK
     * @param mismatches number of molecules with different fragments (without stereochemistry)
     * @param stereoMismatches number of molecules with different fragments when stereochemistry is compared
     * @param fragmentationMilliseconds time spent in the fragmentation of all molecules
     * @param report description of the first mismatches
     */
    private record ComparisonResult(int molecules, int skipped, int mismatches, int stereoMismatches,
                                    long fragmentationMilliseconds, String report) {
    }
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public RECAPFragmenterReferenceTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Compares the fragments of all molecules of the BIOFACQUIM data set with the reference fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimReferenceTest() throws Exception {
        this.assertBiofacquimReference(RECAPFragmenterReferenceTest.BIOFACQUIM_REFERENCE_FILE, false);
    }
    //
    /**
     * Compares the fragments of all molecules of the BIOFACQUIM data set with strict functional group rules with the
     * reference fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimStrictReferenceTest() throws Exception {
        this.assertBiofacquimReference(RECAPFragmenterReferenceTest.BIOFACQUIM_STRICT_REFERENCE_FILE, true);
    }
    //
    /**
     * Checks for all molecules of the BIOFACQUIM data set that the leaves of the RECAP hierarchy are the fragments of
     * the single-pass cleavage (with stereochemistry), and that all hierarchy fragments are returned with the
     * respective output option.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimHierarchyTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        RECAPFragmenter tmpHierarchyFragmenter = new RECAPFragmenter();
        tmpHierarchyFragmenter.setFragmentOutputSetting(RECAPFragmenter.FragmentOutputOption.ALL_HIERARCHY_FRAGMENTS);
        int tmpMolecules = 0;
        int tmpNodes = 0;
        int tmpMaximumNodes = 0;
        StringBuilder tmpMismatches = new StringBuilder();
        try (InputStream tmpStream = RECAPFragmenterReferenceTest.class.getResourceAsStream(
                RECAPFragmenterReferenceTest.BIOFACQUIM_REFERENCE_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t", -1);
                tmpMolecules++;
                IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(tmpColumns[1]);
                List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                RECAPFragmenter.RECAPHierarchyNode tmpRoot = tmpFragmenter.buildHierarchy(tmpMolecule);
                List<RECAPFragmenter.RECAPHierarchyNode> tmpDescendants = tmpRoot.getAllDescendants();
                tmpNodes += 1 + tmpDescendants.size();
                tmpMaximumNodes = Math.max(tmpMaximumNodes, 1 + tmpDescendants.size());
                if (tmpRoot.isLeaf()) {
                    //no cleavable bond, the single pass returns the molecule
                    Assertions.assertEquals(1, tmpFragments.size(), tmpColumns[1]);
                    continue;
                }
                List<String> tmpExpected = RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragments, true);
                List<String> tmpLeaves = RECAPFragmenterTest.toStereoSmilesOfFragments(
                        RECAPFragmenterTest.fragmentsOf(tmpRoot.getLeaves()), true);
                if (!tmpExpected.equals(tmpLeaves)) {
                    tmpMismatches.append(String.format("%s %s%n  single pass: %s%n  leaves:      %s%n", tmpColumns[0],
                            tmpColumns[1], tmpExpected, tmpLeaves));
                }
                Assertions.assertEquals(RECAPFragmenterTest.toCanonicalSmiles(RECAPFragmenterTest.fragmentsOf(tmpDescendants), false),
                        RECAPFragmenterTest.toCanonicalSmiles(tmpHierarchyFragmenter.fragmentMolecule(tmpMolecule), false),
                        tmpColumns[1]);
            }
        }
        System.out.printf("%d molecules, %d hierarchy nodes, at most %d per molecule%n", tmpMolecules, tmpNodes,
                tmpMaximumNodes);
        Assertions.assertTrue(tmpMolecules > 550);
        Assertions.assertEquals("", tmpMismatches.toString());
    }
    //
    /**
     * Checks for all molecules of the BIOFACQUIM data set that the fragments do not depend on the order of the atoms
     * and bonds of the input, with default and with strict functional group rules: every molecule is fragmented once as
     * parsed from the input SMILES and once as parsed from its unique SMILES created by the CDK, whose atom order
     * differs, and the fragments are compared with bond types and without stereochemistry.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void biofacquimAtomOrderTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        RECAPFragmenter tmpStrictFragmenter = new RECAPFragmenter();
        tmpStrictFragmenter.setStrictFunctionalGroupsSetting(true);
        int tmpMolecules = 0;
        int tmpReorderedMolecules = 0;
        StringBuilder tmpMismatches = new StringBuilder();
        try (InputStream tmpStream = RECAPFragmenterReferenceTest.class.getResourceAsStream(
                RECAPFragmenterReferenceTest.BIOFACQUIM_REFERENCE_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t", -1);
                tmpMolecules++;
                IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(tmpColumns[1]);
                //the output order of the atoms, i.e. their indices in the reordered molecule
                int[] tmpOrder = new int[tmpMolecule.getAtomCount()];
                String tmpReorderedSmiles = SmilesGenerator.create(tmpMolecule, SmiFlavor.Unique | SmiFlavor.Stereo,
                        tmpOrder);
                IAtomContainer tmpReorderedMolecule = RECAPFragmenterTest.parse(tmpReorderedSmiles);
                for (int i = 0; i < tmpOrder.length; i++) {
                    if (tmpOrder[i] != i) {
                        tmpReorderedMolecules++;
                        break;
                    }
                }
                for (RECAPFragmenter tmpCurrent : new RECAPFragmenter[] {tmpFragmenter, tmpStrictFragmenter}) {
                    List<String> tmpFragments = RECAPFragmenterTest.toCanonicalSmiles(
                            tmpCurrent.fragmentMolecule(tmpMolecule), false);
                    List<String> tmpReorderedFragments = RECAPFragmenterTest.toCanonicalSmiles(
                            tmpCurrent.fragmentMolecule(tmpReorderedMolecule), false);
                    if (!tmpFragments.equals(tmpReorderedFragments)) {
                        tmpMismatches.append(String.format(
                                "%s %s (strict %s)%n  input order: %s%n  reordered %s: %s%n", tmpColumns[0],
                                tmpColumns[1], tmpCurrent.getStrictFunctionalGroupsSetting(), tmpFragments,
                                tmpReorderedSmiles, tmpReorderedFragments));
                    }
                }
            }
        }
        Assertions.assertTrue(tmpMolecules > 550);
        //the canonical atom order of the CDK differs from the atom order of the input SMILES, written by the RDKit
        Assertions.assertTrue(tmpReorderedMolecules > tmpMolecules / 2);
        Assertions.assertEquals("", tmpMismatches.toString());
    }
    //
    /**
     * Compares the fragments of all molecules of the BIOFACQUIM data set with the given reference file, which must match
     * exactly.
     *
     * @param aFileName name of the reference file in the test resources
     * @param isStrict whether the strict functional group rules are applied
     * @throws Exception if anything goes wrong
     */
    private void assertBiofacquimReference(String aFileName, boolean isStrict) throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        tmpFragmenter.setStrictFunctionalGroupsSetting(isStrict);
        ComparisonResult tmpResult;
        try (InputStream tmpStream = RECAPFragmenterReferenceTest.class.getResourceAsStream(aFileName)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            tmpResult = this.compare(new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8)),
                    tmpFragmenter);
        }
        System.out.println(tmpResult);
        Assertions.assertTrue(tmpResult.molecules() > 550);
        Assertions.assertEquals(0, tmpResult.skipped(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.mismatches(), tmpResult.report());
        Assertions.assertEquals(0, tmpResult.stereoMismatches(), tmpResult.report());
    }
    //
    /**
     * Compares the fragments with the reference file given by the environment variable
     * {@value #REFERENCE_FILE_ENVIRONMENT_VARIABLE}; skipped if it is not set. Large and diverse data sets contain a
     * few molecules for which the toolkits differ (e.g. in the aromaticity perception), so up to one mismatch per
     * thousand molecules (with stereochemistry two) is tolerated; the mismatches are printed. On a sample of every
     * hundredth molecule of ChEMBL 37 (28,978 molecules), 6 molecules have different fragments (tellurium-containing
     * and imide rings that only one of the toolkits perceives as aromatic, and sulfoximine and phosphazene rings with
     * different Kekulé structures) and 11 further molecules have different stereochemistry (e.g. stereo double bonds
     * in small rings that only the CDK keeps, and sulfoxides).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void additionalReferenceTest() throws Exception {
        this.assertAdditionalReference(RECAPFragmenterReferenceTest.REFERENCE_FILE_ENVIRONMENT_VARIABLE, false);
    }
    //
    /**
     * Compares the fragments with strict functional group rules with the reference file given by the environment
     * variable {@value #STRICT_REFERENCE_FILE_ENVIRONMENT_VARIABLE} (generated with the option --strict); skipped if it
     * is not set. The same tolerances as in {@link #additionalReferenceTest()} apply.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void additionalStrictReferenceTest() throws Exception {
        this.assertAdditionalReference(RECAPFragmenterReferenceTest.STRICT_REFERENCE_FILE_ENVIRONMENT_VARIABLE, true);
    }
    //
    /**
     * Compares the fragments with the reference file given by the given environment variable, with the tolerances of
     * {@link #additionalReferenceTest()}; skipped if the variable is not set.
     *
     * @param anEnvironmentVariable name of the environment variable with the path of the reference file
     * @param isStrict whether the strict functional group rules are applied
     * @throws Exception if anything goes wrong
     */
    private void assertAdditionalReference(String anEnvironmentVariable, boolean isStrict) throws Exception {
        String tmpPath = System.getenv(anEnvironmentVariable);
        Assumptions.assumeTrue(tmpPath != null && !tmpPath.isBlank(), "No additional reference file given.");
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        tmpFragmenter.setStrictFunctionalGroupsSetting(isStrict);
        ComparisonResult tmpResult;
        try (BufferedReader tmpReader = Files.newBufferedReader(Path.of(tmpPath), StandardCharsets.UTF_8)) {
            tmpResult = this.compare(tmpReader, tmpFragmenter);
        }
        System.out.println(tmpResult);
        Assertions.assertTrue(tmpResult.mismatches() <= tmpResult.molecules() / 1000, tmpResult.report());
        Assertions.assertTrue(tmpResult.stereoMismatches() <= tmpResult.molecules() / 500, tmpResult.report());
    }
    //
    /**
     * Fragments the input molecules of the given reference file and compares the fragments to the reference fragments.
     *
     * @param aReader reader of the reference file (tab-separated: name, input SMILES, fragments separated by '.',
     *                cleaved bonds)
     * @param aFragmenter fragmenter with the settings of the reference file
     * @return comparison result
     * @throws IOException if the file cannot be read
     */
    private ComparisonResult compare(BufferedReader aReader, RECAPFragmenter aFragmenter) throws IOException {
        RECAPFragmenter tmpFragmenter = aFragmenter;
        int tmpMolecules = 0;
        int tmpSkipped = 0;
        int tmpMismatches = 0;
        int tmpStereoMismatches = 0;
        long tmpFragmentationNanoseconds = 0;
        StringBuilder tmpReport = new StringBuilder();
        String tmpLine;
        while ((tmpLine = aReader.readLine()) != null) {
            if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                continue;
            }
            String[] tmpColumns = tmpLine.split("\t", -1);
            tmpMolecules++;
            List<String> tmpExpected;
            List<String> tmpExpectedStereo;
            List<IAtomContainer> tmpFragments;
            List<String> tmpActual;
            List<String> tmpActualStereo;
            try {
                List<String> tmpReferenceFragments = Arrays.asList(tmpColumns[2].split("\\."));
                tmpExpected = RECAPFragmenterTest.canonical(tmpReferenceFragments, false);
                tmpExpectedStereo = RECAPFragmenterTest.toStereoSmiles(tmpReferenceFragments, true);
                IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(tmpColumns[1]);
                long tmpStart = System.nanoTime();
                tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                tmpFragmentationNanoseconds += System.nanoTime() - tmpStart;
                tmpActual = RECAPFragmenterTest.toCanonicalSmiles(tmpFragments, false);
                tmpActualStereo = RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragments, true);
            } catch (CDKException | CloneNotSupportedException | IllegalArgumentException anException) {
                tmpSkipped++;
                this.report(tmpReport, tmpSkipped + tmpMismatches + tmpStereoMismatches,
                        String.format("%s: skipped, %s%n", tmpColumns[0], anException));
                continue;
            }
            this.checkAtomConservation(tmpColumns[1], tmpFragments, tmpExpected.equals(tmpActual) ? tmpColumns[3] : null);
            if (!tmpExpected.equals(tmpActual)) {
                tmpMismatches++;
                this.report(tmpReport, tmpSkipped + tmpMismatches + tmpStereoMismatches,
                        String.format("%s %s%n  reference bonds: %s%n  actual bonds:    %s%n  expected: %s%n  actual:   %s%n",
                                tmpColumns[0], tmpColumns[1], tmpColumns[3], this.describeCleavedBonds(tmpFragmenter, tmpColumns[1]),
                                tmpExpected, tmpActual));
            } else if (!tmpExpectedStereo.equals(tmpActualStereo)) {
                tmpStereoMismatches++;
                this.report(tmpReport, tmpSkipped + tmpMismatches + tmpStereoMismatches,
                        String.format("%s %s (stereo)%n  expected: %s%n  actual:   %s%n",
                                tmpColumns[0], tmpColumns[1], tmpExpectedStereo, tmpActualStereo));
            }
        }
        return new ComparisonResult(tmpMolecules, tmpSkipped, tmpMismatches, tmpStereoMismatches,
                tmpFragmentationNanoseconds / 1_000_000L, tmpReport.toString());
    }
    //
    /**
     * Checks that every heavy atom of the molecule is contained in exactly one fragment and, if the reference bonds are
     * given (i.e. the fragments match the reference), that two link atoms are created per cleaved bond.
     *
     * @param aSmiles SMILES of the molecule
     * @param aFragmentList fragments
     * @param aReferenceBonds cleaved bonds of the reference, separated by spaces, or null if not checked
     */
    private void checkAtomConservation(String aSmiles, List<IAtomContainer> aFragmentList, String aReferenceBonds) {
        int tmpHeavyAtoms = 0;
        int tmpLinkAtoms = 0;
        for (IAtomContainer tmpFragment : aFragmentList) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom && tmpAtom.getProperty(RECAPFragmenter.BOND_TYPE_PROPERTY_KEY) != null) {
                    tmpLinkAtoms++;
                } else if (tmpAtom.getAtomicNumber() != 1) {
                    tmpHeavyAtoms++;
                }
            }
        }
        try {
            IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(aSmiles);
            int tmpExpectedHeavyAtoms = 0;
            for (IAtom tmpAtom : tmpMolecule.atoms()) {
                if (tmpAtom.getAtomicNumber() == null || tmpAtom.getAtomicNumber() != 1) {
                    tmpExpectedHeavyAtoms++;
                }
            }
            Assertions.assertEquals(tmpExpectedHeavyAtoms, tmpHeavyAtoms, aSmiles);
        } catch (CDKException anException) {
            Assertions.fail(anException);
        }
        if (aReferenceBonds != null) {
            int tmpCleavedBonds = aReferenceBonds.isBlank() ? 0 : aReferenceBonds.trim().split(" ").length;
            Assertions.assertEquals(2 * tmpCleavedBonds, tmpLinkAtoms, aSmiles);
        }
    }
    //
    /**
     * Describes the bonds cleaved by the fragmenter as bond type and atom indices, for the mismatch report.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @return description
     */
    private String describeCleavedBonds(RECAPFragmenter aFragmenter, String aSmiles) {
        try {
            IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(aSmiles);
            aFragmenter.prepareMolecule(tmpMolecule);
            List<String> tmpDescriptions = new ArrayList<>();
            for (RECAPFragmenter.RECAPBond tmpBond : aFragmenter.selectBondsToCleave(tmpMolecule,
                    aFragmenter.findRECAPBonds(tmpMolecule))) {
                tmpDescriptions.add(String.format("%d:%d-%d", tmpBond.bondType(), tmpBond.bond().getBegin().getIndex(),
                        tmpBond.bond().getEnd().getIndex()));
            }
            return String.join(" ", tmpDescriptions);
        } catch (CDKException anException) {
            return anException.toString();
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
    private void report(StringBuilder aReport, int aCount, String aText) {
        if (aCount <= RECAPFragmenterReferenceTest.MAXIMUM_REPORTED_MISMATCHES) {
            aReport.append(aText);
        }
    }
}
