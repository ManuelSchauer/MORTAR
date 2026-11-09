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
 * Cross-validates the {@link BRICSFragmenter} against an independent reference implementation of the original BRICS
 * rules that uses the substructure matcher of the RDKit (script BRICS_Reference_Generator.py in the test resources).
 * Both implementations use the same published SMARTS, so the comparison tests the SMARTS semantics of the CDK, the
 * rule order, the simultaneous cleavage, and the assignment of the link types. The fragments are compared once with
 * link types and without stereochemistry, and once with stereochemistry as MORTAR identifies fragments, i.e. with
 * all link atoms written as '*'; for the latter, stereo elements of atoms that are not stereogenic according to the
 * CDK are removed on both sides, because the RDKit keeps them.
 *
 * <p>The reference fragments of the BIOFACQUIM data set (natural products from Mexico, CC BY 4.0, see
 * Tutorial/BIOFACQUIM-description-licence.txt) are part of the test resources. A further reference file, e.g. for a
 * sample of ChEMBL, can be given with the environment variable {@value #REFERENCE_FILE_ENVIRONMENT_VARIABLE}.</p>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class BRICSFragmenterReferenceTest {
    /**
     * Name of the reference file of the BIOFACQUIM data set in the test resources.
     */
    private static final String BIOFACQUIM_REFERENCE_FILE = "BRICS_Reference_Fragments_BIOFACQUIM.tsv";
    /**
     * Environment variable with the path of an additional reference file.
     */
    private static final String REFERENCE_FILE_ENVIRONMENT_VARIABLE = "MORTAR_BRICS_REFERENCE_FILE";
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
    public BRICSFragmenterReferenceTest() {
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
        ComparisonResult tmpResult;
        try (InputStream tmpStream = BRICSFragmenterReferenceTest.class.getResourceAsStream(
                BRICSFragmenterReferenceTest.BIOFACQUIM_REFERENCE_FILE)) {
            Objects.requireNonNull(tmpStream, "Reference file not found.");
            tmpResult = this.compare(new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8)));
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
     * few molecules for which the toolkits differ, so up to one mismatch per thousand molecules (with stereochemistry
     * two) is tolerated; the mismatches are printed. On a sample of every hundredth molecule of ChEMBL 37 (28,978
     * molecules), 5 molecules have different fragments (tellurium-containing rings that only the RDKit perceives as
     * aromatic, and sulfoximine and phosphazene rings with different Kekulé structures) and 17 further molecules have
     * different stereochemistry (e.g. stereo double bonds in small rings that only the CDK keeps, and sulfoxides).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void additionalReferenceTest() throws Exception {
        String tmpPath = System.getenv(BRICSFragmenterReferenceTest.REFERENCE_FILE_ENVIRONMENT_VARIABLE);
        Assumptions.assumeTrue(tmpPath != null && !tmpPath.isBlank(), "No additional reference file given.");
        ComparisonResult tmpResult;
        try (BufferedReader tmpReader = Files.newBufferedReader(Path.of(tmpPath), StandardCharsets.UTF_8)) {
            tmpResult = this.compare(tmpReader);
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
     * @return comparison result
     * @throws IOException if the file cannot be read
     */
    private ComparisonResult compare(BufferedReader aReader) throws IOException {
        BRICSFragmenter tmpFragmenter = new BRICSFragmenter();
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
                tmpExpected = BRICSFragmenterTest.canonical(tmpReferenceFragments, false);
                tmpExpectedStereo = BRICSFragmenterTest.toStereoSmiles(tmpReferenceFragments, true);
                IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(tmpColumns[1]);
                long tmpStart = System.nanoTime();
                tmpFragments = tmpFragmenter.fragmentMolecule(tmpMolecule);
                tmpFragmentationNanoseconds += System.nanoTime() - tmpStart;
                tmpActual = BRICSFragmenterTest.toCanonicalSmiles(tmpFragments, false);
                tmpActualStereo = BRICSFragmenterTest.toStereoSmilesOfFragments(tmpFragments, true);
            } catch (CDKException | CloneNotSupportedException | IllegalArgumentException anException) {
                tmpSkipped++;
                this.report(tmpReport, tmpSkipped + tmpMismatches + tmpStereoMismatches,
                        String.format("%s: skipped, %s%n", tmpColumns[0], anException));
                continue;
            }
            this.checkAtomConservation(tmpColumns[1], tmpFragments, tmpColumns[3]);
            if (!tmpExpected.equals(tmpActual)) {
                tmpMismatches++;
                this.report(tmpReport, tmpSkipped + tmpMismatches + tmpStereoMismatches,
                        String.format("%s %s%n  reference bonds: %s%n  actual bonds:    %s%n  expected: %s%n  actual:   %s%n",
                                tmpColumns[0], tmpColumns[1], tmpColumns[3], this.describeBRICSBonds(tmpFragmenter, tmpColumns[1]),
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
     * Checks that every heavy atom of the molecule is contained in exactly one fragment and that two link atoms are
     * created per cleaved bond.
     *
     * @param aSmiles SMILES of the molecule
     * @param aFragmentList fragments
     * @param aReferenceBonds cleaved bonds of the reference, separated by spaces
     */
    private void checkAtomConservation(String aSmiles, List<IAtomContainer> aFragmentList, String aReferenceBonds) {
        int tmpHeavyAtoms = 0;
        int tmpLinkAtoms = 0;
        for (IAtomContainer tmpFragment : aFragmentList) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom && tmpAtom.getProperty(BRICSFragmenter.LINK_TYPE_PROPERTY_KEY) != null) {
                    tmpLinkAtoms++;
                } else if (tmpAtom.getAtomicNumber() != 1) {
                    tmpHeavyAtoms++;
                }
            }
        }
        try {
            IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(aSmiles);
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
        int tmpCleavedBonds = aReferenceBonds.isBlank() ? 0 : aReferenceBonds.trim().split(" ").length;
        Assertions.assertEquals(2 * tmpCleavedBonds, tmpLinkAtoms, aSmiles);
    }
    //
    /**
     * Describes the bonds cleaved by the fragmenter as rule name and atom indices, for the mismatch report.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @return description
     */
    private String describeBRICSBonds(BRICSFragmenter aFragmenter, String aSmiles) {
        try {
            IAtomContainer tmpMolecule = BRICSFragmenterTest.parse(aSmiles);
            aFragmenter.prepareMolecule(tmpMolecule);
            List<String> tmpDescriptions = new ArrayList<>();
            for (BRICSFragmenter.BRICSBond tmpBond : aFragmenter.findBRICSBonds(tmpMolecule)) {
                tmpDescriptions.add(String.format("%s:%d-%d", tmpBond.ruleName(), tmpBond.bond().getBegin().getIndex(),
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
        if (aCount <= BRICSFragmenterReferenceTest.MAXIMUM_REPORTED_MISMATCHES) {
            aReport.append(aText);
        }
    }
}
