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
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

/**
 * Tests the RECAP rules of the {@link RECAPFragmenter} with the test molecules of the individual rules of the RECAP
 * re-implementation for the CDK by Jonas Schaub (file RECAP_Rule_Cases.tsv in the test resources). For every molecule,
 * the file gives the RECAP bond type of the tested rule, whether that implementation's rule matches, whether a bond of
 * this type is cleaved by the RECAPFragmenter with default and with strict functional group settings, and a category
 * that explains the differences: e.g. terminal groups that the CDK implementation excludes in its SMARTS are kept by
 * the small-group rule here, pseudo atoms of the input count as attachment points here, aryl-nitrogen and aryl-oxygen
 * bonds are not cleaved (publication), and most of the extended functional groups that the CDK implementation excludes
 * are excluded by the strict setting.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class RECAPFragmenterRuleCaseTest {
    /**
     * Name of the file with the rule cases in the test resources.
     */
    private static final String RULE_CASES_FILE = "RECAP_Rule_Cases.tsv";
    /**
     * Number of rule cases in the file.
     */
    private static final int RULE_CASES_COUNT = 279;
    //
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public RECAPFragmenterRuleCaseTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Checks for every rule case whether a bond of the given bond type is cleaved with default and with strict
     * functional group settings, as expected.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleCasesTest() throws Exception {
        RECAPFragmenter tmpDefaultFragmenter = new RECAPFragmenter();
        RECAPFragmenter tmpStrictFragmenter = new RECAPFragmenter();
        tmpStrictFragmenter.setStrictFunctionalGroupsSetting(true);
        int tmpCases = 0;
        StringBuilder tmpMismatches = new StringBuilder();
        try (InputStream tmpStream = RECAPFragmenterRuleCaseTest.class.getResourceAsStream(
                RECAPFragmenterRuleCaseTest.RULE_CASES_FILE)) {
            Objects.requireNonNull(tmpStream, "Rule cases file not found.");
            BufferedReader tmpReader = new BufferedReader(new InputStreamReader(tmpStream, StandardCharsets.UTF_8));
            String tmpLine;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.isBlank() || tmpLine.startsWith("#")) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t", -1);
                tmpCases++;
                String tmpSmiles = tmpColumns[1];
                int tmpBondType = Integer.parseInt(tmpColumns[2]);
                boolean tmpIsCleavedDefault = this.isCleaved(tmpDefaultFragmenter, tmpSmiles, tmpBondType);
                boolean tmpIsCleavedStrict = this.isCleaved(tmpStrictFragmenter, tmpSmiles, tmpBondType);
                if (tmpIsCleavedDefault != this.parseExpectation(tmpColumns[4])
                        || tmpIsCleavedStrict != this.parseExpectation(tmpColumns[5])) {
                    tmpMismatches.append(String.format("%s %s (bond type %d, %s): expected default %s, strict %s; actual "
                                    + "default %s, strict %s%n", tmpColumns[0], tmpSmiles, tmpBondType, tmpColumns[6],
                            tmpColumns[4], tmpColumns[5], tmpIsCleavedDefault, tmpIsCleavedStrict));
                }
            }
        }
        Assertions.assertEquals(RECAPFragmenterRuleCaseTest.RULE_CASES_COUNT, tmpCases);
        Assertions.assertEquals("", tmpMismatches.toString());
    }
    //
    /**
     * Returns whether the given fragmenter cleaves a bond of the given RECAP bond type in the given molecule.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param aBondType RECAP bond type
     * @return true if a bond of this type is cleaved
     * @throws Exception if anything goes wrong
     */
    private boolean isCleaved(RECAPFragmenter aFragmenter, String aSmiles, int aBondType) throws Exception {
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(aSmiles);
        aFragmenter.prepareMolecule(tmpMolecule);
        List<RECAPFragmenter.RECAPBond> tmpCleavedBonds = aFragmenter.selectBondsToCleave(tmpMolecule,
                aFragmenter.findRECAPBonds(tmpMolecule));
        for (RECAPFragmenter.RECAPBond tmpBond : tmpCleavedBonds) {
            if (tmpBond.bondType() == aBondType) {
                return true;
            }
        }
        return false;
    }
    //
    /**
     * Parses an expectation of the rule cases file.
     *
     * @param anExpectation "cleaved" or "not cleaved"
     * @return true for "cleaved"
     */
    private boolean parseExpectation(String anExpectation) {
        return switch (anExpectation) {
            case "cleaved" -> true;
            case "not cleaved" -> false;
            default -> throw new IllegalArgumentException("Unknown expectation: " + anExpectation);
        };
    }
}
