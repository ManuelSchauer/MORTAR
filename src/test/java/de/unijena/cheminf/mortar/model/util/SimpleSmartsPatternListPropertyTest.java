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

package de.unijena.cheminf.mortar.model.util;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Locale;

/**
 * Test class for the {@link SimpleSmartsPatternListProperty} and its line format.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SimpleSmartsPatternListPropertyTest {
    /**
     * Constructor setting the default locale, which is needed for the language-specific error messages.
     */
    public SimpleSmartsPatternListPropertyTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests that comments, blank lines, and surrounding whitespace are ignored and that the name is the rest of the
     * line after the first whitespace, including inner whitespace.
     */
    @Test
    public void parseEntriesTest() {
        String tmpText = String.join("\n",
                "# comment line with a SMARTS-like C(=O)O",
                "",
                "   ",
                "C(=O)[OH] carboxylic acid",
                "  [#7;!$(N=O)]   amine nitrogen   ",
                "\t# indented comment",
                "c1ccccc1",
                "[OX2H]\thydroxyl\tgroup",
                "C#N nitrile # not a comment");
        List<SimpleSmartsPatternListProperty.SmartsPatternEntry> tmpEntries = SimpleSmartsPatternListProperty.parseEntries(tmpText);
        Assertions.assertEquals(List.of(
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("C(=O)[OH]", "carboxylic acid"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[#7;!$(N=O)]", "amine nitrogen"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("c1ccccc1", ""),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[OX2H]", "hydroxyl\tgroup"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("C#N", "nitrile # not a comment")),
                tmpEntries);
        Assertions.assertEquals("carboxylic acid", tmpEntries.get(0).getLabel());
        Assertions.assertFalse(tmpEntries.get(2).hasName());
        Assertions.assertEquals("c1ccccc1", tmpEntries.get(2).getLabel());
        Assertions.assertTrue(SimpleSmartsPatternListProperty.validate(tmpText).isEmpty());
        Assertions.assertEquals(5, SimpleSmartsPatternListProperty.countEntries(tmpText));
        //line numbers refer to the given text, starting at 1
        List<SimpleSmartsPatternListProperty.ParsedLine> tmpParsedLines = SimpleSmartsPatternListProperty.parseLines(tmpText);
        Assertions.assertEquals(List.of(4, 5, 7, 8, 9), tmpParsedLines.stream().map(SimpleSmartsPatternListProperty.ParsedLine::lineNumber).toList());
        Assertions.assertEquals("[#7;!$(N=O)]   amine nitrogen", tmpParsedLines.get(1).line());
    }
    //
    /**
     * Tests that Windows and old Mac line breaks and a UTF-8 byte order mark at the beginning are handled.
     */
    @Test
    public void lineBreaksAndByteOrderMarkTest() {
        List<SimpleSmartsPatternListProperty.SmartsPatternEntry> tmpExpected = List.of(
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[OH]", "hydroxyl"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("C=O", ""),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[NH2]", "amino"));
        Assertions.assertEquals(tmpExpected, SimpleSmartsPatternListProperty.parseEntries("[OH] hydroxyl\r\nC=O\r\n\r\n[NH2] amino\r\n"));
        Assertions.assertEquals(tmpExpected, SimpleSmartsPatternListProperty.parseEntries("[OH] hydroxyl\rC=O\r[NH2] amino"));
        Assertions.assertEquals(tmpExpected, SimpleSmartsPatternListProperty.parseEntries("\uFEFF[OH] hydroxyl\nC=O\n[NH2] amino\n"));
        Assertions.assertEquals(3, SimpleSmartsPatternListProperty.countEntries("\uFEFF[OH] hydroxyl\r\nC=O\r\n\r\n[NH2] amino\r\n"));
        //line numbers count empty lines, also with Windows line breaks
        List<SimpleSmartsPatternListProperty.ParsedLine> tmpInvalidLines =
                SimpleSmartsPatternListProperty.validate("C\r\n\r\nC(C broken\r\nO");
        Assertions.assertEquals(1, tmpInvalidLines.size());
        Assertions.assertEquals(3, tmpInvalidLines.getFirst().lineNumber());
    }
    //
    /**
     * Tests that invalid lines are reported with line number, content, and the message of the CDK SMARTS parser, while
     * the valid lines are still parsed.
     */
    @Test
    public void invalidLinesTest() {
        String tmpText = String.join("\n",
                "C(=O)O acid",
                "C(C unclosed branch",
                "# comment",
                "[C;H2 unclosed bracket",
                "c1ccc ring not closed",
                "[OH] hydroxyl");
        List<SimpleSmartsPatternListProperty.ParsedLine> tmpInvalidLines = SimpleSmartsPatternListProperty.validate(tmpText);
        Assertions.assertEquals(List.of(2, 4, 5), tmpInvalidLines.stream().map(SimpleSmartsPatternListProperty.ParsedLine::lineNumber).toList());
        for (SimpleSmartsPatternListProperty.ParsedLine tmpInvalidLine : tmpInvalidLines) {
            Assertions.assertFalse(tmpInvalidLine.isValid());
            Assertions.assertNull(tmpInvalidLine.entry());
            Assertions.assertTrue(tmpInvalidLine.errorMessage().startsWith("Invalid SMARTS pattern (CDK SMARTS parser, position"),
                    tmpInvalidLine.errorMessage());
        }
        Assertions.assertEquals("C(C unclosed branch", tmpInvalidLines.getFirst().line());
        Assertions.assertEquals(List.of(
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("C(=O)O", "acid"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[OH]", "hydroxyl")),
                SimpleSmartsPatternListProperty.parseEntries(tmpText));
        String tmpDescription = SimpleSmartsPatternListProperty.describeInvalidLines(tmpInvalidLines);
        Assertions.assertTrue(tmpDescription.startsWith("Line 2 (C(C unclosed branch): Invalid SMARTS pattern"), tmpDescription);
        Assertions.assertEquals(3, tmpDescription.lines().count());
    }
    //
    /**
     * Tests the validation of single SMARTS patterns and names.
     */
    @Test
    public void validateSmartsAndNameTest() {
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateSmarts("[$(C(=O)[OH]),$(C(=O)[O-])]"));
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateSmarts("[OH].[OH]"));
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateSmarts("([OH].[OH])"));
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateSmarts("C/C=C/C"));
        Assertions.assertEquals("SMARTS pattern is missing", SimpleSmartsPatternListProperty.validateSmarts(""));
        Assertions.assertEquals("SMARTS pattern must not contain whitespace", SimpleSmartsPatternListProperty.validateSmarts("C C"));
        Assertions.assertEquals("SMARTS pattern must not contain whitespace", SimpleSmartsPatternListProperty.validateSmarts("C\tC"));
        Assertions.assertNotNull(SimpleSmartsPatternListProperty.validateSmarts("C)"));
        Assertions.assertNotNull(SimpleSmartsPatternListProperty.validateSmarts("[C"));
        Assertions.assertNotNull(SimpleSmartsPatternListProperty.validateSmarts("[$(C]"));
        Assertions.assertEquals("Reaction SMARTS patterns ('>') are not supported",
                SimpleSmartsPatternListProperty.validateSmarts("[C:1][OH]>>[C:1]=O"));
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateName(""));
        Assertions.assertNull(SimpleSmartsPatternListProperty.validateName("name with spaces, \u00e4\u00f6\u00fc and #"));
        Assertions.assertEquals("Name must not contain line breaks", SimpleSmartsPatternListProperty.validateName("two\nlines"));
        Assertions.assertEquals("Name must not contain line breaks", SimpleSmartsPatternListProperty.validateName("two\rlines"));
    }
    //
    /**
     * Tests that the validation does not throw for any incomplete SMARTS pattern, as it occurs while typing in the GUI;
     * the CDK SMARTS parser throws exceptions for some of them, e.g. "[".
     */
    @Test
    public void incompleteSmartsTest() {
        String[] tmpSmartsPatterns = new String[] {
                "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1-&!@[C@@H](O)/C=C/[13CH3]",
                "c%10ccccc%10.([OH]).([OH])",
                "[C;!$(C=O)]-[OX2H1]>>[C:1][O:2]",
                "[NX3;H2,H1;!$(NC=O)]-[CX4@@H](-[CH3])C(=O)[O;H1,-1]"};
        Assertions.assertNotNull(SimpleSmartsPatternListProperty.validateSmarts("["));
        for (String tmpSmarts : tmpSmartsPatterns) {
            for (int i = 1; i <= tmpSmarts.length(); i++) {
                String tmpPrefix = tmpSmarts.substring(0, i);
                Assertions.assertDoesNotThrow(() -> SimpleSmartsPatternListProperty.validateSmarts(tmpPrefix), tmpPrefix);
            }
        }
        Assertions.assertTrue(SimpleSmartsPatternListProperty.validate("[\n[C@\n[$(").size() == 3);
    }
    //
    /**
     * Tests that the property only accepts valid text and keeps its value otherwise, and that the text created from
     * entries is parsed to the same entries.
     */
    @Test
    public void propertyTest() {
        SimpleSmartsPatternListProperty tmpProperty = new SimpleSmartsPatternListProperty(this, "Test property");
        Assertions.assertEquals("", tmpProperty.get());
        Assertions.assertTrue(tmpProperty.getEntries().isEmpty());
        String tmpValidText = "# my patterns\nC(=O)O acid\n\n[OH]";
        tmpProperty.set(tmpValidText);
        Assertions.assertEquals(tmpValidText, tmpProperty.get());
        Assertions.assertEquals(2, tmpProperty.getEntries().size());
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpProperty.set("C(=O)O acid\nC(C"));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpProperty.setValue("[C"));
        Assertions.assertThrows(NullPointerException.class, () -> tmpProperty.set(null));
        Assertions.assertEquals(tmpValidText, tmpProperty.get());
        Assertions.assertThrows(IllegalArgumentException.class,
                () -> new SimpleSmartsPatternListProperty(this, "Test property", "C(C"));
        //text created from entries
        List<SimpleSmartsPatternListProperty.SmartsPatternEntry> tmpEntries = List.of(
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("[#6]~[#8]", "carbon oxygen pair"),
                new SimpleSmartsPatternListProperty.SmartsPatternEntry("c1ccccc1", ""));
        String tmpText = SimpleSmartsPatternListProperty.toText(tmpEntries);
        Assertions.assertEquals("[#6]~[#8] carbon oxygen pair\nc1ccccc1", tmpText);
        Assertions.assertEquals(tmpEntries, SimpleSmartsPatternListProperty.parseEntries(tmpText));
        Assertions.assertEquals("", SimpleSmartsPatternListProperty.toText(List.of()));
    }
}
