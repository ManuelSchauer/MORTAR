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

import de.unijena.cheminf.mortar.message.Message;

import javafx.beans.property.SimpleStringProperty;

import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smarts.SmartsResult;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * A JavaFX property for an ordered list of SMARTS patterns with optional names. The value is plain text with one
 * entry per line in the format {@code SMARTS[<whitespace>name]}, i.e. the SMARTS pattern, optionally followed by
 * whitespace and a name that may contain further whitespace. Blank lines and lines starting with
 * {@value #COMMENT_LINE_PREFIX} are ignored, leading and trailing whitespace of a line is ignored, and line breaks of
 * all platforms are accepted. This is also the format of the text files the patterns can be imported from, so that
 * file content, GUI state, and persisted setting share one representation.
 * <br>The value is a string on purpose: the settings dialogs store the value of every setting when they are opened
 * and restore it when the changes are cancelled, which requires an immutable value.
 * <br>Only text in which every non-comment line holds a valid SMARTS pattern can be set; otherwise, an
 * {@link IllegalArgumentException} is thrown. The SMARTS patterns are validated with the CDK SMARTS parser.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SimpleSmartsPatternListProperty extends SimpleStringProperty {
    //<editor-fold desc="Record SmartsPatternEntry">
    /**
     * One entry of the pattern list, i.e. a SMARTS pattern and its optional name.
     *
     * @param smarts the SMARTS pattern, not empty and without whitespace
     * @param name the name of the pattern, an empty string if the pattern is unnamed; never null
     */
    public record SmartsPatternEntry(String smarts, String name) {
        /**
         * Constructor that checks the parameters.
         *
         * @param smarts the SMARTS pattern
         * @param name the name of the pattern, an empty string if the pattern is unnamed
         * @throws NullPointerException if a parameter is null
         */
        public SmartsPatternEntry {
            Objects.requireNonNull(smarts, "Given SMARTS pattern is null.");
            Objects.requireNonNull(name, "Given name is null.");
        }
        //
        /**
         * Returns whether the pattern has a name.
         *
         * @return true if the name is not empty
         */
        public boolean hasName() {
            return !this.name.isEmpty();
        }
        //
        /**
         * Returns the label of the pattern, i.e. its name or, if it is unnamed, the SMARTS pattern itself.
         *
         * @return name or SMARTS pattern
         */
        public String getLabel() {
            return this.hasName() ? this.name : this.smarts;
        }
        //
        /**
         * Returns the line that represents this entry in the text format of the property.
         *
         * @return the SMARTS pattern, followed by a space and the name if the pattern is named
         */
        public String toLine() {
            return this.hasName() ? this.smarts + " " + this.name : this.smarts;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Record ParsedLine">
    /**
     * One parsed line of the text format that is neither blank nor a comment. Exactly one of entry and error message
     * is not null.
     *
     * @param lineNumber the line number in the parsed text, starting at 1
     * @param line the content of the line without leading and trailing whitespace
     * @param entry the entry of the line if it is valid, null otherwise
     * @param errorMessage the reason why the line is invalid, null if it is valid
     */
    public record ParsedLine(int lineNumber, String line, SmartsPatternEntry entry, String errorMessage) {
        /**
         * Returns whether the line holds a valid entry.
         *
         * @return true if the line is valid
         */
        public boolean isValid() {
            return Objects.isNull(this.errorMessage);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Lines starting with this prefix (after leading whitespace) are comments and ignored.
     */
    public static final String COMMENT_LINE_PREFIX = "#";

    /**
     * Line separator used when the text is created from entries.
     */
    public static final String LINE_SEPARATOR = "\n";
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Byte order mark that some editors write at the beginning of UTF-8 text files.
     */
    private static final String BYTE_ORDER_MARK = "\uFEFF";
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * Constructor with all parameters.
     *
     * @param bean the bean of this property
     * @param name the name of this property
     * @param initialValue the initial pattern list text
     * @throws NullPointerException if the initial value is null
     * @throws IllegalArgumentException if the initial value contains an invalid line
     */
    public SimpleSmartsPatternListProperty(Object bean, String name, String initialValue)
            throws NullPointerException, IllegalArgumentException {
        super(bean, name, SimpleSmartsPatternListProperty.checkText(initialValue));
    }
    //
    /**
     * Constructor with an empty pattern list as initial value.
     *
     * @param bean the bean of this property
     * @param name the name of this property
     */
    public SimpleSmartsPatternListProperty(Object bean, String name) {
        super(bean, name, "");
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get/set">
    /**
     * Sets the pattern list text.
     *
     * @param newValue the new pattern list text
     * @throws NullPointerException if the given value is null
     * @throws IllegalArgumentException if the given text contains a line with an invalid SMARTS pattern
     */
    @Override
    public void set(String newValue) throws NullPointerException, IllegalArgumentException {
        super.set(SimpleSmartsPatternListProperty.checkText(newValue));
    }

    /**
     * Sets the pattern list text.
     *
     * @param v the new pattern list text
     * @throws NullPointerException if the given value is null
     * @throws IllegalArgumentException if the given text contains a line with an invalid SMARTS pattern
     */
    @Override
    public void setValue(String v) throws NullPointerException, IllegalArgumentException {
        this.set(v);
    }

    /**
     * Returns the entries of the current pattern list text in their order.
     *
     * @return list of entries, empty if the text holds no patterns
     */
    public List<SmartsPatternEntry> getEntries() {
        return SimpleSmartsPatternListProperty.parseEntries(this.get());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Parses all lines of the given text that are neither blank nor comments and validates their SMARTS patterns.
     * Leading and trailing whitespace of a line and a byte order mark at the beginning of the text are ignored. The
     * SMARTS pattern of a line ends at the first whitespace character, the rest of the line is the name.
     *
     * @param aText the pattern list text
     * @return the parsed lines in their order, valid and invalid ones
     * @throws NullPointerException if the given text is null
     */
    public static List<ParsedLine> parseLines(String aText) throws NullPointerException {
        Objects.requireNonNull(aText, "Given text is null.");
        String[] tmpLines = SimpleSmartsPatternListProperty.removeByteOrderMark(aText).split("\\R", -1);
        List<ParsedLine> tmpParsedLines = new ArrayList<>(tmpLines.length);
        for (int i = 0; i < tmpLines.length; i++) {
            String tmpLine = tmpLines[i].strip();
            if (tmpLine.isEmpty() || tmpLine.startsWith(SimpleSmartsPatternListProperty.COMMENT_LINE_PREFIX)) {
                continue;
            }
            String[] tmpParts = tmpLine.split("\\s+", 2);
            String tmpSmarts = tmpParts[0];
            String tmpName = tmpParts.length > 1 ? tmpParts[1].strip() : "";
            String tmpError = SimpleSmartsPatternListProperty.validateSmarts(tmpSmarts);
            if (Objects.isNull(tmpError)) {
                tmpParsedLines.add(new ParsedLine(i + 1, tmpLine, new SmartsPatternEntry(tmpSmarts, tmpName), null));
            } else {
                tmpParsedLines.add(new ParsedLine(i + 1, tmpLine, null, tmpError));
            }
        }
        return tmpParsedLines;
    }

    /**
     * Returns the valid entries of the given text in their order; invalid lines are skipped.
     *
     * @param aText the pattern list text
     * @return list of valid entries
     * @throws NullPointerException if the given text is null
     */
    public static List<SmartsPatternEntry> parseEntries(String aText) throws NullPointerException {
        List<ParsedLine> tmpParsedLines = SimpleSmartsPatternListProperty.parseLines(aText);
        List<SmartsPatternEntry> tmpEntries = new ArrayList<>(tmpParsedLines.size());
        for (ParsedLine tmpParsedLine : tmpParsedLines) {
            if (tmpParsedLine.isValid()) {
                tmpEntries.add(tmpParsedLine.entry());
            }
        }
        return tmpEntries;
    }

    /**
     * Validates the given text and returns its invalid lines with the reasons.
     *
     * @param aText the pattern list text
     * @return invalid lines, empty if the text is valid
     * @throws NullPointerException if the given text is null
     */
    public static List<ParsedLine> validate(String aText) throws NullPointerException {
        List<ParsedLine> tmpParsedLines = SimpleSmartsPatternListProperty.parseLines(aText);
        List<ParsedLine> tmpInvalidLines = new ArrayList<>(0);
        for (ParsedLine tmpParsedLine : tmpParsedLines) {
            if (!tmpParsedLine.isValid()) {
                tmpInvalidLines.add(tmpParsedLine);
            }
        }
        return tmpInvalidLines;
    }

    /**
     * Counts the entries of the given text, i.e. its lines that are neither blank nor comments, without validating
     * them. For a valid text, this is the number of SMARTS patterns.
     *
     * @param aText the pattern list text
     * @return number of entries
     * @throws NullPointerException if the given text is null
     */
    public static int countEntries(String aText) throws NullPointerException {
        Objects.requireNonNull(aText, "Given text is null.");
        int tmpCount = 0;
        for (String tmpLine : SimpleSmartsPatternListProperty.removeByteOrderMark(aText).split("\\R")) {
            String tmpStrippedLine = tmpLine.strip();
            if (!tmpStrippedLine.isEmpty() && !tmpStrippedLine.startsWith(SimpleSmartsPatternListProperty.COMMENT_LINE_PREFIX)) {
                tmpCount++;
            }
        }
        return tmpCount;
    }

    /**
     * Creates the pattern list text of the given entries, one line per entry.
     *
     * @param anEntryList the entries
     * @return pattern list text, an empty string if no entries are given
     * @throws NullPointerException if the given list or an entry is null
     */
    public static String toText(List<SmartsPatternEntry> anEntryList) throws NullPointerException {
        Objects.requireNonNull(anEntryList, "Given entry list is null.");
        List<String> tmpLines = new ArrayList<>(anEntryList.size());
        for (SmartsPatternEntry tmpEntry : anEntryList) {
            tmpLines.add(Objects.requireNonNull(tmpEntry, "Given entry is null.").toLine());
        }
        return String.join(SimpleSmartsPatternListProperty.LINE_SEPARATOR, tmpLines);
    }

    /**
     * Validates a single SMARTS pattern with the CDK SMARTS parser.
     *
     * @param aSmarts the SMARTS pattern
     * @return null if the SMARTS pattern is valid, otherwise a language-specific description of the problem,
     * including the message of the CDK SMARTS parser if it rejects the pattern
     * @throws NullPointerException if the given SMARTS pattern is null
     */
    public static String validateSmarts(String aSmarts) throws NullPointerException {
        Objects.requireNonNull(aSmarts, "Given SMARTS pattern is null.");
        if (aSmarts.isEmpty()) {
            return Message.get("SimpleSmartsPatternListProperty.Error.missingSmarts");
        }
        for (int i = 0; i < aSmarts.length(); i++) {
            if (Character.isWhitespace(aSmarts.charAt(i))) {
                return Message.get("SimpleSmartsPatternListProperty.Error.whitespaceInSmarts");
            }
        }
        if (aSmarts.indexOf('>') >= 0) {
            //the parser accepts reaction SMARTS, whose role-specific atoms do not match molecules as expected
            return Message.get("SimpleSmartsPatternListProperty.Error.reactionSmarts");
        }
        SmartsResult tmpResult;
        try {
            //thread-safe variant of Smarts.parse(), whose error message is stored in a static variable
            tmpResult = Smarts.parseToResult(new QueryAtomContainer(SilentChemObjectBuilder.getInstance()), aSmarts);
        } catch (RuntimeException anException) {
            //the parser throws instead of reporting an error for some incomplete patterns, e.g. "["
            return String.format(Message.get("SimpleSmartsPatternListProperty.Error.parserException"), anException.toString());
        }
        if (!tmpResult.ok()) {
            return String.format(Message.get("SimpleSmartsPatternListProperty.Error.invalidSmarts"),
                    tmpResult.getPosition(), tmpResult.getMessage());
        }
        try {
            //the fragmenter compiles the patterns like this, so make sure that it works
            SmartsPattern.create(aSmarts, SilentChemObjectBuilder.getInstance());
        } catch (RuntimeException anException) {
            return String.format(Message.get("SimpleSmartsPatternListProperty.Error.notCompilable"), anException.getMessage());
        }
        return null;
    }

    /**
     * Validates the name of a pattern, which must not contain line breaks.
     *
     * @param aName the name, may be empty
     * @return null if the name is valid, otherwise a language-specific description of the problem
     * @throws NullPointerException if the given name is null
     */
    public static String validateName(String aName) throws NullPointerException {
        Objects.requireNonNull(aName, "Given name is null.");
        if (aName.indexOf('\n') >= 0 || aName.indexOf('\r') >= 0) {
            return Message.get("SimpleSmartsPatternListProperty.Error.lineBreakInName");
        }
        return null;
    }

    /**
     * Creates a language-specific description of the given invalid lines, one line per invalid line.
     *
     * @param anInvalidLineList invalid lines, e.g. as returned by {@link #validate(String)}
     * @return description of the invalid lines
     * @throws NullPointerException if the given list is null
     */
    public static String describeInvalidLines(List<ParsedLine> anInvalidLineList) throws NullPointerException {
        Objects.requireNonNull(anInvalidLineList, "Given list is null.");
        List<String> tmpDescriptions = new ArrayList<>(anInvalidLineList.size());
        for (ParsedLine tmpInvalidLine : anInvalidLineList) {
            tmpDescriptions.add(String.format(Message.get("SimpleSmartsPatternListProperty.invalidLine"),
                    tmpInvalidLine.lineNumber(), tmpInvalidLine.line(), tmpInvalidLine.errorMessage()));
        }
        return String.join(SimpleSmartsPatternListProperty.LINE_SEPARATOR, tmpDescriptions);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static methods">
    /**
     * Removes a byte order mark at the beginning of the given text, which some editors write into UTF-8 text files.
     *
     * @param aText the text
     * @return the text without a leading byte order mark
     */
    private static String removeByteOrderMark(String aText) {
        return aText.startsWith(SimpleSmartsPatternListProperty.BYTE_ORDER_MARK)
                ? aText.substring(SimpleSmartsPatternListProperty.BYTE_ORDER_MARK.length()) : aText;
    }

    /**
     * Checks that the given text is not null and contains no invalid lines.
     *
     * @param aText the pattern list text to check
     * @return the given text
     * @throws NullPointerException if the given text is null
     * @throws IllegalArgumentException if the given text contains invalid lines
     */
    private static String checkText(String aText) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aText, "Given pattern list text is null.");
        List<ParsedLine> tmpInvalidLines = SimpleSmartsPatternListProperty.validate(aText);
        if (!tmpInvalidLines.isEmpty()) {
            throw new IllegalArgumentException("The given SMARTS pattern list contains invalid lines:"
                    + SimpleSmartsPatternListProperty.LINE_SEPARATOR
                    + SimpleSmartsPatternListProperty.describeInvalidLines(tmpInvalidLines));
        }
        return aText;
    }
    //</editor-fold>
}
