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

package de.unijena.cheminf.mortar.preference;

import de.unijena.cheminf.mortar.configuration.Configuration;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.SmartsFragmenter;
import de.unijena.cheminf.mortar.model.util.SimpleSmartsPatternListProperty;

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Base64;
import java.util.List;
import java.util.Locale;

/**
 * Test class for the {@link PreferenceUtil} translation between JavaFX properties and preferences, here for the
 * SMARTS pattern list property, whose text contains characters that preferences do not allow.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class PreferenceUtilTest {
    /**
     * Pattern list text with SMARTS special characters, non-ASCII characters in names, comments, blank lines, and
     * trailing whitespace.
     */
    private static final String SPECIAL_CHARACTERS_PATTERN_LIST = String.join("\n",
            "# patterns with = # $ ! & @ + * ~ % and ; , characters",
            "[#6;R1]~[$([OX2H]),$([O-])] ring carbon with OH or O-",
            "*!@[N+](=O)[O-] nitro group, non-ring bond",
            "",
            "c%10ccccc%10 benzene (two-digit ring closure)",
            "[C;!$(C=O)]-[OX2H1]-&!@* alcohol ether",
            "C=,#C\tdouble or triple \u00e4\u00f6\u00fc \u00df \u20ac",
            "[CH2]~[NH2+] trailing whitespace   ",
            "C/C=C/C");
    //
    /**
     * Constructor to initialise locale and configuration.
     *
     * @throws Exception if anything goes wrong
     */
    public PreferenceUtilTest() throws Exception {
        Locale.setDefault(Locale.of("en", "GB"));
        Configuration.getInstance();
    }
    //
    /**
     * Tests that a pattern list with special characters and an empty pattern list survive the translation into a
     * preference container, writing it to a file, reading it again, and updating fresh fragmenter settings from it.
     *
     * @param aTempDir temporary directory for the preference container file
     * @throws Exception if anything goes wrong
     */
    @Test
    public void smartsPatternListRoundTripTest(@TempDir Path aTempDir) throws Exception {
        for (String tmpText : new String[] {PreferenceUtilTest.SPECIAL_CHARACTERS_PATTERN_LIST, ""}) {
            SmartsFragmenter tmpFragmenter = new SmartsFragmenter();
            tmpFragmenter.setPatternsSetting(tmpText);
            tmpFragmenter.setFilterOverlappingMatchesSetting(false);
            tmpFragmenter.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
            Assertions.assertDoesNotThrow(() -> PreferenceUtil.checkPropertiesForPreferenceRestrictions(tmpFragmenter.settingsProperties()));
            File tmpFile = aTempDir.resolve("SmartsFragmenter" + (tmpText.isEmpty() ? "Empty" : "") + ".txt").toFile();
            PreferenceContainer tmpContainer = PreferenceUtil.translateJavaFxPropertiesToPreferences(
                    tmpFragmenter.settingsProperties(), tmpFile.getPath());
            //all settings are translated, none is skipped because of invalid content
            Assertions.assertEquals(tmpFragmenter.settingsProperties().size(), tmpContainer.getSize());
            String tmpContent = tmpContainer.getPreferences(tmpFragmenter.patternsSettingProperty().getName())[0].getContentRepresentative();
            Assertions.assertTrue(SingleTermPreference.isValidContent(tmpContent), tmpContent);
            tmpContainer.writeRepresentation();
            PreferenceContainer tmpReloadedContainer = new PreferenceContainer(tmpFile);
            SmartsFragmenter tmpReloadedFragmenter = new SmartsFragmenter();
            //a value different from the persisted one, which must be overwritten
            tmpReloadedFragmenter.setPatternsSetting("[OH]");
            PreferenceUtil.updatePropertiesFromPreferences(tmpReloadedFragmenter.settingsProperties(), tmpReloadedContainer);
            Assertions.assertEquals(tmpText, tmpReloadedFragmenter.getPatternsSetting());
            Assertions.assertFalse(tmpReloadedFragmenter.getFilterOverlappingMatchesSetting());
            Assertions.assertEquals(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION,
                    tmpReloadedFragmenter.getFragmentSaturationSetting());
        }
    }
    //
    /**
     * Tests that a persisted pattern list with an invalid SMARTS pattern, e.g. because a newer CDK version rejects it,
     * is restored without the invalid line, and that a persisted value that cannot be decoded leaves the setting
     * unchanged.
     *
     * @param aTempDir temporary directory for the preference container file
     */
    @Test
    public void invalidPersistedSmartsPatternListTest(@TempDir Path aTempDir) {
        SimpleSmartsPatternListProperty tmpProperty = new SimpleSmartsPatternListProperty(this, "SMARTS patterns setting");
        List<Property<?>> tmpProperties = List.of(tmpProperty);
        PreferenceContainer tmpContainer = new PreferenceContainer(aTempDir.resolve("Invalid.txt").toFile().getPath());
        String tmpPersistedText = "[OH] hydroxyl\nC(C broken\nC=O carbonyl";
        tmpContainer.add(new SingleTermPreference(tmpProperty.getName(), "base64:"
                + Base64.getUrlEncoder().withoutPadding().encodeToString(tmpPersistedText.getBytes(StandardCharsets.UTF_8))));
        PreferenceUtil.updatePropertiesFromPreferences(tmpProperties, tmpContainer);
        Assertions.assertEquals("[OH] hydroxyl\nC=O carbonyl", tmpProperty.get());
        //content without the expected prefix, the setting remains unchanged
        PreferenceContainer tmpBrokenContainer = new PreferenceContainer(aTempDir.resolve("Broken.txt").toFile().getPath());
        tmpBrokenContainer.add(new SingleTermPreference(tmpProperty.getName(), "plain text without prefix"));
        PreferenceUtil.updatePropertiesFromPreferences(tmpProperties, tmpBrokenContainer);
        Assertions.assertEquals("[OH] hydroxyl\nC=O carbonyl", tmpProperty.get());
        //invalid Base64 code, the setting remains unchanged
        PreferenceContainer tmpInvalidBase64Container = new PreferenceContainer(aTempDir.resolve("Base64.txt").toFile().getPath());
        tmpInvalidBase64Container.add(new SingleTermPreference(tmpProperty.getName(), "base64:a"));
        PreferenceUtil.updatePropertiesFromPreferences(tmpProperties, tmpInvalidBase64Container);
        Assertions.assertEquals("[OH] hydroxyl\nC=O carbonyl", tmpProperty.get());
    }
}
