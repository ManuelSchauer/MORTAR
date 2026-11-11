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

package de.unijena.cheminf.mortar.model.io;

import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.filter.MoleculeFilterReason;
import de.unijena.cheminf.mortar.model.filter.MoleculeFilterSettings;
import de.unijena.cheminf.mortar.model.filter.MoleculeFilterSettingsSnapshot;
import de.unijena.cheminf.mortar.model.filter.MoleculeFilterSummary;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.File;
import java.net.URL;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * End-to-end tests of the molecule filters at import through the public
 * {@link Importer#importMoleculeFile(File, boolean, boolean, boolean, MoleculeFilterSettingsSnapshot)} method, with a
 * small synthetic dirty SD file and SMILES file (resources "MoleculeFiltersTestFile.sdf" and
 * "MoleculeFiltersTestFile.smi" in this package). Both files contain the same twelve molecules in the same order:
 * <ol start="0">
 *     <li>ethanol (CCO)</li>
 *     <li>ethanol_duplicate (OCC, duplicate of 0)</li>
 *     <li>L-alanine</li>
 *     <li>D-alanine (stereoisomer of 2)</li>
 *     <li>triethylamine_hydrochloride (salt, CCN(CC)CC.Cl)</li>
 *     <li>triethylamine (parent compound of 4)</li>
 *     <li>sodium_acetate (salt with a disallowed counter-ion)</li>
 *     <li>paracetamol_hydrate (hydrate)</li>
 *     <li>cisplatin (connected metal complex)</li>
 *     <li>chloride_ion (lone ion, one heavy atom)</li>
 *     <li>R_group_molecule (*CCO, pseudo atom)</li>
 *     <li>triacontane (30 heavy atoms)</li>
 * </ol>
 * The SD file additionally contains a 13th record without atoms ("empty_record"), for which no SMILES string can be
 * generated. Every SD record has the data fields "Comment" and "Source". The tests run headless and do not touch the
 * MORTAR settings directory.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ImporterMoleculeFilterTest {
    //<editor-fold desc="Static initializer">
    /**
     * Sets the default locale to British English, required for the message texts of the filter enums.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Name of the SD test file resource.
     */
    private static final String SD_FILE_NAME = "MoleculeFiltersTestFile.sdf";

    /**
     * Name of the SMILES test file resource.
     */
    private static final String SMILES_FILE_NAME = "MoleculeFiltersTestFile.smi";

    /**
     * Names of the twelve molecules contained in both test files, in file order.
     */
    private static final List<String> MOLECULE_NAMES = List.of("ethanol", "ethanol_duplicate", "L-alanine",
            "D-alanine", "triethylamine_hydrochloride", "triethylamine", "sodium_acetate", "paracetamol_hydrate",
            "cisplatin", "chloride_ion", "R_group_molecule", "triacontane");

    /**
     * Name of the additional record without atoms in the SD file.
     */
    private static final String EMPTY_RECORD_NAME = "empty_record";
    //</editor-fold>
    //
    //<editor-fold desc="Tests">
    /**
     * Tests that an import with the default (inactive) filter settings gives exactly the same molecules (names and
     * SMILES strings) as the old method without filter settings, for both files, all combinations of the
     * stereochemistry, Kekule, and keep atom container settings, and also with a snapshot of a new settings container.
     * Also checks the neutral summary and that no molecule is logged as removed (the empty SD record is not imported,
     * as before, but logged as without filters).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDefaultSnapshotEqualsImportWithoutFilters() throws Exception {
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            for (boolean tmpIsKeepAtomContainer : new boolean[] {false, true}) {
                for (boolean tmpIsRegardStereo : new boolean[] {true, false}) {
                    for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
                        Importer tmpOldImporter = this.createImporter(tmpIsKeepAtomContainer);
                        List<MoleculeDataModel> tmpOldResult = tmpOldImporter.importMoleculeFile(this.getResourceFile(tmpFileName),
                                tmpIsRegardStereo, true, tmpIsKekulizationEnforced);
                        Importer tmpNewImporter = this.createImporter(tmpIsKeepAtomContainer);
                        List<String> tmpLogLines = new ArrayList<>();
                        List<MoleculeDataModel> tmpNewResult = this.importAndCaptureLog(tmpNewImporter, tmpFileName,
                                tmpIsRegardStereo, tmpIsKekulizationEnforced, MoleculeFilterSettingsSnapshot.createDefault(),
                                tmpLogLines);
                        List<MoleculeDataModel> tmpSettingsResult = this.createImporter(tmpIsKeepAtomContainer).importMoleculeFile(
                                this.getResourceFile(tmpFileName), tmpIsRegardStereo, true, tmpIsKekulizationEnforced,
                                new MoleculeFilterSettings().createSnapshot());
                        Assertions.assertEquals(ImporterMoleculeFilterTest.MOLECULE_NAMES, this.getNames(tmpOldResult));
                        Assertions.assertEquals(this.getNames(tmpOldResult), this.getNames(tmpNewResult));
                        Assertions.assertEquals(this.getSmiles(tmpOldResult), this.getSmiles(tmpNewResult));
                        Assertions.assertEquals(this.getNames(tmpOldResult), this.getNames(tmpSettingsResult));
                        Assertions.assertEquals(this.getSmiles(tmpOldResult), this.getSmiles(tmpSettingsResult));
                        for (MoleculeDataModel tmpMolecule : tmpNewResult) {
                            Assertions.assertEquals(tmpIsKeepAtomContainer, tmpMolecule.isKeepAtomContainer());
                        }
                        boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
                        int tmpReadCount = tmpIsSdFile ? 13 : 12;
                        Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = tmpIsSdFile
                                ? Map.of(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1) : Map.of();
                        this.assertSummary(tmpNewImporter.getMoleculeFilterSummary(), tmpReadCount, 12, 0,
                                tmpExpectedRemovedCounts);
                        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(),
                                tmpOldImporter.getMoleculeFilterSummary().getSettingsSnapshot());
                        Assertions.assertFalse(tmpNewImporter.getMoleculeFilterSummary().isDataSetChangedByFilters());
                        Assertions.assertEquals(0, tmpNewImporter.getMoleculeFilterSummary().getRemovedByFiltersCount());
                        Assertions.assertEquals(List.of(), this.getRemovalLines(tmpLogLines));
                        Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines,
                                "Successfully imported 12 molecules from file: " + tmpFileName));
                        Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines, "Molecule filters ("));
                    }
                }
            }
        }
    }

    /**
     * Tests the duplicate filter alone: with stereochemistry, only the second ethanol is removed (first occurrence
     * kept); without stereochemistry, D-alanine is a duplicate of L-alanine, too. The salt and its parent compound are
     * no duplicates because the salt is not reduced. With an active filter, the empty SD record is logged as removed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDuplicateFilter() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 0, 0, true);
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            //with stereochemistry
            Importer tmpImporter = this.createImporter(false);
            List<String> tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                    tmpSnapshot, tmpLogLines);
            Assertions.assertEquals(this.withoutNames("ethanol_duplicate"), this.getNames(tmpResult));
            List<String> tmpSmiles = this.getSmiles(tmpResult);
            Assertions.assertNotEquals(tmpSmiles.get(1), tmpSmiles.get(2), "L- and D-alanine must differ with stereo");
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.DUPLICATE, 1);
            List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                    this.createRemovalLine("ethanol_duplicate", MoleculeFilterReason.DUPLICATE)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME,
                        MoleculeFilterReason.SMILES_GENERATION_FAILED));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 11, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
            Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines, "Successfully imported 11 molecules"));
            //without stereochemistry
            tmpImporter = this.createImporter(false);
            tmpResult = tmpImporter.importMoleculeFile(this.getResourceFile(tmpFileName), false, true, false, tmpSnapshot);
            Assertions.assertEquals(this.withoutNames("ethanol_duplicate", "D-alanine"), this.getNames(tmpResult));
            Assertions.assertEquals(tmpIsSdFile ? 13 : 12, tmpImporter.getMoleculeFilterSummary().getReadCount());
            Assertions.assertEquals(10, tmpImporter.getMoleculeFilterSummary().getKeptCount());
            Assertions.assertEquals(2, tmpImporter.getMoleculeFilterSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
            Assertions.assertEquals(0, tmpImporter.getMoleculeFilterSummary().getPendingCount());
        }
    }

    /**
     * Tests the disconnected structures setting "remove molecule" alone: the salts and the hydrate are removed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testRemoveDisconnectedStructures() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE, false, false, 0, 0, false);
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            Importer tmpImporter = this.createImporter(false);
            List<String> tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                    tmpSnapshot, tmpLogLines);
            Assertions.assertEquals(this.withoutNames("triethylamine_hydrochloride", "sodium_acetate", "paracetamol_hydrate"),
                    this.getNames(tmpResult));
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.DISCONNECTED_STRUCTURE, 3);
            List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                    this.createRemovalLine("triethylamine_hydrochloride", MoleculeFilterReason.DISCONNECTED_STRUCTURE),
                    this.createRemovalLine("sodium_acetate", MoleculeFilterReason.DISCONNECTED_STRUCTURE),
                    this.createRemovalLine("paracetamol_hydrate", MoleculeFilterReason.DISCONNECTED_STRUCTURE)));
            if (tmpIsSdFile) {
                //rejected as such because a structural filter is active
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME,
                        MoleculeFilterReason.EMPTY_MOLECULE));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 9, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
            for (String tmpSmiles : this.getSmiles(tmpResult)) {
                Assertions.assertFalse(tmpSmiles.contains("."), tmpSmiles);
            }
        }
    }

    /**
     * Tests the disconnected structures setting "keep largest component" alone and together with the duplicate filter,
     * with both values of the keep atom container setting: the salts and the hydrate are reduced but keep their names
     * and SD data fields; the reduced salt gets the SMILES string of its parent compound, which then becomes a
     * duplicate of it (first occurrence kept, i.e. the entry named like the salt).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testKeepLargestComponent() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, false, false, 0, 0, false);
        MoleculeFilterSettingsSnapshot tmpSnapshotWithDuplicates = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, false, false, 0, 0, true);
        for (boolean tmpIsKeepAtomContainer : new boolean[] {false, true}) {
            //SD file with data fields
            Importer tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
            List<MoleculeDataModel> tmpResult = tmpImporter.importMoleculeFile(
                    this.getResourceFile(ImporterMoleculeFilterTest.SD_FILE_NAME), true, true, false, tmpSnapshot);
            Assertions.assertEquals(ImporterMoleculeFilterTest.MOLECULE_NAMES, this.getNames(tmpResult));
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), 13, 12, 3,
                    Map.of(MoleculeFilterReason.EMPTY_MOLECULE, 1));
            MoleculeDataModel tmpReducedSalt = tmpResult.get(4);
            MoleculeDataModel tmpParent = tmpResult.get(5);
            Assertions.assertEquals("triethylamine_hydrochloride", tmpReducedSalt.getName());
            Assertions.assertEquals(tmpParent.getUniqueSmiles(), tmpReducedSalt.getUniqueSmiles());
            Assertions.assertEquals("salt of triethylamine", tmpReducedSalt.getProperties().get("Comment"));
            Assertions.assertEquals("synthetic MORTAR import filter test set", tmpReducedSalt.getProperties().get("Source"));
            Assertions.assertEquals("triethylamine_hydrochloride",
                    tmpReducedSalt.getProperties().get(Importer.MOLECULE_NAME_PROPERTY_KEY));
            Assertions.assertEquals(tmpIsKeepAtomContainer, tmpReducedSalt.isKeepAtomContainer());
            IAtomContainer tmpReducedSaltAtomContainer = tmpReducedSalt.getAtomContainer();
            //seven heavy atoms of triethylamine, hydrogen atoms are implicit
            Assertions.assertEquals(7, tmpReducedSaltAtomContainer.getAtomCount());
            if (tmpIsKeepAtomContainer) {
                Assertions.assertEquals("salt of triethylamine", tmpReducedSaltAtomContainer.getProperty("Comment"));
                Assertions.assertEquals("triethylamine_hydrochloride", tmpReducedSaltAtomContainer.getTitle());
            }
            Assertions.assertEquals("salt with a disallowed counter-ion", tmpResult.get(6).getProperties().get("Comment"));
            Assertions.assertEquals("hydrate", tmpResult.get(7).getProperties().get("Comment"));
            for (int i : new int[] {4, 6, 7}) {
                Assertions.assertFalse(tmpResult.get(i).getUniqueSmiles().contains("."), tmpResult.get(i).getUniqueSmiles());
            }
            Assertions.assertEquals("OCC", tmpResult.getFirst().getUniqueSmiles());
            //SD file, additionally with the duplicate filter
            tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
            List<String> tmpLogLines = new ArrayList<>();
            tmpResult = this.importAndCaptureLog(tmpImporter, ImporterMoleculeFilterTest.SD_FILE_NAME, true, false,
                    tmpSnapshotWithDuplicates, tmpLogLines);
            Assertions.assertEquals(this.withoutNames("ethanol_duplicate", "triethylamine"), this.getNames(tmpResult));
            Assertions.assertEquals("salt of triethylamine", tmpResult.get(3).getProperties().get("Comment"));
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), 13, 10, 3, Map.of(
                    MoleculeFilterReason.DUPLICATE, 2,
                    MoleculeFilterReason.EMPTY_MOLECULE, 1));
            Assertions.assertEquals(List.of(
                    this.createRemovalLine("ethanol_duplicate", MoleculeFilterReason.DUPLICATE),
                    this.createRemovalLine("triethylamine", MoleculeFilterReason.DUPLICATE),
                    this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, MoleculeFilterReason.EMPTY_MOLECULE)),
                    this.getRemovalLines(tmpLogLines));
            //SMILES file
            tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
            tmpResult = tmpImporter.importMoleculeFile(this.getResourceFile(ImporterMoleculeFilterTest.SMILES_FILE_NAME),
                    true, true, false, tmpSnapshotWithDuplicates);
            Assertions.assertEquals(this.withoutNames("ethanol_duplicate", "triethylamine"), this.getNames(tmpResult));
            Assertions.assertEquals("triethylamine_hydrochloride", tmpResult.get(3).getName());
            //unique SMILES strings of triethylamine and acetate
            Assertions.assertEquals("N(CC)(CC)CC", tmpResult.get(3).getUniqueSmiles());
            Assertions.assertEquals("O=C([O-])C", tmpResult.get(4).getUniqueSmiles());
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), 12, 10, 3, Map.of(MoleculeFilterReason.DUPLICATE, 2));
        }
    }

    /**
     * Tests the element filter alone and after the reduction to the largest component: the sodium salt is removed
     * because of its counter-ion only if it is not reduced; the metal complex and the molecule with an R group are
     * always removed; chlorine (hydrochloride, chloride ion) is allowed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testElementFilter() throws Exception {
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            Importer tmpImporter = this.createImporter(false);
            List<String> tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                    this.createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, true, false, 0, 0, false),
                    tmpLogLines);
            Assertions.assertEquals(this.withoutNames("sodium_acetate", "cisplatin", "R_group_molecule"), this.getNames(tmpResult));
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.DISALLOWED_ELEMENT, 3);
            List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                    this.createRemovalLine("sodium_acetate", MoleculeFilterReason.DISALLOWED_ELEMENT),
                    this.createRemovalLine("cisplatin", MoleculeFilterReason.DISALLOWED_ELEMENT),
                    this.createRemovalLine("R_group_molecule", MoleculeFilterReason.DISALLOWED_ELEMENT)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME,
                        MoleculeFilterReason.EMPTY_MOLECULE));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 9, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
            tmpImporter = this.createImporter(false);
            tmpResult = tmpImporter.importMoleculeFile(this.getResourceFile(tmpFileName), true, true, false,
                    this.createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, false, 0, 0, false));
            Assertions.assertEquals(this.withoutNames("cisplatin", "R_group_molecule"), this.getNames(tmpResult));
            Assertions.assertEquals(2, tmpImporter.getMoleculeFilterSummary().getRemovedCount(MoleculeFilterReason.DISALLOWED_ELEMENT));
            Assertions.assertEquals(3, tmpImporter.getMoleculeFilterSummary().getReducedToLargestComponentCount());
        }
    }

    /**
     * Tests the heavy atom count filter alone with inclusive bounds: [3, 29] removes the chloride ion (1) and
     * triacontane (30) and keeps ethanol (3) and the R group molecule (3, pseudo atoms are not counted); [2, 30] keeps
     * triacontane.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testHeavyAtomCountFilter() throws Exception {
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            Importer tmpImporter = this.createImporter(false);
            List<String> tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                    this.createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 3, 29, false),
                    tmpLogLines);
            Assertions.assertEquals(this.withoutNames("chloride_ion", "triacontane"), this.getNames(tmpResult));
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE, 2);
            List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                    this.createRemovalLine("chloride_ion", MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE),
                    this.createRemovalLine("triacontane", MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME,
                        MoleculeFilterReason.EMPTY_MOLECULE));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 10, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
            tmpResult = this.createImporter(false).importMoleculeFile(this.getResourceFile(tmpFileName), true, true, false,
                    this.createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 2, 30, false));
            Assertions.assertEquals(this.withoutNames("chloride_ion"), this.getNames(tmpResult));
        }
    }

    /**
     * Tests all filters together (largest component, elements, heavy atom count [2, 20], duplicates) on both files:
     * kept molecules, their SMILES strings in the SMILES file, the summary, the display text for the summary dialog,
     * and exactly one log line per removed molecule with its name and reason, in file order.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testAllFiltersTogether() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, true, 2, 20, true);
        List<String> tmpExpectedNames = List.of("ethanol", "L-alanine", "D-alanine", "triethylamine_hydrochloride",
                "sodium_acetate", "paracetamol_hydrate");
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            for (boolean tmpIsKeepAtomContainer : new boolean[] {false, true}) {
                boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
                Importer tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
                List<String> tmpLogLines = new ArrayList<>();
                List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                        tmpSnapshot, tmpLogLines);
                Assertions.assertEquals(tmpExpectedNames, this.getNames(tmpResult));
                if (!tmpIsSdFile) {
                    //unique SMILES strings of ethanol, L-alanine, D-alanine, triethylamine, acetate, and paracetamol
                    Assertions.assertEquals(List.of("OCC", "C[C@@H](C(=O)O)N", "C[C@H](C(=O)O)N", "N(CC)(CC)CC",
                            "O=C([O-])C", "O=C(Nc1ccc(O)cc1)C"), this.getSmiles(tmpResult));
                }
                Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.DUPLICATE, 2);
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.DISALLOWED_ELEMENT, 2);
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE, 2);
                List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                        this.createRemovalLine("ethanol_duplicate", MoleculeFilterReason.DUPLICATE),
                        this.createRemovalLine("triethylamine", MoleculeFilterReason.DUPLICATE),
                        this.createRemovalLine("cisplatin", MoleculeFilterReason.DISALLOWED_ELEMENT),
                        this.createRemovalLine("chloride_ion", MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE),
                        this.createRemovalLine("R_group_molecule", MoleculeFilterReason.DISALLOWED_ELEMENT),
                        this.createRemovalLine("triacontane", MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE)));
                if (tmpIsSdFile) {
                    tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                    tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME,
                            MoleculeFilterReason.EMPTY_MOLECULE));
                }
                MoleculeFilterSummary tmpSummary = tmpImporter.getMoleculeFilterSummary();
                this.assertSummary(tmpSummary, tmpIsSdFile ? 13 : 12, 6, 3, tmpExpectedRemovedCounts);
                Assertions.assertEquals(tmpSnapshot, tmpSummary.getSettingsSnapshot());
                Assertions.assertTrue(tmpSummary.isDataSetChangedByFilters());
                Assertions.assertEquals(tmpIsSdFile ? 7 : 6, tmpSummary.getRemovedByFiltersCount());
                //exactly one line per removed molecule, with its name and reason
                Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
                Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines,
                        "Successfully imported 6 molecules from file: " + tmpFileName));
                Assertions.assertEquals(List.of(tmpSummary.toString()), this.getLinesStartingWith(tmpLogLines, "Molecule filters ("));
                //text of the summary dialog, only reasons with removed molecules
                String tmpExpectedDisplayText = "Molecules read: " + (tmpIsSdFile ? 13 : 12) + "\n" +
                        "Molecules imported: 6\n" +
                        "Molecules reduced to their largest component: 3\n" +
                        "Molecules removed: " + (tmpIsSdFile ? 7 : 6) + "\n" +
                        (tmpIsSdFile ? "  - empty molecule: 1\n" : "") +
                        "  - disallowed element: 2\n" +
                        "  - heavy atom count out of range: 2\n" +
                        "  - duplicate: 2";
                Assertions.assertEquals(tmpExpectedDisplayText, tmpSummary.getDisplayText());
            }
        }
    }

    /**
     * Tests an import in which the filters remove every molecule: an empty list (not null) is returned, and the
     * summary reports all molecules as removed by the filters, so that the GUI can tell the user. Also checks that the
     * message texts of the GUI dialog exist.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testAllMoleculesRemoved() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 100, 200, false);
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            Importer tmpImporter = this.createImporter(false);
            List<String> tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, true, false,
                    tmpSnapshot, tmpLogLines);
            Assertions.assertNotNull(tmpResult);
            Assertions.assertTrue(tmpResult.isEmpty());
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE, 12);
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
            }
            MoleculeFilterSummary tmpSummary = tmpImporter.getMoleculeFilterSummary();
            this.assertSummary(tmpSummary, tmpIsSdFile ? 13 : 12, 0, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpSummary.getReadCount(), tmpSummary.getRemovedByFiltersCount());
            Assertions.assertTrue(tmpSummary.isDataSetChangedByFilters());
            Assertions.assertEquals(tmpIsSdFile ? 13 : 12, this.getRemovalLines(tmpLogLines).size());
            Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines, "Successfully imported 0 molecules"));
        }
        //the texts of the summary dialog of the main view controller exist (the dialog itself needs the GUI)
        for (String tmpKey : new String[] {"MainViewController.ImportFilterSummary.Title",
                "MainViewController.ImportFilterSummary.Header", "MainViewController.ImportFilterSummary.AllRemoved.Header",
                "MainViewController.ImportFilterSummary.AllRemoved.Hint", "MainViewController.ImportFilterSummary.LogFileHint"}) {
            String tmpText = Message.get(tmpKey);
            Assertions.assertFalse(tmpText.isBlank(), tmpKey);
            Assertions.assertFalse(tmpText.contains("not found"), tmpText);
        }
    }

    /**
     * Tests the summary handling of the importer: it is null before the first import and after an import of an
     * unsupported file type, and the old import method stores a summary of the default settings.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testSummaryAvailability() throws Exception {
        Importer tmpImporter = this.createImporter(false);
        Assertions.assertNull(tmpImporter.getMoleculeFilterSummary());
        tmpImporter.importMoleculeFile(this.getResourceFile(ImporterMoleculeFilterTest.SMILES_FILE_NAME), true, true, false,
                this.createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE, false, false, 0, 0, false));
        Assertions.assertEquals(9, tmpImporter.getMoleculeFilterSummary().getKeptCount());
        //unsupported file type, null is returned and the summary of the previous import is discarded
        Assertions.assertNull(tmpImporter.importMoleculeFile(new File("unsupported.xyz"), true, true, false,
                MoleculeFilterSettingsSnapshot.createDefault()));
        Assertions.assertNull(tmpImporter.getMoleculeFilterSummary());
        Assertions.assertThrows(NullPointerException.class, () -> tmpImporter.importMoleculeFile(
                this.getResourceFile(ImporterMoleculeFilterTest.SMILES_FILE_NAME), true, true, false, null));
        tmpImporter.importMoleculeFile(this.getResourceFile(ImporterMoleculeFilterTest.SMILES_FILE_NAME), true, true, false);
        Assertions.assertEquals(MoleculeFilterSettingsSnapshot.createDefault(),
                tmpImporter.getMoleculeFilterSummary().getSettingsSnapshot());
        Assertions.assertEquals(12, tmpImporter.getMoleculeFilterSummary().getKeptCount());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Creates an importer with default global settings except for the given keep atom container setting. Since the
     * getter of this deprecated setting always returns false in the settings container, it is overridden here to reach
     * the respective branch of the importer.
     *
     * @param isKeepAtomContainer value of the keep atom container in data model setting
     * @return new importer
     */
    @SuppressWarnings("deprecation")
    private Importer createImporter(boolean isKeepAtomContainer) {
        SettingsContainer tmpSettingsContainer = new SettingsContainer() {
            @Override
            public boolean getKeepAtomContainerInDataModelSetting() {
                return isKeepAtomContainer;
            }
        };
        Assertions.assertEquals(isKeepAtomContainer, tmpSettingsContainer.getKeepAtomContainerInDataModelSetting());
        return new Importer(tmpSettingsContainer);
    }

    /**
     * Returns the test resource file with the given name from the package of this class.
     *
     * @param aFileName name of the resource file
     * @return the file
     * @throws Exception if the resource cannot be found
     */
    private File getResourceFile(String aFileName) throws Exception {
        URL tmpURL = this.getClass().getResource(aFileName);
        Assertions.assertNotNull(tmpURL, aFileName);
        return Paths.get(tmpURL.toURI()).toFile();
    }

    /**
     * Imports the given test file with filling of open valences and the given settings and collects the formatted
     * messages logged by the importer meanwhile.
     *
     * @param anImporter the importer to use
     * @param aFileName name of the resource file
     * @param isRegardStereo stereochemistry setting
     * @param isKekulizationEnforced Kekule setting
     * @param aSnapshot molecule filter settings
     * @param aLogLinesList list to add the logged messages to
     * @return imported molecules
     * @throws Exception if anything goes wrong
     */
    private List<MoleculeDataModel> importAndCaptureLog(Importer anImporter, String aFileName, boolean isRegardStereo,
                                                        boolean isKekulizationEnforced, MoleculeFilterSettingsSnapshot aSnapshot,
                                                        List<String> aLogLinesList) throws Exception {
        Logger tmpLogger = Logger.getLogger(Importer.class.getName());
        SimpleFormatter tmpFormatter = new SimpleFormatter();
        Handler tmpHandler = new Handler() {
            @Override
            public void publish(LogRecord aRecord) {
                if (aRecord.getLevel().intValue() >= Level.INFO.intValue()) {
                    aLogLinesList.add(tmpFormatter.formatMessage(aRecord));
                }
            }
            @Override
            public void flush() {
                //nothing to do
            }
            @Override
            public void close() {
                //nothing to do
            }
        };
        tmpLogger.addHandler(tmpHandler);
        try {
            return anImporter.importMoleculeFile(this.getResourceFile(aFileName), isRegardStereo, true,
                    isKekulizationEnforced, aSnapshot);
        } finally {
            tmpLogger.removeHandler(tmpHandler);
        }
    }

    /**
     * Creates a settings snapshot with the given values.
     *
     * @param anOption disconnected structures option
     * @param isElementFilterActive element filter
     * @param isHeavyAtomCountFilterActive heavy atom count filter
     * @param aMinimum minimum heavy atom count
     * @param aMaximum maximum heavy atom count
     * @param isDuplicateFilterActive duplicate filter
     * @return snapshot
     */
    private MoleculeFilterSettingsSnapshot createSnapshot(MoleculeFilterSettings.DisconnectedStructuresOption anOption,
                                                          boolean isElementFilterActive,
                                                          boolean isHeavyAtomCountFilterActive,
                                                          int aMinimum,
                                                          int aMaximum,
                                                          boolean isDuplicateFilterActive) {
        return new MoleculeFilterSettingsSnapshot(anOption, isElementFilterActive, isHeavyAtomCountFilterActive,
                aMinimum, aMaximum, isDuplicateFilterActive);
    }

    /**
     * Returns the names of the given molecules.
     *
     * @param aMoleculeList molecules
     * @return names in the same order
     */
    private List<String> getNames(List<MoleculeDataModel> aMoleculeList) {
        List<String> tmpNames = new ArrayList<>(aMoleculeList.size());
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            tmpNames.add(tmpMolecule.getName());
        }
        return tmpNames;
    }

    /**
     * Returns the unique SMILES strings of the given molecules.
     *
     * @param aMoleculeList molecules
     * @return SMILES strings in the same order
     */
    private List<String> getSmiles(List<MoleculeDataModel> aMoleculeList) {
        List<String> tmpSmiles = new ArrayList<>(aMoleculeList.size());
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            tmpSmiles.add(tmpMolecule.getUniqueSmiles());
        }
        return tmpSmiles;
    }

    /**
     * Returns the names of the twelve test molecules without the given ones, in file order.
     *
     * @param aRemovedNames names to leave out
     * @return remaining names
     */
    private List<String> withoutNames(String... aRemovedNames) {
        List<String> tmpNames = new ArrayList<>(ImporterMoleculeFilterTest.MOLECULE_NAMES);
        for (String tmpName : aRemovedNames) {
            Assertions.assertTrue(tmpNames.remove(tmpName), tmpName);
        }
        return Collections.unmodifiableList(tmpNames);
    }

    /**
     * Creates the log line the importer writes for a removed molecule.
     *
     * @param aName molecule name
     * @param aReason removal reason
     * @return expected log line
     */
    private String createRemovalLine(String aName, MoleculeFilterReason aReason) {
        return "Molecule " + aName + " was not imported, reason: " + aReason.name() + ".";
    }

    /**
     * Returns the log lines about removed molecules.
     *
     * @param aLogLinesList all captured log lines
     * @return lines about removed molecules, in logging order
     */
    private List<String> getRemovalLines(List<String> aLogLinesList) {
        List<String> tmpLines = new ArrayList<>();
        for (String tmpLine : aLogLinesList) {
            if (tmpLine.contains(" was not imported, reason: ")) {
                tmpLines.add(tmpLine);
            }
        }
        return tmpLines;
    }

    /**
     * Returns the log lines starting with the given prefix.
     *
     * @param aLogLinesList all captured log lines
     * @param aPrefix prefix
     * @return matching lines
     */
    private List<String> getLinesStartingWith(List<String> aLogLinesList, String aPrefix) {
        List<String> tmpLines = new ArrayList<>();
        for (String tmpLine : aLogLinesList) {
            if (tmpLine.startsWith(aPrefix)) {
                tmpLines.add(tmpLine);
            }
        }
        return tmpLines;
    }

    /**
     * Counts the log lines starting with the given prefix.
     *
     * @param aLogLinesList all captured log lines
     * @param aPrefix prefix
     * @return number of matching lines
     */
    private int countLinesStartingWith(List<String> aLogLinesList, String aPrefix) {
        return this.getLinesStartingWith(aLogLinesList, aPrefix).size();
    }

    /**
     * Checks the counts of the given summary; reasons not in the given map must have the count 0, and the counts must
     * add up (no pending molecules).
     *
     * @param aSummary the summary
     * @param aReadCount expected number of read molecules
     * @param aKeptCount expected number of kept molecules
     * @param aReducedCount expected number of reduced molecules
     * @param anExpectedRemovedCountsMap expected numbers of removed molecules per reason
     */
    private void assertSummary(MoleculeFilterSummary aSummary, int aReadCount, int aKeptCount, int aReducedCount,
                               Map<MoleculeFilterReason, Integer> anExpectedRemovedCountsMap) {
        Assertions.assertNotNull(aSummary);
        Assertions.assertEquals(aReadCount, aSummary.getReadCount(), aSummary.toString());
        Assertions.assertEquals(aKeptCount, aSummary.getKeptCount(), aSummary.toString());
        Assertions.assertEquals(aReducedCount, aSummary.getReducedToLargestComponentCount(), aSummary.toString());
        for (MoleculeFilterReason tmpReason : MoleculeFilterReason.values()) {
            Assertions.assertEquals(anExpectedRemovedCountsMap.getOrDefault(tmpReason, 0), aSummary.getRemovedCount(tmpReason),
                    tmpReason + ": " + aSummary);
        }
        Assertions.assertEquals(0, aSummary.getPendingCount(), aSummary.toString());
    }
    //</editor-fold>
}
