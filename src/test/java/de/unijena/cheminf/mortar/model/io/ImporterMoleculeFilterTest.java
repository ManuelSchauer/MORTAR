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
import org.junit.jupiter.api.io.TempDir;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * End-to-end tests of the molecule filters at import through the public
 * {@link Importer#importMoleculeFile(File, boolean, boolean, boolean, MoleculeFilterSettingsSnapshot)} method, with
 * small synthetic test files in this package:
 * <ul>
 *     <li>"MoleculeFiltersTestFile.sdf" and "MoleculeFiltersTestFile.smi", a dirty data set; both contain the same
 *     twelve molecules in the same order (position in brackets): ethanol (1), ethanol_duplicate (2, OCC, duplicate of
 *     1), L-alanine (3), D-alanine (4, stereoisomer of 3), triethylamine_hydrochloride (5, salt, CCN(CC)CC.Cl),
 *     triethylamine (6, parent compound of 5), sodium_acetate (7, salt with a disallowed counter-ion),
 *     paracetamol_hydrate (8), cisplatin (9, connected metal complex), chloride_ion (10, lone ion, one heavy atom),
 *     R_group_molecule (11, *CCO, pseudo atom), triacontane (12, 30 heavy atoms); the SD file additionally contains a
 *     13th record without atoms ("empty_record"), for which no SMILES string can be generated, and every SD record has
 *     the data fields "Comment" and "Source";</li>
 *     <li>"MoleculeFiltersDuplicatesTestFile.smi", toluene, pyridine, and naphthalene in aromatic and different Kekule
 *     notations and the alanine enantiomers;</li>
 *     <li>"MoleculeFiltersDuplicatesTestFile.sdf", toluene twice with different atom order and shifted double bonds,
 *     and paracetamol as Kekule structure and with the aromatic bond type (which sets aromaticity flags) in reversed
 *     atom order.</li>
 * </ul>
 * The tests run headless and do not touch the MORTAR settings directory.
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
     * Name of the SD test file resource with the dirty data set.
     */
    private static final String SD_FILE_NAME = "MoleculeFiltersTestFile.sdf";

    /**
     * Name of the SMILES test file resource with the dirty data set.
     */
    private static final String SMILES_FILE_NAME = "MoleculeFiltersTestFile.smi";

    /**
     * Name of the SD test file resource with duplicates in different notations.
     */
    private static final String DUPLICATES_SD_FILE_NAME = "MoleculeFiltersDuplicatesTestFile.sdf";

    /**
     * Name of the SMILES test file resource with duplicates in different notations.
     */
    private static final String DUPLICATES_SMILES_FILE_NAME = "MoleculeFiltersDuplicatesTestFile.smi";

    /**
     * Names of the twelve molecules contained in both dirty test files, in file order.
     */
    private static final List<String> MOLECULE_NAMES = List.of("ethanol", "ethanol_duplicate", "L-alanine",
            "D-alanine", "triethylamine_hydrochloride", "triethylamine", "sodium_acetate", "paracetamol_hydrate",
            "cisplatin", "chloride_ion", "R_group_molecule", "triacontane");

    /**
     * Expected SMILES strings of the twelve molecules of the dirty SD file imported with the default settings (filters
     * inactive, stereochemistry regarded, aromaticity encoded), identical to an import without molecule filters.
     */
    private static final List<String> DEFAULT_SD_FILE_SMILES = List.of("OCC", "OCC", "C[C@@H](C(=O)O)N",
            "C[C@H](C(=O)O)N", "Cl.N(CC)(CC)CC", "N(CC)(CC)CC", "[Na+].O=C([O-])C", "O=C(NC1=CC=C(O)C=C1)C.O",
            "Cl[Pt](Cl)(N)N", "[Cl-]", "*CCO", "C".repeat(30));

    /**
     * Expected SMILES strings of the twelve molecules of the dirty SMILES file imported with the default settings; the
     * same as for the SD file except for paracetamol, which has aromaticity flags from the SMILES notation.
     */
    private static final List<String> DEFAULT_SMILES_FILE_SMILES = List.of("OCC", "OCC", "C[C@@H](C(=O)O)N",
            "C[C@H](C(=O)O)N", "Cl.N(CC)(CC)CC", "N(CC)(CC)CC", "[Na+].O=C([O-])C", "O=C(Nc1ccc(O)cc1)C.O",
            "Cl[Pt](Cl)(N)N", "[Cl-]", "*CCO", "C".repeat(30));

    /**
     * Name of the additional record without atoms in the SD file.
     */
    private static final String EMPTY_RECORD_NAME = "empty_record";
    //</editor-fold>
    //
    //<editor-fold desc="Tests">
    /**
     * Tests that an import with the default (inactive) filter settings gives the expected molecules (names and SMILES
     * strings pinned literally, equal to an import without molecule filters), with the old method without filter
     * settings, the default snapshot, and the snapshot of a new settings container, for both files and both values of
     * the keep atom container setting. Also checks the neutral summary and that no molecule is logged as removed (the
     * empty SD record is not imported, as before, but not logged individually either).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDefaultSettingsResult() throws Exception {
        for (String tmpFileName : new String[] {ImporterMoleculeFilterTest.SD_FILE_NAME, ImporterMoleculeFilterTest.SMILES_FILE_NAME}) {
            boolean tmpIsSdFile = tmpFileName.equals(ImporterMoleculeFilterTest.SD_FILE_NAME);
            List<String> tmpExpectedSmiles = tmpIsSdFile ? ImporterMoleculeFilterTest.DEFAULT_SD_FILE_SMILES
                    : ImporterMoleculeFilterTest.DEFAULT_SMILES_FILE_SMILES;
            for (boolean tmpIsKeepAtomContainer : new boolean[] {false, true}) {
                Importer tmpOldImporter = this.createImporter(tmpIsKeepAtomContainer);
                List<MoleculeDataModel> tmpOldResult = tmpOldImporter.importMoleculeFile(this.getResourceFile(tmpFileName),
                        true, true, false);
                Importer tmpNewImporter = this.createImporter(tmpIsKeepAtomContainer);
                List<String> tmpLogLines = new ArrayList<>();
                List<MoleculeDataModel> tmpNewResult = this.importAndCaptureLog(tmpNewImporter, tmpFileName,
                        true, false, MoleculeFilterSettingsSnapshot.createDefault(), tmpLogLines);
                List<MoleculeDataModel> tmpSettingsResult = this.createImporter(tmpIsKeepAtomContainer).importMoleculeFile(
                        this.getResourceFile(tmpFileName), true, true, false, new MoleculeFilterSettings().createSnapshot());
                for (List<MoleculeDataModel> tmpResult : List.of(tmpOldResult, tmpNewResult, tmpSettingsResult)) {
                    Assertions.assertEquals(ImporterMoleculeFilterTest.MOLECULE_NAMES, this.getNames(tmpResult));
                    Assertions.assertEquals(tmpExpectedSmiles, this.getSmiles(tmpResult));
                    for (MoleculeDataModel tmpMolecule : tmpResult) {
                        Assertions.assertEquals(tmpIsKeepAtomContainer, tmpMolecule.isKeepAtomContainer());
                    }
                }
                Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = tmpIsSdFile
                        ? Map.of(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1) : Map.of();
                this.assertSummary(tmpNewImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 12, 0,
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

    /**
     * Tests the duplicate filter alone on the dirty data set: with stereochemistry, only the second ethanol is removed
     * (first occurrence kept, referenced in the log line); without stereochemistry, D-alanine is a duplicate of
     * L-alanine, too; the result does not depend on the Kekule setting. The salt and its parent compound are no
     * duplicates because the salt is not reduced. The SMILES strings of the kept molecules are the same as without the
     * filter. With an active filter, the empty SD record is logged as removed.
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
            List<String> tmpExpectedSmiles = new ArrayList<>(tmpIsSdFile ? ImporterMoleculeFilterTest.DEFAULT_SD_FILE_SMILES
                    : ImporterMoleculeFilterTest.DEFAULT_SMILES_FILE_SMILES);
            tmpExpectedSmiles.remove(1);
            Assertions.assertEquals(tmpExpectedSmiles, this.getSmiles(tmpResult));
            Map<MoleculeFilterReason, Integer> tmpExpectedRemovedCounts = new EnumMap<>(MoleculeFilterReason.class);
            tmpExpectedRemovedCounts.put(MoleculeFilterReason.DUPLICATE, 1);
            List<String> tmpExpectedRemovalLines = new ArrayList<>(List.of(
                    this.createDuplicateLine("ethanol_duplicate", 2, "ethanol", 1)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13,
                        MoleculeFilterReason.SMILES_GENERATION_FAILED));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 11, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(tmpExpectedRemovalLines, this.getRemovalLines(tmpLogLines));
            Assertions.assertEquals(1, this.countLinesStartingWith(tmpLogLines, "Successfully imported 11 molecules"));
            //Kekule structures enforced: the same molecules are duplicates
            tmpResult = this.createImporter(false).importMoleculeFile(this.getResourceFile(tmpFileName), true, true, true, tmpSnapshot);
            Assertions.assertEquals(this.withoutNames("ethanol_duplicate"), this.getNames(tmpResult));
            //without stereochemistry, with both Kekule settings
            for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
                tmpImporter = this.createImporter(false);
                tmpLogLines = new ArrayList<>();
                tmpResult = this.importAndCaptureLog(tmpImporter, tmpFileName, false, tmpIsKekulizationEnforced,
                        tmpSnapshot, tmpLogLines);
                Assertions.assertEquals(this.withoutNames("ethanol_duplicate", "D-alanine"), this.getNames(tmpResult));
                Assertions.assertEquals(tmpIsSdFile ? 13 : 12, tmpImporter.getMoleculeFilterSummary().getReadCount());
                Assertions.assertEquals(10, tmpImporter.getMoleculeFilterSummary().getKeptCount());
                Assertions.assertEquals(2, tmpImporter.getMoleculeFilterSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
                Assertions.assertEquals(0, tmpImporter.getMoleculeFilterSummary().getPendingCount());
                Assertions.assertEquals(List.of(this.createDuplicateLine("ethanol_duplicate", 2, "ethanol", 1),
                        this.createDuplicateLine("D-alanine", 4, "L-alanine", 3)),
                        this.getRemovalLines(tmpLogLines).subList(0, 2));
            }
        }
    }

    /**
     * Tests that the duplicate filter recognises one compound in different notations as duplicates, independently of
     * the atom order and of aromatic or Kekule notation and of the Kekule setting of the import, although the stored
     * SMILES strings of these notations differ by default: toluene, pyridine, and naphthalene in aromatic and different
     * Kekule notations in a SMILES file; toluene with permuted atom order and shifted double bonds and paracetamol as
     * Kekule structure and with aromatic bonds (aromaticity flags set) in an SD file. Stereoisomers stay distinct with
     * stereochemistry and become duplicates without. The stored SMILES strings of the kept molecules are the same as
     * without the filter.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDuplicateFilterIsIndependentOfNotation() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 0, 0, true);
        //SMILES file, without the filter: all imported, the stored SMILES strings of the toluene notations differ
        List<MoleculeDataModel> tmpUnfilteredResult = this.createImporter(false).importMoleculeFile(
                this.getResourceFile(ImporterMoleculeFilterTest.DUPLICATES_SMILES_FILE_NAME), true, true, false,
                MoleculeFilterSettingsSnapshot.createDefault());
        Assertions.assertEquals(9, tmpUnfilteredResult.size());
        Set<String> tmpTolueneSmiles = new HashSet<>(this.getSmiles(tmpUnfilteredResult.subList(0, 3)));
        Assertions.assertTrue(tmpTolueneSmiles.size() > 1, tmpTolueneSmiles.toString());
        //with the filter
        for (boolean tmpIsKeepAtomContainer : new boolean[] {false, true}) {
            for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
                Importer tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
                List<String> tmpLogLines = new ArrayList<>();
                List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter,
                        ImporterMoleculeFilterTest.DUPLICATES_SMILES_FILE_NAME, true, tmpIsKekulizationEnforced,
                        tmpSnapshot, tmpLogLines);
                Assertions.assertEquals(List.of("toluene_kekule", "pyridine_aromatic", "naphthalene_kekule", "L-alanine",
                        "D-alanine"), this.getNames(tmpResult));
                Assertions.assertEquals(List.of(
                        this.createDuplicateLine("toluene_kekule_shifted", 2, "toluene_kekule", 1),
                        this.createDuplicateLine("toluene_aromatic", 3, "toluene_kekule", 1),
                        this.createDuplicateLine("pyridine_kekule", 5, "pyridine_aromatic", 4),
                        this.createDuplicateLine("naphthalene_aromatic", 7, "naphthalene_kekule", 6)),
                        this.getRemovalLines(tmpLogLines));
                this.assertSummary(tmpImporter.getMoleculeFilterSummary(), 9, 5, 0, Map.of(MoleculeFilterReason.DUPLICATE, 4));
                if (!tmpIsKekulizationEnforced) {
                    //the stored SMILES strings are not affected by the filter
                    Assertions.assertEquals(List.of(tmpUnfilteredResult.get(0).getUniqueSmiles(),
                            tmpUnfilteredResult.get(3).getUniqueSmiles(), tmpUnfilteredResult.get(5).getUniqueSmiles(),
                            tmpUnfilteredResult.get(7).getUniqueSmiles(), tmpUnfilteredResult.get(8).getUniqueSmiles()),
                            this.getSmiles(tmpResult));
                }
            }
        }
        //without stereochemistry, the enantiomers are duplicates, too
        Importer tmpImporter = this.createImporter(false);
        List<String> tmpLogLines = new ArrayList<>();
        List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter,
                ImporterMoleculeFilterTest.DUPLICATES_SMILES_FILE_NAME, false, false, tmpSnapshot, tmpLogLines);
        Assertions.assertEquals(List.of("toluene_kekule", "pyridine_aromatic", "naphthalene_kekule", "L-alanine"),
                this.getNames(tmpResult));
        Assertions.assertEquals(this.createDuplicateLine("D-alanine", 9, "L-alanine", 8), this.getRemovalLines(tmpLogLines).getLast());
        //SD file, without the filter: all imported, the stored SMILES strings of the paracetamol records differ
        tmpUnfilteredResult = this.createImporter(false).importMoleculeFile(
                this.getResourceFile(ImporterMoleculeFilterTest.DUPLICATES_SD_FILE_NAME), true, true, false,
                MoleculeFilterSettingsSnapshot.createDefault());
        Assertions.assertEquals(List.of("toluene", "toluene_permuted", "paracetamol_kekule", "paracetamol_aromatic_bonds"),
                this.getNames(tmpUnfilteredResult));
        Assertions.assertNotEquals(tmpUnfilteredResult.get(2).getUniqueSmiles(), tmpUnfilteredResult.get(3).getUniqueSmiles());
        //with the filter
        for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
            tmpImporter = this.createImporter(false);
            tmpLogLines = new ArrayList<>();
            tmpResult = this.importAndCaptureLog(tmpImporter, ImporterMoleculeFilterTest.DUPLICATES_SD_FILE_NAME, true,
                    tmpIsKekulizationEnforced, tmpSnapshot, tmpLogLines);
            Assertions.assertEquals(List.of("toluene", "paracetamol_kekule"), this.getNames(tmpResult));
            Assertions.assertEquals(List.of(
                    this.createDuplicateLine("toluene_permuted", 2, "toluene", 1),
                    this.createDuplicateLine("paracetamol_aromatic_bonds", 4, "paracetamol_kekule", 3)),
                    this.getRemovalLines(tmpLogLines));
            Assertions.assertEquals("paracetamol, Kekule structure", tmpResult.get(1).getProperties().get("Comment"));
            if (!tmpIsKekulizationEnforced) {
                Assertions.assertEquals(List.of(tmpUnfilteredResult.get(0).getUniqueSmiles(),
                        tmpUnfilteredResult.get(2).getUniqueSmiles()), this.getSmiles(tmpResult));
            }
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
                    this.createRemovalLine("triethylamine_hydrochloride", 5, MoleculeFilterReason.DISCONNECTED_STRUCTURE),
                    this.createRemovalLine("sodium_acetate", 7, MoleculeFilterReason.DISCONNECTED_STRUCTURE),
                    this.createRemovalLine("paracetamol_hydrate", 8, MoleculeFilterReason.DISCONNECTED_STRUCTURE)));
            if (tmpIsSdFile) {
                //rejected as such because a structural filter is active
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13,
                        MoleculeFilterReason.EMPTY_MOLECULE));
            }
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), tmpIsSdFile ? 13 : 12, 9, 0, tmpExpectedRemovedCounts);
            Assertions.assertEquals(3, tmpImporter.getMoleculeFilterSummary().getRemovedByFiltersCount());
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
     * duplicate of it (first occurrence kept, i.e. the entry named like the salt). An import that only reduced
     * molecules (the empty SD record does not count as removed by the filters) has changed the data set.
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
            MoleculeFilterSummary tmpSummary = tmpImporter.getMoleculeFilterSummary();
            this.assertSummary(tmpSummary, 13, 12, 3, Map.of(MoleculeFilterReason.EMPTY_MOLECULE, 1));
            //only reduced, nothing removed by the filters: the data set was changed nevertheless
            Assertions.assertEquals(0, tmpSummary.getRemovedByFiltersCount());
            Assertions.assertTrue(tmpSummary.isDataSetChangedByFilters());
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
            //SMILES file, reduction only
            tmpImporter = this.createImporter(tmpIsKeepAtomContainer);
            tmpImporter.importMoleculeFile(this.getResourceFile(ImporterMoleculeFilterTest.SMILES_FILE_NAME), true, true,
                    false, tmpSnapshot);
            this.assertSummary(tmpImporter.getMoleculeFilterSummary(), 12, 12, 3, Map.of());
            Assertions.assertTrue(tmpImporter.getMoleculeFilterSummary().isDataSetChangedByFilters());
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
                    this.createDuplicateLine("ethanol_duplicate", 2, "ethanol", 1),
                    this.createDuplicateLine("triethylamine", 6, "triethylamine_hydrochloride", 5),
                    this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13, MoleculeFilterReason.EMPTY_MOLECULE)),
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
                    this.createRemovalLine("sodium_acetate", 7, MoleculeFilterReason.DISALLOWED_ELEMENT),
                    this.createRemovalLine("cisplatin", 9, MoleculeFilterReason.DISALLOWED_ELEMENT),
                    this.createRemovalLine("R_group_molecule", 11, MoleculeFilterReason.DISALLOWED_ELEMENT)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13,
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
                    this.createRemovalLine("chloride_ion", 10, MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE),
                    this.createRemovalLine("triacontane", 12, MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE)));
            if (tmpIsSdFile) {
                tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13,
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
     * kept molecules, their SMILES strings in the SMILES file, the summary, the display text for the summary dialog
     * (the empty SD record is listed although it does not count as removed by the filters), and exactly one log line
     * per removed molecule with its name, position, and reason, in file order.
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
                        this.createDuplicateLine("ethanol_duplicate", 2, "ethanol", 1),
                        this.createDuplicateLine("triethylamine", 6, "triethylamine_hydrochloride", 5),
                        this.createRemovalLine("cisplatin", 9, MoleculeFilterReason.DISALLOWED_ELEMENT),
                        this.createRemovalLine("chloride_ion", 10, MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE),
                        this.createRemovalLine("R_group_molecule", 11, MoleculeFilterReason.DISALLOWED_ELEMENT),
                        this.createRemovalLine("triacontane", 12, MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE)));
                if (tmpIsSdFile) {
                    tmpExpectedRemovedCounts.put(MoleculeFilterReason.EMPTY_MOLECULE, 1);
                    tmpExpectedRemovalLines.add(this.createRemovalLine(ImporterMoleculeFilterTest.EMPTY_RECORD_NAME, 13,
                            MoleculeFilterReason.EMPTY_MOLECULE));
                }
                MoleculeFilterSummary tmpSummary = tmpImporter.getMoleculeFilterSummary();
                this.assertSummary(tmpSummary, tmpIsSdFile ? 13 : 12, 6, 3, tmpExpectedRemovedCounts);
                Assertions.assertEquals(tmpSnapshot, tmpSummary.getSettingsSnapshot());
                Assertions.assertTrue(tmpSummary.isDataSetChangedByFilters());
                //the empty record does not count as removed by the filters
                Assertions.assertEquals(6, tmpSummary.getRemovedByFiltersCount());
                //exactly one line per removed molecule, with its name, position, and reason
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
     * summary reports all molecules (except the empty SD record) as removed by the filters, so that the GUI can tell
     * the user. Also checks that the message texts of the GUI dialog exist.
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
            Assertions.assertEquals(12, tmpSummary.getRemovedByFiltersCount());
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
     * Tests that empty records do not count as removed by the filters: an SD file with only empty records, imported
     * with an active structural filter, gives an empty list, and the summary lists the empty molecules but reports no
     * change by the filters, so that the GUI shows the usual message for an empty import instead of the one for an
     * import emptied by the filters.
     *
     * @param aTempDir temporary directory provided by JUnit
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testEmptyRecordsAreNoChangeByFilters(@TempDir Path aTempDir) throws Exception {
        String tmpEmptyRecord = "empty\n  MORTAR test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n";
        File tmpFile = aTempDir.resolve("EmptyRecords.sdf").toFile();
        Files.writeString(tmpFile.toPath(), tmpEmptyRecord + tmpEmptyRecord, StandardCharsets.US_ASCII);
        Importer tmpImporter = this.createImporter(false);
        List<MoleculeDataModel> tmpResult = tmpImporter.importMoleculeFile(tmpFile, true, true, false, this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE, true, true, 2, 100, true));
        Assertions.assertTrue(tmpResult.isEmpty());
        MoleculeFilterSummary tmpSummary = tmpImporter.getMoleculeFilterSummary();
        this.assertSummary(tmpSummary, 2, 0, 0, Map.of(MoleculeFilterReason.EMPTY_MOLECULE, 2));
        Assertions.assertEquals(0, tmpSummary.getRemovedByFiltersCount());
        Assertions.assertFalse(tmpSummary.isDataSetChangedByFilters());
        Assertions.assertTrue(tmpSummary.getDisplayText().endsWith("  - empty molecule: 2"), tmpSummary.getDisplayText());
    }

    /**
     * Tests the log lines of removed molecules: molecules with the same name are distinguished by their position, a
     * duplicate refers to name and position of the first occurrence, and at most
     * {@link Importer#MAXIMUM_NUMBER_OF_INDIVIDUALLY_LOGGED_REMOVALS} removals are logged individually per import,
     * followed by one line with the number of further removals; with exactly this number of removals, no such line is
     * written.
     *
     * @param aTempDir temporary directory provided by JUnit
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testRemovalLogLines(@TempDir Path aTempDir) throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 0, 0, true);
        //equal names
        File tmpSameNamesFile = aTempDir.resolve("SameNames.smi").toFile();
        Files.writeString(tmpSameNamesFile.toPath(), "CCO\tnan\nCCN\tnan\nOCC\tnan\nNCC\tnan\n", StandardCharsets.US_ASCII);
        List<String> tmpLogLines = new ArrayList<>();
        this.importAndCaptureLog(this.createImporter(false), tmpSameNamesFile, true, false, tmpSnapshot, tmpLogLines);
        Assertions.assertEquals(List.of(
                "Molecule nan (no. 3 of the read molecules) was not imported, reason: DUPLICATE of molecule nan (no. 1).",
                "Molecule nan (no. 4 of the read molecules) was not imported, reason: DUPLICATE of molecule nan (no. 2)."),
                this.getRemovalLines(tmpLogLines));
        //the cap
        int tmpMaximum = Importer.MAXIMUM_NUMBER_OF_INDIVIDUALLY_LOGGED_REMOVALS;
        for (int tmpRemovedCount : new int[] {tmpMaximum, tmpMaximum + 2}) {
            StringBuilder tmpContent = new StringBuilder();
            for (int i = 1; i <= tmpRemovedCount + 1; i++) {
                tmpContent.append("CCO\tmol").append(i).append("\n");
            }
            File tmpFile = aTempDir.resolve("Ethanol" + tmpRemovedCount + ".smi").toFile();
            Files.writeString(tmpFile.toPath(), tmpContent.toString(), StandardCharsets.US_ASCII);
            Importer tmpImporter = this.createImporter(false);
            tmpLogLines = new ArrayList<>();
            List<MoleculeDataModel> tmpResult = this.importAndCaptureLog(tmpImporter, tmpFile, true, false, tmpSnapshot, tmpLogLines);
            Assertions.assertEquals(List.of("mol1"), this.getNames(tmpResult));
            Assertions.assertEquals(tmpRemovedCount, tmpImporter.getMoleculeFilterSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
            List<String> tmpRemovalLines = this.getRemovalLines(tmpLogLines);
            Assertions.assertEquals(tmpMaximum, tmpRemovalLines.size());
            Assertions.assertEquals(this.createDuplicateLine("mol2", 2, "mol1", 1), tmpRemovalLines.getFirst());
            Assertions.assertEquals(this.createDuplicateLine("mol" + (tmpMaximum + 1), tmpMaximum + 1, "mol1", 1),
                    tmpRemovalLines.getLast());
            List<String> tmpFurtherLines = new ArrayList<>();
            for (String tmpLine : tmpLogLines) {
                if (tmpLine.contains("further molecules were not imported")) {
                    tmpFurtherLines.add(tmpLine);
                }
            }
            if (tmpRemovedCount > tmpMaximum) {
                Assertions.assertEquals(List.of("2 further molecules were not imported but are not listed individually "
                        + "(at most " + tmpMaximum + " per import); see the summary of the molecule filters below."), tmpFurtherLines);
            } else {
                Assertions.assertEquals(List.of(), tmpFurtherLines);
            }
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
     * Imports the given test resource file with filling of open valences and the given settings and collects the
     * formatted messages logged by the importer meanwhile.
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
        return this.importAndCaptureLog(anImporter, this.getResourceFile(aFileName), isRegardStereo,
                isKekulizationEnforced, aSnapshot, aLogLinesList);
    }

    /**
     * Imports the given file with filling of open valences and the given settings and collects the formatted messages
     * logged by the importer meanwhile.
     *
     * @param anImporter the importer to use
     * @param aFile the file
     * @param isRegardStereo stereochemistry setting
     * @param isKekulizationEnforced Kekule setting
     * @param aSnapshot molecule filter settings
     * @param aLogLinesList list to add the logged messages to
     * @return imported molecules
     * @throws Exception if anything goes wrong
     */
    private List<MoleculeDataModel> importAndCaptureLog(Importer anImporter, File aFile, boolean isRegardStereo,
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
            return anImporter.importMoleculeFile(aFile, isRegardStereo, true, isKekulizationEnforced, aSnapshot);
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
     * Returns the names of the twelve molecules of the dirty test files without the given ones, in file order.
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
     * Creates the log line the importer writes for a molecule removed for a reason other than being a duplicate.
     *
     * @param aName molecule name
     * @param aPosition 1-based position of the molecule among the read molecules
     * @param aReason removal reason
     * @return expected log line
     */
    private String createRemovalLine(String aName, int aPosition, MoleculeFilterReason aReason) {
        return "Molecule " + aName + " (no. " + aPosition + " of the read molecules) was not imported, reason: "
                + aReason.name() + ".";
    }

    /**
     * Creates the log line the importer writes for a molecule removed as duplicate.
     *
     * @param aName molecule name
     * @param aPosition 1-based position of the molecule among the read molecules
     * @param aFirstOccurrenceName name of the first occurrence
     * @param aFirstOccurrencePosition position of the first occurrence
     * @return expected log line
     */
    private String createDuplicateLine(String aName, int aPosition, String aFirstOccurrenceName, int aFirstOccurrencePosition) {
        return "Molecule " + aName + " (no. " + aPosition + " of the read molecules) was not imported, reason: DUPLICATE "
                + "of molecule " + aFirstOccurrenceName + " (no. " + aFirstOccurrencePosition + ").";
    }

    /**
     * Returns the log lines about individual removed molecules.
     *
     * @param aLogLinesList all captured log lines
     * @return lines about removed molecules, in logging order
     */
    private List<String> getRemovalLines(List<String> aLogLinesList) {
        List<String> tmpLines = new ArrayList<>();
        for (String tmpLine : aLogLinesList) {
            if (tmpLine.contains(" of the read molecules) was not imported, reason: ")) {
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
