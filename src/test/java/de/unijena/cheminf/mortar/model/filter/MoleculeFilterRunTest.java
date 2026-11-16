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

package de.unijena.cheminf.mortar.model.filter;

import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Tests the {@link MoleculeFilterRun} on a small synthetic "dirty" data set that is run through the importer's
 * preprocessing and SMILES generation in the same way as the importer is going to use the filters: neutral defaults,
 * each filter alone, the fixed order of application and its consequences, the duplicate check with the stereochemistry
 * and Kekule settings, odd input, and the summary counts.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class MoleculeFilterRunTest {
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
     * Synthetic dirty data set as SMILES, in the order of a hypothetical input file; the indices are referred to in the
     * tests.
     */
    private static final String[] DIRTY_SET = {
            "CCO", //0 ethanol
            "OCC", //1 ethanol in another notation, duplicate of 0
            "CC(=O)[O-].[Na+]", //2 sodium acetate
            "CC(=O)O", //3 acetic acid
            "CC(=O)[O-].[K+]", //4 potassium acetate
            "CCO.O", //5 ethanol hydrate
            "c1ccccc1.C1CCCCC1", //6 1:1 mixture of benzene and cyclohexane (tie, 6 heavy atoms each)
            "[Na+].[Cl-]", //7 sodium chloride (tie, 1 heavy atom each)
            "[Na+]", //8 lone ion
            "Cl[Pt](Cl)(N)N", //9 cisplatin, connected metal complex
            "*CCO", //10 R group (pseudo atom)
            "[2H]C([2H])([2H])O", //11 deuterated methanol, explicit deuterium atoms
            "[13CH3]CO", //12 carbon-13 labelled ethanol
            "C[C@H](N)C(=O)O", //13 L-alanine
            "C[C@@H](N)C(=O)O", //14 D-alanine
            "c1ccccc1O", //15 phenol, aromatic notation
            "C1=CC=CC=C1O", //16 phenol, Kekule notation
            "C".repeat(120), //17 huge molecule, 120 heavy atoms
            "", //18 empty molecule
            "C", //19 methane, 1 heavy atom
            "C[N+](C)(C)C.Cc1ccc(cc1)S(=O)(=O)[O-]", //20 tetramethylammonium tosylate, counter-ion is larger
            "OB(O)c1ccccc1", //21 phenylboronic acid
            "C[Se]C", //22 dimethyl selenide
            "C[As](C)C" //23 trimethylarsine, metalloid that is not allowed
    };
    //</editor-fold>
    //
    //<editor-fold desc="Tests">
    /**
     * Tests that with the default settings, every molecule is passed through unchanged (same container object, not
     * reduced) and the import result is identical to an import without molecule filters, for all combinations of the
     * stereochemistry and Kekule settings; the summary counts only the molecule without SMILES (the empty one).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDefaultSettingsAreNeutral() throws Exception {
        for (boolean tmpIsRegardStereo : new boolean[] {true, false}) {
            for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
                List<IAtomContainer> tmpMolecules = this.preprocessDirtySet();
                List<String> tmpExpectedSmiles = new ArrayList<>(tmpMolecules.size());
                for (IAtomContainer tmpMolecule : tmpMolecules) {
                    //the importer's parse loop without filters
                    String tmpSmiles = ChemUtil.createUniqueSmiles(tmpMolecule, tmpIsRegardStereo, !tmpIsKekulizationEnforced);
                    if (tmpSmiles != null && !tmpSmiles.isBlank()) {
                        tmpExpectedSmiles.add(tmpSmiles);
                    }
                }
                MoleculeFilterRun tmpRun = new MoleculeFilterRun(new MoleculeFilterSettings());
                for (IAtomContainer tmpMolecule : tmpMolecules) {
                    MoleculeFilterResult tmpResult = new MoleculeFilterRun(MoleculeFilterSettingsSnapshot.createDefault())
                            .applyStructuralFilters(tmpMolecule);
                    Assertions.assertTrue(tmpResult.isKept());
                    Assertions.assertSame(tmpMolecule, tmpResult.getMolecule());
                    Assertions.assertFalse(tmpResult.isReducedToLargestComponent());
                    Assertions.assertNull(tmpResult.getRejectionReason());
                }
                List<String> tmpKeptSmiles = this.runImport(tmpRun, tmpMolecules, tmpIsRegardStereo, tmpIsKekulizationEnforced);
                Assertions.assertEquals(tmpExpectedSmiles, tmpKeptSmiles);
                MoleculeFilterSummary tmpSummary = tmpRun.getSummary();
                Assertions.assertEquals(MoleculeFilterRunTest.DIRTY_SET.length, tmpSummary.getReadCount());
                Assertions.assertEquals(MoleculeFilterRunTest.DIRTY_SET.length - 1, tmpSummary.getKeptCount());
                Assertions.assertEquals(0, tmpSummary.getReducedToLargestComponentCount());
                this.assertRemovedCounts(tmpSummary, Map.of(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1));
            }
        }
    }

    /**
     * Tests the disconnected structures setting "remove molecule" alone: salts, hydrates, mixtures, and ion pairs are
     * removed, all other molecules pass. The empty molecule is rejected as such because a structural filter is active.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testRemoveDisconnectedStructures() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE, false, false, 0, 0, false));
        List<IAtomContainer> tmpMolecules = this.preprocessDirtySet();
        for (int tmpIndex : new int[] {2, 4, 5, 6, 7, 20}) {
            MoleculeFilterResult tmpResult = tmpRun.applyStructuralFilters(tmpMolecules.get(tmpIndex));
            Assertions.assertFalse(tmpResult.isKept());
            Assertions.assertNull(tmpResult.getMolecule());
            Assertions.assertEquals(MoleculeFilterReason.DISCONNECTED_STRUCTURE, tmpResult.getRejectionReason());
        }
        for (int tmpIndex : new int[] {0, 3, 8, 9, 10, 17, 19}) {
            Assertions.assertSame(tmpMolecules.get(tmpIndex), tmpRun.applyStructuralFilters(tmpMolecules.get(tmpIndex)).getMolecule());
        }
        MoleculeFilterRun tmpFullRun = new MoleculeFilterRun(tmpRun.getSettingsSnapshot());
        List<String> tmpKeptSmiles = this.runImport(tmpFullRun, this.preprocessDirtySet(), true, false);
        Assertions.assertEquals(17, tmpKeptSmiles.size());
        this.assertSummary(tmpFullRun.getSummary(), 24, 17, 0, Map.of(
                MoleculeFilterReason.DISCONNECTED_STRUCTURE, 6,
                MoleculeFilterReason.EMPTY_MOLECULE, 1));
    }

    /**
     * Tests the disconnected structures setting "keep largest component" alone: salts, hydrates, the 1:1 mixture (tie),
     * the ion pair (tie), and the salt with the larger counter-ion are reduced; the reduced molecules keep their names.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testKeepLargestComponent() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, false, false, 0, 0, false));
        List<IAtomContainer> tmpMolecules = this.preprocessDirtySet();
        Map<Integer, String> tmpExpectedComponents = Map.of(
                2, "CC(=O)[O-]",
                4, "CC(=O)[O-]",
                5, "CCO",
                6, "C1CCCCC1",
                7, "[Cl-]",
                //limitation: the larger counter-ion is chosen
                20, "Cc1ccc(cc1)S(=O)(=O)[O-]");
        for (Map.Entry<Integer, String> tmpEntry : tmpExpectedComponents.entrySet()) {
            IAtomContainer tmpMolecule = tmpMolecules.get(tmpEntry.getKey());
            MoleculeFilterResult tmpResult = tmpRun.applyStructuralFilters(tmpMolecule);
            Assertions.assertTrue(tmpResult.isKept());
            Assertions.assertTrue(tmpResult.isReducedToLargestComponent());
            Assertions.assertNotSame(tmpMolecule, tmpResult.getMolecule());
            Assertions.assertEquals("mol" + tmpEntry.getKey(), tmpResult.getMolecule().getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY));
            Assertions.assertEquals(this.uniqueSmiles(new PreprocessingTestImporter().parseAndPreprocess(tmpEntry.getValue()).getFirst()),
                    this.uniqueSmiles(tmpResult.getMolecule()));
        }
        MoleculeFilterResult tmpConnectedResult = tmpRun.applyStructuralFilters(tmpMolecules.get(0));
        Assertions.assertSame(tmpMolecules.get(0), tmpConnectedResult.getMolecule());
        Assertions.assertFalse(tmpConnectedResult.isReducedToLargestComponent());
        MoleculeFilterRun tmpFullRun = new MoleculeFilterRun(tmpRun.getSettingsSnapshot());
        List<String> tmpKeptSmiles = this.runImport(tmpFullRun, this.preprocessDirtySet(), true, false);
        Assertions.assertEquals(23, tmpKeptSmiles.size());
        this.assertSummary(tmpFullRun.getSummary(), 24, 23, 6, Map.of(MoleculeFilterReason.EMPTY_MOLECULE, 1));
    }

    /**
     * Tests the element filter alone (disconnected structures unchanged): metals, the metal complex, the pseudo atom,
     * and arsenic are removed, also the salts because the sodium and potassium ions are still part of them;
     * deuterium, carbon-13, boron, and selenium pass.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testElementFilter() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, true, false, 0, 0, false));
        List<IAtomContainer> tmpMolecules = this.preprocessDirtySet();
        for (int tmpIndex = 0; tmpIndex < tmpMolecules.size(); tmpIndex++) {
            MoleculeFilterResult tmpResult = tmpRun.applyStructuralFilters(tmpMolecules.get(tmpIndex));
            switch (tmpIndex) {
                case 2, 4, 7, 8, 9, 10, 23 -> Assertions.assertEquals(MoleculeFilterReason.DISALLOWED_ELEMENT,
                        tmpResult.getRejectionReason(), MoleculeFilterRunTest.DIRTY_SET[tmpIndex]);
                case 18 -> Assertions.assertEquals(MoleculeFilterReason.EMPTY_MOLECULE, tmpResult.getRejectionReason());
                default -> Assertions.assertSame(tmpMolecules.get(tmpIndex), tmpResult.getMolecule(),
                        MoleculeFilterRunTest.DIRTY_SET[tmpIndex]);
            }
        }
        //explicit (non-isotopic) hydrogen atoms that are not removed by the preprocessing, e.g. in H2
        IAtomContainer tmpHydrogen = new PreprocessingTestImporter().parseAndPreprocess("[H][H]").getFirst();
        Assertions.assertEquals(2, tmpHydrogen.getAtomCount());
        Assertions.assertTrue(tmpRun.applyStructuralFilters(tmpHydrogen).isKept());
        this.assertSummary(tmpRun.getSummary(), 25, 0, 0, Map.of(
                MoleculeFilterReason.DISALLOWED_ELEMENT, 7,
                MoleculeFilterReason.EMPTY_MOLECULE, 1));
    }

    /**
     * Tests the heavy atom count filter alone: inclusive bounds, lone ions and molecules without heavy atoms, explicit
     * deuterium atoms and pseudo atoms not counted, disconnected molecules counted as a whole if not reduced, and the
     * default range when the filter is activated.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testHeavyAtomCountFilter() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 2, 3, false));
        String[] tmpKeptSmiles = {"CC", "CCO", "[2H]C([2H])([2H])O", "*CCO", "[2H][2H].CC"};
        String[] tmpRemovedSmiles = {"C", "CCCC", "[Na+]", "[H+]", "[H][H]", "CC(=O)[O-].[Na+]", "CCO.O"};
        for (IAtomContainer tmpMolecule : new PreprocessingTestImporter().parseAndPreprocess(tmpKeptSmiles)) {
            Assertions.assertTrue(tmpRun.applyStructuralFilters(tmpMolecule).isKept());
        }
        for (IAtomContainer tmpMolecule : new PreprocessingTestImporter().parseAndPreprocess(tmpRemovedSmiles)) {
            Assertions.assertEquals(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE,
                    tmpRun.applyStructuralFilters(tmpMolecule).getRejectionReason());
        }
        this.assertSummary(tmpRun.getSummary(), 12, 0, 0, Map.of(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE, 7));
        //lower bound 0 keeps molecules without heavy atoms
        MoleculeFilterRun tmpZeroRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 0, 0, false));
        List<IAtomContainer> tmpSmallMolecules = new PreprocessingTestImporter().parseAndPreprocess("[H+]", "[H][H]", "C");
        Assertions.assertTrue(tmpZeroRun.applyStructuralFilters(tmpSmallMolecules.get(0)).isKept());
        Assertions.assertTrue(tmpZeroRun.applyStructuralFilters(tmpSmallMolecules.get(1)).isKept());
        Assertions.assertFalse(tmpZeroRun.applyStructuralFilters(tmpSmallMolecules.get(2)).isKept());
        //the default range [2, 100] with the filter activated
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        tmpSettings.setHeavyAtomCountFilterSetting(true);
        MoleculeFilterRun tmpDefaultRangeRun = new MoleculeFilterRun(tmpSettings);
        List<IAtomContainer> tmpMolecules = new PreprocessingTestImporter().parseAndPreprocess(
                "C".repeat(100), "C".repeat(101), "[Na+]", "CC");
        Assertions.assertTrue(tmpDefaultRangeRun.applyStructuralFilters(tmpMolecules.get(0)).isKept());
        Assertions.assertFalse(tmpDefaultRangeRun.applyStructuralFilters(tmpMolecules.get(1)).isKept());
        Assertions.assertFalse(tmpDefaultRangeRun.applyStructuralFilters(tmpMolecules.get(2)).isKept());
        Assertions.assertTrue(tmpDefaultRangeRun.applyStructuralFilters(tmpMolecules.get(3)).isKept());
    }

    /**
     * Tests the duplicate filter alone on the dirty set, with the duplicate key the importer uses: different notations
     * of one compound are duplicates, also the aromatic and the Kekule notation of phenol, independently of the Kekule
     * setting; the first occurrence is kept; stereoisomers are only distinct if stereochemistry is regarded; mass
     * numbers are not encoded, so the carbon-13 labelled ethanol is a duplicate of ethanol, while deuterated methanol
     * keeps its explicit hydrogen atoms and differs from methanol. Salts are no duplicates of each other without
     * reduction. Also tests the first occurrence returned by the run.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testDuplicateFilter() throws Exception {
        MoleculeFilterSettingsSnapshot tmpSnapshot = this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 0, 0, true);
        for (boolean tmpIsKekulizationEnforced : new boolean[] {false, true}) {
            //stereo regarded (MORTAR default)
            MoleculeFilterRun tmpRun = new MoleculeFilterRun(tmpSnapshot);
            List<String> tmpKeptSmiles = this.runImport(tmpRun, this.preprocessDirtySet(), true, tmpIsKekulizationEnforced);
            //duplicates: 1 (of 0), 12 (of 0), and 16 (of 15); the empty molecule has no SMILES
            this.assertSummary(tmpRun.getSummary(), 24, 20, 0, Map.of(
                    MoleculeFilterReason.DUPLICATE, 3,
                    MoleculeFilterReason.SMILES_GENERATION_FAILED, 1));
            Assertions.assertEquals(tmpKeptSmiles.size(), Set.copyOf(tmpKeptSmiles).size());
            //first occurrence kept: ethanol is the first kept molecule
            Assertions.assertEquals(ChemUtil.createUniqueSmiles(this.preprocessDirtySet().getFirst(), true, !tmpIsKekulizationEnforced),
                    tmpKeptSmiles.getFirst());
            //stereo ignored: D-alanine becomes a duplicate of L-alanine
            tmpRun = new MoleculeFilterRun(tmpSnapshot);
            this.runImport(tmpRun, this.preprocessDirtySet(), false, tmpIsKekulizationEnforced);
            Assertions.assertEquals(4, tmpRun.getSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
        }
        //the check itself is a plain string comparison and remembers the first occurrence
        MoleculeFilterRun tmpStringRun = new MoleculeFilterRun(tmpSnapshot);
        Assertions.assertNull(tmpStringRun.checkDuplicate("CCO", "ethanol", 1));
        Assertions.assertNull(tmpStringRun.checkDuplicate("OCC", "ethanol in another notation", 2));
        Assertions.assertEquals(new MoleculeFilterRun.FirstOccurrence("ethanol", 1), tmpStringRun.checkDuplicate("CCO", "nan", 3));
        Assertions.assertEquals(new MoleculeFilterRun.FirstOccurrence("ethanol", 1), tmpStringRun.checkDuplicate("CCO", null, 0));
        Assertions.assertTrue(tmpStringRun.isDuplicate("CCO"));
        Assertions.assertFalse(tmpStringRun.isDuplicate("CCN"));
        Assertions.assertEquals(new MoleculeFilterRun.FirstOccurrence(null, 0), tmpStringRun.checkDuplicate("CCN", "ethylamine", 5));
        Assertions.assertThrows(NullPointerException.class, () -> tmpStringRun.isDuplicate(null));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpStringRun.isDuplicate(" "));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpStringRun.checkDuplicate("CCC", "propane", -1));
        Assertions.assertEquals(3, tmpStringRun.getSummary().getKeptCount());
        Assertions.assertEquals(4, tmpStringRun.getSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
        //inactive filter: never a duplicate, nothing registered
        MoleculeFilterRun tmpInactiveRun = new MoleculeFilterRun(MoleculeFilterSettingsSnapshot.createDefault());
        Assertions.assertFalse(tmpInactiveRun.isDuplicate("CCO"));
        Assertions.assertNull(tmpInactiveRun.checkDuplicate("CCO", "ethanol", 2));
        Assertions.assertEquals(2, tmpInactiveRun.getSummary().getKeptCount());
    }

    /**
     * Tests the consequences of the fixed order of application: with "keep largest component", sodium acetate is
     * judged by the element and heavy atom count filters on the acetate only, and two salt forms of acetic acid
     * become duplicates of each other; without reduction, the sodium ion leads to removal by the element filter, counts
     * as heavy atom, and the salt forms are distinct.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testOrderOfApplication() throws Exception {
        PreprocessingTestImporter tmpImporter = new PreprocessingTestImporter();
        String[] tmpSmiles = {"CC(=O)[O-].[Na+]", "CC(=O)O", "CC(=O)[O-].[K+]"};
        //reduction before the element filter
        MoleculeFilterRun tmpReducingRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, false, 0, 0, false));
        MoleculeFilterRun tmpUnchangedRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, true, false, 0, 0, false));
        List<IAtomContainer> tmpMolecules = tmpImporter.parseAndPreprocess(tmpSmiles);
        Assertions.assertTrue(tmpReducingRun.applyStructuralFilters(tmpMolecules.get(0)).isKept());
        Assertions.assertTrue(tmpReducingRun.applyStructuralFilters(tmpMolecules.get(1)).isKept());
        tmpMolecules = tmpImporter.parseAndPreprocess(tmpSmiles);
        Assertions.assertEquals(MoleculeFilterReason.DISALLOWED_ELEMENT,
                tmpUnchangedRun.applyStructuralFilters(tmpMolecules.get(0)).getRejectionReason());
        Assertions.assertTrue(tmpUnchangedRun.applyStructuralFilters(tmpMolecules.get(1)).isKept());
        //reduction before the heavy atom count filter: acetate has 4, sodium acetate 5 heavy atoms
        tmpReducingRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, false, true, 5, 5, false));
        tmpUnchangedRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, true, 5, 5, false));
        Assertions.assertEquals(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE,
                tmpReducingRun.applyStructuralFilters(tmpImporter.parseAndPreprocess(tmpSmiles[0]).getFirst()).getRejectionReason());
        Assertions.assertTrue(tmpUnchangedRun.applyStructuralFilters(tmpImporter.parseAndPreprocess(tmpSmiles[0]).getFirst()).isKept());
        //reduction before the duplicate filter
        tmpReducingRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, false, false, 0, 0, true));
        tmpUnchangedRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED, false, false, 0, 0, true));
        List<String> tmpReducedKept = this.runImport(tmpReducingRun, tmpImporter.parseAndPreprocess(tmpSmiles), true, false);
        List<String> tmpUnchangedKept = this.runImport(tmpUnchangedRun, tmpImporter.parseAndPreprocess(tmpSmiles), true, false);
        Assertions.assertEquals(2, tmpReducedKept.size());
        Assertions.assertEquals(1, tmpReducingRun.getSummary().getRemovedCount(MoleculeFilterReason.DUPLICATE));
        Assertions.assertEquals(3, tmpUnchangedKept.size());
        //a molecule can be reduced and then removed by a later filter
        MoleculeFilterRun tmpAllRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, true, 2, 100, true));
        MoleculeFilterResult tmpResult = tmpAllRun.applyStructuralFilters(
                tmpImporter.parseAndPreprocess("CC[Sn](CC)(CC)CC.O").getFirst());
        Assertions.assertEquals(MoleculeFilterReason.DISALLOWED_ELEMENT, tmpResult.getRejectionReason());
        Assertions.assertTrue(tmpResult.isReducedToLargestComponent());
        Assertions.assertEquals(1, tmpAllRun.getSummary().getReducedToLargestComponentCount());
    }

    /**
     * Tests all filters together on the dirty set and that the summary counts add up (read = kept + removed). Expected
     * per molecule: kept 0, 2 (as acetate), 3, 6 (as cyclohexane), 11, 13, 14, 15, 20 (as tosylate), 21, 22;
     * duplicates 1, 4, 5 (as ethanol), 12, 16 (Kekule notation of 15); disallowed elements 8, 9, 10, 23; heavy atom
     * count 7 (as chloride), 17, 19; empty 18; reduced 2, 4, 5, 6, 7, 20.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testAllFiltersAndSummary() throws Exception {
        MoleculeFilterSettings tmpSettings = new MoleculeFilterSettings();
        tmpSettings.setDisconnectedStructuresSetting(MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT);
        tmpSettings.setElementFilterSetting(true);
        tmpSettings.setHeavyAtomCountFilterSetting(true);
        tmpSettings.setDuplicateFilterSetting(true);
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(tmpSettings);
        //a settings change during the run has no effect
        tmpSettings.restoreDefaultSettings();
        List<String> tmpKeptSmiles = this.runImport(tmpRun, this.preprocessDirtySet(), true, false);
        MoleculeFilterSummary tmpSummary = tmpRun.getSummary();
        this.assertSummary(tmpSummary, 24, 11, 6, Map.of(
                MoleculeFilterReason.DUPLICATE, 5,
                MoleculeFilterReason.DISALLOWED_ELEMENT, 4,
                MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE, 3,
                MoleculeFilterReason.EMPTY_MOLECULE, 1));
        Assertions.assertEquals(11, tmpKeptSmiles.size());
        List<IAtomContainer> tmpExpectedMolecules = new PreprocessingTestImporter().parseAndPreprocess(
                "CCO", "CC(=O)[O-]", "CC(=O)O", "C1CCCCC1", "[2H]C([2H])([2H])O", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O",
                "c1ccccc1O", "Cc1ccc(cc1)S(=O)(=O)[O-]", "OB(O)c1ccccc1", "C[Se]C");
        List<String> tmpExpectedSmiles = new ArrayList<>(tmpExpectedMolecules.size());
        for (IAtomContainer tmpMolecule : tmpExpectedMolecules) {
            tmpExpectedSmiles.add(this.uniqueSmiles(tmpMolecule));
        }
        Assertions.assertEquals(tmpExpectedSmiles, tmpKeptSmiles);
        String tmpLogMessage = tmpSummary.toString();
        Assertions.assertTrue(tmpLogMessage.contains("24 molecules read, 11 kept, 13 removed"), tmpLogMessage);
        //the empty molecule does not count as removed by the filters
        Assertions.assertEquals(12, tmpSummary.getRemovedByFiltersCount());
        Assertions.assertTrue(tmpLogMessage.contains("6 reduced to their largest component"), tmpLogMessage);
        Assertions.assertTrue(tmpLogMessage.contains("heavy atom count filter: on [2, 100]"), tmpLogMessage);
        Assertions.assertFalse(tmpLogMessage.contains("run incomplete"), tmpLogMessage);
    }

    /**
     * Tests odd input: empty molecules with and without active structural filters, a molecule for which a filter check
     * throws an exception (removed and counted, nothing propagated), a molecule whose preprocessing failed, a null
     * molecule, and an incomplete run (pending count).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testOddInput() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, true, 2, 100, true));
        IAtomContainer tmpEmpty = SilentChemObjectBuilder.getInstance().newAtomContainer();
        Assertions.assertEquals(MoleculeFilterReason.EMPTY_MOLECULE, tmpRun.applyStructuralFilters(tmpEmpty).getRejectionReason());
        MoleculeFilterRun tmpDefaultRun = new MoleculeFilterRun(MoleculeFilterSettingsSnapshot.createDefault());
        Assertions.assertSame(tmpEmpty, tmpDefaultRun.applyStructuralFilters(tmpEmpty).getMolecule());
        //a molecule for which the element check throws
        IAtomContainer tmpFailingMolecule = new AtomContainer() {
            @Override
            public Iterable<IAtom> atoms() {
                throw new IllegalStateException("Test exception thrown on purpose.");
            }
        };
        tmpFailingMolecule.addAtom(SilentChemObjectBuilder.getInstance().newInstance(IAtom.class, "C"));
        MoleculeFilterResult tmpFailedResult = Assertions.assertDoesNotThrow(() -> tmpRun.applyStructuralFilters(tmpFailingMolecule));
        Assertions.assertEquals(MoleculeFilterReason.FILTER_CHECK_FAILED, tmpFailedResult.getRejectionReason());
        //cyclopentadienyl without hydrogen cannot be kekulized in the preprocessing, the filters still work
        IAtomContainer tmpUnkekulizable = new PreprocessingTestImporter().parseAndPreprocess("c1cccc1").getFirst();
        Assertions.assertTrue(tmpRun.applyStructuralFilters(tmpUnkekulizable).isKept());
        Assertions.assertThrows(NullPointerException.class, () -> tmpRun.applyStructuralFilters(null));
        MoleculeFilterSummary tmpSummary = tmpRun.getSummary();
        Assertions.assertEquals(3, tmpSummary.getReadCount());
        Assertions.assertEquals(1, tmpSummary.getRemovedCount(MoleculeFilterReason.EMPTY_MOLECULE));
        Assertions.assertEquals(1, tmpSummary.getRemovedCount(MoleculeFilterReason.FILTER_CHECK_FAILED));
        //the kept molecule has not been finalised with isDuplicate() yet
        Assertions.assertEquals(1, tmpSummary.getPendingCount());
        Assertions.assertTrue(tmpSummary.toString().contains("run incomplete"));
        tmpRun.registerSmilesGenerationFailure();
        Assertions.assertEquals(0, tmpRun.getSummary().getPendingCount());
        Assertions.assertEquals(1, tmpRun.getSummary().getRemovedCount(MoleculeFilterReason.SMILES_GENERATION_FAILED));
        Assertions.assertThrows(NullPointerException.class, () -> new MoleculeFilterRun((MoleculeFilterSettingsSnapshot) null));
    }

    /**
     * Tests that the structural filtering can be called from several threads at once and the counts still add up.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testConcurrentStructuralFiltering() throws Exception {
        MoleculeFilterRun tmpRun = new MoleculeFilterRun(this.createSnapshot(
                MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT, true, true, 2, 100, false));
        int tmpNumberOfTasks = 4;
        int tmpRepetitions = 10;
        ExecutorService tmpExecutor = Executors.newFixedThreadPool(tmpNumberOfTasks);
        try {
            List<Future<?>> tmpFutures = new ArrayList<>(tmpNumberOfTasks);
            for (int i = 0; i < tmpNumberOfTasks; i++) {
                tmpFutures.add(tmpExecutor.submit(() -> {
                    for (int j = 0; j < tmpRepetitions; j++) {
                        for (IAtomContainer tmpMolecule : this.preprocessDirtySet()) {
                            if (tmpRun.applyStructuralFilters(tmpMolecule).isKept()) {
                                tmpRun.isDuplicate("C");
                            }
                        }
                    }
                    return null;
                }));
            }
            for (Future<?> tmpFuture : tmpFutures) {
                tmpFuture.get();
            }
        } finally {
            tmpExecutor.shutdown();
        }
        MoleculeFilterSummary tmpSummary = tmpRun.getSummary();
        int tmpRuns = tmpNumberOfTasks * tmpRepetitions;
        //per dirty set: 4 disallowed elements, 3 heavy atom count, 1 empty, 6 reduced, 16 kept
        Assertions.assertEquals(24 * tmpRuns, tmpSummary.getReadCount());
        Assertions.assertEquals(16 * tmpRuns, tmpSummary.getKeptCount());
        Assertions.assertEquals(6 * tmpRuns, tmpSummary.getReducedToLargestComponentCount());
        Assertions.assertEquals(8 * tmpRuns, tmpSummary.getTotalRemovedCount());
        Assertions.assertEquals(0, tmpSummary.getPendingCount());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Parses and preprocesses the dirty data set like the importer.
     *
     * @return preprocessed molecules, named "mol" + index
     * @throws Exception if anything goes wrong
     */
    private List<IAtomContainer> preprocessDirtySet() throws Exception {
        return new PreprocessingTestImporter().parseAndPreprocess(MoleculeFilterRunTest.DIRTY_SET);
    }

    /**
     * Runs the given molecules through the filters the way the importer uses them (including the duplicate key) and
     * returns the stored SMILES strings of the kept molecules in input order.
     *
     * @param aRun the filter run
     * @param aMoleculeList preprocessed molecules
     * @param isRegardStereo the importer's stereochemistry setting
     * @param isKekulizationEnforced the importer's Kekule setting
     * @return SMILES strings of the kept molecules
     */
    private List<String> runImport(MoleculeFilterRun aRun, List<IAtomContainer> aMoleculeList, boolean isRegardStereo,
                                   boolean isKekulizationEnforced) {
        List<String> tmpKeptSmiles = new ArrayList<>(aMoleculeList.size());
        for (IAtomContainer tmpMolecule : aMoleculeList) {
            MoleculeFilterResult tmpResult = aRun.applyStructuralFilters(tmpMolecule);
            if (!tmpResult.isKept()) {
                continue;
            }
            String tmpSmiles = ChemUtil.createUniqueSmiles(tmpResult.getMolecule(), isRegardStereo, !isKekulizationEnforced);
            if (tmpSmiles == null || tmpSmiles.isBlank()) {
                aRun.registerSmilesGenerationFailure();
                continue;
            }
            String tmpDuplicateKey = aRun.getSettingsSnapshot().isDuplicateFilterActive()
                    ? MoleculeFilterUtil.createDuplicateKey(tmpResult.getMolecule(), isRegardStereo, tmpSmiles, !isKekulizationEnforced)
                    : tmpSmiles;
            if (aRun.isDuplicate(tmpDuplicateKey)) {
                continue;
            }
            tmpKeptSmiles.add(tmpSmiles);
        }
        return tmpKeptSmiles;
    }

    /**
     * Creates a settings snapshot.
     *
     * @param anOption disconnected structures option
     * @param isElementFilterActive element filter
     * @param isHeavyAtomCountFilterActive heavy atom count filter
     * @param aMinimum minimum heavy atom count
     * @param aMaximum maximum heavy atom count
     * @param isDuplicateFilterActive duplicate filter
     * @return the snapshot
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
     * Returns MORTAR's unique SMILES with stereochemistry and aromaticity encoded (the importer's default).
     *
     * @param aMolecule the molecule
     * @return unique SMILES
     */
    private String uniqueSmiles(IAtomContainer aMolecule) {
        return ChemUtil.createUniqueSmiles(aMolecule, true, true);
    }

    /**
     * Asserts the counts of a summary and that they add up.
     *
     * @param aSummary the summary
     * @param aReadCount expected read count
     * @param aKeptCount expected kept count; ignored for the pending check if no molecule was finalised
     * @param aReducedCount expected reduced count
     * @param anExpectedRemovedCountsMap expected removed counts, missing reasons are expected to be 0
     */
    private void assertSummary(MoleculeFilterSummary aSummary, int aReadCount, int aKeptCount, int aReducedCount,
                               Map<MoleculeFilterReason, Integer> anExpectedRemovedCountsMap) {
        Assertions.assertEquals(aReadCount, aSummary.getReadCount());
        Assertions.assertEquals(aKeptCount, aSummary.getKeptCount());
        Assertions.assertEquals(aReducedCount, aSummary.getReducedToLargestComponentCount());
        this.assertRemovedCounts(aSummary, anExpectedRemovedCountsMap);
        if (aKeptCount > 0) {
            Assertions.assertEquals(0, aSummary.getPendingCount());
            Assertions.assertEquals(aSummary.getReadCount(), aSummary.getKeptCount() + aSummary.getTotalRemovedCount());
        }
    }

    /**
     * Asserts the removed counts per reason of a summary.
     *
     * @param aSummary the summary
     * @param anExpectedRemovedCountsMap expected removed counts, missing reasons are expected to be 0
     */
    private void assertRemovedCounts(MoleculeFilterSummary aSummary, Map<MoleculeFilterReason, Integer> anExpectedRemovedCountsMap) {
        Map<MoleculeFilterReason, Integer> tmpActual = new EnumMap<>(MoleculeFilterReason.class);
        Map<MoleculeFilterReason, Integer> tmpExpected = new EnumMap<>(MoleculeFilterReason.class);
        int tmpExpectedTotal = 0;
        for (MoleculeFilterReason tmpReason : MoleculeFilterReason.values()) {
            tmpActual.put(tmpReason, aSummary.getRemovedCount(tmpReason));
            int tmpExpectedCount = anExpectedRemovedCountsMap.getOrDefault(tmpReason, 0);
            tmpExpected.put(tmpReason, tmpExpectedCount);
            tmpExpectedTotal += tmpExpectedCount;
        }
        Assertions.assertEquals(tmpExpected, tmpActual);
        Assertions.assertEquals(tmpExpectedTotal, aSummary.getTotalRemovedCount());
    }
    //</editor-fold>
}
