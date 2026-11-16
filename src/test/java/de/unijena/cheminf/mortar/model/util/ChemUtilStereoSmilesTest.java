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

import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.AtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.TetrahedralChirality;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Tests that the stereo SMILES codes created by {@link ChemUtil#createUniqueSmiles(IAtomContainer, boolean)} do not
 * contain configurations of non-stereogenic atoms, so that the same compound always gets the same SMILES code, while
 * real stereocentres, including pseudo-asymmetric ones, are kept. The reference data were created with RDKit by the
 * scripts generate_stereo_smiles_reference.py and generate_biofacquim_stereo_reference.py in the test resources.
 *
 * @author Manuel Schauer
 */
class ChemUtilStereoSmilesTest {
    /**
     * Sets the default locale to British English, required for the settings container of the importer.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Two writings of the same chiral 1,5-anhydropentitol that differ only in the configuration given for the
     * non-stereogenic middle carbon atom C3 must give the same stereo SMILES code, without a configuration for C3
     * (RDKit writes C1[C@@H](C([C@H](CO1)O)O)O for both). The two meso 1,5-anhydropentitols, in which C3 is
     * pseudo-asymmetric, are different compounds and must give different SMILES codes. The given atom containers
     * must not be altered.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testAnhydropentitolMiddleCarbonAtom() throws Exception {
        String tmpChiral1 = this.createStereoSmiles("O[C@H]1COC[C@H](O)[C@H]1O");
        String tmpChiral2 = this.createStereoSmiles("O[C@H]1COC[C@H](O)[C@@H]1O");
        Assertions.assertEquals("C1[C@@H](C([C@H](CO1)O)O)O", tmpChiral1);
        Assertions.assertEquals(tmpChiral1, tmpChiral2);
        String tmpMeso1 = this.createStereoSmiles("O[C@@H]1COC[C@H](O)[C@H]1O");
        String tmpMeso2 = this.createStereoSmiles("O[C@@H]1COC[C@H](O)[C@@H]1O");
        Assertions.assertNotEquals(tmpMeso1, tmpMeso2);
        Assertions.assertNotEquals(tmpChiral1, tmpMeso1);
        Assertions.assertNotEquals(tmpChiral1, tmpMeso2);
        //the two writings in which the CRUSH fragmenter produced this fragment from different parent molecules
        Assertions.assertEquals(this.createStereoSmiles("C1[C@@H]([C@H]([C@H](CO1)O)O)O"),
                this.createStereoSmiles("C1[C@@H]([C@@H]([C@H](CO1)O)O)O"));
        //the given atom container is not altered
        IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer("O[C@H]1COC[C@H](O)[C@H]1O", false, false);
        ChemUtil.createUniqueSmiles(tmpMolecule, true);
        Assertions.assertEquals(3, this.countStereoElements(tmpMolecule));
    }
    //
    /**
     * A molecule whose only stereo element is removed as non-stereogenic must get the same SMILES code as the molecule
     * without any stereo element, i.e. not one in the canonical atom order of stereo SMILES codes.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testOnlyNonStereogenicElement() throws Exception {
        Assertions.assertEquals(this.createStereoSmiles("*C(*)O"), this.createStereoSmiles("*[C@H](*)O"));
        Assertions.assertEquals(this.createStereoSmiles("*C(*)O"), this.createStereoSmiles("*[C@@H](*)O"));
        Assertions.assertEquals(this.createStereoSmiles("CC(C)=CC"), this.createStereoSmiles("C/C(C)=C/C"));
    }
    //
    /**
     * Tetrahedral stereo elements with repeated carriers, as they occur in some alkyl fragments (e.g. on the methyl
     * carbon atoms of an isobutane fragment, with the central atom given three times), must not lead to a stereo
     * SMILES code for isobutane.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testDegenerateTetrahedralElements() throws Exception {
        IAtomContainer tmpIsobutane = ChemUtil.parseSmilesToAtomContainer("CC(C)C", false, false);
        IAtom tmpCentre = tmpIsobutane.getAtom(1);
        tmpIsobutane.addStereoElement(new TetrahedralChirality(tmpCentre,
                new IAtom[] {tmpIsobutane.getAtom(0), tmpIsobutane.getAtom(2), tmpIsobutane.getAtom(3), tmpCentre},
                ITetrahedralChirality.Stereo.CLOCKWISE));
        for (int i : new int[] {0, 2, 3}) {
            IAtom tmpMethyl = tmpIsobutane.getAtom(i);
            tmpIsobutane.addStereoElement(new TetrahedralChirality(tmpMethyl,
                    new IAtom[] {tmpCentre, tmpMethyl, tmpMethyl, tmpMethyl}, ITetrahedralChirality.Stereo.CLOCKWISE));
        }
        Assertions.assertEquals("CC(C)C", ChemUtil.createUniqueSmiles(tmpIsobutane, true));
        Assertions.assertEquals(4, this.countStereoElements(tmpIsobutane));
    }
    //
    /**
     * Hand-picked reference compounds (anhydropentitols, pentitols, disubstituted cyclohexanes, inositols, tartaric
     * and trihydroxyglutaric acids, the same with link atoms, normal stereocentres and double bonds), each in several
     * writings: all configurations of the marked atoms and bonds, and random atom orders. The stereo SMILES codes must
     * be equal for all writings of one compound (same RDKit canonical isomeric SMILES) and different for different
     * compounds.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testStereoSmilesReference() throws Exception {
        Map<String, Set<String>> tmpSmilesByCompound = new LinkedHashMap<>();
        Map<String, Set<String>> tmpCompoundsBySmiles = new HashMap<>();
        int tmpWritingCount = 0;
        for (String[] tmpRow : this.readTsv("StereoSmilesReference.tsv")) {
            //columns: name, input SMILES, RDKit canonical SMILES, InChIKey
            IAtomContainer tmpMolecule = ChemUtil.parseSmilesToAtomContainer(tmpRow[1], false, false);
            String tmpSmiles = ChemUtil.createUniqueSmiles(tmpMolecule, true);
            Assertions.assertNotNull(tmpSmiles, tmpRow[1]);
            tmpSmilesByCompound.computeIfAbsent(tmpRow[2], aKey -> new TreeSet<>()).add(tmpSmiles);
            tmpCompoundsBySmiles.computeIfAbsent(tmpSmiles, aKey -> new TreeSet<>()).add(tmpRow[2]);
            tmpWritingCount++;
        }
        Assertions.assertEquals(725, tmpWritingCount);
        Assertions.assertEquals(149, tmpSmilesByCompound.size());
        List<String> tmpSplitCompounds = new ArrayList<>();
        tmpSmilesByCompound.forEach((aCompound, aSmilesSet) -> {
            if (aSmilesSet.size() > 1) {
                tmpSplitCompounds.add(aCompound + " -> " + aSmilesSet);
            }
        });
        List<String> tmpMergedCompounds = new ArrayList<>();
        tmpCompoundsBySmiles.forEach((aSmiles, aCompoundSet) -> {
            if (aCompoundSet.size() > 1) {
                tmpMergedCompounds.add(aSmiles + " <- " + aCompoundSet);
            }
        });
        Assertions.assertTrue(tmpSplitCompounds.isEmpty(), "compounds with more than one SMILES: " + tmpSplitCompounds);
        Assertions.assertTrue(tmpMergedCompounds.isEmpty(), "SMILES of more than one compound: " + tmpMergedCompounds);
    }
    //
    /**
     * Data set check on BIOFACQUIM (Tutorial/biofacquim-06-2025.sdf): for every molecule, the stereo SMILES code of
     * the molecule as imported from the SD file must be the same as the ones of three RDKit writings of it with random
     * atom orders, and molecules that RDKit regards as different must not get the same SMILES code. The SD file is
     * imported with the preprocessing of the MORTAR importer.
     * <br>Excluded difference unrelated to stereochemistry: MORTAR kekulizes imported molecules without perceiving
     * aromaticity, and the unique SMILES code of a Kekule structure depends on its Kekule form, which differs between
     * the SD file and the RDKit writings (e.g. for record 6). Aromaticity is therefore perceived on both sides (Daylight
     * model) before the SMILES codes are created. Writings that RDKit itself does not parse back to the same compound
     * (some writings of record 456) were already discarded by the generator script.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testBiofacquimWritings() throws Exception {
        IAtomContainerSet tmpSet = new AtomContainerSet();
        try (IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(
                Paths.get("Tutorial", "biofacquim-06-2025.sdf").toFile()), SilentChemObjectBuilder.getInstance())) {
            while (tmpReader.hasNext()) {
                tmpSet.addAtomContainer(tmpReader.next());
            }
        }
        PreprocessingImporter tmpImporter = new PreprocessingImporter();
        tmpImporter.preprocess(tmpSet);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));
        List<String[]> tmpRows = this.readTsv("BIOFACQUIMStereoSmilesReference.tsv");
        Assertions.assertEquals(583, tmpSet.getAtomContainerCount());
        Assertions.assertEquals(583, tmpRows.size());
        List<String> tmpInconsistentRecords = new ArrayList<>();
        Map<String, String> tmpCompoundBySmiles = new HashMap<>();
        List<String> tmpMergedRecords = new ArrayList<>();
        for (int i = 0; i < tmpRows.size(); i++) {
            //columns: record, RDKit canonical SMILES, writings separated by spaces
            String[] tmpRow = tmpRows.get(i);
            Assertions.assertEquals(String.valueOf(i + 1), tmpRow[0]);
            IAtomContainer tmpMolecule = tmpSet.getAtomContainer(i);
            tmpAromaticity.apply(tmpMolecule);
            String tmpSmiles = ChemUtil.createUniqueSmiles(tmpMolecule, true);
            Assertions.assertNotNull(tmpSmiles, tmpRow[0]);
            for (String tmpWriting : tmpRow[2].split(" ")) {
                IAtomContainerSet tmpWritingSet = new AtomContainerSet();
                tmpWritingSet.addAtomContainer(ChemUtil.parseSmilesToAtomContainer(tmpWriting));
                tmpImporter.preprocess(tmpWritingSet);
                IAtomContainer tmpWritingMolecule = tmpWritingSet.getAtomContainer(0);
                tmpAromaticity.apply(tmpWritingMolecule);
                String tmpWritingSmiles = ChemUtil.createUniqueSmiles(tmpWritingMolecule, true);
                if (!tmpSmiles.equals(tmpWritingSmiles)) {
                    tmpInconsistentRecords.add(tmpRow[0] + ": " + tmpSmiles + " vs. " + tmpWritingSmiles);
                }
            }
            String tmpPreviousCompound = tmpCompoundBySmiles.put(tmpSmiles, tmpRow[1]);
            if (tmpPreviousCompound != null && !tmpPreviousCompound.equals(tmpRow[1])) {
                tmpMergedRecords.add(tmpRow[0] + ": " + tmpSmiles);
            }
        }
        Assertions.assertTrue(tmpInconsistentRecords.isEmpty(), "inconsistent records: " + tmpInconsistentRecords);
        Assertions.assertTrue(tmpMergedRecords.isEmpty(), "merged records: " + tmpMergedRecords);
    }
    //
    /**
     * Parses the given SMILES code without kekulization and atom typing and creates its stereo SMILES code.
     *
     * @param aSmiles SMILES code
     * @return unique stereo SMILES code
     * @throws Exception if anything goes wrong
     */
    private String createStereoSmiles(String aSmiles) throws Exception {
        return ChemUtil.createUniqueSmiles(ChemUtil.parseSmilesToAtomContainer(aSmiles, false, false), true);
    }
    //
    /**
     * Counts the stereo elements of the given atom container.
     *
     * @param anAtomContainer atom container
     * @return number of stereo elements
     */
    @SuppressWarnings("rawtypes")
    private int countStereoElements(IAtomContainer anAtomContainer) {
        int tmpCount = 0;
        for (IStereoElement tmpElement : anAtomContainer.stereoElements()) {
            tmpCount++;
        }
        return tmpCount;
    }
    //
    /**
     * Reads the given tab-separated test resource without its comment lines (starting with '#') and its header line.
     *
     * @param aResourceName name of the resource next to this class
     * @return rows split at tabs
     * @throws Exception if anything goes wrong
     */
    private List<String[]> readTsv(String aResourceName) throws Exception {
        List<String[]> tmpRows = new ArrayList<>();
        try (BufferedReader tmpReader = new BufferedReader(new InputStreamReader(
                Objects.requireNonNull(this.getClass().getResourceAsStream(aResourceName)), StandardCharsets.UTF_8))) {
            String tmpLine;
            boolean tmpIsHeaderSkipped = false;
            while ((tmpLine = tmpReader.readLine()) != null) {
                if (tmpLine.startsWith("#")) {
                    continue;
                }
                if (!tmpIsHeaderSkipped) {
                    tmpIsHeaderSkipped = true;
                    continue;
                }
                tmpRows.add(tmpLine.split("\t", -1));
            }
        }
        return tmpRows;
    }
    //
    /**
     * Gives access to the preprocessing of the MORTAR importer (atom typing, implicit hydrogen atoms, kekulization).
     */
    private static class PreprocessingImporter extends Importer {
        /**
         * Constructor.
         */
        PreprocessingImporter() {
            super(new SettingsContainer());
        }

        /**
         * Preprocesses the given molecules like the importer does, with open valences filled with implicit hydrogen
         * atoms (default setting).
         *
         * @param aMoleculeSet molecules
         */
        void preprocess(IAtomContainerSet aMoleculeSet) {
            this.preprocessMoleculeSet(aMoleculeSet, true);
        }
    }
}
