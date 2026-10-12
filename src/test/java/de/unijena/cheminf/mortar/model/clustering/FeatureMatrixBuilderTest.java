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

package de.unijena.cheminf.mortar.model.clustering;

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Tests for {@link FeatureMatrixBuilder}.
 *
 * @author Manuel Schauer
 */
class FeatureMatrixBuilderTest {
    /**
     * Name of the test fragmentation.
     */
    private static final String FRAGMENTATION_NAME = "TestFragmentation";
    //
    /**
     * Constructor that sets the default locale to British English, which is needed for the MORTAR message bundle
     * that is used by the enums of IMoleculeFragmenter.
     */
    FeatureMatrixBuilderTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Creates a molecule data model with the given fragment frequencies for the test fragmentation.
     *
     * @param aSmiles SMILES of the molecule
     * @param aName name of the molecule
     * @param aFrequencies fragment frequencies, may be null if the molecule was not fragmented
     * @return molecule data model
     */
    private static MoleculeDataModel createMolecule(String aSmiles, String aName, Map<String, Integer> aFrequencies) {
        MoleculeDataModel tmpMolecule = new MoleculeDataModel(aSmiles, aName, new HashMap<>());
        if (aFrequencies != null) {
            tmpMolecule.getFragmentFrequencies().put(FeatureMatrixBuilderTest.FRAGMENTATION_NAME, aFrequencies);
        }
        return tmpMolecule;
    }
    //
    /**
     * Tests the fragment count and bit fingerprint blocks, the column order and the minimum molecule frequency.
     */
    @Test
    void testBuildFragmentFingerprintMatrix() {
        List<MoleculeDataModel> tmpMolecules = List.of(
                FeatureMatrixBuilderTest.createMolecule("CCO", "m1", Map.of("C", 2, "O", 1)),
                FeatureMatrixBuilderTest.createMolecule("CC", "m2", Map.of("C", 2)),
                FeatureMatrixBuilderTest.createMolecule("N", "m3", null));
        FeatureMatrix tmpCountMatrix = FeatureMatrixBuilder.buildFragmentFingerprintMatrix(
                tmpMolecules, FeatureMatrixBuilderTest.FRAGMENTATION_NAME, true, 1, 0);
        //"C" occurs in two molecules and comes first, "O" in one molecule
        Assertions.assertEquals(2, tmpCountMatrix.getColumnCount());
        Assertions.assertEquals("C", tmpCountMatrix.getFeatureName(0));
        Assertions.assertEquals("O", tmpCountMatrix.getFeatureName(1));
        Assertions.assertEquals(FeatureType.COUNT, tmpCountMatrix.getFeatureType(0));
        Assertions.assertEquals(2.0f, tmpCountMatrix.getValue(0, 0));
        Assertions.assertEquals(1.0f, tmpCountMatrix.getValue(0, 1));
        Assertions.assertEquals(0.0f, tmpCountMatrix.getValue(2, 0));
        FeatureMatrix tmpBitMatrix = FeatureMatrixBuilder.buildFragmentFingerprintMatrix(
                tmpMolecules, FeatureMatrixBuilderTest.FRAGMENTATION_NAME, false, 2, 0);
        Assertions.assertEquals(1, tmpBitMatrix.getColumnCount());
        Assertions.assertEquals(FeatureType.BIT, tmpBitMatrix.getFeatureType(0));
        Assertions.assertEquals(1.0f, tmpBitMatrix.getValue(0, 0));
        Assertions.assertEquals("m2", tmpBitMatrix.getRowName(1));
    }
    //
    /**
     * Tests that the maximum number of fragments keeps the fragments that occur in the most molecules.
     */
    @Test
    void testMaximumNumberOfFragments() {
        List<MoleculeDataModel> tmpMolecules = List.of(
                FeatureMatrixBuilderTest.createMolecule("CCO", "m1", Map.of("C", 2, "O", 1, "N", 1)),
                FeatureMatrixBuilderTest.createMolecule("CC", "m2", Map.of("C", 2, "O", 1)),
                FeatureMatrixBuilderTest.createMolecule("CN", "m3", Map.of("C", 1)));
        FeatureMatrix tmpMatrix = FeatureMatrixBuilder.buildFragmentFingerprintMatrix(
                tmpMolecules, FeatureMatrixBuilderTest.FRAGMENTATION_NAME, true, 1, 2);
        Assertions.assertEquals(2, tmpMatrix.getColumnCount());
        Assertions.assertEquals("C", tmpMatrix.getFeatureName(0));
        Assertions.assertEquals("O", tmpMatrix.getFeatureName(1));
        Assertions.assertThrows(IllegalArgumentException.class, () -> FeatureMatrixBuilder.buildFragmentFingerprintMatrix(
                tmpMolecules, FeatureMatrixBuilderTest.FRAGMENTATION_NAME, true, 1, -1));
    }
    //
    /**
     * Tests that the descriptor calculation applies the given aromaticity model (MORTAR default Daylight with CDK
     * aromatic set) and does not change the atom container of the data model.
     *
     * @throws Exception if the calculation fails
     */
    @Test
    void testBuildDescriptorMatrixAppliesAromaticity() throws Exception {
        MoleculeDataModel tmpBenzene = FeatureMatrixBuilderTest.createMolecule("C1=CC=CC=C1", "benzene", null);
        tmpBenzene.setKeepAtomContainer(true);
        IAtomContainer tmpOriginal = tmpBenzene.getAtomContainer();
        Aromaticity tmpAromaticity = ChemUtil.createAromaticityModel(
                IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT, IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET);
        FeatureMatrix tmpMatrix = FeatureMatrixBuilder.buildDescriptorMatrix(List.of(tmpBenzene),
                new Descriptor[] {Descriptor.AROMATIC_ATOMS_COUNT, Descriptor.ATOM_COUNT_C}, tmpAromaticity, false);
        Assertions.assertEquals(2, tmpMatrix.getColumnCount());
        Assertions.assertEquals(6.0f, tmpMatrix.getValue(0, 0));
        Assertions.assertEquals(6.0f, tmpMatrix.getValue(0, 1));
        Assertions.assertEquals(FeatureType.CONTINUOUS, tmpMatrix.getFeatureType(0));
        Assertions.assertEquals(Descriptor.AROMATIC_ATOMS_COUNT.name(), tmpMatrix.getFeatureGroupName(0));
        //the cached atom container of the data model must not have been aromatised
        for (IAtom tmpAtom : tmpOriginal.atoms()) {
            Assertions.assertFalse(tmpAtom.isAromatic());
        }
    }
    //
    /**
     * Tests the horizontal concatenation of the fragment and descriptor blocks.
     *
     * @throws Exception if the calculation fails
     */
    @Test
    void testConcatenate() throws Exception {
        List<MoleculeDataModel> tmpMolecules = List.of(
                FeatureMatrixBuilderTest.createMolecule("CCO", "m1", Map.of("C", 2, "O", 1)),
                FeatureMatrixBuilderTest.createMolecule("CC", "m2", Map.of("C", 2)));
        FeatureMatrix tmpFragments = FeatureMatrixBuilder.buildFragmentFingerprintMatrix(
                tmpMolecules, FeatureMatrixBuilderTest.FRAGMENTATION_NAME, false, 1, 0);
        FeatureMatrix tmpDescriptors = FeatureMatrixBuilder.buildDescriptorMatrix(tmpMolecules,
                new Descriptor[] {Descriptor.ATOM_COUNT_C},
                ChemUtil.createAromaticityModel(IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT,
                        IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET), true);
        FeatureMatrix tmpCombined = FeatureMatrix.concatenate(tmpFragments, tmpDescriptors);
        Assertions.assertEquals(3, tmpCombined.getColumnCount());
        Assertions.assertEquals(FeatureType.BIT, tmpCombined.getFeatureType(1));
        Assertions.assertEquals(FeatureType.CONTINUOUS, tmpCombined.getFeatureType(2));
        Assertions.assertEquals(2.0f, tmpCombined.getValue(1, 2));
    }
}
