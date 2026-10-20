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

import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.stereo.TetrahedralChirality;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Integration tests for {@link AlkylStructureFragmenter}.
 * <p>
 *     Instead of comparing against hand-written expected fragments, these tests run the complete fragmentation
 *     pipeline (filter, fragmentation, saturation, stereo remap) on a broad set of hydrocarbons under every
 *     combination of settings and check structural invariants that must hold for every result:
 *     every atom of the input ends up in exactly one fragment, every fragment bond is an input bond with unchanged
 *     order and is used only once, fragments are connected and contain no dangling bonds, carbon valences are not
 *     exceeded, linear chains respect the maximum chain length, and the result does not depend on the order in which
 *     atoms and bonds are stored in the input container.
 * </p>
 *
 * @author Maximilian Rottmann, Manuel Schauer
 * @version 1.0.0.0
 */
public class AlkylStructureFragmenterIntegrationTest {
    /**
     * Static setter for Locale, see {@link AlkylStructureFragmenterTest}.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    /**
     * Hydrocarbon test structures covering chains, branching, rings, ring systems, (conjugated) multi bonds,
     * allenes, aromatics and combinations of these, including several natural products.
     */
    private static final String[] TEST_SMILES = {
            //linear and branched alkanes
            "C", "CC", "CCCCCCCCCCCCCCCC", "CC(C)C", "CC(C)(C)C", "CC(C)(C)C(C)(C)C", "CC(C)C(C)C",
            "CCC(C)CC(C)(C)CCC", "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C", "CCCCCCCCCC(CC)(CCCCCCC)CCCCCCCC",
            "CC(C)(C)CCCC(C)C", "CC(C)CC(C)(C)CC(C)(C)CC(C)C",
            //multi bonds and conjugated systems
            "C=C", "C#C", "CC=CC", "CC#CC", "C=CC=C", "C#CC#C", "C=C=C", "C=C=C=C", "CC=C=CC", "C=CC(C)=C",
            "CC(C)=CCCC(=C)C=C", "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C", "C=CC=CC=C", "C#CC#CC",
            "CCC#CCC", "CC(C)C#CC(C)C", "C=CC=CC(C)(C)C", "C=CCC=C", "C=CCCC#C", "CC(=C)C(C)C", "CC(C)(C)C=C",
            "CC(C)(C)C#C", "CC(C)(C=C)C#C", "C=CC(C)(C)C=C",
            //rings and ring systems
            "C1CC1", "CC1CC1", "C1CCCCC1", "CC1CCCCC1", "CC1(C)CCCCC1(C)C", "C1CCC2(CC1)CCCC2", "C1CCC2CCCCC2C1",
            "C1C2CC3CC1CC(C2)C3", "CC1CCCCCCCCC1", "C1CCC(CC1)C1CCCCC1", "C1CCC(CC1)CC1CCCCC1",
            "CC(C)CC(C)C1CCCCC1", "CC(C)(C)C1CCCCC1", "CC(C)(C)CC1CCCCC1", "CC(C)(C1CCCCC1)C1CCCCC1",
            "C=C1CCCCC1", "CC(C)=C1CCCCC1", "C=CC1CCC=CC1", "C=C1CCCC2CC12", "C1CC=CC1", "C1=CCC=C1",
            "CC1=CCC(CC1)C(=C)C", "CC1=CCC2CC1C2(C)C", "CC1(C)C2CCC(=C)C1C2", "C[C@@H]1C[C@@H](C)[C@@H](C)C1",
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C",
            "CC1CCC(C(C)C)C(C)C1", "C1CCC(=C2CCCCC2)CC1",
            //aromatics
            "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "CC(C)c1ccccc1", "CC(C)(C)c1ccccc1", "CC(C)(C)Cc1ccccc1",
            "Cc1ccc(cc1)C(C)C", "c1ccc(cc1)-c1ccccc1", "c1ccc2ccccc2c1", "c1ccc2cccc2cc1", "C=Cc1ccccc1",
            "C(=C/c1ccccc1)\\c1ccccc1", "C(#Cc1ccccc1)c1ccccc1", "C(c1ccccc1)(c1ccccc1)c1ccccc1",
            "CC(C)/C=C/C=C/c1ccccc1", "C=CC1=CC=C2C=CC=CC2=C1", "C=C1C=CC=C1", "C1=CC=CC=CC1",
            "C12=CC(=CC=C1C(=CC(=C2)C3CCCCC3)CC(C)(C)C/C=C/CCC)/C=C\\C(C)C",
            "CC(CC1C2CC2C(=C)CC1c1ccc(cc1C\\C=C/c1ccccc1)C(C)C)(CCCCCCCCC)CC#CC",
            "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C",
            "CC(C)C1=CC=C(C)CC1", "CC(C)c1ccc(C)cc1C(C)(C)C", "C(Cc1ccccc1)Cc1ccccc1",
            //disconnected input
            "CCCC.CCCCCCCC", "CC(C)C.C1CCCCC1C1CCCCC1"
    };
    /**
     * Natural product library (PubChem natural products, 02/2026) in the test resources, one isomeric SMILES and the
     * COCONUT identifier per line, tab separated. The SMILES were generated with the CDK from the original SDF file
     * (stereo chemistry perceived from its 2D wedge bonds).
     */
    private static final String NATURAL_PRODUCT_LIBRARY = "pubchem-nps-02-2026.smi";
    /**
     * Aromaticity model used to normalize fragment SMILES.
     */
    private static final Aromaticity AROMATICITY = new Aromaticity(ElectronDonation.daylight(), Cycles.all());
    /**
     * Statistics of the stereo chemistry checks (see checkStereoChemistry()), collected over all tests of this class.
     */
    private static final Map<String, Integer> STEREO_STATISTICS = new java.util.concurrent.ConcurrentHashMap<>();
    /**
     * Maximum chain length setting values to test.
     */
    private static final int[] MAX_CHAIN_LENGTHS = {1, 2, 3, 6};
    /**
     * Number of random atom/bond order permutations per molecule in the order invariance test.
     */
    private static final int PERMUTATION_COUNT = 5;

    /**
     * Runs all molecules of the test SMILES list and the substructure library through all combinations of settings
     * and checks the structural invariants of each fragmentation result.
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void fragmentationInvariantsTest() throws Exception {
        List<IAtomContainer> tmpMolecules = this.loadTestMolecules();
        List<String> tmpViolations = new ArrayList<>();
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (IAtomContainer tmpMolecule : tmpMolecules) {
                tmpViolations.addAll(this.checkInvariants(tmpASF, tmpMolecule));
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Checks that the fragmentation result does not depend on the order of atoms and bonds in the input container
     * and on the begin/end orientation of the bonds.
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void atomOrderInvarianceTest() throws Exception {
        List<IAtomContainer> tmpMolecules = this.loadTestMolecules();
        List<String> tmpViolations = new ArrayList<>();
        Random tmpRandom = new Random(42);
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (IAtomContainer tmpMolecule : tmpMolecules) {
                IAtomContainer tmpReference = this.copyWithoutStereo(tmpMolecule, null);
                List<String> tmpExpected = this.fragmentToSortedSmilesOrException(tmpASF, tmpReference);
                for (int i = 0; i < AlkylStructureFragmenterIntegrationTest.PERMUTATION_COUNT; i++) {
                    IAtomContainer tmpPermuted = this.copyWithoutStereo(tmpMolecule, tmpRandom);
                    List<String> tmpActual = this.fragmentToSortedSmilesOrException(tmpASF, tmpPermuted);
                    if (!tmpExpected.equals(tmpActual)) {
                        tmpViolations.add(String.format("%s [%s]: %s vs. permuted %s",
                                ChemUtil.createUniqueSmiles(tmpMolecule, false), this.describeSettings(tmpASF),
                                tmpExpected, tmpActual));
                        break;
                    }
                }
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Dedicated test of the stereo chemistry transfer (see checkStereoChemistry()) on stereo test structures, the
     * natural product skeletons (stereo chemistry of the library's isomeric SMILES) and
     * all setting combinations, with stereo chemistry preservation switched on and off. Fails if no stereo element was
     * checked at all, so that the test cannot pass vacuously.
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void stereoChemistryTest() throws Exception {
        String[] tmpStereoSmiles = {
                "C[C@@H]1C[C@@H](C)[C@@H](C)C1", "C[C@H](CC)CCC", "CC[C@@H](C)[C@H](C)CC", "C/C=C/C", "C/C=C\\C",
                "C/C=C/C=C/C", "CC/C=C(/C)CC", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C",
                "CCCC[C@H](C)[C@H]1CC[C@]2(C)[C@H](CC[C@@H]2[C@H]1C)[C@H](C)CCCC(C)C", "C[C@@H](c1ccccc1)CC",
                "C(=C/c1ccccc1)\\c1ccccc1", "C[C@H]1CC[C@@H](C(C)C)CC1", "C[C@@H]1CC[C@H](C)C1/C=C/C"};
        List<IAtomContainer> tmpMolecules = new ArrayList<>();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        for (String tmpSmiles : tmpStereoSmiles) {
            tmpMolecules.add(tmpParser.parseSmiles(tmpSmiles));
        }
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            IAtomContainer tmpSkeleton = this.createCarbonSkeleton(tmpMolecule);
            if (!tmpSkeleton.isEmpty() && tmpSkeleton.stereoElements().iterator().hasNext()) {
                tmpMolecules.add(tmpSkeleton);
            }
        }
        STEREO_STATISTICS.clear();
        List<String> tmpViolations = new ArrayList<>();
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (boolean tmpPreserveStereo : new boolean[]{true, false}) {
                tmpASF.setPreserveStereoChemistrySetting(tmpPreserveStereo);
                for (IAtomContainer tmpMolecule : tmpMolecules) {
                    tmpViolations.addAll(this.checkInvariants(tmpASF, tmpMolecule));
                }
            }
        }
        STEREO_STATISTICS.forEach((tmpKey, tmpValue) -> System.out.println("STATISTICS " + tmpKey + ": " + tmpValue));
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(STEREO_STATISTICS.getOrDefault("tetrahedral centers expected in fragments", 0) > 0);
        Assertions.assertTrue(STEREO_STATISTICS.getOrDefault("double bond stereo expected in fragments", 0) > 0);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Checks that explicit hydrogen atoms in the input do not change the fragmentation result.
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void explicitHydrogenInvarianceTest() throws Exception {
        List<String> tmpViolations = new ArrayList<>();
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (IAtomContainer tmpMolecule : this.loadTestMolecules()) {
                IAtomContainer tmpImplicit = this.copyWithoutStereo(tmpMolecule, null);
                IAtomContainer tmpExplicit = this.copyWithoutStereo(tmpMolecule, null);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpExplicit);
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpExplicit);
                List<String> tmpExpected = this.fragmentToSortedSmilesOrException(tmpASF, tmpImplicit);
                List<String> tmpActual = this.fragmentToSortedSmilesOrException(tmpASF, tmpExplicit);
                if (!tmpExpected.equals(tmpActual)) {
                    tmpViolations.add(String.format("%s [%s]: %s vs. explicit hydrogens %s",
                            ChemUtil.createUniqueSmiles(tmpMolecule, false), this.describeSettings(tmpASF),
                            tmpExpected, tmpActual));
                }
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Checks that fragmenting the same molecule twice with the same instance gives identical results and that the
     * input molecule is not altered by the fragmentation.
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void repeatabilityTest() throws Exception {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (IAtomContainer tmpMolecule : this.loadTestMolecules()) {
            String tmpInputSmiles = ChemUtil.createUniqueSmiles(tmpMolecule, true);
            List<String> tmpFirst = this.fragmentToSortedSmiles(tmpASF, tmpMolecule);
            List<String> tmpSecond = this.fragmentToSortedSmiles(tmpASF, tmpMolecule);
            Assertions.assertEquals(tmpFirst, tmpSecond, tmpInputSmiles);
            Assertions.assertEquals(tmpInputSmiles, ChemUtil.createUniqueSmiles(tmpMolecule, true));
        }
    }

    /**
     * Runs all molecules of the natural product library (PubChem natural products) unchanged through the pipeline.
     * As most natural products contain hetero atoms, this mainly tests the filter: pure hydrocarbons must pass and be
     * fragmented without violating invariants, all other molecules must be filtered, or returned unchanged as a
     * single fragment if non-fragmentable molecules are kept in the pipeline.
     *
     * @throws Exception if the library cannot be read
     */
    @Test
    public void naturalProductLibraryTest() throws Exception {
        List<String> tmpViolations = new ArrayList<>();
        //non-fragmentable molecules are only filtered if they are not kept in the pipeline
        AlkylStructureFragmenter tmpDefaultASF = new AlkylStructureFragmenter();
        tmpDefaultASF.setKeepNonFragmentableMoleculesSetting(false);
        AlkylStructureFragmenter tmpKeepingASF = new AlkylStructureFragmenter();
        Assertions.assertTrue(tmpKeepingASF.getKeepNonFragmentableMoleculesSetting(), "Non-fragmentable molecules are kept by default.");
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            String tmpContext = tmpMolecule.getTitle() + " " + ChemUtil.createUniqueSmiles(tmpMolecule, false);
            boolean tmpIsHydrocarbon = this.isHydrocarbon(tmpMolecule);
            try {
                if (tmpDefaultASF.shouldBeFiltered(tmpMolecule) == tmpIsHydrocarbon) {
                    tmpViolations.add(tmpContext + ": filter result does not match hydrocarbon check");
                }
                if (tmpIsHydrocarbon) {
                    tmpViolations.addAll(this.checkInvariants(tmpDefaultASF, tmpMolecule));
                    continue;
                }
                if (tmpKeepingASF.shouldBeFiltered(tmpMolecule)) {
                    tmpViolations.add(tmpContext + ": filtered although non-fragmentable molecules are kept");
                }
                List<IAtomContainer> tmpFragments = tmpKeepingASF.fragmentMolecule(tmpMolecule);
                if (tmpFragments.size() != 1 || tmpFragments.getFirst() != tmpMolecule) {
                    tmpViolations.add(tmpContext + ": kept non-fragmentable molecule was altered");
                }
            } catch (Exception anException) {
                tmpViolations.add(tmpContext + ": exception " + anException);
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Derives a hydrocarbon skeleton from every molecule of the natural product library (see
     * {@link #createCarbonSkeleton(IAtomContainer)}) and checks the fragmentation invariants and the atom order
     * invariance for all setting combinations. This gives a large set of realistic carbon frameworks (fused and
     * bridged ring systems, highly substituted centers, stereo centers) that hand-written test structures do not cover.
     *
     * @throws Exception if the library cannot be read
     */
    @Test
    public void naturalProductSkeletonTest() throws Exception {
        List<IAtomContainer> tmpSkeletons = new ArrayList<>();
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            IAtomContainer tmpSkeleton = this.createCarbonSkeleton(tmpMolecule);
            if (!tmpSkeleton.isEmpty()) {
                tmpSkeletons.add(tmpSkeleton);
            }
        }
        List<String> tmpViolations = new ArrayList<>();
        Random tmpRandom = new Random(42);
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (IAtomContainer tmpSkeleton : tmpSkeletons) {
                tmpViolations.addAll(this.checkInvariants(tmpASF, tmpSkeleton));
                List<String> tmpExpected = this.fragmentToSortedSmilesOrException(tmpASF, this.copyWithoutStereo(tmpSkeleton, null));
                List<String> tmpActual = this.fragmentToSortedSmilesOrException(tmpASF, this.copyWithoutStereo(tmpSkeleton, tmpRandom));
                if (!tmpExpected.equals(tmpActual)) {
                    tmpViolations.add(String.format("%s %s [%s]: %s vs. permuted %s", tmpSkeleton.getTitle(),
                            ChemUtil.createUniqueSmiles(tmpSkeleton, false), this.describeSettings(tmpASF),
                            tmpExpected, tmpActual));
                }
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Simulates the intended pipeline use: the natural products are fragmented with the
     * {@link ErtlFunctionalGroupsFinderFragmenter} (functional groups and alkane fragments) and the resulting fragments
     * are passed on to the alkyl fragmenter. As in the MORTAR pipeline, every fragment is stored as unique SMILES and
     * parsed again (with a fallback without kekulization, which gives aromatic bonds without bond order for remnants of
     * aromatic rings). Every alkane fragment must be accepted and fragmented without violating any invariant or
     * chemistry rule under all setting combinations; functional groups must be returned unchanged if non-fragmentable
     * molecules are kept in the pipeline.
     *
     * @throws Exception if the library cannot be read
     */
    @Test
    public void ertlPipelineTest() throws Exception {
        ErtlFunctionalGroupsFinderFragmenter tmpErtlFragmenter = new ErtlFunctionalGroupsFinderFragmenter();
        tmpErtlFragmenter.setReturnedFragmentsSetting(ErtlFunctionalGroupsFinderFragmenter.EFGFFragmenterReturnedFragmentsOption.ALL_FRAGMENTS);
        Map<String, IAtomContainer> tmpAlkaneFragmentsMap = new HashMap<>();
        Map<String, IAtomContainer> tmpFunctionalGroupsMap = new HashMap<>();
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            IAtomContainer tmpParsedMolecule = this.roundTripAsInPipeline(tmpMolecule);
            if (tmpErtlFragmenter.shouldBeFiltered(tmpParsedMolecule)) {
                continue;
            }
            for (IAtomContainer tmpFragment : tmpErtlFragmenter.fragmentMolecule(tmpParsedMolecule)) {
                String tmpSmiles = ChemUtil.createUniqueSmiles(tmpFragment, false);
                if (tmpSmiles == null) {
                    continue;
                }
                IAtomContainer tmpParsedFragment = this.roundTripAsInPipeline(tmpFragment);
                tmpParsedFragment.setTitle(tmpMolecule.getTitle() + " Ertl fragment");
                if (ErtlFunctionalGroupsFinderFragmenter.FRAGMENT_CATEGORY_ALKANE_VALUE.equals(
                        tmpFragment.getProperty(IMoleculeFragmenter.FRAGMENT_CATEGORY_PROPERTY_KEY))) {
                    tmpAlkaneFragmentsMap.putIfAbsent(tmpSmiles, tmpParsedFragment);
                } else {
                    tmpFunctionalGroupsMap.putIfAbsent(tmpSmiles, tmpParsedFragment);
                }
            }
        }
        Assertions.assertFalse(tmpAlkaneFragmentsMap.isEmpty());
        List<String> tmpViolations = new ArrayList<>();
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            for (Map.Entry<String, IAtomContainer> tmpEntry : tmpAlkaneFragmentsMap.entrySet()) {
                if (tmpASF.shouldBeFiltered(tmpEntry.getValue())) {
                    tmpViolations.add(tmpEntry.getKey() + ": Ertl alkane fragment filtered by alkyl fragmenter");
                    continue;
                }
                tmpViolations.addAll(this.checkInvariants(tmpASF, tmpEntry.getValue()));
            }
        }
        AlkylStructureFragmenter tmpKeepingASF = new AlkylStructureFragmenter();
        tmpKeepingASF.setKeepNonFragmentableMoleculesSetting(true);
        for (Map.Entry<String, IAtomContainer> tmpEntry : tmpFunctionalGroupsMap.entrySet()) {
            if (this.isHydrocarbon(tmpEntry.getValue())) {
                //e.g. C=C or C#C are functional groups in Ertl's definition and fragmented as usual
                continue;
            }
            try {
                if (tmpKeepingASF.shouldBeFiltered(tmpEntry.getValue())) {
                    tmpViolations.add(tmpEntry.getKey() + ": functional group filtered although non-fragmentable molecules are kept");
                    continue;
                }
                List<IAtomContainer> tmpResult = tmpKeepingASF.fragmentMolecule(tmpEntry.getValue());
                if (tmpResult.size() != 1 || tmpResult.getFirst() != tmpEntry.getValue()) {
                    tmpViolations.add(tmpEntry.getKey() + ": kept functional group was altered");
                }
            } catch (Exception anException) {
                tmpViolations.add(tmpEntry.getKey() + ": exception " + anException);
            }
        }
        tmpViolations.forEach(System.out::println);
        Assertions.assertTrue(tmpViolations.isEmpty(), () -> tmpViolations.size() + " violations:\n"
                + String.join("\n", tmpViolations.subList(0, Math.min(tmpViolations.size(), 200))));
    }

    /**
     * Differential test of the rule-based {@link AlkylStructureFragmenter} against the previous implementation
     * ({@link LegacyAlkylStructureFragmenter}, test reference only): both are applied with all setting combinations to
     * the hand-picked test molecules, the substructure library, the natural product skeletons and the Ertl alkane
     * fragments of the natural products, and their fragments (unique stereo SMILES) are compared.
     * <p>
     * Two differences were decided to be intended (see AlkylStructureFragmenterTest.specificTest10 and specificTest11);
     * they are recognized structurally and counted, every other difference fails the test:
     * </p>
     * <ul>
     *     <li>A: isolation and separation from rings active, molecule with a tertiary/quaternary carbon bonded to a ring:
     *     the neighbors of this carbon are ordinary chain atoms (legacy cut off their first chain atom).</li>
     *     <li>B: isolation active, separation from rings inactive, molecule with a tertiary/quaternary carbon bonded to a
     *     ring and to a conjugated non-ring atom: the bond to the conjugated system is kept (legacy cut it).</li>
     * </ul>
     *
     * @throws Exception if the test structures cannot be read
     */
    @Test
    public void differentialTest() throws Exception {
        List<IAtomContainer> tmpMolecules = new ArrayList<>(this.loadTestMolecules());
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            IAtomContainer tmpSkeleton = this.createCarbonSkeleton(tmpMolecule);
            if (!tmpSkeleton.isEmpty()) {
                tmpMolecules.add(tmpSkeleton);
            }
        }
        ErtlFunctionalGroupsFinderFragmenter tmpErtlFragmenter = new ErtlFunctionalGroupsFinderFragmenter();
        tmpErtlFragmenter.setReturnedFragmentsSetting(ErtlFunctionalGroupsFinderFragmenter.EFGFFragmenterReturnedFragmentsOption.ONLY_ALKANE_FRAGMENTS);
        Set<String> tmpSeenErtlFragments = new HashSet<>();
        for (IAtomContainer tmpMolecule : this.loadNaturalProductLibrary()) {
            IAtomContainer tmpParsedMolecule = this.roundTripAsInPipeline(tmpMolecule);
            if (tmpErtlFragmenter.shouldBeFiltered(tmpParsedMolecule)) {
                continue;
            }
            for (IAtomContainer tmpFragment : tmpErtlFragmenter.fragmentMolecule(tmpParsedMolecule)) {
                String tmpSmiles = ChemUtil.createUniqueSmiles(tmpFragment, false);
                if (tmpSmiles != null && tmpSeenErtlFragments.add(tmpSmiles)) {
                    tmpMolecules.add(this.roundTripAsInPipeline(tmpFragment));
                }
            }
        }
        List<String> tmpDifferences = new ArrayList<>();
        Map<String, Integer> tmpIntendedDifferenceCounts = new HashMap<>();
        for (AlkylStructureFragmenter tmpASF : this.createAllSettingCombinations()) {
            LegacyAlkylStructureFragmenter tmpLegacyASF = new LegacyAlkylStructureFragmenter();
            tmpLegacyASF.setFragmentSideChainsSetting(tmpASF.getFragmentSideChainsSetting());
            tmpLegacyASF.setMaxChainLengthSetting(tmpASF.getMaxChainLengthSetting());
            tmpLegacyASF.setIsolateTertQuatCarbonsSetting(tmpASF.getIsolateTertQuatCarbonSetting());
            tmpLegacyASF.setSeparateTertQuatCarbonFromRingSetting(tmpASF.getSeparateTertQuatCarbonFromRingSetting());
            for (IAtomContainer tmpMolecule : tmpMolecules) {
                List<String> tmpNew = this.fragmentToSortedStereoSmiles(tmpASF, tmpMolecule);
                List<String> tmpLegacy = this.fragmentToSortedStereoSmiles(tmpLegacyASF, tmpMolecule);
                if (!tmpNew.equals(tmpLegacy)) {
                    String tmpIntendedDifference = this.getIntendedDifference(tmpASF, tmpMolecule);
                    if (tmpIntendedDifference != null) {
                        tmpIntendedDifferenceCounts.merge(tmpIntendedDifference, 1, Integer::sum);
                        continue;
                    }
                    tmpDifferences.add(String.format("%s [%s]: legacy %s vs. new %s", ChemUtil.createUniqueSmiles(tmpMolecule, true),
                            this.describeSettings(tmpASF), tmpLegacy, tmpNew));
                }
            }
        }
        tmpIntendedDifferenceCounts.forEach((tmpKey, tmpValue) -> System.out.println("INTENDED " + tmpKey + ": " + tmpValue));
        tmpDifferences.forEach(tmpDifference -> System.out.println("DIFF " + tmpDifference));
        Assertions.assertTrue(tmpDifferences.isEmpty(), () -> tmpDifferences.size() + " differences:\n"
                + String.join("\n", tmpDifferences.subList(0, Math.min(tmpDifferences.size(), 200))));
    }

    /**
     * Checks whether the given molecule and settings fall under one of the two intended differences to the previous
     * implementation (see {@link #differentialTest()}).
     *
     * @param anASF configured fragmenter
     * @param aMolecule molecule
     * @return name of the intended difference, null if none applies
     * @throws CloneNotSupportedException if cloning fails
     */
    private String getIntendedDifference(AlkylStructureFragmenter anASF, IAtomContainer aMolecule) throws CloneNotSupportedException {
        if (!anASF.getIsolateTertQuatCarbonSetting()) {
            return null;
        }
        IAtomContainer tmpPrepared = AtomContainerManipulator.suppressHydrogens(aMolecule.clone());
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = anASF.classify(tmpPrepared);
        for (int i = 0; i < tmpPrepared.getAtomCount(); i++) {
            if (!tmpClassification.isTertQuatCarbonBondedToRing[i]) {
                continue;
            }
            if (anASF.getSeparateTertQuatCarbonFromRingSetting()) {
                return "A (neighbors of tertiary/quaternary carbons separated from rings are chain atoms)";
            }
            for (IAtom tmpNeighbor : tmpPrepared.getAtom(i).neighbors()) {
                int tmpNeighborIndex = tmpPrepared.indexOf(tmpNeighbor);
                if (tmpClassification.isConjugatedAtom[tmpNeighborIndex] && !tmpClassification.isRingAtom[tmpNeighborIndex]) {
                    return "B (tertiary/quaternary carbons kept at rings keep their bond to conjugated systems)";
                }
            }
        }
        return null;
    }

    /**
     * Fragments a copy of the given molecule with the given fragmenter and returns the sorted unique stereo SMILES of
     * the fragments, or the exception class name if the fragmentation fails.
     *
     * @param aFragmenter fragmenter
     * @param aMolecule molecule to fragment
     * @return sorted fragment SMILES
     */
    private List<String> fragmentToSortedStereoSmiles(IMoleculeFragmenter aFragmenter, IAtomContainer aMolecule) {
        try {
            IAtomContainer tmpCopy = aMolecule.clone();
            aFragmenter.shouldBeFiltered(tmpCopy);
            List<String> tmpSmilesList = new ArrayList<>();
            for (IAtomContainer tmpFragment : aFragmenter.fragmentMolecule(tmpCopy)) {
                tmpSmilesList.add(ChemUtil.createUniqueSmiles(tmpFragment, true, true));
            }
            Collections.sort(tmpSmilesList);
            return tmpSmilesList;
        } catch (Exception anException) {
            return List.of(anException.getClass().getSimpleName());
        }
    }

    /**
     * Converts the given molecule to unique SMILES and parses it again the same way MORTAR does for fragments passed
     * on in a pipeline (see MoleculeDataModel.getAtomContainer()): with kekulization and atom type perception, or,
     * if that fails, without both.
     *
     * @param aMolecule molecule to convert
     * @return the parsed molecule
     * @throws Exception if the SMILES cannot be parsed at all
     */
    private IAtomContainer roundTripAsInPipeline(IAtomContainer aMolecule) throws Exception {
        String tmpSmiles = ChemUtil.createUniqueSmiles(aMolecule, false);
        try {
            return ChemUtil.parseSmilesToAtomContainer(tmpSmiles, true, true);
        } catch (CDKException anException) {
            return ChemUtil.parseSmilesToAtomContainer(tmpSmiles, false, false);
        }
    }

    /**
     * Loads the natural product library.
     *
     * @return list of natural products, titled by their COCONUT identifier
     * @throws Exception if reading or parsing fails
     */
    private List<IAtomContainer> loadNaturalProductLibrary() throws Exception {
        List<IAtomContainer> tmpMolecules = new ArrayList<>();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        try (InputStream tmpStream = this.getClass().getResourceAsStream(AlkylStructureFragmenterIntegrationTest.NATURAL_PRODUCT_LIBRARY)) {
            Assertions.assertNotNull(tmpStream, "Natural product library not found.");
            for (String tmpLine : new String(tmpStream.readAllBytes(), StandardCharsets.UTF_8).split("\n")) {
                if (tmpLine.isBlank()) {
                    continue;
                }
                String[] tmpColumns = tmpLine.split("\t");
                IAtomContainer tmpMolecule = tmpParser.parseSmiles(tmpColumns[0]);
                tmpMolecule.setTitle(tmpColumns.length > 1 ? tmpColumns[1] : "NP#" + (tmpMolecules.size() + 1));
                tmpMolecules.add(tmpMolecule);
            }
        }
        Assertions.assertFalse(tmpMolecules.isEmpty(), "Natural product library could not be read.");
        return tmpMolecules;
    }

    /**
     * Creates a hydrocarbon skeleton of the given molecule: every hetero atom whose bond order sum allows it is
     * replaced by an uncharged carbon atom (e.g. a pyranose ring becomes a cyclohexane ring), all other hetero atoms
     * are removed together with their bonds. Implicit hydrogen counts are recalculated.
     *
     * @param aMolecule molecule to derive the skeleton from
     * @return the hydrocarbon skeleton, possibly disconnected or empty
     * @throws Exception if cloning or hydrogen assignment fails
     */
    private IAtomContainer createCarbonSkeleton(IAtomContainer aMolecule) throws Exception {
        IAtomContainer tmpSkeleton = aMolecule.clone();
        List<IAtom> tmpAtomsToRemove = new ArrayList<>();
        for (IAtom tmpAtom : tmpSkeleton.atoms()) {
            int tmpAtomicNumber = tmpAtom.getAtomicNumber() == null ? 0 : tmpAtom.getAtomicNumber();
            if (tmpAtomicNumber == 1) {
                continue;
            }
            int tmpBondOrderSum = 0;
            for (IBond tmpBond : tmpAtom.bonds()) {
                tmpBondOrderSum += tmpBond.getOrder().numeric();
            }
            if (tmpAtomicNumber == 0 || tmpBondOrderSum > 4) {
                tmpAtomsToRemove.add(tmpAtom);
            } else {
                tmpAtom.setAtomicNumber(6);
                tmpAtom.setSymbol("C");
                tmpAtom.setFormalCharge(0);
                tmpAtom.setAtomTypeName(null);
                tmpAtom.setImplicitHydrogenCount(null);
            }
        }
        for (IAtom tmpAtom : tmpAtomsToRemove) {
            tmpSkeleton.removeAtom(tmpAtom);
        }
        //hydrogens bound to removed hetero atoms are left without neighbors and removed as well
        List<IAtom> tmpIsolatedHydrogens = new ArrayList<>();
        for (IAtom tmpAtom : tmpSkeleton.atoms()) {
            if (tmpAtom.getAtomicNumber() == 1 && tmpAtom.getBondCount() == 0) {
                tmpIsolatedHydrogens.add(tmpAtom);
            }
        }
        tmpIsolatedHydrogens.forEach(tmpSkeleton::removeAtom);
        for (IAtom tmpAtom : tmpSkeleton.atoms()) {
            tmpAtom.setImplicitHydrogenCount(null);
        }
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpSkeleton);
        CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpSkeleton);
        tmpSkeleton.setTitle(aMolecule.getTitle() + "-skeleton");
        return tmpSkeleton;
    }

    /**
     * Checks whether the given molecule consists of carbon and hydrogen atoms only.
     *
     * @param aMolecule molecule to check
     * @return true if hydrocarbon
     */
    private boolean isHydrocarbon(IAtomContainer aMolecule) {
        boolean tmpContainsCarbon = false;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            Integer tmpAtomicNumber = tmpAtom.getAtomicNumber();
            if (tmpAtomicNumber == null || (tmpAtomicNumber != 6 && tmpAtomicNumber != 1)) {
                return false;
            }
            tmpContainsCarbon |= tmpAtomicNumber == 6;
        }
        return tmpContainsCarbon;
    }

    /**
     * Checks invariants of one fragmentation result.
     *
     * @param anASF configured fragmenter
     * @param aMolecule input molecule
     * @return list of violation descriptions, empty if all invariants hold
     */
    private List<String> checkInvariants(AlkylStructureFragmenter anASF, IAtomContainer aMolecule) {
        List<String> tmpViolations = new ArrayList<>();
        String tmpContext = ChemUtil.createUniqueSmiles(aMolecule, false) + " [" + this.describeSettings(anASF) + "]";
        List<IAtomContainer> tmpFragments;
        try {
            Assertions.assertFalse(anASF.shouldBeFiltered(aMolecule), tmpContext);
            tmpFragments = anASF.fragmentMolecule(aMolecule);
        } catch (Exception anException) {
            tmpViolations.add(tmpContext + ": exception " + anException);
            return tmpViolations;
        }
        //index of input atoms and bonds as assigned by the fragmenter (iteration order of the input container)
        Map<Long, IBond> tmpInputBondByAtomPair = new HashMap<>();
        for (IBond tmpBond : aMolecule.bonds()) {
            tmpInputBondByAtomPair.put(this.atomPairKey(aMolecule.indexOf(tmpBond.getBegin()),
                    aMolecule.indexOf(tmpBond.getEnd())), tmpBond);
        }
        int[] tmpAtomOccurrences = new int[aMolecule.getAtomCount()];
        Set<Long> tmpUsedBonds = new HashSet<>();
        List<String> tmpFragmentSmilesList = new ArrayList<>();
        for (IAtomContainer tmpFragment : tmpFragments) {
            String tmpFragmentSmiles = ChemUtil.createUniqueSmiles(tmpFragment, false);
            tmpFragmentSmilesList.add(tmpFragmentSmiles);
            String tmpFragmentContext = tmpContext + " fragment " + tmpFragmentSmiles;
            if (tmpFragment.isEmpty()) {
                tmpViolations.add(tmpFragmentContext + ": empty fragment");
                continue;
            }
            if (!ConnectivityChecker.isConnected(tmpFragment)) {
                tmpViolations.add(tmpFragmentContext + ": fragment is not connected");
            }
            int tmpRealAtomCount = 0;
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (this.isPseudo(tmpAtom)) {
                    continue;
                }
                tmpRealAtomCount++;
                Integer tmpIndex = tmpAtom.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                if (tmpIndex == null) {
                    tmpViolations.add(tmpFragmentContext + ": atom without index");
                    continue;
                }
                tmpAtomOccurrences[tmpIndex]++;
                IAtom tmpInputAtom = aMolecule.getAtom(tmpIndex);
                if (!tmpInputAtom.getAtomicNumber().equals(tmpAtom.getAtomicNumber())) {
                    tmpViolations.add(tmpFragmentContext + ": element changed for atom " + tmpIndex);
                }
                //an atom can never have more (real + pseudo) neighbors in a fragment than in the input
                if (tmpFragment.getConnectedBondsCount(tmpAtom) > tmpInputAtom.getBondCount()) {
                    tmpViolations.add(tmpFragmentContext + ": atom " + tmpIndex + " has more bonds than in input");
                }
                //valence check, implicit hydrogens were (re)assigned by the saturation step
                if (tmpAtom.getAtomicNumber() == 6) {
                    int tmpValence = (tmpAtom.getImplicitHydrogenCount() == null ? 0 : tmpAtom.getImplicitHydrogenCount());
                    for (IBond tmpBond : tmpFragment.getConnectedBondsList(tmpAtom)) {
                        tmpValence += tmpBond.getOrder().numeric();
                    }
                    if (tmpValence > 4) {
                        tmpViolations.add(tmpFragmentContext + ": carbon " + tmpIndex + " has valence " + tmpValence);
                    }
                }
            }
            if (tmpRealAtomCount == 0) {
                tmpViolations.add(tmpFragmentContext + ": fragment consists only of pseudo atoms");
            }
            for (IBond tmpBond : tmpFragment.bonds()) {
                if (tmpBond.getBegin() == null || tmpBond.getEnd() == null) {
                    tmpViolations.add(tmpFragmentContext + ": bond with null atom");
                    continue;
                }
                if (!tmpFragment.contains(tmpBond.getBegin()) || !tmpFragment.contains(tmpBond.getEnd())) {
                    tmpViolations.add(tmpFragmentContext + ": bond atom not contained in fragment");
                }
                if (this.isPseudo(tmpBond.getBegin()) || this.isPseudo(tmpBond.getEnd())) {
                    continue;
                }
                int tmpBeginIndex = tmpBond.getBegin().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                int tmpEndIndex = tmpBond.getEnd().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                long tmpKey = this.atomPairKey(tmpBeginIndex, tmpEndIndex);
                IBond tmpInputBond = tmpInputBondByAtomPair.get(tmpKey);
                if (tmpInputBond == null) {
                    tmpViolations.add(tmpFragmentContext + ": bond " + tmpBeginIndex + "-" + tmpEndIndex + " not in input");
                    continue;
                }
                if (tmpInputBond.getOrder() != tmpBond.getOrder()) {
                    tmpViolations.add(tmpFragmentContext + ": bond order changed for " + tmpBeginIndex + "-" + tmpEndIndex);
                }
                if (!tmpUsedBonds.add(tmpKey)) {
                    tmpViolations.add(tmpFragmentContext + ": bond " + tmpBeginIndex + "-" + tmpEndIndex + " in several fragments");
                }
            }
            //linear, saturated, acyclic chains without pseudo atoms are restricted by the maximum chain length,
            // remnants of tertiary/quaternary systems (contain an atom with three or more bonds in the input) excluded
            if (anASF.getFragmentSideChainsSetting() && this.isLinearAlkane(tmpFragment)
                    && !this.containsBranchingAtomOfInput(tmpFragment, aMolecule)
                    && tmpFragment.getAtomCount() > anASF.getMaxChainLengthSetting()) {
                tmpViolations.add(tmpFragmentContext + ": linear chain longer than " + anASF.getMaxChainLengthSetting());
            }
        }
        for (int i = 0; i < tmpAtomOccurrences.length; i++) {
            if (tmpAtomOccurrences[i] != 1) {
                tmpViolations.add(tmpContext + ": atom " + i + " occurs " + tmpAtomOccurrences[i] + " times " + tmpFragmentSmilesList);
            }
        }
        //the chemistry rules rely on a consistent atom mapping, so they are only checked if the basic invariants hold
        if (tmpViolations.isEmpty()) {
            tmpViolations.addAll(this.checkChemistryRules(anASF, aMolecule, tmpFragments, tmpUsedBonds,
                    tmpContext + " " + tmpFragmentSmilesList));
            tmpViolations.addAll(this.checkStereoChemistry(anASF, aMolecule, tmpFragments, tmpUsedBonds,
                    tmpContext + " " + tmpFragmentSmilesList));
        }
        return tmpViolations;
    }

    /**
     * Checks the transfer of the input molecule's stereo chemistry onto the fragments, independent of the fragmenter's
     * remapping implementation:
     * <ul>
     *     <li>if stereo chemistry is not to be preserved, no fragment carries stereo elements;</li>
     *     <li>a tetrahedral stereo center whose focus atom and all neighbors (carriers) remain bonded in the same
     *     fragment is present in that fragment with an equivalent configuration (compared via the parity of the carrier
     *     permutation);</li>
     *     <li>a double bond stereo element whose double bond and both carrier bonds remain in the same fragment is present
     *     with an equivalent configuration;</li>
     *     <li>every stereo element of a fragment corresponds to such an input stereo element and only references atoms and
     *     bonds of its fragment;</li>
     *     <li>the stereo SMILES of every fragment can be generated.</li>
     * </ul>
     * Statistics about the checked stereo elements are added to {@link #STEREO_STATISTICS}.
     *
     * @param anASF configured fragmenter
     * @param aMolecule input molecule
     * @param aFragmentList fragments of the input molecule
     * @param aUsedBondKeys atom pair keys of all input bonds present in the fragments
     * @param aContext description of the molecule and settings for violation messages
     * @return list of violation descriptions, empty if all rules hold
     */
    private List<String> checkStereoChemistry(AlkylStructureFragmenter anASF, IAtomContainer aMolecule,
                                              List<IAtomContainer> aFragmentList, Set<Long> aUsedBondKeys, String aContext) {
        List<String> tmpViolations = new ArrayList<>();
        if (!anASF.getPreserveStereoChemistrySetting()) {
            for (IAtomContainer tmpFragment : aFragmentList) {
                if (tmpFragment.stereoElements().iterator().hasNext()) {
                    tmpViolations.add(aContext + ": stereo element present although stereo chemistry is not preserved");
                }
            }
            return tmpViolations;
        }
        //fragment index of every input atom
        int[] tmpFragmentIndex = new int[aMolecule.getAtomCount()];
        for (int f = 0; f < aFragmentList.size(); f++) {
            for (IAtom tmpAtom : aFragmentList.get(f).atoms()) {
                if (!this.isPseudo(tmpAtom)) {
                    tmpFragmentIndex[(int) tmpAtom.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY)] = f;
                }
            }
        }
        //validity of the fragments' stereo elements
        Set<IStereoElement> tmpUnmatchedFragmentElements = new HashSet<>();
        for (IAtomContainer tmpFragment : aFragmentList) {
            for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                tmpUnmatchedFragmentElements.add(tmpElement);
                boolean tmpIsContained = tmpElement.getFocus() instanceof IAtom tmpFocusAtom ? tmpFragment.contains(tmpFocusAtom)
                        : tmpFragment.contains((IBond) tmpElement.getFocus());
                for (Object tmpCarrier : tmpElement.getCarriers()) {
                    tmpIsContained &= tmpCarrier instanceof IAtom tmpAtom ? tmpFragment.contains(tmpAtom) : tmpFragment.contains((IBond) tmpCarrier);
                }
                if (!tmpIsContained) {
                    tmpViolations.add(aContext + ": stereo element references atoms or bonds outside of its fragment");
                }
            }
            if (tmpFragment.stereoElements().iterator().hasNext() && ChemUtil.createUniqueSmiles(tmpFragment, true) == null) {
                tmpViolations.add(aContext + ": stereo SMILES of fragment cannot be generated");
            }
        }
        for (IStereoElement tmpInputElement : aMolecule.stereoElements()) {
            if (tmpInputElement instanceof TetrahedralChirality tmpChirality) {
                int tmpFocus = aMolecule.indexOf(tmpChirality.getChiralAtom());
                List<Integer> tmpInputCarriers = new ArrayList<>(4);
                boolean tmpIsComplete = tmpFocus >= 0;
                for (IAtom tmpCarrier : tmpChirality.getLigands()) {
                    int tmpCarrierIndex = aMolecule.indexOf(tmpCarrier);
                    tmpInputCarriers.add(tmpCarrierIndex);
                    tmpIsComplete &= tmpCarrierIndex >= 0 && (tmpCarrierIndex == tmpFocus
                            || aUsedBondKeys.contains(this.atomPairKey(tmpFocus, tmpCarrierIndex)));
                }
                if (!tmpIsComplete) {
                    STEREO_STATISTICS.merge("tetrahedral centers with cut neighbors (not transferable)", 1, Integer::sum);
                    continue;
                }
                IAtomContainer tmpFragment = aFragmentList.get(tmpFragmentIndex[tmpFocus]);
                if (!this.isStereocenterInFragment(tmpFragment, tmpFocus)) {
                    //substituents became identical by fragmentation, the stereo element has to be removed
                    STEREO_STATISTICS.merge("tetrahedral centers no longer stereocenters in fragment", 1, Integer::sum);
                    continue;
                }
                STEREO_STATISTICS.merge("tetrahedral centers expected in fragments", 1, Integer::sum);
                TetrahedralChirality tmpFragmentChirality = null;
                for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                    if (tmpElement instanceof TetrahedralChirality tmpCandidate
                            && Integer.valueOf(tmpFocus).equals(tmpCandidate.getChiralAtom().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY))) {
                        tmpFragmentChirality = tmpCandidate;
                    }
                }
                if (tmpFragmentChirality == null) {
                    tmpViolations.add(String.format("%s: tetrahedral stereo center %d lost", aContext, tmpFocus));
                    continue;
                }
                tmpUnmatchedFragmentElements.remove(tmpFragmentChirality);
                List<Integer> tmpFragmentCarriers = new ArrayList<>(4);
                for (IAtom tmpCarrier : tmpFragmentChirality.getLigands()) {
                    tmpFragmentCarriers.add(tmpCarrier.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY));
                }
                Integer tmpParity = this.permutationParity(tmpInputCarriers, tmpFragmentCarriers);
                if (tmpParity == null) {
                    tmpViolations.add(String.format("%s: stereo center %d has different neighbors in fragment", aContext, tmpFocus));
                } else if ((tmpChirality.getConfigOrder() == tmpFragmentChirality.getConfigOrder()) != (tmpParity == 0)) {
                    tmpViolations.add(String.format("%s: configuration of stereo center %d inverted", aContext, tmpFocus));
                } else {
                    STEREO_STATISTICS.merge("tetrahedral centers transferred correctly", 1, Integer::sum);
                }
            } else if (tmpInputElement instanceof DoubleBondStereochemistry tmpDoubleBondStereo) {
                IBond tmpFocusBond = tmpDoubleBondStereo.getStereoBond();
                int tmpBegin = aMolecule.indexOf(tmpFocusBond.getBegin());
                int tmpEnd = aMolecule.indexOf(tmpFocusBond.getEnd());
                IBond[] tmpCarrierBonds = tmpDoubleBondStereo.getBonds();
                boolean tmpIsComplete = tmpBegin >= 0 && tmpEnd >= 0 && aUsedBondKeys.contains(this.atomPairKey(tmpBegin, tmpEnd));
                int[] tmpSubstituents = new int[2];
                for (int i = 0; i < 2; i++) {
                    IAtom tmpSideAtom = i == 0 ? tmpFocusBond.getBegin() : tmpFocusBond.getEnd();
                    int tmpSideIndex = aMolecule.indexOf(tmpSideAtom);
                    tmpSubstituents[i] = aMolecule.indexOf(tmpCarrierBonds[i].getOther(tmpSideAtom));
                    tmpIsComplete &= tmpSubstituents[i] >= 0 && aUsedBondKeys.contains(this.atomPairKey(tmpSideIndex, tmpSubstituents[i]));
                }
                if (!tmpIsComplete) {
                    STEREO_STATISTICS.merge("double bond stereo with cut substituents (not transferable)", 1, Integer::sum);
                    continue;
                }
                IAtomContainer tmpFragment = aFragmentList.get(tmpFragmentIndex[tmpBegin]);
                if (!this.isStereocenterInFragment(tmpFragment, tmpBegin) || !this.isStereocenterInFragment(tmpFragment, tmpEnd)) {
                    STEREO_STATISTICS.merge("double bonds no longer stereogenic in fragment", 1, Integer::sum);
                    continue;
                }
                STEREO_STATISTICS.merge("double bond stereo expected in fragments", 1, Integer::sum);
                DoubleBondStereochemistry tmpFragmentStereo = null;
                for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                    if (tmpElement instanceof DoubleBondStereochemistry tmpCandidate) {
                        int tmpCandidateBegin = tmpCandidate.getStereoBond().getBegin().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                        int tmpCandidateEnd = tmpCandidate.getStereoBond().getEnd().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                        if (this.atomPairKey(tmpCandidateBegin, tmpCandidateEnd) == this.atomPairKey(tmpBegin, tmpEnd)) {
                            tmpFragmentStereo = tmpCandidate;
                        }
                    }
                }
                if (tmpFragmentStereo == null) {
                    tmpViolations.add(String.format("%s: double bond stereo %d=%d lost", aContext, tmpBegin, tmpEnd));
                    continue;
                }
                tmpUnmatchedFragmentElements.remove(tmpFragmentStereo);
                //substituents of the fragment element, oriented like the input double bond
                IBond tmpFragmentFocus = tmpFragmentStereo.getStereoBond();
                boolean tmpIsSameOrientation = Integer.valueOf(tmpBegin).equals(
                        tmpFragmentFocus.getBegin().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY));
                int[] tmpFragmentSubstituents = new int[2];
                for (int i = 0; i < 2; i++) {
                    IAtom tmpSideAtom = i == 0 ? tmpFragmentFocus.getBegin() : tmpFragmentFocus.getEnd();
                    IAtom tmpOther = tmpFragmentStereo.getBonds()[i].getOther(tmpSideAtom);
                    Integer tmpOtherIndex = tmpOther == null ? null : tmpOther.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
                    tmpFragmentSubstituents[tmpIsSameOrientation ? i : 1 - i] = tmpOtherIndex == null ? -1 : tmpOtherIndex;
                }
                int tmpFlips = (tmpFragmentSubstituents[0] == tmpSubstituents[0] ? 0 : 1) + (tmpFragmentSubstituents[1] == tmpSubstituents[1] ? 0 : 1);
                if ((tmpDoubleBondStereo.getConfigOrder() == tmpFragmentStereo.getConfigOrder()) != (tmpFlips % 2 == 0)) {
                    tmpViolations.add(String.format("%s: configuration of double bond %d=%d inverted", aContext, tmpBegin, tmpEnd));
                } else {
                    STEREO_STATISTICS.merge("double bond stereo transferred correctly", 1, Integer::sum);
                }
            } else {
                STEREO_STATISTICS.merge("other stereo element types (not transferred by design): " + tmpInputElement.getClass().getSimpleName(), 1, Integer::sum);
            }
        }
        for (IStereoElement tmpElement : tmpUnmatchedFragmentElements) {
            tmpViolations.add(aContext + ": unexpected stereo element in fragment: " + tmpElement.getClass().getSimpleName());
        }
        return tmpViolations;
    }

    /**
     * Checks with the symmetry based stereocenter detection of the CDK whether the input atom with the given index is a
     * stereocenter in the given fragment (on a copy with recalculated implicit hydrogen counts).
     *
     * @param aFragment fragment containing the atom
     * @param anInputAtomIndex index of the atom in the input molecule
     * @return true if the atom is a stereocenter in the fragment
     */
    private boolean isStereocenterInFragment(IAtomContainer aFragment, int anInputAtomIndex) {
        for (int i = 0; i < aFragment.getAtomCount(); i++) {
            if (Integer.valueOf(anInputAtomIndex).equals(aFragment.getAtom(i).getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY))) {
                return Stereocenters.of(aFragment).isStereocenter(i);
            }
        }
        return false;
    }

    /**
     * Determines the parity of the permutation transforming the first into the second list.
     *
     * @param aFirstList first list
     * @param aSecondList second list with the same elements
     * @return 0 for an even, 1 for an odd permutation, null if the lists do not contain the same elements
     */
    private Integer permutationParity(List<Integer> aFirstList, List<Integer> aSecondList) {
        if (aFirstList.size() != aSecondList.size() || !new HashSet<>(aFirstList).equals(new HashSet<>(aSecondList))
                || new HashSet<>(aFirstList).size() != aFirstList.size()) {
            return null;
        }
        List<Integer> tmpWorkList = new ArrayList<>(aSecondList);
        int tmpSwaps = 0;
        for (int i = 0; i < aFirstList.size(); i++) {
            int j = tmpWorkList.indexOf(aFirstList.get(i));
            if (j != i) {
                Collections.swap(tmpWorkList, i, j);
                tmpSwaps++;
            }
        }
        return tmpSwaps % 2;
    }

    /**
     * Checks the chemical rules of the alkyl fragmentation on one fragmentation result. All structural properties of
     * the input (rings, multi bonds, conjugation, tertiary/quaternary carbons) are determined independently of the
     * fragmenter's internal markers:
     * <ul>
     *     <li>ring bonds, double/triple bonds, bonds between aromatic atoms (of any order, also without bond order) and
     *     single bonds of conjugated systems (D-S-D, T-S-T, aromatic atoms counting like double bond atoms) are never
     *     cut; aromatic atoms are never tertiary/quaternary carbons;</li>
     *     <li>with isolation of tertiary/quaternary carbons, such a carbon is returned as single carbon saturated with
     *     pseudo atoms, unless it is bonded to a ring and not to be separated from it;</li>
     *     <li>without isolation (or when kept at a ring), a tertiary/quaternary carbon keeps its saturated non-ring
     *     neighbors; every cut bond of such a carbon (e.g. to an isolated multi bond, an isolated
     *     tertiary/quaternary carbon or a ring it is separated from) is marked with a pseudo atom;</li>
     *     <li>a tertiary/quaternary carbon bonded to a ring stays bonded to it, or is not part of any fragment containing
     *     ring atoms if it is to be separated from rings;</li>
     *     <li>every fragment belongs to exactly one substructure class: ring/conjugated system (only with attached
     *     multi bond atoms or tertiary/quaternary systems), tertiary/quaternary system, isolated multi bond or linear
     *     chain;</li>
     *     <li>pseudo atoms only saturate tertiary/quaternary carbons.</li>
     * </ul>
     *
     * @param anASF configured fragmenter
     * @param aMolecule input molecule
     * @param aFragmentList fragments of the input molecule
     * @param aUsedBondKeys atom pair keys of all input bonds present in the fragments
     * @param aContext description of the molecule and settings for violation messages
     * @return list of violation descriptions, empty if all rules hold
     */
    private List<String> checkChemistryRules(AlkylStructureFragmenter anASF, IAtomContainer aMolecule,
                                             List<IAtomContainer> aFragmentList, Set<Long> aUsedBondKeys, String aContext) {
        List<String> tmpViolations = new ArrayList<>();
        int tmpAtomCount = aMolecule.getAtomCount();
        //<editor-fold desc="Input properties">
        IAtomContainer tmpRingMolecule;
        try {
            tmpRingMolecule = aMolecule.clone();
        } catch (CloneNotSupportedException anException) {
            tmpViolations.add(aContext + ": clone failed");
            return tmpViolations;
        }
        Cycles.markRingAtomsAndBonds(tmpRingMolecule);
        boolean[] tmpIsRing = new boolean[tmpAtomCount];
        boolean[] tmpHasDouble = new boolean[tmpAtomCount];
        boolean[] tmpHasTriple = new boolean[tmpAtomCount];
        //aromatic bonds without bond order (non-kekulizable aromatic remnants, e.g. from a previous Ertl fragmentation)
        boolean[] tmpHasAromaticWithoutOrder = new boolean[tmpAtomCount];
        boolean[] tmpIsTertQuat = new boolean[tmpAtomCount];
        for (int i = 0; i < tmpAtomCount; i++) {
            IAtom tmpAtom = aMolecule.getAtom(i);
            tmpIsRing[i] = tmpRingMolecule.getAtom(i).isInRing();
            int tmpHeavyDegree = 0;
            boolean tmpOnlySingleBonds = true;
            for (IBond tmpBond : tmpAtom.bonds()) {
                tmpHasDouble[i] |= tmpBond.getOrder() == IBond.Order.DOUBLE;
                tmpHasTriple[i] |= tmpBond.getOrder() == IBond.Order.TRIPLE;
                tmpHasAromaticWithoutOrder[i] |= this.isAromaticWithoutOrder(tmpBond);
                tmpOnlySingleBonds &= tmpBond.getOrder() == IBond.Order.SINGLE;
                if (tmpBond.getOther(tmpAtom).getAtomicNumber() != 1) {
                    tmpHeavyDegree++;
                }
            }
            //aromatic atoms (also remnants without aromatic bonds) are never tertiary/quaternary carbons
            tmpIsTertQuat[i] = !tmpIsRing[i] && !tmpAtom.isAromatic() && tmpOnlySingleBonds && (tmpHeavyDegree == 3 || tmpHeavyDegree == 4);
        }
        //conjugated single bonds and atoms of conjugated systems (extended along the multi bonds)
        Set<Long> tmpConjugatedSingleBonds = new HashSet<>();
        boolean[] tmpIsConjugated = new boolean[tmpAtomCount];
        //aromatic atoms are always conjugated, bonds between aromatic atoms (of any order) are never cut
        Set<Long> tmpAromaticBonds = new HashSet<>();
        for (int i = 0; i < tmpAtomCount; i++) {
            tmpIsConjugated[i] = aMolecule.getAtom(i).isAromatic();
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            boolean tmpIsBeginPi = tmpHasDouble[tmpBegin] || tmpHasAromaticWithoutOrder[tmpBegin] || tmpBond.getBegin().isAromatic();
            boolean tmpIsEndPi = tmpHasDouble[tmpEnd] || tmpHasAromaticWithoutOrder[tmpEnd] || tmpBond.getEnd().isAromatic();
            if (this.isAromaticWithoutOrder(tmpBond) || (tmpBond.getBegin().isAromatic() && tmpBond.getEnd().isAromatic())) {
                tmpAromaticBonds.add(this.atomPairKey(tmpBegin, tmpEnd));
                tmpIsConjugated[tmpBegin] = true;
                tmpIsConjugated[tmpEnd] = true;
            } else if (tmpBond.getOrder() == IBond.Order.SINGLE
                    && ((tmpIsBeginPi && tmpIsEndPi) || (tmpHasTriple[tmpBegin] && tmpHasTriple[tmpEnd]))) {
                tmpConjugatedSingleBonds.add(this.atomPairKey(tmpBegin, tmpEnd));
                tmpIsConjugated[tmpBegin] = true;
                tmpIsConjugated[tmpEnd] = true;
            }
        }
        for (IBond tmpBond : aMolecule.bonds()) {
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            if (tmpBond.getOrder() != IBond.Order.SINGLE && (tmpIsConjugated[tmpBegin] || tmpIsConjugated[tmpEnd])) {
                tmpIsConjugated[tmpBegin] = true;
                tmpIsConjugated[tmpEnd] = true;
            }
        }
        //fragment index of every input atom
        int[] tmpFragmentIndex = new int[tmpAtomCount];
        for (int f = 0; f < aFragmentList.size(); f++) {
            for (IAtom tmpAtom : aFragmentList.get(f).atoms()) {
                if (!this.isPseudo(tmpAtom)) {
                    tmpFragmentIndex[(int) tmpAtom.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY)] = f;
                }
            }
        }
        //</editor-fold>
        //<editor-fold desc="Bond rules">
        for (IBond tmpBond : aMolecule.bonds()) {
            int tmpBegin = aMolecule.indexOf(tmpBond.getBegin());
            int tmpEnd = aMolecule.indexOf(tmpBond.getEnd());
            long tmpKey = this.atomPairKey(tmpBegin, tmpEnd);
            if (aUsedBondKeys.contains(tmpKey)) {
                continue;
            }
            if (tmpRingMolecule.getBond(tmpRingMolecule.getAtom(tmpBegin), tmpRingMolecule.getAtom(tmpEnd)).isInRing()) {
                tmpViolations.add(String.format("%s: ring bond %d-%d cut", aContext, tmpBegin, tmpEnd));
            } else if (tmpAromaticBonds.contains(tmpKey)) {
                tmpViolations.add(String.format("%s: aromatic bond %d-%d cut", aContext, tmpBegin, tmpEnd));
            } else if (tmpBond.getOrder() != IBond.Order.SINGLE) {
                tmpViolations.add(String.format("%s: multi bond %d-%d cut", aContext, tmpBegin, tmpEnd));
            } else if (tmpConjugatedSingleBonds.contains(tmpKey)) {
                tmpViolations.add(String.format("%s: conjugated system cut at bond %d-%d", aContext, tmpBegin, tmpEnd));
            }
        }
        //</editor-fold>
        //<editor-fold desc="Tertiary/quaternary carbon rules">
        boolean tmpIsolate = anASF.getIsolateTertQuatCarbonSetting();
        boolean tmpSeparateFromRing = anASF.getSeparateTertQuatCarbonFromRingSetting();
        //a tertiary/quaternary carbon is isolated if the setting is active, unless it is bonded to a ring and not to be
        // separated from it
        boolean[] tmpIsIsolatedTertQuat = new boolean[tmpAtomCount];
        for (int i = 0; i < tmpAtomCount; i++) {
            boolean tmpIsBondedToRing = false;
            for (IAtom tmpNeighbor : aMolecule.getAtom(i).neighbors()) {
                tmpIsBondedToRing |= tmpIsRing[aMolecule.indexOf(tmpNeighbor)];
            }
            tmpIsIsolatedTertQuat[i] = tmpIsTertQuat[i] && tmpIsolate && !(tmpIsBondedToRing && !tmpSeparateFromRing);
        }
        for (int i = 0; i < tmpAtomCount; i++) {
            if (!tmpIsTertQuat[i]) {
                continue;
            }
            IAtom tmpAtom = aMolecule.getAtom(i);
            IAtomContainer tmpFragment = aFragmentList.get(tmpFragmentIndex[i]);
            IAtom tmpFragmentAtom = null;
            int tmpPseudoCount = 0;
            for (IAtom tmpCandidate : tmpFragment.atoms()) {
                if (this.isPseudo(tmpCandidate)) {
                    tmpPseudoCount++;
                } else if ((int) tmpCandidate.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY) == i) {
                    tmpFragmentAtom = tmpCandidate;
                }
            }
            if (tmpIsIsolatedTertQuat[i]) {
                if (tmpFragment.getAtomCount() - tmpPseudoCount != 1 || tmpPseudoCount != tmpAtom.getBondCount()) {
                    tmpViolations.add(String.format("%s: tertiary/quaternary carbon %d not isolated", aContext, i));
                }
                continue;
            }
            //every cut bond is marked by a pseudo atom, the carbon keeps its original number of bonds
            int tmpKeptBondCount = 0;
            for (IAtom tmpNeighbor : tmpAtom.neighbors()) {
                int j = aMolecule.indexOf(tmpNeighbor);
                boolean tmpIsBondKept = aUsedBondKeys.contains(this.atomPairKey(i, j));
                tmpKeptBondCount += tmpIsBondKept ? 1 : 0;
                if (tmpIsRing[j]) {
                    if (tmpSeparateFromRing == tmpIsBondKept) {
                        tmpViolations.add(String.format("%s: tertiary/quaternary carbon %d %s ring atom %d", aContext, i,
                                tmpSeparateFromRing ? "not separated from" : "separated from", j));
                    }
                } else if (tmpIsIsolatedTertQuat[j]) {
                    if (tmpIsBondKept) {
                        tmpViolations.add(String.format("%s: tertiary/quaternary carbon %d bonded to isolated carbon %d", aContext, i, j));
                    }
                } else if (!tmpHasDouble[j] && !tmpHasTriple[j] && !tmpIsConjugated[j] && !tmpIsBondKept) {
                    tmpViolations.add(String.format("%s: tertiary/quaternary carbon %d lost its neighbor %d", aContext, i, j));
                }
            }
            int tmpExpectedPseudoNeighbors = tmpAtom.getBondCount() - tmpKeptBondCount;
            int tmpActualPseudoNeighbors = 0;
            for (IAtom tmpNeighbor : tmpFragment.getConnectedAtomsList(tmpFragmentAtom)) {
                tmpActualPseudoNeighbors += this.isPseudo(tmpNeighbor) ? 1 : 0;
            }
            if (tmpActualPseudoNeighbors != tmpExpectedPseudoNeighbors) {
                tmpViolations.add(String.format("%s: tertiary/quaternary carbon %d has %d pseudo atoms, expected %d",
                        aContext, i, tmpActualPseudoNeighbors, tmpExpectedPseudoNeighbors));
            }
        }
        //</editor-fold>
        //<editor-fold desc="Fragment class and pseudo atom rules">
        for (IAtomContainer tmpFragment : aFragmentList) {
            String tmpFragmentContext = aContext + " fragment " + ChemUtil.createUniqueSmiles(tmpFragment, false);
            Set<Integer> tmpRealAtoms = new HashSet<>();
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (this.isPseudo(tmpAtom)) {
                    List<IAtom> tmpPseudoNeighbors = tmpFragment.getConnectedAtomsList(tmpAtom);
                    if (tmpPseudoNeighbors.size() != 1 || this.isPseudo(tmpPseudoNeighbors.getFirst())
                            || !tmpIsTertQuat[(int) tmpPseudoNeighbors.getFirst().getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY)]) {
                        tmpViolations.add(tmpFragmentContext + ": pseudo atom not bonded to a tertiary/quaternary carbon");
                    }
                } else {
                    tmpRealAtoms.add(tmpAtom.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY));
                }
            }
            boolean tmpContainsRingOrConjugated = false;
            boolean tmpContainsTertQuat = false;
            boolean tmpContainsMultiBondAtom = false;
            for (int i : tmpRealAtoms) {
                tmpContainsRingOrConjugated |= tmpIsRing[i] || tmpIsConjugated[i];
                tmpContainsTertQuat |= tmpIsTertQuat[i];
                tmpContainsMultiBondAtom |= tmpHasDouble[i] || tmpHasTriple[i];
            }
            for (int i : tmpRealAtoms) {
                boolean tmpHasMulti = tmpHasDouble[i] || tmpHasTriple[i];
                boolean tmpIsTertQuatSystemAtom = tmpIsTertQuat[i];
                for (IAtom tmpNeighbor : aMolecule.getAtom(i).neighbors()) {
                    int j = aMolecule.indexOf(tmpNeighbor);
                    tmpIsTertQuatSystemAtom |= tmpIsTertQuat[j] && tmpRealAtoms.contains(j);
                }
                String tmpViolation = null;
                if (tmpContainsRingOrConjugated) {
                    if (!tmpIsRing[i] && !tmpIsConjugated[i] && !tmpHasMulti && !tmpIsTertQuatSystemAtom) {
                        tmpViolation = "alkyl chain atom %d attached to ring/conjugated system";
                    }
                } else if (tmpContainsTertQuat) {
                    if (!tmpIsTertQuatSystemAtom) {
                        tmpViolation = "atom %d not part of the tertiary/quaternary system";
                    }
                } else if (tmpContainsMultiBondAtom) {
                    if (!tmpHasMulti) {
                        tmpViolation = "saturated atom %d attached to isolated multi bond";
                    }
                }
                if (tmpViolation != null) {
                    tmpViolations.add(tmpFragmentContext + ": " + String.format(tmpViolation, i));
                }
            }
        }
        //</editor-fold>
        return tmpViolations;
    }

    /**
     * Loads test SMILES and the structures of the substructure library.
     *
     * @return list of test molecules
     * @throws Exception if parsing or reading fails
     */
    private List<IAtomContainer> loadTestMolecules() throws Exception {
        List<IAtomContainer> tmpMolecules = new ArrayList<>();
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        for (String tmpSmiles : AlkylStructureFragmenterIntegrationTest.TEST_SMILES) {
            IAtomContainer tmpMolecule = tmpParser.parseSmiles(tmpSmiles);
            tmpMolecule.setTitle(tmpSmiles);
            tmpMolecules.add(tmpMolecule);
        }
        try (InputStream tmpStream = this.getClass().getResourceAsStream("ASF_Test_Substructure_Library.sdf");
             IteratingSDFReader tmpReader = new IteratingSDFReader(tmpStream, SilentChemObjectBuilder.getInstance())) {
            while (tmpReader.hasNext()) {
                tmpMolecules.add(tmpReader.next());
            }
        }
        return tmpMolecules;
    }

    /**
     * Creates fragmenter instances for all combinations of the boolean settings and the tested maximum chain lengths.
     *
     * @return list of configured fragmenters
     */
    private List<AlkylStructureFragmenter> createAllSettingCombinations() {
        List<AlkylStructureFragmenter> tmpFragmenters = new ArrayList<>();
        for (int tmpMask = 0; tmpMask < 8; tmpMask++) {
            for (int tmpMaxChainLength : AlkylStructureFragmenterIntegrationTest.MAX_CHAIN_LENGTHS) {
                boolean tmpFragmentSideChains = (tmpMask & 1) != 0;
                if (!tmpFragmentSideChains && tmpMaxChainLength != AlkylStructureFragmenter.MAX_CHAIN_LENGTH_SETTING_DEFAULT) {
                    //max chain length has no effect without side chain fragmentation
                    continue;
                }
                AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
                tmpASF.setFragmentSideChainsSetting(tmpFragmentSideChains);
                tmpASF.setIsolateTertQuatCarbonsSetting((tmpMask & 2) != 0);
                tmpASF.setSeparateTertQuatCarbonFromRingSetting((tmpMask & 4) != 0);
                tmpASF.setMaxChainLengthSetting(tmpMaxChainLength);
                tmpFragmenters.add(tmpASF);
            }
        }
        return tmpFragmenters;
    }

    /**
     * Fragments the given molecule and returns the sorted unique SMILES of the fragments.
     *
     * @param anASF configured fragmenter
     * @param aMolecule molecule to fragment
     * @return sorted fragment SMILES
     * @throws Exception if fragmentation fails
     */
    private List<String> fragmentToSortedSmiles(AlkylStructureFragmenter anASF, IAtomContainer aMolecule) throws Exception {
        anASF.shouldBeFiltered(aMolecule);
        List<String> tmpSmilesList = new ArrayList<>();
        for (IAtomContainer tmpFragment : anASF.fragmentMolecule(aMolecule)) {
            //aromaticity is perceived so that different Kekulé structures of the same fragment give the same SMILES
            IAtomContainer tmpClone = tmpFragment.clone();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpClone);
            AlkylStructureFragmenterIntegrationTest.AROMATICITY.apply(tmpClone);
            tmpSmilesList.add(ChemUtil.createUniqueSmiles(tmpClone, false, true));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }

    /**
     * Same as {@link #fragmentToSortedSmiles(AlkylStructureFragmenter, IAtomContainer)} but returns the exception
     * class name as single list element instead of throwing it.
     *
     * @param anASF configured fragmenter
     * @param aMolecule molecule to fragment
     * @return sorted fragment SMILES or exception name
     */
    private List<String> fragmentToSortedSmilesOrException(AlkylStructureFragmenter anASF, IAtomContainer aMolecule) {
        try {
            return this.fragmentToSortedSmiles(anASF, aMolecule);
        } catch (Exception anException) {
            return List.of(anException.getClass().getSimpleName());
        }
    }

    /**
     * Creates a copy of the given molecule without stereo elements, optionally with randomly permuted atoms and
     * bonds and randomly swapped bond orientations.
     *
     * @param aMolecule molecule to copy
     * @param aRandom random source for the permutation, no permutation if null
     * @return the copy
     * @throws CloneNotSupportedException if cloning fails
     */
    private IAtomContainer copyWithoutStereo(IAtomContainer aMolecule, Random aRandom) throws CloneNotSupportedException {
        IAtomContainer tmpClone = aMolecule.clone();
        tmpClone.setStereoElements(new ArrayList<>());
        if (aRandom == null) {
            return tmpClone;
        }
        List<IAtom> tmpAtoms = new ArrayList<>();
        tmpClone.atoms().forEach(tmpAtoms::add);
        List<IBond> tmpBonds = new ArrayList<>();
        tmpClone.bonds().forEach(tmpBonds::add);
        Collections.shuffle(tmpAtoms, aRandom);
        Collections.shuffle(tmpBonds, aRandom);
        IAtomContainer tmpPermuted = SilentChemObjectBuilder.getInstance().newAtomContainer();
        for (IAtom tmpAtom : tmpAtoms) {
            tmpPermuted.addAtom(tmpAtom);
        }
        for (IBond tmpBond : tmpBonds) {
            IAtom tmpBegin = tmpPermuted.getAtom(tmpAtoms.indexOf(tmpBond.getBegin()));
            IAtom tmpEnd = tmpPermuted.getAtom(tmpAtoms.indexOf(tmpBond.getEnd()));
            IBond tmpNewBond = aRandom.nextBoolean()
                    ? SilentChemObjectBuilder.getInstance().newInstance(IBond.class, tmpEnd, tmpBegin, tmpBond.getOrder())
                    : SilentChemObjectBuilder.getInstance().newInstance(IBond.class, tmpBegin, tmpEnd, tmpBond.getOrder());
            tmpNewBond.setIsAromatic(tmpBond.isAromatic());
            tmpPermuted.addBond(tmpNewBond);
        }
        return tmpPermuted;
    }

    /**
     * Checks whether the fragment is an unbranched, acyclic chain of single bonded carbons without pseudo atoms.
     *
     * @param aFragment fragment to check
     * @return true if linear alkane
     */
    private boolean isLinearAlkane(IAtomContainer aFragment) {
        if (aFragment.getBondCount() != aFragment.getAtomCount() - 1) {
            return false;
        }
        for (IAtom tmpAtom : aFragment.atoms()) {
            if (this.isPseudo(tmpAtom) || aFragment.getConnectedBondsCount(tmpAtom) > 2) {
                return false;
            }
        }
        for (IBond tmpBond : aFragment.bonds()) {
            if (tmpBond.getOrder() != IBond.Order.SINGLE) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks whether the fragment contains an atom that has three or more bonds in the input molecule.
     *
     * @param aFragment fragment to check
     * @param aMolecule input molecule
     * @return true if a branching atom of the input is part of the fragment
     */
    private boolean containsBranchingAtomOfInput(IAtomContainer aFragment, IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aFragment.atoms()) {
            Integer tmpIndex = tmpAtom.getProperty(AlkylStructureFragmenter.INTERNAL_ASF_ATOM_INDEX_PROPERTY_KEY);
            if (tmpIndex != null && aMolecule.getAtom(tmpIndex).getBondCount() >= 3) {
                return true;
            }
        }
        return false;
    }

    private boolean isAromaticWithoutOrder(IBond aBond) {
        return aBond.isAromatic() && (aBond.getOrder() == null || aBond.getOrder() == IBond.Order.UNSET);
    }

    private boolean isPseudo(IAtom anAtom) {
        return anAtom instanceof IPseudoAtom || anAtom.getAtomicNumber() == null || anAtom.getAtomicNumber() == 0;
    }

    private long atomPairKey(int anIndex1, int anIndex2) {
        return ((long) Math.min(anIndex1, anIndex2) << 32) | Math.max(anIndex1, anIndex2);
    }

    private String describeSettings(AlkylStructureFragmenter anASF) {
        return String.format("sideChains=%b, maxLength=%d, isolateTertQuat=%b, separateTertQuatFromRing=%b",
                anASF.getFragmentSideChainsSetting(), anASF.getMaxChainLengthSetting(),
                anASF.getIsolateTertQuatCarbonSetting(), anASF.getSeparateTertQuatCarbonFromRingSetting());
    }
}
