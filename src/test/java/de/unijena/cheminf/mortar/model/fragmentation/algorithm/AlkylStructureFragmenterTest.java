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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Unit tests for {@link AlkylStructureFragmenter}: the filter and settings API, the classification of atoms and bonds,
 * and the expected fragments of small example molecules for every fragmentation rule. The rules themselves are
 * additionally checked on a large set of generated molecules in {@link AlkylStructureFragmenterIntegrationTest}.
 * <p>
 * Most fragmentation tests are given as tables with the columns: input SMILES, isolation of tertiary/quaternary
 * carbons, separation of tertiary/quaternary carbons from rings, maximum chain length (0 for no dissection of linear
 * chains) and the expected fragments as SMILES separated by spaces. Expected and actual fragments are compared as
 * sorted unique SMILES.
 * </p>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class AlkylStructureFragmenterTest {
    /**
     * Static setter for Locale.
     * If Locale is set in Test constructor, a MissingResourceException is thrown when user locale diverts from set Locale.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    //<editor-fold desc="Example usage">
    /**
     * Showcases the usage of the AlkylStructureFragmenter: a pure hydrocarbon is checked with the filter and
     * preprocessing methods and fragmented with the default settings.
     *
     * @throws Exception if the SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void exampleUsageTest() throws Exception {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC(C)(C)CCCCCCCC1CCCCC1");
        //constructs an AlkylStructureFragmenter with default settings
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        //pure hydrocarbons are not filtered and need no preprocessing
        Assertions.assertFalse(tmpASF.shouldBeFiltered(tmpMolecule));
        Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpMolecule));
        Assertions.assertTrue(tmpASF.canBeFragmented(tmpMolecule));
        List<IAtomContainer> tmpFragments = tmpASF.fragmentMolecule(tmpMolecule);
        //isolated quaternary carbon with its methyl groups, the chain of seven carbons dissected into chains of at most
        // six carbons, and the ring
        Assertions.assertEquals(this.toSortedUniqueSmiles("*C(*)(*)* C C C CCCCCC C C1CCCCC1", false),
                this.toSortedUniqueSmiles(tmpFragments, false));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Filter and settings">
    /**
     * Tests which molecules are filtered, with and without keeping non-fragmentable molecules in the pipeline:
     * molecules with hetero atoms, pseudo atoms or without carbon atoms are non-fragmentable, explicit hydrogen atoms
     * are allowed.
     *
     * @throws InvalidSmilesException if a SMILES cannot be parsed
     */
    @Test
    public void filterTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //columns: SMILES, fragmentable
        Object[][] tmpCases = {
                {"CCC", true},
                {"[H]C([H])([H])[H]", true},
                {"C=CC#C", true},
                {"CN", false},
                {"C*", false},
                {"[H][H]", false}};
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (boolean tmpKeepNonFragmentable : new boolean[]{true, false}) {
            tmpASF.setKeepNonFragmentableMoleculesSetting(tmpKeepNonFragmentable);
            for (Object[] tmpCase : tmpCases) {
                IAtomContainer tmpMolecule = tmpParser.parseSmiles((String) tmpCase[0]);
                boolean tmpIsFragmentable = (boolean) tmpCase[1];
                String tmpMessage = tmpCase[0] + ", keep non-fragmentable molecules: " + tmpKeepNonFragmentable;
                Assertions.assertEquals(!tmpIsFragmentable && !tmpKeepNonFragmentable, tmpASF.shouldBeFiltered(tmpMolecule), tmpMessage);
                Assertions.assertEquals(tmpIsFragmentable, tmpASF.canBeFragmented(tmpMolecule), tmpMessage);
                Assertions.assertFalse(tmpASF.shouldBePreprocessed(tmpMolecule), tmpMessage);
            }
        }
        Assertions.assertTrue(tmpASF.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
    }

    /**
     * Tests the handling of non-fragmentable molecules: they are returned unchanged as single fragment if they are to
     * be kept in the pipeline (default), otherwise their fragmentation is rejected.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void nonFragmentableMoleculesTest() throws Exception {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (String tmpSmiles : new String[]{"CCO", "C*", "c1ccncc1"}) {
            IAtomContainer tmpMolecule = tmpParser.parseSmiles(tmpSmiles);
            tmpASF.setKeepNonFragmentableMoleculesSetting(true);
            List<IAtomContainer> tmpFragments = tmpASF.fragmentMolecule(tmpMolecule);
            Assertions.assertEquals(1, tmpFragments.size(), tmpSmiles);
            Assertions.assertEquals(ChemUtil.createUniqueSmiles(tmpMolecule, false), ChemUtil.createUniqueSmiles(tmpFragments.getFirst(), false));
            tmpASF.setKeepNonFragmentableMoleculesSetting(false);
            Assertions.assertThrows(IllegalArgumentException.class, () -> tmpASF.fragmentMolecule(tmpMolecule), tmpSmiles);
        }
    }

    /**
     * Test for the handling of molecules containing null atoms. Since the CDK does not allow the addition of null atoms
     * to existing atom container instances, reflection is used in this test to manipulate the inner fields of the atom
     * container. Such molecules are always filtered, independent of the setting for keeping non-fragmentable molecules.
     *
     * @throws InvalidSmilesException if parser cannot parse the given SMILES
     * @throws NoSuchFieldException if either the atoms or bonds field can not be found
     * @throws IllegalAccessException if reflection can not be done
     */
    @Test
    public void nullAtomInputTest() throws InvalidSmilesException, NoSuchFieldException, IllegalAccessException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCCCCCCCCCCC");
        //reflection is used to manually insert/inflict a null value
        Field tmpAtomsField = AtomContainer.class.getDeclaredField("atoms");
        Field tmpBondsField = AtomContainer.class.getDeclaredField("bonds");
        tmpAtomsField.setAccessible(true);
        tmpBondsField.setAccessible(true);
        //"getter" for BaseAtomRef and BondRef
        Object tmpAtoms = tmpAtomsField.get(tmpMolecule);
        Object[] tmpNewAtoms = (Object[]) Array.newInstance(tmpAtoms.getClass().getComponentType(), Array.getLength(tmpAtoms));
        Object tmpBonds = tmpBondsField.get(tmpMolecule);
        Object[] tmpNewBonds = (Object[]) Array.newInstance(tmpBonds.getClass().getComponentType(), Array.getLength(tmpBonds));
        //fill tmpNewAtoms with actual atoms
        for (int i = 0; i < 12; i++) {
            tmpNewAtoms[i] = tmpMolecule.getAtom(i);
            if (i < 11) {
                tmpNewBonds[i] = tmpMolecule.getBond(i);
            }
        }
        //set position 5 null
        tmpNewAtoms[5] = null;
        tmpNewBonds[5] = null;
        //insert the fields back into tmpMolecule
        tmpAtomsField.set(tmpMolecule, tmpNewAtoms);
        tmpBondsField.set(tmpMolecule, tmpNewBonds);
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        for (boolean tmpKeepNonFragmentable : new boolean[]{false, true}) {
            tmpASF.setKeepNonFragmentableMoleculesSetting(tmpKeepNonFragmentable);
            Assertions.assertTrue(tmpASF.shouldBeFiltered(tmpMolecule));
            Assertions.assertFalse(tmpASF.canBeFragmented(tmpMolecule));
        }
    }

    /**
     * Tests the settings: every setting has a display name and a tooltip, copies keep all settings, the defaults can be
     * restored and invalid values are rejected.
     */
    @Test
    public void settingsTest() {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        Assertions.assertEquals(7, tmpASF.settingsProperties().size());
        tmpASF.settingsProperties().forEach(tmpSetting -> {
            Assertions.assertNotNull(tmpASF.getSettingNameToDisplayNameMap().get(tmpSetting.getName()), tmpSetting.getName());
            Assertions.assertNotNull(tmpASF.getSettingNameToTooltipTextMap().get(tmpSetting.getName()), tmpSetting.getName());
        });
        tmpASF.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        tmpASF.setKeepNonFragmentableMoleculesSetting(false);
        tmpASF.setFragmentSideChainsSetting(false);
        tmpASF.setMaxChainLengthSetting(3);
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpASF.setPreserveStereoChemistrySetting(false);
        AlkylStructureFragmenter tmpCopy = (AlkylStructureFragmenter) tmpASF.copy();
        Assertions.assertEquals(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION, tmpCopy.getFragmentSaturationSetting());
        Assertions.assertFalse(tmpCopy.getKeepNonFragmentableMoleculesSetting());
        Assertions.assertFalse(tmpCopy.getFragmentSideChainsSetting());
        Assertions.assertEquals(3, tmpCopy.getMaxChainLengthSetting());
        Assertions.assertFalse(tmpCopy.getIsolateTertQuatCarbonsSetting());
        Assertions.assertTrue(tmpCopy.getSeparateTertQuatCarbonFromRingSetting());
        Assertions.assertFalse(tmpCopy.getPreserveStereoChemistrySetting());
        tmpASF.restoreDefaultSettings();
        Assertions.assertEquals(IMoleculeFragmenter.FRAGMENT_SATURATION_OPTION_DEFAULT, tmpASF.getFragmentSaturationSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.KEEP_NON_FRAGMENTABLE_MOLECULES_SETTING_DEFAULT, tmpASF.getKeepNonFragmentableMoleculesSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.FRAGMENT_SIDE_CHAINS_SETTING_DEFAULT, tmpASF.getFragmentSideChainsSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.MAX_CHAIN_LENGTH_SETTING_DEFAULT, tmpASF.getMaxChainLengthSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.ISOLATE_TERT_QUAT_CARBONS_SETTING_DEFAULT, tmpASF.getIsolateTertQuatCarbonsSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.SEPARATE_TERT_QUAT_CARBON_FROM_RING_SETTING_DEFAULT, tmpASF.getSeparateTertQuatCarbonFromRingSetting());
        Assertions.assertEquals(AlkylStructureFragmenter.PRESERVE_STEREO_CHEMISTRY_SETTING_DEFAULT, tmpASF.getPreserveStereoChemistrySetting());
        //the copy is independent of the original
        Assertions.assertEquals(3, tmpCopy.getMaxChainLengthSetting());
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpASF.setMaxChainLengthSetting(0));
        Assertions.assertThrows(NullPointerException.class, () -> tmpASF.setFragmentSaturationSetting(null));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Classification">
    /**
     * Test for correct classification of tertiary and quaternary carbon atoms and of the neighbors forming a unit with
     * them if they are not isolated. Atom types are deliberately not perceived, the classification has to work on the
     * bonds of the molecule alone.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void tertQuatCarbonClassificationTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setIsolateTertQuatCarbonsSetting(false);
        //isobutane: atom 1 is tertiary, atoms 0, 2, 3 form a unit with it
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon, 1);
        this.assertAtoms(tmpClassification.isQuaternaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2, 3);
        //neopentyl chain: atom 1 is quaternary, atom 5 is not part of its unit
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)(C)CC"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isQuaternaryCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2, 3, 4);
        //branching atom with a double bond (isobutene) is no tertiary carbon
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)=C"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember);
        //isolated tertiary carbons (default) do not form a unit with their neighbors
        tmpClassification = new AlkylStructureFragmenter().classify(tmpParser.parseSmiles("CC(C)C"));
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember);
    }

    /**
     * Test for correct classification of tertiary/quaternary carbon atoms directly bonded to a ring, which are not
     * isolated unless they are to be separated from rings.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void tertQuatCarbonAtRingClassificationTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        //isopropylcyclopropane: tertiary atom 1 is bonded to ring atom 3 and kept at the ring by default
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C1CC1"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon, 1);
        this.assertAtoms(tmpClassification.isTertQuatCarbonBondedToRing, 1);
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon);
        this.assertAtoms(tmpClassification.isTertQuatUnitMember, 0, 2);
        tmpASF.setSeparateTertQuatCarbonFromRingSetting(true);
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC(C)C1CC1"));
        this.assertAtoms(tmpClassification.isIsolatedTertQuatCarbon, 1);
        //ring atoms themselves are never tertiary carbons
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("CC1CC1"));
        this.assertAtoms(tmpClassification.isTertiaryCarbon);
        this.assertAtoms(tmpClassification.isTertQuatCarbonBondedToRing);
    }

    /**
     * Test for correct classification of ring atoms and bonds.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void ringClassificationTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //methylcyclopropane: methyl atom 0 and its bond are not part of the ring
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1CC1");
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = new AlkylStructureFragmenter().classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isRingAtom, 1, 2, 3);
        this.assertBonds(tmpMolecule, tmpClassification.isRingBond, new int[][]{{1, 2}, {2, 3}, {3, 1}});
    }

    /**
     * Test for correct classification of double and triple bond atoms.
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void multipleBondClassificationTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = new AlkylStructureFragmenter().classify(tmpParser.parseSmiles("CC=CC#C"));
        this.assertAtoms(tmpClassification.hasDoubleBond, 1, 2);
        this.assertAtoms(tmpClassification.hasTripleBond, 3, 4);
        this.assertAtoms(tmpClassification.isChainAtom, 0);
    }

    /**
     * Test for the detection of conjugated pi bond systems.
     * Hexatriene has to be classified completely. Double and triple bonds separated by a single bond (D-S-T) are not
     * regarded as conjugated, also not if the double bond is part of an aromatic ring (phenylacetylene).
     *
     * @throws InvalidSmilesException if SMILES cannot be parsed
     */
    @Test
    public void conjugatedPiSystemClassificationTest() throws InvalidSmilesException {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C=CC=CC=C");
        AlkylStructureFragmenter.MoleculeClassification tmpClassification = tmpASF.classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isConjugatedAtom, 0, 1, 2, 3, 4, 5);
        this.assertBonds(tmpMolecule, tmpClassification.isConjugatedBond, new int[][]{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
        //D-S-T is not conjugated
        tmpClassification = tmpASF.classify(tmpParser.parseSmiles("C=CC#C"));
        this.assertAtoms(tmpClassification.isConjugatedAtom);
        //phenylacetylene: only the ring is conjugated
        tmpMolecule = tmpParser.parseSmiles("C#CC1=CC=CC=C1");
        tmpClassification = tmpASF.classify(tmpMolecule);
        this.assertAtoms(tmpClassification.isConjugatedAtom, 2, 3, 4, 5, 6, 7);
        this.assertBonds(tmpMolecule, tmpClassification.isConjugatedBond, new int[][]{{2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 2}});
    }
    //</editor-fold>
    //
    //<editor-fold desc="Fragmentation rules">
    /**
     * Checks all fragmentation examples given in the class documentation of {@link AlkylStructureFragmenter}, so that
     * the documented rules and the implementation cannot diverge.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void documentationExamplesTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                {"C1CCC(CC1)C1CCCCC1", true, false, 6, "C1CCCCC1 C1CCCCC1"},
                {"C=Cc1ccccc1", true, false, 6, "C=Cc1ccccc1"},
                {"C=CC#C", true, false, 6, "C=C C#C"},
                {"cCCC", true, false, 6, "c CCC"},
                {"C=C1CCCCC1", true, false, 6, "C=C1CCCCC1"},
                {"CC=CC", true, false, 6, "C=C C C"},
                {"CC(C)(C)CCCC(C)C", true, false, 6, "*C(*)(*)* *C(*)* C C C C C CCC"},
                {"c1ccc(cc1)C(C)CC", true, true, 6, "c1ccccc1 *C(*)* C CC"},
                {"CC(C)(C)CCCC(C)C", false, false, 6, "CC(C)(C)C CC(C)C C"},
                {"CC(C)CC(C)C", false, false, 6, "CC(C)CC(C)C"},
                {"CC(C)CCC(C)C", false, false, 6, "CC(C)C CC(C)C"},
                {"CC(C)c1ccccc1", true, false, 6, "CC(C)c1ccccc1"},
                {"CC(C1CCCCC1)C=CC=C", true, false, 6, "C=CC=CC(C)C1CCCCC1"},
                {"C=CC(C)(C)C", false, false, 6, "C=C *C(C)(C)C"},
                {"CC(C)c1ccccc1", false, true, 6, "*C(C)C c1ccccc1"},
                {"CC(C)C(C)C1CCCCC1", true, false, 6, "*C(*)* *C(C)C1CCCCC1 C C"},
                {"CC1CCCCC1", true, false, 6, "C1CCCCC1 C"},
                {"CCCCCCCC", true, false, 3, "CCC CCC CC"}});
    }

    /**
     * Tests rule 1: ring bonds are never cut, fused, bridged and spiro ring systems stay one fragment, rings connected
     * by a non-conjugated single bond are separated (independent of the separation of tertiary/quaternary carbons from
     * rings), ring atoms are never tertiary/quaternary carbons.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void ringTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                {"C1CCCCC1", true, false, 6, "C1CCCCC1"},
                {"C1=CCCCC1", true, false, 6, "C1=CCCCC1"},
                {"c1ccccc1", true, false, 6, "c1ccccc1"},
                {"C1CC2CCCC2C1", true, false, 6, "C1CC2CCCC2C1"},
                {"C1CCC2(CCCCC2)CC1", true, false, 6, "C1CCC2(CCCCC2)CC1"},
                {"C1C2CC3CC1CC(C2)C3", true, false, 6, "C1C2CC3CC1CC(C2)C3"},
                {"C1CC(CCC1)C1CCCCC1", true, true, 6, "C1CCCCC1 C1CCCCC1"},
                {"CC1(C)CCCCC1", true, false, 6, "C1CCCCC1 C C"}});
    }

    /**
     * Tests rules 2 and 3: conjugated pi systems (including aromatic rings connected by a single bond) are never cut,
     * other multiple bonds form isolated fragments, also cumulated double bonds (allenes).
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void multipleBondTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                {"C=CCC", true, false, 6, "C=C CC"},
                {"CC#CC", true, false, 6, "C#C C C"},
                {"CC=C=CC", true, false, 6, "C=C=C C C"},
                {"C=CCC=C", true, false, 6, "C=C C C=C"},
                {"C=CC=CC=C", true, false, 6, "C=CC=CC=C"},
                {"C#CC#C", true, false, 6, "C#CC#C"},
                {"C=CC(C)=C", true, false, 6, "C=CC=C C"},
                {"c1ccc(cc1)-c1ccccc1", true, false, 6, "c1ccc(cc1)-c1ccccc1"}});
    }

    /**
     * Tests rules 4 and 5: isolated tertiary/quaternary carbons, non-isolated units (also of two quaternary carbons
     * sharing a neighbor), carbons kept at or separated from rings, and pseudo atoms for cut bonds of
     * tertiary/quaternary carbons, e.g. to an isolated multiple bond.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void tertQuatCarbonTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                {"CCCC(C)(C)CC(C)(C)CCC", true, false, 6, "*C(*)(*)* *C(*)(*)* C C C C C CCC CCC"},
                {"CCCC(C)(C)CC(C)(C)CCC", false, false, 6, "CC(C)(C)CC(C)(C)C CC CC"},
                {"CC(C)(C)CC1CCC(=C)C2CC21", true, false, 6, "*C(*)(*)* C C C C C=C1CCCC2CC12"},
                {"CC(C)(C)CC1CCC(=C)C2CC21", false, false, 6, "CC(C)(C)C C=C1CCCC2CC12"},
                {"CC(C)(C)C1CCCCC1", true, false, 6, "CC(C)(C)C1CCCCC1"},
                {"CC(C)(C)C1CCCCC1", true, true, 6, "*C(*)(*)* C C C C1CCCCC1"},
                {"CC(C)(C)C1CCCCC1", false, true, 6, "*C(C)(C)C C1CCCCC1"},
                {"CC(C)(C)C(C)(C)C1CCCCC1", true, false, 6, "*C(*)(*)* *C(C)(C)C1CCCCC1 C C C"},
                {"c1ccc(cc1)C(CC)CC", true, true, 6, "c1ccccc1 *C(*)* CC CC"},
                {"CC(C)(CCC)C1CCCCC1", true, true, 6, "C1CCCCC1 *C(*)(*)* C C CCC"},
                {"CC(C1CCCCC1)C(=C)C1=CC=CC1", true, false, 6, "C=C(C1=CC=CC1)C(C)C1CCCCC1"},
                {"CC(C1CCCCC1)C(=C)C1=CC=CC1", false, false, 6, "C=C(C1=CC=CC1)C(C)C1CCCCC1"},
                {"CC(C1CCCCC1)C=CC=C", false, false, 6, "C=CC=CC(C)C1CCCCC1"},
                {"C#CC(C)(C)C", false, false, 6, "C#C *C(C)(C)C"},
                {"C=CC(C)C", false, false, 6, "C=C *C(C)C"},
                {"CC(C)C.C1CCCCC1C1CCCCC1", true, false, 6, "*C(*)* C C C C1CCCCC1 C1CCCCC1"},
                {"CC(C)C.C1CCCCC1C1CCCCC1", false, false, 6, "CC(C)C C1CCCCC1 C1CCCCC1"}});
    }

    /**
     * Tests rule 6: linear chains are dissected into chains of the set maximum length (or not at all), walked from the
     * end bonded to the rest of the molecule.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void linearChainTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                {"CCCCCCCCCCCCCC", true, false, 6, "CCCCCC CCCCCC CC"},
                {"CCCCCCCCCCCCCC", true, false, 7, "CCCCCCC CCCCCCC"},
                {"CCCCCCCCCCCCCC", true, false, 1, "C C C C C C C C C C C C C C"},
                {"CCCCCCCCCCCCCC", true, false, 0, "CCCCCCCCCCCCCC"},
                {"CCCCCCCC1CCCCC1", true, false, 3, "C1CCCCC1 CCC CCC C"},
                {"CC(C)(C)CCCCCCCCC", true, false, 3, "*C(*)(*)* C C C CCC CCC CCC"}});
    }

    /**
     * Tests the transfer of the stereo chemistry to the fragments: complete stereo elements are kept (if the setting is
     * active), stereo elements of atoms that are no stereocenters anymore after the fragmentation are removed.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void stereoChemistryTest() throws Exception {
        //columns: input SMILES, preserve stereo chemistry, expected fragments
        Object[][] tmpCases = {
                {"C=C\\C=C/C=C", true, "C=C/C=C\\C=C"},
                {"C=C\\C=C/C=C", false, "C=CC=CC=C"},
                {"c1ccccc1[C@@H](C1CCC2CCCC12)C1CCCCC1", true, "c1ccc(cc1)[C@H](C2CCCCC2)C3CCC4CCCC43"},
                {"c1ccccc1[C@@H](C1CCC2CCCC12)C1CCCCC1", false, "c1ccc(cc1)C(C2CCCCC2)C3CCC4CCCC43"},
                //steroid skeleton: the two former stereocenters next to the ring bear two methyl groups after the
                // fragmentation and are therefore no stereocenters anymore, only the three ring stereocenters are kept
                {"CCCC[C@H](C)[C@H]1CC[C@]2(C)[C@H](CC[C@@H]2[C@H]1C)[C@H](C)CCCC(C)C", true,
                        "*C(*)* CC(C)[C@H]1CCC2[C@H](CC[C@@H]2C(C)C)C1 CCC C C CC C C"}};
        for (Object[] tmpCase : tmpCases) {
            AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
            tmpASF.setPreserveStereoChemistrySetting((boolean) tmpCase[1]);
            this.assertFragments(tmpASF, this.parse((String) tmpCase[0]), (String) tmpCase[2], true, Arrays.toString(tmpCase));
        }
    }

    /**
     * Tests hydrogen atoms that remain explicit after hydrogen suppression (isotopes, molecular hydrogen). They are no
     * neighbors in the sense of the tertiary/quaternary carbon detection, stay bonded to their atom and do not form
     * fragments of their own; saturation does not create molecular hydrogen.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void remainingExplicitHydrogenTest() throws Exception {
        this.assertExpectedFragments(new Object[][]{
                //secondary carbon with deuterium, no tertiary carbon
                {"[2H]C(C)C", true, false, 6, "[2H]C(C)C"},
                //tertiary carbon with deuterium, the deuterium stays on the carbon
                {"[2H]C(C)(C)C", true, false, 6, "[2H]C(*)(*)* C C C"},
                //deuterium in a chain, the carbon is secondary and therefore cut from the ring
                {"[2H]C(CC)C1CCCCC1", true, false, 6, "[2H]CCC C1CCCCC1"},
                //molecular hydrogen is no fragment of its own after saturation
                {"CC.[H][H]", true, false, 6, "CC [H][H]"}});
    }

    /**
     * Tests the fragmentation of aromatic remnants as produced by a previous fragmentation in a pipeline, e.g. alkane
     * fragments of the ErtlFunctionalGroupsFinderFragmenter whose aromatic ring was opened by the removal of a hetero
     * atom. These cannot be kekulized, so their aromatic bonds have no bond order. They are treated as conjugated pi
     * systems and kept intact, alkyl substituents are fragmented as usual. Cut bonds of the aromatic atoms are
     * saturated with hydrogen, so identical aromatic remnants result in identical fragments.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void aromaticRemnantTest() throws Exception {
        Object[][] tmpCases = {
                {"ccc(c1ccc(cc1)C)CC", true, false, 6, "cccc1ccccc1 C CC"},
                {"CC(C)(C)ccc", true, false, 6, "*C(*)(*)* ccc C C C"},
                {"cc(C)cc(CCCCCCCC)c", true, false, 6, "ccccc C CCCCCC CC"},
                {"cc(CC)ccc", true, false, 6, "ccccc CC"},
                //single aromatic remnant atoms (no aromatic bond) stay conjugated with neighboring aromatic atoms and
                // are neither part of alkyl chains nor tertiary/quaternary carbons
                {"c-c1ccccc1", true, false, 6, "c-c1ccccc1"},
                {"cc-c", true, false, 6, "cc-c"},
                {"cCCC", true, false, 6, "c CCC"},
                {"c(C)(C)C", true, false, 6, "c C C C"}};
        for (Object[] tmpCase : tmpCases) {
            //parsed as in the MORTAR pipeline, kekulization of these remnants fails
            Assertions.assertThrows(CDKException.class, () -> ChemUtil.parseSmilesToAtomContainer((String) tmpCase[0], true, true));
        }
        this.assertExpectedFragments(tmpCases);
        //aromatic remnants are only accepted with aromatic flags, bonds without order and flag are still filtered
        IAtomContainer tmpUnsetBondAC = ChemUtil.parseSmilesToAtomContainer("CCCC", false, false);
        tmpUnsetBondAC.getBond(1).setOrder(IBond.Order.UNSET);
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setKeepNonFragmentableMoleculesSetting(false);
        Assertions.assertTrue(tmpASF.shouldBeFiltered(tmpUnsetBondAC));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Robustness">
    /**
     * Tests that missing implicit hydrogen counts are only set for the atoms without a count, so that given counts (here
     * of a radical as in aromatic remnants of the Ertl fragmenter) are kept.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void missingImplicitHydrogenCountTest() throws Exception {
        SmilesParser tmpParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("[CH2]CCC");
        tmpMolecule.getAtom(3).setImplicitHydrogenCount(null);
        this.assertFragments(new AlkylStructureFragmenter(), tmpMolecule, "[CH2]CCC", false, "[CH2]CCC");
    }

    /**
     * Tests that the dissection of linear chains does not depend on the atom order without hydrogen saturation. The
     * chain is walked from its end that is bonded to the rest of the molecule, so the unsaturated (radical) ends of
     * the chain segments are always at the same positions.
     *
     * @throws Exception if a SMILES cannot be parsed or the fragmentation fails
     */
    @Test
    public void chainDissectionAtomOrderTest() throws Exception {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        tmpASF.setFragmentSaturationSetting(IMoleculeFragmenter.FragmentSaturationOption.NO_SATURATION);
        tmpASF.setMaxChainLengthSetting(3);
        //chain atoms before and after the ring atoms, attachment point at the start, in the middle and at the end
        for (String tmpSmiles : new String[]{"CCCCCC1CCCCC1", "C1CCCCC1CCCCC", "C(C1CCCCC1)CCCC", "C1CCCC(C1)CCCCC"}) {
            this.assertFragments(tmpASF, this.parse(tmpSmiles), "[CH2]C[CH2] [CH2]C [CH]1CCCCC1", false, tmpSmiles);
        }
    }

    /**
     * Test that the returned fragments consist of new atom and bond instances: changing a fragment must not change the
     * input molecule.
     *
     * @throws Exception if SMILES cannot be parsed or fragmentation fails
     */
    @Test
    public void fragmentsAreIndependentCopiesTest() throws Exception {
        AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
        IAtomContainer tmpButeneContainer = this.parse("C=CCC");
        String tmpButeneSMILES = ChemUtil.createUniqueSmiles(tmpButeneContainer, false);
        List<IAtomContainer> tmpFragments = tmpASF.fragmentMolecule(tmpButeneContainer);
        Assertions.assertEquals(this.toSortedUniqueSmiles("C=C CC", false), this.toSortedUniqueSmiles(tmpFragments, false));
        for (IAtomContainer tmpFragment : tmpFragments) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                for (IAtom tmpInputAtom : tmpButeneContainer.atoms()) {
                    Assertions.assertNotSame(tmpInputAtom, tmpAtom);
                }
                tmpAtom.setImplicitHydrogenCount(0);
                tmpAtom.setAtomicNumber(14);
            }
            for (IBond tmpBond : tmpFragment.bonds()) {
                tmpBond.setOrder(IBond.Order.TRIPLE);
            }
        }
        Assertions.assertEquals(tmpButeneSMILES, ChemUtil.createUniqueSmiles(tmpButeneContainer, false));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private utility methods">
    /**
     * Fragments the input of every given case with the given settings and compares the fragments with the expected
     * ones (see class documentation for the columns).
     *
     * @param aCaseArray cases: input SMILES, isolate, separate from ring, maximum chain length (0 for no dissection),
     *                   expected fragments separated by spaces
     * @throws CloneNotSupportedException if the input molecule cannot be cloned
     */
    private void assertExpectedFragments(Object[][] aCaseArray) throws CloneNotSupportedException {
        for (Object[] tmpCase : aCaseArray) {
            AlkylStructureFragmenter tmpASF = new AlkylStructureFragmenter();
            tmpASF.setIsolateTertQuatCarbonsSetting((boolean) tmpCase[1]);
            tmpASF.setSeparateTertQuatCarbonFromRingSetting((boolean) tmpCase[2]);
            int tmpMaxChainLength = (int) tmpCase[3];
            tmpASF.setFragmentSideChainsSetting(tmpMaxChainLength > 0);
            if (tmpMaxChainLength > 0) {
                tmpASF.setMaxChainLengthSetting(tmpMaxChainLength);
            }
            IAtomContainer tmpMolecule = this.parse((String) tmpCase[0]);
            Assertions.assertTrue(tmpASF.canBeFragmented(tmpMolecule), () -> Arrays.toString(tmpCase));
            this.assertFragments(tmpASF, tmpMolecule, (String) tmpCase[4], false, Arrays.toString(tmpCase));
        }
    }

    /**
     * Fragments the given molecule and asserts that the fragments contain all atoms of the molecule and are equal to the
     * expected fragments (compared as sorted unique SMILES).
     *
     * @param anASF configured fragmenter
     * @param aMolecule molecule to fragment
     * @param anExpectedFragments expected fragments as SMILES separated by spaces
     * @param anIsStereoChemEncoded whether the stereo chemistry is compared
     * @param aMessage message describing the case
     * @throws CloneNotSupportedException if the molecule cannot be cloned
     */
    private void assertFragments(AlkylStructureFragmenter anASF, IAtomContainer aMolecule, String anExpectedFragments,
                                 boolean anIsStereoChemEncoded, String aMessage) throws CloneNotSupportedException {
        List<IAtomContainer> tmpFragments = anASF.fragmentMolecule(aMolecule);
        int tmpFragmentAtomCount = 0;
        for (IAtomContainer tmpFragment : tmpFragments) {
            tmpFragmentAtomCount += this.countAtoms(tmpFragment);
        }
        Assertions.assertEquals(this.countAtoms(aMolecule), tmpFragmentAtomCount, () -> aMessage + ": number of atoms");
        Assertions.assertEquals(this.toSortedUniqueSmiles(anExpectedFragments, anIsStereoChemEncoded),
                this.toSortedUniqueSmiles(tmpFragments, anIsStereoChemEncoded), aMessage);
    }

    /**
     * Parses the given SMILES like MORTAR does for fragments in a pipeline (see
     * {@link AlkylTestMoleculeGenerator#parseAsInPipeline(String)}).
     *
     * @param aSmiles SMILES to parse
     * @return the parsed molecule
     */
    private IAtomContainer parse(String aSmiles) {
        IAtomContainer tmpMolecule = AlkylTestMoleculeGenerator.parseAsInPipeline(aSmiles);
        Assertions.assertNotNull(tmpMolecule, () -> "Invalid SMILES " + aSmiles);
        return tmpMolecule;
    }

    /**
     * Returns the sorted unique SMILES of the given molecules.
     *
     * @param aMoleculeList molecules, e.g. fragments
     * @param anIsStereoChemEncoded whether the stereo chemistry is encoded
     * @return sorted unique SMILES
     */
    private List<String> toSortedUniqueSmiles(List<IAtomContainer> aMoleculeList, boolean anIsStereoChemEncoded) {
        List<String> tmpSmilesList = new ArrayList<>(aMoleculeList.size());
        for (IAtomContainer tmpMolecule : aMoleculeList) {
            tmpSmilesList.add(ChemUtil.createUniqueSmiles(tmpMolecule, anIsStereoChemEncoded, true));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }

    /**
     * Returns the sorted unique SMILES of the molecules given as SMILES separated by spaces.
     *
     * @param aSmilesList SMILES separated by spaces
     * @param anIsStereoChemEncoded whether the stereo chemistry is encoded
     * @return sorted unique SMILES
     */
    private List<String> toSortedUniqueSmiles(String aSmilesList, boolean anIsStereoChemEncoded) {
        List<IAtomContainer> tmpMolecules = new ArrayList<>();
        for (String tmpSmiles : aSmilesList.split(" ")) {
            tmpMolecules.add(this.parse(tmpSmiles));
        }
        return this.toSortedUniqueSmiles(tmpMolecules, anIsStereoChemEncoded);
    }

    /**
     * Counts the atoms of the given atom container that are no pseudo atoms.
     *
     * @param anAtomContainer atom container to count the atoms of
     * @return number of non-pseudo atoms
     */
    private int countAtoms(IAtomContainer anAtomContainer) {
        int tmpCount = 0;
        for (IAtom tmpAtom : anAtomContainer.atoms()) {
            if (!(tmpAtom instanceof IPseudoAtom) && tmpAtom.getAtomicNumber() != null && tmpAtom.getAtomicNumber() != 0) {
                tmpCount++;
            }
        }
        return tmpCount;
    }

    /**
     * Asserts that exactly the atoms with the given indices are classified as true in the given classification array.
     *
     * @param aClassificationArray classification array indexed by atom index
     * @param anExpectedAtomIndices indices of the atoms expected to be classified as true
     */
    private void assertAtoms(boolean[] aClassificationArray, int... anExpectedAtomIndices) {
        Set<Integer> tmpExpectedIndices = new HashSet<>();
        for (int tmpIndex : anExpectedAtomIndices) {
            tmpExpectedIndices.add(tmpIndex);
        }
        for (int i = 0; i < aClassificationArray.length; i++) {
            Assertions.assertEquals(tmpExpectedIndices.contains(i), aClassificationArray[i], String.format("Classification of atom %d", i));
        }
    }

    /**
     * Asserts that exactly the bonds between the given atom index pairs are classified as true in the given
     * classification array.
     *
     * @param aMolecule the classified molecule
     * @param aClassificationArray classification array indexed by bond index
     * @param anExpectedBonds atom index pairs of the bonds expected to be classified as true
     */
    private void assertBonds(IAtomContainer aMolecule, boolean[] aClassificationArray, int[][] anExpectedBonds) {
        Set<Integer> tmpExpectedBondIndices = new HashSet<>();
        for (int[] tmpAtomPair : anExpectedBonds) {
            IBond tmpBond = aMolecule.getBond(aMolecule.getAtom(tmpAtomPair[0]), aMolecule.getAtom(tmpAtomPair[1]));
            Assertions.assertNotNull(tmpBond, "No bond between atoms " + tmpAtomPair[0] + " and " + tmpAtomPair[1]);
            tmpExpectedBondIndices.add(aMolecule.indexOf(tmpBond));
        }
        for (int i = 0; i < aClassificationArray.length; i++) {
            IBond tmpBond = aMolecule.getBond(i);
            Assertions.assertEquals(tmpExpectedBondIndices.contains(i), aClassificationArray[i],
                    String.format("Classification of bond %d-%d", aMolecule.indexOf(tmpBond.getBegin()), aMolecule.indexOf(tmpBond.getEnd())));
        }
    }
    //</editor-fold>
}
