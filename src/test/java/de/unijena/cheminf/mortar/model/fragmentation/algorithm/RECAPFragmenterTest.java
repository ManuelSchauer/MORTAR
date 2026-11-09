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

import javafx.beans.property.Property;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Tests the {@link RECAPFragmenter}. The expected fragments were derived by hand from the RECAP rules; the bond type
 * that cleaves a bond is given in the comments. Fragments are compared as canonical SMILES with the bond types encoded
 * as mass numbers of the link atoms, e.g. [5*]Oc1ccc(F)cc1.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class RECAPFragmenterTest {
    /**
     * Constructor that sets the default locale to British English, which is important for the correct functioning of the
     * fragmenter because the settings tooltips are imported from the message.properties file.
     */
    public RECAPFragmenterTest() {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //
    /**
     * Tests the settings, the rule set, copying, and the restoration of the default settings.
     */
    @Test
    public void basicSettingsTest() {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        Assertions.assertEquals(RECAPFragmenter.ALGORITHM_NAME, tmpFragmenter.getFragmentationAlgorithmName());
        Assertions.assertFalse(tmpFragmenter.getFragmentationAlgorithmDisplayName().isBlank());
        Assertions.assertEquals(3, tmpFragmenter.settingsProperties().size());
        for (Property<?> tmpSetting : tmpFragmenter.settingsProperties()) {
            Assertions.assertTrue(tmpFragmenter.getSettingNameToDisplayNameMap().containsKey(tmpSetting.getName()));
            Assertions.assertTrue(tmpFragmenter.getSettingNameToTooltipTextMap().containsKey(tmpSetting.getName()));
        }
        Assertions.assertEquals(RECAPFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(RECAPFragmenter.MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT,
                tmpFragmenter.getMinimumTerminalFragmentSizeSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        tmpFragmenter.setMinimumTerminalFragmentSizeSetting(2);
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        RECAPFragmenter tmpCopy = (RECAPFragmenter) tmpFragmenter.copy();
        Assertions.assertEquals(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION, tmpCopy.getSaturationSetting());
        Assertions.assertEquals(2, tmpCopy.getMinimumTerminalFragmentSizeSetting());
        Assertions.assertFalse(tmpCopy.getReturnUnfragmentedMoleculesSetting());
        tmpFragmenter.restoreDefaultSettings();
        Assertions.assertEquals(RECAPFragmenter.SATURATION_OPTION_DEFAULT, tmpFragmenter.getSaturationSetting());
        Assertions.assertEquals(RECAPFragmenter.MINIMUM_TERMINAL_FRAGMENT_SIZE_DEFAULT,
                tmpFragmenter.getMinimumTerminalFragmentSizeSetting());
        Assertions.assertTrue(tmpFragmenter.getReturnUnfragmentedMoleculesSetting());
        Assertions.assertThrows(NullPointerException.class, () -> tmpFragmenter.setSaturationSetting(null));
        Assertions.assertThrows(IllegalArgumentException.class, () -> tmpFragmenter.setMinimumTerminalFragmentSizeSetting(0));
        List<String> tmpBondTypeNames = RECAPFragmenter.getBondTypeNames();
        Assertions.assertEquals(11, tmpBondTypeNames.size());
        Assertions.assertEquals(11, new HashSet<>(tmpBondTypeNames).size());
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(null));
        Assertions.assertTrue(tmpFragmenter.shouldBeFiltered(SilentChemObjectBuilder.getInstance().newAtomContainer()));
    }
    //
    /**
     * Tests the fragmentation of cisapride, the example of the authors' presentation, which shows the fragments
     * [1]C(=O)-aryl, [1]NH-piperidine-[3]N, [3]CH2-CH2-[5]CH2, and [5]O-fluorophenyl: the amide (1), the amine (3), and the
     * ether (5) at the propyl linker are cleaved, but not the methoxy and amino groups, and the oxygen atom of the aryl
     * ether stays at the aromatic ring.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void cisaprideTest() throws Exception {
        this.assertFragments("COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC",
                "[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]NC1CCN([3*])CC1OC", "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1");
        //the stereocentres of the piperidine ring are kept
        this.assertFragments(new RECAPFragmenter(), "COc1cc(N)c(Cl)cc1C(=O)N[C@@H]1CCN(CCCOc2ccc(F)cc2)C[C@@H]1OC", true,
                "[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]N[C@@H]1CCN([3*])C[C@@H]1OC", "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1");
    }
    //
    /**
     * Tests examples for all eleven bond types.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bondTypesTest() throws Exception {
        //1 amide
        this.assertFragments("CCCCCC(=O)NCCCCC", "[1*]C(=O)CCCCC", "[1*]NCCCCC");
        //2 ester, the oxygen atom stays at the alcohol part
        this.assertFragments("CCCCC(=O)OCCCCC", "[2*]C(=O)CCCC", "[2*]OCCCCC");
        //3 amine, all three bonds of a tertiary amine
        this.assertFragments("CCCCCN(CCCCC)CCCCC", "[3*]N([3*])[3*]", "[3*]CCCCC", "[3*]CCCCC", "[3*]CCCCC");
        //4 urea, both bonds
        this.assertFragments("CCCCCNC(=O)NCCCCC", "[4*]C([4*])=O", "[4*]NCCCCC", "[4*]NCCCCC");
        //5 ether, both bonds of a dialkyl ether
        this.assertFragments("CCCCCOCCCCC", "[5*]O[5*]", "[5*]CCCCC", "[5*]CCCCC");
        //5 ether, aryl alkyl ether only at the sp3 carbon atom
        this.assertFragments("CCCCCOc1ccccc1", "[5*]CCCCC", "[5*]Oc1ccccc1");
        //6 olefin
        this.assertFragments("CCCCCC=CCCCCC", "[6*]=CCCCCC", "[6*]=CCCCCC");
        //7 quaternary nitrogen, all four bonds
        this.assertFragments("CCCCC[N+](CCCCC)(CCCCC)CCCCC", "[7*][N+]([7*])([7*])[7*]", "[7*]CCCCC", "[7*]CCCCC",
                "[7*]CCCCC", "[7*]CCCCC");
        //8 aromatic nitrogen - aliphatic carbon
        this.assertFragments("CCCCCn1ccnc1", "[8*]CCCCC", "[8*]n1ccnc1");
        //9 lactam nitrogen - aliphatic carbon
        this.assertFragments("CCCCCN1CCCC1=O", "[9*]CCCCC", "[9*]N1CCCC1=O");
        //10 aromatic carbon - aromatic carbon
        this.assertFragments("c1ccccc1-c1ccccn1", "[10*]c1ccccc1", "[10*]c1ccccn1");
        //11 sulphonamide
        this.assertFragments("CCCCCNS(=O)(=O)c1ccccc1", "[11*]NCCCCC", "[11*]S(=O)(=O)c1ccccc1");
        //carbamate: amide (1) and ester (2) bonds, the carbonyl group links both parts
        this.assertFragments("CC(C)(C)OC(=O)NCCCCC", "[1*]C([2*])=O", "[1*]NCCCCC", "[2*]OC(C)(C)C");
    }
    //
    /**
     * Tests that ring bonds and the bonds of small terminal groups are not cleaved, and that no bond outside the
     * eleven types is cleaved.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void notCleavedTest() throws Exception {
        String[] tmpUnfragmentedMolecules = {
                //small terminal groups (fewer than five heavy atoms): methyl, ethyl, propyl, butyl, methoxy, acetyl,
                //dimethylamino, amino, hydroxy, trimethylammonium
                "Cc1ccccc1", "CCOc1ccccc1", "COc1ccccc1", "CCCCN(CCCC)CCCC", "CC(=O)NCCCCC", "CN(C)Cc1ccccc1",
                "NCc1ccccc1", "OCC1CCCCC1", "CCCCCC[N+](C)(C)C", "NS(=O)(=O)c1ccccc1", "OC(=O)c1ccccc1",
                //no RECAP bond types: aryl ethers and amines, ketones, thioethers, alkyl-aryl bonds, ring bonds
                "c1ccccc1Oc1ccccc1", "c1ccccc1N1CCCCC1", "O=C(c1ccccc1)c1ccccc1", "CCCCCSCCCCC", "CCCCCc1ccccc1",
                "O=C1CCCCN1", "C1CCOC1", "c1ccc2ccccc2c1"};
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        for (String tmpSmiles : tmpUnfragmentedMolecules) {
            List<IAtomContainer> tmpFragments = tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(tmpSmiles));
            Assertions.assertEquals(1, tmpFragments.size(), tmpSmiles);
            Assertions.assertEquals(RECAPFragmenterTest.canonical(tmpSmiles, false),
                    RECAPFragmenterTest.toCanonicalSmiles(tmpFragments.getFirst(), false), tmpSmiles);
        }
    }
    //
    /**
     * Tests the setting for the minimum size of terminal fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void minimumTerminalFragmentSizeTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        tmpFragmenter.setMinimumTerminalFragmentSizeSetting(4);
        this.assertFragments(tmpFragmenter, "CCCCN(CCCC)CCCC", false, "[3*]N([3*])[3*]", "[3*]CCCC", "[3*]CCCC", "[3*]CCCC");
        //the propyl groups are too small, but the pentyl group can be cleaved off the dipropylamino part
        this.assertFragments(tmpFragmenter, "CCCCCN(CCC)CCC", false, "[3*]CCCCC", "[3*]N(CCC)CCC");
        tmpFragmenter.setMinimumTerminalFragmentSizeSetting(1);
        this.assertFragments(tmpFragmenter, "COc1ccccc1", false, "[5*]C", "[5*]Oc1ccccc1");
        //a linking fragment is allowed independent of its size, here the single oxygen atom
        tmpFragmenter.setMinimumTerminalFragmentSizeSetting(8);
        this.assertFragments(tmpFragmenter, "CCCCCCCCOCCCCCCCC", false, "[5*]O[5*]", "[5*]CCCCCCCC", "[5*]CCCCCCCC");
        tmpFragmenter.setMinimumTerminalFragmentSizeSetting(9);
        this.assertFragments(tmpFragmenter, "CCCCCCCCOCCCCCCCC", false, "CCCCCCCCOCCCCCCCC");
    }
    //
    /**
     * Tests the saturation of the open valences with hydrogen atoms.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void hydrogenSaturationTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        this.assertFragments(tmpFragmenter, "CCCCCOCCCCC", false, "O", "CCCCC", "CCCCC");
        //cleaved double bonds are saturated with two hydrogen atoms
        this.assertFragments(tmpFragmenter, "CCCCCC=CCCCCC", false, "CCCCCC", "CCCCCC");
        this.assertFragments(tmpFragmenter, "COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC", false,
                "O=Cc1cc(Cl)c(N)cc1OC", "NC1CCNCC1OC", "CCC", "Oc1ccc(F)cc1");
        for (IAtomContainer tmpFragment : tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("CCCCCOCCCCC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                Assertions.assertFalse(tmpAtom instanceof IPseudoAtom);
                Assertions.assertNotEquals(1, tmpAtom.getAtomicNumber());
            }
        }
    }
    //
    /**
     * Tests the link atoms: pseudo atoms with the bond type as mass number and atom property.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void linkAtomTest() throws Exception {
        Set<Integer> tmpBondTypes = new HashSet<>();
        for (IAtomContainer tmpFragment : new RECAPFragmenter().fragmentMolecule(
                RECAPFragmenterTest.parse("COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC"))) {
            for (IAtom tmpAtom : tmpFragment.atoms()) {
                if (tmpAtom instanceof IPseudoAtom) {
                    Integer tmpBondType = tmpAtom.getProperty(RECAPFragmenter.BOND_TYPE_PROPERTY_KEY);
                    Assertions.assertEquals(tmpBondType, tmpAtom.getMassNumber());
                    Assertions.assertEquals(0, tmpAtom.getImplicitHydrogenCount());
                    Assertions.assertEquals(1, tmpFragment.getConnectedBondsCount(tmpAtom));
                    tmpBondTypes.add(tmpBondType);
                }
            }
        }
        Assertions.assertEquals(Set.of(1, 3, 5), tmpBondTypes);
    }
    //
    /**
     * Tests the setting for molecules without cleavable bonds and the handling of multi-component molecules,
     * explicit hydrogen atoms, and pseudo atoms of preceding fragmentations in a pipeline.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void inputHandlingTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        Assertions.assertEquals(1, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("Cc1ccccc1")).size());
        tmpFragmenter.setReturnUnfragmentedMoleculesSetting(false);
        Assertions.assertTrue(tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("Cc1ccccc1")).isEmpty());
        Assertions.assertEquals(3, tmpFragmenter.fragmentMolecule(RECAPFragmenterTest.parse("CCCCCOCCCCC")).size());
        //counter ions are returned separately, also for molecules without cleaved bonds
        this.assertFragments("CCCCCOCCCCC.[Na+].[Cl-]", "[5*]O[5*]", "[5*]CCCCC", "[5*]CCCCC", "[Na+]", "[Cl-]");
        this.assertFragments("Cc1ccccc1.[Na+]", "Cc1ccccc1", "[Na+]");
        //explicit hydrogen atoms do not change the result
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse("COc1cc(N)c(Cl)cc1C(=O)NC1CCN(CCCOc2ccc(F)cc2)CC1OC");
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        Assertions.assertEquals(
                RECAPFragmenterTest.canonical(List.of("[1*]C(=O)c1cc(Cl)c(N)cc1OC", "[1*]NC1CCN([3*])CC1OC",
                        "[3*]CCC[5*]", "[5*]Oc1ccc(F)cc1"), false),
                RECAPFragmenterTest.toCanonicalSmiles(new RECAPFragmenter().fragmentMolecule(tmpMolecule), false));
        //pseudo atoms of a preceding fragmentation count as attachment points, so *CCC is a linking fragment
        this.assertFragments("*CCCOCCCCC", "[5*]O[5*]", "[5*]CCC*", "[5*]CCCCC");
        //the input molecule is not altered
        IAtomContainer tmpInput = RECAPFragmenterTest.parse("CCCCCOCCCCC");
        new RECAPFragmenter().fragmentMolecule(tmpInput);
        Assertions.assertEquals(11, tmpInput.getAtomCount());
        Assertions.assertEquals(10, tmpInput.getBondCount());
    }
    //
    /**
     * Tests the transfer of tetrahedral and double bond stereo elements to the fragments.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void stereoTest() throws Exception {
        RECAPFragmenter tmpFragmenter = new RECAPFragmenter();
        //the ether oxygen atom of the stereocentre is replaced by the link atom at the same position
        this.assertFragments(tmpFragmenter, "CCCCC[C@@H](C)OCCCCC", true, "CCCCC[C@@H](C)[5*]", "[5*]O[5*]", "[5*]CCCCC");
        this.assertFragments(tmpFragmenter, "CCCCC[C@H](C)OCCCCC", true, "CCCCC[C@H](C)[5*]", "[5*]O[5*]", "[5*]CCCCC");
        //two link atoms at a former stereocentre, the stereo element is removed
        this.assertFragments(tmpFragmenter, "CCCCCN[C@@H](CCCCC)OCCCCC", true, "[3*]C([5*])CCCCC", "[5*]O[5*]",
                "[3*]N[3*]", "[5*]CCCCC", "[3*]CCCCC");
        //the lactam nitrogen atom is a ligand of the stereo double bond (9); the double bond itself is not cleaved (6),
        //because the methyl group is too small
        this.assertFragments(tmpFragmenter, "O=C1CCCN1/C(CCCCC)=C/C", true, "[9*]/C(CCCCC)=C/C", "[9*]N1CCCC1=O");
        //cleaved stereo double bond
        this.assertFragments(tmpFragmenter, "CCCCC/C=C/CCCCC", true, "[6*]=CCCCCC", "[6*]=CCCCCC");
        tmpFragmenter.setSaturationSetting(RECAPFragmenter.SaturationOption.HYDROGEN_SATURATION);
        //with hydrogen saturation, the former stereocentre has two hydrogen atoms
        this.assertFragments(tmpFragmenter, "CCCCC[C@@H](C)OCCCCC", true, "CCCCCCC", "O", "CCCCC");
        //with hydrogen saturation, the configuration refers to the remaining substituent
        this.assertFragments(tmpFragmenter, "O=C1CCCN1/C(CCCCC)=C/C", true, "CCCCC/C=C\\C", "O=C1CCCN1");
    }
    //
    /**
     * Fragments the given molecule with a fragmenter with default settings and compares the fragments to the given
     * expected fragments, without stereochemistry.
     *
     * @param aSmiles SMILES of the molecule
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(String aSmiles, String... anExpectedFragments) throws Exception {
        this.assertFragments(new RECAPFragmenter(), aSmiles, false, anExpectedFragments);
    }
    //
    /**
     * Fragments the given molecule and compares the fragments to the given expected fragments, once with bond types
     * and without stereochemistry, and optionally with stereochemistry as MORTAR encodes it.
     *
     * @param aFragmenter fragmenter
     * @param aSmiles SMILES of the molecule
     * @param isStereo whether the stereochemistry is compared
     * @param anExpectedFragments SMILES of the expected fragments
     * @throws Exception if anything goes wrong
     */
    private void assertFragments(RECAPFragmenter aFragmenter, String aSmiles, boolean isStereo,
                                 String... anExpectedFragments) throws Exception {
        List<IAtomContainer> tmpFragments = aFragmenter.fragmentMolecule(RECAPFragmenterTest.parse(aSmiles));
        Assertions.assertEquals(RECAPFragmenterTest.canonical(List.of(anExpectedFragments), false),
                RECAPFragmenterTest.toCanonicalSmiles(tmpFragments, false), aSmiles);
        if (isStereo) {
            Assertions.assertEquals(RECAPFragmenterTest.toStereoSmiles(List.of(anExpectedFragments), false),
                    RECAPFragmenterTest.toStereoSmilesOfFragments(tmpFragments, false), aSmiles);
        }
    }
    //
    /**
     * Parses the given SMILES code with kekulization, as MORTAR does.
     *
     * @param aSmiles SMILES code
     * @return molecule
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static IAtomContainer parse(String aSmiles) throws CDKException {
        return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(aSmiles);
    }
    //
    /**
     * Returns the sorted canonical SMILES codes of the given fragments, see {@link #canonical(String, boolean)}.
     *
     * @param aFragmentList fragments
     * @param isStereo whether stereochemistry is encoded
     * @return sorted list of canonical SMILES codes
     * @throws CDKException if a SMILES code cannot be created
     */
    static List<String> toCanonicalSmiles(List<IAtomContainer> aFragmentList, boolean isStereo) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(RECAPFragmenterTest.toCanonicalSmiles(tmpFragment, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the canonical SMILES code of the given fragment, see {@link #canonical(String, boolean)}.
     *
     * @param aFragment fragment
     * @param isStereo whether stereochemistry is encoded
     * @return canonical SMILES code
     * @throws CDKException if a SMILES code cannot be created
     */
    static String toCanonicalSmiles(IAtomContainer aFragment, boolean isStereo) throws CDKException {
        if (isStereo) {
            return RECAPFragmenterTest.toStereoSmiles(aFragment, false);
        }
        String tmpSmiles = SmilesGenerator.create(aFragment, SmiFlavor.Unique | SmiFlavor.AtomicMass,
                new int[aFragment.getAtomCount()]);
        return RECAPFragmenterTest.canonical(tmpSmiles, false);
    }
    //
    /**
     * Returns the unique SMILES code with stereochemistry of the given fragment as MORTAR creates it, i.e. with all link
     * atoms written as '*' (the bond types are not encoded). The CDK creates canonical SMILES codes with
     * stereochemistry based on InChI numbers, which do not support pseudo atoms with mass numbers.
     *
     * @param aFragment fragment, it is not altered
     * @param isNonStereogenicElementRemoved whether stereo elements of atoms that are not stereogenic according to the
     *                                       CDK are removed first
     * @return unique SMILES code with stereochemistry
     * @throws CDKException if the fragment cannot be cloned
     */
    @SuppressWarnings("rawtypes")
    static String toStereoSmiles(IAtomContainer aFragment, boolean isNonStereogenicElementRemoved) throws CDKException {
        IAtomContainer tmpFragment;
        try {
            tmpFragment = aFragment.clone();
        } catch (CloneNotSupportedException anException) {
            throw new CDKException(anException.toString(), anException);
        }
        for (IAtom tmpAtom : tmpFragment.atoms()) {
            if (tmpAtom instanceof IPseudoAtom) {
                tmpAtom.setMassNumber(null);
            }
        }
        //the aromaticity is perceived on the fragment for a representation independent of the parent molecule
        SmartsPattern.prepare(tmpFragment);
        if (isNonStereogenicElementRemoved && tmpFragment.stereoElements().iterator().hasNext()) {
            Stereocenters tmpStereocenters = Stereocenters.of(tmpFragment);
            List<IStereoElement> tmpElements = new ArrayList<>();
            for (IStereoElement tmpElement : tmpFragment.stereoElements()) {
                if (tmpElement.getConfigClass() == IStereoElement.TH
                        && !tmpStereocenters.isStereocenter(tmpFragment.indexOf((IAtom) tmpElement.getFocus()))) {
                    continue;
                }
                if (tmpElement.getConfigClass() == IStereoElement.CT
                        && (!tmpStereocenters.isStereocenter(tmpFragment.indexOf(((IBond) tmpElement.getFocus()).getBegin()))
                        || !tmpStereocenters.isStereocenter(tmpFragment.indexOf(((IBond) tmpElement.getFocus()).getEnd())))) {
                    continue;
                }
                tmpElements.add(tmpElement);
            }
            tmpFragment.setStereoElements(tmpElements);
        }
        return ChemUtil.createUniqueSmiles(tmpFragment, true);
    }
    //
    /**
     * Returns the sorted unique SMILES codes with stereochemistry of the given SMILES codes, see
     * {@link #toStereoSmiles(IAtomContainer, boolean)}.
     *
     * @param aSmilesList SMILES codes
     * @param isNonStereogenicElementRemoved whether stereo elements of non-stereogenic atoms are removed first
     * @return sorted list of unique SMILES codes with stereochemistry
     * @throws CDKException if a SMILES code cannot be parsed
     */
    static List<String> toStereoSmiles(List<String> aSmilesList, boolean isNonStereogenicElementRemoved) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesList.size());
        for (String tmpSmiles : aSmilesList) {
            tmpSmilesList.add(RECAPFragmenterTest.toStereoSmiles(RECAPFragmenterTest.parse(tmpSmiles),
                    isNonStereogenicElementRemoved));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the sorted unique SMILES codes with stereochemistry of the given fragments, see
     * {@link #toStereoSmiles(IAtomContainer, boolean)}.
     *
     * @param aFragmentList fragments
     * @param isNonStereogenicElementRemoved whether stereo elements of non-stereogenic atoms are removed first
     * @return sorted list of unique SMILES codes with stereochemistry
     * @throws CDKException if a fragment cannot be cloned
     */
    static List<String> toStereoSmilesOfFragments(List<IAtomContainer> aFragmentList, boolean isNonStereogenicElementRemoved)
            throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aFragmentList.size());
        for (IAtomContainer tmpFragment : aFragmentList) {
            tmpSmilesList.add(RECAPFragmenterTest.toStereoSmiles(tmpFragment, isNonStereogenicElementRemoved));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns the sorted canonical SMILES codes of the given SMILES codes, see {@link #canonical(String, boolean)}.
     *
     * @param aSmilesList SMILES codes
     * @param isStereo whether stereochemistry is encoded
     * @return sorted list of canonical SMILES codes
     * @throws CDKException if a SMILES code cannot be parsed
     */
    static List<String> canonical(List<String> aSmilesList, boolean isStereo) throws CDKException {
        List<String> tmpSmilesList = new ArrayList<>(aSmilesList.size());
        for (String tmpSmiles : aSmilesList) {
            tmpSmilesList.add(RECAPFragmenterTest.canonical(tmpSmiles, isStereo));
        }
        Collections.sort(tmpSmilesList);
        return tmpSmilesList;
    }
    //
    /**
     * Returns a canonical SMILES code of the given SMILES code that is independent of the aromaticity notation of the
     * input: the SMILES is parsed and kekulized, aromaticity is perceived with the Daylight model (as by the
     * fragmenter), and a unique SMILES with mass numbers (i.e. bond types) and aromatic symbols is created.
     *
     * @param aSmiles SMILES code
     * @param isStereo whether stereochemistry is encoded (as MORTAR encodes it, without bond types)
     * @return canonical SMILES code
     * @throws CDKException if the SMILES code cannot be parsed
     */
    static String canonical(String aSmiles, boolean isStereo) throws CDKException {
        IAtomContainer tmpMolecule = RECAPFragmenterTest.parse(aSmiles);
        SmartsPattern.prepare(tmpMolecule);
        if (isStereo) {
            return RECAPFragmenterTest.toStereoSmiles(tmpMolecule, false);
        }
        return SmilesGenerator.create(tmpMolecule, SmiFlavor.Unique | SmiFlavor.AtomicMass | SmiFlavor.UseAromaticSymbols,
                new int[tmpMolecule.getAtomCount()]);
    }
}
