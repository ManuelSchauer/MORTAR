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

import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Stateless utility methods that implement the structural molecule filters (connectivity, largest component, elements,
 * heavy atom count) on CDK atom containers and create the key of the duplicate filter. The methods neither call GUI
 * code nor depend on settings; the order of application and the counting are done by {@link MoleculeFilterRun}.
 * <br>Definitions used throughout: the atomic number of a pseudo atom ({@link IPseudoAtom}, e.g. an R group or a '*'
 * attachment point) is 0. A heavy atom is an atom with an atomic number larger than 1, i.e. hydrogen atoms (including
 * deuterium and tritium, explicit or implicit) and pseudo atoms are not heavy atoms.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class MoleculeFilterUtil {
    //<editor-fold desc="Public static final constants">
    /**
     * Unmodifiable set of the atomic numbers of the elements allowed by the element filter: H, B, C, N, O, F, Si, P, S,
     * Cl, Se, Br, I. Isotopes of these elements are allowed as well since only the atomic number is checked.
     */
    public static final Set<Integer> ALLOWED_ATOMIC_NUMBERS = Set.of(
            IElement.H, IElement.B, IElement.C, IElement.N, IElement.O, IElement.F, IElement.Si, IElement.P,
            IElement.S, IElement.Cl, IElement.Se, IElement.Br, IElement.I);
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(MoleculeFilterUtil.class.getName());

    /**
     * SMILES flavour of the first tie-break key for the largest component: canonical Kekule SMILES without
     * stereochemistry. Without aromatic symbols, CDK normalises the Kekule structure of canonical SMILES, so the key
     * does not depend on the aromatic or Kekule notation of the input (with aromatic symbols, a molecule without
     * aromaticity flags would be written with the bond orders of the input).
     */
    private static final int TIE_BREAK_CONSTITUTION_SMILES_FLAVOUR = SmiFlavor.Canonical;

    /**
     * SMILES flavour of the second tie-break key for the largest component: like the first one plus stereochemistry.
     */
    private static final int TIE_BREAK_STEREO_SMILES_FLAVOUR = MoleculeFilterUtil.TIE_BREAK_CONSTITUTION_SMILES_FLAVOUR
            | SmiFlavor.Stereo;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private MoleculeFilterUtil() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Returns the atomic number of the given atom. Pseudo atoms ({@link IPseudoAtom}) always have the atomic number 0,
     * and so have atoms with an unset atomic number (in CDK 2.13, the element symbol is derived from the atomic number,
     * so such atoms have no symbol either), i.e. they are treated as unknown elements.
     *
     * @param anAtom the atom
     * @return the atomic number, 0 for pseudo atoms and unknown elements
     * @throws NullPointerException if the given atom is null
     */
    public static int getAtomicNumber(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is null.");
        if (anAtom instanceof IPseudoAtom) {
            return 0;
        }
        Integer tmpAtomicNumber = anAtom.getAtomicNumber();
        return Objects.isNull(tmpAtomicNumber) ? 0 : tmpAtomicNumber;
    }

    /**
     * Returns whether the given atom is a heavy atom, i.e. its atomic number is larger than 1. Hydrogen atoms of all
     * isotopes and pseudo atoms are no heavy atoms.
     *
     * @param anAtom the atom
     * @return true if the atom is a heavy atom
     * @throws NullPointerException if the given atom is null
     */
    public static boolean isHeavyAtom(IAtom anAtom) throws NullPointerException {
        return MoleculeFilterUtil.getAtomicNumber(anAtom) > IElement.H;
    }

    /**
     * Counts the heavy atoms of the given molecule (see {@link #isHeavyAtom(IAtom)}). Implicit hydrogen counts are not
     * considered and explicit hydrogen atoms are not counted.
     *
     * @param aMolecule the molecule
     * @return number of heavy atoms, 0 for an empty molecule
     * @throws NullPointerException if the given molecule is null
     */
    public static int countHeavyAtoms(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        int tmpCount = 0;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (MoleculeFilterUtil.isHeavyAtom(tmpAtom)) {
                tmpCount++;
            }
        }
        return tmpCount;
    }

    /**
     * Returns whether the given atom is of an element allowed by the element filter
     * ({@link #ALLOWED_ATOMIC_NUMBERS}). Pseudo atoms are never allowed, isotopes of allowed elements are.
     *
     * @param anAtom the atom
     * @return true if the atom is of an allowed element
     * @throws NullPointerException if the given atom is null
     */
    public static boolean isAllowedElement(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is null.");
        if (anAtom instanceof IPseudoAtom) {
            return false;
        }
        return MoleculeFilterUtil.ALLOWED_ATOMIC_NUMBERS.contains(MoleculeFilterUtil.getAtomicNumber(anAtom));
    }

    /**
     * Returns whether all atoms of the given molecule are of an element allowed by the element filter (see
     * {@link #isAllowedElement(IAtom)}). Implicit hydrogen atoms are always allowed. An empty molecule consists of
     * allowed elements only.
     *
     * @param aMolecule the molecule
     * @return true if every atom is of an allowed element
     * @throws NullPointerException if the given molecule is null
     */
    public static boolean consistsOfAllowedElementsOnly(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (!MoleculeFilterUtil.isAllowedElement(tmpAtom)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns whether the given molecule is connected, i.e. consists of a single component. Molecules with fewer than
     * two atoms (including empty ones) are connected by definition (as in CDK's ConnectivityChecker). Every unconnected
     * atom counts as an own component, also an isolated hydrogen atom or pseudo atom.
     *
     * @param aMolecule the molecule
     * @return true if the molecule is connected
     * @throws NullPointerException if the given molecule is null
     */
    public static boolean isConnected(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return ConnectivityChecker.isConnected(aMolecule);
    }

    /**
     * Returns the largest connected component of the given molecule. If the molecule is connected (or has fewer than
     * two atoms), the given container itself is returned unchanged. Otherwise, a new container is returned.
     * <p>
     * Largest means the most heavy atoms (see {@link #isHeavyAtom(IAtom)}); this is a size criterion only, so if a
     * counter-ion has more heavy atoms than the molecule of interest (e.g. tosylate or pamoate salts of small
     * compounds), the counter-ion is returned. Ties are broken deterministically and independently of the order of the
     * components in the input: (1) the lexicographically smallest canonical SMILES without stereochemistry, (2) the
     * lexicographically smallest canonical SMILES with stereochemistry, (3) the component that comes first in the atom
     * order of the input (only reached if both SMILES are equal, i.e. for identical components, or cannot be
     * generated, e.g. for components whose Kekule structure could not be assigned in the preprocessing; a component
     * with a SMILES string wins over one without). The SMILES strings are canonical Kekule SMILES, which CDK normalises,
     * so the tie-break does not depend on the aromatic or Kekule notation of the input either.
     * <p>
     * The components are created by {@link ConnectivityChecker#partitionIntoMolecules(IAtomContainer)}, which in CDK
     * 2.13 (tested in MoleculeFilterUtilTest) keeps the atoms (the same atom objects as in the input container, so
     * atom properties, flags, charges, isotopes, implicit hydrogen counts, and coordinates are kept), the bonds, single
     * electrons, lone pairs, the stereo elements whose atoms all belong to the component (re-wrapped for the new
     * container), and the Sgroups that lie completely in the component, but drops all container-level data. This
     * method therefore copies the container properties (including the title, the MORTAR molecule name, and all SDF data
     * fields; except the Sgroup property, which CDK has already split), the container ID, and the container flags from
     * the input to the returned component. Because the atom and bond objects are shared, the input container should
     * be discarded and not be changed afterwards.
     *
     * @param aMolecule the molecule
     * @return the given molecule if it is connected, otherwise a new container holding its largest component
     * @throws NullPointerException if the given molecule is null
     */
    public static IAtomContainer extractLargestComponent(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (MoleculeFilterUtil.isConnected(aMolecule)) {
            return aMolecule;
        }
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        if (tmpComponents.getAtomContainerCount() < 2) {
            //not expected since the molecule is not connected
            return aMolecule;
        }
        //the components are ordered by their first atom in the input
        List<IAtomContainer> tmpLargestComponents = new ArrayList<>(tmpComponents.getAtomContainerCount());
        int tmpMaximumHeavyAtomCount = -1;
        for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
            int tmpHeavyAtomCount = MoleculeFilterUtil.countHeavyAtoms(tmpComponent);
            if (tmpHeavyAtomCount > tmpMaximumHeavyAtomCount) {
                tmpMaximumHeavyAtomCount = tmpHeavyAtomCount;
                tmpLargestComponents.clear();
                tmpLargestComponents.add(tmpComponent);
            } else if (tmpHeavyAtomCount == tmpMaximumHeavyAtomCount) {
                tmpLargestComponents.add(tmpComponent);
            }
        }
        IAtomContainer tmpSelectedComponent;
        if (tmpLargestComponents.size() == 1) {
            tmpSelectedComponent = tmpLargestComponents.getFirst();
        } else {
            tmpSelectedComponent = MoleculeFilterUtil.breakLargestComponentTie(tmpLargestComponents);
        }
        MoleculeFilterUtil.copyContainerLevelData(aMolecule, tmpSelectedComponent);
        return tmpSelectedComponent;
    }

    /**
     * Returns the key under which the duplicate filter compares the given molecule: MORTAR's unique SMILES string
     * without aromaticity encoding, i.e. {@link ChemUtil#createUniqueSmiles(IAtomContainer, boolean, boolean)} with
     * aromaticity encoding off, with stereochemistry if requested and present. In this flavour, CDK normalises the
     * Kekule structure, so the key does not depend on the atom order or on the aromatic or Kekule notation of the input
     * (e.g. Cc1ccccc1, CC1=CC=CC=C1, and CC1C=CC=CC=1 give the same key). The unique SMILES string with aromatic
     * symbols, which the importer stores by default, does not have this property, because a molecule without
     * aromaticity flags (e.g. from an SD file) is written with the bond orders of the input.
     * <br>The given stored SMILES string (the one the importer generated for the data model) is returned as key without
     * generating a new one if it was generated without aromaticity encoding (import of aromatics as Kekule structures),
     * and it is the fallback if the key cannot be generated (e.g. if the Kekule structure of the molecule could not be
     * assigned in the preprocessing).
     *
     * @param aMolecule the molecule, as passed to the SMILES generation of the importer
     * @param isRegardStereo whether stereochemistry is encoded, as in the stored SMILES string
     * @param aStoredUniqueSmiles the unique SMILES string the importer generated for the molecule
     * @param isStoredSmilesAromaticityEncoded whether the stored SMILES string was generated with aromaticity encoding
     * @return the duplicate key, never null or blank
     * @throws NullPointerException if the given molecule or SMILES string is null
     * @throws IllegalArgumentException if the given SMILES string is blank
     */
    public static String createDuplicateKey(IAtomContainer aMolecule,
                                            boolean isRegardStereo,
                                            String aStoredUniqueSmiles,
                                            boolean isStoredSmilesAromaticityEncoded)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        Objects.requireNonNull(aStoredUniqueSmiles, "Given SMILES string is null.");
        if (aStoredUniqueSmiles.isBlank()) {
            throw new IllegalArgumentException("Given SMILES string is blank.");
        }
        if (!isStoredSmilesAromaticityEncoded) {
            return aStoredUniqueSmiles;
        }
        //returns null and logs the reason if no SMILES string could be created
        String tmpKey = ChemUtil.createUniqueSmiles(aMolecule, isRegardStereo, false);
        if (Objects.isNull(tmpKey) || tmpKey.isBlank()) {
            return aStoredUniqueSmiles;
        }
        return tmpKey;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static methods">
    /**
     * Selects one of several components with the same heavy atom count, see
     * {@link #extractLargestComponent(IAtomContainer)} for the criteria.
     *
     * @param aComponentList the tied components, in the atom order of the input
     * @return the selected component
     */
    private static IAtomContainer breakLargestComponentTie(List<IAtomContainer> aComponentList) {
        IAtomContainer tmpSelectedComponent = null;
        String tmpSelectedConstitutionSmiles = null;
        String tmpSelectedStereoSmiles = null;
        for (IAtomContainer tmpComponent : aComponentList) {
            String tmpConstitutionSmiles = MoleculeFilterUtil.createSmilesOrNull(tmpComponent,
                    MoleculeFilterUtil.TIE_BREAK_CONSTITUTION_SMILES_FLAVOUR);
            String tmpStereoSmiles = MoleculeFilterUtil.createSmilesOrNull(tmpComponent,
                    MoleculeFilterUtil.TIE_BREAK_STEREO_SMILES_FLAVOUR);
            if (Objects.isNull(tmpSelectedComponent)) {
                tmpSelectedComponent = tmpComponent;
                tmpSelectedConstitutionSmiles = tmpConstitutionSmiles;
                tmpSelectedStereoSmiles = tmpStereoSmiles;
                continue;
            }
            int tmpComparison = MoleculeFilterUtil.compareNullLast(tmpConstitutionSmiles, tmpSelectedConstitutionSmiles);
            if (tmpComparison == 0) {
                tmpComparison = MoleculeFilterUtil.compareNullLast(tmpStereoSmiles, tmpSelectedStereoSmiles);
            }
            //only a strictly smaller key replaces the selection, so the first component in atom order wins ties
            if (tmpComparison < 0) {
                tmpSelectedComponent = tmpComponent;
                tmpSelectedConstitutionSmiles = tmpConstitutionSmiles;
                tmpSelectedStereoSmiles = tmpStereoSmiles;
            }
        }
        return tmpSelectedComponent;
    }

    /**
     * Compares two strings lexicographically with null being larger than every string and equal to null.
     *
     * @param aString first string, may be null
     * @param anotherString second string, may be null
     * @return negative, zero, or positive value as the first string is smaller than, equal to, or larger than the second
     */
    private static int compareNullLast(String aString, String anotherString) {
        if (Objects.isNull(aString)) {
            return Objects.isNull(anotherString) ? 0 : 1;
        }
        if (Objects.isNull(anotherString)) {
            return -1;
        }
        return aString.compareTo(anotherString);
    }

    /**
     * Creates a SMILES string of the given molecule with the given flavour or returns null if this fails. Only used
     * for the tie-break of the largest component.
     *
     * @param aMolecule the molecule
     * @param aFlavour the SMILES flavour, see {@link SmiFlavor}
     * @return the SMILES string or null
     */
    private static String createSmilesOrNull(IAtomContainer aMolecule, int aFlavour) {
        try {
            return new SmilesGenerator(aFlavour).create(aMolecule);
        } catch (CDKException | RuntimeException anException) {
            MoleculeFilterUtil.LOGGER.log(Level.FINE, String.format("SMILES generation for the tie-break of the largest " +
                    "component failed: %s", anException.toString()), anException);
            return null;
        }
    }

    /**
     * Copies the container-level data that {@link ConnectivityChecker#partitionIntoMolecules(IAtomContainer)} drops
     * from the original molecule to one of its components: all properties except the Sgroups (which CDK has already
     * split between the components and which refer to atoms), the ID, and the flags.
     *
     * @param anOriginal the molecule the component was created from
     * @param aComponent the component to copy the data to
     */
    private static void copyContainerLevelData(IAtomContainer anOriginal, IAtomContainer aComponent) {
        for (Map.Entry<Object, Object> tmpEntry : anOriginal.getProperties().entrySet()) {
            if (CDKConstants.CTAB_SGROUPS.equals(tmpEntry.getKey())) {
                continue;
            }
            aComponent.setProperty(tmpEntry.getKey(), tmpEntry.getValue());
        }
        aComponent.setID(anOriginal.getID());
        aComponent.set(anOriginal.flags());
    }
    //</editor-fold>
}
