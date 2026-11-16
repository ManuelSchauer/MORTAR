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

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.io.Importer;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Chemistry utility.
 *
 * @author Samuel Behr
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public final class ChemUtil {
    //<editor-fold desc="Public static final class constants">
    /**
     * Maximum number of different tautomers with explicit Hydrogen atoms added to aromatic Nitrogen atoms generated to
     * get a structure that can be kekulized in {@link #fixAromaticNitrogenAndCreateSMILES(IAtomContainer)}.
     */
    public static final int MAX_TAUTOMER_COMBINATIONS = 1000;

    /**
     * Pattern to match both 'n' and '[nH]' in SMILES (but not 'n' in Rn, Sn, In, Cn, Zn, Mn)
     */
    public static final String AROMATIC_N_REGEX = "\\[nH]|(?<!\\[[RSICZM])n";
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private static final class constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ChemUtil.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private ChemUtil() {
    }
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Public static methods">
    /**
     * Creates a unique SMILES string out of the given atom container or returns null, if the creation was not possible.
     * If the SMILES could not be created in the first place, it is retried with a kekulized clone of the given atom
     * container. Aromaticity information is encoded in the returned SMILES string, if there is any given. Unique SMILES
     * codes do NOT encode stereochemistry by default! This can be turned on with the second parameter. Stereo elements
     * of units that are not stereogenic are not encoded, see {@link #createUniqueSmiles(IAtomContainer, boolean, boolean)}.
     *
     * @param anAtomContainer atom container the unique SMILES should be created of
     * @param isStereoChemEncoded whether stereochemistry should be encoded
     * @return unique SMILES of the given atom container or 'null' if no creation was possible
     */
    public static String createUniqueSmiles(IAtomContainer anAtomContainer, boolean isStereoChemEncoded) {
        return ChemUtil.createUniqueSmiles(anAtomContainer, isStereoChemEncoded, true);
    }

    /**
     * Creates a unique SMILES string out of the given atom container or returns null, if the creation was not possible.
     * If the SMILES could not be created in the first place, it is retried with a kekulized clone of the given atom
     * container. Unique SMILES codes do NOT encode stereochemistry or aromaticity by default! This can be turned on
     * with the parameters.
     * <br>If stereochemistry is encoded, tetrahedral and double bond stereo elements of units that are not stereogenic
     * (e.g. the middle carbon atom of a chiral 1,5-anhydropentitol, whose two constitutionally equivalent ring
     * branches have the same configuration) are not encoded, so that the same compound always gets the same SMILES
     * code, regardless of how it was drawn or from which parent molecule it was cut. Real stereocentres, including
     * pseudo-asymmetric ones like the middle carbon atom of a meso 1,5-anhydropentitol, are kept. The given atom container is not
     * altered by this; if such stereo elements are present, the SMILES code is created for a clone of it.
     *
     * @param anAtomContainer atom container the unique SMILES should be created of
     * @param isStereoChemEncoded whether stereochemistry should be encoded
     * @param isAromaticityEncoded whether aromaticity should be encoded
     * @return unique SMILES of the given atom container or 'null' if no creation was possible
     */
    public static String createUniqueSmiles(IAtomContainer anAtomContainer, boolean isStereoChemEncoded, boolean isAromaticityEncoded) {
        int tmpFlavor = SmiFlavor.Unique;
        if (isAromaticityEncoded) {
            tmpFlavor = tmpFlavor | SmiFlavor.UseAromaticSymbols;
        }
        //the given atom container itself or a clone of it without the stereo elements of non-stereogenic units
        IAtomContainer tmpMoleculeToEncode = anAtomContainer;
        if (isStereoChemEncoded && anAtomContainer.stereoElements().iterator().hasNext()) {
            tmpMoleculeToEncode = ChemUtil.removeStereoElementsOfNonStereogenicUnits(anAtomContainer);
            /* checked again because all stereo elements might have been removed; the canonical atom order of SMILES
               with stereo flavor is different from the one without, so the same compound would otherwise get two
               different SMILES codes, depending on whether its atom container had redundant stereo elements or none */
            if (tmpMoleculeToEncode.stereoElements().iterator().hasNext()) {
                tmpFlavor = tmpFlavor | SmiFlavor.Stereo;
            }
        }
        String tmpSmiles = null;
        try {
            try {
                tmpSmiles = SmilesGenerator.create(tmpMoleculeToEncode, tmpFlavor, new int[tmpMoleculeToEncode.getAtomCount()]);
            } catch (CDKException anException) {
                String tmpFixedSmiles = ChemUtil.fixAromaticNitrogenAndCreateSMILES(tmpMoleculeToEncode);
                if (tmpFixedSmiles != null) {
                    IAtomContainer tmpFixedMol = ChemUtil.parseSmilesToAtomContainer(tmpFixedSmiles, true, false);
                    tmpSmiles = SmilesGenerator.create(tmpFixedMol, tmpFlavor, new int[tmpFixedMol.getAtomCount()]);
                    ChemUtil.LOGGER.log(Level.INFO, String.format("Added hydrogen atom to aromatic nitrogen atom to produce valid Kekule structure of molecule %s", (String) anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)));
                } else {
                    IAtomContainer tmpAtomContainer = tmpMoleculeToEncode.clone();
                    Kekulization.kekulize(tmpAtomContainer);
                    tmpSmiles = SmilesGenerator.create(tmpAtomContainer, tmpFlavor, new int[anAtomContainer.getAtomCount()]);
                    ChemUtil.LOGGER.log(Level.INFO, String.format("Kekulized molecule %s", anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)));
                }
            }
        } catch (CDKException | NullPointerException | IllegalArgumentException | CloneNotSupportedException | ArrayIndexOutOfBoundsException anException){
            ChemUtil.LOGGER.log(Level.SEVERE, String.format("%s; molecule name: %s", anException.toString(), anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)), anException);
        }
        return tmpSmiles;
    }

    /**
     * Returns an IAtomContainer instance which represents the molecule parsed from the SMILES string. Bond types and
     * atom types are assigned to it (the former through kekulization) if required. Aromaticity flags are set only if
     * there is aromaticity information present in the SMILES code, no aromaticity perception is performed here.
     *
     * @param aSmilesCode SMILES representation
     * @param shouldBeKekulized whether explicit bond orders should be assigned or "aromatic bond" can be used if present;
     *                          does not affect aromaticity flags
     * @param shouldAtomTypesBePerceived whether atom types should be perceived and configured
     * @return IAtomContainer atom container of the molecule
     * @throws CDKException if the given SMILES is invalid or if kekulization or atom type matching fails
     */
    public static IAtomContainer parseSmilesToAtomContainer(String aSmilesCode, boolean shouldBeKekulized, boolean shouldAtomTypesBePerceived)
            throws CDKException {
        //no checks because .parseSmiles() checks and throws InvalidSmilesException (subclass of CDKException) if the SMILES cannot be parsed
        IAtomContainer tmpAtomContainer;
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpSmiPar.kekulise(false);
        //throws InvalidSmilesException
        tmpAtomContainer = tmpSmiPar.parseSmiles(aSmilesCode);
        if (shouldBeKekulized) {
            //throws CDKException
            Kekulization.kekulize(tmpAtomContainer);
        }
        if (shouldAtomTypesBePerceived) {
            //throws CDKException
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpAtomContainer);
        }
        return tmpAtomContainer;
    }

    /**
     * Call to {@link #parseSmilesToAtomContainer(String, boolean, boolean)} with kekulization and atom type perception
     * set to true.
     *
     * @param aSmilesCode SMILES representation
     * @return IAtomContainer atom container of the molecule
     * @throws CDKException if the given SMILES is invalid or if kekulization or atom type matching fails
     */
    public static IAtomContainer parseSmilesToAtomContainer(String aSmilesCode) throws CDKException {
        return ChemUtil.parseSmilesToAtomContainer(aSmilesCode, true, true);
    }

    /**
     * Generates the molecular formula of a given atom container.
     *
     * @param anAtomContainer AtomContainer to generate the molecular formula of
     * @return String containing the molecular formula of the given atom container
     * @author Samuel Behr, Jonas Schaub
     */
    public static String generateMolecularFormula(IAtomContainer anAtomContainer) {
        /* the doc of MolecularFormulaManipulator.getMolecularFormula() says "The hydrogens must be implicit." but the
           CDK tests indicate that both implicit and explicit hydrogens are working fine */
        IMolecularFormula tmpMolecularFormula = MolecularFormulaManipulator.getMolecularFormula(anAtomContainer);
        return MolecularFormulaManipulator.getString(tmpMolecularFormula);
    }

    /**
     * Checks whether 3D coordinates are set for all atoms in the given molecule data model.
     *
     * @param aMolecule to check for 3D coordinates
     * @return true if 3D coordinates are set for ALL atoms in the given molecule
     */
    public static boolean has3DCoordinates(MoleculeDataModel aMolecule) {
        IAtomContainer tmpFragment;
        try {
            tmpFragment = aMolecule.getAtomContainer();
        } catch(CDKException anException){
            ChemUtil.LOGGER.log(Level.SEVERE, String.format("%s molecule name: %s", anException.toString(), aMolecule.getName()), anException);
            return false;
        }
        boolean tmpHas3DCoords = true;
        for (IAtom tmpAtom : tmpFragment.atoms()) {
            if (tmpAtom.getPoint3d() == null) {
                tmpHas3DCoords = false;
                break;
            }
        }
        return tmpHas3DCoords;
    }

    /**
     * Checks whether 2D coordinates are set for all atoms in the given molecule data model.
     *
     * @param aMolecule to check for 2D coordinates
     * @return true if 2D coordinates are set for ALL atoms in the given molecule
     */
    public static boolean has2DCoordinates(MoleculeDataModel aMolecule) {
        IAtomContainer tmpFragment;
        try {
            tmpFragment = aMolecule.getAtomContainer();
        } catch(CDKException anException){
            ChemUtil.LOGGER.log(Level.SEVERE, String.format("%s molecule name: %s", anException.toString(), aMolecule.getName()), anException);
            return false;
        }
        boolean tmpHas2DCoords = true;
        for (IAtom tmpAtom : tmpFragment.atoms()) {
            if (tmpAtom.getPoint2d() == null) {
                tmpHas2DCoords = false;
                break;
            }
        }
        return tmpHas2DCoords;
    }

    /**
     * Checks all atoms in all molecules in the given list for 3D or 2D coordinates. If only one atom has neither,
     * false is returned.
     *
     * @param aListOfMolecules to check for 2D or 3D coordinates
     * @return true if 2D or 3D coordinates are set for EVERY atom in EVERY molecule in the given list
     */
    public static boolean checkMoleculeListForCoordinates(List<MoleculeDataModel> aListOfMolecules) {
        if (aListOfMolecules == null || aListOfMolecules.isEmpty()) {
            return false;
        }
        boolean tmpHasCoords = true;
        for (MoleculeDataModel tmpMolecule : aListOfMolecules) {
            if (!ChemUtil.has3DCoordinates(tmpMolecule) && !ChemUtil.has2DCoordinates(tmpMolecule)) {
                tmpHasCoords = false;
                break;
            }
        }
        return tmpHasCoords;
    }

    /**
     * Generates 2D coordinates for the atoms of the given molecule using the CDK StructureDiagramGenerator. Therefore,
     * these coordinates are originally intended for layout! Note that the given atm container is directly manipulated,
     * not cloned. And note that the z coordinate is not set to 0 but left undefined.
     *
     * @param aMolecule the molecule to generate coordinates for
     * @throws NullPointerException if given molecule is null
     * @throws CDKException if coordinates generation fails
     */
    public static void generate2DCoordinates(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            return;
        }
        StructureDiagramGenerator tmpStructureDiagramGenerator = new StructureDiagramGenerator();
        tmpStructureDiagramGenerator.generateCoordinates(aMolecule);
    }

    /**
     * Generates pseudo-3D coordinates for the atoms of the given molecule by retrieving their 2D coordinates and
     * setting z=0. Note that the given atm container is directly manipulated, not cloned.
     *
     * @param aMolecule the molecule to generate coordinates for
     * @throws NullPointerException if given molecule is null
     * @throws IllegalArgumentException if at least one atom of the molecule has no 2D coordinates defined
     */
    public static void generatePseudo3Dfrom2DCoordinates(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            return;
        }
        for (IAtom tmpAtom : aMolecule.atoms()) {
            Point2d tmpAtom2Dcoords = tmpAtom.getPoint2d();
            if (Objects.isNull(tmpAtom2Dcoords)) {
                throw new IllegalArgumentException("At least one atom has no 2D coordinates.");
            }
            Point3d tmpPoint3d = new Point3d(tmpAtom2Dcoords.x, tmpAtom2Dcoords.y, 0.0);
            tmpAtom.setPoint3d(tmpPoint3d);
        }
    }

    /**
     * Sets the 3D coordinates of all atoms in the given molecule to (0,0,0).
     * Note that the given atm container is directly manipulated, not cloned.
     * Use this with caution, if a molecular simulation with these molecules is attempted, it might crash!
     *
     * @param aMolecule the molecule to set (0,0,0) coordinates for
     * @throws NullPointerException if given molecule is null
     */
    public static void generateZero3DCoordinates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            return;
        }
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpAtom.setPoint3d(new Point3d(0.0, 0.0, 0.0));
        }
    }

    /**
     * Saturates free valences in the given molecule or molecular fragment with implicit hydrogen atoms.
     *
     * @param aMolecule to saturate
     * @throws NullPointerException if the given molecule is null
     * @throws CDKException if atom types cannot be assigned to all atoms of the molecule
     */
    public static void saturateWithHydrogen(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            //fragments might be empty on purpose, e.g. when there is no aglycone in a molecule, so throw no exception
            return;
        }
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(aMolecule);
        //just a precaution, you never know what can happen with fragment atom containers
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getImplicitHydrogenCount() == CDKConstants.UNSET) {
                tmpAtom.setImplicitHydrogenCount(0);
            }
        }
    }

    /**
     * Checks whether atoms in the given molecule have free atom pairs correctly assigned if chemically needed. If not,
     * they are added. Note: it does not work really well on most fragments...
     *
     * @param aMolecule to check for missing free electron pairs
     * @throws NullPointerException if the given molecule is null
     * @throws CDKException if atom types cannot be assigned to all atoms of the molecule
     */
    public static void checkAndCorrectElectronConfiguration(IAtomContainer aMolecule)
            throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            //fragments might be empty on purpose, e.g. when there is no aglycon in a molecule, so throw no exception
            return;
        }
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        LonePairElectronChecker tmpElectronChecker = new LonePairElectronChecker();
        try {
            tmpElectronChecker.saturate(aMolecule);
        } catch (NullPointerException aNullPointerException) {
            /* Too many fragment molecules cause this exception, the log file would be overflowing
            ChemUtil.LOGGER.log(Level.WARNING, "Saturation with lone electron pairs not possible, molecule name: "
                    + aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), aNullPointerException);
             */
        }
    }

    /**
     * Fixes radical atoms in a molecule by replacing the single electron with an implicit hydrogen atom.
     * This method converts radical centers to closed-shell configurations by:
     * <ol>
     * <li>Adding an implicit hydrogen to each radical center</li>
     * <li>Removing all single electrons from the molecule</li>
     * <li>Reperceiving atom types to ensure correct configuration</li>
     * </ol>
     * This process effectively caps each radical site with a hydrogen atom, resulting in
     * a more stable molecular representation suitable for further processing.
     * <br><br><b>Important note: expects atom types to be perceived beforehand!</b>
     *
     * @param aMolecule The molecule to process; will be modified in place
     * @throws NullPointerException If the provided molecule is null
     * @throws CDKException If atom type perception fails or other CDK operations encounter problems
     */
    public static void fixRadicals(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            return;
        }
        if (aMolecule.getSingleElectronCount() > 0) {
            //fix properties of the atoms that are radicals
            for (ISingleElectron tmpSingleElectron : aMolecule.singleElectrons()) {
                IAtom tmpAtom = tmpSingleElectron.getAtom();
                //setting to null now, will be re-detected correctly by AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule) below
                tmpAtom.setHybridization(null);
                tmpAtom.setValency(null);
                tmpAtom.setFormalNeighbourCount(null);
                Integer tmpHCount = tmpAtom.getImplicitHydrogenCount();
                if (tmpHCount == null) {
                    tmpHCount = 0;
                }
                tmpAtom.setImplicitHydrogenCount(tmpHCount + 1);
            }
            //remove all single electrons from the molecule
            int tmpSingleElectronCount = aMolecule.getSingleElectronCount();
            // the electron array is re-ordered after a removal, so we need to remove them in reverse order
            for (int i = tmpSingleElectronCount - 1; i >= 0; i--) {
                aMolecule.removeSingleElectron(i);
            }
            //needs to be redone now to set the correct atom types
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
            ChemUtil.LOGGER.log(Level.INFO, "{0}", String.format("Fixed %d radicals in molecule with name %s.",
                    tmpSingleElectronCount, aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)));
        } //else: do nothing
    }

    /**
     * For cases where an atom container cannot be kekulized (which is among other things necessary for absolute SMILES code generation)
     * because of insufficient valence saturation of aromatic nitrogen atoms (i.e. the SMILES code that was originally
     * parsed was missing an "[nH]" instead of "n"), this method tries to fix this by adding such explicit hydrogen atoms
     * to aromatic nitrogen atoms in all possible combinations. If successful, a valid (*non-unique!*) SMILES code will
     * be returned that encodes aromaticity and stereochemistry. If this does not fix the non-kekulizable molecule or
     * the molecule does not contain aromatic nitrogen atoms, null is returned. This method does not check whether the
     * input molecule can be kekulized or not. If the molecule cannot be fixed with the first 1,000 possible solutions
     * generated, the routine is aborted and null is returned. This limit is imposed to prevent excessive computation time
     * (consider the exponential scaling of possible solutions) and most molecules should be fixed within this margin.
     * Note: we generate a SMILES code first, operate on it to generate all the possible solutions, and check for validity
     * by checking whether it can be parsed again by a SmilesParser that is kekulizing. This round-trip may not be optimal
     * but operating on the atom container directly would require making multiple copies of it. This solution here
     * is therefore less memory-intensive. Plus, since the problem originated from a faulty SMILES code at import in most
     * cases, it is also easiest to correct it this way.
     *
     * @param aMolecule the molecule to attempt to fix
     * @return a (*non-unique!*) SMILES code if successful, null otherwise
     * @throws NullPointerException if the given molecule is null
     */
    public static String fixAromaticNitrogenAndCreateSMILES(IAtomContainer aMolecule) throws NullPointerException{
        //checks:
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        if (aMolecule.isEmpty()) {
            return null;
        }
        //count aromatic nitrogen atoms to choose appropriate initial collection sizes below and reject molecules without any
        int tmpAromaticNCount = 0;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getAtomicNumber().equals(IElement.N) && tmpAtom.isAromatic()) {
                tmpAromaticNCount++;
            }
        }
        if (tmpAromaticNCount == 0) {
            return null;
        }
        //end of checks
        //generate SMILES code of molecule to operate on it; encodes aromaticity and stereochemistry but is not(!) unique
        String tmpSmiles;
        try {
            tmpSmiles = SmilesGenerator.create(aMolecule, SmiFlavor.UseAromaticSymbols | SmiFlavor.Stereo, new int[aMolecule.getAtomCount()]);
        } catch (CDKException anException) {
            return null;
        }
        // Pattern to match both 'n' and '[nH]' in SMILES (but not 'n' in Rn, Sn, In, Cn, Zn, Mn)
        Pattern tmpNPattern = Pattern.compile(ChemUtil.AROMATIC_N_REGEX);
        Matcher tmpNMatcher = tmpNPattern.matcher(tmpSmiles);
        // Find all aromatic nitrogen positions
        List<Integer> tmpAromaticNPositions = new ArrayList<>(tmpAromaticNCount);
        List<String> tmpAromaticNTypes = new ArrayList<>(tmpAromaticNCount);
        while (tmpNMatcher.find()) {
            tmpAromaticNPositions.add(tmpNMatcher.start());
            tmpAromaticNTypes.add(tmpNMatcher.group());
        }
        if (tmpAromaticNPositions.isEmpty()) {
            return null;
        }
        //SmilesParser for validation
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //needs to fail if kekulization is not possible
        tmpSmiPar.kekulise(true);
        // Generate all 2^n combinations but already validate while generating and return the first valid solution
        int tmpNCount = tmpAromaticNPositions.size();
        int tmpNrOfTotalCombinations = Math.min((int) Math.pow(2, tmpNCount), ChemUtil.MAX_TAUTOMER_COMBINATIONS);
        tautomerLoop:
        for (int i = 0; i < tmpNrOfTotalCombinations; i++) {
            StringBuilder tmpTautomerBuilder = new StringBuilder(tmpSmiles);
            int tmpOffset = 0;
            for (int j = 0; j < tmpNCount; j++) {
                int tmpPos = tmpAromaticNPositions.get(j) + tmpOffset;
                String tmpCurrentType = tmpAromaticNTypes.get(j);
                boolean tmpIsUseNH = ((i >> j) & 1) == 1;
                String tmpReplacement = tmpIsUseNH ? "[nH]" : "n";
                int tmpOldLength = tmpCurrentType.length();
                int tmpNewLength = tmpReplacement.length();
                tmpTautomerBuilder.replace(tmpPos, tmpPos + tmpOldLength, tmpReplacement);
                tmpOffset += (tmpNewLength - tmpOldLength);
            }
            String tmpTautomerSMILESCode = tmpTautomerBuilder.toString();
            //validation:
            IAtomContainer tmpTautomerAtomContainer;
            try {
                tmpTautomerAtomContainer = tmpSmiPar.parseSmiles(tmpTautomerSMILESCode);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpTautomerAtomContainer);
            } catch (CDKException | NullPointerException anException) {
                continue;
            }
            for (IAtom tmpAtom : tmpTautomerAtomContainer.atoms()) {
                if (tmpAtom.getAtomicNumber().equals(IElement.N) && tmpAtom.isAromatic()) {
                    Integer tmpImplicitH = tmpAtom.getImplicitHydrogenCount();
                    int tmpValence = tmpTautomerAtomContainer.getConnectedBondsCount(tmpAtom) + (tmpImplicitH == null ? 0 : tmpImplicitH);
                    int tmpCharge = tmpAtom.getFormalCharge();
                    if (tmpValence == 4 && tmpCharge == 0) {
                        //the routine created a Nitrogen atom with valence 4 but no charge - invalid!
                        continue tautomerLoop;
                    }
                }
            }
            //if it can be parsed with kekulization, it is valid and therefore returned
            return tmpTautomerSMILESCode;
        }
        //no tautomer was valid, return null
        return null;
    }
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private static methods">
    /**
     * Returns the given atom container itself, if all its tetrahedral and double bond stereo elements describe
     * stereogenic units, or otherwise a clone of it without the stereo elements of the units that are not stereogenic.
     * The given atom container is not altered.
     * <br>A stereo element is regarded as describing a non-stereogenic unit if the molecule with only this element
     * inverted is identical to the molecule as given, including all its stereo elements, i.e. if there is an
     * automorphism of the constitution that maps the one onto the other. Example: the middle carbon atom of a chiral
     * 1,5-anhydropentitol, whose two ring branches have the same configuration; its configuration is arbitrary and
     * would otherwise lead to two different SMILES codes for one compound. In a meso 1,5-anhydropentitol, the two
     * branches are mirror images of each other, the inversion of the middle carbon atom leads to a different compound
     * and its pseudo-asymmetric configuration is kept. The test is applied only to elements whose central atom (or one
     * of whose double bond atoms) has at least two constitutionally equivalent ligands, since all other elements are
     * stereogenic. Each element is tested against the molecule with all its stereo elements as given, so the result
     * does not depend on the order of the elements. Pseudo atoms are equivalent to each other irrespective of their
     * labels and mass numbers, isotopes are ignored in general, as in the unique SMILES codes created in
     * {@link #createUniqueSmiles(IAtomContainer, boolean, boolean)}.
     * <br>Tetrahedral stereo elements with repeated carriers (e.g. the central atom given repeatedly instead of missing
     * ligands, as in some alkyl fragments) can never describe a stereocentre; they are removed without a test and are
     * not regarded in the tests of the other elements.
     * <br>If the molecule has stereo elements of other types (e.g. extended tetrahedral or atropisomeric), which the
     * stereo-aware isomorphism check of the CDK does not regard, or if the test fails or exceeds its search budget,
     * the given atom container is returned unchanged.
     *
     * @param anAtomContainer molecule with at least one stereo element
     * @return the given atom container or a clone of it without the stereo elements of non-stereogenic units
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private static IAtomContainer removeStereoElementsOfNonStereogenicUnits(IAtomContainer anAtomContainer) {
        List<IStereoElement> tmpElements = new ArrayList<>();
        boolean tmpHasGroupInfo = false;
        for (IStereoElement tmpElement : anAtomContainer.stereoElements()) {
            int tmpConfigClass = tmpElement.getConfigClass();
            if (tmpConfigClass != IStereoElement.TH && tmpConfigClass != IStereoElement.CT) {
                return anAtomContainer;
            }
            if (tmpElement.getGroupInfo() != 0) {
                tmpHasGroupInfo = true;
            }
            tmpElements.add(tmpElement);
        }
        int[][] tmpGraph;
        long[] tmpSymmetryClasses;
        try {
            tmpGraph = GraphUtil.toAdjList(anAtomContainer);
            //pseudo atoms all get atomic number 0, mass numbers are not regarded by default
            tmpSymmetryClasses = Canon.symmetry(anAtomContainer, tmpGraph);
        } catch (NullPointerException | IllegalArgumentException anException) {
            ChemUtil.LOGGER.log(Level.FINE, String.format("Symmetry classes could not be determined for molecule %s: %s",
                    anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException));
            return anAtomContainer;
        }
        boolean[] tmpIsDegenerateArray = new boolean[tmpElements.size()];
        boolean tmpHasDegenerateElements = false;
        List<Integer> tmpCandidateIndices = new ArrayList<>(tmpElements.size());
        for (int i = 0; i < tmpElements.size(); i++) {
            IStereoElement tmpElement = tmpElements.get(i);
            if (tmpElement.getConfigClass() == IStereoElement.TH
                    && new HashSet<>(tmpElement.getCarriers()).size() < tmpElement.getCarriers().size()) {
                //e.g. the central atom given repeatedly for missing ligands, can never describe a stereocentre
                tmpIsDegenerateArray[i] = true;
                tmpHasDegenerateElements = true;
            } else if (ChemUtil.hasConstitutionallyEquivalentLigands(tmpElement, anAtomContainer, tmpGraph, tmpSymmetryClasses)) {
                tmpCandidateIndices.add(i);
            }
        }
        if (tmpCandidateIndices.isEmpty() && !tmpHasDegenerateElements) {
            return anAtomContainer;
        }
        //slow path: remove degenerate elements and test the candidates by a stereo-aware automorphism search
        IAtomContainer tmpQuery;
        IAtomContainer tmpTarget;
        try {
            tmpQuery = anAtomContainer.clone();
            tmpTarget = anAtomContainer;
            if (tmpHasGroupInfo || (tmpHasDegenerateElements && !tmpCandidateIndices.isEmpty())) {
                tmpTarget = anAtomContainer.clone();
            }
        } catch (CloneNotSupportedException anException) {
            ChemUtil.LOGGER.log(Level.FINE, String.format("Molecule %s could not be cloned: %s",
                    anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException));
            return anAtomContainer;
        }
        List<IStereoElement> tmpQueryElements = ChemUtil.removeDegenerateStereoElements(tmpQuery, tmpIsDegenerateArray);
        if (tmpQueryElements == null) {
            return anAtomContainer;
        }
        if (tmpTarget != anAtomContainer) {
            List<IStereoElement> tmpTargetElements = ChemUtil.removeDegenerateStereoElements(tmpTarget, tmpIsDegenerateArray);
            if (tmpTargetElements == null) {
                return anAtomContainer;
            }
            //enhanced stereo groups of the target would allow the matching of a group as a whole in inverted form
            for (IStereoElement tmpElement : tmpTargetElements) {
                if (tmpElement != null) {
                    tmpElement.setGroupInfo(0);
                }
            }
        }
        boolean[] tmpIsRedundantArray = new boolean[tmpQueryElements.size()];
        boolean tmpIsAnyElementRedundant = false;
        if (!tmpCandidateIndices.isEmpty()) {
            Map<IAtom, Long> tmpSymmetryClassMap = new IdentityHashMap<>(4 * anAtomContainer.getAtomCount());
            for (int i = 0; i < anAtomContainer.getAtomCount(); i++) {
                tmpSymmetryClassMap.put(anAtomContainer.getAtom(i), tmpSymmetryClasses[i]);
                tmpSymmetryClassMap.put(tmpQuery.getAtom(i), tmpSymmetryClasses[i]);
                if (tmpTarget != anAtomContainer) {
                    tmpSymmetryClassMap.put(tmpTarget.getAtom(i), tmpSymmetryClasses[i]);
                }
            }
            SymmetryClassAtomMatcher tmpAtomMatcher = new SymmetryClassAtomMatcher(tmpSymmetryClassMap);
            org.openscience.cdk.isomorphism.Pattern tmpPattern = VentoFoggia.findIdentical(tmpQuery, tmpAtomMatcher, BondMatcher.forOrder());
            for (int tmpIndex : tmpCandidateIndices) {
                IStereoElement tmpQueryElement = tmpQueryElements.get(tmpIndex);
                int tmpConfigOrder = tmpQueryElement.getConfigOrder();
                if (tmpQueryElement.getConfigClass() != tmpElements.get(tmpIndex).getConfigClass()
                        || (tmpConfigOrder != IStereoElement.LEFT && tmpConfigOrder != IStereoElement.RIGHT)) {
                    continue;
                }
                //LEFT/RIGHT for tetrahedral and OPPOSITE/TOGETHER for double bond stereo elements have the same values
                tmpQueryElement.setConfigOrder(tmpConfigOrder == IStereoElement.LEFT ? IStereoElement.RIGHT : IStereoElement.LEFT);
                tmpAtomMatcher.resetBudget();
                boolean tmpIsInversionIdentical;
                try {
                    tmpIsInversionIdentical = tmpPattern.matches(tmpTarget);
                } catch (RuntimeException anException) {
                    ChemUtil.LOGGER.log(Level.FINE, String.format("Stereo-aware automorphism search failed for molecule %s: %s",
                            anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), anException));
                    tmpIsInversionIdentical = false;
                } finally {
                    tmpQueryElement.setConfigOrder(tmpConfigOrder);
                }
                if (tmpAtomMatcher.isBudgetExhausted()) {
                    ChemUtil.LOGGER.log(Level.INFO, String.format("Search budget for the stereogenicity test exhausted, stereo element kept in molecule %s",
                            anAtomContainer.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)));
                    continue;
                }
                if (tmpIsInversionIdentical) {
                    tmpIsRedundantArray[tmpIndex] = true;
                    tmpIsAnyElementRedundant = true;
                }
            }
        }
        if (!tmpIsAnyElementRedundant && !tmpHasDegenerateElements) {
            return anAtomContainer;
        }
        List<IStereoElement> tmpRemainingElements = new ArrayList<>(tmpQueryElements.size());
        for (int i = 0; i < tmpQueryElements.size(); i++) {
            if (!tmpIsRedundantArray[i] && !tmpIsDegenerateArray[i]) {
                tmpRemainingElements.add(tmpQueryElements.get(i));
            }
        }
        tmpQuery.setStereoElements(tmpRemainingElements);
        return tmpQuery;
    }

    /**
     * Removes the stereo elements marked as degenerate from the given clone and returns its stereo elements in a list
     * that is aligned with the stereo elements of the original atom container, i.e. that has the same size and order,
     * with null at the positions of the removed elements. The elements are read from the clone after the removal,
     * since setting the stereo elements of an atom container may store copies of them. Returns null if the number of
     * stereo elements of the clone does not match the given array.
     *
     * @param aClone clone of the atom container whose stereo elements were marked
     * @param anIsDegenerateArray marks for the stereo elements of the original atom container, in their order
     * @return aligned list of the stereo elements of the clone or null
     */
    @SuppressWarnings("rawtypes")
    private static List<IStereoElement> removeDegenerateStereoElements(IAtomContainer aClone, boolean[] anIsDegenerateArray) {
        List<IStereoElement> tmpAllElements = new ArrayList<>(anIsDegenerateArray.length);
        for (IStereoElement tmpElement : aClone.stereoElements()) {
            tmpAllElements.add(tmpElement);
        }
        if (tmpAllElements.size() != anIsDegenerateArray.length) {
            return null;
        }
        List<IStereoElement> tmpValidElements = new ArrayList<>(tmpAllElements.size());
        for (int i = 0; i < tmpAllElements.size(); i++) {
            if (!anIsDegenerateArray[i]) {
                tmpValidElements.add(tmpAllElements.get(i));
            }
        }
        if (tmpValidElements.size() == tmpAllElements.size()) {
            return tmpAllElements;
        }
        aClone.setStereoElements(tmpValidElements);
        List<IStereoElement> tmpAlignedElements = new ArrayList<>(anIsDegenerateArray.length);
        int tmpValidIndex = 0;
        List<IStereoElement> tmpElementsAfterRemoval = new ArrayList<>(tmpValidElements.size());
        for (IStereoElement tmpElement : aClone.stereoElements()) {
            tmpElementsAfterRemoval.add(tmpElement);
        }
        if (tmpElementsAfterRemoval.size() != tmpValidElements.size()) {
            return null;
        }
        for (boolean tmpIsDegenerate : anIsDegenerateArray) {
            tmpAlignedElements.add(tmpIsDegenerate ? null : tmpElementsAfterRemoval.get(tmpValidIndex++));
        }
        return tmpAlignedElements;
    }

    /**
     * Checks whether the central atom of the given tetrahedral stereo element has at least two constitutionally
     * equivalent ligands or whether one of the two atoms of the given double bond stereo element has two
     * constitutionally equivalent substituents. Only such elements can describe non-stereogenic units. An implicit
     * hydrogen atom, given as the central atom itself in the carriers of a tetrahedral element, is not compared.
     *
     * @param aStereoElement tetrahedral or double bond stereo element
     * @param anAtomContainer molecule the stereo element belongs to
     * @param aGraph adjacency list of the molecule
     * @param aSymmetryClassArray constitutional symmetry classes of the atoms of the molecule
     * @return true if the element has constitutionally equivalent ligands or substituents
     */
    @SuppressWarnings("rawtypes")
    private static boolean hasConstitutionallyEquivalentLigands(IStereoElement aStereoElement,
                                                                IAtomContainer anAtomContainer,
                                                                int[][] aGraph,
                                                                long[] aSymmetryClassArray) {
        if (aStereoElement.getConfigClass() == IStereoElement.TH) {
            IAtom tmpFocus = (IAtom) aStereoElement.getFocus();
            List<Long> tmpLigandClasses = new ArrayList<>(4);
            for (Object tmpCarrier : aStereoElement.getCarriers()) {
                IAtom tmpLigand = (IAtom) tmpCarrier;
                if (tmpLigand.equals(tmpFocus)) {
                    continue;
                }
                int tmpLigandIndex = anAtomContainer.indexOf(tmpLigand);
                if (tmpLigandIndex < 0) {
                    return false;
                }
                Long tmpLigandClass = aSymmetryClassArray[tmpLigandIndex];
                if (tmpLigandClasses.contains(tmpLigandClass)) {
                    return true;
                }
                tmpLigandClasses.add(tmpLigandClass);
            }
            return false;
        } else if (aStereoElement.getConfigClass() == IStereoElement.CT) {
            IBond tmpDoubleBond = (IBond) aStereoElement.getFocus();
            int tmpBeginIndex = anAtomContainer.indexOf(tmpDoubleBond.getBegin());
            int tmpEndIndex = anAtomContainer.indexOf(tmpDoubleBond.getEnd());
            if (tmpBeginIndex < 0 || tmpEndIndex < 0) {
                return false;
            }
            int[][] tmpAtomPairs = {{tmpBeginIndex, tmpEndIndex}, {tmpEndIndex, tmpBeginIndex}};
            for (int[] tmpAtomPair : tmpAtomPairs) {
                long tmpFirstSubstituentClass = -1;
                int tmpSubstituentCount = 0;
                for (int tmpNeighbourIndex : aGraph[tmpAtomPair[0]]) {
                    if (tmpNeighbourIndex == tmpAtomPair[1]) {
                        continue;
                    }
                    tmpSubstituentCount++;
                    if (tmpSubstituentCount == 1) {
                        tmpFirstSubstituentClass = aSymmetryClassArray[tmpNeighbourIndex];
                    } else if (tmpSubstituentCount == 2 && tmpFirstSubstituentClass == aSymmetryClassArray[tmpNeighbourIndex]) {
                        return true;
                    }
                }
            }
            return false;
        }
        return false;
    }
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private static class">
    /**
     * Atom matcher for the stereo-aware automorphism search in
     * {@link #removeStereoElementsOfNonStereogenicUnits(IAtomContainer)}: two atoms match if they are in the same
     * constitutional symmetry class. This restricts the search to automorphisms of the constitution. The number of
     * comparisons per search is limited, since the stereo configurations are checked only for complete mappings and
     * molecules with many equivalent groups (e.g. many tert-butyl or phosphate groups) have very many automorphisms.
     * If the budget is exhausted, no further atoms match and the search ends without a match.
     */
    private static final class SymmetryClassAtomMatcher extends AtomMatcher {
        /**
         * Maximum number of atom comparisons per search.
         */
        private static final int MAX_COMPARISONS_PER_SEARCH = 1_000_000;

        /**
         * Symmetry classes of the atoms of the query and the target.
         */
        private final Map<IAtom, Long> symmetryClassMap;

        /**
         * Remaining number of atom comparisons of the current search.
         */
        private int remainingComparisons;

        /**
         * Constructor.
         *
         * @param aSymmetryClassMap symmetry classes of the atoms of the query and the target
         */
        SymmetryClassAtomMatcher(Map<IAtom, Long> aSymmetryClassMap) {
            this.symmetryClassMap = aSymmetryClassMap;
            this.remainingComparisons = SymmetryClassAtomMatcher.MAX_COMPARISONS_PER_SEARCH;
        }

        /**
         * Resets the number of remaining atom comparisons for a new search.
         */
        void resetBudget() {
            this.remainingComparisons = SymmetryClassAtomMatcher.MAX_COMPARISONS_PER_SEARCH;
        }

        /**
         * Returns whether the budget of atom comparisons of the current search is exhausted.
         *
         * @return true if the budget is exhausted
         */
        boolean isBudgetExhausted() {
            return this.remainingComparisons <= 0;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public boolean matches(IAtom anAtom1, IAtom anAtom2) {
            if (this.remainingComparisons <= 0) {
                return false;
            }
            this.remainingComparisons--;
            Long tmpClass1 = this.symmetryClassMap.get(anAtom1);
            return tmpClass1 != null && tmpClass1.equals(this.symmetryClassMap.get(anAtom2));
        }
    }
    //</editor-fold>
}
