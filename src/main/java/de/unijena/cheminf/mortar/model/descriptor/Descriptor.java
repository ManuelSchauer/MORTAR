/*
 * Copyright (C) 2026 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package de.unijena.cheminf.mortar.model.descriptor;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SpiroAtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

/**
 * Catalog of the CDK descriptors and fingerprints that can be calculated with the {@link DescriptorCalculator}.
 * <p>
 *     Each enum constant carries seven fields that describe its characteristics and can be queried via the
 *     corresponding getter methods:
 *     <ul>
 *         <li><b>isFast</b> ({@code boolean}) – whether this descriptor can be calculated quickly.
 *             The classification is based on a performance benchmark: descriptors that complete computation
 *             for 10,000 molecules in under one second are considered fast.
 *             Slow descriptors (e.g., {@link #CHI_CHAIN}, {@link #CHI_CLUSTER}, {@link #CHI_PATH_CLUSTER},
 *             {@link #CHI_PATH}, {@link #WEIGHTED_PATH}) require more intensive computation.
 *             Use {@link #isFast()} to query this field.</li>
 *         <li><b>isSafe</b> ({@code boolean}) – whether this descriptor reliably produces a valid result
 *             for a valid molecule input. Unsafe descriptors (e.g., {@link #BCUT}, {@link #KAPPA_SHAPE_INDICES})
 *             may produce {@code NaN} values for certain molecular inputs they are not developed for.
 *             Use {@link #isSafe()} to query this field.</li>
 *         <li><b>isFingerprint</b> ({@code boolean}) – whether this descriptor generates a binary
 *             fingerprint bit vector rather than a set of numerical descriptor values. Fingerprint
 *             descriptors (e.g., {@link #PUBCHEM_FINGERPRINTER}, {@link #MACCS_FINGERPRINTER},
 *             all {@code CIRCULAR_FINGERPRINTER_*} variants) are handled via a thread-safe pool of
 *             {@link IFingerprinter} instances internally. By default, fingerprints are calculated as bit
 *             fingerprints (components 0 or 1); circular fingerprints can be calculated as folded count fingerprints
 *             instead, see {@link DescriptorCalculator#DescriptorCalculator(int, int, boolean)}.
 *             Use {@link #isFingerprint()} to query this field.</li>
 *         <li><b>needsExplicitHydrogens</b> ({@code boolean}) – whether this descriptor requires
 *             explicit hydrogen atoms to be present in the molecule before calculation.
 *             Examples include {@link #A_LOG_P}, {@link #X_LOG_P}, {@link #KAPPA_SHAPE_INDICES},
 *             and {@link #HYBRIDIZATION_RATIO}.
 *             Use {@link #needsExplicitHydrogens()} to query this field.</li>
 *         <li><b>requires3DCoordinates</b> ({@code boolean}) – whether this descriptor requires
 *             three-dimensional coordinates to be present in the molecule before calculation.
 *             Use {@link #requires3DCoordinates()} to query this field.</li>
 *         <li><b>descriptorComponentNumber</b> ({@code int}) – the number of individual numerical
 *             values produced by this descriptor. Most descriptors return a single value (e.g.,
 *             {@link #MOLECULAR_WEIGHT}, {@link #TPSA}), while multi-value descriptors return arrays
 *             (e.g., {@link #BCUT} returns 6, {@link #KAPPA_SHAPE_INDICES} returns 3,
 *             {@link #CHI_PATH} returns 16, {@link #KIER_HALL_SMARTS} returns 79,
 *             {@link #PUBCHEM_FINGERPRINTER} returns 881).
 *             Use {@link #getDescriptorComponentNumber()} to query this field.</li>
 *         <li><b>componentType</b> ({@link DescriptorComponentType}) – the type of the values of all components:
 *             {@code BINARY} for fingerprint bits, {@code COUNT} for non-negative integer counts (e.g.
 *             {@link #ATOM_COUNT_C}, {@link #SMALL_RING}) and {@code CONTINUOUS} for real values (e.g.
 *             {@link #MOLECULAR_WEIGHT}, {@link #TPSA}). Circular fingerprints are {@code BINARY} by default and
 *             {@code COUNT} if the calculator calculates them as count fingerprints.
 *             Use {@link #getComponentType()} to query this field.</li>
 *         <li><b>name</b> ({@code String}) – a human-readable display name for this descriptor,
 *             intended for labeling output columns or user-facing representations (e.g.,
 *             {@code "Molecular Weight"}, {@code "TPSA"}, {@code "ALogP"}).
 *             Use {@link #getName()} to query this field.</li>
 *     </ul>
 * <p>
 *     The calculation of descriptor values, e.g. for filling a data matrix of molecules, is done by the
 *     {@link DescriptorCalculator} class.
 * </p>
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
public enum Descriptor {
    /*
     * Note for developers: for adding a new descriptor,
     * 1. add the enum constant here, at "Add new descriptor information here!", including the type of its component
     *    values (DescriptorComponentType) and the recipe for its calculation, i.e. molecular(...) for a CDK
     *    molecular descriptor, fingerprint(...) for a CDK fingerprinter or circular(...) for a CDK circular
     *    fingerprint,
     * 2. add tests in the test class DescriptorCalculatorTest, at "Add new descriptor tests here!" and
     *    "Add new descriptor information here!".
     */
    /**
     * Molecular weight, adds up the natural masses (weighted average of all known
     * isotopes of the particular element based on their natural abundances) of every atom
     * in the given molecule, it is NOT the exact mass.
     *
     * @see WeightDescriptor
     */
    MOLECULAR_WEIGHT(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Molecular Weight",
            molecular(WeightDescriptor::new)),
    /**
     * Wiener number, returns Wiener path number and Wiener polarity number.
     * Path number: sum of the distances between any two atoms in the molecule.<br>
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.<br>
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     *
     * @see WienerNumbersDescriptor
     */
    WIENER_NUMBER(true, true, false, false, false, 2, DescriptorComponentType.CONTINUOUS, "Wiener Number",
            molecular(WienerNumbersDescriptor::new)),
    /**
     * Atom count, counts the number of all atoms in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count",
            molecular(AtomCountDescriptor::new)),
    /**
     * Atom count heavy, counts the number of all heavy atoms in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_HEAVY(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count Heavy",
            molecular(() -> configured(new AtomCountDescriptor(), "#"))),
    /**
     * Atom count C, counts the number of all carbon atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_C(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count C",
            molecular(() -> configured(new AtomCountDescriptor(), "C"))),
    /**
     * Atom count H, counts the number of all hydrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_H(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count H",
            molecular(() -> configured(new AtomCountDescriptor(), "H"))),
    /**
     * Atom count N counts the number of all nitrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_N(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count N",
            molecular(() -> configured(new AtomCountDescriptor(), "N"))),
    /**
     * Atom count O, counts the number of all oxygen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_O(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count O",
            molecular(() -> configured(new AtomCountDescriptor(), "O"))),
    /**
     * Atom count S, counts the number of all sulfur atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_S(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count S",
            molecular(() -> configured(new AtomCountDescriptor(), "S"))),
    /**
     * Atom count P, counts the number of all phosphorus atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_P(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count P",
            molecular(() -> configured(new AtomCountDescriptor(), "P"))),
    /**
     * Atom count F, counts the number of all fluorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_F(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count F",
            molecular(() -> configured(new AtomCountDescriptor(), "F"))),
    /**
     * Atom count Br, counts the number of all bromine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_BR(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count Br",
            molecular(() -> configured(new AtomCountDescriptor(), "Br"))),
    /**
     * Atom count Cl, counts the number of all chlorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_CL(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count Cl",
            molecular(() -> configured(new AtomCountDescriptor(), "Cl"))),
    /**
     * Atom count I, counts the number of all iodine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_I(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Atom Count I",
            molecular(() -> configured(new AtomCountDescriptor(), "I"))),
    /**
     * Total bond count, counts the number of all bonds in a molecule, neglecting the order.
     * Double and triple bonds are counted as one bond.
     * Bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_ALL(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Bond Count All",
            molecular(BondCountDescriptor::new)),
    /**
     * Bond count single, counts the number of single bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_SINGLE(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Bond Count Single",
            molecular(() -> configured(new BondCountDescriptor(), "s"))),
    /**
     * Bond count double, counts the number of double bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_DOUBLE(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Bond Count Double",
            molecular(() -> configured(new BondCountDescriptor(), "d"))),
    /**
     * Bond count triple, counts the number of triple bonds in a molecule.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_TRIPLE(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Bond Count Triple",
            molecular(() -> configured(new BondCountDescriptor(), "t"))),
    /**
     * H bond acceptor count, counts hydrogen bond acceptors based on a simplified PHACIR scheme.
     * It includes: Oxygen atoms with formal charge ≤ 0 (excluding: Aromatic ether oxygens and oxygens adjacent to nitrogen)
     * and nitrogen atoms with formal charge ≤ 0 (excluding: nitrogens adjacent to oxygen).
     *
     * @see HBondAcceptorCountDescriptor
     */
    H_BOND_ACCEPTOR_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "H-Bond Acceptor Count",
            molecular(HBondAcceptorCountDescriptor::new)),
    /**
     * H bond donor count, counts hydrogen bond donors based on a simplified PHACIR classification.
     * It includes: OH groups where the oxygen has a formal charge ≥ 0 and NH groups where the nitrogen has a formal charge ≥ 0.
     *
     * @see HBondDonorCountDescriptor
     */
    H_BOND_DONOR_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "H-Bond Donor Count",
            molecular(HBondDonorCountDescriptor::new)),
    /**
     * Aromatic atoms count, counts the number of aromatic atoms in a molecule.
     * Note: Requires that aromatic atoms in the molecule have already been detected and marked.
     *
     * @see AromaticAtomsCountDescriptor
     */
    AROMATIC_ATOMS_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Aromatic Atoms Count",
            molecular(AromaticAtomsCountDescriptor::new)),
    /**
     * Aromatic bonds count, counts the number of aromatic bonds in a molecule.
     * Note: Requires that aromatic bonds in the molecule have already been detected and marked.
     *
     * @see AromaticBondsCountDescriptor
     */
    AROMATIC_BONDS_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Aromatic Bonds Count",
            molecular(AromaticBondsCountDescriptor::new)),
    /**
     * Rotatable bonds count, counts the number of rotatable bonds in a molecule.
     * A rotatable bond is defined as any single non-ring bond, where atoms on both sides
     * have at least two heavy-atom neighbors. Excluding terminal bonds.
     *
     * @see RotatableBondsCountDescriptor
     */
    ROTATABLE_BONDS_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Rotatable Bonds Count",
            molecular(RotatableBondsCountDescriptor::new)),
    /**
     * Basic group count, returns the number of basic groups in a molecule.
     *
     * @see BasicGroupCountDescriptor
     */
    BASIC_GROUP_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Basic Group Count",
            molecular(() -> initialised(new BasicGroupCountDescriptor()))),
    /**
     * Acidic group count, returns the number of acidic groups in a molecule.
     *
     * @see AcidicGroupCountDescriptor
     */
    ACIDIC_GROUP_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Acidic Group Count",
            molecular(() -> initialised(new AcidicGroupCountDescriptor()))),
    /**
     * TPSA descriptor, calculates the topological polar surface area (TPSA) of a molecule.
     * TPSA is the sum of the surface areas of polar atoms (typically oxygen and nitrogen)
     * and their attached hydrogens, based on a topological approximation (2D structure only).
     *
     * @see TPSADescriptor
     */
    TPSA(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "TPSA",
            molecular(TPSADescriptor::new)),
    /**
     * Largest chain descriptor, calculates the number of atoms in the longest chain in the molecule.
     * This is a simple topological descriptor that provides a measure of molecular linearity.
     *
     * @see LargestChainDescriptor
     */
    LARGEST_CHAIN(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Largest Chain",
            molecular(LargestChainDescriptor::new)),
    /**
     * Longest aliphatic chain descriptor, calculates the number of atoms in the longest aliphatic chain.
     * This descriptor provides information about the maximum linear extent of non-aromatic
     * portions of the molecular structure, which relates to molecular shape properties.
     *
     * @see LongestAliphaticChainDescriptor
     */
    LONGEST_ALIPHATIC_CHAIN(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Longest Aliphatic Chain",
            molecular(LongestAliphaticChainDescriptor::new)),
    /**
     * BCUT descriptor, calculates Burden matrix modified eigenvalues with different weighting schemes. Returns 6 values:<br>
     * 1. BCUTw-1l, BCUTw-2l ... - nlow lowest atom weighted BCUTS<br>
     * 2. BCUTw-1h, BCUTw-2h ... - nhigh highest atom weighted BCUTS<br>
     * 3. BCUTc-1l, BCUTc-2l ... - nlow lowest partial charge weighted BCUTS<br>
     * 4. BCUTc-1h, BCUTc-2h ... - nhigh highest partial charge weighted BCUTS<br>
     * 5. BCUTp-1l, BCUTp-2l ... - nlow lowest polarizability weighted BCUTS<br>
     * 6. BCUTp-1h, BCUTp-2h ... - nhigh highest polarizability weighted BCUTS<br>
     * Note: No array for one parameter is returned, just the highest and lowest numbers. (Default Parameters: nhigh = 1 and nlow = 1)
     *
     * @see BCUTDescriptor
     */
    BCUT(false, false, false, false, false, 6, DescriptorComponentType.CONTINUOUS, "BCUT",
            molecular(() -> configured(new BCUTDescriptor(), 1, 1, false))),
    /**
     * Bond polarizability descriptor.
     * The BPolDescriptor calculates the bond polarizability of a molecule.
     * Bond polarizability is a simple sum of polarizability contributions from all bonds, based on bond types and involved atoms.
     * It provides a rough estimate of how easily the electron cloud in a molecule can be distorted,
     * which relates to intermolecular interactions, polarizability and refractive behavior.
     *
     * @see BPolDescriptor
     */
    B_POL(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "BPol",
            molecular(BPolDescriptor::new)),
    /**
     * Rule of five descriptor, calculates the number of failures of Lipinski's Rule of Five.
     * The descriptor returns the number of violations (0-4).
     *
     * @see RuleOfFiveDescriptor
     */
    RULE_OF_FIVE(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Rule of Five",
            molecular(RuleOfFiveDescriptor::new)),
    /**
     * FMF (Framework Match Fraction) descriptor, calculates the ratio of heavy atoms in
     * the framework to the total number of heavy atoms in the molecule.
     * This provides an indication of the proportion of the molecule that is part of the
     * scaffold or core structure, versus the proportion that is in side chains.
     *
     * @see FMFDescriptor
     */
    FMF(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "FMF",
            molecular(FMFDescriptor::new)),
    /**
     * Fractional C SP3 descriptor, characterizes the non-flatness of a molecule by calculating
     * the fraction of sp3 hybridized carbon atoms over the total carbon count.
     * This provides information about the three-dimensionality and complexity
     * of a molecule, which relates to drug-likeness properties.
     *
     * @see FractionalCSP3Descriptor
     */
    FRACTIONAL_CSP3(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Fractional CSP3",
            molecular(FractionalCSP3Descriptor::new)),
    /**
     * Hybridization ratio descriptor, calculates the ratio of sp3 carbons to sp2 carbons.
     * This provides valuable information about the three-dimensionality and flatness
     * of a molecule, which can be useful for predicting drug-like properties and
     * comparing structural characteristics.
     *
     * @see HybridizationRatioDescriptor
     */
    HYBRIDIZATION_RATIO(true, true, false, true, false, 1, DescriptorComponentType.CONTINUOUS, "Hybridization Ratio",
            molecular(HybridizationRatioDescriptor::new)),
    /**
     * Kappa shape indices descriptor, calculates Kier and Hall kappa molecular shape indices.
     * These indices compare the molecular graph with minimal and maximal molecular graphs. Returns 3 values:<br>
     * 1. Kier1 - First kappa shape index<br>
     * 2. Kier2 - Second kappa shape index<br>
     * 3. Kier3 - Third kappa shape index<br>
     * Note: Hydrogens are ignored in the calculation.
     *
     * @see KappaShapeIndicesDescriptor
     */
    KAPPA_SHAPE_INDICES(false, true, false, true, false, 3, DescriptorComponentType.CONTINUOUS, "Kappa Shape Indices",
            molecular(KappaShapeIndicesDescriptor::new)),
    /**
     * Petitjean number descriptor, calculates an index characterizing molecular graph topology.
     * This topological descriptor is based on the calculation of the graph eccentricity
     * and provides information about the molecular shape and branching pattern.
     *
     * @see PetitjeanNumberDescriptor
     */
    PETITJEAN_NUMBER(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Petitjean Number",
            molecular(PetitjeanNumberDescriptor::new)),
    /**
     * Spiro atom count descriptor, calculates the number of spiro atoms in a molecule.
     *
     * @see SpiroAtomCountDescriptor
     */
    SPIRO_ATOM_COUNT(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Spiro Atom Count",
            molecular(SpiroAtomCountDescriptor::new)),
    /**
     * VAdjMa descriptor, calculates the Vertex adjacency information (magnitude).
     * This is calculated as 1 + log2 m, where m is the number of heavy-heavy bonds.
     * If m is zero, then zero is returned.
     * This descriptor characterizes molecular complexity in terms of edge connectivity.
     *
     * @see VAdjMaDescriptor
     */
    V_ADJ_MAT(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "VAdjMa",
            molecular(VAdjMaDescriptor::new)),
    /**
     * Weighted path descriptor, evaluates the weighted path descriptors for a molecule.
     * Returns 5 values:<br>
     * 1. WTPT1 - molecular ID<br>
     * 2. WTPT2 - molecular ID / number of atoms<br>
     * 3. WTPT3 - sum of path lengths starting from heteroatoms<br>
     * 4. WTPT4 - sum of path lengths starting from oxygens<br>
     * 5. WTPT5 - sum of path lengths starting from nitrogens<br>
     * <p>
     * Note: This descriptor computes all paths which is an NP-hard problem, do not use it for complex molecules.
     * @see WeightedPathDescriptor
     */
    WEIGHTED_PATH(false, true, false, false, false, 5, DescriptorComponentType.CONTINUOUS, "Weighted Path",
            molecular(WeightedPathDescriptor::new)),
    /**
     * Zagreb index descriptor, calculates the Zagreb index of a molecule.
     * The Zagreb index is the sum of the squares of atom degrees over all heavy atoms,
     * which provides information about the molecular complexity and topological structure.
     *
     * @see ZagrebIndexDescriptor
     */
    ZAGREB_INDEX(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Zagreb Index",
            molecular(ZagrebIndexDescriptor::new)),
    /**
     * CarbonTypes descriptor, calculates the frequency of occurrence of 9 different types of carbon atoms. Returns 9 values:<br>
     * 1. C1SP1 - triply bound carbon bound to one other carbon<br>
     * 2. C2SP1 - triply bound carbon bound to two other carbons<br>
     * 3. C1SP2 - doubly bound carbon bound to one other carbon<br>
     * 4. C2SP2 - doubly bound carbon bound to two other carbons<br>
     * 5. C3SP2 - doubly bound carbon bound to three other carbons<br>
     * 6. C1SP3 - singly bound carbon bound to one other carbon<br>
     * 7. C2SP3 - singly bound carbon bound to two other carbons<br>
     * 8. C3SP3 - singly bound carbon bound to three other carbons<br>
     * 9. C4SP3 - singly bound carbon bound to four other carbons
     *
     * @see CarbonTypesDescriptor
     */
    CARBON_TYPES(true, true, false, false, false, 9, DescriptorComponentType.COUNT, "Carbon Types",
            molecular(CarbonTypesDescriptor::new)),
    /**
     * ALogP descriptor, calculates Ghose-Crippen LogP values, molar refractivity values
     * and ALogP squared values. Returns 3 values:<br>
     * 1. ALogP (logP value) is the Ghose-Crippen octanol-water partition coefficient.<br>
     * 2. ALogP² is the squared ALogP value.<br>
     * 3. Molar Refractivity (MR) measures the volume occupied by an atom or group of atoms.<br>
     *
     * @see ALOGPDescriptor
     */
    A_LOG_P(true, true, false, true, false, 3, DescriptorComponentType.CONTINUOUS, "ALogP",
            molecular(ALOGPDescriptor::new)),
    /**
     * XLogP descriptor, calculates logP based on the atom-type method called XLogP.
     * Requires all hydrogen's to be explicit.
     *
     * @see XLogPDescriptor
     */
    X_LOG_P(true, true, false, true, false, 1, DescriptorComponentType.CONTINUOUS, "XLogP",
            molecular(XLogPDescriptor::new)),
    /**
     * JPlogP descriptor, calculates the octanol-water partition coefficient based on an atom contribution model.
     *
     * @see JPlogPDescriptor
     */
    JP_LOG_P(true, false, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "JPlogP",
            molecular(JPlogPDescriptor::new)),
    /**
     * Mannhold LogP descriptor, calculates the octanol-water partition coefficient (logP) using the Mannhold method.
     * LogP describes the hydrophilicity or lipophilicity of a compound and is crucial for
     * predicting solubility, permeability, and bioavailability.
     *
     * @see MannholdLogPDescriptor
     */
    MANNHOLD_LOGP(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Mannhold LogP",
            molecular(MannholdLogPDescriptor::new)),
    /**
     * APol descriptor, calculates the sum of the atomic polarizabilities (including implicit hydrogens).
     *
     * @see APolDescriptor
     */
    A_POL(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "APol",
            molecular(APolDescriptor::new)),
    /**
     * Autocorrelation charge descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to charge distribution across the molecular structure.
     * This descriptor correlates atomic partial charges along the molecular topology
     * to characterize charge-related structural patterns in the molecule.
     * Returns 5 values representing charge autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorCharge
     */
    AUTOCORRELATION_CHARGE(true, false, false, false, false, 5, DescriptorComponentType.CONTINUOUS, "Autocorrelation Charge",
            molecular(AutocorrelationDescriptorCharge::new)),
    /**
     * Autocorrelation mass descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to atomic mass distribution across the molecular structure.
     * This descriptor correlates atomic masses along the molecular topology
     * to characterize mass-related structural patterns in the molecule.
     * Returns 5 values representing mass autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorMass
     */
    AUTOCORRELATION_MASS(true, true, false, false, false, 5, DescriptorComponentType.CONTINUOUS, "Autocorrelation Mass",
            molecular(AutocorrelationDescriptorMass::new)),
    /**
     * Autocorrelation polarizability descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to polarizability distribution across the molecular structure.
     * This descriptor correlates atomic polarizabilities along the molecular topology
     * to characterize polarizability-related structural patterns in the molecule.
     * Returns 5 values representing polarizability autocorrelation at different topological distances.
     * NOTE: Method is not validated in the CDK so not validated in this implementation as well
     *
     * @see AutocorrelationDescriptorPolarizability
     */
    AUTOCORRELATION_POLARIZABILITY(true, true, false, false, false, 5, DescriptorComponentType.CONTINUOUS, "Autocorrelation Polarizability",
            molecular(AutocorrelationDescriptorPolarizability::new)),
    /**
     * Fragment complexity descriptor, calculates the complexity of a molecular system.
     * The complexity is defined as [Nilakantan, R. et al. Journal of chemical information and modeling. 2006. 46]:
     * C = abs(B^2 - A^2 + A) + H/100
     * where:
     * (C = complexity,
     * A = number of non-hydrogen atoms,
     * B = number of bonds,
     * H = number of heteroatoms,).
     * This provides a measure of structural complexity that correlates with synthetic accessibility.
     *
     * @see FragmentComplexityDescriptor
     */
    FRAGMENT_COMPLEXITY(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Fragment Complexity",
            molecular(FragmentComplexityDescriptor::new)),
    /**
     * Chi chain descriptor, calculates the Kier + Hall chi chain indices of orders 3 through 7.
     * These values characterize a molecular graph based on its chain subgraphs.
     * Returns 10 values:<br>
     * 1. SCH-3 - Simple chain, order 3<br>
     * 2. SCH-4 - Simple chain, order 4<br>
     * 3. SCH-5 - Simple chain, order 5<br>
     * 4. SCH-6 - Simple chain, order 6<br>
     * 5. SCH-7 - Simple chain, order 7<br>
     * 6. VCH-3 - Valence chain, order 3<br>
     * 7. VCH-4 - Valence chain, order 4<br>
     * 8. VCH-5 - Valence chain, order 5<br>
     * 9. VCH-6 - Valence chain, order 6<br>
     * 10. VCH-7 - Valence chain, order 7
     *
     * @see ChiChainDescriptor
     */
    CHI_CHAIN(false, true, false, false, false, 10, DescriptorComponentType.CONTINUOUS, "Chi Chain",
            molecular(ChiChainDescriptor::new)),
    /**
     * Chi cluster descriptor, calculates Kier + Hall chi cluster indices of orders 3 through 6.
     * These values characterize a molecular graph based on its cluster subgraphs.
     * Returns 8 values:<br>
     * 1. SC-3 - Simple cluster, order 3<br>
     * 2. SC-4 - Simple cluster, order 4<br>
     * 3. SC-5 - Simple cluster, order 5<br>
     * 4. SC-6 - Simple cluster, order 6<br>
     * 5. VC-3 - Valence cluster, order 3<br>
     * 6. VC-4 - Valence cluster, order 4<br>
     * 7. VC-5 - Valence cluster, order 5<br>
     * 8. VC-6 - Valence cluster, order 6
     *
     * @see ChiClusterDescriptor
     */
    CHI_CLUSTER(false, true, false, false, false, 8, DescriptorComponentType.CONTINUOUS, "Chi Cluster",
            molecular(ChiClusterDescriptor::new)),
    /**
     * Chi path cluster descriptor, calculates Kier + Hall chi path cluster indices of orders 4 through 6.
     * These values characterize a molecular graph based on its path cluster subgraphs.
     * Returns 6 values:<br>
     * 1. SPC-4 - Simple path cluster, order 4<br>
     * 2. SPC-5 - Simple path cluster, order 5<br>
     * 3. SPC-6 - Simple path cluster, order 6<br>
     * 4. VPC-4 - Valence path cluster, order 4<br>
     * 5. VPC-5 - Valence path cluster, order 5<br>
     * 6. VPC-6 - Valence path cluster, order 6
     *
     * @see ChiPathClusterDescriptor
     */
    CHI_PATH_CLUSTER(false, true, false, false, false, 6, DescriptorComponentType.CONTINUOUS, "Chi Path Cluster",
            molecular(ChiPathClusterDescriptor::new)),
    /**
     * Chi path descriptor, calculates Kier + Hall chi path indices of orders 0 through 7.
     * These values characterize a molecular graph based on its path subgraphs.
     * Returns 16 values:<br>
     * 1.  SP-0 - Simple path, order 0<br>
     * 2.  SP-1 - Simple path, order 1<br>
     * 3.  SP-2 - Simple path, order 2<br>
     * 4.  SP-3 - Simple path, order 3<br>
     * 5.  SP-4 - Simple path, order 4<br>
     * 6.  SP-5 - Simple path, order 5<br>
     * 7.  SP-6 - Simple path, order 6<br>
     * 8.  SP-7 - Simple path, order 7<br>
     * 9.  VP-0 - Valence path, order 0<br>
     * 10. VP-1 - Valence path, order 1<br>
     * 11. VP-2 - Valence path, order 2<br>
     * 12. VP-3 - Valence path, order 3<br>
     * 13. VP-4 - Valence path, order 4<br>
     * 14. VP-5 - Valence path, order 5<br>
     * 15. VP-6 - Valence path, order 6<br>
     * 16. VP-7 - Valence path, order 7<br>
     *
     * @see ChiPathDescriptor
     */
    CHI_PATH(false, true, false, false, false, 16, DescriptorComponentType.CONTINUOUS, "Chi Path",
            molecular(ChiPathDescriptor::new)),
    /**
     * Fractional PSA descriptor, calculates the ratio of polar surface area to molecular weight.
     * This descriptor provides the polar surface area efficiency, which is the TPSADescriptor value divided by the
     * molecular weight, measured in square Angstroms per Dalton.
     *
     * @see FractionalPSADescriptor
     */
    FRACTIONAL_PSA(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Fractional PSA",
            molecular(FractionalPSADescriptor::new)),
    /**
     * Largest pi system descriptor, calculates the number of atoms in the largest pi system.
     * This descriptor identifies the largest conjugated pi system within a molecule and
     * returns the count of atoms participating in it.
     *
     * @see LargestPiSystemDescriptor
     */
    LARGEST_PI_SYSTEM(true, true, false, false, false, 1, DescriptorComponentType.COUNT, "Largest Pi System",
            molecular(() -> configured(new LargestPiSystemDescriptor(), false))),
    /**
     * Descriptor that calculates small ring information.
     * Returns 11 values:<br>
     * 1. nSmallRings - total number of small rings (of size 3 through 9)<br>
     * 2. nAromRings - total number of small aromatic rings<br>
     * 3. nRingBlocks - total number of distinct ring blocks<br>
     * 4. nAromBlocks - total number of aromatically connected components<br>
     * 5. nRings3 - total number of 3-membered rings<br>
     * 6. nRings4 - total number of 4-membered rings<br>
     * 7. nRings5 - total number of 5-membered rings<br>
     * 8. nRings6 - total number of 6-membered rings<br>
     * 9. nRings7 - total number of 7-membered rings<br>
     * 10. nRings8 - total number of 8-membered rings<br>
     * 11. nRings9 - total number of 9-membered rings<br>
     *
     * @see SmallRingDescriptor
     */
    SMALL_RING(true, true, false, false, false, 11, DescriptorComponentType.COUNT, "Small Ring",
            molecular(SmallRingDescriptor::new)),
    /**
     * Amino acid count descriptor, calculates the number of each amino acid in a molecule.
     * Returns 20 values, one for each of the 20 standard amino acids:
     * Alanine, Arginine, Asparagine, Aspartic acid, Cysteine, Glutamic acid, Glutamine,
     * Glycine, Histidine, Isoleucine, Leucine, Lysine, Methionine, Phenylalanine,
     * Proline, Serine, Threonine, Tryptophan, Tyrosine, and Valine.
     * This descriptor helps identify and quantify amino acid composition in peptides and proteins.
     *
     * @see AminoAcidCountDescriptor
     */
    AMINO_ACID_COUNT(false, true, false, false, false, 20, DescriptorComponentType.COUNT, "Amino Acid Count",
            molecular("AA.", AminoAcidCountDescriptor::new)),
    /**
     * Kier-Hall SMARTS descriptor that calculates counts of functional groups and substructures
     * based on the Kier and Hall SMARTS patterns, used for QSAR modeling and molecular characterization.
     * Note: This descriptor provides 79 values representing different molecular fragments.
     *
     * @see KierHallSmartsDescriptor
     */
    KIER_HALL_SMARTS(true, true, false, false, false, 79, DescriptorComponentType.COUNT, "Kier Hall SMARTS",
            molecular(KierHallSmartsDescriptor::new)),
    /**
     * Eccentric connectivity index descriptor, calculates a topological descriptor that combines
     * distance and adjacency information.
     * It is defined as the sum of the products of eccentricity and vertex degree for each atom.
     * This index provides information about the distribution of atoms in the molecular structure
     * and helps characterize molecular complexity, branching, and overall shape.
     *
     * @see EccentricConnectivityIndexDescriptor
     */
    ECCENTRIC_CONNECTIVITY_INDEX(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "Eccentric Connectivity Index",
            molecular(EccentricConnectivityIndexDescriptor::new)),
    /**
     * MDE descriptor, calculates molecular distance edge descriptors for carbon, oxygen and nitrogen atoms.
     * These descriptors encode information about the connectivity and distance of atoms of specific types
     * and hybridization states in the molecular graph. Returns 19 values representing various molecular
     * distance edge counts between different types of carbon, oxygen, and nitrogen atoms:<br>
     * MDEC-11<br>
     * MDEC-12<br>
     * MDEC-13<br>
     * MDEC-14<br>
     * MDEC-22<br>
     * MDEC-23<br>
     * MDEC-24<br>
     * MDEC-33<br>
     * MDEC-34<br>
     * MDEC-44<br>
     * MDEO-11<br>
     * MDEO-12<br>
     * MDEO-22<br>
     * MDEN-11<br>
     * MDEN-12<br>
     * MDEN-13<br>
     * MDEN-22<br>
     * MDEN-23<br>
     * MDEN-33<br>
     *
     * @see MDEDescriptor
     */
    MDE(true, true, false, false, false, 19, DescriptorComponentType.CONTINUOUS, "MDE",
            molecular(MDEDescriptor::new)),
    /**
     * VABC descriptor, calculates the volume descriptor using the van der Waals volume calculation approach.
     * This descriptor estimates molecular volume based on atom contributions, considering bond types
     * and atomic properties, providing insights into molecular size and steric properties.
     *
     * @see VABCDescriptor
     */
    VABC(true, true, false, false, false, 1, DescriptorComponentType.CONTINUOUS, "VABC",
            molecular(VABCDescriptor::new)),
    /**
     * PubChem fingerprinter, generates a 881-bit binary fingerprint based on PubChem's substructure keys.
     * This fingerprint encodes the presence or absence of specific substructural features
     * and is useful for similarity searching and chemical space analysis.
     *
     * @see PubchemFingerprinter
     */
    PUBCHEM_FINGERPRINTER(false, true, true, false, false, 881, DescriptorComponentType.BINARY, "PubChem Fingerprinter",
            fingerprint("PubChem", () -> new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()))),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_0(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter ECFP-0",
            circular(CircularFingerprinter.CLASS_ECFP0)),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_0(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter FCFP-0",
            circular(CircularFingerprinter.CLASS_FCFP0)),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_2(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter ECFP-2",
            circular(CircularFingerprinter.CLASS_ECFP2)),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_2(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter FCFP-2",
            circular(CircularFingerprinter.CLASS_FCFP2)),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_4(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter ECFP-4",
            circular(CircularFingerprinter.CLASS_ECFP4)),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_4(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter FCFP-4",
            circular(CircularFingerprinter.CLASS_FCFP4)),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_6(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter ECFP-6",
            circular(CircularFingerprinter.CLASS_ECFP6)),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_6(true, true, true, false, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
            DescriptorComponentType.BINARY, "Circular Fingerprinter FCFP-6",
            circular(CircularFingerprinter.CLASS_FCFP6)),
    /**
     * MACCS fingerprinter, generates a 166-bit binary fingerprint based on the MACCS structural keys.
     *
     * @see MACCSFingerprinter
     */
    MACCS_FINGERPRINTER(true, true, true, false, false, 166, DescriptorComponentType.BINARY, "MACCS Fingerprinter",
            fingerprint("MACCS", () -> new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()))),
    /**
     * Charged Partial Surface Area (CPSA) descriptor, calculates 29 surface and charge descriptors.
     *
     * @see CPSADescriptor
     */
    CPSA(false, true, false, false, true, 29, DescriptorComponentType.CONTINUOUS, "CPSA",
            molecular(CPSADescriptor::new)),
    /**
     * Gravitational Index descriptor, calculates 9 indices based on mass and 3D distances.
     *
     * @see GravitationalIndexDescriptor
     */
    GRAVITATIONAL_INDEX(true, true, false, false, true, 9, DescriptorComponentType.CONTINUOUS, "Gravitational Index",
            molecular(GravitationalIndexDescriptor::new)),
    /**
     * Moment of Inertia descriptor, calculates 7 components from the principal moments of inertia.
     *
     * @see MomentOfInertiaDescriptor
     */
    MOMENT_OF_INERTIA(true, true, false, false, true, 7, DescriptorComponentType.CONTINUOUS, "Moment of Inertia",
            molecular(MomentOfInertiaDescriptor::new)),
    /**
     * WHIM (Weighted Holistic Invariant Molecular) descriptor, calculates 17 directional descriptors.
     *
     * @see WHIMDescriptor
     */
    WHIM(true, false, false, false, true, 17, DescriptorComponentType.CONTINUOUS, "WHIM",
            molecular(WHIMDescriptor::new)),
    /**
     * Length Over Breadth descriptor, calculates maximum and minimum length-to-breadth ratios (LOBMAX, LOBMIN).
     *
     * @see LengthOverBreadthDescriptor
     */
    LENGTH_OVER_BREADTH(true, true, false, false, true, 2, DescriptorComponentType.CONTINUOUS, "Length Over Breadth",
            molecular(LengthOverBreadthDescriptor::new)),
    /**
     * Petitjean Shape Index descriptor, calculates topological and geometric shape indices (topoShape, geomShape).
     *
     * @see PetitjeanShapeIndexDescriptor
     */
    PETITJEAN_SHAPE_INDEX(true, true, false, false, true, 2, DescriptorComponentType.CONTINUOUS, "Petitjean Shape Index",
            molecular(PetitjeanShapeIndexDescriptor::new));

    // Add new descriptor information here!

    /**
     * Indicates whether this descriptor is quickly calculable.
     * Fast descriptors have lower computational complexity and can be calculated efficiently,
     * while slow descriptors require more intensive calculations (e.g., CHI descriptors).
     */
    private final boolean isFast;

    /**
     * Indicates whether this descriptor is safe to calculate.
     * Safe descriptors are reliable and validated, while unsafe descriptors may produce
     * inconsistent results or NaN values under certain conditions.
     */
    private final boolean isSafe;

    /**
     * Indicates whether the enum entry is a fingerprint or not.
     */
    private final boolean isFingerprint;

    /**
     * Indicates whether this descriptor needs explicit hydrogens.
     */
    private final boolean needsExplicitHydrogens;

    /**
     * Indicates whether this descriptor requires 3D coordinates.
     */
    private final boolean requires3DCoordinates;

    /**
     * The number of components calculated by this descriptor (for circular fingerprints the default size).
     */
    private final int descriptorComponentNumber;

    /**
     * The type of the values of all components of this descriptor (for circular fingerprints the default type).
     */
    private final DescriptorComponentType componentType;

    /**
     * The human-readable name of this descriptor for output purposes.
     */
    private final String name;

    /**
     * The recipe for the calculation of this descriptor, used by {@link DescriptorCalculator}. Note: No CDK
     * descriptor or fingerprinter instance is created when the enum is loaded, only when the factory is called.
     */
    private final DescriptorCalculationFactory calculationFactory;

    /**
     * Constructs a Descriptor with the given field values.
     *
     * @param isFast true if the descriptor is quickly calculable, false if it requires
     *               more intensive computation
     * @param isSafe true if the descriptor is safe and reliable, false if it may produce
     *               inconsistent results or NaN values
     * @param isFingerprint true if the "descriptor" is actually a fingerprint, false otherwise
     * @param needsExplicitHydrogens true if the descriptor requires explicit hydrogens, false otherwise
     * @param requires3DCoordinates true if the descriptor requires 3D coordinates, false otherwise
     * @param descriptorComponentNumber the number of components calculated by this descriptor
     * @param componentType the type of the values of all components of this descriptor
     * @param name the human-readable name of this descriptor for output purposes
     * @param calculationFactory the recipe for the calculation of this descriptor
     */
    Descriptor(boolean isFast, boolean isSafe, boolean isFingerprint, boolean needsExplicitHydrogens, boolean requires3DCoordinates, int descriptorComponentNumber, DescriptorComponentType componentType, String name, DescriptorCalculationFactory calculationFactory) {
        this.isFast = isFast;
        this.isSafe = isSafe;
        this.isFingerprint = isFingerprint;
        this.needsExplicitHydrogens = needsExplicitHydrogens;
        this.requires3DCoordinates = requires3DCoordinates;
        this.descriptorComponentNumber = descriptorComponentNumber;
        this.componentType = componentType;
        this.name = name;
        this.calculationFactory = calculationFactory;
    }

    /**
     * Returns whether this descriptor can be calculated quickly. This classification is based on a performance benchmark
     * of 10,000 molecules: descriptors computed in under 1 second for 10,000 molecules are considered fast.
     *
     * @return true if the descriptor has low computational complexity and can be
     *         calculated efficiently, false if it requires intensive computation
     */
    public boolean isFast() {
        return this.isFast;
    }

    /**
     * Returns whether this descriptor is safe to calculate.
     *
     * @return true if the descriptor is reliable and validated, false if it may
     *         produce inconsistent results or NaN values under certain conditions
     */
    public boolean isSafe() {
        return this.isSafe;
    }

    /**
     * Returns whether this descriptor calculates a fingerprint.
     *
     * @return true if the enum entry calculates a fingerprint, false otherwise
     */
    public boolean isFingerprint() {
        return this.isFingerprint;
    }

    /**
     * Returns whether this descriptor needs explicit hydrogens.
     *
     * @return true if the descriptor requires explicit hydrogens, false otherwise
     */
    public boolean needsExplicitHydrogens() {
        return this.needsExplicitHydrogens;
    }

    /**
     * Returns whether this descriptor requires 3D coordinates.
     *
     * @return true if the descriptor requires 3D coordinates, false otherwise
     */
    public boolean requires3DCoordinates() {
        return this.requires3DCoordinates;
    }

    /**
     * Returns the number of components calculated by this descriptor.
     * Note: For circular fingerprints (CIRCULAR_FINGERPRINTER_*), this is the default size
     * (1024); the actual number of components depends on the configuration of
     * the {@link DescriptorCalculator}, see {@link DescriptorCalculator#getNumberOfComponents(Descriptor)}.
     *
     * @return the number of components calculated by this descriptor
     */
    public int getDescriptorComponentNumber() {
        return this.descriptorComponentNumber;
    }

    /**
     * Returns the type of the values of all components of this descriptor, e.g. for choosing a suitable preprocessing
     * before machine learning.
     * Note: For circular fingerprints (CIRCULAR_FINGERPRINTER_*), this is the default type
     * ({@link DescriptorComponentType#BINARY}); the actual type depends on the configuration of the
     * {@link DescriptorCalculator}, see {@link DescriptorCalculator#getComponentType(Descriptor)}.
     *
     * @return the type of the values of all components of this descriptor
     */
    public DescriptorComponentType getComponentType() {
        return this.componentType;
    }

    /**
     * Returns the human-readable name of this descriptor.
     *
     * @return the human-readable name of this descriptor
     */
    public String getName() {
        return this.name;
    }

    /**
     * Returns the recipe for the calculation of this descriptor. Used by {@link DescriptorCalculator} only.
     *
     * @return the calculation factory of this descriptor
     */
    DescriptorCalculationFactory getCalculationFactory() {
        return this.calculationFactory;
    }

    /**
     * Creates a (configured) CDK molecular descriptor instance. Unlike a Supplier, it may throw a CDKException,
     * because some CDK descriptor constructors and setParameters() declare it.
     */
    @FunctionalInterface
    private interface CdkDescriptorSupplier {
        /**
         * Creates the CDK descriptor instance.
         *
         * @return the CDK descriptor instance
         * @throws CDKException if the CDK descriptor cannot be created or configured
         */
        IMolecularDescriptor get() throws CDKException;
    }

    /**
     * Returns the recipe for a CDK molecular descriptor: one CDK descriptor instance that is shared by all
     * calculation threads.
     *
     * @param cdkDescriptorSupplier creates the (configured) CDK descriptor instance
     * @return the calculation factory
     */
    private static DescriptorCalculationFactory molecular(CdkDescriptorSupplier cdkDescriptorSupplier) {
        return Descriptor.molecular("", cdkDescriptorSupplier);
    }

    /**
     * Returns the recipe for a CDK molecular descriptor whose component names get a prefix, because its CDK
     * descriptor names are not unique among all descriptors (e.g. "AA." for the amino acid counts, whose CDK names
     * like "nC" for cysteine collide with the atom counts like "nC" for carbon).
     *
     * @param namePrefix prefix of the component names
     * @param cdkDescriptorSupplier creates the (configured) CDK descriptor instance
     * @return the calculation factory
     */
    private static DescriptorCalculationFactory molecular(String namePrefix, CdkDescriptorSupplier cdkDescriptorSupplier) {
        return settings -> {
            try {
                return new MolecularCalculation(cdkDescriptorSupplier.get(), namePrefix);
            } catch (CDKException exception) {
                throw new IllegalStateException("Failed to create the CDK descriptor, this should never happen.", exception);
            }
        };
    }

    /**
     * Returns the recipe for a CDK fingerprinter (bit fingerprint): a pool of fingerprinter instances, because
     * fingerprinters are not thread-safe.
     *
     * @param namePrefix prefix of the component names, e.g. "PubChem" for "PubChem_0" ... "PubChem_880"
     * @param fingerprinterFactory creates a fingerprinter instance
     * @return the calculation factory
     */
    private static DescriptorCalculationFactory fingerprint(String namePrefix, Supplier<IFingerprinter> fingerprinterFactory) {
        return settings -> new FingerprintCalculation(fingerprinterFactory, settings.getFingerprintPoolSize(), false, namePrefix);
    }

    /**
     * Returns the recipe for a CDK circular fingerprint: a pool of CircularFingerprinter instances with the
     * configured circular fingerprint size, calculated as bit or folded count fingerprint depending on the
     * configuration (see {@link DescriptorCalculator#DescriptorCalculator(int, int, boolean)}).
     *
     * @param circularFingerprinterClass the class of the circular fingerprint, e.g. CircularFingerprinter.CLASS_ECFP4
     * @return the calculation factory
     */
    private static DescriptorCalculationFactory circular(int circularFingerprinterClass) {
        return settings -> new FingerprintCalculation(
                () -> new CircularFingerprinter(circularFingerprinterClass, settings.getCircularFingerprintSize()),
                settings.getFingerprintPoolSize(),
                settings.isCircularFingerprintCounts(),
                Descriptor.circularFingerprintName(circularFingerprinterClass));
    }

    /**
     * Returns the name of a circular fingerprint class, used as prefix of its component names, e.g. "ECFP4".
     *
     * @param circularFingerprinterClass the class of the circular fingerprint, e.g. CircularFingerprinter.CLASS_ECFP4
     * @return the name of the class
     * @throws IllegalArgumentException if the class is unknown
     */
    private static String circularFingerprintName(int circularFingerprinterClass) {
        switch (circularFingerprinterClass) {
            case CircularFingerprinter.CLASS_ECFP0: return "ECFP0";
            case CircularFingerprinter.CLASS_ECFP2: return "ECFP2";
            case CircularFingerprinter.CLASS_ECFP4: return "ECFP4";
            case CircularFingerprinter.CLASS_ECFP6: return "ECFP6";
            case CircularFingerprinter.CLASS_FCFP0: return "FCFP0";
            case CircularFingerprinter.CLASS_FCFP2: return "FCFP2";
            case CircularFingerprinter.CLASS_FCFP4: return "FCFP4";
            case CircularFingerprinter.CLASS_FCFP6: return "FCFP6";
            default: throw new IllegalArgumentException("Unknown circular fingerprint class: " + circularFingerprinterClass);
        }
    }

    /**
     * Sets the parameters of a CDK molecular descriptor.
     *
     * @param cdkDescriptor the CDK descriptor to configure (IS CHANGED)
     * @param parameters the parameters, see the documentation of the respective CDK descriptor
     * @return the configured CDK descriptor
     * @throws CDKException if the parameters are invalid, which should never happen
     */
    private static IMolecularDescriptor configured(IMolecularDescriptor cdkDescriptor, Object... parameters) throws CDKException {
        cdkDescriptor.setParameters(parameters);
        return cdkDescriptor;
    }

    /**
     * Initialises a CDK molecular descriptor with the silent chem object builder.
     *
     * @param cdkDescriptor the CDK descriptor to initialise (IS CHANGED)
     * @return the initialised CDK descriptor
     */
    private static IMolecularDescriptor initialised(IMolecularDescriptor cdkDescriptor) {
        cdkDescriptor.initialise(SilentChemObjectBuilder.getInstance());
        return cdkDescriptor;
    }

    /**
     * Default size used for all circular fingerprint "descriptors".
     * Note: This constant must remain static final because it is required during
     * the static initialization of the enum constants above. Enum constants are
     * instantiated before static variables are initialized, so using a non-final
     * variable would result in a default value of 0.
     */
    static final int CIRCULAR_FINGERPRINT_DEFAULT_SIZE = 1024;

    /**
     * Returns all available descriptors.
     *
     * @return All available descriptors
     */
    public static Descriptor[] getAllDescriptors(){
        return Descriptor.values();
    }

    /**
     * Returns all available fingerprint descriptors.
     * Note: Descriptors that require 3D coordinates are excluded.
     *
     * @return All available fingerprint descriptors
     */
    public static Descriptor[] getAllFingerprints() {
        return Arrays.stream(Descriptor.values())
                .filter(d -> d.isFingerprint() && !d.requires3DCoordinates())
                .toArray(Descriptor[]::new);
    }

    /**
     * Returns specified available descriptors.
     * Note: Fast, safe, 2D/1D, and non-fingerprint descriptors are automatically included by default.
     *
     * @param isSlowlyCalculableDescriptorInclusion True: Slowly calculable descriptors are included in the result, false: Otherwise.
     * @param isUnsafeDescriptorInclusion True: Unsafe descriptors are included in the result, false: Unsafe descriptors
     *                                     are excluded from the result, this does not mean no NaN's can be produced.
     * @param isFingerprintAsDescriptorInclusion True: Fingerprint descriptors are included in the result, false: Fingerprint descriptors are excluded.
     * @param is3dDescriptorInclusion True: 3D descriptors are included in the result, false: 3D descriptors are excluded.
     * @return Specified descriptors
     */
    public static Descriptor[] getSpecifiedDescriptors(
            boolean isSlowlyCalculableDescriptorInclusion,
            boolean isUnsafeDescriptorInclusion,
            boolean isFingerprintAsDescriptorInclusion,
            boolean is3dDescriptorInclusion
    ) {
        // Initialize ArrayList with maximum possible capacity to avoid internal resizing during element addition
        List<Descriptor> result = new ArrayList<>(Descriptor.values().length);

        for (Descriptor descriptor : Descriptor.values()) {
            boolean includeDescriptor = true;

            // Exclude slow descriptors if not requested
            if (!descriptor.isFast() && !isSlowlyCalculableDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude unsafe descriptors if not requested
            if (!descriptor.isSafe() && !isUnsafeDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude fingerprint descriptors if not requested
            if (descriptor.isFingerprint() && !isFingerprintAsDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude 3D descriptors if not requested
            if (descriptor.requires3DCoordinates() && !is3dDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Only fast, safe, 2D/1D, and non-fingerprint descriptors remain if all flags are false
            if (includeDescriptor) {
                result.add(descriptor);
            }
        }

        return result.toArray(new Descriptor[0]);
    }

    /**
     * Returns sum of number of calculated components of an array of defined descriptors.
     * Note: The default size is used for circular fingerprints; for the number of matrix columns of a
     * {@link DescriptorCalculator} with another circular fingerprint size, use
     * {@link DescriptorCalculator#getNumberOfComponents(Descriptor[])}.
     *
     * @param descriptors Array of descriptors (IS NOT CHANGED); if it is empty, 0 is returned;
     *                     may not be null (the array or one of its elements)
     * @return Sum of number of calculated components of array of descriptors
     */
    public static int getNumberOfComponents(
        Descriptor[] descriptors
    ) {
        // Checks
        if (descriptors == null) {
            throw new NullPointerException("Descriptor.getNumberOfComponents: descriptor is null.");
        }
        if (descriptors.length == 0) {
            return 0;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor == null) {
                throw new NullPointerException("Descriptor.getNumberOfComponents: At least one descriptor in descriptors is null.");
            }
        }

        int totalNumberOfComponents = 0;
        for (Descriptor descriptor : descriptors) {
            totalNumberOfComponents += descriptor.getDescriptorComponentNumber();
        }
        return totalNumberOfComponents;
    }


    /**
     * Returns the descriptor and component index for a given position in a descriptor array.
     * Note: The default size is used for circular fingerprints; for a {@link DescriptorCalculator} with another
     * circular fingerprint size, use {@link DescriptorCalculator#getDescriptorAndComponentIndex(Descriptor[], int)}.
     * <p>
     * Example: Given descriptors [MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)]:
     * <ul>
     *   <li>Index 0 → MOLECULAR_WEIGHT, component 0</li>
     *   <li>Index 1 → BCUT, component 0</li>
     *   <li>...</li>
     *   <li>Index 6 → BCUT, component 5</li>
     *   <li>Index 7 → WIENER_NUMBER, component 0</li>
     *   <li>Index 8 → WIENER_NUMBER, component 1</li>
     * </ul>
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A two-element int array: [descriptor index in descriptors, component index within that descriptor]
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static int[] getDescriptorAndComponentIndex(Descriptor[] descriptors, int anIndex)
            throws IllegalArgumentException {
        // Checks
        if (descriptors == null) {
            throw new NullPointerException(
                    "Descriptor.getDescriptorAndComponentIndex: descriptors must not be null."
            );
        }
        if (descriptors.length == 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: descriptors must not be null or empty."
            );
        }
        if (anIndex < 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: anIndex must be >= 0."
            );
        }
        if (anIndex >= Descriptor.getNumberOfComponents(descriptors)) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: anIndex (" + anIndex +
                            ") exceeds total component count (" + Descriptor.getNumberOfComponents(descriptors) + ")."
            );
        }

        // Calculate cumulative component counts and find the target descriptor
        int cumulativeCount = 0;
        int[] result = new int[2];
        for (int i = 0; i < descriptors.length; i++) {
            if (descriptors[i] == null) {
                throw new NullPointerException(
                        "Descriptor.getDescriptorAndComponentIndex: descriptors contains null element at index " + i
                );
            }

            int componentCount = descriptors[i].getDescriptorComponentNumber();
            int nextCumulativeCount = cumulativeCount + componentCount;

            // Check if anIndex falls within this descriptor's range
            if (anIndex < nextCumulativeCount) {
                result[0] = i;
                result[1] = anIndex - cumulativeCount;
                break;
            }

            cumulativeCount = nextCumulativeCount;
        }
        return result;
    }

    /**
     * Returns the descriptor and its component name for a given position in a descriptor array.
     * Note: The default size is used for circular fingerprints; for a {@link DescriptorCalculator} with another
     * circular fingerprint size, use {@link DescriptorCalculator#getDescriptorAndComponentInfo(Descriptor[], int)},
     * which additionally returns the column name of the component.
     * <p>
     * This is a convenience method that extends {@link #getDescriptorAndComponentIndex} by also
     * providing human-readable names for the descriptor and component.
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A String array: [descriptor name, component index as string]
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static String[] getDescriptorAndComponentInfo(Descriptor[] descriptors, int anIndex)
            throws IllegalArgumentException {
        int[] indices = Descriptor.getDescriptorAndComponentIndex(descriptors, anIndex);
        Descriptor descriptor = descriptors[indices[0]];

        return new String[] {
                descriptor.getName(),
                String.valueOf(indices[1])
        };
    }
}
