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
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.HydrogenState;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

/**
 * Descriptor related calculations based on the CDK for the enrichment of data vectors.
 * <p>
 *     The available descriptors and fingerprints, including their characteristics (e.g. the number of
 *     calculated components), are defined by the {@link Descriptor} enum.
 * </p>
 * <p>
 *     There are 4 different "public boolean setDescriptorsForMolecules...()" methods with different forms of
 *     (parallelized) calculation. There is single molecule processing or batch processing for large datastructures.
 *     Single and Batch processing can be used with either an IAtomContainer array or with a String Array of SMILES codes.
 * </p>
 * <p>
 *     The four primary public calculation methods differ along two independent axes – the <i>input type</i>
 *     (pre-parsed {@link IAtomContainer} objects vs. raw SMILES strings) and the
 *     <i>parallelization strategy</i> (molecule-level parallelization vs. batch-level parallelization) – and are
 *     summarized in the table below. Every method accepts a {@code boolean isParallelCalculation} flag; passing
 *     {@code false} disables parallelization and the method runs sequentially instead.
 *     <br><br>
 *     <table border="1">
 *         <caption>Overview of the four public batch-calculation methods</caption>
 *         <tr>
 *             <th>Method</th>
 *             <th>Input type</th>
 *             <th>Parallelisation unit</th>
 *             <th>Extra parameter</th>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesByMoleculeParallelization}</td>
 *             <td>{@link IAtomContainer IAtomContainer[]}</td>
 *             <td>One thread per <b>molecule</b> (Java parallel stream over all molecule indices)</td>
 *             <td>–</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesByBatchParallelization}</td>
 *             <td>{@link IAtomContainer IAtomContainer[]}</td>
 *             <td>One thread per <b>batch</b> of molecules (parallel stream over batch indices)</td>
 *             <td>{@code batchSize} – number of molecules per batch</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesBySmilesStringParallelization}</td>
 *             <td>{@code String[]} (SMILES)</td>
 *             <td>One thread per <b>molecule</b> (Java parallel stream over all molecule indices)</td>
 *             <td>{@code electronDonationModel} – aromaticity model applied during SMILES parsing</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculeBySmilesStringsBatchParallelization}</td>
 *             <td>{@code String[]} (SMILES)</td>
 *             <td>One thread per <b>batch</b> of molecules (parallel stream over batch indices)</td>
 *             <td>{@code batchSize}, {@code electronDonationModel}</td>
 *         </tr>
 *     </table>
 *     <br>
 *     <b>Molecule-level parallelization</b> ({@code ...ByMoleculeParallelization}) dispatches every molecule as an
 *     independent work item to Java's common fork-join pool via {@link IntStream#parallel()}.
 *     This provides fine-grained load balancing and is generally the preferred choice when the molecule set is
 *     large and individual descriptor calculations vary in cost.
 *     <br><br>
 *     <b>Batch-level parallelization</b> ({@code ...ByBatchParallelization}) first divides the molecule array into
 *     fixed-size batches of {@code batchSize} molecules and then submits one batch per fork-join task. Within
 *     each batch, molecules are processed sequentially. This strategy can reduce fork-join overhead for very fast
 *     descriptors or when the overhead of spawning one task per molecule would dominate.
 *     <br><br>
 *     <b>SMILES-based methods</b> additionally parse each SMILES string into an
 *     {@link IAtomContainer} and apply aromaticity perception (using the supplied
 *     {@link ElectronDonation} model, or {@code Daylight} by default if
 *     {@code null} is passed) before descriptor calculation. Pre-parsed {@code IAtomContainer} methods skip this
 *     step; aromaticity must have been perceived on the containers beforehand (e.g. via
 *     {@link #setAromaticity(IAtomContainer, ElectronDonation)}).
 *     <br><br>
 *     All four methods share the same result convention: they return {@code true} if every descriptor value was
 *     computed without producing {@link Float#NaN}, and {@code false} if at least one NaN was encountered or if the
 *     calculation was aborted before it started, i.e. if the given descriptor or molecule array is empty or a
 *     descriptor requires 3D coordinates that are not available (always the case for SMILES input); in the latter
 *     cases, the matrix is not changed. NaN
 *     positions are additionally recorded in the caller-supplied {@code nanPositionsList} list as
 *     {@code int[]{moleculeIndex, componentIndex}} pairs, sorted by molecule index and component index. Any
 *     modifiable list can be used, also for parallel executions.
 *     <br><br>
 *     Thread safety for fingerprints is handled internally via a pool of
 *     {@link IFingerprinter} instances (one pool per fingerprint descriptor),
 *     because CDK fingerprinter objects are not thread-safe. The pool size defaults to 4 and can be adjusted
 *     with {@link #DescriptorCalculator(int, int, boolean)}.
 * <p>
 *     A calculator object holds its configuration (fingerprint pool size, circular fingerprint size and bit or count
 *     mode for circular fingerprints) and its CDK descriptor and fingerprinter instances. The configuration is fixed
 *     at construction, so one calculator can be used by any number of threads at the same time, and calculators
 *     with different configurations can be used side by side. Because creating a calculator creates all CDK
 *     descriptor and fingerprinter instances, a calculator should be reused for many calculations.
 * <p>
 *     Example usage of the DescriptorCalculator class:
 * </p>
 * <pre>{@code
 * // Create a calculator with the default configuration (or e.g. new DescriptorCalculator(8, 2048, true))
 * DescriptorCalculator calculator = new DescriptorCalculator();
 *
 * // Preprocessing
 * // A: Molecule preprocessing:
 * // 1: Create a molecule
 * IAtomContainer molecule = ...;
 *
 * // 2. IMPORTANT: Perform aromaticity perception before descriptor calculation
 * DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
 *
 * // 3. Create a molecule array
 * IAtomContainer[] molecules = new IAtomContainer[] { molecule };
 *
 * // B: String preprocessing
 * // 1: Create a String array of SMILES codes
 * String[] moleculeSmilesStrings = new String[] {"CCO", "CCC(O)O"};
 *
 * // Define parameters
 * // 1: Define descriptor array
 * Descriptor[] descriptors = new Descriptor[] {
 *     Descriptor.MOLECULAR_WEIGHT,
 *     Descriptor.TPSA,
 *     Descriptor.A_LOG_P
 * };
 * // or use a method to get specific descriptors instead, e.g.:
 * // Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(false, false, false, false);
 *
 * // 2. Determine total number of components needed (because some descriptors return multiple values, it is not equal to the number of descriptors)
 * int componentCount = calculator.getNumberOfComponents(descriptors);
 *
 * // 3. Determine total number of molecules
 * int moleculeCount = molecules.length;
 *
 * // 4. Create data matrix for results with appropriate size (molecules as rows, descriptors as columns)
 * float[][] matrix = new float[moleculeCount][componentCount];
 * float[][] smilesMatrix = new float[moleculeSmilesStrings.length][componentCount]; // for SMILES input
 *
 * // 5. Set startIndex to 0 for filling the matrix with descriptor values starting at the 0th column
 * int startIndex = 0; // if the matrix is part of a larger data structure, set startIndex accordingly
 *
 * // 6. Choose whether to use parallel calculation
 * boolean isParallelCalculation = true; // false for sequential calculation
 *
 * // 7. Create a list for keeping track of NaN positions (any modifiable list, also for parallel calculations)
 * List<int[]> nanPositions = new ArrayList<>();
 *
 * // For Batch Processing:
 * // Define batch size
 * int batchSize = 100;
 *
 * // For SMILES Processing
 * // Define anElectronDonation model for aromaticity handling
 * ElectronDonation electronDonationModel = Aromaticity.Model.Daylight;
 *
 * // Calculate descriptors (choose one of the following methods)
 *
 * // Option A: High-Performance-Parallel Descriptor Calculation with IAtomContainer array
 * boolean isSuccessful = calculator.setDescriptorsForMoleculesByMoleculeParallelization(
 *                      descriptors,
 *                      molecules,
 *                      matrix,
 *                      startIndex,
 *                      isParallelCalculation,
 *                      nanPositions
 * );
 * // isSuccessful is false if at least one NaN was produced (see nanPositions) or the calculation was aborted
 * // you can also use calculator.setDescriptorsForMoleculesByBatchParallelization(...) -> batchSize needs to be specified accordingly
 *
 * // Option B: High-Performance-Parallel Descriptor Calculation with SMILES string array
 * boolean isSuccessfulSmiles = calculator.setDescriptorsForMoleculesBySmilesStringParallelization(
 *                      descriptors,
 *                      moleculeSmilesStrings,
 *                      smilesMatrix,
 *                      startIndex,
 *                      electronDonationModel,
 *                      isParallelCalculation,
 *                      nanPositions
 * );
 * // you can also use calculator.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(...) -> batchSize needs to be specified accordingly
 * }</pre>
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
public final class DescriptorCalculator {
    /**
     * Logger of this class.
     */
    static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(DescriptorCalculator.class);

    /**
     * SMILES Parser for SMILES String batch processing.
     */
    private static final SmilesParser SMILES_PARSER = new SmilesParser(SilentChemObjectBuilder.getInstance());

    /**
     * Default pool size for fingerprinter instances, which should be sufficient for regular users.
     */
    public static final int DEFAULT_FINGERPRINT_POOL_SIZE = 4;

    /**
     * Pool size for fingerprinter instances, i.e. number of fingerprinter instances per fingerprint "descriptor".
     */
    private final int fingerprintPoolSize;

    /**
     * Size used for all circular fingerprint "descriptors", i.e. their number of components.
     */
    private final int circularFingerprintSize;

    /**
     * True: circular fingerprints are calculated as folded count fingerprints; false: as bit fingerprints.
     */
    private final boolean isCircularFingerprintCounts;

    /**
     * EnumMap that maps every descriptor enum constant to its calculation, created from the descriptor's
     * {@link DescriptorCalculationFactory} with the configuration of this calculator (fingerprint pool size, circular
     * fingerprint size and circular fingerprint counts). The map is never modified after construction.
     */
    private final EnumMap<Descriptor, DescriptorCalculation> calculations;

    /**
     * Creates a calculator with the default configuration: a fingerprint pool size of
     * {@link #DEFAULT_FINGERPRINT_POOL_SIZE}, circular fingerprints with the default size of the {@link Descriptor}
     * enum (1024), calculated as bit fingerprints.
     */
    public DescriptorCalculator() {
        this(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, false);
    }

    /**
     * Creates a calculator with the given configuration. The configuration cannot be changed afterwards, so one
     * calculator can be used by any number of threads at the same time. Creating a calculator creates all CDK
     * descriptor and fingerprinter instances, so a calculator should be reused for many calculations.
     * <p>
     * The fingerprint pool size is the number of fingerprinter instances per fingerprint "descriptor", because CDK
     * fingerprinter instances are not thread-safe. Increasing the pool size can improve performance in highly
     * parallel environments but increases memory usage.
     * <p>
     * In count mode, all circular fingerprint "descriptors" (CIRCULAR_FINGERPRINTER_*) are calculated as folded count
     * fingerprints instead of bit fingerprints. The count fingerprint of the CDK CircularFingerprinter is not folded,
     * it enumerates hash codes and their counts. To obtain a fixed number of components (the circular fingerprint
     * size), every hash code is folded in exactly the same way as for the bit fingerprint, and the counts of all hash
     * codes folded into the same component are added up. Therefore, a component of the folded count fingerprint is
     * greater than 0 if and only if the corresponding bit of the bit fingerprint is set, and the number of components
     * does not change. Note: In count mode, the components of circular fingerprints are counts (0, 1, 2, ...) instead
     * of 0 or 1, although {@link Descriptor#isFingerprint()} is still true for them. PubChem and MACCS fingerprints are
     * always bit fingerprints.
     *
     * @param aFingerprintPoolSize        number of fingerprinter instances per fingerprint "descriptor" (must be
     *                                    greater than 0); default {@link #DEFAULT_FINGERPRINT_POOL_SIZE}
     * @param aCircularFingerprintSize    size, i.e. number of components, of all circular fingerprint "descriptors"
     *                                    (must be greater than 0); default 1024
     * @param isCircularFingerprintCounts true: circular fingerprints are calculated as folded count fingerprints;
     *                                    false: as bit fingerprints (default)
     * @throws IllegalArgumentException if aFingerprintPoolSize or aCircularFingerprintSize is less than or equal to 0
     */
    public DescriptorCalculator(int aFingerprintPoolSize, int aCircularFingerprintSize, boolean isCircularFingerprintCounts) {
        if (aFingerprintPoolSize <= 0) {
            throw new IllegalArgumentException("DescriptorCalculator: aFingerprintPoolSize must be greater than 0.");
        }
        if (aCircularFingerprintSize <= 0) {
            throw new IllegalArgumentException("DescriptorCalculator: aCircularFingerprintSize must be greater than 0.");
        }
        this.fingerprintPoolSize = aFingerprintPoolSize;
        this.circularFingerprintSize = aCircularFingerprintSize;
        this.isCircularFingerprintCounts = isCircularFingerprintCounts;
        this.calculations = this.createCalculations();
    }

    /**
     * Calculates descriptor or fingerprint values for a molecule and stores them in the result vector.
     * <p>
     * This is the core calculation method. It delegates to the {@link DescriptorCalculation} of the descriptor, i.e.
     * a {@link MolecularCalculation} (one shared CDK descriptor instance) or a {@link FingerprintCalculation} (a pool
     * of fingerprinter instances, because fingerprinters are not thread-safe).
     * <p>
     * <b>Important:</b> This method does NOT perform input validation. All necessary checks
     * must be performed by the calling public methods before invoking this method.
     *
     * @param descriptor    the descriptor to calculate (IS NOT CHANGED)
     * @param atomContainer the molecule to calculate descriptors for (IS NOT CHANGED);
     *                        must have aromaticity already perceived if required by the descriptor
     * @param vector         the result vector to store calculated values (MAY BE CHANGED);
     *                        values are stored starting at startIndex
     * @param startIndex     the starting index in vector where results should be written;
     *                        subsequent values are written to startIndex + 1, startIndex + 2, etc.
     * @return Returns true if the calculation was successful and results were stored in the vector, false otherwise.
     * @throws InterruptedException if the current thread is interrupted while waiting for a pooled fingerprinter
     */
    private boolean calculate(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex) throws InterruptedException {
        try {
            return this.calculations.get(descriptor).calculate(descriptor, atomContainer, vector, startIndex);
        } catch (InterruptedException exception) {
            LOGGER.warn("Interrupted while waiting for fingerprinter instance from pool. This should not happen.", exception);
            throw exception;
        } catch (Exception exception) {
            LOGGER.error("Descriptor ", descriptor.name(), " failed to calculate: ", exception.getMessage());
            return false;
        }
    }

    /**
     * Creates the calculations of all descriptors from their {@link DescriptorCalculationFactory} with the
     * configuration of this calculator (fingerprint pool size, circular fingerprint size and circular fingerprint
     * counts).
     *
     * @return a new map with the calculation of every descriptor
     * @throws UnsupportedOperationException if a calculation cannot be created, which should never happen
     */
    private EnumMap<Descriptor, DescriptorCalculation> createCalculations() {
        try {
            CalculationSettings settings = new CalculationSettings(this.fingerprintPoolSize,
                    this.circularFingerprintSize, this.isCircularFingerprintCounts);
            EnumMap<Descriptor, DescriptorCalculation> newCalculations = new EnumMap<>(Descriptor.class);
            for (Descriptor descriptor : Descriptor.values()) {
                newCalculations.put(descriptor, descriptor.getCalculationFactory().create(settings));
            }
            return newCalculations;
        } catch (Exception exception) {
            throw new UnsupportedOperationException("Failed to initialize descriptors, this should never happen. ", exception);
        }
    }

    /**
     * Returns the pool size for fingerprinter instances of this calculator.
     *
     * @return fingerprint pool size
     */
    public int getFingerprintPoolSize() {
        return this.fingerprintPoolSize;
    }

    /**
     * Returns the size used for circular fingerprints by this calculator.
     *
     * @return circular fingerprint size
     */
    public int getCircularFingerprintSize() {
        return this.circularFingerprintSize;
    }

    /**
     * Returns whether circular fingerprints are calculated as folded count fingerprints by this calculator.
     *
     * @return true: folded count fingerprints; false: bit fingerprints (default)
     */
    public boolean isCircularFingerprintCounts() {
        return this.isCircularFingerprintCounts;
    }

    /**
     * Returns the number of components of the given descriptor for this calculator, i.e. its number of matrix
     * columns: the circular fingerprint size of this calculator for circular fingerprints (CIRCULAR_FINGERPRINTER_*),
     * otherwise the (fixed) number of components of the descriptor, see
     * {@link Descriptor#getDescriptorComponentNumber()}.
     *
     * @param descriptor the descriptor (IS NOT CHANGED); may not be null
     * @return the number of components of the descriptor
     * @throws NullPointerException if descriptor is null
     */
    public int getNumberOfComponents(Descriptor descriptor) {
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.getNumberOfComponents: descriptor is null.");
        }
        if (descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_")) {
            return this.circularFingerprintSize;
        }
        return descriptor.getDescriptorComponentNumber();
    }

    /**
     * Returns the sum of the numbers of components of the given descriptors for this calculator, i.e. the number of
     * matrix columns needed for them (see {@link #getNumberOfComponents(Descriptor)}).
     *
     * @param descriptors Array of descriptors (IS NOT CHANGED); if it is empty, 0 is returned;
     *                     may not be null (the array or one of its elements)
     * @return Sum of number of calculated components of array of descriptors
     * @throws NullPointerException if descriptors or one of its elements is null
     */
    public int getNumberOfComponents(Descriptor[] descriptors) {
        if (descriptors == null) {
            throw new NullPointerException("DescriptorCalculator.getNumberOfComponents: descriptors is null.");
        }
        int totalNumberOfComponents = 0;
        for (Descriptor descriptor : descriptors) {
            if (descriptor == null) {
                throw new NullPointerException("DescriptorCalculator.getNumberOfComponents: At least one descriptor in descriptors is null.");
            }
            totalNumberOfComponents += this.getNumberOfComponents(descriptor);
        }
        return totalNumberOfComponents;
    }

    /**
     * Returns the type of the values of all components of the given descriptor for this calculator: for circular
     * fingerprints (CIRCULAR_FINGERPRINTER_*), {@link DescriptorComponentType#COUNT} if this calculator calculates
     * them as folded count fingerprints, otherwise the type defined in the enum, see
     * {@link Descriptor#getComponentType()} (e.g. {@link DescriptorComponentType#BINARY} for circular fingerprints in
     * the default bit mode).
     *
     * @param descriptor the descriptor (IS NOT CHANGED); may not be null
     * @return the type of the values of all components of the descriptor
     * @throws NullPointerException if descriptor is null
     */
    public DescriptorComponentType getComponentType(Descriptor descriptor) {
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.getComponentType: descriptor is null.");
        }
        if (this.isCircularFingerprintCounts && descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_")) {
            return DescriptorComponentType.COUNT;
        }
        return descriptor.getComponentType();
    }

    /**
     * Returns the type of the values of every matrix column for the given descriptors, in the order of the columns,
     * i.e. the i-th type belongs to column startIndex + i of a matrix calculated by this calculator with the same
     * descriptor array (see {@link #getComponentType(Descriptor)}). The types can be used together with the column
     * names (see {@link #getColumnNames(Descriptor[])}), e.g. for choosing a suitable preprocessing of each column
     * before machine learning. The number of types is equal to {@link #getNumberOfComponents(Descriptor[])}.
     *
     * @param descriptors Array of descriptors (IS NOT CHANGED); if it is empty, an empty array is returned;
     *                     may not be null (the array or one of its elements)
     * @return the component types, one per column
     * @throws NullPointerException if descriptors or one of its elements is null
     */
    public DescriptorComponentType[] getComponentTypes(Descriptor[] descriptors) {
        DescriptorComponentType[] componentTypes = new DescriptorComponentType[this.getNumberOfComponents(descriptors)];
        int startIndex = 0;
        for (Descriptor descriptor : descriptors) {
            int endIndex = startIndex + this.getNumberOfComponents(descriptor);
            Arrays.fill(componentTypes, startIndex, endIndex, this.getComponentType(descriptor));
            startIndex = endIndex;
        }
        return componentTypes;
    }

    /**
     * Returns the descriptor and component index for a given column index of a matrix calculated by this calculator
     * with the given descriptor array (see {@link Descriptor#getDescriptorAndComponentIndex(Descriptor[], int)}, which
     * uses the default circular fingerprint size).
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The column index relative to the start index of the descriptors (must be &gt;= 0 and &lt; total
     *                component count)
     * @return A two-element int array: [descriptor index in descriptors, component index within that descriptor]
     * @throws NullPointerException if descriptors or one of its elements is null
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public int[] getDescriptorAndComponentIndex(Descriptor[] descriptors, int anIndex) throws IllegalArgumentException {
        int numberOfComponents = this.getNumberOfComponents(descriptors);
        if (descriptors.length == 0) {
            throw new IllegalArgumentException("DescriptorCalculator.getDescriptorAndComponentIndex: descriptors must not be empty.");
        }
        if (anIndex < 0) {
            throw new IllegalArgumentException("DescriptorCalculator.getDescriptorAndComponentIndex: anIndex must be >= 0.");
        }
        if (anIndex >= numberOfComponents) {
            throw new IllegalArgumentException("DescriptorCalculator.getDescriptorAndComponentIndex: anIndex (" + anIndex
                    + ") exceeds total component count (" + numberOfComponents + ").");
        }
        int cumulativeCount = 0;
        int[] result = new int[2];
        for (int i = 0; i < descriptors.length; i++) {
            int nextCumulativeCount = cumulativeCount + this.getNumberOfComponents(descriptors[i]);
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
     * Returns the descriptor name, the component index and the component name (column name, see
     * {@link #getColumnNames(Descriptor[])}) for a given column index of a matrix calculated by this calculator with
     * the given descriptor array, e.g. ["Wiener Numbers", "1", "WPOL"].
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The column index relative to the start index of the descriptors (must be &gt;= 0 and &lt; total
     *                component count)
     * @return A String array: [descriptor name, component index as string, component name]
     * @throws NullPointerException if descriptors or one of its elements is null
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public String[] getDescriptorAndComponentInfo(Descriptor[] descriptors, int anIndex) throws IllegalArgumentException {
        int[] indices = this.getDescriptorAndComponentIndex(descriptors, anIndex);
        Descriptor descriptor = descriptors[indices[0]];
        return new String[] {
                descriptor.getName(),
                String.valueOf(indices[1]),
                this.calculations.get(descriptor).getComponentNames()[indices[1]]
        };
    }

    /**
     * Returns the name of every matrix column for the given descriptors, in the order of the columns, i.e. the i-th
     * name belongs to column startIndex + i of a matrix calculated by this calculator with the same descriptor array.
     * For molecular descriptors, these are the descriptor names of the CDK (e.g. "WPATH" and "WPOL" for the Wiener
     * numbers); for fingerprints, which have no bit names in the CDK, they are generated as prefix and bit index (e.g.
     * "PubChem_0" or "ECFP4_1023"). The number of names is equal to {@link #getNumberOfComponents(Descriptor[])}.
     *
     * @param descriptors Array of descriptors (IS NOT CHANGED); if it is empty, an empty array is returned;
     *                     may not be null (the array or one of its elements)
     * @return the column names
     * @throws NullPointerException if descriptors or one of its elements is null
     */
    public String[] getColumnNames(Descriptor[] descriptors) {
        if (descriptors == null) {
            throw new NullPointerException("DescriptorCalculator.getColumnNames: descriptors is null.");
        }
        List<String> columnNames = new ArrayList<>();
        for (Descriptor descriptor : descriptors) {
            if (descriptor == null) {
                throw new NullPointerException("DescriptorCalculator.getColumnNames: At least one descriptor in descriptors is null.");
            }
            columnNames.addAll(Arrays.asList(this.calculations.get(descriptor).getComponentNames()));
        }
        return columnNames.toArray(new String[0]);
    }


    /**
     * Uses a specified aromaticity model to modify molecule.
     * The method performs a complete workflow of:
     * <ol>
     *     <li>Suppresses explicit hydrogens.</li>
     *     <li>Perceives atom types and configures atoms.</li>
     *     <li>Clears existing aromaticity flags.</li>
     *     <li>Detects rings.</li>
     *     <li>Applies the specified {@link ElectronDonation} model.</li>
     * </ol>
     * Note: This method changes the input molecule by applying the specified aromaticity model and reconfiguring its hydrogens.
     *
     * @param molecule Molecule that will be modified (IS CHANGED)
     * @param electronDonationModel The aromaticity model that will be used for aromaticity detection
     * @throws CDKException if adding implicit hydrogen atoms or detection of atom types fails
     */
    public static void setAromaticity(
            IAtomContainer molecule,
            ElectronDonation electronDonationModel
    ) throws CDKException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (molecule.isEmpty()) {
            return;
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("Aromaticity model must not be null");
        }
        AtomContainerManipulator.normalizeHydrogens(molecule, HydrogenState.Minimal);
        // Needed for VABC Descriptor and CDK_AtomTypes
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        // Clears all aromatic flags before applying the aromaticity model.
        Aromaticity.clear(molecule);
        Cycles.markRingAtomsAndBonds(molecule);
        Aromaticity.apply(electronDonationModel, molecule);
    }

    /**
     * Calculates a single descriptor for the given molecule and returns the result as a float array.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for
     * single descriptor calculations. If the calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     * <p>
     * <b>Important:</b> The molecule must have aromaticity already perceived if required by the descriptor.
     * Use {@link #setAromaticity(IAtomContainer, ElectronDonation)} before calling this method if needed.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param molecule The molecule to calculate the descriptor for (must not be null);
     *                  aromaticity must be perceived beforehand if required; if it is empty, an empty array is returned
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getNumberOfComponents(Descriptor)}. Contains NaN values if calculation fails.
     */
    public float[] calculateDescriptor(Descriptor descriptor, IAtomContainer molecule)  {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }
        if (descriptor.requires3DCoordinates() && !GeometryUtil.has3DCoordinates(molecule)) {
            float[] result = new float[this.getNumberOfComponents(descriptor)];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("DescriptorCalculator.calculateDescriptor: Descriptor " + descriptor.getName() + " requires 3D coordinates, but the molecule does not have 3D coordinates.");
            return result;
        }

        float[] result = new float[this.getNumberOfComponents(descriptor)];
        try {
            boolean success;
            if (descriptor.needsExplicitHydrogens()) {
                IAtomContainer moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(molecule);
                success = this.calculate(descriptor, moleculeWithExplicitHydrogens, result, 0);
            } else {
                success = this.calculate(descriptor, molecule, result, 0);
            }
            if (!success) {
                Arrays.fill(result, Float.NaN);
            }
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("Interrupted while calculating descriptor " + descriptor.getName(), exception);
        } catch (Exception exception) {
            // Fill with NaN on failure and log the error
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
        }
        return result;
    }

    /**
     * Calculates a single descriptor for the given molecule with automatic aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically perceives aromaticity using the specified electron donation model before calculating the descriptor.
     * If aromaticity perception or calculation fails, all values in the result array will be set to {@link Float#NaN}.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param molecule The molecule to calculate the descriptor for (must not be null); if it is empty, an empty array is returned
     * @param electronDonationModel The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getNumberOfComponents(Descriptor)}. Contains NaN values if calculation fails.
     */
    public float[] calculateDescriptor(Descriptor descriptor,
                                              IAtomContainer molecule,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: electronDonationModel must not be null.");
        }

        try {
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
            return this.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[this.getNumberOfComponents(descriptor)];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
            return result;
        }
    }

    /**
     * Calculates a single descriptor for the given SMILES string with automatic parsing and aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically parses the SMILES string into a molecule,
     * perceives aromaticity using the specified electron donation model, and then calculates the descriptor.
     * If parsing, aromaticity perception, or calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param smilesString The SMILES string representing the molecule (must not be null); if it is blank, an empty array is returned
     * @param electronDonationModel The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getNumberOfComponents(Descriptor)}. Contains NaN values if calculation fails.
     */
    public float[] calculateDescriptor(Descriptor descriptor,
                                              String smilesString,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (smilesString == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: smilesString must not be null.");
        }
        if (smilesString.trim().isEmpty()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: electronDonationModel must not be null.");
        }
        if (descriptor.requires3DCoordinates()) {
            float[] result = new float[this.getNumberOfComponents(descriptor)];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("DescriptorCalculator.calculateDescriptor: Descriptor " + descriptor.getName() + " requires 3D coordinates, which cannot be computed from a SMILES string.");
            return result;
        }

        try {
            IAtomContainer molecule = DescriptorCalculator.SMILES_PARSER.parseSmiles(smilesString);
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
            return this.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[this.getNumberOfComponents(descriptor)];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
            return result;
        }
    }


    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecule batches. If parallel computation is used, the atom
     * container batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size is set in the constructor, see {@link #DescriptorCalculator(int, int, boolean)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainerArray Array of molecules. Note: atomContainerArray[i] corresponds to matrix[i] data
     *                             vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED); every
     *                             molecule object may only occur once in the array (use copies for duplicates)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                atomContainerArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param batchSize Number of molecules to process in each batch
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to which the NaN positions are added as [moleculeIndex, componentIndex] pairs
     *                      (MAY BE CHANGED); any modifiable list can be used, also for parallel calculations. The
     *                      positions are added after the calculation, sorted by molecule index and component index.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN, a global exception occurred, or the calculation was aborted before it started
     *         because the given descriptor or molecule array is empty or a descriptor requires 3D coordinates that are
     *         not available (in this case, the matrix is not changed)
     * @throws IllegalArgumentException if the matrix dimensions are invalid, if batch size is &lt;= 0 or if a
     *                                  molecule object occurs more than once in atomContainerArray
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected batch is aborted silently
     */
    public boolean setDescriptorsForMoleculesByBatchParallelization(
            Descriptor[] descriptors,
            IAtomContainer[] atomContainerArray,
            float[][] matrix,
            int startIndex,
            int batchSize,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesByBatchParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return false;
        }
        if (!DescriptorCalculator.validateAtomContainerArray(atomContainerArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given atom container array is empty, calculation aborted.");
            return false;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    if (!GeometryUtil.has3DCoordinates(atomContainerArray[i])) {
                        DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, but molecule " + i + " does not have 3D coordinates. Calculation aborted.");
                        return false;
                    }
                }
                break;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        this.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

        int numberOfMolecules = atomContainerArray.length;
        int numberOfBatches = (int) Math.ceil((double) numberOfMolecules / batchSize);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += this.getNumberOfComponents(descriptors[i]);
        }

        // NaN positions are collected in a thread-safe queue and handed over to nanPositionsList at the end
        Collection<int[]> collectedNanPositions = new ConcurrentLinkedQueue<>();
        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                IntStream.range(0, numberOfBatches).parallel().forEach(batchIndex -> {
                    int batchStart = batchIndex * batchSize;
                    int batchEnd = Math.min(batchStart + batchSize, numberOfMolecules);

                    for (int i = batchStart; i < batchEnd; i++) {
                        try {
                            boolean success = this.setDescriptorsForSingleMolecule(
                                    descriptors,
                                    atomContainerArray[i],
                                    matrix[i],
                                    startIndices,
                                    i,
                                    collectedNanPositions
                            );
                            if (!success) {
                                hasNaN.set(true);
                            }
                        } catch (InterruptedException exception) {
                            Thread.currentThread().interrupt(); // Preserve interrupt status
                            return;
                        } catch (Exception exception) {
                            hasNaN.set(true);
                            DescriptorCalculator.LOGGER.warn(
                                    String.format("DescriptorCalculator.setDescriptorsForMoleculesByBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!this.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, collectedNanPositions)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("DescriptorCalculator.setDescriptorsForMoleculesByBatchParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "DescriptorCalculator.setDescriptorsForMoleculesByBatchParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        } finally {
            // hand over in the calling thread, sorted by molecule index and component index
            DescriptorCalculator.addSortedNanPositions(collectedNanPositions, nanPositionsList);
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecule batches. If parallel computation is used, the smiles string
     * batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size is set in the constructor, see {@link #DescriptorCalculator(int, int, boolean)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesStringArray Array of molecule smiles strings. Note: moleculeSmilesStringArray[i] corresponds to matrix[i] data
     *                                   vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                moleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param batchSize Number of molecules to process in each batch
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                                Aromaticity.Model.Daylight is used as default.
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to which the NaN positions are added as [moleculeIndex, componentIndex] pairs
     *                      (MAY BE CHANGED); any modifiable list can be used, also for parallel calculations. The
     *                      positions are added after the calculation, sorted by molecule index and component index.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN, a global exception occurred, or the calculation was aborted before it started
     *         because the given descriptor or molecule array is empty or a descriptor requires 3D coordinates that are
     *         not available (in this case, the matrix is not changed)
     * @throws IllegalArgumentException if the matrix dimensions are invalid or if batch size is &lt;= 0
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected batch is aborted silently
     */
    public boolean setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
            Descriptor[] descriptors,
            String[] moleculeSmilesStringArray,
            float[][] matrix,
            int startIndex,
            int batchSize,
            ElectronDonation electronDonationModel,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculeBySmilesStringsBatchParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return false;
        }
        if (!DescriptorCalculator.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given SMILES string array is empty, calculation aborted.");
            return false;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, which cannot be computed from SMILES strings. Calculation aborted.");
                return false;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        this.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

        int numberOfMolecules = moleculeSmilesStringArray.length;
        int numberOfBatches = (int) Math.ceil((double) numberOfMolecules / batchSize);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += this.getNumberOfComponents(descriptors[i]);
        }

        // NaN positions are collected in a thread-safe queue and handed over to nanPositionsList at the end
        Collection<int[]> collectedNanPositions = new ConcurrentLinkedQueue<>();
        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                IntStream.range(0, numberOfBatches).parallel().forEach(batchIndex -> {
                    int batchStart = batchIndex * batchSize;
                    int batchEnd = Math.min(batchStart + batchSize, numberOfMolecules);

                    for (int i = batchStart; i < batchEnd; i++) {
                        try {
                            boolean success = this.setDescriptorsForSingleMoleculeSmilesString(
                                    descriptors,
                                    moleculeSmilesStringArray[i],
                                    matrix[i],
                                    startIndices,
                                    i,
                                    electronDonationModel,
                                    collectedNanPositions
                            );
                            if (!success) {
                                hasNaN.set(true);
                            }
                        } catch (InterruptedException exception) {
                            Thread.currentThread().interrupt(); // Preserve interrupt status
                            return;
                        } catch (Exception exception) {
                            hasNaN.set(true);
                            DescriptorCalculator.LOGGER.warn(
                                    String.format("DescriptorCalculator.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!this.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, collectedNanPositions)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("DescriptorCalculator.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "DescriptorCalculator.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        } finally {
            // hand over in the calling thread, sorted by molecule index and component index
            DescriptorCalculator.addSortedNanPositions(collectedNanPositions, nanPositionsList);
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size is set in the constructor, see {@link #DescriptorCalculator(int, int, boolean)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainerArray Array of molecules. Note: atomContainerArray[i] corresponds to matrix[i] data
     *                             vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED); every
     *                             molecule object may only occur once in the array (use copies for duplicates)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                atomContainerArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to which the NaN positions are added as [moleculeIndex, componentIndex] pairs
     *                      (MAY BE CHANGED); any modifiable list can be used, also for parallel calculations. The
     *                      positions are added after the calculation, sorted by molecule index and component index.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN, a global exception occurred, or the calculation was aborted before it started
     *         because the given descriptor or molecule array is empty or a descriptor requires 3D coordinates that are
     *         not available (in this case, the matrix is not changed)
     * @throws IllegalArgumentException if the matrix dimensions are invalid or if a molecule object occurs
     *                                  more than once in atomContainerArray
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected molecule is skipped silently
     */
    public boolean setDescriptorsForMoleculesByMoleculeParallelization(
            Descriptor[] descriptors,
            IAtomContainer[] atomContainerArray,
            float[][] matrix,
            int startIndex,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesByMoleculeParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return false;
        }
        if (!DescriptorCalculator.validateAtomContainerArray(atomContainerArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given atom container array is empty, calculation aborted.");
            return false;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    if (!GeometryUtil.has3DCoordinates(atomContainerArray[i])) {
                        DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, but molecule " + i + " does not have 3D coordinates. Calculation aborted.");
                        return false;
                    }
                }
                break;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        this.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += this.getNumberOfComponents(descriptors[i]);
        }

        // NaN positions are collected in a thread-safe queue and handed over to nanPositionsList at the end
        Collection<int[]> collectedNanPositions = new ConcurrentLinkedQueue<>();
        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                // Advise by Oracle: Parallel streams should use the common Fork-join pool
                IntStream.range(0, atomContainerArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                boolean success = this.setDescriptorsForSingleMolecule(
                                        descriptors,
                                        atomContainerArray[i],
                                        matrix[i],
                                        startIndices,
                                        i,
                                        collectedNanPositions
                                );
                                if (!success) {
                                    hasNaN.set(true);
                                }
                            } catch (InterruptedException exception) {
                                Thread.currentThread().interrupt(); // Preserve interrupt status
                            } catch (Exception exception) {
                                hasNaN.set(true);
                                DescriptorCalculator.LOGGER.warn(
                                        String.format("DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelization: One descriptor calculation caused an exception, molecule index: %d.", i),
                                        exception
                                );
                            }
                        }
                );
            } else {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    try {
                        if (!this.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, collectedNanPositions)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelization: Global exception occurred in descriptor calculation: ",
                    exception
            );
            return false;
        } finally {
            // hand over in the calling thread, sorted by molecule index and component index
            DescriptorCalculator.addSortedNanPositions(collectedNanPositions, nanPositionsList);
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecules. If parallel computation is used, the smiles strings
     * (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size is set in the constructor, see {@link #DescriptorCalculator(int, int, boolean)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesStringArray Array of molecule smiles strings. Note: moleculeSmilesStringArray[i] corresponds to matrix[i] data
     *                                   vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                moleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                                Aromaticity.Model.Daylight is used as default.
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to which the NaN positions are added as [moleculeIndex, componentIndex] pairs
     *                      (MAY BE CHANGED); any modifiable list can be used, also for parallel calculations. The
     *                      positions are added after the calculation, sorted by molecule index and component index.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN, a global exception occurred, or the calculation was aborted before it started
     *         because the given descriptor or molecule array is empty or a descriptor requires 3D coordinates that are
     *         not available (in this case, the matrix is not changed)
     * @throws IllegalArgumentException if the matrix dimensions are invalid
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected molecule is skipped silently
     */
    public boolean setDescriptorsForMoleculesBySmilesStringParallelization(
            Descriptor[] descriptors,
            String[] moleculeSmilesStringArray,
            float[][] matrix,
            int startIndex,
            ElectronDonation electronDonationModel,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesBySmilesStringParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return false;
        }
        if (!DescriptorCalculator.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given SMILES string array is empty, calculation aborted.");
            return false;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, which cannot be computed from SMILES strings. Calculation aborted.");
                return false;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        this.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

        int numberOfMolecules = moleculeSmilesStringArray.length;
        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += this.getNumberOfComponents(descriptors[i]);
        }
        // NaN positions are collected in a thread-safe queue and handed over to nanPositionsList at the end
        Collection<int[]> collectedNanPositions = new ConcurrentLinkedQueue<>();
        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                // Advise by Oracle: Parallel streams should use the common Fork-join pool
                IntStream.range(0, moleculeSmilesStringArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                boolean success = this.setDescriptorsForSingleMoleculeSmilesString(
                                        descriptors,
                                        moleculeSmilesStringArray[i],
                                        matrix[i],
                                        startIndices,
                                        i,
                                        electronDonationModel,
                                        collectedNanPositions
                                );
                                if (!success) {
                                    hasNaN.set(true);
                                }
                            } catch (InterruptedException exception) {
                                Thread.currentThread().interrupt(); // Preserve interrupt status
                            } catch (Exception exception) {
                                hasNaN.set(true);
                                DescriptorCalculator.LOGGER.warn(
                                        String.format("DescriptorCalculator.setDescriptorsForMoleculesBySmilesStringParallelization: Exception in molecule index: %d", i),
                                        exception
                                );
                            }

                        }
                );
            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!this.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, collectedNanPositions)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("DescriptorCalculator.setDescriptorsForMoleculesBySmilesStringParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "DescriptorCalculator.setDescriptorsForMoleculesBySmilesStringParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        } finally {
            // hand over in the calling thread, sorted by molecule index and component index
            DescriptorCalculator.addSortedNanPositions(collectedNanPositions, nanPositionsList);
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vector (that corresponds to atomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in vector to be filled with calculated components of descriptors (each
     *                      descriptor in descriptors has its dedicated start index here)
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws CloneNotSupportedException if copying the molecule fails
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #setDescriptor})
     */
    private boolean setDescriptorsForSingleMolecule (
            Descriptor[] descriptors,
            IAtomContainer atomContainer,
            float[] vector,
            int[] aStartIndices,
            int moleculeIndex,
            Collection<int[]> nanPositionsList
    ) throws CloneNotSupportedException, InterruptedException {
        IAtomContainer moleculeWithExplicitHydrogens = null;
        for (Descriptor descriptor : descriptors) {
            if (descriptor.needsExplicitHydrogens()){
                moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(atomContainer);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : atomContainer;
            if (!this.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
                isSuccessful = false;
            }
        }
        return isSuccessful;
    }

    /**
     * Sets calculated descriptor components in vector (that corresponds to atomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesString Molecule SMILES String (IS NOT CHANGED)
     * @param vector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in vector to be filled with calculated components of descriptors (each
     *                      descriptor in descriptors has its dedicated start index here)
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     * Aromaticity.Model.Daylight is used as default.
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws CloneNotSupportedException if copying the molecule fails
     * @throws CDKException if SMILES parsing or aromaticity detection fails
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #setDescriptor})
     */
    private boolean setDescriptorsForSingleMoleculeSmilesString (
            Descriptor[] descriptors,
            String moleculeSmilesString,
            float[] vector,
            int[] aStartIndices,
            int moleculeIndex,
            ElectronDonation electronDonationModel,
            Collection<int[]> nanPositionsList
    ) throws CDKException, CloneNotSupportedException, InterruptedException {
        IAtomContainer molecule = DescriptorCalculator.SMILES_PARSER.parseSmiles(moleculeSmilesString);
        if (electronDonationModel == null) {
            DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
        } else {
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
        }

        IAtomContainer moleculeWithExplicitHydrogens = null;
        for (Descriptor descriptor : descriptors) {
            if (descriptor.needsExplicitHydrogens()){
                moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(molecule);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : molecule;
            if (!this.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
                isSuccessful = false;
            }
        }
        return isSuccessful;
    }

    /**
     * Sets component values of descriptor for atomContainer in vector beginning with startIndex.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Vector of molecule (row of data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param startIndex Start index in vector to be filled with calculated components of descriptors
     * @param moleculeIndex Index of the current molecule being processed (row index in data matrix)
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #calculate}); all other exceptions are caught,
     *                              logged, and converted to NaN values in vector
     */
    private boolean setDescriptor(
            Descriptor descriptor,
            IAtomContainer atomContainer,
            float[] vector,
            int startIndex,
            int moleculeIndex,
            Collection<int[]> nanPositionsList
    ) throws InterruptedException {
        try {
            boolean success = this.calculate(descriptor, atomContainer, vector, startIndex);
            if (!success) {
                int numComponents = this.getNumberOfComponents(descriptor);
                for (int i = 0; i < numComponents; i++) {
                    vector[startIndex + i] = Float.NaN;
                    if (nanPositionsList != null) {
                        nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                    }
                }
                return false;
            }

            // Check for NaN values in the calculated result and track them
            int numComponents = this.getNumberOfComponents(descriptor);
            return !DescriptorCalculator.checkAndTrackNaNValues(vector, startIndex, numComponents, moleculeIndex, nanPositionsList);
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            int numComponents = this.getNumberOfComponents(descriptor);
            for (int i = 0; i < numComponents; i++) {
                vector[startIndex + i] = Float.NaN;
                // Track NaN position if nanPositionsList is provided
                if (nanPositionsList != null) {
                    nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                }
            }
            DescriptorCalculator.LOGGER.warn(
                    String.format("DescriptorCalculator.setDescriptor: An exception occurred while calculating descriptor %s for molecule index %d.",
                            descriptor,
                            moleculeIndex),
                    exception
            );
            return false;
        }
    }

    /**
     * Helper method to check for NaN values in a calculated descriptor result and track their positions.
     * This method is thread-safe when used with a thread-safe collection (the public methods use a
     * ConcurrentLinkedQueue). Note that not the entire data vector is checked but only the positions
     * [startIndex, startIndex + numComponents].
     *
     * @param vector The vector containing calculated descriptor values
     * @param startIndex The start index in the vector for this descriptor
     * @param numComponents The number of components for this descriptor
     * @param moleculeIndex The index of the current molecule
     * @param nanPositionsList Collection to track NaN positions (can be null if NaN positions should not be tracked);
     *                      must be thread-safe for parallel access; will be filled with int[] {[moleculeIndex, componentIndex]}
     *                      pairs of NaN value positions (i.e. row and column index in the data matrix)
     * @return true if any NaN values were found, false otherwise
     */
    static boolean checkAndTrackNaNValues(
            float[] vector,
            int startIndex,
            int numComponents,
            int moleculeIndex,
            Collection<int[]> nanPositionsList
    ) {
        boolean foundNaN = false;
        for (int i = 0; i < numComponents; i++) {
            if (Float.isNaN(vector[startIndex + i])) {
                foundNaN = true;
                // Track NaN position if nanPositionsList is provided
                // Note: nanPositionsList must be thread-safe for parallel access
                if (nanPositionsList != null) {
                    nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                }
            }
        }
        return foundNaN;
    }

    /**
     * Adds the collected NaN positions to the given list, sorted by molecule index and then by component index, so
     * that the order is reproducible and independent of thread scheduling.
     *
     * @param collectedNanPositions NaN positions collected during the calculation (IS NOT CHANGED)
     * @param nanPositionsList list to which the sorted NaN positions are added (MAY BE CHANGED)
     */
    private static void addSortedNanPositions(Collection<int[]> collectedNanPositions, List<int[]> nanPositionsList) {
        List<int[]> sortedNanPositions = new ArrayList<>(collectedNanPositions);
        sortedNanPositions.sort(Comparator.<int[]>comparingInt(position -> position[0]).thenComparingInt(position -> position[1]));
        nanPositionsList.addAll(sortedNanPositions);
    }

    /**
     * Validates descriptor array input. Throws NullPointerExceptions if the array or one of its elements is null.
     * Returns false if the array is empty.
     *
     * @param descriptorArray the descriptor array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if descriptor array or one of its elements is null
     * @return false if descriptor array is empty; true otherwise
     */
    static boolean validateDescriptors(Descriptor[] descriptorArray, String methodName) {
        if (descriptorArray == null) {
            throw new NullPointerException(methodName + ": Given descriptor array is null");
        }
        if (descriptorArray.length == 0) {
            return false;
        }
        for (Descriptor descriptor : descriptorArray) {
            if (descriptor == null) {
                throw new NullPointerException(methodName + ": A single descriptor in descriptors is null.");
            }
        }
        return true;
    }

    /**
     * Validates an atom container array input and its contents. Throws NullPointerExceptions if the array or one of
     * its elements is null and an IllegalArgumentException if the same molecule object occurs more than once in the
     * array (the CDK descriptors may change a molecule during calculation, e.g. by setting properties, so one molecule
     * object must not be calculated by several threads at the same time). Returns false if the array is empty. Note
     * that the molecules in the array are allowed to be empty.
     *
     * @param atomContainerArray the atom container array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if molecule array or one of its elements is null
     * @throws IllegalArgumentException if the same molecule object occurs more than once in the array
     * @return false if the atom container array is empty; true otherwise
     */
    static boolean validateAtomContainerArray(IAtomContainer[] atomContainerArray, String methodName) {
        if (atomContainerArray == null) {
            throw new NullPointerException(methodName + ": atomContainerArray is null.");
        }
        if (atomContainerArray.length == 0) {
            return false;
        }
        Set<IAtomContainer> molecules = Collections.newSetFromMap(new IdentityHashMap<>(atomContainerArray.length));
        for (int i = 0; i < atomContainerArray.length; i++) {
            if (atomContainerArray[i] == null) {
                throw new NullPointerException(methodName + ": A single molecule in atomContainerArray is null.");
            }
            if (!molecules.add(atomContainerArray[i])) {
                throw new IllegalArgumentException(methodName + ": The molecule object at index " + i
                        + " occurs more than once in atomContainerArray, use a copy instead.");
            }
            //do nothing if a molecule is empty
        }
        return true;
    }

    /**
     * Validates a SMILES string array input and its content. Throws NullPointerExceptions if the array or one of
     * its elements is null. Returns false if the array is empty.
     *
     * @param moleculeSmilesStringArray the molecule SMILES string array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if the SMILES string array or one of its elements is null
     * @return false if the SMILES string array is empty; true otherwise
     */
    private static boolean validateSmilesStringArray(String[] moleculeSmilesStringArray, String methodName) {
        if (moleculeSmilesStringArray == null) {
            throw new NullPointerException(methodName + ": moleculeSmilesStringArray is null.");
        }
        if (moleculeSmilesStringArray.length == 0) {
            return false;
        }
        for (String molecule : moleculeSmilesStringArray) {
            if (molecule == null) {
                throw new NullPointerException(methodName + ": A single SMILES string in moleculeSmilesStringArray is null.");
            }
            //do nothing if a String is empty
        }
        return true;
    }

    /**
     * Validates matrix dimensions and start index. Throws a NullPointerException if the matrix or one of its rows is
     * null. Returns false if the matrix is empty. Note that some parameters like {@code descriptors} have their own validation methods.
     *
     * @param matrix the matrix to validate
     * @param descriptors the descriptors intended to be calculated for validating the matrix dimensions (min. nr. of columns)
     * @param moleculeArray the input molecule array (either an {@code IAtomContainer[]} or a {@code String[]}) for
     *                       validating the matrix dimensions (exact(!) nr. of rows)
     * @param startIndex the start index in the matrix (starting there, enough columns must be left to fit all descriptor results)
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if the matrix or one of its rows is null
     * @throws IllegalArgumentException if matrix dimensions are invalid
     */
    void validateMatrix(
            float[][] matrix,
            Descriptor[] descriptors,
            Object[] moleculeArray,
            int startIndex,
            String methodName
    ) {
        if (matrix == null) {
            throw new NullPointerException(methodName + ": matrix is null.");
        }
        if (matrix.length == 0) {
            throw new IllegalArgumentException(methodName + ": matrix is of length 0.");
        }
        if (matrix.length != moleculeArray.length) {
            throw new IllegalArgumentException(methodName + ": matrix and moleculeArray must have the same length.");
        }
        if (!(moleculeArray instanceof IAtomContainer[]) && !(moleculeArray instanceof String[])) {
            throw new IllegalArgumentException(methodName + ": moleculeArray must be an IAtomContainer[] or a String[].");
        }
        int colNum = -1;
        for (float[] vector : matrix) {
            if (vector == null) {
                throw new NullPointerException(methodName + ": A vector in matrix is null.");
            }
            if (vector.length == 0) {
                throw new IllegalArgumentException(methodName + ": A vector in matrix has length 0.");
            }
            if (colNum == -1) {
                colNum = vector.length;
            }
            if (vector.length != colNum) {
                throw new IllegalArgumentException(methodName + ": All vectors in matrix must have the same length.");
            }
            if (startIndex > vector.length) {
                throw new IllegalArgumentException(methodName + ": startIndex is greater than or equal to vector length.");
            }
            int numberOfComponents = this.getNumberOfComponents(descriptors);
            if (startIndex + numberOfComponents > vector.length) {
                throw new IllegalArgumentException(methodName + ": Not enough space in vector for descriptors.");
            }
        }
    }

    /**
     * Creates a deep copy of the input molecule.
     * <p>
     * Notes:
     * <ul>
     *     <li>The copy is built with the same {@link org.openscience.cdk.interfaces.IChemObjectBuilder}
     *         as the input, using the in-container factory methods
     *         {@code newAtom(int, int)} / {@code newBond(IAtom, IAtom, IBond.Order)} which are
     *         faster than the dynamic {@code newInstance(...)} dispatch.</li>
     *     <li>Implicit hydrogen counts are preserved exactly (including a possible {@code null}).</li>
     *     <li>Bond endpoints are resolved through an explicit original-to-copy atom map, so no
     *         assumption is made about identical atom indices in both containers.</li>
     *     <li>If needed, atom types must be perceived and configured manually after creation.</li>
     * </ul>
     *
     * @param molecule Source molecule to be copied (NOT MODIFIED)
     * @return New independent instance of the molecule
     * @throws NullPointerException     if {@code molecule} is {@code null}
     * @throws CloneNotSupportedException if the molecule cannot be properly copied
     */
    static IAtomContainer copyMolecule(IAtomContainer molecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        try {
            IAtomContainer moleculeCopy = molecule.getBuilder().newInstance(IAtomContainer.class);
            // Molecule-level properties and identifier
            if (molecule.getProperties() != null) {
                moleculeCopy.addProperties(new HashMap<>(molecule.getProperties()));
            }
            // Explicit original -> copy maps; used for bonds, stereo elements, electrons and lone pairs.
            int atomCapacity = Math.max(16, (int) (molecule.getAtomCount() / 0.75f) + 1);
            int bondCapacity = Math.max(16, (int) (molecule.getBondCount() / 0.75f) + 1);
            Map<IAtom, IAtom> atomMap = new HashMap<>(atomCapacity);
            Map<IBond, IBond> bondMap = new HashMap<>(bondCapacity);

            // ---- Atoms ----
            for (IAtom atom : molecule.atoms()) {
                Integer atomicNumber = atom.getAtomicNumber();
                Integer implicitHCount = atom.getImplicitHydrogenCount();
                // fast in-container creation (avoids reflective newInstance dispatch)
                IAtom newAtom = moleculeCopy.newAtom(
                        atomicNumber != null ? atomicNumber : 0,          // 0 == wildcard element
                        implicitHCount != null ? implicitHCount : 0);
                newAtom.setSymbol(atom.getSymbol());
                newAtom.setImplicitHydrogenCount(implicitHCount);         // preserves a possible null
                newAtom.setMassNumber(atom.getMassNumber());
                newAtom.setFormalCharge(atom.getFormalCharge());
                newAtom.setCharge(atom.getCharge());
                newAtom.setValency(atom.getValency());
                newAtom.setHybridization(atom.getHybridization());
                newAtom.setAtomTypeName(atom.getAtomTypeName());
                // Coordinates (needed for stereo perception / preservation)
                Point2d p2d = atom.getPoint2d();
                if (p2d != null) {
                    newAtom.setPoint2d(new Point2d(p2d));
                }
                Point3d p3d = atom.getPoint3d();
                if (p3d != null) {
                    newAtom.setPoint3d(new Point3d(p3d));
                }
                Point3d fp3d = atom.getFractionalPoint3d();
                if (fp3d != null) {
                    newAtom.setFractionalPoint3d(new Point3d(fp3d));
                }
                // Flags
                newAtom.setIsAromatic(atom.isAromatic());
                newAtom.setIsInRing(atom.isInRing());
                // Full generic property map (carries the SRU unique atom index + detection markers)
                if (atom.getProperties() != null) {
                    newAtom.addProperties(new HashMap<>(atom.getProperties()));
                }
                atomMap.put(atom, newAtom);
            }

            // ---- Bonds ----
            for (IBond bond : molecule.bonds()) {
                IAtom beginCopy = atomMap.get(bond.getBegin());
                IAtom endCopy = atomMap.get(bond.getEnd());
                IBond newBond = moleculeCopy.newBond(beginCopy, endCopy, bond.getOrder());
                newBond.setDisplay(bond.getDisplay());
                newBond.setIsAromatic(bond.isAromatic());
                newBond.setIsInRing(bond.isInRing());
                if (bond.getProperties() != null) {
                    newBond.addProperties(new HashMap<>(bond.getProperties()));
                }
                bondMap.put(bond, newBond);
            }

            // ---- Single electrons ----
            for (ISingleElectron singleElectron : molecule.singleElectrons()) {
                ISingleElectron newSingleElectron = molecule.getBuilder().newInstance(
                        ISingleElectron.class, atomMap.get(singleElectron.getAtom()));
                if (singleElectron.getProperties() != null) {
                    newSingleElectron.addProperties(new HashMap<>(singleElectron.getProperties()));
                }
                moleculeCopy.addSingleElectron(newSingleElectron);
            }

            // ---- Stereo elements (preserve stereochemistry, re-mapped onto the copy) ----
            for (IStereoElement<?, ?> stereoElement : molecule.stereoElements()) {
                moleculeCopy.addStereoElement(stereoElement.map(atomMap, bondMap));
            }

            return moleculeCopy;
        } catch (Exception exception) {
            throw new CloneNotSupportedException("Could not clone molecule: " + exception.getMessage());
        }
    }

    /**
     * Creates a molecule copy with all implicit hydrogen atoms made explicit.
     *
     * @param molecule Source molecule with implicit hydrogens (NOT MODIFIED)
     * @return New molecule (copy of the parameter) with all hydrogens made explicit
     * @throws CloneNotSupportedException If the molecule cannot be properly processed
     */
    static IAtomContainer createMoleculeWithExplicitHydrogens(
            IAtomContainer molecule
    ) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (molecule.isEmpty()) {
            //returns empty atom container
            return DescriptorCalculator.copyMolecule(molecule);
        }
        try {
            // Create a deep copy of the molecule first
            IAtomContainer moleculeCopy = DescriptorCalculator.copyMolecule(molecule);
            // Add explicit hydrogen atoms
            AtomContainerManipulator.normalizeHydrogens(moleculeCopy, HydrogenState.Explicit);
            return moleculeCopy;
        } catch (Exception exception) {
            throw new CloneNotSupportedException("Could not create molecule with explicit hydrogens: " + exception.getMessage());
        }
    }
}
