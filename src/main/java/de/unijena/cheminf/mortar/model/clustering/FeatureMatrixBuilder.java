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

import de.unijena.cheminf.fragment.fingerprint.FragmentFingerprinter;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorCalculator;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorComponentType;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Builds {@link FeatureMatrix} instances for a list of molecules: the fragment fingerprint block (bits or counts of
 * the fragments of a given fragmentation) and the descriptor block (CDK descriptors and fingerprints calculated on the
 * molecules via a {@link DescriptorCalculator}, see {@link #getDescriptorCalculator(int)}).
 * Both blocks can be combined with {@link FeatureMatrix#concatenate(FeatureMatrix...)}.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class FeatureMatrixBuilder {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(FeatureMatrixBuilder.class.getName());
    /**
     * Descriptors by enum constant name, which is the feature group name of their columns.
     */
    private static final Map<String, Descriptor> DESCRIPTOR_OF_GROUP_NAME_MAP = FeatureMatrixBuilder.createDescriptorOfGroupNameMap();
    /**
     * Component index of the Wiener polarity number ("WPOL") in {@link Descriptor#WIENER_NUMBER}. The descriptor engine
     * defines one type per descriptor, CONTINUOUS for the Wiener numbers because of the Wiener path number ("WPATH"), but
     * the Wiener polarity number is the number of atom pairs with a distance of three bonds, i.e. a count; this single
     * exception is handled here instead of in the engine.
     */
    private static final int WIENER_POLARITY_NUMBER_COMPONENT_INDEX = 1;
    //</editor-fold>
    //
    //<editor-fold desc="private static variables" defaultstate="collapsed">
    /**
     * Descriptor calculator of the last requested circular fingerprint size; reused, because creating a calculator
     * creates all CDK descriptor and fingerprinter instances, and replaced when another size is requested. A calculator
     * is immutable and thread-safe, so it can be shared.
     */
    private static volatile DescriptorCalculator descriptorCalculator;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor" defaultstate="collapsed">
    /**
     * Private parameter-less constructor.
     * Introduced because javadoc build complained about classes without declared default constructor.
     */
    private FeatureMatrixBuilder() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Builds the fragment fingerprint block for the given molecules. The fragment dictionary consists of all fragments
     * of the given fragmentation that occur in at least the given number of molecules. The fragments are sorted by
     * descending molecule frequency and then by unique SMILES so that the column order is reproducible.
     * Molecules that did not undergo the given fragmentation get a row of zeros.
     *
     * @param aMoleculeList molecules (rows of the matrix)
     * @param aFragmentationName name of the fragmentation whose fragment frequencies are used
     * @param anIsCountFingerprint true: fragment frequencies are used (type COUNT); false: presence bits (type BIT)
     * @param aMinimumMoleculeFrequency minimum number of molecules a fragment must occur in to become a feature, at least 1
     * @param aMaximumNumberOfFragments maximum number of fragment features; if more fragments fulfil the minimum
     *                                  frequency, the ones occurring in the most molecules are kept (ties broken by
     *                                  unique SMILES); 0 means no limit
     * @return feature matrix of the fragment fingerprint block, may have zero columns
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the minimum molecule frequency is smaller than 1 or the maximum is negative
     */
    public static FeatureMatrix buildFragmentFingerprintMatrix(
            List<MoleculeDataModel> aMoleculeList,
            String aFragmentationName,
            boolean anIsCountFingerprint,
            int aMinimumMoleculeFrequency,
            int aMaximumNumberOfFragments)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aMoleculeList, "Given molecule list is null.");
        Objects.requireNonNull(aFragmentationName, "Given fragmentation name is null.");
        if (aMinimumMoleculeFrequency < 1) {
            throw new IllegalArgumentException("Minimum molecule frequency must be at least 1.");
        }
        if (aMaximumNumberOfFragments < 0) {
            throw new IllegalArgumentException("Maximum number of fragments must not be negative.");
        }
        //count in how many molecules every fragment occurs
        Map<String, Integer> tmpMoleculeFrequencyMap = new HashMap<>();
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            Map<String, Integer> tmpFrequencies = FeatureMatrixBuilder.getFragmentFrequencies(tmpMolecule, aFragmentationName);
            for (String tmpSmiles : tmpFrequencies.keySet()) {
                tmpMoleculeFrequencyMap.merge(tmpSmiles, 1, Integer::sum);
            }
        }
        List<String> tmpDictionary = new ArrayList<>();
        for (Map.Entry<String, Integer> tmpEntry : tmpMoleculeFrequencyMap.entrySet()) {
            if (tmpEntry.getValue() >= aMinimumMoleculeFrequency) {
                tmpDictionary.add(tmpEntry.getKey());
            }
        }
        tmpDictionary.sort(Comparator.comparing((String aSmiles) -> -tmpMoleculeFrequencyMap.get(aSmiles))
                .thenComparing(Comparator.naturalOrder()));
        //the dictionary is sorted by descending molecule frequency, so the most frequent fragments are kept
        if (aMaximumNumberOfFragments > 0 && tmpDictionary.size() > aMaximumNumberOfFragments) {
            tmpDictionary = new ArrayList<>(tmpDictionary.subList(0, aMaximumNumberOfFragments));
        }
        int tmpColumnCount = tmpDictionary.size();
        float[][] tmpData = new float[aMoleculeList.size()][tmpColumnCount];
        if (tmpColumnCount > 0) {
            //the fragment fingerprinter maps the fragment frequencies of a molecule onto the dictionary
            FragmentFingerprinter tmpFingerprinter = new FragmentFingerprinter(tmpDictionary);
            for (int i = 0; i < aMoleculeList.size(); i++) {
                Map<String, Integer> tmpFrequencies = FeatureMatrixBuilder.getFragmentFrequencies(aMoleculeList.get(i), aFragmentationName);
                if (tmpFrequencies.isEmpty()) {
                    continue;
                }
                int[] tmpValues = anIsCountFingerprint
                        ? tmpFingerprinter.getCountArray(tmpFrequencies)
                        : tmpFingerprinter.getBitArray(tmpFrequencies);
                for (int j = 0; j < tmpColumnCount; j++) {
                    tmpData[i][j] = tmpValues[j];
                }
            }
        }
        FeatureType[] tmpTypes = new FeatureType[tmpColumnCount];
        Arrays.fill(tmpTypes, anIsCountFingerprint ? FeatureType.COUNT : FeatureType.BIT);
        String[] tmpGroupNames = new String[tmpColumnCount];
        Arrays.fill(tmpGroupNames, aFragmentationName);
        return new FeatureMatrix(tmpData, FeatureMatrixBuilder.getRowNames(aMoleculeList),
                tmpDictionary.toArray(new String[0]), tmpTypes, tmpGroupNames);
    }
    //
    /**
     * Builds the descriptor block for the given molecules. Every molecule is copied, its atom types are perceived and
     * the given aromaticity model is applied to the copy before the descriptors are calculated, so the atom containers
     * of the data models are never changed. Molecules whose preparation fails get a row of NaN values; failed descriptor
     * calculations also result in NaN values.
     * The columns get the component names of the CDK as feature names (e.g. "MW", "WPATH" or "ECFP4_1023"), the
     * component types of the calculator as feature types (BINARY: BIT, COUNT: COUNT, CONTINUOUS: CONTINUOUS, see
     * {@link #toFeatureType(DescriptorComponentType)}), so that every column gets the transformation of its type, and
     * the descriptor enum constant names as feature group names. Exception: the Wiener polarity number ("WPOL") of
     * {@link Descriptor#WIENER_NUMBER} gets the type COUNT, while the Wiener path number stays CONTINUOUS.
     *
     * @param aMoleculeList molecules (rows of the matrix)
     * @param aDescriptorArray descriptors to calculate
     * @param anAromaticityModel aromaticity model to apply, e.g. created via
     *                           {@link de.unijena.cheminf.mortar.model.util.ChemUtil#createAromaticityModel}
     * @param aDescriptorCalculator calculator, which defines e.g. the circular fingerprint size, see
     *                              {@link #getDescriptorCalculator(int)}
     * @param anIsParallelCalculation true: molecules are processed in parallel
     * @return feature matrix of the descriptor block
     * @throws NullPointerException if a parameter is null
     * @throws InterruptedException if the calculation is interrupted
     */
    public static FeatureMatrix buildDescriptorMatrix(
            List<MoleculeDataModel> aMoleculeList,
            Descriptor[] aDescriptorArray,
            Aromaticity anAromaticityModel,
            DescriptorCalculator aDescriptorCalculator,
            boolean anIsParallelCalculation)
            throws NullPointerException, InterruptedException {
        Objects.requireNonNull(aMoleculeList, "Given molecule list is null.");
        Objects.requireNonNull(aDescriptorArray, "Given descriptor array is null.");
        Objects.requireNonNull(anAromaticityModel, "Given aromaticity model is null.");
        Objects.requireNonNull(aDescriptorCalculator, "Given descriptor calculator is null.");
        int tmpColumnCount = aDescriptorCalculator.getNumberOfComponents(aDescriptorArray);
        float[][] tmpData = new float[aMoleculeList.size()][tmpColumnCount];
        //prepare copies of the molecules, failed molecules are excluded from the calculation
        List<IAtomContainer> tmpPreparedList = new ArrayList<>(aMoleculeList.size());
        List<Integer> tmpPreparedRowIndices = new ArrayList<>(aMoleculeList.size());
        for (int i = 0; i < aMoleculeList.size(); i++) {
            try {
                IAtomContainer tmpCopy = aMoleculeList.get(i).getAtomContainer().clone();
                FeatureMatrixBuilder.prepareMolecule(tmpCopy, anAromaticityModel);
                tmpPreparedList.add(tmpCopy);
                tmpPreparedRowIndices.add(i);
            } catch (Exception anException) {
                FeatureMatrixBuilder.LOGGER.log(Level.WARNING, String.format(
                        "Molecule %s could not be prepared for descriptor calculation.", aMoleculeList.get(i).getName()), anException);
                Arrays.fill(tmpData[i], Float.NaN);
            }
        }
        if (tmpColumnCount > 0 && !tmpPreparedList.isEmpty()) {
            float[][] tmpPreparedData = new float[tmpPreparedList.size()][tmpColumnCount];
            //NaN positions are tracked by the engine, the NaN values themselves are kept in the matrix
            List<int[]> tmpNaNPositionList = new ArrayList<>();
            aDescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelization(aDescriptorArray,
                    tmpPreparedList.toArray(new IAtomContainer[0]), tmpPreparedData, 0, anIsParallelCalculation, tmpNaNPositionList);
            for (int k = 0; k < tmpPreparedRowIndices.size(); k++) {
                tmpData[tmpPreparedRowIndices.get(k)] = tmpPreparedData[k];
            }
        }
        //CDK component names and component types, one per column
        String[] tmpNames = aDescriptorCalculator.getColumnNames(aDescriptorArray);
        DescriptorComponentType[] tmpComponentTypes = aDescriptorCalculator.getComponentTypes(aDescriptorArray);
        FeatureType[] tmpTypes = new FeatureType[tmpColumnCount];
        String[] tmpGroupNames = new String[tmpColumnCount];
        int tmpColumnIndex = 0;
        for (Descriptor tmpDescriptor : aDescriptorArray) {
            int tmpComponentNumber = aDescriptorCalculator.getNumberOfComponents(tmpDescriptor);
            for (int c = 0; c < tmpComponentNumber; c++) {
                //the Wiener polarity number is a count, unlike the Wiener path number, see the constant
                tmpTypes[tmpColumnIndex] = tmpDescriptor == Descriptor.WIENER_NUMBER
                        && c == FeatureMatrixBuilder.WIENER_POLARITY_NUMBER_COMPONENT_INDEX
                        ? FeatureType.COUNT
                        : FeatureMatrixBuilder.toFeatureType(tmpComponentTypes[tmpColumnIndex]);
                tmpGroupNames[tmpColumnIndex] = tmpDescriptor.name();
                tmpColumnIndex++;
            }
        }
        return new FeatureMatrix(tmpData, FeatureMatrixBuilder.getRowNames(aMoleculeList), tmpNames, tmpTypes, tmpGroupNames);
    }
    //
    /**
     * Builds the raw feature vector of a single molecule with exactly the columns of the raw matrix of the given clustering
     * run: the fragment columns (presence or frequency of the run's fragments in the molecule, taken from the molecule's
     * fragment frequencies of the run's fragmentation) followed by the run's descriptors calculated on the molecule with
     * the run's aromaticity model. Used to search molecules similar to a molecule that was not clustered.
     *
     * @param aMolecule molecule; its fragment frequencies for the run's fragmentation must be set (a molecule that was
     *                  not fragmented gets zeros for all fragment columns)
     * @param aRunResult clustering run that defines the columns
     * @return raw feature vector, may contain NaN for failed descriptor calculations
     * @throws NullPointerException if a parameter is null
     * @throws InterruptedException if the descriptor calculation is interrupted
     */
    public static float[] buildRawFeatureVector(MoleculeDataModel aMolecule, ClusteringRunResult aRunResult)
            throws NullPointerException, InterruptedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        ClusteringParameters tmpParameters = aRunResult.parameters();
        FeatureMatrix tmpRaw = aRunResult.rawMatrix();
        float[] tmpVector = new float[tmpRaw.getColumnCount()];
        //the fragment columns come first and belong to the group of the fragmentation name
        int tmpFragmentColumnCount = 0;
        while (tmpFragmentColumnCount < tmpRaw.getColumnCount()
                && tmpRaw.getFeatureGroupName(tmpFragmentColumnCount).equals(tmpParameters.fragmentationName())
                && !FeatureMatrixBuilder.isMolecularDescriptorColumn(tmpRaw, tmpFragmentColumnCount)) {
            tmpFragmentColumnCount++;
        }
        Map<String, Integer> tmpFrequencies = FeatureMatrixBuilder.getFragmentFrequencies(aMolecule, tmpParameters.fragmentationName());
        for (int j = 0; j < tmpFragmentColumnCount; j++) {
            Integer tmpFrequency = tmpFrequencies.get(tmpRaw.getFeatureName(j));
            int tmpCount = tmpFrequency == null ? 0 : tmpFrequency;
            tmpVector[j] = tmpParameters.isCountFingerprint() ? tmpCount : (tmpCount > 0 ? 1.0f : 0.0f);
        }
        //descriptors with the circular fingerprint size of the run
        FeatureMatrix tmpDescriptors = FeatureMatrixBuilder.buildDescriptorMatrix(List.of(aMolecule), tmpParameters.descriptors(),
                ChemUtil.createAromaticityModel(tmpParameters.electronDonationModelOption(), tmpParameters.cycleFinderOption()),
                FeatureMatrixBuilder.getDescriptorCalculator(tmpParameters.circularFingerprintSize()), false);
        if (tmpFragmentColumnCount + tmpDescriptors.getColumnCount() != tmpVector.length) {
            throw new IllegalStateException("The features of the molecule do not match the columns of the clustering run.");
        }
        for (int j = 0; j < tmpDescriptors.getColumnCount(); j++) {
            tmpVector[tmpFragmentColumnCount + j] = tmpDescriptors.getValue(0, j);
        }
        return tmpVector;
    }
    //
    /**
     * Prepares a molecule for descriptor calculation: perceives atom types, clears existing aromaticity flags, marks
     * ring atoms and bonds and applies the given aromaticity model. The given molecule is changed.
     *
     * @param aMolecule molecule to prepare (IS CHANGED)
     * @param anAromaticityModel aromaticity model to apply
     * @throws Exception if atom type perception or aromaticity detection fails
     */
    public static void prepareMolecule(IAtomContainer aMolecule, Aromaticity anAromaticityModel) throws Exception {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        Aromaticity.clear(aMolecule);
        Cycles.markRingAtomsAndBonds(aMolecule);
        anAromaticityModel.apply(aMolecule);
    }
    //
    /**
     * Returns a descriptor calculator with the given circular fingerprint size; circular fingerprints are calculated as
     * bit fingerprints. The calculator of the last request is reused if it has the same size, otherwise a new one is
     * created and kept.
     *
     * @param aCircularFingerprintSize size (number of components) of circular fingerprints
     * @return descriptor calculator
     * @throws IllegalArgumentException if the size is not positive
     */
    public static DescriptorCalculator getDescriptorCalculator(int aCircularFingerprintSize) throws IllegalArgumentException {
        DescriptorCalculator tmpCalculator = FeatureMatrixBuilder.descriptorCalculator;
        if (tmpCalculator == null || tmpCalculator.getCircularFingerprintSize() != aCircularFingerprintSize) {
            tmpCalculator = new DescriptorCalculator(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, aCircularFingerprintSize, false);
            FeatureMatrixBuilder.descriptorCalculator = tmpCalculator;
        }
        return tmpCalculator;
    }
    //
    /**
     * Returns the feature type of the given descriptor component type: BINARY gives BIT, COUNT gives COUNT and
     * CONTINUOUS gives CONTINUOUS.
     *
     * @param aComponentType component type of a descriptor
     * @return feature type
     * @throws NullPointerException if the component type is null
     */
    public static FeatureType toFeatureType(DescriptorComponentType aComponentType) throws NullPointerException {
        Objects.requireNonNull(aComponentType, "Given component type is null.");
        return switch (aComponentType) {
            case BINARY -> FeatureType.BIT;
            case COUNT -> FeatureType.COUNT;
            case CONTINUOUS -> FeatureType.CONTINUOUS;
        };
    }
    //
    /**
     * Returns whether the given column of the given matrix is a component of a molecular descriptor, i.e. of a
     * descriptor that is no fingerprint (e.g. "MW" or "nC"), regardless of its type (COUNT or CONTINUOUS). Fragment
     * and fingerprint columns are no descriptor columns.
     *
     * @param aMatrix feature matrix whose group names are descriptor enum constant names or fragmentation names
     * @param aColumnIndex column index
     * @return true if the column belongs to a molecular descriptor
     */
    public static boolean isMolecularDescriptorColumn(FeatureMatrix aMatrix, int aColumnIndex) {
        Descriptor tmpDescriptor = FeatureMatrixBuilder.DESCRIPTOR_OF_GROUP_NAME_MAP.get(aMatrix.getFeatureGroupName(aColumnIndex));
        return tmpDescriptor != null && !tmpDescriptor.isFingerprint();
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Creates the map of the descriptors by enum constant name.
     *
     * @return unmodifiable map from enum constant name to descriptor
     */
    private static Map<String, Descriptor> createDescriptorOfGroupNameMap() {
        Map<String, Descriptor> tmpMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(Descriptor.values().length));
        for (Descriptor tmpDescriptor : Descriptor.values()) {
            tmpMap.put(tmpDescriptor.name(), tmpDescriptor);
        }
        return Map.copyOf(tmpMap);
    }
    //
    /**
     * Returns the fragment frequencies of the given molecule for the given fragmentation or an empty map if the molecule
     * did not undergo this fragmentation.
     *
     * @param aMolecule molecule
     * @param aFragmentationName fragmentation name
     * @return map of fragment unique SMILES to frequency, never null
     */
    private static Map<String, Integer> getFragmentFrequencies(MoleculeDataModel aMolecule, String aFragmentationName) {
        Map<String, Integer> tmpFrequencies = aMolecule.getFragmentFrequencyOfSpecificFragmentation(aFragmentationName);
        return tmpFrequencies == null ? Map.of() : tmpFrequencies;
    }
    //
    /**
     * Returns the names of the given molecules as array; package-private, also used by {@link DescriptorMatrixCache}.
     *
     * @param aMoleculeList molecules
     * @return array of molecule names
     */
    static String[] getRowNames(List<MoleculeDataModel> aMoleculeList) {
        String[] tmpRowNames = new String[aMoleculeList.size()];
        for (int i = 0; i < aMoleculeList.size(); i++) {
            tmpRowNames[i] = aMoleculeList.get(i).getName();
        }
        return tmpRowNames;
    }
    //</editor-fold>
}
