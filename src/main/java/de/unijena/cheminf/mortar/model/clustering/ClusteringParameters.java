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

import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * Immutable snapshot of all parameters of one clustering run. It is created from the settings of the
 * {@link ClusteringService} when a run is started, so that changing the settings during a run has no effect, and it is
 * written to the export to make the run reproducible.
 *
 * @param fragmentationName name of the fragmentation whose fragments form the fragment fingerprint
 * @param isCountFingerprint true: fragment counts, false: fragment presence bits
 * @param minimumFragmentMoleculeFrequency minimum number of molecules a fragment must occur in to become a feature
 * @param maximumNumberOfFragmentFeatures maximum number of fragment features (the most frequent fragments are kept), 0 = no limit
 * @param descriptors descriptors calculated on the molecules, may be empty
 * @param electronDonationModelOption electron donation model of the aromaticity model
 * @param cycleFinderOption cycle finder of the aromaticity model
 * @param countTransformation transformation of COUNT features
 * @param continuousTransformation transformation of CONTINUOUS features
 * @param maximumNaNFraction maximum fraction of NaN values per feature
 * @param isCorrelationFilterEnabled true if the correlation filter is applied
 * @param correlationThreshold absolute correlation above which features are regarded as redundant (only used if the
 *                             correlation filter is enabled)
 * @param correlationMethod correlation coefficient of the correlation filter
 * @param algorithm ART-2A kernel
 * @param vigilance vigilance parameter
 * @param maximumNumberOfClusters maximum number of clusters
 * @param maximumNumberOfEpochs maximum number of epochs
 * @param convergenceThreshold convergence threshold of the kernel
 * @param learningParameter learning parameter
 * @param offsetForContrastEnhancement offset for the contrast enhancement threshold
 * @param randomSeed seed of the random number generator
 * @param circularFingerprintSize length of the circular fingerprint in bits (only relevant if one is selected)
 * @param isParallelCalculation true: descriptor calculation and winner search are parallelized
 * @param isFragmentFingerprintUsed true: the fragment fingerprint is part of the features; false: only the descriptors
 *                                  and fingerprints calculated on the molecules are used (the fragment settings are
 *                                  ignored)
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record ClusteringParameters(
        String fragmentationName,
        boolean isCountFingerprint,
        int minimumFragmentMoleculeFrequency,
        int maximumNumberOfFragmentFeatures,
        Descriptor[] descriptors,
        IMoleculeFragmenter.ElectronDonationModelOption electronDonationModelOption,
        IMoleculeFragmenter.CycleFinderOption cycleFinderOption,
        FeatureTransformation countTransformation,
        FeatureTransformation continuousTransformation,
        double maximumNaNFraction,
        boolean isCorrelationFilterEnabled,
        double correlationThreshold,
        FeatureFilter.CorrelationMethod correlationMethod,
        Art2aClusterer.Art2aAlgorithm algorithm,
        float vigilance,
        int maximumNumberOfClusters,
        int maximumNumberOfEpochs,
        float convergenceThreshold,
        float learningParameter,
        float offsetForContrastEnhancement,
        long randomSeed,
        int circularFingerprintSize,
        boolean isParallelCalculation,
        boolean isFragmentFingerprintUsed) {
    /**
     * Compact constructor that checks the parameters for null and copies the descriptor array.
     *
     * @throws NullPointerException if an object parameter is null
     */
    public ClusteringParameters {
        Objects.requireNonNull(fragmentationName, "Given fragmentation name is null.");
        Objects.requireNonNull(descriptors, "Given descriptor array is null.");
        Objects.requireNonNull(electronDonationModelOption, "Given electron donation model option is null.");
        Objects.requireNonNull(cycleFinderOption, "Given cycle finder option is null.");
        Objects.requireNonNull(countTransformation, "Given count transformation is null.");
        Objects.requireNonNull(continuousTransformation, "Given continuous transformation is null.");
        Objects.requireNonNull(correlationMethod, "Given correlation method is null.");
        Objects.requireNonNull(algorithm, "Given algorithm is null.");
        descriptors = descriptors.clone();
    }
    //
    /**
     * Constructor with the fragment fingerprint as part of the features, see the record description for the parameters.
     *
     * @param fragmentationName name of the fragmentation whose fragments form the fragment fingerprint
     * @param isCountFingerprint true: fragment counts, false: fragment presence bits
     * @param minimumFragmentMoleculeFrequency minimum number of molecules a fragment must occur in to become a feature
     * @param maximumNumberOfFragmentFeatures maximum number of fragment features, 0 = no limit
     * @param descriptors descriptors calculated on the molecules, may be empty
     * @param electronDonationModelOption electron donation model of the aromaticity model
     * @param cycleFinderOption cycle finder of the aromaticity model
     * @param countTransformation transformation of COUNT features
     * @param continuousTransformation transformation of CONTINUOUS features
     * @param maximumNaNFraction maximum fraction of NaN values per feature
     * @param isCorrelationFilterEnabled true if the correlation filter is applied
     * @param correlationThreshold absolute correlation above which features are regarded as redundant
     * @param correlationMethod correlation coefficient of the correlation filter
     * @param algorithm ART-2A kernel
     * @param vigilance vigilance parameter
     * @param maximumNumberOfClusters maximum number of clusters
     * @param maximumNumberOfEpochs maximum number of epochs
     * @param convergenceThreshold convergence threshold of the kernel
     * @param learningParameter learning parameter
     * @param offsetForContrastEnhancement offset for the contrast enhancement threshold
     * @param randomSeed seed of the random number generator
     * @param circularFingerprintSize length of the circular fingerprint in bits
     * @param isParallelCalculation true: descriptor calculation and winner search are parallelized
     * @throws NullPointerException if an object parameter is null
     */
    public ClusteringParameters(String fragmentationName, boolean isCountFingerprint, int minimumFragmentMoleculeFrequency,
                                int maximumNumberOfFragmentFeatures, Descriptor[] descriptors,
                                IMoleculeFragmenter.ElectronDonationModelOption electronDonationModelOption,
                                IMoleculeFragmenter.CycleFinderOption cycleFinderOption, FeatureTransformation countTransformation,
                                FeatureTransformation continuousTransformation, double maximumNaNFraction, boolean isCorrelationFilterEnabled,
                                double correlationThreshold, FeatureFilter.CorrelationMethod correlationMethod,
                                Art2aClusterer.Art2aAlgorithm algorithm, float vigilance, int maximumNumberOfClusters,
                                int maximumNumberOfEpochs, float convergenceThreshold, float learningParameter,
                                float offsetForContrastEnhancement, long randomSeed, int circularFingerprintSize,
                                boolean isParallelCalculation) throws NullPointerException {
        this(fragmentationName, isCountFingerprint, minimumFragmentMoleculeFrequency, maximumNumberOfFragmentFeatures, descriptors,
                electronDonationModelOption, cycleFinderOption, countTransformation, continuousTransformation, maximumNaNFraction,
                isCorrelationFilterEnabled, correlationThreshold, correlationMethod, algorithm, vigilance, maximumNumberOfClusters,
                maximumNumberOfEpochs, convergenceThreshold, learningParameter, offsetForContrastEnhancement, randomSeed,
                circularFingerprintSize, isParallelCalculation, true);
    }
    //
    /**
     * Returns a copy of the descriptor array.
     *
     * @return descriptors calculated on the molecules
     */
    @Override
    public Descriptor[] descriptors() {
        return this.descriptors.clone();
    }
    //
    /**
     * Returns all parameters as "name=value" lines, e.g. for the parameter file of the export.
     *
     * @return list of parameter lines
     */
    public List<String> getParameterLines() {
        List<String> tmpLines = new ArrayList<>(20);
        tmpLines.add("fragmentationName=" + this.fragmentationName);
        tmpLines.add("isFragmentFingerprintUsed=" + this.isFragmentFingerprintUsed);
        tmpLines.add("isCountFingerprint=" + this.isCountFingerprint);
        tmpLines.add("minimumFragmentMoleculeFrequency=" + this.minimumFragmentMoleculeFrequency);
        tmpLines.add("maximumNumberOfFragmentFeatures=" + this.maximumNumberOfFragmentFeatures);
        tmpLines.add("descriptors=" + String.join(",", Arrays.stream(this.descriptors).map(Descriptor::name).toList()));
        tmpLines.add("electronDonationModel=" + this.electronDonationModelOption.name());
        tmpLines.add("cycleFinder=" + this.cycleFinderOption.name());
        tmpLines.add("countTransformation=" + this.countTransformation.name());
        tmpLines.add("continuousTransformation=" + this.continuousTransformation.name());
        tmpLines.add("maximumNaNFraction=" + this.maximumNaNFraction);
        tmpLines.add("isCorrelationFilterEnabled=" + this.isCorrelationFilterEnabled);
        tmpLines.add("correlationThreshold=" + this.correlationThreshold);
        tmpLines.add("correlationMethod=" + this.correlationMethod.name());
        tmpLines.add("algorithm=" + this.algorithm.name());
        tmpLines.add("vigilance=" + this.vigilance);
        tmpLines.add("maximumNumberOfClusters=" + this.maximumNumberOfClusters);
        tmpLines.add("maximumNumberOfEpochs=" + this.maximumNumberOfEpochs);
        tmpLines.add("convergenceThreshold=" + this.convergenceThreshold);
        tmpLines.add("learningParameter=" + this.learningParameter);
        tmpLines.add("offsetForContrastEnhancement=" + this.offsetForContrastEnhancement);
        tmpLines.add("randomSeed=" + this.randomSeed);
        tmpLines.add("circularFingerprintSize=" + this.circularFingerprintSize);
        return tmpLines;
    }
}
