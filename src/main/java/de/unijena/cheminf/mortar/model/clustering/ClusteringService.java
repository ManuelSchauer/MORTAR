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

import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorCalculator;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.model.util.FileUtil;
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;
import de.unijena.cheminf.mortar.model.util.SimpleIDisplayEnumConstantProperty;
import de.unijena.cheminf.mortar.preference.PreferenceContainer;
import de.unijena.cheminf.mortar.preference.PreferenceUtil;
import de.unijena.cheminf.mortar.preference.SingleTermPreference;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;

import org.openscience.cdk.aromaticity.Aromaticity;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Properties;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Service for the clustering of molecules by their fragment fingerprint and descriptors. It holds all clustering
 * settings as JavaFX properties (including the selection of the descriptors), so that they can be displayed and edited
 * in the {@link de.unijena.cheminf.mortar.gui.views.SettingsView} like the fragmenter settings, persists and reloads
 * them, and creates {@link ClusteringTask}s with an immutable snapshot of the settings.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringService {
    //<editor-fold desc="Enum FragmentFingerprintTypeOption" defaultstate="collapsed">
    /**
     * Enum for the type of the fragment fingerprint.
     */
    public enum FragmentFingerprintTypeOption implements IDisplayEnum {
        /**
         * Presence bits of the fragments.
         */
        BIT(Message.get("ClusteringService.FragmentFingerprintTypeOption.bit.displayName"),
                Message.get("ClusteringService.FragmentFingerprintTypeOption.bit.tooltip")),
        /**
         * Frequencies of the fragments in the molecule.
         */
        COUNT(Message.get("ClusteringService.FragmentFingerprintTypeOption.count.displayName"),
                Message.get("ClusteringService.FragmentFingerprintTypeOption.count.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor setting the display name and tooltip.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        FragmentFingerprintTypeOption(String aDisplayName, String aTooltip) {
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum CountTransformationOption" defaultstate="collapsed">
    /**
     * Enum for the transformations available for COUNT features.
     */
    public enum CountTransformationOption implements IDisplayEnum {
        /**
         * No transformation.
         */
        NONE(FeatureTransformation.NONE, Message.get("ClusteringService.TransformationOption.none.displayName"),
                Message.get("ClusteringService.TransformationOption.none.tooltip")),
        /**
         * ln(1 + x).
         */
        LOG1P(FeatureTransformation.LOG1P, Message.get("ClusteringService.TransformationOption.log1p.displayName"),
                Message.get("ClusteringService.TransformationOption.log1p.tooltip")),
        /**
         * sqrt(x).
         */
        SQRT(FeatureTransformation.SQRT, Message.get("ClusteringService.TransformationOption.sqrt.displayName"),
                Message.get("ClusteringService.TransformationOption.sqrt.tooltip")),
        /**
         * Presence bit.
         */
        BINARIZE(FeatureTransformation.BINARIZE, Message.get("ClusteringService.TransformationOption.binarize.displayName"),
                Message.get("ClusteringService.TransformationOption.binarize.tooltip"));
        /**
         * Associated feature transformation.
         */
        private final FeatureTransformation transformation;
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param aTransformation associated feature transformation
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        CountTransformationOption(FeatureTransformation aTransformation, String aDisplayName, String aTooltip) {
            this.transformation = aTransformation;
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        /**
         * Returns the associated feature transformation.
         *
         * @return feature transformation
         */
        public FeatureTransformation getTransformation() {
            return this.transformation;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum ContinuousTransformationOption" defaultstate="collapsed">
    /**
     * Enum for the transformations available for CONTINUOUS features.
     */
    public enum ContinuousTransformationOption implements IDisplayEnum {
        /**
         * No transformation.
         */
        NONE(FeatureTransformation.NONE, Message.get("ClusteringService.TransformationOption.none.displayName"),
                Message.get("ClusteringService.TransformationOption.none.tooltip")),
        /**
         * Clipping at the 1st and 99th percentile.
         */
        WINSORIZE(FeatureTransformation.WINSORIZE, Message.get("ClusteringService.TransformationOption.winsorize.displayName"),
                Message.get("ClusteringService.TransformationOption.winsorize.tooltip"));
        /**
         * Associated feature transformation.
         */
        private final FeatureTransformation transformation;
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param aTransformation associated feature transformation
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        ContinuousTransformationOption(FeatureTransformation aTransformation, String aDisplayName, String aTooltip) {
            this.transformation = aTransformation;
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        /**
         * Returns the associated feature transformation.
         *
         * @return feature transformation
         */
        public FeatureTransformation getTransformation() {
            return this.transformation;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum CorrelationMethodOption" defaultstate="collapsed">
    /**
     * Enum for the correlation coefficient of the correlation filter.
     */
    public enum CorrelationMethodOption implements IDisplayEnum {
        /**
         * Pearson correlation.
         */
        PEARSON(FeatureFilter.CorrelationMethod.PEARSON, Message.get("ClusteringService.CorrelationMethodOption.pearson.displayName"),
                Message.get("ClusteringService.CorrelationMethodOption.pearson.tooltip")),
        /**
         * Spearman rank correlation.
         */
        SPEARMAN(FeatureFilter.CorrelationMethod.SPEARMAN, Message.get("ClusteringService.CorrelationMethodOption.spearman.displayName"),
                Message.get("ClusteringService.CorrelationMethodOption.spearman.tooltip"));
        /**
         * Associated correlation method.
         */
        private final FeatureFilter.CorrelationMethod method;
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param aMethod associated correlation method
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        CorrelationMethodOption(FeatureFilter.CorrelationMethod aMethod, String aDisplayName, String aTooltip) {
            this.method = aMethod;
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        /**
         * Returns the associated correlation method.
         *
         * @return correlation method
         */
        public FeatureFilter.CorrelationMethod getMethod() {
            return this.method;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum AlgorithmOption" defaultstate="collapsed">
    /**
     * Enum for the ART-2A kernel.
     */
    public enum AlgorithmOption implements IDisplayEnum {
        /**
         * Classic ART-2A with cosine similarity.
         */
        ART_2A_COSINE(Art2aClusterer.Art2aAlgorithm.COSINE, Art2aClusterer.CONVERGENCE_THRESHOLD_COSINE_DEFAULT,
                Message.get("ClusteringService.AlgorithmOption.cosine.displayName"),
                Message.get("ClusteringService.AlgorithmOption.cosine.tooltip")),
        /**
         * ART-2A-Euclid with squared Euclidean distance.
         */
        ART_2A_EUCLID(Art2aClusterer.Art2aAlgorithm.EUCLID, Art2aClusterer.CONVERGENCE_THRESHOLD_EUCLID_DEFAULT,
                Message.get("ClusteringService.AlgorithmOption.euclid.displayName"),
                Message.get("ClusteringService.AlgorithmOption.euclid.tooltip"));
        /**
         * Associated ART-2A algorithm.
         */
        private final Art2aClusterer.Art2aAlgorithm algorithm;
        /**
         * Default convergence threshold of the kernel.
         */
        private final double convergenceThresholdDefault;
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor.
         *
         * @param anAlgorithm associated ART-2A algorithm
         * @param aConvergenceThresholdDefault default convergence threshold of the kernel
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        AlgorithmOption(Art2aClusterer.Art2aAlgorithm anAlgorithm, float aConvergenceThresholdDefault,
                        String aDisplayName, String aTooltip) {
            this.algorithm = anAlgorithm;
            //converted via the decimal representation, so that e.g. 0.99f is shown as 0.99 and not as 0.9900000095367432
            this.convergenceThresholdDefault = ClusteringService.toDecimalDouble(aConvergenceThresholdDefault);
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        /**
         * Returns the associated ART-2A algorithm.
         *
         * @return ART-2A algorithm
         */
        public Art2aClusterer.Art2aAlgorithm getAlgorithm() {
            return this.algorithm;
        }
        /**
         * Returns the default convergence threshold of the kernel.
         *
         * @return default convergence threshold
         */
        public double getConvergenceThresholdDefault() {
            return this.convergenceThresholdDefault;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum DescriptorPresetOption" defaultstate="collapsed">
    /**
     * Enum for the built-in descriptor presets. Choosing a preset sets the check boxes of the descriptor list; changing a
     * check box afterwards switches the selection to custom, which can be saved as a user-defined preset. Fingerprints
     * are not part of the presets, they are chosen separately.
     */
    public enum DescriptorPresetOption implements IDisplayEnum {
        /**
         * No descriptors, the molecules are clustered by their fragment fingerprint (and selected fingerprints) only.
         */
        NONE(Message.get("ClusteringService.DescriptorPresetOption.none.displayName"),
                Message.get("ClusteringService.DescriptorPresetOption.none.tooltip")),
        /**
         * A small set of well interpretable physicochemical descriptors, see {@link #BASIC_DESCRIPTORS}.
         */
        BASIC(Message.get("ClusteringService.DescriptorPresetOption.basic.displayName"),
                Message.get("ClusteringService.DescriptorPresetOption.basic.tooltip")),
        /**
         * All fast and safe descriptors that do not need 3D coordinates.
         */
        FAST(Message.get("ClusteringService.DescriptorPresetOption.fast.displayName"),
                Message.get("ClusteringService.DescriptorPresetOption.fast.tooltip")),
        /**
         * All descriptors that do not need 3D coordinates, including slow and unsafe ones.
         */
        ALL_2D(Message.get("ClusteringService.DescriptorPresetOption.all2d.displayName"),
                Message.get("ClusteringService.DescriptorPresetOption.all2d.tooltip")),
        /**
         * All descriptors including the ones that need 3D coordinates (NaN for molecules without 3D coordinates).
         */
        ALL(Message.get("ClusteringService.DescriptorPresetOption.all.displayName"),
                Message.get("ClusteringService.DescriptorPresetOption.all.tooltip"));
        /**
         * Descriptors of the BASIC preset: size, lipophilicity, polarity, hydrogen bonding, flexibility, aromaticity,
         * saturation and complexity (proposal, to be agreed on).
         */
        public static final Set<Descriptor> BASIC_DESCRIPTORS = Collections.unmodifiableSet(EnumSet.of(
                Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT_HEAVY, Descriptor.X_LOG_P, Descriptor.TPSA,
                Descriptor.H_BOND_ACCEPTOR_COUNT, Descriptor.H_BOND_DONOR_COUNT, Descriptor.ROTATABLE_BONDS_COUNT,
                Descriptor.AROMATIC_ATOMS_COUNT, Descriptor.FRACTIONAL_CSP3, Descriptor.FRAGMENT_COMPLEXITY));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor setting the display name and tooltip.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        DescriptorPresetOption(String aDisplayName, String aTooltip) {
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        /**
         * Returns whether the given descriptor belongs to this preset. Fingerprints never belong to a preset.
         *
         * @param aDescriptor descriptor
         * @return true if the descriptor belongs to the preset
         */
        public boolean includes(Descriptor aDescriptor) {
            if (aDescriptor.isFingerprint()) {
                return false;
            }
            return switch (this) {
                case NONE -> false;
                case BASIC -> DescriptorPresetOption.BASIC_DESCRIPTORS.contains(aDescriptor);
                case FAST -> aDescriptor.isFast() && aDescriptor.isSafe() && !aDescriptor.requires3DCoordinates();
                case ALL_2D -> !aDescriptor.requires3DCoordinates();
                case ALL -> true;
            };
        }
        /**
         * Returns this built-in preset as {@link DescriptorPreset}; its key is the name of the enum constant.
         *
         * @return descriptor preset
         */
        public DescriptorPreset toDescriptorPreset() {
            Set<Descriptor> tmpDescriptors = EnumSet.noneOf(Descriptor.class);
            for (Descriptor tmpDescriptor : Descriptor.getAllDescriptors()) {
                if (this.includes(tmpDescriptor)) {
                    tmpDescriptors.add(tmpDescriptor);
                }
            }
            return new DescriptorPreset(this.name(), this.displayName, this.tooltip, tmpDescriptors, false);
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum CircularFingerprintTypeOption" defaultstate="collapsed">
    /**
     * Enum for the type of the circular fingerprint.
     */
    public enum CircularFingerprintTypeOption implements IDisplayEnum {
        /**
         * Extended-connectivity fingerprint, atoms are described by their element and connectivity.
         */
        ECFP(Message.get("ClusteringService.CircularFingerprintTypeOption.ecfp.displayName"),
                Message.get("ClusteringService.CircularFingerprintTypeOption.ecfp.tooltip")),
        /**
         * Functional-class fingerprint, atoms are described by their pharmacophoric role.
         */
        FCFP(Message.get("ClusteringService.CircularFingerprintTypeOption.fcfp.displayName"),
                Message.get("ClusteringService.CircularFingerprintTypeOption.fcfp.tooltip"));
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Language-specific tooltip text for display in GUI.
         */
        private final String tooltip;
        /**
         * Constructor setting the display name and tooltip.
         *
         * @param aDisplayName display name
         * @param aTooltip tooltip text
         */
        CircularFingerprintTypeOption(String aDisplayName, String aTooltip) {
            this.displayName = aDisplayName;
            this.tooltip = aTooltip;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return this.tooltip;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Enum CircularFingerprintRadiusOption" defaultstate="collapsed">
    /**
     * Enum for the radius of the circular fingerprint. The number in the common names (ECFP4, FCFP6, ...) is the
     * diameter, i.e. twice the radius.
     */
    public enum CircularFingerprintRadiusOption implements IDisplayEnum {
        /**
         * Radius 0 (ECFP0/FCFP0): only the atoms themselves.
         */
        RADIUS_0(0, Message.get("ClusteringService.CircularFingerprintRadiusOption.radius0.displayName")),
        /**
         * Radius 1 (ECFP2/FCFP2): atoms and their direct neighbours.
         */
        RADIUS_1(1, Message.get("ClusteringService.CircularFingerprintRadiusOption.radius1.displayName")),
        /**
         * Radius 2 (ECFP4/FCFP4): environments up to two bonds, the most common choice.
         */
        RADIUS_2(2, Message.get("ClusteringService.CircularFingerprintRadiusOption.radius2.displayName")),
        /**
         * Radius 3 (ECFP6/FCFP6): environments up to three bonds.
         */
        RADIUS_3(3, Message.get("ClusteringService.CircularFingerprintRadiusOption.radius3.displayName"));
        /**
         * Radius in bonds.
         */
        private final int radius;
        /**
         * Language-specific name for display in GUI.
         */
        private final String displayName;
        /**
         * Constructor.
         *
         * @param aRadius radius in bonds
         * @param aDisplayName display name
         */
        CircularFingerprintRadiusOption(int aRadius, String aDisplayName) {
            this.radius = aRadius;
            this.displayName = aDisplayName;
        }
        /**
         * Returns the radius in bonds.
         *
         * @return radius
         */
        public int getRadius() {
            return this.radius;
        }
        /**
         * Returns the diameter, i.e. the number in the names ECFP4, FCFP6 etc.
         *
         * @return diameter (twice the radius)
         */
        public int getDiameter() {
            return 2 * this.radius;
        }
        @Override
        public String getDisplayName() {
            return this.displayName;
        }
        @Override
        public String getTooltipText() {
            return Message.get("ClusteringService.CircularFingerprintRadiusOption.tooltip");
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Name of the subfolder of the settings directory in which the clustering settings are persisted.
     */
    public static final String CLUSTERING_SETTINGS_SUBFOLDER_NAME = "Clustering_Settings";
    /**
     * File name (without extension) of the persisted clustering settings.
     */
    public static final String CLUSTERING_SETTINGS_FILE_NAME = "ClusteringSettings";
    /**
     * File name (without extension) of the persisted descriptor selection.
     */
    public static final String DESCRIPTOR_SELECTION_FILE_NAME = "DescriptorSelection";
    /**
     * Default descriptor preset.
     */
    public static final DescriptorPresetOption DESCRIPTOR_PRESET_OPTION_DEFAULT = DescriptorPresetOption.FAST;
    /**
     * Key of the selected preset if the descriptor selection does not correspond to a preset.
     */
    public static final String CUSTOM_PRESET_KEY = "CUSTOM";
    /**
     * File name of the user-defined descriptor presets in the clustering settings directory.
     */
    public static final String USER_DESCRIPTOR_PRESETS_FILE_NAME = "DescriptorPresets.properties";
    /**
     * Default descriptor columns shown with their cluster medians in the clustering result tabs, as CDK component names
     * (molecular weight, XLogP and TPSA); only the ones that were calculated are shown.
     */
    public static final List<String> DISPLAYED_DESCRIPTOR_COLUMNS_DEFAULT = List.of("MW", "XLogP", "TopoPSA");
    /**
     * Maximum number of component names listed in the tooltip of a descriptor.
     */
    private static final int MAXIMUM_NUMBER_OF_COMPONENT_NAMES_IN_TOOLTIP = 10;
    /**
     * Separator of the displayed descriptor column names in their setting; also terminates the list, so that the setting
     * is never empty, which a preference does not allow.
     */
    public static final String DISPLAYED_DESCRIPTOR_COLUMNS_SEPARATOR = ";";
    /**
     * Default maximum number of fragment features, 0 means no limit.
     */
    public static final int MAXIMUM_NUMBER_OF_FRAGMENT_FEATURES_DEFAULT = 0;
    /**
     * Default of whether a circular fingerprint is calculated.
     */
    public static final boolean CIRCULAR_FINGERPRINT_DEFAULT = false;
    /**
     * Default type of the circular fingerprint.
     */
    public static final CircularFingerprintTypeOption CIRCULAR_FINGERPRINT_TYPE_OPTION_DEFAULT = CircularFingerprintTypeOption.ECFP;
    /**
     * Default radius of the circular fingerprint (ECFP4).
     */
    public static final CircularFingerprintRadiusOption CIRCULAR_FINGERPRINT_RADIUS_OPTION_DEFAULT = CircularFingerprintRadiusOption.RADIUS_2;
    /**
     * Default length of the circular fingerprint in bits.
     */
    public static final int CIRCULAR_FINGERPRINT_SIZE_DEFAULT = 1024;
    /**
     * Default of whether the MACCS fingerprint is calculated.
     */
    public static final boolean MACCS_FINGERPRINT_DEFAULT = false;
    /**
     * Default of whether the PubChem fingerprint is calculated.
     */
    public static final boolean PUBCHEM_FINGERPRINT_DEFAULT = false;
    /**
     * Default fragment fingerprint type.
     */
    public static final FragmentFingerprintTypeOption FRAGMENT_FINGERPRINT_TYPE_OPTION_DEFAULT = FragmentFingerprintTypeOption.COUNT;
    /**
     * Default minimum number of molecules a fragment must occur in to become a feature.
     */
    public static final int MINIMUM_FRAGMENT_MOLECULE_FREQUENCY_DEFAULT = 2;
    /**
     * Default electron donation model, MORTAR's standard model.
     */
    public static final IMoleculeFragmenter.ElectronDonationModelOption ELECTRON_DONATION_MODEL_OPTION_DEFAULT =
            IMoleculeFragmenter.ElectronDonationModelOption.DAYLIGHT;
    /**
     * Default cycle finder, MORTAR's standard cycle finder.
     */
    public static final IMoleculeFragmenter.CycleFinderOption CYCLE_FINDER_OPTION_DEFAULT =
            IMoleculeFragmenter.CycleFinderOption.CDK_AROMATIC_SET;
    /**
     * Default transformation of COUNT features.
     */
    public static final CountTransformationOption COUNT_TRANSFORMATION_OPTION_DEFAULT = CountTransformationOption.LOG1P;
    /**
     * Default transformation of CONTINUOUS features.
     */
    public static final ContinuousTransformationOption CONTINUOUS_TRANSFORMATION_OPTION_DEFAULT = ContinuousTransformationOption.WINSORIZE;
    /**
     * Default maximum fraction of NaN values per feature.
     */
    public static final double MAXIMUM_NAN_FRACTION_DEFAULT = 0.1;
    /**
     * Default number of similar molecules returned by the similarity search.
     */
    public static final int NUMBER_OF_SIMILAR_MOLECULES_DEFAULT = 10;
    /**
     * Default of whether the correlation filter is applied.
     */
    public static final boolean CORRELATION_FILTER_DEFAULT = false;
    /**
     * Default of whether the fragment fingerprint is part of the features.
     */
    public static final boolean FRAGMENT_FINGERPRINT_DEFAULT = true;
    /**
     * Name used instead of a fragmentation name for clusterings without fragmentation, i.e. only by the descriptors and
     * fingerprints calculated on the molecules (e.g. started from the molecules tab); appears in titles and file names.
     */
    public static final String WITHOUT_FRAGMENTATION_NAME = "Descriptors only";
    /**
     * Default of whether the selected descriptors are part of the features of a clustering.
     */
    public static final boolean DESCRIPTORS_DEFAULT = true;
    /**
     * Default of whether the selected fingerprints are part of the features of a clustering.
     */
    public static final boolean FINGERPRINTS_DEFAULT = true;
    /**
     * Default correlation threshold.
     */
    public static final double CORRELATION_THRESHOLD_DEFAULT = 0.95;
    /**
     * Default correlation method.
     */
    public static final CorrelationMethodOption CORRELATION_METHOD_OPTION_DEFAULT = CorrelationMethodOption.SPEARMAN;
    /**
     * Default ART-2A kernel.
     */
    public static final AlgorithmOption ALGORITHM_OPTION_DEFAULT = AlgorithmOption.ART_2A_COSINE;
    /**
     * Default learning parameter, the library default as exact decimal double.
     */
    public static final double LEARNING_PARAMETER_DEFAULT = ClusteringService.toDecimalDouble(Art2aClusterer.LEARNING_PARAMETER_DEFAULT);
    /**
     * Default offset for the contrast enhancement, the library default as exact decimal double.
     */
    public static final double OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT =
            ClusteringService.toDecimalDouble(Art2aClusterer.OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT);
    /**
     * Default vigilance parameter.
     */
    public static final double VIGILANCE_DEFAULT = 0.5;
    /**
     * Default maximum number of clusters.
     */
    public static final int MAXIMUM_NUMBER_OF_CLUSTERS_DEFAULT = 1000;
    /**
     * Smallest allowed value for parameters that must lie in the open interval (0,1).
     */
    public static final double OPEN_INTERVAL_MINIMUM = 1.0e-4;
    /**
     * Largest allowed value for parameters that must lie in the open interval (0,1).
     */
    public static final double OPEN_INTERVAL_MAXIMUM = 1.0 - 1.0e-4;
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusteringService.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * All general clustering settings in display order.
     */
    private final List<Property<?>> settings;
    /**
     * Settings of the fragment fingerprint section of the clustering tab.
     */
    private final List<Property<?>> fragmentFingerprintSettings;
    /**
     * Setting whether the fragment fingerprint is used, in its own list, so that the other fragment fingerprint settings
     * can be hidden while it is switched off.
     */
    private final List<Property<?>> fragmentFingerprintSwitchSettings;
    /**
     * Setting whether the fragment fingerprint is part of the features.
     */
    private final SimpleBooleanProperty fragmentFingerprintSetting;
    /**
     * Settings of the aromaticity section of the clustering tab.
     */
    private final List<Property<?>> aromaticitySettings;
    /**
     * Settings of the transformation section of the clustering tab.
     */
    private final List<Property<?>> transformationSettings;
    /**
     * Settings of the filter section of the clustering tab, including the correlation filter switch.
     */
    private final List<Property<?>> filterSettings;
    /**
     * Settings of the correlation filter, only relevant (and shown) if the correlation filter is switched on.
     */
    private final List<Property<?>> correlationSettings;
    /**
     * Settings of the ART-2A section of the clustering tab.
     */
    private final List<Property<?>> art2aSettings;
    /**
     * Setting whether the correlation filter is applied.
     */
    private final SimpleBooleanProperty correlationFilterSetting;
    /**
     * Settings of the similarity search section of the clustering tab.
     */
    private final List<Property<?>> similaritySearchSettings;
    /**
     * Setting for the number of similar molecules returned by the similarity search.
     */
    private final SimpleIntegerProperty numberOfSimilarMoleculesSetting;
    /**
     * Keys of all descriptor presets for the selection in the GUI (own presets first, then built-in, then custom); shared
     * by all preset combo boxes so that they are updated when presets are saved or deleted.
     */
    private final ObservableList<String> descriptorPresetKeyList;
    /**
     * Descriptor selection settings in display order.
     */
    private final List<Property<?>> descriptorSettings;
    /**
     * Map of descriptor to its selection setting.
     */
    private final Map<Descriptor, SimpleBooleanProperty> descriptorToSettingMap;
    /**
     * Map of setting name to display name.
     */
    private final Map<String, String> settingNameToDisplayNameMap;
    /**
     * Map of setting name to tooltip text.
     */
    private final Map<String, String> settingNameToTooltipTextMap;
    /**
     * Descriptor preset setting.
     */
    private final SimpleStringProperty descriptorPresetSetting;
    /**
     * Settings of the fingerprint section of the descriptor tab.
     */
    private final List<Property<?>> fingerprintSettings;
    /**
     * Setting whether the selected descriptors are used for the clustering, in its own list for the settings view.
     */
    private final List<Property<?>> descriptorSwitchSettings;
    /**
     * Setting whether the selected fingerprints are used for the clustering, in its own list for the settings view.
     */
    private final List<Property<?>> fingerprintSwitchSettings;
    /**
     * Setting whether the selected descriptors are part of the features of a clustering; the selection is kept while it
     * is switched off.
     */
    private final SimpleBooleanProperty descriptorsSetting;
    /**
     * Setting whether the selected fingerprints are part of the features of a clustering; the selection is kept while it
     * is switched off.
     */
    private final SimpleBooleanProperty fingerprintsSetting;
    /**
     * Advanced settings of the fingerprint section of the descriptor tab (hidden by default).
     */
    private final List<Property<?>> advancedFingerprintSettings;
    /**
     * Check box settings of the descriptor list in display order.
     */
    private final List<Property<?>> descriptorSelectionSettings;
    /**
     * User-defined descriptor presets, sorted by name (case-insensitive).
     */
    private final Map<String, Set<Descriptor>> userDescriptorPresetMap;
    /**
     * Maximum number of fragment features setting, 0 means no limit.
     */
    private final SimpleIntegerProperty maximumNumberOfFragmentFeaturesSetting;
    /**
     * Setting whether a circular fingerprint is calculated.
     */
    private final SimpleBooleanProperty circularFingerprintSetting;
    /**
     * Circular fingerprint type setting.
     */
    private final SimpleIDisplayEnumConstantProperty circularFingerprintTypeSetting;
    /**
     * Circular fingerprint radius setting.
     */
    private final SimpleIDisplayEnumConstantProperty circularFingerprintRadiusSetting;
    /**
     * Circular fingerprint size setting.
     */
    private final SimpleIntegerProperty circularFingerprintSizeSetting;
    /**
     * Setting whether the MACCS fingerprint is calculated.
     */
    private final SimpleBooleanProperty maccsFingerprintSetting;
    /**
     * Setting whether the PubChem fingerprint is calculated.
     */
    private final SimpleBooleanProperty pubChemFingerprintSetting;
    /**
     * Fragment fingerprint type setting.
     */
    private final SimpleIDisplayEnumConstantProperty fragmentFingerprintTypeSetting;
    /**
     * Minimum fragment molecule frequency setting.
     */
    private final SimpleIntegerProperty minimumFragmentMoleculeFrequencySetting;
    /**
     * Electron donation model setting of the aromaticity model.
     */
    private final SimpleIDisplayEnumConstantProperty electronDonationModelSetting;
    /**
     * Cycle finder setting of the aromaticity model.
     */
    private final SimpleIDisplayEnumConstantProperty cycleFinderSetting;
    /**
     * Count transformation setting.
     */
    private final SimpleIDisplayEnumConstantProperty countTransformationSetting;
    /**
     * Continuous transformation setting.
     */
    private final SimpleIDisplayEnumConstantProperty continuousTransformationSetting;
    /**
     * Maximum NaN fraction setting.
     */
    private final SimpleDoubleProperty maximumNaNFractionSetting;
    /**
     * Correlation threshold setting.
     */
    private final SimpleDoubleProperty correlationThresholdSetting;
    /**
     * Correlation method setting.
     */
    private final SimpleIDisplayEnumConstantProperty correlationMethodSetting;
    /**
     * Algorithm setting.
     */
    private final SimpleIDisplayEnumConstantProperty algorithmSetting;
    /**
     * Names of the descriptor columns shown with their cluster medians in the clustering result tabs, each followed by
     * the separator; not shown in the settings view, chosen in the result tabs.
     */
    private final SimpleStringProperty displayedDescriptorColumnsSetting;
    /**
     * Cache of calculated descriptor values, shared by the clusterings, vigilance scans and descriptor calculations, so
     * that descriptors are not calculated twice for the same molecules.
     */
    private final DescriptorMatrixCache descriptorMatrixCache = new DescriptorMatrixCache();
    /**
     * Vigilance setting.
     */
    private final SimpleDoubleProperty vigilanceSetting;
    /**
     * Maximum number of clusters setting.
     */
    private final SimpleIntegerProperty maximumNumberOfClustersSetting;
    /**
     * Maximum number of epochs setting.
     */
    private final SimpleIntegerProperty maximumNumberOfEpochsSetting;
    /**
     * Convergence threshold setting.
     */
    private final SimpleDoubleProperty convergenceThresholdSetting;
    /**
     * Learning parameter setting.
     */
    private final SimpleDoubleProperty learningParameterSetting;
    /**
     * Offset for contrast enhancement setting.
     */
    private final SimpleDoubleProperty offsetForContrastEnhancementSetting;
    /**
     * Random seed setting.
     */
    private final SimpleIntegerProperty randomSeedSetting;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Priority table for the correlation filter, loaded on first use.
     */
    private DescriptorPriorityTable priorityTable;
    /**
     * Directory in which the clustering settings and the user-defined presets are persisted, ending with a file separator.
     */
    private final String settingsDirectoryPath;
    /**
     * True while a preset is applied to the descriptor check boxes, so that the changes do not switch the preset to custom.
     */
    private boolean isPresetBeingApplied;
    /**
     * Display name of the preset that was applied last, used to propose a name when the modified selection is saved.
     */
    private String lastAppliedPresetDisplayName;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, creates all settings with their default values. The settings are persisted in the subfolder
     * {@link #CLUSTERING_SETTINGS_SUBFOLDER_NAME} of the MORTAR settings directory.
     */
    public ClusteringService() {
        this(FileUtil.getSettingsDirPath() + ClusteringService.CLUSTERING_SETTINGS_SUBFOLDER_NAME + File.separator);
    }
    //
    /**
     * Constructor, creates all settings with their default values. The settings are persisted in the given directory,
     * which is e.g. used by tests to not touch the MORTAR settings directory.
     *
     * @param aSettingsDirectoryPath directory for the persisted settings and presets
     * @throws NullPointerException if the path is null
     */
    public ClusteringService(String aSettingsDirectoryPath) throws NullPointerException {
        Objects.requireNonNull(aSettingsDirectoryPath, "Given settings directory path is null.");
        this.settingsDirectoryPath = aSettingsDirectoryPath.endsWith(File.separator)
                ? aSettingsDirectoryPath
                : aSettingsDirectoryPath + File.separator;
        this.settings = new ArrayList<>(24);
        this.fragmentFingerprintSettings = new ArrayList<>(3);
        this.fragmentFingerprintSwitchSettings = new ArrayList<>(1);
        this.aromaticitySettings = new ArrayList<>(2);
        this.transformationSettings = new ArrayList<>(2);
        this.filterSettings = new ArrayList<>(2);
        this.correlationSettings = new ArrayList<>(2);
        this.art2aSettings = new ArrayList<>(8);
        this.similaritySearchSettings = new ArrayList<>(1);
        this.descriptorSettings = new ArrayList<>(Descriptor.getAllDescriptors().length);
        this.settingNameToDisplayNameMap = new HashMap<>(64);
        this.settingNameToTooltipTextMap = new HashMap<>(64);
        //fragment fingerprint
        this.fragmentFingerprintSetting = this.createBooleanSetting(this.fragmentFingerprintSwitchSettings, "Fragment fingerprint setting",
                ClusteringService.FRAGMENT_FINGERPRINT_DEFAULT, "fragmentFingerprintSetting");
        this.fragmentFingerprintTypeSetting = this.createEnumSetting(this.fragmentFingerprintSettings, "Fragment fingerprint type setting",
                ClusteringService.FRAGMENT_FINGERPRINT_TYPE_OPTION_DEFAULT, FragmentFingerprintTypeOption.class,
                "fragmentFingerprintTypeSetting");
        this.minimumFragmentMoleculeFrequencySetting = this.createIntegerSetting(this.fragmentFingerprintSettings, "Minimum fragment molecule frequency setting",
                ClusteringService.MINIMUM_FRAGMENT_MOLECULE_FREQUENCY_DEFAULT, 1, Integer.MAX_VALUE,
                "minimumFragmentMoleculeFrequencySetting");
        this.maximumNumberOfFragmentFeaturesSetting = this.createIntegerSetting(this.fragmentFingerprintSettings, "Maximum number of fragment features setting",
                ClusteringService.MAXIMUM_NUMBER_OF_FRAGMENT_FEATURES_DEFAULT, 0, Integer.MAX_VALUE,
                "maximumNumberOfFragmentFeaturesSetting");
        //aromaticity model for the descriptor calculation
        this.electronDonationModelSetting = this.createEnumSetting(this.aromaticitySettings, "Electron donation model setting",
                ClusteringService.ELECTRON_DONATION_MODEL_OPTION_DEFAULT, IMoleculeFragmenter.ElectronDonationModelOption.class,
                "electronDonationModelSetting");
        this.cycleFinderSetting = this.createEnumSetting(this.aromaticitySettings, "Cycle finder algorithm setting",
                ClusteringService.CYCLE_FINDER_OPTION_DEFAULT, IMoleculeFragmenter.CycleFinderOption.class,
                "cycleFinderSetting");
        //transformation and filter
        this.countTransformationSetting = this.createEnumSetting(this.transformationSettings, "Count transformation setting",
                ClusteringService.COUNT_TRANSFORMATION_OPTION_DEFAULT, CountTransformationOption.class,
                "countTransformationSetting");
        this.continuousTransformationSetting = this.createEnumSetting(this.transformationSettings, "Continuous transformation setting",
                ClusteringService.CONTINUOUS_TRANSFORMATION_OPTION_DEFAULT, ContinuousTransformationOption.class,
                "continuousTransformationSetting");
        this.maximumNaNFractionSetting = this.createDoubleSetting(this.filterSettings, "Maximum NaN fraction setting",
                ClusteringService.MAXIMUM_NAN_FRACTION_DEFAULT, 0.0, 1.0, "maximumNaNFractionSetting");
        this.correlationFilterSetting = this.createBooleanSetting(this.filterSettings, "Correlation filter setting",
                ClusteringService.CORRELATION_FILTER_DEFAULT, "correlationFilterSetting");
        this.correlationThresholdSetting = this.createDoubleSetting(this.correlationSettings, "Correlation threshold setting",
                ClusteringService.CORRELATION_THRESHOLD_DEFAULT, 0.0, 1.0, "correlationThresholdSetting");
        this.correlationMethodSetting = this.createEnumSetting(this.correlationSettings, "Correlation method setting",
                ClusteringService.CORRELATION_METHOD_OPTION_DEFAULT, CorrelationMethodOption.class,
                "correlationMethodSetting");
        //ART-2A
        this.algorithmSetting = this.createEnumSetting(this.art2aSettings, "Clustering algorithm setting",
                ClusteringService.ALGORITHM_OPTION_DEFAULT, AlgorithmOption.class, "algorithmSetting");
        this.vigilanceSetting = this.createDoubleSetting(this.art2aSettings, "Vigilance setting", ClusteringService.VIGILANCE_DEFAULT,
                ClusteringService.OPEN_INTERVAL_MINIMUM, ClusteringService.OPEN_INTERVAL_MAXIMUM, "vigilanceSetting");
        this.maximumNumberOfClustersSetting = this.createIntegerSetting(this.art2aSettings, "Maximum number of clusters setting",
                ClusteringService.MAXIMUM_NUMBER_OF_CLUSTERS_DEFAULT, 2, Integer.MAX_VALUE, "maximumNumberOfClustersSetting");
        this.maximumNumberOfEpochsSetting = this.createIntegerSetting(this.art2aSettings, "Maximum number of epochs setting",
                Art2aClusterer.MAXIMUM_NUMBER_OF_EPOCHS_DEFAULT, 1, Integer.MAX_VALUE, "maximumNumberOfEpochsSetting");
        this.convergenceThresholdSetting = this.createDoubleSetting(this.art2aSettings, "Convergence threshold setting",
                ClusteringService.ALGORITHM_OPTION_DEFAULT.getConvergenceThresholdDefault(),
                ClusteringService.OPEN_INTERVAL_MINIMUM, Double.MAX_VALUE, "convergenceThresholdSetting");
        this.learningParameterSetting = this.createDoubleSetting(this.art2aSettings, "Learning parameter setting",
                ClusteringService.LEARNING_PARAMETER_DEFAULT, ClusteringService.OPEN_INTERVAL_MINIMUM,
                ClusteringService.OPEN_INTERVAL_MAXIMUM, "learningParameterSetting");
        this.offsetForContrastEnhancementSetting = this.createDoubleSetting(this.art2aSettings, "Offset for contrast enhancement setting",
                ClusteringService.OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT, ClusteringService.OPEN_INTERVAL_MINIMUM,
                Double.MAX_VALUE, "offsetForContrastEnhancementSetting");
        this.randomSeedSetting = this.createIntegerSetting(this.art2aSettings, "Random seed setting",
                (int) Art2aClusterer.RANDOM_SEED_DEFAULT, 1, Integer.MAX_VALUE, "randomSeedSetting");
        this.numberOfSimilarMoleculesSetting = this.createIntegerSetting(this.similaritySearchSettings,
                "Number of similar molecules setting", ClusteringService.NUMBER_OF_SIMILAR_MOLECULES_DEFAULT, 1,
                Integer.MAX_VALUE, "numberOfSimilarMoleculesSetting");
        //all general settings in persistence order; the algorithm comes before the convergence threshold, so that a
        // reload sets the persisted threshold after the algorithm listener has set the default of the kernel
        this.settings.addAll(this.fragmentFingerprintSwitchSettings);
        this.settings.addAll(this.fragmentFingerprintSettings);
        this.settings.addAll(this.aromaticitySettings);
        this.settings.addAll(this.transformationSettings);
        this.settings.addAll(this.filterSettings);
        this.settings.addAll(this.correlationSettings);
        this.settings.addAll(this.art2aSettings);
        this.settings.addAll(this.similaritySearchSettings);
        //persisted with the general settings, but chosen in the result tabs
        this.displayedDescriptorColumnsSetting = new SimpleStringProperty(this, "Displayed descriptor columns setting",
                ClusteringService.joinDisplayedDescriptorColumnNames(ClusteringService.DISPLAYED_DESCRIPTOR_COLUMNS_DEFAULT));
        this.registerSetting(this.settings, this.displayedDescriptorColumnsSetting, "displayedDescriptorColumnsSetting");
        //the meaning of the convergence threshold depends on the kernel, so its default is set when the kernel changes
        this.algorithmSetting.addListener((anObservable, anOldValue, aNewValue) -> {
            if (aNewValue != null && aNewValue != anOldValue) {
                this.convergenceThresholdSetting.set(((AlgorithmOption) aNewValue).getConvergenceThresholdDefault());
            }
        });
        //descriptor tab: preset, fingerprint section (with advanced settings) and descriptor list
        this.isPresetBeingApplied = false;
        this.userDescriptorPresetMap = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);
        this.descriptorPresetKeyList = FXCollections.observableArrayList();
        this.updateDescriptorPresetKeyList();
        this.fingerprintSettings = new ArrayList<>(4);
        this.descriptorSwitchSettings = new ArrayList<>(1);
        this.fingerprintSwitchSettings = new ArrayList<>(1);
        this.descriptorsSetting = this.createBooleanSetting(this.descriptorSwitchSettings, "Descriptors setting",
                ClusteringService.DESCRIPTORS_DEFAULT, "descriptorsSetting");
        this.fingerprintsSetting = this.createBooleanSetting(this.fingerprintSwitchSettings, "Fingerprints setting",
                ClusteringService.FINGERPRINTS_DEFAULT, "fingerprintsSetting");
        this.advancedFingerprintSettings = new ArrayList<>(2);
        this.descriptorSelectionSettings = new ArrayList<>(Descriptor.getAllDescriptors().length);
        //the preset is a string setting (key of a built-in or user-defined preset) and has its own control in the view
        this.descriptorPresetSetting = new SimpleStringProperty(this, "Descriptor preset setting",
                ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.name());
        this.registerSetting(this.descriptorSettings, this.descriptorPresetSetting, "descriptorPresetSetting");
        this.lastAppliedPresetDisplayName = ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.getDisplayName();
        this.circularFingerprintSetting = this.createBooleanSetting(this.fingerprintSettings, "Circular fingerprint setting",
                ClusteringService.CIRCULAR_FINGERPRINT_DEFAULT, "circularFingerprintSetting");
        this.circularFingerprintTypeSetting = this.createEnumSetting(this.fingerprintSettings, "Circular fingerprint type setting",
                ClusteringService.CIRCULAR_FINGERPRINT_TYPE_OPTION_DEFAULT, CircularFingerprintTypeOption.class,
                "circularFingerprintTypeSetting");
        this.maccsFingerprintSetting = this.createBooleanSetting(this.fingerprintSettings, "MACCS fingerprint setting",
                ClusteringService.MACCS_FINGERPRINT_DEFAULT, "maccsFingerprintSetting");
        this.pubChemFingerprintSetting = this.createBooleanSetting(this.fingerprintSettings, "PubChem fingerprint setting",
                ClusteringService.PUBCHEM_FINGERPRINT_DEFAULT, "pubChemFingerprintSetting");
        this.circularFingerprintRadiusSetting = this.createEnumSetting(this.advancedFingerprintSettings, "Circular fingerprint radius setting",
                ClusteringService.CIRCULAR_FINGERPRINT_RADIUS_OPTION_DEFAULT, CircularFingerprintRadiusOption.class,
                "circularFingerprintRadiusSetting");
        this.circularFingerprintSizeSetting = this.createIntegerSetting(this.advancedFingerprintSettings, "Circular fingerprint size setting",
                ClusteringService.CIRCULAR_FINGERPRINT_SIZE_DEFAULT, 1, Integer.MAX_VALUE, "circularFingerprintSizeSetting");
        //the list only contains real descriptors, the fingerprints are chosen via the settings above; the calculator gives
        // the CDK names of their components for the tooltips
        DescriptorCalculator tmpCalculator = FeatureMatrixBuilder.getDescriptorCalculator(this.circularFingerprintSizeSetting.get());
        this.descriptorToSettingMap = new EnumMap<>(Descriptor.class);
        for (Descriptor tmpDescriptor : Descriptor.getAllDescriptors()) {
            if (tmpDescriptor.isFingerprint()) {
                continue;
            }
            SimpleBooleanProperty tmpSetting = new SimpleBooleanProperty(this,
                    ClusteringService.getDescriptorSettingName(tmpDescriptor),
                    ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.includes(tmpDescriptor));
            this.descriptorSelectionSettings.add(tmpSetting);
            this.descriptorToSettingMap.put(tmpDescriptor, tmpSetting);
            this.settingNameToDisplayNameMap.put(tmpSetting.getName(), ClusteringService.getDescriptorDisplayName(tmpDescriptor));
            String[] tmpComponentNames = tmpCalculator.getColumnNames(new Descriptor[] {tmpDescriptor});
            String tmpComponentNamesText = String.join(", ", Arrays.asList(tmpComponentNames).subList(0,
                    Math.min(tmpComponentNames.length, ClusteringService.MAXIMUM_NUMBER_OF_COMPONENT_NAMES_IN_TOOLTIP)))
                    + (tmpComponentNames.length > ClusteringService.MAXIMUM_NUMBER_OF_COMPONENT_NAMES_IN_TOOLTIP ? ", ..." : "");
            this.settingNameToTooltipTextMap.put(tmpSetting.getName(), String.format(
                    Message.get("ClusteringService.descriptorSetting.tooltip"),
                    tmpComponentNames.length, tmpComponentNamesText,
                    //the Wiener polarity number is used as count, see FeatureMatrixBuilder
                    tmpDescriptor == Descriptor.WIENER_NUMBER ? tmpDescriptor.getComponentType() + " (WPOL: COUNT)"
                            : tmpDescriptor.getComponentType().name(), tmpDescriptor.isFast(),
                    tmpDescriptor.isSafe(), tmpDescriptor.requires3DCoordinates()));
            //a manual change of a descriptor that deviates from the selected preset switches the selection to custom
            tmpSetting.addListener((anObservable, anOldValue, aNewValue) -> {
                DescriptorPreset tmpPreset = this.getDescriptorPreset(this.descriptorPresetSetting.get());
                if (!this.isPresetBeingApplied && tmpPreset != null && aNewValue != tmpPreset.includes(tmpDescriptor)) {
                    this.descriptorPresetSetting.set(ClusteringService.CUSTOM_PRESET_KEY);
                }
            });
        }
        //all settings of the descriptor tab in persistence order: the preset first, so that a reload applies it first
        this.descriptorSettings.addAll(this.descriptorSwitchSettings);
        this.descriptorSettings.addAll(this.fingerprintSwitchSettings);
        this.descriptorSettings.addAll(this.fingerprintSettings);
        this.descriptorSettings.addAll(this.advancedFingerprintSettings);
        this.descriptorSettings.addAll(this.descriptorSelectionSettings);
        //choosing a preset sets the descriptor check boxes accordingly
        this.descriptorPresetSetting.addListener((anObservable, anOldValue, aNewValue) -> {
            DescriptorPreset tmpPreset = this.getDescriptorPreset(aNewValue);
            if (tmpPreset != null) {
                this.applyDescriptorPreset(tmpPreset);
            } else if (aNewValue != null && !ClusteringService.CUSTOM_PRESET_KEY.equals(aNewValue)) {
                //unknown key, e.g. a deleted user preset in the persisted settings
                this.descriptorPresetSetting.set(ClusteringService.CUSTOM_PRESET_KEY);
            }
        });
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the general clustering settings in display order.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> settingsProperties() {
        return Collections.unmodifiableList(this.settings);
    }
    //
    /**
     * Returns the descriptor selection settings in display order.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> descriptorSettingsProperties() {
        return Collections.unmodifiableList(this.descriptorSettings);
    }
    //
    /**
     * Returns the map of setting names to display names (general and descriptor settings).
     *
     * @return unmodifiable map
     */
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return Collections.unmodifiableMap(this.settingNameToDisplayNameMap);
    }
    //
    /**
     * Returns the map of setting names to tooltip texts (general and descriptor settings).
     *
     * @return unmodifiable map
     */
    public Map<String, String> getSettingNameToTooltipTextMap() {
        return Collections.unmodifiableMap(this.settingNameToTooltipTextMap);
    }
    //
    /**
     * Returns the descriptors used for a clustering: the selected descriptors in enum order (only if the descriptors are
     * switched on), followed by the selected fingerprints (circular, MACCS, PubChem; only if the fingerprints are
     * switched on).
     *
     * @return descriptors and fingerprints of a clustering, may be empty
     */
    public Descriptor[] getSelectedDescriptors() {
        List<Descriptor> tmpSelected = new ArrayList<>();
        if (this.descriptorsSetting.get()) {
            tmpSelected.addAll(Arrays.asList(this.getSelectedNonFingerprintDescriptors()));
        }
        if (!this.fingerprintsSetting.get()) {
            return tmpSelected.toArray(new Descriptor[0]);
        }
        if (this.circularFingerprintSetting.get()) {
            tmpSelected.add(ClusteringService.getCircularFingerprintDescriptor(
                    (CircularFingerprintTypeOption) this.circularFingerprintTypeSetting.get(),
                    (CircularFingerprintRadiusOption) this.circularFingerprintRadiusSetting.get()));
        }
        if (this.maccsFingerprintSetting.get()) {
            tmpSelected.add(Descriptor.MACCS_FINGERPRINTER);
        }
        if (this.pubChemFingerprintSetting.get()) {
            tmpSelected.add(Descriptor.PUBCHEM_FINGERPRINTER);
        }
        return tmpSelected.toArray(new Descriptor[0]);
    }
    //
    /**
     * Returns the currently selected descriptors without the fingerprints, e.g. for a table of descriptor values in which
     * fingerprint bits would only add hundreds of bit columns.
     *
     * @return selected descriptors that are no fingerprints, in enum order, may be empty
     */
    public Descriptor[] getSelectedNonFingerprintDescriptors() {
        //independent of the descriptor switch, which only concerns the clustering
        List<Descriptor> tmpSelected = new ArrayList<>();
        for (Map.Entry<Descriptor, SimpleBooleanProperty> tmpEntry : this.descriptorToSettingMap.entrySet()) {
            if (tmpEntry.getValue().get()) {
                tmpSelected.add(tmpEntry.getKey());
            }
        }
        return tmpSelected.toArray(new Descriptor[0]);
    }
    //
    /**
     * Creates the aromaticity model of the current aromaticity settings, which is applied to the molecules before the
     * descriptor calculation.
     *
     * @return aromaticity model
     */
    public Aromaticity createAromaticityModel() {
        return ChemUtil.createAromaticityModel((IMoleculeFragmenter.ElectronDonationModelOption) this.electronDonationModelSetting.get(),
                (IMoleculeFragmenter.CycleFinderOption) this.cycleFinderSetting.get());
    }
    //
    /**
     * Selects or deselects the given descriptor; a deviation from the current preset switches it to custom.
     *
     * @param aDescriptor descriptor, must not be a fingerprint
     * @param anIsSelected true to select
     * @throws NullPointerException if the descriptor is null
     * @throws IllegalArgumentException if the descriptor is a fingerprint
     */
    public void setDescriptorSelected(Descriptor aDescriptor, boolean anIsSelected) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aDescriptor, "Given descriptor is null.");
        if (aDescriptor.isFingerprint()) {
            throw new IllegalArgumentException("Fingerprints are selected via the fingerprint settings, not via the descriptor list.");
        }
        this.descriptorToSettingMap.get(aDescriptor).set(anIsSelected);
    }
    //
    /**
     * Sets the aromaticity settings (electron donation model and cycle finder) to MORTAR's standard model, e.g. for the
     * default button of the descriptor calculation view.
     */
    public void restoreDefaultAromaticitySettings() {
        this.electronDonationModelSetting.set(ClusteringService.ELECTRON_DONATION_MODEL_OPTION_DEFAULT);
        this.cycleFinderSetting.set(ClusteringService.CYCLE_FINDER_OPTION_DEFAULT);
    }
    //
    /**
     * Restores the default values of the fragment fingerprint settings (switch, type, minimum frequency, maximum number of
     * features), which are shown in the descriptor tab of the settings view.
     */
    public void restoreDefaultFragmentFingerprintSettings() {
        this.fragmentFingerprintSetting.set(ClusteringService.FRAGMENT_FINGERPRINT_DEFAULT);
        this.fragmentFingerprintTypeSetting.set(ClusteringService.FRAGMENT_FINGERPRINT_TYPE_OPTION_DEFAULT);
        this.minimumFragmentMoleculeFrequencySetting.set(ClusteringService.MINIMUM_FRAGMENT_MOLECULE_FREQUENCY_DEFAULT);
        this.maximumNumberOfFragmentFeaturesSetting.set(ClusteringService.MAXIMUM_NUMBER_OF_FRAGMENT_FEATURES_DEFAULT);
    }
    //
    /**
     * Restores the default values of the general clustering settings of the clustering tab of the settings view; the
     * fragment fingerprint settings (descriptor tab) are restored by {@link #restoreDefaultFragmentFingerprintSettings()}.
     */
    public void restoreDefaultSettings() {
        this.restoreDefaultAromaticitySettings();
        this.countTransformationSetting.set(ClusteringService.COUNT_TRANSFORMATION_OPTION_DEFAULT);
        this.continuousTransformationSetting.set(ClusteringService.CONTINUOUS_TRANSFORMATION_OPTION_DEFAULT);
        this.maximumNaNFractionSetting.set(ClusteringService.MAXIMUM_NAN_FRACTION_DEFAULT);
        this.correlationFilterSetting.set(ClusteringService.CORRELATION_FILTER_DEFAULT);
        this.numberOfSimilarMoleculesSetting.set(ClusteringService.NUMBER_OF_SIMILAR_MOLECULES_DEFAULT);
        this.correlationThresholdSetting.set(ClusteringService.CORRELATION_THRESHOLD_DEFAULT);
        this.correlationMethodSetting.set(ClusteringService.CORRELATION_METHOD_OPTION_DEFAULT);
        this.algorithmSetting.set(ClusteringService.ALGORITHM_OPTION_DEFAULT);
        this.vigilanceSetting.set(ClusteringService.VIGILANCE_DEFAULT);
        this.maximumNumberOfClustersSetting.set(ClusteringService.MAXIMUM_NUMBER_OF_CLUSTERS_DEFAULT);
        this.maximumNumberOfEpochsSetting.set(Art2aClusterer.MAXIMUM_NUMBER_OF_EPOCHS_DEFAULT);
        this.convergenceThresholdSetting.set(ClusteringService.ALGORITHM_OPTION_DEFAULT.getConvergenceThresholdDefault());
        this.learningParameterSetting.set(ClusteringService.LEARNING_PARAMETER_DEFAULT);
        this.offsetForContrastEnhancementSetting.set(ClusteringService.OFFSET_FOR_CONTRAST_ENHANCEMENT_DEFAULT);
        this.randomSeedSetting.set((int) Art2aClusterer.RANDOM_SEED_DEFAULT);
    }
    //
    /**
     * Restores the default descriptor preset and fingerprint settings, which also resets the descriptor check boxes.
     * User-defined presets are kept.
     */
    public void restoreDefaultDescriptorSelection() {
        this.descriptorsSetting.set(ClusteringService.DESCRIPTORS_DEFAULT);
        this.fingerprintsSetting.set(ClusteringService.FINGERPRINTS_DEFAULT);
        this.selectDescriptorPreset(ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.name());
        this.circularFingerprintSetting.set(ClusteringService.CIRCULAR_FINGERPRINT_DEFAULT);
        this.circularFingerprintTypeSetting.set(ClusteringService.CIRCULAR_FINGERPRINT_TYPE_OPTION_DEFAULT);
        this.circularFingerprintRadiusSetting.set(ClusteringService.CIRCULAR_FINGERPRINT_RADIUS_OPTION_DEFAULT);
        this.circularFingerprintSizeSetting.set(ClusteringService.CIRCULAR_FINGERPRINT_SIZE_DEFAULT);
        this.maccsFingerprintSetting.set(ClusteringService.MACCS_FINGERPRINT_DEFAULT);
        this.pubChemFingerprintSetting.set(ClusteringService.PUBCHEM_FINGERPRINT_DEFAULT);
    }
    //
    /**
     * Returns all descriptor presets for the selection in the GUI: the user-defined presets first (sorted by name),
     * then the built-in presets.
     *
     * @return list of presets
     */
    public List<DescriptorPreset> getDescriptorPresets() {
        List<DescriptorPreset> tmpPresets = new ArrayList<>(this.userDescriptorPresetMap.size() + DescriptorPresetOption.values().length);
        for (Map.Entry<String, Set<Descriptor>> tmpEntry : this.userDescriptorPresetMap.entrySet()) {
            tmpPresets.add(ClusteringService.createUserPreset(tmpEntry.getKey(), tmpEntry.getValue()));
        }
        for (DescriptorPresetOption tmpOption : DescriptorPresetOption.values()) {
            tmpPresets.add(tmpOption.toDescriptorPreset());
        }
        return tmpPresets;
    }
    //
    /**
     * Returns the preset with the given key.
     *
     * @param aKey preset key
     * @return preset or null if there is no preset with this key (e.g. for {@link #CUSTOM_PRESET_KEY})
     */
    public DescriptorPreset getDescriptorPreset(String aKey) {
        if (aKey == null) {
            return null;
        }
        if (aKey.startsWith(DescriptorPreset.USER_PRESET_KEY_PREFIX)) {
            String tmpName = aKey.substring(DescriptorPreset.USER_PRESET_KEY_PREFIX.length());
            Set<Descriptor> tmpDescriptors = this.userDescriptorPresetMap.get(tmpName);
            return tmpDescriptors == null ? null : ClusteringService.createUserPreset(tmpName, tmpDescriptors);
        }
        for (DescriptorPresetOption tmpOption : DescriptorPresetOption.values()) {
            if (tmpOption.name().equals(aKey)) {
                return tmpOption.toDescriptorPreset();
            }
        }
        return null;
    }
    //
    /**
     * Returns the setting holding the key of the selected preset ({@link #CUSTOM_PRESET_KEY} if the check boxes were
     * changed after the preset was applied).
     *
     * @return preset key setting
     */
    public SimpleStringProperty descriptorPresetKeyProperty() {
        return this.descriptorPresetSetting;
    }
    //
    /**
     * Selects the preset with the given key and sets the descriptor check boxes accordingly, also if the preset was
     * already selected but the check boxes were changed in the meantime.
     *
     * @param aKey preset key
     * @throws IllegalArgumentException if there is no preset with this key
     */
    public void selectDescriptorPreset(String aKey) throws IllegalArgumentException {
        DescriptorPreset tmpPreset = this.getDescriptorPreset(aKey);
        if (tmpPreset == null) {
            throw new IllegalArgumentException("Unknown descriptor preset: " + aKey);
        }
        this.descriptorPresetSetting.set(aKey);
        this.applyDescriptorPreset(tmpPreset);
    }
    //
    /**
     * Returns a proposal for the name under which the current (modified) selection could be saved, e.g. "Fast (custom)".
     *
     * @return proposed preset name
     */
    public String getProposedUserPresetName() {
        return this.lastAppliedPresetDisplayName + " " + Message.get("ClusteringService.userPreset.proposedNameSuffix");
    }
    //
    /**
     * Saves the currently checked descriptors as user-defined preset with the given name, selects it and writes all
     * user-defined presets to the settings directory. An existing user-defined preset with the same name is replaced.
     *
     * @param aName preset name; must be a valid preference term and must not be the name of a built-in preset
     * @return key of the saved preset
     * @throws NullPointerException if the name is null
     * @throws IllegalArgumentException if the name is invalid
     * @throws IOException if the presets cannot be written
     */
    public String saveCurrentDescriptorSelectionAsPreset(String aName) throws NullPointerException, IllegalArgumentException, IOException {
        String tmpName = Objects.requireNonNull(aName, "Given preset name is null.").trim();
        if (!ClusteringService.isValidUserPresetName(tmpName)) {
            throw new IllegalArgumentException(Message.get("ClusteringService.userPreset.invalidName"));
        }
        Set<Descriptor> tmpDescriptors = EnumSet.noneOf(Descriptor.class);
        for (Map.Entry<Descriptor, SimpleBooleanProperty> tmpEntry : this.descriptorToSettingMap.entrySet()) {
            if (tmpEntry.getValue().get()) {
                tmpDescriptors.add(tmpEntry.getKey());
            }
        }
        this.userDescriptorPresetMap.put(tmpName, tmpDescriptors);
        this.updateDescriptorPresetKeyList();
        String tmpKey = DescriptorPreset.getUserPresetKey(tmpName);
        this.descriptorPresetSetting.set(tmpKey);
        this.lastAppliedPresetDisplayName = tmpName;
        this.persistUserDescriptorPresets();
        return tmpKey;
    }
    //
    /**
     * Deletes the user-defined preset with the given key and writes the remaining user-defined presets. If it was
     * selected, the selection becomes custom (the check boxes are kept).
     *
     * @param aKey key of a user-defined preset
     * @throws IllegalArgumentException if the key does not belong to a user-defined preset
     * @throws IOException if the presets cannot be written
     */
    public void deleteUserDescriptorPreset(String aKey) throws IllegalArgumentException, IOException {
        DescriptorPreset tmpPreset = this.getDescriptorPreset(aKey);
        if (tmpPreset == null || !tmpPreset.isUserDefined()) {
            throw new IllegalArgumentException("No user-defined descriptor preset: " + aKey);
        }
        this.userDescriptorPresetMap.remove(tmpPreset.displayName());
        this.updateDescriptorPresetKeyList();
        if (aKey.equals(this.descriptorPresetSetting.get())) {
            this.descriptorPresetSetting.set(ClusteringService.CUSTOM_PRESET_KEY);
        }
        this.persistUserDescriptorPresets();
    }
    //
    /**
     * Sets the descriptor check boxes according to the given preset without switching the selection to custom.
     *
     * @param aPreset preset to apply
     */
    public void applyDescriptorPreset(DescriptorPreset aPreset) {
        this.isPresetBeingApplied = true;
        try {
            for (Map.Entry<Descriptor, SimpleBooleanProperty> tmpEntry : this.descriptorToSettingMap.entrySet()) {
                tmpEntry.getValue().set(aPreset.includes(tmpEntry.getKey()));
            }
            this.lastAppliedPresetDisplayName = aPreset.displayName();
        } finally {
            this.isPresetBeingApplied = false;
        }
    }
    //
    /**
     * Sets the descriptor check boxes according to the given built-in preset, see {@link #applyDescriptorPreset(DescriptorPreset)}.
     *
     * @param aPresetOption built-in preset to apply
     */
    public void applyDescriptorPreset(DescriptorPresetOption aPresetOption) {
        this.applyDescriptorPreset(aPresetOption.toDescriptorPreset());
    }
    //
    /**
     * Returns the observable list of all preset keys (own presets first, then built-in presets, then the custom key).
     * It is updated when presets are saved, deleted or reloaded.
     *
     * @return unmodifiable observable list of preset keys
     */
    public ObservableList<String> getDescriptorPresetKeyList() {
        return FXCollections.unmodifiableObservableList(this.descriptorPresetKeyList);
    }
    //
    /**
     * Returns the settings of the fragment fingerprint section of the clustering tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> fragmentFingerprintSettingsProperties() {
        return Collections.unmodifiableList(this.fragmentFingerprintSettings);
    }
    //
    /**
     * Returns the settings of the aromaticity section of the clustering tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> aromaticitySettingsProperties() {
        return Collections.unmodifiableList(this.aromaticitySettings);
    }
    //
    /**
     * Returns the settings of the transformation section of the clustering tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> transformationSettingsProperties() {
        return Collections.unmodifiableList(this.transformationSettings);
    }
    //
    /**
     * Returns the settings of the filter section of the clustering tab, including the correlation filter switch.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> filterSettingsProperties() {
        return Collections.unmodifiableList(this.filterSettings);
    }
    //
    /**
     * Returns the settings of the correlation filter (threshold and coefficient).
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> correlationSettingsProperties() {
        return Collections.unmodifiableList(this.correlationSettings);
    }
    //
    /**
     * Returns the settings of the ART-2A section of the clustering tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> art2aSettingsProperties() {
        return Collections.unmodifiableList(this.art2aSettings);
    }
    //
    /**
     * Returns the settings of the similarity search section of the clustering tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> similaritySearchSettingsProperties() {
        return Collections.unmodifiableList(this.similaritySearchSettings);
    }
    //
    /**
     * Sets the vigilance setting, e.g. to a value chosen from a vigilance scan.
     *
     * @param aVigilance vigilance in (0,1)
     * @throws IllegalArgumentException if the value is outside the allowed range
     */
    public void setVigilance(double aVigilance) throws IllegalArgumentException {
        this.vigilanceSetting.set(aVigilance);
    }
    //
    /**
     * Returns the vigilance setting, e.g. to bind the vigilance field of the clustering bars to it.
     *
     * @return vigilance setting
     */
    public SimpleDoubleProperty vigilanceProperty() {
        return this.vigilanceSetting;
    }
    //
    /**
     * Returns the names of the descriptor columns shown with their cluster medians in the clustering result tabs.
     *
     * @return list of descriptor column names in display order
     */
    public List<String> getDisplayedDescriptorColumnNames() {
        List<String> tmpNames = new ArrayList<>();
        for (String tmpName : this.displayedDescriptorColumnsSetting.get().split(ClusteringService.DISPLAYED_DESCRIPTOR_COLUMNS_SEPARATOR)) {
            if (!tmpName.isBlank()) {
                tmpNames.add(tmpName);
            }
        }
        return tmpNames;
    }
    //
    /**
     * Sets the names of the descriptor columns shown with their cluster medians in the clustering result tabs. Names that
     * cannot be persisted (characters not allowed in a preference or the separator) are left out.
     *
     * @param aNameList descriptor column names in display order
     * @throws NullPointerException if the list is null
     */
    public void setDisplayedDescriptorColumnNames(List<String> aNameList) throws NullPointerException {
        Objects.requireNonNull(aNameList, "Given name list is null.");
        this.displayedDescriptorColumnsSetting.set(ClusteringService.joinDisplayedDescriptorColumnNames(aNameList));
    }
    //
    /**
     * Returns the setting of the displayed descriptor columns, e.g. to update all result tabs when it changes; use
     * {@link #getDisplayedDescriptorColumnNames()} to read it.
     *
     * @return displayed descriptor columns setting
     */
    public SimpleStringProperty displayedDescriptorColumnsProperty() {
        return this.displayedDescriptorColumnsSetting;
    }
    //
    /**
     * Returns the number of similar molecules the similarity search returns.
     *
     * @return number of similar molecules
     */
    public int getNumberOfSimilarMolecules() {
        return this.numberOfSimilarMoleculesSetting.get();
    }
    //
    /**
     * Returns the setting whether the selected descriptors are used for the clustering, in its own list for the settings
     * view.
     *
     * @return unmodifiable list with the descriptor switch
     */
    public List<Property<?>> descriptorSwitchSettingsProperties() {
        return Collections.unmodifiableList(this.descriptorSwitchSettings);
    }
    //
    /**
     * Returns the setting whether the selected fingerprints are used for the clustering, in its own list for the settings
     * view.
     *
     * @return unmodifiable list with the fingerprint switch
     */
    public List<Property<?>> fingerprintSwitchSettingsProperties() {
        return Collections.unmodifiableList(this.fingerprintSwitchSettings);
    }
    //
    /**
     * Returns the setting whether the selected descriptors are used for the clustering, e.g. to show the preset and the
     * descriptor list only while it is switched on.
     *
     * @return descriptor switch setting
     */
    public SimpleBooleanProperty descriptorsProperty() {
        return this.descriptorsSetting;
    }
    //
    /**
     * Returns the setting whether the selected fingerprints are used for the clustering, e.g. to show the fingerprint
     * settings only while it is switched on.
     *
     * @return fingerprint switch setting
     */
    public SimpleBooleanProperty fingerprintsProperty() {
        return this.fingerprintsSetting;
    }
    //
    /**
     * Returns the setting whether the fragment fingerprint is used, in its own list, e.g. for a grid of the settings view
     * above the other fragment fingerprint settings.
     *
     * @return unmodifiable list with the fragment fingerprint switch
     */
    public List<Property<?>> fragmentFingerprintSwitchSettingsProperties() {
        return Collections.unmodifiableList(this.fragmentFingerprintSwitchSettings);
    }
    //
    /**
     * Returns the setting whether the fragment fingerprint is part of the features, e.g. to show the other fragment
     * fingerprint settings only while it is switched on.
     *
     * @return fragment fingerprint setting
     */
    public SimpleBooleanProperty fragmentFingerprintProperty() {
        return this.fragmentFingerprintSetting;
    }
    //
    /**
     * Returns the setting whether the correlation filter is applied.
     *
     * @return correlation filter setting
     */
    public SimpleBooleanProperty correlationFilterProperty() {
        return this.correlationFilterSetting;
    }
    //
    /**
     * Returns the settings of the fingerprint section of the descriptor tab.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> fingerprintSettingsProperties() {
        return Collections.unmodifiableList(this.fingerprintSettings);
    }
    //
    /**
     * Returns the advanced settings of the fingerprint section (radius and length of the circular fingerprint).
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> advancedFingerprintSettingsProperties() {
        return Collections.unmodifiableList(this.advancedFingerprintSettings);
    }
    //
    /**
     * Returns the check box settings of the descriptor list.
     *
     * @return unmodifiable list of settings
     */
    public List<Property<?>> descriptorSelectionSettingsProperties() {
        return Collections.unmodifiableList(this.descriptorSelectionSettings);
    }
    //
    /**
     * Writes the user-defined presets to {@link #USER_DESCRIPTOR_PRESETS_FILE_NAME} in the clustering settings directory
     * (preset name as key, comma-separated descriptor constant names as value).
     *
     * @throws IOException if the file cannot be written
     */
    public void persistUserDescriptorPresets() throws IOException {
        File tmpDirectory = new File(this.settingsDirectoryPath);
        if (!tmpDirectory.exists() && !tmpDirectory.mkdirs()) {
            throw new IOException("Clustering settings directory could not be created.");
        }
        Properties tmpProperties = new Properties();
        for (Map.Entry<String, Set<Descriptor>> tmpEntry : this.userDescriptorPresetMap.entrySet()) {
            tmpProperties.setProperty(tmpEntry.getKey(),
                    String.join(",", tmpEntry.getValue().stream().map(Descriptor::name).toList()));
        }
        try (Writer tmpWriter = Files.newBufferedWriter(Paths.get(this.settingsDirectoryPath
                + ClusteringService.USER_DESCRIPTOR_PRESETS_FILE_NAME), StandardCharsets.UTF_8)) {
            tmpProperties.store(tmpWriter, "MORTAR user-defined descriptor presets");
        }
    }
    //
    /**
     * Creates an immutable snapshot of the current settings for a clustering run.
     *
     * @param aFragmentationName name of the fragmentation whose fragments form the fragment fingerprint
     * @param anIsParallelCalculation true: descriptor calculation and winner search are parallelized
     * @return clustering parameters
     * @throws NullPointerException if the fragmentation name is null
     */
    public ClusteringParameters createClusteringParameters(String aFragmentationName, boolean anIsParallelCalculation)
            throws NullPointerException {
        return new ClusteringParameters(
                aFragmentationName,
                this.fragmentFingerprintTypeSetting.get() == FragmentFingerprintTypeOption.COUNT,
                this.minimumFragmentMoleculeFrequencySetting.get(),
                this.maximumNumberOfFragmentFeaturesSetting.get(),
                this.getSelectedDescriptors(),
                (IMoleculeFragmenter.ElectronDonationModelOption) this.electronDonationModelSetting.get(),
                (IMoleculeFragmenter.CycleFinderOption) this.cycleFinderSetting.get(),
                ((CountTransformationOption) this.countTransformationSetting.get()).getTransformation(),
                ((ContinuousTransformationOption) this.continuousTransformationSetting.get()).getTransformation(),
                this.maximumNaNFractionSetting.get(),
                this.correlationFilterSetting.get(),
                this.correlationThresholdSetting.get(),
                ((CorrelationMethodOption) this.correlationMethodSetting.get()).getMethod(),
                ((AlgorithmOption) this.algorithmSetting.get()).getAlgorithm(),
                (float) this.vigilanceSetting.get(),
                this.maximumNumberOfClustersSetting.get(),
                this.maximumNumberOfEpochsSetting.get(),
                (float) this.convergenceThresholdSetting.get(),
                (float) this.learningParameterSetting.get(),
                (float) this.offsetForContrastEnhancementSetting.get(),
                this.randomSeedSetting.get(),
                this.circularFingerprintSizeSetting.get(),
                anIsParallelCalculation,
                this.fragmentFingerprintSetting.get());
    }
    //
    /**
     * Creates a clustering task for the given molecules with the current settings.
     *
     * @param aMoleculeList molecules to cluster
     * @param aFragmentationName name of the fragmentation whose fragments form the fragment fingerprint
     * @param anIsParallelCalculation true: descriptor calculation and winner search are parallelized
     * @return new clustering task, not started
     * @throws NullPointerException if a parameter is null
     * @throws IOException if the priority table cannot be loaded
     */
    public ClusteringTask createClusteringTask(List<MoleculeDataModel> aMoleculeList, String aFragmentationName,
                                               boolean anIsParallelCalculation) throws NullPointerException, IOException {
        return this.createClusteringTask(aMoleculeList, this.createClusteringParameters(aFragmentationName, anIsParallelCalculation));
    }
    //
    /**
     * Creates a clustering task for the given molecules with the given parameters; descriptors already calculated for
     * these molecules (in the descriptor cache) are not calculated again.
     *
     * @param aMoleculeList molecules to cluster
     * @param aParameters parameters of the run
     * @return new clustering task, not started
     * @throws NullPointerException if a parameter is null
     * @throws IOException if the priority table cannot be loaded
     */
    public ClusteringTask createClusteringTask(List<MoleculeDataModel> aMoleculeList, ClusteringParameters aParameters)
            throws NullPointerException, IOException {
        return new ClusteringTask(aMoleculeList, aParameters, this.getPriorityTable(), this.descriptorMatrixCache);
    }
    //
    /**
     * Creates the clustering parameters of the current settings, but with exactly the given descriptors and fingerprints
     * instead of the selected ones (e.g. the features calculated in a descriptor tab, which are then taken from the
     * descriptor cache).
     *
     * @param aFragmentationName name of the fragmentation whose fragments form the fragment fingerprint
     * @param aDescriptors descriptors and fingerprints to use instead of the selected ones, in this order
     * @param anIsParallelCalculation true: descriptor calculation and winner search are parallelized
     * @return clustering parameters
     * @throws NullPointerException if a parameter is null
     */
    public ClusteringParameters createClusteringParameters(String aFragmentationName, Descriptor[] aDescriptors,
                                                           boolean anIsParallelCalculation) throws NullPointerException {
        Objects.requireNonNull(aDescriptors, "Given descriptor array is null.");
        ClusteringParameters tmpParameters = this.createClusteringParameters(aFragmentationName, anIsParallelCalculation);
        return new ClusteringParameters(tmpParameters.fragmentationName(), tmpParameters.isCountFingerprint(),
                tmpParameters.minimumFragmentMoleculeFrequency(), tmpParameters.maximumNumberOfFragmentFeatures(),
                Arrays.stream(aDescriptors).distinct().toArray(Descriptor[]::new), tmpParameters.electronDonationModelOption(), tmpParameters.cycleFinderOption(),
                tmpParameters.countTransformation(), tmpParameters.continuousTransformation(), tmpParameters.maximumNaNFraction(),
                tmpParameters.isCorrelationFilterEnabled(), tmpParameters.correlationThreshold(), tmpParameters.correlationMethod(),
                tmpParameters.algorithm(), tmpParameters.vigilance(), tmpParameters.maximumNumberOfClusters(),
                tmpParameters.maximumNumberOfEpochs(), tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(),
                tmpParameters.offsetForContrastEnhancement(), tmpParameters.randomSeed(), tmpParameters.circularFingerprintSize(),
                tmpParameters.isParallelCalculation(), tmpParameters.isFragmentFingerprintUsed());
    }
    //
    /**
     * Creates the clustering parameters of the current settings for a clustering without fragmentation: the fragment
     * fingerprint is switched off, the fragmentation name is {@link #WITHOUT_FRAGMENTATION_NAME}, and exactly the given
     * descriptors and fingerprints (or the selected ones) are the features.
     *
     * @param aDescriptors descriptors and fingerprints to use instead of the selected ones; null for the selected ones
     * @param anIsParallelCalculation true: descriptor calculation and winner search are parallelized
     * @return clustering parameters without fragment fingerprint
     */
    public ClusteringParameters createClusteringParametersWithoutFragmentation(Descriptor[] aDescriptors, boolean anIsParallelCalculation) {
        ClusteringParameters tmpParameters = aDescriptors == null
                ? this.createClusteringParameters(ClusteringService.WITHOUT_FRAGMENTATION_NAME, anIsParallelCalculation)
                : this.createClusteringParameters(ClusteringService.WITHOUT_FRAGMENTATION_NAME, aDescriptors, anIsParallelCalculation);
        return new ClusteringParameters(tmpParameters.fragmentationName(), tmpParameters.isCountFingerprint(),
                tmpParameters.minimumFragmentMoleculeFrequency(), tmpParameters.maximumNumberOfFragmentFeatures(),
                tmpParameters.descriptors(), tmpParameters.electronDonationModelOption(), tmpParameters.cycleFinderOption(),
                tmpParameters.countTransformation(), tmpParameters.continuousTransformation(), tmpParameters.maximumNaNFraction(),
                tmpParameters.isCorrelationFilterEnabled(), tmpParameters.correlationThreshold(), tmpParameters.correlationMethod(),
                tmpParameters.algorithm(), tmpParameters.vigilance(), tmpParameters.maximumNumberOfClusters(),
                tmpParameters.maximumNumberOfEpochs(), tmpParameters.convergenceThreshold(), tmpParameters.learningParameter(),
                tmpParameters.offsetForContrastEnhancement(), tmpParameters.randomSeed(), tmpParameters.circularFingerprintSize(),
                tmpParameters.isParallelCalculation(), false);
    }
    //
    /**
     * Calculates the given descriptors of the given molecules with the aromaticity model of the current settings; cached
     * values are reused and the calculated ones are added to the descriptor cache, so that a following clustering with
     * these descriptors does not calculate them again.
     *
     * @param aMoleculeList molecules (rows of the matrix)
     * @param aDescriptorArray descriptors and fingerprints to calculate
     * @param anIsParallelCalculation true: molecules are processed in parallel
     * @return descriptor matrix with the columns in the order of the given descriptors
     * @throws NullPointerException if a parameter is null
     * @throws InterruptedException if the calculation is interrupted
     */
    public FeatureMatrix calculateDescriptorMatrix(List<MoleculeDataModel> aMoleculeList, Descriptor[] aDescriptorArray,
                                                   boolean anIsParallelCalculation) throws NullPointerException, InterruptedException {
        return this.descriptorMatrixCache.getDescriptorMatrix(aMoleculeList, aDescriptorArray,
                (IMoleculeFragmenter.ElectronDonationModelOption) this.electronDonationModelSetting.get(),
                (IMoleculeFragmenter.CycleFinderOption) this.cycleFinderSetting.get(),
                FeatureMatrixBuilder.getDescriptorCalculator(this.circularFingerprintSizeSetting.get()), anIsParallelCalculation);
    }
    //
    /**
     * Returns the cache of calculated descriptor values, e.g. to pass it to a pipeline run of a vigilance scan.
     *
     * @return descriptor cache
     */
    public DescriptorMatrixCache getDescriptorMatrixCache() {
        return this.descriptorMatrixCache;
    }
    //
    /**
     * Discards all cached descriptor values, e.g. when new molecules are imported.
     */
    public void clearDescriptorMatrixCache() {
        this.descriptorMatrixCache.clear();
    }
    //
    /**
     * Returns the priority table for the correlation filter; the default table is loaded on first use.
     *
     * @return priority table
     * @throws IOException if the default table cannot be loaded
     */
    public DescriptorPriorityTable getPriorityTable() throws IOException {
        if (this.priorityTable == null) {
            this.priorityTable = DescriptorPriorityTable.loadDefaultTable();
        }
        return this.priorityTable;
    }
    //
    /**
     * Persists the general settings and the descriptor selection in preference container files in a subfolder of the
     * settings directory. Failures are logged and shown to the user.
     */
    public void persistClusteringSettings() {
        String tmpDirectoryPath = this.settingsDirectoryPath;
        File tmpDirectory = new File(tmpDirectoryPath);
        if (!tmpDirectory.exists() && !tmpDirectory.mkdirs() || !tmpDirectory.canWrite()) {
            GuiUtil.guiMessageAlert(javafx.scene.control.Alert.AlertType.ERROR, Message.get("Error.ExceptionAlert.Title"),
                    Message.get("Error.ExceptionAlert.Header"),
                    Message.get("ClusteringService.Error.settingsPersistence"));
            return;
        }
        try {
            PreferenceUtil.translateJavaFxPropertiesToPreferences(this.settings, tmpDirectoryPath
                    + ClusteringService.CLUSTERING_SETTINGS_FILE_NAME + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION)
                    .writeRepresentation();
            PreferenceUtil.translateJavaFxPropertiesToPreferences(this.descriptorSettings, tmpDirectoryPath
                    + ClusteringService.DESCRIPTOR_SELECTION_FILE_NAME + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION)
                    .writeRepresentation();
        } catch (NullPointerException | IllegalArgumentException | IOException | SecurityException anException) {
            ClusteringService.LOGGER.log(Level.WARNING, String.format("Clustering settings persistence went wrong, exception: %s",
                    anException.toString()), anException);
            GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                    Message.get("Error.ExceptionAlert.Header"),
                    Message.get("ClusteringService.Error.settingsPersistence"),
                    anException);
        }
    }
    //
    /**
     * Reloads the user-defined descriptor presets, the persisted general settings and the descriptor selection. If a file is missing or cannot be read, the
     * respective settings remain in default and the problem is logged.
     */
    public void reloadClusteringSettings() {
        String tmpDirectoryPath = this.settingsDirectoryPath;
        //user presets first, so that a persisted selection of a user preset can be resolved
        this.reloadUserDescriptorPresets(new File(tmpDirectoryPath + ClusteringService.USER_DESCRIPTOR_PRESETS_FILE_NAME));
        ClusteringService.reloadSettings(this.settings, new File(tmpDirectoryPath
                + ClusteringService.CLUSTERING_SETTINGS_FILE_NAME + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION));
        ClusteringService.reloadSettings(this.descriptorSettings, new File(tmpDirectoryPath
                + ClusteringService.DESCRIPTOR_SELECTION_FILE_NAME + BasicDefinitions.PREFERENCE_CONTAINER_FILE_EXTENSION));
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Returns whether the given descriptor is selected by default, i.e. by the default preset.
     *
     * @param aDescriptor descriptor
     * @return true if selected by default
     */
    public static boolean isDescriptorSelectedByDefault(Descriptor aDescriptor) {
        return ClusteringService.DESCRIPTOR_PRESET_OPTION_DEFAULT.includes(aDescriptor);
    }
    //
    /**
     * Converts the given float into the double with the same decimal representation, e.g. 0.99f into 0.99 instead of
     * 0.9900000095367432 (the exact value of 0.99f). Used for the library defaults, which are floats, because the
     * settings are double properties and their values are shown in the GUI.
     *
     * @param aValue float value
     * @return double value with the same decimal representation
     */
    public static double toDecimalDouble(float aValue) {
        return Double.parseDouble(Float.toString(aValue));
    }
    //
    /**
     * Returns whether the given name can be used for a user-defined preset: it must not be blank, must be a valid
     * preference term (letters, digits, spaces and simple punctuation) and must not equal the name or display name of
     * a built-in preset (case-insensitive).
     *
     * @param aName preset name
     * @return true if the name is valid
     */
    public static boolean isValidUserPresetName(String aName) {
        if (aName == null || aName.isBlank() || !SingleTermPreference.isValidContent(aName)
                || ClusteringService.CUSTOM_PRESET_KEY.equalsIgnoreCase(aName.trim())) {
            return false;
        }
        for (DescriptorPresetOption tmpOption : DescriptorPresetOption.values()) {
            if (tmpOption.name().equalsIgnoreCase(aName.trim()) || tmpOption.getDisplayName().equalsIgnoreCase(aName.trim())) {
                return false;
            }
        }
        return true;
    }
    //
    /**
     * Returns the circular fingerprint descriptor of the given type and radius, e.g. CIRCULAR_FINGERPRINTER_ECFP_4 for
     * ECFP with radius 2.
     *
     * @param aType fingerprint type
     * @param aRadius fingerprint radius
     * @return circular fingerprint descriptor
     */
    public static Descriptor getCircularFingerprintDescriptor(CircularFingerprintTypeOption aType, CircularFingerprintRadiusOption aRadius) {
        return Descriptor.valueOf("CIRCULAR_FINGERPRINTER_" + aType.name() + "_" + aRadius.getDiameter());
    }
    //
    /**
     * Returns the display name of a descriptor in the descriptor list, with hints for slow and unsafe descriptors and
     * descriptors that need 3D coordinates, e.g. "Chi Path (slow)".
     *
     * @param aDescriptor descriptor
     * @return display name
     */
    public static String getDescriptorDisplayName(Descriptor aDescriptor) {
        List<String> tmpHints = new ArrayList<>(3);
        if (!aDescriptor.isFast()) {
            tmpHints.add(Message.get("ClusteringService.descriptorSetting.slowHint"));
        }
        if (!aDescriptor.isSafe()) {
            tmpHints.add(Message.get("ClusteringService.descriptorSetting.unsafeHint"));
        }
        if (aDescriptor.requires3DCoordinates()) {
            tmpHints.add(Message.get("ClusteringService.descriptorSetting.threeDimensionalHint"));
        }
        return tmpHints.isEmpty() ? aDescriptor.getName() : aDescriptor.getName() + " (" + String.join(", ", tmpHints) + ")";
    }
    //
    /**
     * Returns the name of the selection setting of the given descriptor. Underscores are replaced by spaces because
     * preference names must not contain them.
     *
     * @param aDescriptor descriptor
     * @return setting name
     */
    public static String getDescriptorSettingName(Descriptor aDescriptor) {
        return "Descriptor " + aDescriptor.name().replace('_', ' ');
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Reads the user-defined presets from the given file if it exists; unknown descriptor names are skipped.
     *
     * @param aFile presets file
     */
    private void reloadUserDescriptorPresets(File aFile) {
        if (!aFile.isFile()) {
            return;
        }
        Properties tmpProperties = new Properties();
        try (Reader tmpReader = Files.newBufferedReader(aFile.toPath(), StandardCharsets.UTF_8)) {
            tmpProperties.load(tmpReader);
        } catch (IOException anException) {
            ClusteringService.LOGGER.log(Level.WARNING, "Unable to reload the user-defined descriptor presets.", anException);
            return;
        }
        for (String tmpName : tmpProperties.stringPropertyNames()) {
            if (!ClusteringService.isValidUserPresetName(tmpName)) {
                continue;
            }
            Set<Descriptor> tmpDescriptors = EnumSet.noneOf(Descriptor.class);
            for (String tmpDescriptorName : tmpProperties.getProperty(tmpName).split(",")) {
                try {
                    Descriptor tmpDescriptor = Descriptor.valueOf(tmpDescriptorName.trim());
                    if (!tmpDescriptor.isFingerprint()) {
                        tmpDescriptors.add(tmpDescriptor);
                    }
                } catch (IllegalArgumentException anException) {
                    //descriptor no longer available or empty preset
                    ClusteringService.LOGGER.log(Level.INFO, "Unknown descriptor {0} in user preset skipped.", tmpDescriptorName);
                }
            }
            this.userDescriptorPresetMap.put(tmpName, tmpDescriptors);
        }
        this.updateDescriptorPresetKeyList();
    }
    //
    /**
     * Refills the observable list of preset keys: own presets first, then the built-in presets, then the custom key.
     * The currently selected key is kept, because a combo box bound to it may reset it while its items are replaced.
     */
    private void updateDescriptorPresetKeyList() {
        List<String> tmpKeys = new ArrayList<>();
        for (DescriptorPreset tmpPreset : this.getDescriptorPresets()) {
            tmpKeys.add(tmpPreset.key());
        }
        tmpKeys.add(ClusteringService.CUSTOM_PRESET_KEY);
        String tmpSelectedKey = this.descriptorPresetSetting == null ? null : this.descriptorPresetSetting.get();
        this.descriptorPresetKeyList.setAll(tmpKeys);
        if (tmpSelectedKey != null && !tmpSelectedKey.equals(this.descriptorPresetSetting.get())) {
            this.descriptorPresetSetting.set(tmpSelectedKey);
        }
    }
    //
    /**
     * Creates an enum setting, adds it to the settings list and registers its display name and tooltip.
     *
     * @param aTargetList settings list the new setting is added to
     * @param aName setting name (must be a valid preference name)
     * @param aDefault default value
     * @param anEnumClass enum class of the options
     * @param aMessageKey key of the setting in the message properties (without prefix and suffix)
     * @return new setting
     */
    private SimpleIDisplayEnumConstantProperty createEnumSetting(List<Property<?>> aTargetList, String aName, IDisplayEnum aDefault,
                                                                 Class<?> anEnumClass, String aMessageKey) {
        SimpleIDisplayEnumConstantProperty tmpSetting = new SimpleIDisplayEnumConstantProperty(this, aName, aDefault, anEnumClass) {
            @Override
            public void set(IDisplayEnum newValue) throws NullPointerException, IllegalArgumentException {
                try {
                    super.set(newValue);
                } catch (NullPointerException | IllegalArgumentException anException) {
                    ClusteringService.handleIllegalSettingValue(anException);
                    //re-throws the exception to properly reset the binding
                    throw anException;
                }
            }
        };
        this.registerSetting(aTargetList, tmpSetting, aMessageKey);
        return tmpSetting;
    }
    //
    /**
     * Creates a boolean setting, adds it to the given settings list and registers its display name and tooltip.
     *
     * @param aTargetList settings list the new setting is added to
     * @param aName setting name (must be a valid preference name)
     * @param aDefault default value
     * @param aMessageKey key of the setting in the message properties (without prefix and suffix)
     * @return new setting
     */
    private SimpleBooleanProperty createBooleanSetting(List<Property<?>> aTargetList, String aName, boolean aDefault, String aMessageKey) {
        SimpleBooleanProperty tmpSetting = new SimpleBooleanProperty(this, aName, aDefault);
        this.registerSetting(aTargetList, tmpSetting, aMessageKey);
        return tmpSetting;
    }
    //
    /**
     * Creates an integer setting with the given inclusive range, adds it to the settings list and registers its display
     * name and tooltip.
     *
     * @param aTargetList settings list the new setting is added to
     * @param aName setting name (must be a valid preference name)
     * @param aDefault default value
     * @param aMinimum minimum allowed value
     * @param aMaximum maximum allowed value
     * @param aMessageKey key of the setting in the message properties (without prefix and suffix)
     * @return new setting
     */
    private SimpleIntegerProperty createIntegerSetting(List<Property<?>> aTargetList, String aName, int aDefault, int aMinimum,
                                                       int aMaximum, String aMessageKey) {
        SimpleIntegerProperty tmpSetting = new SimpleIntegerProperty(this, aName, aDefault) {
            @Override
            public void set(int newValue) throws IllegalArgumentException {
                if (newValue < aMinimum || newValue > aMaximum) {
                    IllegalArgumentException tmpException = new IllegalArgumentException(String.format(
                            "%s must be in [%d, %d] but is %d.", aName, aMinimum, aMaximum, newValue));
                    ClusteringService.handleIllegalSettingValue(tmpException);
                    throw tmpException;
                }
                super.set(newValue);
            }
        };
        this.registerSetting(aTargetList, tmpSetting, aMessageKey);
        return tmpSetting;
    }
    //
    /**
     * Creates a double setting with the given inclusive range, adds it to the settings list and registers its display
     * name and tooltip.
     *
     * @param aTargetList settings list the new setting is added to
     * @param aName setting name (must be a valid preference name)
     * @param aDefault default value
     * @param aMinimum minimum allowed value
     * @param aMaximum maximum allowed value
     * @param aMessageKey key of the setting in the message properties (without prefix and suffix)
     * @return new setting
     */
    private SimpleDoubleProperty createDoubleSetting(List<Property<?>> aTargetList, String aName, double aDefault, double aMinimum,
                                                     double aMaximum, String aMessageKey) {
        SimpleDoubleProperty tmpSetting = new SimpleDoubleProperty(this, aName, aDefault) {
            @Override
            public void set(double newValue) throws IllegalArgumentException {
                if (!(newValue >= aMinimum && newValue <= aMaximum)) {
                    IllegalArgumentException tmpException = new IllegalArgumentException(String.format(
                            "%s must be in [%s, %s] but is %s.", aName, aMinimum, aMaximum, newValue));
                    ClusteringService.handleIllegalSettingValue(tmpException);
                    throw tmpException;
                }
                super.set(newValue);
            }
        };
        this.registerSetting(aTargetList, tmpSetting, aMessageKey);
        return tmpSetting;
    }
    //
    /**
     * Adds the given setting to the given settings list and registers its display name and tooltip from the message
     * properties.
     *
     * @param aTargetList settings list the setting is added to
     * @param aSetting setting
     * @param aMessageKey key of the setting in the message properties (without prefix and suffix)
     */
    private void registerSetting(List<Property<?>> aTargetList, Property<?> aSetting, String aMessageKey) {
        aTargetList.add(aSetting);
        this.settingNameToDisplayNameMap.put(aSetting.getName(), Message.get("ClusteringService." + aMessageKey + ".displayName"));
        this.settingNameToTooltipTextMap.put(aSetting.getName(), Message.get("ClusteringService." + aMessageKey + ".tooltip"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Joins the given descriptor column names to the value of the displayed descriptor columns setting; every name is
     * followed by the separator, names that cannot be persisted are left out.
     *
     * @param aNameList descriptor column names
     * @return setting value, the separator alone for an empty list
     */
    private static String joinDisplayedDescriptorColumnNames(List<String> aNameList) {
        StringBuilder tmpBuilder = new StringBuilder();
        for (String tmpName : aNameList) {
            if (tmpName == null || tmpName.isBlank() || tmpName.contains(ClusteringService.DISPLAYED_DESCRIPTOR_COLUMNS_SEPARATOR)
                    || !SingleTermPreference.isValidContent(tmpName)) {
                ClusteringService.LOGGER.log(Level.INFO, "Descriptor column {0} cannot be persisted and is not displayed.", tmpName);
                continue;
            }
            tmpBuilder.append(tmpName).append(ClusteringService.DISPLAYED_DESCRIPTOR_COLUMNS_SEPARATOR);
        }
        return tmpBuilder.isEmpty() ? ClusteringService.DISPLAYED_DESCRIPTOR_COLUMNS_SEPARATOR : tmpBuilder.toString();
    }
    //
    /**
     * Creates the descriptor preset object of a user-defined preset.
     *
     * @param aName preset name
     * @param aDescriptors descriptors of the preset
     * @return descriptor preset
     */
    private static DescriptorPreset createUserPreset(String aName, Set<Descriptor> aDescriptors) {
        return new DescriptorPreset(DescriptorPreset.getUserPresetKey(aName), aName,
                String.format(Message.get("ClusteringService.userPreset.tooltip"), aDescriptors.size()), aDescriptors, true);
    }
    //
    /**
     * Updates the given settings from the given preference container file if it exists.
     *
     * @param aSettingsList settings to update
     * @param aFile preference container file
     */
    private static void reloadSettings(List<Property<?>> aSettingsList, File aFile) {
        if (!aFile.isFile() || !aFile.canRead()) {
            ClusteringService.LOGGER.log(Level.INFO, "No persisted clustering settings available in {0}.", aFile.getName());
            return;
        }
        try {
            PreferenceUtil.updatePropertiesFromPreferences(aSettingsList, new PreferenceContainer(aFile));
        } catch (IllegalArgumentException | IOException anException) {
            ClusteringService.LOGGER.log(Level.WARNING, String.format("Unable to reload clustering settings from %s: %s",
                    aFile.getName(), anException.toString()), anException);
        }
    }
    //
    /**
     * Logs an illegal setting value and shows it to the user, as done for the fragmenter settings.
     *
     * @param anException exception describing the illegal value
     */
    private static void handleIllegalSettingValue(RuntimeException anException) {
        ClusteringService.LOGGER.log(Level.WARNING, anException.toString(), anException);
        GuiUtil.guiExceptionAlert(Message.get("Fragmenter.IllegalSettingValue.Title"),
                Message.get("Fragmenter.IllegalSettingValue.Header"),
                anException.toString(),
                anException);
    }
    //</editor-fold>
}
