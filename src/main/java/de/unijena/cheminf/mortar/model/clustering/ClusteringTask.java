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

import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorCalculator;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import javafx.concurrent.Task;

import java.util.List;
import java.util.Objects;
import java.util.concurrent.CancellationException;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;

/**
 * JavaFX task that runs the complete clustering pipeline in the background: building the feature matrix (fragment
 * fingerprint and descriptors), transformation, filtering and ART-2A clustering. Progress and status messages are
 * reported via the task properties. Cancellation is checked between the steps; a running descriptor calculation
 * or ART-2A clustering itself cannot be interrupted by the libraries, so cancellation takes effect after that step.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringTask extends Task<ClusteringRunResult> {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Number of pipeline steps, used for the progress.
     */
    private static final int NUMBER_OF_STEPS = 5;
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Molecules to cluster.
     */
    private final List<MoleculeDataModel> moleculeList;
    /**
     * Parameters of the run.
     */
    private final ClusteringParameters parameters;
    /**
     * Priority table for the correlation filter.
     */
    private final DescriptorPriorityTable priorityTable;
    /**
     * Cache of calculated descriptor values, null to calculate all descriptors.
     */
    private final DescriptorMatrixCache descriptorMatrixCache;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aMoleculeList molecules to cluster
     * @param aParameters parameters of the run
     * @param aPriorityTable priority table for the correlation filter
     * @throws NullPointerException if a parameter is null
     */
    public ClusteringTask(List<MoleculeDataModel> aMoleculeList, ClusteringParameters aParameters,
                          DescriptorPriorityTable aPriorityTable) throws NullPointerException {
        this(aMoleculeList, aParameters, aPriorityTable, null);
    }
    //
    /**
     * Constructor with a cache of calculated descriptor values, from which cached descriptors are taken instead of being
     * calculated again.
     *
     * @param aMoleculeList molecules to cluster
     * @param aParameters parameters of the run
     * @param aPriorityTable priority table for the correlation filter
     * @param aDescriptorMatrixCache cache of descriptor values, null to calculate all descriptors
     * @throws NullPointerException if the molecules, parameters or priority table are null
     */
    public ClusteringTask(List<MoleculeDataModel> aMoleculeList, ClusteringParameters aParameters,
                          DescriptorPriorityTable aPriorityTable, DescriptorMatrixCache aDescriptorMatrixCache) throws NullPointerException {
        this.moleculeList = List.copyOf(Objects.requireNonNull(aMoleculeList, "Given molecule list is null."));
        this.parameters = Objects.requireNonNull(aParameters, "Given parameters are null.");
        this.priorityTable = Objects.requireNonNull(aPriorityTable, "Given priority table is null.");
        this.descriptorMatrixCache = aDescriptorMatrixCache;
    }
    //</editor-fold>
    //
    //<editor-fold desc="overridden methods" defaultstate="collapsed">
    /**
     * Runs the clustering pipeline.
     *
     * @return result of the run
     * @throws Exception if a step fails or the task was cancelled
     */
    @Override
    protected ClusteringRunResult call() throws Exception {
        return ClusteringTask.runPipeline(this.moleculeList, this.parameters, this.priorityTable, this.descriptorMatrixCache, this::isCancelled,
                (aStep, aMessage) -> {
                    this.updateMessage(aMessage);
                    this.updateProgress(aStep, ClusteringTask.NUMBER_OF_STEPS);
                });
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Runs the clustering pipeline without JavaFX dependencies, so it can be used and tested outside the GUI.
     *
     * @param aMoleculeList molecules to cluster, at least one
     * @param aParameters parameters of the run
     * @param aPriorityTable priority table for the correlation filter
     * @param anIsCancelledSupplier returns true if the run should be stopped
     * @param aProgressConsumer receives the number of the step that starts and a status message
     * @return result of the run
     * @throws IllegalArgumentException if no molecule is given
     * @throws IllegalStateException if fewer than two features remain after filtering
     * @throws CancellationException if the run was cancelled
     * @throws Exception if the descriptor calculation or the clustering fails
     */
    public static ClusteringRunResult runPipeline(
            List<MoleculeDataModel> aMoleculeList,
            ClusteringParameters aParameters,
            DescriptorPriorityTable aPriorityTable,
            BooleanSupplier anIsCancelledSupplier,
            BiConsumer<Integer, String> aProgressConsumer)
            throws IllegalArgumentException, IllegalStateException, CancellationException, Exception {
        return ClusteringTask.runPipeline(aMoleculeList, aParameters, aPriorityTable, null, anIsCancelledSupplier, aProgressConsumer);
    }
    //
    /**
     * Runs the clustering pipeline like {@link #runPipeline(List, ClusteringParameters, DescriptorPriorityTable,
     * BooleanSupplier, BiConsumer)}, taking descriptors that are already in the given cache from it instead of
     * calculating them again; the calculated ones are added to the cache.
     *
     * @param aMoleculeList molecules to cluster, at least one
     * @param aParameters parameters of the run
     * @param aPriorityTable priority table for the correlation filter
     * @param aDescriptorMatrixCache cache of descriptor values, null to calculate all descriptors
     * @param anIsCancelledSupplier returns true if the run should be stopped
     * @param aProgressConsumer receives the number of the step that starts and a status message
     * @return result of the run
     * @throws IllegalArgumentException if no molecule is given
     * @throws IllegalStateException if fewer than two features remain after filtering
     * @throws CancellationException if the run was cancelled
     * @throws Exception if the descriptor calculation or the clustering fails
     */
    public static ClusteringRunResult runPipeline(
            List<MoleculeDataModel> aMoleculeList,
            ClusteringParameters aParameters,
            DescriptorPriorityTable aPriorityTable,
            DescriptorMatrixCache aDescriptorMatrixCache,
            BooleanSupplier anIsCancelledSupplier,
            BiConsumer<Integer, String> aProgressConsumer)
            throws IllegalArgumentException, IllegalStateException, CancellationException, Exception {
        if (aMoleculeList.isEmpty()) {
            throw new IllegalArgumentException("No molecules to cluster.");
        }
        //step 1: fragment fingerprint block
        aProgressConsumer.accept(0, "Building fragment fingerprints...");
        FeatureMatrix tmpFragmentMatrix = FeatureMatrixBuilder.buildFragmentFingerprintMatrix(aMoleculeList,
                aParameters.fragmentationName(), aParameters.isCountFingerprint(), aParameters.minimumFragmentMoleculeFrequency(),
                aParameters.maximumNumberOfFragmentFeatures());
        ClusteringTask.checkCancellation(anIsCancelledSupplier);
        //step 2: descriptor block, calculated on the molecules
        aProgressConsumer.accept(1, "Calculating descriptors...");
        //note: the circular fingerprint size is a global setting of the descriptor engine
        if (DescriptorCalculator.getCircularFingerprintSize() != aParameters.circularFingerprintSize()) {
            DescriptorCalculator.setCircularFingerprintSize(aParameters.circularFingerprintSize());
        }
        FeatureMatrix tmpDescriptorMatrix = aDescriptorMatrixCache == null
                ? FeatureMatrixBuilder.buildDescriptorMatrix(aMoleculeList, aParameters.descriptors(),
                        ChemUtil.createAromaticityModel(aParameters.electronDonationModelOption(), aParameters.cycleFinderOption()),
                        aParameters.isParallelCalculation())
                : aDescriptorMatrixCache.getDescriptorMatrix(aMoleculeList, aParameters.descriptors(),
                        aParameters.electronDonationModelOption(), aParameters.cycleFinderOption(),
                        aParameters.circularFingerprintSize(), aParameters.isParallelCalculation());
        FeatureMatrix tmpRawMatrix = FeatureMatrix.concatenate(tmpFragmentMatrix, tmpDescriptorMatrix);
        ClusteringTask.checkCancellation(anIsCancelledSupplier);
        //step 3: transformation
        aProgressConsumer.accept(2, "Transforming features...");
        FeatureMatrix tmpTransformedMatrix = FeatureTransformer.transform(tmpRawMatrix,
                aParameters.countTransformation(), aParameters.continuousTransformation());
        ClusteringTask.checkCancellation(anIsCancelledSupplier);
        //step 4: filtering
        aProgressConsumer.accept(3, "Filtering features...");
        FeatureFilter.FilterResult tmpFilterResult = FeatureFilter.filter(tmpTransformedMatrix, aPriorityTable,
                aParameters.maximumNaNFraction(),
                //a threshold of 1 disables the correlation filter
                aParameters.isCorrelationFilterEnabled() ? aParameters.correlationThreshold() : 1.0,
                aParameters.correlationMethod());
        FeatureMatrix tmpProcessedMatrix = tmpFilterResult.filteredMatrix();
        if (tmpProcessedMatrix.getColumnCount() < 2) {
            throw new IllegalStateException(String.format(
                    "Only %d feature(s) left after filtering, at least 2 are needed for the clustering.",
                    tmpProcessedMatrix.getColumnCount()));
        }
        ClusteringTask.checkCancellation(anIsCancelledSupplier);
        //step 5: clustering
        aProgressConsumer.accept(4, "Clustering...");
        ClusteringResult tmpClusteringResult = Art2aClusterer.cluster(tmpProcessedMatrix, aParameters.algorithm(),
                aParameters.vigilance(), aParameters.maximumNumberOfClusters(), aParameters.maximumNumberOfEpochs(),
                aParameters.convergenceThreshold(), aParameters.learningParameter(),
                aParameters.offsetForContrastEnhancement(), aParameters.randomSeed(), aParameters.isParallelCalculation());
        ClusteringTask.checkCancellation(anIsCancelledSupplier);
        aProgressConsumer.accept(ClusteringTask.NUMBER_OF_STEPS, "Clustering finished.");
        return new ClusteringRunResult(aParameters, aMoleculeList, tmpRawMatrix, tmpTransformedMatrix,
                tmpProcessedMatrix, tmpFilterResult.removedFeatures(), tmpClusteringResult);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Throws a cancellation exception if the given supplier returns true.
     *
     * @param anIsCancelledSupplier returns true if the run should be stopped
     * @throws CancellationException if the run was cancelled
     */
    private static void checkCancellation(BooleanSupplier anIsCancelledSupplier) throws CancellationException {
        if (anIsCancelledSupplier.getAsBoolean()) {
            throw new CancellationException("Clustering was cancelled.");
        }
    }
    //</editor-fold>
}
