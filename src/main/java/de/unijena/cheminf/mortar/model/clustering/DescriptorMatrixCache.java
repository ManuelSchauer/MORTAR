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
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.descriptor.DescriptorCalculator;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.util.ChemUtil;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Cache of calculated descriptor values, so that descriptors are not calculated again, e.g. when the molecules are
 * clustered after a descriptor calculation with the same descriptors or clustered again with another vigilance. The
 * cache holds the values of one molecule set calculated with one aromaticity model, per descriptor as a block of
 * columns. A request takes the cached blocks and calculates only the missing descriptors.
 * <ul>
 *     <li>A request for a subset of the cached molecules (identical data model objects) reuses the cached rows; if
 *     descriptors are missing for such a subset, or other molecules or another aromaticity model are requested, the
 *     cache is reset to the requested molecules.</li>
 *     <li>A changed length of the circular fingerprint discards the cached fingerprint blocks.</li>
 * </ul>
 * The methods are synchronized, since the calculations run in background threads.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class DescriptorMatrixCache {
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Descriptor values of the cached molecules, one block of columns per descriptor.
     */
    private final Map<Descriptor, FeatureMatrix> blockOfDescriptorMap;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Cached molecules in the order of the rows of the blocks; empty if nothing is cached.
     */
    private List<MoleculeDataModel> moleculeList;
    /**
     * Row index of every cached molecule, by identity of the data model.
     */
    private Map<MoleculeDataModel, Integer> rowIndexOfMoleculeMap;
    /**
     * Electron donation model of the aromaticity model the cached values were calculated with.
     */
    private IMoleculeFragmenter.ElectronDonationModelOption electronDonationModelOption;
    /**
     * Cycle finder of the aromaticity model the cached values were calculated with.
     */
    private IMoleculeFragmenter.CycleFinderOption cycleFinderOption;
    /**
     * Size of the circular fingerprints the cached fingerprint blocks were calculated with.
     */
    private int circularFingerprintSize;
    /**
     * True if the cached circular fingerprint blocks were calculated as count fingerprints.
     */
    private boolean isCircularFingerprintCounts;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, creates an empty cache.
     */
    public DescriptorMatrixCache() {
        this.blockOfDescriptorMap = new EnumMap<>(Descriptor.class);
        this.moleculeList = List.of();
        this.rowIndexOfMoleculeMap = new IdentityHashMap<>();
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the descriptor matrix of the given molecules and descriptors, with the columns in the order of the given
     * descriptors like {@link FeatureMatrixBuilder#buildDescriptorMatrix}. Cached descriptors are reused, only the
     * missing ones are calculated and added to the cache.
     *
     * @param aMoleculeList molecules (rows of the matrix)
     * @param aDescriptorArray descriptors to calculate
     * @param anElectronDonationModelOption electron donation model of the aromaticity model
     * @param aCycleFinderOption cycle finder of the aromaticity model
     * @param aDescriptorCalculator calculator for the missing descriptors; the cached fingerprint blocks are discarded
     *                              if its circular fingerprint configuration (size, bits or counts) differs from the
     *                              one they were calculated with
     * @param anIsParallelCalculation true: molecules are processed in parallel
     * @return feature matrix of the descriptor block
     * @throws NullPointerException if a parameter is null
     * @throws InterruptedException if the calculation is interrupted
     */
    public synchronized FeatureMatrix getDescriptorMatrix(List<MoleculeDataModel> aMoleculeList, Descriptor[] aDescriptorArray,
                                                          IMoleculeFragmenter.ElectronDonationModelOption anElectronDonationModelOption,
                                                          IMoleculeFragmenter.CycleFinderOption aCycleFinderOption,
                                                          DescriptorCalculator aDescriptorCalculator,
                                                          boolean anIsParallelCalculation)
            throws NullPointerException, InterruptedException {
        Objects.requireNonNull(aMoleculeList, "Given molecule list is null.");
        Objects.requireNonNull(aDescriptorArray, "Given descriptor array is null.");
        Objects.requireNonNull(aDescriptorCalculator, "Given descriptor calculator is null.");
        Objects.requireNonNull(anElectronDonationModelOption, "Given electron donation model is null.");
        Objects.requireNonNull(aCycleFinderOption, "Given cycle finder is null.");
        if (aDescriptorArray.length == 0) {
            return new FeatureMatrix(new float[aMoleculeList.size()][0], FeatureMatrixBuilder.getRowNames(aMoleculeList),
                    new String[0], new FeatureType[0], new String[0]);
        }
        if (anElectronDonationModelOption != this.electronDonationModelOption || aCycleFinderOption != this.cycleFinderOption
                || !this.containsAll(aMoleculeList)) {
            this.reset(aMoleculeList, anElectronDonationModelOption, aCycleFinderOption);
        }
        if (aDescriptorCalculator.getCircularFingerprintSize() != this.circularFingerprintSize
                || aDescriptorCalculator.isCircularFingerprintCounts() != this.isCircularFingerprintCounts) {
            this.blockOfDescriptorMap.keySet().removeIf(Descriptor::isFingerprint);
            this.circularFingerprintSize = aDescriptorCalculator.getCircularFingerprintSize();
            this.isCircularFingerprintCounts = aDescriptorCalculator.isCircularFingerprintCounts();
        }
        List<Descriptor> tmpMissingDescriptors = new ArrayList<>();
        for (Descriptor tmpDescriptor : aDescriptorArray) {
            if (!this.blockOfDescriptorMap.containsKey(tmpDescriptor) && !tmpMissingDescriptors.contains(tmpDescriptor)) {
                tmpMissingDescriptors.add(tmpDescriptor);
            }
        }
        //missing descriptors are calculated for the cached molecules only if exactly these are requested, so that no
        // values are calculated for molecules that are not needed
        if (!tmpMissingDescriptors.isEmpty() && aMoleculeList.size() != this.moleculeList.size()) {
            this.reset(aMoleculeList, anElectronDonationModelOption, aCycleFinderOption);
            tmpMissingDescriptors.clear();
            for (Descriptor tmpDescriptor : aDescriptorArray) {
                if (!tmpMissingDescriptors.contains(tmpDescriptor)) {
                    tmpMissingDescriptors.add(tmpDescriptor);
                }
            }
        }
        if (!tmpMissingDescriptors.isEmpty()) {
            FeatureMatrix tmpCalculated = FeatureMatrixBuilder.buildDescriptorMatrix(this.moleculeList,
                    tmpMissingDescriptors.toArray(new Descriptor[0]),
                    ChemUtil.createAromaticityModel(anElectronDonationModelOption, aCycleFinderOption), aDescriptorCalculator,
                    anIsParallelCalculation);
            this.addBlocks(tmpCalculated, tmpMissingDescriptors, aDescriptorCalculator);
        }
        //the requested rows of every block, in the order of the requested descriptors
        int[] tmpRowIndices = new int[aMoleculeList.size()];
        boolean tmpIsIdentity = aMoleculeList.size() == this.moleculeList.size();
        for (int i = 0; i < tmpRowIndices.length; i++) {
            tmpRowIndices[i] = this.rowIndexOfMoleculeMap.get(aMoleculeList.get(i));
            tmpIsIdentity = tmpIsIdentity && tmpRowIndices[i] == i;
        }
        FeatureMatrix[] tmpBlocks = new FeatureMatrix[aDescriptorArray.length];
        for (int d = 0; d < aDescriptorArray.length; d++) {
            FeatureMatrix tmpBlock = this.blockOfDescriptorMap.get(aDescriptorArray[d]);
            tmpBlocks[d] = tmpIsIdentity ? tmpBlock : tmpBlock.selectRows(tmpRowIndices);
        }
        return FeatureMatrix.concatenate(tmpBlocks);
    }
    //
    /**
     * Returns whether the values of the given descriptor are cached for all given molecules with the given aromaticity
     * model, i.e. whether a request would not calculate it.
     *
     * @param aMoleculeList molecules
     * @param aDescriptor descriptor
     * @param anElectronDonationModelOption electron donation model of the aromaticity model
     * @param aCycleFinderOption cycle finder of the aromaticity model
     * @return true if cached
     */
    public synchronized boolean isCached(List<MoleculeDataModel> aMoleculeList, Descriptor aDescriptor,
                                         IMoleculeFragmenter.ElectronDonationModelOption anElectronDonationModelOption,
                                         IMoleculeFragmenter.CycleFinderOption aCycleFinderOption) {
        return anElectronDonationModelOption == this.electronDonationModelOption && aCycleFinderOption == this.cycleFinderOption
                && this.containsAll(aMoleculeList) && this.blockOfDescriptorMap.containsKey(aDescriptor);
    }
    //
    /**
     * Discards all cached values, e.g. when new molecules are imported.
     */
    public synchronized void clear() {
        this.reset(List.of(), null, null);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns whether all given molecules are cached (by identity).
     *
     * @param aMoleculeList molecules
     * @return true if every molecule is cached
     */
    private boolean containsAll(List<MoleculeDataModel> aMoleculeList) {
        if (this.moleculeList.isEmpty()) {
            return aMoleculeList.isEmpty();
        }
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            if (!this.rowIndexOfMoleculeMap.containsKey(tmpMolecule)) {
                return false;
            }
        }
        return true;
    }
    //
    /**
     * Discards all blocks and caches the given molecules and aromaticity model.
     *
     * @param aMoleculeList molecules to cache values of
     * @param anElectronDonationModelOption electron donation model, null for an empty cache
     * @param aCycleFinderOption cycle finder, null for an empty cache
     */
    private void reset(List<MoleculeDataModel> aMoleculeList, IMoleculeFragmenter.ElectronDonationModelOption anElectronDonationModelOption,
                       IMoleculeFragmenter.CycleFinderOption aCycleFinderOption) {
        this.blockOfDescriptorMap.clear();
        this.moleculeList = List.copyOf(aMoleculeList);
        this.rowIndexOfMoleculeMap = new IdentityHashMap<>(aMoleculeList.size());
        for (int i = 0; i < this.moleculeList.size(); i++) {
            this.rowIndexOfMoleculeMap.put(this.moleculeList.get(i), i);
        }
        this.electronDonationModelOption = anElectronDonationModelOption;
        this.cycleFinderOption = aCycleFinderOption;
    }
    //
    /**
     * Splits the given matrix into one block per descriptor (by the feature group name, which is the descriptor enum
     * constant name) and caches the blocks.
     *
     * @param aMatrix calculated descriptor matrix of the cached molecules
     * @param aDescriptorList descriptors of the matrix in column order
     * @param aDescriptorCalculator calculator the matrix was calculated with, which gives the number of columns of every
     *                              descriptor
     */
    private void addBlocks(FeatureMatrix aMatrix, List<Descriptor> aDescriptorList, DescriptorCalculator aDescriptorCalculator) {
        int tmpColumn = 0;
        for (Descriptor tmpDescriptor : aDescriptorList) {
            int tmpNumberOfComponents = aDescriptorCalculator.getNumberOfComponents(tmpDescriptor);
            int[] tmpColumnIndices = new int[tmpNumberOfComponents];
            for (int c = 0; c < tmpNumberOfComponents; c++) {
                tmpColumnIndices[c] = tmpColumn++;
            }
            this.blockOfDescriptorMap.put(tmpDescriptor, aMatrix.selectColumns(tmpColumnIndices));
        }
    }
    //</editor-fold>
}
