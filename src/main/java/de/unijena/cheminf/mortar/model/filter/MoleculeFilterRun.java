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

import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.EnumMap;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Applies the molecule filters to the molecules of one import run and counts the outcomes. Create one instance per
 * import run, with a settings snapshot taken at its start, and use it for all molecules of that run. The filters are
 * applied in a fixed order:
 * <ol>
 *     <li>{@link #applyStructuralFilters(IAtomContainer)}: disconnected structures (keep, reduce to the largest
 *     component, or remove), then elements, then heavy atom count; the latter two judge the possibly reduced
 *     molecule;</li>
 *     <li>the importer generates the unique SMILES string of the molecule returned by step 1;</li>
 *     <li>{@link #isDuplicate(String)} with that SMILES string (or {@link #registerSmilesGenerationFailure()} if none
 *     could be generated); the first occurrence of a SMILES string is kept.</li>
 * </ol>
 * Consequently, with "reduce to the largest component", a counter-ion does not cause removal by the element filter,
 * and different salt forms of one compound become duplicates of each other.
 * <br>With the default settings (all filters inactive), the molecules are passed through unchanged and uninspected,
 * so an import behaves exactly as without molecule filters; only the counting takes place.
 * <br>Odd input: an empty molecule (no atoms) is rejected with {@link MoleculeFilterReason#EMPTY_MOLECULE} if at least
 * one structural filter is active, otherwise it is passed on unchanged as before. If an active structural filter
 * throws an exception (or a stack overflow error) for a molecule, the molecule is removed (fail-closed, because it is
 * unknown whether it meets the criteria), counted with {@link MoleculeFilterReason#FILTER_CHECK_FAILED}, and the
 * exception is logged; no exception is propagated.
 * <br>Thread safety: the structural filtering is stateless and may run in parallel; the counting, the duplicate check,
 * and {@link #getSummary()} are synchronised, so the summary may also be retrieved by another thread (e.g. the GUI
 * thread) during the run. Note that with parallel callers, "first occurrence" means the first call of
 * {@link #isDuplicate(String)}, not the first molecule in the file. No GUI code is called.
 * <br>Removed molecules are not logged individually by this class, so that every removal is logged in one place only:
 * the importer logs one line per removed molecule with its name and the reason (see
 * {@link Importer#importMoleculeFile(java.io.File, boolean, boolean, boolean, MoleculeFilterSettingsSnapshot)}); the
 * summary can be logged via {@link MoleculeFilterSummary#toString()}. Only failed filter checks are logged here, with
 * the exception (see above).
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class MoleculeFilterRun {
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(MoleculeFilterRun.class.getName());

    /**
     * Initial capacity of the set of registered SMILES strings.
     */
    private static final int INITIAL_SMILES_SET_CAPACITY = CollectionUtil.calculateInitialHashCollectionCapacity(
            1000, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Settings of this run.
     */
    private final MoleculeFilterSettingsSnapshot settingsSnapshot;

    /**
     * Unique SMILES strings of the molecules kept so far; only filled if the duplicate filter is active. Guarded by
     * this.
     */
    private final Set<String> registeredSmilesSet;

    /**
     * Number of removed molecules per reason. Guarded by this.
     */
    private final EnumMap<MoleculeFilterReason, Integer> removedCountsMap;
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * Number of molecules passed to the structural filters. Guarded by this.
     */
    private int readCount;

    /**
     * Number of molecules that passed all filters. Guarded by this.
     */
    private int keptCount;

    /**
     * Number of molecules reduced to their largest component. Guarded by this.
     */
    private int reducedToLargestComponentCount;
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * Constructor for a run with the given settings snapshot.
     *
     * @param aSettingsSnapshot the settings to apply during the whole run
     * @throws NullPointerException if the given snapshot is null
     */
    public MoleculeFilterRun(MoleculeFilterSettingsSnapshot aSettingsSnapshot) throws NullPointerException {
        Objects.requireNonNull(aSettingsSnapshot, "Given settings snapshot is null.");
        this.settingsSnapshot = aSettingsSnapshot;
        this.registeredSmilesSet = new HashSet<>(MoleculeFilterRun.INITIAL_SMILES_SET_CAPACITY,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.removedCountsMap = new EnumMap<>(MoleculeFilterReason.class);
        for (MoleculeFilterReason tmpReason : MoleculeFilterReason.values()) {
            this.removedCountsMap.put(tmpReason, 0);
        }
        this.readCount = 0;
        this.keptCount = 0;
        this.reducedToLargestComponentCount = 0;
    }

    /**
     * Constructor for a run with a snapshot of the given settings, taken now.
     *
     * @param aSettings the settings to take the snapshot of
     * @throws NullPointerException if the given settings are null
     */
    public MoleculeFilterRun(MoleculeFilterSettings aSettings) throws NullPointerException {
        this(Objects.requireNonNull(aSettings, "Given settings are null.").createSnapshot());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the settings snapshot of this run.
     *
     * @return settings snapshot
     */
    public MoleculeFilterSettingsSnapshot getSettingsSnapshot() {
        return this.settingsSnapshot;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Applies the structural filters to the given molecule in the order disconnected structures, elements, heavy atom
     * count, and counts it as read. Call this exactly once per molecule of the import, after the importer's
     * preprocessing. If the returned result is kept, continue with {@link MoleculeFilterResult#getMolecule()} (which
     * may be a new, reduced container) and finally call {@link #isDuplicate(String)} or
     * {@link #registerSmilesGenerationFailure()} for it. If none of the structural filters is active, the given
     * container is returned unchanged without inspection.
     *
     * @param aMolecule the preprocessed molecule
     * @return the result, never null
     * @throws NullPointerException if the given molecule is null
     */
    public MoleculeFilterResult applyStructuralFilters(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        MoleculeFilterResult tmpResult = this.determineStructuralFilterResult(aMolecule);
        synchronized (this) {
            this.readCount++;
            if (tmpResult.isReducedToLargestComponent()) {
                this.reducedToLargestComponentCount++;
            }
            if (!tmpResult.isKept()) {
                this.removedCountsMap.merge(tmpResult.getRejectionReason(), 1, Integer::sum);
            }
        }
        return tmpResult;
    }

    /**
     * Duplicate check, the last filter step for a molecule that passed the structural filters and for which the
     * importer generated the given unique SMILES string. If the duplicate filter is active and the SMILES string was
     * registered before in this run, the molecule is counted as duplicate and true is returned; otherwise, the SMILES
     * string is registered (only if the filter is active), the molecule is counted as kept, and false is returned.
     * Call this exactly once per such molecule. The comparison is a plain string comparison, so whether stereoisomers
     * or aromatic and Kekule notations of one compound count as duplicates depends on how the importer generated the
     * SMILES strings (stereochemistry and Kekule settings).
     *
     * @param aUniqueSmiles the unique SMILES string generated by the importer
     * @return true if the molecule is a duplicate and must not be taken over; false if it is kept
     * @throws NullPointerException if the given SMILES string is null; use {@link #registerSmilesGenerationFailure()}
     * in that case
     * @throws IllegalArgumentException if the given SMILES string is blank; use
     * {@link #registerSmilesGenerationFailure()} in that case
     */
    public synchronized boolean isDuplicate(String aUniqueSmiles) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aUniqueSmiles, "Given SMILES string is null.");
        if (aUniqueSmiles.isBlank()) {
            throw new IllegalArgumentException("Given SMILES string is blank.");
        }
        if (this.settingsSnapshot.isDuplicateFilterActive() && !this.registeredSmilesSet.add(aUniqueSmiles)) {
            this.removedCountsMap.merge(MoleculeFilterReason.DUPLICATE, 1, Integer::sum);
            return true;
        }
        this.keptCount++;
        return false;
    }

    /**
     * Registers that the importer could not generate a SMILES string for a molecule that passed the structural
     * filters, so that the summary is complete. Call this instead of {@link #isDuplicate(String)} for such a molecule.
     */
    public synchronized void registerSmilesGenerationFailure() {
        this.removedCountsMap.merge(MoleculeFilterReason.SMILES_GENERATION_FAILED, 1, Integer::sum);
    }

    /**
     * Returns an immutable summary of the counts so far.
     *
     * @return summary of this run
     */
    public synchronized MoleculeFilterSummary getSummary() {
        return new MoleculeFilterSummary(this.settingsSnapshot, this.readCount, this.keptCount,
                this.reducedToLargestComponentCount, this.removedCountsMap);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Applies the structural filters without counting, see {@link #applyStructuralFilters(IAtomContainer)}.
     *
     * @param aMolecule the molecule
     * @return the result
     */
    private MoleculeFilterResult determineStructuralFilterResult(IAtomContainer aMolecule) {
        if (!this.settingsSnapshot.isAnyStructuralFilterActive()) {
            return MoleculeFilterResult.createKeptResult(aMolecule, false);
        }
        if (aMolecule.isEmpty()) {
            return MoleculeFilterResult.createRejectedResult(MoleculeFilterReason.EMPTY_MOLECULE, false);
        }
        IAtomContainer tmpMolecule = aMolecule;
        boolean tmpIsReduced = false;
        try {
            switch (this.settingsSnapshot.disconnectedStructuresOption()) {
                case MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED -> {
                    //nothing to do
                }
                case MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_LARGEST_COMPONENT -> {
                    IAtomContainer tmpLargestComponent = MoleculeFilterUtil.extractLargestComponent(tmpMolecule);
                    if (tmpLargestComponent != tmpMolecule) {
                        tmpMolecule = tmpLargestComponent;
                        tmpIsReduced = true;
                    }
                }
                case MoleculeFilterSettings.DisconnectedStructuresOption.REMOVE_MOLECULE -> {
                    if (!MoleculeFilterUtil.isConnected(tmpMolecule)) {
                        return MoleculeFilterResult.createRejectedResult(MoleculeFilterReason.DISCONNECTED_STRUCTURE, false);
                    }
                }
            }
            if (this.settingsSnapshot.isElementFilterActive()
                    && !MoleculeFilterUtil.consistsOfAllowedElementsOnly(tmpMolecule)) {
                return MoleculeFilterResult.createRejectedResult(MoleculeFilterReason.DISALLOWED_ELEMENT, tmpIsReduced);
            }
            if (this.settingsSnapshot.isHeavyAtomCountFilterActive()) {
                int tmpHeavyAtomCount = MoleculeFilterUtil.countHeavyAtoms(tmpMolecule);
                if (tmpHeavyAtomCount < this.settingsSnapshot.minimumHeavyAtomCount()
                        || tmpHeavyAtomCount > this.settingsSnapshot.maximumHeavyAtomCount()) {
                    return MoleculeFilterResult.createRejectedResult(MoleculeFilterReason.HEAVY_ATOM_COUNT_OUT_OF_RANGE,
                            tmpIsReduced);
                }
            }
        } catch (Exception | StackOverflowError anError) {
            //stack overflow errors are caught as well because CDK's connected components search is recursive and might
            // overflow for huge molecules; the molecule is removed because it is unknown whether it meets the criteria
            MoleculeFilterRun.LOGGER.log(Level.WARNING, String.format("A molecule filter check failed for molecule %s, " +
                    "the molecule is removed: %s", aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY),
                    anError.toString()), anError);
            return MoleculeFilterResult.createRejectedResult(MoleculeFilterReason.FILTER_CHECK_FAILED, tmpIsReduced);
        }
        return MoleculeFilterResult.createKeptResult(tmpMolecule, tmpIsReduced);
    }
    //</editor-fold>
}
