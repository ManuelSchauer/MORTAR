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

import de.unijena.cheminf.mortar.message.Message;

import java.util.EnumMap;
import java.util.Map;
import java.util.Objects;

/**
 * Immutable summary of one import run with molecule filters, created by {@link MoleculeFilterRun#getSummary()}. For a
 * completed run, the counts add up: read = kept + total removed, i.e. {@link #getPendingCount()} is 0. The count of
 * molecules reduced to their largest component is informative only; these molecules are part of the kept or removed
 * counts as well.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class MoleculeFilterSummary {
    //<editor-fold desc="Private final variables">
    /**
     * Settings the run was performed with.
     */
    private final MoleculeFilterSettingsSnapshot settingsSnapshot;

    /**
     * Number of molecules passed to the structural filters.
     */
    private final int readCount;

    /**
     * Number of molecules that passed all filters and were taken over into the data set.
     */
    private final int keptCount;

    /**
     * Number of molecules that were reduced to their largest component.
     */
    private final int reducedToLargestComponentCount;

    /**
     * Number of removed molecules per reason; contains all reasons.
     */
    private final EnumMap<MoleculeFilterReason, Integer> removedCountsMap;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, only used by {@link MoleculeFilterRun}.
     *
     * @param aSettingsSnapshot settings the run was performed with
     * @param aReadCount number of molecules passed to the structural filters
     * @param aKeptCount number of molecules taken over into the data set
     * @param aReducedCount number of molecules reduced to their largest component
     * @param aRemovedCountsMap number of removed molecules per reason; missing reasons are counted as 0; the map is copied
     * @throws NullPointerException if the settings snapshot or the map is null
     */
    MoleculeFilterSummary(MoleculeFilterSettingsSnapshot aSettingsSnapshot,
                          int aReadCount,
                          int aKeptCount,
                          int aReducedCount,
                          Map<MoleculeFilterReason, Integer> aRemovedCountsMap) throws NullPointerException {
        Objects.requireNonNull(aSettingsSnapshot, "Given settings snapshot is null.");
        Objects.requireNonNull(aRemovedCountsMap, "Given map of removed counts is null.");
        this.settingsSnapshot = aSettingsSnapshot;
        this.readCount = aReadCount;
        this.keptCount = aKeptCount;
        this.reducedToLargestComponentCount = aReducedCount;
        this.removedCountsMap = new EnumMap<>(MoleculeFilterReason.class);
        for (MoleculeFilterReason tmpReason : MoleculeFilterReason.values()) {
            this.removedCountsMap.put(tmpReason, aRemovedCountsMap.getOrDefault(tmpReason, 0));
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the settings the run was performed with.
     *
     * @return settings snapshot
     */
    public MoleculeFilterSettingsSnapshot getSettingsSnapshot() {
        return this.settingsSnapshot;
    }

    /**
     * Returns the number of molecules passed to the structural filters, i.e. the number of molecules read from the
     * file that reached the filtering step.
     *
     * @return number of read molecules
     */
    public int getReadCount() {
        return this.readCount;
    }

    /**
     * Returns the number of molecules that passed all filters and were taken over into the data set.
     *
     * @return number of kept molecules
     */
    public int getKeptCount() {
        return this.keptCount;
    }

    /**
     * Returns the number of molecules that were reduced to their largest component. These molecules are also counted
     * as kept or removed, depending on the later filters.
     *
     * @return number of reduced molecules
     */
    public int getReducedToLargestComponentCount() {
        return this.reducedToLargestComponentCount;
    }

    /**
     * Returns the number of molecules removed for the given reason.
     *
     * @param aReason the reason
     * @return number of molecules removed for this reason
     * @throws NullPointerException if the given reason is null
     */
    public int getRemovedCount(MoleculeFilterReason aReason) throws NullPointerException {
        Objects.requireNonNull(aReason, "Given reason is null.");
        return this.removedCountsMap.get(aReason);
    }

    /**
     * Returns the number of removed molecules over all reasons, including the molecules for which the importer could
     * not generate a SMILES string.
     *
     * @return total number of removed molecules
     */
    public int getTotalRemovedCount() {
        int tmpSum = 0;
        for (int tmpCount : this.removedCountsMap.values()) {
            tmpSum += tmpCount;
        }
        return tmpSum;
    }

    /**
     * Returns the number of molecules that passed the structural filters but whose final outcome (kept, duplicate, or
     * SMILES generation failed) was not registered with the run, i.e. read - kept - total removed. This is 0 for a
     * completed run; a positive value means that the run was interrupted or that the importer missed a call.
     *
     * @return number of molecules without a final outcome
     */
    public int getPendingCount() {
        return this.readCount - this.keptCount - this.getTotalRemovedCount();
    }

    /**
     * Returns the number of molecules removed by the molecule filters, i.e. the total number of removed molecules
     * without the empty molecules and the molecules for which the importer could not generate a SMILES string. These
     * are not imported without molecule filters either (an empty molecule has no SMILES string), so they are listed
     * in the counts per reason but do not count as changed by the filters. Always 0 if all filters are inactive.
     *
     * @return number of molecules removed by the filters
     */
    public int getRemovedByFiltersCount() {
        return this.getTotalRemovedCount()
                - this.getRemovedCount(MoleculeFilterReason.EMPTY_MOLECULE)
                - this.getRemovedCount(MoleculeFilterReason.SMILES_GENERATION_FAILED);
    }

    /**
     * Returns whether the molecule filters changed the imported data set, i.e. whether at least one molecule was removed
     * by a filter (see {@link #getRemovedByFiltersCount()}) or reduced to its largest component. Always false if all
     * filters are inactive, so that, e.g., no summary is shown to the user for an import without molecule filters, and
     * also false if only empty molecules or molecules without SMILES string were not imported.
     *
     * @return true if the filters removed or reduced at least one molecule
     */
    public boolean isDataSetChangedByFilters() {
        return this.getRemovedByFiltersCount() > 0 || this.reducedToLargestComponentCount > 0;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Returns a language-specific, multi-line description of the counts for display in the GUI: the numbers of read,
     * imported, and reduced molecules, the total number of removed molecules, and one indented line per reason with
     * at least one removed molecule, in the order of {@link MoleculeFilterReason}, using its display name. The settings
     * are not included.
     *
     * @return summary for display
     */
    public String getDisplayText() {
        StringBuilder tmpBuilder = new StringBuilder(400);
        tmpBuilder.append(String.format(Message.get("MoleculeFilterSummary.DisplayText.read"), this.readCount))
                .append("\n")
                .append(String.format(Message.get("MoleculeFilterSummary.DisplayText.kept"), this.keptCount))
                .append("\n")
                .append(String.format(Message.get("MoleculeFilterSummary.DisplayText.reduced"),
                        this.reducedToLargestComponentCount))
                .append("\n")
                .append(String.format(Message.get("MoleculeFilterSummary.DisplayText.removed"),
                        this.getTotalRemovedCount()));
        for (Map.Entry<MoleculeFilterReason, Integer> tmpEntry : this.removedCountsMap.entrySet()) {
            if (tmpEntry.getValue() > 0) {
                tmpBuilder.append("\n  - ").append(tmpEntry.getKey().getDisplayName()).append(": ")
                        .append(tmpEntry.getValue());
            }
        }
        return tmpBuilder.toString();
    }

    /**
     * Returns an English one-line description of the summary for the log file, including the settings.
     *
     * @return summary for logging
     */
    @Override
    public String toString() {
        StringBuilder tmpBuilder = new StringBuilder(400);
        tmpBuilder.append("Molecule filters (disconnected structures: ")
                .append(this.settingsSnapshot.disconnectedStructuresOption().name())
                .append(", element filter: ").append(this.settingsSnapshot.isElementFilterActive() ? "on" : "off")
                .append(", heavy atom count filter: ");
        if (this.settingsSnapshot.isHeavyAtomCountFilterActive()) {
            tmpBuilder.append("on [").append(this.settingsSnapshot.minimumHeavyAtomCount()).append(", ")
                    .append(this.settingsSnapshot.maximumHeavyAtomCount()).append("]");
        } else {
            tmpBuilder.append("off");
        }
        tmpBuilder.append(", duplicate filter: ").append(this.settingsSnapshot.isDuplicateFilterActive() ? "on" : "off")
                .append("): ").append(this.readCount).append(" molecules read, ")
                .append(this.keptCount).append(" kept, ")
                .append(this.getTotalRemovedCount()).append(" removed (");
        boolean tmpIsFirst = true;
        for (Map.Entry<MoleculeFilterReason, Integer> tmpEntry : this.removedCountsMap.entrySet()) {
            if (!tmpIsFirst) {
                tmpBuilder.append(", ");
            }
            tmpBuilder.append(tmpEntry.getKey().name()).append(": ").append(tmpEntry.getValue());
            tmpIsFirst = false;
        }
        tmpBuilder.append("), ").append(this.reducedToLargestComponentCount).append(" reduced to their largest component");
        int tmpPendingCount = this.getPendingCount();
        if (tmpPendingCount != 0) {
            tmpBuilder.append(", ").append(tmpPendingCount).append(" without final outcome (run incomplete)");
        }
        tmpBuilder.append(".");
        return tmpBuilder.toString();
    }
    //</editor-fold>
}
