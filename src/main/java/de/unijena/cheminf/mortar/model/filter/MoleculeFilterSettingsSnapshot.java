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

import java.util.Objects;

/**
 * Immutable snapshot of the molecule filter settings, taken at the start of an import run (see
 * {@link MoleculeFilterSettings#createSnapshot()}) so that a settings change during a running import cannot produce a
 * mixed result. All filters described by a snapshot are applied by a {@link MoleculeFilterRun}.
 *
 * @param disconnectedStructuresOption how molecules with more than one unconnected component are treated
 * @param isElementFilterActive whether only molecules consisting of the allowed elements
 *                              ({@link MoleculeFilterUtil#ALLOWED_ATOMIC_NUMBERS}) are taken over
 * @param isHeavyAtomCountFilterActive whether only molecules with a heavy atom count within the inclusive range
 *                                     [minimumHeavyAtomCount, maximumHeavyAtomCount] are taken over
 * @param minimumHeavyAtomCount inclusive lower bound of the heavy atom count range, not negative
 * @param maximumHeavyAtomCount inclusive upper bound of the heavy atom count range, not smaller than the lower bound
 * @param isDuplicateFilterActive whether molecules with a structure (duplicate key, see
 *                                {@link MoleculeFilterUtil#createDuplicateKey(org.openscience.cdk.interfaces.IAtomContainer,
 *                                boolean, String, boolean)}) that occurred before in the same import run are removed
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record MoleculeFilterSettingsSnapshot(
        MoleculeFilterSettings.DisconnectedStructuresOption disconnectedStructuresOption,
        boolean isElementFilterActive,
        boolean isHeavyAtomCountFilterActive,
        int minimumHeavyAtomCount,
        int maximumHeavyAtomCount,
        boolean isDuplicateFilterActive
) {
    //<editor-fold desc="Constructor">
    /**
     * Compact canonical constructor that checks the given values.
     *
     * @throws NullPointerException if the disconnected structures option is null
     * @throws IllegalArgumentException if the minimum heavy atom count is negative or larger than the maximum heavy
     * atom count
     */
    public MoleculeFilterSettingsSnapshot {
        Objects.requireNonNull(disconnectedStructuresOption, "Given disconnected structures option is null.");
        if (minimumHeavyAtomCount < 0) {
            throw new IllegalArgumentException("The minimum heavy atom count must not be negative but is "
                    + minimumHeavyAtomCount + ".");
        }
        if (minimumHeavyAtomCount > maximumHeavyAtomCount) {
            throw new IllegalArgumentException("The minimum heavy atom count " + minimumHeavyAtomCount
                    + " is larger than the maximum heavy atom count " + maximumHeavyAtomCount + ".");
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Creates a snapshot holding the default values of all settings as declared in the public constants of
     * {@link MoleculeFilterSettings}. With these, all filters are inactive and an import behaves exactly as without
     * molecule filters.
     *
     * @return snapshot of the default settings
     */
    public static MoleculeFilterSettingsSnapshot createDefault() {
        return new MoleculeFilterSettingsSnapshot(
                MoleculeFilterSettings.DISCONNECTED_STRUCTURES_SETTING_DEFAULT,
                MoleculeFilterSettings.ELEMENT_FILTER_SETTING_DEFAULT,
                MoleculeFilterSettings.HEAVY_ATOM_COUNT_FILTER_SETTING_DEFAULT,
                MoleculeFilterSettings.MINIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                MoleculeFilterSettings.MAXIMUM_HEAVY_ATOM_COUNT_SETTING_DEFAULT,
                MoleculeFilterSettings.DUPLICATE_FILTER_SETTING_DEFAULT);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Returns whether at least one of the structural filters, i.e. the filters that inspect the atom container
     * (disconnected structures, elements, heavy atom count), is active. If none is active, the structural filtering of
     * a molecule leaves it untouched.
     *
     * @return true if a structural filter is active
     */
    public boolean isAnyStructuralFilterActive() {
        return this.disconnectedStructuresOption != MoleculeFilterSettings.DisconnectedStructuresOption.KEEP_UNCHANGED
                || this.isElementFilterActive
                || this.isHeavyAtomCountFilterActive;
    }

    /**
     * Returns whether at least one filter is active, i.e. whether the import result can differ from an import without
     * molecule filters.
     *
     * @return true if any filter is active
     */
    public boolean isAnyFilterActive() {
        return this.isAnyStructuralFilterActive() || this.isDuplicateFilterActive;
    }
    //</editor-fold>
}
