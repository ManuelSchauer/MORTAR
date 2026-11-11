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
import de.unijena.cheminf.mortar.model.util.IDisplayEnum;

/**
 * Reasons for which a molecule is not taken over into the data set during an import with molecule filters, see
 * {@link MoleculeFilterRun}. The order of the constants is the order in which the reasons should be listed in a
 * summary; it does not equal the order in which the filters are applied, which is: empty molecule check, disconnected
 * structures, elements, heavy atom count, (SMILES generation by the importer), duplicates.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public enum MoleculeFilterReason implements IDisplayEnum {
    /**
     * The molecule contains no atoms. Only checked if at least one structural filter (disconnected structures,
     * elements, heavy atom count) is active; otherwise, empty molecules are passed on unchanged as before (the importer
     * then fails to generate a SMILES string for them).
     */
    EMPTY_MOLECULE(
            Message.get("MoleculeFilterReason.EmptyMolecule.displayName"),
            Message.get("MoleculeFilterReason.EmptyMolecule.tooltip")),
    /**
     * The molecule consists of more than one unconnected component and the disconnected structures setting is set to
     * {@link MoleculeFilterSettings.DisconnectedStructuresOption#REMOVE_MOLECULE}.
     */
    DISCONNECTED_STRUCTURE(
            Message.get("MoleculeFilterReason.DisconnectedStructure.displayName"),
            Message.get("MoleculeFilterReason.DisconnectedStructure.tooltip")),
    /**
     * The molecule contains an element that is not in {@link MoleculeFilterUtil#ALLOWED_ATOMIC_NUMBERS} or a pseudo
     * atom.
     */
    DISALLOWED_ELEMENT(
            Message.get("MoleculeFilterReason.DisallowedElement.displayName"),
            Message.get("MoleculeFilterReason.DisallowedElement.tooltip")),
    /**
     * The heavy atom count of the molecule lies outside the set inclusive range.
     */
    HEAVY_ATOM_COUNT_OUT_OF_RANGE(
            Message.get("MoleculeFilterReason.HeavyAtomCountOutOfRange.displayName"),
            Message.get("MoleculeFilterReason.HeavyAtomCountOutOfRange.tooltip")),
    /**
     * A molecule with the same unique SMILES string was taken over before in the same import run.
     */
    DUPLICATE(
            Message.get("MoleculeFilterReason.Duplicate.displayName"),
            Message.get("MoleculeFilterReason.Duplicate.tooltip")),
    /**
     * An active structural filter threw an exception for the molecule, so it could not be decided whether the molecule
     * meets the filter criteria. Such molecules are removed (fail-closed) and the exception is logged.
     */
    FILTER_CHECK_FAILED(
            Message.get("MoleculeFilterReason.FilterCheckFailed.displayName"),
            Message.get("MoleculeFilterReason.FilterCheckFailed.tooltip")),
    /**
     * The importer could not generate a SMILES string for the molecule; this is independent of the filters and was
     * also the case before they were introduced, it is only counted here to make the summary complete.
     */
    SMILES_GENERATION_FAILED(
            Message.get("MoleculeFilterReason.SmilesGenerationFailed.displayName"),
            Message.get("MoleculeFilterReason.SmilesGenerationFailed.tooltip"));
    //
    /**
     * Language-specific name for display in GUI.
     */
    private final String displayName;
    /**
     * Language-specific tooltip text for display in GUI.
     */
    private final String tooltip;
    //
    /**
     * Constructor.
     *
     * @param aDisplayName display name
     * @param aTooltip tooltip text
     */
    private MoleculeFilterReason(String aDisplayName, String aTooltip) {
        this.displayName = aDisplayName;
        this.tooltip = aTooltip;
    }
    //
    @Override
    public String getDisplayName() {
        return this.displayName;
    }
    //
    @Override
    public String getTooltipText() {
        return this.tooltip;
    }
}
