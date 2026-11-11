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

import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Objects;

/**
 * Immutable outcome of the structural molecule filters for one molecule, see
 * {@link MoleculeFilterRun#applyStructuralFilters(IAtomContainer)}. Either the molecule is kept, then
 * {@link #getMolecule()} returns the container to continue with (the input container or, if it was reduced to its
 * largest component, a new container), or it is rejected, then {@link #getRejectionReason()} returns the reason.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public final class MoleculeFilterResult {
    //<editor-fold desc="Private final variables">
    /**
     * The molecule to continue with, null if the molecule was rejected.
     */
    private final IAtomContainer molecule;

    /**
     * The reason for the rejection, null if the molecule was kept.
     */
    private final MoleculeFilterReason rejectionReason;

    /**
     * Whether the molecule was reduced to its largest component (also if it was rejected by a later filter).
     */
    private final boolean isReducedToLargestComponent;
    //</editor-fold>
    //
    //<editor-fold desc="Private constructor">
    /**
     * Private constructor, use the static factory methods.
     *
     * @param aMolecule the molecule to continue with or null
     * @param aRejectionReason the rejection reason or null
     * @param isReduced whether the molecule was reduced to its largest component
     */
    private MoleculeFilterResult(IAtomContainer aMolecule, MoleculeFilterReason aRejectionReason, boolean isReduced) {
        this.molecule = aMolecule;
        this.rejectionReason = aRejectionReason;
        this.isReducedToLargestComponent = isReduced;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Creates the result for a molecule that passed the structural filters.
     *
     * @param aMolecule the molecule to continue with
     * @param isReduced whether the given molecule is the largest component of the original input
     * @return result of a kept molecule
     * @throws NullPointerException if the given molecule is null
     */
    public static MoleculeFilterResult createKeptResult(IAtomContainer aMolecule, boolean isReduced) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        return new MoleculeFilterResult(aMolecule, null, isReduced);
    }

    /**
     * Creates the result for a molecule that was rejected by a filter.
     *
     * @param aReason the rejection reason
     * @param isReduced whether the molecule had been reduced to its largest component before it was rejected
     * @return result of a rejected molecule
     * @throws NullPointerException if the given reason is null
     */
    public static MoleculeFilterResult createRejectedResult(MoleculeFilterReason aReason, boolean isReduced) throws NullPointerException {
        Objects.requireNonNull(aReason, "Given rejection reason is null.");
        return new MoleculeFilterResult(null, aReason, isReduced);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns whether the molecule passed the structural filters.
     *
     * @return true if the molecule is kept
     */
    public boolean isKept() {
        return Objects.isNull(this.rejectionReason);
    }

    /**
     * Returns the molecule to continue with, i.e. to generate the SMILES string and the data model of. This is the
     * input container itself or, if it was reduced to its largest component, a new container. Null if the molecule
     * was rejected.
     *
     * @return the kept molecule or null
     */
    public IAtomContainer getMolecule() {
        return this.molecule;
    }

    /**
     * Returns the reason why the molecule was rejected, null if it was kept.
     *
     * @return rejection reason or null
     */
    public MoleculeFilterReason getRejectionReason() {
        return this.rejectionReason;
    }

    /**
     * Returns whether the molecule was reduced to its largest component. This is also true if the reduced molecule was
     * rejected by a later filter.
     *
     * @return true if the molecule was reduced
     */
    public boolean isReducedToLargestComponent() {
        return this.isReducedToLargestComponent;
    }
    //</editor-fold>
}
