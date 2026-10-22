/*
 * Copyright (C) 2026 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package de.unijena.cheminf.mortar.model.descriptor;

/**
 * Immutable snapshot of the configuration of the {@link DescriptorCalculator} that is needed to create the
 * {@link DescriptorCalculation} of a descriptor (see {@link DescriptorCalculationFactory}). New settings can be added
 * here without changing the signature of the factory.
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
final class CalculationSettings {

    /**
     * Number of fingerprinter instances in the pool of a fingerprint calculation.
     */
    private final int fingerprintPoolSize;

    /**
     * Size (number of components) of circular fingerprints.
     */
    private final int circularFingerprintSize;

    /**
     * True: circular fingerprints are calculated as folded count fingerprints; false: as bit fingerprints.
     */
    private final boolean isCircularFingerprintCounts;

    /**
     * Constructor.
     *
     * @param fingerprintPoolSize         number of fingerprinter instances in the pool of a fingerprint calculation
     * @param circularFingerprintSize     size (number of components) of circular fingerprints
     * @param isCircularFingerprintCounts true: circular fingerprints are calculated as folded count fingerprints;
     *                                    false: as bit fingerprints
     */
    CalculationSettings(int fingerprintPoolSize, int circularFingerprintSize, boolean isCircularFingerprintCounts) {
        this.fingerprintPoolSize = fingerprintPoolSize;
        this.circularFingerprintSize = circularFingerprintSize;
        this.isCircularFingerprintCounts = isCircularFingerprintCounts;
    }

    /**
     * Returns the number of fingerprinter instances in the pool of a fingerprint calculation.
     *
     * @return the fingerprint pool size
     */
    int getFingerprintPoolSize() {
        return this.fingerprintPoolSize;
    }

    /**
     * Returns the size (number of components) of circular fingerprints.
     *
     * @return the circular fingerprint size
     */
    int getCircularFingerprintSize() {
        return this.circularFingerprintSize;
    }

    /**
     * Returns whether circular fingerprints are calculated as folded count fingerprints.
     *
     * @return true: folded count fingerprints; false: bit fingerprints
     */
    boolean isCircularFingerprintCounts() {
        return this.isCircularFingerprintCounts;
    }
}
