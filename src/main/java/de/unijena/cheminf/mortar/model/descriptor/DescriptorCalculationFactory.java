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
 * Recipe for the calculation of one {@link Descriptor}: creates its {@link DescriptorCalculation} for the current
 * configuration of the {@link DescriptorCalculator}. Every descriptor enum constant carries its factory, so a new
 * descriptor only has to be added in one place.
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
interface DescriptorCalculationFactory {

    /**
     * Creates the calculation of the descriptor.
     *
     * @param fingerprintPoolSize     number of fingerprinter instances in the pool of a fingerprint calculation
     *                                (ignored for molecular descriptors)
     * @param circularFingerprintSize size of circular fingerprints (ignored for all other descriptors)
     * @return a new calculation
     */
    DescriptorCalculation create(int fingerprintPoolSize, int circularFingerprintSize);
}
