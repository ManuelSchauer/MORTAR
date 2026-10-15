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

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Calculates the components of one {@link Descriptor} for one molecule and writes them into a data vector.
 * Instances are created by the {@link DescriptorCalculationFactory} of the descriptor and are shared by all
 * calculation threads of the {@link DescriptorCalculator}, so implementations must be thread-safe.
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
interface DescriptorCalculation {

    /**
     * Calculates the components of the given descriptor for the given molecule and writes them into the vector.
     * Note: No input validation is performed here, this is done by the public methods of the
     * {@link DescriptorCalculator}.
     *
     * @param descriptor    the descriptor to calculate; provides the number of components and the name for logging
     * @param atomContainer the molecule (IS NOT CHANGED); aromaticity must already be perceived if required
     * @param vector        the data vector (row of the data matrix) to write the components to (MAY BE CHANGED);
     *                      the components are written from startIndex on
     * @param startIndex    the index in vector of the first component
     * @return true if all components were written to the vector, false if the calculation failed
     * @throws InterruptedException if the current thread is interrupted while waiting for a pooled instance
     */
    boolean calculate(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex)
            throws InterruptedException;
}
