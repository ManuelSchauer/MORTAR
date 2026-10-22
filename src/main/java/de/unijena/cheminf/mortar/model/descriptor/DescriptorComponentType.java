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
 * Type of the values of a descriptor component, i.e. of a column of a data matrix calculated by the
 * {@link DescriptorCalculator}. The type tells how the values of a column are to be interpreted, e.g. for choosing a
 * suitable preprocessing before machine learning: binary values are usually used as they are, counts are often
 * transformed (e.g. with log(1 + x)) and continuous values are usually standardised. The type does not prescribe a
 * preprocessing; this depends on the model.
 * <p>
 * The type of every descriptor is defined in the {@link Descriptor} enum, see {@link Descriptor#getComponentType()}.
 * Note: The type of the circular fingerprints depends on the configuration of the calculator, see
 * {@link DescriptorCalculator#getComponentType(Descriptor)}.
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
public enum DescriptorComponentType {
    /**
     * Binary values, 0 or 1, e.g. the bits of a fingerprint.
     */
    BINARY,
    /**
     * Non-negative integer counts, e.g. an atom count or a component of a folded count fingerprint.
     */
    COUNT,
    /**
     * Real values, e.g. the molecular weight or a logP value. Integer-valued indices that are not counts (e.g. the
     * Zagreb index) are continuous, too.
     */
    CONTINUOUS
}
