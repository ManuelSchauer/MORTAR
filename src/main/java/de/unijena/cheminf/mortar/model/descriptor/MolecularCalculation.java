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
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * Calculation of a CDK molecular descriptor. One CDK descriptor instance is shared by all calculation threads
 * (the thread-safety of the used CDK descriptors has been checked).
 * <p>
 * The calculation handles four CDK result types:
 * <ul>
 *   <li>{@link DoubleResult} - Single double value (e.g., molecular weight)</li>
 *   <li>{@link IntegerResult} - Single integer value (e.g., atom count)</li>
 *   <li>{@link DoubleArrayResult} - Array of doubles (e.g., BCUT returns 6 values)</li>
 *   <li>{@link IntegerArrayResult} - Array of integers (e.g., amino acid counts)</li>
 * </ul>
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
final class MolecularCalculation implements DescriptorCalculation {

    /**
     * Logger of this class.
     */
    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(MolecularCalculation.class);

    /**
     * The CDK descriptor instance, shared by all calculation threads.
     */
    private final IMolecularDescriptor cdkDescriptor;

    /**
     * Prefix of the component names, empty if the CDK descriptor names are used unchanged.
     */
    private final String namePrefix;

    /**
     * Constructor.
     *
     * @param cdkDescriptor the (configured) CDK descriptor instance
     */
    MolecularCalculation(IMolecularDescriptor cdkDescriptor) {
        this(cdkDescriptor, "");
    }

    /**
     * Constructor with a prefix for the component names, used where the CDK descriptor names are not unique among all
     * descriptors (e.g. "AA." for the amino acid counts, whose CDK names like "nC" collide with the atom counts).
     *
     * @param cdkDescriptor the (configured) CDK descriptor instance
     * @param namePrefix    prefix of the component names (may be empty)
     */
    MolecularCalculation(IMolecularDescriptor cdkDescriptor, String namePrefix) {
        this.cdkDescriptor = cdkDescriptor;
        this.namePrefix = namePrefix;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns the descriptor names of the configured CDK descriptor instance (e.g. "WPATH" and "WPOL" for the Wiener
     * numbers, or "nC" for an atom count descriptor configured to count carbon atoms), with the name prefix if one
     * is set (e.g. "AA.nC" for the cysteine count of the amino acid count descriptor).
     */
    @Override
    public String[] getComponentNames() {
        String[] names = this.cdkDescriptor.getDescriptorNames().clone();
        for (int i = 0; i < names.length; i++) {
            names[i] = this.namePrefix + names[i];
        }
        return names;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean calculate(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex) {
        try {
            IDescriptorResult result = this.cdkDescriptor.calculate(atomContainer).getValue();
            //note: we ignore all the Sonar suggestions here to keep the code compatible with Java 8
            if (result instanceof DoubleResult) {
                DoubleResult doubleResult = (DoubleResult) result;
                vector[startIndex] = (float) doubleResult.doubleValue();
            } else if (result instanceof IntegerResult) {
                IntegerResult integerResult = (IntegerResult) result;
                vector[startIndex] = (float) integerResult.intValue();
            } else if (result instanceof DoubleArrayResult) {
                DoubleArrayResult arrayResult = (DoubleArrayResult) result;
                for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                    vector[startIndex + i] = (float) arrayResult.get(i);
                }
            } else if (result instanceof IntegerArrayResult) {
                IntegerArrayResult arrayResult = (IntegerArrayResult) result;
                for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                    vector[startIndex + i] = (float) arrayResult.get(i);
                }
            } else {
                LOGGER.warn("Unexpected result type ", result.getClass().getName(), " for descriptor ", descriptor.name(), ". Expected DoubleResult, IntegerResult, DoubleArrayResult or IntegerArrayResult.");
                return false;
            }
            return true;
        } catch (Exception exception) {
            LOGGER.error("Descriptor ", descriptor.name(), " failed to calculate: ", exception.getMessage());
            return false;
        }
    }
}
