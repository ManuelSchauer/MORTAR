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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.util.Arrays;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.Supplier;

/**
 * Calculation of a CDK fingerprint. Fingerprinter instances are not thread-safe, so this calculation owns a pool of
 * fingerprinter instances: every calculation takes one instance from the pool and returns it to the same pool
 * afterwards.
 * <p>
 * The fingerprint is written as bit fingerprint (1.0f for a set bit, 0.0f otherwise) or, for circular fingerprints, as
 * folded count fingerprint: the count fingerprint of the CDK is not folded (it enumerates hash codes and their counts),
 * so every hash code is folded into the fingerprint size in exactly the same way as
 * {@code CircularFingerprinter.getBitFingerprint()} does it, and the counts of all hash codes folded into the same
 * position are added up. Therefore, a component of the folded count fingerprint is greater than 0 if and only if the
 * corresponding bit of the bit fingerprint is set.
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 */
final class FingerprintCalculation implements DescriptorCalculation {

    /**
     * Logger of this class.
     */
    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(FingerprintCalculation.class);

    /**
     * Pool of fingerprinter instances owned by this calculation.
     */
    private final BlockingQueue<IFingerprinter> pool;

    /**
     * True: the fingerprint is written as folded count fingerprint; false: as bit fingerprint.
     */
    private final boolean isCounts;

    /**
     * Prefix of the component names, e.g. "PubChem" or "ECFP4".
     */
    private final String namePrefix;

    /**
     * Size of the fingerprint, i.e. number of components, as given by the fingerprinter.
     */
    private final int size;

    /**
     * Constructor, creates and fills the pool of fingerprinter instances.
     *
     * @param fingerprinterFactory creates a (configured) fingerprinter instance
     * @param poolSize             number of fingerprinter instances in the pool (must be greater than 0)
     * @param isCounts             true: the fingerprint is written as folded count fingerprint (only for circular
     *                             fingerprints, see the class documentation); false: as bit fingerprint
     * @param namePrefix           prefix of the component names, e.g. "PubChem" or "ECFP4"
     */
    FingerprintCalculation(Supplier<IFingerprinter> fingerprinterFactory, int poolSize, boolean isCounts, String namePrefix) {
        this.pool = new LinkedBlockingQueue<>(poolSize);
        for (int i = 0; i < poolSize; i++) {
            this.pool.add(fingerprinterFactory.get());
        }
        this.isCounts = isCounts;
        this.namePrefix = namePrefix;
        this.size = this.pool.peek().getSize();
    }

    /**
     * {@inheritDoc}
     * <p>
     * The CDK fingerprinters do not name their bits, so the names are generated as prefix and bit index, e.g.
     * "PubChem_0" ... "PubChem_880" or "ECFP4_0" ... "ECFP4_1023".
     */
    @Override
    public String[] getComponentNames() {
        String[] names = new String[this.size];
        for (int i = 0; i < this.size; i++) {
            names[i] = this.namePrefix + "_" + i;
        }
        return names;
    }

    /**
     * Returns whether the fingerprint is written as folded count fingerprint.
     *
     * @return true: folded count fingerprint; false: bit fingerprint
     */
    boolean isCounts() {
        return this.isCounts;
    }

    /**
     * Folds a hash code of a circular fingerprint into a fingerprint of the given size, in exactly the same way as
     * {@code CircularFingerprinter.getBitFingerprint()} does it (the hash code is interpreted as unsigned 32-bit value).
     *
     * @param hashCode the hash code
     * @param size     the fingerprint size
     * @return the position in the folded fingerprint, between 0 and size - 1
     */
    static int fold(int hashCode, int size) {
        long unsignedHashCode = hashCode >= 0 ? hashCode : ((hashCode & 0x7FFFFFFF) | (1L << 31));
        return (int) (unsignedHashCode % size);
    }

    /**
     * {@inheritDoc}
     *
     * @throws InterruptedException if the current thread is interrupted while waiting for a fingerprinter instance
     *                              from the pool (or is already interrupted when called)
     */
    @Override
    public boolean calculate(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex)
            throws InterruptedException {
        // take() blocks until a fingerprinter is available; InterruptedException propagates to the caller
        IFingerprinter fingerprinter = this.pool.take();
        try {
            if (this.isCounts) {
                ICountFingerprint countFingerprint = fingerprinter.getCountFingerprint(atomContainer);
                Arrays.fill(vector, startIndex, startIndex + this.size, 0.0f);
                for (int i = 0; i < countFingerprint.numOfPopulatedbins(); i++) {
                    vector[startIndex + FingerprintCalculation.fold(countFingerprint.getHash(i), this.size)] += countFingerprint.getCount(i);
                }
            } else {
                IBitFingerprint bitFingerprint = fingerprinter.getBitFingerprint(atomContainer);
                for (int i = 0; i < this.size; i++) {
                    vector[startIndex + i] = bitFingerprint.get(i) ? 1.0f : 0.0f;
                }
            }
            return true;
        } catch (CDKException | RuntimeException exception) {
            LOGGER.warn("Failed to calculate: " + descriptor.name(), exception);
            return false;
        } finally {
            // always back to the pool it was taken from; the pool cannot be full here
            this.pool.offer(fingerprinter);
        }
    }
}
