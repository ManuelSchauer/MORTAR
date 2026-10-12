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
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.IntFunction;

/**
 * Calculation of a CDK fingerprint. Fingerprinter instances are not thread-safe, so this calculation owns a pool of
 * fingerprinter instances: every calculation takes one instance from the pool and returns it to the same pool
 * afterwards.
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
     * Constructor, creates and fills the pool of fingerprinter instances.
     *
     * @param fingerprinterFactory    creates a fingerprinter instance; the argument is the circular fingerprint size
     * @param poolSize                number of fingerprinter instances in the pool (must be greater than 0)
     * @param circularFingerprintSize size of circular fingerprints, passed to the fingerprinter factory
     */
    FingerprintCalculation(IntFunction<IFingerprinter> fingerprinterFactory, int poolSize, int circularFingerprintSize) {
        this.pool = new LinkedBlockingQueue<>(poolSize);
        for (int i = 0; i < poolSize; i++) {
            this.pool.add(fingerprinterFactory.apply(circularFingerprintSize));
        }
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
            IBitFingerprint bitFingerprint = fingerprinter.getBitFingerprint(atomContainer);
            for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                vector[startIndex + i] = bitFingerprint.get(i) ? 1.0f : 0.0f;
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
