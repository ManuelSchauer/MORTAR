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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.EnumMap;

/**
 * Test class for the {@link Descriptor} enum, i.e. the descriptor metadata and the static query methods of
 * the descriptor catalog. The calculation of descriptor values is tested in {@link DescriptorCalculatorTest}.
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorTest {
    /**
     * Test method for getDescriptorAndComponentIndex and getDescriptorAndComponentInfo methods.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_getDescriptorAndComponentIndexAndInfo() throws Exception {
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.BCUT, Descriptor.WIENER_NUMBER};
        // MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)

        int[] resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 0);
        Assertions.assertArrayEquals(new int[]{0, 0}, resultIndex);
        String[] resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 0);
        Assertions.assertArrayEquals(new String[]{Descriptor.MOLECULAR_WEIGHT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 1);
        Assertions.assertArrayEquals(new int[]{1, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 1);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 6);
        Assertions.assertArrayEquals(new int[]{1, 5}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 6);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "5"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 7);
        Assertions.assertArrayEquals(new int[]{2, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 7);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 8);
        Assertions.assertArrayEquals(new int[]{2, 1}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 8);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "1"}, resultInfo);

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(descriptors, -1);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(descriptors, -1);
        });

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(descriptors, 9);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(descriptors, 9);
        });
    }

    @Test
    void test_GetAllFingerprints_Excludes3D() {
        Descriptor[] fingerprints = Descriptor.getAllFingerprints();
        Assertions.assertTrue(fingerprints.length > 0);
        for (Descriptor fp : fingerprints) {
            Assertions.assertTrue(fp.isFingerprint());
            Assertions.assertFalse(fp.requires3DCoordinates());
        }
    }

    @Test
    void test_GetSpecifiedDescriptors_3DInclusion() {
        Descriptor[] without3D = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        for (Descriptor d : without3D) {
            Assertions.assertFalse(d.requires3DCoordinates());
        }

        Descriptor[] with3D = Descriptor.getSpecifiedDescriptors(true, true, true, true);
        Assertions.assertEquals(Descriptor.values().length, with3D.length);

        Descriptor[] defaultSelection = Descriptor.getSpecifiedDescriptors(false, false, false, false);
        for (Descriptor d : defaultSelection) {
            Assertions.assertTrue(d.isFast());
            Assertions.assertTrue(d.isSafe());
            Assertions.assertFalse(d.isFingerprint());
            Assertions.assertFalse(d.requires3DCoordinates());
        }
    }

    @Test
    void test_Requires3DCoordinates_Getter() {
        for (Descriptor descriptor : Descriptor.values()) {
            boolean is3D = descriptor == Descriptor.CPSA
                    || descriptor == Descriptor.GRAVITATIONAL_INDEX
                    || descriptor == Descriptor.MOMENT_OF_INERTIA
                    || descriptor == Descriptor.WHIM
                    || descriptor == Descriptor.LENGTH_OVER_BREADTH
                    || descriptor == Descriptor.PETITJEAN_SHAPE_INDEX;
            Assertions.assertEquals(is3D, descriptor.requires3DCoordinates(), "Mismatch for " + descriptor);
        }
    }

    /**
     * Tests the initialization of all enum constants, so that an incorrectly implemented new enum constant is
     * detected immediately:
     * <ul>
     *     <li>every descriptor has a non-empty name, at least one component and a calculation recipe (factory);</li>
     *     <li>the calculation of every descriptor can be created from its recipe in the same way as
     *         DescriptorCalculator.createCalculations() does it, i.e. the CDK descriptor or fingerprinter can be
     *         instantiated and configured (e.g. valid parameters) without an exception, both for bit and for count
     *         circular fingerprints;</li>
     *     <li>the number of fingerprints and their sizes are as expected.</li>
     * </ul>
     */
    @Test
    void test_Initialization() {
        Assertions.assertEquals(Descriptor.values().length, Descriptor.getAllDescriptors().length);
        int numberOfFingerprints = 0;
        for (Descriptor descriptor : Descriptor.values()) {
            Assertions.assertNotNull(descriptor.getName(), descriptor.name());
            Assertions.assertFalse(descriptor.getName().trim().isEmpty(), descriptor.name());
            Assertions.assertTrue(descriptor.getDescriptorComponentNumber() > 0, descriptor.name());
            Assertions.assertNotNull(descriptor.getCalculationFactory(), descriptor.name());
            if (descriptor.isFingerprint()) {
                numberOfFingerprints++;
            }
        }

        // create all calculations in the same way as DescriptorCalculator.createCalculations(); the recipes are called
        // directly instead of creating a DescriptorCalculator (which creates all calculations, too), so that the
        // message names the incorrect enum constant (both for bit and for count circular fingerprints)
        for (boolean isCircularFingerprintCounts : new boolean[]{false, true}) {
            CalculationSettings settings = new CalculationSettings(1, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, isCircularFingerprintCounts);
            EnumMap<Descriptor, DescriptorCalculation> calculations = new EnumMap<>(Descriptor.class);
            for (Descriptor descriptor : Descriptor.values()) {
                DescriptorCalculation calculation = Assertions.assertDoesNotThrow(
                        () -> descriptor.getCalculationFactory().create(settings),
                        "Calculation of " + descriptor.name() + " cannot be created from its recipe (counts="
                                + isCircularFingerprintCounts + ")");
                Assertions.assertNotNull(calculation, descriptor.name());
                calculations.put(descriptor, calculation);
            }
            Assertions.assertEquals(Descriptor.values().length, calculations.size());
        }

        Assertions.assertEquals(10, numberOfFingerprints);
        Assertions.assertEquals(881, Descriptor.PUBCHEM_FINGERPRINTER.getDescriptorComponentNumber());
        Assertions.assertEquals(166, Descriptor.MACCS_FINGERPRINTER.getDescriptorComponentNumber());
        for (Descriptor descriptor : Descriptor.values()) {
            if (descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_")) {
                Assertions.assertEquals(Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE,
                        descriptor.getDescriptorComponentNumber(), descriptor.name());
            }
        }
    }
}
