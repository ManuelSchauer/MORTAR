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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Test class for the column names and component types of the {@link DescriptorCalculator}, taken from the test class
 * of the CDK descriptor calculation engine (the tests of the descriptor values need test resources of the CDK and stay
 * there).
 */
class DescriptorCalculatorTest {

    /**
     * Calculator with the default configuration.
     */
    private static final DescriptorCalculator CALCULATOR = new DescriptorCalculator();

    /**
     * Tests getColumnNames(): for every descriptor, the number of names is equal to its number of components and no
     * name is empty; the names of all descriptors are unique (the amino acid counts get the prefix "AA." because their
     * CDK names like "nC" collide with the atom counts); for several descriptors, the names are returned in the order
     * of the matrix columns (example with the CDK descriptor names "MW", "WPATH", "WPOL" and "nC"); fingerprints get
     * generated names (prefix and bit index); empty and null input.
     */
    @Test
    void test_getColumnNames() {
        for (Descriptor descriptor : Descriptor.values()) {
            Descriptor[] singleDescriptor = new Descriptor[]{descriptor};
            String[] names = CALCULATOR.getColumnNames(singleDescriptor);
            Assertions.assertEquals(Descriptor.getNumberOfComponents(singleDescriptor), names.length, descriptor.name());
            for (String name : names) {
                Assertions.assertNotNull(name, descriptor.name());
                Assertions.assertFalse(name.trim().isEmpty(), descriptor.name());
            }
        }
        Descriptor[] allDescriptors = Descriptor.values();
        String[] allNames = CALCULATOR.getColumnNames(allDescriptors);
        Assertions.assertEquals(Descriptor.getNumberOfComponents(allDescriptors), allNames.length);
        Set<String> uniqueNames = new HashSet<>();
        for (String name : allNames) {
            Assertions.assertTrue(uniqueNames.add(name), "Duplicate column name: " + name);
        }

        // amino acid counts: prefixed CDK names, atom counts unchanged
        String[] aminoAcidNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.AMINO_ACID_COUNT});
        Assertions.assertEquals("AA.nA", aminoAcidNames[0]);
        Assertions.assertTrue(Arrays.asList(aminoAcidNames).contains("AA.nC"));
        Assertions.assertArrayEquals(new String[]{"nC"}, CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.ATOM_COUNT_C}));

        // names in the order of the matrix columns
        Assertions.assertArrayEquals(new String[]{"MW", "WPATH", "WPOL", "nC"}, CALCULATOR.getColumnNames(
                new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.WIENER_NUMBER, Descriptor.ATOM_COUNT_C}));
        Assertions.assertArrayEquals(new String[]{"nC", "MW"}, CALCULATOR.getColumnNames(
                new Descriptor[]{Descriptor.ATOM_COUNT_C, Descriptor.MOLECULAR_WEIGHT}));

        // fingerprints: generated names
        String[] pubChemNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER});
        Assertions.assertEquals("PubChem_0", pubChemNames[0]);
        Assertions.assertEquals("PubChem_880", pubChemNames[880]);
        String[] maccsNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.MACCS_FINGERPRINTER});
        Assertions.assertEquals("MACCS_0", maccsNames[0]);
        Assertions.assertEquals("MACCS_165", maccsNames[165]);
        int circularSize = CALCULATOR.getCircularFingerprintSize();
        String[] ecfp4Names = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4});
        Assertions.assertEquals(circularSize, ecfp4Names.length);
        Assertions.assertEquals("ECFP4_0", ecfp4Names[0]);
        Assertions.assertEquals("ECFP4_" + (circularSize - 1), ecfp4Names[circularSize - 1]);
        Assertions.assertEquals("FCFP6_0", CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6})[0]);

        // empty and null input
        Assertions.assertEquals(0, CALCULATOR.getColumnNames(new Descriptor[0]).length);
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getColumnNames(null));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, null}));
    }

    /**
     * Tests getComponentType() and getComponentTypes(): every descriptor has a type, all fingerprints are binary by
     * default (also the circular fingerprints), and only the circular fingerprints become counts in count mode; the
     * types are returned per column in the order of the columns, i.e. with the same length as the column names, also
     * for a non-default circular fingerprint size; empty and null input.
     */
    @Test
    void test_getComponentTypes() {
        for (Descriptor descriptor : Descriptor.values()) {
            Assertions.assertNotNull(descriptor.getComponentType(), descriptor.name());
            Assertions.assertEquals(descriptor.getComponentType(), CALCULATOR.getComponentType(descriptor), descriptor.name());
            Assertions.assertEquals(descriptor.isFingerprint(), descriptor.getComponentType() == DescriptorComponentType.BINARY, descriptor.name());
        }
        Assertions.assertEquals(DescriptorComponentType.CONTINUOUS, Descriptor.MOLECULAR_WEIGHT.getComponentType());
        Assertions.assertEquals(DescriptorComponentType.COUNT, Descriptor.ATOM_COUNT_C.getComponentType());
        Assertions.assertEquals(DescriptorComponentType.BINARY, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4.getComponentType());

        // one type per column, in the order of the columns
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT_C, Descriptor.MACCS_FINGERPRINTER, Descriptor.WIENER_NUMBER};
        DescriptorComponentType[] types = CALCULATOR.getComponentTypes(descriptors);
        Assertions.assertEquals(CALCULATOR.getNumberOfComponents(descriptors), types.length);
        Assertions.assertEquals(CALCULATOR.getColumnNames(descriptors).length, types.length);
        Assertions.assertEquals(DescriptorComponentType.CONTINUOUS, types[0]);
        Assertions.assertEquals(DescriptorComponentType.COUNT, types[1]);
        for (int i = 2; i < 2 + 166; i++) {
            Assertions.assertEquals(DescriptorComponentType.BINARY, types[i], "column " + i);
        }
        Assertions.assertEquals(DescriptorComponentType.CONTINUOUS, types[168]);
        Assertions.assertEquals(DescriptorComponentType.CONTINUOUS, types[169]);
        Assertions.assertEquals(CALCULATOR.getNumberOfComponents(Descriptor.values()), CALCULATOR.getComponentTypes(Descriptor.values()).length);

        // count mode: only the circular fingerprints become counts; a non-default size changes the number of types
        DescriptorCalculator countCalculator = new DescriptorCalculator(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, 2048, true);
        for (Descriptor descriptor : Descriptor.values()) {
            DescriptorComponentType expected = descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_")
                    ? DescriptorComponentType.COUNT : descriptor.getComponentType();
            Assertions.assertEquals(expected, countCalculator.getComponentType(descriptor), descriptor.name());
        }
        Descriptor[] fingerprints = new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4};
        DescriptorComponentType[] countTypes = countCalculator.getComponentTypes(fingerprints);
        Assertions.assertEquals(881 + 2048, countTypes.length);
        Assertions.assertEquals(DescriptorComponentType.BINARY, countTypes[880]);
        Assertions.assertEquals(DescriptorComponentType.COUNT, countTypes[881]);
        Assertions.assertEquals(DescriptorComponentType.COUNT, countTypes[881 + 2047]);

        // empty and null input
        Assertions.assertEquals(0, CALCULATOR.getComponentTypes(new Descriptor[0]).length);
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getComponentType(null));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getComponentTypes(null));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getComponentTypes(new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, null}));
    }
}
