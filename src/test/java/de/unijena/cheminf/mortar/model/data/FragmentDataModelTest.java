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

package de.unijena.cheminf.mortar.model.data;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Test class for the display labels of {@link FragmentDataModel}.
 */
class FragmentDataModelTest {
    /**
     * Tests that the display labels are stored per parent molecule and combined in alphabetical order over all parent
     * molecules, also if the first parent molecule or a later one has no display label for the fragment, and that the
     * frequencies of the labels are only returned for molecules with more than one label, in the given order.
     */
    @Test
    void displayLabelsPerParentMoleculeTest() {
        FragmentDataModel tmpFragment = new FragmentDataModel("O", null, new HashMap<>());
        MoleculeDataModel tmpFirst = new MoleculeDataModel("COC", "first", new HashMap<>());
        MoleculeDataModel tmpSecond = new MoleculeDataModel("CCO", "second", new HashMap<>());
        MoleculeDataModel tmpThird = new MoleculeDataModel("COCCO", "third", new HashMap<>());
        MoleculeDataModel tmpFourth = new MoleculeDataModel("CCOCC", "fourth", new HashMap<>());
        MoleculeDataModel tmpFifth = new MoleculeDataModel("CCCO", "fifth", new HashMap<>());
        Assertions.assertNull(tmpFragment.getDisplayLabelText());
        Assertions.assertNull(tmpFragment.getDisplayLabelText(tmpFirst));
        Assertions.assertNull(tmpFragment.getDisplayLabelFrequencies(tmpFirst));
        tmpFragment.addDisplayLabels(tmpFirst, Map.of("Ether", 1));
        tmpFragment.addDisplayLabels(tmpSecond, Map.of("Hydroxyl", 2));
        Map<String, Integer> tmpThirdLabels = new LinkedHashMap<>();
        tmpThirdLabels.put("Hydroxyl", 1);
        tmpThirdLabels.put("Ether", 3);
        tmpFragment.addDisplayLabels(tmpThird, tmpThirdLabels);
        tmpFragment.addDisplayLabels(tmpFourth, Map.of("Ether", 1));
        tmpFragment.addDisplayLabels(tmpFifth, null);
        Assertions.assertEquals("Ether, Hydroxyl", tmpFragment.getDisplayLabelText());
        Assertions.assertEquals("Ether", tmpFragment.getDisplayLabelText(tmpFirst));
        Assertions.assertEquals("Hydroxyl", tmpFragment.getDisplayLabelText(tmpSecond));
        Assertions.assertEquals("Ether, Hydroxyl", tmpFragment.getDisplayLabelText(tmpThird));
        Assertions.assertEquals("Ether", tmpFragment.getDisplayLabelText(tmpFourth));
        Assertions.assertNull(tmpFragment.getDisplayLabelText(tmpFifth));
        Assertions.assertNull(tmpFragment.getDisplayLabelFrequencies(tmpFirst));
        Assertions.assertNull(tmpFragment.getDisplayLabelFrequencies(tmpSecond));
        Assertions.assertEquals(List.of(Map.entry("Hydroxyl", 1), Map.entry("Ether", 3)),
                List.copyOf(tmpFragment.getDisplayLabelFrequencies(tmpThird).entrySet()));
        Assertions.assertNull(tmpFragment.getDisplayLabelFrequencies(tmpFifth));
        //altering the given map does not alter the stored frequencies
        tmpThirdLabels.put("Epoxide", 1);
        Assertions.assertEquals(2, tmpFragment.getDisplayLabelFrequencies(tmpThird).size());
        //first parent molecule without display label
        FragmentDataModel tmpOtherFragment = new FragmentDataModel("N", null, new HashMap<>());
        tmpOtherFragment.addDisplayLabels(tmpFirst, Map.of());
        tmpOtherFragment.addDisplayLabels(tmpSecond, Map.of("Amine", 1));
        Assertions.assertEquals("Amine", tmpOtherFragment.getDisplayLabelText());
        Assertions.assertNull(tmpOtherFragment.getDisplayLabelText(tmpFirst));
        Assertions.assertEquals("Amine", tmpOtherFragment.getDisplayLabelText(tmpSecond));
        Assertions.assertThrows(NullPointerException.class, () -> tmpOtherFragment.addDisplayLabels(null, Map.of("Amine", 1)));
    }
}
