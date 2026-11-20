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

package de.unijena.cheminf.mortar.model.util;

import de.unijena.cheminf.mortar.model.data.FragmentDataModel;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Test class for CollectionUtil routines.
 */
public class CollectionUtilTest {
    /**
     * Tests whether a few examples of number of elements and load factor pairs generate suitable initial capacities that
     * lead to rehash thresholds higher than the number of elements using the specific calculation method from CollectionUtil.
     */
    @Test
    public void calculateInitialHashMapSizeTest() {
        int[] tmpNumberOfElements = new int[]{10, 100, 4353456, 30};
        float[] tmpLoadFactor = new float[]{0.75f, 0.75f, 0.6f, 0.75f};
        int[] tmpExpectedInitialCapacity = new int[]{15, 135, 7255762, 42};
        for (int i = 0; i < tmpNumberOfElements.length; i++) {
            int tmpCalculatedInitialHashMapCapacity = CollectionUtil.calculateInitialHashCollectionCapacity(tmpNumberOfElements[i], tmpLoadFactor[i]);
            float tmpRehashThreshold = tmpCalculatedInitialHashMapCapacity * tmpLoadFactor[i];
            Assertions.assertTrue(tmpNumberOfElements[i] < tmpRehashThreshold);
            Assertions.assertEquals(tmpExpectedInitialCapacity[i], tmpCalculatedInitialHashMapCapacity);
        }
    }
    //
    /**
     * Tests that fragments with the same display label are sorted below each other, groups by summed frequency and
     * fragments within a group by frequency, and that fragments without display label are listed last.
     */
    @Test
    public void sortGivenFragmentListByDisplayLabelGroupsTest() {
        List<FragmentDataModel> tmpList = new ArrayList<>();
        //label, SMILES, frequency
        Object[][] tmpData = {{"Ethyl", "CC", 5}, {"Acetyl", "CCO", 4}, {null, "N", 9}, {"Ethyl", "C=C", 7},
                {"Acetyl", "CC=O", 3}, {"Methyl", "C", 1}};
        for (Object[] tmpEntry : tmpData) {
            FragmentDataModel tmpFragment = new FragmentDataModel((String) tmpEntry[1], null, new HashMap<>());
            if (tmpEntry[0] != null) {
                tmpFragment.addDisplayLabel((String) tmpEntry[0]);
            }
            for (int i = 0; i < (int) tmpEntry[2]; i++) {
                tmpFragment.incrementAbsoluteFrequency();
            }
            tmpList.add(tmpFragment);
        }
        CollectionUtil.sortGivenFragmentListByDisplayLabelGroups(tmpList);
        Assertions.assertEquals(List.of("C=C", "CC", "CCO", "CC=O", "C", "N"),
                tmpList.stream().map(FragmentDataModel::getUniqueSmiles).toList());
    }
}
