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

package de.unijena.cheminf.mortar.model.clustering;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Properties;

/**
 * Tests for {@link FeatureFilter}.
 *
 * @author Manuel Schauer
 */
class FeatureFilterTest {
    /**
     * Creates a matrix with the given columns (column-wise input), all of type CONTINUOUS, with the given group names.
     *
     * @param aGroupNames group names of the columns
     * @param aColumns column values
     * @return feature matrix
     */
    private static FeatureMatrix createMatrix(String[] aGroupNames, float[]... aColumns) {
        int tmpRowCount = aColumns[0].length;
        float[][] tmpData = new float[tmpRowCount][aColumns.length];
        String[] tmpRowNames = new String[tmpRowCount];
        for (int i = 0; i < tmpRowCount; i++) {
            tmpRowNames[i] = "m" + i;
            for (int j = 0; j < aColumns.length; j++) {
                tmpData[i][j] = aColumns[j][i];
            }
        }
        FeatureType[] tmpTypes = new FeatureType[aColumns.length];
        java.util.Arrays.fill(tmpTypes, FeatureType.CONTINUOUS);
        return new FeatureMatrix(tmpData, tmpRowNames, aGroupNames.clone(), tmpTypes, aGroupNames);
    }
    //
    /**
     * Creates a priority table from the given key-value pairs.
     *
     * @param aKeyValuePairs alternating group names and priorities
     * @return priority table
     */
    private static DescriptorPriorityTable createTable(String... aKeyValuePairs) {
        Properties tmpProperties = new Properties();
        for (int i = 0; i < aKeyValuePairs.length; i += 2) {
            tmpProperties.setProperty(aKeyValuePairs[i], aKeyValuePairs[i + 1]);
        }
        return new DescriptorPriorityTable(tmpProperties);
    }
    //
    /**
     * Tests that of two perfectly correlated columns the one with the higher priority is kept, independent of the
     * column order.
     */
    @Test
    void testCorrelatedColumnWithLowerPriorityIsRemoved() {
        float[] tmpA = new float[] {1, 2, 3, 4, 5};
        float[] tmpB = new float[] {2, 4, 6, 8, 10};
        float[] tmpC = new float[] {5, 1, 4, 2, 3};
        FeatureMatrix tmpMatrix = FeatureFilterTest.createMatrix(new String[] {"A", "B", "C"}, tmpA, tmpB, tmpC);
        FeatureFilter.FilterResult tmpResult = FeatureFilter.filter(tmpMatrix,
                FeatureFilterTest.createTable("B", "10"), 0.0, 0.9, FeatureFilter.CorrelationMethod.PEARSON);
        Assertions.assertEquals(2, tmpResult.filteredMatrix().getColumnCount());
        Assertions.assertEquals("B", tmpResult.filteredMatrix().getFeatureName(0));
        Assertions.assertEquals("C", tmpResult.filteredMatrix().getFeatureName(1));
        Assertions.assertEquals(1, tmpResult.removedFeatures().size());
        RemovedFeature tmpRemoved = tmpResult.removedFeatures().get(0);
        Assertions.assertEquals("A", tmpRemoved.featureName());
        Assertions.assertEquals(RemovedFeature.RemovalReason.CORRELATION, tmpRemoved.reason());
        Assertions.assertEquals("B", tmpRemoved.partnerFeatureName());
        Assertions.assertEquals(1.0, tmpRemoved.correlation(), 1.0e-9);
    }
    //
    /**
     * Tests that for equal priority the first column is kept.
     */
    @Test
    void testEqualPriorityKeepsFirstColumn() {
        float[] tmpA = new float[] {1, 2, 3, 4, 5};
        float[] tmpB = new float[] {5, 4, 3, 2, 1};
        FeatureMatrix tmpMatrix = FeatureFilterTest.createMatrix(new String[] {"A", "B"}, tmpA, tmpB);
        FeatureFilter.FilterResult tmpResult = FeatureFilter.filter(tmpMatrix,
                FeatureFilterTest.createTable(), 0.0, 0.9, FeatureFilter.CorrelationMethod.PEARSON);
        Assertions.assertEquals(1, tmpResult.filteredMatrix().getColumnCount());
        Assertions.assertEquals("A", tmpResult.filteredMatrix().getFeatureName(0));
        //negative correlation counts as well
        Assertions.assertEquals(-1.0, tmpResult.removedFeatures().get(0).correlation(), 1.0e-9);
    }
    //
    /**
     * Tests that Spearman detects a monotone but non-linear dependence that Pearson does not exceed the threshold for.
     */
    @Test
    void testSpearmanDetectsMonotoneDependence() {
        float[] tmpA = new float[] {1, 2, 3, 4, 5, 6};
        float[] tmpB = new float[] {1, 2, 3, 4, 5, 1000};
        Assertions.assertTrue(FeatureFilter.getCorrelation(tmpA, tmpB, FeatureFilter.CorrelationMethod.PEARSON) < 0.9);
        Assertions.assertEquals(1.0, FeatureFilter.getCorrelation(tmpA, tmpB, FeatureFilter.CorrelationMethod.SPEARMAN), 1.0e-9);
    }
    //
    /**
     * Tests the removal of columns with too many NaN values, the median imputation and the removal of constant columns.
     */
    @Test
    void testMissingValuesAndConstantColumns() {
        float[] tmpManyNaN = new float[] {Float.NaN, Float.NaN, 1, 2};
        float[] tmpOneNaN = new float[] {1, Float.NaN, 3, 10};
        float[] tmpConstant = new float[] {7, 7, 7, 7};
        FeatureMatrix tmpMatrix = FeatureFilterTest.createMatrix(new String[] {"N", "O", "K"}, tmpManyNaN, tmpOneNaN, tmpConstant);
        FeatureFilter.FilterResult tmpResult = FeatureFilter.filter(tmpMatrix,
                FeatureFilterTest.createTable(), 0.25, 1.0, FeatureFilter.CorrelationMethod.PEARSON);
        Assertions.assertEquals(1, tmpResult.filteredMatrix().getColumnCount());
        Assertions.assertEquals("O", tmpResult.filteredMatrix().getFeatureName(0));
        //median of 1, 3, 10
        Assertions.assertEquals(3.0f, tmpResult.filteredMatrix().getValue(1, 0));
        Assertions.assertEquals(RemovedFeature.RemovalReason.MISSING_VALUES, tmpResult.removedFeatures().get(0).reason());
        Assertions.assertEquals(RemovedFeature.RemovalReason.CONSTANT, tmpResult.removedFeatures().get(1).reason());
    }
    //
    /**
     * Tests that BIT columns do not take part in the correlation filter.
     */
    @Test
    void testBitColumnsAreNotCorrelationFiltered() {
        float[][] tmpData = new float[][] {{1, 1}, {0, 0}, {1, 1}};
        FeatureMatrix tmpMatrix = new FeatureMatrix(tmpData, new String[] {"m1", "m2", "m3"}, new String[] {"b1", "b2"},
                new FeatureType[] {FeatureType.BIT, FeatureType.BIT}, new String[] {"F", "F"});
        FeatureFilter.FilterResult tmpResult = FeatureFilter.filter(tmpMatrix,
                FeatureFilterTest.createTable(), 0.0, 0.5, FeatureFilter.CorrelationMethod.PEARSON);
        Assertions.assertEquals(2, tmpResult.filteredMatrix().getColumnCount());
        Assertions.assertTrue(tmpResult.removedFeatures().isEmpty());
    }
    //
    /**
     * Tests that the default priority table can be loaded and returns the default priority for unknown groups.
     *
     * @throws Exception if the table cannot be loaded
     */
    @Test
    void testLoadDefaultPriorityTable() throws Exception {
        DescriptorPriorityTable tmpTable = DescriptorPriorityTable.loadDefaultTable();
        Assertions.assertEquals(DescriptorPriorityTable.DEFAULT_PRIORITY, tmpTable.getPriority("UNKNOWN_GROUP"));
    }
}
