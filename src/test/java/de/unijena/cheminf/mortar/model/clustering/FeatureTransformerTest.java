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

/**
 * Tests for {@link FeatureTransformer}.
 *
 * @author Manuel Schauer
 */
class FeatureTransformerTest {
    /**
     * Creates a feature matrix with one BIT, one COUNT and one CONTINUOUS column.
     *
     * @return test matrix
     */
    private static FeatureMatrix createTestMatrix() {
        float[][] tmpData = new float[][] {
                {1.0f, 0.0f, -2.0f},
                {0.0f, 3.0f, 1.0f},
                {1.0f, Float.NaN, 100.0f}
        };
        return new FeatureMatrix(tmpData, new String[] {"m1", "m2", "m3"},
                new String[] {"bit", "count", "continuous"},
                new FeatureType[] {FeatureType.BIT, FeatureType.COUNT, FeatureType.CONTINUOUS},
                new String[] {"G", "G", "G"});
    }
    //
    /**
     * Tests that LOG1P is applied to COUNT columns only, that BIT columns stay unchanged and NaN stays NaN.
     */
    @Test
    void testLog1pAppliedToCountColumnsOnly() {
        FeatureMatrix tmpMatrix = FeatureTransformerTest.createTestMatrix();
        FeatureMatrix tmpResult = FeatureTransformer.transform(tmpMatrix, FeatureTransformation.LOG1P, FeatureTransformation.NONE);
        Assertions.assertEquals(1.0f, tmpResult.getValue(0, 0));
        Assertions.assertEquals((float) Math.log1p(3.0), tmpResult.getValue(1, 1), 1.0e-6f);
        Assertions.assertTrue(Float.isNaN(tmpResult.getValue(2, 1)));
        Assertions.assertEquals(-2.0f, tmpResult.getValue(0, 2));
        //the input matrix is not changed
        Assertions.assertEquals(3.0f, tmpMatrix.getValue(1, 1));
    }
    //
    /**
     * Tests that no copy is created if no column is transformed and that a copy is created otherwise.
     */
    @Test
    void testNoCopyWithoutTransformation() {
        FeatureMatrix tmpMatrix = FeatureTransformerTest.createTestMatrix();
        Assertions.assertFalse(FeatureTransformer.isTransforming(tmpMatrix, FeatureTransformation.NONE, FeatureTransformation.NONE));
        Assertions.assertSame(tmpMatrix, FeatureTransformer.transform(tmpMatrix, FeatureTransformation.NONE, FeatureTransformation.NONE));
        Assertions.assertTrue(FeatureTransformer.isTransforming(tmpMatrix, FeatureTransformation.NONE, FeatureTransformation.WINSORIZE));
        Assertions.assertNotSame(tmpMatrix, FeatureTransformer.transform(tmpMatrix, FeatureTransformation.NONE, FeatureTransformation.WINSORIZE));
        //only BIT columns: a count transformation changes nothing
        FeatureMatrix tmpBitMatrix = tmpMatrix.selectColumns(new int[] {0});
        Assertions.assertSame(tmpBitMatrix, FeatureTransformer.transform(tmpBitMatrix, FeatureTransformation.LOG1P, FeatureTransformation.NONE));
    }
    //
    /**
     * Tests that BINARIZE turns counts into presence bits.
     */
    @Test
    void testBinarize() {
        FeatureMatrix tmpResult = FeatureTransformer.transform(FeatureTransformerTest.createTestMatrix(),
                FeatureTransformation.BINARIZE, FeatureTransformation.NONE);
        Assertions.assertEquals(0.0f, tmpResult.getValue(0, 1));
        Assertions.assertEquals(1.0f, tmpResult.getValue(1, 1));
    }
    //
    /**
     * Tests that LOG1P on a column with negative values throws an exception.
     */
    @Test
    void testLog1pOnNegativeValuesThrows() {
        Assertions.assertThrows(IllegalArgumentException.class, () -> FeatureTransformer.transform(
                FeatureTransformerTest.createTestMatrix(), FeatureTransformation.NONE, FeatureTransformation.LOG1P));
    }
    //
    /**
     * Tests that WINSORIZE clips an outlier to the upper percentile.
     */
    @Test
    void testWinsorizeClipsOutlier() {
        float[] tmpColumn = new float[101];
        for (int i = 0; i < 100; i++) {
            tmpColumn[i] = i;
        }
        tmpColumn[100] = 10000.0f;
        FeatureTransformer.transformColumn(tmpColumn, FeatureTransformation.WINSORIZE, "test");
        Assertions.assertEquals(99.0f, tmpColumn[100], 1.0e-4f);
        Assertions.assertEquals(1.0f, tmpColumn[0], 1.0e-4f);
        Assertions.assertEquals(50.0f, tmpColumn[50], 1.0e-4f);
    }
    //
    /**
     * Tests the percentile calculation with linear interpolation.
     */
    @Test
    void testGetPercentile() {
        float[] tmpSorted = new float[] {1.0f, 2.0f, 3.0f, 4.0f};
        Assertions.assertEquals(2.5f, FeatureTransformer.getPercentile(tmpSorted, 0.5), 1.0e-6f);
        Assertions.assertEquals(1.0f, FeatureTransformer.getPercentile(tmpSorted, 0.0), 1.0e-6f);
        Assertions.assertEquals(4.0f, FeatureTransformer.getPercentile(tmpSorted, 1.0), 1.0e-6f);
    }
}
