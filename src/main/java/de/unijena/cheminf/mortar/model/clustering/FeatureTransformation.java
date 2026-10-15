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

/**
 * Enumeration of the non-linear transformations that can be applied to the features before clustering.
 * <p>
 * Note: The ART-2A kernels always min-max scale every column to [0,1] internally. Any positive affine scaling done
 * beforehand (z-score, min-max, multiplication with a weight) therefore has no effect on the clustering result. Only
 * non-linear transformations like the ones defined here change the clustering.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public enum FeatureTransformation {
    /**
     * Values are not changed.
     */
    NONE,
    /**
     * x is replaced by ln(1 + x); only applicable to non-negative values (counts), dampens high counts.
     */
    LOG1P,
    /**
     * x is replaced by sqrt(x); only applicable to non-negative values (counts), dampens high counts less than LOG1P.
     */
    SQRT,
    /**
     * x is replaced by 1 if x is greater than 0, else 0; turns counts into presence bits.
     */
    BINARIZE,
    /**
     * Values below the lower or above the upper percentile of the column are clipped to the respective percentile,
     * see {@link FeatureTransformer#WINSORIZE_LOWER_PERCENTILE} and {@link FeatureTransformer#WINSORIZE_UPPER_PERCENTILE}.
     * Reduces the influence of outliers on the internal min-max scaling of the kernels.
     */
    WINSORIZE;
}
