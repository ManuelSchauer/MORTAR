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
 * Describes a feature that was removed by the {@link FeatureFilter} and why. Used for the filter report in the export.
 *
 * @param featureName name of the removed feature
 * @param reason reason of the removal
 * @param partnerFeatureName name of the kept feature that the removed one is correlated with; null for other reasons
 * @param correlation correlation coefficient with the partner feature; NaN for other reasons
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record RemovedFeature(String featureName, RemovalReason reason, String partnerFeatureName, double correlation) {
    /**
     * Enumeration of the reasons why a feature can be removed.
     */
    public enum RemovalReason {
        /**
         * The fraction of missing (NaN) values of the feature exceeds the allowed maximum.
         */
        MISSING_VALUES,
        /**
         * The feature has the same value for all molecules.
         */
        CONSTANT,
        /**
         * The feature is highly correlated with a kept feature of higher (or equal) priority.
         */
        CORRELATION;
    }
}
