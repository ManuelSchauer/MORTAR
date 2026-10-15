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

import de.unijena.cheminf.mortar.model.descriptor.Descriptor;

import java.util.Collections;
import java.util.EnumSet;
import java.util.Objects;
import java.util.Set;

/**
 * A named set of descriptors that can be selected at once in the clustering settings. Presets are either built in
 * (see {@link ClusteringService.DescriptorPresetOption}) or defined and saved by the user. Fingerprints are never part
 * of a preset, they are chosen separately.
 *
 * @param key unique key of the preset, used to persist the selected preset
 * @param displayName name for display in the GUI
 * @param tooltip tooltip text for display in the GUI
 * @param descriptors descriptors of the preset
 * @param isUserDefined true for presets saved by the user, false for built-in presets
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public record DescriptorPreset(String key, String displayName, String tooltip, Set<Descriptor> descriptors, boolean isUserDefined) {
    /**
     * Prefix of the keys of user-defined presets, followed by the preset name.
     */
    public static final String USER_PRESET_KEY_PREFIX = "USER:";
    //
    /**
     * Compact constructor that checks the parameters and makes the descriptor set unmodifiable.
     *
     * @throws NullPointerException if a parameter is null
     */
    public DescriptorPreset {
        Objects.requireNonNull(key, "Given key is null.");
        Objects.requireNonNull(displayName, "Given display name is null.");
        Objects.requireNonNull(tooltip, "Given tooltip is null.");
        Objects.requireNonNull(descriptors, "Given descriptor set is null.");
        descriptors = descriptors.isEmpty()
                ? Collections.emptySet()
                : Collections.unmodifiableSet(EnumSet.copyOf(descriptors));
    }
    //
    /**
     * Returns whether the given descriptor belongs to this preset.
     *
     * @param aDescriptor descriptor
     * @return true if the descriptor belongs to the preset
     */
    public boolean includes(Descriptor aDescriptor) {
        return this.descriptors.contains(aDescriptor);
    }
    //
    /**
     * Returns the key of a user-defined preset with the given name.
     *
     * @param aName preset name
     * @return preset key
     */
    public static String getUserPresetKey(String aName) {
        return DescriptorPreset.USER_PRESET_KEY_PREFIX + aName;
    }
}
