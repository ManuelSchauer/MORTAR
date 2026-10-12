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

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Properties;

/**
 * Lookup table that assigns a priority to every feature group (e.g. a descriptor enum constant name). It is used by the
 * correlation filter of {@link FeatureFilter} to decide which of two highly correlated features is kept: the feature with
 * the higher priority is kept, the other one is removed. Groups that are not listed get {@link #DEFAULT_PRIORITY}.
 * <p>
 * The default table is read from {@link #DEFAULT_TABLE_RESOURCE_PATH}; custom tables can be created from any
 * {@link Properties} object (key: group name, value: integer priority).
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class DescriptorPriorityTable {
    //<editor-fold desc="public static final class constants" defaultstate="collapsed">
    /**
     * Class path of the default priority table.
     */
    public static final String DEFAULT_TABLE_RESOURCE_PATH = "de/unijena/cheminf/mortar/clustering/DescriptorPriorities.properties";
    /**
     * Priority of feature groups that are not listed in the table.
     */
    public static final int DEFAULT_PRIORITY = 0;
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Map of feature group name to priority.
     */
    private final Map<String, Integer> priorityMap;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor. Creates the table from the given properties.
     *
     * @param aProperties properties with group names as keys and integer priorities as values
     * @throws NullPointerException if the given properties are null
     * @throws IllegalArgumentException if a value is not an integer
     */
    public DescriptorPriorityTable(Properties aProperties) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aProperties, "Given properties are null.");
        this.priorityMap = new HashMap<>(aProperties.size() * 2);
        for (String tmpKey : aProperties.stringPropertyNames()) {
            String tmpValue = aProperties.getProperty(tmpKey).trim();
            try {
                this.priorityMap.put(tmpKey.trim(), Integer.parseInt(tmpValue));
            } catch (NumberFormatException anException) {
                throw new IllegalArgumentException("Priority of " + tmpKey + " is not an integer: " + tmpValue, anException);
            }
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the priority of the given feature group.
     *
     * @param aFeatureGroupName feature group name, e.g. a descriptor enum constant name
     * @return priority, {@link #DEFAULT_PRIORITY} if the group is not listed
     */
    public int getPriority(String aFeatureGroupName) {
        return this.priorityMap.getOrDefault(aFeatureGroupName, DescriptorPriorityTable.DEFAULT_PRIORITY);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public static methods" defaultstate="collapsed">
    /**
     * Loads the default priority table from the class path.
     *
     * @return default priority table
     * @throws IOException if the resource cannot be read
     */
    public static DescriptorPriorityTable loadDefaultTable() throws IOException {
        try (InputStream tmpInputStream = DescriptorPriorityTable.class.getClassLoader()
                .getResourceAsStream(DescriptorPriorityTable.DEFAULT_TABLE_RESOURCE_PATH)) {
            if (tmpInputStream == null) {
                throw new IOException("Default descriptor priority table not found: " + DescriptorPriorityTable.DEFAULT_TABLE_RESOURCE_PATH);
            }
            Properties tmpProperties = new Properties();
            tmpProperties.load(tmpInputStream);
            return new DescriptorPriorityTable(tmpProperties);
        }
    }
    //</editor-fold>
}
