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

package de.unijena.cheminf.mortar.gui.controls;

import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;
import de.unijena.cheminf.mortar.model.clustering.DescriptorPreset;

import javafx.scene.control.ComboBox;
import javafx.scene.control.ListCell;

import java.util.Objects;

/**
 * Combo box for the descriptor presets of a {@link ClusteringService}. It shows the own presets first, then the built-in
 * presets and the custom entry, and is bound to the selected preset of the service, so that all instances (in the
 * clustering settings and in the bottom bars of the fragments and items tabs) show the same selection. Choosing a preset
 * sets the descriptor check boxes of the service, also if the preset was already selected.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class DescriptorPresetComboBox extends ComboBox<String> {
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Clustering service whose presets are shown.
     */
    private final ClusteringService clusteringService;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aClusteringService clustering service whose presets are shown
     * @throws NullPointerException if the service is null
     */
    public DescriptorPresetComboBox(ClusteringService aClusteringService) throws NullPointerException {
        super(Objects.requireNonNull(aClusteringService, "Given clustering service is null.").getDescriptorPresetKeyList());
        this.clusteringService = aClusteringService;
        this.setPrefWidth(GuiDefinitions.GUI_SETTING_COMBO_BOX_PREF_WIDTH_VALUE);
        this.setTooltip(GuiUtil.createTooltip(Message.get("ClusteringService.descriptorPresetSetting.tooltip")));
        this.setCellFactory(aListView -> this.createPresetListCell());
        this.setButtonCell(this.createPresetListCell());
        this.valueProperty().bindBidirectional(aClusteringService.descriptorPresetKeyProperty());
        //choosing a preset (also the one already selected) sets the check boxes; custom keeps them
        this.setOnAction(anEvent -> {
            String tmpKey = this.getValue();
            if (this.clusteringService.getDescriptorPreset(tmpKey) != null) {
                this.clusteringService.selectDescriptorPreset(tmpKey);
            }
        });
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Creates a list cell that shows the display name and tooltip of a preset key.
     *
     * @return list cell
     */
    private ListCell<String> createPresetListCell() {
        return new ListCell<>() {
            @Override
            protected void updateItem(String aKey, boolean anIsEmpty) {
                super.updateItem(aKey, anIsEmpty);
                if (anIsEmpty || aKey == null) {
                    this.setText(null);
                    this.setTooltip(null);
                    return;
                }
                DescriptorPreset tmpPreset = DescriptorPresetComboBox.this.clusteringService.getDescriptorPreset(aKey);
                if (tmpPreset == null) {
                    this.setText(Message.get("ClusteringService.DescriptorPresetOption.custom.displayName"));
                    this.setTooltip(GuiUtil.createTooltip(Message.get("ClusteringService.DescriptorPresetOption.custom.tooltip")));
                } else {
                    this.setText(tmpPreset.displayName());
                    this.setTooltip(GuiUtil.createTooltip(tmpPreset.tooltip()));
                }
            }
        };
    }
    //</editor-fold>
}
