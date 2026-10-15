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

package de.unijena.cheminf.mortar.controller;

import de.unijena.cheminf.mortar.configuration.IConfiguration;
import de.unijena.cheminf.mortar.gui.controls.DescriptorPresetComboBox;
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.gui.views.SettingsView;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;
import de.unijena.cheminf.mortar.model.clustering.DescriptorPreset;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import javafx.beans.binding.Bindings;
import javafx.beans.property.Property;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.TextInputDialog;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.Modality;
import javafx.stage.Stage;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Controls a {@link SettingsView} for the clustering settings. Both tabs consist of visually separated sections whose
 * rows are created from the properties of the {@link ClusteringService} like the fragmentation settings:
 * the clustering tab (fragment fingerprint, aromaticity, transformation, filter with the correlation settings only shown
 * if the correlation filter is switched on, ART-2A) and the descriptor tab (preset with saving and deleting of own
 * presets, fingerprints with hidden advanced settings, descriptor list). A button at the bottom applies the settings
 * and starts the clustering. In the descriptor calculation mode, only the aromaticity, preset and descriptor sections
 * are shown (fingerprints are not part of a descriptor table) and the button starts the descriptor calculation; the
 * selection and the aromaticity model are the same as for the clustering.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringSettingsViewController {
    //<editor-fold desc="private static final class constants">
    /**
     * Style of the section titles of the descriptor tab.
     */
    private static final String SECTION_TITLE_STYLE = "-fx-font-weight: bold; -fx-font-size: 14px;";
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusteringSettingsViewController.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="private and private final class variables">
    /**
     * Main stage object of the application.
     */
    private final Stage mainStage;
    /**
     * Stage for the SettingsView.
     */
    private Stage clusteringSettingsViewStage;
    /**
     * SettingsView.
     */
    private SettingsView settingsView;
    /**
     * Clustering service whose settings are displayed.
     */
    private final ClusteringService clusteringService;
    /**
     * Values of the general settings when the view was opened, used to restore them on cancel.
     */
    private final Map<String, Object> recentSettingsMap;
    /**
     * Values of all settings of the descriptor tab when the view was opened, used to restore them on cancel.
     */
    private final Map<String, Object> recentDescriptorSettingsMap;
    /**
     * Combo box of the descriptor presets, holds the preset keys.
     */
    private DescriptorPresetComboBox presetComboBox;
    /**
     * Action that starts the clustering, may be null if no clustering can be started.
     */
    private final Runnable startClusteringAction;
    /**
     * True if the view is opened to choose the descriptors of a descriptor calculation instead of the clustering settings.
     */
    private final boolean isDescriptorCalculationMode;
    /**
     * Configuration class to read resource file paths from.
     */
    private final IConfiguration configuration;
    //</editor-fold>
    //
    /**
     * Constructor, opens the view.
     *
     * @param aStage main stage
     * @param aClusteringService clustering service whose settings are displayed
     * @param aConfiguration configuration instance to read resource file paths from
     * @param aStartClusteringAction action that starts the clustering with the applied settings; null disables the start
     *                               button, e.g. if there is no fragmentation yet
     * @throws NullPointerException if the stage, service or configuration is null
     */
    public ClusteringSettingsViewController(Stage aStage, ClusteringService aClusteringService, IConfiguration aConfiguration,
                                            Runnable aStartClusteringAction) throws NullPointerException {
        this(aStage, aClusteringService, aConfiguration, aStartClusteringAction, false);
    }
    //
    /**
     * Constructor, opens the view either with the clustering settings or for the choice of the descriptors of a
     * descriptor calculation.
     *
     * @param aStage main stage
     * @param aClusteringService clustering service whose settings are displayed
     * @param aConfiguration configuration instance to read resource file paths from
     * @param aStartAction action that starts the clustering or the descriptor calculation with the applied settings;
     *                     null disables the start button
     * @param anIsDescriptorCalculationMode true to show only the preset and descriptor sections and to start the
     *                                      descriptor calculation
     * @throws NullPointerException if the stage, service or configuration is null
     */
    public ClusteringSettingsViewController(Stage aStage, ClusteringService aClusteringService, IConfiguration aConfiguration,
                                            Runnable aStartAction, boolean anIsDescriptorCalculationMode) throws NullPointerException {
        this.startClusteringAction = aStartAction;
        this.isDescriptorCalculationMode = anIsDescriptorCalculationMode;
        this.mainStage = Objects.requireNonNull(aStage, "aStage (instance of Stage) is null");
        this.clusteringService = Objects.requireNonNull(aClusteringService, "aClusteringService is null");
        this.configuration = Objects.requireNonNull(aConfiguration, "aConfiguration is null");
        this.recentSettingsMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                aClusteringService.settingsProperties().size()));
        this.recentDescriptorSettingsMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                aClusteringService.descriptorSettingsProperties().size()));
        this.openClusteringSettingsView();
    }
    //
    /**
     * Initialises and opens the settings view.
     */
    private void openClusteringSettingsView() {
        this.settingsView = new SettingsView();
        this.clusteringSettingsViewStage = new Stage();
        Scene tmpScene = new Scene(this.settingsView, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE, GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        this.clusteringSettingsViewStage.setScene(tmpScene);
        this.clusteringSettingsViewStage.initModality(Modality.WINDOW_MODAL);
        this.clusteringSettingsViewStage.initOwner(this.mainStage);
        this.clusteringSettingsViewStage.setTitle(Message.get(this.isDescriptorCalculationMode
                ? "ClusteringSettingsView.descriptorCalculation.title.text" : "ClusteringSettingsView.title.text"));
        this.clusteringSettingsViewStage.setMinHeight(GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        this.clusteringSettingsViewStage.setMinWidth(GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE);
        String tmpIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder")
                        + this.configuration.getProperty("mortar.logo.icon.name")).toExternalForm();
        this.clusteringSettingsViewStage.getIcons().add(new Image(tmpIconURL));
        //the recent values are recorded here because the tabs consist of several sections
        for (Property<?> tmpProperty : this.clusteringService.settingsProperties()) {
            this.recentSettingsMap.put(tmpProperty.getName(), tmpProperty.getValue());
        }
        for (Property<?> tmpProperty : this.clusteringService.descriptorSettingsProperties()) {
            this.recentDescriptorSettingsMap.put(tmpProperty.getName(), tmpProperty.getValue());
        }
        if (!this.isDescriptorCalculationMode) {
            this.settingsView.addTab(Message.get("ClusteringSettingsView.generalTab.text"), this.createGeneralTabContent());
        }
        this.settingsView.addTab(Message.get("ClusteringSettingsView.descriptorTab.text"), this.createDescriptorTabContent());
        this.addListener();
        this.clusteringSettingsViewStage.show();
    }
    //
    /**
     * Creates the content of the clustering tab with the sections fragment fingerprint, aromaticity, transformation,
     * filter (correlation settings only visible if the correlation filter is switched on) and ART-2A.
     *
     * @return content of the clustering tab
     */
    private VBox createGeneralTabContent() {
        VBox tmpContent = new VBox(GuiDefinitions.GUI_SPACING_VALUE);
        tmpContent.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        tmpContent.getChildren().addAll(
                this.createSectionTitle("ClusteringSettingsView.fragmentFingerprintSection.title"),
                this.createGrid(this.clusteringService.fragmentFingerprintSettingsProperties()),
                this.createSectionTitle("ClusteringSettingsView.aromaticitySection.title"),
                this.createGrid(this.clusteringService.aromaticitySettingsProperties()),
                this.createSectionTitle("ClusteringSettingsView.transformationSection.title"),
                this.createGrid(this.clusteringService.transformationSettingsProperties()),
                this.createSectionTitle("ClusteringSettingsView.filterSection.title"),
                this.createGrid(this.clusteringService.filterSettingsProperties()));
        //the correlation settings are only shown if the correlation filter is switched on
        GridPane tmpCorrelationGrid = this.createGrid(this.clusteringService.correlationSettingsProperties());
        tmpCorrelationGrid.visibleProperty().bind(this.clusteringService.correlationFilterProperty());
        tmpCorrelationGrid.managedProperty().bind(this.clusteringService.correlationFilterProperty());
        tmpContent.getChildren().addAll(tmpCorrelationGrid,
                this.createSectionTitle("ClusteringSettingsView.art2aSection.title"),
                this.createGrid(this.clusteringService.art2aSettingsProperties()),
                this.createSectionTitle("ClusteringSettingsView.similaritySearchSection.title"),
                this.createGrid(this.clusteringService.similaritySearchSettingsProperties()));
        return tmpContent;
    }
    //
    /**
     * Creates a grid pane with one row per given setting via {@link SettingsView#createPropertyGridPane}.
     *
     * @param aSettingsList settings
     * @return grid pane
     */
    private GridPane createGrid(List<Property<?>> aSettingsList) {
        //the recent values are already recorded in the constructor
        return this.settingsView.createPropertyGridPane(aSettingsList,
                this.clusteringService.getSettingNameToDisplayNameMap(),
                this.clusteringService.getSettingNameToTooltipTextMap(),
                new HashMap<>());
    }
    //
    /**
     * Creates the content of the descriptor tab with the sections aromaticity model (only in the descriptor calculation
     * mode, otherwise in the clustering tab), preset, fingerprints (not in the descriptor calculation mode) and
     * descriptors.
     *
     * @return content of the descriptor tab
     */
    private VBox createDescriptorTabContent() {
        VBox tmpContent = new VBox(GuiDefinitions.GUI_SPACING_VALUE);
        tmpContent.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        //the aromaticity model is applied before the descriptor calculation, so it is offered here without the clustering tab
        if (this.isDescriptorCalculationMode) {
            tmpContent.getChildren().addAll(this.createSectionTitle("ClusteringSettingsView.aromaticitySection.title"),
                    this.createGrid(this.clusteringService.aromaticitySettingsProperties()));
        }
        //<editor-fold desc="preset section">
        this.presetComboBox = new DescriptorPresetComboBox(this.clusteringService);
        this.presetComboBox.setPrefWidth(GuiDefinitions.GUI_SETTING_COMBO_BOX_PREF_WIDTH_VALUE * 1.5);
        Button tmpSavePresetButton = new Button(Message.get("ClusteringSettingsView.savePresetButton.text"));
        tmpSavePresetButton.setTooltip(new Tooltip(Message.get("ClusteringSettingsView.savePresetButton.tooltip")));
        tmpSavePresetButton.setOnAction(anEvent -> this.saveDescriptorSelectionAsPreset());
        Button tmpDeletePresetButton = new Button(Message.get("ClusteringSettingsView.deletePresetButton.text"));
        tmpDeletePresetButton.setTooltip(new Tooltip(Message.get("ClusteringSettingsView.deletePresetButton.tooltip")));
        tmpDeletePresetButton.setOnAction(anEvent -> this.deleteSelectedPreset());
        //only own presets can be deleted
        tmpDeletePresetButton.disableProperty().bind(Bindings.createBooleanBinding(() -> {
            DescriptorPreset tmpPreset = this.clusteringService.getDescriptorPreset(this.presetComboBox.getValue());
            return tmpPreset == null || !tmpPreset.isUserDefined();
        }, this.presetComboBox.valueProperty()));
        HBox tmpPresetBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, this.presetComboBox, tmpSavePresetButton, tmpDeletePresetButton);
        tmpPresetBox.setAlignment(Pos.CENTER_LEFT);
        tmpPresetBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        Label tmpPresetLabel = new Label(Message.get("ClusteringService.descriptorPresetSetting.tooltip"));
        tmpPresetLabel.setWrapText(true);
        tmpContent.getChildren().addAll(this.createSectionTitle("ClusteringSettingsView.presetSection.title"),
                tmpPresetLabel, tmpPresetBox);
        //</editor-fold>
        //<editor-fold desc="fingerprint section">
        GridPane tmpFingerprintGrid = this.createGrid(this.clusteringService.fingerprintSettingsProperties());
        //radius and length of the circular fingerprint are hidden until the user asks for them
        GridPane tmpAdvancedGrid = this.createGrid(this.clusteringService.advancedFingerprintSettingsProperties());
        CheckBox tmpAdvancedCheckBox = new CheckBox(Message.get("ClusteringSettingsView.advancedFingerprintSettings.text"));
        tmpAdvancedCheckBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        tmpAdvancedGrid.visibleProperty().bind(tmpAdvancedCheckBox.selectedProperty());
        tmpAdvancedGrid.managedProperty().bind(tmpAdvancedCheckBox.selectedProperty());
        if (!this.isDescriptorCalculationMode) {
            tmpContent.getChildren().addAll(this.createSectionTitle("ClusteringSettingsView.fingerprintSection.title"),
                    tmpFingerprintGrid, tmpAdvancedCheckBox, tmpAdvancedGrid);
        }
        //</editor-fold>
        //<editor-fold desc="descriptor section">
        GridPane tmpDescriptorGrid = this.createGrid(this.clusteringService.descriptorSelectionSettingsProperties());
        tmpContent.getChildren().addAll(this.createSectionTitle("ClusteringSettingsView.descriptorSection.title"),
                tmpDescriptorGrid);
        //</editor-fold>
        return tmpContent;
    }
    //
    /**
     * Creates a section title with a separator below.
     *
     * @param aMessageKey message key of the title text
     * @return box with title and separator
     */
    private VBox createSectionTitle(String aMessageKey) {
        Label tmpTitle = new Label(Message.get(aMessageKey));
        tmpTitle.setStyle(ClusteringSettingsViewController.SECTION_TITLE_STYLE);
        VBox tmpBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpTitle, new Separator());
        tmpBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, 0, 0, 0));
        return tmpBox;
    }
    //
    /**
     * Asks for a name and saves the currently checked descriptors as own preset.
     */
    private void saveDescriptorSelectionAsPreset() {
        TextInputDialog tmpDialog = new TextInputDialog(this.clusteringService.getProposedUserPresetName());
        tmpDialog.initOwner(this.clusteringSettingsViewStage);
        tmpDialog.setTitle(Message.get("ClusteringSettingsView.savePresetDialog.title"));
        tmpDialog.setHeaderText(Message.get("ClusteringSettingsView.savePresetDialog.header"));
        Optional<String> tmpName = tmpDialog.showAndWait();
        if (tmpName.isEmpty()) {
            return;
        }
        try {
            String tmpKey = this.clusteringService.saveCurrentDescriptorSelectionAsPreset(tmpName.get());
            this.clusteringService.descriptorPresetKeyProperty().set(tmpKey);
        } catch (IllegalArgumentException | IOException anException) {
            ClusteringSettingsViewController.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                    Message.get("ClusteringSettingsView.presetError.header"), anException.getMessage(), anException);
        }
    }
    //
    /**
     * Deletes the selected own preset; the check boxes are kept and the selection becomes custom.
     */
    private void deleteSelectedPreset() {
        try {
            this.clusteringService.deleteUserDescriptorPreset(this.presetComboBox.getValue());
        } catch (IllegalArgumentException | IOException anException) {
            ClusteringSettingsViewController.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                    Message.get("ClusteringSettingsView.presetError.header"), anException.getMessage(), anException);
        }
    }
    //
    /**
     * Adds the listeners of the buttons and the close request.
     */
    private void addListener() {
        //close request discards the changes, like the cancel button
        this.clusteringSettingsViewStage.setOnCloseRequest(event -> this.restoreRecentSettings());
        this.settingsView.getApplyButton().setOnAction(event -> this.clusteringSettingsViewStage.close());
        //applies the settings and starts the clustering or the descriptor calculation
        Button tmpStartButton = new Button(Message.get(this.isDescriptorCalculationMode
                ? "ClusteringSettingsView.calculateDescriptorsButton.text" : "ClusteringSettingsView.startClusteringButton.text"));
        tmpStartButton.setTooltip(GuiUtil.createTooltip(Message.get(this.isDescriptorCalculationMode
                ? "ClusteringSettingsView.calculateDescriptorsButton.tooltip" : "ClusteringSettingsView.startClusteringButton.tooltip")));
        tmpStartButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpStartButton.setDisable(this.startClusteringAction == null);
        tmpStartButton.setOnAction(event -> {
            this.clusteringSettingsViewStage.close();
            this.startClusteringAction.run();
        });
        this.settingsView.addButtonToRightSide(tmpStartButton);
        this.settingsView.getCancelButton().setOnAction(event -> {
            this.restoreRecentSettings();
            this.clusteringSettingsViewStage.close();
        });
        //the default button restores the defaults of the currently selected tab
        this.settingsView.getDefaultButton().setOnAction(event -> {
            if (this.isDescriptorCalculationMode) {
                this.clusteringService.restoreDefaultAromaticitySettings();
                this.clusteringService.restoreDefaultDescriptorSelection();
            } else if (this.settingsView.getSelectionModel().getSelectedIndex() == 0) {
                this.clusteringService.restoreDefaultSettings();
            } else {
                this.clusteringService.restoreDefaultDescriptorSelection();
            }
        });
    }
    //
    /**
     * Sets all settings back to the values they had when the view was opened. Saved or deleted own presets are not
     * restored, they are written immediately.
     */
    private void restoreRecentSettings() {
        ClusteringSettingsViewController.setRecentProperties(this.clusteringService.settingsProperties(), this.recentSettingsMap);
        //the preset is restored first, then the check boxes (see order of descriptorSettingsProperties())
        ClusteringSettingsViewController.setRecentProperties(this.clusteringService.descriptorSettingsProperties(),
                this.recentDescriptorSettingsMap);
    }
    //
    /**
     * Sets the given properties to the values of the given map.
     *
     * @param aPropertyList properties
     * @param aRecentPropertiesMap map of property name to recent value
     */
    private static void setRecentProperties(List<Property<?>> aPropertyList, Map<String, Object> aRecentPropertiesMap) {
        for (Property tmpProperty : aPropertyList) {
            if (aRecentPropertiesMap.containsKey(tmpProperty.getName())) {
                tmpProperty.setValue(aRecentPropertiesMap.get(tmpProperty.getName()));
            }
        }
    }
}
