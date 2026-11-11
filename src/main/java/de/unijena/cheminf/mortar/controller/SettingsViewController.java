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
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.views.SettingsView;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.filter.MoleculeFilterSettings;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import javafx.application.Platform;
import javafx.beans.property.Property;
import javafx.scene.Scene;
import javafx.scene.control.Tab;
import javafx.scene.image.Image;
import javafx.stage.Modality;
import javafx.stage.Stage;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * SettingsViewController controls {@link SettingsView} for {@link SettingsContainer} (first tab) and
 * {@link MoleculeFilterSettings} (second tab, import filters).
 *
 * @author Felix Baensch
 * @version 1.0.0.0
 */
public class SettingsViewController {
    //<editor-fold desc="private and private final class variables">
    /**
     * SettingsContainer.
     */
    private final SettingsContainer settingsContainer;
    /**
     * Settings of the molecule filters applied at import, shown in a second tab.
     */
    private final MoleculeFilterSettings moleculeFilterSettings;
    /**
     * Main stage object of the application.
     */
    private final Stage mainStage;
    /**
     * Configuration class to read resource file paths from.
     */
    private final IConfiguration configuration;
    /**
     * Stage for the SettingsView.
     */
    private Stage settingsViewStage;
    /**
     * SettingsView.
     */
    private SettingsView settingsView;
    /**
     * Map to hold the initial settings properties.
     */
    private final Map<String, Object> recentProperties;
    /**
     * Map to hold the initial values of the molecule filter settings.
     */
    private final Map<String, Object> recentMoleculeFilterProperties;
    /**
     * Tab of the molecule filter settings; null until the tabs are added.
     */
    private Tab moleculeFilterSettingsTab;
    /**
     * Boolean value to check if the rowsPerPage property has changed.
     */
    private boolean hasRowsPerPageChanged;
    /**
     * Boolean value to check if the keepAtomContainerInDataModel property has changed.
     */
    private boolean hasKeepAtomContainerInDataModelChanged;
    //</editor-fold>
    //
    /**
     * Constructor.
     *
     * @param aStage Parent stage
     * @param aSettingsContainer SettingsContainer
     * @param aMoleculeFilterSettings settings of the molecule filters applied at import
     * @param aConfiguration configuration instance to read resource file paths from
     */
    public SettingsViewController(Stage aStage, SettingsContainer aSettingsContainer,
                                  MoleculeFilterSettings aMoleculeFilterSettings, IConfiguration aConfiguration) {
        this.mainStage = aStage;
        this.settingsContainer = aSettingsContainer;
        this.moleculeFilterSettings = aMoleculeFilterSettings;
        this.configuration = aConfiguration;
        this.recentProperties = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(this.settingsContainer.settingsProperties().size()));
        this.recentMoleculeFilterProperties = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(
                this.moleculeFilterSettings.settingsProperties().size()));
        this.showSettingsView();
    }
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Initialises and opens settingsView.
     */
    private void showSettingsView() {
        if (this.settingsView == null) {
            this.settingsView = new SettingsView();
        }
        this.settingsViewStage = new Stage();
        Scene tmpScene = new Scene(this.settingsView, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE, GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        this.settingsViewStage.setScene(tmpScene);
        this.settingsViewStage.initModality(Modality.WINDOW_MODAL);
        this.settingsViewStage.initOwner(this.mainStage);
        this.settingsViewStage.setTitle(Message.get("SettingsView.title.default.text"));
        this.settingsViewStage.setMinHeight(GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        this.settingsViewStage.setMinWidth(GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE);
        String tmpIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder")
                        + this.configuration.getProperty("mortar.logo.icon.name")).toExternalForm();
        this.settingsViewStage.getIcons().add(new Image(tmpIconURL));
        Platform.runLater(()->{
            this.addListeners();
            this.settingsView.addTab(Message.get("GlobalSettingsView.title.text"),
                    this.settingsContainer.settingsProperties(), this.settingsContainer.getSettingNameToDisplayNameMap(),
                    this.settingsContainer.getSettingNameToTooltipTextMap(), this.recentProperties);
            this.moleculeFilterSettingsTab = this.settingsView.addTab(this.moleculeFilterSettings.getDisplayName(),
                    this.moleculeFilterSettings.settingsProperties(), this.moleculeFilterSettings.getSettingNameToDisplayNameMap(),
                    this.moleculeFilterSettings.getSettingNameToTooltipTextMap(), this.recentMoleculeFilterProperties);
        });
        this.settingsViewStage.showAndWait();
    }
    //
    /**
     * Adds listeners and event handlers to the buttons of the settings view.
     */
    private void addListeners() {
        //stage close request
        this.settingsViewStage.setOnCloseRequest(event -> {
            this.setRecentProperties();
            this.settingsViewStage.close();
        });
        //apply button
        this.settingsView.getApplyButton().setOnAction(event -> {
            this.hasRowsPerPageChanged = this.settingsContainer.rowsPerPageSettingProperty().getValue()
                    != (int) this.recentProperties.get(this.settingsContainer.rowsPerPageSettingProperty().getName());
            this.hasKeepAtomContainerInDataModelChanged = this.settingsContainer.keepAtomContainerInDataModelSettingProperty().getValue()
                    != this.recentProperties.get(this.settingsContainer.keepAtomContainerInDataModelSettingProperty().getName());
            this.settingsViewStage.close();
        });
        //cancel button
        this.settingsView.getCancelButton().setOnAction(event -> {
            this.setRecentProperties();
            this.settingsViewStage.close();
        });
        //default button, restores only the settings of the selected tab (like the fragmentation settings view)
        this.settingsView.getDefaultButton().setOnAction(event -> {
            if (this.moleculeFilterSettingsTab != null
                    && this.settingsView.getSelectionModel().getSelectedItem() == this.moleculeFilterSettingsTab) {
                this.moleculeFilterSettings.restoreDefaultSettings();
            } else {
                this.settingsContainer.restoreDefaultSettings();
            }
        });
    }
    //
    /**
     * Sets the properties of both tabs to the values of the respective recent properties map.
     */
    private void setRecentProperties() {
        Platform.runLater(()->{
            SettingsViewController.restoreRecentValues(this.settingsContainer.settingsProperties(), this.recentProperties);
            //note: the heavy atom count bounds can be restored in any order, see MoleculeFilterSettings
            SettingsViewController.restoreRecentValues(this.moleculeFilterSettings.settingsProperties(),
                    this.recentMoleculeFilterProperties);
        });
    }
    //
    /**
     * Sets the given properties to their values in the given recent properties map, if contained.
     *
     * @param aPropertiesList properties to restore
     * @param aRecentPropertiesMap map of setting names and recent values
     */
    private static void restoreRecentValues(List<Property<?>> aPropertiesList, Map<String, Object> aRecentPropertiesMap) {
        for (Property tmpProperty : aPropertiesList) {
            if (aRecentPropertiesMap.containsKey(tmpProperty.getName())){
                tmpProperty.setValue(aRecentPropertiesMap.get(tmpProperty.getName()));
            }
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="public properties" defaultstate="collapsed">
    /**
     * Returns boolean value whether if rowsPerPage property has changed or not.
     *
     * @return hasRowsPerPageChanged
     */
    public boolean hasRowsPerPageChanged() {
        return this.hasRowsPerPageChanged;
    }
    //
    /**
     * Returns boolean value whether if keepAtomContainerInDataModel property has changed or not.
     *
     * @return hasKeepAtomContainerInDataModelChanged
     */
    public boolean hasKeepAtomContainerInDataModelChanged() {
        return this.hasKeepAtomContainerInDataModelChanged;
    }
    //</editor-fold>
}
