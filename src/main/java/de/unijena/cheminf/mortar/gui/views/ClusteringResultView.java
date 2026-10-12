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

package de.unijena.cheminf.mortar.gui.views;

import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusterProfiler;
import de.unijena.cheminf.mortar.model.clustering.ClusteringResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;

import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.image.ImageView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;

import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

/**
 * View that shows the result of a clustering run: a summary, a warning in case of a cluster overflow and a table with
 * one row per cluster ranked from the largest to the smallest cluster (rank, cluster index, size, name and structure of
 * the representative). The buttons open the overview of all cluster representatives, open the selected cluster in the
 * overview and export the clustering data. A similarity search finds the molecules most similar to a molecule given by
 * name or SMILES (search field) or to a cluster representative (context menu of the table). The actions are set by the
 * controller.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringResultView extends BorderPane {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Number of fragments and descriptors shown in the profile of the selected cluster.
     */
    private static final int NUMBER_OF_PROFILE_FEATURES = 10;
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Result that is displayed.
     */
    private final ClusteringRunResult runResult;
    /**
     * Table of the clusters, ranked by size.
     */
    private final TableView<ClusteringRunResult.ClusterSummary> clusterTableView;
    /**
     * Button to open the overview of the representatives of all clusters.
     */
    private final Button representativesOverviewButton;
    /**
     * Button to open the selected cluster in the overview view.
     */
    private final Button overviewButton;
    /**
     * Button to export the clustering data.
     */
    private final Button exportButton;
    /**
     * Button to create a training/test partitioning based on the clusters.
     */
    private final Button trainTestSplitButton;
    /**
     * Button to start a vigilance scan with the features of this clustering.
     */
    private final Button vigilanceScanButton;
    /**
     * Button to show the window of the last vigilance scan again; only visible after a scan.
     */
    private final Button showVigilanceScanButton;
    /**
     * Button to select a diverse subset (molecules per cluster).
     */
    private final Button diverseSubsetButton;
    /**
     * Text field for the name or SMILES of the query of the similarity search.
     */
    private final TextField similaritySearchTextField;
    /**
     * Button that starts the similarity search for the query of the text field.
     */
    private final Button similaritySearchButton;
    /**
     * Context menu item of the table that starts the similarity search for the representative of the selected cluster.
     */
    private final MenuItem similaritySearchMenuItem;
    /**
     * Context menu item of the table that selects typical or diverse molecules from the selected cluster.
     */
    private final MenuItem clusterMoleculeSelectionMenuItem;
    /**
     * Box at the bottom holding the button row and, if set, the clustering bar below it.
     */
    private final VBox bottomBox;
    /**
     * Profiler of the clusters, computes the profile of the selected cluster.
     */
    private final ClusterProfiler clusterProfiler;
    /**
     * Title of the profile of the selected cluster.
     */
    private final Label profileTitleLabel;
    /**
     * Table of the most enriched fragments of the selected cluster.
     */
    private final TableView<ClusterProfiler.FeatureProfile> fragmentProfileTableView;
    /**
     * Table of the descriptors of the selected cluster that deviate most from the whole data set.
     */
    private final TableView<ClusterProfiler.FeatureProfile> descriptorProfileTableView;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, builds the view for the given result.
     *
     * @param aRunResult clustering result to display
     * @throws NullPointerException if the result is null
     */
    public ClusteringResultView(ClusteringRunResult aRunResult) throws NullPointerException {
        super();
        this.runResult = Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        ClusteringResult tmpResult = aRunResult.clusteringResult();
        //clusters ranked by size (largest first); the not clusterable molecules are listed last
        List<ClusteringRunResult.ClusterSummary> tmpRankedClusters = aRunResult.getClusterSummariesRankedBySize();
        //<editor-fold desc="summary and overflow warning">
        VBox tmpTopBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE);
        tmpTopBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        Label tmpSummaryLabel = new Label(String.format(Message.get("ClusteringResultView.summaryLabel.text"),
                aRunResult.molecules().size(),
                tmpResult.numberOfClusters(),
                tmpResult.getNumberOfUnclusteredRows(),
                aRunResult.processedMatrix().getColumnCount(),
                aRunResult.rawMatrix().getColumnCount(),
                aRunResult.parameters().algorithm().name(),
                tmpResult.vigilance(),
                tmpResult.isConverged()));
        tmpSummaryLabel.setWrapText(true);
        tmpTopBox.getChildren().add(tmpSummaryLabel);
        if (tmpResult.isOverflow()) {
            Label tmpOverflowLabel = new Label(Message.get("ClusteringResultView.overflowLabel.text"));
            tmpOverflowLabel.setWrapText(true);
            tmpOverflowLabel.setStyle("-fx-text-fill: darkred;");
            tmpTopBox.getChildren().add(tmpOverflowLabel);
        }
        this.setTop(tmpTopBox);
        //</editor-fold>
        this.clusterTableView = this.createClusterTableView(tmpRankedClusters);
        this.representativesOverviewButton = new Button(Message.get("ClusteringResultView.representativesOverviewButton.text"));
        this.representativesOverviewButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.representativesOverviewButton.tooltip")));
        this.overviewButton = new Button(Message.get("ClusteringResultView.overviewButton.text"));
        this.overviewButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.overviewButton.tooltip")));
        this.overviewButton.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        this.exportButton = new Button(Message.get("ClusteringResultView.exportButton.text"));
        this.exportButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.exportButton.tooltip")));
        this.trainTestSplitButton = new Button(Message.get("ClusteringResultView.trainTestSplitButton.text"));
        this.trainTestSplitButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.trainTestSplitButton.tooltip")));
        //the vigilance scan button is shown in the clustering bar next to the start button, see setClusteringBar()
        this.vigilanceScanButton = new Button(Message.get("ClusteringResultView.vigilanceScanButton.text"));
        this.vigilanceScanButton.setMinWidth(Region.USE_PREF_SIZE);
        this.showVigilanceScanButton = new Button(Message.get("ClusteringResultView.showVigilanceScanButton.text"));
        this.showVigilanceScanButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.showVigilanceScanButton.tooltip")));
        this.showVigilanceScanButton.setMinWidth(Region.USE_PREF_SIZE);
        //hidden until a scan of this result exists
        this.showVigilanceScanButton.setVisible(false);
        this.showVigilanceScanButton.managedProperty().bind(this.showVigilanceScanButton.visibleProperty());
        this.vigilanceScanButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.vigilanceScanButton.tooltip")));
        this.diverseSubsetButton = new Button(Message.get("ClusteringResultView.diverseSubsetButton.text"));
        this.diverseSubsetButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.diverseSubsetButton.tooltip")));
        //similarity search on the left, view and export buttons on the right
        this.similaritySearchTextField = new TextField();
        this.similaritySearchTextField.setPromptText(Message.get("ClusteringSimilaritySearch.textField.prompt"));
        this.similaritySearchTextField.setTooltip(new Tooltip(Message.get("ClusteringSimilaritySearch.textField.tooltip")));
        HBox.setHgrow(this.similaritySearchTextField, Priority.ALWAYS);
        this.similaritySearchButton = new Button(Message.get("ClusteringSimilaritySearch.button.text"));
        this.similaritySearchButton.setMinWidth(Region.USE_PREF_SIZE);
        this.similaritySearchButton.disableProperty().bind(this.similaritySearchTextField.textProperty().isEmpty());
        //enter in the text field starts the search
        this.similaritySearchTextField.setOnAction(anEvent -> this.similaritySearchButton.fire());
        this.similaritySearchMenuItem = new MenuItem(Message.get("ClusteringSimilaritySearch.contextMenuItem.text"));
        this.clusterMoleculeSelectionMenuItem = new MenuItem(Message.get("ClusterMoleculeSelection.tableMenuItem.text"));
        this.clusterMoleculeSelectionMenuItem.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        this.clusterTableView.setContextMenu(new ContextMenu(this.similaritySearchMenuItem, this.clusterMoleculeSelectionMenuItem));
        this.similaritySearchMenuItem.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        HBox tmpButtonBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE,
                new Label(Message.get("ClusteringSimilaritySearch.label.text")), this.similaritySearchTextField,
                this.similaritySearchButton, new Separator(Orientation.VERTICAL),
                this.representativesOverviewButton, this.overviewButton, this.exportButton, this.trainTestSplitButton,
                this.diverseSubsetButton);
        tmpButtonBox.setAlignment(Pos.CENTER_LEFT);
        tmpButtonBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.bottomBox = new VBox(tmpButtonBox);
        this.setBottom(this.bottomBox);
        //profile of the selected cluster on the right
        this.clusterProfiler = new ClusterProfiler(aRunResult);
        this.profileTitleLabel = new Label(Message.get("ClusterProfile.noSelection.text"));
        this.profileTitleLabel.setStyle("-fx-font-weight: bold;");
        this.profileTitleLabel.setWrapText(true);
        this.fragmentProfileTableView = this.createFragmentProfileTableView();
        this.descriptorProfileTableView = this.createDescriptorProfileTableView();
        VBox tmpProfileBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE,
                this.profileTitleLabel,
                new Label(Message.get("ClusterProfile.fragmentSection.text")), this.fragmentProfileTableView,
                new Label(Message.get("ClusterProfile.descriptorSection.text")), this.descriptorProfileTableView);
        tmpProfileBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        VBox.setVgrow(this.fragmentProfileTableView, Priority.ALWAYS);
        VBox.setVgrow(this.descriptorProfileTableView, Priority.ALWAYS);
        this.clusterTableView.getSelectionModel().selectedItemProperty().addListener(
                (anObservable, anOldValue, aNewValue) -> this.showProfile(aNewValue));
        SplitPane tmpSplitPane = new SplitPane(this.clusterTableView, tmpProfileBox);
        tmpSplitPane.setDividerPositions(0.55);
        this.setCenter(tmpSplitPane);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the displayed clustering result.
     *
     * @return clustering run result
     */
    public ClusteringRunResult getRunResult() {
        return this.runResult;
    }
    //
    /**
     * Returns the table of the clusters.
     *
     * @return cluster table view
     */
    public TableView<ClusteringRunResult.ClusterSummary> getClusterTableView() {
        return this.clusterTableView;
    }
    //
    /**
     * Sets the bar with the clustering controls (descriptor preset, settings, start and cancel), shown as the last row at
     * the bottom of the view like in the fragments and items tabs, so that one can cluster again directly from here. The
     * vigilance scan buttons of this view are inserted directly after the given start button of the bar.
     *
     * @param aClusteringBar bar with the clustering controls
     * @param aStartButton start button contained in the bar, the vigilance scan button is placed next to it
     */
    public void setClusteringBar(HBox aClusteringBar, Node aStartButton) {
        if (this.bottomBox.getChildren().size() > 1) {
            this.bottomBox.getChildren().remove(1);
        }
        int tmpStartIndex = aClusteringBar.getChildren().indexOf(aStartButton);
        aClusteringBar.getChildren().addAll(tmpStartIndex < 0 ? aClusteringBar.getChildren().size() : tmpStartIndex + 1,
                List.of(this.vigilanceScanButton, this.showVigilanceScanButton));
        this.bottomBox.getChildren().add(aClusteringBar);
    }
    //
    /**
     * Returns the text field for the name or SMILES of the query of the similarity search.
     *
     * @return similarity search text field
     */
    public TextField getSimilaritySearchTextField() {
        return this.similaritySearchTextField;
    }
    //
    /**
     * Returns the button that starts the similarity search for the query of the text field.
     *
     * @return similarity search button
     */
    public Button getSimilaritySearchButton() {
        return this.similaritySearchButton;
    }
    //
    /**
     * Returns the context menu item of the table that starts the similarity search for the representative of the
     * selected cluster.
     *
     * @return similarity search menu item
     */
    public MenuItem getSimilaritySearchMenuItem() {
        return this.similaritySearchMenuItem;
    }
    //
    /**
     * Returns the context menu item of the table that selects typical or diverse molecules from the selected cluster.
     *
     * @return cluster molecule selection menu item
     */
    public MenuItem getClusterMoleculeSelectionMenuItem() {
        return this.clusterMoleculeSelectionMenuItem;
    }
    //
    /**
     * Returns the button that opens the overview of the representatives of all clusters.
     *
     * @return representatives overview button
     */
    public Button getRepresentativesOverviewButton() {
        return this.representativesOverviewButton;
    }
    //
    /**
     * Returns the button that opens the selected cluster in the overview view.
     *
     * @return overview button
     */
    public Button getOverviewButton() {
        return this.overviewButton;
    }
    //
    /**
     * Returns the button that starts a vigilance scan with the features of this clustering.
     *
     * @return vigilance scan button
     */
    public Button getVigilanceScanButton() {
        return this.vigilanceScanButton;
    }
    //
    /**
     * Returns the button that shows the window of the last vigilance scan again; it is only visible after a scan.
     *
     * @return show vigilance scan button
     */
    public Button getShowVigilanceScanButton() {
        return this.showVigilanceScanButton;
    }
    //
    /**
     * Returns the button that selects a diverse subset.
     *
     * @return diverse subset button
     */
    public Button getDiverseSubsetButton() {
        return this.diverseSubsetButton;
    }
    //
    /**
     * Returns the button that creates a training/test partitioning based on the clusters.
     *
     * @return train/test split button
     */
    public Button getTrainTestSplitButton() {
        return this.trainTestSplitButton;
    }
    //
    /**
     * Returns the button that exports the clustering data.
     *
     * @return export button
     */
    public Button getExportButton() {
        return this.exportButton;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Shows the profile of the given cluster: the fragments that are most enriched in it and the descriptors that deviate
     * most from the whole data set. For the not clusterable molecules no profile is shown.
     *
     * @param aCluster selected cluster, may be null
     */
    private void showProfile(ClusteringRunResult.ClusterSummary aCluster) {
        if (aCluster == null || aCluster.clusterIndex() < 0) {
            this.profileTitleLabel.setText(Message.get("ClusterProfile.noSelection.text"));
            this.fragmentProfileTableView.getItems().clear();
            this.descriptorProfileTableView.getItems().clear();
            return;
        }
        ClusterProfiler.ClusterProfile tmpProfile = this.clusterProfiler.getProfile(aCluster.clusterIndex());
        this.profileTitleLabel.setText(String.format(Message.get("ClusterProfile.title.text"), aCluster.clusterIndex(), tmpProfile.size()));
        this.fragmentProfileTableView.setItems(FXCollections.observableArrayList(
                tmpProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, ClusteringResultView.NUMBER_OF_PROFILE_FEATURES, true)));
        this.descriptorProfileTableView.setItems(FXCollections.observableArrayList(
                tmpProfile.getMostDistinctFeatures(ClusterProfiler.FeatureCategory.DESCRIPTOR, ClusteringResultView.NUMBER_OF_PROFILE_FEATURES)));
    }
    //
    /**
     * Creates the table of the enriched fragments: structure, share in the cluster, share in the whole data set and the
     * difference in percentage points.
     *
     * @return fragment profile table
     */
    private TableView<ClusterProfiler.FeatureProfile> createFragmentProfileTableView() {
        TableView<ClusterProfiler.FeatureProfile> tmpTableView = new TableView<>();
        tmpTableView.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        tmpTableView.setPlaceholder(new Label(Message.get("ClusterProfile.noFragments.text")));
        TableColumn<ClusterProfiler.FeatureProfile, ImageView> tmpStructureColumn =
                new TableColumn<>(Message.get("ClusterProfile.fragmentColumn.header"));
        //the feature name of a fragment is its unique SMILES; the structure is depicted when the row becomes visible
        tmpStructureColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(
                new MoleculeDataModel(aCellData.getValue().featureName(), aCellData.getValue().featureName(), new HashMap<>()).getStructure()));
        tmpStructureColumn.setSortable(false);
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpClusterColumn =
                new TableColumn<>(Message.get("ClusterProfile.clusterShareColumn.header"));
        tmpClusterColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(ClusteringResultView.formatPercentage(aCellData.getValue().clusterValue())));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpOverallColumn =
                new TableColumn<>(Message.get("ClusterProfile.overallShareColumn.header"));
        tmpOverallColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(ClusteringResultView.formatPercentage(aCellData.getValue().overallValue())));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpDifferenceColumn =
                new TableColumn<>(Message.get("ClusterProfile.differenceColumn.header"));
        tmpDifferenceColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(
                String.format(Locale.getDefault(), "%+.0f", 100.0 * aCellData.getValue().score())));
        tmpTableView.getColumns().addAll(List.of(tmpStructureColumn, tmpClusterColumn, tmpOverallColumn, tmpDifferenceColumn));
        return tmpTableView;
    }
    //
    /**
     * Creates the table of the most distinct descriptors: name, mean in the cluster, mean in the whole data set and the
     * deviation in standard deviations of the whole data set.
     *
     * @return descriptor profile table
     */
    private TableView<ClusterProfiler.FeatureProfile> createDescriptorProfileTableView() {
        TableView<ClusterProfiler.FeatureProfile> tmpTableView = new TableView<>();
        tmpTableView.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        tmpTableView.setPlaceholder(new Label(Message.get("ClusterProfile.noDescriptors.text")));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpNameColumn =
                new TableColumn<>(Message.get("ClusterProfile.descriptorColumn.header"));
        tmpNameColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(aCellData.getValue().featureName()));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpClusterColumn =
                new TableColumn<>(Message.get("ClusterProfile.clusterMeanColumn.header"));
        tmpClusterColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(ClusteringResultView.formatNumber(aCellData.getValue().clusterValue())));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpOverallColumn =
                new TableColumn<>(Message.get("ClusterProfile.overallMeanColumn.header"));
        tmpOverallColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(ClusteringResultView.formatNumber(aCellData.getValue().overallValue())));
        TableColumn<ClusterProfiler.FeatureProfile, String> tmpDeviationColumn =
                new TableColumn<>(Message.get("ClusterProfile.deviationColumn.header"));
        tmpDeviationColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(
                String.format(Locale.getDefault(), "%+.2f", aCellData.getValue().score())));
        tmpTableView.getColumns().addAll(List.of(tmpNameColumn, tmpClusterColumn, tmpOverallColumn, tmpDeviationColumn));
        return tmpTableView;
    }
    //
    /**
     * Formats a share as percentage with one decimal.
     *
     * @param aShare share in [0,1], may be NaN
     * @return formatted percentage, "-" for NaN
     */
    private static String formatPercentage(double aShare) {
        return Double.isNaN(aShare) ? "-" : String.format(Locale.getDefault(), "%.1f %%", 100.0 * aShare);
    }
    //
    /**
     * Formats a descriptor value with two decimals.
     *
     * @param aValue value, may be NaN
     * @return formatted value, "-" for NaN
     */
    private static String formatNumber(double aValue) {
        return Double.isNaN(aValue) ? "-" : String.format(Locale.getDefault(), "%.2f", aValue);
    }
    //
    /**
     * Creates the table of the clusters in the given order with the columns rank, cluster, size, representative name and
     * representative structure.
     *
     * @param aRankedClusters clusters ranked by size
     * @return cluster table view
     */
    private TableView<ClusteringRunResult.ClusterSummary> createClusterTableView(List<ClusteringRunResult.ClusterSummary> aRankedClusters) {
        TableView<ClusteringRunResult.ClusterSummary> tmpTableView = new TableView<>();
        tmpTableView.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        TableColumn<ClusteringRunResult.ClusterSummary, Integer> tmpRankColumn =
                new TableColumn<>(Message.get("ClusteringResultView.rankColumn.header"));
        //the rank is the position in the size ranking, the not clusterable molecules have no rank
        tmpRankColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().clusterIndex() < 0
                ? null : aRankedClusters.indexOf(aCellData.getValue()) + 1));
        TableColumn<ClusteringRunResult.ClusterSummary, String> tmpClusterColumn =
                new TableColumn<>(Message.get("ClusteringResultView.clusterColumn.header"));
        tmpClusterColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(
                aCellData.getValue().clusterIndex() < 0
                        ? Message.get("ClusteringResultView.notClusterable.text")
                        : Integer.toString(aCellData.getValue().clusterIndex())));
        TableColumn<ClusteringRunResult.ClusterSummary, Integer> tmpSizeColumn =
                new TableColumn<>(Message.get("ClusteringResultView.sizeColumn.header"));
        tmpSizeColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().getSize()));
        TableColumn<ClusteringRunResult.ClusterSummary, String> tmpNameColumn =
                new TableColumn<>(Message.get("ClusteringResultView.representativeNameColumn.header"));
        tmpNameColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(
                aCellData.getValue().representative() == null ? "" : aCellData.getValue().representative().getName()));
        TableColumn<ClusteringRunResult.ClusterSummary, ImageView> tmpStructureColumn =
                new TableColumn<>(Message.get("ClusteringResultView.representativeStructureColumn.header"));
        //structures are depicted lazily when a row becomes visible
        tmpStructureColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(
                aCellData.getValue().representative() == null ? null : aCellData.getValue().representative().getStructure()));
        tmpStructureColumn.setSortable(false);
        tmpTableView.getColumns().addAll(List.of(tmpRankColumn, tmpClusterColumn, tmpSizeColumn, tmpNameColumn, tmpStructureColumn));
        tmpTableView.setItems(FXCollections.observableArrayList(aRankedClusters));
        return tmpTableView;
    }
    //</editor-fold>
}
