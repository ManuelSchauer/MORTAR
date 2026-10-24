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
import de.unijena.cheminf.mortar.model.clustering.ClusterChemistryCheck;
import de.unijena.cheminf.mortar.model.clustering.ClusterDescriptorStatistics;
import de.unijena.cheminf.mortar.model.clustering.ClusterProfiler;
import de.unijena.cheminf.mortar.model.clustering.ClusterQuality;
import de.unijena.cheminf.mortar.model.clustering.ClusterStability;
import de.unijena.cheminf.mortar.model.clustering.ClusteringResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.depict.DepictionUtil;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import javafx.beans.Observable;
import javafx.beans.binding.Bindings;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ObjectProperty;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleObjectProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.collections.transformation.SortedList;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableRow;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;

import org.openscience.cdk.exception.CDKException;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * View that shows the result of a clustering run: a summary, a warning in case of a cluster overflow and a table with
 * one row per cluster ranked from the largest to the smallest cluster (rank, cluster index, size, name and structure of
 * the representative, median of the chosen descriptors). Clusters can be marked by check box or by the range filter on
 * the descriptor medians; marked clusters are listed first and, if any are marked, only they are used for the export,
 * the diverse subset and the training/test partitioning. The buttons open the overview of all cluster representatives,
 * open the selected cluster in the overview and export the clustering data. A similarity search finds the molecules most
 * similar to a molecule given by name or SMILES (search field) or to a cluster representative (context menu of the
 * table). On the right, the profile of the selected cluster shows its most enriched fragments and its most distinct
 * descriptors. The actions of the buttons are set by the controller.
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
    /**
     * Fixed row height of the descriptor table of the profile.
     */
    private static final double DESCRIPTOR_ROW_HEIGHT = 22.0;
    /**
     * Height of the header of the descriptor table of the profile, including the borders.
     */
    private static final double DESCRIPTOR_TABLE_HEADER_HEIGHT = 30.0;
    /**
     * Background style of the rows of marked clusters, a light blue like the bars of the cluster overview.
     */
    private static final String MARKED_ROW_STYLE = "-fx-background-color: #e3f1f8;";
    /**
     * Width of the structure columns, the default width of the structure images plus a margin.
     */
    public static final double STRUCTURE_COLUMN_WIDTH = BasicDefinitions.DEFAULT_IMAGE_WIDTH_DEFAULT + 10.0;
    /**
     * Minimum width of a descriptor column; wider if its header needs more space.
     */
    private static final double DESCRIPTOR_COLUMN_MIN_WIDTH = 70.0;
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusteringResultView.class.getName());
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
     * Button to select a diverse subset (molecules per cluster).
     */
    private final Button diverseSubsetButton;
    /**
     * Button to open the interactive PCA of this clustering in a new tab.
     */
    private final Button pcaButton;
    /**
     * Button to open the quality of all clusters and the distances between them.
     */
    private final Button clusterQualityButton;
    /**
     * Quality line of the selected cluster below the profile title: silhouette, neighbour purity and nearest cluster.
     */
    private final Label profileQualityLabel;
    /**
     * Quality of the clusters, set when its calculation in the background has finished; null before.
     */
    private final ObjectProperty<ClusterQuality> clusterQualityProperty;
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
     * Context menu item of the table that shows the descriptor values of the molecules of the selected cluster.
     */
    private final MenuItem descriptorValuesMenuItem;
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
    /**
     * Medians of the descriptors per cluster, from the raw (untransformed) values.
     */
    private final ClusterDescriptorStatistics descriptorStatistics;
    /**
     * Marked state of every cluster, by cluster index; the not clusterable molecules cannot be marked.
     */
    private final Map<Integer, BooleanProperty> markedPropertyOfClusterMap;
    /**
     * Indices of the marked clusters in ascending order.
     */
    private final Set<Integer> markedClusterIndices;
    /**
     * Label that tells how many clusters are marked and used.
     */
    private final Label markedStatusLabel;
    /**
     * Column group of the cluster table with one median column per displayed descriptor.
     */
    private final TableColumn<ClusteringRunResult.ClusterSummary, Double> descriptorMedianColumnGroup;
    /**
     * Button to choose the displayed descriptors.
     */
    private final Button chooseDescriptorsButton;
    /**
     * Button to mark the clusters by the range filter on the descriptor medians.
     */
    private final Button rangeFilterButton;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Displayed descriptor columns that exist in this result, in display order.
     */
    private List<String> displayedDescriptorNames;
    /**
     * Conditions of the last applied range filter, proposed when the filter is opened again.
     */
    private List<ClusterDescriptorStatistics.RangeCondition> lastRangeConditions;
    /**
     * Result of the optional stability check, null if it was not run.
     */
    private ClusterStability.StabilityResult stabilityResult;
    /**
     * Result of the optional chemical check, null if it was not run.
     */
    private ClusterChemistryCheck.ChemistryCheckResult chemistryCheckResult;
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
        //marks and descriptor medians of the cluster table
        this.descriptorStatistics = new ClusterDescriptorStatistics(aRunResult);
        this.displayedDescriptorNames = List.of();
        this.lastRangeConditions = List.of();
        this.markedClusterIndices = new TreeSet<>();
        this.markedPropertyOfClusterMap = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(tmpResult.numberOfClusters()));
        this.markedStatusLabel = new Label();
        for (int k = 0; k < tmpResult.numberOfClusters(); k++) {
            int tmpClusterIndex = k;
            BooleanProperty tmpMarkedProperty = new SimpleBooleanProperty(false);
            tmpMarkedProperty.addListener((anObservable, anOldValue, aNewValue) -> {
                if (aNewValue) {
                    this.markedClusterIndices.add(tmpClusterIndex);
                } else {
                    this.markedClusterIndices.remove(tmpClusterIndex);
                }
                this.updateMarkedStatusLabel();
            });
            this.markedPropertyOfClusterMap.put(k, tmpMarkedProperty);
        }
        this.updateMarkedStatusLabel();
        this.descriptorMedianColumnGroup = new TableColumn<>(Message.get("ClusteringResultView.medianColumnGroup.header"));
        this.clusterTableView = this.createClusterTableView(tmpRankedClusters);
        //marking controls above the cluster table
        Button tmpMarkAllButton = new Button(Message.get("ClusteringResultView.markAllButton.text"));
        tmpMarkAllButton.setOnAction(anEvent -> this.markedPropertyOfClusterMap.values().forEach(aProperty -> aProperty.set(true)));
        Button tmpMarkNoneButton = new Button(Message.get("ClusteringResultView.markNoneButton.text"));
        tmpMarkNoneButton.setOnAction(anEvent -> this.markedPropertyOfClusterMap.values().forEach(aProperty -> aProperty.set(false)));
        Button tmpInvertMarksButton = new Button(Message.get("ClusteringResultView.invertMarksButton.text"));
        tmpInvertMarksButton.setOnAction(anEvent -> this.markedPropertyOfClusterMap.values().forEach(aProperty -> aProperty.set(!aProperty.get())));
        this.rangeFilterButton = new Button(Message.get("ClusteringResultView.rangeFilterButton.text"));
        this.rangeFilterButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.rangeFilterButton.tooltip")));
        this.chooseDescriptorsButton = new Button(Message.get("ClusteringResultView.chooseDescriptorsButton.text"));
        this.chooseDescriptorsButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.chooseDescriptorsButton.tooltip")));
        //no descriptors, no filter
        boolean tmpHasNoDescriptors = this.descriptorStatistics.getDescriptorColumnNames().isEmpty();
        this.rangeFilterButton.setDisable(tmpHasNoDescriptors);
        this.chooseDescriptorsButton.setDisable(tmpHasNoDescriptors);
        for (Button tmpButton : List.of(tmpMarkAllButton, tmpMarkNoneButton, tmpInvertMarksButton, this.rangeFilterButton, this.chooseDescriptorsButton)) {
            tmpButton.setMinWidth(Region.USE_PREF_SIZE);
        }
        tmpMarkAllButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.markAllButton.tooltip")));
        HBox tmpMarkBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, new Label(Message.get("ClusteringResultView.markLabel.text")),
                tmpMarkAllButton, tmpMarkNoneButton, tmpInvertMarksButton, this.rangeFilterButton,
                new Separator(Orientation.VERTICAL), this.chooseDescriptorsButton, this.markedStatusLabel);
        tmpMarkBox.setAlignment(Pos.CENTER_LEFT);
        tmpMarkBox.setPadding(new Insets(0.0, 0.0, GuiDefinitions.GUI_INSETS_VALUE, 0.0));
        BorderPane tmpClusterPane = new BorderPane(this.clusterTableView);
        tmpClusterPane.setTop(tmpMarkBox);
        this.representativesOverviewButton = new Button(Message.get("ClusteringResultView.representativesOverviewButton.text"));
        this.representativesOverviewButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.representativesOverviewButton.tooltip")));
        this.overviewButton = new Button(Message.get("ClusteringResultView.overviewButton.text"));
        this.overviewButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.overviewButton.tooltip")));
        this.overviewButton.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        this.exportButton = new Button(Message.get("ClusteringResultView.exportButton.text"));
        this.exportButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.exportButton.tooltip")));
        this.trainTestSplitButton = new Button(Message.get("ClusteringResultView.trainTestSplitButton.text"));
        this.trainTestSplitButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.trainTestSplitButton.tooltip")));
        this.diverseSubsetButton = new Button(Message.get("ClusteringResultView.diverseSubsetButton.text"));
        this.diverseSubsetButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.diverseSubsetButton.tooltip")));
        this.pcaButton = new Button(Message.get("ClusteringResultView.pcaButton.text"));
        this.pcaButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.pcaButton.tooltip")));
        this.pcaButton.setMinWidth(Region.USE_PREF_SIZE);
        this.clusterQualityButton = new Button(Message.get("ClusteringResultView.clusterQualityButton.text"));
        this.clusterQualityButton.setTooltip(new Tooltip(Message.get("ClusteringResultView.clusterQualityButton.tooltip")));
        this.clusterQualityButton.setMinWidth(Region.USE_PREF_SIZE);
        //the quality is calculated in the background, the button is enabled when it is available
        this.clusterQualityProperty = new SimpleObjectProperty<>(null);
        this.clusterQualityButton.disableProperty().bind(this.clusterQualityProperty.isNull());
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
        this.descriptorValuesMenuItem = new MenuItem(Message.get("ClusteringResultView.descriptorValuesMenuItem.text"));
        this.descriptorValuesMenuItem.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        //the descriptor values are set apart from the other entries
        this.clusterTableView.setContextMenu(new ContextMenu(this.similaritySearchMenuItem, this.clusterMoleculeSelectionMenuItem,
                new SeparatorMenuItem(), this.descriptorValuesMenuItem));
        this.similaritySearchMenuItem.disableProperty().bind(this.clusterTableView.getSelectionModel().selectedItemProperty().isNull());
        HBox tmpButtonBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE,
                new Label(Message.get("ClusteringSimilaritySearch.label.text")), this.similaritySearchTextField,
                this.similaritySearchButton, new Separator(Orientation.VERTICAL),
                this.representativesOverviewButton, this.overviewButton, this.exportButton, this.trainTestSplitButton,
                this.diverseSubsetButton, new Separator(Orientation.VERTICAL), this.clusterQualityButton, this.pcaButton);
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
        this.profileQualityLabel = new Label();
        this.profileQualityLabel.setWrapText(true);
        this.profileQualityLabel.setTooltip(new Tooltip(Message.get("ClusterProfile.quality.tooltip")));
        VBox tmpProfileBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE,
                this.profileTitleLabel, this.profileQualityLabel,
                new Label(Message.get("ClusterProfile.fragmentSection.text")), this.fragmentProfileTableView,
                new Label(Message.get("ClusterProfile.descriptorSection.text")), this.descriptorProfileTableView);
        tmpProfileBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        //the descriptor table is as high as the shown descriptors need, the fragments get the remaining space
        this.descriptorProfileTableView.setFixedCellSize(ClusteringResultView.DESCRIPTOR_ROW_HEIGHT);
        double tmpDescriptorTableHeight = ClusteringResultView.DESCRIPTOR_TABLE_HEADER_HEIGHT
                + ClusteringResultView.NUMBER_OF_PROFILE_FEATURES * ClusteringResultView.DESCRIPTOR_ROW_HEIGHT;
        this.descriptorProfileTableView.setPrefHeight(tmpDescriptorTableHeight);
        this.descriptorProfileTableView.setMinHeight(tmpDescriptorTableHeight);
        this.descriptorProfileTableView.setMaxHeight(tmpDescriptorTableHeight);
        VBox.setVgrow(this.fragmentProfileTableView, Priority.ALWAYS);
        VBox.setVgrow(this.descriptorProfileTableView, Priority.NEVER);
        this.clusterTableView.getSelectionModel().selectedItemProperty().addListener(
                (anObservable, anOldValue, aNewValue) -> this.showProfile(aNewValue));
        tmpClusterPane.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        SplitPane tmpSplitPane = new SplitPane(tmpClusterPane, tmpProfileBox);
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
     * Returns the medians of the descriptors per cluster.
     *
     * @return descriptor statistics of this result
     */
    public ClusterDescriptorStatistics getDescriptorStatistics() {
        return this.descriptorStatistics;
    }
    //
    /**
     * Returns the indices of the marked clusters; if none is marked, all clusters are used.
     *
     * @return unmodifiable copy of the marked cluster indices in ascending order, empty if none is marked
     */
    public Set<Integer> getMarkedClusterIndices() {
        return Collections.unmodifiableSet(new TreeSet<>(this.markedClusterIndices));
    }
    //
    /**
     * Marks exactly the given clusters, e.g. the ones matching the range filter; they are listed first.
     *
     * @param aClusterIndices indices of the clusters to mark
     * @throws NullPointerException if the set is null
     */
    public void setMarkedClusterIndices(Set<Integer> aClusterIndices) throws NullPointerException {
        Objects.requireNonNull(aClusterIndices, "Given cluster index set is null.");
        this.markedPropertyOfClusterMap.forEach((aClusterIndex, aProperty) -> aProperty.set(aClusterIndices.contains(aClusterIndex)));
        this.clusterTableView.scrollTo(0);
    }
    //
    /**
     * Sets the descriptors whose medians are shown in the cluster table; names of descriptors that were
     * not calculated in this result are ignored.
     *
     * @param aNameList descriptor column names in display order
     * @throws NullPointerException if the list is null
     */
    public void setDisplayedDescriptorNames(List<String> aNameList) throws NullPointerException {
        Objects.requireNonNull(aNameList, "Given name list is null.");
        this.displayedDescriptorNames = aNameList.stream().filter(this.descriptorStatistics::hasDescriptorColumn).toList();
        this.descriptorMedianColumnGroup.getColumns().clear();
        for (String tmpName : this.displayedDescriptorNames) {
            TableColumn<ClusteringRunResult.ClusterSummary, Double> tmpColumn = new TableColumn<>(tmpName);
            ClusteringResultView.setFixedWidthForHeader(tmpColumn, tmpName);
            //null instead of NaN, so that clusters without value are sorted consistently
            tmpColumn.setCellValueFactory(aCellData -> {
                double tmpMedian = this.descriptorStatistics.getMedian(aCellData.getValue().clusterIndex(), tmpName);
                return new ReadOnlyObjectWrapper<>(Double.isNaN(tmpMedian) ? null : tmpMedian);
            });
            tmpColumn.setCellFactory(aColumn -> new TableCell<>() {
                @Override
                protected void updateItem(Double aValue, boolean anIsEmpty) {
                    super.updateItem(aValue, anIsEmpty);
                    this.setText(anIsEmpty || aValue == null ? null : ClusteringResultView.formatNumber(aValue));
                }
            });
            this.descriptorMedianColumnGroup.getColumns().add(tmpColumn);
        }
        //the group without columns is hidden
        if (this.displayedDescriptorNames.isEmpty()) {
            this.clusterTableView.getColumns().remove(this.descriptorMedianColumnGroup);
        } else if (!this.clusterTableView.getColumns().contains(this.descriptorMedianColumnGroup)) {
            this.clusterTableView.getColumns().add(this.descriptorMedianColumnGroup);
        }
    }
    //
    /**
     * Returns the descriptors whose medians are shown, i.e. the chosen ones that exist in this result.
     *
     * @return descriptor column names in display order
     */
    public List<String> getDisplayedDescriptorNames() {
        return this.displayedDescriptorNames;
    }
    //
    /**
     * Returns the conditions of the last applied range filter.
     *
     * @return range conditions, empty if the filter was not applied yet
     */
    public List<ClusterDescriptorStatistics.RangeCondition> getLastRangeConditions() {
        return this.lastRangeConditions;
    }
    //
    /**
     * Sets the conditions of the last applied range filter, which are proposed when the filter is opened again.
     *
     * @param aConditionList range conditions
     * @throws NullPointerException if the list is null
     */
    public void setLastRangeConditions(List<ClusterDescriptorStatistics.RangeCondition> aConditionList) throws NullPointerException {
        this.lastRangeConditions = List.copyOf(Objects.requireNonNull(aConditionList, "Given condition list is null."));
    }
    //
    /**
     * Returns the button that opens the choice of the displayed descriptors.
     *
     * @return choose descriptors button
     */
    public Button getChooseDescriptorsButton() {
        return this.chooseDescriptorsButton;
    }
    //
    /**
     * Returns the button that opens the range filter on the descriptor medians.
     *
     * @return range filter button
     */
    public Button getRangeFilterButton() {
        return this.rangeFilterButton;
    }
    //
    /**
     * Sets the bar with the clustering controls (descriptor preset, settings, vigilance, start, vigilance scan and
     * cancel), shown as the last row at the bottom of the view like in the fragments and items tabs, so that one can
     * cluster or scan again directly from here.
     *
     * @param aClusteringBar bar with the clustering controls
     */
    public void setClusteringBar(HBox aClusteringBar) {
        if (this.bottomBox.getChildren().size() > 1) {
            this.bottomBox.getChildren().remove(1);
        }
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
     * Returns the context menu item of the table that shows the descriptor values of the molecules of the selected
     * cluster.
     *
     * @return descriptor values menu item
     */
    public MenuItem getDescriptorValuesMenuItem() {
        return this.descriptorValuesMenuItem;
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
     * Returns the quality of the clusters and the distances between them.
     *
     * @return cluster quality of this result, null while it is calculated in the background
     */
    public ClusterQuality getClusterQuality() {
        return this.clusterQualityProperty.get();
    }
    //
    /**
     * Sets the quality of the clusters when its calculation has finished; the profile of the selected cluster is updated.
     *
     * @param aClusterQuality cluster quality of this result
     */
    public void setClusterQuality(ClusterQuality aClusterQuality) {
        this.clusterQualityProperty.set(aClusterQuality);
        this.showProfile(this.clusterTableView.getSelectionModel().getSelectedItem());
    }
    //
    /**
     * Returns the result of the optional stability check.
     *
     * @return stability result, null if the check was not run
     */
    public ClusterStability.StabilityResult getStabilityResult() {
        return this.stabilityResult;
    }
    //
    /**
     * Sets the result of the optional stability check, so that it is shown again when the quality window is reopened.
     *
     * @param aStabilityResult stability result
     */
    public void setStabilityResult(ClusterStability.StabilityResult aStabilityResult) {
        this.stabilityResult = aStabilityResult;
    }
    //
    /**
     * Returns the result of the optional chemical check with an independent fingerprint and the scaffolds.
     *
     * @return chemical check result, null if the check was not run
     */
    public ClusterChemistryCheck.ChemistryCheckResult getChemistryCheckResult() {
        return this.chemistryCheckResult;
    }
    //
    /**
     * Sets the result of the optional chemical check, so that it is shown again when the quality window is reopened.
     *
     * @param aChemistryCheckResult chemical check result
     */
    public void setChemistryCheckResult(ClusterChemistryCheck.ChemistryCheckResult aChemistryCheckResult) {
        this.chemistryCheckResult = aChemistryCheckResult;
    }
    //
    /**
     * Returns the button that opens the quality of all clusters and the distances between them.
     *
     * @return cluster quality button
     */
    public Button getClusterQualityButton() {
        return this.clusterQualityButton;
    }
    //
    /**
     * Returns the button that opens the interactive PCA of this clustering.
     *
     * @return PCA button
     */
    public Button getPcaButton() {
        return this.pcaButton;
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
            this.profileQualityLabel.setText("");
            this.fragmentProfileTableView.getItems().clear();
            this.descriptorProfileTableView.getItems().clear();
            return;
        }
        ClusterProfiler.ClusterProfile tmpProfile = this.clusterProfiler.getProfile(aCluster.clusterIndex());
        this.profileTitleLabel.setText(String.format(Message.get("ClusterProfile.title.text"), aCluster.clusterIndex(), tmpProfile.size()));
        ClusterQuality tmpClusterQuality = this.clusterQualityProperty.get();
        if (tmpClusterQuality == null) {
            this.profileQualityLabel.setText(Message.get("ClusterProfile.quality.calculating.text"));
        } else {
            ClusterQuality.ClusterQualityEntry tmpQuality = tmpClusterQuality.getEntry(aCluster.clusterIndex());
            if (tmpQuality.nearestClusterIndex() < 0) {
                this.profileQualityLabel.setText(String.format(Locale.getDefault(), Message.get("ClusterProfile.quality.single.text"),
                        tmpQuality.meanSimilarityToOwnCentre()));
            } else if (tmpQuality.numberOfNeighbours() == 0) {
                //a single molecule has no neighbours in its cluster and no distance within the cluster
                this.profileQualityLabel.setText(String.format(Locale.getDefault(), Message.get("ClusterProfile.quality.singleton.text"),
                        tmpQuality.nearestClusterIndex(), tmpQuality.meanDistanceToNearestCluster()));
            } else {
                this.profileQualityLabel.setText(String.format(Locale.getDefault(), Message.get("ClusterProfile.quality.text"), tmpQuality.silhouette(),
                            100.0 * tmpQuality.neighbourPurity(), tmpQuality.numberOfNegativeSilhouettes(), tmpQuality.nearestClusterIndex(),
                            tmpQuality.meanDistanceToNearestCluster(), tmpQuality.meanDistanceWithinCluster(), tmpQuality.numberOfNeighbours()));
            }
        }
        this.fragmentProfileTableView.setItems(FXCollections.observableArrayList(
                tmpProfile.getTopFeatures(ClusterProfiler.FeatureCategory.FRAGMENT, ClusteringResultView.NUMBER_OF_PROFILE_FEATURES, true)));
        this.descriptorProfileTableView.setItems(FXCollections.observableArrayList(
                tmpProfile.getMostDistinctFeatures(ClusterProfiler.FeatureCategory.DESCRIPTOR, ClusteringResultView.NUMBER_OF_PROFILE_FEATURES)));
    }
    //
    /**
     * Shows how many clusters are marked and whether all clusters or only the marked ones are used.
     */
    private void updateMarkedStatusLabel() {
        int tmpNumberOfClusters = this.markedPropertyOfClusterMap.size();
        this.markedStatusLabel.setText(this.markedClusterIndices.isEmpty()
                ? String.format(Message.get("ClusteringResultView.markedStatus.none.text"), tmpNumberOfClusters)
                : String.format(Message.get("ClusteringResultView.markedStatus.some.text"), this.markedClusterIndices.size(), tmpNumberOfClusters));
    }
    //
    /**
     * Returns whether the given cluster is marked; the not clusterable molecules never are.
     *
     * @param aCluster cluster
     * @return true if marked
     */
    private boolean isMarked(ClusteringRunResult.ClusterSummary aCluster) {
        return aCluster != null && this.markedClusterIndices.contains(aCluster.clusterIndex());
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
        //numeric columns, so that they sort by value; shares and difference in percent
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpClusterColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.clusterShareColumn.header"), aFeature -> 100.0 * aFeature.clusterValue(), "%.1f %%");
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpOverallColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.overallShareColumn.header"), aFeature -> 100.0 * aFeature.overallValue(), "%.1f %%");
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpDifferenceColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.differenceColumn.header"), aFeature -> 100.0 * aFeature.score(), "%+.0f");
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
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpClusterColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.clusterMeanColumn.header"), ClusterProfiler.FeatureProfile::clusterValue, "%.2f");
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpOverallColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.overallMeanColumn.header"), ClusterProfiler.FeatureProfile::overallValue, "%.2f");
        //shown with sign, sorted by absolute value: the most deviating descriptors first, in either direction
        TableColumn<ClusterProfiler.FeatureProfile, Double> tmpDeviationColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterProfile.deviationColumn.header"), ClusterProfiler.FeatureProfile::score, "%+.2f");
        tmpDeviationColumn.setComparator(Comparator.nullsFirst(Comparator.comparingDouble(Math::abs)));
        tmpDeviationColumn.setSortType(TableColumn.SortType.DESCENDING);
        tmpTableView.getColumns().addAll(List.of(tmpNameColumn, tmpClusterColumn, tmpOverallColumn, tmpDeviationColumn));
        tmpTableView.getSortOrder().add(tmpDeviationColumn);
        //new items are sorted again by the current sort order, by default by absolute deviation
        tmpTableView.itemsProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpTableView.sort());
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
     * Returns the number of a text of the cluster column for sorting: the cluster index, or the largest integer for the
     * text of the not clusterable molecules, so that they are sorted after all clusters in ascending order.
     *
     * @param aText text of the cluster column
     * @return cluster index or {@link Integer#MAX_VALUE}
     */
    private static int getClusterNumberForSorting(String aText) {
        try {
            return Integer.parseInt(aText);
        } catch (NumberFormatException anException) {
            return Integer.MAX_VALUE;
        }
    }
    //
    /**
     * Creates a numeric column, which sorts by value and shows the value in the given format; NaN is shown as "-" and
     * sorted like a missing value.
     *
     * @param aHeader column header
     * @param aValueFunction function returning the value of a row
     * @param aFormat format of the value, e.g. "%.2f"
     * @param <T> type of the rows
     * @return table column
     */
    public static <T> TableColumn<T, Double> createNumberColumn(String aHeader, ToDoubleFunction<T> aValueFunction, String aFormat) {
        TableColumn<T, Double> tmpColumn = new TableColumn<>(aHeader);
        tmpColumn.setCellValueFactory(aCellData -> {
            double tmpValue = aValueFunction.applyAsDouble(aCellData.getValue());
            return new ReadOnlyObjectWrapper<>(Double.isNaN(tmpValue) ? null : tmpValue);
        });
        tmpColumn.setCellFactory(aColumn -> new TableCell<>() {
            @Override
            protected void updateItem(Double aValue, boolean anIsEmpty) {
                super.updateItem(aValue, anIsEmpty);
                this.setText(anIsEmpty ? null : aValue == null ? "-" : String.format(Locale.getDefault(), aFormat, aValue));
            }
        });
        return tmpColumn;
    }
    //
    /**
     * Creates a cell that depicts the structure of its molecule when it becomes visible, centred and always in MORTAR's
     * default image size. The image size stored in a molecule is set by the other tables of MORTAR (e.g. very wide and
     * low in the molecules tab) and is therefore not used here.
     *
     * @param <S> type of the table rows
     * @return structure cell
     */
    public static <S> TableCell<S, MoleculeDataModel> createStructureCell() {
        TableCell<S, MoleculeDataModel> tmpCell = new TableCell<>() {
            @Override
            protected void updateItem(MoleculeDataModel aMolecule, boolean anIsEmpty) {
                super.updateItem(aMolecule, anIsEmpty);
                this.setText(null);
                if (anIsEmpty || aMolecule == null) {
                    this.setGraphic(null);
                    return;
                }
                Image tmpImage;
                try {
                    tmpImage = DepictionUtil.depictImageWithTransparentBackground(aMolecule.getAtomContainer(), 1.0,
                            BasicDefinitions.DEFAULT_IMAGE_WIDTH_DEFAULT, BasicDefinitions.DEFAULT_IMAGE_HEIGHT_DEFAULT, false);
                } catch (CDKException anException) {
                    ClusteringResultView.LOGGER.log(Level.WARNING, anException.toString(), anException);
                    tmpImage = DepictionUtil.depictErrorImage(anException.getMessage(),
                            (int) BasicDefinitions.DEFAULT_IMAGE_WIDTH_DEFAULT, (int) BasicDefinitions.DEFAULT_IMAGE_HEIGHT_DEFAULT);
                }
                this.setGraphic(new ImageView(tmpImage));
            }
        };
        tmpCell.setAlignment(Pos.CENTER);
        return tmpCell;
    }
    //
    /**
     * Gives the column a fixed width that shows the whole given header text, at least the minimum descriptor column width.
     *
     * @param aColumn table column
     * @param aHeader header text of the column
     */
    public static void setFixedWidthForHeader(TableColumn<?, ?> aColumn, String aHeader) {
        double tmpTextWidth = new Text(aHeader).getLayoutBounds().getWidth();
        //margin for the cell padding and the sort arrow
        aColumn.setPrefWidth(Math.max(ClusteringResultView.DESCRIPTOR_COLUMN_MIN_WIDTH, Math.ceil(tmpTextWidth) + 30.0));
    }
    //
    /**
     * Creates the table of the clusters with the columns mark, rank, cluster, size, representative name, representative
     * structure and the medians of the displayed descriptors (see {@link #setDisplayedDescriptorNames(List)}). The marked
     * clusters are always listed first (highlighted), within both groups the clusters are in the given order or sorted by
     * the column chosen by the user.
     *
     * @param aRankedClusters clusters ranked by size
     * @return cluster table view
     */
    private TableView<ClusteringRunResult.ClusterSummary> createClusterTableView(List<ClusteringRunResult.ClusterSummary> aRankedClusters) {
        TableView<ClusteringRunResult.ClusterSummary> tmpTableView = new TableView<>();
        //fixed column widths, with many descriptor columns the table scrolls horizontally
        tmpTableView.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
        TableColumn<ClusteringRunResult.ClusterSummary, Boolean> tmpMarkColumn = new TableColumn<>();
        tmpMarkColumn.setGraphic(new Label("\u2713"));
        tmpMarkColumn.setSortable(false);
        tmpMarkColumn.setResizable(false);
        tmpMarkColumn.setPrefWidth(30.0);
        tmpMarkColumn.setCellFactory(aColumn -> this.createMarkCell());
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
        //numeric order instead of text order (2 before 11), the not clusterable molecules after the clusters
        tmpClusterColumn.setComparator(Comparator.comparingInt(ClusteringResultView::getClusterNumberForSorting));
        TableColumn<ClusteringRunResult.ClusterSummary, Integer> tmpSizeColumn =
                new TableColumn<>(Message.get("ClusteringResultView.sizeColumn.header"));
        tmpSizeColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().getSize()));
        TableColumn<ClusteringRunResult.ClusterSummary, String> tmpNameColumn =
                new TableColumn<>(Message.get("ClusteringResultView.representativeNameColumn.header"));
        tmpNameColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(
                aCellData.getValue().representative() == null ? "" : aCellData.getValue().representative().getName()));
        TableColumn<ClusteringRunResult.ClusterSummary, MoleculeDataModel> tmpStructureColumn =
                new TableColumn<>(Message.get("ClusteringResultView.representativeStructureColumn.header"));
        tmpStructureColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().representative()));
        tmpStructureColumn.setCellFactory(aColumn -> ClusteringResultView.createStructureCell());
        tmpStructureColumn.setSortable(false);
        tmpStructureColumn.setPrefWidth(ClusteringResultView.STRUCTURE_COLUMN_WIDTH);
        tmpRankColumn.setPrefWidth(50.0);
        tmpClusterColumn.setPrefWidth(110.0);
        tmpSizeColumn.setPrefWidth(60.0);
        tmpNameColumn.setPrefWidth(150.0);
        tmpTableView.getColumns().addAll(List.of(tmpMarkColumn, tmpRankColumn, tmpClusterColumn, tmpSizeColumn, tmpNameColumn, tmpStructureColumn));
        //a change of a mark updates the list, so that the sorted list moves the cluster into its group
        ObservableList<ClusteringRunResult.ClusterSummary> tmpItems = FXCollections.observableArrayList(aCluster ->
                aCluster.clusterIndex() < 0 ? new Observable[0] : new Observable[] {this.markedPropertyOfClusterMap.get(aCluster.clusterIndex())});
        tmpItems.addAll(aRankedClusters);
        SortedList<ClusteringRunResult.ClusterSummary> tmpSortedItems = new SortedList<>(tmpItems);
        //marked first, then the column chosen by the user or the size ranking
        Comparator<ClusteringRunResult.ClusterSummary> tmpRankingComparator = Comparator.comparingInt(aRankedClusters::indexOf);
        tmpSortedItems.comparatorProperty().bind(Bindings.createObjectBinding(() -> {
            Comparator<ClusteringRunResult.ClusterSummary> tmpTableComparator = tmpTableView.getComparator();
            return Comparator.comparing((ClusteringRunResult.ClusterSummary aCluster) -> !this.isMarked(aCluster))
                    .thenComparing(tmpTableComparator == null ? tmpRankingComparator : tmpTableComparator.thenComparing(tmpRankingComparator));
        }, tmpTableView.comparatorProperty()));
        //the sorting is done by the sorted list, whose comparator includes the one of the table
        tmpTableView.setSortPolicy(aTable -> true);
        tmpTableView.setItems(tmpSortedItems);
        tmpTableView.setRowFactory(aTable -> new TableRow<>() {
            @Override
            protected void updateItem(ClusteringRunResult.ClusterSummary aCluster, boolean anIsEmpty) {
                super.updateItem(aCluster, anIsEmpty);
                this.updateMarkedStyle();
            }
            //
            @Override
            public void updateSelected(boolean anIsSelected) {
                super.updateSelected(anIsSelected);
                this.updateMarkedStyle();
            }
            //
            /**
             * Highlights the row of a marked cluster unless it is selected, so that the selection stays visible.
             */
            private void updateMarkedStyle() {
                this.setStyle(!this.isEmpty() && !this.isSelected() && ClusteringResultView.this.isMarked(this.getItem())
                        ? ClusteringResultView.MARKED_ROW_STYLE : "");
            }
        });
        return tmpTableView;
    }
    //
    /**
     * Creates a cell with a check box that marks its cluster; the row of the not clusterable molecules has none.
     *
     * @return mark cell
     */
    private TableCell<ClusteringRunResult.ClusterSummary, Boolean> createMarkCell() {
        return new TableCell<>() {
            /**
             * Check box of the cell.
             */
            private final CheckBox checkBox = new CheckBox();
            /**
             * Marked property the check box is currently bound to.
             */
            private BooleanProperty boundProperty;
            //
            @Override
            protected void updateItem(Boolean anItem, boolean anIsEmpty) {
                super.updateItem(anItem, anIsEmpty);
                if (this.boundProperty != null) {
                    this.checkBox.selectedProperty().unbindBidirectional(this.boundProperty);
                    this.boundProperty = null;
                }
                //the cluster of the row, read via the index, since the column has no cell value
                int tmpIndex = this.getIndex();
                ClusteringRunResult.ClusterSummary tmpCluster = anIsEmpty || tmpIndex < 0 || tmpIndex >= this.getTableView().getItems().size()
                        ? null : this.getTableView().getItems().get(tmpIndex);
                if (tmpCluster == null || tmpCluster.clusterIndex() < 0) {
                    this.setGraphic(null);
                    return;
                }
                this.boundProperty = ClusteringResultView.this.markedPropertyOfClusterMap.get(tmpCluster.clusterIndex());
                this.checkBox.selectedProperty().bindBidirectional(this.boundProperty);
                this.setGraphic(this.checkBox);
                this.setAlignment(Pos.CENTER);
            }
        };
    }
    //</editor-fold>
}
