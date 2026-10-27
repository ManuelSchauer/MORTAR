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
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.gui.views.ClusteringResultView;
import de.unijena.cheminf.mortar.gui.views.PcaView;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusterChemistryCheck;
import de.unijena.cheminf.mortar.model.clustering.ClusterDescriptorStatistics;
import de.unijena.cheminf.mortar.model.clustering.ClusterQuality;
import de.unijena.cheminf.mortar.model.clustering.ClusterStability;
import de.unijena.cheminf.mortar.model.clustering.ClusteringParameters;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;
import de.unijena.cheminf.mortar.model.clustering.DiverseSubsetSelector;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrixBuilder;
import de.unijena.cheminf.mortar.model.clustering.PrincipalComponentAnalysis;
import de.unijena.cheminf.mortar.model.clustering.SimilaritySearch;
import de.unijena.cheminf.mortar.model.clustering.TrainTestSplitter;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.FragmentationService;
import de.unijena.cheminf.mortar.model.io.Exporter;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.LogUtil;

import javafx.application.Platform;
import javafx.beans.binding.Bindings;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.value.ChangeListener;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.collections.transformation.FilteredList;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.ListView;
import javafx.scene.control.Spinner;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.cell.CheckBoxListCell;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;
import javafx.util.StringConverter;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Controls the tab of one clustering result ({@link ClusteringResultView}) and everything that is opened from it: the
 * overviews of the clusters and the cluster representatives, the similarity search, the descriptor values, the choice
 * of the displayed descriptors, the range filter, the cluster quality with the stability check, the PCA, the diverse
 * subset, the selection from a cluster, the training/test partitioning and the export. The controller is created by
 * the {@link ClusteringViewController}, which adds the tab, and is kept as user data of the tab; its background tasks
 * use the clustering controls of the {@link ClusteringViewController}.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringResultViewController {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusteringResultViewController.class.getName());
    /**
     * Training/test ratio proposed in the dialog of the cluster-based training/test partitioning.
     */
    private static final String TRAIN_TEST_RATIO_DEFAULT = "80:20";
    /**
     * Number of molecules proposed in the dialog for selecting molecules from a cluster.
     */
    private static final int CLUSTER_MOLECULE_SELECTION_NUMBER_DEFAULT = 10;
    /**
     * Width of the column with the most frequent scaffold in the quality window.
     */
    private static final double SCAFFOLD_COLUMN_WIDTH = 220.0;
    //</editor-fold>
    //
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * Controller of the clustering, which runs the background tasks and creates the clustering bars.
     */
    private final ClusteringViewController clusteringViewController;
    /**
     * Primary stage of the application, owner of the dialogs and windows.
     */
    private final Stage primaryStage;
    /**
     * Tab pane of the main view, which gets the PCA tabs.
     */
    private final TabPane mainTabPane;
    /**
     * Configuration to read resource file paths from.
     */
    private final IConfiguration configuration;
    /**
     * Global settings, e.g. the CSV separator.
     */
    private final SettingsContainer settingsContainer;
    /**
     * Clustering service, which holds the displayed descriptors and the similarity search settings.
     */
    private final ClusteringService clusteringService;
    /**
     * Fragmentation service, used to fragment the query of a similarity search.
     */
    private final FragmentationService fragmentationService;
    /**
     * Imported molecules, searched by name or SMILES for the query of a similarity search.
     */
    private final ObservableList<MoleculeDataModel> moleculeDataModelList;
    /**
     * Result of the clustering run shown in the tab.
     */
    private final ClusteringRunResult runResult;
    /**
     * View of the clustering result.
     */
    private final ClusteringResultView view;
    /**
     * Tab that shows the view; its user data is this controller.
     */
    private final Tab tab;
    /**
     * Calculation of the cluster quality in the background, cancelled when the tab is closed.
     */
    private final Task<ClusterQuality> qualityTask;
    /**
     * Updates the displayed medians when the chosen descriptors change; removed when the tab is closed.
     */
    private final ChangeListener<String> displayedDescriptorsListener;
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    /**
     * Similarity search in the feature space of the clustering run, created on first use.
     */
    private SimilaritySearch similaritySearch;
    //</editor-fold>
    //
    //<editor-fold desc="constructor" defaultstate="collapsed">
    /**
     * Constructor. Creates the view and the tab for the given clustering result, connects the buttons of the view and
     * starts the calculation of the cluster quality in the background. The clustering bar of the view and the closing of
     * the tab are set by the {@link ClusteringViewController}, which also adds the tab.
     *
     * @param aClusteringViewController controller of the clustering
     * @param aPrimaryStage primary stage of the application
     * @param aMainTabPane tab pane of the main view
     * @param aConfiguration configuration to read resource file paths from
     * @param aSettingsContainer global settings
     * @param aClusteringService clustering service
     * @param aFragmentationService fragmentation service
     * @param aMoleculeDataModelList imported molecules
     * @param aRunResult result of the clustering run
     * @param aTabTitle title of the tab
     * @throws NullPointerException if a parameter is null
     */
    ClusteringResultViewController(ClusteringViewController aClusteringViewController, Stage aPrimaryStage,
                                   TabPane aMainTabPane, IConfiguration aConfiguration, SettingsContainer aSettingsContainer,
                                   ClusteringService aClusteringService, FragmentationService aFragmentationService,
                                   ObservableList<MoleculeDataModel> aMoleculeDataModelList, ClusteringRunResult aRunResult,
                                   String aTabTitle) throws NullPointerException {
        this.clusteringViewController = Objects.requireNonNull(aClusteringViewController, "aClusteringViewController is null");
        this.primaryStage = Objects.requireNonNull(aPrimaryStage, "aPrimaryStage is null");
        this.mainTabPane = Objects.requireNonNull(aMainTabPane, "aMainTabPane is null");
        this.configuration = Objects.requireNonNull(aConfiguration, "aConfiguration is null");
        this.settingsContainer = Objects.requireNonNull(aSettingsContainer, "aSettingsContainer is null");
        this.clusteringService = Objects.requireNonNull(aClusteringService, "aClusteringService is null");
        this.fragmentationService = Objects.requireNonNull(aFragmentationService, "aFragmentationService is null");
        this.moleculeDataModelList = Objects.requireNonNull(aMoleculeDataModelList, "aMoleculeDataModelList is null");
        this.runResult = Objects.requireNonNull(aRunResult, "aRunResult is null");
        Objects.requireNonNull(aTabTitle, "aTabTitle is null");
        this.view = new ClusteringResultView(this.runResult);
        this.tab = new Tab(aTabTitle);
        this.tab.setId(TabNames.CLUSTERING.toString());
        this.tab.setClosable(true);
        this.tab.setContent(this.view);
        this.tab.setUserData(this);
        //medians of the chosen descriptors, updated in all result tabs when the choice changes
        this.view.setDisplayedDescriptorNames(this.clusteringService.getDisplayedDescriptorColumnNames());
        this.displayedDescriptorsListener = (anObservable, anOldValue, aNewValue) ->
                this.view.setDisplayedDescriptorNames(this.clusteringService.getDisplayedDescriptorColumnNames());
        this.clusteringService.displayedDescriptorColumnsProperty().addListener(this.displayedDescriptorsListener);
        this.addListener();
        //the cluster quality (all pairwise distances) is calculated in the background, the tab can be used meanwhile
        this.qualityTask = new Task<>() {
            @Override
            protected ClusterQuality call() {
                return new ClusterQuality(ClusteringResultViewController.this.runResult, this::isCancelled);
            }
        };
        this.qualityTask.setOnSucceeded(event -> this.view.setClusterQuality(this.qualityTask.getValue()));
        this.qualityTask.setOnFailed(event -> ClusteringResultViewController.LOGGER.log(Level.WARNING,
                "Cluster quality could not be calculated.", event.getSource().getException()));
        Thread tmpQualityThread = new Thread(this.qualityTask);
        tmpQualityThread.setName(MainViewController.ThreadType.CLUSTERING_THREAD.getThreadName());
        tmpQualityThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        tmpQualityThread.setDaemon(true);
        tmpQualityThread.start();
    }
    //</editor-fold>
    //
    //<editor-fold desc="package-private methods" defaultstate="collapsed">
    /**
     * Returns the view of the clustering result.
     *
     * @return clustering result view
     */
    ClusteringResultView getView() {
        return this.view;
    }
    //
    /**
     * Returns the tab that shows the view.
     *
     * @return tab of the clustering result
     */
    Tab getTab() {
        return this.tab;
    }
    //
    /**
     * Releases what the tab holds when it is closed: cancels the quality calculation and stops following the choice of
     * the displayed descriptors.
     */
    void close() {
        this.qualityTask.cancel(true);
        this.clusteringService.displayedDescriptorColumnsProperty().removeListener(this.displayedDescriptorsListener);
    }
    //
    /**
     * Opens an overview view with the representative of every cluster, ranked from the largest to the smallest cluster,
     * with the cluster size and its share of the clustered molecules below each structure, like the fragment frequencies
     * in the fragments overview. A double click on a representative opens all molecules of its cluster in another
     * overview view. The molecules that could not be clustered are not shown.
     */
    void openClusterRepresentativesOverview() {
        List<ClusteringRunResult.ClusterSummary> tmpClusters = this.runResult.getClusterSummariesRankedBySize().stream()
                .filter(aCluster -> aCluster.clusterIndex() >= 0).toList();
        if (tmpClusters.isEmpty()) {
            return;
        }
        List<MoleculeDataModel> tmpRepresentatives = new ArrayList<>(tmpClusters.size());
        List<String> tmpLabels = new ArrayList<>(tmpClusters.size());
        List<Double> tmpFractions = new ArrayList<>(tmpClusters.size());
        DecimalFormat tmpPercentageFormat = new DecimalFormat("#.#", new DecimalFormatSymbols(Locale.getDefault()));
        int tmpNumberOfMolecules = this.runResult.molecules().size();
        for (ClusteringRunResult.ClusterSummary tmpCluster : tmpClusters) {
            tmpRepresentatives.add(tmpCluster.representative());
            double tmpFraction = (double) tmpCluster.getSize() / tmpNumberOfMolecules;
            tmpLabels.add(tmpCluster.getSize() + " (" + tmpPercentageFormat.format(100.0 * tmpFraction) + "%)");
            tmpFractions.add(tmpFraction);
        }
        OverviewViewController tmpOverviewController = this.clusteringViewController.getClusterRepresentativesOverviewViewController();
        tmpOverviewController.setSimilaritySearchAction(this::searchSimilarMolecules);
        //a representative stands for its cluster, so the descriptor values of the whole cluster are shown
        Map<MoleculeDataModel, ClusteringRunResult.ClusterSummary> tmpClusterOfRepresentativeMap = new IdentityHashMap<>(tmpClusters.size());
        for (ClusteringRunResult.ClusterSummary tmpCluster : tmpClusters) {
            tmpClusterOfRepresentativeMap.put(tmpCluster.representative(), tmpCluster);
        }
        tmpOverviewController.setDescriptorValuesOfMoleculeAction(aMolecule -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpClusterOfRepresentativeMap.get(aMolecule);
            if (tmpCluster != null) {
                this.showDescriptorValuesOfCluster(tmpCluster);
            }
        }, Message.get("OverviewView.contextMenu.descriptorValuesOfClusterMenuItem"));
        tmpOverviewController.initializeAndShowOverviewView(this.primaryStage,
                OverviewViewController.DataSources.CLUSTER_REPRESENTATIVES, this.tab.getText(), tmpRepresentatives, tmpLabels,
                tmpFractions,
                //opened after the mouse event has been processed and owned by the representatives overview, so that the
                // remaining mouse events do not reach a blocked window when the cluster overview is closed
                anIndex -> Platform.runLater(() -> {
                    Stage tmpOwnerStage = tmpOverviewController.getOverviewViewStage();
                    this.openClusterOverview(tmpClusters.get(anIndex),
                            tmpOwnerStage != null && tmpOwnerStage.isShowing() ? tmpOwnerStage : this.primaryStage);
                }));
    }
    //
    /**
     * Searches the molecules most similar to the given molecule in the feature space of the clustering run and
     * shows them in an overview view, the query first and the similarity below every molecule (also as bar). A clustered
     * molecule keeps its own cluster; another molecule is assigned to the cluster with the most similar cluster vector.
     * A molecule that was not fragmented by the run's fragmentation is fragmented with the same fragmenter first, which
     * is only possible for single-algorithm fragmentations.
     *
     * @param aQuery query molecule
     */
    void searchSimilarMolecules(MoleculeDataModel aQuery) {
        String tmpFragmentationName = this.runResult.parameters().fragmentationName();
        try {
            if (this.similaritySearch == null) {
                this.similaritySearch = new SimilaritySearch(this.runResult);
            }
            int tmpNumberOfMolecules = this.clusteringService.getNumberOfSimilarMolecules();
            int tmpRowIndex = -1;
            for (int i = 0; i < this.runResult.molecules().size(); i++) {
                if (this.runResult.molecules().get(i) == aQuery) {
                    tmpRowIndex = i;
                    break;
                }
            }
            SimilaritySearch.SearchResult tmpResult;
            if (tmpRowIndex >= 0) {
                tmpResult = this.similaritySearch.searchByRowIndex(tmpRowIndex, tmpNumberOfMolecules);
            } else {
                MoleculeDataModel tmpFragmentedQuery = aQuery;
                //without fragment fingerprint, the query needs no fragmentation
                if (this.runResult.parameters().isFragmentFingerprintUsed()
                        && aQuery.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName) == null) {
                    if (!this.fragmentationService.canFragmentSingleMolecule(tmpFragmentationName)) {
                        GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringSimilaritySearch.error.title"),
                                Message.get("ClusteringSimilaritySearch.cannotFragment.header"),
                                String.format(Message.get("ClusteringSimilaritySearch.cannotFragment.content"), tmpFragmentationName));
                        return;
                    }
                    //a copy, so that the fragments of the query are not added to an imported molecule
                    tmpFragmentedQuery = new MoleculeDataModel(aQuery.getUniqueSmiles(), aQuery.getName(), new HashMap<>());
                    this.fragmentationService.fragmentSingleMolecule(tmpFragmentedQuery, tmpFragmentationName);
                }
                tmpResult = this.similaritySearch.searchByRawFeatureVector(
                        FeatureMatrixBuilder.buildRawFeatureVector(tmpFragmentedQuery, this.runResult), tmpNumberOfMolecules);
            }
            this.showSimilaritySearchResult(aQuery, tmpResult);
        } catch (Exception anException) {
            ClusteringResultViewController.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("ClusteringSimilaritySearch.error.title"),
                    Message.get("ClusteringSimilaritySearch.error.header"), anException.getMessage(), anException);
        }
    }
    //
    /**
     * Shows the descriptor values of the given molecule, see {@link #showDescriptorValues(List, String)}.
     *
     * @param aMolecule molecule
     */
    void showDescriptorValuesOfMolecule(MoleculeDataModel aMolecule) {
        this.showDescriptorValues(List.of(aMolecule), aMolecule.getName());
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Connects the buttons and menu items of the view.
     */
    private void addListener() {
        //opens the selected cluster in the overview view, the representative first
        Runnable tmpOpenOverview = () -> {
            ClusteringRunResult.ClusterSummary tmpCluster = this.view.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster == null) {
                return;
            }
            this.openClusterOverview(tmpCluster, this.primaryStage);
        };
        this.view.getClusterMoleculeSelectionMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = this.view.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null && tmpCluster.clusterIndex() >= 0) {
                this.selectMoleculesFromCluster(tmpCluster.clusterIndex());
            }
        });
        this.view.getOverviewButton().setOnAction(event -> tmpOpenOverview.run());
        this.view.getClusterTableView().setOnMouseClicked(event -> {
            if (event.getClickCount() == 2) {
                tmpOpenOverview.run();
            }
        });
        //export, training/test partitioning and diverse subset use the marked clusters, all if none is marked
        this.view.getExportButton().setOnAction(event ->
                this.clusteringViewController.exportClusteringData(this.runResult, this.view.getMarkedClusterIndices()));
        this.view.getTrainTestSplitButton().setOnAction(event -> this.createTrainTestSplit(this.view.getMarkedClusterIndices()));
        this.view.getChooseDescriptorsButton().setOnAction(event -> this.chooseDisplayedDescriptors());
        this.view.getDescriptorValuesMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = this.view.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null) {
                this.showDescriptorValuesOfCluster(tmpCluster);
            }
        });
        this.view.getRangeFilterButton().setOnAction(event -> this.applyClusterRangeFilter());
        this.view.getDiverseSubsetButton().setOnAction(event -> this.createDiverseSubset(this.view.getMarkedClusterIndices()));
        this.view.getPcaButton().setOnAction(event -> this.openPcaTab(this.view.getMarkedClusterIndices()));
        this.view.getClusterQualityButton().setOnAction(event -> this.showClusterQuality());
        //similarity search for a name or SMILES and for the representative of the selected cluster
        this.view.getSimilaritySearchButton().setOnAction(event ->
                this.searchSimilarMolecules(this.view.getSimilaritySearchTextField().getText()));
        this.view.getSimilaritySearchMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = this.view.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null && tmpCluster.representative() != null) {
                this.searchSimilarMolecules(tmpCluster.representative());
            } else if (tmpCluster != null && !tmpCluster.members().isEmpty()) {
                this.searchSimilarMolecules(tmpCluster.members().getFirst());
            }
        });
        this.view.getRepresentativesOverviewButton().setOnAction(event -> this.openClusterRepresentativesOverview());
    }
    //
    /**
     * Returns the focused window of the application or the main window, see
     * {@link ClusteringViewController#getFocusedStage(Stage)}.
     *
     * @return focused stage or the main stage
     */
    private Stage getFocusedStage() {
        return ClusteringViewController.getFocusedStage(this.primaryStage);
    }
    //
    /**
     * Creates an overview view controller for molecules of this clustering result, whose context menu offers the
     * similarity search and the descriptor values of a molecule.
     *
     * @return new overview view controller
     */
    private OverviewViewController createOverviewViewController() {
        OverviewViewController tmpController = new OverviewViewController(this.configuration, this.settingsContainer);
        tmpController.setSimilaritySearchAction(this::searchSimilarMolecules);
        tmpController.setDescriptorValuesOfMoleculeAction(this::showDescriptorValuesOfMolecule);
        return tmpController;
    }
    //
    //<editor-fold desc="clusters and similarity search" defaultstate="collapsed">
    /**
     * Opens the molecules of the given cluster in a new overview view, the representative first. Its context menu offers
     * the similarity search, the descriptor values and the selection of typical or diverse molecules from the cluster.
     *
     * @param aCluster cluster to show
     * @param anOwnerStage owner of the overview window
     */
    private void openClusterOverview(ClusteringRunResult.ClusterSummary aCluster, Stage anOwnerStage) {
        String tmpClusterName = aCluster.clusterIndex() < 0
                ? Message.get("ClusteringResultView.notClusterable.text")
                : Message.get("OverviewView.titleOfDataSource.cluster") + " " + aCluster.clusterIndex();
        //a new controller per cluster, so that its context menu entries refer to this clustering result
        OverviewViewController tmpController = this.createOverviewViewController();
        if (aCluster.clusterIndex() >= 0) {
            tmpController.setClusterMoleculeSelectionAction(aMolecule -> this.selectMoleculesFromCluster(aCluster.clusterIndex()));
        }
        tmpController.setDescriptorValuesOfAllMoleculesAction(() -> this.showDescriptorValuesOfCluster(aCluster));
        tmpController.initializeAndShowOverviewView(anOwnerStage, OverviewViewController.DataSources.CLUSTER, tmpClusterName, aCluster.members());
    }
    //
    /**
     * Searches the molecules most similar to the molecule given by the text: a clustered or imported molecule with this
     * name or unique SMILES, otherwise the text is parsed as SMILES.
     *
     * @param aText name or SMILES of the query
     */
    private void searchSimilarMolecules(String aText) {
        String tmpText = aText == null ? "" : aText.trim();
        if (tmpText.isEmpty()) {
            return;
        }
        MoleculeDataModel tmpQuery = null;
        //clustered molecules first, then all imported molecules
        for (List<MoleculeDataModel> tmpList : List.of(this.runResult.molecules(), this.moleculeDataModelList)) {
            for (MoleculeDataModel tmpMolecule : tmpList) {
                if (tmpText.equalsIgnoreCase(tmpMolecule.getName()) || tmpText.equals(tmpMolecule.getUniqueSmiles())) {
                    tmpQuery = tmpMolecule;
                    break;
                }
            }
            if (tmpQuery != null) {
                break;
            }
        }
        if (tmpQuery == null) {
            try {
                IAtomContainer tmpAtomContainer = ChemUtil.parseSmilesToAtomContainer(tmpText);
                tmpQuery = new MoleculeDataModel(tmpAtomContainer, this.settingsContainer.getRegardStereochemistrySetting());
                tmpQuery.setName(tmpText);
            } catch (Exception anException) {
                ClusteringResultViewController.LOGGER.log(Level.INFO, "Similarity search query is neither a molecule name nor a valid SMILES: " + tmpText, anException);
                GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringSimilaritySearch.error.title"),
                        Message.get("ClusteringSimilaritySearch.invalidQuery.header"),
                        String.format(Message.get("ClusteringSimilaritySearch.invalidQuery.content"), tmpText));
                return;
            }
        }
        this.searchSimilarMolecules(tmpQuery);
    }
    //
    /**
     * Shows the result of a similarity search in a new overview view: the query first (highlighted), then the most
     * similar molecules with their similarity as text and light blue bar. The title names the query and its cluster.
     *
     * @param aQuery query molecule
     * @param aResult result of the search
     */
    private void showSimilaritySearchResult(MoleculeDataModel aQuery, SimilaritySearch.SearchResult aResult) {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(aResult.neighbors().size() + 1);
        List<String> tmpLabels = new ArrayList<>(aResult.neighbors().size() + 1);
        List<Double> tmpFractions = new ArrayList<>(aResult.neighbors().size() + 1);
        //the query is highlighted as first structure, its label is not shown
        tmpMolecules.add(aQuery);
        tmpLabels.add("");
        tmpFractions.add(1.0);
        DecimalFormat tmpSimilarityFormat = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.getDefault()));
        for (SimilaritySearch.Neighbor tmpNeighbor : aResult.neighbors()) {
            tmpMolecules.add(this.runResult.molecules().get(tmpNeighbor.rowIndex()));
            tmpLabels.add(tmpSimilarityFormat.format(tmpNeighbor.similarity()));
            tmpFractions.add(Math.max(0.0, Math.min(1.0, tmpNeighbor.similarity())));
        }
        String tmpTitle;
        if (aResult.clusterIndex() < 0) {
            tmpTitle = String.format(Message.get("ClusteringSimilaritySearch.title.noCluster"), aQuery.getName());
        } else if (aResult.isWithinVigilance()) {
            tmpTitle = String.format(Message.get("ClusteringSimilaritySearch.title.cluster"), aQuery.getName(), aResult.clusterIndex());
        } else {
            tmpTitle = String.format(Message.get("ClusteringSimilaritySearch.title.newCluster"), aQuery.getName(), aResult.clusterIndex());
        }
        //a new controller per result, so that a search started from a result view opens another view on top of it
        OverviewViewController tmpController = this.createOverviewViewController();
        Stage tmpOwner = this.getFocusedStage();
        tmpController.initializeAndShowOverviewView(tmpOwner, OverviewViewController.DataSources.SIMILARITY_SEARCH, tmpTitle,
                tmpMolecules, tmpLabels, tmpFractions, null);
    }
    //
    /**
     * Asks for the number of molecules and the selection mode and shows the selected molecules of the given cluster in a
     * new overview view, with their selection rank below every structure.
     *
     * @param aClusterIndex index of the cluster
     */
    private void selectMoleculesFromCluster(int aClusterIndex) {
        //the window the action came from (e.g. the cluster overview) owns the dialog and the result window; the main
        // window is blocked by the modal overview, and focusing it on closing would cause a system beep
        Stage tmpOwner = this.getFocusedStage();
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(tmpOwner);
        tmpDialog.setTitle(Message.get("ClusterMoleculeSelection.dialog.title"));
        tmpDialog.setHeaderText(String.format(Message.get("ClusterMoleculeSelection.dialog.header"), aClusterIndex));
        Spinner<Integer> tmpNumberSpinner = new Spinner<>(1, Integer.MAX_VALUE, ClusteringResultViewController.CLUSTER_MOLECULE_SELECTION_NUMBER_DEFAULT);
        tmpNumberSpinner.setEditable(true);
        ComboBox<DiverseSubsetSelector.SelectionMode> tmpModeComboBox = ClusteringResultViewController.createSelectionModeComboBox();
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("ClusterMoleculeSelection.number.label")), tmpNumberSpinner);
        tmpGrid.addRow(1, new Label(Message.get("DiverseSubset.selectionMode.label")), tmpModeComboBox);
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        List<DiverseSubsetSelector.SelectedMolecule> tmpSelection = DiverseSubsetSelector.selectFromCluster(this.runResult,
                aClusterIndex, tmpNumberSpinner.getValue(), tmpModeComboBox.getValue());
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(tmpSelection.size());
        List<String> tmpLabels = new ArrayList<>(tmpSelection.size());
        for (DiverseSubsetSelector.SelectedMolecule tmpSelected : tmpSelection) {
            tmpMolecules.add(this.runResult.molecules().get(tmpSelected.rowIndex()));
            tmpLabels.add("#" + tmpSelected.rankInCluster());
        }
        OverviewViewController tmpController = this.createOverviewViewController();
        tmpController.initializeAndShowOverviewView(tmpOwner, OverviewViewController.DataSources.DIVERSE_SUBSET,
                String.format(Message.get("ClusterMoleculeSelection.overview.title"), aClusterIndex, tmpSelection.size(),
                        Message.get("DiverseSubset.selectionMode." + tmpModeComboBox.getValue().name())),
                tmpMolecules, tmpLabels, null, null);
    }
    //</editor-fold>
    //
    //<editor-fold desc="descriptors" defaultstate="collapsed">
    /**
     * Shows the descriptor values of all molecules of the given cluster, the representative first, see
     * {@link #showDescriptorValues(List, String)}.
     *
     * @param aCluster cluster whose molecules are shown
     */
    private void showDescriptorValuesOfCluster(ClusteringRunResult.ClusterSummary aCluster) {
        String tmpClusterName = aCluster.clusterIndex() < 0
                ? Message.get("ClusteringResultView.notClusterable.text")
                : Message.get("OverviewView.titleOfDataSource.cluster") + " " + aCluster.clusterIndex();
        this.showDescriptorValues(aCluster.members(), tmpClusterName);
    }
    //
    /**
     * Shows the untransformed values of the displayed descriptors for the given molecules in a new window: one row per
     * molecule with name, cluster, structure and one sortable column per descriptor with fixed width; the table scrolls
     * horizontally. If all molecules belong to the same cluster, the column headers show its median as reference.
     * Descriptors removed by the filter are marked with an asterisk. Molecules that are not part of the clustering are
     * left out.
     *
     * @param aMoleculeList molecules, e.g. the members of a cluster or a single molecule
     * @param aName name of the molecules shown in the window title, e.g. the cluster or the molecule name
     */
    private void showDescriptorValues(List<MoleculeDataModel> aMoleculeList, String aName) {
        ClusterDescriptorStatistics tmpStatistics = new ClusterDescriptorStatistics(this.runResult);
        List<String> tmpNames = this.clusteringService.getDisplayedDescriptorColumnNames().stream()
                .filter(tmpStatistics::hasDescriptorColumn).toList();
        if (tmpNames.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.noDescriptors.header"), null);
            return;
        }
        Map<MoleculeDataModel, Integer> tmpRowIndexMap = ClusteringResultViewController.getRowIndexMap(this.runResult);
        List<MoleculeDataModel> tmpMolecules = aMoleculeList.stream().filter(tmpRowIndexMap::containsKey).toList();
        if (tmpMolecules.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.notClustered.header"), null);
            return;
        }
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        //the cluster all molecules belong to, -1 if they belong to different clusters or could not be clustered
        int tmpCommonCluster = tmpLabels[tmpRowIndexMap.get(tmpMolecules.getFirst())];
        for (MoleculeDataModel tmpMolecule : tmpMolecules) {
            if (tmpLabels[tmpRowIndexMap.get(tmpMolecule)] != tmpCommonCluster) {
                tmpCommonCluster = -1;
                break;
            }
        }
        TableView<MoleculeDataModel> tmpTableView = new TableView<>(FXCollections.observableArrayList(tmpMolecules));
        //fixed column widths, with many descriptors the table scrolls horizontally
        tmpTableView.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
        TableColumn<MoleculeDataModel, String> tmpNameColumn = new TableColumn<>(Message.get("ClusteringResultView.descriptorValues.nameColumn.header"));
        tmpNameColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(aCellData.getValue().getName()));
        tmpNameColumn.setPrefWidth(150.0);
        TableColumn<MoleculeDataModel, Integer> tmpClusterColumn = new TableColumn<>(Message.get("ClusteringResultView.clusterColumn.header"));
        tmpClusterColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(tmpLabels[tmpRowIndexMap.get(aCellData.getValue())]));
        tmpClusterColumn.setPrefWidth(70.0);
        TableColumn<MoleculeDataModel, MoleculeDataModel> tmpStructureColumn = new TableColumn<>(Message.get("ClusteringResultView.descriptorValues.structureColumn.header"));
        tmpStructureColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue()));
        //the structure of the molecule of the row is depicted by the cell itself when it becomes visible, centred and
        // fitted to the column
        tmpStructureColumn.setCellFactory(aColumn -> ClusteringResultView.createStructureCell());
        tmpStructureColumn.setSortable(false);
        tmpStructureColumn.setPrefWidth(ClusteringResultView.STRUCTURE_COLUMN_WIDTH);
        tmpTableView.getColumns().add(tmpNameColumn);
        tmpTableView.getColumns().add(tmpClusterColumn);
        tmpTableView.getColumns().add(tmpStructureColumn);
        for (String tmpName : tmpNames) {
            String tmpHeader = tmpStatistics.isUsedForClustering(tmpName) ? tmpName : tmpName + " *";
            //the median of the common cluster is shown in the header, as reference for the values of the molecules
            String tmpMedianLine = tmpCommonCluster < 0 ? "" : String.format(Message.get("ClusteringResultView.descriptorValues.median.text"),
                    String.format(Locale.getDefault(), "%.2f", tmpStatistics.getMedian(tmpCommonCluster, tmpName)));
            TableColumn<MoleculeDataModel, Double> tmpColumn = new TableColumn<>(tmpMedianLine.isEmpty() ? tmpHeader
                    : tmpHeader + System.lineSeparator() + tmpMedianLine);
            ClusteringResultView.setFixedWidthForHeader(tmpColumn, tmpHeader.length() >= tmpMedianLine.length() ? tmpHeader : tmpMedianLine);
            //null instead of NaN, so that molecules without value are sorted consistently
            tmpColumn.setCellValueFactory(aCellData -> {
                double tmpValue = tmpStatistics.getValue(tmpRowIndexMap.get(aCellData.getValue()), tmpName);
                return new ReadOnlyObjectWrapper<>(Double.isNaN(tmpValue) ? null : tmpValue);
            });
            tmpColumn.setCellFactory(aColumn -> new TableCell<>() {
                @Override
                protected void updateItem(Double aValue, boolean anIsEmpty) {
                    super.updateItem(aValue, anIsEmpty);
                    this.setText(anIsEmpty ? null : aValue == null ? "-" : String.format(Locale.getDefault(), "%.2f", aValue));
                }
            });
            tmpTableView.getColumns().add(tmpColumn);
        }
        Label tmpHintLabel = new Label(Message.get("ClusteringResultView.descriptorValues.hint.text"));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        BorderPane tmpPane = new BorderPane(tmpTableView);
        tmpPane.setTop(tmpHintLabel);
        Stage tmpStage = new Stage();
        tmpStage.initOwner(this.getFocusedStage());
        tmpStage.setTitle(String.format(Message.get("ClusteringResultView.descriptorValues.window.title"), aName, tmpMolecules.size()));
        //a single molecule needs only a small window
        tmpStage.setScene(new Scene(tmpPane, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE,
                tmpMolecules.size() == 1 ? GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE / 2.0 : GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE));
        tmpStage.show();
    }
    //
    /**
     * Opens a dialog to choose the descriptors whose medians are shown in the clustering result tabs: a searchable list of
     * all descriptor columns of this result with check boxes. The choice is stored in the clustering settings and
     * applied to all result tabs; chosen descriptors that were not calculated in this result are kept.
     */
    private void chooseDisplayedDescriptors() {
        List<String> tmpAvailableNames = this.view.getDescriptorStatistics().getDescriptorColumnNames();
        List<String> tmpChosenNames = this.clusteringService.getDisplayedDescriptorColumnNames();
        Map<String, BooleanProperty> tmpCheckedPropertyMap = new LinkedHashMap<>(
                CollectionUtil.calculateInitialHashCollectionCapacity(tmpAvailableNames.size()));
        for (String tmpName : tmpAvailableNames) {
            tmpCheckedPropertyMap.put(tmpName, new SimpleBooleanProperty(tmpChosenNames.contains(tmpName)));
        }
        TextField tmpSearchField = new TextField();
        tmpSearchField.setPromptText(Message.get("ClusteringResultView.chooseDescriptors.search.prompt"));
        FilteredList<String> tmpFilteredNames = new FilteredList<>(FXCollections.observableArrayList(tmpAvailableNames));
        tmpSearchField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpFilteredNames.setPredicate(aName ->
                aNewValue == null || aNewValue.isBlank() || aName.toLowerCase(Locale.ROOT).contains(aNewValue.trim().toLowerCase(Locale.ROOT))));
        ListView<String> tmpListView = new ListView<>(tmpFilteredNames);
        tmpListView.setCellFactory(CheckBoxListCell.forListView(tmpCheckedPropertyMap::get));
        tmpListView.setPrefHeight(400);
        Button tmpNoneButton = new Button(Message.get("ClusteringResultView.markNoneButton.text"));
        tmpNoneButton.setOnAction(anEvent -> tmpCheckedPropertyMap.values().forEach(aProperty -> aProperty.set(false)));
        HBox tmpSearchBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpSearchField, tmpNoneButton);
        HBox.setHgrow(tmpSearchField, Priority.ALWAYS);
        VBox tmpContentBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpSearchBox, tmpListView);
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.getFocusedStage());
        tmpDialog.setTitle(Message.get("ClusteringResultView.chooseDescriptors.title"));
        tmpDialog.setHeaderText(Message.get("ClusteringResultView.chooseDescriptors.header"));
        tmpDialog.getDialogPane().setContent(tmpContentBox);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        //the previous order is kept, newly chosen descriptors are appended in column order
        List<String> tmpNewNames = new ArrayList<>(tmpChosenNames.size());
        for (String tmpName : tmpChosenNames) {
            BooleanProperty tmpCheckedProperty = tmpCheckedPropertyMap.get(tmpName);
            if (tmpCheckedProperty == null || tmpCheckedProperty.get()) {
                tmpNewNames.add(tmpName);
            }
        }
        tmpCheckedPropertyMap.forEach((aName, aProperty) -> {
            if (aProperty.get() && !tmpNewNames.contains(aName)) {
                tmpNewNames.add(aName);
            }
        });
        this.clusteringService.setDisplayedDescriptorColumnNames(tmpNewNames);
    }
    //
    /**
     * Opens the range filter for the clusters of this result: one or more conditions, each a descriptor with an
     * optional lower and upper bound on its median in the cluster. The number of matching clusters is updated while the
     * conditions are edited. On OK, exactly the clusters that satisfy all conditions are marked and the descriptors of
     * the conditions are added to the displayed descriptors.
     */
    private void applyClusterRangeFilter() {
        ClusterDescriptorStatistics tmpStatistics = this.view.getDescriptorStatistics();
        List<String> tmpNames = tmpStatistics.getDescriptorColumnNames();
        if (tmpNames.isEmpty()) {
            return;
        }
        GridPane tmpConditionGrid = new GridPane();
        tmpConditionGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpConditionGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpConditionGrid.addRow(0, new Label(Message.get("ClusteringResultView.rangeFilter.descriptor.header")),
                new Label(Message.get("ClusteringResultView.rangeFilter.minimum.header")),
                new Label(Message.get("ClusteringResultView.rangeFilter.maximum.header")),
                new Label(Message.get("ClusteringResultView.rangeFilter.overallMedian.header")));
        List<ComboBox<String>> tmpDescriptorComboBoxes = new ArrayList<>();
        List<TextField> tmpMinimumFields = new ArrayList<>();
        List<TextField> tmpMaximumFields = new ArrayList<>();
        Label tmpMatchLabel = new Label();
        //counts the matching clusters for the current input
        Runnable tmpUpdateMatchLabel = () -> {
            List<ClusterDescriptorStatistics.RangeCondition> tmpConditions = ClusteringResultViewController.readRangeConditions(
                    tmpDescriptorComboBoxes, tmpMinimumFields, tmpMaximumFields);
            tmpMatchLabel.setText(tmpConditions == null
                    ? Message.get("ClusteringResultView.rangeFilter.invalidInput.text")
                    : String.format(Message.get("ClusteringResultView.rangeFilter.matching.text"),
                            tmpStatistics.getMatchingClusters(tmpConditions).size(), this.runResult.clusteringResult().numberOfClusters()));
        };
        //adds a condition row, filled with the given condition if not null
        Consumer<ClusterDescriptorStatistics.RangeCondition> tmpAddRow = aCondition -> {
            ComboBox<String> tmpComboBox = new ComboBox<>(FXCollections.observableArrayList(tmpNames));
            tmpComboBox.setValue(aCondition != null && tmpNames.contains(aCondition.columnName()) ? aCondition.columnName()
                    : (this.view.getDisplayedDescriptorNames().isEmpty() ? tmpNames.getFirst() : this.view.getDisplayedDescriptorNames().getFirst()));
            TextField tmpMinimumField = new TextField(aCondition == null || Double.isNaN(aCondition.lowerBound()) ? "" : Double.toString(aCondition.lowerBound()));
            TextField tmpMaximumField = new TextField(aCondition == null || Double.isNaN(aCondition.upperBound()) ? "" : Double.toString(aCondition.upperBound()));
            tmpMinimumField.setPromptText(Message.get("ClusteringResultView.rangeFilter.noBound.prompt"));
            tmpMaximumField.setPromptText(Message.get("ClusteringResultView.rangeFilter.noBound.prompt"));
            tmpMinimumField.setPrefColumnCount(8);
            tmpMaximumField.setPrefColumnCount(8);
            Label tmpOverallMedianLabel = new Label();
            //the overall median helps to choose the bounds
            Runnable tmpUpdateOverallMedian = () -> tmpOverallMedianLabel.setText(tmpComboBox.getValue() == null ? ""
                    : String.format(Locale.getDefault(), "%.2f", tmpStatistics.getOverallMedian(tmpComboBox.getValue())));
            tmpUpdateOverallMedian.run();
            tmpComboBox.valueProperty().addListener((anObservable, anOldValue, aNewValue) -> {
                tmpUpdateOverallMedian.run();
                tmpUpdateMatchLabel.run();
            });
            tmpMinimumField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpUpdateMatchLabel.run());
            tmpMaximumField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpUpdateMatchLabel.run());
            tmpDescriptorComboBoxes.add(tmpComboBox);
            tmpMinimumFields.add(tmpMinimumField);
            tmpMaximumFields.add(tmpMaximumField);
            tmpConditionGrid.addRow(tmpDescriptorComboBoxes.size(), tmpComboBox, tmpMinimumField, tmpMaximumField, tmpOverallMedianLabel);
            tmpUpdateMatchLabel.run();
        };
        if (this.view.getLastRangeConditions().isEmpty()) {
            tmpAddRow.accept(null);
        } else {
            this.view.getLastRangeConditions().forEach(tmpAddRow);
        }
        Button tmpAddButton = new Button(Message.get("ClusteringResultView.rangeFilter.addCondition.text"));
        tmpAddButton.setOnAction(anEvent -> {
            tmpAddRow.accept(null);
            tmpConditionGrid.getScene().getWindow().sizeToScene();
        });
        Label tmpHintLabel = new Label(Message.get("ClusteringResultView.rangeFilter.hint.text"));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setMaxWidth(500);
        VBox tmpContentBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpConditionGrid, tmpAddButton, tmpMatchLabel, tmpHintLabel);
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.getFocusedStage());
        tmpDialog.setTitle(Message.get("ClusteringResultView.rangeFilter.title"));
        tmpDialog.setHeaderText(Message.get("ClusteringResultView.rangeFilter.header"));
        tmpDialog.getDialogPane().setContent(tmpContentBox);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        //OK only with valid bounds
        tmpDialog.getDialogPane().lookupButton(ButtonType.OK).disableProperty().bind(Bindings.createBooleanBinding(
                () -> tmpMatchLabel.getText().equals(Message.get("ClusteringResultView.rangeFilter.invalidInput.text")),
                tmpMatchLabel.textProperty()));
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        List<ClusterDescriptorStatistics.RangeCondition> tmpConditions = ClusteringResultViewController.readRangeConditions(
                tmpDescriptorComboBoxes, tmpMinimumFields, tmpMaximumFields);
        if (tmpConditions == null) {
            return;
        }
        this.view.setLastRangeConditions(tmpConditions);
        this.view.setMarkedClusterIndices(tmpStatistics.getMatchingClusters(tmpConditions));
        //the medians the filter refers to are shown in the table
        List<String> tmpDisplayedNames = new ArrayList<>(this.clusteringService.getDisplayedDescriptorColumnNames());
        for (ClusterDescriptorStatistics.RangeCondition tmpCondition : tmpConditions) {
            if (!tmpDisplayedNames.contains(tmpCondition.columnName())) {
                tmpDisplayedNames.add(tmpCondition.columnName());
            }
        }
        this.clusteringService.setDisplayedDescriptorColumnNames(tmpDisplayedNames);
    }
    //</editor-fold>
    //
    //<editor-fold desc="quality, stability and PCA" defaultstate="collapsed">
    /**
     * Shows the quality of all clusters of this result in a new window (see {@link ClusterQuality}): overall
     * silhouette and neighbour purity at the top, a sortable table with one row per cluster and, for the selected
     * cluster, its distances to all other clusters (mean distance of the molecules, distance of the centres, share of
     * neighbour contacts), nearest first. The optional stability check clusters the features again with other seeds and
     * adds the stability of every cluster; the quality can be exported as CSV. A double click on a row opens the cluster.
     */
    private void showClusterQuality() {
        ClusterQuality tmpQuality = this.view.getClusterQuality();
        if (tmpQuality == null) {
            return;
        }
        Map<Integer, ClusteringRunResult.ClusterSummary> tmpClusterMap = this.getClusterMap();
        //table of all clusters
        TableView<ClusterQuality.ClusterQualityEntry> tmpQualityTable = new TableView<>(FXCollections.observableArrayList(tmpQuality.getEntries()));
        tmpQualityTable.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpClusterColumn = new TableColumn<>(Message.get("ClusteringResultView.clusterColumn.header"));
        tmpClusterColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().clusterIndex()));
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpSizeColumn = new TableColumn<>(Message.get("ClusteringResultView.sizeColumn.header"));
        tmpSizeColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().size()));
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpNegativeColumn = new TableColumn<>(Message.get("ClusterQuality.negativeColumn.header"));
        tmpNegativeColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().numberOfNegativeSilhouettes()));
        //number of neighbours used for the purity, min(10, size - 1), so that small clusters can reach 100 %
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpNeighbourCountColumn = new TableColumn<>(Message.get("ClusterQuality.neighbourCountColumn.header"));
        tmpNeighbourCountColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(
                aCellData.getValue().numberOfNeighbours() == 0 ? null : aCellData.getValue().numberOfNeighbours()));
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpNearestColumn = new TableColumn<>(Message.get("ClusterQuality.nearestColumn.header"));
        tmpNearestColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(
                aCellData.getValue().nearestClusterIndex() < 0 ? null : aCellData.getValue().nearestClusterIndex()));
        TableColumn<ClusterQuality.ClusterQualityEntry, Integer> tmpContactColumn = new TableColumn<>(Message.get("ClusterQuality.contactColumn.header"));
        tmpContactColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(
                aCellData.getValue().mostContactsClusterIndex() < 0 ? null : aCellData.getValue().mostContactsClusterIndex()));
        TableColumn<ClusterQuality.ClusterQualityEntry, Double> tmpStabilityColumn = ClusteringResultView.createNumberColumn(
                Message.get("ClusterQuality.stabilityColumn.header"), anEntry -> this.view.getStabilityResult() == null ? Double.NaN
                        : this.view.getStabilityResult().getClusterJaccardIndex(anEntry.clusterIndex()), "%.2f");
        //values of the optional chemical check, "-" until it was run
        TableColumn<ClusterQuality.ClusterQualityEntry, String> tmpScaffoldColumn = new TableColumn<>(Message.get("ClusterQuality.scaffoldColumn.header"));
        tmpScaffoldColumn.setCellValueFactory(aCellData -> {
            ClusterChemistryCheck.ChemistryCheckResult tmpResult = this.view.getChemistryCheckResult();
            String tmpScaffold = tmpResult == null ? null : tmpResult.getMostFrequentScaffold(aCellData.getValue().clusterIndex());
            return new ReadOnlyStringWrapper(tmpScaffold == null ? "-"
                    : tmpScaffold.isEmpty() ? Message.get("ClusterQuality.acyclicScaffold.text") : tmpScaffold);
        });
        List<TableColumn<ClusterQuality.ClusterQualityEntry, ?>> tmpColumns = List.of(tmpClusterColumn, tmpSizeColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.silhouetteColumn.header"),
                        ClusterQuality.ClusterQualityEntry::silhouette, "%+.3f"),
                tmpNegativeColumn,
                tmpNeighbourCountColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.purityColumn.header"),
                        anEntry -> 100.0 * anEntry.neighbourPurity(), "%.0f %%"),
                //a large cluster has a high purity by chance alone
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.randomPurityColumn.header"),
                        anEntry -> 100.0 * anEntry.randomNeighbourPurity(), "%.0f %%"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.withinColumn.header"),
                        ClusterQuality.ClusterQualityEntry::meanDistanceWithinCluster, "%.3f"),
                tmpNearestColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.nearestDistanceColumn.header"),
                        ClusterQuality.ClusterQualityEntry::meanDistanceToNearestCluster, "%.3f"),
                tmpContactColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.contactShareColumn.header"),
                        anEntry -> 100.0 * anEntry.contactShareWithMostContacts(), "%.0f %%"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.cohesionColumn.header"),
                        ClusterQuality.ClusterQualityEntry::meanSimilarityToOwnCentre, "%.3f"),
                tmpStabilityColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.tanimotoColumn.header"),
                        anEntry -> this.view.getChemistryCheckResult() == null ? Double.NaN
                                : this.view.getChemistryCheckResult().getMeanTanimoto(anEntry.clusterIndex()), "%.2f"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.randomTanimotoColumn.header"),
                        anEntry -> this.view.getChemistryCheckResult() == null ? Double.NaN
                                : this.view.getChemistryCheckResult().getRandomMeanTanimoto(anEntry.clusterIndex()), "%.2f"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.scaffoldShareColumn.header"),
                        anEntry -> this.view.getChemistryCheckResult() == null ? Double.NaN
                                : 100.0 * this.view.getChemistryCheckResult().getScaffoldShare(anEntry.clusterIndex()), "%.0f %%"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.randomScaffoldShareColumn.header"),
                        anEntry -> this.view.getChemistryCheckResult() == null ? Double.NaN
                                : 100.0 * this.view.getChemistryCheckResult().getRandomScaffoldShare(anEntry.clusterIndex()), "%.0f %%"),
                tmpScaffoldColumn);
        for (TableColumn<ClusterQuality.ClusterQualityEntry, ?> tmpColumn : tmpColumns) {
            ClusteringResultView.setFixedWidthForHeader(tmpColumn, tmpColumn.getText());
            tmpQualityTable.getColumns().add(tmpColumn);
        }
        //the scaffold SMILES need more room than the header
        tmpScaffoldColumn.setPrefWidth(ClusteringResultViewController.SCAFFOLD_COLUMN_WIDTH);
        //distances of the selected cluster to all other clusters
        TableView<ClusterQuality.ClusterDistance> tmpDistanceTable = new TableView<>();
        tmpDistanceTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        tmpDistanceTable.setPlaceholder(new Label(Message.get("ClusterQuality.distances.noSelection.text")));
        TableColumn<ClusterQuality.ClusterDistance, Integer> tmpOtherColumn = new TableColumn<>(Message.get("ClusteringResultView.clusterColumn.header"));
        tmpOtherColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().clusterIndex()));
        TableColumn<ClusterQuality.ClusterDistance, Integer> tmpOtherSizeColumn = new TableColumn<>(Message.get("ClusteringResultView.sizeColumn.header"));
        tmpOtherSizeColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(tmpQuality.getEntry(aCellData.getValue().clusterIndex()).size()));
        tmpDistanceTable.getColumns().addAll(List.of(tmpOtherColumn, tmpOtherSizeColumn,
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.meanDistanceColumn.header"),
                        ClusterQuality.ClusterDistance::meanDistance, "%.3f"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.centreDistanceColumn.header"),
                        ClusterQuality.ClusterDistance::centreDistance, "%.3f"),
                ClusteringResultView.createNumberColumn(Message.get("ClusterQuality.neighbourContactColumn.header"),
                        aDistance -> 100.0 * aDistance.neighbourContactShare(), "%.1f %%")));
        Label tmpDistanceLabel = new Label(Message.get("ClusterQuality.distances.noSelection.text"));
        tmpDistanceLabel.setWrapText(true);
        tmpQualityTable.getSelectionModel().selectedItemProperty().addListener((anObservable, anOldValue, aNewValue) -> {
            if (aNewValue == null) {
                tmpDistanceTable.getItems().clear();
                tmpDistanceLabel.setText(Message.get("ClusterQuality.distances.noSelection.text"));
                return;
            }
            tmpDistanceTable.setItems(FXCollections.observableArrayList(tmpQuality.getDistancesToOtherClusters(aNewValue.clusterIndex())));
            tmpDistanceLabel.setText(String.format(Message.get("ClusterQuality.distances.title.text"), aNewValue.clusterIndex()));
        });
        //a double click opens the cluster in the overview
        tmpQualityTable.setOnMouseClicked(event -> {
            ClusterQuality.ClusterQualityEntry tmpEntry = tmpQualityTable.getSelectionModel().getSelectedItem();
            if (event.getClickCount() == 2 && tmpEntry != null && tmpClusterMap.containsKey(tmpEntry.clusterIndex())) {
                this.openClusterOverview(tmpClusterMap.get(tmpEntry.clusterIndex()), this.getFocusedStage());
            }
        });
        tmpDistanceTable.setOnMouseClicked(event -> {
            ClusterQuality.ClusterDistance tmpDistance = tmpDistanceTable.getSelectionModel().getSelectedItem();
            if (event.getClickCount() == 2 && tmpDistance != null && tmpClusterMap.containsKey(tmpDistance.clusterIndex())) {
                this.openClusterOverview(tmpClusterMap.get(tmpDistance.clusterIndex()), this.getFocusedStage());
            }
        });
        //overall values, the stability once it was checked
        Label tmpOverallLabel = new Label(String.format(Locale.getDefault(), Message.get("ClusterQuality.overall.text"),
                tmpQuality.getOverallSilhouette(), 100.0 * tmpQuality.getOverallNeighbourPurity(), ClusterQuality.NUMBER_OF_NEIGHBOURS,
                ClusterQuality.NUMBER_OF_NEIGHBOURS + 1)
                + (tmpQuality.isSampled() ? " " + String.format(Message.get("ClusterQuality.sampled.text"), tmpQuality.getNumberOfEvaluatedMolecules()) : ""));
        tmpOverallLabel.setStyle("-fx-font-weight: bold;");
        tmpOverallLabel.setWrapText(true);
        //random baseline as reference for the overall values
        Label tmpRandomLabel = new Label(String.format(Locale.getDefault(), Message.get("ClusterQuality.random.text"),
                ClusterQuality.NUMBER_OF_RANDOM_ASSIGNMENTS, tmpQuality.getRandomOverallSilhouette(),
                tmpQuality.getRandomOverallSilhouetteDeviation(), 100.0 * tmpQuality.getRandomOverallNeighbourPurity(),
                100.0 * tmpQuality.getRandomOverallNeighbourPurityDeviation()));
        tmpRandomLabel.setWrapText(true);
        Label tmpStabilityLabel = new Label(this.view.getStabilityResult() == null ? Message.get("ClusterQuality.stability.notChecked.text")
                : ClusteringResultViewController.formatStabilityText(this.view.getStabilityResult()));
        tmpStabilityLabel.setWrapText(true);
        Button tmpStabilityButton = new Button(Message.get("ClusterQuality.stabilityButton.text"));
        tmpStabilityButton.setTooltip(GuiUtil.createTooltip(Message.get("ClusterQuality.stabilityButton.tooltip")));
        tmpStabilityButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpStabilityButton.disableProperty().bind(this.clusteringViewController.isClusteringRunningProperty());
        tmpStabilityButton.setOnAction(event -> this.checkClusterStability(tmpQualityTable, tmpStabilityLabel));
        //optional chemical check with a fingerprint independent of the features and the scaffolds
        Label tmpChemistryLabel = new Label(this.view.getChemistryCheckResult() == null ? Message.get("ClusterQuality.chemistry.notChecked.text")
                : this.formatChemistryCheckText(this.view.getChemistryCheckResult()));
        tmpChemistryLabel.setWrapText(true);
        Button tmpChemistryButton = new Button(Message.get("ClusterQuality.chemistryButton.text"));
        tmpChemistryButton.setTooltip(GuiUtil.createTooltip(Message.get("ClusterQuality.chemistryButton.tooltip")));
        tmpChemistryButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpChemistryButton.disableProperty().bind(this.clusteringViewController.isClusteringRunningProperty());
        tmpChemistryButton.setOnAction(event -> this.checkClusterChemistry(tmpQualityTable, tmpChemistryLabel));
        Button tmpExportButton = new Button(Message.get("ClusterQuality.exportButton.text"));
        tmpExportButton.setMinWidth(Region.USE_PREF_SIZE);
        Label tmpHintLabel = new Label(Message.get("ClusterQuality.hint.text"));
        tmpHintLabel.setWrapText(true);
        HBox tmpStabilityBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpStabilityButton, tmpStabilityLabel);
        tmpStabilityBox.setAlignment(Pos.CENTER_LEFT);
        HBox tmpChemistryBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpChemistryButton, tmpChemistryLabel);
        tmpChemistryBox.setAlignment(Pos.CENTER_LEFT);
        HBox.setHgrow(tmpChemistryLabel, Priority.ALWAYS);
        VBox tmpTopBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpOverallLabel, tmpRandomLabel, tmpStabilityBox, tmpChemistryBox, tmpHintLabel);
        tmpTopBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        VBox tmpDistanceBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpDistanceLabel, tmpDistanceTable);
        tmpDistanceBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        VBox.setVgrow(tmpDistanceTable, Priority.ALWAYS);
        SplitPane tmpSplitPane = new SplitPane(tmpQualityTable, tmpDistanceBox);
        tmpSplitPane.setDividerPositions(0.68);
        HBox tmpBottomBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpExportButton);
        tmpBottomBox.setAlignment(Pos.CENTER_RIGHT);
        tmpBottomBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        BorderPane tmpPane = new BorderPane(tmpSplitPane);
        tmpPane.setTop(tmpTopBox);
        tmpPane.setBottom(tmpBottomBox);
        //the cluster selected in the result tab is preselected
        ClusteringRunResult.ClusterSummary tmpSelectedCluster = this.view.getClusterTableView().getSelectionModel().getSelectedItem();
        if (tmpSelectedCluster != null && tmpSelectedCluster.clusterIndex() >= 0) {
            tmpQualityTable.getSelectionModel().select(tmpQuality.getEntry(tmpSelectedCluster.clusterIndex()));
            tmpQualityTable.scrollTo(tmpQuality.getEntry(tmpSelectedCluster.clusterIndex()));
        }
        Stage tmpStage = new Stage();
        tmpStage.initOwner(this.primaryStage);
        tmpStage.setTitle(String.format(Message.get("ClusterQuality.window.title"), this.tab.getText()));
        tmpStage.setScene(new Scene(tmpPane, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE, GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE));
        tmpExportButton.setOnAction(event -> {
            Exporter tmpExporter = new Exporter(this.settingsContainer);
            File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(tmpStage, Exporter.ExportTypes.CLUSTER_QUALITY_CSV_FILE,
                    this.runResult.parameters().fragmentationName());
            if (tmpDirectory != null) {
                ClusteringViewController.showExportOutcome(tmpExporter.exportClusterQualityCsvFile(tmpDirectory, this.runResult, tmpQuality,
                        this.view.getStabilityResult(), this.view.getChemistryCheckResult(),
                        this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
            }
        });
        tmpStage.show();
    }
    //
    /**
     * Runs the optional stability check of this clustering result in the background (see {@link ClusterStability})
     * and shows its result in the given quality table and label; the result is kept in the result view.
     *
     * @param aQualityTable quality table whose stability column is updated
     * @param aStabilityLabel label that shows the mean adjusted Rand index
     */
    private void checkClusterStability(TableView<ClusterQuality.ClusterQualityEntry> aQualityTable, Label aStabilityLabel) {
        if (this.clusteringViewController.isClusteringRunning()) {
            return;
        }
        //mode and number of runs
        ComboBox<ClusterStability.StabilityMode> tmpModeComboBox = new ComboBox<>(FXCollections.observableArrayList(ClusterStability.StabilityMode.values()));
        tmpModeComboBox.setConverter(new StringConverter<>() {
            @Override
            public String toString(ClusterStability.StabilityMode aMode) {
                return aMode == null ? "" : String.format(Message.get("ClusterStability.mode." + aMode.name()), 100.0 * ClusterStability.SUBSAMPLE_FRACTION);
            }
            @Override
            public ClusterStability.StabilityMode fromString(String aString) {
                return null;
            }
        });
        tmpModeComboBox.setValue(this.view.getStabilityResult() == null ? ClusterStability.StabilityMode.SEEDS : this.view.getStabilityResult().mode());
        Spinner<Integer> tmpRunsSpinner = new Spinner<>(1, ClusterStability.MAXIMUM_NUMBER_OF_RUNS, ClusterStability.NUMBER_OF_RUNS_DEFAULT);
        tmpRunsSpinner.setEditable(true);
        Label tmpHintLabel = new Label(String.format(Message.get("ClusterStability.dialog.hint.text"), 100.0 * ClusterStability.SUBSAMPLE_FRACTION));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setMaxWidth(480);
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("ClusterStability.mode.label")), tmpModeComboBox);
        tmpGrid.addRow(1, new Label(Message.get("ClusterStability.runs.label")), tmpRunsSpinner);
        tmpGrid.add(tmpHintLabel, 0, 2, 2, 1);
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.getFocusedStage());
        tmpDialog.setTitle(Message.get("ClusterStability.dialog.title"));
        tmpDialog.setHeaderText(Message.get("ClusterQuality.stabilityButton.tooltip"));
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        ClusterStability.StabilityMode tmpMode = tmpModeComboBox.getValue();
        int tmpNumberOfRuns;
        try {
            //a typed value is only taken over by the spinner with enter, so it is read from the editor
            tmpNumberOfRuns = Integer.parseInt(tmpRunsSpinner.getEditor().getText().trim());
        } catch (NumberFormatException anException) {
            tmpNumberOfRuns = tmpRunsSpinner.getValue();
        }
        int tmpFinalNumberOfRuns = Math.max(1, Math.min(ClusterStability.MAXIMUM_NUMBER_OF_RUNS, tmpNumberOfRuns));
        Task<ClusterStability.StabilityResult> tmpTask = new Task<>() {
            @Override
            protected ClusterStability.StabilityResult call() throws Exception {
                return ClusterStability.check(ClusteringResultViewController.this.runResult, tmpMode, tmpFinalNumberOfRuns,
                        this::isCancelled, (aDone, aTotal) -> this.updateProgress(aDone, aTotal));
            }
        };
        String tmpPreviousText = aStabilityLabel.getText();
        aStabilityLabel.setText(Message.get("ClusterQuality.stability.running.text"));
        this.clusteringViewController.startBackgroundTask(tmpTask, Message.get("Status.stability"),
                Message.get("ClusterQuality.stabilityButton.tooltip"), Message.get("MainViewController.ClusteringError.Header"),
                aResult -> {
                    this.view.setStabilityResult(aResult);
                    aStabilityLabel.setText(ClusteringResultViewController.formatStabilityText(aResult));
                    aQualityTable.refresh();
                },
                //a cancelled or failed check keeps the previous result
                () -> aStabilityLabel.setText(tmpPreviousText));
    }
    //
    /**
     * Returns the text of the overall value of the given stability check, naming its mode.
     *
     * @param aResult result of the stability check
     * @return text for the quality window
     */
    private static String formatStabilityText(ClusterStability.StabilityResult aResult) {
        return String.format(Locale.getDefault(), Message.get("ClusterQuality.stability." + aResult.mode().name() + ".text"),
                aResult.meanAdjustedRandIndex(), aResult.numberOfRuns(), 100.0 * ClusterStability.SUBSAMPLE_FRACTION);
    }
    //
    /**
     * Asks for the fingerprint of the chemical check and runs the check in the background (see
     * {@link ClusterChemistryCheck}); its result is shown in the given quality table and label and kept in the result
     * view.
     *
     * @param aQualityTable quality table whose chemistry columns are updated
     * @param aChemistryLabel label that shows the overall values
     */
    private void checkClusterChemistry(TableView<ClusterQuality.ClusterQualityEntry> aQualityTable, Label aChemistryLabel) {
        if (this.clusteringViewController.isClusteringRunning()) {
            return;
        }
        Descriptor tmpFingerprint = this.chooseChemistryCheckFingerprint();
        if (tmpFingerprint == null) {
            return;
        }
        Task<ClusterChemistryCheck.ChemistryCheckResult> tmpTask = new Task<>() {
            @Override
            protected ClusterChemistryCheck.ChemistryCheckResult call() throws Exception {
                return ClusterChemistryCheck.check(ClusteringResultViewController.this.runResult, tmpFingerprint,
                        ClusterChemistryCheck.NUMBER_OF_RANDOM_RUNS_DEFAULT, this::isCancelled,
                        (aDone, aTotal) -> this.updateProgress(aDone, aTotal));
            }
        };
        String tmpPreviousText = aChemistryLabel.getText();
        aChemistryLabel.setText(Message.get("ClusterQuality.chemistry.running.text"));
        this.clusteringViewController.startBackgroundTask(tmpTask, Message.get("Status.chemistryCheck"),
                Message.get("ClusterChemistryCheck.dialog.title"), Message.get("ClusterQuality.chemistry.error.header"),
                aResult -> {
                    this.view.setChemistryCheckResult(aResult);
                    aChemistryLabel.setText(this.formatChemistryCheckText(aResult));
                    aQualityTable.refresh();
                },
                //a cancelled or failed check keeps the previous result
                () -> aChemistryLabel.setText(tmpPreviousText));
    }
    //
    /**
     * Opens the dialog for the fingerprint of the chemical check: a choice of the bit fingerprints, preset to one that
     * was not used for the clustering, a hint why an independent fingerprint is better and a warning if the chosen one
     * was used for the clustering or is related to a fingerprint used for it.
     *
     * @return chosen fingerprint or null if the dialog was cancelled
     */
    private Descriptor chooseChemistryCheckFingerprint() {
        ClusteringParameters tmpParameters = this.runResult.parameters();
        ComboBox<Descriptor> tmpComboBox = new ComboBox<>(FXCollections.observableArrayList(ClusterChemistryCheck.getBitFingerprints()));
        tmpComboBox.setConverter(new StringConverter<>() {
            @Override
            public String toString(Descriptor aDescriptor) {
                return aDescriptor == null ? "" : aDescriptor.getName();
            }
            @Override
            public Descriptor fromString(String aString) {
                return null;
            }
        });
        Label tmpHintLabel = new Label(Message.get("ClusterChemistryCheck.hint.text"));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setMaxWidth(480);
        Label tmpIndependenceLabel = new Label();
        tmpIndependenceLabel.setWrapText(true);
        tmpIndependenceLabel.setMaxWidth(480);
        //the warning names the problem of the chosen fingerprint
        Runnable tmpUpdateIndependence = () -> {
            Descriptor tmpFingerprint = tmpComboBox.getValue();
            if (tmpFingerprint == null) {
                tmpIndependenceLabel.setText("");
                return;
            }
            ClusterChemistryCheck.Independence tmpIndependence = ClusterChemistryCheck.getIndependence(tmpFingerprint, tmpParameters);
            tmpIndependenceLabel.setText(String.format(Message.get("ClusterChemistryCheck." + tmpIndependence.name().toLowerCase(Locale.ROOT) + ".text"),
                    tmpFingerprint.getName()));
            tmpIndependenceLabel.setStyle(tmpIndependence == ClusterChemistryCheck.Independence.INDEPENDENT ? "" : "-fx-text-fill: darkorange;");
        };
        tmpComboBox.valueProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpUpdateIndependence.run());
        tmpComboBox.setValue(ClusterChemistryCheck.getDefaultFingerprint(tmpParameters));
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("ClusterChemistryCheck.fingerprint.label")), tmpComboBox);
        tmpGrid.add(tmpIndependenceLabel, 0, 1, 2, 1);
        tmpGrid.add(tmpHintLabel, 0, 2, 2, 1);
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.getFocusedStage());
        tmpDialog.setTitle(Message.get("ClusterChemistryCheck.dialog.title"));
        tmpDialog.setHeaderText(Message.get("ClusterChemistryCheck.dialog.header"));
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return null;
        }
        return tmpComboBox.getValue();
    }
    //
    /**
     * Returns the text of the overall values of the given chemical check, with a note if the fingerprint was not
     * independent of the features and the number of molecules left out.
     *
     * @param aResult result of the chemical check
     * @return text for the quality window
     */
    private String formatChemistryCheckText(ClusterChemistryCheck.ChemistryCheckResult aResult) {
        String tmpText = String.format(Locale.getDefault(), Message.get("ClusterQuality.chemistry.text"),
                aResult.fingerprint().getName(), aResult.overallMeanTanimoto(), aResult.overallRandomMeanTanimoto(),
                aResult.overallRandomMeanTanimotoDeviation(), 100.0 * aResult.overallScaffoldShare(),
                100.0 * aResult.overallRandomScaffoldShare(), 100.0 * aResult.overallRandomScaffoldShareDeviation(),
                aResult.numberOfRandomRuns());
        if (ClusterChemistryCheck.getIndependence(aResult.fingerprint(), this.runResult.parameters()) != ClusterChemistryCheck.Independence.INDEPENDENT) {
            tmpText += " " + Message.get("ClusterQuality.chemistry.notIndependent.text");
        }
        if (aResult.numberOfFailedMolecules() > 0) {
            tmpText += " " + String.format(Message.get("ClusterQuality.chemistry.failed.text"), aResult.numberOfFailedMolecules());
        }
        return tmpText;
    }
    //
    /**
     * Calculates the PCA of this clustering result in the background (on the vectors the ART-2A clustering compared)
     * and shows it in a new tab as interactive scatter plot; a double click on a point or cluster centre opens the cluster.
     *
     * @param aVisibleClusterIndices clusters shown at the start, e.g. the marked ones; empty to show all
     */
    private void openPcaTab(Set<Integer> aVisibleClusterIndices) {
        if (this.clusteringViewController.isClusteringRunning()) {
            return;
        }
        if (this.runResult.molecules().size() < 2) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("PcaView.title.text"),
                    Message.get("PcaView.tooFewMolecules.header"), null);
            return;
        }
        Task<PrincipalComponentAnalysis.PcaResult> tmpTask = new Task<>() {
            @Override
            protected PrincipalComponentAnalysis.PcaResult call() {
                return PrincipalComponentAnalysis.calculate(ClusteringResultViewController.this.runResult,
                        PrincipalComponentAnalysis.NUMBER_OF_COMPONENTS_DEFAULT, this::isCancelled);
            }
        };
        this.clusteringViewController.startBackgroundTask(tmpTask, Message.get("Status.pca"), Message.get("PcaView.title.text"),
                Message.get("PcaView.error.header"), aPcaResult -> {
                    PcaView tmpView = new PcaView(this.runResult, aPcaResult, aVisibleClusterIndices);
                    Map<Integer, ClusteringRunResult.ClusterSummary> tmpClusterMap = this.getClusterMap();
                    tmpView.setClusterDoubleClickAction(aClusterIndex -> {
                        ClusteringRunResult.ClusterSummary tmpCluster = tmpClusterMap.get(aClusterIndex);
                        if (tmpCluster != null) {
                            this.openClusterOverview(tmpCluster, this.primaryStage);
                        }
                    });
                    Tab tmpTab = new Tab(Message.get("PcaView.tabTitle") + " - " + this.tab.getText());
                    tmpTab.setId(TabNames.PCA.toString());
                    tmpTab.setClosable(true);
                    tmpTab.setContent(tmpView);
                    this.mainTabPane.getTabs().add(tmpTab);
                    this.mainTabPane.getSelectionModel().select(tmpTab);
                }, null);
    }
    //
    /**
     * Returns the clusters of this result by cluster index.
     *
     * @return map from cluster index to cluster
     */
    private Map<Integer, ClusteringRunResult.ClusterSummary> getClusterMap() {
        Map<Integer, ClusteringRunResult.ClusterSummary> tmpClusterMap = new HashMap<>();
        this.runResult.getClusterSummaries().forEach(aCluster -> tmpClusterMap.put(aCluster.clusterIndex(), aCluster));
        return tmpClusterMap;
    }
    //</editor-fold>
    //
    //<editor-fold desc="subset and partitioning" defaultstate="collapsed">
    /**
     * Asks for the number of molecules per cluster, selects the diverse subset of this clustering result, shows it
     * in an overview view (cluster and rank below every structure) and optionally exports it as CSV.
     *
     * @param aClusterIndices indices of the marked clusters to select from; empty for all clusters
     */
    private void createDiverseSubset(Set<Integer> aClusterIndices) {
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.primaryStage);
        tmpDialog.setTitle(Message.get("DiverseSubset.dialog.title"));
        tmpDialog.setHeaderText(ClusteringResultViewController.appendMarkedClustersNote(Message.get("DiverseSubset.dialog.header"), aClusterIndices));
        Spinner<Integer> tmpNumberSpinner = new Spinner<>(1, Integer.MAX_VALUE, 1);
        tmpNumberSpinner.setEditable(true);
        ComboBox<DiverseSubsetSelector.SelectionMode> tmpModeComboBox = ClusteringResultViewController.createSelectionModeComboBox();
        CheckBox tmpUnclusteredCheckBox = new CheckBox(Message.get("DiverseSubset.unclustered.label"));
        //the not clusterable molecules belong to no marked cluster
        tmpUnclusteredCheckBox.setDisable(!aClusterIndices.isEmpty());
        CheckBox tmpExportCheckBox = new CheckBox(Message.get("DiverseSubset.export.label"));
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("DiverseSubset.number.label")), tmpNumberSpinner);
        tmpGrid.addRow(1, new Label(Message.get("DiverseSubset.selectionMode.label")), tmpModeComboBox);
        tmpGrid.add(tmpUnclusteredCheckBox, 0, 2, 2, 1);
        tmpGrid.add(tmpExportCheckBox, 0, 3, 2, 1);
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        List<DiverseSubsetSelector.SelectedMolecule> tmpSelection = DiverseSubsetSelector.select(this.runResult,
                tmpNumberSpinner.getValue(), tmpUnclusteredCheckBox.isSelected(), tmpModeComboBox.getValue(), aClusterIndices);
        if (tmpExportCheckBox.isSelected()) {
            Exporter tmpExporter = new Exporter(this.settingsContainer);
            File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage,
                    Exporter.ExportTypes.DIVERSE_SUBSET_CSV_FILE, this.runResult.parameters().fragmentationName());
            if (tmpDirectory != null) {
                ClusteringViewController.showExportOutcome(tmpExporter.exportDiverseSubsetCsvFile(tmpDirectory, this.runResult,
                        tmpSelection, this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
            }
        }
        if (tmpSelection.isEmpty()) {
            return;
        }
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(tmpSelection.size());
        List<String> tmpLabels = new ArrayList<>(tmpSelection.size());
        for (DiverseSubsetSelector.SelectedMolecule tmpSelected : tmpSelection) {
            tmpMolecules.add(this.runResult.molecules().get(tmpSelected.rowIndex()));
            tmpLabels.add(tmpSelected.clusterIndex() < 0
                    ? Message.get("ClusteringResultView.notClusterable.text")
                    : String.format(Message.get("DiverseSubset.structureLabel.text"), tmpSelected.clusterIndex(), tmpSelected.rankInCluster()));
        }
        OverviewViewController tmpController = this.createOverviewViewController();
        tmpController.initializeAndShowOverviewView(this.primaryStage, OverviewViewController.DataSources.DIVERSE_SUBSET,
                String.format(Message.get("DiverseSubset.overview.title"), tmpSelection.size()), tmpMolecules, tmpLabels, null, null);
    }
    //
    /**
     * Asks for the training/test ratio and the selection mode, creates the cluster-based training/test partitioning of
     * this clustering result and exports it as CSV files into a directory chosen by the user.
     *
     * @param aClusterIndices indices of the marked clusters whose molecules are partitioned; empty for all molecules
     */
    private void createTrainTestSplit(Set<Integer> aClusterIndices) {
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.primaryStage);
        tmpDialog.setTitle(Message.get("TrainTestSplit.dialog.title"));
        tmpDialog.setHeaderText(ClusteringResultViewController.appendMarkedClustersNote(Message.get("TrainTestSplit.dialog.header"), aClusterIndices));
        TextField tmpRatioTextField = new TextField(ClusteringResultViewController.TRAIN_TEST_RATIO_DEFAULT);
        tmpRatioTextField.setPromptText(ClusteringResultViewController.TRAIN_TEST_RATIO_DEFAULT);
        ComboBox<TrainTestSplitter.SelectionMode> tmpModeComboBox = new ComboBox<>(
                FXCollections.observableArrayList(TrainTestSplitter.SelectionMode.values()));
        tmpModeComboBox.setValue(TrainTestSplitter.SelectionMode.CLOSEST_TO_CLUSTER_CENTER);
        tmpModeComboBox.setConverter(new StringConverter<>() {
            @Override
            public String toString(TrainTestSplitter.SelectionMode aMode) {
                return aMode == null ? "" : Message.get("TrainTestSplit.selectionMode." + aMode.name());
            }
            @Override
            public TrainTestSplitter.SelectionMode fromString(String aString) {
                return null;
            }
        });
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("TrainTestSplit.ratio.label")), tmpRatioTextField);
        tmpGrid.addRow(1, new Label(Message.get("TrainTestSplit.selectionMode.label")), tmpModeComboBox);
        Label tmpRulesLabel = new Label(Message.get("TrainTestSplit.rules.text"));
        tmpRulesLabel.setWrapText(true);
        tmpRulesLabel.setMaxWidth(450);
        tmpGrid.add(tmpRulesLabel, 0, 2, 2, 1);
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        TrainTestSplitter.TrainTestSplit tmpSplit;
        try {
            tmpSplit = TrainTestSplitter.split(this.runResult, TrainTestSplitter.parseTrainingFraction(tmpRatioTextField.getText()),
                    tmpModeComboBox.getValue(), this.runResult.parameters().randomSeed(), aClusterIndices);
        } catch (IllegalArgumentException anException) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("TrainTestSplit.dialog.title"),
                    Message.get("TrainTestSplit.invalidRatio.header"), anException.getMessage());
            return;
        }
        Exporter tmpExporter = new Exporter(this.settingsContainer);
        File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage,
                Exporter.ExportTypes.TRAIN_TEST_SPLIT_CSV_FILES, this.runResult.parameters().fragmentationName());
        if (tmpDirectory == null) {
            return;
        }
        List<String> tmpFailedFiles = tmpExporter.exportTrainTestSplitCsvFiles(tmpDirectory, this.runResult, tmpSplit,
                this.settingsContainer.getCsvExportSeparatorSettingCharacter());
        if (tmpFailedFiles == null || !tmpFailedFiles.isEmpty()) {
            ClusteringViewController.showExportOutcome(tmpFailedFiles);
            return;
        }
        DecimalFormat tmpPercentFormat = new DecimalFormat("0.0", new DecimalFormatSymbols(Locale.getDefault()));
        GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("TrainTestSplit.dialog.title"),
                Message.get("TrainTestSplit.finished.header"),
                String.format(Message.get("TrainTestSplit.finished.content"), tmpSplit.trainingIndices().length,
                        tmpSplit.testIndices().length, tmpPercentFormat.format(100 * tmpSplit.getActualTrainingFraction()),
                        tmpPercentFormat.format(100 * tmpSplit.requestedTrainingFraction()), tmpDirectory.getAbsolutePath()));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns the row index of every molecule of the given clustering result, by identity of the molecule.
     *
     * @param aRunResult clustering result
     * @return map from molecule to its row index in the clustered data
     */
    private static Map<MoleculeDataModel, Integer> getRowIndexMap(ClusteringRunResult aRunResult) {
        Map<MoleculeDataModel, Integer> tmpRowIndexMap = new IdentityHashMap<>(aRunResult.molecules().size());
        for (int i = 0; i < aRunResult.molecules().size(); i++) {
            tmpRowIndexMap.put(aRunResult.molecules().get(i), i);
        }
        return tmpRowIndexMap;
    }
    //
    /**
     * Reads the conditions of the range filter from the given input fields; a blank bound means no bound, decimal point
     * and comma are accepted.
     *
     * @param aDescriptorComboBoxList descriptor of every condition
     * @param aMinimumFieldList lower bound of every condition
     * @param aMaximumFieldList upper bound of every condition
     * @return conditions or null if a bound is no number or a lower bound is greater than its upper bound
     */
    private static List<ClusterDescriptorStatistics.RangeCondition> readRangeConditions(List<ComboBox<String>> aDescriptorComboBoxList,
                                                                                        List<TextField> aMinimumFieldList,
                                                                                        List<TextField> aMaximumFieldList) {
        List<ClusterDescriptorStatistics.RangeCondition> tmpConditions = new ArrayList<>(aDescriptorComboBoxList.size());
        for (int i = 0; i < aDescriptorComboBoxList.size(); i++) {
            String tmpName = aDescriptorComboBoxList.get(i).getValue();
            if (tmpName == null) {
                continue;
            }
            double tmpMinimum;
            double tmpMaximum;
            try {
                tmpMinimum = ClusteringResultViewController.parseBound(aMinimumFieldList.get(i).getText());
                tmpMaximum = ClusteringResultViewController.parseBound(aMaximumFieldList.get(i).getText());
            } catch (NumberFormatException anException) {
                return null;
            }
            if (!Double.isNaN(tmpMinimum) && !Double.isNaN(tmpMaximum) && tmpMinimum > tmpMaximum) {
                return null;
            }
            tmpConditions.add(new ClusterDescriptorStatistics.RangeCondition(tmpName, tmpMinimum, tmpMaximum));
        }
        return tmpConditions;
    }
    //
    /**
     * Parses a bound of the range filter.
     *
     * @param aText text of the bound
     * @return bound, NaN if the text is blank
     * @throws NumberFormatException if the text is no number
     */
    private static double parseBound(String aText) throws NumberFormatException {
        if (aText == null || aText.isBlank()) {
            return Double.NaN;
        }
        double tmpValue = Double.parseDouble(aText.trim().replace(',', '.'));
        if (Double.isNaN(tmpValue)) {
            throw new NumberFormatException("NaN is no valid bound.");
        }
        return tmpValue;
    }
    //
    /**
     * Appends the note that only the marked clusters are used to the given dialog header, if clusters are marked.
     *
     * @param aHeader dialog header
     * @param aClusterIndices indices of the marked clusters, empty if none is marked
     * @return header with or without note
     */
    private static String appendMarkedClustersNote(String aHeader, Set<Integer> aClusterIndices) {
        return aClusterIndices.isEmpty() ? aHeader
                : aHeader + System.lineSeparator() + String.format(Message.get("ClusteringResultView.markedClustersUsed.text"), aClusterIndices.size());
    }
    //
    /**
     * Creates the combo box for the selection mode within a cluster (typical or diverse), preset to the default.
     *
     * @return combo box
     */
    private static ComboBox<DiverseSubsetSelector.SelectionMode> createSelectionModeComboBox() {
        ComboBox<DiverseSubsetSelector.SelectionMode> tmpComboBox = new ComboBox<>(
                FXCollections.observableArrayList(DiverseSubsetSelector.SelectionMode.values()));
        tmpComboBox.setValue(DiverseSubsetSelector.SELECTION_MODE_DEFAULT);
        tmpComboBox.setConverter(new StringConverter<>() {
            @Override
            public String toString(DiverseSubsetSelector.SelectionMode aMode) {
                return aMode == null ? "" : Message.get("DiverseSubset.selectionMode." + aMode.name());
            }
            @Override
            public DiverseSubsetSelector.SelectionMode fromString(String aString) {
                return null;
            }
        });
        return tmpComboBox;
    }
    //</editor-fold>
    //</editor-fold>
}
