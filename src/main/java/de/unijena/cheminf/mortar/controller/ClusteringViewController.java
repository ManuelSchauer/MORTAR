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
import de.unijena.cheminf.mortar.gui.controls.GridTabForTableView;
import de.unijena.cheminf.mortar.gui.controls.VigilanceSpinner;
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.gui.views.DescriptorTableView;
import de.unijena.cheminf.mortar.gui.views.MoleculesDataTableView;
import de.unijena.cheminf.mortar.gui.views.VigilanceScanView;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusteringParameters;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;
import de.unijena.cheminf.mortar.model.clustering.ClusteringTask;
import de.unijena.cheminf.mortar.model.clustering.DescriptorPriorityTable;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.VigilanceScan;
import de.unijena.cheminf.mortar.model.data.FragmentDataModel;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.descriptor.Descriptor;
import de.unijena.cheminf.mortar.model.fragmentation.FragmentationService;
import de.unijena.cheminf.mortar.model.io.Exporter;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;
import de.unijena.cheminf.mortar.model.util.LogUtil;

import javafx.application.Platform;
import javafx.beans.binding.Bindings;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ReadOnlyBooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.collections.ObservableMap;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.stage.Stage;
import javafx.stage.Window;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Controls the clustering in the main view: the clustering bars of the fragments, items and clustering result tabs,
 * the start and cancellation of a clustering, the vigilance scan, the descriptor calculation with its descriptor tabs
 * and the clustering settings. Every clustering result is shown in a tab that is controlled by its own
 * {@link ClusteringResultViewController}. All background tasks of the clustering (clustering, scan, descriptor
 * calculation, PCA, stability check) use the same controls: only one of them runs at a time, the start controls are
 * disabled and the cancel controls are visible while it runs.
 * <br>The main view controller keeps the menu, the status bar, the export of files and the clearing of the data and
 * creates this controller with the objects it shares.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class ClusteringViewController {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(ClusteringViewController.class.getName());
    /**
     * Number of clusters proposed in the dialog of the vigilance scan for the search of a target number of clusters.
     */
    private static final int TARGET_NUMBER_OF_CLUSTERS_DEFAULT = 20;
    /**
     * Prefix of the key of a vigilance scan started in the fragments or items tab of a fragmentation, followed by the
     * fragmentation name.
     */
    private static final String VIGILANCE_SCAN_KEY_FRAGMENTATION_PREFIX = "Fragmentation:";
    /**
     * Prefix of the key of a vigilance scan started in a clustering result tab, followed by the number of the tab.
     */
    private static final String VIGILANCE_SCAN_KEY_RESULT_PREFIX = "Result:";
    //</editor-fold>
    //
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * Primary stage of the application, owner of the dialogs and windows.
     */
    private final Stage primaryStage;
    /**
     * Tab pane of the main view, which gets the clustering result, descriptor and PCA tabs.
     */
    private final TabPane mainTabPane;
    /**
     * Configuration to read resource file paths from.
     */
    private final IConfiguration configuration;
    /**
     * Global settings, e.g. the number of tasks and the CSV separator.
     */
    private final SettingsContainer settingsContainer;
    /**
     * Service for the clustering of the molecules, holds the clustering settings.
     */
    private final ClusteringService clusteringService;
    /**
     * Fragmentation service, used to fragment the query of a similarity search and to get the last fragmentation.
     */
    private final FragmentationService fragmentationService;
    /**
     * View tools manager, whose overview views offer the similarity search and the descriptor values once there is a
     * clustering result.
     */
    private final ViewToolsManager viewToolsManager;
    /**
     * Imported molecules; the selected ones are clustered.
     */
    private final ObservableList<MoleculeDataModel> moleculeDataModelList;
    /**
     * Fragments of every fragmentation by fragmentation name, only read to check which fragmentations exist.
     */
    private final Map<String, ObservableList<FragmentDataModel>> mapOfFragmentDataModelLists;
    /**
     * Updates the status bar of the main view with the given thread and message.
     */
    private final BiConsumer<Thread, String> statusBarUpdater;
    /**
     * Exports the data of a clustering result (marked clusters, all if empty) with the export task of the main view.
     */
    private final BiConsumer<ClusteringRunResult, Set<Integer>> clusteringDataExporter;
    /**
     * True while a background task of the clustering is running; the start and cancel controls of the menu and the
     * tab bars are bound to it.
     */
    private final BooleanProperty isClusteringRunningProperty;
    /**
     * Window of the last vigilance scan of every clustering bar, by scan key (see the key prefix constants); a scan of a
     * fragmentation is kept until it is discarded or the data is cleared, a scan of a clustering result at the latest
     * until the result tab is closed. The show scan buttons of the bars observe this map.
     */
    private final ObservableMap<String, Stage> vigilanceScanStageMap;
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    /**
     * Running background task of the clustering, cancelled by the cancel controls; null if none was started yet.
     */
    private Task<?> runningTask;
    /**
     * Number of clustering result tabs created so far, used to give every tab a unique title.
     */
    private int clusteringRunCounter;
    /**
     * Number of descriptor tabs created so far, used to give every tab a unique title.
     */
    private int descriptorTabCounter;
    /**
     * Controller of the last clustering result, used by the similarity search and the descriptor values if no
     * clustering tab is selected; null if there was no clustering yet.
     */
    private ClusteringResultViewController lastResultViewController;
    /**
     * Overview view controller for the representatives of the clusters, shared by all clustering results; separate from
     * the one of the view tools manager, which opens a cluster on top of the representatives overview.
     */
    private OverviewViewController clusterRepresentativesOverviewViewController;
    //</editor-fold>
    //
    //<editor-fold desc="constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aPrimaryStage primary stage of the application
     * @param aMainTabPane tab pane of the main view
     * @param aConfiguration configuration to read resource file paths from
     * @param aSettingsContainer global settings
     * @param aClusteringService clustering service
     * @param aFragmentationService fragmentation service
     * @param aViewToolsManager view tools manager
     * @param aMoleculeDataModelList imported molecules
     * @param aMapOfFragmentDataModelLists fragments of every fragmentation by fragmentation name
     * @param aStatusBarUpdater updates the status bar with a thread and message
     * @param aClusteringDataExporter exports the data of a clustering result with the export task of the main view
     * @throws NullPointerException if a parameter is null
     */
    public ClusteringViewController(Stage aPrimaryStage, TabPane aMainTabPane, IConfiguration aConfiguration,
                                    SettingsContainer aSettingsContainer, ClusteringService aClusteringService,
                                    FragmentationService aFragmentationService, ViewToolsManager aViewToolsManager,
                                    ObservableList<MoleculeDataModel> aMoleculeDataModelList,
                                    Map<String, ObservableList<FragmentDataModel>> aMapOfFragmentDataModelLists,
                                    BiConsumer<Thread, String> aStatusBarUpdater,
                                    BiConsumer<ClusteringRunResult, Set<Integer>> aClusteringDataExporter)
            throws NullPointerException {
        this.primaryStage = Objects.requireNonNull(aPrimaryStage, "aPrimaryStage is null");
        this.mainTabPane = Objects.requireNonNull(aMainTabPane, "aMainTabPane is null");
        this.configuration = Objects.requireNonNull(aConfiguration, "aConfiguration is null");
        this.settingsContainer = Objects.requireNonNull(aSettingsContainer, "aSettingsContainer is null");
        this.clusteringService = Objects.requireNonNull(aClusteringService, "aClusteringService is null");
        this.fragmentationService = Objects.requireNonNull(aFragmentationService, "aFragmentationService is null");
        this.viewToolsManager = Objects.requireNonNull(aViewToolsManager, "aViewToolsManager is null");
        this.moleculeDataModelList = Objects.requireNonNull(aMoleculeDataModelList, "aMoleculeDataModelList is null");
        this.mapOfFragmentDataModelLists = Objects.requireNonNull(aMapOfFragmentDataModelLists, "aMapOfFragmentDataModelLists is null");
        this.statusBarUpdater = Objects.requireNonNull(aStatusBarUpdater, "aStatusBarUpdater is null");
        this.clusteringDataExporter = Objects.requireNonNull(aClusteringDataExporter, "aClusteringDataExporter is null");
        this.isClusteringRunningProperty = new SimpleBooleanProperty(false);
        this.vigilanceScanStageMap = FXCollections.observableHashMap();
        this.clusteringRunCounter = 0;
        this.descriptorTabCounter = 0;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the property that is true while a background task of the clustering is running.
     *
     * @return read-only running property
     */
    public ReadOnlyBooleanProperty isClusteringRunningProperty() {
        return this.isClusteringRunningProperty;
    }
    //
    /**
     * Returns whether a background task of the clustering is running.
     *
     * @return true if a clustering, scan, descriptor calculation, PCA or stability check is running
     */
    public boolean isClusteringRunning() {
        return this.isClusteringRunningProperty.get();
    }
    //
    /**
     * Starts the clustering of the selected molecules by the fragments of the fragmentation of the selected fragments or
     * items tab (or of the last fragmentation if another tab is selected) and the selected descriptors. The result is
     * shown in a new tab.
     */
    public void startClustering() {
        this.startClustering(null);
    }
    //
    /**
     * Cancels the running background task of the clustering. The current step (e.g. descriptor calculation or ART-2A
     * clustering) cannot be interrupted, the task stops after it and its result is discarded.
     */
    public void interruptClustering() {
        if (this.runningTask != null) {
            this.runningTask.cancel(true);
        }
        this.setClusteringRunningState(false);
    }
    //
    /**
     * Opens the clustering settings; their start button starts the clustering if a fragmentation is available.
     */
    public void openClusteringSettingsView() {
        this.openClusteringSettingsView(null);
    }
    //
    /**
     * Creates the clustering bar for an additional row below the pagination of a fragments or items tab, so that a
     * clustering or a vigilance scan can be started without opening the settings. A vigilance scan started here
     * calculates the features of the selected molecules once and is shared by the fragments and items tab of the
     * fragmentation.
     *
     * @param aFragmentationName name of the fragmentation of the tab
     * @return box with the clustering controls
     */
    public HBox createClusteringBar(String aFragmentationName) {
        return this.createClusteringBar(aFragmentationName, null,
                ClusteringViewController.VIGILANCE_SCAN_KEY_FRAGMENTATION_PREFIX + aFragmentationName);
    }
    //
    /**
     * Adds the clustering controls to the molecules tab: a bottom row with a button that calculates descriptors of the
     * selected molecules without fragmentation (the descriptor tab opened then clusters without fragment fingerprint)
     * and the context menu items for the similarity search and the descriptor values of the selected molecule, which
     * are available once there is a clustering result.
     *
     * @param aMoleculesTab molecules tab
     * @param aMoleculesDataTableView table view of the molecules tab
     */
    public void addClusteringControlsToMoleculesTab(GridTabForTableView aMoleculesTab, MoleculesDataTableView aMoleculesDataTableView) {
        aMoleculesTab.addNodeAsBottomRow(this.createMoleculesDescriptorBar());
        aMoleculesDataTableView.getContextMenu().setOnShowing(event -> {
            boolean tmpIsDisabled = this.getActiveResultViewController() == null
                    || aMoleculesDataTableView.getSelectionModel().getSelectedItem() == null;
            aMoleculesDataTableView.getSimilaritySearchMenuItem().setDisable(tmpIsDisabled);
            aMoleculesDataTableView.getDescriptorValuesMenuItem().setDisable(tmpIsDisabled);
        });
        aMoleculesDataTableView.getDescriptorValuesMenuItem().setOnAction(event -> {
            Object tmpSelectedItem = aMoleculesDataTableView.getSelectionModel().getSelectedItem();
            if (tmpSelectedItem instanceof MoleculeDataModel tmpMolecule) {
                this.showDescriptorValuesOfMoleculeInActiveResult(tmpMolecule);
            }
        });
        aMoleculesDataTableView.getSimilaritySearchMenuItem().setOnAction(event -> {
            Object tmpSelectedItem = aMoleculesDataTableView.getSelectionModel().getSelectedItem();
            if (tmpSelectedItem instanceof MoleculeDataModel tmpMolecule) {
                this.searchSimilarMoleculesInActiveResult(tmpMolecule);
            }
        });
    }
    //
    /**
     * Discards everything that belongs to the current data before new molecules are imported: cancels a running
     * background task, closes the kept vigilance scan windows (removing all tabs does not fire their closed events),
     * clears the cached descriptor values and forgets the last clustering result.
     */
    public void clear() {
        if (this.isClusteringRunning()) {
            this.interruptClustering();
        }
        this.vigilanceScanStageMap.values().forEach(Stage::close);
        this.vigilanceScanStageMap.clear();
        //the cached descriptor values belong to the molecules that are discarded now
        this.clusteringService.clearDescriptorMatrixCache();
        this.lastResultViewController = null;
    }
    //</editor-fold>
    //
    //<editor-fold desc="package-private methods" defaultstate="collapsed">
    /**
     * Runs the given background task of the clustering in a daemon thread with the clustering controls (start disabled,
     * cancel visible) and the status bar. Only one such task runs at a time, the callers check
     * {@link #isClusteringRunning()} before.
     *
     * @param aTask task to run
     * @param aStatusMessage message of the status bar while the task runs
     * @param anErrorTitle title of the alert if the task fails
     * @param anErrorHeader header of the alert if the task fails
     * @param aSucceededAction action with the value of the task when it has succeeded
     * @param aStoppedAction action when the task is cancelled or has failed, e.g. to reset a label; may be null
     * @param <T> type of the value of the task
     */
    <T> void startBackgroundTask(Task<T> aTask, String aStatusMessage, String anErrorTitle, String anErrorHeader,
                                 Consumer<T> aSucceededAction, Runnable aStoppedAction) {
        Thread tmpThread = new Thread(aTask);
        tmpThread.setName(MainViewController.ThreadType.CLUSTERING_THREAD.getThreadName());
        tmpThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        tmpThread.setDaemon(true);
        tmpThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number, as for the fragmentation
        aTask.setOnSucceeded(event -> {
            this.statusBarUpdater.accept(tmpThread, Message.get("Status.finished"));
            this.setClusteringRunningState(false);
            aSucceededAction.accept(aTask.getValue());
        });
        aTask.setOnCancelled(event -> {
            this.statusBarUpdater.accept(tmpThread, Message.get("Status.canceled"));
            this.setClusteringRunningState(false);
            if (aStoppedAction != null) {
                aStoppedAction.run();
            }
            ClusteringViewController.LOGGER.info("Background task of the clustering cancelled by user.");
        });
        aTask.setOnFailed(event -> {
            this.statusBarUpdater.accept(tmpThread, Message.get("Status.failed"));
            this.setClusteringRunningState(false);
            if (aStoppedAction != null) {
                aStoppedAction.run();
            }
            Throwable tmpCause = event.getSource().getException();
            ClusteringViewController.LOGGER.log(Level.SEVERE, tmpCause.toString(), tmpCause);
            GuiUtil.guiExceptionAlert(anErrorTitle, anErrorHeader, tmpCause.getMessage(),
                    tmpCause instanceof Exception tmpException ? tmpException : new Exception(tmpCause));
        });
        this.runningTask = aTask;
        this.setClusteringRunningState(true);
        this.statusBarUpdater.accept(tmpThread, aStatusMessage);
        tmpThread.start();
    }
    //
    /**
     * Exports the data of the given clustering result as CSV files with the export task of the main view.
     *
     * @param aRunResult clustering result
     * @param aClusterIndices indices of the marked clusters to export; empty for all molecules
     */
    void exportClusteringData(ClusteringRunResult aRunResult, Set<Integer> aClusterIndices) {
        this.clusteringDataExporter.accept(aRunResult, aClusterIndices);
    }
    //
    /**
     * Returns the overview view controller for the representatives of the clusters, shared by all clustering results;
     * created on first use.
     *
     * @return overview view controller of the cluster representatives
     */
    OverviewViewController getClusterRepresentativesOverviewViewController() {
        if (this.clusterRepresentativesOverviewViewController == null) {
            this.clusterRepresentativesOverviewViewController = new OverviewViewController(this.configuration, this.settingsContainer);
        }
        return this.clusterRepresentativesOverviewViewController;
    }
    //
    /**
     * Returns the focused window of the application, e.g. an overview view from whose context menu an action was
     * started, to use it as owner of dialogs and further windows; the given stage if no other window is focused.
     *
     * @param aDefaultStage stage returned if no window is focused, e.g. the main stage
     * @return focused stage or the given stage
     */
    static Stage getFocusedStage(Stage aDefaultStage) {
        for (Window tmpWindow : Window.getWindows()) {
            if (tmpWindow instanceof Stage tmpStage && tmpStage.isShowing() && tmpStage.isFocused()) {
                return tmpStage;
            }
        }
        return aDefaultStage;
    }
    //
    /**
     * Shows a warning if an export returned failed files or an invalid directory; otherwise nothing.
     *
     * @param aFailedFileList result of an export method of the exporter
     */
    static void showExportOutcome(List<String> aFailedFileList) {
        if (aFailedFileList == null || !aFailedFileList.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                    aFailedFileList == null ? null : String.join("\n", aFailedFileList));
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    //<editor-fold desc="start and settings" defaultstate="collapsed">
    /**
     * Starts the clustering of the selected molecules by the fragments of the given fragmentation and the selected
     * descriptors. The result is shown in a new tab.
     *
     * @param aFragmentationName name of the fragmentation to use; null to use the fragmentation of the selected fragments
     *                           or items tab (or the last fragmentation if another tab is selected)
     */
    private void startClustering(String aFragmentationName) {
        this.startClustering(aFragmentationName, null, null);
    }
    //
    /**
     * Starts the clustering of the given molecules by the fragments of the given fragmentation and the given
     * descriptors, e.g. those of a descriptor tab; descriptors already calculated for these molecules are taken from the
     * descriptor cache. The result is shown in a new tab.
     *
     * @param aFragmentationName name of the fragmentation to use; null to use the fragmentation of the selected fragments
     *                           or items tab (or the last fragmentation if another tab is selected);
     *                           {@link ClusteringService#WITHOUT_FRAGMENTATION_NAME} to cluster without fragment fingerprint
     * @param aMoleculeList molecules to cluster; null for the selected molecules
     * @param aNonFingerprintDescriptors descriptors to use instead of the selected ones (the selected fingerprints are
     *                                   added); null for the selected descriptors
     */
    private void startClustering(String aFragmentationName, List<MoleculeDataModel> aMoleculeList, Descriptor[] aNonFingerprintDescriptors) {
        if (this.isClusteringRunning()) {
            return;
        }
        //without fragmentation, e.g. from the molecules tab, only descriptors and fingerprints are the features
        boolean tmpIsWithoutFragmentation = ClusteringService.WITHOUT_FRAGMENTATION_NAME.equals(aFragmentationName);
        String tmpFragmentationName = tmpIsWithoutFragmentation ? aFragmentationName : this.getFragmentationNameToCluster(aFragmentationName);
        if (tmpFragmentationName == null) {
            return;
        }
        List<MoleculeDataModel> tmpSelectedMolecules = aMoleculeList == null ? this.getMoleculesToCluster() : aMoleculeList;
        if (tmpSelectedMolecules == null) {
            return;
        }
        boolean tmpIsParallel = this.settingsContainer.getNumberOfTasksForFragmentationSetting() > 1;
        ClusteringTask tmpClusteringTask;
        try {
            ClusteringParameters tmpParameters;
            if (tmpIsWithoutFragmentation) {
                tmpParameters = this.clusteringService.createClusteringParametersWithoutFragmentation(aNonFingerprintDescriptors, tmpIsParallel);
            } else if (aNonFingerprintDescriptors == null) {
                tmpParameters = this.clusteringService.createClusteringParameters(tmpFragmentationName, tmpIsParallel);
            } else {
                tmpParameters = this.clusteringService.createClusteringParameters(tmpFragmentationName, aNonFingerprintDescriptors, tmpIsParallel);
            }
            tmpClusteringTask = this.clusteringService.createClusteringTask(tmpSelectedMolecules, tmpParameters);
        } catch (IOException anException) {
            ClusteringViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("MainViewController.ClusteringError.Title"),
                    Message.get("MainViewController.ClusteringError.Header"),
                    Message.get("MainViewController.ClusteringError.Content"), anException);
            return;
        }
        this.startBackgroundTask(tmpClusteringTask, Message.get("Status.clustering"),
                Message.get("MainViewController.ClusteringError.Title"), Message.get("MainViewController.ClusteringError.Header"),
                this::addClusteringResultTab, null);
    }
    //
    /**
     * Returns the name of the fragmentation to cluster with: the given name or, if null, the fragmentation of the
     * selected fragments or items tab or the last fragmentation. If there is no such fragmentation, the user is informed.
     *
     * @param aFragmentationName name of the fragmentation, null to determine it from the selected tab
     * @return name of the fragmentation or null if there is none
     */
    private String getFragmentationNameToCluster(String aFragmentationName) {
        String tmpFragmentationName = aFragmentationName;
        Tab tmpSelectedTab = this.mainTabPane.getSelectionModel().getSelectedItem();
        if (tmpFragmentationName == null && tmpSelectedTab instanceof GridTabForTableView tmpGridTab
                && !TabNames.MOLECULES.toString().equals(tmpGridTab.getId())) {
            tmpFragmentationName = tmpGridTab.getFragmentationNameOutOfTitle();
        }
        if (tmpFragmentationName == null) {
            tmpFragmentationName = this.fragmentationService.getCurrentFragmentationName();
        }
        if (tmpFragmentationName == null || !this.mapOfFragmentDataModelLists.containsKey(tmpFragmentationName)) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("MainViewController.ClusteringError.Title"),
                    Message.get("MainViewController.Clustering.NoFragmentation.Header"),
                    Message.get("MainViewController.Clustering.NoFragmentation.Content"));
            return null;
        }
        return tmpFragmentationName;
    }
    //
    /**
     * Returns the selected molecules, which are clustered. If no molecule is selected, the user is informed.
     *
     * @return selected molecules or null if no molecule is selected
     */
    private List<MoleculeDataModel> getMoleculesToCluster() {
        List<MoleculeDataModel> tmpSelectedMolecules = this.moleculeDataModelList.stream().filter(MoleculeDataModel::isSelected).toList();
        if (tmpSelectedMolecules.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("MainViewController.ClusteringError.Title"),
                    Message.get("MainViewController.Clustering.NoMolecules.Header"), null);
            return null;
        }
        return tmpSelectedMolecules;
    }
    //
    /**
     * Sets whether a background task of the clustering is running; the start and cancel controls are bound to it.
     *
     * @param anIsRunning true if a background task of the clustering is running
     */
    private void setClusteringRunningState(boolean anIsRunning) {
        this.isClusteringRunningProperty.set(anIsRunning);
    }
    //
    /**
     * Opens the clustering settings; their start button starts the clustering with the given fragmentation.
     *
     * @param aFragmentationName name of the fragmentation to cluster with, null for the fragmentation of the selected tab
     */
    private void openClusteringSettingsView(String aFragmentationName) {
        //without fragmentation (clustering by descriptors only, or no fragmentation done yet), only the fragment fingerprint
        // settings are greyed out; fingerprints and descriptors stay available
        boolean tmpIsWithoutFragmentation = ClusteringService.WITHOUT_FRAGMENTATION_NAME.equals(aFragmentationName);
        boolean tmpIsFragmentFingerprintAvailable = !tmpIsWithoutFragmentation && !this.mapOfFragmentDataModelLists.isEmpty();
        new ClusteringSettingsViewController(this.primaryStage, this.clusteringService, this.configuration,
                (this.mapOfFragmentDataModelLists.isEmpty() && !tmpIsWithoutFragmentation) || this.isClusteringRunning()
                        ? null
                        : () -> this.startClustering(aFragmentationName), false, tmpIsFragmentFingerprintAvailable);
    }
    //
    /**
     * Creates the settings button with the gear icon that opens the clustering settings for the given fragmentation.
     *
     * @param aFragmentationName name of the fragmentation to cluster with
     * @return settings button
     */
    private Button createSettingsButton(String aFragmentationName) {
        Button tmpSettingsButton = new Button();
        String tmpIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder")
                        + this.configuration.getProperty("mortar.icon.gear.name")).toExternalForm();
        tmpSettingsButton.setGraphic(new ImageView(new Image(tmpIconURL)));
        tmpSettingsButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.settingsButton.tooltip")));
        tmpSettingsButton.setOnAction(anEvent -> this.openClusteringSettingsView(aFragmentationName));
        return tmpSettingsButton;
    }
    //
    /**
     * Creates a clustering bar that clusters with the given fragmentation: descriptor preset, settings button,
     * vigilance, start button, vigilance scan button, button to show the last scan again (only visible if there is
     * one) and cancel button; in the fragments and items tabs also a button on the right to calculate the descriptors
     * of the molecules without clustering.
     *
     * @param aFragmentationName name of the fragmentation to cluster with
     * @param aBaseRunResult clustering result whose features a vigilance scan reuses (in a clustering result tab); null
     *                       to calculate the features for the scan (in a fragments or items tab)
     * @param aScanKey key of the vigilance scan window of this bar in the scan window map
     * @return box with the clustering controls
     */
    private HBox createClusteringBar(String aFragmentationName, ClusteringRunResult aBaseRunResult, String aScanKey) {
        HBox tmpClusteringBar = new HBox();
        tmpClusteringBar.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        tmpClusteringBar.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpClusteringBar.setAlignment(Pos.CENTER_LEFT);
        //separated from the pagination above by a line
        tmpClusteringBar.setStyle("-fx-border-color: lightgrey; -fx-border-width: 1 0 0 0;");
        Label tmpPresetLabel = new Label(Message.get("MainTabPane.clusteringBar.presetLabel.text"));
        DescriptorPresetComboBox tmpPresetComboBox = new DescriptorPresetComboBox(this.clusteringService);
        Button tmpSettingsButton = this.createSettingsButton(aFragmentationName);
        //the buttons get the width of their text instead of the standard width, so that the text is not cut
        Button tmpStartButton = new Button(Message.get("MainTabPane.clusteringBar.startButton.text"));
        tmpStartButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpStartButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        tmpStartButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.startButton.tooltip")));
        tmpStartButton.disableProperty().bind(this.isClusteringRunningProperty);
        tmpStartButton.setOnAction(anEvent -> this.startClustering(aFragmentationName));
        //vigilance next to the start button, the scan helps to choose it
        Label tmpVigilanceLabel = new Label(Message.get("MainTabPane.clusteringBar.vigilanceLabel.text"));
        VigilanceSpinner tmpVigilanceSpinner = new VigilanceSpinner(this.clusteringService);
        Button tmpScanButton = new Button(Message.get("MainTabPane.clusteringBar.vigilanceScanButton.text"));
        tmpScanButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpScanButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        tmpScanButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.vigilanceScanButton.tooltip")));
        tmpScanButton.disableProperty().bind(this.isClusteringRunningProperty);
        tmpScanButton.setOnAction(anEvent -> this.startVigilanceScan(aFragmentationName, aBaseRunResult, aScanKey));
        Button tmpShowScanButton = new Button(Message.get("MainTabPane.clusteringBar.showVigilanceScanButton.text"));
        tmpShowScanButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpShowScanButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        tmpShowScanButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.showVigilanceScanButton.tooltip")));
        //only visible while a scan of this bar is kept
        tmpShowScanButton.visibleProperty().bind(Bindings.createBooleanBinding(
                () -> this.vigilanceScanStageMap.containsKey(aScanKey), this.vigilanceScanStageMap));
        tmpShowScanButton.managedProperty().bind(tmpShowScanButton.visibleProperty());
        tmpShowScanButton.setOnAction(anEvent -> {
            Stage tmpScanStage = this.vigilanceScanStageMap.get(aScanKey);
            if (tmpScanStage != null) {
                tmpScanStage.show();
                tmpScanStage.toFront();
            }
        });
        Button tmpCancelButton = new Button(Message.get("MainTabPane.clusteringBar.cancelButton.text"));
        tmpCancelButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpCancelButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        tmpCancelButton.visibleProperty().bind(this.isClusteringRunningProperty);
        tmpCancelButton.managedProperty().bind(this.isClusteringRunningProperty);
        tmpCancelButton.setOnAction(anEvent -> this.interruptClustering());
        tmpClusteringBar.getChildren().addAll(tmpPresetLabel, tmpPresetComboBox, tmpSettingsButton, tmpVigilanceLabel,
                tmpVigilanceSpinner, tmpStartButton, tmpScanButton, tmpShowScanButton, tmpCancelButton);
        //in the fragments and items tabs, the descriptors of the molecules can be calculated without clustering, on the right
        if (aBaseRunResult == null) {
            Region tmpSpacer = new Region();
            HBox.setHgrow(tmpSpacer, Priority.ALWAYS);
            Button tmpCalculateDescriptorsButton = new Button(Message.get("MainTabPane.clusteringBar.calculateDescriptorsButton.text"));
            tmpCalculateDescriptorsButton.setMinWidth(Region.USE_PREF_SIZE);
            tmpCalculateDescriptorsButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
            tmpCalculateDescriptorsButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.calculateDescriptorsButton.tooltip")));
            tmpCalculateDescriptorsButton.disableProperty().bind(this.isClusteringRunningProperty);
            tmpCalculateDescriptorsButton.setOnAction(anEvent -> this.openDescriptorCalculationView(aFragmentationName));
            tmpClusteringBar.getChildren().addAll(tmpSpacer, tmpCalculateDescriptorsButton);
        }
        return tmpClusteringBar;
    }
    //</editor-fold>
    //
    //<editor-fold desc="clustering results" defaultstate="collapsed">
    /**
     * Adds a tab that shows the given clustering result, controlled by a new {@link ClusteringResultViewController},
     * with a clustering bar to cluster again with this fragmentation and changed settings (its scan reuses the
     * features). The overview views offer the similarity search and the descriptor values from now on, and the
     * overview of the cluster representatives opens directly.
     *
     * @param aRunResult result of the clustering run
     */
    private void addClusteringResultTab(ClusteringRunResult aRunResult) {
        this.clusteringRunCounter++;
        String tmpTitle = Message.get("ClusteringResultView.tabTitle") + " " + this.clusteringRunCounter
                + " - " + aRunResult.parameters().fragmentationName();
        ClusteringResultViewController tmpController = new ClusteringResultViewController(this, this.primaryStage,
                this.mainTabPane, this.configuration, this.settingsContainer, this.clusteringService,
                this.fragmentationService, this.moleculeDataModelList, aRunResult, tmpTitle);
        //the scan window of this result is discarded at the latest when the result tab is closed
        String tmpScanKey = ClusteringViewController.VIGILANCE_SCAN_KEY_RESULT_PREFIX + this.clusteringRunCounter;
        tmpController.getView().setClusteringBar(this.createClusteringBar(aRunResult.parameters().fragmentationName(), aRunResult, tmpScanKey));
        Tab tmpTab = tmpController.getTab();
        tmpTab.setOnClosed(event -> {
            tmpController.close();
            Stage tmpScanStage = this.vigilanceScanStageMap.remove(tmpScanKey);
            if (tmpScanStage != null) {
                tmpScanStage.close();
            }
        });
        //the overview views offer the similarity search in their context menu from now on
        this.lastResultViewController = tmpController;
        this.viewToolsManager.setOverviewSimilaritySearchAction(this::searchSimilarMoleculesInActiveResult);
        this.viewToolsManager.setOverviewDescriptorValuesAction(this::showDescriptorValuesOfMoleculeInActiveResult);
        this.mainTabPane.getTabs().add(tmpTab);
        this.mainTabPane.getSelectionModel().select(tmpTab);
        //the overview of the cluster representatives opens directly, after the success handler has finished
        Platform.runLater(tmpController::openClusterRepresentativesOverview);
        if (aRunResult.clusteringResult().isOverflow()) {
            ClusteringViewController.LOGGER.warning("Clustering finished with cluster overflow.");
        }
    }
    //
    /**
     * Returns the controller of the clustering result the similarity search and the descriptor values refer to: the
     * one of the selected clustering tab, otherwise the one of the last clustering run.
     *
     * @return controller of the clustering result or null if there was no clustering yet
     */
    private ClusteringResultViewController getActiveResultViewController() {
        Tab tmpSelectedTab = this.mainTabPane.getSelectionModel().getSelectedItem();
        if (tmpSelectedTab != null && tmpSelectedTab.getUserData() instanceof ClusteringResultViewController tmpController) {
            return tmpController;
        }
        return this.lastResultViewController;
    }
    //
    /**
     * Searches the molecules most similar to the given molecule in the active clustering result, see
     * {@link #getActiveResultViewController()}; if there is none, the user is informed.
     *
     * @param aQuery query molecule
     */
    private void searchSimilarMoleculesInActiveResult(MoleculeDataModel aQuery) {
        ClusteringResultViewController tmpController = this.getActiveResultViewController();
        if (tmpController == null) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringSimilaritySearch.error.title"),
                    Message.get("ClusteringSimilaritySearch.noClustering.header"), null);
            return;
        }
        tmpController.searchSimilarMolecules(aQuery);
    }
    //
    /**
     * Shows the descriptor values of the given molecule in the active clustering result, see
     * {@link #getActiveResultViewController()}; if there is none, the user is informed.
     *
     * @param aMolecule molecule
     */
    private void showDescriptorValuesOfMoleculeInActiveResult(MoleculeDataModel aMolecule) {
        ClusteringResultViewController tmpController = this.getActiveResultViewController();
        if (tmpController == null) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.noClustering.header"), null);
            return;
        }
        tmpController.showDescriptorValuesOfMolecule(aMolecule);
    }
    //</editor-fold>
    //
    //<editor-fold desc="descriptor calculation" defaultstate="collapsed">
    /**
     * Creates the bar for an additional row below the pagination of the molecules tab with a button on the right that
     * calculates descriptors of the selected molecules without fragmentation; the descriptor tab opened then clusters
     * without fragment fingerprint.
     *
     * @return box with the descriptor calculation button
     */
    private HBox createMoleculesDescriptorBar() {
        HBox tmpBar = new HBox();
        tmpBar.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        tmpBar.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpBar.setAlignment(Pos.CENTER_RIGHT);
        //separated from the pagination above by a line, like the clustering bars
        tmpBar.setStyle("-fx-border-color: lightgrey; -fx-border-width: 1 0 0 0;");
        Button tmpCalculateDescriptorsButton = new Button(Message.get("MainTabPane.clusteringBar.calculateDescriptorsButton.text"));
        tmpCalculateDescriptorsButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpCalculateDescriptorsButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        tmpCalculateDescriptorsButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.moleculesTab.calculateDescriptorsButton.tooltip")));
        tmpCalculateDescriptorsButton.disableProperty().bind(this.isClusteringRunningProperty);
        tmpCalculateDescriptorsButton.setOnAction(anEvent -> this.openDescriptorCalculationView(ClusteringService.WITHOUT_FRAGMENTATION_NAME));
        tmpBar.getChildren().add(tmpCalculateDescriptorsButton);
        return tmpBar;
    }
    //
    /**
     * Opens the descriptor tab of the clustering settings for the choice of the descriptors (presets and single
     * descriptors, the same selection as for the clustering); its button starts the descriptor calculation.
     *
     * @param aFragmentationName fragmentation of the tab the calculation was started from, used by the start clustering
     *                           button of the descriptor tab
     */
    private void openDescriptorCalculationView(String aFragmentationName) {
        new ClusteringSettingsViewController(this.primaryStage, this.clusteringService, this.configuration,
                this.isClusteringRunning() ? null : () -> this.startDescriptorCalculation(aFragmentationName), true);
    }
    //
    /**
     * Calculates the selected descriptors (without fingerprints) of the selected molecules in the background, with the
     * aromaticity model of the clustering settings, and shows them in a new descriptor tab. The values are cached, so that
     * a following clustering of these molecules does not calculate them again.
     *
     * @param aFragmentationName fragmentation of the tab the calculation was started from, used by the start clustering
     *                           button of the descriptor tab
     */
    private void startDescriptorCalculation(String aFragmentationName) {
        if (this.isClusteringRunning()) {
            return;
        }
        List<MoleculeDataModel> tmpMolecules = this.getMoleculesToCluster();
        if (tmpMolecules == null) {
            return;
        }
        //the switched on descriptors and fingerprints, e.g. only fingerprints
        Descriptor[] tmpDescriptors = this.clusteringService.getSelectedDescriptors();
        if (tmpDescriptors.length == 0) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("DescriptorTableView.title.text"),
                    Message.get("DescriptorTableView.noDescriptors.header"), null);
            return;
        }
        boolean tmpIsParallel = this.settingsContainer.getNumberOfTasksForFragmentationSetting() > 1;
        Task<FeatureMatrix> tmpTask = new Task<>() {
            @Override
            protected FeatureMatrix call() throws Exception {
                //cached values are reused and the calculated ones are cached for a following clustering
                return ClusteringViewController.this.clusteringService.calculateDescriptorMatrix(tmpMolecules, tmpDescriptors, tmpIsParallel);
            }
        };
        this.startBackgroundTask(tmpTask, Message.get("Status.descriptorCalculation"), Message.get("DescriptorTableView.title.text"),
                Message.get("DescriptorTableView.error.header"),
                aMatrix -> this.addDescriptorTab(aFragmentationName, tmpMolecules, tmpDescriptors, aMatrix), null);
    }
    //
    /**
     * Adds a tab that shows the given descriptor values of the given molecules and connects its export button. Its
     * clustering controls (vigilance and start button) cluster exactly these molecules with exactly these descriptors
     * and the given fragmentation; the descriptors are taken from the cache and not calculated again.
     *
     * @param aFragmentationName fragmentation to cluster with
     * @param aMoleculeList molecules in the order of the matrix rows
     * @param aDescriptors calculated descriptors
     * @param aDescriptorMatrix descriptor values, one row per molecule
     */
    private void addDescriptorTab(String aFragmentationName, List<MoleculeDataModel> aMoleculeList, Descriptor[] aDescriptors,
                                  FeatureMatrix aDescriptorMatrix) {
        this.descriptorTabCounter++;
        int tmpTabNumber = this.descriptorTabCounter;
        String tmpTabName = Message.get("DescriptorTableView.tabTitle") + " " + tmpTabNumber;
        DescriptorTableView tmpView = new DescriptorTableView(aMoleculeList, aDescriptorMatrix);
        //clustering of exactly these molecules and descriptors
        Button tmpStartButton = new Button(Message.get("DescriptorTableView.startClusteringButton.text"));
        tmpStartButton.setMinWidth(Region.USE_PREF_SIZE);
        tmpStartButton.setTooltip(GuiUtil.createTooltip(ClusteringService.WITHOUT_FRAGMENTATION_NAME.equals(aFragmentationName)
                ? Message.get("DescriptorTableView.startClusteringButton.withoutFragmentation.tooltip")
                : String.format(Message.get("DescriptorTableView.startClusteringButton.tooltip"), aFragmentationName)));
        tmpStartButton.disableProperty().bind(this.isClusteringRunningProperty);
        tmpStartButton.setOnAction(event -> this.startClustering(aFragmentationName, aMoleculeList, aDescriptors));
        //settings, e.g. to choose the fingerprints that are calculated additionally for the clustering
        HBox tmpClusteringBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, this.createSettingsButton(aFragmentationName),
                new Label(Message.get("MainTabPane.clusteringBar.vigilanceLabel.text")), new VigilanceSpinner(this.clusteringService),
                tmpStartButton);
        tmpClusteringBox.setAlignment(Pos.CENTER_LEFT);
        tmpView.setClusteringControls(tmpClusteringBox);
        Tab tmpTab = new Tab(tmpTabName);
        tmpTab.setId(TabNames.DESCRIPTORS.toString());
        tmpTab.setClosable(true);
        tmpTab.setContent(tmpView);
        tmpView.getExportButton().setOnAction(event -> {
            Exporter tmpExporter = new Exporter(this.settingsContainer);
            File tmpFile = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage, Exporter.ExportTypes.DESCRIPTOR_CSV_FILE,
                    Integer.toString(tmpTabNumber));
            if (tmpFile != null) {
                ClusteringViewController.showExportOutcome(tmpExporter.exportDescriptorCsvFile(tmpFile, tmpView.getMoleculeList(),
                        tmpView.getDescriptorMatrix(), this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
            }
        });
        this.mainTabPane.getTabs().add(tmpTab);
        this.mainTabPane.getSelectionModel().select(tmpTab);
    }
    //</editor-fold>
    //
    //<editor-fold desc="vigilance scan" defaultstate="collapsed">
    /**
     * Asks for the vigilance range or a target number of clusters and starts a vigilance scan or the search for the target
     * (see {@link VigilanceScan#searchForNumberOfClusters}) in the background; the result is shown in a new window, for a
     * search with the best clustering highlighted. With a given clustering result, its features and parameters are
     * reused and only the clustering is repeated. Without one (scan before the first clustering), the features of the
     * selected molecules are calculated once with the current clustering settings.
     *
     * @param aFragmentationName name of the fragmentation to cluster with
     * @param aBaseRunResult clustering result whose features are reused; null to calculate the features
     * @param aScanKey key of the scan window in the scan window map; a kept scan with this key is replaced
     */
    private void startVigilanceScan(String aFragmentationName, ClusteringRunResult aBaseRunResult, String aScanKey) {
        if (this.isClusteringRunning()) {
            return;
        }
        List<MoleculeDataModel> tmpMolecules = null;
        ClusteringParameters tmpParameters = null;
        DescriptorPriorityTable tmpPriorityTable = null;
        if (aBaseRunResult == null) {
            if (this.getFragmentationNameToCluster(aFragmentationName) == null) {
                return;
            }
            tmpMolecules = this.getMoleculesToCluster();
            if (tmpMolecules == null) {
                return;
            }
        }
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.primaryStage);
        tmpDialog.setTitle(Message.get("VigilanceScan.dialog.title"));
        tmpDialog.setHeaderText(Message.get(aBaseRunResult == null ? "VigilanceScan.dialog.header.new" : "VigilanceScan.dialog.header"));
        //either a range of vigilance values or a target number of clusters
        ToggleGroup tmpModeGroup = new ToggleGroup();
        RadioButton tmpRangeRadioButton = new RadioButton(Message.get("VigilanceScan.rangeMode.label"));
        tmpRangeRadioButton.setToggleGroup(tmpModeGroup);
        tmpRangeRadioButton.setSelected(true);
        RadioButton tmpTargetRadioButton = new RadioButton(Message.get("VigilanceScan.targetMode.label"));
        tmpTargetRadioButton.setToggleGroup(tmpModeGroup);
        tmpTargetRadioButton.setTooltip(GuiUtil.createTooltip(Message.get("VigilanceScan.targetMode.tooltip")));
        TextField tmpStartField = new TextField("0.1");
        TextField tmpEndField = new TextField("0.9");
        TextField tmpStepField = new TextField("0.1");
        int tmpNumberOfMolecules = aBaseRunResult == null ? tmpMolecules.size() : aBaseRunResult.molecules().size();
        Spinner<Integer> tmpTargetSpinner = new Spinner<>(1, Math.max(1, tmpNumberOfMolecules),
                Math.min(ClusteringViewController.TARGET_NUMBER_OF_CLUSTERS_DEFAULT, Math.max(1, tmpNumberOfMolecules)));
        tmpTargetSpinner.setEditable(true);
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.add(tmpRangeRadioButton, 0, 0, 2, 1);
        Label tmpStartLabel = new Label(Message.get("VigilanceScan.start.label"));
        Label tmpEndLabel = new Label(Message.get("VigilanceScan.end.label"));
        Label tmpStepLabel = new Label(Message.get("VigilanceScan.step.label"));
        tmpGrid.addRow(1, tmpStartLabel, tmpStartField);
        tmpGrid.addRow(2, tmpEndLabel, tmpEndField);
        tmpGrid.addRow(3, tmpStepLabel, tmpStepField);
        tmpGrid.add(tmpTargetRadioButton, 0, 4, 2, 1);
        Label tmpTargetLabel = new Label(Message.get("VigilanceScan.target.label"));
        tmpGrid.addRow(5, tmpTargetLabel, tmpTargetSpinner);
        //only the fields of the chosen mode are enabled
        for (Node tmpNode : List.of(tmpStartLabel, tmpStartField, tmpEndLabel, tmpEndField, tmpStepLabel, tmpStepField)) {
            tmpNode.disableProperty().bind(tmpTargetRadioButton.selectedProperty());
        }
        tmpTargetLabel.disableProperty().bind(tmpRangeRadioButton.selectedProperty());
        tmpTargetSpinner.disableProperty().bind(tmpRangeRadioButton.selectedProperty());
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        boolean tmpIsTargetSearch = tmpTargetRadioButton.isSelected();
        float[] tmpVigilanceValues = null;
        int tmpTargetNumberOfClusters = 0;
        try {
            if (tmpIsTargetSearch) {
                //a typed value is only taken over by the spinner with enter, so it is read from the editor
                tmpTargetNumberOfClusters = Integer.parseInt(tmpTargetSpinner.getEditor().getText().trim());
                if (tmpTargetNumberOfClusters < 1) {
                    throw new IllegalArgumentException("The target number of clusters must be at least 1.");
                }
            } else {
                tmpVigilanceValues = VigilanceScan.getVigilanceValues(
                        Double.parseDouble(tmpStartField.getText().trim().replace(',', '.')),
                        Double.parseDouble(tmpEndField.getText().trim().replace(',', '.')),
                        Double.parseDouble(tmpStepField.getText().trim().replace(',', '.')));
            }
        } catch (IllegalArgumentException anException) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("VigilanceScan.dialog.title"),
                    Message.get("VigilanceScan.invalidInput.header"), anException.getMessage());
            return;
        }
        if (aBaseRunResult == null) {
            //the settings are read now, so that changes during the scan do not affect it
            tmpParameters = this.clusteringService.createClusteringParameters(aFragmentationName,
                    this.settingsContainer.getNumberOfTasksForFragmentationSetting() > 1);
            try {
                tmpPriorityTable = this.clusteringService.getPriorityTable();
            } catch (IOException anException) {
                ClusteringViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
                GuiUtil.guiExceptionAlert(Message.get("VigilanceScan.dialog.title"), Message.get("MainViewController.ClusteringError.Header"),
                        Message.get("MainViewController.ClusteringError.Content"), anException);
                return;
            }
        }
        List<MoleculeDataModel> tmpFinalMolecules = tmpMolecules;
        ClusteringParameters tmpFinalParameters = tmpParameters;
        DescriptorPriorityTable tmpFinalPriorityTable = tmpPriorityTable;
        float[] tmpFinalVigilanceValues = tmpVigilanceValues;
        int tmpFinalTargetNumberOfClusters = tmpTargetNumberOfClusters;
        //holds the clustering result the scan is based on, calculated in the task if not given
        ClusteringRunResult[] tmpBaseRunResult = new ClusteringRunResult[] {aBaseRunResult};
        //holds the result of a target search, null for a range scan
        VigilanceScan.TargetSearchResult[] tmpSearchResult = new VigilanceScan.TargetSearchResult[1];
        Task<List<VigilanceScan.ScanEntry>> tmpScanTask = new Task<>() {
            @Override
            protected List<VigilanceScan.ScanEntry> call() throws Exception {
                if (tmpBaseRunResult[0] == null) {
                    //features once, clustered with the current vigilance; only the clustering is repeated in the scan
                    tmpBaseRunResult[0] = ClusteringTask.runPipeline(tmpFinalMolecules, tmpFinalParameters, tmpFinalPriorityTable,
                            ClusteringViewController.this.clusteringService.getDescriptorMatrixCache(),
                            this::isCancelled, (aStep, aMessage) -> this.updateMessage(aMessage));
                }
                if (tmpIsTargetSearch) {
                    tmpSearchResult[0] = VigilanceScan.searchForNumberOfClusters(tmpBaseRunResult[0], tmpFinalTargetNumberOfClusters,
                            this::isCancelled, (aDone, aTotal) -> this.updateProgress(aDone, aTotal));
                    return tmpSearchResult[0].entries();
                }
                return VigilanceScan.scan(tmpBaseRunResult[0], tmpFinalVigilanceValues, this::isCancelled,
                        (aDone, aTotal) -> this.updateProgress(aDone, aTotal));
            }
        };
        this.startBackgroundTask(tmpScanTask, Message.get("Status.clustering"), Message.get("VigilanceScan.dialog.title"),
                Message.get("MainViewController.ClusteringError.Header"),
                anEntryList -> this.showVigilanceScanResult(aScanKey, tmpBaseRunResult[0], anEntryList, aBaseRunResult != null,
                        tmpSearchResult[0]), null);
    }
    //
    /**
     * Shows the result of a vigilance scan in a new window with chart and table. "Use This Clustering" (or a double
     * click on a row) opens the clustering of the selected vigilance as new result tab and uses its vigilance for the
     * next clusterings. "Close" and the window close button only hide the window, the scan is kept and can be shown
     * again with the show scan button of the clustering bar; "Discard Scan" closes the window and discards the scan.
     *
     * @param aScanKey key of the scan window in the scan window map; a kept window with this key is replaced
     * @param aBaseRunResult clustering result the scan is based on
     * @param anEntryList entries of the scan
     * @param anIsSelectingBaseVigilance true to select the row of the vigilance of the base result, i.e. of the
     *                                   clustering shown in the result tab the scan was started from
     * @param aTargetSearchResult result of a search for a target number of clusters, whose best entry is highlighted and
     *                            selected instead; null for a range scan
     */
    private void showVigilanceScanResult(String aScanKey, ClusteringRunResult aBaseRunResult, List<VigilanceScan.ScanEntry> anEntryList,
                                         boolean anIsSelectingBaseVigilance, VigilanceScan.TargetSearchResult aTargetSearchResult) {
        VigilanceScanView tmpView = new VigilanceScanView(anEntryList);
        if (aTargetSearchResult != null) {
            tmpView.showTargetSearchResult(aTargetSearchResult);
        } else if (anIsSelectingBaseVigilance) {
            tmpView.selectVigilance(aBaseRunResult.clusteringResult().vigilance());
        }
        Stage tmpOldStage = this.vigilanceScanStageMap.remove(aScanKey);
        if (tmpOldStage != null) {
            tmpOldStage.close();
        }
        Stage tmpStage = new Stage();
        tmpStage.initOwner(this.primaryStage);
        tmpStage.setTitle(String.format(Message.get("VigilanceScan.window.title"), aBaseRunResult.parameters().fragmentationName()));
        tmpStage.setScene(new Scene(tmpView, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE, GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE));
        Runnable tmpAdoptAction = () -> {
            VigilanceScan.ScanEntry tmpEntry = tmpView.getScanTableView().getSelectionModel().getSelectedItem();
            if (tmpEntry == null) {
                return;
            }
            //the float is converted via its decimal representation, e.g. 0.3 instead of 0.30000001192092896
            this.clusteringService.setVigilance(ClusteringService.toDecimalDouble(tmpEntry.vigilance()));
            this.addClusteringResultTab(tmpEntry.runResult());
        };
        tmpView.getAdoptButton().setOnAction(event -> tmpAdoptAction.run());
        tmpView.getScanTableView().setOnMouseClicked(event -> {
            if (event.getClickCount() == 2) {
                tmpAdoptAction.run();
            }
        });
        tmpView.getExportButton().setOnAction(event -> {
            Exporter tmpExporter = new Exporter(this.settingsContainer);
            File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(tmpStage, Exporter.ExportTypes.VIGILANCE_SCAN_CSV_FILE,
                    aBaseRunResult.parameters().fragmentationName());
            if (tmpDirectory != null) {
                ClusteringViewController.showExportOutcome(tmpExporter.exportVigilanceScanCsvFile(tmpDirectory, aBaseRunResult,
                        anEntryList, this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
            }
        });
        //closing only hides the window, the scan is kept
        tmpView.getCloseButton().setOnAction(event -> tmpStage.hide());
        tmpView.getDiscardButton().setOnAction(event -> {
            //only removed if it is still this window, not a newer scan with the same key
            this.vigilanceScanStageMap.remove(aScanKey, tmpStage);
            tmpStage.close();
        });
        this.vigilanceScanStageMap.put(aScanKey, tmpStage);
        tmpStage.show();
    }
    //</editor-fold>
    //</editor-fold>
}
