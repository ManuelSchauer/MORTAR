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
import de.unijena.cheminf.mortar.gui.controls.CustomPaginationSkin;
import de.unijena.cheminf.mortar.gui.controls.DescriptorPresetComboBox;
import de.unijena.cheminf.mortar.gui.controls.GridTabForTableView;
import de.unijena.cheminf.mortar.gui.controls.VigilanceSpinner;
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.gui.views.ClusteringResultView;
import de.unijena.cheminf.mortar.gui.views.FragmentsDataTableView;
import de.unijena.cheminf.mortar.gui.views.IDataTableView;
import de.unijena.cheminf.mortar.gui.views.ItemizationDataTableView;
import de.unijena.cheminf.mortar.gui.views.MainView;
import de.unijena.cheminf.mortar.gui.views.MoleculesDataTableView;
import de.unijena.cheminf.mortar.gui.views.VigilanceScanView;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusterDescriptorStatistics;
import de.unijena.cheminf.mortar.model.clustering.ClusteringParameters;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;
import de.unijena.cheminf.mortar.model.clustering.ClusteringTask;
import de.unijena.cheminf.mortar.model.clustering.DescriptorPriorityTable;
import de.unijena.cheminf.mortar.model.clustering.DiverseSubsetSelector;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrixBuilder;
import de.unijena.cheminf.mortar.model.clustering.SimilaritySearch;
import de.unijena.cheminf.mortar.model.clustering.TrainTestSplitter;
import de.unijena.cheminf.mortar.model.clustering.VigilanceScan;
import de.unijena.cheminf.mortar.model.data.FragmentDataModel;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.fragmentation.FragmentationService;
import de.unijena.cheminf.mortar.model.fragmentation.algorithm.IMoleculeFragmenter;
import de.unijena.cheminf.mortar.model.io.ChemFileTypes;
import de.unijena.cheminf.mortar.model.io.Exporter;
import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.ChemUtil;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;
import de.unijena.cheminf.mortar.model.util.FileUtil;
import de.unijena.cheminf.mortar.model.util.LogUtil;

import javafx.application.Platform;
import javafx.beans.Observable;
import javafx.beans.binding.Bindings;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.value.ChangeListener;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.collections.ObservableMap;
import javafx.collections.transformation.FilteredList;
import javafx.concurrent.Task;
import javafx.event.EventHandler;
import javafx.event.EventType;
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
import javafx.scene.control.Pagination;
import javafx.scene.control.RadioMenuItem;
import javafx.scene.control.SortEvent;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.Tooltip;
import javafx.scene.control.cell.CheckBoxListCell;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.Dragboard;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;
import javafx.scene.input.TransferMode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import javafx.stage.Window;
import javafx.stage.WindowEvent;
import javafx.util.StringConverter;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.File;
import java.io.IOException;
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
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ExecutionException;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * MainViewController controls  {@link MainView}.
 *
 * @author Felix Baensch
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class MainViewController {
    //<editor-fold desc="private and private final class variables" defaultstate="collapsed">
    /**
     * Primary Stage.
     */
    private final Stage primaryStage;
    /**
     * MainView.
     */
    private final MainView mainView;
    /**
     * Scene.
     */
    private final Scene scene;
    /**
     * TabPane which holds the different tabs.
     */
    private final TabPane mainTabPane;
    /**
     * ObservableList to hold MoleculeDataModels for visualisation in MoleculesDataTableView.
     */
    private final ObservableList<MoleculeDataModel> moleculeDataModelList;
    /**
     * MoleculesDataTableView to show imported molecules.
     */
    private MoleculesDataTableView moleculesDataTableView;
    /**
     * SettingsContainer.
     */
    private final SettingsContainer settingsContainer;
    /**
     * FragmentationService.
     */
    private final FragmentationService fragmentationService;
    /**
     * ViewToolsManager.
     */
    private final ViewToolsManager viewToolsManager;
    /**
     * Button to start single algorithm fragmentation.
     */
    private Button fragmentationButton;
    /**
     * Button to cancel running fragmentation.
     */
    private Button cancelFragmentationButton;
    /**
     * HashMap to hold Lists of FragmentDataModels for each fragmentation.
     */
    private final HashMap<String, ObservableList<FragmentDataModel>> mapOfFragmentDataModelLists;
    /**
     * Boolean value whether fragmentation is running.
     */
    private boolean isFragmentationRunning;
    /**
     * Task for parallel fragmentation.
     */
    private Task<Void> parallelFragmentationMainTask;
    /**
     * Thread for task for parallel fragmentation.
     */
    private Thread fragmentationThread;
    /**
     * Thread for molecule imports, so GUI thread is always responsive.
     */
    private Thread importerThread;
    /**
     * Task for molecule file import.
     */
    private Task<List<MoleculeDataModel>> importTask;
    /**
     * Storing the name of the last imported file.
     */
    private String importedFileName;
    /**
     * Thread for molecule exports, so GUI thread is always responsive.
     */
    private Thread exporterThread;
    /**
     * Task for molecule file export.
     */
    private Task<List<String>> exportTask;
    /**
     * Service for the clustering of the molecules, holds the clustering settings.
     */
    private final ClusteringService clusteringService;
    /**
     * Task of a running clustering.
     */
    private ClusteringTask clusteringTask;
    /**
     * Task of a running vigilance scan, cancelled by the cancel controls of the clustering.
     */
    private Task<List<VigilanceScan.ScanEntry>> vigilanceScanTask;
    /**
     * Thread of a running clustering.
     */
    private Thread clusteringThread;
    /**
     * Number of clustering result tabs created so far, used to give every tab a unique title.
     */
    private int clusteringRunCounter;
    /**
     * True while a clustering is running; the start and cancel controls of the menu and the tab bars are bound to it.
     */
    private final BooleanProperty isClusteringRunningProperty;
    /**
     * Overview view controller for the representatives of the clusters; separate from the one of the view tools manager,
     * which opens a cluster on top of the representatives overview.
     */
    private OverviewViewController clusterRepresentativesOverviewViewController;
    /**
     * Result of the last clustering run, used by the similarity search if no clustering tab is selected; null if there
     * was no clustering yet.
     */
    private ClusteringRunResult lastClusteringRunResult;
    /**
     * Similarity searches of the clustering runs, created on first use.
     */
    private final Map<ClusteringRunResult, SimilaritySearch> similaritySearchMap = new IdentityHashMap<>();
    /**
     * Window of the last vigilance scan of every clustering bar, by scan key (see the key prefix constants); a scan of a
     * fragmentation is kept until it is discarded or the data is cleared, a scan of a clustering result at the latest
     * until the result tab is closed. The show scan buttons of the bars observe this map.
     */
    private final ObservableMap<String, Stage> vigilanceScanStageMap = FXCollections.observableHashMap();
    /**
     * Training/test ratio proposed in the dialog of the cluster-based training/test partitioning.
     */
    private static final String TRAIN_TEST_RATIO_DEFAULT = "80:20";
    /**
     * Number of molecules proposed in the dialog for selecting molecules from a cluster.
     */
    private static final int CLUSTER_MOLECULE_SELECTION_NUMBER_DEFAULT = 10;
    /**
     * Prefix of the key of a vigilance scan started in the fragments or items tab of a fragmentation, followed by the
     * fragmentation name.
     */
    private static final String VIGILANCE_SCAN_KEY_FRAGMENTATION_PREFIX = "Fragmentation:";
    /**
     * Prefix of the key of a vigilance scan started in a clustering result tab, followed by the number of the tab.
     */
    private static final String VIGILANCE_SCAN_KEY_RESULT_PREFIX = "Result:";
    /**
     * BooleanProperty whether import is running.
     */
    private final BooleanProperty isImportRunningProperty;
    /**
     * BooleanProperty whether export is running.
     */
    private final BooleanProperty isExportRunningProperty;
    /**
     * Thread safe list to hold running threads to update StatusBar.
     */
    private final CopyOnWriteArrayList<Thread> threadList;
    /**
     * Configuration class to read resource file paths from.
     */
    private final IConfiguration configuration;
    //</editor-fold>
    //
    //<editor-fold desc="private static final variables" defaultstate="collapsed">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(MainViewController.class.getName());
    //</editor-fold>
    //
    /**
     * Constructor. Starts the application.
     *
     * @param aStage    Stage
     * @param aMainView MainView
     * @param anAppDir  String path to app dir
     * @param aConfiguration configuration class reading from properties file
     * @throws IllegalArgumentException given application directory is either no directory or does not exist
     * @throws NullPointerException if one param is null
     */
    public MainViewController(Stage aStage, MainView aMainView, String anAppDir, IConfiguration aConfiguration)
            throws IllegalArgumentException, NullPointerException {
        //<editor-fold desc="checks" defaultstate="collapsed">
        Objects.requireNonNull(aStage, "aStage (instance of Stage) is null");
        Objects.requireNonNull(aMainView, "aMainView (instance of MainView) is null");
        Objects.requireNonNull(aMainView, "anAppDir (instance of String) is null");
        Objects.requireNonNull(aConfiguration, "aConfiguration (instance of IConfiguration) is null");
        File tmpAppDirFile = new File(anAppDir);
        if (!tmpAppDirFile.isDirectory() || !tmpAppDirFile.exists()) {
            throw new IllegalArgumentException("The given application directory is either no directory or does not exist");
        }
        //</editor-fold>
        this.configuration = aConfiguration;
        this.moleculeDataModelList = FXCollections.observableArrayList(param -> new Observable[]{param.selectionProperty()});
        this.primaryStage = aStage;
        this.mainView = aMainView;
        this.settingsContainer = new SettingsContainer();
        this.settingsContainer.reloadGlobalSettings();
        this.fragmentationService = new FragmentationService();
        this.fragmentationService.reloadFragmenterSettings();
        this.fragmentationService.reloadActiveFragmenterAndPipeline();
        this.viewToolsManager = new ViewToolsManager(this.configuration, this.settingsContainer);
        this.viewToolsManager.reloadViewToolsSettings();
        this.clusteringService = new ClusteringService();
        this.clusteringService.reloadClusteringSettings();
        this.clusteringRunCounter = 0;
        //<editor-fold desc="show MainView inside primaryStage" defaultstate="collapsed">
        this.mainTabPane = new TabPane();
        this.mainView.getMainCenterPane().getChildren().add(this.mainTabPane);
        GuiUtil.guiBindControlSizeToParentPane(this.mainView.getMainCenterPane(), this.mainTabPane);
        this.scene = new Scene(this.mainView, GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE, GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        String tmpStyleSheetURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.styleFolder")
                        + this.configuration.getProperty("mortar.stylesheet.name")).toExternalForm();
        this.scene.getStylesheets().add(tmpStyleSheetURL);
        this.primaryStage.setTitle(Message.get("Title.text"));
        this.primaryStage.setScene(this.scene);
        this.primaryStage.show();
        this.primaryStage.setMinHeight(GuiDefinitions.GUI_MAIN_VIEW_HEIGHT_VALUE);
        this.primaryStage.setMinWidth(GuiDefinitions.GUI_MAIN_VIEW_WIDTH_VALUE);
        String tmpIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder") + this.configuration.getProperty("mortar.logo.icon.name")).toExternalForm();
        this.primaryStage.getIcons().add(new Image(tmpIconURL));
        //</editor-fold>
        this.isImportRunningProperty = new SimpleBooleanProperty(false);
        this.isExportRunningProperty = new SimpleBooleanProperty(false);
        this.isClusteringRunningProperty = new SimpleBooleanProperty(false);
        this.mapOfFragmentDataModelLists = new HashMap<>(CollectionUtil.calculateInitialHashCollectionCapacity(5));
        this.threadList = new CopyOnWriteArrayList<>();
        this.addListener();
        this.addFragmentationAlgorithmCheckMenuItems();
    }
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Adds listeners and event handlers to control elements etc.
     */
    private void addListener() {
        this.mainView.getMainMenuBar().getExitMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.closeApplication(0)
        );
        this.mainView.getMainMenuBar().getOpenMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.chooseAndImportMoleculeFile(this.primaryStage)
        );
        this.mainView.getMainMenuBar().getCancelImportMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.interruptImport()
        );
        this.mainView.getMainMenuBar().getCancelImportMenuItem().visibleProperty().bind(this.isImportRunningProperty);
        this.mainView.getMainMenuBar().getCancelExportMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.interruptExport()
        );
        this.mainView.getMainMenuBar().getCancelExportMenuItem().visibleProperty().bind(this.isExportRunningProperty);
        //<editor-fold desc="export">
        //fragments export to CSV
        this.mainView.getMainMenuBar().getFragmentsExportToCSVMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.FRAGMENT_CSV_FILE));
        //fragments export to PDB
        this.mainView.getMainMenuBar().getFragmentsExportToPDBMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.FRAGMENT_PDB_FILE));
        //fragments export to PDF
        this.mainView.getMainMenuBar().getFragmentsExportToPDFMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.FRAGMENT_PDF_FILE));
        //fragments export to single SDF
        this.mainView.getMainMenuBar().getFragmentsExportToSingleSDFMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.FRAGMENT_SINGLE_SD_FILE));
        //fragments export to separate SDFs
        this.mainView.getMainMenuBar().getFragmentsExportToSeparateSDFsMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.FRAGMENT_MULTIPLE_SD_FILES));
        //items export to CSV
        this.mainView.getMainMenuBar().getItemsExportToCSVMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.ITEM_CSV_FILE));
        //items export to PDF
        this.mainView.getMainMenuBar().getItemsExportToPDFMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.exportFile(Exporter.ExportTypes.ITEM_PDF_FILE));
        //</editor-fold>
        this.mainView.getMainMenuBar().getFragmentationSettingsMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.openFragmentationSettingsView()
        );
        //<editor-fold desc="clustering">
        this.mainView.getMainMenuBar().getClusteringSettingsMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.openClusteringSettingsView()
        );
        this.mainView.getMainMenuBar().getStartClusteringMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.startClustering()
        );
        this.mainView.getMainMenuBar().getCancelClusteringMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.interruptClustering()
        );
        this.mainView.getMainMenuBar().getCancelClusteringMenuItem().visibleProperty().bind(this.isClusteringRunningProperty);
        //</editor-fold>
        this.mainView.getMainMenuBar().getGlobalSettingsMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.openGlobalSettingsView()
        );
        this.mainView.getMainMenuBar().getPipelineSettingsMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.openPipelineSettingsView()
        );
        this.mainView.getMainMenuBar().getHistogramViewerMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> this.openHistogramView()
        );
        this.mainView.getMainMenuBar().getOverviewViewMenuItem().addEventHandler(
                EventType.ROOT,
                anEvent -> {
                    if (this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.MOLECULES.toString())) {
                        this.openOverviewView(OverviewViewController.DataSources.MOLECULES_TAB);
                    } else if (this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.FRAGMENTS.toString())) {
                        this.openOverviewView(OverviewViewController.DataSources.FRAGMENTS_TAB);
                    } else if (this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.ITEMIZATION.toString())) {
                        //should not happen, since menu item should be disabled if items tab is active
                        throw new IllegalStateException();
                    }
                }
        );
        this.primaryStage.addEventFilter(WindowEvent.WINDOW_CLOSE_REQUEST, (this::closeWindowEvent));
        this.mainView.getMainMenuBar().getAboutViewMenuItem().setOnAction(actionEvent -> new AboutViewController(this.primaryStage, this.configuration));
        this.scene.addEventFilter(KeyEvent.KEY_PRESSED, keyEvent -> {
            Tab tmpSelectedTab = this.mainTabPane.getSelectionModel().getSelectedItem();
            if (tmpSelectedTab == null) {
                keyEvent.consume();
                return;
            }
            //clustering result tabs have no pagination; the key event is left to their table
            if (!(tmpSelectedTab instanceof GridTabForTableView tmpGrid)) {
                return;
            }
            if (GuiDefinitions.KEY_CODE_LAST_PAGE.match(keyEvent) || keyEvent.getCode() == KeyCode.END) {
                tmpGrid.getPagination().setCurrentPageIndex(tmpGrid.getPagination().getPageCount() - 1);
                keyEvent.consume();
            }
            else if (GuiDefinitions.KEY_CODE_FIRST_PAGE.match(keyEvent) || keyEvent.getCode() == KeyCode.HOME) {
                tmpGrid.getPagination().setCurrentPageIndex(0);
                keyEvent.consume();
            }
            else if (keyEvent.getCode() == KeyCode.RIGHT || keyEvent.getCode() == KeyCode.PAGE_UP) {
                tmpGrid.getPagination().setCurrentPageIndex(tmpGrid.getPagination().getCurrentPageIndex() + 1);
                keyEvent.consume();
            }
            else if (keyEvent.getCode() == KeyCode.LEFT || keyEvent.getCode() == KeyCode.PAGE_DOWN) {
                tmpGrid.getPagination().setCurrentPageIndex(tmpGrid.getPagination().getCurrentPageIndex() - 1);
                keyEvent.consume();
            }
        });
        this.mainTabPane.getSelectionModel().selectedItemProperty().addListener((observable, oldValue, newValue) -> Platform.runLater(() -> {
            if (newValue == null) {
                return;
            }
            //a clustering result tab has its own overview button and no histogram
            boolean tmpIsClusteringTab = TabNames.CLUSTERING.toString().equals(newValue.getId());
            this.mainView.getMainMenuBar().getHistogramViewerMenuItem().setDisable(tmpIsClusteringTab
                    || TabNames.MOLECULES.toString().equals(newValue.getId()));
            this.mainView.getMainMenuBar().getOverviewViewMenuItem().setDisable(tmpIsClusteringTab
                    || TabNames.ITEMIZATION.toString().equals(newValue.getId()));
        }));
        this.mainView.getMainCenterPane().setOnDragOver(aDragEvent -> {
            if (aDragEvent.getGestureSource() != this.mainView.getMainCenterPane() && aDragEvent.getDragboard().hasFiles()){
                if (
                  aDragEvent.getDragboard().getFiles().size() > 1 ||
                  !Importer.VALID_IMPORT_FILE_EXTENSIONS_SET.contains(FileUtil.getFileExtension(aDragEvent.getDragboard().getFiles().getFirst().getName()))
                ) {
                    aDragEvent.consume();
                    return;
                }
                aDragEvent.acceptTransferModes(TransferMode.COPY);
            }
            aDragEvent.consume();
        });
        this.mainView.getMainCenterPane().setOnDragDropped(aDragEvent -> {
            Dragboard tmpDragboard = aDragEvent.getDragboard();
            boolean tmpSucceeded = false;
            if (tmpDragboard.hasFiles()) {
                this.importMoleculeFile(tmpDragboard.getFiles().getFirst());
                tmpSucceeded = true;
            }
            aDragEvent.setDropCompleted(tmpSucceeded);
            aDragEvent.consume();
        });
    }
    //
    /**
     * Closes application.
     *
     * @param aStatus the status to use for calling System.exit(); a nonzero status code indicates abnormal termination
     */
    private void closeApplication(int aStatus) {
        if (!moleculeDataModelList.isEmpty() && (!this.isFragmentationStopAndDataLossConfirmed())) {
            return;
        }
        this.settingsContainer.preserveSettings();
        this.viewToolsManager.persistViewToolsSettings();
        this.fragmentationService.persistFragmenterSettings();
        this.fragmentationService.persistSelectedFragmenterAndPipeline();
        this.clusteringService.persistClusteringSettings();
        if (this.isFragmentationRunning) {
            this.interruptFragmentation();
        }
        if (this.clusteringTask != null && this.clusteringTask.isRunning()) {
            this.interruptClustering();
        }
        MainViewController.LOGGER.info(BasicDefinitions.MORTAR_SESSION_END);
        Platform.exit();
        System.exit(aStatus);
    }
    //
    /**
     * Opens a dialog to warn the user of possible data loss and stopping a running fragmentation, e.g. when a new
     * molecule set should be imported or the application shut down. Returns true if "OK" was clicked, "false" for cancel
     * button.
     *
     * @return true if "OK" was clicked, false for "Cancel"
     */
    private boolean isFragmentationStopAndDataLossConfirmed() {
        ButtonType tmpConfirmationResult;
        if (this.isFragmentationRunning) {
            tmpConfirmationResult = GuiUtil.guiConfirmationAlert(
                    Message.get("MainViewController.Warning.FragmentationRunning.Title"),
                    Message.get("MainViewController.Warning.FragmentationRunning.Header"),
                    Message.get("MainViewController.Warning.FragmentationRunning.Content"));
        } else {
            tmpConfirmationResult = GuiUtil.guiConfirmationAlert(
                    Message.get("MainViewController.Warning.DataLoss.Title"),
                    Message.get("MainViewController.Warning.DataLoss.Header"),
                    Message.get("MainViewController.Warning.DataLoss.Content"));
        }
        return tmpConfirmationResult == ButtonType.OK;
    }
    //
    /**
     * Closes the application via closeApplication method when close window event was fired.
     *
     * @param anEvent WindowEvent
     */
    private void closeWindowEvent(WindowEvent anEvent) {
        this.closeApplication(0);
        anEvent.consume();
    }
    //

    /**
     * Opens a file choose, loads the chosen file and opens molecules tab
     *
     * @param aParentStage Stage where to open the file chooser dialog
     */
    private void chooseAndImportMoleculeFile(Stage aParentStage) {
        Importer tmpImporter = new Importer(this.settingsContainer);
        File tmpFile = tmpImporter.openFile(aParentStage);
        if (tmpFile == null)
            return;
        this.importMoleculeFile(tmpFile, tmpImporter);
    }
    //
    /**
     * Loads molecule file and opens molecules tab.
     * Convenient method to avoid using a null parameter for the importer.
     *
     * @param aFile File that contains molecular data
     */
    private void importMoleculeFile(File aFile) {
        this.importMoleculeFile(aFile, new Importer(this.settingsContainer));
    }
    //
    /**
     * Loads molecule file and opens molecules tab.
     *
     * @param aFile File that contains molecular data
     */
    private void importMoleculeFile(File aFile, Importer anImporter) {
        if (Objects.isNull(aFile)) {
            return;
        }
        if (!this.moleculeDataModelList.isEmpty()) {
            if (!this.isFragmentationStopAndDataLossConfirmed()) {
                return;
            }
            this.fragmentationService.clearCache();
        }
        Importer tmpImporter = Objects.requireNonNullElseGet(anImporter, () -> new Importer(this.settingsContainer));
        if (this.isFragmentationRunning) {
            this.interruptFragmentation();
        }
        if (this.isImportRunningProperty.get()) {
            this.interruptImport();
        }
        if (this.isExportRunningProperty.get()) {
            this.interruptExport();
        }
        this.clearGuiAndCollections();
        boolean tmpIsRegardStereo = this.settingsContainer.getRegardStereochemistrySetting();
        boolean tmpIsFillOpenValences = this.settingsContainer.getAddImplicitHydrogensAtImportSetting();
        boolean tmpIsKekulizationEnforced = this.settingsContainer.getImportAromaticsAsKekuleStructuresSetting();
        this.importTask = new Task<>() {
            @Override
            protected List<MoleculeDataModel> call() throws Exception {
                List<MoleculeDataModel> tmpSet = tmpImporter.importMoleculeFile(aFile, tmpIsRegardStereo, tmpIsFillOpenValences, tmpIsKekulizationEnforced);
                return tmpSet;
            }
        };
        this.importTask.setOnSucceeded(event ->
            //note: setOnSucceeded() takes place in the JavaFX GUI thread again but still runLater() is necessary to wait
            // for the thread to be free for the update
            Platform.runLater(() -> {
                List<MoleculeDataModel> tmpImportedMoleculeDataModels = null;
                try {
                    tmpImportedMoleculeDataModels = this.importTask.get();
                } catch (InterruptedException | ExecutionException anException) {
                    MainViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
                    GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                            Message.get("Importer.FileImportExceptionAlert.Header"),
                            Message.get("Importer.FileImportExceptionAlert.Text"),
                            anException);
                    this.updateStatusBar(this.importerThread, Message.get("Status.importFailed"));
                }
                this.moleculeDataModelList.addAll(tmpImportedMoleculeDataModels);
                if (tmpImportedMoleculeDataModels == null || tmpImportedMoleculeDataModels.isEmpty() || this.moleculeDataModelList.isEmpty()) {
                    MainViewController.LOGGER.log(Level.WARNING, "Import failed, set of imported molecules is null or empty");
                    this.updateStatusBar(this.importerThread, Message.get("Status.importFailed"));
                    this.isImportRunningProperty.setValue(false);
                    Platform.runLater(() -> {
                        GuiUtil.guiMessageAlert(Alert.AlertType.WARNING,
                                Message.get("Error.ExceptionAlert.Title"),
                                Message.get("Importer.FileImportEmptyAlert.Header"),
                                Message.get("Importer.FileImportEmptyAlert.Content"));
                    });
                    return;
                }
                this.mainView.getMainMenuBar().getExportMenu().setDisable(true);
                this.mainView.getMainMenuBar().getHistogramViewerMenuItem().setDisable(true);
                this.mainView.getMainMenuBar().getOverviewViewMenuItem().setDisable(false);
                this.primaryStage.setTitle(Message.get("Title.text") + " - " + tmpImporter.getFileName() + " - " + tmpImportedMoleculeDataModels.size() +
                        " " + Message.get((tmpImportedMoleculeDataModels.size() == 1 ? "Title.molecule" : "Title.molecules")));
                this.importedFileName = tmpImporter.getFileName();
                this.updateStatusBar(this.importerThread, Message.get("Status.imported"));
                this.isImportRunningProperty.setValue(false);
                this.mainView.getMainCenterPane().setStyle("-fx-background-image: none");
                this.openMoleculesTab();
            })
        );
        this.importTask.setOnCancelled(event -> {
            this.updateStatusBar(this.importerThread, Message.get("Status.canceled"));
            this.isImportRunningProperty.setValue(false);
        });
        this.importTask.setOnFailed(event -> {
            Exception tmpCause = (Exception) event.getSource().getException();
            MainViewController.LOGGER.log(Level.SEVERE, tmpCause.toString(), tmpCause);
            this.updateStatusBar(this.importerThread, Message.get("Status.importFailed"));
            this.isImportRunningProperty.setValue(false);
            Platform.runLater(() -> {
                GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                        Message.get("Importer.FileImportExceptionAlert.Header"),
                        Message.get("Importer.FileImportExceptionAlert.Text"),
                        tmpCause);
            });
        });
        this.importerThread = new Thread(importTask);
        this.importerThread.setName(ThreadType.IMPORT_THREAD.getThreadName());
        this.importerThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        this.importerThread.setDaemon(false);
        this.importerThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number
        this.isImportRunningProperty.setValue(true);
        this.updateStatusBar(this.importerThread, Message.get("Status.importing"));
        this.importerThread.start();
    }
    //
    /**
     * Exports the given type of file.
     *
     * @param anExportType Enum to specify what type of file to export
     */
    private void exportFile(Exporter.ExportTypes anExportType) {
        //a clustering result tab has no fragments or items, its data is exported via the button in the tab
        if (!(this.mainTabPane.getSelectionModel().getSelectedItem() instanceof GridTabForTableView)) {
            GuiUtil.guiMessageAlert(
                    Alert.AlertType.INFORMATION,
                    Message.get("Exporter.MessageAlert.NoDataAvailable.title"),
                    Message.get("Exporter.MessageAlert.NoDataAvailable.header"),
                    null
            );
            return;
        }
        if ((this.mainTabPane.getSelectionModel().getSelectedItem()).getId().equals(TabNames.MOLECULES.toString())) {
            GuiUtil.guiConfirmationAlert(Message.get("Exporter.confirmationAlert.moleculesTabSelected.title"),
                    Message.get("Exporter.confirmationAlert.moleculesTabSelected.header"),
                    Message.get("Exporter.confirmationAlert.moleculesTabSelected.text"));
            return;
        }
        switch (anExportType) {
            case Exporter.ExportTypes.FRAGMENT_CSV_FILE, Exporter.ExportTypes.FRAGMENT_PDB_FILE, Exporter.ExportTypes.FRAGMENT_PDF_FILE, Exporter.ExportTypes.FRAGMENT_SINGLE_SD_FILE, FRAGMENT_MULTIPLE_SD_FILES:
                if (this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS) == null ||
                        this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS).isEmpty() ||
                        ((GridTabForTableView) mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle() == null) {
                    GuiUtil.guiMessageAlert(
                            Alert.AlertType.INFORMATION,
                            Message.get("Exporter.MessageAlert.NoDataAvailable.title"),
                            Message.get("Exporter.MessageAlert.NoDataAvailable.header"),
                            null
                    );
                    return;
                }
                break;
            case Exporter.ExportTypes.ITEM_CSV_FILE, Exporter.ExportTypes.ITEM_PDF_FILE:
                if (this.getItemsListOfSelectedFragmentationByTabId(TabNames.ITEMIZATION) == null ||
                        this.getItemsListOfSelectedFragmentationByTabId(TabNames.ITEMIZATION).isEmpty() ||
                        this.moleculeDataModelList == null || this.moleculeDataModelList.isEmpty() ||
                        ((GridTabForTableView) mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle() == null) {
                    GuiUtil.guiMessageAlert(
                            Alert.AlertType.INFORMATION,
                            Message.get("Exporter.MessageAlert.NoDataAvailable.title"),
                            Message.get("Exporter.MessageAlert.NoDataAvailable.header"),
                            null
                    );
                    return;
                }
                break;
        }
        Exporter tmpExporter = new Exporter(this.settingsContainer);
        if (this.isExportRunningProperty.get()) {
            this.interruptExport();
        }
        //returns null if file chooser dialog was cancelled
        File tmpExportFile = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage, anExportType,
                ((GridTabForTableView) this.mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle());
        if (tmpExportFile == null) {
            return;
        }
        boolean tmpGenerate2dAtomCoordinates = false;
        if ((anExportType.equals(Exporter.ExportTypes.FRAGMENT_PDB_FILE)
                || anExportType.equals(Exporter.ExportTypes.FRAGMENT_SINGLE_SD_FILE)
                || anExportType.equals(Exporter.ExportTypes.FRAGMENT_MULTIPLE_SD_FILES))
                && (!ChemUtil.checkMoleculeListForCoordinates(this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS)))) {
            ButtonType tmpConfirmationResult = GuiUtil.guiYesNoCancelConfirmationAlert(
                    Message.get("Exporter.FragmentsTab.ConfirmationAlert.No3dInformationAvailable.title"),
                    Message.get("Exporter.FragmentsTab.ConfirmationAlert.No3dInformationAvailable.header"),
                    Message.get("Exporter.FragmentsTab.ConfirmationAlert.No3dInformationAvailable.text")
            );
            /*
            yes -> generate coordinates
            no -> export but do not generate coordinates
            cancel -> abort export
             */
            if (tmpConfirmationResult == ButtonType.CANCEL) {
                return;
            }
            tmpGenerate2dAtomCoordinates = tmpConfirmationResult == ButtonType.YES;
        }
        //reassigned because variable needs to be effectively final to be used in the inner classes below
        boolean tmpGenerate2dAtomCoordinatesFinal = tmpGenerate2dAtomCoordinates;
        this.exportTask = new Task<>() {
            @Override
            protected List<String> call() throws Exception {
                return switch (anExportType) {
                    case Exporter.ExportTypes.FRAGMENT_CSV_FILE -> tmpExporter.exportCsvFile(
                            tmpExportFile,
                            MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                            ((GridTabForTableView) MainViewController.this.mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle(),
                            MainViewController.this.settingsContainer.getCsvExportSeparatorSettingCharacter(),
                            TabNames.FRAGMENTS
                    );
                    case Exporter.ExportTypes.FRAGMENT_PDB_FILE ->
                            tmpExporter.exportFragmentsAsChemicalFile(
                                    tmpExportFile,
                                    MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                                    ChemFileTypes.PDB,
                                    tmpGenerate2dAtomCoordinatesFinal
                            );
                    case Exporter.ExportTypes.FRAGMENT_PDF_FILE -> tmpExporter.exportPdfFile(
                            tmpExportFile,
                            MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                            MainViewController.this.moleculeDataModelList,
                            ((GridTabForTableView) MainViewController.this.mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle(),
                            MainViewController.this.importedFileName,
                            TabNames.FRAGMENTS
                    );
                    case Exporter.ExportTypes.FRAGMENT_SINGLE_SD_FILE ->
                            tmpExporter.exportFragmentsAsChemicalFile(
                                    tmpExportFile,
                                    MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                                    ChemFileTypes.SDF,
                                    tmpGenerate2dAtomCoordinatesFinal,
                                    true
                            );
                    case Exporter.ExportTypes.FRAGMENT_MULTIPLE_SD_FILES ->
                            tmpExporter.exportFragmentsAsChemicalFile(
                                    tmpExportFile,
                                    MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                                    ChemFileTypes.SDF,
                                    tmpGenerate2dAtomCoordinatesFinal,
                                    false
                            );
                    case Exporter.ExportTypes.ITEM_CSV_FILE -> tmpExporter.exportCsvFile(
                            tmpExportFile,
                            MainViewController.this.moleculeDataModelList,
                            ((GridTabForTableView) MainViewController.this.mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle(),
                            MainViewController.this.settingsContainer.getCsvExportSeparatorSettingCharacter(),
                            TabNames.ITEMIZATION
                    );
                    case Exporter.ExportTypes.ITEM_PDF_FILE -> tmpExporter.exportPdfFile(
                            tmpExportFile,
                            MainViewController.this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS),
                            MainViewController.this.moleculeDataModelList,
                            ((GridTabForTableView) MainViewController.this.mainTabPane.getSelectionModel().getSelectedItem()).getFragmentationNameOutOfTitle(),
                            MainViewController.this.importedFileName,
                            TabNames.ITEMIZATION
                    );
                    default -> throw new UnsupportedOperationException("Unknown export type.");
                };
            }
        };
        this.exportTask.setOnSucceeded(event -> {
            this.isExportRunningProperty.setValue(false);
            this.updateStatusBar(this.exporterThread, Message.get("Status.finished"));
            List<String> tmpFailedExportFragments = this.exportTask.getValue();
            if (tmpFailedExportFragments == null) {
                GuiUtil.guiMessageAlert(Alert.AlertType.WARNING,
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                        null);
                return;
            }
            if (!tmpFailedExportFragments.isEmpty()) {
                StringBuilder tmpStringBuilder = new StringBuilder();
                for (String tmpFragmentName : tmpFailedExportFragments) {
                    tmpStringBuilder.append(tmpFragmentName).append("\n");
                }
                GuiUtil.guiExpandableAlert(
                        Alert.AlertType.WARNING.toString(),
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                        tmpFailedExportFragments.size() + " " + Message.get("Exporter.FragmentsTab.ExportNotPossible.label"),
                        tmpStringBuilder.toString()
                );
            }
        });
        this.exportTask.setOnCancelled(event -> {
            this.isExportRunningProperty.setValue(false);
            MainViewController.LOGGER.log(Level.SEVERE, "Export canceled");
            this.updateStatusBar(this.exporterThread, Message.get("Status.canceled"));
        });
        this.exportTask.setOnFailed(event -> {
            this.isExportRunningProperty.setValue(false);
            MainViewController.LOGGER.log(Level.WARNING, event.getSource().getException().toString(), event.getSource().getException());
            this.updateStatusBar(this.exporterThread, Message.get("Status.failed"));
            GuiUtil.guiMessageAlert(
                    Alert.AlertType.WARNING,
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                    null);

        });
        this.exporterThread = new Thread(this.exportTask);
        this.exporterThread.setName(ThreadType.EXPORT_THREAD.getThreadName());
        this.exporterThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        this.exporterThread.setDaemon(false);
        this.exporterThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number
        this.isExportRunningProperty.setValue(true);
        this.updateStatusBar(this.exporterThread, Message.get("Status.exporting"));
        this.exporterThread.start();
    }
    //
    /**
     * Opens settings view for fragmentation settings.
     */
    private void openFragmentationSettingsView() {
        new FragmentationSettingsViewController(this.primaryStage,
                this.fragmentationService.getFragmenters(),
                this.fragmentationService.getSelectedFragmenter().getFragmentationAlgorithmDisplayName(),
                this.configuration);
    }
    //
    /**
     * Opens PipelineSettingsView.
     */
    private void openPipelineSettingsView() {
        PipelineSettingsViewController tmpPipelineSettingsViewController =
                new PipelineSettingsViewController(this.primaryStage, this.fragmentationService, !this.moleculeDataModelList.isEmpty(), this.isFragmentationRunning, this.configuration);
        if (tmpPipelineSettingsViewController.isFragmentationStarted()) {
            this.startFragmentation(tmpPipelineSettingsViewController.isFragmentationStarted());
        }
    }
    //
    /**
     * Opens HistogramView.
     */
    private void openHistogramView()  {
        List<MoleculeDataModel> tmpMoleculesList = this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS);
        List<FragmentDataModel> tmpFragmentsList = new ArrayList<>(tmpMoleculesList.size());
        for (MoleculeDataModel tmpMolecule : tmpMoleculesList) {
            tmpFragmentsList.add((FragmentDataModel) tmpMolecule);
        }
        this.viewToolsManager.openHistogramView(this.primaryStage, tmpFragmentsList);
    }
    //
    /**
     * Adds CheckMenuItems for fragmentation algorithms to MainMenuBar.
     */
    private void addFragmentationAlgorithmCheckMenuItems() {
        ToggleGroup tmpToggleGroup = new ToggleGroup();
        for (IMoleculeFragmenter tmpFragmenter : this.fragmentationService.getFragmenters()) {
            RadioMenuItem tmpRadioMenuItem = new RadioMenuItem(tmpFragmenter.getFragmentationAlgorithmDisplayName());
            tmpRadioMenuItem.setToggleGroup(tmpToggleGroup);
            this.mainView.getMainMenuBar().getFragmentationAlgorithmMenu().getItems().add(tmpRadioMenuItem);
            if (!Objects.isNull(this.fragmentationService.getSelectedFragmenter())
                    && tmpFragmenter.getFragmentationAlgorithmDisplayName()
                        .equals(this.fragmentationService.getSelectedFragmenter().getFragmentationAlgorithmDisplayName())) {
                tmpToggleGroup.selectToggle(tmpRadioMenuItem);
            }
        }
        tmpToggleGroup.selectedToggleProperty().addListener((observableValue, oldValue, newValue) -> {
            if (tmpToggleGroup.getSelectedToggle() != null) {
                this.fragmentationService.setSelectedFragmenter(((RadioMenuItem) newValue).getText());
                this.fragmentationService.setSelectedFragmenterDisplayName(((RadioMenuItem) newValue).getText());
            }
        });
    }
    //
    /**
     * Opens settings view for global settings.
     */
    private void openGlobalSettingsView() {
        SettingsViewController tmpSettingsViewController = new SettingsViewController(this.primaryStage, this.settingsContainer, this.configuration);
        Platform.runLater(() -> {
            if (tmpSettingsViewController.hasRowsPerPageChanged()) {
                for (Tab tmpTab : this.mainTabPane.getTabs()) {
                    //clustering result tabs have no pagination
                    if (!(tmpTab instanceof GridTabForTableView)) {
                        continue;
                    }
                    // type of generic not given because it does not matter here, only the size of the items list
                    TableView<?> tmpTableView = ((GridTabForTableView) tmpTab).getTableView();
                    int tmpListSize = ((IDataTableView) tmpTableView).getItemsList().size();
                    int tmpPageIndex = ((GridTabForTableView) tmpTab).getPagination().getCurrentPageIndex();
                    int tmpRowsPerPage = this.settingsContainer.getRowsPerPageSetting();
                    int tmpPageCount = tmpListSize / tmpRowsPerPage;
                    if (tmpListSize % tmpRowsPerPage > 0) {
                        tmpPageCount++;
                    }
                    if (tmpPageIndex > tmpPageCount) {
                        tmpPageIndex = tmpPageCount;
                    }
                    /*
                    the following might cause "javafx.scene.control.skin.VirtualFlow addTrailingCells
                    INFO: index exceeds maxCellCount. Check size calculations for class javafx.scene.control.TableRow"
                    when the new rows per page value is smaller than the older one, but it is not a real problem;
                    the refreshed GUI just needs to "scroll" to a different position
                    */
                    ((GridTabForTableView) tmpTab).getPagination().setPageCount(tmpPageCount);
                    ((GridTabForTableView) tmpTab).getPagination().setCurrentPageIndex(tmpPageIndex);
                    ((GridTabForTableView) tmpTab).getTableView().refresh();
                    GuiUtil.setImageStructureHeight(((GridTabForTableView) tmpTab).getTableView(), ((GridTabForTableView) tmpTab).getTableView().getHeight(), this.settingsContainer.getRowsPerPageSetting());
                    ((GridTabForTableView) tmpTab).getTableView().refresh();
                }
            }
            if (tmpSettingsViewController.hasKeepAtomContainerInDataModelChanged()) {
                for (MoleculeDataModel tmpMoleculeDataModel : this.moleculeDataModelList) {
                    tmpMoleculeDataModel.setKeepAtomContainer(this.settingsContainer.getKeepAtomContainerInDataModelSetting());
                }
                for (ObservableList<FragmentDataModel> tmpFragmentDataModelList : this.mapOfFragmentDataModelLists.values()) {
                    for (FragmentDataModel tmpFragmentDataModel : tmpFragmentDataModelList) {
                        tmpFragmentDataModel.setKeepAtomContainer(this.settingsContainer.getKeepAtomContainerInDataModelSetting());
                    }
                }
            }
        });
    }
    //
    /**
     * Opens OverviewView.
     *
     * @param aDataSource Source of the data to be shown in the overview view
     */
    private void openOverviewView(OverviewViewController.DataSources aDataSource) {
        try {
            switch (aDataSource) {
                case OverviewViewController.DataSources.MOLECULES_TAB -> {
                    if (!(this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.MOLECULES.toString())))
                        //should not happen
                        throw new IllegalStateException();
                    this.viewToolsManager.openOverviewView(
                            this.primaryStage,
                            aDataSource,
                            ((GridTabForTableView) mainTabPane.getSelectionModel().getSelectedItem()).getTitle(),
                            getItemsListOfSelectedFragmentationByTabId(TabNames.MOLECULES)
                    );
                }
                case OverviewViewController.DataSources.FRAGMENTS_TAB -> {
                    if (!(this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.FRAGMENTS.toString())))
                        //should not happen
                        throw new IllegalStateException();
                    this.viewToolsManager.openOverviewView(
                            this.primaryStage,
                            aDataSource,
                            ((GridTabForTableView) mainTabPane.getSelectionModel().getSelectedItem()).getTitle(),
                            this.getItemsListOfSelectedFragmentationByTabId(TabNames.FRAGMENTS)
                    );
                }
                case OverviewViewController.DataSources.PARENT_MOLECULES_SAMPLE -> {
                    if (!(this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.FRAGMENTS.toString())))
                        //should not happen
                        throw new IllegalStateException();
                    //Parent-Molecules of the Fragments-Tab (showing all fragments of one molecule in the overview view)
                    GridTabForTableView tmpSelectedTab = (GridTabForTableView) this.mainTabPane.getSelectionModel().getSelectedItem();
                    //IllegalStateException if there is more than one cell selected
                    if (((TableView<?>) tmpSelectedTab.getTableView()).getSelectionModel().getSelectedCells().size() > 1) {
                        //should not happen
                        throw new IllegalStateException();
                    }
                    //getting the data for the overview view
                    List<MoleculeDataModel> tmpDataForOverviewView = new ArrayList<>();
                    int tmpSelectedRowIndex = ((TableView<?>) tmpSelectedTab.getTableView()).getSelectionModel().getSelectedCells().getFirst().getRow();
                    int tmpIndexInDataList = tmpSelectedTab.getPagination().getCurrentPageIndex() * this.settingsContainer.getRowsPerPageSetting() + tmpSelectedRowIndex;
                    //adding the fragment itself
                    tmpDataForOverviewView.add(((IDataTableView) tmpSelectedTab.getTableView()).getItemsList().get(tmpIndexInDataList));
                    //adding the sample of parent molecules
                    tmpDataForOverviewView.addAll(((FragmentDataModel) ((IDataTableView) tmpSelectedTab.getTableView()).getItemsList().get(tmpIndexInDataList)).getParentMolecules());
                    this.viewToolsManager.openOverviewView(
                            this.primaryStage,
                            OverviewViewController.DataSources.PARENT_MOLECULES_SAMPLE,
                            null,
                            tmpDataForOverviewView
                    );
                }
                case OverviewViewController.DataSources.ITEM_WITH_FRAGMENTS_SAMPLE -> {
                    if (!(this.mainTabPane.getSelectionModel().getSelectedItem().getId().equals(TabNames.ITEMIZATION.toString())))
                        //should not happen
                        throw new IllegalStateException();
                    //Items-Tab (showing all fragments of one molecule in the overview view)
                    GridTabForTableView tmpSelectedTab = (GridTabForTableView) this.mainTabPane.getSelectionModel().getSelectedItem();
                    //IllegalStateException if there is more than one cell selected
                    if (((TableView<?>) tmpSelectedTab.getTableView()).getSelectionModel().getSelectedCells().size() > 1) {
                        //should not happen
                        throw new IllegalStateException();
                    }
                    //getting the data for the overview view
                    List<MoleculeDataModel> tmpDataForOverviewView = new ArrayList<>();
                    int tmpSelectedRowIndex = ((TableView<?>) tmpSelectedTab.getTableView()).getSelectionModel().getSelectedCells().getFirst().getRow();
                    int tmpIndexInDataList = tmpSelectedTab.getPagination().getCurrentPageIndex() * this.settingsContainer.getRowsPerPageSetting() + tmpSelectedRowIndex;
                    //adding the item itself
                    tmpDataForOverviewView.add(((IDataTableView) tmpSelectedTab.getTableView()).getItemsList().get(tmpIndexInDataList));
                    //adding the sample of fragments
                    tmpDataForOverviewView.addAll(((IDataTableView) tmpSelectedTab.getTableView()).getItemsList().get(tmpIndexInDataList).getFragmentsOfSpecificFragmentation(tmpSelectedTab.getFragmentationNameOutOfTitle()));
                    this.viewToolsManager.openOverviewView(
                            this.primaryStage,
                            OverviewViewController.DataSources.ITEM_WITH_FRAGMENTS_SAMPLE,
                            null,
                            tmpDataForOverviewView
                    );
                }
                default -> {
                    //should not happen
                    throw new IllegalStateException();
                }
            }
        } catch (IllegalStateException anIllegalStateException) {
            MainViewController.LOGGER.log(Level.SEVERE, anIllegalStateException.toString(), anIllegalStateException);
            return;
        }
        //
        int tmpIndexOfMoleculeDataModelToReturnTo = this.viewToolsManager.getCachedIndexOfStructureInMoleculeDataModelList();
        //since -1 is returned, if no specific structure should be shown
        if (tmpIndexOfMoleculeDataModelToReturnTo >= 0) {
            //go to page showing the structure of the MoleculeDataModel with the given index
            int tmpNewPageIndex = tmpIndexOfMoleculeDataModelToReturnTo / this.settingsContainer.getRowsPerPageSetting();
            ((GridTabForTableView) this.mainTabPane.getSelectionModel().getSelectedItem()).getPagination()
                    .setCurrentPageIndex(tmpNewPageIndex);
            // unnecessary to provide generic type
            TableView tmpSelectedTabTableView = ((GridTabForTableView) this.mainTabPane.getSelectionModel()
                    .getSelectedItem()).getTableView();
            if (tmpSelectedTabTableView.getClass() == MoleculesDataTableView.class) {
                //select structure cell
                int tmpRowIndexOfStructure = tmpIndexOfMoleculeDataModelToReturnTo
                        % this.settingsContainer.getRowsPerPageSetting();
                tmpSelectedTabTableView.getSelectionModel().clearSelection();
                tmpSelectedTabTableView.getSelectionModel().select(tmpRowIndexOfStructure,
                        ((MoleculesDataTableView) tmpSelectedTabTableView).getStructureColumn());
            } else if (tmpSelectedTabTableView.getClass() == FragmentsDataTableView.class) {
                //select structure cell
                int tmpRowIndexOfStructure = tmpIndexOfMoleculeDataModelToReturnTo
                        % this.settingsContainer.getRowsPerPageSetting();
                tmpSelectedTabTableView.getSelectionModel().clearSelection();
                tmpSelectedTabTableView.getSelectionModel().select(tmpRowIndexOfStructure,
                        ((FragmentsDataTableView) tmpSelectedTabTableView).getStructureColumn());
            }
        }
        this.viewToolsManager.resetCachedIndexOfStructureInMoleculeDataModelList();
    }
    //
    /**
     * Opens molecules tab.
     */
    private void openMoleculesTab() {
        this.moleculesDataTableView = new MoleculesDataTableView(this.configuration);
        this.moleculesDataTableView.setItemsList(this.moleculeDataModelList);
        GridTabForTableView tmpMoleculesTab = new GridTabForTableView(Message.get("MainTabPane.moleculesTab.title"), TabNames.MOLECULES.name(), this.moleculesDataTableView);
        this.mainTabPane.getTabs().add(tmpMoleculesTab);
        Pagination tmpPagination = this.createPaginationWithSuitablePageCount(this.moleculeDataModelList.size());
        tmpPagination.setPageFactory(pageIndex -> this.moleculesDataTableView.createMoleculeTableViewPage(pageIndex, this.settingsContainer));
        tmpMoleculesTab.addPaginationToGridPane(tmpPagination);
        HBox tmpFragmentationButtonsHBox = new HBox();
        tmpFragmentationButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpFragmentationButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpFragmentationButtonsHBox.setAlignment(Pos.CENTER_LEFT);
        this.fragmentationButton = new Button();
        this.fragmentationButton.textProperty().bind(this.fragmentationService.selectedFragmenterDisplayNameProperty());
        Tooltip tmpTooltip = GuiUtil.createTooltip("");
        tmpTooltip.textProperty().bind(Bindings.format(Message.get("MainTabPane.moleculesTab.fragmentButton.text"), this.fragmentationService.selectedFragmenterDisplayNameProperty()));
        this.fragmentationButton.setTooltip(tmpTooltip);
        double tmpTextWidth = new Text(this.fragmentationService.getSelectedFragmenterDisplayName()).getLayoutBounds().getWidth() + 20;
        this.fragmentationButton.setPrefWidth(tmpTextWidth);
        this.fragmentationButton.setMinWidth(tmpTextWidth);
        this.fragmentationButton.setMaxWidth(tmpTextWidth);
        this.fragmentationButton.setPrefHeight(GuiDefinitions.GUI_BUTTON_HEIGHT_VALUE);
        this.fragmentationService.selectedFragmenterDisplayNameProperty().addListener((observable, oldValue, newValue) -> {
            double tmpTextWidthChange = new Text(newValue).getLayoutBounds().getWidth() + 20;
            this.fragmentationButton.setPrefWidth(tmpTextWidthChange);
            this.fragmentationButton.setMinWidth(tmpTextWidthChange);
            this.fragmentationButton.setMaxWidth(tmpTextWidthChange);
        });
        tmpFragmentationButtonsHBox.getChildren().add(this.fragmentationButton);
        this.cancelFragmentationButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.moleculesTab.cancelFragmentationButton.text"));
        this.cancelFragmentationButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.moleculesTab.cancelFragmentationButton.tooltip")));
        this.cancelFragmentationButton.setVisible(false);
        tmpFragmentationButtonsHBox.getChildren().add(this.cancelFragmentationButton);
        tmpMoleculesTab.addNodeToGridPane(tmpFragmentationButtonsHBox, 0, 1, 1, 1);
        this.fragmentationButton.setOnAction(event -> this.startFragmentation());
        this.cancelFragmentationButton.setOnAction(event -> this.interruptFragmentation());
        HBox tmpViewButtonsHBox = new HBox();
        tmpViewButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpViewButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpViewButtonsHBox.setAlignment(Pos.CENTER_RIGHT);
        tmpViewButtonsHBox.setMaxWidth(GuiDefinitions.GUI_GRIDPANE_FOR_NODE_ALIGNMENT_THIRD_COL_WIDTH);
        Button tmpOpenOverviewViewButton = GuiUtil.getButtonOfStandardSize(Message.get("MainView.showOverviewViewButton.text"));
        tmpOpenOverviewViewButton.setTooltip(GuiUtil.createTooltip(Message.get("MainView.showOverviewViewButton.tooltip")));
        tmpViewButtonsHBox.getChildren().add(tmpOpenOverviewViewButton);
        tmpMoleculesTab.addNodeToGridPane(tmpViewButtonsHBox, 2, 1, 1, 1);
        tmpOpenOverviewViewButton.setOnAction(event -> this.openOverviewView(OverviewViewController.DataSources.MOLECULES_TAB));
        this.moleculesDataTableView.addTableViewHeightListener(this.settingsContainer);
        this.moleculesDataTableView.getCopyMenuItem().setOnAction(event -> GuiUtil.copySelectedTableViewCellsToClipboard(this.moleculesDataTableView));
        //similarity search for the selected molecule, available once there is a clustering result
        this.moleculesDataTableView.getContextMenu().setOnShowing(event -> {
            boolean tmpIsDisabled = this.getActiveClusteringRunResult() == null
                    || this.moleculesDataTableView.getSelectionModel().getSelectedItem() == null;
            this.moleculesDataTableView.getSimilaritySearchMenuItem().setDisable(tmpIsDisabled);
            this.moleculesDataTableView.getDescriptorValuesMenuItem().setDisable(tmpIsDisabled);
        });
        this.moleculesDataTableView.getDescriptorValuesMenuItem().setOnAction(event -> {
            Object tmpSelectedItem = this.moleculesDataTableView.getSelectionModel().getSelectedItem();
            if (tmpSelectedItem instanceof MoleculeDataModel tmpMolecule) {
                this.showDescriptorValuesOfMolecule(tmpMolecule, this.getActiveClusteringRunResult());
            }
        });
        this.moleculesDataTableView.getSimilaritySearchMenuItem().setOnAction(event -> {
            Object tmpSelectedItem = this.moleculesDataTableView.getSelectionModel().getSelectedItem();
            if (tmpSelectedItem instanceof MoleculeDataModel tmpMolecule) {
                this.searchSimilarMolecules(tmpMolecule, this.getActiveClusteringRunResult());
            }
        });
        this.moleculesDataTableView.setOnKeyPressed(event -> {
            if (GuiDefinitions.KEY_CODE_COPY.match(event)) {
                GuiUtil.copySelectedTableViewCellsToClipboard(this.moleculesDataTableView);
            }
        });
        int tmpRowsPerPage = this.settingsContainer.getRowsPerPageSetting();
        this.moleculesDataTableView.setOnSort((EventHandler<SortEvent<TableView>>) event ->
                GuiUtil.sortTableViewGlobally(event, tmpPagination, tmpRowsPerPage));
        this.moleculesDataTableView.widthProperty().addListener((observable, oldValue, newValue) -> {
            for(Object tmpObject : this.moleculesDataTableView.getItems()) {
                ((MoleculeDataModel) tmpObject).setStructureImageWidth(this.moleculesDataTableView.getStructureColumn().getWidth());
            }
        });
    }
    //
    /**
     * Creates a new JavaFx pagination control that is configured with a suitable page count for the given number
     * of molecules/fragments taking into account the rows per page setting. Also sets the MORTAR custom pagination skin
     * as skin of the new pagination instance and configures its growth behavior. The page factory is *NOT* set.
     *
     * @param aListSize number of molecules/fragments to display
     * @return configured pagination control instance
     */
    private Pagination createPaginationWithSuitablePageCount(int aListSize) {
        int tmpRowsPerPage = this.settingsContainer.getRowsPerPageSetting();
        int tmpPageCount = aListSize / tmpRowsPerPage;
        if (aListSize % tmpRowsPerPage > 0) {
            tmpPageCount++;
        }
        if (aListSize == 0) {
            tmpPageCount = 1;
        }
        Pagination tmpPagination = new Pagination(tmpPageCount, 0);
        tmpPagination.setSkin(new CustomPaginationSkin(tmpPagination));
        VBox.setVgrow(tmpPagination, Priority.ALWAYS);
        HBox.setHgrow(tmpPagination, Priority.ALWAYS);
        return tmpPagination;
    }
    //
    /**
     * Cancels import task and interrupts the corresponding thread.
     */
    private void interruptImport() {
        this.importTask.cancel();
        this.importerThread.interrupt();
    }
    //
    /**
     * Cancels export task and interrupts the corresponding thread.
     */
    private void interruptExport() {
        this.exportTask.cancel();
        this.exporterThread.interrupt();
    }
    //
    /**
     * Gets called by the cancel fragmentation button.
     */
    private void interruptFragmentation() {
        //cancel() of the task was overridden to shut down the executor service in FragmentationService
        this.parallelFragmentationMainTask.cancel(true);
        this.cancelFragmentationButton.setVisible(false);
        this.fragmentationButton.setDisable(false);
    }
    //
    /**
     * Starts fragmentation for only one algorithm.
     */
    private void startFragmentation() {
        this.startFragmentation(false);
    }
    //
    /**
     * Starts fragmentation task and opens fragment and itemization tabs.
     */
    private void startFragmentation(boolean isPipelining) {
        long tmpStartTime = System.nanoTime();
        MainViewController.LOGGER.info("Start of method startFragmentation");
        List<MoleculeDataModel> tmpSelectedMolecules = this.moleculeDataModelList.stream().filter(MoleculeDataModel::isSelected).toList();
        int tmpNumberOfCores = this.settingsContainer.getNumberOfTasksForFragmentationSetting();
        boolean tmpIsKeepLastFragmentSetting = this.settingsContainer.isKeepLastFragmentSetting();
        boolean tmpIsStereoChemRegarded = this.settingsContainer.getRegardStereochemistrySetting();
        try {
            this.fragmentationButton.setDisable(true);
            this.cancelFragmentationButton.setVisible(true);
            this.parallelFragmentationMainTask = new Task<>() {
                @Override
                protected Void call() throws Exception {
                    if (isPipelining) {
                        MainViewController.this.fragmentationService.startPipelineFragmentation(tmpSelectedMolecules,
                                tmpNumberOfCores, tmpIsStereoChemRegarded, tmpIsKeepLastFragmentSetting);
//                        fragmentationService.startPipelineFragmentationMolByMol(tmpSelectedMolecules, tmpNumberOfCores);
                    } else {
                        MainViewController.this.fragmentationService.startSingleFragmentation(tmpSelectedMolecules,
                                tmpNumberOfCores, tmpIsStereoChemRegarded);
                    }
                    return null;
                }
                //
                @Override
                public boolean cancel(boolean anInterruptThread) {
                    MainViewController.this.fragmentationService.abortExecutor();
                    return super.cancel(anInterruptThread);
                }
            };
            this.parallelFragmentationMainTask.setOnSucceeded(event ->
                //note: setOnSucceeded() takes place in the JavaFX GUI thread again but still runLater() is necessary to wait
                // for the thread to be free for the update
                Platform.runLater(() -> {
                    try {
                        ObservableList<FragmentDataModel> tmpObservableFragments = FXCollections.observableArrayList();
                        Set<String> tmpKeys = this.fragmentationService.getFragments().keySet();
                        for (String tmpKey : tmpKeys) {
                            tmpObservableFragments.add(this.fragmentationService.getFragments().get(tmpKey));
                        }
                        this.mapOfFragmentDataModelLists.put(this.fragmentationService.getCurrentFragmentationName(), tmpObservableFragments);
                        this.addFragmentationResultTabs(this.fragmentationService.getCurrentFragmentationName());
                        this.updateStatusBar(this.fragmentationThread, Message.get("Status.finished"));
                        this.mainView.getMainMenuBar().getExportMenu().setDisable(false);
                        this.mainView.getMainMenuBar().getHistogramViewerMenuItem().setDisable(false);
                        this.mainView.getMainMenuBar().getStartClusteringMenuItem().setDisable(false);
                        this.fragmentationButton.setDisable(false);
                        this.cancelFragmentationButton.setVisible(false);
                        this.isFragmentationRunning = false;
                        long tmpEndTime = System.nanoTime();
                        MainViewController.LOGGER.info("End of method startFragmentation after " + (tmpEndTime - tmpStartTime) / 1000000000.0 + " seconds.");
                    } catch (Exception anException) {
                        MainViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
                    }
                })
            );
            this.parallelFragmentationMainTask.setOnCancelled(event -> {
                this.updateStatusBar(this.fragmentationThread, Message.get("Status.canceled"));
                //note: the export functions give a warning to the user should the fragment set be empty
                this.mainView.getMainMenuBar().getExportMenu().setDisable(false);
                this.fragmentationButton.setDisable(false);
                this.cancelFragmentationButton.setVisible(false);
                this.isFragmentationRunning = false;
                MainViewController.LOGGER.info("Fragmentation cancelled by user.");
            });
            // see also how exceptions thrown in the task's sub-tasks are handled in startFragmentation()
            this.parallelFragmentationMainTask.setOnFailed(event -> {
                this.updateStatusBar(this.fragmentationThread, Message.get("Status.failed"));
                //note: the export functions give a warning to the user should the fragment set be empty
                this.mainView.getMainMenuBar().getExportMenu().setDisable(false);
                this.fragmentationButton.setDisable(false);
                this.cancelFragmentationButton.setVisible(false);
                this.isFragmentationRunning = false;
                Exception tmpCause = (Exception) event.getSource().getException();
                MainViewController.LOGGER.log(Level.SEVERE, tmpCause.toString(), tmpCause);
                Platform.runLater(() -> {
                    GuiUtil.guiExceptionAlert(Message.get("MainViewController.FragmentationError.Title"),
                            Message.get("MainViewController.FragmentationError.Header"),
                            Message.get("MainViewController.FragmentationError.Content"),
                            tmpCause);
                });
            });
            this.fragmentationThread = new Thread(this.parallelFragmentationMainTask);
            this.fragmentationThread.setName(ThreadType.FRAGMENTATION_THREAD.getThreadName());
            this.fragmentationThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
            this.fragmentationThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number, do not touch
            this.updateStatusBar(this.fragmentationThread, Message.get("Status.running"));
            this.isFragmentationRunning = true;
            this.fragmentationThread.start();
        } catch (Exception anException) {
            MainViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("MainViewController.FragmentationError.Title"),
                    Message.get("MainViewController.FragmentationError.Header"),
                    Message.get("MainViewController.FragmentationError.Content"),
                    anException);
        }
    }
    //
    //<editor-fold desc="clustering" defaultstate="collapsed">
    /**
     * Starts the clustering of the selected molecules by the fragments of the fragmentation of the selected fragments or
     * items tab (or of the last fragmentation if another tab is selected) and the selected descriptors. The result is
     * shown in a new tab.
     */
    private void startClustering() {
        this.startClustering(null);
    }
    //
    /**
     * Starts the clustering of the selected molecules by the fragments of the given fragmentation and the selected
     * descriptors. The result is shown in a new tab.
     *
     * @param aFragmentationName name of the fragmentation to use; null to use the fragmentation of the selected fragments
     *                           or items tab (or the last fragmentation if another tab is selected)
     */
    private void startClustering(String aFragmentationName) {
        if (this.isClusteringRunningProperty.get()) {
            return;
        }
        String tmpFragmentationName = this.getFragmentationNameToCluster(aFragmentationName);
        if (tmpFragmentationName == null) {
            return;
        }
        List<MoleculeDataModel> tmpSelectedMolecules = this.getMoleculesToCluster();
        if (tmpSelectedMolecules == null) {
            return;
        }
        try {
            this.clusteringTask = this.clusteringService.createClusteringTask(tmpSelectedMolecules, tmpFragmentationName,
                    this.settingsContainer.getNumberOfTasksForFragmentationSetting() > 1);
        } catch (IOException anException) {
            MainViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("MainViewController.ClusteringError.Title"),
                    Message.get("MainViewController.ClusteringError.Header"),
                    Message.get("MainViewController.ClusteringError.Content"), anException);
            return;
        }
        this.clusteringTask.setOnSucceeded(event -> {
            this.updateStatusBar(this.clusteringThread, Message.get("Status.finished"));
            this.setClusteringRunningState(false);
            this.addClusteringResultTab(this.clusteringTask.getValue());
        });
        this.clusteringTask.setOnCancelled(event -> {
            this.updateStatusBar(this.clusteringThread, Message.get("Status.canceled"));
            this.setClusteringRunningState(false);
            MainViewController.LOGGER.info("Clustering cancelled by user.");
        });
        this.clusteringTask.setOnFailed(event -> {
            this.updateStatusBar(this.clusteringThread, Message.get("Status.failed"));
            this.setClusteringRunningState(false);
            Throwable tmpCause = event.getSource().getException();
            MainViewController.LOGGER.log(Level.SEVERE, tmpCause.toString(), tmpCause);
            GuiUtil.guiExceptionAlert(Message.get("MainViewController.ClusteringError.Title"),
                    Message.get("MainViewController.ClusteringError.Header"),
                    tmpCause.getMessage(),
                    tmpCause instanceof Exception tmpException ? tmpException : new Exception(tmpCause));
        });
        this.clusteringThread = new Thread(this.clusteringTask);
        this.clusteringThread.setName(ThreadType.CLUSTERING_THREAD.getThreadName());
        this.clusteringThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        this.clusteringThread.setDaemon(true);
        this.clusteringThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number, as for the fragmentation
        this.setClusteringRunningState(true);
        this.updateStatusBar(this.clusteringThread, Message.get("Status.clustering"));
        this.clusteringThread.start();
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
     * Cancels a running clustering. The current step (descriptor calculation or ART-2A clustering) cannot be interrupted,
     * the run stops after it and its result is discarded.
     */
    private void interruptClustering() {
        if (this.clusteringTask != null) {
            this.clusteringTask.cancel(true);
        }
        if (this.vigilanceScanTask != null) {
            this.vigilanceScanTask.cancel(true);
        }
        this.setClusteringRunningState(false);
    }
    //
    /**
     * Enables or disables the clustering menu items according to whether a clustering is running.
     *
     * @param anIsRunning true if a clustering is running
     */
    private void setClusteringRunningState(boolean anIsRunning) {
        this.isClusteringRunningProperty.set(anIsRunning);
        this.mainView.getMainMenuBar().getStartClusteringMenuItem().setDisable(anIsRunning || this.mapOfFragmentDataModelLists.isEmpty());
    }
    //
    /**
     * Opens the clustering settings; their start button starts the clustering if a fragmentation is available.
     */
    private void openClusteringSettingsView() {
        this.openClusteringSettingsView(null);
    }
    //
    /**
     * Opens the clustering settings; their start button starts the clustering with the given fragmentation.
     *
     * @param aFragmentationName name of the fragmentation to cluster with, null for the fragmentation of the selected tab
     */
    private void openClusteringSettingsView(String aFragmentationName) {
        new ClusteringSettingsViewController(this.primaryStage, this.clusteringService, this.configuration,
                this.mapOfFragmentDataModelLists.isEmpty() || this.isClusteringRunningProperty.get()
                        ? null
                        : () -> this.startClustering(aFragmentationName));
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
    private HBox createClusteringBar(String aFragmentationName) {
        return this.createClusteringBar(aFragmentationName, null,
                MainViewController.VIGILANCE_SCAN_KEY_FRAGMENTATION_PREFIX + aFragmentationName);
    }
    //
    /**
     * Creates a clustering bar that clusters with the given fragmentation: descriptor preset, settings button,
     * vigilance, start button, vigilance scan button, button to show the last scan again (only visible if there is
     * one) and cancel button.
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
        Button tmpSettingsButton = new Button();
        String tmpIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder")
                        + this.configuration.getProperty("mortar.icon.gear.name")).toExternalForm();
        tmpSettingsButton.setGraphic(new ImageView(new Image(tmpIconURL)));
        tmpSettingsButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.clusteringBar.settingsButton.tooltip")));
        tmpSettingsButton.setOnAction(anEvent -> this.openClusteringSettingsView(aFragmentationName));
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
        return tmpClusteringBar;
    }
    //
    /**
     * Adds a tab that shows the given clustering result and connects its buttons to the overview view and the export.
     *
     * @param aRunResult result of the clustering run
     */
    private void addClusteringResultTab(ClusteringRunResult aRunResult) {
        this.clusteringRunCounter++;
        ClusteringResultView tmpView = new ClusteringResultView(aRunResult);
        Tab tmpTab = new Tab(Message.get("ClusteringResultView.tabTitle") + " " + this.clusteringRunCounter
                + " - " + aRunResult.parameters().fragmentationName());
        tmpTab.setId(TabNames.CLUSTERING.toString());
        tmpTab.setClosable(true);
        tmpTab.setContent(tmpView);
        //opens the selected cluster in the overview view, the representative first
        Runnable tmpOpenOverview = () -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpView.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster == null) {
                return;
            }
            this.openClusterOverview(aRunResult, tmpCluster, this.primaryStage);
        };
        tmpView.getClusterMoleculeSelectionMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpView.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null && tmpCluster.clusterIndex() >= 0) {
                this.selectMoleculesFromCluster(aRunResult, tmpCluster.clusterIndex());
            }
        });
        tmpView.getOverviewButton().setOnAction(event -> tmpOpenOverview.run());
        tmpView.getClusterTableView().setOnMouseClicked(event -> {
            if (event.getClickCount() == 2) {
                tmpOpenOverview.run();
            }
        });
        //export, training/test partitioning and diverse subset use the marked clusters, all if none is marked
        tmpView.getExportButton().setOnAction(event -> this.exportClusteringData(aRunResult, tmpView.getMarkedClusterIndices()));
        tmpView.getTrainTestSplitButton().setOnAction(event -> this.createTrainTestSplit(aRunResult, tmpView.getMarkedClusterIndices()));
        //medians of the chosen descriptors, updated in all result tabs when the choice changes
        tmpView.setDisplayedDescriptorNames(this.clusteringService.getDisplayedDescriptorColumnNames());
        ChangeListener<String> tmpDisplayedDescriptorsListener = (anObservable, anOldValue, aNewValue) ->
                tmpView.setDisplayedDescriptorNames(this.clusteringService.getDisplayedDescriptorColumnNames());
        this.clusteringService.displayedDescriptorColumnsProperty().addListener(tmpDisplayedDescriptorsListener);
        tmpView.getChooseDescriptorsButton().setOnAction(event -> this.chooseDisplayedDescriptors(tmpView));
        tmpView.getDescriptorValuesMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpView.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null) {
                this.showDescriptorValuesOfCluster(aRunResult, tmpCluster);
            }
        });
        tmpView.getRangeFilterButton().setOnAction(event -> this.applyClusterRangeFilter(tmpView));
        //the scan window of this result is discarded at the latest when the result tab is closed
        String tmpScanKey = MainViewController.VIGILANCE_SCAN_KEY_RESULT_PREFIX + this.clusteringRunCounter;
        tmpTab.setOnClosed(event -> {
            this.clusteringService.displayedDescriptorColumnsProperty().removeListener(tmpDisplayedDescriptorsListener);
            Stage tmpScanStage = this.vigilanceScanStageMap.remove(tmpScanKey);
            if (tmpScanStage != null) {
                tmpScanStage.close();
            }
        });
        tmpView.getDiverseSubsetButton().setOnAction(event -> this.createDiverseSubset(aRunResult, tmpView.getMarkedClusterIndices()));
        //similarity search for a name or SMILES and for the representative of the selected cluster
        tmpView.getSimilaritySearchButton().setOnAction(event ->
                this.searchSimilarMolecules(tmpView.getSimilaritySearchTextField().getText(), aRunResult));
        tmpView.getSimilaritySearchMenuItem().setOnAction(event -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpView.getClusterTableView().getSelectionModel().getSelectedItem();
            if (tmpCluster != null && tmpCluster.representative() != null) {
                this.searchSimilarMolecules(tmpCluster.representative(), aRunResult);
            } else if (tmpCluster != null && !tmpCluster.members().isEmpty()) {
                this.searchSimilarMolecules(tmpCluster.members().getFirst(), aRunResult);
            }
        });
        //the overview views offer the similarity search in their context menu from now on
        this.lastClusteringRunResult = aRunResult;
        this.viewToolsManager.setOverviewSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, this.getActiveClusteringRunResult()));
        this.viewToolsManager.setOverviewDescriptorValuesAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, this.getActiveClusteringRunResult()));
        //clustering bar to cluster again with this fragmentation and changed settings; its scan reuses the features
        tmpView.setClusteringBar(this.createClusteringBar(aRunResult.parameters().fragmentationName(), aRunResult, tmpScanKey));
        tmpView.getRepresentativesOverviewButton().setOnAction(event -> this.openClusterRepresentativesOverview(aRunResult, tmpTab.getText()));
        this.mainTabPane.getTabs().add(tmpTab);
        this.mainTabPane.getSelectionModel().select(tmpTab);
        //the overview of the cluster representatives opens directly, after the success handler has finished
        Platform.runLater(() -> this.openClusterRepresentativesOverview(aRunResult, tmpTab.getText()));
        if (aRunResult.clusteringResult().isOverflow()) {
            MainViewController.LOGGER.warning("Clustering finished with cluster overflow.");
        }
    }
    //
    /**
     * Opens an overview view with the representative of every cluster, ranked from the largest to the smallest cluster,
     * with the cluster size and its share of the clustered molecules below each structure, like the fragment frequencies
     * in the fragments overview. A double click on a representative opens all molecules of its cluster in another
     * overview view. The molecules that could not be clustered are not shown.
     *
     * @param aRunResult result of the clustering run
     * @param aTitle title of the clustering result, used in the window title
     */
    private void openClusterRepresentativesOverview(ClusteringRunResult aRunResult, String aTitle) {
        List<ClusteringRunResult.ClusterSummary> tmpClusters = aRunResult.getClusterSummariesRankedBySize().stream()
                .filter(aCluster -> aCluster.clusterIndex() >= 0).toList();
        if (tmpClusters.isEmpty()) {
            return;
        }
        List<MoleculeDataModel> tmpRepresentatives = new ArrayList<>(tmpClusters.size());
        List<String> tmpLabels = new ArrayList<>(tmpClusters.size());
        List<Double> tmpFractions = new ArrayList<>(tmpClusters.size());
        DecimalFormat tmpPercentageFormat = new DecimalFormat("#.#", new DecimalFormatSymbols(Locale.getDefault()));
        int tmpNumberOfMolecules = aRunResult.molecules().size();
        for (ClusteringRunResult.ClusterSummary tmpCluster : tmpClusters) {
            tmpRepresentatives.add(tmpCluster.representative());
            double tmpFraction = (double) tmpCluster.getSize() / tmpNumberOfMolecules;
            tmpLabels.add(tmpCluster.getSize() + " (" + tmpPercentageFormat.format(100.0 * tmpFraction) + "%)");
            tmpFractions.add(tmpFraction);
        }
        if (this.clusterRepresentativesOverviewViewController == null) {
            this.clusterRepresentativesOverviewViewController = new OverviewViewController(this.configuration, this.settingsContainer);
        }
        this.clusterRepresentativesOverviewViewController.setSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, aRunResult));
        //a representative stands for its cluster, so the descriptor values of the whole cluster are shown
        Map<MoleculeDataModel, ClusteringRunResult.ClusterSummary> tmpClusterOfRepresentativeMap = new IdentityHashMap<>(tmpClusters.size());
        for (ClusteringRunResult.ClusterSummary tmpCluster : tmpClusters) {
            tmpClusterOfRepresentativeMap.put(tmpCluster.representative(), tmpCluster);
        }
        this.clusterRepresentativesOverviewViewController.setDescriptorValuesOfMoleculeAction(aMolecule -> {
            ClusteringRunResult.ClusterSummary tmpCluster = tmpClusterOfRepresentativeMap.get(aMolecule);
            if (tmpCluster != null) {
                this.showDescriptorValuesOfCluster(aRunResult, tmpCluster);
            }
        }, Message.get("OverviewView.contextMenu.descriptorValuesOfClusterMenuItem"));
        this.clusterRepresentativesOverviewViewController.initializeAndShowOverviewView(this.primaryStage,
                OverviewViewController.DataSources.CLUSTER_REPRESENTATIVES, aTitle, tmpRepresentatives, tmpLabels,
                tmpFractions,
                //opened after the mouse event has been processed and owned by the representatives overview, so that the
                // remaining mouse events do not reach a blocked window when the cluster overview is closed
                anIndex -> Platform.runLater(() -> {
                    Stage tmpOwnerStage = this.clusterRepresentativesOverviewViewController.getOverviewViewStage();
                    this.openClusterOverview(aRunResult, tmpClusters.get(anIndex),
                            tmpOwnerStage != null && tmpOwnerStage.isShowing() ? tmpOwnerStage : this.primaryStage);
                }));
    }
    //
    /**
     * Returns the clustering result the similarity search refers to: the one of the selected clustering tab, otherwise
     * the result of the last clustering run.
     *
     * @return clustering run result or null if there was no clustering yet
     */
    private ClusteringRunResult getActiveClusteringRunResult() {
        Tab tmpSelectedTab = this.mainTabPane.getSelectionModel().getSelectedItem();
        if (tmpSelectedTab != null && tmpSelectedTab.getContent() instanceof ClusteringResultView tmpView) {
            return tmpView.getRunResult();
        }
        return this.lastClusteringRunResult;
    }
    //
    /**
     * Searches the molecules most similar to the molecule given by the text: a clustered or imported molecule with this
     * name or unique SMILES, otherwise the text is parsed as SMILES.
     *
     * @param aText name or SMILES of the query
     * @param aRunResult clustering run whose feature space is searched
     */
    private void searchSimilarMolecules(String aText, ClusteringRunResult aRunResult) {
        String tmpText = aText == null ? "" : aText.trim();
        if (tmpText.isEmpty()) {
            return;
        }
        MoleculeDataModel tmpQuery = null;
        //clustered molecules first, then all imported molecules
        for (List<MoleculeDataModel> tmpList : List.of(aRunResult.molecules(), this.moleculeDataModelList)) {
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
                MainViewController.LOGGER.log(Level.INFO, "Similarity search query is neither a molecule name nor a valid SMILES: " + tmpText, anException);
                GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringSimilaritySearch.error.title"),
                        Message.get("ClusteringSimilaritySearch.invalidQuery.header"),
                        String.format(Message.get("ClusteringSimilaritySearch.invalidQuery.content"), tmpText));
                return;
            }
        }
        this.searchSimilarMolecules(tmpQuery, aRunResult);
    }
    //
    /**
     * Searches the molecules most similar to the given molecule in the feature space of the given clustering run and
     * shows them in an overview view, the query first and the similarity below every molecule (also as bar). A clustered
     * molecule keeps its own cluster; another molecule is assigned to the cluster with the most similar cluster vector.
     * A molecule that was not fragmented by the run's fragmentation is fragmented with the same fragmenter first, which
     * is only possible for single-algorithm fragmentations.
     *
     * @param aQuery query molecule
     * @param aRunResult clustering run whose feature space is searched; null if there is no clustering result
     */
    private void searchSimilarMolecules(MoleculeDataModel aQuery, ClusteringRunResult aRunResult) {
        if (aRunResult == null) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringSimilaritySearch.error.title"),
                    Message.get("ClusteringSimilaritySearch.noClustering.header"), null);
            return;
        }
        String tmpFragmentationName = aRunResult.parameters().fragmentationName();
        try {
            SimilaritySearch tmpSearch = this.similaritySearchMap.computeIfAbsent(aRunResult, SimilaritySearch::new);
            int tmpNumberOfMolecules = this.clusteringService.getNumberOfSimilarMolecules();
            int tmpRowIndex = -1;
            for (int i = 0; i < aRunResult.molecules().size(); i++) {
                if (aRunResult.molecules().get(i) == aQuery) {
                    tmpRowIndex = i;
                    break;
                }
            }
            SimilaritySearch.SearchResult tmpResult;
            if (tmpRowIndex >= 0) {
                tmpResult = tmpSearch.searchByRowIndex(tmpRowIndex, tmpNumberOfMolecules);
            } else {
                MoleculeDataModel tmpFragmentedQuery = aQuery;
                if (aQuery.getFragmentFrequencyOfSpecificFragmentation(tmpFragmentationName) == null) {
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
                tmpResult = tmpSearch.searchByRawFeatureVector(
                        FeatureMatrixBuilder.buildRawFeatureVector(tmpFragmentedQuery, aRunResult), tmpNumberOfMolecules);
            }
            this.showSimilaritySearchResult(aQuery, tmpResult, aRunResult);
        } catch (Exception anException) {
            MainViewController.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("ClusteringSimilaritySearch.error.title"),
                    Message.get("ClusteringSimilaritySearch.error.header"), anException.getMessage(), anException);
        }
    }
    //
    /**
     * Shows the result of a similarity search in a new overview view: the query first (highlighted), then the most
     * similar molecules with their similarity as text and light blue bar. The title names the query and its cluster.
     *
     * @param aQuery query molecule
     * @param aResult result of the search
     * @param aRunResult clustering run that was searched
     */
    private void showSimilaritySearchResult(MoleculeDataModel aQuery, SimilaritySearch.SearchResult aResult, ClusteringRunResult aRunResult) {
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(aResult.neighbors().size() + 1);
        List<String> tmpLabels = new ArrayList<>(aResult.neighbors().size() + 1);
        List<Double> tmpFractions = new ArrayList<>(aResult.neighbors().size() + 1);
        //the query is highlighted as first structure, its label is not shown
        tmpMolecules.add(aQuery);
        tmpLabels.add("");
        tmpFractions.add(1.0);
        DecimalFormat tmpSimilarityFormat = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.getDefault()));
        for (SimilaritySearch.Neighbor tmpNeighbor : aResult.neighbors()) {
            tmpMolecules.add(aRunResult.molecules().get(tmpNeighbor.rowIndex()));
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
        OverviewViewController tmpController = new OverviewViewController(this.configuration, this.settingsContainer);
        tmpController.setSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, aRunResult));
        tmpController.setDescriptorValuesOfMoleculeAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, aRunResult));
        Stage tmpOwner = this.getFocusedStage();
        tmpController.initializeAndShowOverviewView(tmpOwner, OverviewViewController.DataSources.SIMILARITY_SEARCH, tmpTitle,
                tmpMolecules, tmpLabels, tmpFractions, null);
    }
    //
    /**
     * Asks for the vigilance range and starts a vigilance scan in the background; the result is shown in a new window.
     * With a given clustering result, its features and parameters are reused and only the clustering is repeated. Without
     * one (scan before the first clustering), the features of the selected molecules are calculated once with the
     * current clustering settings.
     *
     * @param aFragmentationName name of the fragmentation to cluster with
     * @param aBaseRunResult clustering result whose features are reused; null to calculate the features
     * @param aScanKey key of the scan window in the scan window map; a kept scan with this key is replaced
     */
    private void startVigilanceScan(String aFragmentationName, ClusteringRunResult aBaseRunResult, String aScanKey) {
        if (this.isClusteringRunningProperty.get()) {
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
        TextField tmpStartField = new TextField("0.1");
        TextField tmpEndField = new TextField("0.9");
        TextField tmpStepField = new TextField("0.1");
        GridPane tmpGrid = new GridPane();
        tmpGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpGrid.addRow(0, new Label(Message.get("VigilanceScan.start.label")), tmpStartField);
        tmpGrid.addRow(1, new Label(Message.get("VigilanceScan.end.label")), tmpEndField);
        tmpGrid.addRow(2, new Label(Message.get("VigilanceScan.step.label")), tmpStepField);
        tmpDialog.getDialogPane().setContent(tmpGrid);
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        float[] tmpVigilanceValues;
        try {
            tmpVigilanceValues = VigilanceScan.getVigilanceValues(
                    Double.parseDouble(tmpStartField.getText().trim().replace(',', '.')),
                    Double.parseDouble(tmpEndField.getText().trim().replace(',', '.')),
                    Double.parseDouble(tmpStepField.getText().trim().replace(',', '.')));
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
                MainViewController.LOGGER.log(Level.SEVERE, anException.toString(), anException);
                GuiUtil.guiExceptionAlert(Message.get("VigilanceScan.dialog.title"), Message.get("MainViewController.ClusteringError.Header"),
                        Message.get("MainViewController.ClusteringError.Content"), anException);
                return;
            }
        }
        List<MoleculeDataModel> tmpFinalMolecules = tmpMolecules;
        ClusteringParameters tmpFinalParameters = tmpParameters;
        DescriptorPriorityTable tmpFinalPriorityTable = tmpPriorityTable;
        //holds the clustering result the scan is based on, calculated in the task if not given
        ClusteringRunResult[] tmpBaseRunResult = new ClusteringRunResult[] {aBaseRunResult};
        Task<List<VigilanceScan.ScanEntry>> tmpScanTask = new Task<>() {
            @Override
            protected List<VigilanceScan.ScanEntry> call() throws Exception {
                if (tmpBaseRunResult[0] == null) {
                    //features once, clustered with the current vigilance; only the clustering is repeated in the scan
                    tmpBaseRunResult[0] = ClusteringTask.runPipeline(tmpFinalMolecules, tmpFinalParameters, tmpFinalPriorityTable,
                            this::isCancelled, (aStep, aMessage) -> this.updateMessage(aMessage));
                }
                return VigilanceScan.scan(tmpBaseRunResult[0], tmpVigilanceValues, this::isCancelled,
                        (aDone, aTotal) -> this.updateProgress(aDone, aTotal));
            }
        };
        //the scan uses the clustering controls (start disabled, cancel visible) like a clustering
        Thread tmpScanThread = new Thread(tmpScanTask);
        tmpScanThread.setName(ThreadType.CLUSTERING_THREAD.getThreadName());
        tmpScanThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        tmpScanThread.setDaemon(true);
        tmpScanTask.setOnSucceeded(event -> {
            this.updateStatusBar(tmpScanThread, Message.get("Status.finished"));
            this.setClusteringRunningState(false);
            this.showVigilanceScanResult(aScanKey, tmpBaseRunResult[0], tmpScanTask.getValue(), aBaseRunResult != null);
        });
        tmpScanTask.setOnCancelled(event -> {
            this.updateStatusBar(tmpScanThread, Message.get("Status.canceled"));
            this.setClusteringRunningState(false);
        });
        tmpScanTask.setOnFailed(event -> {
            this.updateStatusBar(tmpScanThread, Message.get("Status.failed"));
            this.setClusteringRunningState(false);
            Throwable tmpCause = event.getSource().getException();
            MainViewController.LOGGER.log(Level.SEVERE, tmpCause.toString(), tmpCause);
            GuiUtil.guiExceptionAlert(Message.get("VigilanceScan.dialog.title"), Message.get("MainViewController.ClusteringError.Header"),
                    tmpCause.getMessage(), tmpCause instanceof Exception tmpException ? tmpException : new Exception(tmpCause));
        });
        this.vigilanceScanTask = tmpScanTask;
        this.setClusteringRunningState(true);
        this.updateStatusBar(tmpScanThread, Message.get("Status.clustering"));
        tmpScanThread.start();
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
     */
    private void showVigilanceScanResult(String aScanKey, ClusteringRunResult aBaseRunResult, List<VigilanceScan.ScanEntry> anEntryList,
                                         boolean anIsSelectingBaseVigilance) {
        VigilanceScanView tmpView = new VigilanceScanView(anEntryList);
        if (anIsSelectingBaseVigilance) {
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
                this.showExportOutcome(tmpExporter.exportVigilanceScanCsvFile(tmpDirectory, aBaseRunResult, anEntryList,
                        this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
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
    //
    /**
     * Opens a dialog to choose the descriptors whose medians are shown in the clustering result tabs: a searchable list of
     * all descriptor columns of the given result with check boxes. The choice is stored in the clustering settings and
     * applied to all result tabs; chosen descriptors that were not calculated in the given result are kept.
     *
     * @param aResultView result view whose descriptors are offered
     */
    private void chooseDisplayedDescriptors(ClusteringResultView aResultView) {
        List<String> tmpAvailableNames = aResultView.getDescriptorStatistics().getDescriptorColumnNames();
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
     * Opens the range filter for the clusters of the given result view: one or more conditions, each a descriptor with an
     * optional lower and upper bound on its median in the cluster. The number of matching clusters is updated while the
     * conditions are edited. On OK, exactly the clusters that satisfy all conditions are marked and the descriptors of
     * the conditions are added to the displayed descriptors.
     *
     * @param aResultView result view whose clusters are filtered
     */
    private void applyClusterRangeFilter(ClusteringResultView aResultView) {
        ClusterDescriptorStatistics tmpStatistics = aResultView.getDescriptorStatistics();
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
            List<ClusterDescriptorStatistics.RangeCondition> tmpConditions = MainViewController.readRangeConditions(
                    tmpDescriptorComboBoxes, tmpMinimumFields, tmpMaximumFields);
            tmpMatchLabel.setText(tmpConditions == null
                    ? Message.get("ClusteringResultView.rangeFilter.invalidInput.text")
                    : String.format(Message.get("ClusteringResultView.rangeFilter.matching.text"),
                            tmpStatistics.getMatchingClusters(tmpConditions).size(), aResultView.getRunResult().clusteringResult().numberOfClusters()));
        };
        //adds a condition row, filled with the given condition if not null
        Consumer<ClusterDescriptorStatistics.RangeCondition> tmpAddRow = aCondition -> {
            ComboBox<String> tmpComboBox = new ComboBox<>(FXCollections.observableArrayList(tmpNames));
            tmpComboBox.setValue(aCondition != null && tmpNames.contains(aCondition.columnName()) ? aCondition.columnName()
                    : (aResultView.getDisplayedDescriptorNames().isEmpty() ? tmpNames.getFirst() : aResultView.getDisplayedDescriptorNames().getFirst()));
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
        if (aResultView.getLastRangeConditions().isEmpty()) {
            tmpAddRow.accept(null);
        } else {
            aResultView.getLastRangeConditions().forEach(tmpAddRow);
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
        List<ClusterDescriptorStatistics.RangeCondition> tmpConditions = MainViewController.readRangeConditions(
                tmpDescriptorComboBoxes, tmpMinimumFields, tmpMaximumFields);
        if (tmpConditions == null) {
            return;
        }
        aResultView.setLastRangeConditions(tmpConditions);
        aResultView.setMarkedClusterIndices(tmpStatistics.getMatchingClusters(tmpConditions));
        //the medians the filter refers to are shown in the table
        List<String> tmpDisplayedNames = new ArrayList<>(this.clusteringService.getDisplayedDescriptorColumnNames());
        for (ClusterDescriptorStatistics.RangeCondition tmpCondition : tmpConditions) {
            if (!tmpDisplayedNames.contains(tmpCondition.columnName())) {
                tmpDisplayedNames.add(tmpCondition.columnName());
            }
        }
        this.clusteringService.setDisplayedDescriptorColumnNames(tmpDisplayedNames);
    }
    //
    /**
     * Asks for the number of molecules per cluster, selects the diverse subset of the given clustering result, shows it
     * in an overview view (cluster and rank below every structure) and optionally exports it as CSV.
     *
     * @param aRunResult clustering result
     * @param aClusterIndices indices of the marked clusters to select from; empty for all clusters
     */
    private void createDiverseSubset(ClusteringRunResult aRunResult, Set<Integer> aClusterIndices) {
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.primaryStage);
        tmpDialog.setTitle(Message.get("DiverseSubset.dialog.title"));
        tmpDialog.setHeaderText(MainViewController.appendMarkedClustersNote(Message.get("DiverseSubset.dialog.header"), aClusterIndices));
        Spinner<Integer> tmpNumberSpinner = new Spinner<>(1, Integer.MAX_VALUE, 1);
        tmpNumberSpinner.setEditable(true);
        ComboBox<DiverseSubsetSelector.SelectionMode> tmpModeComboBox = MainViewController.createSelectionModeComboBox();
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
        List<DiverseSubsetSelector.SelectedMolecule> tmpSelection = DiverseSubsetSelector.select(aRunResult,
                tmpNumberSpinner.getValue(), tmpUnclusteredCheckBox.isSelected(), tmpModeComboBox.getValue(), aClusterIndices);
        if (tmpExportCheckBox.isSelected()) {
            Exporter tmpExporter = new Exporter(this.settingsContainer);
            File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage,
                    Exporter.ExportTypes.DIVERSE_SUBSET_CSV_FILE, aRunResult.parameters().fragmentationName());
            if (tmpDirectory != null) {
                this.showExportOutcome(tmpExporter.exportDiverseSubsetCsvFile(tmpDirectory, aRunResult, tmpSelection,
                        this.settingsContainer.getCsvExportSeparatorSettingCharacter()));
            }
        }
        if (tmpSelection.isEmpty()) {
            return;
        }
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(tmpSelection.size());
        List<String> tmpLabels = new ArrayList<>(tmpSelection.size());
        for (DiverseSubsetSelector.SelectedMolecule tmpSelected : tmpSelection) {
            tmpMolecules.add(aRunResult.molecules().get(tmpSelected.rowIndex()));
            tmpLabels.add(tmpSelected.clusterIndex() < 0
                    ? Message.get("ClusteringResultView.notClusterable.text")
                    : String.format(Message.get("DiverseSubset.structureLabel.text"), tmpSelected.clusterIndex(), tmpSelected.rankInCluster()));
        }
        OverviewViewController tmpController = new OverviewViewController(this.configuration, this.settingsContainer);
        tmpController.setSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, aRunResult));
        tmpController.setDescriptorValuesOfMoleculeAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, aRunResult));
        tmpController.initializeAndShowOverviewView(this.primaryStage, OverviewViewController.DataSources.DIVERSE_SUBSET,
                String.format(Message.get("DiverseSubset.overview.title"), tmpSelection.size()), tmpMolecules, tmpLabels, null, null);
    }
    //
    /**
     * Opens the molecules of the given cluster in a new overview view, the representative first. Its context menu offers
     * the similarity search and the selection of typical or diverse molecules from the cluster.
     *
     * @param aRunResult clustering result
     * @param aCluster cluster to show
     * @param anOwnerStage owner of the overview window
     */
    private void openClusterOverview(ClusteringRunResult aRunResult, ClusteringRunResult.ClusterSummary aCluster, Stage anOwnerStage) {
        String tmpClusterName = aCluster.clusterIndex() < 0
                ? Message.get("ClusteringResultView.notClusterable.text")
                : Message.get("OverviewView.titleOfDataSource.cluster") + " " + aCluster.clusterIndex();
        //a new controller per cluster, so that its context menu entries refer to this clustering result
        OverviewViewController tmpController = new OverviewViewController(this.configuration, this.settingsContainer);
        tmpController.setSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, aRunResult));
        tmpController.setDescriptorValuesOfMoleculeAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, aRunResult));
        if (aCluster.clusterIndex() >= 0) {
            tmpController.setClusterMoleculeSelectionAction(aMolecule -> this.selectMoleculesFromCluster(aRunResult, aCluster.clusterIndex()));
        }
        tmpController.setDescriptorValuesOfMoleculeAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, aRunResult));
        tmpController.setDescriptorValuesOfAllMoleculesAction(() -> this.showDescriptorValuesOfCluster(aRunResult, aCluster));
        tmpController.initializeAndShowOverviewView(anOwnerStage, OverviewViewController.DataSources.CLUSTER, tmpClusterName, aCluster.members());
    }
    //
    /**
     * Shows the descriptor values of all molecules of the given cluster, the representative first, see
     * {@link #showDescriptorValues(ClusteringRunResult, List, String)}.
     *
     * @param aRunResult clustering result
     * @param aCluster cluster whose molecules are shown
     */
    private void showDescriptorValuesOfCluster(ClusteringRunResult aRunResult, ClusteringRunResult.ClusterSummary aCluster) {
        String tmpClusterName = aCluster.clusterIndex() < 0
                ? Message.get("ClusteringResultView.notClusterable.text")
                : Message.get("OverviewView.titleOfDataSource.cluster") + " " + aCluster.clusterIndex();
        this.showDescriptorValues(aRunResult, aCluster.members(), tmpClusterName);
    }
    //
    /**
     * Shows the descriptor values of the given molecule, see {@link #showDescriptorValues(ClusteringRunResult, List, String)}.
     *
     * @param aMolecule molecule
     * @param aRunResult clustering result whose values are shown; if null, the user is informed that there is none
     */
    private void showDescriptorValuesOfMolecule(MoleculeDataModel aMolecule, ClusteringRunResult aRunResult) {
        if (aRunResult == null) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.noClustering.header"), null);
            return;
        }
        this.showDescriptorValues(aRunResult, List.of(aMolecule), aMolecule.getName());
    }
    //
    /**
     * Shows the untransformed values of the displayed descriptors for the given molecules in a new window: one row per
     * molecule with name, cluster, structure and one sortable column per descriptor with fixed width; the table scrolls
     * horizontally. If all molecules belong to the same cluster, the column headers show its median as reference.
     * Descriptors removed by the filter are marked with an asterisk. Molecules that are not part of the clustering are
     * left out.
     *
     * @param aRunResult clustering result whose values are shown
     * @param aMoleculeList molecules, e.g. the members of a cluster or a single molecule
     * @param aName name of the molecules shown in the window title, e.g. the cluster or the molecule name
     */
    private void showDescriptorValues(ClusteringRunResult aRunResult, List<MoleculeDataModel> aMoleculeList, String aName) {
        ClusterDescriptorStatistics tmpStatistics = new ClusterDescriptorStatistics(aRunResult);
        List<String> tmpNames = this.clusteringService.getDisplayedDescriptorColumnNames().stream()
                .filter(tmpStatistics::hasDescriptorColumn).toList();
        if (tmpNames.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.noDescriptors.header"), null);
            return;
        }
        Map<MoleculeDataModel, Integer> tmpRowIndexMap = MainViewController.getRowIndexMap(aRunResult);
        List<MoleculeDataModel> tmpMolecules = aMoleculeList.stream().filter(tmpRowIndexMap::containsKey).toList();
        if (tmpMolecules.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("ClusteringResultView.descriptorValues.title"),
                    Message.get("ClusteringResultView.descriptorValues.notClustered.header"), null);
            return;
        }
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
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
        //the structure of the molecule of the row is depicted by the cell itself when it becomes visible, so that every
        // cell has its own image
        tmpStructureColumn.setCellFactory(aColumn -> new TableCell<>() {
            @Override
            protected void updateItem(MoleculeDataModel aMolecule, boolean anIsEmpty) {
                super.updateItem(aMolecule, anIsEmpty);
                this.setText(null);
                this.setGraphic(anIsEmpty || aMolecule == null ? null : aMolecule.getStructure());
            }
        });
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
     * Asks for the number of molecules and the selection mode and shows the selected molecules of the given cluster in a
     * new overview view, with their selection rank below every structure.
     *
     * @param aRunResult clustering result
     * @param aClusterIndex index of the cluster
     */
    private void selectMoleculesFromCluster(ClusteringRunResult aRunResult, int aClusterIndex) {
        //the window the action came from (e.g. the cluster overview) owns the dialog and the result window; the main
        // window is blocked by the modal overview, and focusing it on closing would cause a system beep
        Stage tmpOwner = this.getFocusedStage();
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(tmpOwner);
        tmpDialog.setTitle(Message.get("ClusterMoleculeSelection.dialog.title"));
        tmpDialog.setHeaderText(String.format(Message.get("ClusterMoleculeSelection.dialog.header"), aClusterIndex));
        Spinner<Integer> tmpNumberSpinner = new Spinner<>(1, Integer.MAX_VALUE, MainViewController.CLUSTER_MOLECULE_SELECTION_NUMBER_DEFAULT);
        tmpNumberSpinner.setEditable(true);
        ComboBox<DiverseSubsetSelector.SelectionMode> tmpModeComboBox = MainViewController.createSelectionModeComboBox();
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
        List<DiverseSubsetSelector.SelectedMolecule> tmpSelection = DiverseSubsetSelector.selectFromCluster(aRunResult,
                aClusterIndex, tmpNumberSpinner.getValue(), tmpModeComboBox.getValue());
        List<MoleculeDataModel> tmpMolecules = new ArrayList<>(tmpSelection.size());
        List<String> tmpLabels = new ArrayList<>(tmpSelection.size());
        for (DiverseSubsetSelector.SelectedMolecule tmpSelected : tmpSelection) {
            tmpMolecules.add(aRunResult.molecules().get(tmpSelected.rowIndex()));
            tmpLabels.add("#" + tmpSelected.rankInCluster());
        }
        OverviewViewController tmpController = new OverviewViewController(this.configuration, this.settingsContainer);
        tmpController.setSimilaritySearchAction(aMolecule -> this.searchSimilarMolecules(aMolecule, aRunResult));
        tmpController.setDescriptorValuesOfMoleculeAction(aMolecule -> this.showDescriptorValuesOfMolecule(aMolecule, aRunResult));
        tmpController.initializeAndShowOverviewView(tmpOwner, OverviewViewController.DataSources.DIVERSE_SUBSET,
                String.format(Message.get("ClusterMoleculeSelection.overview.title"), aClusterIndex, tmpSelection.size(),
                        Message.get("DiverseSubset.selectionMode." + tmpModeComboBox.getValue().name())),
                tmpMolecules, tmpLabels, null, null);
    }
    //
    /**
     * Returns the focused window of the application, e.g. an overview view from whose context menu an action was
     * started, to use it as owner of dialogs and further windows; the main window if no other window is focused.
     *
     * @return focused stage or the main stage
     */
    private Stage getFocusedStage() {
        for (Window tmpWindow : Window.getWindows()) {
            if (tmpWindow instanceof Stage tmpStage && tmpStage.isShowing() && tmpStage.isFocused()) {
                return tmpStage;
            }
        }
        return this.primaryStage;
    }
    //
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
                tmpMinimum = MainViewController.parseBound(aMinimumFieldList.get(i).getText());
                tmpMaximum = MainViewController.parseBound(aMaximumFieldList.get(i).getText());
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
    //
    /**
     * Shows a warning if an export returned failed files or an invalid directory; otherwise nothing.
     *
     * @param aFailedFileList result of an export method of the exporter
     */
    private void showExportOutcome(List<String> aFailedFileList) {
        if (aFailedFileList == null || !aFailedFileList.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                    aFailedFileList == null ? null : String.join("\n", aFailedFileList));
        }
    }
    //
    /**
     * Asks for the training/test ratio and the selection mode, creates the cluster-based training/test partitioning of
     * the given clustering result and exports it as CSV files into a directory chosen by the user.
     *
     * @param aRunResult clustering result the partitioning is based on
     * @param aClusterIndices indices of the marked clusters whose molecules are partitioned; empty for all molecules
     */
    private void createTrainTestSplit(ClusteringRunResult aRunResult, Set<Integer> aClusterIndices) {
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        tmpDialog.initOwner(this.primaryStage);
        tmpDialog.setTitle(Message.get("TrainTestSplit.dialog.title"));
        tmpDialog.setHeaderText(MainViewController.appendMarkedClustersNote(Message.get("TrainTestSplit.dialog.header"), aClusterIndices));
        TextField tmpRatioTextField = new TextField(MainViewController.TRAIN_TEST_RATIO_DEFAULT);
        tmpRatioTextField.setPromptText(MainViewController.TRAIN_TEST_RATIO_DEFAULT);
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
            tmpSplit = TrainTestSplitter.split(aRunResult, TrainTestSplitter.parseTrainingFraction(tmpRatioTextField.getText()),
                    tmpModeComboBox.getValue(), aRunResult.parameters().randomSeed(), aClusterIndices);
        } catch (IllegalArgumentException anException) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("TrainTestSplit.dialog.title"),
                    Message.get("TrainTestSplit.invalidRatio.header"), anException.getMessage());
            return;
        }
        Exporter tmpExporter = new Exporter(this.settingsContainer);
        File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage,
                Exporter.ExportTypes.TRAIN_TEST_SPLIT_CSV_FILES, aRunResult.parameters().fragmentationName());
        if (tmpDirectory == null) {
            return;
        }
        List<String> tmpFailedFiles = tmpExporter.exportTrainTestSplitCsvFiles(tmpDirectory, aRunResult, tmpSplit,
                this.settingsContainer.getCsvExportSeparatorSettingCharacter());
        if (tmpFailedFiles == null || !tmpFailedFiles.isEmpty()) {
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING, Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                    tmpFailedFiles == null ? null : String.join("\n", tmpFailedFiles));
            return;
        }
        DecimalFormat tmpPercentFormat = new DecimalFormat("0.0", new DecimalFormatSymbols(Locale.getDefault()));
        GuiUtil.guiMessageAlert(Alert.AlertType.INFORMATION, Message.get("TrainTestSplit.dialog.title"),
                Message.get("TrainTestSplit.finished.header"),
                String.format(Message.get("TrainTestSplit.finished.content"), tmpSplit.trainingIndices().length,
                        tmpSplit.testIndices().length, tmpPercentFormat.format(100 * tmpSplit.getActualTrainingFraction()),
                        tmpPercentFormat.format(100 * tmpSplit.requestedTrainingFraction()), tmpDirectory.getAbsolutePath()));
    }
    //
    /**
     * Exports the data of the given clustering run as CSV files into a directory chosen by the user.
     *
     * @param aRunResult result of the clustering run
     * @param aClusterIndices indices of the marked clusters to export; empty for all molecules
     */
    private void exportClusteringData(ClusteringRunResult aRunResult, Set<Integer> aClusterIndices) {
        Exporter tmpExporter = new Exporter(this.settingsContainer);
        if (this.isExportRunningProperty.get()) {
            this.interruptExport();
        }
        //returns null if the directory chooser dialog was cancelled
        File tmpDirectory = tmpExporter.openFileChooserForExportFileOrDir(this.primaryStage,
                Exporter.ExportTypes.CLUSTERING_CSV_FILES, aRunResult.parameters().fragmentationName());
        if (tmpDirectory == null) {
            return;
        }
        char tmpSeparator = this.settingsContainer.getCsvExportSeparatorSettingCharacter();
        this.exportTask = new Task<>() {
            @Override
            protected List<String> call() {
                return tmpExporter.exportClusteringCsvFiles(tmpDirectory, aRunResult, aClusterIndices, tmpSeparator);
            }
        };
        this.exportTask.setOnSucceeded(event -> {
            this.isExportRunningProperty.setValue(false);
            this.updateStatusBar(this.exporterThread, Message.get("Status.finished"));
            List<String> tmpFailedFiles = this.exportTask.getValue();
            if (tmpFailedFiles == null || !tmpFailedFiles.isEmpty()) {
                GuiUtil.guiMessageAlert(Alert.AlertType.WARNING,
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                        Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                        tmpFailedFiles == null ? null : String.join("\n", tmpFailedFiles));
            }
        });
        this.exportTask.setOnFailed(event -> {
            this.isExportRunningProperty.setValue(false);
            MainViewController.LOGGER.log(Level.WARNING, event.getSource().getException().toString(), event.getSource().getException());
            this.updateStatusBar(this.exporterThread, Message.get("Status.failed"));
            GuiUtil.guiMessageAlert(Alert.AlertType.WARNING,
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.title"),
                    Message.get("Exporter.FragmentsTab.ExportNotPossible.header"),
                    null);
        });
        this.exporterThread = new Thread(this.exportTask);
        this.exporterThread.setName(ThreadType.EXPORT_THREAD.getThreadName());
        this.exporterThread.setUncaughtExceptionHandler(LogUtil.getUncaughtExceptionHandler());
        this.exporterThread.setPriority(Thread.currentThread().getPriority() - 2); //magic number
        this.isExportRunningProperty.setValue(true);
        this.updateStatusBar(this.exporterThread, Message.get("Status.exporting"));
        this.exporterThread.start();
    }
    //</editor-fold>
    //
    /**
     * Adds a tab for fragments and a tab for items (results of fragmentation).
     *
     * @param aFragmentationName name of the fragmentation process
     */
    private void addFragmentationResultTabs(String aFragmentationName) {
        //fragments tab
        Tab tmpFragmentsTab = this.createFragmentsTab(aFragmentationName);
        //itemization tab
        Tab tmpItemsTab = this.createItemsTab(aFragmentationName);
        this.mainTabPane.getSelectionModel().select(tmpFragmentsTab);
    }
    //
    /**
     * Creates and returns a tab, which visualizes the resulting fragments of the fragmentation with given name.
     *
     * @param aFragmentationName String, unique name for fragmentation job
     * @return Tab
     */
    private Tab createFragmentsTab(String aFragmentationName){
        FragmentsDataTableView tmpFragmentsDataTableView = new FragmentsDataTableView(this.configuration);
        GridTabForTableView tmpFragmentsTab = new GridTabForTableView(Message.get("MainTabPane.fragmentsTab.title") + " - " + aFragmentationName, TabNames.FRAGMENTS.name(), tmpFragmentsDataTableView);
        this.mainTabPane.getTabs().add(tmpFragmentsTab);
        ObservableList<MoleculeDataModel> tmpList = FXCollections.observableArrayList(this.mapOfFragmentDataModelLists.get(aFragmentationName));
        for (MoleculeDataModel tmpMoleculeDataModel : tmpList) {
            tmpMoleculeDataModel.setStructureImageWidth(tmpFragmentsDataTableView.getStructureColumn().getWidth());
        }
        tmpFragmentsDataTableView.setItemsList(tmpList);
        Pagination tmpPagination = this.createPaginationWithSuitablePageCount(tmpList.size());
        tmpPagination.setPageFactory(pageIndex -> tmpFragmentsDataTableView.createFragmentsTableViewPage(pageIndex, this.settingsContainer));
        tmpFragmentsTab.addPaginationToGridPane(tmpPagination);
        Button tmpExportCsvButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.fragments.buttonCSV.txt"));
        tmpExportCsvButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.fragments.buttonCSV.tooltip")));
        Button tmpExportPdfButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.fragments.buttonPDF.txt"));
        tmpExportPdfButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.fragments.buttonPDF.tooltip")));
        Button tmpCancelExportButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.fragments.buttonCancelExport.txt"));
        tmpCancelExportButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.fragments.buttonCancelExport.tooltip")));
        tmpCancelExportButton.visibleProperty().bind(this.isExportRunningProperty);
        HBox tmpExportButtonsHBox = new HBox();
        tmpExportButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpExportButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpExportButtonsHBox.setAlignment(Pos.CENTER_LEFT);
        tmpExportButtonsHBox.getChildren().addAll(tmpExportCsvButton, tmpExportPdfButton, tmpCancelExportButton);
        tmpFragmentsTab.addNodeToGridPane(tmpExportButtonsHBox, 0, 1, 1, 1);
        tmpFragmentsTab.addNodeAsBottomRow(this.createClusteringBar(aFragmentationName));
        tmpExportPdfButton.setOnAction(event -> this.exportFile(Exporter.ExportTypes.FRAGMENT_PDF_FILE));
        tmpExportCsvButton.setOnAction(event -> this.exportFile(Exporter.ExportTypes.FRAGMENT_CSV_FILE));
        tmpCancelExportButton.setOnAction(event -> this.interruptExport());
        HBox tmpViewButtonsHBox = new HBox();
        tmpViewButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpViewButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpViewButtonsHBox.setAlignment(Pos.CENTER_RIGHT);
        tmpViewButtonsHBox.setMaxWidth(GuiDefinitions.GUI_GRIDPANE_FOR_NODE_ALIGNMENT_THIRD_COL_WIDTH);
        Button tmpOpenOverviewViewButton = GuiUtil.getButtonOfStandardSize(Message.get("MainView.showOverviewViewButton.text"));
        tmpOpenOverviewViewButton.setTooltip(GuiUtil.createTooltip(Message.get("MainView.showOverviewViewButton.tooltip")));
        Button tmpOpenHistogramViewButton = GuiUtil.getButtonOfStandardSize(Message.get("MainView.showHistogramViewButton.text"));
        tmpOpenHistogramViewButton.setTooltip(GuiUtil.createTooltip(Message.get("MainView.showHistogramViewButton.tooltip")));
        tmpViewButtonsHBox.getChildren().addAll(tmpOpenOverviewViewButton, tmpOpenHistogramViewButton);
        tmpFragmentsTab.addNodeToGridPane(tmpViewButtonsHBox, 2, 1, 1, 1);
        tmpOpenOverviewViewButton.setOnAction(event -> this.openOverviewView(OverviewViewController.DataSources.FRAGMENTS_TAB));
        tmpOpenHistogramViewButton.setOnAction(event -> this.openHistogramView());
        if (tmpList.isEmpty()) {
            tmpOpenOverviewViewButton.setDisable(true);
            tmpOpenHistogramViewButton.setDisable(true);
        }
        int tmpRowsPerPage = this.settingsContainer.getRowsPerPageSetting();
        tmpFragmentsDataTableView.setOnSort((EventHandler<SortEvent<TableView>>) event -> GuiUtil.sortTableViewGlobally(event, tmpPagination, tmpRowsPerPage));
        tmpFragmentsDataTableView.widthProperty().addListener((observable, oldValue, newValue) -> {
            for (Object tmpObject : tmpFragmentsDataTableView.getItems()) {
                ((MoleculeDataModel) tmpObject).setStructureImageWidth(tmpFragmentsDataTableView.getStructureColumn().getWidth());
                ((FragmentDataModel) tmpObject).getFirstParentMolecule().setStructureImageWidth(tmpFragmentsDataTableView.getParentMolColumn().getWidth());
            }
        });
        tmpFragmentsDataTableView.addTableViewHeightListener(this.settingsContainer);
        tmpFragmentsDataTableView.getCopyMenuItem().setOnAction(event -> GuiUtil.copySelectedTableViewCellsToClipboard(tmpFragmentsDataTableView));
        tmpFragmentsDataTableView.getOverviewViewMenuItem().setOnAction(event -> this.openOverviewView(OverviewViewController.DataSources.PARENT_MOLECULES_SAMPLE));
        tmpFragmentsDataTableView.setOnKeyPressed(event -> {
            if (GuiDefinitions.KEY_CODE_COPY.match(event)) {
                GuiUtil.copySelectedTableViewCellsToClipboard(tmpFragmentsDataTableView);
            }
        });
        return tmpFragmentsTab;
    }
    //
    /**
     * Creates and returns a tab which visualizes the resulting fragments of each molecule that has undergone the
     * fragmentation with the given name.
     *
     * @param aFragmentationName String, unique name for the fragmentation job
     * @return Tab
     */
    private Tab createItemsTab(String aFragmentationName){
        ItemizationDataTableView tmpItemizationDataTableView = new ItemizationDataTableView(aFragmentationName, this.configuration);
        tmpItemizationDataTableView.setItemsList(
                //developers note: a modifiable list is needed for sorting, so don't let SonarCloud tell you that the Collectors are not needed here!
                this.moleculeDataModelList.stream().filter(x -> x.hasMoleculeUndergoneSpecificFragmentation(aFragmentationName)).collect(Collectors.toList()));
        GridTabForTableView tmpItemizationTab = new GridTabForTableView(Message.get("MainTabPane.itemizationTab.title") + " - " + aFragmentationName, TabNames.ITEMIZATION.name(), tmpItemizationDataTableView);
        this.mainTabPane.getTabs().add(tmpItemizationTab);
        Pagination tmpPagination = this.createPaginationWithSuitablePageCount(tmpItemizationDataTableView.getItemsList().size());
        tmpPagination.setPageFactory(pageIndex -> tmpItemizationDataTableView.createItemizationTableViewPage(pageIndex, aFragmentationName, this.settingsContainer));
        tmpItemizationTab.addPaginationToGridPane(tmpPagination);
        Button tmpItemizationTabExportPDfButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.itemizationTab.pdfButton.txt"));
        tmpItemizationTabExportPDfButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.itemizationTab.pdfButton.tooltip")));
        Button tmpItemizationExportCsvButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.itemizationTab.csvButton.txt"));
        tmpItemizationExportCsvButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.itemizationTab.csvButton.tooltip")));
        Button tmpCancelExportButton = GuiUtil.getButtonOfStandardSize(Message.get("MainTabPane.fragments.buttonCancelExport.txt"));
        tmpCancelExportButton.setTooltip(GuiUtil.createTooltip(Message.get("MainTabPane.fragments.buttonCancelExport.tooltip")));
        tmpCancelExportButton.visibleProperty().bind(this.isExportRunningProperty);
        tmpItemizationExportCsvButton.setOnAction(event -> this.exportFile(Exporter.ExportTypes.ITEM_CSV_FILE));
        tmpItemizationTabExportPDfButton.setOnAction(event -> this.exportFile(Exporter.ExportTypes.ITEM_PDF_FILE));
        tmpCancelExportButton.setOnAction(event -> this.interruptExport());
        HBox tmpExportButtonsHBox = new HBox();
        tmpExportButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpExportButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpExportButtonsHBox.setAlignment(Pos.CENTER_LEFT);
        tmpExportButtonsHBox.getChildren().addAll(tmpItemizationExportCsvButton, tmpItemizationTabExportPDfButton, tmpCancelExportButton);
        tmpItemizationTab.addNodeToGridPane(tmpExportButtonsHBox, 0, 1, 1, 1);
        tmpItemizationTab.addNodeAsBottomRow(this.createClusteringBar(aFragmentationName));
        HBox tmpViewButtonsHBox = new HBox();
        tmpViewButtonsHBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE));
        tmpViewButtonsHBox.setSpacing(GuiDefinitions.GUI_SPACING_VALUE);
        tmpViewButtonsHBox.setAlignment(Pos.CENTER_RIGHT);
        tmpViewButtonsHBox.setMaxWidth(GuiDefinitions.GUI_GRIDPANE_FOR_NODE_ALIGNMENT_THIRD_COL_WIDTH);
        Button tmpOpenHistogramViewButton = GuiUtil.getButtonOfStandardSize(Message.get("MainView.showHistogramViewButton.text"));
        tmpOpenHistogramViewButton.setTooltip(GuiUtil.createTooltip(Message.get("MainView.showHistogramViewButton.tooltip")));
        tmpViewButtonsHBox.getChildren().add(tmpOpenHistogramViewButton);
        tmpItemizationTab.addNodeToGridPane(tmpViewButtonsHBox, 2, 1, 1, 1);
        tmpOpenHistogramViewButton.setOnAction(event -> this.openHistogramView());
        int tmpRowsPerPage = this.settingsContainer.getRowsPerPageSetting();
        tmpItemizationDataTableView.setOnSort((EventHandler<SortEvent<TableView>>) event -> GuiUtil.sortTableViewGlobally(event, tmpPagination, tmpRowsPerPage));
        tmpItemizationDataTableView.widthProperty().addListener((observable, oldValue, newValue) -> {
            for (Object tmpObject : tmpItemizationDataTableView.getItems()) {
                ((MoleculeDataModel) tmpObject).setStructureImageWidth(tmpItemizationDataTableView.getMoleculeStructureColumn().getWidth());
            }
        });
        tmpItemizationDataTableView.addTableViewHeightListener(this.settingsContainer);
        tmpItemizationDataTableView.getCopyMenuItem().setOnAction(event -> GuiUtil.copySelectedTableViewCellsToClipboard(tmpItemizationDataTableView));
        tmpItemizationDataTableView.getOverviewViewMenuItem().setOnAction(event -> this.openOverviewView(OverviewViewController.DataSources.ITEM_WITH_FRAGMENTS_SAMPLE));
        tmpItemizationDataTableView.setOnKeyPressed(event -> {
            if (GuiDefinitions.KEY_CODE_COPY.match(event)) {
                GuiUtil.copySelectedTableViewCellsToClipboard(tmpItemizationDataTableView);
            }
        });
        if (this.mapOfFragmentDataModelLists.get(aFragmentationName).isEmpty()) {
            tmpOpenHistogramViewButton.setDisable(true);
        }
        return tmpItemizationTab;
    }
    //
    /**
     * Clears the gui and all collections.
     */
    private void clearGuiAndCollections() {
        if (this.clusteringTask != null && this.clusteringTask.isRunning()) {
            this.interruptClustering();
        }
        this.mainView.getMainMenuBar().getStartClusteringMenuItem().setDisable(true);
        //removing all tabs does not fire their closed events, so the kept scan windows and searches are discarded here
        this.vigilanceScanStageMap.values().forEach(Stage::close);
        this.vigilanceScanStageMap.clear();
        this.similaritySearchMap.clear();
        this.lastClusteringRunResult = null;
        this.moleculeDataModelList.clear();
        this.mapOfFragmentDataModelLists.clear();
        this.moleculesDataTableView = null;
        this.mainTabPane.getTabs().clear();
    }
    //
    /**
     * Returns the items list of the table view of the selected tab.
     *
     * @param aTabName Enum which specifies which kind of tab
     * @return List {@literal <}MoleculeDataModel{@literal >}
     */
    private List<MoleculeDataModel> getItemsListOfSelectedFragmentationByTabId(TabNames aTabName) {
        //only table tabs belong to a fragmentation, clustering result tabs are skipped
        if (!(this.mainTabPane.getSelectionModel().getSelectedItem() instanceof GridTabForTableView tmpSelectedGridTab)) {
            return new ArrayList<>();
        }
        GridTabForTableView tmpSelectedTab =  (GridTabForTableView) (this.mainTabPane.getTabs().stream().filter(tab ->
                tab instanceof GridTabForTableView tmpGridTab
                        && tmpSelectedGridTab.getFragmentationNameOutOfTitle().equals(tmpGridTab.getFragmentationNameOutOfTitle())
                        && aTabName.name().equals(tab.getId())
        ).findFirst().orElse(null));
        if (tmpSelectedTab == null) {
            return new ArrayList<>();
        } else {
            return ((IDataTableView) tmpSelectedTab.getTableView()).getItemsList();
        }
    }
    //
    /**
     * Updates StatusBar.
     *
     * @param aThread  Thread which was started or ended
     * @param aMessage String message to display in StatusBar
     */
    private void updateStatusBar(Thread aThread, String aMessage) {
        if (!this.threadList.contains(aThread)) {
            this.threadList.add(aThread);
            this.mainView.getStatusBar().getStatusLabel().setText(aMessage);
            this.mainView.getStatusBar().getStatusLabel().setVisible(true);
            this.mainView.getStatusBar().getProgressBar().setVisible(true);
            //return;
        } else {
            this.threadList.remove(aThread);
            if (this.threadList.isEmpty()) {
                this.mainView.getStatusBar().getProgressBar().setVisible(false);
                this.mainView.getStatusBar().getStatusLabel().setText(aMessage);
                return;
            }
            this.mainView.getStatusBar().getStatusLabel().setText(
                    this.getStatusMessageByThreadType(
                            Objects.requireNonNull(ThreadType.get(
                                    this.threadList.getLast().getName()
                            ))
                    )
            );
        }
    }
    //
    /**
     * Returns status message as string by given ThreadType.
     *
     * @param aThreadType ThreadType
     * @return String status message
     */
    private String getStatusMessageByThreadType(ThreadType aThreadType) {
        return switch (aThreadType) {
            case FRAGMENTATION_THREAD -> Message.get("Status.running");
            case IMPORT_THREAD -> Message.get("Status.importing");
            case EXPORT_THREAD -> Message.get("Status.exporting");
            case CLUSTERING_THREAD -> Message.get("Status.clustering");
            default -> "Could not find message";
        };
    }
    //</editor-fold>
    //
    //<editor-fold desc="public enum" defaultstate="collapsed">
    /**
     * Enum for different thread types, set as thread name.
     */
    public enum ThreadType {
        /**
         * enum value for fragmentation thread
         */
        FRAGMENTATION_THREAD("Fragmentation_Thread"),
        /**
         * enum value for clustering thread
         */
        CLUSTERING_THREAD("Clustering_Thread"),
        /**
         * enum value for import thread
         */
        IMPORT_THREAD("Import_Thread"),
        /**
         * enum value for export thread
         */
        EXPORT_THREAD("Export_Thread");

        private final String threadName;

        ThreadType(String aThreadName) {
            this.threadName = aThreadName;
        }
        //
        /**
         * Returns the name of this thread type
         *
         * @return String
         */
        public String getThreadName() {
            return this.threadName;
        }
        //
        /**
         * Reverse lookup
         * Returns ThreadType by given thread name;
         * returns null if thread name does not correspond to any ThreadType
         *
         * @param aThreadName String
         * @return ThreadType
         */
        public static ThreadType get(String aThreadName) {
            for (ThreadType aType : values()) {
                if (aType.threadName.equals(aThreadName)) {
                    return aType;
                }
            }
            return null;
        }
    }
    //</editor-fold>
}
