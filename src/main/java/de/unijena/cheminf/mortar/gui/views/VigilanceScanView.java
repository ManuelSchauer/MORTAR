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
import de.unijena.cheminf.mortar.model.clustering.VigilanceScan;

import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Pos;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableRow;
import javafx.scene.control.TableView;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;

import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.ToIntFunction;

/**
 * View that shows the result of a vigilance scan: a line chart of the number of clusters and of single-molecule clusters
 * against the vigilance and a table with all values per vigilance (clusters, singletons, largest cluster, not clusterable
 * molecules, mean similarity to the cluster vector, overflow, convergence). The buttons take over the clustering of the
 * selected vigilance (also by double click on a row), export the scan, discard it or close the window, which keeps the
 * scan; their actions are set by the controller. The result of a search for a target number of clusters is shown in
 * the same way, with the best entry highlighted and selected, a summary at the top and the target as line in the chart.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class VigilanceScanView extends BorderPane {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Tolerance for comparing vigilance values, which are rounded to 4 decimals in the scan.
     */
    private static final float VIGILANCE_TOLERANCE = 1.0e-5f;
    /**
     * Style of the row of the best entry of a target search.
     */
    private static final String HIGHLIGHTED_ROW_STYLE = "-fx-background-color: #ffe9a8; -fx-font-weight: bold;";
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Table of the scan entries.
     */
    private final TableView<VigilanceScan.ScanEntry> scanTableView;
    /**
     * Chart of the number of clusters and singletons.
     */
    private final LineChart<Number, Number> chart;
    /**
     * Hint at the top, replaced by the summary of a target search.
     */
    private final Label hintLabel;
    /**
     * Button to take over the clustering of the selected vigilance.
     */
    private final Button adoptButton;
    /**
     * Button to export the scan as CSV.
     */
    private final Button exportButton;
    /**
     * Button to discard the scan and close the window.
     */
    private final Button discardButton;
    /**
     * Button to close (hide) the window; the scan is kept.
     */
    private final Button closeButton;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Highlighted entry, the best entry of a target search; null if none.
     */
    private VigilanceScan.ScanEntry highlightedEntry;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, builds the view for the given scan entries.
     *
     * @param anEntryList entries of the scan in ascending order of the vigilance
     * @throws NullPointerException if the list is null
     */
    public VigilanceScanView(List<VigilanceScan.ScanEntry> anEntryList) throws NullPointerException {
        super();
        Objects.requireNonNull(anEntryList, "Given entry list is null.");
        //chart: number of clusters and of singletons against the vigilance
        NumberAxis tmpXAxis = new NumberAxis();
        tmpXAxis.setLabel(Message.get("VigilanceScan.chart.xAxis.label"));
        tmpXAxis.setForceZeroInRange(false);
        NumberAxis tmpYAxis = new NumberAxis();
        tmpYAxis.setLabel(Message.get("VigilanceScan.chart.yAxis.label"));
        LineChart<Number, Number> tmpChart = new LineChart<>(tmpXAxis, tmpYAxis);
        this.chart = tmpChart;
        tmpChart.setTitle(Message.get("VigilanceScan.chart.title"));
        tmpChart.setAnimated(false);
        XYChart.Series<Number, Number> tmpClusterSeries = new XYChart.Series<>();
        tmpClusterSeries.setName(Message.get("VigilanceScan.clustersColumn.header"));
        XYChart.Series<Number, Number> tmpSingletonSeries = new XYChart.Series<>();
        tmpSingletonSeries.setName(Message.get("VigilanceScan.singletonsColumn.header"));
        for (VigilanceScan.ScanEntry tmpEntry : anEntryList) {
            tmpClusterSeries.getData().add(new XYChart.Data<>(tmpEntry.vigilance(), tmpEntry.numberOfClusters()));
            tmpSingletonSeries.getData().add(new XYChart.Data<>(tmpEntry.vigilance(), tmpEntry.numberOfSingletons()));
        }
        tmpChart.getData().addAll(List.of(tmpClusterSeries, tmpSingletonSeries));
        //table with all values
        this.scanTableView = new TableView<>(FXCollections.observableArrayList(anEntryList));
        this.scanTableView.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
        this.scanTableView.getColumns().addAll(List.of(
                VigilanceScanView.createVigilanceColumn(),
                VigilanceScanView.createIntegerColumn("VigilanceScan.clustersColumn.header", VigilanceScan.ScanEntry::numberOfClusters),
                VigilanceScanView.createIntegerColumn("VigilanceScan.singletonsColumn.header", VigilanceScan.ScanEntry::numberOfSingletons),
                VigilanceScanView.createIntegerColumn("VigilanceScan.largestColumn.header", VigilanceScan.ScanEntry::largestClusterSize),
                VigilanceScanView.createIntegerColumn("VigilanceScan.unclusteredColumn.header", VigilanceScan.ScanEntry::numberOfUnclusteredMolecules),
                VigilanceScanView.createColumn("VigilanceScan.similarityColumn.header",
                        anEntry -> String.format(Locale.getDefault(), "%.3f", anEntry.meanSimilarityToClusterVector())),
                VigilanceScanView.createColumn("VigilanceScan.overflowColumn.header",
                        anEntry -> anEntry.isOverflow() ? Message.get("VigilanceScan.yes.text") : ""),
                VigilanceScanView.createColumn("VigilanceScan.convergedColumn.header",
                        anEntry -> anEntry.isConverged() ? Message.get("VigilanceScan.yes.text") : Message.get("VigilanceScan.no.text"))));
        SplitPane tmpSplitPane = new SplitPane(tmpChart, this.scanTableView);
        tmpSplitPane.setOrientation(Orientation.VERTICAL);
        tmpSplitPane.setDividerPositions(0.5);
        this.setCenter(tmpSplitPane);
        this.hintLabel = new Label(Message.get("VigilanceScan.hint.text"));
        this.hintLabel.setWrapText(true);
        this.hintLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.setTop(this.hintLabel);
        //the best entry of a target search is highlighted, unless it is selected
        this.scanTableView.setRowFactory(aTable -> new TableRow<>() {
            @Override
            protected void updateItem(VigilanceScan.ScanEntry anEntry, boolean anIsEmpty) {
                super.updateItem(anEntry, anIsEmpty);
                this.setStyle(!anIsEmpty && anEntry != null && anEntry == VigilanceScanView.this.highlightedEntry
                        ? VigilanceScanView.HIGHLIGHTED_ROW_STYLE : "");
            }
        });
        //buttons for the selected vigilance on the left, for the window on the right
        this.adoptButton = new Button(Message.get("VigilanceScan.adoptButton.text"));
        this.adoptButton.setTooltip(new Tooltip(Message.get("VigilanceScan.adoptButton.tooltip")));
        this.adoptButton.setMinWidth(Region.USE_PREF_SIZE);
        this.adoptButton.setDefaultButton(true);
        this.adoptButton.disableProperty().bind(this.scanTableView.getSelectionModel().selectedItemProperty().isNull());
        this.exportButton = new Button(Message.get("VigilanceScan.exportButton.text"));
        this.exportButton.setMinWidth(Region.USE_PREF_SIZE);
        this.discardButton = new Button(Message.get("VigilanceScan.discardButton.text"));
        this.discardButton.setTooltip(new Tooltip(Message.get("VigilanceScan.discardButton.tooltip")));
        this.discardButton.setMinWidth(Region.USE_PREF_SIZE);
        this.closeButton = new Button(Message.get("VigilanceScan.closeButton.text"));
        this.closeButton.setTooltip(new Tooltip(Message.get("VigilanceScan.closeButton.tooltip")));
        this.closeButton.setMinWidth(Region.USE_PREF_SIZE);
        this.closeButton.setPrefWidth(GuiDefinitions.GUI_BUTTON_WIDTH_VALUE);
        Region tmpSpacer = new Region();
        HBox.setHgrow(tmpSpacer, Priority.ALWAYS);
        HBox tmpButtonBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, this.adoptButton, this.exportButton, tmpSpacer,
                this.discardButton, this.closeButton);
        tmpButtonBox.setAlignment(Pos.CENTER_LEFT);
        tmpButtonBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.setBottom(tmpButtonBox);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the table of the scan entries.
     *
     * @return scan table view
     */
    public TableView<VigilanceScan.ScanEntry> getScanTableView() {
        return this.scanTableView;
    }
    //
    /**
     * Returns the button that takes over the clustering of the selected vigilance.
     *
     * @return adopt button
     */
    public Button getAdoptButton() {
        return this.adoptButton;
    }
    //
    /**
     * Returns the button that exports the scan as CSV.
     *
     * @return export button
     */
    public Button getExportButton() {
        return this.exportButton;
    }
    //
    /**
     * Returns the button that discards the scan and closes the window.
     *
     * @return discard button
     */
    public Button getDiscardButton() {
        return this.discardButton;
    }
    //
    /**
     * Returns the button that closes (hides) the window and keeps the scan.
     *
     * @return close button
     */
    public Button getCloseButton() {
        return this.closeButton;
    }
    //
    /**
     * Shows the result of a search for a target number of clusters: the best entry is highlighted and selected, the
     * summary replaces the hint at the top and the target is drawn as horizontal line in the chart.
     *
     * @param aSearchResult result of the target search, whose entries this view shows
     * @throws NullPointerException if the result is null
     */
    public void showTargetSearchResult(VigilanceScan.TargetSearchResult aSearchResult) throws NullPointerException {
        Objects.requireNonNull(aSearchResult, "Given search result is null.");
        VigilanceScan.ScanEntry tmpBest = aSearchResult.bestEntry();
        this.highlightedEntry = tmpBest;
        this.scanTableView.refresh();
        this.scanTableView.getSelectionModel().select(tmpBest);
        this.scanTableView.scrollTo(tmpBest);
        String tmpVigilance = VigilanceScan.formatVigilance(tmpBest.vigilance());
        this.hintLabel.setText(String.format(Message.get(aSearchResult.isExact()
                        ? "VigilanceScan.targetSearch.exact.text" : "VigilanceScan.targetSearch.closest.text"),
                aSearchResult.targetNumberOfClusters(), tmpVigilance, tmpBest.numberOfClusters(), aSearchResult.entries().size()));
        this.hintLabel.setStyle("-fx-font-weight: bold;");
        //the target as horizontal line over the searched range
        XYChart.Series<Number, Number> tmpTargetSeries = new XYChart.Series<>();
        tmpTargetSeries.setName(Message.get("VigilanceScan.targetSeries.name"));
        List<VigilanceScan.ScanEntry> tmpEntries = aSearchResult.entries();
        tmpTargetSeries.getData().add(new XYChart.Data<>(tmpEntries.getFirst().vigilance(), aSearchResult.targetNumberOfClusters()));
        tmpTargetSeries.getData().add(new XYChart.Data<>(tmpEntries.getLast().vigilance(), aSearchResult.targetNumberOfClusters()));
        this.chart.getData().add(tmpTargetSeries);
    }
    //
    /**
     * Selects the row of the given vigilance, e.g. the vigilance of the clustering the scan was started from; nothing is
     * selected if the scan does not contain it.
     *
     * @param aVigilance vigilance to select
     */
    public void selectVigilance(float aVigilance) {
        for (VigilanceScan.ScanEntry tmpEntry : this.scanTableView.getItems()) {
            if (Math.abs(tmpEntry.vigilance() - aVigilance) < VigilanceScanView.VIGILANCE_TOLERANCE) {
                this.scanTableView.getSelectionModel().select(tmpEntry);
                this.scanTableView.scrollTo(tmpEntry);
                return;
            }
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Creates a text column.
     *
     * @param aHeaderKey message key of the header
     * @param aValueFunction function returning the text of an entry
     * @return table column
     */
    private static TableColumn<VigilanceScan.ScanEntry, String> createColumn(String aHeaderKey,
                                                                              Function<VigilanceScan.ScanEntry, String> aValueFunction) {
        TableColumn<VigilanceScan.ScanEntry, String> tmpColumn = new TableColumn<>(Message.get(aHeaderKey));
        tmpColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(aValueFunction.apply(aCellData.getValue())));
        return tmpColumn;
    }
    //
    /**
     * Creates the vigilance column: sorted by value, shown with as many decimals as the vigilance has (at most four).
     *
     * @return vigilance column
     */
    private static TableColumn<VigilanceScan.ScanEntry, Float> createVigilanceColumn() {
        TableColumn<VigilanceScan.ScanEntry, Float> tmpColumn = new TableColumn<>(Message.get("VigilanceScan.vigilanceColumn.header"));
        tmpColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aCellData.getValue().vigilance()));
        tmpColumn.setCellFactory(aColumn -> new TableCell<>() {
            @Override
            protected void updateItem(Float aVigilance, boolean anIsEmpty) {
                super.updateItem(aVigilance, anIsEmpty);
                this.setText(anIsEmpty || aVigilance == null ? null : VigilanceScan.formatVigilance(aVigilance));
            }
        });
        return tmpColumn;
    }
    //
    /**
     * Creates an integer column, which sorts numerically.
     *
     * @param aHeaderKey message key of the header
     * @param aValueFunction function returning the value of an entry
     * @return table column
     */
    private static TableColumn<VigilanceScan.ScanEntry, Integer> createIntegerColumn(String aHeaderKey,
                                                                                      ToIntFunction<VigilanceScan.ScanEntry> aValueFunction) {
        TableColumn<VigilanceScan.ScanEntry, Integer> tmpColumn = new TableColumn<>(Message.get(aHeaderKey));
        tmpColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(aValueFunction.applyAsInt(aCellData.getValue())));
        return tmpColumn;
    }
    //</editor-fold>
}
