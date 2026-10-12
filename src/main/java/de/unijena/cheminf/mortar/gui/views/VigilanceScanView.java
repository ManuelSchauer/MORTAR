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
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;

import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.ToIntFunction;

/**
 * View that shows the result of a vigilance scan: a line chart of the number of clusters and of single-molecule clusters
 * against the vigilance and a table with all values per vigilance (clusters, singletons, largest cluster, not clusterable
 * molecules, mean similarity to the cluster vector, overflow, convergence). The buttons open the clustering of the
 * selected vigilance as result tab, take over its vigilance into the clustering settings and export the scan; their
 * actions are set by the controller.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class VigilanceScanView extends BorderPane {
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Table of the scan entries.
     */
    private final TableView<VigilanceScan.ScanEntry> scanTableView;
    /**
     * Button to open the clustering of the selected vigilance as result tab.
     */
    private final Button openResultButton;
    /**
     * Button to take over the selected vigilance into the clustering settings.
     */
    private final Button useVigilanceButton;
    /**
     * Button to export the scan as CSV.
     */
    private final Button exportButton;
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
                VigilanceScanView.createColumn("VigilanceScan.vigilanceColumn.header", anEntry -> String.format(Locale.getDefault(), "%.4f", anEntry.vigilance())),
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
        Label tmpHintLabel = new Label(Message.get("VigilanceScan.hint.text"));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.setTop(tmpHintLabel);
        //buttons for the selected vigilance
        this.openResultButton = new Button(Message.get("VigilanceScan.openResultButton.text"));
        this.openResultButton.setTooltip(new Tooltip(Message.get("VigilanceScan.openResultButton.tooltip")));
        this.openResultButton.disableProperty().bind(this.scanTableView.getSelectionModel().selectedItemProperty().isNull());
        this.useVigilanceButton = new Button(Message.get("VigilanceScan.useVigilanceButton.text"));
        this.useVigilanceButton.setTooltip(new Tooltip(Message.get("VigilanceScan.useVigilanceButton.tooltip")));
        this.useVigilanceButton.disableProperty().bind(this.scanTableView.getSelectionModel().selectedItemProperty().isNull());
        this.exportButton = new Button(Message.get("VigilanceScan.exportButton.text"));
        HBox tmpButtonBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, this.openResultButton, this.useVigilanceButton, this.exportButton);
        tmpButtonBox.setAlignment(Pos.CENTER_RIGHT);
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
     * Returns the button that opens the clustering of the selected vigilance as result tab.
     *
     * @return open result button
     */
    public Button getOpenResultButton() {
        return this.openResultButton;
    }
    //
    /**
     * Returns the button that takes over the selected vigilance into the clustering settings.
     *
     * @return use vigilance button
     */
    public Button getUseVigilanceButton() {
        return this.useVigilanceButton;
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
