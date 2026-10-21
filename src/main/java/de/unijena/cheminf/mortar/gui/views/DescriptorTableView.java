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
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.FeatureType;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;

import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;

import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.stream.IntStream;

/**
 * View of a tab with calculated descriptor values: one row per molecule with name, structure and one sortable column
 * per descriptor component (untransformed values). Calculated fingerprint bits are not shown as columns (only their
 * number), but are part of the matrix, the export and a clustering started from the tab. The columns have fixed widths, so that the table scrolls
 * horizontally with many descriptors. At the bottom, the clustering controls set by the controller (e.g. to cluster
 * these molecules with these descriptors) are shown on the left and the export button on the right; their actions are
 * set by the controller.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class DescriptorTableView extends BorderPane {
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Molecules in the order of the matrix rows.
     */
    private final List<MoleculeDataModel> moleculeList;
    /**
     * Descriptor values, one row per molecule.
     */
    private final FeatureMatrix descriptorMatrix;
    /**
     * Table with one row per molecule, the items are the row indices of the matrix.
     */
    private final TableView<Integer> tableView;
    /**
     * Button to export the table as CSV.
     */
    private final Button exportButton;
    /**
     * Bottom bar with the clustering controls on the left and the export button on the right.
     */
    private final BorderPane bottomPane;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, builds the table for the given molecules and descriptor values.
     *
     * @param aMoleculeList molecules in the order of the matrix rows
     * @param aDescriptorMatrix descriptor values, one row per molecule
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the number of molecules does not match the number of rows
     */
    public DescriptorTableView(List<MoleculeDataModel> aMoleculeList, FeatureMatrix aDescriptorMatrix)
            throws NullPointerException, IllegalArgumentException {
        super();
        this.moleculeList = List.copyOf(Objects.requireNonNull(aMoleculeList, "Given molecule list is null."));
        this.descriptorMatrix = Objects.requireNonNull(aDescriptorMatrix, "Given descriptor matrix is null.");
        if (this.moleculeList.size() != aDescriptorMatrix.getRowCount()) {
            throw new IllegalArgumentException("Number of molecules does not match the number of rows.");
        }
        //fingerprint bits (BIT columns) are calculated and exported, but not shown as hundreds of columns
        int tmpBitCount = 0;
        for (int j = 0; j < aDescriptorMatrix.getColumnCount(); j++) {
            if (aDescriptorMatrix.getFeatureType(j) == FeatureType.BIT) {
                tmpBitCount++;
            }
        }
        Label tmpSummaryLabel = new Label(String.format(Message.get("DescriptorTableView.summary.text"),
                this.moleculeList.size(), aDescriptorMatrix.getColumnCount() - tmpBitCount, tmpBitCount));
        tmpSummaryLabel.setWrapText(true);
        tmpSummaryLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.setTop(tmpSummaryLabel);
        this.tableView = new TableView<>(FXCollections.observableArrayList(
                IntStream.range(0, this.moleculeList.size()).boxed().toList()));
        //fixed column widths, with many descriptors the table scrolls horizontally
        this.tableView.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
        TableColumn<Integer, String> tmpNameColumn = new TableColumn<>(Message.get("ClusteringResultView.descriptorValues.nameColumn.header"));
        tmpNameColumn.setCellValueFactory(aCellData -> new ReadOnlyStringWrapper(this.moleculeList.get(aCellData.getValue()).getName()));
        tmpNameColumn.setPrefWidth(150.0);
        TableColumn<Integer, MoleculeDataModel> tmpStructureColumn = new TableColumn<>(
                Message.get("ClusteringResultView.descriptorValues.structureColumn.header"));
        tmpStructureColumn.setCellValueFactory(aCellData -> new ReadOnlyObjectWrapper<>(this.moleculeList.get(aCellData.getValue())));
        //the structure is depicted by the cell itself when it becomes visible, centred and fitted to the column
        tmpStructureColumn.setCellFactory(aColumn -> ClusteringResultView.createStructureCell());
        tmpStructureColumn.setSortable(false);
        tmpStructureColumn.setPrefWidth(ClusteringResultView.STRUCTURE_COLUMN_WIDTH);
        this.tableView.getColumns().add(tmpNameColumn);
        this.tableView.getColumns().add(tmpStructureColumn);
        for (int j = 0; j < aDescriptorMatrix.getColumnCount(); j++) {
            if (aDescriptorMatrix.getFeatureType(j) != FeatureType.BIT) {
                this.tableView.getColumns().add(this.createDescriptorColumn(j));
            }
        }
        this.setCenter(this.tableView);
        this.exportButton = new Button(Message.get("DescriptorTableView.exportButton.text"));
        this.exportButton.setTooltip(new Tooltip(Message.get("DescriptorTableView.exportButton.tooltip")));
        this.bottomPane = new BorderPane();
        this.bottomPane.setRight(this.exportButton);
        BorderPane.setAlignment(this.exportButton, Pos.CENTER_RIGHT);
        this.bottomPane.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.setBottom(this.bottomPane);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the molecules in the order of the matrix rows.
     *
     * @return unmodifiable molecule list
     */
    public List<MoleculeDataModel> getMoleculeList() {
        return this.moleculeList;
    }
    //
    /**
     * Returns the descriptor values.
     *
     * @return descriptor matrix, one row per molecule
     */
    public FeatureMatrix getDescriptorMatrix() {
        return this.descriptorMatrix;
    }
    //
    /**
     * Returns the button that exports the table as CSV.
     *
     * @return export button
     */
    public Button getExportButton() {
        return this.exportButton;
    }
    //
    /**
     * Sets the clustering controls shown on the left of the bottom bar, e.g. vigilance and a button to cluster these
     * molecules with these descriptors.
     *
     * @param aClusteringControls clustering controls, null to remove them
     */
    public void setClusteringControls(HBox aClusteringControls) {
        this.bottomPane.setLeft(aClusteringControls);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Creates the numeric column of the given matrix column, with a fixed width fitting its header; NaN (calculation
     * failed) is shown as "-" and sorted like a missing value.
     *
     * @param aColumnIndex matrix column index
     * @return table column
     */
    private TableColumn<Integer, Double> createDescriptorColumn(int aColumnIndex) {
        String tmpName = this.descriptorMatrix.getFeatureName(aColumnIndex);
        TableColumn<Integer, Double> tmpColumn = new TableColumn<>(tmpName);
        ClusteringResultView.setFixedWidthForHeader(tmpColumn, tmpName);
        tmpColumn.setCellValueFactory(aCellData -> {
            float tmpValue = this.descriptorMatrix.getValue(aCellData.getValue(), aColumnIndex);
            return new ReadOnlyObjectWrapper<>(Float.isNaN(tmpValue) ? null : (double) tmpValue);
        });
        tmpColumn.setCellFactory(aColumn -> new TableCell<>() {
            @Override
            protected void updateItem(Double aValue, boolean anIsEmpty) {
                super.updateItem(aValue, anIsEmpty);
                this.setText(anIsEmpty ? null : aValue == null ? "-" : String.format(Locale.getDefault(), "%.2f", aValue));
            }
        });
        return tmpColumn;
    }
    //</editor-fold>
}
