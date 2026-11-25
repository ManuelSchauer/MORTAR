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

import de.unijena.cheminf.mortar.configuration.IConfiguration;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.data.DataModelPropertiesForTableView;
import de.unijena.cheminf.mortar.model.data.FragmentDataModel;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.depict.DepictionUtil;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import javafx.beans.binding.Bindings;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.VBox;
import javafx.scene.shape.Rectangle;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Custom table view for the itemization table view.
 *
 * @author Felix Baensch
 * @version 1.0.0.0
 *
 */
public class ItemizationDataTableView extends TableView implements IDataTableView {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Height of the header bar above fragments with a display label that marks the group of the fragment, e.g. the
     * substructures matched by the same SMARTS pattern.
     */
    private static final double GROUP_HEADER_HEIGHT = 20.0;
    /**
     * Style of the header bar marking a fragment group, similar to a table column header. The border widths are
     * appended depending on the position of the cell in the group.
     */
    private static final String GROUP_HEADER_STYLE = "-fx-background-color: linear-gradient(#f4f4f4, #e2e2e2);"
            + " -fx-border-color: #c8c8c8; -fx-font-size: 11px; -fx-padding: 0 0 0 6;"
            //the text fill is set explicitly because selected table cells would otherwise turn the label text white
            + " -fx-text-fill: black;";
    /**
     * Horizontal distance between the header bars of two neighbouring fragment groups.
     */
    private static final double GROUP_HEADER_GAP = 3.0;
    /**
     * Style of the fragment cells if fragments with display labels are shown. The cell padding is removed and only the
     * right border (the vertical grid line) is kept, so the cell content starts directly at the left cell border. The
     * content is as wide as the whole cell and covers the grid line in the header area, so that the header bar of a
     * group is continuous over all its cells while the grid lines remain visible below it.
     */
    private static final String LABELLED_FRAGMENT_CELL_STYLE = "-fx-alignment: TOP_LEFT; -fx-padding: 0;"
            + " -fx-border-width: 0 1 0 0;";
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    /**
     * TableColumn for name of the molecule.
     */
    private final TableColumn<MoleculeDataModel, String> nameColumn;
    /**
     * TableColumn for 2D structure of the molecule.
     */
    private final TableColumn<MoleculeDataModel, Image> moleculeStructureColumn;
    /**
     * TableColumn for 2D structure of the fragments of the molecule with labelling of how often
     * the respective fragment occurs in the molecule. BorderPane is used to align text and image.
     */
    private final TableColumn<MoleculeDataModel, BorderPane> fragmentStructureColumn;
    /**
     * Name of the fragmentation algorithm used.
     */
    private final String fragmentationName;
    /**
     * List which contains all items to be shown in this tableview not only the displayed ones for this page (Pagination).
     */
    private List<MoleculeDataModel> itemsList;
    /**
     * ContextMenu ot the TableView.
     */
    private final ContextMenu contextMenu;
    /**
     * MenuItem of ContextMenu to copy selected cell to clipboard.
     */
    private final MenuItem copyMenuItem;
    /**
     * MenuItem of ContextMenu to open an overview view with the item and its fragments.
     */
    private final MenuItem overviewViewMenuItem;
    /**
     * Configuration class to read resource file paths from.
     */
    private final IConfiguration configuration;
    /**
     * Cached {@link FontMetrics} instance (created once in the constructor using
     * {@link DepictionUtil#getGraphicsInstanceWithStandardFont(int, int)}) that is reused by the
     * fragment-column cell-value factories to avoid allocating a new {@link java.awt.image.BufferedImage}
     * and {@link Graphics2D} on every cell render.
     */
    private final FontMetrics cachedFontMetrics;
    //</editor-fold>
    //
    /**
     * Constructor.
     *
     * @param aFragmentationName String of fragmentation name used as title
     * @param aConfiguration configuration instance to read resource file paths from
     */
    public ItemizationDataTableView(String aFragmentationName, IConfiguration aConfiguration) {
        super();
        this.configuration = aConfiguration;
        // Build FontMetrics once; the backing Graphics2D is disposed immediately after extraction.
        // The FontMetrics object itself is self-contained and can be reused across all cell renders.
        Graphics2D tmpGraphics2D = DepictionUtil.getGraphicsInstanceWithStandardFont(1, 1);
        this.cachedFontMetrics = tmpGraphics2D.getFontMetrics();
        tmpGraphics2D.dispose();
        this.setEditable(false);
        this.fragmentationName = aFragmentationName;
        this.getSelectionModel().setCellSelectionEnabled(true);
        //-nameColumn
        this.nameColumn = new TableColumn<>(Message.get("MainTabPane.itemizationTab.tableView.nameColumn.header"));
        this.nameColumn.setMinWidth(100); //magic number
        this.nameColumn.setResizable(true);
        this.nameColumn.setEditable(false);
        this.nameColumn.setSortable(true);
        this.nameColumn.setCellValueFactory(new PropertyValueFactory<>(DataModelPropertiesForTableView.NAME.getText()));
        this.nameColumn.setCellFactory(TextFieldTableCell.forTableColumn());
        this.nameColumn.setStyle("-fx-alignment: CENTER");
        this.nameColumn.prefWidthProperty().bind(
                this.widthProperty().multiply(0.15) //magic number
        );
        //-moleculeStructureColumn
        this.moleculeStructureColumn = new TableColumn<>(Message.get("MainTabPane.itemizationTab.tableView.moleculeStructureColumn.header"));
        this.moleculeStructureColumn.setMinWidth(300); //magic number
        this.moleculeStructureColumn.setResizable(true);
        this.moleculeStructureColumn.setEditable(false);
        this.moleculeStructureColumn.setSortable(false);
        this.moleculeStructureColumn.setCellValueFactory(new PropertyValueFactory(DataModelPropertiesForTableView.STRUCTURE.getText()));
        this.moleculeStructureColumn.setStyle("-fx-alignment: CENTER");
        this.moleculeStructureColumn.prefWidthProperty().bind(
                this.widthProperty().multiply(0.3) //magic number
        );
        //-fragmentStructureColumn
        this.fragmentStructureColumn = new TableColumn<>(Message.get("MainTabPane.itemizationTab.tableView.fragmentsColumn.header"));
        this.fragmentStructureColumn.setResizable(true);
        this.fragmentStructureColumn.setEditable(false);
        this.fragmentStructureColumn.setSortable(false);
        this.fragmentStructureColumn.setStyle("-fx-alignment: CENTER");
        //
        this.getColumns().addAll(this.nameColumn, this.moleculeStructureColumn, this.fragmentStructureColumn);
        //context menu
        this.contextMenu = new ContextMenu();
        this.setContextMenu(this.contextMenu);
        //-copyMenuItem
        this.copyMenuItem = new MenuItem(Message.get("TableView.contextMenu.copyMenuItem"));
        String tmpCopyIconURL = this.getClass().getClassLoader().getResource(
                this.configuration.getProperty("mortar.imagesFolder")
                        + this.configuration.getProperty("mortar.icon.copy.name")).toExternalForm();
        this.copyMenuItem.setGraphic(new ImageView(new Image(tmpCopyIconURL)));
        this.contextMenu.getItems().add(this.copyMenuItem);
        //-separatorMenuItem
        this.contextMenu.getItems().add(new SeparatorMenuItem());
        //-overviewViewMenuItem
        this.overviewViewMenuItem = new MenuItem(Message.get("TableView.contextMenu.itemsTab.overviewViewMenuItem"));
        this.contextMenu.getItems().add(this.overviewViewMenuItem);
    }
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Creates and returns an itemization tableview page.
     *
     * @param aPageIndex integer value for the page index
     * @param aFragmentationName String for unique name of fragmentation job
     * @param aSettingsContainer SettingsContainer
     * @return Node BorderPane which holds TableView as page for Pagination
     */
    public Node createItemizationTableViewPage(int aPageIndex, String aFragmentationName, SettingsContainer aSettingsContainer) {
        int tmpRowsPerPage = aSettingsContainer.getRowsPerPageSetting();
        int fromIndex = aPageIndex * tmpRowsPerPage;
        int toIndex = Math.min(fromIndex + tmpRowsPerPage, this.itemsList.size());
        List<MoleculeDataModel> tmpList = this.itemsList.subList(fromIndex, toIndex);
        boolean tmpHasLabelledFragments = this.hasLabelledFragments(tmpList);
        //a fragment produced with several display labels in a molecule takes one column per label
        int tmpItemAmount = tmpHasLabelledFragments ? this.getLargestNumberOfFragmentCells(tmpList)
                : GuiUtil.getLargestNumberOfFragmentsForGivenMoleculeListAndFragmentationName(tmpList, aFragmentationName);
        this.resetFragmentStructureColumns(tmpItemAmount, tmpHasLabelledFragments);
        for (MoleculeDataModel tmpMoleculeDataModel : tmpList) {
            tmpMoleculeDataModel.setStructureImageWidth(this.moleculeStructureColumn.getWidth());
            //width of depictions of fragments does not need to be set, the fragment-sub-columns will adjust and have as much space as they need
        }
        //also sets heights of fragment depictions for items tab
        GuiUtil.setImageStructureHeight(this, this.getHeight(), aSettingsContainer.getRowsPerPageSetting());
        this.setItems(FXCollections.observableArrayList(tmpList));
        this.scrollTo(0);
        return new BorderPane(this);
    }
    //
    /**
     * Adds a change listener to the height property of table view which sets the height for structure images to
     * each MoleculeDataModel object of the items list and refreshes the table view.
     * If image height is too small it will be set to GuiDefinitions.GUI_STRUCTURE_IMAGE_MIN_HEIGHT (50.0).
     *
     * @param aSettingsContainer SettingsContainer
     */
    public void addTableViewHeightListener(SettingsContainer aSettingsContainer) {
        this.heightProperty().addListener((observable, oldValue, newValue) -> {
            GuiUtil.setImageStructureHeight(this, newValue.doubleValue(), aSettingsContainer.getRowsPerPageSetting());
            this.refresh();
        });
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns true if at least one fragment of the given molecules has a display label in its molecule (see
     * {@link FragmentDataModel#getDisplayLabelText(MoleculeDataModel)}).
     *
     * @param aMoleculeList molecules displayed on the current page
     * @return true if labelled fragments are present
     */
    private boolean hasLabelledFragments(List<MoleculeDataModel> aMoleculeList) {
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            if (!tmpMolecule.hasMoleculeUndergoneSpecificFragmentation(this.fragmentationName)) {
                continue;
            }
            for (FragmentDataModel tmpFragment : tmpMolecule.getFragmentsOfSpecificFragmentation(this.fragmentationName)) {
                if (tmpFragment.getDisplayLabelText(tmpMolecule) != null) {
                    return true;
                }
            }
        }
        return false;
    }
    //
    /**
     * Returns the contents of the fragment cells of the given molecule in the order of the fragment columns: one cell
     * per fragment with the frequency of the fragment in the molecule, but one cell per display label for a fragment
     * that was produced with several display labels in this molecule (see
     * {@link FragmentDataModel#getDisplayLabelFrequencies(MoleculeDataModel)}), e.g. the single-atom fragment "O"
     * matched by different Biosynfoni substructures, each with the frequency of the fragment with this label. In this
     * case, the cells are ordered by the first occurrence of their display label, so that the cells with the same
     * label are neighbours and form one group; otherwise, the order of the fragments is kept.
     *
     * @param aMolecule molecule that has undergone the fragmentation of this view
     * @return contents of the fragment cells
     */
    private List<FragmentCell> getFragmentCells(MoleculeDataModel aMolecule) {
        List<FragmentDataModel> tmpFragmentList = aMolecule.getFragmentsOfSpecificFragmentation(this.fragmentationName);
        Map<String, Integer> tmpFrequencyMap = aMolecule.getFragmentFrequencyOfSpecificFragmentation(this.fragmentationName);
        List<FragmentCell> tmpCellList = new ArrayList<>(tmpFragmentList.size());
        boolean tmpHasSplitFragments = false;
        for (FragmentDataModel tmpFragment : tmpFragmentList) {
            Map<String, Integer> tmpLabelToFrequencyMap = tmpFragment.getDisplayLabelFrequencies(aMolecule);
            if (tmpLabelToFrequencyMap == null) {
                tmpCellList.add(new FragmentCell(tmpFragment, tmpFragment.getDisplayLabelText(aMolecule),
                        tmpFrequencyMap.get(tmpFragment.getUniqueSmiles())));
            } else {
                tmpHasSplitFragments = true;
                for (Map.Entry<String, Integer> tmpEntry : tmpLabelToFrequencyMap.entrySet()) {
                    tmpCellList.add(new FragmentCell(tmpFragment, tmpEntry.getKey(), tmpEntry.getValue()));
                }
            }
        }
        if (tmpHasSplitFragments) {
            //the cells of a split fragment are placed at the position of the fragment, which can separate them from the
            // other cells with the same label; the sort is stable, so the order within a group is kept
            Map<String, Integer> tmpLabelToGroupIndexMap = new HashMap<>();
            for (FragmentCell tmpCell : tmpCellList) {
                tmpLabelToGroupIndexMap.putIfAbsent(tmpCell.displayLabel(), tmpLabelToGroupIndexMap.size());
            }
            tmpCellList.sort(Comparator.comparingInt(aCell -> tmpLabelToGroupIndexMap.get(aCell.displayLabel())));
        }
        return tmpCellList;
    }
    //
    /**
     * Returns the largest number of fragment cells of one molecule of the given list, see
     * {@link #getFragmentCells(MoleculeDataModel)}.
     *
     * @param aMoleculeList molecules displayed on the current page
     * @return largest number of fragment cells
     */
    private int getLargestNumberOfFragmentCells(List<MoleculeDataModel> aMoleculeList) {
        int tmpLargestNumber = 0;
        for (MoleculeDataModel tmpMolecule : aMoleculeList) {
            if (tmpMolecule.hasMoleculeUndergoneSpecificFragmentation(this.fragmentationName)) {
                tmpLargestNumber = Math.max(tmpLargestNumber, this.getFragmentCells(tmpMolecule).size());
            }
        }
        return tmpLargestNumber;
    }
    //
    /**
     * Deletes the old columns of the fragment structures and adds new ones for the given number.
     *
     * @param anItemAmount int value for the number of columns to add
     * @param aHasLabelledFragments whether fragments with display labels are shown; if true, the cells are styled
     *                              without padding and only with the right grid line so that the group header bars
     *                              are continuous
     */
    private void resetFragmentStructureColumns(int anItemAmount, boolean aHasLabelledFragments) {
        this.fragmentStructureColumn.getColumns().clear();
        for (int i = 0; i < anItemAmount; i++) {
            int tmpIndex = i;
            TableColumn<MoleculeDataModel, Node> tmpColumn = new TableColumn<>("Fragment " + (i + 1)); //+1 to avoid 0 in GUI
            tmpColumn.setResizable(true);
            tmpColumn.setEditable(false);
            tmpColumn.setSortable(false);
            tmpColumn.setStyle(aHasLabelledFragments ? ItemizationDataTableView.LABELLED_FRAGMENT_CELL_STYLE : "-fx-alignment: CENTER;");
            tmpColumn.setCellValueFactory(cellData -> Bindings.createObjectBinding(() -> {
                if (!cellData.getValue().hasMoleculeUndergoneSpecificFragmentation(this.fragmentationName)) {
                    return null;
                }
                List<FragmentCell> tmpCellList = this.getFragmentCells(cellData.getValue());
                if (tmpIndex >= tmpCellList.size()) {
                    return null;
                }
                FragmentCell tmpCell = tmpCellList.get(tmpIndex);
                FragmentDataModel tmpFragment = tmpCell.fragment();
                int tmpFrequency = tmpCell.frequency();
                String tmpFrequencyString = String.valueOf(tmpFrequency);
                // just a precaution; it is highly unlikely that we get that many fragments of the same type for one(!) molecule
                // cachedFontMetrics is reused here to avoid allocating a BufferedImage/Graphics2D per cell render
                if (!DepictionUtil.isTextNarrowerThanImage(tmpFragment.getStructureImageWidth(), tmpFrequencyString, this.cachedFontMetrics)) {
                    tmpFrequencyString = DepictionUtil.fitIntegerDisplayToImageWidth(tmpFragment.getStructureImageWidth(), tmpFrequency, this.cachedFontMetrics);
                }
                // fragments with display labels, e.g. the name of the substructure pattern that produced them, are
                // grouped by a header bar with the label above them; the frequency is displayed below the structure;
                // the labels are those with which the fragment was produced in this molecule
                if (tmpCell.displayLabel() != null) {
                    return this.createLabelledFragmentCell(tmpCellList, tmpIndex, tmpFrequencyString, tmpColumn);
                }
                return tmpFragment.getStructureWithText(tmpFrequencyString);
            }));
            tmpColumn.setMinWidth(300);
            this.fragmentStructureColumn.getColumns().add(tmpColumn);
        }
    }
    //
    /**
     * Creates the content of a fragment cell for a fragment with a display label: a header bar with the label above
     * the structure depiction with the frequency below it. The cells of a molecule that share a display label are
     * neighbours (see {@link #getFragmentCells(MoleculeDataModel)}). Therefore, the label is only written in the first
     * cell of a group, and the bar has side borders only at the start and end of the group.
     *
     * @param aCellList contents of the fragment cells of the molecule of the row, in the order of the fragment columns
     * @param anIndex index of the cell in the given list; its display label must not be null
     * @param aFrequencyText text displayed below the structure, i.e. the frequency of the fragment with this label in
     *                       the molecule
     * @param aColumn fragment column the cell belongs to, whose width the content takes
     * @return box with the header bar and the structure depiction
     */
    private VBox createLabelledFragmentCell(List<FragmentCell> aCellList, int anIndex, String aFrequencyText,
                                            TableColumn<MoleculeDataModel, Node> aColumn) {
        FragmentDataModel tmpFragment = aCellList.get(anIndex).fragment();
        String tmpLabelText = aCellList.get(anIndex).displayLabel();
        boolean tmpIsFirstOfGroup = anIndex == 0
                || !tmpLabelText.equals(aCellList.get(anIndex - 1).displayLabel());
        boolean tmpIsLastOfGroup = anIndex == aCellList.size() - 1
                || !tmpLabelText.equals(aCellList.get(anIndex + 1).displayLabel());
        Label tmpHeader = new Label(tmpIsFirstOfGroup ? tmpLabelText : "");
        tmpHeader.setAlignment(Pos.CENTER_LEFT);
        tmpHeader.setMaxWidth(Double.MAX_VALUE);
        tmpHeader.setMinHeight(ItemizationDataTableView.GROUP_HEADER_HEIGHT);
        tmpHeader.setPrefHeight(ItemizationDataTableView.GROUP_HEADER_HEIGHT);
        tmpHeader.setStyle(ItemizationDataTableView.GROUP_HEADER_STYLE + String.format(" -fx-border-width: 1 %d 1 %d;",
                tmpIsLastOfGroup ? 1 : 0, tmpIsFirstOfGroup ? 1 : 0));
        VBox.setMargin(tmpHeader, new Insets(0,
                tmpIsLastOfGroup ? ItemizationDataTableView.GROUP_HEADER_GAP : 0,
                0,
                tmpIsFirstOfGroup ? ItemizationDataTableView.GROUP_HEADER_GAP : 0));
        double tmpImageHeight = Math.max(tmpFragment.getStructureImageHeight() - ItemizationDataTableView.GROUP_HEADER_HEIGHT, 1.0);
        ImageView tmpImage = tmpFragment.getStructureWithText(aFrequencyText, tmpFragment.getStructureImageWidth(), tmpImageHeight);
        VBox tmpCellBox = new VBox(tmpHeader, tmpImage);
        tmpCellBox.setAlignment(Pos.TOP_CENTER);
        //the box has exactly the cell width, independent of the image width, and is anchored at the left cell
        // border, so that the header bars of neighbouring cells touch and every group starts with the same
        // distance to the left cell border; an image wider than the cell is centered and clipped
        tmpCellBox.minWidthProperty().bind(aColumn.widthProperty());
        tmpCellBox.prefWidthProperty().bind(aColumn.widthProperty());
        tmpCellBox.maxWidthProperty().bind(aColumn.widthProperty());
        Rectangle tmpClip = new Rectangle();
        tmpClip.widthProperty().bind(tmpCellBox.widthProperty());
        tmpClip.heightProperty().bind(tmpCellBox.heightProperty());
        tmpCellBox.setClip(tmpClip);
        return tmpCellBox;
    }
    //</editor-fold>
    //
    //<editor-fold desc="public properties" defaultstate="collapsed">
    /**
     * Returns TableColumn for the 2D structure of the molecule.
     *
     * @return TableColumn for 2D structure
     */
    public TableColumn<MoleculeDataModel, Image> getMoleculeStructureColumn() {
        return this.moleculeStructureColumn;
    }
    //
    /**
     * Returns the fragment displayed in the given fragment column of the row of the given molecule. The fragment
     * columns can differ from the fragment list of the molecule, since a fragment produced with several display labels
     * is shown once per label, see {@link #getFragmentCells(MoleculeDataModel)}.
     *
     * @param aMolecule molecule of the row
     * @param aFragmentColumnIndex index of the fragment column, starting at 0
     * @return fragment displayed in this cell
     * @throws IndexOutOfBoundsException if the molecule has no fragment in this column
     */
    public FragmentDataModel getFragmentOfCell(MoleculeDataModel aMolecule, int aFragmentColumnIndex) throws IndexOutOfBoundsException {
        return this.getFragmentCells(aMolecule).get(aFragmentColumnIndex).fragment();
    }
    //
    /**
     * Returns name of fragmentation.
     *
     * @return String
     */
    public String getFragmentationName() {
        return this.fragmentationName;
    }
    //
    /**
     * Returns list of items shown in TableView.
     *
     * @return List {@literal <}MoleculeDataModel {@literal >}
     */
    public List<MoleculeDataModel> getItemsList() {
        return this.itemsList;
    }
    //
    /**
     * Returns MenuItem to copy selected cell.
     *
     * @return MenuItem
     */
    public MenuItem getCopyMenuItem() {
        return this.copyMenuItem;
    }
    //
    /**
     * Returns MenuItem to open an overview view of an item and its fragments.
     *
     * @return MenuItem
     */
    public MenuItem getOverviewViewMenuItem() {
        return overviewViewMenuItem;
    }
    //
    /**
     * Sets items to list.
     *
     * @param aListOfFragments List {@literal <}MoleculeDataModel {@literal >}
     */
    public void setItemsList(List<MoleculeDataModel> aListOfFragments) {
        this.itemsList = aListOfFragments;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private record" defaultstate="collapsed">
    /**
     * Content of a fragment cell: a fragment of the molecule of the row, the display label with which it was produced
     * in this molecule (null if it has none), and how often it was produced with this label in this molecule.
     *
     * @param fragment fragment displayed in the cell
     * @param displayLabel display label of the cell or null
     * @param frequency frequency displayed below the structure
     */
    private record FragmentCell(FragmentDataModel fragment, String displayLabel, int frequency) {
    }
    //</editor-fold>
}
