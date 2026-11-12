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

import de.unijena.cheminf.mortar.configuration.Configuration;
import de.unijena.cheminf.mortar.configuration.IConfiguration;
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;
import de.unijena.cheminf.mortar.model.util.SimpleSmartsPatternListProperty;

import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.StringProperty;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonBar;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Modality;
import javafx.stage.Stage;
import javafx.stage.Window;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Modal dialog for editing a SMARTS pattern list (see {@link SimpleSmartsPatternListProperty}), e.g. the patterns of
 * the custom SMARTS fragmenter. Every pattern is shown in a row with a field for the SMARTS pattern, a field for its
 * optional name, and buttons to move or remove it. Invalid SMARTS patterns are highlighted while typing, with the
 * message of the CDK SMARTS parser as tooltip, and the changes can only be applied while all non-empty rows are valid.
 * Patterns can be appended from a text file in the format of the property, and the patterns can be exported in this
 * format. The dialog returns the new pattern list text if the changes are applied and no result if they are cancelled;
 * empty rows are dropped and comments are not kept.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class SmartsPatternListEditorView extends Dialog<String> {
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(SmartsPatternListEditorView.class.getName());

    /**
     * Style of a text field with an invalid entry.
     */
    private static final String INVALID_FIELD_STYLE = "-fx-text-box-border: red; -fx-focus-color: red; -fx-control-inner-background: #ffecec;";

    /**
     * Preferred width of the SMARTS fields.
     */
    private static final double SMARTS_FIELD_PREF_WIDTH = 380.0;

    /**
     * Preferred width of the name fields.
     */
    private static final double NAME_FIELD_PREF_WIDTH = 220.0;

    /**
     * Width of the buttons of a row.
     */
    private static final double ROW_BUTTON_WIDTH = 32.0;

    /**
     * Preferred width of the dialog content.
     */
    private static final double CONTENT_PREF_WIDTH = 820.0;

    /**
     * Preferred height of the dialog content.
     */
    private static final double CONTENT_PREF_HEIGHT = 480.0;

    /**
     * File name patterns of the text files offered for import and export.
     */
    private static final String[] TEXT_FILE_PATTERNS = new String[] {"*.txt", "*.smarts", "*.sma"};
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Container of the pattern rows, in pattern order.
     */
    private final VBox rowsBox;

    /**
     * Whether at least one non-empty row is invalid.
     */
    private final SimpleBooleanProperty hasInvalidRows;

    /**
     * Recent directory path setting of the application used for the file choosers, may be null.
     */
    private final StringProperty recentDirectoryPathProperty;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, creates the dialog with one row per pattern of the given text, or one empty row if the text holds no
     * patterns. Call {@link #showAndWait()} to show it.
     *
     * @param anOwner owner window of the dialog, may be null
     * @param aPatternListText the current pattern list text
     * @param aRecentDirectoryPathProperty recent directory path setting that the file choosers start in and update; may
     *                                     be null, then the user home directory is used
     * @throws NullPointerException if the given pattern list text is null
     */
    public SmartsPatternListEditorView(Window anOwner, String aPatternListText, StringProperty aRecentDirectoryPathProperty)
            throws NullPointerException {
        super();
        Objects.requireNonNull(aPatternListText, "Given pattern list text is null.");
        this.recentDirectoryPathProperty = aRecentDirectoryPathProperty;
        this.hasInvalidRows = new SimpleBooleanProperty(false);
        if (!Objects.isNull(anOwner)) {
            this.initOwner(anOwner);
        }
        this.initModality(Modality.WINDOW_MODAL);
        this.setResizable(true);
        this.setTitle(Message.get("SmartsPatternListEditorView.title"));
        this.setHeaderText(Message.get("SmartsPatternListEditorView.headerText"));
        DialogPane tmpDialogPane = this.getDialogPane();
        //rows
        this.rowsBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE);
        this.rowsBox.setPadding(new Insets(GuiDefinitions.GUI_SPACING_VALUE));
        ScrollPane tmpScrollPane = new ScrollPane(this.rowsBox);
        tmpScrollPane.setFitToWidth(true);
        tmpScrollPane.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);
        //column headers, aligned with the fields of the rows
        Label tmpSmartsHeader = new Label(Message.get("SmartsPatternListEditorView.smartsColumn.text"));
        tmpSmartsHeader.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(tmpSmartsHeader, Priority.ALWAYS);
        Label tmpNameHeader = new Label(Message.get("SmartsPatternListEditorView.nameColumn.text"));
        tmpNameHeader.setPrefWidth(SmartsPatternListEditorView.NAME_FIELD_PREF_WIDTH);
        tmpNameHeader.setMinWidth(SmartsPatternListEditorView.NAME_FIELD_PREF_WIDTH);
        Region tmpButtonsSpacer = new Region();
        tmpButtonsSpacer.setMinWidth(3 * SmartsPatternListEditorView.ROW_BUTTON_WIDTH + 3 * GuiDefinitions.GUI_SPACING_VALUE);
        HBox tmpHeaderBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpSmartsHeader, tmpNameHeader, tmpButtonsSpacer);
        //spacing of the rows box and scroll pane border
        tmpHeaderBox.setPadding(new Insets(0, GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_SPACING_VALUE,
                GuiDefinitions.GUI_SPACING_VALUE + 1));
        //tool buttons
        Button tmpAddRowButton = new Button(Message.get("SmartsPatternListEditorView.addRowButton.text"));
        tmpAddRowButton.setTooltip(GuiUtil.createTooltip(Message.get("SmartsPatternListEditorView.addRowButton.tooltip")));
        tmpAddRowButton.setStyle("-fx-font-weight: bold");
        tmpAddRowButton.setMinWidth(GuiDefinitions.GUI_PIPELINE_SETTINGS_VIEW_BUTTON_WIDTH_VALUE);
        tmpAddRowButton.setOnAction(anEvent -> {
            PatternRow tmpRow = this.addRow("", "");
            this.updateValidity();
            tmpScrollPane.layout();
            tmpScrollPane.setVvalue(tmpScrollPane.getVmax());
            tmpRow.smartsField.requestFocus();
        });
        Button tmpImportButton = new Button(Message.get("SmartsPatternListEditorView.importButton.text"));
        tmpImportButton.setTooltip(GuiUtil.createTooltip(Message.get("SmartsPatternListEditorView.importButton.tooltip")));
        tmpImportButton.setOnAction(anEvent -> this.importPatterns());
        Button tmpExportButton = new Button(Message.get("SmartsPatternListEditorView.exportButton.text"));
        tmpExportButton.setTooltip(GuiUtil.createTooltip(Message.get("SmartsPatternListEditorView.exportButton.tooltip")));
        tmpExportButton.setOnAction(anEvent -> this.exportPatterns());
        tmpExportButton.disableProperty().bind(this.hasInvalidRows);
        Button tmpRemoveAllButton = new Button(Message.get("SmartsPatternListEditorView.removeAllButton.text"));
        tmpRemoveAllButton.setTooltip(GuiUtil.createTooltip(Message.get("SmartsPatternListEditorView.removeAllButton.tooltip")));
        tmpRemoveAllButton.setOnAction(anEvent -> {
            this.rowsBox.getChildren().clear();
            this.updateValidity();
        });
        HBox tmpToolBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpAddRowButton, tmpImportButton, tmpExportButton,
                tmpRemoveAllButton);
        tmpToolBox.setAlignment(Pos.CENTER_LEFT);
        Label tmpInvalidRowsHint = new Label(Message.get("SmartsPatternListEditorView.invalidRowsHint.text"));
        tmpInvalidRowsHint.setStyle("-fx-text-fill: red;");
        tmpInvalidRowsHint.visibleProperty().bind(this.hasInvalidRows);
        tmpInvalidRowsHint.managedProperty().bind(this.hasInvalidRows);
        VBox tmpBottomBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpToolBox, tmpInvalidRowsHint);
        tmpBottomBox.setPadding(new Insets(GuiDefinitions.GUI_SPACING_VALUE, 0, 0, 0));
        BorderPane tmpContentPane = new BorderPane(tmpScrollPane, tmpHeaderBox, null, tmpBottomBox, null);
        tmpContentPane.setPrefSize(SmartsPatternListEditorView.CONTENT_PREF_WIDTH, SmartsPatternListEditorView.CONTENT_PREF_HEIGHT);
        tmpDialogPane.setContent(tmpContentPane);
        //apply and cancel buttons
        ButtonType tmpApplyButtonType = new ButtonType(Message.get("SmartsPatternListEditorView.applyButton.text"),
                ButtonBar.ButtonData.OK_DONE);
        ButtonType tmpCancelButtonType = new ButtonType(Message.get("SmartsPatternListEditorView.cancelButton.text"),
                ButtonBar.ButtonData.CANCEL_CLOSE);
        tmpDialogPane.getButtonTypes().addAll(tmpApplyButtonType, tmpCancelButtonType);
        tmpDialogPane.lookupButton(tmpApplyButtonType).disableProperty().bind(this.hasInvalidRows);
        this.setResultConverter(aButtonType -> aButtonType == tmpApplyButtonType ? this.getPatternListText() : null);
        SmartsPatternListEditorView.addIcon((Stage) tmpDialogPane.getScene().getWindow());
        //initial rows
        for (SimpleSmartsPatternListProperty.SmartsPatternEntry tmpEntry : SimpleSmartsPatternListProperty.parseEntries(aPatternListText)) {
            this.addRow(tmpEntry.smarts(), tmpEntry.name());
        }
        if (this.rowsBox.getChildren().isEmpty()) {
            this.addRow("", "");
        }
        this.updateValidity();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Appends a row with the given SMARTS pattern and name.
     *
     * @param aSmarts SMARTS pattern, may be empty
     * @param aName name, may be empty
     * @return the new row
     */
    private PatternRow addRow(String aSmarts, String aName) {
        PatternRow tmpRow = new PatternRow(aSmarts, aName);
        this.rowsBox.getChildren().add(tmpRow);
        return tmpRow;
    }

    /**
     * Moves the given row up or down by the given offset, if possible.
     *
     * @param aRow the row to move
     * @param anOffset -1 to move it up, 1 to move it down
     */
    private void moveRow(PatternRow aRow, int anOffset) {
        List<Node> tmpRows = this.rowsBox.getChildren();
        int tmpIndex = tmpRows.indexOf(aRow);
        int tmpNewIndex = tmpIndex + anOffset;
        if (tmpIndex < 0 || tmpNewIndex < 0 || tmpNewIndex >= tmpRows.size()) {
            return;
        }
        tmpRows.remove(tmpIndex);
        tmpRows.add(tmpNewIndex, aRow);
    }

    /**
     * Returns the rows in their order.
     *
     * @return rows
     */
    private List<PatternRow> getRows() {
        List<PatternRow> tmpRows = new ArrayList<>(this.rowsBox.getChildren().size());
        for (Node tmpNode : this.rowsBox.getChildren()) {
            tmpRows.add((PatternRow) tmpNode);
        }
        return tmpRows;
    }

    /**
     * Updates whether at least one row is invalid; empty rows are valid.
     */
    private void updateValidity() {
        boolean tmpHasInvalidRows = false;
        for (PatternRow tmpRow : this.getRows()) {
            if (!tmpRow.isValid()) {
                tmpHasInvalidRows = true;
                break;
            }
        }
        this.hasInvalidRows.set(tmpHasInvalidRows);
    }

    /**
     * Returns the pattern list text of all non-empty rows, which must be valid.
     *
     * @return pattern list text
     */
    private String getPatternListText() {
        List<SimpleSmartsPatternListProperty.SmartsPatternEntry> tmpEntries = new ArrayList<>(this.rowsBox.getChildren().size());
        for (PatternRow tmpRow : this.getRows()) {
            if (!tmpRow.isEmpty()) {
                tmpEntries.add(new SimpleSmartsPatternListProperty.SmartsPatternEntry(tmpRow.getSmarts(), tmpRow.getName()));
            }
        }
        return SimpleSmartsPatternListProperty.toText(tmpEntries);
    }

    /**
     * Lets the user choose a text file and appends its valid patterns that are not in the list yet. Afterwards, a
     * summary lists every skipped line with its line number and the reason.
     */
    private void importPatterns() {
        FileChooser tmpFileChooser = this.createFileChooser(Message.get("SmartsPatternListEditorView.fileChooser.importTitle"));
        File tmpFile = tmpFileChooser.showOpenDialog(this.getDialogPane().getScene().getWindow());
        if (Objects.isNull(tmpFile)) {
            return;
        }
        this.updateRecentDirectory(tmpFile);
        String tmpText;
        try {
            tmpText = SmartsPatternListEditorView.readTextFile(tmpFile);
        } catch (IOException | SecurityException anException) {
            SmartsPatternListEditorView.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                    Message.get("SmartsPatternListEditorView.importError.header"),
                    anException.toString(),
                    anException);
            return;
        }
        //empty rows would end up between the existing and the imported patterns
        this.rowsBox.getChildren().removeIf(aNode -> ((PatternRow) aNode).isEmpty());
        Set<String> tmpSmartsInList = new HashSet<>();
        for (PatternRow tmpRow : this.getRows()) {
            tmpSmartsInList.add(tmpRow.getSmarts());
        }
        List<SimpleSmartsPatternListProperty.ParsedLine> tmpSkippedLines = new ArrayList<>(0);
        int tmpImportedCount = 0;
        for (SimpleSmartsPatternListProperty.ParsedLine tmpLine : SimpleSmartsPatternListProperty.parseLines(tmpText)) {
            if (!tmpLine.isValid()) {
                tmpSkippedLines.add(tmpLine);
            } else if (!tmpSmartsInList.add(tmpLine.entry().smarts())) {
                tmpSkippedLines.add(new SimpleSmartsPatternListProperty.ParsedLine(tmpLine.lineNumber(), tmpLine.line(),
                        null, Message.get("SmartsPatternListEditorView.importSummary.duplicate")));
            } else {
                this.addRow(tmpLine.entry().smarts(), tmpLine.entry().name());
                tmpImportedCount++;
            }
        }
        if (this.rowsBox.getChildren().isEmpty()) {
            this.addRow("", "");
        }
        this.updateValidity();
        SmartsPatternListEditorView.LOGGER.log(Level.INFO, "Imported {0} SMARTS patterns from {1}, {2} lines skipped.",
                new Object[] {tmpImportedCount, tmpFile.getAbsolutePath(), tmpSkippedLines.size()});
        this.showImportSummary(tmpFile.getName(), tmpImportedCount, tmpSkippedLines);
    }

    /**
     * Shows a summary of an import with all skipped lines in an expanded details area.
     *
     * @param aFileName name of the imported file
     * @param anImportedCount number of imported patterns
     * @param aSkippedLineList skipped lines with the reasons
     */
    private void showImportSummary(String aFileName, int anImportedCount, List<SimpleSmartsPatternListProperty.ParsedLine> aSkippedLineList) {
        Alert tmpAlert = new Alert(aSkippedLineList.isEmpty() ? Alert.AlertType.INFORMATION : Alert.AlertType.WARNING);
        tmpAlert.initOwner(this.getDialogPane().getScene().getWindow());
        tmpAlert.setTitle(Message.get("SmartsPatternListEditorView.importSummary.title"));
        tmpAlert.setHeaderText(String.format(Message.get("SmartsPatternListEditorView.importSummary.header"),
                anImportedCount, aFileName));
        if (aSkippedLineList.isEmpty()) {
            tmpAlert.setContentText(Message.get("SmartsPatternListEditorView.importSummary.allImported"));
        } else {
            tmpAlert.setContentText(String.format(Message.get("SmartsPatternListEditorView.importSummary.linesSkipped"),
                    aSkippedLineList.size()));
            TextArea tmpTextArea = new TextArea(SimpleSmartsPatternListProperty.describeInvalidLines(aSkippedLineList));
            tmpTextArea.setEditable(false);
            tmpTextArea.setWrapText(false);
            tmpTextArea.setMaxWidth(Double.MAX_VALUE);
            tmpTextArea.setMaxHeight(Double.MAX_VALUE);
            VBox.setVgrow(tmpTextArea, Priority.ALWAYS);
            VBox tmpDetailsBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE,
                    new Label(Message.get("SmartsPatternListEditorView.importSummary.detailsLabel")), tmpTextArea);
            tmpDetailsBox.setMaxWidth(Double.MAX_VALUE);
            tmpAlert.getDialogPane().setExpandableContent(tmpDetailsBox);
            tmpAlert.getDialogPane().setExpanded(true);
            tmpAlert.setResizable(true);
        }
        tmpAlert.getDialogPane().setMinHeight(Region.USE_PREF_SIZE);
        tmpAlert.getDialogPane().setMinWidth(Region.USE_PREF_SIZE);
        SmartsPatternListEditorView.addIcon((Stage) tmpAlert.getDialogPane().getScene().getWindow());
        tmpAlert.showAndWait();
    }

    /**
     * Lets the user choose a file and writes the patterns into it, preceded by a comment line that describes the format.
     */
    private void exportPatterns() {
        FileChooser tmpFileChooser = this.createFileChooser(Message.get("SmartsPatternListEditorView.fileChooser.exportTitle"));
        tmpFileChooser.setInitialFileName(Message.get("SmartsPatternListEditorView.exportFile.defaultName"));
        File tmpFile = tmpFileChooser.showSaveDialog(this.getDialogPane().getScene().getWindow());
        if (Objects.isNull(tmpFile)) {
            return;
        }
        this.updateRecentDirectory(tmpFile);
        String tmpContent = SimpleSmartsPatternListProperty.COMMENT_LINE_PREFIX + " "
                + Message.get("SmartsPatternListEditorView.exportFile.header")
                + SimpleSmartsPatternListProperty.LINE_SEPARATOR
                + this.getPatternListText()
                + SimpleSmartsPatternListProperty.LINE_SEPARATOR;
        try {
            Files.writeString(tmpFile.toPath(), tmpContent, StandardCharsets.UTF_8);
        } catch (IOException | SecurityException anException) {
            SmartsPatternListEditorView.LOGGER.log(Level.WARNING, anException.toString(), anException);
            GuiUtil.guiExceptionAlert(Message.get("Error.ExceptionAlert.Title"),
                    Message.get("SmartsPatternListEditorView.exportError.header"),
                    anException.toString(),
                    anException);
        }
    }

    /**
     * Creates a file chooser for text files that starts in the recent directory.
     *
     * @param aTitle title of the file chooser
     * @return file chooser
     */
    private FileChooser createFileChooser(String aTitle) {
        FileChooser tmpFileChooser = new FileChooser();
        tmpFileChooser.setTitle(aTitle);
        tmpFileChooser.getExtensionFilters().addAll(
                new FileChooser.ExtensionFilter(Message.get("SmartsPatternListEditorView.fileChooser.textFilesDescription"),
                        SmartsPatternListEditorView.TEXT_FILE_PATTERNS),
                new FileChooser.ExtensionFilter(Message.get("SmartsPatternListEditorView.fileChooser.allFilesDescription"),
                        "*.*"));
        File tmpDirectory = null;
        if (!Objects.isNull(this.recentDirectoryPathProperty) && !Objects.isNull(this.recentDirectoryPathProperty.get())) {
            tmpDirectory = new File(this.recentDirectoryPathProperty.get());
        }
        if (Objects.isNull(tmpDirectory) || !tmpDirectory.isDirectory()) {
            tmpDirectory = new File(SettingsContainer.RECENT_DIRECTORY_PATH_SETTING_DEFAULT);
        }
        if (tmpDirectory.isDirectory()) {
            tmpFileChooser.setInitialDirectory(tmpDirectory);
        }
        return tmpFileChooser;
    }

    /**
     * Sets the recent directory path setting to the directory of the given file, if the setting is available and
     * accepts the path.
     *
     * @param aFile the chosen file
     */
    private void updateRecentDirectory(File aFile) {
        if (Objects.isNull(this.recentDirectoryPathProperty) || Objects.isNull(aFile.getParent())) {
            return;
        }
        try {
            this.recentDirectoryPathProperty.set(aFile.getParent() + File.separator);
        } catch (IllegalArgumentException anException) {
            //the setting rejects some paths, e.g. with characters that cannot be persisted; it is logged by the setting
            SmartsPatternListEditorView.LOGGER.log(Level.FINE, anException.toString(), anException);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static methods">
    /**
     * Reads the given text file as UTF-8 or, if it is no valid UTF-8 text, as ISO-8859-1.
     *
     * @param aFile the file
     * @return the content of the file
     * @throws IOException if the file cannot be read
     */
    private static String readTextFile(File aFile) throws IOException {
        byte[] tmpBytes = Files.readAllBytes(aFile.toPath());
        try {
            return StandardCharsets.UTF_8.newDecoder().decode(ByteBuffer.wrap(tmpBytes)).toString();
        } catch (CharacterCodingException anException) {
            SmartsPatternListEditorView.LOGGER.log(Level.INFO, "File {0} is no valid UTF-8 text and is read as ISO-8859-1.",
                    aFile.getAbsolutePath());
            return new String(tmpBytes, StandardCharsets.ISO_8859_1);
        }
    }

    /**
     * Adds the MORTAR icon to the given dialog stage, as done for the alerts of the application.
     *
     * @param aStage the stage of a dialog
     */
    private static void addIcon(Stage aStage) {
        try {
            IConfiguration tmpConfiguration = Configuration.getInstance();
            String tmpIconURL = SmartsPatternListEditorView.class.getClassLoader().getResource(
                    tmpConfiguration.getProperty("mortar.imagesFolder")
                            + tmpConfiguration.getProperty("mortar.logo.icon.name")).toExternalForm();
            aStage.getIcons().add(new Image(tmpIconURL));
        } catch (IOException | NullPointerException anException) {
            SmartsPatternListEditorView.LOGGER.log(Level.WARNING, anException.toString(), anException);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Class PatternRow">
    /**
     * Row of the dialog for one pattern, with fields for the SMARTS pattern and the name and buttons to move the row
     * up or down and to remove it. The row validates its fields whenever they change.
     */
    private final class PatternRow extends HBox {
        /**
         * Field for the SMARTS pattern.
         */
        private final TextField smartsField;
        /**
         * Field for the optional name.
         */
        private final TextField nameField;
        /**
         * Description of the problem of the row, null if the row is valid or empty.
         */
        private String errorMessage;
        /**
         * Constructor.
         *
         * @param aSmarts SMARTS pattern, may be empty
         * @param aName name, may be empty
         */
        private PatternRow(String aSmarts, String aName) {
            super(GuiDefinitions.GUI_SPACING_VALUE);
            this.setAlignment(Pos.CENTER_LEFT);
            this.smartsField = new TextField(aSmarts);
            this.smartsField.setPromptText(Message.get("SmartsPatternListEditorView.smartsField.promptText"));
            this.smartsField.setPrefWidth(SmartsPatternListEditorView.SMARTS_FIELD_PREF_WIDTH);
            HBox.setHgrow(this.smartsField, Priority.ALWAYS);
            this.nameField = new TextField(aName);
            this.nameField.setPromptText(Message.get("SmartsPatternListEditorView.nameField.promptText"));
            this.nameField.setPrefWidth(SmartsPatternListEditorView.NAME_FIELD_PREF_WIDTH);
            this.nameField.setMinWidth(SmartsPatternListEditorView.NAME_FIELD_PREF_WIDTH);
            Button tmpMoveUpButton = this.createButton("SmartsPatternListEditorView.moveUpButton");
            tmpMoveUpButton.setOnAction(anEvent -> SmartsPatternListEditorView.this.moveRow(this, -1));
            Button tmpMoveDownButton = this.createButton("SmartsPatternListEditorView.moveDownButton");
            tmpMoveDownButton.setOnAction(anEvent -> SmartsPatternListEditorView.this.moveRow(this, 1));
            Button tmpRemoveButton = this.createButton("SmartsPatternListEditorView.removeRowButton");
            tmpRemoveButton.setStyle("-fx-font-weight: bold");
            tmpRemoveButton.setOnAction(anEvent -> {
                SmartsPatternListEditorView.this.rowsBox.getChildren().remove(this);
                SmartsPatternListEditorView.this.updateValidity();
            });
            this.getChildren().addAll(this.smartsField, this.nameField, tmpMoveUpButton, tmpMoveDownButton, tmpRemoveButton);
            this.smartsField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> {
                this.validate();
                SmartsPatternListEditorView.this.updateValidity();
            });
            this.nameField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> {
                this.validate();
                SmartsPatternListEditorView.this.updateValidity();
            });
            this.validate();
        }
        //
        /**
         * Returns the SMARTS pattern without leading and trailing whitespace.
         *
         * @return SMARTS pattern, may be empty
         */
        private String getSmarts() {
            return Objects.requireNonNullElse(this.smartsField.getText(), "").strip();
        }
        //
        /**
         * Returns the name without leading and trailing whitespace.
         *
         * @return name, may be empty
         */
        private String getName() {
            return Objects.requireNonNullElse(this.nameField.getText(), "").strip();
        }
        //
        /**
         * Returns whether both fields are empty.
         *
         * @return true if the row is empty
         */
        private boolean isEmpty() {
            return this.getSmarts().isEmpty() && this.getName().isEmpty();
        }
        //
        /**
         * Returns whether the row is valid or empty.
         *
         * @return true if the row has no problem
         */
        private boolean isValid() {
            return Objects.isNull(this.errorMessage);
        }
        //
        /**
         * Validates the fields and marks invalid ones with a style and a tooltip that describes the problem.
         */
        private void validate() {
            String tmpSmartsError = this.isEmpty() ? null : SimpleSmartsPatternListProperty.validateSmarts(this.getSmarts());
            String tmpNameError = SimpleSmartsPatternListProperty.validateName(Objects.requireNonNullElse(this.nameField.getText(), ""));
            this.smartsField.setStyle(Objects.isNull(tmpSmartsError) ? "" : SmartsPatternListEditorView.INVALID_FIELD_STYLE);
            this.smartsField.setTooltip(Objects.isNull(tmpSmartsError) ? null : GuiUtil.createTooltip(tmpSmartsError));
            this.nameField.setStyle(Objects.isNull(tmpNameError) ? "" : SmartsPatternListEditorView.INVALID_FIELD_STYLE);
            this.nameField.setTooltip(GuiUtil.createTooltip(Objects.isNull(tmpNameError)
                    ? Message.get("SmartsPatternListEditorView.nameField.tooltip") : tmpNameError));
            this.errorMessage = Objects.isNull(tmpSmartsError) ? tmpNameError : tmpSmartsError;
        }
        //
        /**
         * Creates a small button of the row with the text and tooltip of the given message key prefix.
         *
         * @param aMessageKeyPrefix message key without ".text" and ".tooltip"
         * @return the button
         */
        private Button createButton(String aMessageKeyPrefix) {
            Button tmpButton = new Button(Message.get(aMessageKeyPrefix + ".text"));
            tmpButton.setTooltip(GuiUtil.createTooltip(Message.get(aMessageKeyPrefix + ".tooltip")));
            tmpButton.setMinWidth(SmartsPatternListEditorView.ROW_BUTTON_WIDTH);
            tmpButton.setPrefWidth(SmartsPatternListEditorView.ROW_BUTTON_WIDTH);
            tmpButton.setMaxWidth(SmartsPatternListEditorView.ROW_BUTTON_WIDTH);
            return tmpButton;
        }
    }
    //</editor-fold>
}
