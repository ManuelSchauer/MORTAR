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

import de.unijena.cheminf.mortar.gui.controls.AlignedNumberAxis;
import de.unijena.cheminf.mortar.gui.util.GuiDefinitions;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusterDescriptorStatistics;
import de.unijena.cheminf.mortar.model.clustering.ClusteringRunResult;
import de.unijena.cheminf.mortar.model.clustering.FeatureMatrix;
import de.unijena.cheminf.mortar.model.clustering.PrincipalComponentAnalysis;
import de.unijena.cheminf.mortar.model.data.MoleculeDataModel;
import de.unijena.cheminf.mortar.model.depict.DepictionUtil;

import javafx.application.ConditionalFeature;
import javafx.application.Platform;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.value.ChangeListener;
import javafx.collections.FXCollections;
import javafx.collections.transformation.FilteredList;
import javafx.geometry.Insets;
import javafx.geometry.Point2D;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.ListCell;
import javafx.scene.control.ListView;
import javafx.scene.control.Separator;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleButton;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.Tooltip;
import javafx.scene.control.cell.CheckBoxListCell;
import javafx.scene.image.ImageView;
import javafx.scene.input.MouseButton;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.paint.CycleMethod;
import javafx.scene.paint.LinearGradient;
import javafx.scene.paint.Stop;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Line;
import javafx.scene.shape.Polygon;
import javafx.scene.shape.Rectangle;
import javafx.scene.text.Text;
import javafx.util.StringConverter;

import org.openscience.cdk.exception.CDKException;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.IntConsumer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * View of an interactive principal component analysis (PCA) of a clustering result: a scatter plot of the molecules on
 * two chosen principal components, coloured by cluster, with the cluster vectors (centres) as diamonds and the
 * representatives as larger points with a black border. The scroll wheel zooms around the mouse position, dragging moves
 * the view, hovering a point shows its molecule and cluster, a click shows its structure on the right and a double click
 * calls the action set by the controller (e.g. opening the cluster). On the right, the components of both axes can be
 * chosen and clusters can be shown, hidden and highlighted (by selecting them in the list). Instead of by cluster, the
 * points can be coloured by a calculated descriptor (raw value, viridis colour scale from the 2nd to the 98th percentile,
 * grey for missing values), to see what the components mean chemically. The tick marks of both axes lie on round
 * values (see {@link AlignedNumberAxis}). What the components stand for is shown below the plot and at the axes: the
 * features with the largest weight (loading) in every component, with the direction in the raw values (sign of the
 * correlation r of the raw feature values with the component scores), e.g. "PC1: nHeavyAtom +, MW +". With "3D", the molecules are shown on the first three components in a rotatable
 * three-dimensional view ({@link Pca3dView}) with the same colours, visibility, highlighting and mouse actions; there,
 * dimmed clusters are drawn in a pale colour instead of transparent.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class PcaView extends BorderPane {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Radius of the point of a molecule.
     */
    private static final double POINT_RADIUS = 3.5;
    /**
     * Radius of the point of a cluster representative.
     */
    private static final double REPRESENTATIVE_RADIUS = 5.5;
    /**
     * Factor by which the radius grows while the point is hovered.
     */
    private static final double HOVER_SCALE = 1.8;
    /**
     * Half size of the diamond of a cluster vector.
     */
    private static final double CENTRE_SIZE = 7.0;
    /**
     * Opacity of the points.
     */
    private static final double POINT_OPACITY = 0.8;
    /**
     * Opacity of the points of clusters that are not highlighted while a cluster is highlighted.
     */
    private static final double DIMMED_OPACITY = 0.1;
    /**
     * Factor of one scroll step of the zoom.
     */
    private static final double ZOOM_FACTOR = 0.85;
    /**
     * Relative margin around the points when the view is reset.
     */
    private static final double VIEW_MARGIN = 0.05;
    /**
     * Size of the structure shown for a clicked point.
     */
    private static final double DETAIL_IMAGE_SIZE = 220.0;
    /**
     * Golden angle in degrees, used to give neighbouring cluster indices clearly different hues.
     */
    private static final double GOLDEN_ANGLE = 137.508;
    /**
     * Anchor colours of the viridis colour scale, from low to high values.
     */
    private static final Color[] VIRIDIS_COLORS = {Color.web("#440154"), Color.web("#3b528b"), Color.web("#21918c"),
            Color.web("#5ec962"), Color.web("#fde725")};
    /**
     * Lower percentile of the descriptor colour scale; smaller values get the lowest colour.
     */
    private static final double COLOUR_SCALE_LOWER_PERCENTILE = 0.02;
    /**
     * Upper percentile of the descriptor colour scale; larger values get the highest colour.
     */
    private static final double COLOUR_SCALE_UPPER_PERCENTILE = 0.98;
    /**
     * Item of the colour combo box that colours the points by cluster.
     */
    private static final String CLUSTER_COLOUR_ITEM = "";
    /**
     * Number of main features shown per component below the plot.
     */
    private static final int NUMBER_OF_MAIN_FEATURES = 5;
    /**
     * Number of main features shown in the axis labels.
     */
    private static final int NUMBER_OF_AXIS_FEATURES = 3;
    /**
     * Maximum length of a feature name in the axis labels and the overview; longer names (e.g. fragment SMILES) are cut.
     */
    private static final int MAXIMUM_FEATURE_NAME_LENGTH = 28;
    /**
     * Share of the smaller score range (largest absolute score of the two shown components) that the longest possible
     * loading vector of all features reaches; all vectors share this scale, so that their lengths can be compared.
     */
    private static final double LOADING_VECTOR_LENGTH_SHARE = 0.9;
    /**
     * Colour of the loading vectors, distinct from the cluster colours and the viridis scale.
     */
    private static final Color LOADING_VECTOR_COLOR = Color.web("#b2182b");
    /**
     * Length of the sides of the arrow head of a loading vector, in pixels.
     */
    private static final double ARROW_HEAD_LENGTH = 9.0;
    /**
     * Half the opening angle of the arrow head of a loading vector, in radians.
     */
    private static final double ARROW_HEAD_ANGLE = Math.toRadians(25.0);
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(PcaView.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Displayed clustering result.
     */
    private final ClusteringRunResult runResult;
    /**
     * PCA of the clustering result.
     */
    private final PrincipalComponentAnalysis.PcaResult pcaResult;
    /**
     * Axis of the component on the x-axis.
     */
    private final AlignedNumberAxis xAxis;
    /**
     * Axis of the component on the y-axis.
     */
    private final AlignedNumberAxis yAxis;
    /**
     * Scatter plot.
     */
    private final ScatterChart<Number, Number> chart;
    /**
     * Data point of every molecule, by row index.
     */
    private final List<XYChart.Data<Number, Number>> moleculeDataList;
    /**
     * Data point of every cluster vector, by cluster index.
     */
    private final List<XYChart.Data<Number, Number>> centreDataList;
    /**
     * Visibility of every cluster in the plot, by cluster index (-1 for the not clusterable molecules).
     */
    private final Map<Integer, BooleanProperty> visiblePropertyOfClusterMap;
    /**
     * Combo box of the component on the x-axis.
     */
    private final ComboBox<Integer> xComponentComboBox;
    /**
     * Combo box of the component on the y-axis.
     */
    private final ComboBox<Integer> yComponentComboBox;
    /**
     * Check box to show the cluster vectors.
     */
    private final CheckBox centresCheckBox;
    /**
     * List of the clusters with visibility check boxes; the selected cluster is highlighted.
     */
    private final ListView<Integer> clusterListView;
    /**
     * Label that describes the hovered point.
     */
    private final Label hoverLabel;
    /**
     * Raw descriptor values of the clustering result, used to colour the points by a descriptor.
     */
    private final ClusterDescriptorStatistics descriptorStatistics;
    /**
     * Combo box to colour the points by cluster (null item) or by a descriptor (its column name).
     */
    private final ComboBox<String> colourComboBox;
    /**
     * Legend of the descriptor colour scale; only visible while the points are coloured by a descriptor.
     */
    private final VBox colourLegendBox;
    /**
     * Lower value of the legend.
     */
    private final Label legendMinimumLabel;
    /**
     * Upper value of the legend.
     */
    private final Label legendMaximumLabel;
    /**
     * Whether every molecule is the representative of its cluster, by row index.
     */
    private final boolean[] isRepresentative;
    /**
     * Axis label of every component: name, explained variance and main features.
     */
    private final String[] axisLabels;
    /**
     * Toggle button of the three-dimensional view.
     */
    private final ToggleButton threeDimensionalToggleButton;
    /**
     * Pane in the centre holding the scatter plot or the three-dimensional view.
     */
    private final StackPane plotPane;
    /**
     * Structure of the clicked point.
     */
    private final ImageView detailImageView;
    /**
     * Description of the clicked point.
     */
    private final Label detailLabel;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Action called with the cluster index of a double-clicked point; null for none.
     */
    private IntConsumer clusterDoubleClickAction;
    /**
     * Scene coordinates and axis bounds at the start of a drag, used to move the view.
     */
    private double[] dragStart;
    /**
     * Lower and upper value of the descriptor colour scale; null while the points are coloured by cluster.
     */
    private double[] colourScaleBounds;
    /**
     * Three-dimensional view, created when it is shown for the first time; null before.
     */
    private Pca3dView pca3dView;
    /**
     * Pane over the scatter plot with the loading vectors (biplot); it does not take mouse events.
     */
    private Pane loadingVectorPane;
    /**
     * Features (columns of the processed matrix) whose loading vectors are shown, in the order they were chosen.
     */
    private final Set<Integer> loadingVectorFeatureIndices = new LinkedHashSet<>();
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, builds the plot of the given PCA.
     *
     * @param aRunResult clustering result
     * @param aPcaResult PCA of the clustering result, see {@link PrincipalComponentAnalysis#calculate(ClusteringRunResult, int,
     *                   java.util.function.BooleanSupplier)}
     * @param aVisibleClusterIndices clusters shown at the start, e.g. the marked ones; empty to show all
     * @throws NullPointerException if a parameter is null
     */
    public PcaView(ClusteringRunResult aRunResult, PrincipalComponentAnalysis.PcaResult aPcaResult, Set<Integer> aVisibleClusterIndices)
            throws NullPointerException {
        super();
        this.runResult = Objects.requireNonNull(aRunResult, "Given clustering run result is null.");
        this.pcaResult = Objects.requireNonNull(aPcaResult, "Given PCA result is null.");
        Objects.requireNonNull(aVisibleClusterIndices, "Given cluster index set is null.");
        int[] tmpLabels = aRunResult.clusteringResult().clusterLabels();
        int[] tmpRepresentatives = aRunResult.clusteringResult().representativeIndices();
        //<editor-fold desc="visibility of the clusters">
        this.visiblePropertyOfClusterMap = new HashMap<>();
        List<ClusteringRunResult.ClusterSummary> tmpRankedClusters = aRunResult.getClusterSummariesRankedBySize();
        for (ClusteringRunResult.ClusterSummary tmpCluster : tmpRankedClusters) {
            BooleanProperty tmpVisibleProperty = new SimpleBooleanProperty(aVisibleClusterIndices.isEmpty()
                    || aVisibleClusterIndices.contains(tmpCluster.clusterIndex()));
            tmpVisibleProperty.addListener((anObservable, anOldValue, aNewValue) -> this.updatePointStyles());
            this.visiblePropertyOfClusterMap.put(tmpCluster.clusterIndex(), tmpVisibleProperty);
        }
        //</editor-fold>
        //<editor-fold desc="scatter plot">
        //what the components stand for, at the axes and below the plot
        this.axisLabels = new String[aPcaResult.getNumberOfComponents()];
        VBox tmpMeaningBox = new VBox(2.0);
        Label tmpMeaningTitle = new Label(Message.get("PcaView.loadings.title"));
        tmpMeaningTitle.setStyle("-fx-font-weight: bold;");
        tmpMeaningTitle.setTooltip(new Tooltip(Message.get("PcaView.loadings.tooltip")));
        tmpMeaningBox.getChildren().add(tmpMeaningTitle);
        for (int c = 0; c < aPcaResult.getNumberOfComponents(); c++) {
            this.axisLabels[c] = PcaView.getComponentName(c, aPcaResult) + ": " + this.describeComponent(c, PcaView.NUMBER_OF_AXIS_FEATURES, false);
            Label tmpComponentLabel = new Label(PcaView.getComponentName(c, aPcaResult) + ": "
                    + this.describeComponent(c, PcaView.NUMBER_OF_MAIN_FEATURES, true));
            tmpComponentLabel.setWrapText(true);
            tmpComponentLabel.setTooltip(new Tooltip(this.describeComponentInFull(c)));
            tmpMeaningBox.getChildren().add(tmpComponentLabel);
        }
        tmpMeaningBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE, GuiDefinitions.GUI_INSETS_VALUE, 0.0, GuiDefinitions.GUI_INSETS_VALUE));
        this.xAxis = new AlignedNumberAxis();
        this.yAxis = new AlignedNumberAxis();
        this.chart = new ScatterChart<>(this.xAxis, this.yAxis);
        this.chart.setAnimated(false);
        this.chart.setLegendVisible(false);
        XYChart.Series<Number, Number> tmpMoleculeSeries = new XYChart.Series<>();
        this.moleculeDataList = new ArrayList<>(tmpLabels.length);
        this.isRepresentative = new boolean[tmpLabels.length];
        for (int i = 0; i < tmpLabels.length; i++) {
            XYChart.Data<Number, Number> tmpData = new XYChart.Data<>(0.0, 0.0);
            boolean tmpIsRepresentative = tmpLabels[i] >= 0 && tmpRepresentatives[tmpLabels[i]] == i;
            tmpData.setNode(this.createMoleculeNode(i, tmpLabels[i], tmpIsRepresentative));
            this.isRepresentative[i] = tmpIsRepresentative;
            this.moleculeDataList.add(tmpData);
        }
        tmpMoleculeSeries.getData().addAll(this.moleculeDataList);
        XYChart.Series<Number, Number> tmpCentreSeries = new XYChart.Series<>();
        this.centreDataList = new ArrayList<>(aPcaResult.clusterVectorScores().length);
        for (int k = 0; k < aPcaResult.clusterVectorScores().length; k++) {
            XYChart.Data<Number, Number> tmpData = new XYChart.Data<>(0.0, 0.0);
            tmpData.setNode(this.createCentreNode(k));
            this.centreDataList.add(tmpData);
        }
        tmpCentreSeries.getData().addAll(this.centreDataList);
        this.chart.getData().addAll(List.of(tmpMoleculeSeries, tmpCentreSeries));
        //zoom with the scroll wheel around the mouse position, move by dragging
        this.chart.setOnScroll(anEvent -> {
            if (anEvent.getDeltaY() == 0.0) {
                return;
            }
            double tmpFactor = anEvent.getDeltaY() > 0.0 ? PcaView.ZOOM_FACTOR : 1.0 / PcaView.ZOOM_FACTOR;
            double tmpX = this.xAxis.getValueForDisplay(this.xAxis.sceneToLocal(anEvent.getSceneX(), anEvent.getSceneY()).getX()).doubleValue();
            double tmpY = this.yAxis.getValueForDisplay(this.yAxis.sceneToLocal(anEvent.getSceneX(), anEvent.getSceneY()).getY()).doubleValue();
            PcaView.setBounds(this.xAxis, tmpX - (tmpX - this.xAxis.getLowerBound()) * tmpFactor,
                    tmpX + (this.xAxis.getUpperBound() - tmpX) * tmpFactor);
            PcaView.setBounds(this.yAxis, tmpY - (tmpY - this.yAxis.getLowerBound()) * tmpFactor,
                    tmpY + (this.yAxis.getUpperBound() - tmpY) * tmpFactor);
            anEvent.consume();
        });
        this.chart.setOnMousePressed(anEvent -> this.dragStart = new double[] {anEvent.getSceneX(), anEvent.getSceneY(),
                this.xAxis.getLowerBound(), this.xAxis.getUpperBound(), this.yAxis.getLowerBound(), this.yAxis.getUpperBound()});
        this.chart.setOnMouseDragged(anEvent -> {
            if (this.dragStart == null || this.xAxis.getWidth() <= 0.0 || this.yAxis.getHeight() <= 0.0) {
                return;
            }
            double tmpDeltaX = (anEvent.getSceneX() - this.dragStart[0]) * (this.dragStart[3] - this.dragStart[2]) / this.xAxis.getWidth();
            double tmpDeltaY = (anEvent.getSceneY() - this.dragStart[1]) * (this.dragStart[5] - this.dragStart[4]) / this.yAxis.getHeight();
            PcaView.setBounds(this.xAxis, this.dragStart[2] - tmpDeltaX, this.dragStart[3] - tmpDeltaX);
            PcaView.setBounds(this.yAxis, this.dragStart[4] + tmpDeltaY, this.dragStart[5] + tmpDeltaY);
        });
        this.chart.setOnMouseReleased(anEvent -> this.dragStart = null);
        //loading vectors over the plot, redrawn whenever the axes have been laid out anew (zoom, move, resize)
        this.loadingVectorPane = new Pane();
        this.loadingVectorPane.setMouseTransparent(true);
        ChangeListener<Boolean> tmpRedrawListener = (anObservable, anOldValue, aNewValue) -> {
            if (!aNewValue) {
                this.updateLoadingVectors();
            }
        };
        this.xAxis.needsLayoutProperty().addListener(tmpRedrawListener);
        this.yAxis.needsLayoutProperty().addListener(tmpRedrawListener);
        this.chart.needsLayoutProperty().addListener(tmpRedrawListener);
        //</editor-fold>
        //<editor-fold desc="controls on the right">
        List<Integer> tmpComponents = new ArrayList<>();
        for (int c = 0; c < aPcaResult.getNumberOfComponents(); c++) {
            tmpComponents.add(c);
        }
        this.xComponentComboBox = this.createComponentComboBox(tmpComponents, 0);
        this.yComponentComboBox = this.createComponentComboBox(tmpComponents, Math.min(1, tmpComponents.size() - 1));
        this.centresCheckBox = new CheckBox(Message.get("PcaView.centresCheckBox.text"));
        this.centresCheckBox.setSelected(true);
        this.centresCheckBox.selectedProperty().addListener((anObservable, anOldValue, aNewValue) -> this.updatePointStyles());
        Button tmpResetViewButton = new Button(Message.get("PcaView.resetViewButton.text"));
        tmpResetViewButton.setOnAction(anEvent -> this.resetView());
        //loading vectors of chosen features in the 2D plot
        Button tmpLoadingVectorsButton = new Button(Message.get("PcaView.loadingVectorsButton.text"));
        tmpLoadingVectorsButton.setTooltip(new Tooltip(Message.get("PcaView.loadingVectorsButton.tooltip")));
        tmpLoadingVectorsButton.setOnAction(anEvent -> this.chooseLoadingVectors());
        GridPane tmpAxisGrid = new GridPane();
        tmpAxisGrid.setHgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpAxisGrid.setVgap(GuiDefinitions.GUI_SPACING_VALUE);
        tmpAxisGrid.addRow(0, new Label(Message.get("PcaView.xAxis.label")), this.xComponentComboBox);
        tmpAxisGrid.addRow(1, new Label(Message.get("PcaView.yAxis.label")), this.yComponentComboBox);
        this.clusterListView = new ListView<>(FXCollections.observableArrayList(
                tmpRankedClusters.stream().map(ClusteringRunResult.ClusterSummary::clusterIndex).toList()));
        Map<Integer, Integer> tmpSizeOfClusterMap = new HashMap<>();
        tmpRankedClusters.forEach(aCluster -> tmpSizeOfClusterMap.put(aCluster.clusterIndex(), aCluster.getSize()));
        this.clusterListView.setCellFactory(aListView -> this.createClusterCell(tmpSizeOfClusterMap));
        //the selected cluster is highlighted, the others are dimmed
        this.clusterListView.getSelectionModel().selectedItemProperty().addListener((anObservable, anOldValue, aNewValue) ->
                this.updatePointStyles());
        VBox.setVgrow(this.clusterListView, Priority.ALWAYS);
        Button tmpAllButton = new Button(Message.get("ClusteringResultView.markAllButton.text"));
        tmpAllButton.setOnAction(anEvent -> this.visiblePropertyOfClusterMap.values().forEach(aProperty -> aProperty.set(true)));
        Button tmpNoneButton = new Button(Message.get("ClusteringResultView.markNoneButton.text"));
        tmpNoneButton.setOnAction(anEvent -> this.visiblePropertyOfClusterMap.values().forEach(aProperty -> aProperty.set(false)));
        Button tmpClearHighlightButton = new Button(Message.get("PcaView.clearHighlightButton.text"));
        tmpClearHighlightButton.setOnAction(anEvent -> this.clusterListView.getSelectionModel().clearSelection());
        HBox tmpListButtonBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, tmpAllButton, tmpNoneButton, tmpClearHighlightButton);
        //colour by cluster or by a calculated descriptor
        this.descriptorStatistics = new ClusterDescriptorStatistics(aRunResult);
        List<String> tmpColourItems = new ArrayList<>();
        tmpColourItems.add(PcaView.CLUSTER_COLOUR_ITEM);
        tmpColourItems.addAll(this.descriptorStatistics.getDescriptorColumnNames());
        this.colourComboBox = new ComboBox<>(FXCollections.observableArrayList(tmpColourItems));
        this.colourComboBox.setCellFactory(aListView -> PcaView.createColourItemCell());
        this.colourComboBox.setButtonCell(PcaView.createColourItemCell());
        this.colourComboBox.setValue(PcaView.CLUSTER_COLOUR_ITEM);
        this.colourComboBox.setMaxWidth(Double.MAX_VALUE);
        this.colourComboBox.setDisable(tmpColourItems.size() < 2);
        this.colourComboBox.valueProperty().addListener((anObservable, anOldValue, aNewValue) -> this.updatePointColours());
        this.legendMinimumLabel = new Label();
        this.legendMaximumLabel = new Label();
        Rectangle tmpGradient = new Rectangle(PcaView.DETAIL_IMAGE_SIZE, 12.0);
        List<Stop> tmpStops = new ArrayList<>();
        for (int i = 0; i < PcaView.VIRIDIS_COLORS.length; i++) {
            tmpStops.add(new Stop((double) i / (PcaView.VIRIDIS_COLORS.length - 1), PcaView.VIRIDIS_COLORS[i]));
        }
        tmpGradient.setFill(new LinearGradient(0.0, 0.0, 1.0, 0.0, true, CycleMethod.NO_CYCLE, tmpStops));
        BorderPane tmpLegendLabels = new BorderPane();
        tmpLegendLabels.setLeft(this.legendMinimumLabel);
        tmpLegendLabels.setRight(this.legendMaximumLabel);
        Label tmpLegendHint = new Label(Message.get("PcaView.colourLegend.hint.text"));
        tmpLegendHint.setWrapText(true);
        this.colourLegendBox = new VBox(2.0, tmpGradient, tmpLegendLabels, tmpLegendHint);
        this.colourLegendBox.setVisible(false);
        this.colourLegendBox.managedProperty().bind(this.colourLegendBox.visibleProperty());
        this.detailImageView = new ImageView();
        this.detailLabel = new Label(Message.get("PcaView.detail.noSelection.text"));
        this.detailLabel.setWrapText(true);
        //2D scatter plot or 3D view on the first three components
        ToggleButton tmpTwoDimensionalToggleButton = new ToggleButton(Message.get("PcaView.twoDimensional.text"));
        this.threeDimensionalToggleButton = new ToggleButton(Message.get("PcaView.threeDimensional.text"));
        this.threeDimensionalToggleButton.setTooltip(new Tooltip(Message.get("PcaView.threeDimensional.tooltip")));
        ToggleGroup tmpDimensionGroup = new ToggleGroup();
        tmpTwoDimensionalToggleButton.setToggleGroup(tmpDimensionGroup);
        this.threeDimensionalToggleButton.setToggleGroup(tmpDimensionGroup);
        tmpTwoDimensionalToggleButton.setSelected(true);
        //one of both is always selected
        tmpDimensionGroup.selectedToggleProperty().addListener((anObservable, anOldValue, aNewValue) -> {
            if (aNewValue == null) {
                anOldValue.setSelected(true);
            }
        });
        if (aPcaResult.getNumberOfComponents() < 3 || !Platform.isSupported(ConditionalFeature.SCENE3D)) {
            this.threeDimensionalToggleButton.setDisable(true);
            this.threeDimensionalToggleButton.setTooltip(new Tooltip(Message.get("PcaView.threeDimensional.unavailable.tooltip")));
        }
        this.threeDimensionalToggleButton.selectedProperty().addListener((anObservable, anOldValue, aNewValue) -> this.showDimension(aNewValue));
        //the axes are fixed to PC1 to PC3 in 3D
        tmpAxisGrid.disableProperty().bind(this.threeDimensionalToggleButton.selectedProperty());
        tmpLoadingVectorsButton.disableProperty().bind(this.threeDimensionalToggleButton.selectedProperty());
        HBox tmpDimensionBox = new HBox(GuiDefinitions.GUI_BUTTON_SPACING_VALUE, tmpTwoDimensionalToggleButton, this.threeDimensionalToggleButton);
        VBox tmpRightBox = new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpDimensionBox, tmpAxisGrid, new Label(Message.get("PcaView.colourBy.label")),
                this.colourComboBox, this.colourLegendBox, this.centresCheckBox, tmpLoadingVectorsButton, tmpResetViewButton,
                new Separator(), new Label(Message.get("PcaView.clusterList.label")), this.clusterListView, tmpListButtonBox,
                new Separator(), this.detailImageView, this.detailLabel);
        tmpRightBox.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        tmpRightBox.setPrefWidth(PcaView.DETAIL_IMAGE_SIZE + 2 * GuiDefinitions.GUI_INSETS_VALUE);
        tmpRightBox.setMinWidth(PcaView.DETAIL_IMAGE_SIZE + 2 * GuiDefinitions.GUI_INSETS_VALUE);
        //</editor-fold>
        //<editor-fold desc="information at the top and bottom">
        StringBuilder tmpVariance = new StringBuilder();
        for (int c = 0; c < aPcaResult.getNumberOfComponents(); c++) {
            if (c > 0) {
                tmpVariance.append(", ");
            }
            tmpVariance.append(PcaView.getComponentName(c, aPcaResult));
        }
        Label tmpInfoLabel = new Label(String.format(Message.get("PcaView.info.text"), tmpLabels.length, tmpVariance));
        tmpInfoLabel.setWrapText(true);
        tmpInfoLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        this.hoverLabel = new Label(Message.get("PcaView.hover.text"));
        this.hoverLabel.setPadding(new Insets(GuiDefinitions.GUI_INSETS_VALUE));
        //</editor-fold>
        this.setTop(tmpInfoLabel);
        this.plotPane = new StackPane(this.chart, this.loadingVectorPane);
        this.setCenter(this.plotPane);
        this.setRight(tmpRightBox);
        this.setBottom(new VBox(tmpMeaningBox, this.hoverLabel));
        this.updateCoordinates();
        this.updatePointStyles();
        this.updatePointColours();
        this.resetView();
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Sets the action called with the cluster index of a double-clicked point, e.g. to open the cluster.
     *
     * @param anAction action, null for none
     */
    public void setClusterDoubleClickAction(IntConsumer anAction) {
        this.clusterDoubleClickAction = anAction;
    }
    //
    /**
     * Sets the axis bounds so that all visible points and, if shown, the cluster vectors fit into the plot.
     */
    public void resetView() {
        if (this.pca3dView != null) {
            this.pca3dView.resetView();
        }
        double tmpMinX = Double.POSITIVE_INFINITY;
        double tmpMaxX = Double.NEGATIVE_INFINITY;
        double tmpMinY = Double.POSITIVE_INFINITY;
        double tmpMaxY = Double.NEGATIVE_INFINITY;
        List<XYChart.Data<Number, Number>> tmpDataList = new ArrayList<>(this.moleculeDataList);
        tmpDataList.addAll(this.centreDataList);
        for (XYChart.Data<Number, Number> tmpData : tmpDataList) {
            if (tmpData.getNode() == null || !tmpData.getNode().isVisible()) {
                continue;
            }
            tmpMinX = Math.min(tmpMinX, tmpData.getXValue().doubleValue());
            tmpMaxX = Math.max(tmpMaxX, tmpData.getXValue().doubleValue());
            tmpMinY = Math.min(tmpMinY, tmpData.getYValue().doubleValue());
            tmpMaxY = Math.max(tmpMaxY, tmpData.getYValue().doubleValue());
        }
        if (tmpMinX > tmpMaxX) {
            //nothing visible
            tmpMinX = -1.0;
            tmpMaxX = 1.0;
            tmpMinY = -1.0;
            tmpMaxY = 1.0;
        }
        double tmpMarginX = Math.max((tmpMaxX - tmpMinX) * PcaView.VIEW_MARGIN, 1.0e-6);
        double tmpMarginY = Math.max((tmpMaxY - tmpMinY) * PcaView.VIEW_MARGIN, 1.0e-6);
        PcaView.setBounds(this.xAxis, tmpMinX - tmpMarginX, tmpMaxX + tmpMarginX);
        PcaView.setBounds(this.yAxis, tmpMinY - tmpMarginY, tmpMaxY + tmpMarginY);
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Creates the point of a molecule: a circle in the colour of its cluster, larger with black border for a
     * representative; hovering enlarges it and describes it, a click shows its structure, a double click calls the
     * cluster action.
     *
     * @param aRowIndex row index of the molecule
     * @param aClusterIndex cluster of the molecule, -1 if not clusterable
     * @param anIsRepresentative true for the representative of its cluster
     * @return point node
     */
    private Node createMoleculeNode(int aRowIndex, int aClusterIndex, boolean anIsRepresentative) {
        Circle tmpCircle = new Circle(anIsRepresentative ? PcaView.REPRESENTATIVE_RADIUS : PcaView.POINT_RADIUS,
                PcaView.getClusterColor(aClusterIndex));
        if (anIsRepresentative) {
            tmpCircle.setStroke(Color.BLACK);
            tmpCircle.setStrokeWidth(1.2);
        }
        tmpCircle.setOnMouseEntered(anEvent -> {
            tmpCircle.setScaleX(PcaView.HOVER_SCALE);
            tmpCircle.setScaleY(PcaView.HOVER_SCALE);
            tmpCircle.toFront();
            this.hoverLabel.setText(this.describeMolecule(aRowIndex));
        });
        tmpCircle.setOnMouseExited(anEvent -> {
            tmpCircle.setScaleX(1.0);
            tmpCircle.setScaleY(1.0);
        });
        tmpCircle.setOnMouseClicked(anEvent -> {
            if (anEvent.getButton() != MouseButton.PRIMARY) {
                return;
            }
            if (anEvent.getClickCount() == 2 && this.clusterDoubleClickAction != null && aClusterIndex >= 0) {
                this.clusterDoubleClickAction.accept(aClusterIndex);
            } else {
                this.showDetails(aRowIndex);
            }
        });
        return tmpCircle;
    }
    //
    /**
     * Creates the diamond of a cluster vector in the colour of its cluster with a black border.
     *
     * @param aClusterIndex cluster index
     * @return diamond node
     */
    private Node createCentreNode(int aClusterIndex) {
        double tmpSize = PcaView.CENTRE_SIZE;
        Polygon tmpDiamond = new Polygon(0.0, -tmpSize, tmpSize, 0.0, 0.0, tmpSize, -tmpSize, 0.0);
        tmpDiamond.setFill(PcaView.getClusterColor(aClusterIndex));
        tmpDiamond.setStroke(Color.BLACK);
        tmpDiamond.setStrokeWidth(1.5);
        tmpDiamond.setOnMouseEntered(anEvent -> this.hoverLabel.setText(String.format(Message.get("PcaView.centre.text"), aClusterIndex,
                this.getClusterSize(aClusterIndex))));
        tmpDiamond.setOnMouseClicked(anEvent -> {
            if (anEvent.getClickCount() == 2 && this.clusterDoubleClickAction != null) {
                this.clusterDoubleClickAction.accept(aClusterIndex);
            }
        });
        return tmpDiamond;
    }
    //
    /**
     * Creates a combo box to choose a component, showing its explained variance.
     *
     * @param aComponentList component indices
     * @param aSelectedComponent initially selected component
     * @return combo box
     */
    private ComboBox<Integer> createComponentComboBox(List<Integer> aComponentList, int aSelectedComponent) {
        ComboBox<Integer> tmpComboBox = new ComboBox<>(FXCollections.observableArrayList(aComponentList));
        tmpComboBox.setCellFactory(aListView -> new ListCell<>() {
            @Override
            protected void updateItem(Integer aComponent, boolean anIsEmpty) {
                super.updateItem(aComponent, anIsEmpty);
                this.setText(anIsEmpty || aComponent == null ? null : PcaView.getComponentName(aComponent, PcaView.this.pcaResult));
            }
        });
        tmpComboBox.setButtonCell(tmpComboBox.getCellFactory().call(null));
        tmpComboBox.setValue(aSelectedComponent);
        tmpComboBox.valueProperty().addListener((anObservable, anOldValue, aNewValue) -> {
            this.updateCoordinates();
            this.resetView();
        });
        return tmpComboBox;
    }
    //
    /**
     * Creates a cell of the cluster list: visibility check box, colour swatch, cluster index and size.
     *
     * @param aSizeOfClusterMap size of every cluster
     * @return list cell
     */
    private ListCell<Integer> createClusterCell(Map<Integer, Integer> aSizeOfClusterMap) {
        return new ListCell<>() {
            /**
             * Visibility check box of the cell.
             */
            private final CheckBox checkBox = new CheckBox();
            /**
             * Visibility property the check box is currently bound to.
             */
            private BooleanProperty boundProperty;
            //
            @Override
            protected void updateItem(Integer aClusterIndex, boolean anIsEmpty) {
                super.updateItem(aClusterIndex, anIsEmpty);
                if (this.boundProperty != null) {
                    this.checkBox.selectedProperty().unbindBidirectional(this.boundProperty);
                    this.boundProperty = null;
                }
                if (anIsEmpty || aClusterIndex == null) {
                    this.setGraphic(null);
                    this.setText(null);
                    return;
                }
                this.boundProperty = PcaView.this.visiblePropertyOfClusterMap.get(aClusterIndex);
                this.checkBox.selectedProperty().bindBidirectional(this.boundProperty);
                Rectangle tmpSwatch = new Rectangle(12.0, 12.0, PcaView.getClusterColor(aClusterIndex));
                HBox tmpGraphic = new HBox(GuiDefinitions.GUI_SPACING_VALUE, this.checkBox, tmpSwatch);
                tmpGraphic.setAlignment(Pos.CENTER_LEFT);
                this.setGraphic(tmpGraphic);
                this.setText(aClusterIndex < 0
                        ? String.format(Message.get("PcaView.notClusterable.text"), aSizeOfClusterMap.get(aClusterIndex))
                        : String.format(Message.get("PcaView.cluster.text"), aClusterIndex, aSizeOfClusterMap.get(aClusterIndex)));
            }
        };
    }
    //
    /**
     * Returns the main features of a component: the features with the largest weight (loading), each with "+" or "-" for
     * the direction in the raw values (sign of the correlation of the raw feature values with the component scores; the
     * sign of the loading if there are no raw values) and optionally the correlation r.
     *
     * @param aComponent component index
     * @param aNumber number of features
     * @param anIsWithCorrelation true to add the correlation r of every feature
     * @return description, e.g. "Heavy Atom Count +, TPSA -"
     */
    private String describeComponent(int aComponent, int aNumber, boolean anIsWithCorrelation) {
        List<String> tmpParts = new ArrayList<>();
        for (PrincipalComponentAnalysis.FeatureLoading tmpLoading : this.pcaResult.getTopLoadings(aComponent, aNumber)) {
            String tmpName = this.runResult.processedMatrix().getFeatureName(tmpLoading.featureIndex());
            double tmpCorrelation = this.getCorrelationWithScores(tmpName, aComponent);
            double tmpDirection = Double.isNaN(tmpCorrelation) ? tmpLoading.loading() : tmpCorrelation;
            String tmpShortName = tmpName.length() > PcaView.MAXIMUM_FEATURE_NAME_LENGTH
                    ? tmpName.substring(0, PcaView.MAXIMUM_FEATURE_NAME_LENGTH - 1) + "\u2026" : tmpName;
            String tmpPart = tmpShortName + (tmpDirection >= 0.0 ? " +" : " \u2212");
            if (anIsWithCorrelation && !Double.isNaN(tmpCorrelation)) {
                tmpPart += String.format(Locale.getDefault(), " (r %.2f)", tmpCorrelation);
            }
            tmpParts.add(tmpPart);
        }
        return String.join(", ", tmpParts);
    }
    //
    /**
     * Returns the main features of a component with full names, loadings and correlations, one per line, e.g. for a
     * tooltip.
     *
     * @param aComponent component index
     * @return description
     */
    private String describeComponentInFull(int aComponent) {
        StringBuilder tmpText = new StringBuilder();
        for (PrincipalComponentAnalysis.FeatureLoading tmpLoading : this.pcaResult.getTopLoadings(aComponent, PcaView.NUMBER_OF_MAIN_FEATURES)) {
            String tmpName = this.runResult.processedMatrix().getFeatureName(tmpLoading.featureIndex());
            double tmpCorrelation = this.getCorrelationWithScores(tmpName, aComponent);
            tmpText.append(String.format(Locale.getDefault(), Message.get("PcaView.loadings.line.text"), tmpName, tmpLoading.loading(),
                    Double.isNaN(tmpCorrelation) ? "-" : String.format(Locale.getDefault(), "%.2f", tmpCorrelation))).append(System.lineSeparator());
        }
        return tmpText.toString().trim();
    }
    //
    /**
     * Returns the Pearson correlation of the raw values of the feature with the given name with the scores of the given
     * component, i.e. the direction of the component in the units of the feature; missing values are left out.
     *
     * @param aFeatureName name of the feature (column of the raw matrix)
     * @param aComponent component index
     * @return correlation in [-1, 1], NaN if the feature has no raw column or no variance
     */
    private double getCorrelationWithScores(String aFeatureName, int aComponent) {
        FeatureMatrix tmpRaw = this.runResult.rawMatrix();
        int tmpColumn = -1;
        for (int j = 0; j < tmpRaw.getColumnCount() && tmpColumn < 0; j++) {
            if (tmpRaw.getFeatureName(j).equals(aFeatureName)) {
                tmpColumn = j;
            }
        }
        if (tmpColumn < 0) {
            return Double.NaN;
        }
        double[][] tmpScores = this.pcaResult.scores();
        double tmpSumX = 0.0;
        double tmpSumY = 0.0;
        int tmpCount = 0;
        for (int i = 0; i < tmpScores.length; i++) {
            float tmpValue = tmpRaw.getValue(i, tmpColumn);
            if (Float.isFinite(tmpValue)) {
                tmpSumX += tmpValue;
                tmpSumY += tmpScores[i][aComponent];
                tmpCount++;
            }
        }
        if (tmpCount < 2) {
            return Double.NaN;
        }
        double tmpMeanX = tmpSumX / tmpCount;
        double tmpMeanY = tmpSumY / tmpCount;
        double tmpCovariance = 0.0;
        double tmpVarianceX = 0.0;
        double tmpVarianceY = 0.0;
        for (int i = 0; i < tmpScores.length; i++) {
            float tmpValue = tmpRaw.getValue(i, tmpColumn);
            if (Float.isFinite(tmpValue)) {
                double tmpDeltaX = tmpValue - tmpMeanX;
                double tmpDeltaY = tmpScores[i][aComponent] - tmpMeanY;
                tmpCovariance += tmpDeltaX * tmpDeltaY;
                tmpVarianceX += tmpDeltaX * tmpDeltaX;
                tmpVarianceY += tmpDeltaY * tmpDeltaY;
            }
        }
        return tmpVarianceX > 0.0 && tmpVarianceY > 0.0 ? tmpCovariance / Math.sqrt(tmpVarianceX * tmpVarianceY) : Double.NaN;
    }
    //
    /**
     * Returns the descriptor the points are coloured by.
     *
     * @return descriptor column name, null if the points are coloured by cluster
     */
    private String getColourDescriptor() {
        String tmpValue = this.colourComboBox == null ? null : this.colourComboBox.getValue();
        return tmpValue == null || PcaView.CLUSTER_COLOUR_ITEM.equals(tmpValue) ? null : tmpValue;
    }
    //
    /**
     * Colours the points by cluster or by the chosen descriptor: the raw value is mapped to the viridis colour scale
     * between the 2nd and the 98th percentile of the descriptor (values outside get the end colours), missing values are
     * grey. The cluster centres are white while the points are coloured by a descriptor, since their colour would be
     * mistaken for a value. The legend shows the bounds of the scale.
     */
    private void updatePointColours() {
        String tmpDescriptor = this.getColourDescriptor();
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        this.colourScaleBounds = null;
        if (tmpDescriptor != null) {
            double[] tmpValues = new double[tmpLabels.length];
            List<Double> tmpFinite = new ArrayList<>(tmpLabels.length);
            for (int i = 0; i < tmpLabels.length; i++) {
                tmpValues[i] = this.descriptorStatistics.getValue(i, tmpDescriptor);
                if (Double.isFinite(tmpValues[i])) {
                    tmpFinite.add(tmpValues[i]);
                }
            }
            tmpFinite.sort(null);
            if (!tmpFinite.isEmpty()) {
                double tmpLower = tmpFinite.get((int) Math.floor(PcaView.COLOUR_SCALE_LOWER_PERCENTILE * (tmpFinite.size() - 1)));
                double tmpUpper = tmpFinite.get((int) Math.ceil(PcaView.COLOUR_SCALE_UPPER_PERCENTILE * (tmpFinite.size() - 1)));
                if (!(tmpUpper > tmpLower)) {
                    //few distinct values: the whole range
                    tmpLower = tmpFinite.getFirst();
                    tmpUpper = tmpFinite.getLast();
                }
                this.colourScaleBounds = new double[] {tmpLower, tmpUpper};
            }
            for (int i = 0; i < this.moleculeDataList.size(); i++) {
                ((Circle) this.moleculeDataList.get(i).getNode()).setFill(this.getDescriptorColor(tmpValues[i]));
            }
        } else {
            for (int i = 0; i < this.moleculeDataList.size(); i++) {
                ((Circle) this.moleculeDataList.get(i).getNode()).setFill(PcaView.getClusterColor(tmpLabels[i]));
            }
        }
        for (int k = 0; k < this.centreDataList.size(); k++) {
            ((Polygon) this.centreDataList.get(k).getNode()).setFill(tmpDescriptor == null ? PcaView.getClusterColor(k) : Color.WHITE);
        }
        this.colourLegendBox.setVisible(tmpDescriptor != null && this.colourScaleBounds != null);
        if (this.colourScaleBounds != null) {
            this.legendMinimumLabel.setText(String.format(Locale.getDefault(), "%.2f", this.colourScaleBounds[0]));
            this.legendMaximumLabel.setText(String.format(Locale.getDefault(), "%.2f", this.colourScaleBounds[1]));
        }
        this.update3dView();
    }
    //
    /**
     * Returns the colour of a descriptor value on the viridis colour scale between the bounds of the current scale.
     *
     * @param aValue raw descriptor value
     * @return colour, grey for a missing value
     */
    private Color getDescriptorColor(double aValue) {
        if (!Double.isFinite(aValue) || this.colourScaleBounds == null) {
            return Color.LIGHTGREY;
        }
        double tmpRange = this.colourScaleBounds[1] - this.colourScaleBounds[0];
        double tmpFraction = tmpRange > 0.0 ? (aValue - this.colourScaleBounds[0]) / tmpRange : 0.5;
        //256 steps are indistinguishable from a continuous scale and keep the palette of the 3D view small
        double tmpStep = Math.round(Math.min(1.0, Math.max(0.0, tmpFraction)) * 255.0) / 255.0;
        return PcaView.getViridisColor(tmpStep);
    }
    //
    /**
     * Sets the coordinates of all points to the scores of the chosen components.
     */
    private void updateCoordinates() {
        int tmpX = this.xComponentComboBox.getValue();
        int tmpY = this.yComponentComboBox.getValue();
        double[][] tmpScores = this.pcaResult.scores();
        for (int i = 0; i < this.moleculeDataList.size(); i++) {
            this.moleculeDataList.get(i).setXValue(tmpScores[i][tmpX]);
            this.moleculeDataList.get(i).setYValue(tmpScores[i][tmpY]);
        }
        double[][] tmpCentreScores = this.pcaResult.clusterVectorScores();
        for (int k = 0; k < this.centreDataList.size(); k++) {
            this.centreDataList.get(k).setXValue(tmpCentreScores[k][tmpX]);
            this.centreDataList.get(k).setYValue(tmpCentreScores[k][tmpY]);
        }
        this.xAxis.setLabel(this.axisLabels[tmpX]);
        this.yAxis.setLabel(this.axisLabels[tmpY]);
        this.updateLoadingVectors();
    }
    //
    /**
     * Opens the dialog to choose the features whose loading vectors are shown: a searchable list of all features of the
     * processed matrix with check boxes, sorted by the length of their vector on the shown components; "Main features"
     * checks the features with the largest weight in both components, "None" unchecks all.
     */
    private void chooseLoadingVectors() {
        int tmpX = this.xComponentComboBox.getValue();
        int tmpY = this.yComponentComboBox.getValue();
        double[][] tmpLoadings = this.pcaResult.loadings();
        FeatureMatrix tmpProcessed = this.runResult.processedMatrix();
        List<Integer> tmpFeatures = new ArrayList<>(tmpProcessed.getColumnCount());
        for (int j = 0; j < tmpProcessed.getColumnCount(); j++) {
            tmpFeatures.add(j);
        }
        tmpFeatures.sort(Comparator.comparingDouble((Integer j) -> -Math.hypot(tmpLoadings[tmpX][j], tmpLoadings[tmpY][j])));
        Map<Integer, BooleanProperty> tmpCheckedPropertyMap = new LinkedHashMap<>();
        for (Integer tmpFeature : tmpFeatures) {
            tmpCheckedPropertyMap.put(tmpFeature, new SimpleBooleanProperty(this.loadingVectorFeatureIndices.contains(tmpFeature)));
        }
        TextField tmpSearchField = new TextField();
        tmpSearchField.setPromptText(Message.get("ClusteringResultView.chooseDescriptors.search.prompt"));
        FilteredList<Integer> tmpFilteredFeatures = new FilteredList<>(FXCollections.observableArrayList(tmpFeatures));
        tmpSearchField.textProperty().addListener((anObservable, anOldValue, aNewValue) -> tmpFilteredFeatures.setPredicate(aFeature ->
                aNewValue == null || aNewValue.isBlank()
                        || tmpProcessed.getFeatureName(aFeature).toLowerCase(Locale.ROOT).contains(aNewValue.trim().toLowerCase(Locale.ROOT))));
        ListView<Integer> tmpListView = new ListView<>(tmpFilteredFeatures);
        //name, group and the loadings on the shown components
        tmpListView.setCellFactory(CheckBoxListCell.forListView(tmpCheckedPropertyMap::get, new StringConverter<>() {
            @Override
            public String toString(Integer aFeature) {
                return aFeature == null ? "" : String.format(Locale.getDefault(), Message.get("PcaView.loadingVectors.item.text"),
                        tmpProcessed.getFeatureName(aFeature), tmpProcessed.getFeatureGroupName(aFeature),
                        PcaView.getComponentLabel(tmpX), tmpLoadings[tmpX][aFeature], PcaView.getComponentLabel(tmpY), tmpLoadings[tmpY][aFeature]);
            }
            @Override
            public Integer fromString(String aString) {
                return null;
            }
        }));
        tmpListView.setPrefHeight(400);
        tmpListView.setPrefWidth(520);
        Button tmpMainButton = new Button(Message.get("PcaView.loadingVectors.mainButton.text"));
        tmpMainButton.setOnAction(anEvent -> {
            for (int tmpComponent : new int[] {tmpX, tmpY}) {
                for (PrincipalComponentAnalysis.FeatureLoading tmpLoading : this.pcaResult.getTopLoadings(tmpComponent, PcaView.NUMBER_OF_AXIS_FEATURES)) {
                    tmpCheckedPropertyMap.get(tmpLoading.featureIndex()).set(true);
                }
            }
        });
        Button tmpNoneButton = new Button(Message.get("ClusteringResultView.markNoneButton.text"));
        tmpNoneButton.setOnAction(anEvent -> tmpCheckedPropertyMap.values().forEach(aProperty -> aProperty.set(false)));
        HBox tmpSearchBox = new HBox(GuiDefinitions.GUI_SPACING_VALUE, tmpSearchField, tmpMainButton, tmpNoneButton);
        HBox.setHgrow(tmpSearchField, Priority.ALWAYS);
        Label tmpHintLabel = new Label(Message.get("PcaView.loadingVectors.hint.text"));
        tmpHintLabel.setWrapText(true);
        tmpHintLabel.setMaxWidth(520);
        Dialog<ButtonType> tmpDialog = new Dialog<>();
        if (this.getScene() != null) {
            tmpDialog.initOwner(this.getScene().getWindow());
        }
        tmpDialog.setTitle(Message.get("PcaView.loadingVectorsButton.text"));
        tmpDialog.setHeaderText(String.format(Message.get("PcaView.loadingVectors.header"), PcaView.getComponentLabel(tmpX),
                PcaView.getComponentLabel(tmpY)));
        tmpDialog.getDialogPane().setContent(new VBox(GuiDefinitions.GUI_SPACING_VALUE, tmpSearchBox, tmpListView, tmpHintLabel));
        tmpDialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        if (tmpDialog.showAndWait().orElse(ButtonType.CANCEL) != ButtonType.OK) {
            return;
        }
        this.loadingVectorFeatureIndices.clear();
        tmpCheckedPropertyMap.forEach((aFeature, aProperty) -> {
            if (aProperty.get()) {
                this.loadingVectorFeatureIndices.add(aFeature);
            }
        });
        this.updateLoadingVectors();
    }
    //
    /**
     * Draws the loading vectors of the chosen features over the 2D plot: an arrow from the origin of the components (the
     * mean of the molecules) in the direction of the feature's weights on the shown components, labelled with the
     * feature name. All vectors share one scale, chosen so that the longest possible vector of all features reaches
     * {@link #LOADING_VECTOR_LENGTH_SHARE} of the smaller score range; so the lengths can be compared, and a feature with
     * a long vector separates the molecules along it. The vectors are clipped to the plot area.
     */
    private void updateLoadingVectors() {
        if (this.loadingVectorPane == null) {
            return;
        }
        this.loadingVectorPane.getChildren().clear();
        if (this.loadingVectorFeatureIndices.isEmpty() || !this.chart.isVisible() || this.xAxis.getWidth() <= 0.0 || this.yAxis.getHeight() <= 0.0) {
            return;
        }
        int tmpX = this.xComponentComboBox.getValue();
        int tmpY = this.yComponentComboBox.getValue();
        double[][] tmpLoadings = this.pcaResult.loadings();
        double tmpMaximumScoreX = 0.0;
        double tmpMaximumScoreY = 0.0;
        for (double[] tmpScore : this.pcaResult.scores()) {
            tmpMaximumScoreX = Math.max(tmpMaximumScoreX, Math.abs(tmpScore[tmpX]));
            tmpMaximumScoreY = Math.max(tmpMaximumScoreY, Math.abs(tmpScore[tmpY]));
        }
        double tmpMaximumLength = 0.0;
        for (int j = 0; j < tmpLoadings[tmpX].length; j++) {
            tmpMaximumLength = Math.max(tmpMaximumLength, Math.hypot(tmpLoadings[tmpX][j], tmpLoadings[tmpY][j]));
        }
        if (tmpMaximumLength <= 0.0) {
            return;
        }
        double tmpScale = PcaView.LOADING_VECTOR_LENGTH_SHARE * Math.min(tmpMaximumScoreX, tmpMaximumScoreY) / tmpMaximumLength;
        //plot area in the coordinates of the pane, the vectors are clipped to it
        Point2D tmpTopLeft = this.loadingVectorPane.sceneToLocal(this.xAxis.localToScene(0.0, 0.0).getX(),
                this.yAxis.localToScene(0.0, 0.0).getY());
        this.loadingVectorPane.setClip(new Rectangle(tmpTopLeft.getX(), tmpTopLeft.getY(), this.xAxis.getWidth(), this.yAxis.getHeight()));
        Point2D tmpOrigin = this.toLoadingVectorPane(0.0, 0.0);
        FeatureMatrix tmpProcessed = this.runResult.processedMatrix();
        for (int tmpFeature : this.loadingVectorFeatureIndices) {
            Point2D tmpTip = this.toLoadingVectorPane(tmpScale * tmpLoadings[tmpX][tmpFeature], tmpScale * tmpLoadings[tmpY][tmpFeature]);
            Line tmpLine = new Line(tmpOrigin.getX(), tmpOrigin.getY(), tmpTip.getX(), tmpTip.getY());
            tmpLine.setStroke(PcaView.LOADING_VECTOR_COLOR);
            tmpLine.setStrokeWidth(1.5);
            this.loadingVectorPane.getChildren().add(tmpLine);
            double tmpAngle = Math.atan2(tmpTip.getY() - tmpOrigin.getY(), tmpTip.getX() - tmpOrigin.getX());
            if (tmpOrigin.distance(tmpTip) > PcaView.ARROW_HEAD_LENGTH) {
                for (double tmpSide : new double[] {-1.0, 1.0}) {
                    double tmpHeadAngle = tmpAngle + Math.PI + tmpSide * PcaView.ARROW_HEAD_ANGLE;
                    Line tmpHead = new Line(tmpTip.getX(), tmpTip.getY(), tmpTip.getX() + PcaView.ARROW_HEAD_LENGTH * Math.cos(tmpHeadAngle),
                            tmpTip.getY() + PcaView.ARROW_HEAD_LENGTH * Math.sin(tmpHeadAngle));
                    tmpHead.setStroke(PcaView.LOADING_VECTOR_COLOR);
                    tmpHead.setStrokeWidth(1.5);
                    this.loadingVectorPane.getChildren().add(tmpHead);
                }
            }
            //the name a little beyond the tip, on the side the arrow points to
            String tmpName = tmpProcessed.getFeatureName(tmpFeature);
            Text tmpLabel = new Text(tmpName.length() > PcaView.MAXIMUM_FEATURE_NAME_LENGTH
                    ? tmpName.substring(0, PcaView.MAXIMUM_FEATURE_NAME_LENGTH - 1) + "\u2026" : tmpName);
            tmpLabel.setFill(PcaView.LOADING_VECTOR_COLOR);
            tmpLabel.setStyle("-fx-font-weight: bold;");
            double tmpWidth = tmpLabel.getLayoutBounds().getWidth();
            tmpLabel.setX(tmpTip.getX() + 6.0 * Math.cos(tmpAngle) - (Math.cos(tmpAngle) < 0.0 ? tmpWidth : 0.0));
            tmpLabel.setY(tmpTip.getY() + 6.0 * Math.sin(tmpAngle) + (Math.sin(tmpAngle) > 0.0 ? 10.0 : 0.0));
            this.loadingVectorPane.getChildren().add(tmpLabel);
        }
    }
    //
    /**
     * Returns the position of the given component values in the coordinates of the loading vector pane.
     *
     * @param aXValue value on the x axis
     * @param aYValue value on the y axis
     * @return position in the pane
     */
    private Point2D toLoadingVectorPane(double aXValue, double aYValue) {
        double tmpSceneX = this.xAxis.localToScene(this.xAxis.getDisplayPosition(aXValue), 0.0).getX();
        double tmpSceneY = this.yAxis.localToScene(0.0, this.yAxis.getDisplayPosition(aYValue)).getY();
        return this.loadingVectorPane.sceneToLocal(tmpSceneX, tmpSceneY);
    }
    //
    /**
     * Sets visibility and opacity of all points: hidden clusters are invisible, and while a cluster is highlighted
     * (selected in the list) the points of the other clusters are dimmed; the cluster vectors follow their clusters and
     * the check box.
     */
    private void updatePointStyles() {
        Integer tmpHighlighted = this.clusterListView == null ? null : this.clusterListView.getSelectionModel().getSelectedItem();
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        for (int i = 0; i < this.moleculeDataList.size(); i++) {
            Node tmpNode = this.moleculeDataList.get(i).getNode();
            tmpNode.setVisible(this.isClusterVisible(tmpLabels[i]));
            tmpNode.setOpacity(tmpHighlighted == null || tmpHighlighted == tmpLabels[i] ? PcaView.POINT_OPACITY : PcaView.DIMMED_OPACITY);
        }
        boolean tmpIsShowingCentres = this.centresCheckBox == null || this.centresCheckBox.isSelected();
        for (int k = 0; k < this.centreDataList.size(); k++) {
            Node tmpNode = this.centreDataList.get(k).getNode();
            tmpNode.setVisible(tmpIsShowingCentres && this.isClusterVisible(k));
            tmpNode.setOpacity(tmpHighlighted == null || tmpHighlighted == k ? 1.0 : PcaView.DIMMED_OPACITY);
            if (tmpNode.isVisible()) {
                tmpNode.toFront();
            }
        }
        this.update3dView();
    }
    //
    /**
     * Shows the 2D scatter plot or the 3D view; the 3D view is created when it is shown for the first time.
     *
     * @param anIsThreeDimensional true for the 3D view
     */
    private void showDimension(boolean anIsThreeDimensional) {
        if (anIsThreeDimensional && this.pca3dView == null) {
            this.pca3dView = new Pca3dView(this.pcaResult, this.axisLabels);
            int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
            this.pca3dView.setMoleculeActions(
                    aRowIndex -> this.hoverLabel.setText(aRowIndex < 0 ? Message.get("PcaView.hover.text") : this.describeMolecule(aRowIndex)),
                    this::showDetails,
                    aRowIndex -> {
                        if (this.clusterDoubleClickAction != null && tmpLabels[aRowIndex] >= 0) {
                            this.clusterDoubleClickAction.accept(tmpLabels[aRowIndex]);
                        }
                    });
            this.pca3dView.setCentreActions(
                    aClusterIndex -> this.hoverLabel.setText(String.format(Message.get("PcaView.centre.text"), aClusterIndex,
                            this.getClusterSize(aClusterIndex))),
                    aClusterIndex -> {
                        if (this.clusterDoubleClickAction != null) {
                            this.clusterDoubleClickAction.accept(aClusterIndex);
                        }
                    });
            this.plotPane.getChildren().add(this.pca3dView);
        }
        this.chart.setVisible(!anIsThreeDimensional);
        //the loading vectors are only drawn in the 2D plot
        this.loadingVectorPane.setVisible(!anIsThreeDimensional);
        if (this.pca3dView != null) {
            this.pca3dView.setVisible(anIsThreeDimensional);
        }
        this.update3dView();
    }
    //
    /**
     * Passes the current colours and visibility of the points and centres to the 3D view, if it is shown; dimmed
     * clusters get a pale colour, since the points of one mesh cannot have different transparency.
     */
    private void update3dView() {
        if (this.pca3dView == null || !this.pca3dView.isVisible()) {
            return;
        }
        Integer tmpHighlighted = this.clusterListView.getSelectionModel().getSelectedItem();
        int[] tmpLabels = this.runResult.clusteringResult().clusterLabels();
        Color[] tmpColors = new Color[this.moleculeDataList.size()];
        boolean[] tmpIsVisible = new boolean[this.moleculeDataList.size()];
        for (int i = 0; i < tmpColors.length; i++) {
            Color tmpColor = (Color) ((Circle) this.moleculeDataList.get(i).getNode()).getFill();
            tmpColors[i] = tmpHighlighted == null || tmpHighlighted == tmpLabels[i] ? tmpColor : PcaView.getPaleColor(tmpColor);
            tmpIsVisible[i] = this.isClusterVisible(tmpLabels[i]);
        }
        boolean tmpIsShowingCentres = this.centresCheckBox.isSelected();
        Color[] tmpCentreColors = new Color[this.centreDataList.size()];
        boolean[] tmpIsCentreVisible = new boolean[this.centreDataList.size()];
        for (int k = 0; k < tmpCentreColors.length; k++) {
            Color tmpColor = (Color) ((Polygon) this.centreDataList.get(k).getNode()).getFill();
            tmpCentreColors[k] = tmpHighlighted == null || tmpHighlighted == k ? tmpColor : PcaView.getPaleColor(tmpColor);
            tmpIsCentreVisible[k] = tmpIsShowingCentres && this.isClusterVisible(k);
        }
        this.pca3dView.update(tmpColors, tmpIsVisible, this.isRepresentative, tmpCentreColors, tmpIsCentreVisible);
    }
    //
    /**
     * Returns whether the given cluster is shown.
     *
     * @param aClusterIndex cluster index, -1 for the not clusterable molecules
     * @return true if shown
     */
    private boolean isClusterVisible(int aClusterIndex) {
        BooleanProperty tmpProperty = this.visiblePropertyOfClusterMap.get(aClusterIndex);
        return tmpProperty == null || tmpProperty.get();
    }
    //
    /**
     * Returns the size of the given cluster.
     *
     * @param aClusterIndex cluster index
     * @return number of molecules
     */
    private int getClusterSize(int aClusterIndex) {
        int tmpSize = 0;
        for (int tmpLabel : this.runResult.clusteringResult().clusterLabels()) {
            if (tmpLabel == aClusterIndex) {
                tmpSize++;
            }
        }
        return tmpSize;
    }
    //
    /**
     * Returns the description of a molecule: name, cluster and coordinates on the chosen components.
     *
     * @param aRowIndex row index of the molecule
     * @return description
     */
    private String describeMolecule(int aRowIndex) {
        int tmpCluster = this.runResult.clusteringResult().clusterLabels()[aRowIndex];
        String tmpClusterText = tmpCluster < 0 ? Message.get("ClusteringResultView.notClusterable.text")
                : Message.get("OverviewView.titleOfDataSource.cluster") + " " + tmpCluster;
        double[] tmpScores = this.pcaResult.scores()[aRowIndex];
        String tmpDescription = String.format(Locale.getDefault(), "%s - %s - %s %.3f, %s %.3f", this.runResult.molecules().get(aRowIndex).getName(),
                tmpClusterText, "PC" + (this.xComponentComboBox.getValue() + 1), tmpScores[this.xComponentComboBox.getValue()],
                "PC" + (this.yComponentComboBox.getValue() + 1), tmpScores[this.yComponentComboBox.getValue()]);
        //the value of the descriptor the points are coloured by
        String tmpDescriptor = this.getColourDescriptor();
        if (tmpDescriptor != null) {
            double tmpValue = this.descriptorStatistics.getValue(aRowIndex, tmpDescriptor);
            tmpDescription += String.format(Locale.getDefault(), " - %s %s", tmpDescriptor,
                    Double.isNaN(tmpValue) ? "-" : String.format(Locale.getDefault(), "%.2f", tmpValue));
        }
        return tmpDescription;
    }
    //
    /**
     * Shows the structure and description of the given molecule on the right.
     *
     * @param aRowIndex row index of the molecule
     */
    private void showDetails(int aRowIndex) {
        MoleculeDataModel tmpMolecule = this.runResult.molecules().get(aRowIndex);
        try {
            this.detailImageView.setImage(DepictionUtil.depictImageWithTransparentBackground(tmpMolecule.getAtomContainer(), 1.0,
                    PcaView.DETAIL_IMAGE_SIZE, PcaView.DETAIL_IMAGE_SIZE, false));
        } catch (CDKException anException) {
            PcaView.LOGGER.log(Level.WARNING, anException.toString(), anException);
            this.detailImageView.setImage(DepictionUtil.depictErrorImage(anException.getMessage(),
                    (int) PcaView.DETAIL_IMAGE_SIZE, (int) PcaView.DETAIL_IMAGE_SIZE));
        }
        this.detailLabel.setText(this.describeMolecule(aRowIndex));
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns a pale version of the given colour for dimmed points of the 3D view.
     *
     * @param aColor colour
     * @return colour mixed with mostly white
     */
    private static Color getPaleColor(Color aColor) {
        return aColor.interpolate(Color.WHITE, 0.85);
    }
    //
    /**
     * Returns the colour of the given position on the viridis colour scale, interpolated between its anchor colours.
     *
     * @param aFraction position in [0, 1], 0 = lowest value
     * @return colour
     */
    private static Color getViridisColor(double aFraction) {
        double tmpPosition = aFraction * (PcaView.VIRIDIS_COLORS.length - 1);
        int tmpIndex = Math.min((int) Math.floor(tmpPosition), PcaView.VIRIDIS_COLORS.length - 2);
        return PcaView.VIRIDIS_COLORS[tmpIndex].interpolate(PcaView.VIRIDIS_COLORS[tmpIndex + 1], tmpPosition - tmpIndex);
    }
    //
    /**
     * Creates a cell of the colour combo box: "Cluster" for the cluster colours, otherwise the descriptor name.
     *
     * @return list cell
     */
    private static ListCell<String> createColourItemCell() {
        return new ListCell<>() {
            @Override
            protected void updateItem(String anItem, boolean anIsEmpty) {
                super.updateItem(anItem, anIsEmpty);
                this.setText(anIsEmpty || anItem == null ? null
                        : PcaView.CLUSTER_COLOUR_ITEM.equals(anItem) ? Message.get("PcaView.colourByCluster.text") : anItem);
            }
        };
    }
    //
    /**
     * Returns the colour of a cluster: hues by the golden angle, so that neighbouring indices differ clearly; grey for
     * the not clusterable molecules.
     *
     * @param aClusterIndex cluster index, -1 for the not clusterable molecules
     * @return colour
     */
    private static Color getClusterColor(int aClusterIndex) {
        if (aClusterIndex < 0) {
            return Color.GREY;
        }
        return Color.hsb((aClusterIndex * PcaView.GOLDEN_ANGLE) % 360.0, 0.75, 0.85);
    }
    //
    /**
     * Returns the short name of a component, e.g. "PC1".
     *
     * @param aComponent component index
     * @return short name
     */
    private static String getComponentLabel(int aComponent) {
        return "PC" + (aComponent + 1);
    }
    //
    /**
     * Returns the name of a component with its share of the explained variance, e.g. "PC1 (18.3 %)".
     *
     * @param aComponent component index
     * @param aPcaResult PCA result
     * @return component name
     */
    private static String getComponentName(int aComponent, PrincipalComponentAnalysis.PcaResult aPcaResult) {
        return String.format(Locale.getDefault(), "PC%d (%.1f %%)", aComponent + 1, 100.0 * aPcaResult.explainedVarianceRatios()[aComponent]);
    }
    //
    /**
     * Sets the bounds of the given axis and a tick unit that gives about ten ticks.
     *
     * @param anAxis axis
     * @param aLowerBound lower bound
     * @param anUpperBound upper bound, greater than the lower bound
     */
    private static void setBounds(AlignedNumberAxis anAxis, double aLowerBound, double anUpperBound) {
        if (!(anUpperBound > aLowerBound)) {
            return;
        }
        anAxis.setLowerBound(aLowerBound);
        anAxis.setUpperBound(anUpperBound);
        anAxis.setTickUnit(PcaView.getNiceTickUnit((anUpperBound - aLowerBound) / 10.0));
    }
    //
    /**
     * Rounds the given tick unit to 1, 2 or 5 times a power of ten.
     *
     * @param aRawUnit raw tick unit, greater than 0
     * @return rounded tick unit
     */
    private static double getNiceTickUnit(double aRawUnit) {
        double tmpPower = Math.pow(10.0, Math.floor(Math.log10(aRawUnit)));
        double tmpFraction = aRawUnit / tmpPower;
        double tmpUnit = (tmpFraction < 1.5 ? 1.0 : tmpFraction < 3.5 ? 2.0 : tmpFraction < 7.5 ? 5.0 : 10.0) * tmpPower;
        //one significant digit, so that no value like 0.30000000000000004 reaches the axis labels
        return new BigDecimal(tmpUnit).round(new MathContext(1)).doubleValue();
    }
    //</editor-fold>
}
