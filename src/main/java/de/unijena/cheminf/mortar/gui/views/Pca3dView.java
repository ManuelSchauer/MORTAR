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

import de.unijena.cheminf.mortar.model.clustering.PrincipalComponentAnalysis;

import javafx.geometry.Bounds;
import javafx.geometry.Point2D;
import javafx.geometry.Point3D;
import javafx.scene.AmbientLight;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.PerspectiveCamera;
import javafx.scene.SceneAntialiasing;
import javafx.scene.SubScene;
import javafx.scene.control.Label;
import javafx.scene.image.PixelWriter;
import javafx.scene.image.WritableImage;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;
import javafx.scene.input.PickResult;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Cylinder;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.function.IntConsumer;

/**
 * Three-dimensional view of a PCA of a clustering result on the first three principal components. All molecule points
 * form one mesh (a small octahedron per molecule, coloured via a palette texture), so that also tens of thousands of
 * molecules can be rotated smoothly; the cluster centres form a second mesh with larger octahedra. Dragging rotates the
 * view, the scroll wheel zooms; the molecule or centre under the mouse is determined from the picked triangle and passed
 * to the actions set by the {@link PcaView}, which also provides colours and visibility (see {@link #update}).
 * <p>
 * PC1 is drawn to the right, PC2 upwards and PC3 towards the viewer; all components are scaled with the same factor,
 * so that the proportions of the projection are kept. The axes have tick marks at round values. Their labels are not
 * part of the 3D scene, where text would be rotated and scaled by the perspective and become blurred, but sharp 2D labels
 * in an overlay that are moved to the projected positions of their anchors after every rotation, zoom and resize.
 * Loading vectors of chosen features can be shown as arrows from the origin (see {@link #setLoadingVectors}), labelled
 * the same way.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class Pca3dView extends Pane {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Half the edge length of the cube the points are scaled into.
     */
    private static final double HALF_SIZE = 200.0;
    /**
     * Radius of the octahedron of a molecule.
     */
    private static final double POINT_RADIUS = 1.6;
    /**
     * Radius of the octahedron of a representative.
     */
    private static final double REPRESENTATIVE_RADIUS = 3.2;
    /**
     * Radius of the octahedron of a cluster centre.
     */
    private static final double CENTRE_RADIUS = 5.0;
    /**
     * Number of triangles of an octahedron.
     */
    private static final int FACES_PER_POINT = 8;
    /**
     * Distance of the camera at the start.
     */
    private static final double CAMERA_DISTANCE = 900.0;
    /**
     * Rotation around the x-axis at the start, in degrees.
     */
    private static final double START_ANGLE_X = -20.0;
    /**
     * Rotation around the y-axis at the start, in degrees.
     */
    private static final double START_ANGLE_Y = -30.0;
    /**
     * Degrees of rotation per pixel of dragging.
     */
    private static final double ROTATION_PER_PIXEL = 0.4;
    /**
     * Factor of one scroll step of the zoom.
     */
    private static final double ZOOM_FACTOR = 0.9;
    /**
     * Radius of the axis cylinders.
     */
    private static final double AXIS_RADIUS = 0.7;
    /**
     * Radius of the tick marks on the axes.
     */
    private static final double TICK_RADIUS = 1.6;
    /**
     * Approximate number of tick marks per half axis.
     */
    private static final int TICKS_PER_HALF_AXIS = 4;
    /**
     * Radius of the shaft of a loading vector.
     */
    private static final double LOADING_VECTOR_RADIUS = 1.2;
    /**
     * Radius of the sphere at the tip of a loading vector.
     */
    private static final double LOADING_VECTOR_TIP_RADIUS = 3.5;
    /**
     * Colour of the loading vectors, the same as in the 2D plot.
     */
    private static final Color LOADING_VECTOR_COLOR = Color.web("#b2182b");
    /**
     * Style of the name of an axis.
     */
    private static final String AXIS_NAME_STYLE = "-fx-font-size: 13px; -fx-font-weight: bold; -fx-text-fill: #404040;";
    /**
     * Style of the value of a tick mark.
     */
    private static final String TICK_VALUE_STYLE = "-fx-font-size: 10px; -fx-text-fill: #606060;";
    /**
     * Style of the name of a loading vector.
     */
    private static final String LOADING_VECTOR_NAME_STYLE = "-fx-font-size: 12px; -fx-font-weight: bold; -fx-text-fill: #b2182b;";
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * PCA whose first three components are shown.
     */
    private final PrincipalComponentAnalysis.PcaResult pcaResult;
    /**
     * Group of all 3D content, rotated by dragging.
     */
    private final Group world;
    /**
     * Rotation around the x-axis.
     */
    private final Rotate rotateX;
    /**
     * Rotation around the y-axis.
     */
    private final Rotate rotateY;
    /**
     * Camera of the sub scene.
     */
    private final PerspectiveCamera camera;
    /**
     * Mesh view of the molecules.
     */
    private final MeshView moleculeMeshView;
    /**
     * Mesh view of the cluster centres.
     */
    private final MeshView centreMeshView;
    /**
     * Factor that scales the scores into the cube.
     */
    private final double scale;
    /**
     * Overlay with the sharp 2D labels of the axes and tick marks.
     */
    private final Pane labelOverlay;
    /**
     * Anchor of every label in the 3D scene, in the order of the labels of the overlay.
     */
    private final List<Node> labelAnchors;
    /**
     * Number of labels of the axes and tick marks; the labels of the loading vectors follow them in the overlay.
     */
    private final int numberOfAxisLabels;
    /**
     * Group of the loading vectors in the 3D scene.
     */
    private final Group loadingVectorGroup;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Row index of the molecule of every point of the molecule mesh, in mesh order.
     */
    private int[] rowOfPoint;
    /**
     * Cluster index of every point of the centre mesh, in mesh order.
     */
    private int[] clusterOfCentre;
    /**
     * Action called with the row index of the hovered molecule, -1 if none.
     */
    private IntConsumer moleculeHoverAction;
    /**
     * Action called with the row index of a clicked molecule.
     */
    private IntConsumer moleculeClickAction;
    /**
     * Action called with the row index of a double-clicked molecule.
     */
    private IntConsumer moleculeDoubleClickAction;
    /**
     * Action called with the cluster index of a hovered centre.
     */
    private IntConsumer centreHoverAction;
    /**
     * Action called with the cluster index of a double-clicked centre.
     */
    private IntConsumer centreDoubleClickAction;
    /**
     * Mouse position and rotation at the start of a drag.
     */
    private double[] dragStart;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor, builds the scene with the axes; the points are added by {@link #update}.
     *
     * @param aPcaResult PCA with at least three components
     * @param anAxisNameArray names of the first three axes, e.g. with explained variance and main features
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the PCA has fewer than three components or fewer than three names are given
     */
    public Pca3dView(PrincipalComponentAnalysis.PcaResult aPcaResult, String[] anAxisNameArray) throws NullPointerException, IllegalArgumentException {
        super();
        this.pcaResult = Objects.requireNonNull(aPcaResult, "Given PCA result is null.");
        Objects.requireNonNull(anAxisNameArray, "Given axis names are null.");
        if (anAxisNameArray.length < 3) {
            throw new IllegalArgumentException("Three axis names are needed.");
        }
        if (aPcaResult.getNumberOfComponents() < 3) {
            throw new IllegalArgumentException("A 3D view needs three components.");
        }
        //one factor for all components, so that the proportions are kept
        double tmpMaximum = 0.0;
        for (double[] tmpScore : aPcaResult.scores()) {
            for (int c = 0; c < 3; c++) {
                tmpMaximum = Math.max(tmpMaximum, Math.abs(tmpScore[c]));
            }
        }
        this.scale = tmpMaximum > 0.0 ? Pca3dView.HALF_SIZE / tmpMaximum : 1.0;
        this.moleculeMeshView = new MeshView();
        this.centreMeshView = new MeshView();
        this.labelOverlay = new Pane();
        this.labelOverlay.setMouseTransparent(true);
        this.labelAnchors = new ArrayList<>();
        this.loadingVectorGroup = new Group();
        this.world = new Group(this.moleculeMeshView, this.centreMeshView, this.loadingVectorGroup);
        this.world.getChildren().addAll(this.createAxes(tmpMaximum, anAxisNameArray));
        this.numberOfAxisLabels = this.labelAnchors.size();
        this.rotateX = new Rotate(Pca3dView.START_ANGLE_X, Rotate.X_AXIS);
        this.rotateY = new Rotate(Pca3dView.START_ANGLE_Y, Rotate.Y_AXIS);
        this.world.getTransforms().addAll(this.rotateX, this.rotateY);
        //flat colours without shading, so that the colours can be compared
        AmbientLight tmpLight = new AmbientLight(Color.WHITE);
        Group tmpRoot = new Group(this.world, tmpLight);
        this.camera = new PerspectiveCamera(true);
        this.camera.setNearClip(0.1);
        this.camera.setFarClip(20000.0);
        this.camera.setTranslateZ(-Pca3dView.CAMERA_DISTANCE);
        SubScene tmpSubScene = new SubScene(tmpRoot, 400.0, 400.0, true, SceneAntialiasing.BALANCED);
        tmpSubScene.setFill(Color.WHITE);
        tmpSubScene.setCamera(this.camera);
        tmpSubScene.widthProperty().bind(this.widthProperty());
        tmpSubScene.heightProperty().bind(this.heightProperty());
        this.getChildren().addAll(tmpSubScene, this.labelOverlay);
        this.setMinSize(0.0, 0.0);
        //the labels follow every rotation and zoom
        this.rotateX.angleProperty().addListener((anObservable, anOldValue, aNewValue) -> this.updateLabelPositions());
        this.rotateY.angleProperty().addListener((anObservable, anOldValue, aNewValue) -> this.updateLabelPositions());
        this.camera.translateZProperty().addListener((anObservable, anOldValue, aNewValue) -> this.updateLabelPositions());
        //rotation by dragging, zoom with the scroll wheel
        tmpSubScene.setOnMousePressed(anEvent -> this.dragStart = new double[] {anEvent.getSceneX(), anEvent.getSceneY(),
                this.rotateX.getAngle(), this.rotateY.getAngle()});
        tmpSubScene.setOnMouseDragged(anEvent -> {
            if (this.dragStart == null) {
                return;
            }
            this.rotateY.setAngle(this.dragStart[3] - (anEvent.getSceneX() - this.dragStart[0]) * Pca3dView.ROTATION_PER_PIXEL);
            this.rotateX.setAngle(this.dragStart[2] + (anEvent.getSceneY() - this.dragStart[1]) * Pca3dView.ROTATION_PER_PIXEL);
        });
        tmpSubScene.setOnMouseReleased(anEvent -> this.dragStart = null);
        tmpSubScene.setOnScroll(anEvent -> {
            if (anEvent.getDeltaY() == 0.0) {
                return;
            }
            double tmpFactor = anEvent.getDeltaY() > 0.0 ? Pca3dView.ZOOM_FACTOR : 1.0 / Pca3dView.ZOOM_FACTOR;
            this.camera.setTranslateZ(Math.min(-10.0, this.camera.getTranslateZ() * tmpFactor));
            anEvent.consume();
        });
        tmpSubScene.setOnMouseMoved(this::handleMouseMoved);
        tmpSubScene.setOnMouseClicked(this::handleMouseClicked);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Builds the meshes of the given visible molecules and centres with the given colours.
     *
     * @param aPointColors colour of every molecule by row index
     * @param anIsVisible visibility of every molecule by row index
     * @param anIsRepresentative whether every molecule is a representative, by row index
     * @param aCentreColors colour of every cluster centre by cluster index
     * @param anIsCentreVisible visibility of every cluster centre by cluster index
     * @throws NullPointerException if a parameter is null
     */
    public void update(Color[] aPointColors, boolean[] anIsVisible, boolean[] anIsRepresentative, Color[] aCentreColors,
                       boolean[] anIsCentreVisible) throws NullPointerException {
        Objects.requireNonNull(aPointColors, "Given colours are null.");
        List<Integer> tmpRows = new ArrayList<>();
        for (int i = 0; i < aPointColors.length; i++) {
            if (anIsVisible[i]) {
                tmpRows.add(i);
            }
        }
        this.rowOfPoint = tmpRows.stream().mapToInt(Integer::intValue).toArray();
        double[][] tmpPositions = new double[this.rowOfPoint.length][];
        double[] tmpRadii = new double[this.rowOfPoint.length];
        Color[] tmpColors = new Color[this.rowOfPoint.length];
        for (int p = 0; p < this.rowOfPoint.length; p++) {
            int tmpRow = this.rowOfPoint[p];
            tmpPositions[p] = this.toScene(this.pcaResult.scores()[tmpRow]);
            tmpRadii[p] = anIsRepresentative[tmpRow] ? Pca3dView.REPRESENTATIVE_RADIUS : Pca3dView.POINT_RADIUS;
            tmpColors[p] = aPointColors[tmpRow];
        }
        Pca3dView.fillMeshView(this.moleculeMeshView, tmpPositions, tmpRadii, tmpColors);
        List<Integer> tmpClusters = new ArrayList<>();
        for (int k = 0; k < aCentreColors.length; k++) {
            if (anIsCentreVisible[k]) {
                tmpClusters.add(k);
            }
        }
        this.clusterOfCentre = tmpClusters.stream().mapToInt(Integer::intValue).toArray();
        double[][] tmpCentrePositions = new double[this.clusterOfCentre.length][];
        double[] tmpCentreRadii = new double[this.clusterOfCentre.length];
        Color[] tmpCentreColors = new Color[this.clusterOfCentre.length];
        for (int c = 0; c < this.clusterOfCentre.length; c++) {
            tmpCentrePositions[c] = this.toScene(this.pcaResult.clusterVectorScores()[this.clusterOfCentre[c]]);
            tmpCentreRadii[c] = Pca3dView.CENTRE_RADIUS;
            tmpCentreColors[c] = aCentreColors[this.clusterOfCentre[c]];
        }
        Pca3dView.fillMeshView(this.centreMeshView, tmpCentrePositions, tmpCentreRadii, tmpCentreColors);
    }
    //
    /**
     * Shows the given loading vectors as arrows from the origin of the components to their tips, labelled with the
     * given names; replaces the vectors shown before. The tips are given in score units of the first three components,
     * scaled by the caller so that their lengths can be compared.
     *
     * @param aTipScoresArray tip of every vector in score units of PC1 to PC3
     * @param aNameArray name of every vector
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if the arrays have different lengths
     */
    public void setLoadingVectors(double[][] aTipScoresArray, String[] aNameArray) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aTipScoresArray, "Given vector tips are null.");
        Objects.requireNonNull(aNameArray, "Given vector names are null.");
        if (aTipScoresArray.length != aNameArray.length) {
            throw new IllegalArgumentException("The number of vector tips and names differ.");
        }
        this.loadingVectorGroup.getChildren().clear();
        //the labels of the vectors shown before follow the labels of the axes
        while (this.labelAnchors.size() > this.numberOfAxisLabels) {
            this.labelAnchors.removeLast();
            this.labelOverlay.getChildren().removeLast();
        }
        PhongMaterial tmpMaterial = new PhongMaterial(Pca3dView.LOADING_VECTOR_COLOR);
        Point3D tmpCylinderAxis = new Point3D(0.0, 1.0, 0.0);
        for (int v = 0; v < aTipScoresArray.length; v++) {
            double[] tmpTip = this.toScene(aTipScoresArray[v]);
            Point3D tmpDirection = new Point3D(tmpTip[0], tmpTip[1], tmpTip[2]);
            double tmpLength = tmpDirection.magnitude();
            if (!(tmpLength > 0.0)) {
                continue;
            }
            //the cylinder lies centred on the y-axis; it is rotated onto the direction and moved to the middle of the vector
            Cylinder tmpShaft = new Cylinder(Pca3dView.LOADING_VECTOR_RADIUS, tmpLength, 12);
            tmpShaft.setMaterial(tmpMaterial);
            Point3D tmpUnit = tmpDirection.normalize();
            Point3D tmpRotationAxis = tmpCylinderAxis.crossProduct(tmpUnit);
            double tmpAngle = Math.toDegrees(Math.acos(Math.max(-1.0, Math.min(1.0, tmpCylinderAxis.dotProduct(tmpUnit)))));
            tmpShaft.getTransforms().add(new Translate(tmpTip[0] / 2.0, tmpTip[1] / 2.0, tmpTip[2] / 2.0));
            if (tmpRotationAxis.magnitude() > 1.0e-9) {
                tmpShaft.getTransforms().add(new Rotate(tmpAngle, tmpRotationAxis));
            } else if (tmpUnit.getY() < 0.0) {
                //exactly opposite to the cylinder axis
                tmpShaft.getTransforms().add(new Rotate(180.0, Rotate.X_AXIS));
            }
            Sphere tmpHead = new Sphere(Pca3dView.LOADING_VECTOR_TIP_RADIUS, 12);
            tmpHead.setMaterial(tmpMaterial);
            tmpHead.setTranslateX(tmpTip[0]);
            tmpHead.setTranslateY(tmpTip[1]);
            tmpHead.setTranslateZ(tmpTip[2]);
            //the name a little beyond the tip
            double tmpLabelFactor = (tmpLength + 2.0 * Pca3dView.LOADING_VECTOR_TIP_RADIUS) / tmpLength;
            this.loadingVectorGroup.getChildren().addAll(tmpShaft, tmpHead, this.addLabel(new double[] {tmpTip[0] * tmpLabelFactor,
                    tmpTip[1] * tmpLabelFactor, tmpTip[2] * tmpLabelFactor}, aNameArray[v], Pca3dView.LOADING_VECTOR_NAME_STYLE));
        }
        this.updateLabelPositions();
    }
    //
    /**
     * Returns the half edge length of the cube the scores are scaled into, in score units, i.e. the largest absolute
     * score of the first three components.
     *
     * @return largest absolute score of PC1 to PC3
     */
    public double getMaximumScore() {
        return Pca3dView.HALF_SIZE / this.scale;
    }
    //
    /**
     * Sets the rotation and the zoom back to the start.
     */
    public void resetView() {
        this.rotateX.setAngle(Pca3dView.START_ANGLE_X);
        this.rotateY.setAngle(Pca3dView.START_ANGLE_Y);
        this.camera.setTranslateZ(-Pca3dView.CAMERA_DISTANCE);
    }
    //
    /**
     * Sets the actions called for molecules under the mouse.
     *
     * @param aHoverAction called with the row index of the hovered molecule, -1 if none
     * @param aClickAction called with the row index of a clicked molecule
     * @param aDoubleClickAction called with the row index of a double-clicked molecule
     */
    public void setMoleculeActions(IntConsumer aHoverAction, IntConsumer aClickAction, IntConsumer aDoubleClickAction) {
        this.moleculeHoverAction = aHoverAction;
        this.moleculeClickAction = aClickAction;
        this.moleculeDoubleClickAction = aDoubleClickAction;
    }
    //
    /**
     * Sets the actions called for cluster centres under the mouse.
     *
     * @param aHoverAction called with the cluster index of the hovered centre
     * @param aDoubleClickAction called with the cluster index of a double-clicked centre
     */
    public void setCentreActions(IntConsumer aHoverAction, IntConsumer aDoubleClickAction) {
        this.centreHoverAction = aHoverAction;
        this.centreDoubleClickAction = aDoubleClickAction;
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Converts the first three scores into scene coordinates: PC1 to the right, PC2 upwards, PC3 towards the viewer.
     *
     * @param aScores scores
     * @return scene coordinates x, y, z
     */
    private double[] toScene(double[] aScores) {
        return new double[] {aScores[0] * this.scale, -aScores[1] * this.scale, -aScores[2] * this.scale};
    }
    //
    /**
     * Creates the three axes through the origin with tick marks at round values; their labels are added to the overlay.
     *
     * @param aMaximumScore largest absolute score of the first three components, the half length of the axes in score
     *                      units
     * @param anAxisNameArray names of the three axes
     * @return axis nodes
     */
    private List<Node> createAxes(double aMaximumScore, String[] anAxisNameArray) {
        List<Node> tmpNodes = new ArrayList<>();
        PhongMaterial tmpMaterial = new PhongMaterial(Color.GREY);
        double tmpLength = 2.2 * Pca3dView.HALF_SIZE;
        for (int tmpComponent = 0; tmpComponent < 3; tmpComponent++) {
            //the cylinder lies along the y-axis and is rotated onto the x- and z-axis
            Cylinder tmpAxis = new Cylinder(Pca3dView.AXIS_RADIUS, tmpLength, 24);
            tmpAxis.setMaterial(tmpMaterial);
            if (tmpComponent == 0) {
                tmpAxis.getTransforms().add(new Rotate(90.0, Rotate.Z_AXIS));
            } else if (tmpComponent == 2) {
                tmpAxis.getTransforms().add(new Rotate(90.0, Rotate.X_AXIS));
            }
            tmpNodes.add(tmpAxis);
            //name of the axis at its positive end
            double[] tmpEnd = this.toScene(Pca3dView.unitScores(tmpComponent, 1.12 * Pca3dView.HALF_SIZE / this.scale));
            tmpNodes.add(this.addLabel(tmpEnd, anAxisNameArray[tmpComponent], Pca3dView.AXIS_NAME_STYLE));
            if (!(aMaximumScore > 0.0)) {
                continue;
            }
            //tick marks at multiples of a round tick unit, except the origin
            double tmpTickUnit = Pca3dView.getNiceTickUnit(aMaximumScore / Pca3dView.TICKS_PER_HALF_AXIS);
            int tmpDecimals = Math.max(0, BigDecimal.valueOf(tmpTickUnit).stripTrailingZeros().scale());
            long tmpLast = (long) Math.floor(aMaximumScore / tmpTickUnit + 1.0e-9);
            PhongMaterial tmpTickMaterial = new PhongMaterial(Color.DIMGREY);
            for (long t = -tmpLast; t <= tmpLast; t++) {
                if (t == 0) {
                    continue;
                }
                double tmpValue = t * tmpTickUnit;
                double[] tmpPosition = this.toScene(Pca3dView.unitScores(tmpComponent, tmpValue));
                Sphere tmpTick = new Sphere(Pca3dView.TICK_RADIUS, 8);
                tmpTick.setMaterial(tmpTickMaterial);
                tmpTick.setTranslateX(tmpPosition[0]);
                tmpTick.setTranslateY(tmpPosition[1]);
                tmpTick.setTranslateZ(tmpPosition[2]);
                tmpNodes.add(tmpTick);
                tmpNodes.add(this.addLabel(tmpPosition, String.format(Locale.getDefault(), "%." + tmpDecimals + "f", tmpValue),
                        Pca3dView.TICK_VALUE_STYLE));
            }
        }
        return tmpNodes;
    }
    //
    /**
     * Adds a sharp label to the overlay and returns the invisible anchor at its position in the 3D scene.
     *
     * @param aPosition scene position of the label
     * @param aText text of the label
     * @param aStyle style of the label, e.g. {@link #AXIS_NAME_STYLE}
     * @return anchor node, to be added to the 3D scene
     */
    private Node addLabel(double[] aPosition, String aText, String aStyle) {
        Group tmpAnchor = new Group();
        tmpAnchor.setTranslateX(aPosition[0]);
        tmpAnchor.setTranslateY(aPosition[1]);
        tmpAnchor.setTranslateZ(aPosition[2]);
        Label tmpLabel = new Label(aText);
        tmpLabel.setStyle(aStyle);
        this.labelAnchors.add(tmpAnchor);
        this.labelOverlay.getChildren().add(tmpLabel);
        return tmpAnchor;
    }
    //
    /**
     * Moves every label of the overlay to the projected position of its anchor; called after every rotation, zoom and
     * layout.
     */
    private void updateLabelPositions() {
        if (this.getScene() == null) {
            return;
        }
        for (int i = 0; i < this.labelAnchors.size(); i++) {
            Node tmpLabel = this.labelOverlay.getChildren().get(i);
            //the projection of the sub scene camera is included with rootScene = true
            Point3D tmpScenePoint = this.labelAnchors.get(i).localToScene(Point3D.ZERO, true);
            Point2D tmpLocal = this.labelOverlay.sceneToLocal(tmpScenePoint.getX(), tmpScenePoint.getY());
            if (tmpLocal == null) {
                continue;
            }
            Bounds tmpBounds = tmpLabel.getLayoutBounds();
            //right of the anchor
            tmpLabel.relocate(tmpLocal.getX() + 4.0, tmpLocal.getY() - tmpBounds.getHeight() / 2.0);
        }
    }
    //
    /**
     * Lays out the sub scene and the overlay and moves the labels to their projected positions.
     */
    @Override
    protected void layoutChildren() {
        super.layoutChildren();
        this.labelOverlay.resizeRelocate(0.0, 0.0, this.getWidth(), this.getHeight());
        this.updateLabelPositions();
    }
    //
    /**
     * Passes the molecule or centre under the mouse to the hover actions.
     *
     * @param anEvent mouse event
     */
    private void handleMouseMoved(MouseEvent anEvent) {
        PickResult tmpPick = anEvent.getPickResult();
        int tmpFace = tmpPick == null ? -1 : tmpPick.getIntersectedFace();
        if (tmpPick != null && tmpPick.getIntersectedNode() == this.moleculeMeshView && tmpFace >= 0) {
            if (this.moleculeHoverAction != null) {
                this.moleculeHoverAction.accept(this.rowOfPoint[tmpFace / Pca3dView.FACES_PER_POINT]);
            }
        } else if (tmpPick != null && tmpPick.getIntersectedNode() == this.centreMeshView && tmpFace >= 0) {
            if (this.centreHoverAction != null) {
                this.centreHoverAction.accept(this.clusterOfCentre[tmpFace / Pca3dView.FACES_PER_POINT]);
            }
        } else if (this.moleculeHoverAction != null) {
            this.moleculeHoverAction.accept(-1);
        }
    }
    //
    /**
     * Passes a clicked molecule or centre to the click actions; clicks at the end of a drag are ignored.
     *
     * @param anEvent mouse event
     */
    private void handleMouseClicked(MouseEvent anEvent) {
        if (anEvent.getButton() != MouseButton.PRIMARY || !anEvent.isStillSincePress()) {
            return;
        }
        PickResult tmpPick = anEvent.getPickResult();
        int tmpFace = tmpPick == null ? -1 : tmpPick.getIntersectedFace();
        if (tmpPick != null && tmpPick.getIntersectedNode() == this.moleculeMeshView && tmpFace >= 0) {
            int tmpRow = this.rowOfPoint[tmpFace / Pca3dView.FACES_PER_POINT];
            IntConsumer tmpAction = anEvent.getClickCount() == 2 ? this.moleculeDoubleClickAction : this.moleculeClickAction;
            if (tmpAction != null) {
                tmpAction.accept(tmpRow);
            }
        } else if (tmpPick != null && tmpPick.getIntersectedNode() == this.centreMeshView && tmpFace >= 0
                && anEvent.getClickCount() == 2 && this.centreDoubleClickAction != null) {
            this.centreDoubleClickAction.accept(this.clusterOfCentre[tmpFace / Pca3dView.FACES_PER_POINT]);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Returns scores that are 0 except for the given component.
     *
     * @param aComponent component index 0 to 2
     * @param aValue score of the component
     * @return scores of the three components
     */
    private static double[] unitScores(int aComponent, double aValue) {
        double[] tmpScores = new double[3];
        tmpScores[aComponent] = aValue;
        return tmpScores;
    }
    //
    /**
     * Rounds the given tick unit to 1, 2 or 5 times a power of ten, with one significant digit.
     *
     * @param aRawUnit raw tick unit, greater than 0
     * @return rounded tick unit
     */
    private static double getNiceTickUnit(double aRawUnit) {
        double tmpPower = Math.pow(10.0, Math.floor(Math.log10(aRawUnit)));
        double tmpFraction = aRawUnit / tmpPower;
        double tmpUnit = (tmpFraction < 1.5 ? 1.0 : tmpFraction < 3.5 ? 2.0 : tmpFraction < 7.5 ? 5.0 : 10.0) * tmpPower;
        return new BigDecimal(tmpUnit).round(new MathContext(1)).doubleValue();
    }
    //
    /**
     * Fills the given mesh view with one octahedron per position. The colours are taken from a palette texture with one
     * pixel per distinct colour, so that one mesh with one material can show any number of colours.
     *
     * @param aMeshView mesh view to fill
     * @param aPositions scene positions of the octahedra
     * @param aRadii radius of every octahedron
     * @param aColors colour of every octahedron
     */
    private static void fillMeshView(MeshView aMeshView, double[][] aPositions, double[] aRadii, Color[] aColors) {
        //palette of the distinct colours
        Map<Color, Integer> tmpPaletteIndexMap = new LinkedHashMap<>();
        for (Color tmpColor : aColors) {
            tmpPaletteIndexMap.putIfAbsent(tmpColor, tmpPaletteIndexMap.size());
        }
        int tmpPaletteSize = Math.max(1, tmpPaletteIndexMap.size());
        WritableImage tmpPalette = new WritableImage(tmpPaletteSize, 1);
        PixelWriter tmpWriter = tmpPalette.getPixelWriter();
        tmpPaletteIndexMap.forEach((aColor, anIndex) -> tmpWriter.setColor(anIndex, 0, aColor));
        TriangleMesh tmpMesh = new TriangleMesh();
        for (int c = 0; c < tmpPaletteSize; c++) {
            tmpMesh.getTexCoords().addAll((float) ((c + 0.5) / tmpPaletteSize), 0.5f);
        }
        float[] tmpPoints = new float[aPositions.length * 6 * 3];
        int[] tmpFaces = new int[aPositions.length * Pca3dView.FACES_PER_POINT * 6];
        //octahedron: vertices +x, -x, +y, -y, +z, -z and eight faces, outward oriented
        int[][] tmpOctahedronFaces = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
        for (int p = 0; p < aPositions.length; p++) {
            double[] tmpPosition = aPositions[p];
            double tmpRadius = aRadii[p];
            double[][] tmpVertices = {{tmpRadius, 0, 0}, {-tmpRadius, 0, 0}, {0, tmpRadius, 0}, {0, -tmpRadius, 0},
                    {0, 0, tmpRadius}, {0, 0, -tmpRadius}};
            for (int v = 0; v < 6; v++) {
                int tmpOffset = (p * 6 + v) * 3;
                tmpPoints[tmpOffset] = (float) (tmpPosition[0] + tmpVertices[v][0]);
                tmpPoints[tmpOffset + 1] = (float) (tmpPosition[1] + tmpVertices[v][1]);
                tmpPoints[tmpOffset + 2] = (float) (tmpPosition[2] + tmpVertices[v][2]);
            }
            int tmpTexture = tmpPaletteIndexMap.getOrDefault(aColors[p], 0);
            for (int f = 0; f < Pca3dView.FACES_PER_POINT; f++) {
                int tmpOffset = (p * Pca3dView.FACES_PER_POINT + f) * 6;
                for (int c = 0; c < 3; c++) {
                    tmpFaces[tmpOffset + 2 * c] = p * 6 + tmpOctahedronFaces[f][c];
                    tmpFaces[tmpOffset + 2 * c + 1] = tmpTexture;
                }
            }
        }
        tmpMesh.getPoints().setAll(tmpPoints);
        tmpMesh.getFaces().setAll(tmpFaces);
        PhongMaterial tmpMaterial = new PhongMaterial(Color.WHITE);
        tmpMaterial.setDiffuseMap(tmpPalette);
        aMeshView.setMesh(tmpMesh);
        aMeshView.setMaterial(tmpMaterial);
    }
    //</editor-fold>
}
