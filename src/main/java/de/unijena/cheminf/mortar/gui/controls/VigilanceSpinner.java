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

package de.unijena.cheminf.mortar.gui.controls;

import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.clustering.ClusteringService;

import javafx.beans.property.DoubleProperty;
import javafx.scene.control.Spinner;
import javafx.scene.control.SpinnerValueFactory;
import javafx.util.StringConverter;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Objects;

/**
 * Editable spinner for the vigilance of a {@link ClusteringService}. It is bound to the vigilance setting of the service,
 * so that all instances (in the bottom bars of the fragments, items and clustering result tabs) and the clustering
 * settings show the same value. The arrows change the value in steps of 0.05; a typed value is taken over with enter or
 * when the field loses the focus, an invalid value is replaced by the current one. Both decimal point and comma are
 * accepted.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class VigilanceSpinner extends Spinner<Double> {
    //<editor-fold desc="private static final class constants" defaultstate="collapsed">
    /**
     * Step size of the spinner arrows.
     */
    private static final double STEP_SIZE = 0.05;
    /**
     * Number of decimals the value is rounded to, as in the vigilance scan.
     */
    private static final int NUMBER_OF_DECIMALS = 4;
    /**
     * Preferred width of the spinner.
     */
    private static final double PREF_WIDTH = 85.0;
    //</editor-fold>
    //
    //<editor-fold desc="private final variables" defaultstate="collapsed">
    /**
     * Vigilance setting of the clustering service.
     */
    private final DoubleProperty vigilanceSetting;
    //</editor-fold>
    //
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * True while the value is synchronised with the setting, to avoid a loop of the two listeners.
     */
    private boolean isUpdating;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aClusteringService clustering service whose vigilance is shown and set
     * @throws NullPointerException if the service is null
     */
    public VigilanceSpinner(ClusteringService aClusteringService) throws NullPointerException {
        super();
        this.vigilanceSetting = Objects.requireNonNull(aClusteringService, "Given clustering service is null.").vigilanceProperty();
        SpinnerValueFactory.DoubleSpinnerValueFactory tmpValueFactory = new SpinnerValueFactory.DoubleSpinnerValueFactory(
                ClusteringService.OPEN_INTERVAL_MINIMUM, ClusteringService.OPEN_INTERVAL_MAXIMUM,
                this.vigilanceSetting.get(), VigilanceSpinner.STEP_SIZE);
        tmpValueFactory.setConverter(this.createConverter());
        this.setValueFactory(tmpValueFactory);
        this.setEditable(true);
        this.setPrefWidth(VigilanceSpinner.PREF_WIDTH);
        this.setTooltip(GuiUtil.createTooltip(Message.get("ClusteringService.vigilanceSetting.tooltip")));
        //value of the spinner -> setting; the arrows can produce values like 0.6000000000000001, so it is rounded
        tmpValueFactory.valueProperty().addListener((anObservable, anOldValue, aNewValue) -> {
            if (this.isUpdating || aNewValue == null) {
                return;
            }
            double tmpRoundedValue = VigilanceSpinner.round(aNewValue);
            this.isUpdating = true;
            try {
                tmpValueFactory.setValue(tmpRoundedValue);
                this.vigilanceSetting.set(tmpRoundedValue);
            } catch (IllegalArgumentException anException) {
                //already shown by the setting, the value of the setting is restored
                tmpValueFactory.setValue(this.vigilanceSetting.get());
            } finally {
                this.isUpdating = false;
            }
        });
        //setting -> value of the spinner, e.g. after choosing a vigilance in a scan or changing it in the settings
        this.vigilanceSetting.addListener((anObservable, anOldValue, aNewValue) -> {
            if (this.isUpdating) {
                return;
            }
            this.isUpdating = true;
            try {
                tmpValueFactory.setValue(aNewValue.doubleValue());
            } finally {
                this.isUpdating = false;
            }
        });
        //a typed value is only taken over with enter by default, here also when the field loses the focus
        this.getEditor().focusedProperty().addListener((anObservable, anOldValue, aNewValue) -> {
            if (!aNewValue) {
                this.commitEditorText();
            }
        });
        this.getEditor().setOnAction(anEvent -> this.commitEditorText());
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Takes over the text of the editor as value; a text that is no number or outside the allowed range is replaced by
     * the current value.
     */
    private void commitEditorText() {
        SpinnerValueFactory<Double> tmpValueFactory = this.getValueFactory();
        double tmpValue = tmpValueFactory.getConverter().fromString(this.getEditor().getText());
        if (tmpValue >= ClusteringService.OPEN_INTERVAL_MINIMUM && tmpValue <= ClusteringService.OPEN_INTERVAL_MAXIMUM) {
            tmpValueFactory.setValue(tmpValue);
        }
        //shows the (possibly rounded or restored) value
        this.getEditor().setText(tmpValueFactory.getConverter().toString(tmpValueFactory.getValue()));
    }
    //
    /**
     * Creates the converter between value and text: up to four decimals without trailing zeros, decimal point or comma
     * when parsing.
     *
     * @return converter, which returns the current vigilance for a text that is no number, so that the spinner never
     *         gets an empty value
     */
    private StringConverter<Double> createConverter() {
        return new StringConverter<>() {
            @Override
            public String toString(Double aValue) {
                if (aValue == null) {
                    return "";
                }
                return BigDecimal.valueOf(VigilanceSpinner.round(aValue)).stripTrailingZeros().toPlainString();
            }
            //
            @Override
            public Double fromString(String aText) {
                if (aText == null || aText.isBlank()) {
                    return VigilanceSpinner.this.vigilanceSetting.get();
                }
                try {
                    return VigilanceSpinner.round(Double.parseDouble(aText.trim().replace(',', '.')));
                } catch (NumberFormatException anException) {
                    return VigilanceSpinner.this.vigilanceSetting.get();
                }
            }
        };
    }
    //</editor-fold>
    //
    //<editor-fold desc="private static methods" defaultstate="collapsed">
    /**
     * Rounds the given value to the number of decimals of the spinner.
     *
     * @param aValue value
     * @return rounded value
     */
    private static double round(double aValue) {
        return BigDecimal.valueOf(aValue).setScale(VigilanceSpinner.NUMBER_OF_DECIMALS, RoundingMode.HALF_UP).doubleValue();
    }
    //</editor-fold>
}
