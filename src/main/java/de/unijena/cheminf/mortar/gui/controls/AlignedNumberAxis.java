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

import javafx.scene.chart.ValueAxis;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Number axis whose tick marks always lie on multiples of the tick unit (e.g. 0.1, 0.2, 0.3) instead of starting at the
 * lower bound like the tick marks of {@link javafx.scene.chart.NumberAxis}, so that the labels stay round numbers also
 * after zooming and moving the view. The bounds and the tick unit are set by the user of the axis; auto-ranging is not
 * supported. The labels show as many decimals as the tick unit needs.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
public class AlignedNumberAxis extends ValueAxis<Number> {
    //<editor-fold desc="private variables" defaultstate="collapsed">
    /**
     * Distance between two major tick marks.
     */
    private double tickUnit;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor with the range [0, 1] and a tick unit of 0.1.
     */
    public AlignedNumberAxis() {
        super(0.0, 1.0);
        this.tickUnit = 0.1;
        this.setAutoRanging(false);
    }
    //</editor-fold>
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     * Returns the distance between two major tick marks.
     *
     * @return tick unit
     */
    public double getTickUnit() {
        return this.tickUnit;
    }
    //
    /**
     * Sets the distance between two major tick marks and lays out the axis again.
     *
     * @param aTickUnit tick unit, greater than 0
     */
    public void setTickUnit(double aTickUnit) {
        if (aTickUnit > 0.0) {
            this.tickUnit = aTickUnit;
            this.invalidateRange();
            this.requestAxisLayout();
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="overridden methods" defaultstate="collapsed">
    /**
     * Returns the multiples of the tick unit within the given range.
     *
     * @param aLength length of the axis in display units
     * @param aRange range as returned by {@link #getRange()}
     * @return tick values in ascending order
     */
    @Override
    protected List<Number> calculateTickValues(double aLength, Object aRange) {
        double[] tmpRange = (double[]) aRange;
        List<Number> tmpTicks = new ArrayList<>();
        if (!(this.tickUnit > 0.0) || !(tmpRange[1] > tmpRange[0])) {
            return tmpTicks;
        }
        long tmpFirst = (long) Math.ceil(tmpRange[0] / this.tickUnit - 1.0e-9);
        long tmpLast = (long) Math.floor(tmpRange[1] / this.tickUnit + 1.0e-9);
        //at most a few hundred ticks, also for an extreme ratio of range and tick unit
        for (long i = tmpFirst; i <= tmpLast && tmpTicks.size() < 500; i++) {
            //rounded to the decimals of the tick unit, so that no value like 0.30000000000000004 appears
            tmpTicks.add(BigDecimal.valueOf(i * this.tickUnit).setScale(this.getNumberOfDecimals() + 2, java.math.RoundingMode.HALF_UP).doubleValue());
        }
        return tmpTicks;
    }
    //
    /**
     * Returns the minor tick marks: four between two major tick marks.
     *
     * @return minor tick values
     */
    @Override
    protected List<Number> calculateMinorTickMarks() {
        List<Number> tmpMinorTicks = new ArrayList<>();
        int tmpCount = Math.max(1, this.getMinorTickCount());
        double tmpMinorUnit = this.tickUnit / tmpCount;
        double tmpStart = Math.ceil(this.getLowerBound() / tmpMinorUnit) * tmpMinorUnit;
        for (double tmpValue = tmpStart; tmpValue <= this.getUpperBound() && tmpMinorTicks.size() < 2500; tmpValue += tmpMinorUnit) {
            tmpMinorTicks.add(tmpValue);
        }
        return tmpMinorTicks;
    }
    //
    /**
     * Sets the given range: lower bound, upper bound and scale.
     *
     * @param aRange range as returned by {@link #getRange()}
     * @param anIsAnimating ignored, the range is set without animation
     */
    @Override
    protected void setRange(Object aRange, boolean anIsAnimating) {
        double[] tmpRange = (double[]) aRange;
        this.setLowerBound(tmpRange[0]);
        this.setUpperBound(tmpRange[1]);
        this.currentLowerBound.set(tmpRange[0]);
        this.setScale(tmpRange[2]);
    }
    //
    /**
     * Returns the current range: lower bound, upper bound and the scale for the current length of the axis.
     *
     * @return range as double array
     */
    @Override
    protected Object getRange() {
        double tmpLength = this.getSide() != null && this.getSide().isVertical() ? this.getHeight() : this.getWidth();
        return new double[] {this.getLowerBound(), this.getUpperBound(),
                this.calculateNewScale(tmpLength, this.getLowerBound(), this.getUpperBound())};
    }
    //
    /**
     * Returns the label of a tick mark with as many decimals as the tick unit needs.
     *
     * @param aValue tick value
     * @return label
     */
    @Override
    protected String getTickMarkLabel(Number aValue) {
        return String.format(Locale.getDefault(), "%." + this.getNumberOfDecimals() + "f", aValue.doubleValue());
    }
    //</editor-fold>
    //
    //<editor-fold desc="private methods" defaultstate="collapsed">
    /**
     * Returns the number of decimals the tick unit needs, e.g. 1 for 0.1 or 0.5, 2 for 0.05.
     *
     * @return number of decimals, at least 0
     */
    private int getNumberOfDecimals() {
        return Math.max(0, BigDecimal.valueOf(this.tickUnit).stripTrailingZeros().scale());
    }
    //</editor-fold>
}
