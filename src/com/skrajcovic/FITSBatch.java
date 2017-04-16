package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.Arrays;
import java.util.HashMap;

public class FITSBatch {
    private HashMap<double[], FITSObject> data;

    public FITSBatch(HashMap<double[], FITSObject> data) {
        this.data = data;
    }

    public void doRegression() {
        SimpleRegression regression = new SimpleRegression();
        for (FITSObject object : data.values()) {
            regression.addData(object.getX(), object.getY());
        }
        System.out.println(regression.getSlope());
        System.out.println(regression.getIntercept());
    }

    public String toString() {
        for (double[] d : data.keySet()) {
            System.out.print(Arrays.toString(d));
        }
        System.out.println();
        return this.data.toString();
    }
}
