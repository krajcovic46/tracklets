package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class FITSBatch {
    private HashMap<double[], FITSObject> data;
    SimpleRegression regression;

    public FITSBatch() {
        data = new HashMap<>();
        regression = new SimpleRegression();
    }

    public FITSBatch(HashMap<double[], FITSObject> data) {
        this.data = data;
    }

    public void getPoints() {
        double m = regression.getSlope();
        double c = regression.getIntercept();

        System.out.println(regression.getRSquare());

        for (Map.Entry<double[], FITSObject> e : data.entrySet()) {
            System.out.println("y = " + e.getKey()[1]);
            System.out.println("x = " + m*e.getKey()[0] + c);
            if (e.getKey()[1] == m*e.getKey()[0] + c) {
                System.out.println("wut");
            }
            System.out.println("-----------------------------");
        }

        System.out.println(m + " " + c);
    }

    public HashMap<double[], FITSObject> getDataStructure() {
        return data;
    }

    public String toString() {
        for (double[] d : data.keySet()) {
            System.out.println(d[0] + " " + d[1]);
        }
        System.out.println();
        return this.data.toString();
    }
}
