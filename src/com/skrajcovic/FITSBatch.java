package com.skrajcovic;

import com.sun.corba.se.impl.encoding.OSFCodeSetRegistry;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSBatch {
    private HashMap<double[], FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<FITSObject[], SimpleRegression> regressions;

    public static final boolean DEBUG = true;

    public FITSBatch() {
        regressions = new HashMap<>();
        data = new HashMap<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

    public FITSBatch(HashMap<double[], FITSObject> data) {
        this.data = data;
    }

    public void doTheThing() {

        findRegressions();

        for (int i = 1; i < 101; i++) {
            double[] tmp = testRegressions(i);
            if (tmp != null) {
                System.out.format("Threshold: %d\nNumber of all points under threshold: %d\nReal points: %d\nSuccess rate: %f",
                        i, (int)tmp[0], (int)tmp[1], tmp[2]); System.out.println("%");
                System.out.println("---------------------------------");
            }
        }
//        testRegressions();
    }

    void findRegressions() {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getX(), obj1.getY());
                sr.addData(obj2.getX(), obj2.getY());
                regressions.put(new FITSObject[]{obj1, obj2}, sr);
            }
        }
    }

    double[] testRegressions(double threshold) {
//        double threshold = 70;
        for (Map.Entry<FITSObject[], SimpleRegression> regression : regressions.entrySet()) {
            Set<FITSObject> result = new HashSet<>();
            int real = 0;
            for (FITSObject fitsObject : data.values()) {
                double x = fitsObject.getX();
                double y = fitsObject.getY();
                double m = regression.getValue().getSlope();
                double c = regression.getValue().getIntercept();
                double b = (y > 0) ? 1 : -1;

                double distance = (-m * x + b * y - c) / Math.sqrt(Math.pow(m, 2) + Math.pow(b, 2));

                if (Math.abs(distance) <= threshold) {
                    if (fitsObject.isReal()) {
                        real++;
                    }
                    result.add(fitsObject);
                }
            }
            if (regression.getKey()[0].isReal() && regression.getKey()[1].isReal()) {
                return new double[]{result.size(), real, (real / (double) result.size())*100};
//                    System.out.println(Arrays.toString(regression.getKey()) + " -- " + result);
//                    System.out.println(result.size());
//                    System.out.println(real);
//                    System.out.println("Success rate: " + (real / (double) result.size()) * 100 + "%");
//                    System.out.println("-------------");
            }
        }
        return null;
    }

    public void mainDataInsert(double[] id, FITSObject object) {
        data.put(id, object);
    }

    public void firstSetInsert(FITSObject fitsObject) {
        fSet.add(fitsObject);
    }

    public void secondSetInsert(FITSObject fitsObject) {
        sSet.add(fitsObject);
    }

    public HashMap<double[], FITSObject> getDataStructure() {
        return data;
    }

    public String toString() {
//        for (double[] d : data.keySet()) {
//            System.out.println(String.valueOf(d[0]).replace(".", ",") + " " + String.valueOf(d[1]).replace(".", ","));
//        }
//        System.out.println();
        return this.data.toString();
    }
}
