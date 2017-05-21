package com.skrajcovic;

import com.sun.corba.se.impl.encoding.OSFCodeSetRegistry;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import sun.java2d.pipe.SpanShapeRenderer;

import java.util.*;

public class FITSBatch {
    private HashMap<double[], FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<FITSObject[], SimpleRegression> regressions;

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

        testRegressions();
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

    void testRegressions() {
        double threshold = 15;
        for (Map.Entry<FITSObject[], SimpleRegression> regression : regressions.entrySet()) {
            Set<FITSObject> result = new HashSet<>();
            int real = 0;
            for (FITSObject fitsObject : data.values()) {
//                if (fitsObject != regression.getKey()[0] && fitsObject != regression.getKey()[1]) {
                    double x = fitsObject.getX();
                    double y = fitsObject.getY();
                    double m = regression.getValue().getSlope();
                    double c = regression.getValue().getIntercept();

                    if (y <= (m * x + c) + threshold && y >= (m * x + c) - threshold) {
                        if (x == 599.61 || x == 598.63 || x == 599.79 || x == 600.07 || x == 599.98 || x == 600.62) {
                            real++;
                        }
                        result.add(fitsObject);
                    }
                }
                if (result.size() > 0) {
                    System.out.println(result);
                    System.out.println(result.size());
                    System.out.println(real);
                    System.out.println("-------------");
                }
            }
//        }
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
