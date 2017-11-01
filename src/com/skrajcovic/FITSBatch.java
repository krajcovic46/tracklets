package com.skrajcovic;

import eap.fits.*;
import eap.fitsbrowser.HDUDisplay;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

public class FITSBatch {
    private List<FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<ArrayList<FITSObject>, SimpleRegression> regressions;
    private Set<FitsData> images;

    public FITSBatch() {
        regressions = new HashMap<>();
        data = new ArrayList<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

    public void doTheThing() {

        findRegressions();
        for (int i = 100; i < 101; i++) {
            double[] tmp = fitPointsToRegressions(i);
            if (tmp != null) {
                System.out.format("Threshold: %d\nNumber of all points under threshold: %d\nReal points: %d\nSuccess rate: %f",
                        i, (int)tmp[0], (int)tmp[1], tmp[2]); System.out.println("%");
                System.out.println("---------------------------------");
            }
        }
    }

    void findRegressions() {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getX(), obj1.getY());
                sr.addData(obj2.getX(), obj2.getY());
                regressions.put(new ArrayList<>(Arrays.asList(obj1, obj2)), sr);
            }
        }
    }

    double[] fitPointsToRegressions(double threshold) {
        for (Map.Entry<ArrayList<FITSObject>, SimpleRegression> regression : regressions.entrySet()) {
            if (regression.getKey().get(0).isReal() && regression.getKey().get(1).isReal()) {
                FITSObject last = null;
                double lastSpeed = Double.MAX_VALUE;
                ArrayList<FITSObject> regressionPoints = regression.getKey();

                double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0));
                double baseHeading = regressionPoints.get(0).getHeading(regressionPoints.get(1));

                int real = 0;
                for (FITSObject fitsObject : data) {
                    if (fitsObject.isWithinLineThreshold(regression.getValue(), threshold)
                            && regressionPoints.get(regressionPoints.size() - 1).isWithinAngleThreshold(fitsObject, baseHeading, 20)
                            ) {

                        if (last != null && !fitsObject.getFileName().equals(last.getFileName()) && !regressionPoints.contains(last)) {
                            regressionPoints.add(last);
                            regression.getValue().addData(last.getX(), last.getY());

                            //cleanup
                            last = null;
                            lastSpeed = Double.MAX_VALUE;
                        }

                        double currentSpeed = regressionPoints.get(regressionPoints.size() - 1).calculateSpeed(fitsObject);
                        if (Math.abs(averageCombinedSpeed - currentSpeed) < Math.abs(averageCombinedSpeed - lastSpeed)) {
                            last = fitsObject;
                            lastSpeed = currentSpeed;
                        }
                    }
                }
                // toto sa stane ak uplne posledny prvok vo forcykle zapada - musi sa pridat
                if (last != null && !regressionPoints.get(regressionPoints.size() - 1).getFileName().equals(last.getFileName())) {
                    regressionPoints.add(last);
                    regression.getValue().addData(last.getX(), last.getY());
                }
                if (regressionPoints.get(0).isReal() && regressionPoints.get(1).isReal()) {
                    System.out.println(regressionPoints);
                    for (FITSObject obj : regressionPoints) {
                        if (obj.isReal()) {
                            real++;
                        }
                    }
                    return new double[]{regressionPoints.size(), real, (real / (double) regressionPoints.size()) * 100};
                }
            }
        }
        return null;
    }

    public void mainDataInsert(FITSObject object) {
        data.add(object);
    }

    public void firstSetInsert(FITSObject fitsObject) {
        fSet.add(fitsObject);
    }

    public void secondSetInsert(FITSObject fitsObject) {
        sSet.add(fitsObject);
    }

    public List<FITSObject> getDataStructure() {
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
