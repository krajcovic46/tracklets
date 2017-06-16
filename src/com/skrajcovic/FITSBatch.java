package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSBatch {
    private List<FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<ArrayList<FITSObject>, SimpleRegression> regressions;

    public static final boolean DEBUG = false;

    public FITSBatch() {
        regressions = new HashMap<>();
        data = new ArrayList<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

    public void doTheThing() {

        findRegressions();
        for (int i = 1; i < 101; i++) {
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
            Set<FITSObject> result = new HashSet<>();
            FITSObject last = null;
            double lastSpeed = Double.MAX_VALUE;
            ArrayList<FITSObject> regressionPoints = regression.getKey();

            double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0),
                    regressionPoints.get(1).calculateDeltaTime(regressionPoints.get(0)));

            int real = 0;
            for (FITSObject fitsObject : data) {

                double deltaTime;

                if (fitsObject.isWithinLineThreshold(regression.getValue(), threshold)) {
                    deltaTime = regressionPoints.get(regressionPoints.size() - 1).calculateDeltaTime(fitsObject);
                    double currentSpeed = regressionPoints.get(1).calculateSpeed(fitsObject, deltaTime);
                    if (Math.abs(averageCombinedSpeed - currentSpeed) < Math.abs(averageCombinedSpeed - lastSpeed)) {
                        last = fitsObject;
                        lastSpeed = currentSpeed;
                    }

                    if (last != null && !fitsObject.getName().equals(last.getName()) &&
                            !regressionPoints.get(regressionPoints.size() - 1).getName().equals(last.getName())) {
                        regressionPoints.add(last);
                        regression.getValue().addData(last.getX(), last.getY());
                        averageCombinedSpeed = (averageCombinedSpeed + lastSpeed) / 2;

                        //cleanup
                        last = null;
                        lastSpeed = Double.MAX_VALUE;
                    }

                    result.add(fitsObject);

                    if (fitsObject.isReal()) {
                        real++;
                    }
                }
            }
            // toto sa stane ak uplne posledny prvok vo forcykle zapada - musi sa pridat
            if (last != null && !regressionPoints.get(regressionPoints.size() - 1).getName().equals(last.getName())) {
                regressionPoints.add(last);
                regression.getValue().addData(last.getX(), last.getY());
            }
            if (regressionPoints.get(0).isReal() && regressionPoints.get(1).isReal()) {
                return new double[]{result.size(), real, (real / (double) result.size())*100};
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
