package com.skrajcovic.algorithms;

import com.skrajcovic.FITSBatch;
import com.skrajcovic.FITSObject;
import com.skrajcovic.FITSTracklet;
import com.skrajcovic.datastructures.Type;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSLinearRegression {

    private static final double distanceThreshold = 50;
    private static final double angleThreshold = 50;
    private static final double speedThreshold = 100;

    public static void perform(FITSBatch batch) {
        FITSLinearRegression.findInitialRegressions((HashSet<FITSObject>) batch.getFirstSet(),
                (HashSet<FITSObject>) batch.getSecondSet(), batch.getRegressions());

        FITSLinearRegression.fitPointsToRegressions2(distanceThreshold, batch);

        System.out.println(FITSBatch.tracklets);

//        System.out.println(batch.getRegressionResults());

    }

    private static void findInitialRegressions(HashSet<FITSObject> fSet, HashSet<FITSObject> sSet, Map<ArrayList<FITSObject>, SimpleRegression> regressions) {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getxComponent(), obj1.getyComponent());
                sr.addData(obj2.getyComponent(), obj2.getyComponent());
//                if (obj1.getType() == Type.H && obj2.getType() == Type.H) {
//                    System.out.println(obj1.getFileName());
//                    System.out.println(obj1.getxComponent());
//                    System.out.println(obj1.getyComponent());
//
//                    System.out.println(obj2.getFileName());
//                    System.out.println(obj2.getxComponent());
//                    System.out.println(obj2.getyComponent());
//
//                    System.out.println("slope: " + sr.getSlope());
//                    System.out.println("intercept: " + sr.getIntercept());
//                }
                regressions.put(new ArrayList<>(Arrays.asList(obj1, obj2)), sr);
            }
        }
    }

    private static void fitPointsToRegressions2(double threshold, FITSBatch batch) {
        Map<ArrayList<FITSObject>, SimpleRegression> regressions = batch.getRegressions();
        ArrayList<FITSObject> data = (ArrayList<FITSObject>) batch.getMainData();

        String name = "";

        for (Map.Entry<ArrayList<FITSObject>, SimpleRegression> mapDataRegression : regressions.entrySet()) {
            FITSTracklet tracklet = new FITSTracklet();
            SimpleRegression regression = mapDataRegression.getValue();
            ArrayList<FITSObject> regressionPoints = mapDataRegression.getKey();

            FITSObject firstObject = regressionPoints.get(0);
            FITSObject secondObject = regressionPoints.get(1);
            FITSObject lastObject = regressionPoints.get(regressionPoints.size() - 1);

            double baselineHeading = secondObject.calculateHeading(firstObject);
            double baselineSpeed = firstObject.calculateSpeed(secondObject);

            tracklet.baselineHeading = baselineHeading;
            tracklet.baselineSpeed = baselineSpeed;

            ArrayList<Map<FITSObject, Double>> firstList = new ArrayList<>();
            Map<FITSObject, Double> firstMap = new HashMap<>(); firstMap.put(firstObject, 0D);
            firstList.add(firstMap);
            tracklet.objects.add(firstList);

            ArrayList<Map<FITSObject, Double>> secondList = new ArrayList<>();
            Map<FITSObject, Double> secondMap = new HashMap<>(); secondMap.put(firstObject, baselineHeading + baselineSpeed);
            secondList.add(secondMap);
            tracklet.objects.add(secondList);

            ArrayList<Map<FITSObject, Double>> nextSet = new ArrayList<>();

            name = secondObject.getFileName();
            int numberOfObjects = 0;

            for (FITSObject fitsObject : data) {
                if (!fitsObject.getFileName().equals(name)) {
                    numberOfObjects = 0;
                    name = fitsObject.getFileName();
                }
                double distanceToLine = fitsObject.calculateDistanceToLine(regression);
                if (fitsObject.isWithinLineThreshold(distanceToLine, threshold)) {
                    double angle = fitsObject.calculateHeading(lastObject);
                    if (fitsObject.isWithinAngleThreshold(angle, baselineHeading, angleThreshold)) {
                        double speed = fitsObject.calculateSpeed(lastObject);
                        if (fitsObject.isWithinSpeedThreshold(speed, baselineSpeed, speedThreshold)) {
                            System.out.println("trigger: " + fitsObject);

                            HashMap<FITSObject, Double> temp = new HashMap<>();
                            temp.put(fitsObject, distanceToLine + angle + speed);

                            nextSet.add(temp);
                        }
                    }
                }
                numberOfObjects++;
//                System.out.print(numberOfObjects + " = ");
//                System.out.println(FITSBatch.objectsCount.get(fitsObject.getFileName()));
                if (numberOfObjects == FITSBatch.objectsCount.get(fitsObject.getFileName())) {
                    tracklet.objects.add(nextSet);
                    tracklet.sortLast();
                    nextSet = new ArrayList<>();
                }
            }
            FITSBatch.tracklets.add(tracklet);
        }
    }

//    private static void fitPointsToRegressions(double threshold, FITSBatch batch) {
//        Map<ArrayList<FITSObject>, SimpleRegression> regressions = batch.getRegressions();
//        ArrayList<FITSObject> data = (ArrayList<FITSObject>) batch.getMainData();
//
//        for (Map.Entry<ArrayList<FITSObject>, SimpleRegression> regression : regressions.entrySet()) {
//            FITSObject last = null;
//            double lastSpeed = Double.MAX_VALUE;
//            ArrayList<FITSObject> regressionPoints = regression.getKey();
//
//            if (regressionPoints.get(0).getType() == Type.H && regressionPoints.get(1).getType() == Type.H) {
//
//                double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0));
//                double baseHeading = regressionPoints.get(0).calculateHeading(regressionPoints.get(1));
//
//                for (FITSObject fitsObject : data) {
//                    if (fitsObject.isWithinLineThreshold(regression.getValue(), threshold)
//                            && regressionPoints.get(regressionPoints.size() - 1).isWithinAngleThreshold(fitsObject, baseHeading, angleThreshold)) {
//
////                        System.out.println(fitsObject.getFileName());
////                        System.out.println(fitsObject.getType());
////                        System.out.println(fitsObject.getX());
////                        System.out.println(fitsObject.getY());
////                        System.out.println("-----------------");
//
//                        if (last != null && !fitsObject.getFileName().equals(last.getFileName()) && !regressionPoints.contains(last)) {
//                            regressionPoints.add(last);
//                            regression.getValue().addData(last.getxComponent(), last.getyComponent());
//
//                            //cleanup
//                            last = null;
//                            lastSpeed = Double.MAX_VALUE;
//                        }
//
//                        double currentSpeed = regressionPoints.get(regressionPoints.size() - 1).calculateSpeed(fitsObject);
//                        if (Math.abs(averageCombinedSpeed - currentSpeed) < Math.abs(averageCombinedSpeed - lastSpeed)) {
//                            last = fitsObject;
//                            lastSpeed = currentSpeed;
//                        }
//                    }
//                }
//                // toto sa stane ak uplne posledny prvok vo forcykle zapada - musi sa pridat
//                if (last != null && !regressionPoints.get(regressionPoints.size() - 1).getFileName().equals(last.getFileName())) {
//                    regressionPoints.add(last);
////                    regression.getValue().addData(last.getxComponent(), last.getyComponent());
//                }
//                batch.setRegressionResults(regressionPoints);
//            }
//        }
//    }
}
