package com.skrajcovic.algorithms;

import com.skrajcovic.SDTBatch;
import com.skrajcovic.SDTObject;
import com.skrajcovic.SDTTracklet;
import com.skrajcovic.datastructures.Type;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class SDTLinearRegression {

    private static final double distanceThreshold = 50;
    private static final double angleThreshold = 50;
    private static final double speedThreshold = 100000;

    public static void perform(SDTBatch batch) {
        SDTLinearRegression.findInitialRegressions((HashSet<SDTObject>) batch.getFirstSet(),
                (HashSet<SDTObject>) batch.getSecondSet(), batch.getRegressions());

        SDTLinearRegression.fitPointsToRegressions2(distanceThreshold, batch);

        System.out.println(SDTBatch.tracklets);
        batch.filterOutEmptyTracklets();
        System.out.println("--------------");
        System.out.println(SDTBatch.tracklets);

//        System.out.println(batch.getRegressionResults());

    }

    private static void findInitialRegressions(HashSet<SDTObject> fSet, HashSet<SDTObject> sSet, Map<ArrayList<SDTObject>, SimpleRegression> regressions) {
        for (SDTObject obj1 : fSet) {
            for (SDTObject obj2 : sSet) {
                if (obj1.isUnidentified() && obj2.isUnidentified()) {
                    SimpleRegression sr = new SimpleRegression();
                    sr.addData(obj1.getxComponent(), obj1.getyComponent());
                    sr.addData(obj2.getxComponent(), obj2.getyComponent());
                    regressions.put(new ArrayList<>(Arrays.asList(obj1, obj2)), sr);
                }
            }
        }
    }

    private static void fitPointsToRegressions2(double threshold, SDTBatch batch) {
        Map<ArrayList<SDTObject>, SimpleRegression> regressions = batch.getRegressions();
        ArrayList<SDTObject> data = (ArrayList<SDTObject>) batch.getMainData();

        String name = "";

        for (Map.Entry<ArrayList<SDTObject>, SimpleRegression> mapDataRegression : regressions.entrySet()) {
            SDTTracklet tracklet = new SDTTracklet();
            SimpleRegression regression = mapDataRegression.getValue();
            ArrayList<SDTObject> regressionPoints = mapDataRegression.getKey();

            SDTObject firstObject = regressionPoints.get(0);
            SDTObject secondObject = regressionPoints.get(1);
            SDTObject lastObject = regressionPoints.get(regressionPoints.size() - 1);

            double baselineHeading = secondObject.calculateHeading(firstObject);
            double baselineSpeed = firstObject.calculateSpeed(secondObject);

            ArrayList<Map<SDTObject, Double>> firstList = new ArrayList<>();
            Map<SDTObject, Double> firstMap = new HashMap<>(); firstMap.put(firstObject, 0D);
            firstList.add(firstMap);
            tracklet.objects.add(firstList);

            ArrayList<Map<SDTObject, Double>> secondList = new ArrayList<>();
            Map<SDTObject, Double> secondMap = new HashMap<>(); secondMap.put(secondObject, baselineHeading + baselineSpeed);
            secondList.add(secondMap);
            tracklet.objects.add(secondList);

            tracklet.valueToCompareWith = baselineHeading + baselineSpeed;

            ArrayList<Map<SDTObject, Double>> nextSet = new ArrayList<>();

            name = secondObject.getFileName();
            int numberOfObjects = 0;

            for (SDTObject SDTObject : data) {
                if (!SDTObject.getFileName().equals(name)) {
                    numberOfObjects = 0;
                    name = SDTObject.getFileName();
                }
                double distanceToLine = SDTObject.calculateDistanceToLine(regression);
                if (SDTObject.isWithinLineThreshold(distanceToLine, threshold)) {
                    if (firstObject.getType() == Type.H && secondObject.getType() == Type.H) {
                        System.out.println("trigger line threshold: " + SDTObject);
                    }
                    double angle = SDTObject.calculateHeading(lastObject);
                    if (SDTObject.isWithinAngleThreshold(angle, baselineHeading, angleThreshold)) {
                        if (firstObject.getType() == Type.H && secondObject.getType() == Type.H) {
                            System.out.println("trigger angle threshold: " + SDTObject);
                        }
                        double speed = SDTObject.calculateSpeed(lastObject);
                        if (SDTObject.isWithinSpeedThreshold(speed, baselineSpeed, speedThreshold)) {
                            if (firstObject.getType() == Type.H && secondObject.getType() == Type.H) {
                                System.out.println("trigger speed threshold: " + SDTObject);
                            }
                            HashMap<SDTObject, Double> temp = new HashMap<>();
                            temp.put(SDTObject, distanceToLine + angle + speed);

                            nextSet.add(temp);
                        }
                    }
                }
                numberOfObjects++;
                if (numberOfObjects == SDTBatch.objectsCount.get(SDTObject.getFileName())) {
                    tracklet.objects.add(nextSet);
                    tracklet.sortLast();
                    nextSet = new ArrayList<>();

                    SDTObject pointToAdd = tracklet.getLastPoint();
                    if (pointToAdd != null) {
                        regression.addData(pointToAdd.getxComponent(), pointToAdd.getyComponent());
                    }
                }
            }
            SDTBatch.tracklets.add(tracklet);
        }
    }

//    private static void fitPointsToRegressions(double threshold, SDTBatch batch) {
//        Map<ArrayList<SDTObject>, SimpleRegression> regressions = batch.getRegressions();
//        ArrayList<SDTObject> data = (ArrayList<SDTObject>) batch.getMainData();
//
//        for (Map.Entry<ArrayList<SDTObject>, SimpleRegression> regression : regressions.entrySet()) {
//            SDTObject last = null;
//            double lastSpeed = Double.MAX_VALUE;
//            ArrayList<SDTObject> regressionPoints = regression.getKey();
//
//            if (regressionPoints.get(0).getType() == Type.H && regressionPoints.get(1).getType() == Type.H) {
//
//                double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0));
//                double baseHeading = regressionPoints.get(0).calculateHeading(regressionPoints.get(1));
//
//                for (SDTObject fitsObject : data) {
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
