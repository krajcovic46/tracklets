package com.skrajcovic.algorithms;

import com.skrajcovic.FITSBatch;
import com.skrajcovic.FITSObject;
import com.skrajcovic.datastructures.Type;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSLinearRegression {

    public static void perform(FITSBatch batch) {
        FITSLinearRegression.findInitialRegressions((HashSet<FITSObject>) batch.getFirstSet(),
                (HashSet<FITSObject>) batch.getSecondSet(), batch.getRegressions());

        FITSLinearRegression.fitPointsToRegressions(50, batch);

        System.out.println(batch.getRegressionResults());

    }

    private static void findInitialRegressions(HashSet<FITSObject> fSet, HashSet<FITSObject> sSet, Map<ArrayList<FITSObject>, SimpleRegression> regressions) {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getxComponent(), obj1.getyComponent());
                sr.addData(obj2.getyComponent(), obj2.getyComponent());
                regressions.put(new ArrayList<>(Arrays.asList(obj1, obj2)), sr);
            }
        }
    }

    //TODO - make the method work on RADEC instead of x/y
    private static void fitPointsToRegressions(double threshold, FITSBatch batch) {
        Map<ArrayList<FITSObject>, SimpleRegression> regressions = batch.getRegressions();
        ArrayList<FITSObject> data = (ArrayList<FITSObject>) batch.getMainData();

        for (Map.Entry<ArrayList<FITSObject>, SimpleRegression> regression : regressions.entrySet()) {
            FITSObject last = null;
            double lastSpeed = Double.MAX_VALUE;
            ArrayList<FITSObject> regressionPoints = regression.getKey();

            /*TODO - figure out a way to determine correct starting points for regression
            this is a quick hack to make sure that the algorithm is correct for the right
            two starting points (both identified as type H  in Astrometrica tool
             */
            if (regressionPoints.get(0).getType() == Type.H && regressionPoints.get(1).getType() == Type.H) {

                double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0));
                double baseHeading = regressionPoints.get(0).getHeading(regressionPoints.get(1));

                for (FITSObject fitsObject : data) {
                    if (fitsObject.isWithinLineThreshold(regression.getValue(), threshold)
                            && regressionPoints.get(regressionPoints.size() - 1).isWithinAngleThreshold(fitsObject, baseHeading, 20)) {

//                        System.out.println(fitsObject.getFileName());
//                        System.out.println(fitsObject.getType());
//                        System.out.println(fitsObject.getX());
//                        System.out.println(fitsObject.getY());
//                        System.out.println("-----------------");

                        if (last != null && !fitsObject.getFileName().equals(last.getFileName()) && !regressionPoints.contains(last)) {
                            regressionPoints.add(last);
                            regression.getValue().addData(last.getxComponent(), last.getyComponent());

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
                    regression.getValue().addData(last.getxComponent(), last.getyComponent());
                }
                batch.setRegressionResults(regressionPoints);
            }
        }
    }
}
