package com.skrajcovic.algorithms;

import com.skrajcovic.FITSBatch;
import com.skrajcovic.FITSObject;
import com.skrajcovic.datastructures.Type;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;

public class FITSLinearRegression {

    public static void perform(FITSBatch batch) {
        FITSLinearRegression.findInitialRegressions((HashSet<FITSObject>) batch.getFirstSet(),
                (HashSet<FITSObject>) batch.getSecondSet(), batch.getRegressions());

        FITSLinearRegression.fitPointsToRegressions(100, batch.getRegressions(),
                (ArrayList<FITSObject>) batch.getMainData());

    }

    private static void findInitialRegressions(HashSet<FITSObject> fSet, HashSet<FITSObject> sSet, Map<ArrayList<FITSObject>, SimpleRegression> regressions) {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getX(), obj1.getY());
                sr.addData(obj2.getX(), obj2.getY());
                regressions.put(new ArrayList<>(Arrays.asList(obj1, obj2)), sr);
            }
        }
    }

    //TODO - make the method work on RADEC instead of x/y
    private static double[] fitPointsToRegressions(double threshold, Map<ArrayList<FITSObject>, SimpleRegression> regressions, ArrayList<FITSObject> data) {
        for (Map.Entry<ArrayList<FITSObject>, SimpleRegression> regression : regressions.entrySet()) {
            FITSObject last = null;
            double lastSpeed = Double.MAX_VALUE;
            ArrayList<FITSObject> regressionPoints = regression.getKey();

            /*TODO - figure out a way to determine correct starting points for regression
            this is a quick hack to make sure that the algorithm is correct for the right
            two starting points (both identified as type H  in Astrometrica tool
             */
            if (regressionPoints.get(0).getType() == Type.H && regressionPoints.get(1).getType() == Type.H) {

//                System.out.println(regressionPoints.get(0).getType());
//                System.out.println(regressionPoints.get(1).getType());

                double averageCombinedSpeed = regressionPoints.get(1).calculateSpeed(regressionPoints.get(0));
                double baseHeading = regressionPoints.get(0).getHeading(regressionPoints.get(1));

//                System.out.println(averageCombinedSpeed);

                int real = 0;
                for (FITSObject fitsObject : data) {
                    if (fitsObject.isWithinLineThreshold(regression.getValue(), threshold)
                            && regressionPoints.get(regressionPoints.size() - 1).isWithinAngleThreshold(fitsObject, baseHeading, 20)) {

                        System.out.println(fitsObject.getFileName());
                        System.out.println(fitsObject.getType());
                        System.out.println(fitsObject.getX());
                        System.out.println(fitsObject.getY());
                        System.out.println("-----------------");

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
                //            for (FITSObject obj : regressionPoints) {
                //                if (obj.isReal()) {
                //                    real++;
                //                }
                //            }
                //            return new double[]{regressionPoints.size(), real, (real / (double) regressionPoints.size()) * 100};
            }
        }
        return null;
    }
}
