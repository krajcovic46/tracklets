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

    public static final boolean DEBUG = false;

    public FITSBatch() {
        regressions = new HashMap<>();
        data = new ArrayList<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

//    public void doTheThing() {
//        System.out.println("fset: " + fSet);
//        System.out.println("sSet: " + sSet);
//        System.out.println("data: " + data);
//        findRegressions();
//        for (int i = 100; i < 101; i++) {
//            double[] tmp = fitPointsToRegressions(i);
//            if (tmp != null) {
//                System.out.format("Threshold: %d\nNumber of all points under threshold: %d\nReal points: %d\nSuccess rate: %f",
//                        i, (int)tmp[0], (int)tmp[1], tmp[2]); System.out.println("%");
//                System.out.println("---------------------------------");
//            }
//        }
//    }

    public void mainDataInsert(FITSObject object) {
        data.add(object);
    }

    public void firstSetInsert(FITSObject fitsObject) {
        fSet.add(fitsObject);
    }

    public void secondSetInsert(FITSObject fitsObject) {
        sSet.add(fitsObject);
    }

    public Set<FITSObject> getFirstSet() {
        return this.fSet;
    }

    public Set<FITSObject> getSecondSet() {
        return this.sSet;
    }

    public List<FITSObject> getMainData() {
        return data;
    }

    public Map<ArrayList<FITSObject>, SimpleRegression> getRegressions() {
        return this.regressions;
    }

    public String toString() {
//        for (double[] d : data.keySet()) {
//            System.out.println(String.valueOf(d[0]).replace(".", ",") + " " + String.valueOf(d[1]).replace(".", ","));
//        }
//        System.out.println();
        return this.data.toString();
    }
}
