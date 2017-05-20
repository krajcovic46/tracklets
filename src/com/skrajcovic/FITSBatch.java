package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSBatch {
    private HashMap<double[], FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Set<SimpleRegression> regressions;

    public FITSBatch() {
        regressions = new HashSet<>();
        data = new HashMap<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

    public FITSBatch(HashMap<double[], FITSObject> data) {
        this.data = data;
    }

    public void getPoints() {

        for (Map.Entry<double[], FITSObject> e : data.entrySet()) {
            //pass for now
//            if (e.getKey()[1] == m*e.getKey()[0] + c) {
                //regression check
//            }
        }

        System.out.println(fSet);

        System.out.println(sSet);

        findreal();
    }

    void findreal() {
        for (FITSObject obj1 : fSet) {
            for (FITSObject obj2 : sSet) {
                SimpleRegression sr = new SimpleRegression();
                sr.addData(obj1.getX(), obj1.getY());
                sr.addData(obj2.getX(), obj2.getY());
                regressions.add(sr);
            }
        }
        System.out.println(regressions.size());
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
        for (double[] d : data.keySet()) {
            System.out.println(String.valueOf(d[0]).replace(".", ",") + " " + String.valueOf(d[1]).replace(".", ","));
        }
        System.out.println();
        return this.data.toString();
    }
}
