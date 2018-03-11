package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSBatch {
    private List<FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<ArrayList<FITSObject>, SimpleRegression> regressions;
    private List<FITSObject> regressionResults;

    public static final boolean DEBUG = false;

    public FITSBatch() {
        regressions = new HashMap<>();
        data = new ArrayList<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
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

    public void setRegressionResults(List<FITSObject> list) {
        regressionResults = list;
    }

    public List<FITSObject> getRegressionResults() {
        return regressionResults;
    }

    public String toString() {
        return this.data.toString();
    }
}
