package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class SDTBatch {
    private List<SDTObject> data;
    private Set<SDTObject> fSet;
    private Set<SDTObject> sSet;
    private Map<ArrayList<SDTObject>, SimpleRegression> regressions;
    private List<SDTObject> regressionResults;
    public static Map<String, Integer> objectsCount = new TreeMap<>();

    public static String dirName;

    public HashSet<SDTTracklet> tracklets = new HashSet<>();

    public static final boolean DEBUG = false;
    public static final boolean RADEC = true;

    private double slope;
    private double intercept;

    public SDTBatch() {
        regressions = new HashMap<>();
        data = new ArrayList<>();
        fSet = new HashSet<>();
        sSet = new HashSet<>();
    }

    public void mainDataInsert(SDTObject object) {
        data.add(object);
    }

    public void firstSetInsert(SDTObject object) {
        fSet.add(object);
    }

    public void secondSetInsert(SDTObject object) {
        sSet.add(object);
    }

    public Set<SDTObject> getFirstSet() {
        return this.fSet;
    }

    public Set<SDTObject> getSecondSet() {
        return this.sSet;
    }

    public List<SDTObject> getMainData() {
        return data;
    }

    public Map<ArrayList<SDTObject>, SimpleRegression> getRegressions() {
        return this.regressions;
    }

    public void setRegressionResults(List<SDTObject> list) {
        regressionResults = list;
    }

    public List<SDTObject> getRegressionResults() {
        return regressionResults;
    }

    public void filterOutEmptyTracklets() {
        Iterator iterator = tracklets.iterator();
        while (iterator.hasNext()) {
            SDTTracklet tracklet = (SDTTracklet) iterator.next();
            int trackletSize = tracklet.objects.size();
            int numberOfEmpty = 0;
            for (ArrayList<HashMap<SDTObject, Double>> list : tracklet.objects) {
                if (list.isEmpty()) {
                    numberOfEmpty++;
                }
            }
            if (trackletSize - numberOfEmpty == 2) {
                iterator.remove();
            }
        }
    }

    public double getSlope() {
        return slope;
    }

    public void setSlope(double slope) {
        this.slope = slope;
    }

    public double getIntercept() {
        return intercept;
    }

    public void setIntercept(double intercept) {
        this.intercept = intercept;
    }

    public String toString() {
        return this.data.toString();
    }
}
