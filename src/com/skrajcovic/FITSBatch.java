package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.lang.reflect.Array;
import java.util.*;

public class FITSBatch {
    private List<FITSObject> data;
    private Set<FITSObject> fSet;
    private Set<FITSObject> sSet;
    private Map<ArrayList<FITSObject>, SimpleRegression> regressions;
    private List<FITSObject> regressionResults;
    public static Map<String, Integer> objectsCount = new TreeMap<>();

    public static HashSet<FITSTracklet> tracklets = new HashSet<>();

    public static final boolean DEBUG = false;
    public static final boolean RADEC = true;

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

    public void filterOutEmptyTracklets() {
        Iterator iterator = tracklets.iterator();
        while (iterator.hasNext()) {
            FITSTracklet tracklet = (FITSTracklet) iterator.next();
            int trackletSize = tracklet.objects.size();
            int numberOfEmpty = 0;
            for (ArrayList<Map<FITSObject, Double>> list : tracklet.objects) {
                if (list.isEmpty()) {
                    numberOfEmpty++;
                }
            }
            if (trackletSize - numberOfEmpty == 2) {
                iterator.remove();
            }
        }
    }

    public String toString() {
        return this.data.toString();
    }
}
