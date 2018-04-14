package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;

public class FITSTracklet {
    public ArrayList<ArrayList<Map<FITSObject, Double>>> objects;
    public double baselineHeading = 0;
    public double baselineSpeed = 0;

    public SimpleRegression simpleRegression = null;

    public FITSTracklet() {
        objects = new ArrayList<>();
    }

    public String toString() {
        return objects.toString() + "\n";
    }

    public void sortLast() {
        int lastIndex = objects.size() - 1;
        ArrayList<Map<FITSObject, Double>> last = objects.get(lastIndex);
        // don't sort if there's less than 2 objects in the array of objects or only the two beginning objects
        if (objects.size() > 2 && last.size() > 1) {
            ArrayList<Map<FITSObject, Double>> newLast;

            newLast = sort(last);
            objects.remove(lastIndex);

            objects.add(newLast);
        }
    }

    private ArrayList<Map<FITSObject, Double>> sort(ArrayList<Map<FITSObject, Double>> list) {
        List<Map.Entry<FITSObject, Double>> l = new LinkedList<>();

        Collections.sort(l, new Comparator<Map.Entry<FITSObject, Double>>() {
            @Override
            public int compare(Map.Entry<FITSObject, Double> o1, Map.Entry<FITSObject, Double> o2) {
                return o1.getValue().compareTo(o2.getValue());
            }
        });

        return new ArrayList<>(list);
    }
}
