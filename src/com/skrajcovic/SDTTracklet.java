package com.skrajcovic;

import java.util.*;

public class SDTTracklet {
    public ArrayList<ArrayList<Map<SDTObject, Double>>> objects;
    public double valueToCompareWith = 0d;

    public SDTTracklet() {
        objects = new ArrayList<>();
    }

    public String toString() {
        return objects.toString() + "\n";
    }

    public SDTObject getLastPoint() {
        ArrayList<Map<SDTObject, Double>> firstLast = objects.get(objects.size() - 1);
        if (firstLast.size() > 0) {
            LinkedHashMap<SDTObject, Double> toReturn = new LinkedHashMap<>(firstLast.get(0));
            Map.Entry<SDTObject, Double> entry = toReturn.entrySet().iterator().next();
            return entry.getKey();
        }
        return null;
    }

    public void sortLast() {
        int lastIndex = objects.size() - 1;
        ArrayList<Map<SDTObject, Double>> last = objects.get(lastIndex);
        // don't sort if there's less than 2 objects in the array of objects or only the two beginning objects
        if (objects.size() > 2 && last.size() > 1) {
            ArrayList<Map<SDTObject, Double>> newLast;

            newLast = sort(last);
            objects.remove(lastIndex);

            objects.add(newLast);
        }
    }

    private ArrayList<Map<SDTObject, Double>> sort(ArrayList<Map<SDTObject, Double>> list) {
        List<Map.Entry<SDTObject, Double>> l = new LinkedList<>();

        Collections.sort(l, new Comparator<Map.Entry<SDTObject, Double>>() {
            @Override
            public int compare(Map.Entry<SDTObject, Double> o1, Map.Entry<SDTObject, Double> o2) {
                Double firstValue = valueToCompareWith - o1.getValue();
                Double secondValue = valueToCompareWith - o2.getValue();
                return secondValue.compareTo(firstValue);
            }
        });

        return new ArrayList<>(list);
    }
}
