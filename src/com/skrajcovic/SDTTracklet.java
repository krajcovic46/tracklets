package com.skrajcovic;

import java.util.*;

public class SDTTracklet {
    public ArrayList<ArrayList<HashMap<SDTObject, Double>>> objects;
    public double valueToCompareWith = 0d;

    private double slope;
    private double intercept;

    public SDTTracklet() {
        objects = new ArrayList<>();
    }

    public String toString() {
        ArrayList<String> retList = new ArrayList<>(objects.size()*2);
        for (ArrayList<HashMap<SDTObject, Double>> object : objects) {
            retList.add("\n" + object.toString());
        }
        return retList.toString();
    }

    public SDTObject getLastPoint() {
        ArrayList<HashMap<SDTObject, Double>> firstLast = objects.get(objects.size() - 1);
        if (firstLast.size() > 0) {
            return getFirstObject(firstLast.get(0));
        }
        return null;
    }

    public ArrayList<HashMap<SDTObject, Double>> getListOfObjectsOnIndex(int index) {
        for (int i = 0; i < objects.size(); i++) {
            if (i == index) {
                return objects.get(i);
            }
        }
        return null;
    }

    public SDTObject getMostProbableObject(ArrayList<HashMap<SDTObject, Double>> list) {
        if (list.size() > 0) {
            return getFirstObject(list.get(0));
        }
        return null;
    }

    private SDTObject getFirstObject(Map<SDTObject, Double> map) {
        LinkedHashMap<SDTObject, Double> toReturn = new LinkedHashMap<>(map);
        Map.Entry<SDTObject, Double> entry = toReturn.entrySet().iterator().next();
        return entry.getKey();
    }

    public void sortLast() {
        int lastIndex = objects.size() - 1;
        ArrayList<HashMap<SDTObject, Double>> last = objects.get(lastIndex);
        // don't sort if there's less than 2 objects in the array of objects or only the two beginning objects
        if (objects.size() > 2 && last.size() > 1) {
            ArrayList<HashMap<SDTObject, Double>> newLast;

            newLast = sort(last);
            objects.remove(lastIndex);

            objects.add(newLast);
        }
    }

    private ArrayList<HashMap<SDTObject, Double>> sort(ArrayList<HashMap<SDTObject, Double>> list) {
//        List<Map.Entry<SDTObject, Double>> l = new LinkedList<>();

//        Collections.sort(l, new Comparator<Map.Entry<SDTObject, Double>>() {
//            @Override
//            public int compare(Map.Entry<SDTObject, Double> o1, Map.Entry<SDTObject, Double> o2) {
//                Double firstValue = valueToCompareWith - o1.getValue();
//                Double secondValue = valueToCompareWith - o2.getValue();
////                if (firstValue < secondValue) {
////                    return -1;
////                } else if (firstValue == secondValue) {
////                    return 0;
////                }
////                return 1;
//                return secondValue.compareTo(firstValue);
//            }
//        });

        Collections.sort(list, new Comparator<HashMap<SDTObject, Double>>() {
            @Override
            public int compare(HashMap<SDTObject, Double> o1, HashMap<SDTObject, Double> o2) {
                Double o1value = o1.entrySet().iterator().next().getValue();
                Double o2value = o2.entrySet().iterator().next().getValue();

                Double firstValue = Math.abs(valueToCompareWith - o1value);
                Double secondValue = Math.abs(valueToCompareWith - o2value);

                return firstValue.compareTo(secondValue);
            }
        });
        return new ArrayList<>(list);
    }

    public void setSlope(double slope) {
        this.slope = slope;
    }

    public void setIntercept(double intercept) {
        this.intercept = intercept;
    }

    public double getSlope() {
        return this.slope;
    }

    public double getIntercept() {
        return this.intercept;
    }
}
