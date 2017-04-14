package com.skrajcovic;

import java.util.ArrayList;

public class FITSFile {
    private ArrayList<FITSObject> fitsObjects;

    public FITSFile() {
        fitsObjects = new ArrayList<>();
    }

    public void addObject(FITSObject fitsObject) {
        fitsObjects.add(fitsObject);
    }

    public String toString() {
        return fitsObjects.toString();
    }
}
