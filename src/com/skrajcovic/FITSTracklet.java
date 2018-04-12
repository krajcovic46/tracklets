package com.skrajcovic;

import java.util.ArrayList;
import java.util.HashSet;

public class FITSTracklet {
    public ArrayList<HashSet<FITSObject>> objects;

    public FITSTracklet() {
        objects = new ArrayList<>();
    }

    public String toString() {
        return objects.toString();
    }
}
