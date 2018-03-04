package com.skrajcovic.algorithms;

import com.skrajcovic.FITSObject;
import com.skrajcovic.orbitdetermination.Observation;

import java.util.Set;

public class FITSOrbitalDetermination {

    private Observation[] constructSuitableInputData(Set<FITSObject> fitsObjects) {
        return (Observation[]) fitsObjects.toArray();
    }
}
