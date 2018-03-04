package com.skrajcovic.algorithms;

import com.skrajcovic.FITSObject;
import com.skrajcovic.orbitdetermination.Observation;
import com.skrajcovic.orbitdetermination.OrbitDetermination;
import com.skrajcovic.orbitdetermination.compute.Kepler;

import java.util.ArrayList;
import java.util.Set;

public class FITSOrbitDetermination {

    public static void perform(Set<FITSObject> fitsObjects) {
        Observation[] observation = constructSuitableInputData(fitsObjects);
        Kepler kepler = new Kepler();

        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < 1; j++) {
                //TODO - change algo so it doesnt output to a file
                OrbitDetermination.getTheOdResults(observation, i, j, kepler, "2000_072B _1");
            }
        }
    }

    private static Observation[] constructSuitableInputData(Set<FITSObject> fitsObjects) {
        ArrayList<Observation> list = new ArrayList<>();
        for (FITSObject object : fitsObjects) {
            Observation observation = new Observation();
            observation.lon = Math.toRadians(17.2740);
            observation.lat = Math.toRadians(48.3733);
            observation.alt = 531.1;
            observation.ra = Math.toRadians(object.getRectascension().getDegValue());
            observation.dec = Math.toRadians(object.getDeclination().getDegValue());
            list.add(observation);
        }
        return (Observation[]) list.toArray();
    }
}
