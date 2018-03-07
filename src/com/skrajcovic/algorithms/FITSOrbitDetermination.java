package com.skrajcovic.algorithms;

import com.skrajcovic.FITSBatch;
import com.skrajcovic.FITSObject;
import com.skrajcovic.orbitdetermination.Observation;
import com.skrajcovic.orbitdetermination.OrbitDetermination;
import com.skrajcovic.orbitdetermination.compute.Kepler;

import java.util.ArrayList;
import java.util.Set;

public class FITSOrbitDetermination {

    private static final double AGOLON = 17.2740;
    private static final double AGOLAT = 48.3733;
    private static final double AGOALT = 531.1;


    public static void perform(FITSBatch batch) {
        //TODO - create suitable structure in batch to hold results
        Set<FITSObject> fitsObjects = batch.getFirstSet();

        Observation[] observation = constructSuitableInputData(fitsObjects);
        Kepler kepler = new Kepler();

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                //TODO - change algo so it doesnt output to a file?
                OrbitDetermination.getTheOdResults(observation, i, j, kepler, "2000_072B _1");
            }
        }
    }

    //TODO len na tri objekty
    private static Observation[] constructSuitableInputData(Set<FITSObject> fitsObjects) {
        ArrayList<Observation> list = new ArrayList<>();
        for (FITSObject object : fitsObjects) {
            Observation observation = new Observation();
            observation.lon = Math.toRadians(AGOLON);
            observation.lat = Math.toRadians(AGOLAT);
            observation.alt = AGOALT;
            observation.ra = Math.toRadians(object.getRectascension().getDegValue());
            observation.dec = Math.toRadians(object.getDeclination().getDegValue());
            list.add(observation);
        }
        return (Observation[]) list.toArray();
    }
}
