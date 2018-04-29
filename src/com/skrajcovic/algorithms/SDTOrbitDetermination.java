package com.skrajcovic.algorithms;

import com.skrajcovic.SDTBatch;
import com.skrajcovic.SDTObject;
import com.skrajcovic.SDTOutput;
import com.skrajcovic.orbitdetermination.Observation;
import com.skrajcovic.orbitdetermination.OrbitDetermination;
import com.skrajcovic.orbitdetermination.compute.Kepler;
import com.skrajcovic.utils.Combinations;

import java.util.ArrayList;
import java.util.List;

import static com.skrajcovic.SDTBatch.dirName;

public class SDTOrbitDetermination {

    private static final double AGOLON = 17.2740;
    private static final double AGOLAT = 48.3733;
    private static final double AGOALT = 531.1;


    public static void perform(SDTBatch batch) {
        List<SDTObject> SDTObjects = batch.getRegressionResults();
        if (SDTObjects == null) {
            throw new NullPointerException("Empty results list");
        }
//        Combinations comb = new Combinations(SDTObjects);

//        List<SDTObject[]> resultCombinations = comb.getCombinations();

//        System.out.println(resultCombinations);

        Kepler kepler = new Kepler();

        Observation[] observation = transformObjectsIntoObservations(SDTObjects);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                OrbitDetermination.getTheOdResults(observation, i, j, kepler, dirName);
            }
        }
    }

    private static Observation[] transformObjectsIntoObservations(List<SDTObject> SDTObjects) {
        Observation[] observations = new Observation[SDTObjects.size()];
        for (int i = 0; i < SDTObjects.size(); i++) {
            Observation observation = new Observation();
            observation.lon = Math.toRadians(AGOLON);
            observation.lat = Math.toRadians(AGOLAT);
            observation.alt = AGOALT;
            observation.ra = Math.toRadians(SDTObjects.get(i).getRectascension().getDegValue());
            observation.dec = Math.toRadians(SDTObjects.get(i).getDeclination().getDegValue());
            observation.timeMjd = SDTObjects.get(i).getMjd();
            observations[i] = observation;
        }

        return observations;
    }

    private static Observation[] constructSuitableInputData(SDTObject[] SDTObjects) {
        ArrayList<Observation> list = new ArrayList<>();
        for (SDTObject object : SDTObjects) {
            Observation observation = new Observation();
            observation.lon = Math.toRadians(AGOLON);
            observation.lat = Math.toRadians(AGOLAT);
            observation.alt = AGOALT;
            observation.ra = Math.toRadians(object.getRectascension().getDegValue());
            observation.dec = Math.toRadians(object.getDeclination().getDegValue());
            list.add(observation);
        }
        Observation[] ret = new Observation[list.size()];
        return list.toArray(ret);
    }
}
