package com.skrajcovic.algorithms;

import com.skrajcovic.FITSBatch;
import com.skrajcovic.FITSObject;
import com.skrajcovic.orbitdetermination.Observation;
import com.skrajcovic.orbitdetermination.OrbitDetermination;
import com.skrajcovic.orbitdetermination.compute.Kepler;
import com.skrajcovic.utils.Combinations;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class FITSOrbitDetermination {

    private static final double AGOLON = 17.2740;
    private static final double AGOLAT = 48.3733;
    private static final double AGOALT = 531.1;


    public static void perform(FITSBatch batch) {
        List<FITSObject> fitsObjects = batch.getRegressionResults();
        Combinations comb = new Combinations(fitsObjects);

        List<FITSObject[]> resultCombinations = comb.getCombinations();

//        System.out.println(resultCombinations);

        Kepler kepler = new Kepler();

        for (FITSObject[] combination : resultCombinations) {
            Observation[] observation = constructSuitableInputData(combination);
            System.out.println("LENGTH " + observation.length);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    OrbitDetermination.getTheOdResults(observation, i, j, kepler, "2000_072B _1");
                }
            }
        }
    }

    private static Observation[] constructSuitableInputData(FITSObject[] fitsObjects) {
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
        Observation[] ret = new Observation[list.size()];
        return list.toArray(ret);
    }
}
