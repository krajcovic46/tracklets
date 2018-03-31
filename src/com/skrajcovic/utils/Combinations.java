package com.skrajcovic.utils;

import com.skrajcovic.FITSObject;
import com.skrajcovic.orbitdetermination.Observation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Combinations {

    List<FITSObject[]> outputList;
    List<FITSObject> inputList;

    public Combinations(List<FITSObject> inputList) {
        this.outputList = new ArrayList<>();
        this.inputList = inputList;
    }

    public List<FITSObject[]> getCombinations() {
        fillCombination(new FITSObject[3], 0, 0);
        return outputList;
    }

    private void fillCombination(FITSObject[] combination, int innerIndex, int outerIndex) {
        if (innerIndex == 3) {
            outputList.add(combination);
            return;
        }

        if (outerIndex >= inputList.size()) {
            return;
        }

        combination[innerIndex] = inputList.get(outerIndex);
        fillCombination(combination, innerIndex+1, outerIndex+1);

        fillCombination(combination, innerIndex, outerIndex+1);
    }
}
