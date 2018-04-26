package com.skrajcovic.utils;

import com.skrajcovic.SDTObject;

import java.util.ArrayList;
import java.util.List;

public class Combinations {

    List<SDTObject[]> outputList;
    List<SDTObject> inputList;

    public Combinations(List<SDTObject> inputList) {
        this.outputList = new ArrayList<>();
        this.inputList = inputList;
    }

    public List<SDTObject[]> getCombinations() {
        fillCombination(new SDTObject[3], 0, 0);
        return outputList;
    }

    private void fillCombination(SDTObject[] combination, int innerIndex, int outerIndex) {
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
