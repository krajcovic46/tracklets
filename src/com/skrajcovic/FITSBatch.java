package com.skrajcovic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.ArrayList;

public class FITSBatch {
    private ArrayList<FITSFile> files;

    public FITSBatch(ArrayList<FITSFile> files) {
        this.files = files;
    }

    public void doRegression() {
        SimpleRegression regression = new SimpleRegression();
        for (FITSFile file : files) {
//            regression.addData();
        }
    }

    public String toString() {
        return this.files.toString();
    }
}
