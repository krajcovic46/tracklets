package com.skrajcovic;


import com.skrajcovic.algorithms.FITSLinearRegression;
import com.skrajcovic.utils.Arguments;
import eap.fitsbrowser.Browser;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class Main {

    public static void main(String[] args) throws Exception {
        List<String> arguments = Arrays.asList(args);
        for (String arg : arguments) {
            if (!Arguments.getArguments().contains(arg)) {
                throw new Exception("Unknown argument \"" + arg + "\"" );
            }
        }

    	FITSBatch batch = new FITSBatch();
//	    FileReaderOld.processFile(batch, "misc//Data_20120305_2.txt");

        FITSFileHandler.readFiles(new File("misc//data//Oct2017//NEA//2017_PR25_R_7"), batch);

        if (arguments.contains("lr") || arguments.contains("linear")) {
            //TODO - linear regression
        }

        if (arguments.contains("lod") || arguments.contains("orbital")) {
            //TODO - orbital determination
        }

        if (arguments.contains("nn") || arguments.contains("neural")) {
            //TODO - neural network - python script probably
        }
        FITSLinearRegression.perform(batch);
	    if (FITSBatch.DEBUG) {
            Browser browser = new Browser();
            try {
                browser.openFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
