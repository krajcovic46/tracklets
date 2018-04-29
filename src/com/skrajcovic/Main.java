package com.skrajcovic;


import com.skrajcovic.algorithms.SDTLinearRegression;
import com.skrajcovic.algorithms.SDTOrbitDetermination;
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

    	SDTBatch batch = new SDTBatch();
//	    FileReaderOld.processFile(batch, "misc//Data_20120305_2.txt");

        String pathName = "misc//data//Oct2017//NEA//";
        String dirName = "2017_PR25_R_7";
        SDTBatch.dirName = dirName;

        SDTOutput output = new SDTOutput(batch);
        output.setWorkingDirectory("ODResults//" + dirName);

        SDTFileHandler.readFiles(new File(pathName + dirName), batch);

        if (arguments.contains("lr") || arguments.contains("linear")) {
            //TODO - linear regression
        }

        if (arguments.contains("lod") || arguments.contains("orbital")) {
            //TODO - orbital determination
        }
        if (arguments.contains("nn") || arguments.contains("neural")) {
            //TODO - neural network - python script probably
        }
        SDTLinearRegression.perform(batch);
        SDTOrbitDetermination.perform(batch);
        output.processGoodResults();
	    if (SDTBatch.DEBUG) {
            Browser browser = new Browser();
            try {
                browser.openFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
