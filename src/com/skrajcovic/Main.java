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

        String pathName = "misc//data//Oct2017//NEA//";
        String dirName = "00040A_R_4";
        SDTBatch.dirName = dirName;

        SDTOutput output = new SDTOutput(batch);
        output.setWorkingDirectory("ODResults//" + dirName);

        SDTFileHandler.readFiles(new File(pathName + dirName), batch);

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
