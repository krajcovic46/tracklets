package com.skrajcovic;


import com.skrajcovic.handlers.FITSFileHandler;
import com.skrajcovic.handlers.ResourceHandler;
import eap.fitsbrowser.Browser;

import java.io.File;
import java.io.IOException;

public class Main {

    public static void main(String[] args) throws Exception {
    	FITSBatch batch = new FITSBatch();
//	    FileReaderOld.processFile(batch, "misc//Data_20120305_2.txt");
        FITSFileHandler.processFiles(new File("misc//data//Oct2017//NEA//2017_PR25_R_7"));
        ResourceHandler.initValues();
        batch.doTheThing();
	    if (ResourceHandler.getValue("debug") == 1) {
            Browser browser = new Browser();
            try {
                browser.openFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
