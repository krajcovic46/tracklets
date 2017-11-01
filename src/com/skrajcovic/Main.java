package com.skrajcovic;


import eap.fitsbrowser.Browser;

import java.io.File;
import java.io.IOException;

public class Main {

    public static void main(String[] args) {
    	FITSBatch batch = new FITSBatch();
//	    FileReaderOld.processFile(batch, "misc//Data_20120305_2.txt");
        try {
            FITSFileHandler.processFiles(new File("misc//data//Oct2017//NEA//2017_PR25_R_7"));
        } catch (Exception e) {
            e.printStackTrace();
        }
        batch.doTheThing();
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
