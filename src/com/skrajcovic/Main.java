package com.skrajcovic;


import eap.fitsbrowser.Browser;

import java.io.File;
import java.io.IOException;

public class Main {

    public static void main(String[] args) {
    	FITSBatch batch = new FITSBatch();
//	    FileReaderOld.processFile(batch, "misc//Data_20120305_2.txt");
        FITSFileHandler.readCATFiles(new File("misc//data"));
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
