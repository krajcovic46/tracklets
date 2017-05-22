package com.skrajcovic;


import eap.fitsbrowser.Browser;

import java.io.IOException;

public class Main {

    public static void main(String[] args) {
    	FITSBatch batch = new FITSBatch();
	    FileReader.processFile(batch, "Data_20120305_2.txt");
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
