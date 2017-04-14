package com.skrajcovic;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Pattern;

public class FileReader {
    public static ArrayList<FITSFile> processFile(String fileLocation) {
        boolean skipFirstLine = true;
        ArrayList<FITSFile> files = new ArrayList<>();

        File file = new File(fileLocation);
        try {
            BufferedReader bf = new BufferedReader(new java.io.FileReader(file));
            String text;
            Pattern pattern = Pattern.compile("\t");

            files.add(new FITSFile());

            while ((text = bf.readLine()) != null) {
                if (!skipFirstLine) {
                    String[] splitLine = pattern.split(text.replace(',', '.'), 0);

                    FITSFile last = files.get(files.size()-1);

                    if (splitLine.length == 0) {
                        files.add(new FITSFile());
                        continue;
                    }

                    FITSObject object = new FITSObject(splitLine);
                    last.addObject(object);
                } else {
                    skipFirstLine = false;
                }
            }

            bf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return files;
    }
}
