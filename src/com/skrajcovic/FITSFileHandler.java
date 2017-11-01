package com.skrajcovic;


import com.skrajcovic.utils.Type;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

public class FITSFileHandler {
    public static void processFile(FITSBatch batch, String fileLocation) {
        boolean skipFirstLine = true;

        File file = new File(fileLocation);
    }

    static void processFiles(File folder) throws Exception {
        File[] files = folder.listFiles();
        if (files == null) {
            throw new NullPointerException();
        }
        ArrayList<File> fitsFilesList = new ArrayList<>(Arrays.asList(files));
        File astrometryFolder = null;
        Iterator iterator = fitsFilesList.iterator();
        while (iterator.hasNext()) {
            File fileEntry = (File) iterator.next();
            if (fileEntry.getName().equals("Astrometry")) {
                astrometryFolder = fileEntry;
                iterator.remove();
            }
        }
        if (astrometryFolder == null) {
            throw new Exception("Missing \"Astrometry\" folder in the provided directory.");
        }

        List<File> catFiles = getCatFiles(astrometryFolder, fitsFilesList.size()*2);

        Map<File, File> mergedFiles = mergeCATWithFITS(catFiles, fitsFilesList);

        insertBatch(mergedFiles);

        System.out.println(mergedFiles);
    }

    private static void insertBatch(Map<File, File> mergedFiles) {
        FITSObject fitsObject = new FITSObject();
        for (Map.Entry<File, File> entry : mergedFiles.entrySet()) {
            readCATFile(entry.getKey());
        }
        System.out.println(mergedFiles.size());
    }

    private static Map<File, File> mergeCATWithFITS(List<File> catFiles, List<File> fitsFiles) throws Exception {
        if (catFiles.size() != fitsFiles.size()) {
            throw new Exception("The number of .cat and .fits files is different.");
        }
        Map<File, File> mergedFiles = new HashMap<>();
        for (File catFileEntry : catFiles) {
            for (File fitsFileEntry : fitsFiles) {
                String catFileNameSplit =  catFileEntry.getName().split("\\.")[0];
                String fitsFileNameSplit =  fitsFileEntry.getName().split("\\.")[0];
                if (catFileNameSplit.equals(fitsFileNameSplit)) {
                    mergedFiles.put(catFileEntry, fitsFileEntry);
                }
            }
        }
        return mergedFiles;
    }

    private static List<File> getCatFiles(File folder, int numberOfFiles) {
        File[] files = folder.listFiles();
        if (files == null) {
            throw new NullPointerException();
        }
        List<File> catFiles = new ArrayList<>(numberOfFiles);
        for (File fileEntry : files) {
            String fileName = fileEntry.getName();
            int dot = fileName.lastIndexOf('.');
            if (fileName.substring(dot + 1).equals("cat")) {
                catFiles.add(fileEntry);
            }
        }
        return catFiles;
    }

    private static void readCATFile(File file) {
        try {
            BufferedReader bf = new BufferedReader(new java.io.FileReader(file));
            String text;
            boolean creatingFITS = false;
            Pattern pattern = Pattern.compile("\\s+");

            List<FITSObject> fitsObjects = new ArrayList<>();

            while ((text = bf.readLine()) != null) {
                text = text.trim();
                String splitLine[] = pattern.split(text);
                if (creatingFITS) {
                    System.out.println(Arrays.toString(splitLine));

                }
                if (splitLine.length == 1 && !splitLine[0].equals("")) {
                    creatingFITS = true;
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
