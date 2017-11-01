package com.skrajcovic.handlers;


import com.skrajcovic.FITSObject;
import com.skrajcovic.datastructures.Declination;
import com.skrajcovic.datastructures.Rectascension;
import com.skrajcovic.datastructures.Type;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

public class FITSFileHandler {

    private static final int RA_INDEX = 1;
    private static final int DA_INDEX = 5;
    private static final int MAG_INDEX = 9;
    private static final int X_INDEX = 11;
    private static final int Y_INDEX = 12;

    private FITSFileHandler() {}

    public static void processFiles(File folder) throws Exception {
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
            try {
                readCATFile(entry.getKey());
            } catch (Exception e) {
                e.printStackTrace();
            }
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

    private static void readCATFile(File file) throws Exception {
        try {
            BufferedReader bf = new BufferedReader(new java.io.FileReader(file));
            String text;
            boolean creatingFITS = false;
            Pattern pattern = Pattern.compile("\\s+");

            List<FITSObject> fitsObjects = new ArrayList<>();

            while ((text = bf.readLine()) != null) {
//                text = text.trim();
//                List<String> splitLine = Arrays.asList(pattern.split(text));
                if (creatingFITS) {
                    FITSObject fitsObject = new FITSObject();
                    int type_i = ResourceHandler.getValue("type_i");
                    String typeString = text.substring(type_i-1, type_i);
                    Type type;
                    switch (typeString) {
                        case "R":
                            type = Type.R;
                            break;
                        case "S":
                            type = Type.S;
                            break;
                        case "?":
                            type = Type.UNKNOWN;
                            break;
                        default:
                            type = Type.UNKNOWN;
                            break;
                        // last two are distinguished for the future
                    }
                    fitsObject.setType(type);
                    int raHoursI = ResourceHandler.getValue("ra_hours_i");
                    int raMinutesI = ResourceHandler.getValue("ra_minutes_i");
                    int raSecondsI = ResourceHandler.getValue("ra_seconds_i");
                    fitsObject.setRectascension(new Rectascension(Integer.valueOf(text.substring(raHoursI-2, raHoursI)),
                            Integer.valueOf(text.substring(raMinutesI-2, raMinutesI)),
                            Double.valueOf(text.substring(raSecondsI-6, raSecondsI))));
                    int decDegreesI = ResourceHandler.getValue("dec_degrees_i");
                    int decMinutesI = ResourceHandler.getValue("dec_minutes_i");
                    int decSecondsI = ResourceHandler.getValue("dec_seconds_i");
                    fitsObject.setDeclination(new Declination(Integer.valueOf(text.substring(decDegreesI-3, decDegreesI)),
                            Integer.valueOf(text.substring(decMinutesI-2, decMinutesI)),
                            Double.valueOf(text.substring(decSecondsI-5, decSecondsI))));
                    int magI = ResourceHandler.getValue("mag_i");
                    fitsObject.setMagnitude(Double.valueOf(text.substring(magI-5, magI)));
                    int xCoordI = ResourceHandler.getValue("x_coord_i");
                    fitsObject.setX(Double.valueOf(text.substring(xCoordI-7, xCoordI).trim()));
                    int yCoordI = ResourceHandler.getValue("x_coord_i");
                    fitsObject.setY(Double.valueOf(text.substring(yCoordI-7, yCoordI).trim()));
                    System.out.println("FITSObject" + fitsObject);
                }
                if (text.contains("----")) {
                    creatingFITS = true;
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
