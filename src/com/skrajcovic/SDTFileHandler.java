package com.skrajcovic;


import com.skrajcovic.datastructures.Declination;
import com.skrajcovic.datastructures.Rectascension;
import com.skrajcovic.datastructures.Type;
import eap.fits.*;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

public class SDTFileHandler {
    static PrintWriter neuralPreProcessFile;

    public static void processEntry(Map.Entry<File, File> entry, SDTBatch batch, String set) throws Exception {
        boolean read = false;

        File catFile = entry.getKey();
        File fitsFile = entry.getValue();

        if (catFile == null) {
            throw new NullPointerException("Missing .cat file in an entry.");
        }
        if (fitsFile == null) {
            throw new NullPointerException("Missing .fits file in an entry.");
        }

        String fileName = catFile.getName();
        fileName = fileName.substring(0, fileName.length() - 4);

        BufferedReader bf = new BufferedReader(new FileReader(catFile));
        String text;
        Pattern pattern = Pattern.compile("\\s+");

        int numberOfObjects = 0;

        while ((text = bf.readLine()) != null) {
            if (read) {
                List<String> textList = Arrays.asList(pattern.split(text.trim()));
                if (textList.get(0).equals("?") || textList.get(0).equals("H") || textList.get(0).equals("S")) {

                    numberOfObjects++;

                    String type = textList.get(0);

                    Integer RAHours = Integer.valueOf(textList.get(1));
                    Integer RAMinutes = Integer.valueOf(textList.get(2));
                    Double RASeconds = Double.valueOf(textList.get(3));
                    Rectascension ra = new Rectascension(RAHours, RAMinutes, RASeconds);

                    Integer DECDegrees = Integer.valueOf(textList.get(4));
                    Integer DECMinutes = Integer.valueOf(textList.get(5));
                    Double DECSeconds = Double.valueOf(textList.get(6));
                    Declination dec = new Declination(DECDegrees, DECMinutes, DECSeconds);

                    Double magnitude = Double.valueOf(textList.get(7));

                    Double x = Double.valueOf(textList.get(8));
                    Double y = Double.valueOf(textList.get(9));

                    FitsHeader fitsHeader = getFitsHeader(fitsFile);
                    FitsCard fitsCardDATEOBS = fitsHeader.card("DATE-OBS");
                    FitsCard fitsCardEXPTIME;

                    try {
                        fitsCardEXPTIME = fitsHeader.card("EXPTIME");
                    } catch (NoSuchFitsCardException e) {
                        fitsCardEXPTIME = fitsHeader.card("EXPOSURE");
                    }

                    SDTObject SDTObject = new SDTObject(fileName, type, ra, dec, magnitude, x, y);
                    SDTObject.setTime(fitsCardDATEOBS, fitsCardEXPTIME);

                    preProcessForNeural(SDTObject);

                    switch (set) {
                        case "firstSet":
                            batch.firstSetInsert(SDTObject);
                            break;
                        case "secondSet":
                            batch.secondSetInsert(SDTObject);
                            break;
                        default:
                            batch.mainDataInsert(SDTObject);
                            break;
                    }
                }
            }

            if (text.contains("-----")) {
                read = true;
            }
        }

        SDTBatch.objectsCount.put(fileName, numberOfObjects);
        bf.close();
    }

    private static FitsHeader getFitsHeader(File fitsFile) throws IOException {
        InputStreamFitsFile inputStreamFitsFile = new InputStreamFitsFile(new FileInputStream(fitsFile));
        FitsHDU fitsHDU = inputStreamFitsFile.getHDU(0);
        return fitsHDU.getHeader();
    }

    static void readFiles(File folder, SDTBatch batch) throws Exception {
        File[] files = folder.listFiles();
        if (files == null) {
            throw new NullPointerException("No files in the provided folder.");
        }
        ArrayList<File> filesList = new ArrayList<>(Arrays.asList(files));
        File astrometryFolder = null;
        Iterator iterator = filesList.iterator();
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

        List<File> catFiles = getCatFiles(astrometryFolder, filesList.size()*2);

        Map<File, File> mergedFiles = mergeCATWithFITS(catFiles, filesList);

        insertBatch(mergedFiles, batch);

        System.out.println(mergedFiles);
    }

    @SuppressWarnings("unchecked")
    private static void insertBatch(Map<File, File> mergedFiles, SDTBatch batch) throws Exception {
        ArrayList<Map.Entry<File, File>> arr = new ArrayList<>(mergedFiles.entrySet());

//        neuralPreProcessFile = new PrintWriter(new FileWriter("/resources/test.csv", true));

        Map.Entry<File, File> firstEntry = arr.get(0);
        Map.Entry<File, File> secondEntry = arr.get(1);
        processEntry(firstEntry, batch, "firstSet");
        processEntry(secondEntry, batch, "secondSet");
        for (int i = 2; i < arr.size(); i++) {
            processEntry(arr.get(i), batch, "mainSet");
        }

//        neuralPreProcessFile.close();

        System.out.println(SDTBatch.objectsCount);
    }

    private static Map<File, File> mergeCATWithFITS(List<File> catFiles, List<File> fitsFiles) throws Exception {
        if (catFiles.size() != fitsFiles.size()) {
            throw new Exception("The number of .cat and .fits files is different.");
        }
        Map<File, File> mergedFiles = new TreeMap<>();
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

    private static void preProcessForNeural(SDTObject object) throws IOException {
//        BufferedWriter bw = new BufferedWriter(neuralPreProcessFile);

        String type = object.getType() == Type.H ? "1" : "0";
        String x = String.valueOf(object.getX());
        String y = String.valueOf(object.getY());
        String time = String.valueOf(object.getMjd());
//        bw.write("type,x,y,time");
//        bw.write(type+","+x+","+y+","+time);

//        bw.close();
    }

    private static Integer countLines(InputStream is) throws IOException {
        try {
            byte[] c = new byte[1024];
            int count = 0;
            int readChars = 0;
            boolean empty = true;
            while ((readChars = is.read(c)) != -1) {
                empty = false;
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;
                    }
                }
            }
            return (count == 0 && !empty) ? 1 : count;
        } finally {
            is.close();
        }
    }
}
