package com.skrajcovic;


import com.skrajcovic.utils.Declination;
import com.skrajcovic.utils.Rectascension;
import eap.fits.FitsCard;
import eap.fits.FitsHDU;
import eap.fits.FitsHeader;
import eap.fits.InputStreamFitsFile;

import java.io.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.regex.Pattern;

public class FITSFileHandler {
    public static void processEntry(Map.Entry<File, File> entry, FITSBatch batch, String set) throws IOException{
        boolean read = false;

        File catFile = entry.getKey();
        File fitsFile = entry.getValue();

        BufferedReader bf = new BufferedReader(new FileReader(catFile));
        String text;
        Pattern pattern = Pattern.compile("\\s+");

        while ((text = bf.readLine()) != null) {
            if (read) {
                List<String> textList = Arrays.asList(pattern.split(text.trim()));
                if (textList.get(0).equals("?")) {
                    FITSObject fitsObject = new FITSObject();

                    Rectascension ra = new Rectascension();
                    ra.setHours(Integer.valueOf(textList.get(1)));
                    ra.setMinutes(Integer.valueOf(textList.get(2)));
                    ra.setSeconds(Double.valueOf(textList.get(3)));

                    Declination dec = new Declination();
                    dec.setDegrees(Integer.valueOf(textList.get(4)));
                    dec.setMinutes(Integer.valueOf(textList.get(5)));
                    dec.setSeconds(Double.valueOf(textList.get(6)));

                    Double magnitude = Double.valueOf(textList.get(7));

                    Double x = Double.valueOf(textList.get(8));
                    Double y = Double.valueOf(textList.get(9));

                    FitsHeader fitsHeader = getFitsMetadata(fitsFile);
                    FitsCard fitsCard = fitsHeader.card("DATE-OBS");

                    DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ss");
                    LocalDateTime ld = LocalDateTime.parse(fitsCard.value().toString(), dtf);

                    fitsObject.setRectascension(ra);
                    fitsObject.setDeclination(dec);
                    fitsObject.setMagnitude(magnitude);
                    fitsObject.setX(x);
                    fitsObject.setY(y);
                    fitsObject.setLocalDateTime(ld);

                    switch (set) {
                        case "firstSet":
                            batch.firstSetInsert(fitsObject);
                            break;
                        case "secondSet":
                            batch.secondSetInsert(fitsObject);
                            break;
                        default:
                            batch.mainDataInsert(fitsObject);
                            break;
                    }
                }
            }

            if (text.contains("-----")) {
                read = true;
            }
        }

    }

    private static FitsHeader getFitsMetadata(File fitsFile) throws IOException {
        InputStreamFitsFile inputStreamFitsFile = new InputStreamFitsFile(new FileInputStream(fitsFile));
        FitsHDU fitsHDU = inputStreamFitsFile.getHDU(0);
        return fitsHDU.getHeader();
    }

    static void readFiles(File folder, FITSBatch batch) throws Exception {
        File[] files = folder.listFiles();
        if (files == null) {
            throw new NullPointerException();
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
    private static void insertBatch(Map<File, File> mergedFiles, FITSBatch batch) throws Exception {
        Object[] arr = mergedFiles.entrySet().toArray();
        Map.Entry<File, File> firstEntry = (Map.Entry<File, File>) arr[0];
        Map.Entry<File, File> secondEntry = (Map.Entry<File, File>) arr[1];
        processEntry(firstEntry, batch, "firstSet");
        processEntry(secondEntry, batch, "secondSet");
        for (int i = 2; i < mergedFiles.size(); i++) {
            processEntry(firstEntry, batch, "mainSet");
        }
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
}
//                try {
//                    BufferedReader bf = new BufferedReader(new java.io.FileReader(fileEntry));
//                    String text;
//                    Pattern pattern = Pattern.compile("\\s+");
//
//                    while ((text = bf.readLine()) != null) {
//                        String splitLine[] = pattern.split(text);
//                        System.out.println(Arrays.toString(splitLine));
//                    }
//
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
