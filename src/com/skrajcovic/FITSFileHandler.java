package com.skrajcovic;


import com.skrajcovic.datastructures.Declination;
import com.skrajcovic.datastructures.Rectascension;
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
    public static void processEntry(Map.Entry<File, File> entry, FITSBatch batch, String set) throws Exception {
        boolean read = false;

        File catFile = entry.getKey();
        File fitsFile = entry.getValue();

        if (catFile == null) {
            throw new NullPointerException("Missing .cat file in an entry.");
        }
        if (fitsFile == null) {
            throw new NullPointerException("Missing .fits file in an entry.");
        }

        BufferedReader bf = new BufferedReader(new FileReader(catFile));
        String text;
        Pattern pattern = Pattern.compile("\\s+");

        while ((text = bf.readLine()) != null) {
            if (read) {
                List<String> textList = Arrays.asList(pattern.split(text.trim()));
                if (textList.get(0).equals("?")) {
                    FITSObject fitsObject = new FITSObject();

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

                    FitsHeader fitsHeader = getFitsMetadata(fitsFile);
                    FitsCard fitsCard = fitsHeader.card("DATE-OBS");

                    DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ss");
                    LocalDateTime ld = LocalDateTime.parse(fitsCard.value().toString(), dtf);

                    fitsObject.setFileName(catFile.getName());
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
//        Object[] arr = mergdFiles.entrySet().toArray();
        ArrayList<Map.Entry<File, File>> arr = new ArrayList<>(mergedFiles.entrySet());
        Map.Entry<File, File> firstEntry = arr.get(0);
        Map.Entry<File, File> secondEntry = arr.get(1);
        processEntry(firstEntry, batch, "firstSet");
        processEntry(secondEntry, batch, "secondSet");
        for (int i = 2; i < arr.size(); i++) {
            processEntry(arr.get(i), batch, "mainSet");
        }
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
}
