package com.skrajcovic;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.Buffer;
import java.util.*;
import java.util.regex.Pattern;

public class FITSFileHandler {
    public static void processEntry(Map.Entry<File, File> entry) throws IOException{
        boolean read = false;

        File catFile = entry.getKey();
        File fitsFile = entry.getValue();

        BufferedReader bf = new BufferedReader(new FileReader(catFile));
        String text;
        Pattern pattern = Pattern.compile("\\s+");

        while ((text = bf.readLine()) != null) {
            if (read) {

            }

            if (text.contains("-----")) {
                read = true;
            }
        }

    }

    static void readFiles(File folder) throws Exception {
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

        System.out.println(mergedFiles);
    }

    @SuppressWarnings("unchecked")
    private static void insertBatch(Map<File, File> mergedFiles, FITSBatch batch) {
        Object[] arr = mergedFiles.entrySet().toArray();
        Map.Entry<File, File> firstEntry = (Map.Entry<File, File>) arr[0];
        Map.Entry<File, File> secondEntry = (Map.Entry<File, File>) arr[1];
        for (int i = 2; i < mergedFiles.size(); i++) {

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
