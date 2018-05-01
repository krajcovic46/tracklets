package com.skrajcovic;

import com.skrajcovic.datastructures.Declination;
import com.skrajcovic.datastructures.Rectascension;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;

public class SDTOutput {

    private File workingDirectory;
    private SDTBatch batch;

    public SDTOutput(SDTBatch batch) {
        this.batch = batch;
    }

    public void setWorkingDirectory(String pathToDirectory) {
        workingDirectory = new File(pathToDirectory);
    }

    public void processGoodResults() throws NullPointerException, IOException {
        if (workingDirectory == null) {
            throw new NullPointerException("No working directory set");
        }

        File workingFile = null;
        File[] files = workingDirectory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.getName().equals("resultsGood.txt")) {
                    workingFile = file;
                    break;
                }
            }
        } else {
            throw new NullPointerException("Can't find resultsGood file");
        }

        HashSet<Integer> objectIndexes = new HashSet<>();
        if (workingFile != null) {
            BufferedReader bf = new BufferedReader(new FileReader(workingFile));
            String text;

            Pattern pattern = Pattern.compile("\\s+");
            while ((text = bf.readLine()) != null) {
                List<String> textList = Arrays.asList(pattern.split(text.trim()));

                Integer first = Integer.valueOf(textList.get(0));
                Integer second = Integer.valueOf(textList.get(2));
                Integer third = Integer.valueOf(textList.get(4));

                objectIndexes.addAll(Arrays.asList(first, second, third));
            }
        }

        List<SDTObject> objects = batch.getRegressionResults();

        System.out.println(objects);
        System.out.println(objectIndexes);

        List<String> results = new ArrayList<>();
        results.add(batch.getSlope() + "\t" + batch.getIntercept() + "\n");

        for (Integer index : objectIndexes) {
            SDTObject object = objects.get(index);

            String result = String.valueOf(object.getRectascension()) +
                    "\t" +
                    object.getDeclination() +
                    "\t" +
                    object.getMagnitude() +
                    "\t" +
                    object.getMjd() +
                    "\t" +
                    object.getX() +
                    "\t" +
                    object.getY();

            results.add(result);
        }

        System.out.println(results);

        Path file = Paths.get("ODResults/myresults.txt");
        Files.write(file, results, Charset.forName("UTF-8"));

        writeIntoBernFormat(objects, "AGO_70");
        writeIntoMPCFormat(objects, null, null, null, "118");
    }

    private void writeIntoBernFormat(List<SDTObject> objects, String stationName) throws IOException {
        List<String> result = new ArrayList<>();

        for (SDTObject object : objects) {
            String objectName = object.getFileName().substring(0, object.getFileName().indexOf('_'));

            Rectascension ra = object.getRectascension();
            String sec = String.valueOf(ra.getSeconds());
            String seconds = sec.substring(0, sec.indexOf('.')) + sec.substring(sec.indexOf('.')+1, sec.length());
            String rectascension = ra.getHours() + "." + ra.getMinutes() +
                    "" + seconds;

            Declination dec = object.getDeclination();
            String dsec = String.valueOf(dec.getSeconds());
            String dseconds = dsec.substring(0, dsec.indexOf('.')) + dsec.substring(dsec.indexOf('.')+1, dsec.length());
            String dmin = String.valueOf(dec.getMinutes());
            String dminutes = dmin.substring(0, dmin.indexOf('.')) + dmin.substring(dmin.indexOf('.')+1, dmin.length());
            String declination = object.getDeclination().getDegrees() + "." + dminutes +
                    "" + dseconds;

            DecimalFormat format = new DecimalFormat("#.00");
            String ret = stationName + "\t" + objectName + "\t" + object.getMjd() + "\t" + rectascension + "\t" +
                    declination + "\t" + format.format(object.getMagnitude()).replace(',', '.') +
                    "\t" + object.getFileName();

            result.add(ret);
        }

        Path file = Paths.get("ODResults/BernOutput.txt");
        Files.write(file, result, Charset.forName("UTF-8"));
    }

    private void writeIntoMPCFormat(List<SDTObject> objects, String designation, String note1, String note2,
                                    String observatoryCode) throws IOException {
        List<String> result = new ArrayList<>();

        if (designation == null) {
            designation = "designation1";
        } else {
            designation = designation.substring(0, 13);
        }

        if (note1 == null) {
            note1 = "K";
        }

        if (note2 == null) {
            note2 = "C";
        }

        DecimalFormat formatRA = new DecimalFormat("##.000");
        DecimalFormat formatDec = new DecimalFormat("##.00");

        for (SDTObject object : objects) {
            String ret = "";

            Rectascension ra = object.getRectascension();
            String rectascension = ra.getHours() + " " + ra.getMinutes() + " " +
                    formatRA.format(ra.getSeconds()).replace(',', '.');

            Declination dec = object.getDeclination();
            String decDegrees = (dec.getDegrees() > 0) ? ("+" + dec.getDegrees()) : ("-" + dec.getDegrees());
            String decSeconds = formatDec.format(dec.getSeconds());
            if (decSeconds.substring(0, decSeconds.indexOf(',')).length() == 1) {
                decSeconds = "0"+decSeconds;
                decSeconds.replace(',', '.');
            }
            String declination = decDegrees + " " + dec.getMinutes() + " " + decSeconds;

            String magnitude = String.valueOf(object.getMagnitude());

            ret += designation + " " + note1 + note2 + object.getMpcDate() + rectascension + declination + "        " +
                    magnitude + "     " + observatoryCode;
            result.add(ret);
        }

        Path file = Paths.get("ODResults/MPCOutput.txt");
        Files.write(file, result, Charset.forName("UTF-8"));
    }
}
