package com.skrajcovic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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

        System.out.println("help");
        System.out.println(objects);
        System.out.println(objectIndexes);

        List<String> results = new ArrayList<>();
        for (Integer index : objectIndexes) {
            System.out.println(index);
            SDTObject object = objects.get(index);

            StringBuilder builder = new StringBuilder();
            builder.append(object.getRectascension());
            builder.append("\t");
            builder.append(object.getDeclination());
            builder.append("\t");
            builder.append(object.getMagnitude());
            builder.append("\t");
            builder.append(object.getMjd());
            builder.append("\t");
            builder.append(object.getX());
            builder.append("\t");
            builder.append(object.getY());
            String result = builder.toString();

            System.out.println("result " + result);

            results.add(result);
        }

        System.out.println(results);

        Path file = Paths.get("ODResults/myresults.txt");
        Files.write(file, results, Charset.forName("UTF-8"));
    }
}
