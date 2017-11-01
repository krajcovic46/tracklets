package com.skrajcovic.handlers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

public class ResourceHandler {

    private static final String filePath = "settings.txt";
    private static Map<String, Integer> resources;

    private ResourceHandler() {
    }

    public static void initValues() {
        BufferedReader bf;
        Pattern pattern = Pattern.compile("=");
        resources = new HashMap<>();
        try {
            bf = new BufferedReader(new java.io.FileReader(new File(filePath)));
            String text;
            while ((text = bf.readLine()) != null) {
                String[] splitLine = pattern.split(text);
                resources.put(splitLine[0], Integer.valueOf(splitLine[1]));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static Integer getValue(String key) throws Exception {
        if (!resources.containsKey(key)) {
            throw new Exception("Invalid resource name.");
        }
        return resources.get(key);
    }
}
