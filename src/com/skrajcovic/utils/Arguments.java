package com.skrajcovic.utils;

import java.util.ArrayList;
import java.util.List;

public class Arguments {

    private static ArrayList<String> arguments;

    public static ArrayList<String> getArguments() {
        if (arguments == null) {
            arguments = new ArrayList<>();

            arguments.add("lr");
            arguments.add("linear");

            arguments.add("lod");
            arguments.add("orbit");

            arguments.add("nn");
            arguments.add("neural");
        }
        return arguments;
    }

    private Arguments() {}
}
