package com.skrajcovic.utils;

public class Declination {
    private int degrees;
    private int minutes;
    private double seconds;

    public Declination(int degrees, int minutes, double seconds) {
        setDegrees(degrees);
        setMinutes(minutes);
        setSeconds(seconds);
    }

    public int getDegrees() {
        return degrees;
    }

    public void setDegrees(int degrees) {
        this.degrees = degrees;
    }

    public int getMinutes() {
        return minutes;
    }

    public void setMinutes(int minutes) {
        this.minutes = minutes;
    }

    public double getSeconds() {
        return seconds;
    }

    public void setSeconds(double seconds) {
        this.seconds = seconds;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(getDegrees());
        sb.append(" ");
        sb.append(getMinutes());
        sb.append(" ");
        sb.append(getSeconds());
        return sb.toString();
    }
}
