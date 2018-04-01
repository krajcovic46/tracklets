package com.skrajcovic.datastructures;

public class Declination {
    private int degrees;
    private double minutes;
    private double seconds;

    private double degValue;

    public Declination(Integer degrees, Integer minutes, Double seconds) {
        setDegrees(degrees);
        setMinutes(minutes);
        setSeconds(seconds);

        convertToDegrees();
    }

    public int getDegrees() {
        return degrees;
    }

    public void setDegrees(int degrees) {
        this.degrees = degrees;
    }

    public double getMinutes() {
        return minutes;
    }

    public void setMinutes(double minutes) {
        this.minutes = minutes;
    }

    public double getSeconds() {
        return seconds;
    }

    public void setSeconds(double seconds) {
        this.seconds = seconds;
    }

    public double getDegValue() {
        return degValue;
    }

    private void convertToDegrees() {
        double temp = 0;

        temp += getDegrees();
        temp += getMinutes() / 60;
        temp += getSeconds() / 3600;

        degValue = temp;
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
