package com.skrajcovic.datastructures;

public class Declination {
    private int degrees;
    private int minutes;
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

    //TODO
    public double add(Declination other) {
        double temp = getDegValue() + other.getDegValue();

        if (temp > 90) {
            temp = -90 + temp;
        }

        return temp;
    }

    //TODO
    public double sub(Declination other) {
        return 0;
    }

    public boolean isSmaller(Declination other) {
        return getDegValue() < other.getDegValue();
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
