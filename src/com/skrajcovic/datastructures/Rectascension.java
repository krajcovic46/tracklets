package com.skrajcovic.datastructures;

public class Rectascension  {
    private int hours;
    private int minutes;
    private double seconds;

    private double degValue;

    public Rectascension(Integer hours, Integer minutes, Double seconds) {
        setHours(hours);
        setMinutes(minutes);
        setSeconds(seconds);

        convertToDegrees();
    }

    public int getHours() {
        return hours;
    }

    public void setHours(int hours) {
        this.hours = hours;
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
        temp += getHours() * 15;
        temp += getMinutes() * 15 / 60;
        temp += getSeconds() * 15 / 3600;

        degValue = temp;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(getHours());
        sb.append(" ");
        sb.append(getMinutes());
        sb.append(" ");
        sb.append(getSeconds());
        return sb.toString();
    }
}
