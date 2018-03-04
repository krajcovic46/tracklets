package com.skrajcovic.datastructures;

public class Rectascension {
    private int hours;
    private int minutes;
    private double seconds;

    public Rectascension(Integer hours, Integer minutes, Double seconds) {
        setHours(hours);
        setMinutes(minutes);
        setSeconds(seconds);
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

    public boolean isSmaller(Rectascension other) {
        if (this.getHours() >= other.getHours()) {
            if (this.getMinutes() >= other.getMinutes()) {
                if (this.getSeconds() >= other.getSeconds()) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isLarger(Rectascension other) {
        return !isSmaller(other);
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
