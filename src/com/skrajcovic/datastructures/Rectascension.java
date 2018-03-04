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
        temp += getMinutes() / 60;
        temp += getSeconds() / 3600;

        degValue = temp;
    }

    public double add(Rectascension other) {
        double temp = this.getDegValue() + other.getDegValue();
        return temp % 24;
    }

    public double sub(Rectascension other) {
        double myDegValue = getDegValue();
        double otherDegValue = other.getDegValue();

        if (myDegValue > otherDegValue) {
            return myDegValue - otherDegValue;
        } else if (myDegValue < otherDegValue) {
            myDegValue = 24 + myDegValue;
            return myDegValue - otherDegValue;
        }
        return myDegValue;
    }


    public boolean isSmaller(Rectascension other) {
        return getDegValue() < other.getDegValue();
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
