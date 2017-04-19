package com.skrajcovic;

public class FITSObject {
    private String name;
    private boolean real;
    private double x;
    private double y;
    private String mjd;
    private double intensity;

    public FITSObject() {}

    public FITSObject(String name, boolean real, String mjd, double x, double y, double intensity) {
        setName(name);
        setReal(real);
        setMjd(mjd);
        setX(x);
        setY(y);
        setIntensity(intensity);
    }

    public FITSObject(String[] data) {
        setName(data[0]);
        setReal(Boolean.parseBoolean(data[1]));
        setMjd(data[2]);
        setX(Double.valueOf(data[3]));
        setY(Double.valueOf(data[4]));
        setIntensity(Double.valueOf(data[5]));
    }

    public String toString() {
        return "[" + getName() + ", " + getX() + ", " + getY() + "]";
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isReal() {
        return real;
    }

    public void setReal(boolean real) {
        this.real = real;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public String getMjd() {
        return mjd;
    }

    public void setMjd(String mjd) {
        this.mjd = mjd;
    }

    public double getIntensity() {
        return intensity;
    }

    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }
}
