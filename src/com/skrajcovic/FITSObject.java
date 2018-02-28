package com.skrajcovic;

import com.skrajcovic.datastructures.Declination;
import com.skrajcovic.datastructures.Rectascension;
import com.skrajcovic.datastructures.Type;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.time.LocalDateTime;

public class FITSObject implements Comparable<FITSObject> {
    private String fileName;
    private Type type;

    private Rectascension rectascension;
    private Declination declination;

    private double magnitude;

    private double x;
    private double y;

    private double mjd;

    private LocalDateTime ldt;

    private boolean real;

    public FITSObject() {}

    public FITSObject(String fileName, boolean real, double mjd, double x, double y, double magnitude) {
        setFileName(fileName);
        setReal(real);
        setMjd(mjd);
        setX(x);
        setY(y);
        setMagnitude(magnitude);
    }

    public boolean isWithinLineThreshold(SimpleRegression regression, double threshold) {
        double x = this.getX();
        double y = this.getY();
        double m = regression.getSlope();
        double c = regression.getIntercept();
        double b = (y > 0) ? 1 : -1;

        double distance = (-m * x + b * y - c) / Math.sqrt(Math.pow(m, 2) + Math.pow(b, 2));

        return Math.abs(distance) <= threshold;
    }

    public double calculateDeltaTime(FITSObject otherObject) {
//        System.out.println("MOJ MJD: " + this.getMjd());
//        System.out.println("MOJ LDT: " + this.getLdt());
//        System.out.println("IHR MJD: " + otherObject.getMjd());
//        System.out.println("IHR LDT: " + otherObject.getLdt());
//        System.out.println("DEL MJD: " + (this.getMjd() - otherObject.getMjd()));
//        System.out.println("----------------------------------------------------");
        return Math.abs(this.getMjd() - otherObject.getMjd());
    }

    public double calculateSpeed(FITSObject otherObject) {
        return Math.sqrt(Math.pow(this.getX() - otherObject.getX(), 2) + Math.pow(this.getY() - otherObject.getY(), 2))
                / calculateDeltaTime(otherObject);
    }

    public double getHeading(FITSObject otherObject) {
        double theta = Math.toDegrees(Math.atan2(otherObject.getX() - getX(), otherObject.getY() - getY()));
        return Math.abs(theta);
    }

    public boolean isWithinAngleThreshold(FITSObject otherObject, double angle, double threshold) {
        double heading = getHeading(otherObject);
        if (otherObject.isReal()) {
            System.out.println(otherObject.getFileName() + " - " + angle +" <= " + heading + " + " + threshold +
            " && " + angle + " >= " + heading + " - " + threshold);
        }
        return angle <= heading + threshold && angle >= heading - threshold;
    }

    @Override
    public int compareTo(FITSObject o) {
        return Double.valueOf(this.getMjd()).compareTo(o.getMjd());
    }




    public String toString() {
        return "[" + getFileName() + ", " + getX() + ", " + getY() + ", " + isReal() +"]";
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String name) {
        this.fileName = name;
    }

    public Type getType() {
        return type;
    }

    public void setType(String type) {
        switch (type) {
            case "R":
                this.type = Type.R;
                break;
            case "S":
                this.type = Type.S;
                break;
            case "H":
                this.type = Type.H;
                break;
            case "?":
                this.type = Type.UNKNOWN;
                break;
            default:
                break;
        }
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

    public double getMjd() {
        return mjd;
    }

    public void setMjd(double mjd) {
        this.mjd = mjd;
    }

    public double getMagnitude() {
        return magnitude;
    }

    public void setMagnitude(double magnitude) {
        this.magnitude = magnitude;
    }

    public Rectascension getRectascension() {
        return rectascension;
    }

    public void setRectascension(Rectascension rectascension) {
        this.rectascension = rectascension;
    }

    public Declination getDeclination() {
        return declination;
    }

    public void setDeclination(Declination declination) {
        this.declination = declination;
    }

    public LocalDateTime getLdt() {
        return ldt;
    }

    public void setLocalDateTime(LocalDateTime ldt) {
        this.ldt = ldt;
        convertLDTtoMJD();
    }

    private void convertLDTtoMJD() {
        long    MjdMidnight;
        double  FracOfDay;
        int     b;


        int month = getLdt().getMonthValue();
        int year = getLdt().getYear();
        int day = getLdt().getDayOfMonth();
        int hour = getLdt().getHour();
        int min = getLdt().getMinute();
        int sec = getLdt().getSecond();

        if (month <= 2) {
            month+=12; --year;
        }

        if ((10000L * year + 100L * month + day) <= 15821004L ) {
            b = -2 + ((year + 4716) / 4) - 1179;     // Julian calendar
        } else {
            b = (year / 400) - (year / 100) + (year / 4);  // Gregorian calendar
        }

        MjdMidnight = 365L * year - 679004L + b + (int)(30.6001 * (month + 1)) + day;
        FracOfDay   = (hour + min / 60.0 + sec / 3600.0) / 24.0;

        setMjd(MjdMidnight + FracOfDay);
    }
}