/*
 * Satellite.java
 *
 * Created on February 27, 2008, 8:48 PM
 *
 *En: Constructor for datas of stars
 */

//package
package SatEph;

/**
 *
 * @author Jiri Silha
 */
public class Star {
    //Identification datas               
    /**
     * Star ID
     **/
    public String name = "No name";
    /**
     * Right Ascension [rad]
     **/
    public double ra = 0.0;
    /**
     * Declination [rad]
     **/
    public double dec = 0.0;            
     /**
     * Magnitude
     **/
    public double mag = 0.0;
    /**
     * Declination [rad]
     **/
    public double elev = 0.0;            
     /**
     * Magnitude
     **/
    public double az = 0.0;    
    
    /** constructor*/
    public Star() {
        
    }
    
}
