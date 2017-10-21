/*
 * SatEphData.java
 *
 * Created on March 19, 2008, 4:52 AM
 *
 * SK: Konstruktor sluzi na ukladanie aktualnych dat telesapre dany cas
 * EN: Constructor is for satellite ephemeris datas
 */

package SatEph;

import Compute.Silha.*;

/**
 *
 * @author Jiri Silha
 */
public class SatEphData {
    
    /** Creates a new instance of SatEphData */
    public SatEphData() {
    }

    /**
     * Element set epoch in MJD
     **/
    public double epochMjd = 0;             

    /**
     * NORAD catalogue number
     */
    public int noradNo = 0;                    
    
    /**
     *  actual horizontal coordinates - azimuth [rad] <0,2PI>
     */
    public double azimuth = 0;          
    /**
     *  actual horizontal coordinates - altitude [rad] <-PI/2,PI/2>
     */
    public double altitude = 0;         
    
    /**
     * actual equatorial coordinates - right ascension [rad] <0,2PI>
     */
    public double ra = 0;
    /**
     * actual equatorial coordinates - declination [rad] <-PI/2,PI/2>
     */
    public double dec = 0;              
    
    /**
     * satellite range from the observer [m]
     */
    public float rangeObs = 0;
    /**
     * satellite range from the center of the Earth [m]
     */
    //public double rangeCenter = 0;
    
    /**
     * moving visual magnitude, magnitude computed from static magnitude
     */
    public float movingVisMag1 = 0;

    /**
     * static visual magnitude, magnitude of object with given RCS, distance, phase angle
     */
    public float staticVisMag1 = 0;

    /**
     * visual magnitude, formula RCS = f(OCS) Bodhwar 1992,[], RCS - Radar Cross Section, OCS - Optical Cross section
     */
    //public double visMag2 = 0;
    /**
     * absolute magnitude of the object, RCS = OCS, [], RCS - Radar Cross Section, OCS - Optical Cross section
     */
    public float absMag1 = 0;

    /**
     * What is the error magnitude
     */
    public float deltaMag1 = 0;
    
    /**
     * absolute magnitude of the object, formula RCS = f(OCS) Bodhwar 1992,[], RCS - Radar Cross Section, OCS - Optical Cross section
     */
    //public double absMag2 = 0;

    /**
     * actual mean anomaly [rad]
     */
    public float actMeanAnom = 0;
    /**
     * actual true anomaly [rad]
     */
    //public double actTrueAnom = 0;
    /**
     * actual geocentric velocity [m/s]
     */
    //public double velocity = 0;
    /**
     * actual phase angle [rad] <0,PI>
     */
    public float phaseAngle = 0;
    /**
     * is the object visible (is not in the shadow of the Earth?
     */
    public boolean isVisible = true;    
    /**
     * is  in the shadow of the Earth
     */
    public boolean isInShadow = false; 
    /**
     * actual angular speed [rad/s]
     */
    public float angularSpeed = 0;      
    /**
     * actual position angle [rad]
     */
    public float positionAngle = 0; 
    
    /**
     *Actual equatorial position X, Y, Z
     */
    public Vector positionVector = new Vector(3);
    
     /**
     *Actual equatorial velocity X, Y, Z
     */
    public Vector velocityVector = new Vector(3);
}
