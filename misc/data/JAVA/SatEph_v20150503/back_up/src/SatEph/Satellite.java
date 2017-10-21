/*
 * Satellite.java
 *
 * Created on February 27, 2008, 8:48 PM
 *
 *Sk: Konstruktor sluzi na ukladanie dat umelych satelitov Zeme
 *En: Constructor for datas of artificial satellites of Earth
 */

//package
package SatEph;

/**
 *
 * @author Jiri Silha
 */
public class Satellite {
    //Identification datas
    /**
     * NORAD catalogue number
     **/
    public int catalogNo = 0;               
    /**
     * Name of object
     **/
    public String name = "No name";
    /**
     * International designator
     **/
    public String intDes = "No id";
    
    //TLE datas
    /**
     * 1st line of TLE datas
     **/
    public String line1 = "N/A";            
     /**
     * 2nd line of TLE datas
     **/
    public String line2 = "N/A";

    /**
     * 1st report line
     */
    public String line1Rep = "N/A";

     /**
     * 2nd report line
     */
    public String line2Rep = "N/A";
    
     /**
     * Element set epoch
     **/
    //public double epoch = 0;
    //POZOR - ine ako double moze sposobit vacsie nepresnosti
     /**
     * Element set epoch in MJD
     **/
    public double epochMjd = 0;             
     /**
     * 1st derivat of mean motion
     **/
    //public float der1mm = 0;
     /**
     * 2nd derivat of mean motion
     **/
    //public float der2mm = 0;
     /**
     * B*Drag term
     **/
    public float bDrag = 0;                 
     /**
     * Element number
     **/
    public int elementNo = 0;               
     /**
     * Revolution Number at Epoch
     **/
    public int revNoEp = 0;                  
     /**
     * Orbit inclination (degrees)
     **/
    public float incl = 0;
     /**
     * Right Ascension of Ascending Node (degrees)
     **/
    public float node = 0;
     /**
     * Argument of Perigee (degrees)
     **/
    public float argument = 0;
     /**
     * Mean Anomaly (degrees)
     **/
    public float meanAnom = 0;
     /**
     * Eccentricity
     **/
    public float ecc = 0;
     /**
     * Mean Motion (revolutions/day)
     **/
    public float meanMot = 0;
     /**
     * Radar Cross Section
     **/
    //other datas
    public float rcs = 0;                   
     /**
     * Organisation
     **/
    public String source = "";
    //public int source = 0;
     /**
     * Date of launch
     **/
    public String launched = "";           
     /**
     * Date of decayed
     **/
    public String decayed = "On orbit";
     /**
     * Is the object on Earth orbit?
     **/
    public boolean isOnOrbit = true;
     /**
     * Is the object orbital debris? (all debris object except R/B)
     * Debris objects are considered to be any object with the string 'DEB',
     * or 'COOLANT' or 'SHROUD' or 'WESTFORD NEEDLES' in the SatCat object common name.
      **/
    public boolean isOrbitalDebris = false;                                                //or 'COOLANT' or 'SHROUD' or 'WESTFORD NEEDLES' in the SatCat object common name.  
     /**
     * Is the object rocket body (R/B)?
     * Rocket Bodies are considered to be any object which contains the strings 'R/B' or
     * 'AKM' (Apogee Kick Motor) or 'PKM' (Perigee Kick Motor) but not 'DEB' in the SatCat object common name.
     **/
    public boolean isRB = false;
     /**
     * above Earth surface [km]
     **/                                           
    //public double apogee = 0;
     /**
     * above Earth surface [km]
     **/
    //public double perigee = 0;
     /**
     *  semi major axis [m]
     **/
    public float smAxis = 0;
     /**
     * period [min]
     **/
    public float period = 0;
     /**
     * array with actual datas for ephermeris
     **/
    //public SatEphData []array = new SatEphData[1];
    public SatEphData array = new SatEphData();

    //Osculating elements
    /**
     * Osculating orbit inclination (degrees)
     **/
    //public float oscIncl = 0;
     /**
     * Osculating right Ascension of Ascending Node (degrees)
     **/
    //public float oscNode = 0;
     /**
     * Osculating argument of Perigee (degrees)
     **/
    //public float oscArgument = 0;
     /**
     * Osculating mean Anomaly (degrees)
     **/
    //public float oscMeanAnom = 0;
     /**
     * Osculating eccentricity
     **/
    //public float oscEcc = 0;
     /**
     * Osculating mean Motion (revolutions/day)
     **/
    //public float oscMeanMot = 0;
    
    
    /** constructor*/
    public Satellite() {
        
    }
    
}
