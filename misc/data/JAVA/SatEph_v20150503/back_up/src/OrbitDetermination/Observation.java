/*
 * Constructor to hold informations about object position on the celestial
 * sphere.
 */

package OrbitDetermination;

/**
 *  Constructor
 * @author Jiri Silha - 22.10.2009
 */
public class Observation {

    //DATA
    /**
     * Time Modified Julian Date [day]
     */
    public double timeMjd;

    /**
     * Observer's position - longitude [rad], from East to West [0,2*PI]
     */
    public double lon;

     /**
     * Observer's position - latitude [rad], from South to North [-PI/2,PI/2]
     */
    public double lat;

     /**
     * Observer's position - altitude [m] 
     */
    public double alt;

     /**
     * Object's equatorial position on celestial sphere - R.A. [rad]
     */
    public double ra;

    /**
     * Object's equatorial position on celestial sphere - declination [rad]
     */
    public double dec;

    /**
     * Constructor to hold the observation informations.
     */
     public Observation(){
     }

     /**
     * Constructor to hold the observation informations.
     */
     public Observation(double timeMjd, double lon, double lat, double alt,
                double ra, double dec){
         this.timeMjd = timeMjd;
         this.lon = lon;
         this.lat = lat;
         this.alt = alt;
         this.ra = ra;
         this.dec = dec;
     }

     /**
     * Constructor to hold the observation informations.
     */
     //public Observation(double timeMjd, Geodetic observer,
     //           double ra, double dec){
     //    this.timeMjd = timeMjd;
     //    this.lon = observer.lon;
     //    this.lat = observer.lat;
     //    this.alt = observe.alt;
     //    this.ra = ra;
     //    this.dec = dec;
     //}
}
