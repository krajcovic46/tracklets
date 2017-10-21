/*
 * Silha - created november 2009.
 *
 * This class is to compute Moon geocentric position.
 *
 * Source is Montenbruck file SAT_Force.cpp 1999-2005
 */

package Compute.Silha;

/**
 *
 * @author Jiri Silha - november 2009
 */
public class Moon {

    static Constants constants = new Constants();

    public double az, h, ra, dec, distance;
    //Moon geocentric position vector
    //public Vector moonPosVector = new Vector(3);

    /** Creates a new instance of MoonPosition */
    public Moon() {
    }

     public Moon(double az, double h, double ra, double dec) {
         this.az = az;
         this.h = h;
         this. ra = ra;
         this.dec = dec;
    }

    /**
     * getMoonPosition()
     *
     * Method to compute geocentric position of the Moon.
     * Source class SAT_Force.cpp.
     *
     * INPUT:
     *  double mjd - modified julian date [day]
     *
     * OUTPUT:
     *  Vector moonPosVec - [m]
     * Montenbruck note:
     * Lunar position vector [m] with respect to the
     * mean equator and equinox of J2000 (EME2000, ICRF).
     */
    
    public static Vector getMoonPosition(double mjd){
        //Moon geocentirc position vector
        Vector moonPosVec = new Vector(3);

        // Constants
        double eps = Math.toRadians(23.43929111);   // Obliquity of J2000 ecliptic
        double mjd_J2000 = Time.getMjd(2000, 1, 1, 12, 0, 0.0);
        //System.out.println("mjd_J2000 " + mjd_J2000);
        double T   = (mjd-mjd_J2000)/36525.0;       // Julian cent. since J2000

        double arcs = 3600*180/Constants.pi;
        
        // Variables
        double  L_0, l,lp, F, D, dL, S, h, N;
        double  L, B, R, cosB;
        //Vector  r_Moon(3);

        // Mean elements of lunar orbit
        L_0 = getFrac ( 0.606433 + 1336.851344*T );     // Mean longitude [rev]
                                                           // w.r.t. J2000 equinox
        l   = Constants.pi2*getFrac ( 0.374897 + 1325.552410*T );     // Moon's mean anomaly [rad]
        lp  = Constants.pi2*getFrac ( 0.993133 +   99.997361*T );     // Sun's mean anomaly [rad]
        D   = Constants.pi2*getFrac ( 0.827361 + 1236.853086*T );     // Diff. long. Moon-Sun [rad]
        F   = Constants.pi2*getFrac ( 0.259086 + 1342.227825*T );     // Argument of latitude

        // Ecliptic longitude (w.r.t. equinox of J2000)
        dL = +22640*Math.sin(l) - 4586*Math.sin(l-2*D) + 2370*Math.sin(2*D) +  769*Math.sin(2*l)
               -668*Math.sin(lp) - 412*Math.sin(2*F) - 212*Math.sin(2*l-2*D) - 206*Math.sin(l+lp-2*D)
               +192*Math.sin(l+2*D) - 165*Math.sin(lp-2*D) - 125*Math.sin(D) - 110*Math.sin(l+lp)
               +148*Math.sin(l-lp) - 55*Math.sin(2*F-2*D);

        L = Constants.pi2 * getFrac( L_0 + dL/1296.0e3 );  // [rad]

        // Ecliptic latitude
        S  = F + (dL+412*Math.sin(2*F)+541*Math.sin(lp)) / arcs;
        h  = F-2*D;
        N  = -526*Math.sin(h) + 44*Math.sin(l+h) - 31*Math.sin(-l+h) - 23*Math.sin(lp+h)
             +11*Math.sin(-lp+h) - 25*Math.sin(-2*l+F) + 21*Math.sin(-l+F);

        B = ( 18520.0*Math.sin(S) + N ) / arcs;   // [rad]

        cosB = Math.cos(B);

        // Distance [m]
        R = 385000e3 - 20905e3*Math.cos(l) - 3699e3*Math.cos(2*D-l) - 2956e3*Math.cos(2*D)
              -570e3*Math.cos(2*l) + 246e3*Math.cos(2*l-2*D) - 205e3*Math.cos(lp-2*D)
              -171e3*Math.cos(l+2*D) - 152e3*Math.cos(l+lp-2*D);

        // Equatorial coordinates
        moonPosVec.v[0] =  R*Math.cos(L)*cosB;
        moonPosVec.v[1] =  R*Math.sin(L)*cosB;
        moonPosVec.v[2] =  R*Math.sin(B);
        Matrix rotationMatrix = new Matrix(3,3);
        rotationMatrix = rotationMatrix.R_x(-eps);
        moonPosVec = rotationMatrix.matrixMultiplyVector(rotationMatrix, moonPosVec);
        //r_Moon = R_x(-eps) * Vector ( R*cos(L)*cosB, R*sin(L)*cosB, R*sin(B) );

        return moonPosVec;
    }

    //method getFrac
    public static double getFrac(double value){
        double value2 = value - Math.floor(value);
	return value2;
    }

    /**
    *
    * getMoonLocalPositions: Computes the Moon's horizontal and equatorial coordinates
    * from observer view
    *
    * Input:   double MJD (days)
    *          Geodetic observer position [rad, rad, m]
    *
    *
    * Output:
    *
    *  MoonPosition [azimuth, altitude, ra, declination]
    */

    public Moon getMoonLocalPositions(double mjd, Geodetic obs){
        Moon mp = new Moon();

        Vector moonVector = new Vector(3);   //equatorial coordinates of the Moon
        Vector horizontal = new Vector(3);  //horizontal coordinates of object
        Vector equatorial = new Vector(3);  //equatorial coordinates of object from observer view
        Vector observatory = new Vector(3); //equatorial coordinates of observer

        moonVector = (new Moon()).getMoonPosition(mjd);

        horizontal=(new Transformation()).getHorizontalCoordinates3(moonVector,obs, mjd);
        //horizontal coordinates
        mp.az = horizontal.v[0];
        mp.h = horizontal.v[1];
        mp.distance = moonVector.getSize(moonVector);

        //equatorial (observer) coordinates
        //observer eq. coordinates
        observatory = (new Transformation()).fromGeodToGeoc(obs.lon, obs.lat,
                   obs.altitude, (new Constants()).R_Earth, (new Constants()).f_Earth);
        observatory = (new Transformation()).fromGeocToEquat(observatory,mjd);
        //equatorial coordinates of object
        equatorial = equatorial.subtractVectors(moonVector, observatory);
        equatorial = equatorial.fromCartesianToPolar(equatorial);
        //FILLING contructor
        mp.ra = equatorial.v[0];
        mp.dec = equatorial.v[1];

        return mp;
    }

    //test method
    public static void main(String args[]){
        Time time = new Time(2009,10,9,0,0,0.00);
        double mjd = time.getMjd(time);

        Vector moonVec = new Vector(3);
        moonVec = getMoonPosition(mjd);

        System.out.println(moonVec.v[0]/Constants.AU);
        System.out.println(moonVec.v[1]/Constants.AU);
        System.out.println(moonVec.v[2]/Constants.AU);
        System.out.println(moonVec.v[0]/1000);
        System.out.println(moonVec.v[1]/1000);
        System.out.println(moonVec.v[2]/1000);
        System.out.println(moonVec.getSize(moonVec)/1000);

    }
}
