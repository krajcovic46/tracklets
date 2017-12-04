/*
 * SunPosition.java
 *
 * Created on March 5, 2008, 3:56 AM
 *
 *  SK: Vypocet pozicie Slnka
 *  EN: Computation of position of Sun
 *
 * NOTES: source Montebruck
 */

package com.skrajcovic.orbitdetermination.compute;

/**
 *
 * @author Jiri Silha
 */

/**Computing actual Sun position, modification of Montenbruck's class SAT_Force.cpp*/
public class Sun {
    
    static Constants constants = new Constants();
    
    public double az, h, ra, dec;
    
    /** Creates a new instance of SunPosition */
    public Sun() {
    }
    
     public Sun(double az, double h, double ra, double dec) {
         this.az = az;
         this.h = h;
         this. ra = ra;
         this.dec = dec;
    }
    
    /*
    *
    * MiniSun: Computes the Sun's unit position vector using a low precision
    *          analytical series
    *
    * Input:
    *
    *   T         Time in Julian centuries since J2000
    *
    * Output:
    *
    *  Vector [x,y,z]	
    */
	public static Vector  getMiniSun (double T)
	{   //
            // Constants
            //
            double eps = Math.toRadians(23.43929111);
            
            //
            // Variables
            //
            double L,M;
            Vector  e_Sun;
            Matrix matrix = new Matrix(3,3);

            // Mean anomaly and ecliptic longitude
            M  = Math.PI*2 * getFrac ( 0.993133 + 99.997361*T);
            L  = Math.PI*2 * getFrac ( 0.7859453 + M/(Math.PI*2) +
                (6893.0*Math.sin(M)+72.0*Math.sin(2.0*M)+6191.2*T) / 1296.0e3);

            // Equatorial coordinates
            Vector vecL = new Vector(3);
            vecL.v[0] = L;

            vecL.v[1] = 0;
            vecL.v[2] = 0;
		
            e_Sun = matrix.matrixMultiplyVector(matrix.R_x(-eps), vecL.fromPolarToCartesian(vecL.v[0],vecL.v[1],1));
            return e_Sun;
	}
	
    /**
    *
    * getSunPosition: Computes the Sun's position vector using
    *
    * Input:
    *
    *   T         Time in MJD [days]
    *
    * Output:
    *
    *  Vector [x,y,z]	
    */
	public static Vector getSunPosition(double T){
            Vector sunPosition = new Vector(3);
            sunPosition = getMiniSun((T-51544.5)/365.25/100);
            sunPosition = sunPosition.multiplyVector(sunPosition, constants.AU);
            //System.out.println("Actual Sun Position " + sunPosition.getSize(sunPosition)/1000 + " km");
            return sunPosition;
        }
        
	//methor getFrac
	public static double getFrac(double value){
		double value2 = value - Math.floor(value);
		return value2;
	}
    
    /**
    *
    * getSunLocalPositions: Computes the Sun's horizontal and equatorial coordinates
    * from observer view    
    *
    * Input:   double MJD (days)
    *          Geodetic observer position [rad, rad, m]
    *    
    *   
    * Output:
    *
    *  SunPosition [azimuth, altitude, ra, declination]	
    */
    
    public Sun getSunLocalPositions(double mjd, Geodetic obs){
        Sun sp = new Sun();
        
        Vector sunVector = new Vector(3);   //equatorial coordinates of the Sun
        Vector horizontal = new Vector(3);  //horizontal coordinates of object
        Vector equatorial = new Vector(3);  //equatorial coordinates of object from observer view
        Vector observatory = new Vector(3); //equatorial coordinates of observer
        
        sunVector = (new Sun()).getSunPosition(mjd);
       
        horizontal=(new Transformation()).getHorizontalCoordinates3(sunVector,obs, mjd);
        //horizontal coordinates
        sp.az = horizontal.v[0];
        sp.h = horizontal.v[1];
        
        //equatorial (observer) coordinates
        //observer eq. coordinates
        observatory = (new Transformation()).fromGeodToGeoc(obs.lon, obs.lat,
                   obs.altitude, (new Constants()).R_Earth, (new Constants()).f_Earth);
        observatory = (new Transformation()).fromGeocToEquat(observatory,mjd);
        //equatorial coordinates of object
        equatorial = equatorial.subtractVectors(sunVector, observatory);
        equatorial = equatorial.fromCartesianToPolar(equatorial);
        //FILLING contructor
        sp.ra = equatorial.v[0];
        sp.dec = equatorial.v[1];
                                
        return sp;
    }

    /**
     * getAntihelionLocalPositions()
     *
     * Method to compute get Vector of Sun's opposition point for given time.
     *
     * INPUT:
     * double MJD (days)
     * Geodetic observer position [rad, rad, m]
     *
     * OUTPUT:
     * Sun sp
     *
     */

    public Sun getAntihelionLocalPositions(double mjd, Geodetic obs){
        Sun sp = new Sun();

        Vector sunVector = new Vector(3);   //equatorial coordinates of the Sun
        Vector horizontal = new Vector(3);  //horizontal coordinates of object
        Vector equatorial = new Vector(3);  //equatorial coordinates of object from observer view
        Vector observatory = new Vector(3); //equatorial coordinates of observer

        sunVector = sunVector.multiplyVector(new Sun().getSunPosition(mjd),-1);

        horizontal=(new Transformation()).getHorizontalCoordinates3(sunVector,obs, mjd);
        //horizontal coordinates
        sp.az = horizontal.v[0];
        sp.h = horizontal.v[1];

        //equatorial (observer) coordinates
        //observer eq. coordinates
        observatory = (new Transformation()).fromGeodToGeoc(obs.lon, obs.lat,
                   obs.altitude, (new Constants()).R_Earth, (new Constants()).f_Earth);
        observatory = (new Transformation()).fromGeocToEquat(observatory,mjd);
        //equatorial coordinates of object
        equatorial = equatorial.subtractVectors(sunVector, observatory);
        equatorial = equatorial.fromCartesianToPolar(equatorial);
        //FILLING contructor
        sp.ra = equatorial.v[0];
        sp.dec = equatorial.v[1];

        return sp;
    }
    
    /**
     * getAntihelionLocalPositionsForGEO()
     *
     * Method to compute get Vector of Sun's opposition point for given time and for GEO region.
     *
     * INPUT:
     * double MJD (days)
     * Geodetic observer position [rad, rad, m]
     *
     * OUTPUT:
     * Sun sp
     *
     */

    public Sun getAntihelionLocalPositionsForGEO(double mjd, Geodetic obs){
        Sun sp = new Sun();

        Vector sunVector = new Vector(3);   //equatorial coordinates of the Sun
        Vector horizontal = new Vector(3);  //horizontal coordinates of object
        Vector equatorial = new Vector(3);  //equatorial coordinates of object from observer view
        Vector observatory = new Vector(3); //equatorial coordinates of observer

        sunVector = sunVector.multiplyVector(new Sun().getSunPosition(mjd),-1);
        sunVector = Vector.multiplyVector(sunVector, (double)1/Vector.getSize(sunVector));
        sunVector = Vector.multiplyVector(sunVector, 42164150.0);

        horizontal=(new Transformation()).getHorizontalCoordinates3(sunVector,obs, mjd);
        //horizontal coordinates
        sp.az = horizontal.v[0];
        sp.h = horizontal.v[1];

        //equatorial (observer) coordinates
        //observer eq. coordinates
        observatory = (new Transformation()).fromGeodToGeoc(obs.lon, obs.lat,
                   obs.altitude, (new Constants()).R_Earth, (new Constants()).f_Earth);
        observatory = (new Transformation()).fromGeocToEquat(observatory,mjd);
        //equatorial coordinates of object
        equatorial = equatorial.subtractVectors(sunVector, observatory);
        equatorial = equatorial.fromCartesianToPolar(equatorial);
        //FILLING contructor
        sp.ra = equatorial.v[0];
        sp.dec = equatorial.v[1];

        return sp;
    }

    /**
     * test method
     */
    public static void main(String args[]){
        //get the antisolar maximum altitude for given night
        //1st of January 2011, 00:00
        double startMjd = 55562.0;
        double helpMjd = startMjd;
        //observatory
        Geodetic observatory = new Geodetic();
        observatory.lon = Math.toRadians(17.274);
        observatory.lat = Math.toRadians(48.373299);
        observatory.altitude = 531.0;

        //time difference
        double nightDuration = 0;
        double mainScanDuration = 0;
        double nightStartsMjd=0, nightEndsMjd=0;
        double nautNightStartsMjd=0, nautNightEndsMjd=0;
        double startScanningMjd=0, endScanningMjd=0;
        double maxOpositionAltitudeMjd = 0;
        double startScanningOposPointAltitude=0, stopScanningOposPointAltitude=0;
        boolean isNight = false;
        boolean isNauticalNight = false;
        //System.out.println(Time.getDateTime(startMjd).day + ", " + Time.getDateTime(startMjd).month + ", " + Time.getDateTime(startMjd).year);
        //System.out.println(Time.getDateTime(startMjd).hour + ", " + Time.getDateTime(startMjd).min + ", " + Time.getDateTime(startMjd).sec);

        //data holder
        Sun sp = new Sun();
        Sun spA = new Sun();
        //highest value of the oposition point altitude above horizont
        double oposPointMax = 0.0;
        for(int i = 0; i<=(365); i++){
            oposPointMax = 0.0;
            isNight = false;
            nautNightStartsMjd=0;
            nautNightEndsMjd=0;
            for(int j = 400; j<1600; j++){
                helpMjd = startMjd + (i+(float)j/(float)1000);
                //get the oposition point
                spA = sp.getAntihelionLocalPositions(helpMjd, observatory);
                if(spA.h>oposPointMax) {
                    oposPointMax = spA.h;
                    maxOpositionAltitudeMjd = helpMjd;
                }
                //get the sunset and sunrise time
                sp = sp.getSunLocalPositions(helpMjd, observatory);
                if((sp.h<0)&&(isNight==false)){
                    nightStartsMjd = helpMjd;
                    isNight = true;
                }
                if((sp.h>0)&&(isNight==true)){
                    nightEndsMjd = helpMjd;
                    isNight = false;
                }
                //get the nautical night start/end MJDs
                if((sp.h<-Math.PI/15)&&(isNauticalNight==false)){
                    nautNightStartsMjd = helpMjd;
                    startScanningMjd = helpMjd + 0.014;
                    //get the oposition altitude, when the main scanning will start, 20 minutes = 1/0.014 day after nautical night starts
                    spA = sp.getAntihelionLocalPositions(nautNightStartsMjd + 0.014, observatory);
                    startScanningOposPointAltitude = spA.h;
                    isNauticalNight = true;
                }
                if((sp.h>-Math.PI/15)&&(isNauticalNight==true)){
                    nautNightEndsMjd = helpMjd;
                    endScanningMjd = helpMjd - 0.014;
                    //get the oposition altitude, when the main scanning will stop, 20 minutes = 1/0.014 day before nautical night stops
                    spA = sp.getAntihelionLocalPositions(nautNightEndsMjd - 0.014, observatory);
                    stopScanningOposPointAltitude = spA.h;
                    isNauticalNight = false;
                }
            }
            //get the night duration in days
            nightDuration = nightEndsMjd - nightStartsMjd;
            mainScanDuration = endScanningMjd - startScanningMjd;
            
            System.out.println("Night of DOY " + (i) + "., " + nightStartsMjd + ", " + nautNightStartsMjd + ", " +
                    startScanningMjd + ", " + maxOpositionAltitudeMjd + ", " + endScanningMjd + ", " +
                    //maxOpositionAltitudeMjd + ", " + nautNightEndsMjd + ", " + nightEndsMjd + ", altitudes, " +
                    nautNightEndsMjd + ", " + nightEndsMjd + ", altitudes, " +
                    Math.toDegrees(startScanningOposPointAltitude) + ", " + Math.toDegrees(oposPointMax) + ", " +
                    Math.toDegrees(stopScanningOposPointAltitude) + ", duration, " + nightDuration*24 + ", " + mainScanDuration*24);
        }

        //TEST
        Sun sp3 = new Sun();
        sp3 = sp3.getAntihelionLocalPositions(55666.80300140381, observatory);
        System.out.println(Math.toDegrees(sp3.h));
    }
}
