/*
*
* Constants.java
*
*Ucel (purpose):
*
*	Trieda je drzitelom matematickych a astronomickych konstant v SI jednotkach. (Definition of astronomical and mathematical constants (in MKS units)
*
*Poznamka (notes):
*
*	Trieda je modifikacia triedy SAT_Const.h od O. Montenbruck, E. Gill (2005/04/14  OMO  Final version (2nd reprint)
*
* 2007/10/09 - Jiri Silha
*
*/

package com.skrajcovic.orbitdetermination.compute;

/*
* Deklaracia triedy (class declaration)
*/
/** Definition of astronomical and mathematical constants (in MKS units)*/
public class Constants{

	//
	//Definicia matematickyck konstant (definition of mathematical constants)
	//
		
	public static double pi     =   3.141592653589793238462643383279;//Math.PI;
	public static double pi2    =   2*3.141592653589793238462643383279;//Math.PI;
	public static double Rad    =   pi/180;	//pocet radianov na jednotku stupena
	public static double Deg    =   180/pi;	//pocet stupnov na jedenotku radian

	//
	//Zakldane konstanty (General constants)	
	//

	public static double MJD_J2000 = 51544.5;	// Modifikovany juliansky datum pre epochu J2000.0 (Modif. Julian Date of J2000.0)
		
	public static double AU        = 149597870000.0;      // Astronomicka jednotka[m] (Astronomical unit [m]); IAU 1976)
	public static double c_light   = 299792458.0;         // Rychlost svetla (Speed of light  [m/s]); IAU 1976)

	//
	// Fyzikalne parametre Zeme, Slnka a Mesiaca (Physical parameters of the Earth, Sun and Moon)
	//

	// Rovnikove polomery a splostenia (Equatorial radius and flattening)

	public static double R_Earth     =   6378.137e3;      // Polomer Zeme [m] (Radius Earth [m]); WGS-84)
	public static double f_Earth     = 1.0/298.257223563; // Splostenie Zeme (Flattening); WGS-84   
	public static double R_Sun       = 696000.0e3;        // Polomer Slnka [m](Radius Sun [m]); Seidelmann 1992
	public static double R_Moon      =   1738.0e3;        // Polomer Mesiaca [m](Radius Moon [m])

	// Rotacia Zeme - Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)

	public static double omega_Earth = 7.2921158553e-5;   // [rad/s]; Aoki 1982, NIMA 1997

	// Gravitacne konstanty (Gravitational coefficient)

	public static double GM_Earth    = 398600.4415e+9;    // [m^3/s^2]; JGM3
	public static double GM_Sun      = 1.32712438e+20;    // [m^3/s^2]; IAU 1976
	public static double GM_Moon     = GM_Earth/81.300587;// [m^3/s^2]; DE200
        public static double G_const     = 6.67259e-11;       // [m^3/kg/s^2]; http://www.physik.uni-wuerzburg.de/~rkritzer/grav.pdf


	// Tlak slnecneho ziarenia vo vzdialenosti 1 AU (Solar radiation pressure at 1 AU )

	public static double P_Sol       = 4.560E-6;          // [N/m^2] (~1367 W/m^2); IERS 96
        
        //Earth polar radius
        public static double Rpolar_Earth = 6356.750e3;      // Najmensi polomer Zeme [m] (medzi polmi)
}
	