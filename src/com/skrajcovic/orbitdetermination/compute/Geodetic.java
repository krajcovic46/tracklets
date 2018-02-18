/*
*
*Geodetic.java
* 
*Ucel (purpose) :
*
*	Trieda sluzi na ukladanie geodetickych suradnic (zem.dlzka, sirka a nadmorska vyska). 
*	Constructor for geodetic coordinates (lon, lat, altitude)
*
* 2007/10/11 - Jiri Silha
*
*/

package com.skrajcovic.orbitdetermination.compute;


/** Observatory position info*/


public class Geodetic{
	
	//deklaracia premennych (variables)
	public double lon, lat, altitude;
	
	/*
         *Konstruktor (constructor)
         **/
	public Geodetic(double lon, double lat, double altitude){
		this.lon=lon;
		this.lat=lat;
		this.altitude=altitude;
	}
        
        public Geodetic(){}
}