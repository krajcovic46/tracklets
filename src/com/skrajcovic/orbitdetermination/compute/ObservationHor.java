/*
*
*Observation.java
* 
*Ucel (purpose) :
*
*	Trieda sluzi na ukladanie udajov pozorovani (azimut, vysku, cas). 
*	Constructor for observations (azimuth, elevation, time)
*
* 2007/10/30 - Jiri Silha
*
*/

package com.skrajcovic.orbitdetermination.compute;

/*
* Deklaracia triedy (class declaration)
*/

/** Observation properties*/
public class ObservationHor{
	
	//deklaracia premennych (variables)
	public double azim, elev;
	public Time t;
	
	//Konstruktor (constructor)
	public ObservationHor(double azim, double elev, Time t){
		this.azim=azim;
		this.elev=elev;
		this.t=t;
	}

        //Konstruktor (constructor)
	public ObservationHor(){
	}
}