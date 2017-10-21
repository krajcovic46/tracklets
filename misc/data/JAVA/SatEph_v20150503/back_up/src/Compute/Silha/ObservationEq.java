/*
*
*Observation.java
* 
*Ucel (purpose) :
*
*	Trieda sluzi na ukladanie udajov pozorovani (R.A, deklinacia, cas).
*	Constructor for observations (right acc., declination, time)
*
* 2007/10/30 - Jiri Silha
*
*/

package Compute.Silha;

/*
* Deklaracia triedy (class declaration)
*/

/** Observation properties*/
public class ObservationEq{
	
	//deklaracia premennych (variables)
	public double ra, dec;
	public Time t = new Time();
	
	//Konstruktor (constructor)
	public ObservationEq(double ra, double dec, Time t){
		this.ra = ra;
		this.dec = dec;
		this.t = t;
	}
}