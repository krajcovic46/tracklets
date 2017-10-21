/*
*
*SatetVector.java
* 
*Ucel (purpose) :
*
*	Trieda sluzi na ulozenie aktualnej pozicie a rychlosti
*	Constructor for position and velocity
*
* 2007/11/09 - Jiri Silha
*
*/

package Compute.Silha;

//
// Deklaracia triedy (class declaration)
//

/** Position a velocity vectors*/
public class StateVector{
	
	public Vector r;// = new Vector(3);
        public Vector v;// = new Vector(3);
        
        
        //state vector - premenna je drzitelom informacii - pozicny vektor a vektor rychlosti  	
	public StateVector(){
        }
	
	//state vector - premenna je drzitelom informacii - pozicny vektor a vektor rychlosti  	
	public StateVector(Vector r, Vector v){
		this.r=r;
		this.v=v;
	}
}
