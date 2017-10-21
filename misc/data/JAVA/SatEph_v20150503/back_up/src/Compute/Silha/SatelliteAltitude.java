/*
 * Methods to get the satellite - observer distance. Based on work J.Silha (2007)
 */

package Compute.Silha;

// 07/08/07 - Jiri Silha

/*Program sluzi na vypocet vysky telesa na zaklade pozorovanej uhlovej rychlosti. Vyuzivaju sa vztahy z prace:
Bakalarska praca: Uhlove rychlosti a pocetnost satelitov na oblohe - Jiri Silha (2007) (BP)
*/

//deklaracia triedy
public class SatelliteAltitude{
	
	/*
	Deaklaracia premennych v programe
	*/
	
	//deklaracia triedy CubicFormula, pomocou jej metody getRoot zistime koren nasej kubickej rovnice
	CubicFormula cb=new CubicFormula();
	
	
	//Konstanty
	
	//gravitacna konstanta G, Montebruck,Gill (2000); [m3/kg/s2]
	double Gconst=6.67259E-11;
	//gravitacna konstanta G, Montebruck,Gill (2000); [km3/kg/s2]
	double Gconst2=6.67259E-20;
	
	//hmotnost Zeme M; Montebruck, Gill (2000), [kg]
	double Mearth=5.9743E24;
	//polomer Zeme, zdroj NASA (rovnik=6378.1 km, poly=6356.8 km),
	double Rearth=6378.137e3;		//[m]
	
	//Vstupujuce premenne
	double time;		//expozicna doba
	
	double beta;		//uhlova vzdialenost, ktoru teleso preslo na snimke, z pohladu pozorovatela [stupne]
	double betarad;		//uhol beta v radianoch [rad]
	double alpha;		//uhlova vzdialenost, ktoru preslo teleso na kruhovej drahe [stupne]
	double altitude;		//vyska telesa nad povrchom [m]
	double altitude2;	//vzdialenost telesa od stredu Zeme; Rradius+altitude
	
	/** getAltitude()
         *  
         * Method to compute satellite altitude - id depends on angular speed of object
         * metoda sluzi na vypocet vysky telesa [m] v zavisloti od prejdenej drahy beta za expozicnu dobu time
         *
         * IN:
         *  double beta - angle [rad]
         *  double time - exp time [s]
         * 
         * OUT:
         *  double altitude - altitude of object [m]
         */
	public double getAltitude(double beta, double time){
		//prevod bety na rad
		betarad=beta;
		//na zaklade vztahu (12) z Bak.prace (dalej iba BP) definujem koeficienty a,b,c,d, (a*x^3 + b*x^2 + c*x + d = 0), hladame x (resp. vysku telesa altitude)
		double a=1;
		double b=-2*Rearth;
		double c=Rearth*Rearth;
		double d=time*time*Gconst*Mearth/4*(Math.cos(betarad)+1)/(Math.cos(betarad)-1); //cot(beta/2) = 1/4 * (cos(beta)-1)/(cos(beta)+1)
		altitude=cb.getRoot(a,b,c,d);
                
                //System.out.println("root " + altitude/1000);
		
		//treba hodnotu altitude osetrit o pripad, kedy je teleso vo vyske pod 2200 km, kedy metoda pre vypocet kubickych polynomov nema riesenie
		if((""+altitude).equals("NaN")){
			altitude=Rearth+getLowAltitude(beta, time);
		}
			
                //System.out.println("alt1 " + altitude/1000);
		return altitude-Rearth;
		
	}
	
        /** getLowAltitude()
         *  
         * Method to compute satellite altitude - id depends on angular speed of object
	 * metoda sluzi na vypocet vysky [m] v pripade, ze je altitude >2200 km, 
         *
         * IN:
         *  double beta - angle [rad]
         *  double time - exp time [s]
         * 
         * OUT:
         *  double altitude - altitude of object [km]
         */
	public double getLowAltitude(double beta, double time){
		//premenna sluzi na zistenie najmensieho rozdielu vypoctitanej bety od pozorovanej bety
		double deltabeta;
		//pomocna deltabeta
		double deltabeta2=2200000;
		//vypocitana beta v zavislosi na vyske
		double betapoz;
		//vysledna beta
		double lowaltitude=0;
		//for cyklus kontroluje, ktora hodnota vypocitanej bety sa najviac priblizuje k nami pozorovanej bete
		//zaciname od vysky 100 k  nad povrchom
		for(double i=100000; i<=2200e3; i=i+1000){
			//double i=2200000;
			//vypocet betapoz pomocou vztahu (9) a (11) z BP
			betapoz=Math.acos((4*i*i*(i+Rearth)-(time*time*Gconst*Mearth))/(4*i*i*(i+Rearth)+time*time*Gconst*Mearth));
			//prevod na stupne
			betapoz=Math.toDegrees(betapoz);
			//porovnavame bety
			deltabeta=Math.abs(beta-betapoz);
			//ak je deltabeta mensia, ako ta pred tym prirad novu betu, tu krora je blizsie
			if(deltabeta<deltabeta2) {
				deltabeta2=deltabeta;
				lowaltitude=i;
			}
				
		}
		//System.out.println("alt2 " + lowaltitude/1000);
		return lowaltitude;
	}
	
}
