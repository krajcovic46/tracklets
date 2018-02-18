/*
*
*Kepler.java
* 
*Ucel (purpose) :
*
*	Vypocet keplerovskych drahovych elementov (Keplers elements computation)
*
*Poznamka (notes):
*
*	SOURCE: Trieda je modifikacia triedy SAT_Kepler.h od O. Montenbruck, E. Gill (2005/04/14  OMO  Final version (2nd reprint))
*
* 2007/10/10 - Jiri Silha
*
*/

package com.skrajcovic.orbitdetermination.compute;

//
// Deklaracia triedy (class declaration)
//
/** Class - Kepler motion */
 
public class Kepler{
	
	//
	//Deklaracia premennych pre konstruktor (variables declaration)
	//
	public double omega,Omega,incl,e,a,M;
        
        /**
         * Actual mean anomaly
         */
        public double M_actual;
	
	//
	//Konstanty (constants)
	//
	public static Constants consta = new Constants();
	
	public static double minimum = 1e-20;

        /**
         * Epoch MJD [days]
         */
        public double epochMJD;
	
        /**
         *Constructor Kepler()
         */
	public Kepler(double omega, double Omega, double incl, double e, double a, double M){
		this.omega=omega;
		this.Omega=Omega;
		this.incl=incl;
		this.e=e;
		this.a=a;
		this.M=M;
	}
	 /**
         *Constructor Kepler()
         */
	public Kepler(){}
	
	/*
	KEPLEROVSKY POHYB (KEPLERIAN ORBIT)
	*/
	
	/**
	*
	*   Vypocet drahovych elementov na zaklade 2 pozicnych vektorov v danych casoch 
	*   (Computing orbital elements from two given position vectors and  associated times )
	*
	* Input/Output:
	*
	*   GM        Gravitational coefficient
	*             (gravitational constant * mass of central body)
	*   Mjd_a     Time t_a (Modified Julian Date)
	*   Mjd_b     Time t_b (Modified Julian Date)
	*   r_a       Position vector at time t_a
	*   r_b       Position vector at time t_b
	*
	*   <return>  Keplerian elements (a,e,i,Omega,omega,M)
	*               a      Semimajor axis 
	*               e      Eccentricity 
	*               i      Inclination [rad]
	*               Omega  Longitude of the ascending node [rad]
	*               omega  Argument of pericenter  [rad]
	*               M      Mean anomaly  [rad]
	*             at time t_a 
	*
	* Notes:
	*
	*   The function cannot be used with state vectors describing a circular
	*   or non-inclined orbit.
	*
	*/
	
	public Kepler getElements(double GM, double Mjd_a, double Mjd_b, Vector r_a, Vector r_b){
		Kepler kepler;//=new Kepler(0,0,0,0,0,0);
		
		// Variables
  
		double  tau, eta, p;
		double  n, nu, E, u;
		double  s_a, s_b, s_0, fac, sinhH;
		double  cos_dnu, sin_dnu, ecos_nu, esin_nu;
		double  a, e, i, Omega, omega, M;
		Vector  e_a, r_0, e_0, W;

		//Vypocet vektora r_0, ktory je kolmy na vektor r_a a vypocet velkosti vektorov r_a,r_b a r_0
		//(Calculate vector r_0 (fraction of r_b perpendicular to r_a)  and the magnitudes of r_a,r_b and r_0)

		s_a = r_a.getSize(r_a);  
		e_a =r_a.multiplyVector(r_a,(double)1/s_a);// new Vector (r_a.x/s_a,r_a.y/s_a,r_a.z/s_a);
		
		s_b = r_b.getSize(r_b); 
		fac = r_b.getScalarProduct(r_b,e_a); 
		r_0 = r_b.subtractVectors(r_b,e_a.multiplyVector(e_a, fac));
		s_0 = r_0.getSize(r_0);  
		e_0 = r_0.multiplyVector(r_0,(double)1/s_0); //new Vector (r_0.x/s_0,r_0.y/s_0,r_0.z/s_0);
		
		// Inclination and ascending node 
		
		W = r_a.getVectorProduct(e_a,e_0);
		Omega = getOmega(W);                     // Long. ascend. node
		i = getInclination(W); // Inclination        
		u=getArgOfLatitude(W,r_a,e_a,i);  
		// Semilatus rectum
		 
		tau = getTau(GM,Mjd_a,Mjd_b);
		eta = getEta( r_a, r_b, tau );
		p   = Math.pow ( s_a*s_0*eta/tau, 2 );   
		
		// Eccentricity, true anomaly and argument of perihelion
		e  = getEccentricity(r_a, r_b, e_a, p);
		nu = getTrueAnomaly(r_a, r_b, e_a, p);
		omega = getArgOfPerigee(u,nu);
		
		// Perihelion distance, semimajor axis and mean motion
		  
		a = getSMAxis(p,e);

		n = getMeanMotion(GM,a);
		// Mean anomaly and time of perihelion passage
		
		//pomocne premenne
		cos_dnu = fac / s_b;    
		sin_dnu = s_0 / s_b;

		ecos_nu = p / s_a - 1.0;  
		esin_nu = ( ecos_nu * cos_dnu - (p/s_b-1.0) ) / sin_dnu;
		if (e<1.0) {
			E = Math.atan2 ( (Math.sqrt((1.0-e)*(1.0+e)) * esin_nu),( ecos_nu + e*e ));
			this.M = ( E - e*Math.sin(E))%consta.pi2;
			while(this.M < 0) {
                            this.M = this.M + consta.pi2;
                        }
		}
		else {
			sinhH = Math.sqrt((e-1.0)*(e+1.0)) * esin_nu / ( e + e * ecos_nu );
			this.M = e * sinhH - Math.log ( sinhH + Math.sqrt(1.0+sinhH*sinhH) );
			while(this.M < 0) {
                            this.M = this.M + consta.pi2;
                            //System.out.println("M " + this.M);
                            if(this.M < 10e3) break;
                        }
		}

                //EMPIRICAL CORRECTION, TEST
                //this.M = this.M + Math.toRadians(0.154);
                this.M = this.M;
		kepler = new Kepler(omega,Omega,i,e,a,this.M);
		//System.out.println("a: "+kepler.a/1000+" km");
		//System.out.println("e: "+kepler.e);
		//System.out.println("i: "+Math.toDegrees(kepler.incl)+" °");
		//System.out.println("Omega: "+Math.toDegrees(kepler.Omega)+" °");
		//System.out.println("omega: "+Math.toDegrees(kepler.omega)+" °");
		//System.out.println("Ma: "+Math.toDegrees(kepler.M)+" °\n");
		
		return kepler;
	}
	
	/**
	*
	*   Vypocet polohoveho vektora a vektora okamzitej rychlosti z drahovych elementov a aktualnej epochy
	*   (Computes  position and velocity vectors from actual orbital elements)
	*
	* Input/Output:
	*
	*   GM        Gravitational coefficient
	*             (gravitational constant * mass of central body)
	*   Kep       Keplerian elements (a,e,i,Omega,omega,M) with
	*               a      Semimajor axis [m]
	*               e      Eccentricity 
	*               i      Inclination [rad]
	*               Omega  Longitude of the ascending node [rad]
	*               omega  Argument of pericenter  [rad]
	*               M      Mean anomaly at epoch [rad]
	*   dt        Time since epoch [s]
	*   <return>  State vector (x,y,z,vx,vy,vz) [m,m.s-1]
	*
	* Notes:
	*
	*   The semimajor axis a=Kep(0), dt and GM must be given in consistent units, 
	*   e.g. [m], [s] and [m^3/s^2]. The resulting units of length and velocity  
	*   are implied by the units of GM, e.g. [m] and [m/s].
	*
	*/
	
	public StateVector getStateVector(double GM, Kepler kepler, double dt){

		// Variables
		double  a,e,i,Omega,omega,M,M0,n;
		double  E,cosE,sinE, fac, R,V;
		Vector  r=new Vector(3);
		Vector  v=new Vector(3);
		Matrix  PQW=new Matrix(3,3);
		StateVector stateVector=new StateVector(r,v);

		// Keplerian elements at epoch
		a = kepler.a;  Omega =kepler.Omega;
		e = kepler.e;  omega = kepler.omega; 
		i = kepler.incl;  M0    = kepler.M;

		// Mean anomaly

		if (dt==0.0) {
			M = M0;
		}
		else {
			n = Math.sqrt (GM/(a*a*a));
			M = M0 +n*dt;
			//System.out.println("\n"+"Actual M: "+Math.toDegrees(M)+" �\n");
		};
                
                //System.out.println("Mena anom " + Math.toDegrees(M));
                this.M_actual = M;
                
		// Eccentric anomaly
		E  = getEccAnom(M,e);               

		cosE = Math.cos(E); 
		sinE = Math.sin(E);

		// Perifocal coordinates
		fac = Math.sqrt ((1.0-e)*(1.0+e));  

		R = a*(1.0-e*cosE);  // Distance
		V = Math.sqrt(GM*a)/R;    // Velocity

		r.v[0]= a*(cosE-e);
		r.v[1]=a*fac*sinE;
		v.v[0]=-V*sinE;
		v.v[1]=+V*fac*cosE; 
		
		// Transformation to reference system (Gaussian vectors)
		PQW = PQW.getMatrixProduct(PQW.getMatrixProduct(PQW.R_z(-Omega),PQW.R_x(-i)),PQW.R_z(-omega));

		stateVector.r = PQW.matrixMultiplyVector(PQW,r);
		stateVector.v = PQW.matrixMultiplyVector(PQW,v);

		// State vector 
		return stateVector;
	}
	
	/*
	VYPOCET ELEMENTOV
	*/
	
	/**
	*
	*  Vypocet dlzky vystupneho uzla ( computing longitude of the ascending node [rad])
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector W (Gaussian vector)	  
	*
	*   <navratna hodnota> (<return>) Omega
	*
	*/
	
	public static double getOmega(Vector W){
		double Omega;
		Omega = Math.atan2 ( W.v[0],(-W.v[1]) );                     // Long. ascend. node 
		//System.out.println("Omega in method1: " + Math.toDegrees(Math.atan2 ( W.v[0],(-W.v[1]) )));
		while(Omega < 0.0) {
                    Omega = Omega + consta.pi2;
                    //System.out.println("Omega in method2a: " + Math.toDegrees(Omega));
                }
		while(Omega > consta.pi2) {
                    Omega = Omega - consta.pi2;
                    //System.out.println("Omega in method2b: " + Math.toDegrees(Omega));
                }
                //System.out.println("Omega in method2: " + Math.toDegrees(Omega));
		return Omega;
	}
	
	/*
	*   Vypocet argumentu sirky pericentra ( computing argument of pericenter  [rad])
	*
	* Vstup/Vystup (Input/output):
	*
	*   Zemepisna sirka (argument of latitude)
	*   Prava anomalia	(true anomaly)  
	*
	*   <navratna hodnota> (<return>) omega
	*
	*/
	/*
	public static double getArgument(double u, double v){
		double omega;
		omega=(u-v)%consta.pi2;
		while(omega <0.0) omega=omega+consta.pi2;
		while(omega >consta.pi2) omega=omega-consta.pi2;
		return omega;
	}
	*/
        
	/**
	*
	*   Vypocet sklonu drahy ( computing inclination  [rad])
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector W (Gaussian vector)	  
	*
	*   <navratna hodnota> (<return>) inclination 
	*
	*/
	
	public static double getInclination(Vector W){
		double incl;
		incl = Math.atan ( Math.sqrt(W.v[0]*W.v[0]+W.v[1]*W.v[1])/(W.v[2])); // Inclination   
		while(incl < 0) incl = incl + consta.pi;
		while(incl > consta.pi) incl = incl - consta.pi;
		return incl;
	}
	
	/**
	*
	*   Vypocet excentricity drahy (vztah 2.116) ( computing eccentricity (formula 2.116))
	*
	* Vstup/Vystup (Input/output):
	*   Vector r_a
	*   Vector r_b
	*   Vector e_a
	*   p - parameter (semi latus rectum)
	*
	*   <navratna hodnota> (<return>) Excentricita (eccentricity)
	*
	*/
	
	public static double getEccentricity(Vector r_a, Vector r_b, Vector e_a, double p){
		double ecc;
		double cos_dnu, sin_dnu, fac, s_a, s_b, s_0, ecos_nu, esin_nu;
		Vector r_0, e_0;
		
		s_a = r_a.getSize(r_a);  
		e_a =r_a.multiplyVector(r_a,(double)1/s_a);// new Vector (r_a.x/s_a,r_a.y/s_a,r_a.z/s_a);
		
		s_b = r_b.getSize(r_b); 
		fac = r_b.getScalarProduct(r_b,e_a); 
		r_0 = r_b.subtractVectors(r_b,e_a.multiplyVector(e_a, fac));
		s_0 = r_0.getSize(r_0);  
		e_0 = r_0.multiplyVector(r_0,(double)1/s_0); //new Vector (r_0.x/s_0,r_0.y/s_0,r_0.z/s_0);
		  
		// Eccentricity
		cos_dnu = fac / s_b;    
		sin_dnu = s_0 / s_b;

		ecos_nu = p / s_a - 1.0;  
		esin_nu = ( ecos_nu * cos_dnu - (p/s_b-1.0) ) / sin_dnu;

		ecc  = Math.sqrt ( ecos_nu*ecos_nu + esin_nu*esin_nu );
		
		return ecc;
	}

	/**
	*
	*   Vypocet pravej anomalie (computing the true anomaly)
	*
	* Vstup/Vystup (Input/output):
	*   Vector r_a
	*   Vector r_b
	*   Vector e_a (jednotkovy vektor z r_a)
	*   p - sprievodic (semi latus rectum)
	*
	*   <navratna hodnota> (<return>) Prava namoalia (true anomaly)
	*
	*/
	
	public static double getTrueAnomaly(Vector r_a, Vector r_b, Vector e_a, double p){
		double nu;
		double cos_dnu, sin_dnu, fac, s_a, s_b, s_0, ecos_nu, esin_nu;
		Vector r_0, e_0;
		
		s_a = r_a.getSize(r_a);  
		e_a =r_a.multiplyVector(r_a,(double)1/s_a);// new Vector (r_a.x/s_a,r_a.y/s_a,r_a.z/s_a);
		s_b = r_b.getSize(r_b); 
		fac = r_b.getScalarProduct(r_b,e_a); 
		r_0 = r_b.subtractVectors(r_b,e_a.multiplyVector(e_a, fac));
		s_0 = r_0.getSize(r_0);  
		e_0 = r_0.multiplyVector(r_0,(double)1/s_0); //new Vector (r_0.x/s_0,r_0.y/s_0,r_0.z/s_0);
		 
		cos_dnu = fac / s_b;    
		sin_dnu = s_0 / s_b;

		ecos_nu = p / s_a - 1.0;  
		esin_nu = ( ecos_nu * cos_dnu - (p/s_b-1.0) ) / sin_dnu;

		nu = Math.atan2(esin_nu,ecos_nu);
		
		while(nu>consta.pi2) nu=nu-consta.pi2;
		while(nu<0) nu=nu+consta.pi2;

		return nu;
	}
	
	/**
	*
	*  Argument sirky perigea (vztah 2.117) (computing the argument of perigee (formula 2.117))
	*
	* Vstup/Vystup (Input/output):
	*   u  -  (argument of latitude)
	*   nu - prava anomalia (true anomaly)
	*
	*   <navratna hodnota> (<return>) Argument sirky perigea (argument of perigee)
	*
	*/
	
	public static double getArgOfPerigee(double u, double nu){
		double omega;
		omega=(u-nu)%consta.pi2;
		while(omega <0.0) omega=omega+consta.pi2;
		while(omega >consta.pi2) omega=omega-consta.pi2;
		return omega;
	}

	/**
	*
	*  Velka poloos (vztah 2.118) (computing the semi-major axis (formula 2.118))
	*
	* Vstup/Vystup (Input/output):
	*   p - sprievodic (semi-latus rectum)
	*   e - excentricita (eccentricity)
	*
	*   <navratna hodnota> (<return>) Velka poloos (semi-major axis)
	*
	*/
        
	public static double getSMAxis(double p, double e){
		double a;
		a=p/(1.0-e*e);
		return  a;
	}
	
	/**
	*   Velka poloos (vztah 2.61) (computing the semi-major axis (formula 2.61))
	*
	* Vstup/Vystup (Input/output):
	*   GM        Gravitational coefficient
	*             (gravitational constant * mass of central body)
	*   n - stredny denny pohyb (mean motion), [rad.s-1]
	*
	*   <navratna hodnota> (<return>) Velka poloos (semi-major axis) 
	*
	*/
	
	public static double getSMAxis2(double GM, double n){
		double a;
		a=Math.pow(GM/(n*n),(double)1/3);
		return  a;
	}
	
	/**
	*
	*   Stredny denny pohyb (computing the mean motion)
	*
	* Vstup/Vystup (Input/output):
	*  GM        Gravitational coefficient
	*             (gravitational constant * mass of central body)
	*  a - velka poloos (semi-major axis)
	*
	*   <navratna hodnota> (<return>) Stredny denny pohyb (mean motion)
	*
	*/
	
	public static double getMeanMotion(double GM, double a){
		double n;
		n = Math.sqrt ( GM / Math.abs(a*a*a) );
		return n;
	}
	
	/*
	POMOCNE METODY (OTHER METHODS)
	*/
	
	/**
	*   Computing argument of latitude [rad])
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector W (Gaussian vector)
	*   Pozicny vektor (position vector)
	*   Jednotkovy vektor (unit vector)
	*   Sklon (inclination)	  
	*
	*   <navratna hodnota> (<return>) Argument zemepisnej sirky (argument of latitude)
	*
	*/
	
	public static double getArgOfLatitude(Vector W, Vector r_a, Vector e_a, double incl){
		double u;
		if (incl==0.0) u=Math.atan2( r_a.v[1],r_a.v[0]);
		else 	u=Math.atan2(e_a.v[2],(-e_a.v[0]*W.v[1]+e_a.v[1]*W.v[0]));
		
		return u;
	}
	
	/**
	*
	*   Vypocet premennej tau (vztah 2.99) ( computing variable tau (formula 2.99))
	*
	* Vstup/Vystup (Input/output):
	*
	*   GM        Gravitational coefficient
	*             (gravitational constant * mass of central body)
	*   Mjd_a     Time t_a (Modified Julian Date)
	*   Mjd_b     Time t_b (Modified Julian Date)
	*
	*   <navratna hodnota> (<return>) Argumnet zemepisnej sirky (argument of latitude)
	*
	*/
	
	public static double getTau(double GM, double Mjd_a, double Mjd_b){
		double tau;
		tau=Math.sqrt(GM) * 86400.0*Math.abs(Mjd_b-Mjd_a); 
		return tau;
	}
	
	/**
	*
	*  Vypocet funkcie f (vztaht 2.106) ( computing function f (formula 2.106))
	*
	* Vstup/Vystup (Input/output):
	*
	*  Eta (2.98)
	*  m (2.101)
	*  l  (2.101)
	*   <navratna hodnota> (<return>) Hodnota funkcie f (amount)
	*
	*/
	
	public static double getFunction(double eta, double m, double l){
		double function;
		
		// Constants
		double eps = 100.0 * minimum;

		// Variables
		double  w,W,a,n,g;
		    
		w = m/(eta*eta)-l; 

		if (Math.abs(w)<0.1) { // Series expansion
			W = a = 4.0/3.0; 
			n = 0.0;
			do {
				n += 1.0;  
				a *= w*(n+2.0)/(n+1.5);  
				W += a; 
			}
			while (Math.abs(a) >= eps);
		}
		
		else {
			if (w > 0.0) {
				g = 2.0*Math.asin(Math.sqrt(w));  
				W = (2.0*g - Math.sin(2.0*g)) / Math.pow(Math.sin(g), 3);
			}
			else {
				g = 2.0*Math.log(Math.sqrt(-w)+Math.sqrt(1.0-w));  // =2.0*arsinh(sqrt(-w))
				//W = (Math.sinh(2.0*g) - 2.0*g) / Math.pow(Math.sinh(g), 3);
                                W = (((Math.exp(Math.toRadians(2.0*g)) - 
					Math.exp(-Math.toRadians(2.0*g)))/2 - 2.0*g)) / 
                                        (Math.pow((Math.exp(Math.toRadians(g)) - 
					Math.exp(-Math.toRadians(g)))/2, 3));
			}
		}
		  
		function = ( 1.0 - eta + (w+l)*W );
		
		return function;
	}
	
	/**
	*
	*   Vypocet premennej eta (vztahty 2.98 - 2.108) ( computing variable eta (formulas 2.98 - 2.108))
	*
	* Vstup/Vystup (Input/output):
	*
	*   r_a   1.pos.vector
	*   r_b   2.pos.vector
	*   tau
	*
	*   <navratna hodnota> (<return>) Eta (formula 2.98)
	*
	*/
	
	public static double getEta(Vector r_a, Vector r_b, double tau){
		
		// Constants
		int maxit = 30;
		double delta = 100.0*minimum;  

		// Variables
		int    i;
		double kappa, m, l, s_a, s_b, eta_min, eta1, eta2, F1, F2, d_eta;

		// Auxiliary quantities

		s_a = r_a.getSize(r_a);  
		s_b = r_b.getSize(r_b);  

		kappa = Math.sqrt ( 2.0*(s_a*s_b+r_a.getScalarProduct(r_a,r_b)) );

		m = tau*tau / Math.pow(kappa,3);   
		l = (s_a+s_b) / (2.0*kappa) - 0.5;

		eta_min = Math.sqrt(m/(l+1.0));

		// Start with Hansen's approximation

		eta2 = ( 12.0 + 10.0*Math.sqrt(1.0+(44.0/9.0)*m /(l+5.0/6.0)) ) / 22.0;
		eta1 = eta2 + 0.1;   

		// Secant method
		 
		F1 = getFunction(eta1, m, l);   
		F2 = getFunction(eta2, m, l);  

		i = 0;

		while (Math.abs(F2-F1) > delta){
			d_eta = -F2*(eta2-eta1)/(F2-F1);  
			eta1 = eta2; 
			F1 = F2; 
			while (eta2+d_eta<=eta_min)  d_eta *= 0.5;
			eta2 += d_eta;  
			F2 = getFunction(eta2,m,l); 
			++i;
			if ( i == maxit ) {
				//System.out.println( "WARNING: Convergence problems in getEta");
			break;
			}
		  
		}
		return eta2;
	}
	
	/**
	*
	*   Computes the eccentric anomaly for elliptic orbits
	*
	* Input/Output:
	*
	*   M         Mean anomaly in [rad]
	*   e         Eccentricity of the orbit [0,1[
	*   <return>  Eccentric anomaly in [rad]
	*
	*/

	public static double getEccAnom (double M, double e){

		  // Constants
		 int maxit = 15;
		 double eps = 100.0*minimum;

		  // Variables
		  int    i=0;
		  double E, f;

		  // Starting value
		  M = M%consta.pi2;   
		  if (e<0.8) E=M; else E=consta.pi;

		  // Iteration
		  do {
		    f = E - e*Math.sin(E) - M;
		    E = E - f / ( 1.0 - e*Math.cos(E) );
		    ++i;
		    if (i==maxit) {
		      //System.out.println( " convergence problems in EccAnom");
		      break;
		    }
		  }
		  while (Math.abs(f) > eps);

		  return E;
	}
        
        /**
         * getActualMeanAnomaly()
         * 
         * Method to get actual mean anomaly
         * 
         * INPUT:
         * OUTPUT:
         * double m - mean anomaly [rad]
         */
        public double getActualMeanAnomaly(){
            return this.M_actual;
        }

	/*
	TEST METHOD
	*/
	//public static void main(String args[]){
	//}

        /**
         * fromTleToKepler()
         *
         * Method to read TLE data and fill Kepler variable.
         *
         * INPUT:
         * String line2 - 2nd line of TLE
         *
         * OUTPUT:
         * Kepler kepler - kepler filled with TLE data
         *
         */
        public Kepler fromTleToKepler(String line2){
            Kepler kepler = new Kepler();

            try{
                //ascending node
                String string = line2.substring(17,25);
                double value = Double.parseDouble(string);
                value = Math.toRadians(value);
                kepler.Omega = value;
                //inclination
                string = line2.substring(8,16);
                value = Double.parseDouble(string);
                value = Math.toRadians(value);
                kepler.incl = value;
                //argument of perigee
                string = line2.substring(34,42);
                value = Double.parseDouble(string);
                value = Math.toRadians(value);
                kepler.omega = value;
                //mean anomaly
                string = line2.substring(43,51);
                value = Double.parseDouble(string);
                value = Math.toRadians(value);
                kepler.M = value;
                //mean motion
                string = line2.substring(52,63);
                value = Double.parseDouble(string);
                kepler.a = Kepler.getSMAxis2(Constants.GM_Earth,Constants.pi2*value/86400);
                //eccentricity
                string = line2.substring(26,33);
                value = Double.parseDouble("0." + string);
                kepler.e = value;
            }
            catch(NumberFormatException nfe){
                //System.out.println("Method fromTleToKepler: " + nfe.getMessage());
            }
            //catch(Exception e){
            //    System.out.println("Method fromTleToKepler: " + e.getMessage());
            //}
            catch(StringIndexOutOfBoundsException e){
                //System.out.println("Method fromTleToKepler: " + e.getMessage());
            }

            return kepler;
        }

        public static void main (String[] args){
            //double a = 71583232.68;
            //double n = getMeanMotion(new Constants().GM_Earth,a);
            //System.out.println("n " + n*86400/(2*Math.PI));
            /*
            Vector vector_1 = new Vector(3);
            Vector vector_2 = new Vector(3);
            Vector vector_3 = new Vector(3);
            
            vector_1.v[0] = 2.5229729375129677E7;
            vector_1.v[1] = 5.107166666690566E7;
            vector_1.v[2] = 1.4632863695981173E7;
            vector_2.v[0] = 2.5020172534840558E7;
            vector_2.v[1] = 5.0826383550155096E7;
            vector_2.v[2] = 1.455563014741832E7;
            vector_3.v[0] = 4.751504556159597E7;
            vector_3.v[1] = 4060392.57911194;
            vector_3.v[2] = 4744947.378274405;

            double distance = Vector.getSize(Vector.subtractVectors(vector_2, vector_1));
            double distance_2 = Vector.getSize(Vector.subtractVectors(vector_3, vector_2));

            System.out.println("distance_1 " + Vector.getSize(vector_1)/1000 + " km");
            System.out.println("distance_2 " + Vector.getSize(vector_2)/1000 + " km");
            System.out.println("distance_3 " + Vector.getSize(vector_3)/1000 + " km");
            System.out.println("distance_12 " + distance/1000 + " km");
            System.out.println("distance_23 " + distance_2/1000 + " km");
            System.out.println("speed " + distance/1000/90 + " km/s");

            //
            double speed = distance/90;
            double gm = new Constants().GM_Earth;
            double a = (gm*Vector.getSize(vector_1)/(speed*speed))/
                        (gm/(speed*speed)*2 - Vector.getSize(vector_1));
            System.out.println("a " + a/1000 + " km");

            double v_esc = Math.sqrt(2*gm/Vector.getSize(vector_1));
            System.out.println("v_esc " + v_esc/1000 + " km");

            double a_2 = 55969815.69; //[m]
            System.out.println("n " + getMeanMotion(gm,a_2)*86400/(2*Math.PI) + " rev/day");

            //test of method
            //getTLE("2009 JL2_3",55562, "00000A",Time.getMjd(2009, 5, 13, 21, 39, 0),155.0919787,333.2148232,123.4509830,357.6689633,26849328.03,0.737075648);

            //getTLE("VFMO090316",55570, "00000A",Time.getMjd(2009, 3, 16, 22, 11, 47.58),23.12265599,124.7120960,233.3682200,190.9894865,55560966.10,0.209129125);
            //getTLE("Bolid 090403_wrong",55580, "00000A",Time.getMjd(2009, 04, 03, 1, 32, 23.911),48.83046129,127.6176216,96.06421442,0.415067505,13954327.34,0.539564567);
            //getTLE("Bolid 090403",55580, "00000A",Time.getMjd(2009, 04, 03, 1, 32, 23.911),48.82879227,129.8897726,94.82299338,0.668592227,14658330.09,0.561826570);
            //getTLE("VFMO100422_Circle",55595, "00000A",Time.getMjd(2010,04,22,19,55,57.38),23.41952211,108.8175303,18.56781979,15.64258,60000000,0.0);
            getTLE("VFMO100422_NewMethod",55600, "00000A",Time.getMjd(2010,04,22,19,55,57.38), 25.223095953592825,168.42436232579848,197.58293028638417,236.27199104428428,49449529.74534309,0.28833411125675823);

            //energy integral
            double sma = 180000 + Constants.R_Earth;
            double vel = Math.sqrt(Constants.GM_Earth/sma);
            System.out.println("v " + vel/1000 + " km/s");
            */

            String string_TLE = getTLE("Reentry 100709",50001, "00000A",Time.getMjd(2010, 7, 9,23, 1, 58.839999),
                    52.01795458065492,231.50458551391998,151.75801385677477,338.8740497880426,
                    11062825.555889494,0.608126999014875, "85848-3");
            string_TLE = getTLE("Reentry 120322",50001, "00000A",Time.getMjd(2012, 03, 22, 19, 28, 10.73),
                    48.26804764,48.26219992,11.78830567,68.08288951,
                    6502513.614,0.051991745, "85848-3");
            //System.out.println(string_TLE);
            /*
            //Cosmos 2251 before colliosion
            //geoc. position vector of parent body [m], this case is Cosmos 2251
            Vector vec_PosParentBody = new Vector(3);
            vec_PosParentBody.v[0] = -1468065.3559454195;
            vec_PosParentBody.v[1] = 1585916.6866229775;
            vec_PosParentBody.v[2] = 6812734.5420169365;
            //geoc. velocity vector of parent body [m], this case is Cosmos 2251
            Vector vec_VelParentBody = new Vector(3);
            vec_VelParentBody.v[0] = -6995.621549549612;
            vec_VelParentBody.v[1] = -2453.961487093762;
            vec_VelParentBody.v[2] = -934.0938610237473;

            Kepler kepler_Test = getElementsFromPosAndVelVec(Constants.GM_Earth,
                                vec_PosParentBody, vec_VelParentBody);

            System.out.println("a: "+kepler_Test.a/1000+" km");
            System.out.println("e: "+kepler_Test.e);
            System.out.println("i: "+Math.toDegrees(kepler_Test.incl)+" °");
            System.out.println("Omega: "+Math.toDegrees(kepler_Test.Omega)+" °");
            System.out.println("omega: "+Math.toDegrees(kepler_Test.omega)+" °");
            System.out.println("Ma: "+Math.toDegrees(kepler_Test.M)+" °\n");
            */
        }

        /**
         * getTLE()
         *
         * Method to generate TLE from given data.
         *
         * INPUT:
         * String noradString - norad no of given object
         * String intIdString - int ID of given object
         * double mjd - date of epoch (mjd - [days])
         * double incl  - inclination [deg]
         * double Omega - Right Ascension of Ascending Node [deg]
         * double omega - Argument of Perigee [deg]
         * double sma - semi major axis [m]
         * double ecc - eccentricity
         * String bDragString - coefficient B*Drag
         * int revNoEpoch - revolution no of epoch
         *
         * String tleString - TLE --> 2 lines
         */
        public static String getTLE(String name, int norad, String intIdString, double mjd, double incl,
                double Omega, double omega, double meanAnom, double sma, double ecc, String bDragString){
            String tleString = "1 nnnnnU dddddd   tttttttttttttt  .00000000  00000-0  00000-0 0  0000" + "\n" +
                               "2 nnnnn iii.iiii rrr.rrrr eeeeeee aaa.aaaa mmm.mmmm v.vvvvvvvv     01";

            String stringHelp = "";
            if(norad < 10) stringHelp = ("0000" + norad + "     ").substring(0, 5);
            else if((norad >= 10)&&(norad < 100)) stringHelp = ("000" + norad + "     ").substring(0, 5);
            else if((norad >= 100)&&(norad < 1000)) stringHelp = ("00" + norad + "     ").substring(0, 5);
            else if((norad >= 1000)&&(norad < 10000)) stringHelp = ("0" + norad + "     ").substring(0, 5);
            else stringHelp = ("" + norad + "     ").substring(0, 5);
            //System.out.println(stringHelp);
            String stringHelp_2 = (intIdString + "       ").substring(0, 7);
            //System.out.println(stringHelp_2);
            //replacing all positions
            //tleString.replaceAll("nnnnn", stringHelp);
            //tleString.replaceFirst("nnnnn", stringHelp);
            //tleString.s("nnnnn", stringHelp);
            //System.out.println(tleString.indexOf("nnnnn"));
            //System.out.println(tleString);

            Time time = new Time();
            time = time.getDateTime(mjd);
            double mjd_2 = time.getMjd(time.year, 1, 1, 0, 0, 0.0);
            double diffMjd = mjd - mjd_2 + 1;//+1 is correction between TLE form and mjd form
            String stringHelp_3 = "";
            if(time.year<2000) stringHelp_3 = (time.year - 1900) + "";
            else  {
                if(time.year == 2000){
                    stringHelp_3 = "00";
                }
                else if((time.year > 2000)&&(time.year < 2010)){
                    stringHelp_3 = "0" + (time.year - 2000);
                }
                else if(time.year >= 2010){
                    stringHelp_3 = "" + (time.year - 2000);
                }
            }
            
            //date
            String stringHelp_4 = "";
            if(diffMjd < 10) stringHelp_4 = stringHelp_3 + "00" + diffMjd + "          ";
            else if((diffMjd >= 10)&&(diffMjd < 100)) stringHelp_4 = stringHelp_3 + "0" + diffMjd + "          ";
            else if(diffMjd >= 100) stringHelp_4 = stringHelp_3 + "" + diffMjd + "          ";

            stringHelp_4 = stringHelp_4.substring(0,14);
            //System.out.println(stringHelp_4);

            String stringHelp_5 = "";
            if(incl < 10) stringHelp_5 = "00" + incl + "000000000000";
            else if((incl >= 10)&&(incl < 100)) stringHelp_5 = "0" + incl + "000000000000";
            else stringHelp_5 = "" + incl + "000000000000";
            stringHelp_5 = stringHelp_5.substring(0,8);
            //System.out.println(stringHelp_5);

            String stringHelp_6 = "";
            if( Omega < 10) stringHelp_6 = "00" + Omega + "000000000000";
            else if(( Omega >= 10)&&( Omega < 100)) stringHelp_6 = "0" + Omega + "000000000000";
            else stringHelp_6 = "" + Omega + "000000000000";
            stringHelp_6 = stringHelp_6.substring(0,8);
            //System.out.println(stringHelp_6);

            String stringHelp_7 = "";
            if(omega < 10) stringHelp_7 = "00" + omega + "000000000000";
            else if((omega >= 10)&&(omega < 100)) stringHelp_7 = "0" + omega + "000000000000";
            else stringHelp_7 = "" + omega + "000000000000";
            stringHelp_7 = stringHelp_7.substring(0,8);
            //System.out.println(stringHelp_7);

            String stringHelp_8 = "";
            if(meanAnom < 10) stringHelp_8 = "00" + meanAnom + "000000000000";
            else if((meanAnom >= 10)&&(meanAnom < 100)) stringHelp_8 = "0" + meanAnom + "000000000000";
            else stringHelp_8 = "" + meanAnom + "000000000000";
            stringHelp_8 = stringHelp_8.substring(0,8);
            //System.out.println(stringHelp_8);
            
            String stringHelp_9 = "";
            //stringHelp_9 = ecc + "0000000000";
            int eccHelp = (int)(ecc * 10e6);
            if((eccHelp < 10e6)&&(eccHelp >= 10e5)) stringHelp_9 = eccHelp + "";
            else if((eccHelp < 10e5)&&(eccHelp >= 10e4)) stringHelp_9 = "0" + eccHelp;
            else if((eccHelp < 10e4)&&(eccHelp >= 10e3)) stringHelp_9 = "00" + eccHelp;
            else if((eccHelp < 10e3)&&(eccHelp >= 10e2)) stringHelp_9 = "000" + eccHelp;
            else if((eccHelp < 10e2)&&(eccHelp >= 10e1)) stringHelp_9 = "0000" + eccHelp;
            else if((eccHelp < 10e1)&&(eccHelp >= 10e0)) stringHelp_9 = "00000" + eccHelp;
            else stringHelp_9 = "000000" + eccHelp;
            
            //stringHelp_9 = stringHelp_9.substring(2,9);
            //System.out.println("ecc " + (float)ecc);
            //System.out.println("stringHelp_9 " + stringHelp_9);

            String stringHelp_10 = "";
            double doubleHelp = getMeanMotion(Constants.GM_Earth,sma)*86400/(2*Math.PI);
            if(doubleHelp<10) stringHelp_10 = "0" + doubleHelp + "0000000000000";
            else stringHelp_10 = "" + doubleHelp + "0000000000000";
            stringHelp_10 = stringHelp_10.substring(0,11);
            //System.out.println(stringHelp_10);

            if(name.equals("")) tleString = "";
            else tleString = name + "\n";
            String stringName = tleString;
            int lengthStringName = tleString.length();
            tleString = //(tleString + "1 " + stringHelp + "U" + " " + stringHelp_2 + "  " + stringHelp_4 +
                        //"  .00000000  00000-0  " + bDragString.substring(0,7) + " 0  0000").substring(lengthStringName,69+lengthStringName) + "\n" +
                        tleString + "1 " + stringHelp + "U" + " " + stringHelp_2 + "  " + stringHelp_4 +
                        "  .00000000  00000-0  " + bDragString.substring(0,7) + " 0  0000" + "\n" +
                        "2 " + stringHelp + " " + stringHelp_5 + " " + stringHelp_6 +
                        " " + stringHelp_9 + " " + stringHelp_7 + " " + stringHelp_8 +
                        " " + stringHelp_10 + "    01";
            //System.out.println(tleString);
            return tleString;
        }
        /**
         * getTLE()
         *
         * Method to generate TLE from given data.
         *
         * INPUT:
         * String noradString - norad no of given object
         * String intIdString - int ID of given object
         * double mjd - date of epoch (mjd - [days])
         * double incl  - inclination [deg]
         * double Omega - Right Ascension of Ascending Node [deg]
         * double omega - Argument of Perigee [deg]
         * double sma - semi major axis [m]
         * double ecc - eccentricity
         * String bDragString - coefficient B*Drag
         * int revNoEpoch - revolution no of epoch
         *
         * String tleString - TLE --> 2 lines
         */
        public static String getTLE(String name, int norad, String intIdString, double mjd, double incl,
                double Omega, double omega, double meanAnom, double sma, double ecc, double bDrag){
            String tleString = "1 nnnnnU dddddd   tttttttttttttt  .00000000  00000-0  00000-0 0  0000" + "\n" +
                               "2 nnnnn iii.iiii rrr.rrrr eeeeeee aaa.aaaa mmm.mmmm v.vvvvvvvv     01";

            String stringHelp = "";
            if(norad < 10) stringHelp = ("0000" + norad + "     ").substring(0, 5);
            else if((norad >= 10)&&(norad < 100)) stringHelp = ("000" + norad + "     ").substring(0, 5);
            else if((norad >= 100)&&(norad < 1000)) stringHelp = ("00" + norad + "     ").substring(0, 5);
            else if((norad >= 1000)&&(norad < 10000)) stringHelp = ("0" + norad + "     ").substring(0, 5);
            else stringHelp = ("" + norad + "     ").substring(0, 5);
            //System.out.println(stringHelp);
            String stringHelp_2 = (intIdString + "       ").substring(0, 7);
            //System.out.println(stringHelp_2);
            //replacing all positions
            //tleString.replaceAll("nnnnn", stringHelp);
            //tleString.replaceFirst("nnnnn", stringHelp);
            //tleString.s("nnnnn", stringHelp);
            //System.out.println(tleString.indexOf("nnnnn"));
            //System.out.println(tleString);

            Time time = new Time();
            time = time.getDateTime(mjd);
            double mjd_2 = time.getMjd(time.year, 1, 1, 0, 0, 0.0);
            double diffMjd = mjd - mjd_2 + 1;//+1 is correction between TLE form and mjd form
            String stringHelp_3 = "";
            if(time.year<2000) stringHelp_3 = (time.year - 1900) + "";
            else  {
                if(time.year == 2000){
                    stringHelp_3 = "00";
                }
                else if((time.year > 2000)&&(time.year < 2010)){
                    stringHelp_3 = "0" + (time.year - 2000);
                }
                else if(time.year >= 2010){
                    stringHelp_3 = "" + (time.year - 2000);
                }
            }

            //date
            String stringHelp_4 = "";
            if(diffMjd < 10) stringHelp_4 = stringHelp_3 + "00" + diffMjd + "          ";
            else if((diffMjd >= 10)&&(diffMjd < 100)) stringHelp_4 = stringHelp_3 + "0" + diffMjd + "          ";
            else if(diffMjd >= 100) stringHelp_4 = stringHelp_3 + "" + diffMjd + "          ";

            stringHelp_4 = stringHelp_4.substring(0,14);
            //System.out.println(stringHelp_4);

            String stringHelp_5 = "";
            if(incl < 10) stringHelp_5 = "00" + incl + "000000000000";
            else if((incl >= 10)&&(incl < 100)) stringHelp_5 = "0" + incl + "000000000000";
            else stringHelp_5 = "" + incl + "000000000000";
            stringHelp_5 = stringHelp_5.substring(0,8);
            //System.out.println(stringHelp_5);

            String stringHelp_6 = "";
            if( Omega < 10) stringHelp_6 = "00" + Omega + "000000000000";
            else if(( Omega >= 10)&&( Omega < 100)) stringHelp_6 = "0" + Omega + "000000000000";
            else stringHelp_6 = "" + Omega + "000000000000";
            stringHelp_6 = stringHelp_6.substring(0,8);
            //System.out.println(stringHelp_6);

            String stringHelp_7 = "";
            if(omega < 10) stringHelp_7 = "00" + omega + "000000000000";
            else if((omega >= 10)&&(omega < 100)) stringHelp_7 = "0" + omega + "000000000000";
            else stringHelp_7 = "" + omega + "000000000000";
            stringHelp_7 = stringHelp_7.substring(0,8);
            //System.out.println(stringHelp_7);

            String stringHelp_8 = "";
            if(meanAnom < 10) stringHelp_8 = "00" + meanAnom + "000000000000";
            else if((meanAnom >= 10)&&(meanAnom < 100)) stringHelp_8 = "0" + meanAnom + "000000000000";
            else stringHelp_8 = "" + meanAnom + "000000000000";
            stringHelp_8 = stringHelp_8.substring(0,8);
            //System.out.println(stringHelp_8);

            String stringHelp_9 = "";
            //stringHelp_9 = ecc + "0000000000";
            int eccHelp = (int)(ecc * 10e6);
            if((eccHelp < 10e6)&&(eccHelp >= 10e5)) stringHelp_9 = eccHelp + "";
            else if((eccHelp < 10e5)&&(eccHelp >= 10e4)) stringHelp_9 = "0" + eccHelp;
            else if((eccHelp < 10e4)&&(eccHelp >= 10e3)) stringHelp_9 = "00" + eccHelp;
            else if((eccHelp < 10e3)&&(eccHelp >= 10e2)) stringHelp_9 = "000" + eccHelp;
            else if((eccHelp < 10e2)&&(eccHelp >= 10e1)) stringHelp_9 = "0000" + eccHelp;
            else if((eccHelp < 10e1)&&(eccHelp >= 10e0)) stringHelp_9 = "00000" + eccHelp;
            else stringHelp_9 = "000000" + eccHelp;

            //stringHelp_9 = stringHelp_9.substring(2,9);
            //System.out.println("ecc " + (float)ecc);
            //System.out.println("stringHelp_9 " + stringHelp_9);

            String stringHelp_10 = "";
            double doubleHelp = getMeanMotion(Constants.GM_Earth,sma)*86400/(2*Math.PI);
            if(doubleHelp<10) stringHelp_10 = "0" + doubleHelp + "0000000000000";
            else stringHelp_10 = "" + doubleHelp + "0000000000000";
            stringHelp_10 = stringHelp_10.substring(0,11);
            //System.out.println(stringHelp_10);

            String stringHelp_11 = "";
            int logBDrag = (int)Math.floor(Math.log10(bDrag));
            double expBDrag = Math.pow(10,logBDrag);
            int bdragInt = (int)Math.round(bDrag/expBDrag*10000);
            if(bDrag>=0.1) stringHelp_11 = bdragInt + "+" + (logBDrag+1);
            else stringHelp_11 = bdragInt + "" + (logBDrag+1);

            if(name.equals("")) tleString = "";
            else tleString = name + "\n";
            String stringName = tleString;
            int lengthStringName = tleString.length();
            tleString = //(tleString + "1 " + stringHelp + "U" + " " + stringHelp_2 + "  " + stringHelp_4 +
                        //"  .00000000  00000-0  " + bDragString.substring(0,7) + " 0  0000").substring(lengthStringName,69+lengthStringName) + "\n" +
                        tleString + "1 " + stringHelp + "U" + " " + stringHelp_2 + "  " + stringHelp_4 +
                        "  .00000000  00000-0  " + stringHelp_11 + " 0  0000" + "\n" +
                        "2 " + stringHelp + " " + stringHelp_5 + " " + stringHelp_6 +
                        " " + stringHelp_9 + " " + stringHelp_7 + " " + stringHelp_8 +
                        " " + stringHelp_10 + "    01";
            //System.out.println(tleString);
            return tleString;
        }


    /**
    *
    * getElementsFromPosAndVelVec():  computes the elements of an elliptical orbit from position
    *            and velocity vectors
    *
    * Comment: Method from APC_Kepler.cpp
    *
    * Input:
    *
    *   GM       Product of gravitational constant and centre mass [m^3/s^2]
    *   r        Geocentric elliptical position in [m]
    *   v        Geocentric elliptical velocity vector in [m/s]
    *
    * Output:
    *
    *   a        Semimajor axis of the orbit in [m]
    *   e        Eccentricity of the orbit
    *   i        Inclination of the orbit to the ecliptic in [rad]
    *   Omega    Longitude of the ascending node of the orbit in [rad]
    *   omega    Argument of perihelion in [rad]
    *   M        Mean anomaly in [rad]
    *
    */
    public static Kepler getElementsFromPosAndVelVec ( double GM, Vector r, Vector v){
          //
          // Variables
          //
          Vector  h = new Vector(3);
          double H, u, R, v2;
          double eCosE, eSinE, e2, E, nu;

          //orbital elements
          double Omega_In, i_In, u_In, a_In, e_In, M_In, omega_In;

          //kepler
          Kepler kepler = new Kepler();

          h = h.getVectorProduct(r, v);//Cross(r,v);          // Areal velocity
          H = h.getSize(h);//Norm(h);

          Omega_In = Math.atan2(h.v[0], -h.v[1]);                               // Long. ascend. node
          i_In     = Math.atan2(Math.sqrt(h.v[0]*h.v[0]+h.v[1]*h.v[1]), h.v[2]);// Inclination
          u     = Math.atan2(r.v[2]*H, -r.v[0]*h.v[1]+r.v[1]*h.v[0]);           // Arg. of latitude

          R  = r.getSize(r);//Norm(r);                                          // Distance
          v2 = v.getScalarProduct(v,v);//Dot(v, v);                               // Velocity squared

          a_In = 1.0 / (2.0/R-v2/GM);                                           // Semi-major axis

          eCosE = 1.0-R/a_In;                                                      // e*cos(E)
          eSinE = v.getScalarProduct(r,v)/Math.sqrt(GM*a_In);                        //Dot(r, v)/sqrt(GM*a);                      // e*sin(E)

          e2 = eCosE*eCosE +eSinE*eSinE;
          e_In  = Math.sqrt(e2);                                                // Eccentricity
          E  = Math.atan2(eSinE,eCosE);                                         // Eccentric anomaly

          M_In  = E - eSinE;                                                    // Mean anomaly

          nu = Math.atan2(Math.sqrt(1.0-e2)*eSinE, eCosE-e2);                   // True anomaly

          omega_In = u - nu;                                                    // Arg. of perihelion

          if (Omega_In<0.0) Omega_In += 2.0*Math.PI;
          if (omega_In<0.0) omega_In += 2.0*Math.PI;
          if (M_In   <0.0) M_In     += 2.0*Math.PI;

          kepler = new Kepler();
          kepler.omega=omega_In;
          kepler.Omega=Omega_In;
	  kepler.incl=i_In;
	  kepler.e=e_In;
          //from [AU] to [m]
          kepler.a=a_In;
	  kepler.M=M_In;

          return kepler;
    };

}