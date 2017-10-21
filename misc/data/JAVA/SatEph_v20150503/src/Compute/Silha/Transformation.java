//------------------------------------------------------------------------------
//
//Transformation.java
// 
//Ucel (purpose) :
//
//	Trieda sluzi na transformacie vektorov z jedneho ref. systemu do druheho. (Transformation of vectors from 1st coordinate system to other one)
//
//Poznamka (notes):
//
//	Trieda je modifikacia triedy SAT_RefSys.h od O. Montenbruck, E. Gill (2005/04/14  OMO  Final version (2nd reprint))
//
// 2007/10/11 - Jiri Silha
//
//------------------------------------------------------------------------------

package Compute.Silha;

//
// Deklaracia triedy (class declaration)
//

/** Transformation operations*/
public class Transformation{
	
	//
	//Konstanta pre iteracie
	//
	
	public static double minimum=1e-15;
	
	public static Constants c = new Constants();
	
	/**
	*
	* fromGeocToGeod
	*
	* Ucel (purpose):
	*
	*   Vypocet geodetickych suradnic z geocentrickych (greenwichskych, kde x smeru do nulteho poludnika).
	*   Vystupna hodnota bude premenna typu Geodetic, ktora bude obsahovat zem. dlzku, sirku a nadmorsku vysku
	*   Geodetic coordinates from geocentric coordinates 
	*
	* Vstup/Vystup (Input/output):
	*
	*   r - position vector [m]	  
	*   R_equ - equator radius 
	*   f - splostenie (flattening) 
	*
	*   <navratna hodnota> (<return>) Geodetic(lon, lat, altitude)
	*
	*
	*/
	public static Geodetic fromGeocToGeod(Vector r, double R_equ, double f){
		Geodetic geodetic=new Geodetic(0,0,0);
		
		double lon, lat, h;
		
		double  eps     = 1.0e3*minimum;    // Convergence criterion 
		double  epsRequ = eps*R_equ;
		double  e2      = f*(2.0-f);        // Square of eccentricity
  
		double  X = r.v[0];                   // Cartesian coordinates
		double  Y = r.v[1];
		double  Z = r.v[2];
		double  rho2 = X*X + Y*Y;           // Square of distance from z-axis
  
		// Check validity of input data
		if (r.getSize(r)==0.0) {
			lon=0.0; lat=0.0; h=-(new Constants()).R_Earth;
		}
		
		// Iteration 
		double  dZ, dZ_new, SinPhi;
		double  ZdZ, Nh, N;

		dZ = e2*Z;
		
		boolean goIteration=true;

                int j = 0;
		while(true) {
			ZdZ    =  Z + dZ;
			Nh     =  Math.sqrt ( rho2 + ZdZ*ZdZ ); 
			SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
			N      =  R_equ / Math.sqrt(1.0-e2*SinPhi*SinPhi); 
			dZ_new =  N*e2*SinPhi;
			if(Math.abs(dZ-dZ_new) < epsRequ ) break;
			dZ = dZ_new;
                        if(j > 1000) break;
                        j++;
		}
    
		// Longitude, latitude, altitude

		lon = Math.atan ( Y/X );
		//oprava o kvadranty
		//2.kvadrant
		if((X<0)&&(Y>0)) lon=lon+Math.PI;
		if((X<0)&&(Y<0)) lon=lon-Math.PI;

		lat = Math.atan ( ZdZ/Math.sqrt(rho2) );

		h   = Nh - N;
		
		geodetic.lon=lon;
		geodetic.lat=lat;
		geodetic.altitude=h;
		return geodetic;
	}
	
	/**
	*
	*  Vypocet equatorialnych suradnic (os x smeruje do jarneho bodu) z geocentrickych.
	*   Vystupna hodnota bude geocentricky vektor vzhladom na jarny bod a rovnik
	*   Equatorial coordinates from geocentric coordinates 
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector - geocentric position vector
	*  MJD	
	*
	*   <navratna hodnota> (<return>) Equatorial position vector for MJD
	*
	*/
	
	public static Vector fromGeocToEquat(Vector geocentric, double mjd){
		Vector equatorial=new Vector(3);
		
		Time time=new Time();
		Matrix u=new Matrix(3,3);
		//inicializacia translacnej matice z ekvat. => geocentrickych, pouzije sa jej transp. matica
		u=u.R_z(time.getGMST(mjd));
			
		//pozicny vektor pozorovatela v ekvatorilanych suradniciach v danom case pozorovanie
		equatorial=u.matrixMultiplyVector(u.transposeMatrix(u),geocentric);
		
		return equatorial;
	}
	
	/**
	*
	* Vypocet z equatorialnych suradnic (os x smeruje do jarneho bodu)  geocentricke suradnice.
	*   Vystupna hodnota bude geocentricky vektor vzhladom na nulty polidnik
	*   Geocentric coordinates from  equatorial coordinates 
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector - equatorial position vector
	*  MJD	
	*
	*   <navratna hodnota> (<return>) Geocentric position vector
	*
	*/
	
	public static Vector fromEquatToGeoc(Vector equatorial, double mjd){
		Vector geocentric = new Vector(3);
		
		Time time=new Time();
		Matrix u = new Matrix(3,3);
		//inicializacia translacnej matice z ekvat. => geocentrickych, pouzije sa jej transp. matica
		u = u.R_z(time.getGMST(mjd));
			
		//pozicny vektor pozorovatela v ekvatorilanych suradniciach v danom case pozorovanie
		geocentric = u.matrixMultiplyVector(u,equatorial);
		
		return geocentric;
	}
	
	/**
	*
	* fromGeodToGeoc
	*
	* Ucel (purpose):
	*
	*   Vypocet geocentrickych suradnic z geodetickych(greenwichskych, kde x smeru do nulteho poludnika).
	*   Vystupna hodnota bude geocentricky pozicny vektor
	*   Geocentric coordinates from geodetic coordinates 
	*
	* Vstup/Vystup (Input/output):
	*
	*   lon - zem.sirka [rad] 	  
	*   lat - zem.dlzka [rad] 
	*   altitude h - vyska [m]
	*   R_equ - earth equator  
	*   f - flattening 
	*
	*   <navratna hodnota> (<return>) Pozicny vektor
	*
	*/
	
	public static Vector fromGeodToGeoc(double lon, double lat, double h, double R_equ, double f){
		Constants c=new Constants();
		while(lon>c.pi) lon=lon-c.pi2;
		while(lon<-c.pi) lon=lon+c.pi2;
		while(lat<-c.pi/2) lat=lat+c.pi;
		while(lat>c.pi/2) lat=lat-c.pi;
					
		Vector position=new Vector(3);
		double  e2=f*(2.0-f);        // Square of eccentricity
	        double  CosLat = Math.cos(lat);  // (Co)sine of geodetic latitude
                double  SinLat = Math.sin(lat);

                double  N;
                
                // Position vector 
                N = R_equ / Math.sqrt(1.0-e2*SinLat*SinLat);
                
                position.v[0] =  (N+h)*CosLat*Math.cos(lon);
                position.v[1] =  (N+h)*CosLat*Math.sin(lon);
                position.v[2] =  ((1.0-e2)*N+h)*SinLat;

                /*
                //mega test - GEOTRANS coordinates of AGO Modra
                position.v[0] =  4053711.0;
                position.v[1] =  1260573.0;
                position.v[2] =  4744958.0;
                */
                return position;
	}
	
	/*
	TRANSLACNE MATICE (TRANSFORMATION MATRIXES)
	*/
	
	/**
	*
	*  Translacna matica z greenwichskych suradnic do lokalnych horizotalnych suradnic
	*  Transformation from Greenwich meridian system to local tangent coordinates
	*
	* Vstup/Vystup (Input/output):
	*
	*   lambda - zemepisna dlzka  (Geodetic East longitude [rad])
	*   phi - zemepisna sirka ( Geodetic latitude [rad])
	*
	*   <Navratova hodnota> <return>  Rotacna matica z geocentrickych suradnic do horizontalnych
	*	Rotation matrix from the Earth equator and Greenwich meridian
	*       to the local tangent (East-North-Zenith) coordinate system
	*
	*/
	
	public static Matrix getLtcMatrix (double lambda, double phi){
  
		Matrix  M = new Matrix(3,3);
		double  Aux;
  
		// Transformation to Zenith-East-North System
		M = M.getMatrixProduct(M.R_y(-phi),M.R_z(lambda));
  
		// Cyclic shift of rows 0,1,2 to 1,2,0 to obtain East-North-Zenith system
		for (int j=0; j<3; j++) {
			Aux=M.matrix[0][j]; M.matrix[0][j]=M.matrix[1][j]; M.matrix[1][j]=M.matrix[2][j]; M.matrix[2][j]= Aux;
		}
	  
		return  M;
	}
	
	/**
	*   Vypocet koeficientu D podla vztahu 2.127 
	*   Compute coefficient D (formula 2.127)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Unit vector e1
	*   Unit vector e2
	*   Unit vector e3
	*
	*    <Navratova hodnota> <return>  double D
	*
	*/
	
	public static double getCoefficient_D(Vector e1, Vector e2, Vector e3){
		double d=e1.getScalarProduct(e1,e1.getVectorProduct(e2,e3));
		return d;
	}
	
	/**
	*   Vypocet koeficientu Dij podla vztahu 2.128 
	*   Compute coefficient Dij (formula 2.128)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector di - auxiliary vector - doplnujuci vektor (formula 2.125)
	*   Vector Rj - position vector
	*
	*    <Navratova hodnota> <return>  double Dij
	*
	*/
	
	public static double getCoefficient_Dij(Vector di, Vector Ri){
		double dij=di.getScalarProduct(di,Ri);
		return dij;
	}
	
	/**
	*   Vypocet doplnujuceho vektora di podla vztahu 2.125	
	*   Compute auxiliary vector di (formula 2.125)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector ej - unit vector (method getUnitVectorE)
	*   Vector ek - unit vector (method getUnitVectorE)
	*
	*    <Navratova hodnota> <return>  double di
	*
	*/
	
	public static Vector getVector_di(Vector ej, Vector ek){
		Vector di=new Vector(3);
		di=di.getVectorProduct(ej,ek);
		return di;
	}
	
	/*
	*   Vypocet koeficientu n1 podla vztahu 2.132 
	*   Compute coefficient n1 (formula 2.132)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Eta_1
	*   Eta_2
	*   Mjd_1    Time t_1 (Modified Julian Date)
	*   Mjd_2    Time t_2 (Modified Julian Date) 
	*   Mjd_3    Time t_3 (Modified Julian Date) 
	*
	*    <Navratova hodnota> <return>  double n1
	*
	*/
	
	public static double getCoefficient_n1(double Eta_1, double Eta_2, double  Mjd_1, double  Mjd_2, double  Mjd_3){
		double n1=(Eta_2/Eta_1)*(Mjd_3-Mjd_2)/(Mjd_3-Mjd_1);
		return n1;
	}
	
	/**
	*   Vypocet koeficientu n3 podla vztahu 2.132 
	*   Compute coefficient n3 (formula 2.132)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Eta_2
	*   Eta_3
	*   Mjd_1    Time t_1 (Modified Julian Date)
	*   Mjd_2    Time t_2 (Modified Julian Date) 
	*   Mjd_3    Time t_3 (Modified Julian Date) 
	*
	*    <Navratova hodnota> <return>  double n3
	*
	*/
	
	public static double getCoefficient_n3(double Eta_2, double Eta_3, double  Mjd_1, double  Mjd_2, double  Mjd_3){
		double n3=(Eta_2/Eta_3)*(Mjd_2-Mjd_1)/(Mjd_3-Mjd_1);
		return n3;
	}
	
	/**
	*
	*   Vypocet koeficientu Ro1 podla vztahu 2.129
	*   Compute coefficient Ro1 (formula 2.129)
	*
	* Vstup/Vystup (Input/output):
	*
	*   coefficient n1
	*   coefficient n3
	*   D - formula (2.127)
	*   D11 - formula (2.128)
	*   D12 - formula (2.128)
	*   D13 - formula (2.128)
	*
	*    <Navratova hodnota> <return>  double Ro1 - station - satellite distance (1st position)
	*
	*/
	
	public static double getCoefficient_Ro1(double n1, double n3, double d, double d11, double d12, double d13){
		double ro1=-1/(n1*d)*(n1*d11-d12+n3*d13);
		return ro1;
	}
	
	/*
	*
	*   Vypocet koeficientu Ro2 podla vztahu 2.129
	*   Compute coefficient Ro2 (formula 2.129)
	*
	* Vstup/Vystup (Input/output):
	*
	*    coefficient n1
	*    coefficient n3
	*    D - formula (2.127)
	*    D11 - formula (2.128)
	*    D12 - formula (2.128)
	*    D13 - formula (2.128)
	*
	*    <Navratova hodnota> <return>  double Ro2 - station - satellite distance (2nd position)
	*
	*/
	
	public static double getCoefficient_Ro2(double n1, double n3, double d, double d21, double d22, double d23){
		double ro2=1/d*(n1*d21-d22+n3*d23);
		return ro2;
	}
	
	/*
	*
	*   Vypocet koeficientu Ro3 podla vztahu 2.129
	*   Compute coefficient Ro3 (formula 2.129)
	*
	* Vstup/Vystup (Input/output):
	*
	*    coefficient n1
	*    coefficient n3
	*    D - formula (2.127)
	*    D11 - formula (2.128)
	*    D12 - formula (2.128)
	*    D13 - formula (2.128)
	*
	*    <Navratova hodnota> <return>  double Ro3 - station - satellite distance (3rd position)
	*
	*/
	
	public static double getCoefficient_Ro3(double n1, double n3, double d, double d31, double d32, double d33){
		double ro3=-1/(n3*d)*(n1*d31-d32+n3*d33);
		return ro3;
	}
	
	/**
	*   Vypocet geoc. pozicneho vektora satelitu (vztah 2.122)
	*   Compute sat.geoc.position vector (formula 2.122)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Vector R - geoc. position vector of ground station
	*   Vector e - unit vector describes the direction of observation
	*   ro - distance sat - station 
	*
	*    <Navratova hodnota> <return>  r - geoc.sat.pos. vector
	*
	*/
	
	public static Vector getSatPosition(Vector R, Vector e, double ro){
		Vector r = new Vector(3);
		r = r.addVectors(R,e.multiplyVector(e,ro));
		return r;
	}
	
	/**
	*
	*   Vypocet horizontalnych suradnic (azimut, vyska, vzdialenost telesa) z pozicneho vektora (os x smeruje do jarneho bodu) 
	*   Computes horizontal coordinates (azimuth, altitude, range) from position vector
	*	
	* Vstup/Vystup (Input/output):
	*
	*   Vector r - position vector (x,y,z) [m,m,m]
	*   Geodetic R -observatory position vector (lon,lat,altitude) [deg,deg,m]
	*   double time - observatotion time [s], MJD
	* 
	*    <Navratova hodnota> <return>  r - (Az, h, range) - horizonthal  coordinates, [rad,rad,m]
	*
        */
	
        public static Vector getHorizontalCoordinates(Vector r, Geodetic R, double time){
		Matrix ltc,mjd=new Matrix(3,3);
		Vector horizontal=new Vector(3);
		//from equatorial => geocentric
		horizontal=fromEquatToGeoc(r,time);
		//from geocentic => topocentric
		ltc=getLtcMatrix(R.lon,R.lat);
		//from center to observatory
		Vector geodeticR=fromGeodToGeoc(R.lon,R.lat,R.altitude, c.R_Earth, c.f_Earth);
		horizontal=horizontal.subtractVectors(horizontal,geodeticR);
		
		horizontal=ltc.matrixMultiplyVector(ltc,horizontal);
                
                return horizontal;
	}
	
	//oprava Az a h
	public static Vector getHorizontalCoordinates2(Vector r){
		Vector vector=new Vector(3);
		//azimut
		vector.v[0]=Math.atan2(r.v[1],r.v[0]);
		vector.v[0]=-vector.v[0]+Math.PI/2;
		while((vector.v[0]<0)) vector.v[0]=vector.v[0]+Math.PI*2;
		while((vector.v[0]>Math.PI*2)) vector.v[0]=vector.v[0]-Math.PI*2;
		//vyska
		vector.v[1]=Math.acos(r.v[2]/(Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1]+r.v[2]*r.v[2])));
		while((vector.v[1]<0)) vector.v[1]=vector.v[1]+Math.PI;
		while((vector.v[1]>Math.PI)) vector.v[1]=vector.v[1]-Math.PI;
		
		//oprava na vysku (-90,90)
		vector.v[1]=Math.PI/2-vector.v[1];
		//vzdialenost
		vector.v[2]=Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1]+r.v[2]*r.v[2]);
		return vector;
	
	}
	
        /**
	*  Method is to get actual horizontal coordinates
	*
	*IN:   Vector r - position vector (x,y,z) [m,m,m] [equatorial coor. - x => vernal equinox]
	*	Geodetic R -observatory position vector (lon,lat,altitude) [rad,rad,m]
	*	double time - observatotion time [days], MJD
	*
	*OUT:  OrbitalState equatorial coordinates [m,m,m],[m/s,m/s,m/s]
	*
	*/
	
	public static Vector getHorizontalCoordinates3(Vector r, Geodetic R, double time){
		Matrix ltc,mjd=new Matrix(3,3);
		Vector horizontal=new Vector(3);
		//from equatorial => geocentric
		horizontal=fromEquatToGeoc(r,time);
		//from geocentic => topocentric
		ltc=getLtcMatrix(R.lon,R.lat);
		//from center to observatory
		Vector geodeticR=fromGeodToGeoc(R.lon,R.lat,R.altitude, c.R_Earth, c.f_Earth);
		horizontal=horizontal.subtractVectors(horizontal,geodeticR);
		
		horizontal=ltc.matrixMultiplyVector(ltc,horizontal);
		
		horizontal = horizontal.fromCartesianToPolar(horizontal);
		//angle horizontal.v[o] is from geoc. x axis to geoc. y axis, it is
		//supposed to get from local meridian (NORTH )to EAST ==> in opossite direction, and from y
		horizontal.v[0] = Math.PI/2 - horizontal.v[0];
		while(horizontal.v[0] < 0) horizontal.v[0] = horizontal.v[0] +2*Math.PI;
		return horizontal;
	
	}
        
	/**
	*   Vypocet sferickych ekvatorialnych suradnic (R.A,, dec, vzdialenost telesa) z pozicneho vektora (os x smeruje do jarneho bodu) 
	*   Computes spherical equat. coordinates (R.A., dec, range) from position vector
	*	
	* Vstup/Vystup (Input/output):
	*
	*   Vector r - position vector (x,y,z) [m,m,m]
	* 
	*    <Navratova hodnota> <return>  r - (dec, R.A., range) - spher. equat.  coordinates, [rad,rad,m]
	*
	*/
	
	public static Vector getEquatSphericalCoordinates(Vector r){
		Vector spher=getEquatSphericalCoordinates2(r);
		//oprava dec
		while((spher.v[1]<0)) spher.v[1]=spher.v[1]+Math.PI*2;
		while((spher.v[1]>Math.PI*2)) spher.v[1]=spher.v[1]-Math.PI*2;
		//oprava na deklinaciu (-90,90)
		spher.v[1]=Math.PI/2-spher.v[1];
		//vzdialenost
		spher.v[2]=Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1]+r.v[2]*r.v[2]);
		
		return spher;
	}
	
	public static Vector getEquatSphericalCoordinates2(Vector r){
		Vector vector=new Vector(3);
		//R.A.
		vector.v[0]=Math.atan2(r.v[1],r.v[0]);
		while((vector.v[0]<0)) vector.v[0]=vector.v[0]+Math.PI*2;
		while((vector.v[0]>Math.PI*2)) vector.v[0]=vector.v[0]-Math.PI*2;
		
		//deklinacia
		vector.v[1]=Math.atan2(r.v[2],(Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1])));
		
		//vzdialenost
		vector.v[2]=Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1]+r.v[2]*r.v[2]);
		return vector;
	}
        
        /**
         *fromEquatToHorizontal
         *
         *  IN: 
         *  R.A. [rad]
         *  dec  [rad]
         *  GMST [rad]     
         *  
         *  lon  [rad]  - obs. pos.
         *  lat  [rad]  - obs. pos.
         *
         *  OUT:
         *  Vector:
         *  v[0] = Az [rad]
         *  v[1] = h  [rad]
         */
        public static Vector fromEquatToHorizontal(double RA, double dec, double gmst,
                        double lon, double lat){
            Vector vector = new Vector(2);
            //hour angle
            double t = gmst +lon - RA;

            vector.v[0] = Math.atan2(Math.sin(t), Math.cos(t)*Math.sin(lat) - Math.tan(dec)*Math.cos(lat)) + c.pi;//Math.PI;
                       
            vector.v[1] = Math.asin(Math.sin(lat)*Math.sin(dec) + Math.cos(lat)*Math.cos(dec)*Math.cos(t));
            return vector;
        }

        /**
         *fromHorizontalToEquat
         *
         *  IN:
         *  Az [rad]
         *  h  [rad]
         *  GMST [rad]
         *
         *  lon  [rad]  - obs. pos.
         *  lat  [rad]  - obs. pos.
         *
         *  OUT:
         *  Vector:
         *  v[0] = R.A. [rad]
         *  v[1] = dec  [rad]
         */
        public static Vector fromHorizontalToEquat(double Az, double h, double gmst,
                        double lon, double lat){
            Vector vector = new Vector(2);
            //hour angle
            //double t;// = gmst +lon - RA;

            //declination
            vector.v[1] = Math.asin(Math.sin(lat) * Math.sin(h) + Math.cos(lat) * Math.cos(h) * Math.cos(Az));
            //hour angle
            double cos_H = (Math.cos(lat) * Math.sin(h) - Math.sin(lat) * Math.cos(h) * Math.cos(Az))/Math.cos(vector.v[1]);
            double sin_H = (-1) * (Math.sin(Az) * Math.cos(h)) / Math.cos(vector.v[1]);
            double H = Math.atan2(sin_H, cos_H);
            vector.v[0] = gmst + lon - H;
                    
            return vector;
        }
        
        /**
         * fromHoursToDegrees
         *
         *  IN: 
         *  angle [hours]
         *
         *  OUT:
         *  angle [rad]
         */
        public static double fromHourToDegrees(double hours){
            double deg = Math.toRadians(hours*15);
            return deg;
        }
        
         /**
         * fromHoursToDegrees
         *
         *  IN: 
         *  hours
         *  minutes 
         *  sec
         *
         *  OUT:
         *  angle [rad]
         */
        public static double fromHourToDegrees(int hours, int min, double sec){
            double deg = Math.toRadians((hours + (double)min/60 + sec/3600)*15);
            return deg;
        }
        
        /**
         * getAngle
         *
         *  IN: 
         *  deg    [deg]
         *  arcmin [arcmin]
         *  arcsec [arcsec]
         *
         *  OUT:
         *  angle [rad]
         */
        public static double getAngle(int deg, int min, double sec){
            double retDeg = Math.toRadians(deg + (double)min/60 + sec/3600);
            return retDeg;
        }
        
	//------------------------------------------------------------------------------
	//
	// getSpherEquat
	//
	// Ucel (purpose):
	//
	//   Vypocet sferickych equatorialnych suradnic  z kartezskych
	//   Computes spherical equatorial coordinates from cartesian
	//	
	// Vstup/Vystup (Input/output):
	//
	//   Vector R - equatorial coordinates (x,y,z) [m,m,m]
	// 
	//    <Navratova hodnota> <return>  r - (Dec, RA, r) - equatorial spherical coordinates, [rad,rad,m]
	//
	//------------------------------------------------------------------------------
	/*
	public static Vector getSpherEquat(Vector R){
		Vector r=new Vector(3);
		//declination
		r.v[0]=Math.atan2(R.v[1],R.v[0]);
		//if(R.v[0]<0) {
		//	if(r.v[0]<-c.pi/2) r.v[0]+=c.pi;
		//	if(r.v[0]>c.pi/2) r.v[0]-=c.pi;
		//}
		//r.v[0]=Math.atan(R.v[1]/R.v[0]);
		//right ascension
		r.v[1]=Math.atan2(R.v[2],(Math.sqrt(R.v[0]*R.v[0]+R.v[1]*R.v[1])));
		if((r.v[1]>c.pi/2)&&(r.v[1]<c.pi)) r.v[1]=Math.abs(r.v[1]-c.pi);
		else if((r.v[1]>c.pi*(double)3/2)&&(r.v[1]<c.pi2)) r.v[1]=-Math.abs(r.v[1]-c.pi2);
		//r.v[1]=Math.atan(R.v[2]/(Math.sqrt(R.v[0]*R.v[0]+R.v[1]*R.v[1])));
		//distance from center
		r.v[2]=Math.sqrt(R.v[0]*R.v[0]+R.v[1]*R.v[1]+R.v[2]*R.v[2]);
		
		return r;
	}
	
	//------------------------------------------------------------------------------
	//
	// getCartesFromSpherEquat
	//
	// Ucel (purpose):
	//
	//   Vypocet kartezskych equatorialnych suradnic z sferickych
	//   Computes cartesian equatorial coordinates from spherical
	//	
	// Vstup/Vystup (Input/output):
	//
	//   Vector r - equatorial coordinates (Dec,RA,r) [rad,rad,m]
	// 
	//    <Navratova hodnota> <return>  R - (x, y, z) - equatorial cartesian coordinates, [m, m, m]
	//
	//------------------------------------------------------------------------------
	
	public static Vector getCartesFromSpherEquat(Vector r){
		Vector R=new Vector(3);
		//x
		R.v[0]=r.v[2]*Math.cos(r.v[0])*Math.cos(r.v[1]);
		//y
		R.v[1]=r.v[2]*Math.cos(r.v[0])*Math.sin(r.v[1]);
		//z
		R.v[2]=r.v[2]*Math.sin(r.v[0]);
		
		return R;
	}
			
	/*
	TEST METHOD
	*/
	
	public static void main(String []args){
            /*
            Time time = new Time(2009,06,17,21,6,1);
            double timeMjd = time.getMjd(time);
            double gmst = time.getGMST(timeMjd);
            
            double lon = Math.toRadians(17.2740);
            double lat = Math.toRadians(48.3733);
            //double RA = fromHourToDegrees(11,56,24.7);
            //double dec = getAngle(-4, 19, 20.1);
            double RA = Math.toRadians(271.04400211);
            double dec = Math.toRadians(40.90840240);
            System.out.println("ra in     " + Math.toDegrees(RA));
            System.out.println("dec in    " + Math.toDegrees(dec));

            Vector vector = epoch_of_date_to_j2000(timeMjd, RA, dec);
            //System.out.println("ra_0 old  " + Math.toDegrees(vector.v[0]));
            //System.out.println("dec_0 old " + Math.toDegrees(vector.v[1]));

            Vector vector2 = epoch_of_date_to_j2000_2(timeMjd, RA, dec);
            //System.out.println("ra_0 old  " + Math.toDegrees(vector2.v[0]));
            //System.out.println("dec_0 old " + Math.toDegrees(vector2.v[1]));
            //System.out.println("d_ra      " + Math.toDegrees(vector.v[0] - vector2.v[0])*3600 + " arcsec");
            //System.out.println("d_dec     " + Math.toDegrees(vector.v[1] - vector2.v[1])*3600 + " arcsec");
            Vector vector4 = epoch_of_date_from_j2000(timeMjd, vector2.v[0], vector2.v[1]);
            //System.out.println("ra 2x red " + Math.toDegrees(vector4.v[0]));
            //System.out.println("dec 2x    " + Math.toDegrees(vector4.v[1]));

            Vector vector5 = epoch_of_date_from_j2000(timeMjd, RA, dec);
            System.out.println("ra corr  " + Math.toDegrees(vector5.v[0]));
            System.out.println("dec corr " + Math.toDegrees(vector5.v[1]));
            
            Vector vector3 = fromEquatToHorizontal(vector5.v[0], vector5.v[1], gmst, lon, lat);
            vector3.v[0] = Math.toRadians(97.119);
            vector3.v[1] = Math.toRadians(66.839);
            //Vector vector2 = fromEquatToHorizontal(vector.v[0], vector.v[1], gmst, lon, lat);
            System.out.println("Az   " + Math.toDegrees(vector3.v[0]));
            System.out.println("h    " + Math.toDegrees(vector3.v[1]));
            //System.out.println("d_Az " + (Math.toDegrees(vector3.v[0]) - 97.119)*3600);
            //System.out.println("d_h  " + (Math.toDegrees(vector3.v[1]) - 66.839)*3600);

            //back to eequatorial
            Vector vector6 = fromHorizontalToEquat(vector3.v[0],vector3.v[1],gmst,lon,lat);
            System.out.println("ra corr  " + Math.toDegrees(vector6.v[0]));
            System.out.println("dec corr " + Math.toDegrees(vector6.v[1]));

            Vector vector7 = epoch_of_date_to_j2000(timeMjd, vector6.v[0], vector6.v[1]);
            System.out.println("ra out   " + Math.toDegrees(vector7.v[0]));
            System.out.println("dec out  " + Math.toDegrees(vector7.v[1]));

            System.out.println("d_ra in out  " + Math.toDegrees(vector7.v[0] - RA) * 3600 + " arcsec");
            System.out.println("d_dec in out " + Math.toDegrees(vector7.v[1] - dec) * 3600 + " arcsec");
            */

            
            //test new methods - trosromation between hor to eq, and back 19.10.2009
            /*
            //test 1. time - UTC
            //Time time_2 = new Time(2009,10,19,17,51,34);
            //test 2. time
            Time time_2 = new Time(2009,10,19,21,0,32);
            double timeMjd_2 = time_2.getMjd(time_2);
            //AGO Modra position
            double lon_2 = Math.toRadians(17.2740);
            double lat_2 = Math.toRadians(48.3733);

            //incoming values
            double ra_In, dec_In, az_In, h_In;
            //outcoming values
            double ra_Out, dec_Out, az_Out, h_Out;
            //diffrence between incoming and outcoming values
            double d_ra, d_dec, d_az, d_h;

            //test 1.
            //ra_In = Math.toRadians(305.0753102);
            //dec_In = Math.toRadians(59.305433);
            //az_In = Math.toRadians(339.086);
            //h_In = Math.toRadians(77.998);

            //test 2.
            ra_In = Math.toRadians(46.20523884);
            dec_In = Math.toRadians(32.3076368);
            az_In = Math.toRadians(98.247);
            h_In = Math.toRadians(52.501);

            //from equatorial to horizontal
            az_Out = fromEquatToHorizCoord_New(ra_In, dec_In, lon_2, lat_2, timeMjd_2).v[0];
            h_Out = fromEquatToHorizCoord_New(ra_In, dec_In, lon_2, lat_2, timeMjd_2).v[1];

            while(az_Out > Math.PI*2) az_Out = az_Out - Math.PI*2;
            while(az_Out < -Math.PI*2) az_Out = az_Out + Math.PI*2;

            d_az = az_In - az_Out;
            //while(d_az > Math.PI*2) d_az = d_az - Math.PI*2;
            //while(d_az < -Math.PI*2) d_az = d_az + Math.PI*2;

            d_h = h_In - h_Out;
            //while(d_h > Math.PI) d_h = d_h - Math.PI;
            //while(d_h < -Math.PI) d_h = d_h + Math.PI;

            System.out.println("In  Az " + Math.toDegrees(az_In));
            System.out.println("Out Az " + Math.toDegrees(az_Out));
            System.out.println("d   Az " + Math.toDegrees(d_az));
            System.out.println("d   Az " + Math.toDegrees(d_az)*3600 + " arsec");
            System.out.println("In  h  " + Math.toDegrees(h_In));
            System.out.println("Out h  " + Math.toDegrees(h_Out));
            System.out.println("d   h  " + Math.toDegrees(d_h));
            System.out.println("d   h  " + Math.toDegrees(d_h)*3600 + " arsec");

            //from horizontal to equatorial
            ra_Out = fromHorizToEquatCoord_New(az_In, h_In, lon_2, lat_2, timeMjd_2).v[0];
            dec_Out = fromHorizToEquatCoord_New(az_In, h_In, lon_2, lat_2, timeMjd_2).v[1];
            while(ra_Out > Math.PI*2) ra_Out = ra_Out - Math.PI*2;
            while(ra_Out < -Math.PI*2) ra_Out = ra_Out + Math.PI*2;

            d_ra = ra_In - ra_Out;
            //while(d_ra > Math.PI*2) d_ra = d_ra - Math.PI*2;
            //while(d_ra < -Math.PI*2) d_ra = d_ra + Math.PI*2;

            d_dec = dec_In - dec_Out;
            //while(d_h > Math.PI) d_h = d_h - Math.PI;
            //while(d_h < -Math.PI) d_h = d_h + Math.PI;

            System.out.println("In   RA " + Math.toDegrees(ra_In));
            System.out.println("Out  RA " + Math.toDegrees(ra_Out));
            System.out.println("d    RA " + Math.toDegrees(d_ra));
            System.out.println("d    RA " + Math.toDegrees(d_ra)*3600 + " arsec");
            System.out.println("In  dec " + Math.toDegrees(dec_In));
            System.out.println("Out dec " + Math.toDegrees(dec_Out));
            System.out.println("d   dec " + Math.toDegrees(d_dec));
            System.out.println("d   dec " + Math.toDegrees(d_dec)*3600 + " arsec");

            //test of accurancy
            for(int i = 0; i < 10; i++){
                System.out.println(" ----------------- " + i + " ------------------ ");
                //from equatorial to horizontal
                az_Out = fromEquatToHorizCoord_New(ra_Out, dec_Out, lon_2, lat_2, timeMjd_2).v[0];
                h_Out = fromEquatToHorizCoord_New(ra_Out, dec_Out, lon_2, lat_2, timeMjd_2).v[1];

                while(az_Out > Math.PI*2) az_Out = az_Out - Math.PI*2;
                while(az_Out < -Math.PI*2) az_Out = az_Out + Math.PI*2;

                d_az = az_In - az_Out;
                //while(d_az > Math.PI*2) d_az = d_az - Math.PI*2;
                //while(d_az < -Math.PI*2) d_az = d_az + Math.PI*2;
                d_h = h_In - h_Out;

                System.out.println("In  Az " + Math.toDegrees(az_In));
                System.out.println("Out Az " + Math.toDegrees(az_Out));
                System.out.println("d   Az " + Math.toDegrees(d_az));
                System.out.println("d   Az " + Math.toDegrees(d_az)*3600 + " arsec");
                System.out.println("In  h  " + Math.toDegrees(h_In));
                System.out.println("Out h  " + Math.toDegrees(h_Out));
                System.out.println("d   h  " + Math.toDegrees(d_h));
                System.out.println("d   h  " + Math.toDegrees(d_h)*3600 + " arsec");

                //from horizontal to equatorial
                ra_Out = fromHorizToEquatCoord_New(az_Out, h_Out, lon_2, lat_2, timeMjd_2).v[0];
                dec_Out = fromHorizToEquatCoord_New(az_Out, h_Out, lon_2, lat_2, timeMjd_2).v[1];

                while(ra_Out > Math.PI*2) ra_Out = ra_Out - Math.PI*2;
                while(ra_Out < -Math.PI*2) ra_Out = ra_Out + Math.PI*2;

                d_ra = ra_In - ra_Out;
                //while(d_ra > Math.PI*2) d_ra = d_ra - Math.PI*2;
                //while(d_ra < -Math.PI*2) d_ra = d_ra + Math.PI*2;
                d_dec = dec_In - dec_Out;

                System.out.println("In   RA " + Math.toDegrees(ra_In));
                System.out.println("Out  RA " + Math.toDegrees(ra_Out));
                System.out.println("d    RA " + Math.toDegrees(d_ra));
                System.out.println("d    RA " + Math.toDegrees(d_ra)*3600 + " arsec");
                System.out.println("In  dec " + Math.toDegrees(dec_In));
                System.out.println("Out dec " + Math.toDegrees(dec_Out));
                System.out.println("d   dec " + Math.toDegrees(d_dec));
                System.out.println("d   dec " + Math.toDegrees(d_dec)*3600 + " arsec");
            }
            */
            /*
            //test of method getMeanObliquityOfEcliptic
            Time time = new Time(2009,10,20,20,0,0);
            double mjd = time.getMjd(time);
            double meanObl = getMeanObliquityOfEcliptic(mjd);
            System.out.println("meanObl " + Math.toDegrees(meanObl));
            */
            //test of method fromGeocentricToHeliocentric\
            //time - source http://www.erh.noaa.gov/box/equinox.html
            //vernal equinox
            //Time time = new Time(2009,3,20,11,44,0);
            //summer solstice
            //Time time = new Time(2009,6,21,5,45,0);
            //autumn equinox
            //Time time = new Time(2009,9,22,21,18,0);
            //winter solstice
            /*
            Time time = new Time(2009,12,21,17,47,0);

            double mjd = time.getMjd(time);
            //center of the Earth in geocentric equatorial coordinates
            Vector vectorPosGeocEq = new Vector(3);
            vectorPosGeocEq.v[0] = 0;
            vectorPosGeocEq.v[1] = 0;
            vectorPosGeocEq.v[2] = 0;
            //heliocentric ecliptical position of Earth for given time
            Vector vectorPosHeliocEcl = fromGeocentricToHeliocentric(vectorPosGeocEq, mjd);

            System.out.println("Xhel " + vectorPosHeliocEcl.v[0]/1000 + " km");
            System.out.println("Yhel " + vectorPosHeliocEcl.v[1]/1000 + " km");
            System.out.println("Zhel " + vectorPosHeliocEcl.v[2]/1000 + " km");
            System.out.println("Size " + vectorPosHeliocEcl.getSize(vectorPosHeliocEcl)/1000 + " km" );
            */

            //test
            Time time = new Time(2011,05,15,17,05,17.2);
            double mjd = Time.getMjd(time);
            //mjd = 55697.072006225586;
            //observatory
            Geodetic observatory = new Geodetic();
            observatory.lon = Math.toRadians(17.274);
            observatory.lat = Math.toRadians(48.373299);
            observatory.altitude = 531.0;
            //sun properties
            Sun sp = new Sun();
            sp = sp.getSunLocalPositions(mjd, observatory);
            //Equatorial coordinates
            Vector equatVec = new Vector(3);
            //RA
            equatVec.v[0] = sp.ra;
            //DEC
            equatVec.v[1] = sp.dec;
            //ecliptic coordinates            
            Vector eclipVec = new Vector(3);
            eclipVec = new Transformation().fromEquatorialToEcliptic2(equatVec,mjd);
            System.out.println("Time:    " + mjd + ", " + Time.getDateTime(mjd).year +  "," + Time.getDateTime(mjd).month +  "," + Time.getDateTime(mjd).day +
                     "," + Time.getDateTime(mjd).hour +  "," + Time.getDateTime(mjd).min +  "," + Time.getDateTime(mjd).sec);
            System.out.println("Equatorial Bef:    " + mjd + ", " + Math.toDegrees(equatVec.v[0]) + ", " + Math.toDegrees(equatVec.v[1]));
            System.out.println("Eliptic:           " + mjd + ", " + Math.toDegrees(eclipVec.v[0]) + ", " + Math.toDegrees(eclipVec.v[1]));
            //correction 60 degrees from Sun
            eclipVec.v[0] = eclipVec.v[0] + Math.PI/2;
            //test
            Vector equatVec2 = new Vector(3);
            equatVec2 = new Transformation().fromEclipticToEquatorial2(eclipVec, mjd);
            System.out.println("Equatorial Aft:    " + mjd + ", " + Math.toDegrees(equatVec2.v[0]) + ", " + Math.toDegrees(equatVec2.v[1]));
            Vector horVec = new Transformation().fromEquatToHorizCoord_New(equatVec2.v[0], equatVec2.v[1], observatory.lon, observatory.lat, mjd);
            System.out.println("Horizontal Aft:    " + mjd + ", " + Math.toDegrees(horVec.v[0]) + ", " + Math.toDegrees(horVec.v[1]));
            
            //Galactic coordinates tests
            Vector eq_vec = new Vector(2);
            Vector gal_vec = new Vector(2);
            
            eq_vec.v[0] = Math.toRadians(180);
            eq_vec.v[1] = Math.toRadians(45);
            
            gal_vec = Transformation.getGalacticCoordinates(eq_vec);
            System.out.println("LON " + Math.toDegrees(gal_vec.v[0]) + " LAT " + Math.toDegrees(gal_vec.v[1]));
                       
        }
         
        
        /**
         * epoch_of_date_to_j2000()
         * 
         * Method to convert equatorial coordinates from epoch? date to J2000?.
         * Modified from code of Mr. Neoklis Kyriazis (sat_code.zip), 
         * file observe.cpp, lines 78 - 88.
         * 
         * IN:
         *  double mjd - Modified Julian date [days]
         *  double ra - R.A. in epoch date??
         *  doubel dec - declination in epoch date
         * 
         *  OUT:
         *  Vector coordinates[2] - v.[0] - R.A. in J2000
         *                        - v.[1] - declination in J2000
         */
        
        public static Vector epoch_of_date_to_j2000( double mjd, double ra, double dec){
           Vector coordinates = new Vector(2);
           double jd = mjd + 2400000.5;
           double t_centuries = (jd - 2451545.) / 36525.;
           double m = (3.07496 + .00186 * t_centuries / 2.) * (new Constants().pi / 180.) / 240.;
           double n = (1.33621 - .00057 * t_centuries / 2.) * (new Constants().pi / 180.) / 240.;
           double ra_rate  = m + n * Math.sin(ra) * Math.tan(dec);
           double dec_rate = n * Math.cos(ra);

           ra -= t_centuries * ra_rate * 100.;
           dec -= t_centuries * dec_rate * 100.;
           
           coordinates.v[0] = ra;
           coordinates.v[1] = dec;
           
           return coordinates;
        }

        /**
         * epoch_of_date_to_j2000_2()
         *
         * Reduction for preccession - approximate formula
         * Method to convert equatorial coordinates from epoch date to J2000.
         * Source "Astronomical Almanac 2005", site B19
         * For reduction to J2000
         *
         * IN:
         *  double mjd - Modified Julian date [days]
         *  double ra - R.A. in epoch date??
         *  doubel dec - declination in epoch date
         *
         *  OUT:
         *  Vector coordinates[2] - v.[0] - R.A. in J2000
         *                        - v.[1] - declination in J2000
         */

        public static Vector epoch_of_date_to_j2000_2( double mjd, double ra, double dec){
           Vector coordinates = new Vector(2);
           double jd = mjd + 2400000.5;
           double t_centuries = (jd - 2451545.) / 36525.;
           //corrected positions
           double ra_0, dec_0;
           //help positions
           double ra_m, dec_m;

           double m = Math.toRadians(1.2812323)*t_centuries + Math.toRadians(0.0003879)*t_centuries*t_centuries +
                   Math.toRadians(0.0000101)*t_centuries*t_centuries*t_centuries;
           double n = Math.toRadians(0.556753)*t_centuries - Math.toRadians(0.0001185)*t_centuries*t_centuries -
                   Math.toRadians(0.0000116)*t_centuries*t_centuries*t_centuries;

           ra_m = ra - 0.5*(m + n*Math.sin(ra)*Math.tan(dec));
           dec_m = dec - 0.5*n*Math.cos(ra_m);

           ra_0 = ra - m - n*Math.sin(ra_m)*Math.tan(dec_m);
           dec_0 = dec - n*Math.cos(ra_m);

           coordinates.v[0] = ra_0;
           coordinates.v[1] = dec_0;

           return coordinates;
        }
        /**
         * epoch_of_date_from_j2000()
         *
         * Reduction for preccession - approximate formula
         * Method to convert equatorial coordinates from J2000 to epoch date.
         * Source "Astronomical Almanac 2005", site B19
         * For reduction from J2000
         *
         * IN:
         *  double mjd - Modified Julian date [days]
         *  double ra - R.A. in epoch date??
         *  doubel dec - declination in epoch date
         *
         *  OUT:
         *  Vector coordinates[2] - v.[0] - R.A. in J2000
         *                        - v.[1] - declination in J2000
         */

        public static Vector epoch_of_date_from_j2000( double mjd, double ra_0, double dec_0){
           Vector coordinates = new Vector(2);
           double jd = mjd + 2400000.5;
           double t_centuries = (jd - 2451545.) / 36525.;
           //corrected positions
           double ra, dec;
           //help positions
           double ra_m, dec_m;

           double m = Math.toRadians(1.2812323)*t_centuries + Math.toRadians(0.0003879)*t_centuries*t_centuries +
                   Math.toRadians(0.0000101)*t_centuries*t_centuries*t_centuries;
           double n = Math.toRadians(0.556753)*t_centuries - Math.toRadians(0.0001185)*t_centuries*t_centuries -
                   Math.toRadians(0.0000116)*t_centuries*t_centuries*t_centuries;

           ra_m = ra_0 + 0.5*(m + n*Math.sin(ra_0)*Math.tan(dec_0));
           dec_m = dec_0 + 0.5*n*Math.cos(ra_m);

           ra = ra_0 + m + n*Math.sin(ra_m)*Math.tan(dec_m);
           dec = dec_0 + n*Math.cos(ra_m);

           coordinates.v[0] = ra;
           coordinates.v[1] = dec;

           return coordinates;
        }

        /**
         * fromHorizToEquatCoord_New
         *
         * Method to transforme horizontal coordinates to equatorial for given
         * time and observer geodetic position. - with precesion correlation
         *
         *  IN:
         *  Az [rad]
         *  h  [rad]
         *  lon  [rad]  - obs. pos. - longitude
         *  lat  [rad]  - obs. pos. - latitude
         *  mjd [days]
         *
         *  OUT:
         *  Vector:
         *  v[0] = R.A. [rad]
         *  v[1] = dec  [rad]
         */
        public static Vector fromHorizToEquatCoord_New(double Az, double h,
                double lon, double lat, double mjd){

            double gmst = new Time().getGMST(mjd);

            //get equatorial coordinates for epoch mjd
            Vector vector = fromHorizontalToEquat(Az,h,gmst,lon,lat);
            //System.out.println("ra corr  " + Math.toDegrees(vector.v[0]));
            //System.out.println("dec corr " + Math.toDegrees(vector.v[1]));

            //equatorial coordinates for epoch 2000
            Vector vector2 = epoch_of_date_to_j2000_2(mjd, vector.v[0], vector.v[1]);
            //System.out.println("ra out   " + Math.toDegrees(vector2.v[0]));
            //System.out.println("dec out  " + Math.toDegrees(vector2.v[1]));

            return vector2;
        }

        /**
         * fromEquatToHorizCoord_New
         *
         * Method to transforme equatorial coordinates to horizontal for given
         * time and observer geodetic position. - with precession correlation
         *
         *  IN:
         *  ra [rad]
         *  dec  [rad]
         *  lon  [rad]  - obs. pos. - longitude
         *  lat  [rad]  - obs. pos. - latitude
         *  mjd [days]
         *
         *  OUT:
         *  Vector:
         *  v[0] = Az [rad]
         *  v[1] = h  [rad]
         */
        public static Vector fromEquatToHorizCoord_New(double ra, double dec,
                double lon, double lat, double mjd){

            double gmst = new Time().getGMST(mjd);

            //get equatorial coordinates to epoch mjd
            Vector vector = epoch_of_date_from_j2000(mjd, ra, dec);
            //System.out.println("ra corr  " + Math.toDegrees(vector.v[0]));
            //System.out.println("dec corr " + Math.toDegrees(vector.v[1]));

            //equatorial coordinates for epoch 2000
            Vector vector2 = fromEquatToHorizontal(vector.v[0], vector.v[1], gmst, lon, lat);
            //Vector vector2 = fromEquatToHorizontal(ra, dec, gmst, lon, lat);
            //System.out.println("Az   " + Math.toDegrees(vector2.v[0]));
            //System.out.println("h    " + Math.toDegrees(vector2.v[1]));

            return vector2;
        }

        /**
         * getMeanObliquityOfEcliptic()
         *
         * Method to compute angle between equatorial and ecliptical planes.
         * See "Explanatory Supplement to the Astronomical Almanac
         * by P. Kenneth Seidelmann, 2006 - formula no. 3.222-1 and 3.222-2
         *
         * IN:
         *      double mjd - [days]
         *
         * OUT:
         *      double mq - angle [rad]
         */
        public static double getMeanObliquityOfEcliptic(double mjd){
            double mq;
            //diffrence between J2000 and current epoch [centuries]
            double t = (mjd - 51544.5)/(double)36525;

            mq = Math.toRadians(23.439291111111111111111111111111) -
                 Math.toRadians(0.013004166666666666666666666666667)*t -
                 Math.toRadians(1.6388888888888888888888888888889e-7)*t*t +
                 Math.toRadians(5.0361111111111111111111111111111e-7)*t*t*t;

            return mq;
        }

        /**
         * fromGeocentricToHeliocentric()
         *
         * Method to trasforme geocentric position from geocentric position system
         * to heliocentric position system and from equatorial to ecliptical coordinates.
         *
         * IN:
         * Vector vectorGeo - geocentric position vector of object
         * Tdouble mjd - [day]
         *
         * OUT:
         * Vector vectorHel - heliocentric position vector
         *
         */
        public static Vector fromGeocentricToHeliocentric(Vector vectorGeo, double mjd){
            //mean obliquity of ecliptic
            double meanObl = getMeanObliquityOfEcliptic(mjd);
            //rotation matrix from equatorial coordinates to ecliptical
            Matrix rotMatEqToEcl = new Matrix(3,3);
            //geocentric position vector of object in ecliptic coordinates
            Vector vectorGeoEcl = new Vector(3);
            //heliocentric position vector of object in ecliptic coordinates
            Vector vectorHeliocEcl = new Vector(3);
            //geocentric position vector of Sun
            Vector vectorSunGeoc = new Vector(3);
            //geocentric ecliptical position vector of Sun
            Vector vectorSunHelioc = new Vector(3);
            //heliocentric ecliptical position vector of Earth
            Vector vectorEarthHelioc = new Vector(3);

            //get the rotation matrix
            rotMatEqToEcl = rotMatEqToEcl.R_x(meanObl);
            //get geocentric position vector of object in ecliptic coordinates
            vectorGeoEcl = rotMatEqToEcl.matrixMultiplyVector(rotMatEqToEcl, vectorGeo);
            //vectorGeoEcl = rotMatEqToEcl.vectorMultiplyMatrix(rotMatEqToEcl, vectorGeo);

            //get Sun geocentric equatorial position vector
            vectorSunGeoc = new Sun().getSunPosition(mjd);
            //System.out.println("Xgeo " + vectorSunGeoc.v[0]/1000 + " km");
            //System.out.println("Ygeo " + vectorSunGeoc.v[1]/1000 + " km");
            //System.out.println("Zgeo " + vectorSunGeoc.v[2]/1000 + " km");

            //get Sun geocentric ecliptical position vector
            vectorSunHelioc = rotMatEqToEcl.matrixMultiplyVector(rotMatEqToEcl, vectorSunGeoc);
            //System.out.println("Xecl " + vectorSunHelioc.v[0]/1000 + " km");
            //System.out.println("Yecl " + vectorSunHelioc.v[1]/1000 + " km");
            //System.out.println("Zecl " + vectorSunHelioc.v[2]/1000 + " km");
            //vectorSunHelioc = rotMatEqToEcl.vectorMultiplyMatrix(rotMatEqToEcl, vectorSunGeoc);
            //get heliocentric ecliptical position of Earth
            vectorEarthHelioc = vectorSunHelioc.multiplyVector(vectorSunHelioc, -1);

            //get heliocentric position vector of object in ecliptic coordinates
            vectorHeliocEcl =  vectorHeliocEcl.subtractVectors(vectorGeoEcl, vectorEarthHelioc);

            return vectorHeliocEcl;
        }
        /**
         * fromEquatorialToEcliptic()
         * ITS BAAAAAAAAAAAD!!!!!!!!!!!!!!!!!!!!!!!!!
         * Method to transform vector from equatorial plane reference system to ecliptic plane.
         *
         * IN:
         * Vector vectorGeo - geocentric position vector of object
         * Tdouble mjd - [day]
         *
         * OUT:
         * Vector vectorHel - heliocentric position vector
         *
         */
        
        public static Vector fromEquatorialToEcliptic(Vector vectorGeo, double mjd){
            //mean obliquity of ecliptic
            double meanObl = getMeanObliquityOfEcliptic(mjd);
            //System.out.println("meanObl " + Math.toDegrees(meanObl));
            //rotation matrix from equatorial coordinates to ecliptical
            Matrix rotMatEqToEcl = new Matrix(3,3);
            //geocentric position vector of object in ecliptic coordinates
            Vector vectorGeoEcl = new Vector(3);
            //heliocentric position vector of object in ecliptic coordinates
            Vector vectorHeliocEcl = new Vector(3);
            //geocentric position vector of Sun
            Vector vectorSunGeoc = new Vector(3);
            //geocentric ecliptical position vector of Sun
            Vector vectorSunHelioc = new Vector(3);
            //heliocentric ecliptical position vector of Earth
            Vector vectorEarthHelioc = new Vector(3);

            //get the rotation matrix
            rotMatEqToEcl = rotMatEqToEcl.R_x(meanObl);
            //get geocentric position vector of object in ecliptic coordinates
            vectorGeoEcl = rotMatEqToEcl.matrixMultiplyVector(rotMatEqToEcl, vectorGeo);
            //vectorGeoEcl = rotMatEqToEcl.vectorMultiplyMatrix(rotMatEqToEcl, vectorGeo);

            //get Sun geocentric equatorial position vector
            //vectorSunGeoc = new Sun().getSunPosition(mjd);
            //System.out.println("Xgeo " + vectorSunGeoc.v[0]/1000 + " km");
            //System.out.println("Ygeo " + vectorSunGeoc.v[1]/1000 + " km");
            //System.out.println("Zgeo " + vectorSunGeoc.v[2]/1000 + " km");

            //get Sun geocentric ecliptical position vector
            //vectorSunHelioc = rotMatEqToEcl.matrixMultiplyVector(rotMatEqToEcl, vectorSunGeoc);
            //System.out.println("Xecl " + vectorSunHelioc.v[0]/1000 + " km");
            //System.out.println("Yecl " + vectorSunHelioc.v[1]/1000 + " km");
            //System.out.println("Zecl " + vectorSunHelioc.v[2]/1000 + " km");
            //vectorSunHelioc = rotMatEqToEcl.vectorMultiplyMatrix(rotMatEqToEcl, vectorSunGeoc);
            //get heliocentric ecliptical position of Earth
            //vectorEarthHelioc = vectorSunHelioc.multiplyVector(vectorSunHelioc, -1);

            //get heliocentric position vector of object in ecliptic coordinates
            //vectorHeliocEcl =  vectorHeliocEcl.subtractVectors(vectorGeoEcl, vectorEarthHelioc);

            return vectorGeoEcl;
        }

        /**
         * fromEquatorialToEcliptic
         *
         * Method to calculate the ecliptic coordinates from equatorial coordinates.
         *
         * INPUT:
         *  Vector equatVec - 0-RA, 1-DEC
         *
         * OUTPUT:
         *  Vector eclipVec - 0-lon, 1-lat
         */
        public Vector fromEquatorialToEcliptic2(Vector equatVec, double mjd){
            Vector unitEqVec = getEqUnitVector(equatVec);
            //mean obliquity of ecliptic
            double meanObl = getMeanObliquityOfEcliptic(mjd);
            //get the rotation matrix
            Matrix rotMatEqToEcl = new Matrix(3,3);
            rotMatEqToEcl = rotMatEqToEcl.R_x(meanObl);
            //get geocentric position vector of object in ecliptic coordinates
            Vector eclipVec = rotMatEqToEcl.matrixMultiplyVector(rotMatEqToEcl, unitEqVec);
            //get the direction values - ecliptic longituded and latitude
            eclipVec = getEquatFromUnitVec(eclipVec);
            while (eclipVec.v[0]<0) eclipVec.v[0] = eclipVec.v[0] + Math.PI*2;
            while (eclipVec.v[0]>=Math.PI*2) eclipVec.v[0] = eclipVec.v[0] - Math.PI*2;
            return eclipVec;
        }
        /**
         * fromEquatorialToEcliptic
         *
         * Method to calculate the ecliptic coordinates from equatorial coordinates.
         *
         * INPUT:
         *  Vector equatVec - 0-RA, 1-DEC
         *
         * OUTPUT:
         *  Vector eclipVec - 0-lon, 1-lat
         */
        public Vector fromEclipticToEquatorial2(Vector eclipVec, double mjd){
            Vector unitEclVec = getEqUnitVector(eclipVec);
            //mean obliquity of ecliptic
            double meanObl = getMeanObliquityOfEcliptic(mjd);
            //get the rotation matrix
            Matrix rotMatEclToEq = new Matrix(3,3);
            rotMatEclToEq = rotMatEclToEq.R_x(-meanObl);
            //get geocentric position vector of object in ecliptic coordinates
            Vector equatVec = rotMatEclToEq.matrixMultiplyVector(rotMatEclToEq, unitEclVec);
            //get the direction values - ecliptic longituded and latitude
            equatVec = getEquatFromUnitVec(equatVec);
            while (equatVec.v[0]<0) equatVec.v[0] = equatVec.v[0] + Math.PI*2;
            while (equatVec.v[0]>=Math.PI*2) equatVec.v[0] = equatVec.v[0] - Math.PI*2;

            return equatVec;
        }

        /**
         * getEquatFromUnitVec() - !!!!!!!!!!!!!!!!!!!!!!!!! CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!quadrants!!!!!
         *
         * Method to calculate the equatorial coordinates from unit vector. Reverse
         * method to the method getEqUnitVector().
         *
         * INPUT:
         *  Vector unitVec - unit vector
         *
         * OUTPUT:
         *  Vector equatVec - equatorial coordinates RA and DEC [rad]
         */
        public Vector getEquatFromUnitVec(Vector unitVec){
            Vector equatVec = new Vector(3);
            //DEC
            equatVec.v[1] = Math.asin(unitVec.v[2]);
            //RA
            //equatVec.v[0] = Math.asin(unitVec.v[1]/Math.cos(equatVec.v[1]));
            equatVec.v[0] = Math.atan2(unitVec.v[1],unitVec.v[0]);
            return equatVec;
        }


        /**
         * getEqUnitVector()
         *
         * Method to calculate the unit vector in direction of the meteor point's
         * right ascension and declination.
         *
         * INPUT:
         *  Vector eqCoorVec(2) - RA, DEC [rad]
         *
         * OUTPUT:
         *  Vector unitVec - Xi, Eta, Zeta [-]
         */
        public Vector getEqUnitVector(Vector eqCoorVec){
            Vector unitVec = new Vector(3);
            //Xi
            unitVec.v[0] = Math.cos(eqCoorVec.v[1])*Math.cos(eqCoorVec.v[0]);
            //Eta
            unitVec.v[1] = Math.cos(eqCoorVec.v[1])*Math.sin(eqCoorVec.v[0]);
            //Zeta
            unitVec.v[2] = Math.sin(eqCoorVec.v[1]);
            return unitVec;
        }

        /**
         * getEclipticInclAndNodeFromGeocPosVectors()
         *
         * Method to compute ecliptical inclination and RAN from 2 geocentric position vectors (1st position is sooner),
         *
         * Input:
         *      Vector geocPosVec_1 - 1st geocentric position vector [m]
         *      Vector geocPosVec_2 - 2nd geocentric position vector [m]
         *      double mjd_1        - get time for 1st position  [days]
         *      double mjd_2        - get time for 1st position  [days]
         *
         * Output:
         *      double elements[] - ecliptical inclination-elements[0] [rad] and RAN elements[1]
         */
        public static double[] getEclipticInclAndNodeFromGeocPosVectors(Vector geocPosVec_1,Vector geocPosVec_2,double mjd_1,double mjd_2){
            //ecliptical elements
            double elements[] = new double[2];
            //conversion from geocentric vectors to heliocentric vectors
            Vector posVecHelio_1 = new Vector(3);
            posVecHelio_1 = Transformation.fromGeocentricToHeliocentric(geocPosVec_1,mjd_1);
            Vector posVecHelio_2 = new Vector(3);
            posVecHelio_2 = Transformation.fromGeocentricToHeliocentric(geocPosVec_2,mjd_2);
            //get amgular momentum unit vector
            Vector angMomUnitVec = new Vector(3);
            angMomUnitVec = Vector.getVectorProduct(posVecHelio_1, posVecHelio_2);
            angMomUnitVec = Vector.multiplyVector(angMomUnitVec, (double)1/Vector.getSize(angMomUnitVec));

            //get elements
            elements[0] = Math.atan(Math.sqrt(angMomUnitVec.v[0]*angMomUnitVec.v[0]+
                                    angMomUnitVec.v[1]*angMomUnitVec.v[1])/angMomUnitVec.v[2]);
            if(elements[0]<0)elements[0] = elements[0]+Math.PI;
            //System.out.println("TEST Ecl i    : " + Math.toDegrees(eclip_i));
            elements[1] = Math.atan2(angMomUnitVec.v[0],-angMomUnitVec.v[1]);
            if(elements[1]<0)elements[1] = elements[1]+2*Math.PI;
            //System.out.println("TEST Ecl O   : " + Math.toDegrees(eclip_Omega));
            return elements;
        }

        /**
         * getEquatorialInclAndNodeFromGeocPosVectors()
         *
         * Method to compute ecliptical inclination and RAN from 2 geocentric position vectors (1st position is sooner),
         *
         * Input:
         *      Vector geocPosVec_1 - 1st geocentric position vector [m]
         *      Vector geocPosVec_2 - 2nd geocentric position vector [m]
         *      double mjd_1        - get time for 1st position  [days]
         *      double mjd_2        - get time for 1st position  [days]
         *
         * Output:
         *      double elements[] - ecliptical inclination-elements[0] [rad] and RAN elements[1]
         */
        public static double[] getEquatorialInclAndNodeFromGeocPosVectors(Vector geocPosVec_1,Vector geocPosVec_2,double mjd_1,double mjd_2){
            //ecliptical elements
            double elements[] = new double[2];
            //get amgular momentum unit vector
            Vector angMomUnitVec = new Vector(3);
            angMomUnitVec = Vector.getVectorProduct(geocPosVec_1, geocPosVec_2);
            angMomUnitVec = Vector.multiplyVector(angMomUnitVec, (double)1/Vector.getSize(angMomUnitVec));

            //get elements
            elements[0] = Math.atan(Math.sqrt(angMomUnitVec.v[0]*angMomUnitVec.v[0]+
                                    angMomUnitVec.v[1]*angMomUnitVec.v[1])/angMomUnitVec.v[2]);
            if(elements[0]<0)elements[0] = elements[0]+Math.PI;
            //System.out.println("TEST Ecl i    : " + Math.toDegrees(eclip_i));
            elements[1] = Math.atan2(angMomUnitVec.v[0],-angMomUnitVec.v[1]);
            if(elements[1]<0)elements[1] = elements[1]+2*Math.PI;
            //System.out.println("TEST Ecl O   : " + Math.toDegrees(eclip_Omega));
            return elements;
        }
        
        /**
         * getGalacticCoordinates()
         * 
         * Method to comvert equatorial to galactic coordinates.
         * Source: http://scienceworld.wolfram.com/astronomy/GalacticCoordinates.html
         * 
         * INPUT:
         *  Vector eqCoorVec(2) - RA, DEC [rad]
         * 
         * Output:
         *  Vector galCoorVec(2) - LON, LAT [rad]
         */
        public static Vector getGalacticCoordinates(Vector eqCoorVec){
            Vector galCoorVec = new Vector(2);
            double dec_p = Math.toRadians(62.6);
            double ra_p = Math.toRadians(282.25);
            
            galCoorVec.v[1] = Math.asin(Math.sin(eqCoorVec.v[1])*Math.cos(dec_p) - 
                    Math.cos(eqCoorVec.v[1])*Math.sin(eqCoorVec.v[0]-ra_p)*Math.sin(dec_p));
            //getting angle sin(l-33deg)
            double sin_glon = Math.asin(Math.sin(eqCoorVec.v[1])*Math.sin(dec_p)+
                    Math.cos(eqCoorVec.v[1])*Math.sin(eqCoorVec.v[0]-ra_p)*Math.cos(dec_p));
            //getting angle cos(l-33deg);
            double cos_glon = Math.acos(Math.cos(eqCoorVec.v[1])*Math.cos(eqCoorVec.v[0]-ra_p));
            //get the galactic longitude
            galCoorVec.v[0] = Math.atan2(sin_glon,cos_glon) + Math.toRadians(33.0);
            
            return galCoorVec;
        }
}
