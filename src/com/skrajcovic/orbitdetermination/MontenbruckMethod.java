/*
*
*OrbitDetermination.java
* 
*Ucel (purpose) :
*
*	Vypocet keplerovskych drahovych elementov (Keplers elements computation) na zaklade 3 pozorovani
*       Determination of keplerian elements from 3 observations
*
*
* 2007/10/30 - Jiri Silha
*
*/

package com.skrajcovic.orbitdetermination;

//

import com.skrajcovic.orbitdetermination.compute.*;

// Deklaracia triedy (class declaration)
//

/** Orbit determination from 3 observations*/
public class MontenbruckMethod{
	
	/**
	*   Vypocet drahovych elementov na zaklade 3 pozorovani v danych casoch, Montenbruckova metoda, 
	*   citaj 'Satellite orbits' , str. 43 - 46
	*   (Computes orbital elements from 3 observations and  associated times ), Montenbruck Method,
	*    read 'Satellite orbits', sites 43 - 46
	*
	* Input/Output:
	*
	*   Geodetic coordinates of observatory   
	*   Observation 1
	*   Observation 2
	*   Observation 3
	*
	*   <return>  Keplerian elements (a,e,i,Omega,omega,M)
	*               a      Semimajor axis 
	*               e      Eccentricity 
	*               i      Inclination [rad]
	*               Omega  Longitude of the ascending node [rad]
	*               omega  Argument of pericenter  [rad]
	*               M      Mean anomaly  [rad] at time t_1
	*
	* Notes:
	*
	*   The function cannot be used with state vectors describing a circular
	*   or non-inclined orbit.
	*
	*/
	
	public static Kepler getElementsMM(Geodetic geodetic, ObservationHor obs1, ObservationHor obs2, ObservationHor obs3){
		//konstanty
		Constants c = new Constants();
		Transformation transf = new Transformation();
		double minimum = 1e-13;
		
		//koeficienty n1, n3 
		double n1, n3;
		
		//Polohovy vektor pozorovatela
		//Vector sta=new Vector(3);
		
		//Geodeticke suradnice - TREBA!!! zadat zemepisnu dlzku, sirku a nadm. vysku observatoria
		Geodetic geodObs = geodetic;//new Geodetic(Math.toRadians(18.6107),Math.toRadians(53.0217),0);
		
		//prevod suradnic pozorovatela z geodetickych do geocentrickych suradnic
		Vector obsGeocPosition = new Vector(3);
		obsGeocPosition = transf.fromGeodToGeoc(geodObs.lon, geodObs.lat, geodObs.altitude, c.R_Earth, c.f_Earth);
		
		//v kazkom case ma pozorovatel ine ekvatorialne suradnice (os x smeruje do jarneho bodu)
		Vector obsPosition[] = new Vector[3];
		obsPosition[0] = new Vector(3);
		obsPosition[1] = new Vector(3);
		obsPosition[2] = new Vector(3);
		
		//POZOROVANIA - CAS, AZIMUT, VYSKA
		//horizontalne suradnice
		double azi[] = new double[3];
		double elev[] = new double[3];
		Time time[] = new Time[3];
		
		//1. Pozorovanie
		time[0] = obs1.t;//new Time(0,2007,10,20,5,35,47.0);
		azi[0] = obs1.azim;//Math.toRadians(295.7);
		elev[0] = obs1.elev;//Math.toRadians(63.8);
		
		//2. Pozorovanie
		time[1] = obs2.t;//new Time(0,2007,10,20,6,35,47.0);
		azi[1] = obs2.azim;//Math.toRadians(194.6);
		elev[1] = obs2.elev;//Math.toRadians(80);
		
		//3. Pozorovanie
		time[2] = obs3.t;//new Time(0,2007,10,20,7,35,47.0);
		azi[2] = obs3.azim;//Math.toRadians(166.2);
		elev[2] = obs3.elev;//Math.toRadians(51.1);
		
		//Matica na prevod do geocentrickych suradnic
		Matrix ltcMat = new Matrix(3,3);
		ltcMat = transf.getLtcMatrix(geodObs.lon,geodObs.lat);
		
		//Matica na prevod z ekvatorialnych do geocentrickych suradnic 
		Matrix[] u = new Matrix[3];
		//vektor, ktory je prevod z horizontalnych do tangencialnych suradnic 
		Vector[] s = new Vector[3];
		
		for(int i=0;i<3;i++){
			//inicializacia translacnej matice z ekvat. => geocentrickych, pouzije sa jej transp. matica
			u[i] = new Matrix(3,3);
			u[i] = u[i].R_z(time[i].getGMST(time[i].getMjd(time[i])));
			
			//pozicny vektor pozorovatela v ekvatorilanych suradniciach v danom case pozorovanie
			obsPosition[i] = u[i].matrixMultiplyVector(u[i].transposeMatrix(u[i]),obsGeocPosition);
		}
		
		//jednotkove vektory pozorovani
		Vector e[] = new Vector[3];
		for(int i=0; i<3;i++){
			s[i] = new Vector(3);
			s[i] = s[i].fromPolarToCartesian(c.pi/2-azi[i],elev[i],1);
			e[i] = ltcMat.matrixMultiplyVector(u[i].transposeMatrix(u[i]),u[i].matrixMultiplyVector(ltcMat.transposeMatrix(ltcMat),s[i]));
		}
		
		//d koeficienty pre vypocet vzdialenosti sat-poz
		double D,d_11,d_12,d_13,d_21,d_22,d_23,d_31,d_32,d_33;
		Vector d[] = new Vector[3];
		d[0] = transf.getVector_di(e[1],e[2]);
		d[1] = transf.getVector_di(e[2],e[0]);
		d[2] = transf.getVector_di(e[0],e[1]);
		
		D = transf.getCoefficient_D(e[0],e[1],e[2]);
		d_11 = transf.getCoefficient_Dij(d[0],obsPosition[0]);
		d_12 = transf.getCoefficient_Dij(d[0],obsPosition[1]);
		d_13 = transf.getCoefficient_Dij(d[0],obsPosition[2]);
		d_21 = transf.getCoefficient_Dij(d[1],obsPosition[0]);
		d_22 = transf.getCoefficient_Dij(d[1],obsPosition[1]);
		d_23 = transf.getCoefficient_Dij(d[1],obsPosition[2]);
		d_31 = transf.getCoefficient_Dij(d[2],obsPosition[0]);
		d_32 = transf.getCoefficient_Dij(d[2],obsPosition[1]);
		d_33 = transf.getCoefficient_Dij(d[2],obsPosition[2]);
		
		//koeficienty ro - vzdialenosti poz. - sat.
		double ro[] = new double[3];
		
		//ekvatorialne suradnice satelitov v danom case
		Vector[] satPosition = new Vector[3];
		
		//deklaracia triedy kepler
		Kepler kepler = new Kepler(0,0,0,0,0,0);

		//pomocne premenne
		int maxit = 0;
		double n1_help;
		double eta_1,eta_2,eta_3;
		double tau_1,tau_2,tau_3;
		
		//pociatocna podmienka pre iteraciu
		n1 = (time[2].getMjd(time[2])-time[1].getMjd(time[1]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));
		n1_help = n1;
		n3 = (time[1].getMjd(time[1])-time[0].getMjd(time[0]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));
		
		//vypocet koeficientov tau_1, tau_2, tau_3
		tau_1 = kepler.getTau(c.GM_Earth,time[1].getMjd(time[1]),time[2].getMjd(time[2]));
		tau_2 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[2].getMjd(time[2]));
		tau_3 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[1].getMjd(time[1]));
		
		//vypocet vzdialenosti pre dane poc. podmienky
		ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
		ro[1] = (1/D)*(n1*d_21-d_22+n3*d_23);
		ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);
		
		//pole aktualnych vektorov - jednotkovy * vzdialenost
		Vector actualPos[] = new Vector[3];
		
		//vypocet polohovych vektorov
		for(int i = 0;i < 3;i++){
			actualPos[i] = new Vector(3);
			actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
			satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
		}			
		
		//pole double pre kontrolu
		double itRo[][] = new double[3][3];
		
		
		while(true){			
			//vypocet koeficientov eta_1,eta_2,eta_3, vsetko o nich je vo vzahoch 2.98-2.108,
			eta_1 = kepler.getEta(satPosition[1],satPosition[2],tau_1);
			eta_2 = kepler.getEta(satPosition[0],satPosition[2],tau_2);
			eta_3 = kepler.getEta(satPosition[0],satPosition[1],tau_3);
			
			//vypocet koeficientov n1 a n3
			n1 = transf.getCoefficient_n1(eta_1,eta_2,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));
			n3 = transf.getCoefficient_n3(eta_2,eta_3,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));
			
			//vypocet vzdialenosti sat. - poz
			ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
			ro[1] = 1/D*(n1*d_21-d_22+n3*d_23);
			ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);
			
			//vypocet polohovych vektorov
			for(int i = 0;i<3;i++){
				actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
				satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
			}	
			
			//ukladame ziskane ro cez iteracie, na kontrolu
			itRo[0][0] = itRo[1][0];
			itRo[0][1] = itRo[1][1];
			itRo[0][2] = itRo[1][2];
			
			itRo[1][0] = itRo[2][0];
			itRo[1][1] = itRo[2][1];
			itRo[1][2] = itRo[2][2];
			
			itRo[2][0] = ro[0];
			itRo[2][1] = ro[1];
			itRo[2][2] = ro[2];
			
			//test iteracie
			if(Math.abs(n1-n1_help)<minimum) {
				//System.out.println("No. of it: "+maxit+"\n");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i = 0;i<3;i++){
						for(int j = 0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j = 0;j<3;j++){
						//System.out.println(ro[j]);
					}
				}
				
				//if((ro[0]<0)||(ro[1]<0)||(ro[2]<0)) System.out.println("Negative distance!!!");
				break; 
			}
			
			if(maxit > 300) {
				//System.out.println("WARNING: Convergence problems in ro!!!");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i=0;i<3;i++){
						for(int j=0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j=0;j<3;j++){
						//System.out.println(ro[j]);
					}
				}
				
				break;
			}
			
			n1_help=n1;
			maxit++;
		}
		
                //WARNING
                //if((ro[0] < 0)||(ro[1] < 0)||(ro[2] < 0)) System.out.println("WARNING: negative ro!!!");
		//uz mame 3 pozicne vektory, treba vyratat elementy
		Constants consta = new Constants();
		kepler = kepler.getElements(consta.GM_Earth,time[0].getMjd(time[0]),time[2].getMjd(time[2]),satPosition[0],satPosition[2]);
		//System.out.println("x " + satPosition[1].v[0]);
		//System.out.println("y " + satPosition[1].v[1]);
		//System.out.println("z " + satPosition[1].v[2]);
		//System.out.println("x " + satPosition[2].v[0]);
		//System.out.println("y " + satPosition[2].v[1]);
		//System.out.println("z " + satPosition[2].v[2]);
		return kepler;
	}

        /**
	*   Vypocet drahovych elementov na zaklade 3 pozorovani v danych casoch, Montenbruckova metoda,
	*   citaj 'Satellite orbits' , str. 43 - 46. Pozorovatel moze byt na roznych miestach.
	*   (Computes orbital elements from 3 observations and  associated times ), Montenbruck Method,
	*    read 'Satellite orbits', sites 43 - 46. Observer can be on 3 different places.
	*
	* Input/Output:
	*
	*   Geodetic coordinates of observatory - array[3]
	*   Observation 1
	*   Observation 2
	*   Observation 3
	*
	*   <return>  Keplerian elements (a,e,i,Omega,omega,M)
	*               a      Semimajor axis
	*               e      Eccentricity
	*               i      Inclination [rad]
	*               Omega  Longitude of the ascending node [rad]
	*               omega  Argument of pericenter  [rad]
	*               M      Mean anomaly  [rad] at time t_1
	*
	* Notes:
	*
	*   The function cannot be used with state vectors describing a circular
	*   or non-inclined orbit.
	*
	*/

	public static Kepler getElementsMM_2(Geodetic geodetic[], ObservationHor obs1, ObservationHor obs2, ObservationHor obs3){
		//konstanty
		Constants c = new Constants();
		Transformation transf = new Transformation();
		double minimum = 1e-15;

		//koeficienty n1, n3
		double n1, n3;

		//Polohovy vektor pozorovatela
		//Vector sta=new Vector(3);

		//Geodeticke suradnice - TREBA!!! zadat zemepisnu dlzku, sirku a nadm. vysku observatoria
		Geodetic geodObs[] = geodetic;//new Geodetic(Math.toRadians(18.6107),Math.toRadians(53.0217),0);

		//prevod suradnic pozorovatela z geodetickych do geocentrickych suradnic
		Vector obsGeocPosition[] = new Vector[3];
                for(int i=0;i<3;i++){
                    obsGeocPosition[i] = new Vector(3);
                    obsGeocPosition[i] = transf.fromGeodToGeoc(geodObs[i].lon, geodObs[i].lat, geodObs[i].altitude, c.R_Earth, c.f_Earth);
                }
		//v kazkom case ma pozorovatel ine ekvatorialne suradnice (os x smeruje do jarneho bodu)
		Vector obsPosition[] = new Vector[3];
		obsPosition[0] = new Vector(3);
		obsPosition[1] = new Vector(3);
		obsPosition[2] = new Vector(3);

		//POZOROVANIA - CAS, AZIMUT, VYSKA
		//horizontalne suradnice
		double azi[] = new double[3];
		double elev[] = new double[3];
		Time time[] = new Time[3];

		//1. Pozorovanie
		time[0] = obs1.t;//new Time(0,2007,10,20,5,35,47.0);
		azi[0] = obs1.azim;//Math.toRadians(295.7);
		elev[0] = obs1.elev;//Math.toRadians(63.8);

		//2. Pozorovanie
		time[1] = obs2.t;//new Time(0,2007,10,20,6,35,47.0);
		azi[1] = obs2.azim;//Math.toRadians(194.6);
		elev[1] = obs2.elev;//Math.toRadians(80);

		//3. Pozorovanie
		time[2] = obs3.t;//new Time(0,2007,10,20,7,35,47.0);
		azi[2] = obs3.azim;//Math.toRadians(166.2);
		elev[2] = obs3.elev;//Math.toRadians(51.1);

		//Matica na prevod do geocentrickych suradnic
		Matrix ltcMat[] = new Matrix[3];
		for(int i=0;i<3;i++){
                    ltcMat[i] = new Matrix(3,3);
                    ltcMat[i] = transf.getLtcMatrix(geodObs[i].lon,geodObs[i].lat);
                }

		//Matica na prevod z ekvatorialnych do geocentrickych suradnic
		Matrix[] u = new Matrix[3];
		//vektor, ktory je prevod z horizontalnych do tangencialnych suradnic
		Vector[] s = new Vector[3];

		for(int i=0;i<3;i++){
			//inicializacia translacnej matice z ekvat. => geocentrickych, pouzije sa jej transp. matica
			u[i] = new Matrix(3,3);
			u[i] = u[i].R_z(time[i].getGMST(time[i].getMjd(time[i])));

			//pozicny vektor pozorovatela v ekvatorilanych suradniciach v danom case pozorovanie
			obsPosition[i] = u[i].matrixMultiplyVector(u[i].transposeMatrix(u[i]),obsGeocPosition[i]);
		}

		//jednotkove vektory pozorovani
		Vector e[] = new Vector[3];
		for(int i=0; i<3;i++){
			s[i] = new Vector(3);
			s[i] = s[i].fromPolarToCartesian(c.pi/2-azi[i],elev[i],1);
			e[i] = ltcMat[i].matrixMultiplyVector(u[i].transposeMatrix(u[i]),u[i].matrixMultiplyVector(ltcMat[i].transposeMatrix(ltcMat[i]),s[i]));
		}

		//d koeficienty pre vypocet vzdialenosti sat-poz
		double D,d_11,d_12,d_13,d_21,d_22,d_23,d_31,d_32,d_33;
		Vector d[] = new Vector[3];
		d[0] = transf.getVector_di(e[1],e[2]);
		d[1] = transf.getVector_di(e[2],e[0]);
		d[2] = transf.getVector_di(e[0],e[1]);

		D = transf.getCoefficient_D(e[0],e[1],e[2]);
		d_11 = transf.getCoefficient_Dij(d[0],obsPosition[0]);
		d_12 = transf.getCoefficient_Dij(d[0],obsPosition[1]);
		d_13 = transf.getCoefficient_Dij(d[0],obsPosition[2]);
		d_21 = transf.getCoefficient_Dij(d[1],obsPosition[0]);
		d_22 = transf.getCoefficient_Dij(d[1],obsPosition[1]);
		d_23 = transf.getCoefficient_Dij(d[1],obsPosition[2]);
		d_31 = transf.getCoefficient_Dij(d[2],obsPosition[0]);
		d_32 = transf.getCoefficient_Dij(d[2],obsPosition[1]);
		d_33 = transf.getCoefficient_Dij(d[2],obsPosition[2]);

		//koeficienty ro - vzdialenosti poz. - sat.
		double ro[] = new double[3];

		//ekvatorialne suradnice satelitov v danom case
		Vector[] satPosition = new Vector[3];

		//deklaracia triedy kepler
		Kepler kepler = new Kepler(0,0,0,0,0,0);

		//pomocne premenne
		int maxit = 0;
		double n1_help;
		double eta_1,eta_2,eta_3;
		double tau_1,tau_2,tau_3;

		//pociatocna podmienka pre iteraciu
		n1 = (time[2].getMjd(time[2])-time[1].getMjd(time[1]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));
		n1_help = n1;
		n3 = (time[1].getMjd(time[1])-time[0].getMjd(time[0]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));

		//vypocet koeficientov tau_1, tau_2, tau_3
		tau_1 = kepler.getTau(c.GM_Earth,time[1].getMjd(time[1]),time[2].getMjd(time[2]));
		tau_2 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[2].getMjd(time[2]));
		tau_3 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[1].getMjd(time[1]));

		//vypocet vzdialenosti pre dane poc. podmienky
		ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
		ro[1] = (1/D)*(n1*d_21-d_22+n3*d_23);
		ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);

		//pole aktualnych vektorov - jednotkovy * vzdialenost
		Vector actualPos[] = new Vector[3];

		//vypocet polohovych vektorov
		for(int i = 0;i < 3;i++){
			actualPos[i] = new Vector(3);
			actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
			satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
		}

		//pole double pre kontrolu
		double itRo[][] = new double[3][3];


		while(true){
			//vypocet koeficientov eta_1,eta_2,eta_3, vsetko o nich je vo vzahoch 2.98-2.108,
			eta_1 = kepler.getEta(satPosition[1],satPosition[2],tau_1);
			eta_2 = kepler.getEta(satPosition[0],satPosition[2],tau_2);
			eta_3 = kepler.getEta(satPosition[0],satPosition[1],tau_3);

			//vypocet koeficientov n1 a n3
			n1 = transf.getCoefficient_n1(eta_1,eta_2,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));
			n3 = transf.getCoefficient_n3(eta_2,eta_3,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));

			//vypocet vzdialenosti sat. - poz
			ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
			ro[1] = 1/D*(n1*d_21-d_22+n3*d_23);
			ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);

			//vypocet polohovych vektorov
			for(int i = 0;i<3;i++){
				actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
				satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
			}

			//ukladame ziskane ro cez iteracie, na kontrolu
			itRo[0][0] = itRo[1][0];
			itRo[0][1] = itRo[1][1];
			itRo[0][2] = itRo[1][2];

			itRo[1][0] = itRo[2][0];
			itRo[1][1] = itRo[2][1];
			itRo[1][2] = itRo[2][2];

			itRo[2][0] = ro[0];
			itRo[2][1] = ro[1];
			itRo[2][2] = ro[2];

			//test iteracie
			if(Math.abs(n1-n1_help)<minimum) {
				//System.out.println("No. of it: "+maxit+"\n");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i = 0;i<3;i++){
						for(int j = 0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j = 0;j<3;j++){
						//System.out.println(ro[j]);
					}
				}

				//if((ro[0]<0)||(ro[1]<0)||(ro[2]<0)) System.out.println("Negative distance!!!");
				break;
			}

			if(maxit > 300) {
				//System.out.println("WARNING: Convergence problems in ro!!!");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i=0;i<3;i++){
						for(int j=0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j=0;j<3;j++){
						//System.out.println(ro[j]);
					}
				}

				break;
			}

			n1_help=n1;
			maxit++;
		}

                //WARNING
                //if((ro[0] < 0)||(ro[1] < 0)||(ro[2] < 0)) System.out.println("WARNING: negative ro!!!");
		//uz mame 3 pozicne vektory, treba vyratat elementy
		Constants consta = new Constants();
		kepler = kepler.getElements(consta.GM_Earth,time[1].getMjd(time[1]),time[2].getMjd(time[2]),satPosition[1],satPosition[2]);
		//System.out.println("x " + satPosition[1].v[0]);
		//System.out.println("y " + satPosition[1].v[1]);
		//System.out.println("z " + satPosition[1].v[2]);
		//System.out.println("x " + satPosition[2].v[0]);
		//System.out.println("y " + satPosition[2].v[1]);
		//System.out.println("z " + satPosition[2].v[2]);
		return kepler;
	}

        /**
         * getPosVectors()
	*   Vypocet polohovych vektorov na zaklade 3 pozorovani v danych casoch (albo polohach pozorovatela),
        *   Montenbruckova metoda, citaj 'Satellite orbits' , str. 43 - 46
	*   (Computes position vectors from 3 observations (or observers) and  associated times ), Montenbruck Method,
	*    read 'Satellite orbits', sites 43 - 46
	*
	* Input/Output:
	*
	*   Geodetic coordinates of observatory - array[3]
	*   Observation 1
	*   Observation 2
	*   Observation 3
	*
	*   <return>  Position vectors [m] for given times, or observer positions 3
	*
	* Notes:
	*
	*   The function cannot be used with state vectors describing a circular
	*   or non-inclined orbit.
	*
	*/

	public static Vector[] getPosVectors(Geodetic geodetic[], ObservationHor obs1, ObservationHor obs2, ObservationHor obs3){
		//konstanty
		Constants c = new Constants();
		Transformation transf = new Transformation();
		double minimum = 1e-15;

		//koeficienty n1, n3
		double n1, n3;

		//Polohovy vektor pozorovatela
		//Vector sta=new Vector(3);

		//Geodeticke suradnice - TREBA!!! zadat zemepisnu dlzku, sirku a nadm. vysku observatoria
		Geodetic geodObs[] = geodetic;//new Geodetic(Math.toRadians(18.6107),Math.toRadians(53.0217),0);

		//prevod suradnic pozorovatela z geodetickych do geocentrickych suradnic
		Vector obsGeocPosition[] = new Vector[3];
                for(int i=0;i<3;i++){
                    obsGeocPosition[i] = new Vector(3);
                    obsGeocPosition[i] = transf.fromGeodToGeoc(geodObs[i].lon, geodObs[i].lat, geodObs[i].altitude, c.R_Earth, c.f_Earth);
                }
		//v kadom case ma pozorovatel ine ekvatorialne suradnice (os x smeruje do jarneho bodu)
		Vector obsPosition[] = new Vector[3];
		obsPosition[0] = new Vector(3);
		obsPosition[1] = new Vector(3);
		obsPosition[2] = new Vector(3);

		//POZOROVANIA - CAS, AZIMUT, VYSKA
		//horizontalne suradnice
		double azi[] = new double[3];
		double elev[] = new double[3];
		Time time[] = new Time[3];

		//1. Pozorovanie
		time[0] = obs1.t;//new Time(0,2007,10,20,5,35,47.0);
		azi[0] = obs1.azim;//Math.toRadians(295.7);
		elev[0] = obs1.elev;//Math.toRadians(63.8);

		//2. Pozorovanie
		time[1] = obs2.t;//new Time(0,2007,10,20,6,35,47.0);
		azi[1] = obs2.azim;//Math.toRadians(194.6);
		elev[1] = obs2.elev;//Math.toRadians(80);

		//3. Pozorovanie
		time[2] = obs3.t;//new Time(0,2007,10,20,7,35,47.0);
		azi[2] = obs3.azim;//Math.toRadians(166.2);
		elev[2] = obs3.elev;//Math.toRadians(51.1);

		//Matica na prevod do geocentrickych suradnic
		Matrix ltcMat[] = new Matrix[3];
		for(int i=0;i<3;i++){
                    ltcMat[i] = new Matrix(3,3);
                    ltcMat[i] = transf.getLtcMatrix(geodObs[i].lon,geodObs[i].lat);
                }

		//Matica na prevod z ekvatorialnych do geocentrickych suradnic
		Matrix[] u = new Matrix[3];
		//vektor, ktory je prevod z horizontalnych do tangencialnych suradnic
		Vector[] s = new Vector[3];

		for(int i=0;i<3;i++){
			//inicializacia translacnej matice z ekvat. => geocentrickych, pouzije sa jej transp. matica
			u[i] = new Matrix(3,3);
			u[i] = u[i].R_z(time[i].getGMST(time[i].getMjd(time[i])));

			//pozicny vektor pozorovatela v ekvatorilanych suradniciach v danom case pozorovanie
			obsPosition[i] = u[i].matrixMultiplyVector(u[i].transposeMatrix(u[i]),obsGeocPosition[i]);
		}

		//jednotkove vektory pozorovani
		Vector e[] = new Vector[3];
		for(int i=0; i<3;i++){
			s[i] = new Vector(3);
			s[i] = s[i].fromPolarToCartesian(c.pi/2-azi[i],elev[i],1);
			e[i] = ltcMat[i].matrixMultiplyVector(u[i].transposeMatrix(u[i]),u[i].matrixMultiplyVector(ltcMat[i].transposeMatrix(ltcMat[i]),s[i]));
		}

		//d koeficienty pre vypocet vzdialenosti sat-poz
		double D,d_11,d_12,d_13,d_21,d_22,d_23,d_31,d_32,d_33;
		Vector d[] = new Vector[3];
		d[0] = transf.getVector_di(e[1],e[2]);
		d[1] = transf.getVector_di(e[2],e[0]);
		d[2] = transf.getVector_di(e[0],e[1]);

		D = transf.getCoefficient_D(e[0],e[1],e[2]);
		d_11 = transf.getCoefficient_Dij(d[0],obsPosition[0]);
		d_12 = transf.getCoefficient_Dij(d[0],obsPosition[1]);
		d_13 = transf.getCoefficient_Dij(d[0],obsPosition[2]);
		d_21 = transf.getCoefficient_Dij(d[1],obsPosition[0]);
		d_22 = transf.getCoefficient_Dij(d[1],obsPosition[1]);
		d_23 = transf.getCoefficient_Dij(d[1],obsPosition[2]);
		d_31 = transf.getCoefficient_Dij(d[2],obsPosition[0]);
		d_32 = transf.getCoefficient_Dij(d[2],obsPosition[1]);
		d_33 = transf.getCoefficient_Dij(d[2],obsPosition[2]);

		//koeficienty ro - vzdialenosti poz. - sat.
		double ro[] = new double[3];

		//ekvatorialne suradnice satelitov v danom case
		Vector[] satPosition = new Vector[3];

		//deklaracia triedy kepler
		Kepler kepler = new Kepler(0,0,0,0,0,0);

		//pomocne premenne
		int maxit = 0;
		double n1_help;
		double eta_1,eta_2,eta_3;
		double tau_1,tau_2,tau_3;

		//pociatocna podmienka pre iteraciu
		n1 = (time[2].getMjd(time[2])-time[1].getMjd(time[1]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));
		n1_help = n1;
		n3 = (time[1].getMjd(time[1])-time[0].getMjd(time[0]))/(time[2].getMjd(time[2])-time[0].getMjd(time[0]));

		//vypocet koeficientov tau_1, tau_2, tau_3
		tau_1 = kepler.getTau(c.GM_Earth,time[1].getMjd(time[1]),time[2].getMjd(time[2]));
		tau_2 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[2].getMjd(time[2]));
		tau_3 = kepler.getTau(c.GM_Earth,time[0].getMjd(time[0]),time[1].getMjd(time[1]));

		//vypocet vzdialenosti pre dane poc. podmienky
		ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
		ro[1] = (1/D)*(n1*d_21-d_22+n3*d_23);
		ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);

		//pole aktualnych vektorov - jednotkovy * vzdialenost
		Vector actualPos[] = new Vector[3];

		//vypocet polohovych vektorov
		for(int i = 0;i < 3;i++){
			actualPos[i] = new Vector(3);
			actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
			satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
		}

		//pole double pre kontrolu
		double itRo[][] = new double[3][3];


		while(true){
			//vypocet koeficientov eta_1,eta_2,eta_3, vsetko o nich je vo vzahoch 2.98-2.108,
			eta_1 = kepler.getEta(satPosition[1],satPosition[2],tau_1);
			eta_2 = kepler.getEta(satPosition[0],satPosition[2],tau_2);
			eta_3 = kepler.getEta(satPosition[0],satPosition[1],tau_3);

			//vypocet koeficientov n1 a n3
			n1 = transf.getCoefficient_n1(eta_1,eta_2,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));
			n3 = transf.getCoefficient_n3(eta_2,eta_3,time[0].getMjd(time[0]),time[1].getMjd(time[1]),time[2].getMjd(time[2]));

			//vypocet vzdialenosti sat. - poz
			ro[0] = -1/(n1*D)*(n1*d_11-d_12+n3*d_13);
			ro[1] = 1/D*(n1*d_21-d_22+n3*d_23);
			ro[2] = -1/(n3*D)*(n1*d_31-d_32+n3*d_33);

			//vypocet polohovych vektorov
			for(int i = 0;i<3;i++){
				actualPos[i] = actualPos[i].multiplyVector(e[i],ro[i]);
				satPosition[i] = satPosition[i].addVectors(obsPosition[i], actualPos[i]);
			}

			//ukladame ziskane ro cez iteracie, na kontrolu
			itRo[0][0] = itRo[1][0];
			itRo[0][1] = itRo[1][1];
			itRo[0][2] = itRo[1][2];

			itRo[1][0] = itRo[2][0];
			itRo[1][1] = itRo[2][1];
			itRo[1][2] = itRo[2][2];

			itRo[2][0] = ro[0];
			itRo[2][1] = ro[1];
			itRo[2][2] = ro[2];

			//test iteracie
			if(Math.abs(n1-n1_help)<minimum) {
				//System.out.println("No. of it: "+maxit+"\n");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i = 0;i<3;i++){
						for(int j = 0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j = 0;j<3;j++){
						//System.out.println(("ro"+j)+ro[j]);
					}
				}

				if((ro[0]<0)||(ro[1]<0)||(ro[2]<0)) {
                                    //System.out.println("Negative distance!!!");
                                    //System.out.println(" TEST Negative distance!!!");
                                    //Vector satPositionVecArray[] = new Vector[3];
                                    //return satPositionVecArray;
                                };
				break;
			}

			if(maxit > 1000) {
				//System.out.println("WARNING: Convergence problems in ro!!!");
				//vypis vysledkov poslednych 3.iteracii, pre kontrolu
				if(maxit>3){
					for(int i=0;i<3;i++){
						for(int j=0;j<3;j++){
							//System.out.println(("ro"+j) + ": " + itRo[i][j]);
						}
						//System.out.println("\n");
					}
				}
				else{
					for(int j=0;j<3;j++){
						//System.out.println(ro[j]);
					}
				}
                                Vector satPositionVecArray[] = new Vector[3];
                                satPositionVecArray[0] = new Vector(3);
                                satPositionVecArray[1] = new Vector(3);
                                satPositionVecArray[2] = new Vector(3);
                                return satPositionVecArray;

				//break;
			}

			n1_help=n1;
			maxit++;
		}

                //WARNING
                if((ro[0] < 0)||(ro[1] < 0)||(ro[2] < 0)) {
                    //System.out.println("WARNING: negative ro!!!");
                    //if ro is negative do not compute elements
                    Vector satPositionVecArray[] = new Vector[3];
                    satPositionVecArray[0] = new Vector(3);
                    satPositionVecArray[1] = new Vector(3);
                    satPositionVecArray[2] = new Vector(3);
                    return satPositionVecArray;
                }
                else {
                    //if (ro[0]>10000) System.out.println("Distance OK!!!");
                }
		//uz mame 3 pozicne vektory, treba vyratat elementy
		Constants consta = new Constants();
		//kepler = kepler.getElements(consta.GM_Earth,time[1].getMjd(time[1]),time[2].getMjd(time[2]),satPosition[1],satPosition[2]);
		
                //System.out.println("double x1 = " + satPosition[0].v[0] + ";");
		//System.out.println("double y1 = " + satPosition[0].v[1] + ";");
		//System.out.println("double z1 = " + satPosition[0].v[2] + ";");
                //System.out.println("double x2 = " + satPosition[1].v[0] + ";");
		//System.out.println("double y2 = " + satPosition[1].v[1] + ";");
		//System.out.println("double z2 = " + satPosition[1].v[2] + ";");
		//System.out.println("double x3 = " + satPosition[2].v[0] + ";");
                //System.out.println("double y3 = " + satPosition[2].v[1] + ";");
		//System.out.println("double z3 = " + satPosition[2].v[2] + ";");

                //System.out.println("double mjdTime_1 = " + Time.getMjd(time[0]) + ";");
                //System.out.println("double mjdTime_2 = " + Time.getMjd(time[1]) + ";");
                //System.out.println("double mjdTime_3 = " + Time.getMjd(time[2]) + ";");
                Transformation.getEclipticInclAndNodeFromGeocPosVectors(satPosition[0],satPosition[2],Time.getMjd(time[0]),Time.getMjd(time[2]));
		//return kepler;
                Vector satPositionVecArray[] = new Vector[3];
                satPositionVecArray[0] = satPosition[0];
                satPositionVecArray[1] = satPosition[1];
                satPositionVecArray[2] = satPosition[2];
                return satPositionVecArray;
	}

        /**
         * getVectors()
	*   Vypocet polohovych vektorov na zaklade 3 pozorovani v danych casoch (albo polohach pozorovatela), vypocet ventora rychlosti.
        *   Montenbruckova metoda, citaj 'Satellite orbits' , str. 43 - 46
	*   (Computes position vectors from 3 observations (or observers) and  associated times, velocity vector fo 2nd time ), Montenbruck Method,
	*    read 'Satellite orbits', sites 43 - 46
	*
	* Input/Output:
	*   Observation observation[] - arrray of Observations - see constructor Observation
	*
	*   OUTPUT:
        *  Vector vector[4] - Array of vectors:
        *               - 1st position vector for 1st time [m]
        *               - 2nd position vector for 2nd time [m]
        *               - 3rd position vector for 3rd time [m]
        *               - velocity Vector for 2nd time     [m/s]
	*
	*/

	public static Vector[] getVectors(Observation observation[]){
            //output vectors - 3 position and 1 velocity vectors
            Vector vectors[] = new Vector[3];
            
            //transformation from equatroail to horizontal coordinates
            //1st observation
            ObservationHor obsHor[] = new ObservationHor[3];
            for(int i = 0; i < 3; i++){
                obsHor[i] = new ObservationHor();
                obsHor[i].t = Time.getDateTime(observation[i].timeMjd);
                //obsHor[i].azim = new Transformation().fromEquatToHorizCoord_New(observation[i].ra, observation[i].dec,
                //        observation[i].lon, observation[i].lat, observation[i].timeMjd).v[0];
                obsHor[i].azim = new Transformation().fromEquatToHorizontal(observation[i].ra, observation[i].dec,
                        Time.getGMST(observation[i].timeMjd), observation[i].lon, observation[i].lat).v[0];
                //obsHor[i].elev = new Transformation().fromEquatToHorizCoord_New(observation[i].ra, observation[i].dec,
                //        observation[i].lon, observation[i].lat, observation[i].timeMjd).v[1];
                obsHor[i].elev = new Transformation().fromEquatToHorizontal(observation[i].ra, observation[i].dec,
                        Time.getGMST(observation[i].timeMjd), observation[i].lon, observation[i].lat).v[1];
            }

            Geodetic geodetic[] = new Geodetic[3];
            for(int i = 0; i<3; i++){
                geodetic[i] = new Geodetic();
                geodetic[i].lon = observation[i].lon;
                geodetic[i].lat = observation[i].lat;
                geodetic[i].altitude = observation[i].alt;
            }
            
            vectors = getPosVectors(geodetic, obsHor[0], obsHor[1], obsHor[2]);

            Vector vectors2[] = new Vector[4];
            vectors2[0] = vectors[0];   //1st position vector
            vectors2[1] = vectors[1];   //2nd position vector
            vectors2[2] = vectors[2];   //3rd position vector

            //fill arrays for compute the velocity vector
            //position vectors - should be in Erath radius for Escobal method
            Vector vectors3[] = new Vector[3];
            //times for given position vectors - time should be Julian dates [day] for Escobal method
            double time3[] = new double[3];
            for(int i = 0; i<3; i++){
                //vectors3[i] = vectors2[i];
                time3[i] = Time.getJdFromMjd(Time.getMjd(obsHor[i].t));
            }
            //position vector should be in Earth radius, convert them from [m]
            for(int i=0; i<3; i++){
                //for(int j=0; j<3; j++){
                    //[m] to Earth radius
                vectors3[i] = new Vector(3);
                vectors3[i] = vectors2[i];
                //}
                vectors3[i] = Vector.multiplyVector(vectors3[i], (double)1/Constants.R_Earth);
            }

            //constants for geocentric motion - Escobal units
            double k2 = 0.07436574;  //[(e.r)^3/2 / min]
            double mi2 = 1.0;        //[e.m.] - Earth mass
            
            vectors2[3] = new RIterationAnglesOnly().getVelVecFrom3PosVec(vectors3, time3, k2, mi2, false);
            //get velocity from Escobal units to [m/s]
            //for(int i =0; i<3; i++){
                  //vectors2[3].v[i] = vectors2[3].v[i]*25936*0.3048;
                  vectors2[3] = Vector.multiplyVector(vectors2[3],25936*0.3048);
                  //System.out.println("velocity_" + Vector.getSize(vectors2[3]));
            //}
            //System.out.println("geocentric vel " + Vector.getSize(vectors2[3]) + " m/s");
            
            //back to SI units
            //for(int i=0; i<3; i++){
            //    for(int j = 0; j < 3; j++){
            //        //Earth radius to [m]
            //        vectors3[i].v[j] = vectors3[i].v[j]*Constants.R_Earth;
            //    }
            //}

            //System.out.println("M vec1.x " + vectors2[0].v[0]);
            //System.out.println("M vec1.y " + vectors2[0].v[1]);
            //System.out.println("M vec1.z " + vectors2[0].v[2]);

            return vectors2;
        }


	/*
	TEST METHOD
	*/
	
	public static void main(String []args){

            char degChar = 176;
            String stringDeg = degChar + "";

            //constants
            Constants c = new Constants();

            //observatory position - AGO Modra
            Geodetic geodObs = new Geodetic(Math.toRadians(17.2740),Math.toRadians(48.3733),531.1);
                        
            //positions
            int amount = 17;
            double azi[]=new double[amount];
            double elev[]=new double[amount];
            Time time[]=new Time[amount];

            double meanAnom[] = new double[amount];

            ObservationHor obser[] = new ObservationHor[amount];

            Kepler kepler = new Kepler();

            //combinations of observations
            int a1,a2,a3;

            /*
            //Tv camera observation from 17/06/2009
            //Observation of object 33777 - Iridium 33 deb
            //1 33777U 97051Q   09168.71409626  .00000428  00000-0  14492-3 0   952
            //2 33777 086.3874 068.0339 0008045 289.9426 070.0880 14.34550544 18143
            
            //1st pozition
            //time[0] = new Time(2009, 06, 17,21, 6, 0.72);
            //azi[0] = Math.toRadians(97.0627);
            //elev[0] = Math.toRadians(66.7976);
            //meanAnom[0] = getAngle(202.618);    //mean anomaly from SatEph

            //2nd pozition
            //time[1] = new Time(2009, 06, 17,21, 6, 17.48);
            //azi[1] = Math.toRadians(115.2003);
            //elev[1] = Math.toRadians(64.1480);
            //meanAnom[1] = getAngle(203.62);    //mean anomaly from SatEph

            //3rd pozition
            //time[2] = new Time(2009, 06, 17,21, 6, 25.86);
            //azi[2] = Math.toRadians(122.6298);
            //elev[2] = Math.toRadians(62.0616);
            //meanAnom[2] = getAngle(204.12);    //mean anomaly from SatEph

            //4th pozition
            //time[3] = new Time(2009, 06, 17,21, 6, 34.26);
            //azi[3] = Math.toRadians(128.8211);
            //elev[3] = Math.toRadians(59.6636);
            //meanAnom[3] = getAngle(204.623);    //mean anomaly from SatEph

            //5th pozition
            //time[4] = new Time(2009, 06, 17,21, 6, 59.5);
            //azi[4] = Math.toRadians(141.7862);
            //elev[4] = Math.toRadians(51.7536);
            //meanAnom[4] = getAngle(206.131);    //mean anomaly from SatEph

            //1st pozition
            time[0] = new Time(2009, 06, 17,21, 6, 0.76);
            azi[0] = Math.toRadians(96.9753);
            elev[0] = Math.toRadians(66.8813);
            meanAnom[0] = getAngle(202.618);    //mean anomaly from SatEph

            //2nd pozition
            time[1] = new Time(2009, 06, 17,21, 6, 25.91);
            azi[1] = Math.toRadians(122.882);
            elev[1] = Math.toRadians(61.9934);
            meanAnom[1] = getAngle(203.62);    //mean anomaly from SatEph

            //3rd pozition
            time[2] = new Time(2009, 06, 17,21, 6, 59.5);
            azi[2] = Math.toRadians(141.8267);
            elev[2] = Math.toRadians(51.7589);
            meanAnom[2] = getAngle(204.12);    //mean anomaly from SatEph
            
            //fill array of observations
            for(int i = 0; i < amount; i++){
                obser[i] = new ObservationHor(azi[i], elev[i], time[i]);
            }

            a1 = 0;
            a2 = 1;
            a3 = 2;

            kepler = getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);

            System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

            System.out.println("Computed a: " + kepler.a/1000/6378.15 + " e.r.");
            System.out.println("Computed a: "+kepler.a/1000+" km");
            System.out.println("Computed e: "+kepler.e);
            System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+ " " + stringDeg);
            System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+ " " + stringDeg);
            System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+ " " + stringDeg);
            System.out.println("cOMPUTED Ma   : "+Math.toDegrees(kepler.M)+ " " + stringDeg + "\n");
            */
            /*
            //Tv camera observation from 16/09/2009
            //Observation of object 25723 - 
            //1 25723U 99022C   09256.64718178 -.00000716  00000-0 -20684-4 0  3832
            //2 25723 048.4407 324.3833 0021806 174.1355 185.9833 15.16115652572720

            Geodetic geodObs_Array[] = new Geodetic[3];
            //ago
            geodObs_Array[0] = new Geodetic(Math.toRadians(17.2740),Math.toRadians(48.3733),531.1);
            geodObs_Array[1] = new Geodetic(Math.toRadians(17.2740),Math.toRadians(48.3733),531.1);
            //arbo
            geodObs_Array[2] = new Geodetic(Math.toRadians(18.3685), Math.toRadians(48.3235),185.0);
           
            //1st pozition - AGO
            time[0] = new Time(2009,9,16,2,28,17.446);
            azi[0] = Math.toRadians(273.9952);
            elev[0] = Math.toRadians(87.7871);
            meanAnom[0] = getAngle(269.762);    //mean anomaly from SatEph

            //2nd pozition - AGO
            time[1] = new Time(2009,9,16,2,28,21.369);
            azi[1] = Math.toRadians(97.3703);
            elev[1] = Math.toRadians(89.0732);
            meanAnom[1] = getAngle(270.199);    //mean anomaly from SatEph

            //3rd pozition - AGO
            time[2] = new Time(2009,9,16,2,28,24.370);
            azi[2] = Math.toRadians(99.7246);
            elev[2] = Math.toRadians(86.5353);
            meanAnom[2] = getAngle(270.199);    //mean anomaly from SatEph

            //1st position - Arbo
            //time[2] = new Time(2009,9,16,2,28,17.058);
            //azi[2] = Math.toRadians(274.8248);
            //elev[2] = Math.toRadians(78.2012);
            //meanAnom[2] = 269.737;    //mean anomaly from SatEph

            //2nd position Arbo
            time[3] = new Time(2009,9,16,2,28,23.699);
            azi[3] = Math.toRadians(271.6158);
            elev[3] = Math.toRadians(83.3391);
            meanAnom[3] = getAngle(270.157);    //mean anomaly from SatEph

            
            //SatEph positions
            //1st pozition - AGO - Sateph
            //time[0] = new Time(2009,9,16,2,28,17.446);
            //azi[0] = Math.toRadians(268.9854);
            //elev[0] = Math.toRadians(87.70859);
            //meanAnom[0] = getAngle(269.762);    //mean anomaly from SatEph

            //2nd pozition - AGO - Sateph
            //time[1] = new Time(2009,9,16,2,28,24.370);
            //azi[1] = Math.toRadians(100.71089);
            //elev[1] = Math.toRadians(86.73897);
            //meanAnom[1] = getAngle(270.199);    //mean anomaly from SatEph
            
            //1st position - Arbo - Sateph
            //time[2] = new Time(2009,9,16,2,28,17.058);
            //azi[2] = Math.toRadians(273.53066);
            //elev[2] = Math.toRadians(77.92056);
            //meanAnom[2] = 269.737;    //mean anomaly from SatEph
            
            //2nd position Arbo - Sateph
            //time[3] = new Time(2009,9,16,2,28,23.699);
            //azi[3] = Math.toRadians(271.1246);
            //elev[3] = Math.toRadians(83.07175);
            //meanAnom[3] = getAngle(270.157);    //mean anomaly from SatEph
            
            //fill array of observations
            for(int i = 0; i < amount; i++){
                obser[i] = new ObservationHor(azi[i], elev[i], time[i]);
            }

            a1 = 0;
            a2 = 1;
            a3 = 3;

            kepler = getElementsMM_2(geodObs_Array, obser[a1], obser[a2], obser[a3]);

            System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

            System.out.println("Computed a: " + kepler.a/1000/6378.15 + " e.r.");
            System.out.println("Computed a: "+kepler.a/1000+" km");
            System.out.println("Computed e: "+kepler.e);
            System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+ " " + stringDeg);
            System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+ " " + stringDeg);
            System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+ " " + stringDeg);
            System.out.println("cOMPUTED Ma   : "+Math.toDegrees(kepler.M)+ " " + stringDeg + "\n");
            */

            //Tv camera observation from 17/06/2009
            //Observation of object 33777 - Iridium 33 deb
            //1 33777U 97051Q   09168.71409626  .00000428  00000-0  14492-3 0   952
            //2 33777 086.3874 068.0339 0008045 289.9426 070.0880 14.34550544 18143

            //1st pozition
            //time[0] = new Time(2009, 06, 17,21, 6, 0.72);
            //azi[0] = Math.toRadians(97.0627);
            //elev[0] = Math.toRadians(66.7976);
            //meanAnom[0] = getAngle(202.618);    //mean anomaly from SatEph

            //2nd pozition
            //time[1] = new Time(2009, 06, 17,21, 6, 17.48);
            //azi[1] = Math.toRadians(115.2003);
            //elev[1] = Math.toRadians(64.1480);
            //meanAnom[1] = getAngle(203.62);    //mean anomaly from SatEph

            //3rd pozition
            //time[2] = new Time(2009, 06, 17,21, 6, 25.86);
            //azi[2] = Math.toRadians(122.6298);
            //elev[2] = Math.toRadians(62.0616);
            //meanAnom[2] = getAngle(204.12);    //mean anomaly from SatEph

            /**
             * Test from 09/06/2009 - Escobal vs Bucerius/Montenbruck
             */

            //Excercise 6, s.291
            //San Fernando, Spain
            /*
            geodObs = new Geodetic(Math.toRadians(353.79486111),Math.toRadians(36.4638333),24);
            
            //array of geodetic position
            Geodetic geodObsArray[] = new Geodetic[3];
            for(int i=0; i<3; i++){
                geodObsArray[i] = new Geodetic();
                geodObsArray[i] = geodObs;
            }

            time[0] = new Time(1959, 9, 26, 21, 25, 37.403);
            Vector eqToHor = Transformation.fromEquatToHorizontal(Math.toRadians(254.673208333), Math.toRadians(13.0927222), Time.getGMST(Time.getMjd(time[0])), geodObs.lon, geodObs.lat);
            //Vector eqToHor = Transformation.fromEquatToHorizontal(254.673, 13.093, Time.getGMST(Time.getMjd(time[0])), geodObs.lon, geodObs.lat);
            azi[0] = eqToHor.v[0];
            elev[0] = eqToHor.v[1];
            meanAnom[0] = getAngle(0);

            time[1] = new Time(1959, 9, 26, 21, 26, 37.862);
            eqToHor = Transformation.fromEquatToHorizontal(Math.toRadians(260.941125), Math.toRadians(13.28569444), Time.getGMST(Time.getMjd(time[1])), geodObs.lon, geodObs.lat);
            //eqToHor = Transformation.fromEquatToHorizontal(260.94, 13.286, Time.getGMST(Time.getMjd(time[1])), geodObs.lon, geodObs.lat);
            azi[1] = eqToHor.v[0];
            elev[1] = eqToHor.v[1];
            meanAnom[1] = getAngle(0);

            time[2] = new Time(1959, 9, 26, 21, 27, 45.919);
            eqToHor = Transformation.fromEquatToHorizontal(Math.toRadians(269.85375), Math.toRadians(13.00425), Time.getGMST(Time.getMjd(time[2])), geodObs.lon, geodObs.lat);
            //eqToHor = Transformation.fromEquatToHorizontal(269.854, 13.004, Time.getGMST(Time.getMjd(time[2])), geodObs.lon, geodObs.lat);
            azi[2] = eqToHor.v[0];
            elev[2] = eqToHor.v[1];
            meanAnom[2] = getAngle(0);

            for(int i = 0; i < amount; i++){
                obser[i] = new ObservationHor(azi[i], elev[i], time[i]);
            }

            a1 = 0;
            a2 = 1;
            a3 = 2;

            //kepler = getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);
            Vector resultVec[] = new Vector[3];
            resultVec = getPosVectors(geodObsArray, obser[a1], obser[a2], obser[a3]);
            //position vectors
            kepler = kepler.getElements(Constants.GM_Earth,time[0].getMjd(time[0]),time[2].getMjd(time[2]),resultVec[0],resultVec[2]);

            System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

            System.out.println("Escobal  a: 9184.536 km");
            System.out.println("Computed a: "+kepler.a/1000+" km");
            System.out.println("Escobal  e: 0.23167");
            System.out.println("Computed e: "+kepler.e);
            System.out.println("Escobal  i: 33.281 deg");
            System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" deg");
            System.out.println("Escobal  Omega: 205.110deg");
            System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" deg");
            System.out.println("Escobal  omega: 161.790�");
            System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" deg");
            System.out.println("Ma        : ???deg");
            System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" deg\n");
            */
            /*
            //Escobal ref orbit no.VII. / not done yet
            geodObs=new Geodetic(Math.toRadians(250.0),Math.toRadians(40.0),5000);
            
            time[0] = new Time(1959, 9, 26, 21, 25, 37.403);
            Vector eqToHor = Transformation.fromEquatToHorizontal(254.673208333, 13.0927222, Time.getGMST(2438314.7916667 - 2400000.5), geodObs.lon, geodObs.lat);
            //Vector eqToHor = Transformation.fromEquatToHorizontal(254.673, 13.093, Time.getGMST(Time.getMjd(time[0])), geodObs.lon, geodObs.lat);
            azi[0] = eqToHor.v[0];
            elev[0] = eqToHor.v[1];
            meanAnom[0] = getAngle(0);

            time[1] = new Time(1959, 9, 26, 21, 26, 37.862);
            eqToHor = Transformation.fromEquatToHorizontal(260.941125, 13.28569444, Time.getGMST(2438314.8055556 - 2400000.5), geodObs.lon, geodObs.lat);
            //eqToHor = Transformation.fromEquatToHorizontal(260.94, 13.286, Time.getGMST(Time.getMjd(time[1])), geodObs.lon, geodObs.lat);
            azi[1] = eqToHor.v[0];
            elev[1] = eqToHor.v[1];
            meanAnom[1] = getAngle(0);

            time[2] = new Time(1959, 9, 26, 21, 27, 45.919);
            eqToHor = Transformation.fromEquatToHorizontal(269.85375, 13.00425, Time.getGMST(2438314.8194444 - 2400000.5), geodObs.lon, geodObs.lat);
            //eqToHor = Transformation.fromEquatToHorizontal(269.854, 13.004, Time.getGMST(Time.getMjd(time[2])), geodObs.lon, geodObs.lat);
            azi[2] = eqToHor.v[0];
            elev[2] = eqToHor.v[1];
            meanAnom[2] = getAngle(0);

            for(int i = 0; i < amount; i++){
                obser[i] = new Observation(azi[i], elev[i], time[i]);
            }

            a1 = 0;
            a2 = 1;
            a3 = 2;

            kepler = getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);

            System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

            System.out.println("TLE      a: 9184.536 km");
            System.out.println("Computed a: "+kepler.a/1000+" km");
            System.out.println("TLE      e: 0.23167");
            System.out.println("Computed e: "+kepler.e);
            System.out.println("TLE      i: 33.281�");
            System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" �");
            System.out.println("TLE      Omega: 205.110�");
            System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" �");
            System.out.println("TLE      omega: 161.790�");
            System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" �");
            System.out.println("TLE/Kepler: ???�");
            System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" �\n");
            */

            
            //Tv camera observation Reentry 2009/07/09

            //1st pozition
            time[0] = new Time(2010, 7, 9,23, 1, 58.119999);
            azi[0] = Math.toRadians(165.2448885);
            elev[0] = Math.toRadians(45.9058595);
            meanAnom[0] = getAngle(202.618);    //mean anomaly from SatEph

            //2nd pozition
            time[1] = new Time(2010, 7, 9,23, 1, 58.439999);
            azi[1] = Math.toRadians(161.2335997);
            elev[1] = Math.toRadians(45.1443069);
            meanAnom[1] = getAngle(203.62);    //mean anomaly from SatEph

            //3rd pozition
            time[2] = new Time(2010, 7, 9,23, 1, 58.839999);
            azi[2] = Math.toRadians(157.6121359);
            elev[2] = Math.toRadians(44.2116381);
            meanAnom[2] = getAngle(204.12);    //mean anomaly from SatEph

            //4th pozition
            time[3] = new Time(2010, 7, 9,23, 1, 59.239999);
            azi[3] = Math.toRadians(154.1043277);
            elev[3] = Math.toRadians(43.1677772);
            meanAnom[3] = getAngle(204.623);    //mean anomaly from SatEph

            //5th pozition
            time[4] = new Time(2010, 7, 9,23, 1, 59.639999);
            azi[4] = Math.toRadians(150.5520087);
            elev[4] = Math.toRadians(42.0166291);
            meanAnom[4] = getAngle(206.131);    //mean anomaly from SatEph

            //6th pozition
            time[5] = new Time(2010, 7, 9,23, 2, 0.039999);
            azi[5] = Math.toRadians(147.0251502);
            elev[5] = Math.toRadians(40.6803367);
            meanAnom[5] = getAngle(202.618);    //mean anomaly from SatEph

            //7th pozition
            time[6] = new Time(2010, 7, 9,23, 2, 0.439999);
            azi[6] = Math.toRadians(143.61804);
            elev[6] = Math.toRadians(39.2443045);
            meanAnom[6] = getAngle(203.62);    //mean anomaly from SatEph

            //8th pozition
            time[7] = new Time(2010, 7, 9,23, 2, 0.839999);
            azi[7] = Math.toRadians(140.4049605);
            elev[7] = Math.toRadians(37.7205573);
            meanAnom[7] = getAngle(204.12);    //mean anomaly from SatEph

            //9th pozition
            time[8] = new Time(2010, 7, 9,23, 2, 1.239999);
            azi[8] = Math.toRadians(137.2629554);
            elev[8] = Math.toRadians(36.0643443);
            meanAnom[8] = getAngle(204.12);

            //10th pozition
            time[9] = new Time(2010, 7, 9,23, 2, 1.639999);
            azi[9] = Math.toRadians(134.449988);
            elev[9] = Math.toRadians(34.454899);
            meanAnom[9] = getAngle(204.12);

            //11th pozition
            time[10] = new Time(2010, 7, 9,23, 2, 2.039999);
            azi[10] = Math.toRadians(131.7474735);
            elev[10] = Math.toRadians(32.8350332);
            meanAnom[10] = getAngle(204.12);

            //12th pozition
            time[11] = new Time(2010, 7, 9,23, 2, 2.439999);
            azi[11] = Math.toRadians(129.2471499);
            elev[11] = Math.toRadians(31.2334717);
            meanAnom[11] = getAngle(204.12);

            //13th pozition
            time[12] = new Time(2010, 7, 9,23, 2, 2.839999);
            azi[12] = Math.toRadians(126.9674203);
            elev[12] = Math.toRadians(29.6441429);
            meanAnom[12] = getAngle(204.12);

            //14th pozition
            time[13] = new Time(2010, 7, 9,23, 2, 3.239999);
            azi[13] = Math.toRadians(124.8676988);
            elev[13] = Math.toRadians(28.0650011);
            meanAnom[13] = getAngle(204.12);

            //15th pozition
            time[14] = new Time(2010, 7, 9,23, 2, 3.639999);
            azi[14] = Math.toRadians(122.9708784);
            elev[14] = Math.toRadians(26.6128552);
            meanAnom[14] = getAngle(204.12);
            
            //16th pozition
            time[15] = new Time(2010, 7, 9,23, 2, 4.039999);
            azi[15] = Math.toRadians(121.240836);
            elev[15] = Math.toRadians(25.2767419);
            meanAnom[15] = getAngle(204.12);

            //17th pozition
            time[16] = new Time(2010, 7, 9,23, 2, 4.319999);
            azi[16] = Math.toRadians(119.8038891);
            elev[16] = Math.toRadians(24.0455916);
            meanAnom[16] = getAngle(204.12);



            //fill array of observations
            for(int i = 0; i < amount; i++){
                obser[i] = new ObservationHor(azi[i], elev[i], time[i]);
            }

            a1 = 2;
            a2 = 6;
            a3 = 9;

            kepler = getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);

            System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

            System.out.println("Computed a: " + kepler.a/1000/6378.15 + " e.r.");
            System.out.println("Computed a: "+kepler.a/1000+" km");
            System.out.println("Computed e: "+kepler.e);
            System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+ " " + stringDeg);
            System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+ " " + stringDeg);
            System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+ " " + stringDeg);
            System.out.println("cOMPUTED Ma   : "+Math.toDegrees(kepler.M)+ " " + stringDeg + "\n");

            //make the arrray of integers for method getCombinations()
            int intArray[] = new int[obser.length];
            for(int i = 0; i < intArray.length; i++){
                intArray[i] = i;
            }
            //making the array of Vectors, which are holding information about combinations
            // v[0] - 1st observation to use to OD
            // v[1] - 2nd observation to use to OD
            // v[2] - 3rd observation to use to OD
            Vector vectorArray[] = new CombinationGenerator().getCombinations(intArray, false);

            //array of observation for the method getVectors()
            Observation observation2[] = new Observation[3];
            Vector vectorsForOd[] = new Vector[4];
            //orbital elements
            Kepler getKepler = new Kepler();
            for(int i = 0; i < vectorArray.length; i++){
                a1 = (int)vectorArray[i].v[0];
                a2 = (int)vectorArray[i].v[1];
                a3 = (int)vectorArray[i].v[2];
                kepler = getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);
                if(kepler.a > 6378000) {
                    System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

                    System.out.println("Computed a: " + kepler.a/1000/6378.15 + " e.r.");
                    System.out.println("Computed a: "+kepler.a/1000+" km");
                    System.out.println("Computed e: "+kepler.e);
                    System.out.println("Computed q: "+(kepler.a*(1-kepler.e))/1000 + " km");
                    System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+ " " + stringDeg);
                    System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+ " " + stringDeg);
                    System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+ " " + stringDeg);
                    System.out.println("cOMPUTED Ma   : "+Math.toDegrees(kepler.M)+ " " + stringDeg + "\n");
                }
            }
        }
	
        //metoda na zjednodusenie vypisocvania textu
	static void printString(String type, String value, String unit){
		System.out.println(type + ": " + value + " " + unit);
	}
        
        /*
         *Method to get the angle 
         *IN: [deg]
         *
         *OUT: [deg]      
         */
	public static double getAngle(double angle){
            while(angle > 360) angle = angle - 360;
            while(angle < 0) angle = angle + 360;
            
            return angle;
        }
}