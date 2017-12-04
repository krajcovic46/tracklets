/**
 * Class is for computing the geocentric position vectors of body with only
 * 3 topocentric positions on the celestial sphere - Az_i,h_i, or R.A._i, dec_i, for i = 1,2,3.
 * Double r-iteration is used. All formulas are from book P.Escobal - Methods of orbit
 * determination (1976). All formulas numbers in this code are identical with
 * formulas numbers used in Escobal's book.
 *
 * Author of Java source is Jiri Silha - 2009
 */

package com.skrajcovic.orbitdetermination;

import com.skrajcovic.orbitdetermination.compute.*;

/**
 *
 * @author Jiri Silha - 20. April 2009
 */
public class RIterationAnglesOnly {    
    
    //Variables
    
    //Incomming variables - 
    
    /**
     * Array with right accessions [rad]
     */
    private double ra[] = new double[3];
    
    /**
     * Array with declinations [rad]
     */
    private double dec[] = new double[3];
    
    /**
     * Array with observation times
     */
    private double time[] = new double[3];
    
    /**
     * Array with observer's position 
     */
    private Geodetic geodetic[] = new Geodetic[3]; 
    
    /**
     * Earth rotation
     */
    private double dTheta_dTime;
    
    /**
     * Earth flattening
     */
    private double flattening;
    
    /**
     * Earth radius (equator)
     */
    private double a_e;
    
    /**
     * ???
     */
    private double mi;
    
    /**
     * ???
     */
    private double k;
    
    /**
     * Is motion retrodrade
     */
    private boolean isMotionRetrograde;
    
    /**
     * 1st approximations of position vector sizes
     */
    private double r_1Appr, r_2Appr;
    
    /**
     * Number of revolutions - Formula 7/258
     */
    private int lambda;

    /**
     * Formulas 7.274 - 7.275
     */
    private double tau_1, tau_3;

    /**
     * Forulas 7.276
     */
    private Vector l_Vec[] = new Vector[3];
    
    /**
     * Formula 7.277
     */
    private double g_1[] = new double[3];
    
    /**
     * Formula 7.277
     */
    private double g_2[] = new double[3];
    
     /**
     * Formula 7.7.278
     */
    static double theta[] = new double[3];
    
    /**
     * Formula 7.278
     */
    static double thetaG_2;
    
    /**
     * Formula 7.279
     */
    static Vector obsR_Vec[] = new Vector[3];
    
    /**
     * Formula 7.280
     */
    static double c_psi[] = new double[3];
    
    /**
     * Before Formula 7.281
     */
    static double r_1, r_2, r_3, r_1g, r_2g;
    
    /**
     * Formula 7.281
     */
    static double ro[] = new double[3];
    
    /**
     * Formula 7.282
     */
    static Vector r_Vec[] = new Vector[3];
    
    /**
     * Formula 7.283
     */
    static Vector w_Vec = new Vector(3);
    
    /**
     * Formulas 7.287 - 7.288
     */
    static double v[] = new double[3];
    static double cos_v2_v1, cos_v2_v2, cos_v3_v1, cos_v3_v2;
    
    /**
     * Formulas 7.288 - 7.289
     */
    static double sin_v2_v1, sin_v2_v2, sin_v3_v1, sin_v3_v2;
    
    /**
     * Formulas 7.290
     */
    static double c_1, c_3;
    
    /**
     * Formula 7.291
     */
    static double p;
    
    /**
     * Formula 7.292
     */
    static double cr_1, cr_3;
    
    /**
     * Formulas 7.294 - 7.296
     */
    static double e_cosv_1, e_cosv_2, e_cosv_3, e_sinv_1, e_sinv_2, e_sinv_3;
    
    /**
     * Formulas 7.297
     */
    static double e;
    
    /**
     * Formula 7.298
     */
    static double a;
    
    /**
     * Formula 7.299
     */
    static double n;
    
    /**
     * Formulas 7.300 - 7.301
     */
    static double s_e, c_e;
    
    /**
     * Formulas 7.302 - 7.303
     */
    static double sin_E3_E2, cos_E3_E2, sin_E2_E1, cos_E2_E1;
    
    /**
     * Formulas 7.304 - 7.305
     */
    static double m3_M2, m1_M2, e3_E2, e2_E1;
    
    /**
     * Formulas 7.305
     */
    static double f_1, f_2, f_3;
    
    /**
     * Formulas 7.306 - 7.307
     */
    static double s_h, c_h;
    
    /**
     * Formulas 7.312 - 7.313
     */
    static double dF_1_DIV_dr_1, dF_2_DIV_dr_1; 
    
    /**
     * Formulas 7.314 - 7.315
     */
    static double dF_1_DIV_dr_2, dF_2_DIV_dr_2;
    
    /**
     * Formulas 7.316 - 7.320
     */
    static double delta, delta_1, delta_2, deltar_1, deltar_2;
    
    /**
     * Formula 7.322
     */
    static double r_1_n, r_1_nPlus1,r_2_n, r_2_nPlus1;
    
    /**
     * Formulas 7.323 - 7.324
     */
    static double fSerie, gSerie;
    
    /**
     * Formulas 7.325
     */
    static Vector r_2_Dot = new Vector(3);
    
    
    //Global variables - increment of r, tolerance for delta_r, number of iteration allowed
    
    /**
     * Condition before formula 7.312
     */
    private double delta_r_PerCent = 0.000001;
    
    /**
     * Formula 7.321
     */
    private double epsilonTolerance = 1.0*10e-8;//a_e*10e-6;;
    
    /**
     * Allowed number of iteration
     */
    private int iterationBoundery = 50;
    
    
    /**
      * getStateVector()
      *  
      * Escobal's double r-iteration method to get 2nd geocentric position and velocity 
      * vectors of body from 3 angles (R.A., dec)
      * 
      * INPUT:
      *     double ra[i] - right accesions of body for time[i], i = 0,1,2
      *     double dec[i] - declinations of body for time[i], i = 0,1,2
      *     double time[i] - times of observation of body, Julian date [min] i = 0,1,2
      *     Geodetic geodetic[i] - geodetic positions of observer, i = 0,1,2
      *     double dTheta_dTime - Earth rotation
      *     double f - flattening of Earth
      *     double a_e - planet radius (equator)
      *     double mi - 
      *     double k  - 
      *     boolean isMotionRetrograde - is motion retrograde 
      *     double r_1Appr - 1st input approximation of distance r_1 [a_e] 
      *     double r_2Appr - 1st input approximation of distance r_2 [a_e] 
      *     int lambda - number of revolutions    
      * 
      * OUTPUT:
      *     StateVector stateVector - stat vector - positiona and velocity vectors 
      */
    public StateVector getStateVector(double ra[], double dec[], double time[], 
                Geodetic geodetic[], double dTheta_dTime, double flattening, double a_e,
                double mi, double k, boolean isMotionRetrograde, double r_1Appr, double r_2Appr,
                int lambda){
        
        this.k = k;
        this.mi = mi;
        this.ra = ra;
        this.dec = dec;
        this.time = time;
        this.geodetic = geodetic;
        this.dTheta_dTime = dTheta_dTime;
        this.flattening = flattening;
        this.a_e = a_e;
        this.isMotionRetrograde = isMotionRetrograde;
        this.r_1Appr = r_1Appr;
        this.r_2Appr = r_2Appr;
        this.lambda = lambda;
        
        //condition 7.321
        double epsilon = a_e*epsilonTolerance;//a_e*10e-6;
        
        //initialization of main variable - info about position and velocity
        StateVector sv = new StateVector();
        
        // Formula 7.274
        tau_1 = k * (time[0] - time[1]);
        
        // Formula 7.275
        tau_3 = k * (time[2] - time[1]);
        
        //System.out.println("tau_1 " + tau_1);
        //System.out.println("tau_3 " + tau_3);
        
        //Greenwich siderial time - 1.27 [s]
        thetaG_2 = getGreenwichSiderialTime(time[1]);
        //System.out.println("thetaG_2 " + Math.toDegrees(thetaG_2) + " deg");
        
        for(int i = 0; i < 3; i++){
            // Formulas 7.276         
            l_Vec[i] = new Vector(3);
            l_Vec[i].v[0] = Math.cos(dec[i])*Math.cos(ra[i]);
            l_Vec[i].v[1] = Math.cos(dec[i])*Math.sin(ra[i]);
            l_Vec[i].v[2] = Math.sin(dec[i]);
            
            //System.out.println("l_Vec[" + i + "] size " + Vector.getSize(l_Vec[i]));
            
            //Formulas 7.277
            g_1[i] = a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) + 
                    geodetic[i].altitude;
            g_2[i] = (1 - flattening)*(1 - flattening)*a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) + 
                    geodetic[i].altitude;
            
            //System.out.println("g_1[0] " + g_1[0]);
            //System.out.println("g_2[0] " + g_2[0]);
            
            //Formula 7.278
            theta[i] = thetaG_2 + dTheta_dTime*(time[i] - time[1]) + geodetic[i].lon;
            //System.out.println("theta[i] " + i + " " + theta[i]);
            
            //Formulas 7.279
            obsR_Vec[i] = new Vector(3);
            obsR_Vec[i].v[0] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.cos(theta[i]);
            obsR_Vec[i].v[1] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.sin(theta[i]);
            obsR_Vec[i].v[2] = -g_2[i] * Math.sin(geodetic[i].lat);
            
            //System.out.println("observer x " + obsR_Vec[i].v[0] + " e.r.");
            //System.out.println("observer y " + obsR_Vec[i].v[1] + " e.r.");
            //System.out.println("observer z " + obsR_Vec[i].v[2] + " e.r.");
            
            //formula 7.280
            c_psi[i] = -2 * Vector.getScalarProduct(l_Vec[i], obsR_Vec[i]);   
            //System.out.println("c_psi[" + i + "] " + c_psi[i]);
        }
        
        //before formula 7.281
        r_1g = r_1Appr;
        r_2g = r_2Appr;
        //r_1g = 2.0*a_e;//1.1*a_e;
        //r_2g = 2.1*a_e;//1.11*a_e;
        r_1 = r_1g;
        r_2 = r_2g;
        //System.out.println("r_1 " + r_1 + " e.r.");
        //System.out.println("r_2 " + r_2 + " e.r.");
        
        //observer position vector size obsR_Vec[i]
        double obsR_VecSize[] = new double[3];
        
        //loop 7.281 - 7.305/7.311
        double f_1Loop[] = new double[3];
        double f_2Loop[] = new double[3];
        
        boolean lookForValue = true;
        
        //before 7.312
        //double delta_r_PerCent = 0.0001;
        
        //condition 7.321
        //epsilon = a_e*10e-8;//a_e*10e-6;
        
        //number of iteration
        int iteration = 0;
        //int iterationBoundery = 50;
        
        //loop 2.281 - 7.322
        while(lookForValue){
            for(int l = 0; l < 3; l++){
                //System.out.println("------------- Case " + l + ". -----------------");
                //1st case - get F1(r_1,r_2) and F2(r_1,r_2)
                if(l == 0){
                    r_1 = r_1g;
                    r_2 = r_2g;
                }
                //2nd case - get F1{r_1+delta_r_1},r_2), F2(r_1+delta_r_1},r_2)
                else if(l == 1) {
                    r_1 = (1.0 + delta_r_PerCent)*r_1g;
                    r_2 = r_2g;
                }
                //3rd case - get F1{r_1,r_2+delta_r_1), F2(r_1,r_2+delta_r_1}
                else if(l == 2) {
                    r_1 = r_1g;
                    r_2 = (1.0 + delta_r_PerCent)*r_2g;                    
                }

                for(int i = 0; i < 2; i++){
                    obsR_VecSize[i] = Vector.getSize(obsR_Vec[i]);
                    //System.out.println("obsR_VecSize[" + i + "] " + obsR_VecSize[i] + " e.r.");
                    //formula 7.281
                    if(i == 0){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 * 
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1))));//r_1 !!!!!!!
                        //System.out.println("r_1 " + r_1);
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 * 
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    else if(i == 1){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 * 
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2))));//r_2 !!!!!!!
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 * 
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    //formula 7.282
                    Vector l_ro_Vec = Vector.multiplyVector(l_Vec[i], ro[i]);
                    //r_Vec[i] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[i], ro[i]),
                    r_Vec[i] = Vector.subtractVectors(l_ro_Vec, obsR_Vec[i]);
                    //System.out.println("l_ro_Vec[" + i + "].v[1]    " + l_ro_Vec.v[1]);
                    //System.out.println("obsR_Vec[" + i + "].v[1] " + obsR_Vec[i].v[1]);
                    //System.out.println("r_Vec[" + i + "].v[1]    " + r_Vec[i].v[1]);
                }

                //formulas 7.283
                //is motion retrograde
                int retrograde = 1;
                if(isMotionRetrograde) retrograde = -1;

                r_1 = Vector.getSize(r_Vec[0]);
                //System.out.println("1. r_1 " + r_1);
                //System.out.println("1. r_1 " + Vector.getSize(r_Vec[0]));
                r_2 = Vector.getSize(r_Vec[1]);
                //System.out.println("1. r_2 " + r_2);
                //System.out.println("1. r_2 " + Vector.getSize(r_Vec[1]));
                
                //w_Vec.v[0] = (r_Vec[0].v[1]*r_Vec[1].v[2] - r_Vec[1].v[1]*r_Vec[0].v[2])/
                //                (r_1*r_2);
                //w_Vec.v[1] = (r_Vec[1].v[0]*r_Vec[0].v[2] - r_Vec[0].v[0]*r_Vec[1].v[2])/
                //                (r_1*r_2); 
                //w_Vec.v[2] = (r_Vec[0].v[0]*r_Vec[1].v[1] - r_Vec[1].v[0]*r_Vec[0].v[1])/
                //                (r_1*r_2);
                
                w_Vec = Vector.getVectorProduct(r_Vec[0], r_Vec[1]);
                w_Vec = Vector.multiplyVector(w_Vec, 1.0/(r_1*r_2));
                
                //System.out.println("w_Vec.v[2] " + w_Vec.v[2]);

                //condition 20
                w_Vec = Vector.multiplyVector(w_Vec, retrograde);

                //formula 7.284
                double r3_w = Vector.getScalarProduct(obsR_Vec[2], w_Vec);
                double l3_w = Vector.getScalarProduct(l_Vec[2], w_Vec);
                ro[2] = r3_w/l3_w;

                //formula 7.285
                r_Vec[2] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[2], ro[2]),
                           obsR_Vec[2]);

                //formula 7.286
                r_3 = Vector.getSize(r_Vec[2]);
                
                //System.out.println("r_1 " + r_1 + " e.r.");
                //System.out.println("r_2 " + r_2 + " e.r.");
                //System.out.println("r_3 " + r_3 + " e.r.");

                // formulas 7.287
                cos_v2_v1 = Vector.getScalarProduct(r_Vec[1], r_Vec[0])/(r_2*r_1);
                //cos_v2_v2 = Vector.getScalarProduct(r_Vec[1], r_Vec[1])/(r_2*r_2);
                cos_v3_v1 = Vector.getScalarProduct(r_Vec[2], r_Vec[0])/(r_3*r_1);
                cos_v3_v2 = Vector.getScalarProduct(r_Vec[2], r_Vec[1])/(r_3*r_2);
                
                //System.out.println("cos_v2_v1 " + cos_v2_v1);
                //System.out.println("cos_v2_v2 " + cos_v2_v2);
                //System.out.println("cos_v3_v1 " + cos_v3_v1);
                //System.out.println("cos_v3_v2 " + cos_v3_v2);

                //formulas 7.288 - 7.289
                sin_v2_v1 =     (r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v2_v1*cos_v2_v1);
                //sin_v2_v2 = 0;//(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])/
                        //Math.abs(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])*
                        //Math.sqrt(1 - cos_v2_v2*cos_v2_v2);
                sin_v3_v1 =     (r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v3_v1*cos_v3_v1);
                sin_v3_v2 =     (r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])/
                        Math.abs(r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])*
                        Math.sqrt(1 - cos_v3_v2*cos_v3_v2);
                
                //System.out.println("sin_v2_v1 " + sin_v2_v1);
                //System.out.println("sin_v2_v2 " + sin_v2_v2);
                //System.out.println("sin_v3_v1 " + sin_v3_v1);
                //System.out.println("sin_v3_v2 " + sin_v3_v2);
                
                //System.out.println("v2 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1)));
                //System.out.println("v3 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1)));
                //System.out.println("v3 - v2 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v2, cos_v3_v2)));
                
                //condition before 7.288
                if(w_Vec.v[2] < 0){
                    //System.out.println("w_Vec.v[2] < 0");
                    //this is probably solved by method getAngleFromSinAndCos(sin, cos)
                    sin_v2_v1 = (-1) * sin_v2_v1;
                    //sin_v2_v2 = (-1) * sin_v2_v2;
                    sin_v3_v1 = (-1) * sin_v3_v1;
                    sin_v3_v2 = (-1) * sin_v3_v2;
                }

                //formula and condition 7.290
                //double v3_v1 = Math.acos(cos_v3_v1);
                //double v3_v1 = Math.acos(cos_v3_v1);
                double v3_v1 = getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1);
                //double v3_v1Sin = Math.asin(sin_v3_v1);
                //if(v3_v1Sin < 0) v3_v1 = 2*Math.PI - v3_v1;
                
                //System.out.println("v3_v1 " + Math.toDegrees(v3_v1) + " deg");
                if(v3_v1 > Math.PI){
                    //System.out.println("v3_v1 > PI !!!!");
                    //formulas 7.290
                    c_1 = (r_2/r_1)*(sin_v3_v2/sin_v3_v1);
                    c_3 = (r_2/r_3)*(sin_v2_v1/sin_v3_v1);
                    //formula 7.291
                    p = (c_1*r_1 + c_3*r_3 - r_2) / (c_1 + c_3 - 1);
                    //System.out.println("(c_1*r_1 + c_3*r_3 - r_2) " + (c_1*r_1 + c_3*r_3 - r_2));
                    //System.out.println("(c_1 + c_3 - 1) " + (c_1 + c_3 - 1));
                }
                //v3 - v1 <= Math.PI
                else{
                    //System.out.println("v3_v1 <= PI !!!!");
                    //formulas 7.292
                    cr_1 = (r_1/r_2)*(sin_v3_v1/sin_v3_v2);
                    cr_3 = (r_1/r_3)*(sin_v2_v1/sin_v3_v2);
                    //System.out.println("(sin_v3_v1/sin_v3_v2) " + (sin_v3_v1/sin_v3_v2));
                    //System.out.println("(sin_v2_v1/sin_v3_v2) " + (sin_v2_v1/sin_v3_v2));
                    //System.out.println("cr_1 " + cr_1 );
                    //System.out.println("cr_3 " + cr_3);
                    //formula 7.293
                    p = (r_1 + cr_3*r_3 - cr_1*r_2) / (1 + cr_3 - cr_1);
                    //System.out.println("(r_1 + cr_3*r_3 - cr_1*r_2) " + (r_1 + cr_3*r_3 - cr_1*r_2));
                    //System.out.println("(1 + cr_3 - cr_1) " + (1 + cr_3 - cr_1));
                }
                //System.out.println("----- Beware!!! p " + p);
                //formulas 7.294
                e_cosv_1 = p/r_1 - 1;
                e_cosv_2 = p/r_2 - 1;
                e_cosv_3 = p/r_3 - 1;

                //formula 7.295
                //double v2_v1 = Math.acos(cos_v2_v1);
                double v2_v1 = this.getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1);
                //condition before 7.295
                if(v2_v1 != Math.PI){
                    //formula 7.295
                    e_sinv_1 = (cos_v2_v1*e_cosv_1 - e_cosv_2)/sin_v2_v1;
                    e_sinv_2 = (-cos_v2_v1*e_cosv_2 + e_cosv_1)/sin_v2_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }        
                else{
                    //formula 7.296
                    e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }
                if(v3_v1!= Math.PI){
                    //e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                    //e_sinv_3 = (-cos_v3_v2*e_cosv_3 + e_cosv_2)/sin_v3_v1;
                } 

                //formula 7.297
                //e = Math.sqrt(e_cosv_1*e_cosv_1 + e_sinv_1*e_sinv_1);
                //System.out.println("e_1 " + e);
                e = Math.sqrt(e_cosv_2*e_cosv_2 + e_sinv_2*e_sinv_2);
                //System.out.println("e_2 " + e);
                //e = Math.sqrt(e_cosv_3*e_cosv_3 + e_sinv_3*e_sinv_3);
                //System.out.println("e_3 " + e);
                //e = 0.16419;
                //formula 7.298
                a = p/(1.0 - e*e);
                //System.out.println("a " + a);
                
                //System.out.println("asin(e_sinv_2/e) " + Math.asin(e_sinv_2/e));
                //System.out.println("acos(e_cosv_2/e) " + Math.acos(e_cosv_2/e));
                //System.out.println("asin + acos " + (Math.asin(e_sinv_2/e)+ Math.acos(e_cosv_2/e)));
                //System.out.println("asin - acos " + (Math.asin(e_sinv_2/e)- Math.acos(e_cosv_2/e)));
                //System.out.println("pi          " + (Math.PI));
                //System.out.println("v_2 " + this.getAngleFromSinAndCos(e_sinv_2/e, e_cosv_2/e));

                //elliptic case
                if(e*e < 1.0){
                    //System.out.println("Elliptic case!!! ------!!!--------");
                    //formula 7.299
                    n = k*Math.sqrt(mi)*Math.pow(a,-1.5);
                    //System.out.println("n " + n);
                    //formula 7.300
                    s_e = r_2/p * Math.sqrt(1 - e*e) * e_sinv_2;
                    //System.out.println("s_e " + s_e);
                    //formula 7.301
                    c_e = r_2/p * (e*e + e_cosv_2);
                    //System.out.println("c_e " + c_e);
                    //formulas 7.302
                    sin_E3_E2 = r_3*sin_v3_v2/Math.sqrt(a*p) - (r_3/p)*(1 - cos_v3_v2)*s_e;
                    cos_E3_E2 = 1 - (r_3*r_2/(a*p))*(1 - cos_v3_v2);
                    //.out.println("sin_E3_E2 " + sin_E3_E2);
                    //System.out.println("cos_E3_E2 " + cos_E3_E2);
                    //formulas 7.303
                    sin_E2_E1 = r_1*sin_v2_v1/Math.sqrt(a*p) + (r_1/p)*(1 - cos_v2_v1)*s_e;
                    cos_E2_E1 = 1 - (r_2*r_1/(a*p))*(1 - cos_v2_v1);
                    //System.out.println("sin_E2_E1 " + sin_E2_E1);
                    //System.out.println("cos_E2_E1 " + cos_E2_E1);
                    
                    e3_E2 = getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.PI*2 - getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.atan(sin_E3_E2 / cos_E3_E2);
                    e2_E1 = getAngleFromSinAndCos(sin_E2_E1, cos_E2_E1);
                    //e2_E1 = Math.atan(sin_E2_E1 / cos_E2_E1);
                    
                    //System.out.println("e3_E2 " + e3_E2);
                    //System.out.println("e2_E1 " + e2_E1);
                    
                    //formulas 7.304                            
                    m3_M2 = e3_E2 + 2*s_e*Math.pow(Math.sin(e3_E2/2),2) - c_e*sin_E3_E2;
                    m1_M2 = -e2_E1 + 2*s_e*Math.pow(Math.sin(e2_E1/2),2) + c_e*sin_E2_E1;
                    
                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);
                    //formulas 7.305
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;
                    
                    //System.out.println("k*(m1_M2/n) " + k*(m1_M2/n));
                    //System.out.println("tau_1 " + tau_1);
                    //System.out.println("k*(m3_M2/n) " + k*(m3_M2/n));
                    //System.out.println("tau_3 " + tau_3);
                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //hyperbolic case
                else{
                    //System.out.println("Hyperbolic case!!! ------!!!--------");
                    //System.out.println("a " + a);
                    //formula 7.306
                    n = k*Math.sqrt(mi)*Math.pow(-a, 1.5);
                    //System.out.println("n " + n);
                    //formulas 7.307
                    s_h = r_2/p * Math.sqrt(e*e - 1)*e_sinv_2;
                    c_h = r_2/p * (e*e + e_cosv_2);
                    
                    //System.out.println("s_h " + s_h);
                    //System.out.println("c_h " + c_h);
                    //formulas 7.308
                    double sinh_F3_F2, sinh_F2_F1;
                    double f3_F2, f2_F1;
                    sinh_F3_F2 = r_3/Math.sqrt(-a*p)*sin_v3_v2 - r_3/p*(1 - cos_v3_v2)*s_h;
                    sinh_F2_F1 = r_1/Math.sqrt(-a*p)*sin_v2_v1 + r_1/p*(1 - cos_v2_v1)*s_h;
                    
                    //System.out.println("sinh_F3_F2 " + sinh_F3_F2);
                    //System.out.println("sinh_F2_F1 " + sinh_F2_F1);
                    
                    //formulas 7.309
                    f3_F2 = Math.log10(sinh_F3_F2 + Math.sqrt(sinh_F3_F2*sinh_F3_F2 + 1));
                    f2_F1 = Math.log10(sinh_F2_F1 + Math.sqrt(sinh_F2_F1*sinh_F2_F1 + 1));
                    
                    //System.out.println("f3_F2 " + f3_F2);
                    //System.out.println("f2_F1 " + f2_F1);
                    
                    //formulas 7.310
                    //m3_M2 = -(f3_F2) + 2*s_h*(Math.sqrt(1 + sinh_F3_F2*sinh_F3_F2) - 1)/2.0 +    
                    //        c_h*sinh_F3_F2;
                    //m1_M2 = (f2_F1) + 2*s_h*(Math.sqrt(1 + sinh_F2_F1*sinh_F2_F1) - 1)/2.0 +    
                    //        c_h*sinh_F2_F1;
                    m3_M2 = -(f3_F2) + 2*s_h*Math.pow(Math.sinh(f3_F2/2),2) + c_h*Math.sinh(f3_F2);
                    m1_M2 =  (f2_F1) + 2*s_h*Math.pow(Math.sinh(f2_F1/2),2) - c_h*Math.sinh(f2_F1);
                    
                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);
                    
                    //formulas 7.311
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;;
                    
                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //1st/2nd/3rd case - l (low L)
                f_1Loop[l] = f_1;            
                f_2Loop[l] = f_2;            
            }
            
            for(int l = 0; l < 3; l++){
                 //System.out.println("f_1 f_1Loop[" + l + "]" + f_1Loop[l]);
                 //System.out.println("f_2 f_2Loop[" + l + "]" + f_2Loop[l]);
            }
            
            //System.out.println("-------------------- Loop 'for' OUT ------------------");
            
            //
            //double delta_r_PerCent = 0.04;
            
            //System.out.println("r_1g " + r_1g);
            //System.out.println("r_2g " + r_2g);
            
            //formula 7.312
            dF_1_DIV_dr_1 = (f_1Loop[1] - f_1Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.313
            dF_2_DIV_dr_1 = (f_2Loop[1] - f_2Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.314
            dF_1_DIV_dr_2 = (f_1Loop[2] - f_1Loop[0])/(delta_r_PerCent*r_2g);
            //formula 7.315
            dF_2_DIV_dr_2 = (f_2Loop[2] - f_2Loop[0])/(delta_r_PerCent*r_2g);
            
            //System.out.println("delta_r_PerCent*r_1g " + delta_r_PerCent*r_1g);
            //System.out.println("delta_r_PerCent*r_2g " + delta_r_PerCent*r_2g);
            
            //System.out.println("dF_1_DIV_dr_1 " + dF_1_DIV_dr_1);
            //System.out.println("dF_2_DIV_dr_1 " + dF_2_DIV_dr_1);
            //System.out.println("dF_1_DIV_dr_2 " + dF_1_DIV_dr_2);
            //System.out.println("dF_2_DIV_dr_2 " + dF_2_DIV_dr_2);

            //formula 7.316
            delta = dF_1_DIV_dr_1*dF_2_DIV_dr_2 - dF_2_DIV_dr_1*dF_1_DIV_dr_2;
            //formula 7.317
            delta_1 = dF_2_DIV_dr_2*f_1Loop[0] - dF_1_DIV_dr_2*f_2Loop[0];
            //System.out.println("dF_2_DIV_dr_2*f_1Loop[0] " + dF_2_DIV_dr_2*f_1Loop[0]);
            //System.out.println("dF_1_DIV_dr_2*f_2Loop[0] " + dF_1_DIV_dr_2*f_2Loop[0]);
            //formula 7.318
            delta_2 = dF_1_DIV_dr_1*f_2Loop[0] - dF_2_DIV_dr_1*f_1Loop[0];
            //System.out.println("dF_1_DIV_dr_1*f_2Loop[0] " + dF_1_DIV_dr_1*f_2Loop[0]);
            //System.out.println("dF_2_DIV_dr_1*f_1Loop[0] " + dF_2_DIV_dr_1*f_1Loop[0]);
            
            //formula 7.319
            deltar_1 = -delta_1/delta;
            //formula 7.320
            deltar_2 = -delta_2/delta;
            
            //System.out.println("delta " + delta);
            //System.out.println("delta_1 " + delta_1);
            //System.out.println("delta_2 " + delta_2);

            //System.out.println("deltar_1: " + deltar_1);
            //System.out.println("deltar_2: " + deltar_2);
            
            r_1 = r_1g + deltar_1;
            //System.out.println("r_2 " + r_2);
            r_2 = r_2g + deltar_2;
            //System.out.println("r_2 " + r_2);
                
            r_1g = r_1;                
            r_2g = r_2;
            
            //System.out.println("iteration " + iteration);
                
            //loop 7.321 - 7.322
            //condition 7.321
            if((Math.abs(deltar_1) < epsilon)&&(Math.abs(deltar_2) < epsilon)){
                lookForValue = false;
                break;
            }
            //formulas 7.322
            else{
                //r_1 = r_1g + deltar_1;
                //System.out.println("r_2 " + r_2);
                //r_2 = r_2g + deltar_2;
                //System.out.println("r_2 " + r_2);
                
                //r_1g = r_1;                
                //r_2g = r_2;                
            }
            
            //System.out.println("r_1: " + r_1);
            //System.out.println("r_2: " + r_2);
            
            //System.out.println(zmazat);
            if(p < 0) {
                //System.out.println("---------------------------- p is negative!!!!");
                break;
            }
            if((r_1 < a_e)||(r_2 < a_e)) {
                //System.out.println("---------------------------- r is too small!!!!");
                break;
            }
            
            iteration++;
            //System.out.println("iteration " + iteration);
            if(iteration > iterationBoundery) lookForValue = false;
        }
        
        //formula 7.323
        fSerie = 1.0 - a/r_2*(1 - cos_E3_E2);
        //formula 7.324
        gSerie = tau_3 - Math.pow(a, 1.5)/Math.sqrt(mi)*(e3_E2 - sin_E3_E2);
        
        //formula 7.325
        r_2_Dot = Vector.subtractVectors(r_Vec[2], Vector.multiplyVector(r_Vec[1], fSerie));
        r_2_Dot = Vector.multiplyVector(r_2_Dot, 1.0/gSerie);
        
        //State vector
        sv.r = r_Vec[1];
        sv.v = r_2_Dot;
        
        //System.out.println("Computed a = " + a);
        //System.out.println("Computed e = " + e);
        
        return sv;
    }
    /**
      * getStateVector_2()
      *
      * Escobal's double r-iteration method to get all position vectors (geocentric) and 2nd time velocity
      * vectors of body from 3 angles (R.A., dec)
      *
      * INPUT:
      *     double ra[i] - right accesions of body for time[i], i = 0,1,2
      *     double dec[i] - declinations of body for time[i], i = 0,1,2
      *     double time[i] - times of observation of body, Julian date [min] i = 0,1,2
      *     Geodetic geodetic[i] - geodetic positions of observer, i = 0,1,2
      *     double dTheta_dTime - Earth rotation
      *     double f - flattening of Earth
      *     double a_e - planet radius (equator)
      *     double mi -
      *     double k  -
      *     boolean isMotionRetrograde - is motion retrograde
      *     double r_1Appr - 1st input approximation of distance r_1 [a_e]
      *     double r_2Appr - 1st input approximation of distance r_2 [a_e]
      *     int lambda - number of revolutions
      *
      * OUTPUT:
      *     Vector array - 3. Position vectors and 2nd velocity vector
      *     vectorArray[0] - 1st position
      *     vectorArray[1] - 2nd position
      *     vectorArray[2] - 3rd position
      *     vectorArray[3] - 2nd velocity
      */
    public Vector[] getStateVector_2(double ra[], double dec[], double time[],
                Geodetic geodetic[], double dTheta_dTime, double flattening, double a_e,
                double mi, double k, boolean isMotionRetrograde, double r_1Appr, double r_2Appr,
                int lambda){

        this.k = k;
        this.mi = mi;
        this.ra = ra;
        this.dec = dec;
        this.time = time;
        this.geodetic = geodetic;
        this.dTheta_dTime = dTheta_dTime;
        this.flattening = flattening;
        this.a_e = a_e;
        this.isMotionRetrograde = isMotionRetrograde;
        this.r_1Appr = r_1Appr;
        this.r_2Appr = r_2Appr;
        this.lambda = lambda;

        //condition 7.321
        double epsilon = a_e*epsilonTolerance;//a_e*10e-6;
        //double epsilon = 10e-6;

        //initialization of main variable - info about position and velocity
        StateVector sv = new StateVector();

        // Formula 7.274
        tau_1 = k * (time[0] - time[1]);

        // Formula 7.275
        tau_3 = k * (time[2] - time[1]);

        //System.out.println("tau_1 " + tau_1);
        //System.out.println("tau_3 " + tau_3);

        //Greenwich siderial time - 1.27 [s]
        thetaG_2 = getGreenwichSiderialTime(time[1]);
        //System.out.println("thetaG_2 " + Math.toDegrees(thetaG_2) + " deg");

        for(int i = 0; i < 3; i++){
            // Formulas 7.276
            l_Vec[i] = new Vector(3);
            l_Vec[i].v[0] = Math.cos(dec[i])*Math.cos(ra[i]);
            l_Vec[i].v[1] = Math.cos(dec[i])*Math.sin(ra[i]);
            l_Vec[i].v[2] = Math.sin(dec[i]);

            //System.out.println("l_Vec[" + i + "] size " + Vector.getSize(l_Vec[i]));

            //Formulas 7.277
            g_1[i] = a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) +
                    geodetic[i].altitude;
            g_2[i] = (1 - flattening)*(1 - flattening)*a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) +
                    geodetic[i].altitude;

            //System.out.println("g_1[0] " + g_1[0]);
            //System.out.println("g_2[0] " + g_2[0]);

            //Formula 7.278
            theta[i] = thetaG_2 + dTheta_dTime*(time[i] - time[1]) + geodetic[i].lon;
            //System.out.println("theta[i] " + i + " " + theta[i]);

            //Formulas 7.279
            obsR_Vec[i] = new Vector(3);
            obsR_Vec[i].v[0] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.cos(theta[i]);
            obsR_Vec[i].v[1] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.sin(theta[i]);
            obsR_Vec[i].v[2] = -g_2[i] * Math.sin(geodetic[i].lat);

            //System.out.println("observer x " + obsR_Vec[i].v[0] + " e.r.");
            //System.out.println("observer y " + obsR_Vec[i].v[1] + " e.r.");
            //System.out.println("observer z " + obsR_Vec[i].v[2] + " e.r.");

            //formula 7.280
            c_psi[i] = -2 * Vector.getScalarProduct(l_Vec[i], obsR_Vec[i]);
            //System.out.println("c_psi[" + i + "] " + c_psi[i]);
        }

        //before formula 7.281
        r_1g = r_1Appr;
        r_2g = r_2Appr;
        //r_1g = 2.0*a_e;//1.1*a_e;
        //r_2g = 2.1*a_e;//1.11*a_e;
        r_1 = r_1g;
        r_2 = r_2g;
        //System.out.println("r_1 " + r_1 + " e.r.");
        //System.out.println("r_2 " + r_2 + " e.r.");

        //observer position vector size obsR_Vec[i]
        double obsR_VecSize[] = new double[3];

        //loop 7.281 - 7.305/7.311
        double f_1Loop[] = new double[3];
        double f_2Loop[] = new double[3];

        boolean lookForValue = true;

        //before 7.312
        //double delta_r_PerCent = 0.0001;

        //condition 7.321
        //epsilon = a_e*10e-8;//a_e*10e-6;

        //number of iteration
        int iteration = 0;
        //int iterationBoundery = 50;

        //loop 2.281 - 7.322
        while(lookForValue){
            for(int l = 0; l < 3; l++){
                //System.out.println("------------- Case " + l + ". -----------------");
                //1st case - get F1(r_1,r_2) and F2(r_1,r_2)
                if(l == 0){
                    r_1 = r_1g;
                    r_2 = r_2g;
                }
                //2nd case - get F1{r_1+delta_r_1},r_2), F2(r_1+delta_r_1},r_2)
                else if(l == 1) {
                    r_1 = (1.0 + delta_r_PerCent)*r_1g;
                    r_2 = r_2g;
                }
                //3rd case - get F1{r_1,r_2+delta_r_1), F2(r_1,r_2+delta_r_1}
                else if(l == 2) {
                    r_1 = r_1g;
                    r_2 = (1.0 + delta_r_PerCent)*r_2g;
                }

                for(int i = 0; i < 2; i++){
                    obsR_VecSize[i] = Vector.getSize(obsR_Vec[i]);
                    //System.out.println("obsR_VecSize[" + i + "] " + obsR_VecSize[i] + " e.r.");
                    //formula 7.281
                    if(i == 0){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 *
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1))));//r_1 !!!!!!!
                        //System.out.println("r_1 " + r_1);
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 *
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    else if(i == 1){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 *
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2))));//r_2 !!!!!!!
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 *
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    //formula 7.282
                    Vector l_ro_Vec = Vector.multiplyVector(l_Vec[i], ro[i]);
                    //r_Vec[i] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[i], ro[i]),
                    r_Vec[i] = Vector.subtractVectors(l_ro_Vec, obsR_Vec[i]);
                    //System.out.println("l_ro_Vec[" + i + "].v[1]    " + l_ro_Vec.v[1]);
                    //System.out.println("obsR_Vec[" + i + "].v[1] " + obsR_Vec[i].v[1]);
                    //System.out.println("r_Vec[" + i + "].v[1]    " + r_Vec[i].v[1]);
                }

                //formulas 7.283
                //is motion retrograde
                int retrograde = 1;
                if(isMotionRetrograde) retrograde = -1;

                r_1 = Vector.getSize(r_Vec[0]);
                //System.out.println("1. r_1 " + r_1);
                //System.out.println("1. r_1 " + Vector.getSize(r_Vec[0]));
                r_2 = Vector.getSize(r_Vec[1]);
                //System.out.println("1. r_2 " + r_2);
                //System.out.println("1. r_2 " + Vector.getSize(r_Vec[1]));

                //w_Vec.v[0] = (r_Vec[0].v[1]*r_Vec[1].v[2] - r_Vec[1].v[1]*r_Vec[0].v[2])/
                //                (r_1*r_2);
                //w_Vec.v[1] = (r_Vec[1].v[0]*r_Vec[0].v[2] - r_Vec[0].v[0]*r_Vec[1].v[2])/
                //                (r_1*r_2);
                //w_Vec.v[2] = (r_Vec[0].v[0]*r_Vec[1].v[1] - r_Vec[1].v[0]*r_Vec[0].v[1])/
                //                (r_1*r_2);

                w_Vec = Vector.getVectorProduct(r_Vec[0], r_Vec[1]);
                w_Vec = Vector.multiplyVector(w_Vec, 1.0/(r_1*r_2));

                //System.out.println("w_Vec.v[2] " + w_Vec.v[2]);

                //condition 20
                w_Vec = Vector.multiplyVector(w_Vec, retrograde);

                //formula 7.284
                double r3_w = Vector.getScalarProduct(obsR_Vec[2], w_Vec);
                double l3_w = Vector.getScalarProduct(l_Vec[2], w_Vec);
                ro[2] = r3_w/l3_w;

                //formula 7.285
                r_Vec[2] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[2], ro[2]),
                           obsR_Vec[2]);

                //formula 7.286
                r_3 = Vector.getSize(r_Vec[2]);

                //System.out.println("r_1 " + r_1 + " e.r.");
                //System.out.println("r_2 " + r_2 + " e.r.");
                //System.out.println("r_3 " + r_3 + " e.r.");

                // formulas 7.287
                cos_v2_v1 = Vector.getScalarProduct(r_Vec[1], r_Vec[0])/(r_2*r_1);
                //cos_v2_v2 = Vector.getScalarProduct(r_Vec[1], r_Vec[1])/(r_2*r_2);
                cos_v3_v1 = Vector.getScalarProduct(r_Vec[2], r_Vec[0])/(r_3*r_1);
                cos_v3_v2 = Vector.getScalarProduct(r_Vec[2], r_Vec[1])/(r_3*r_2);

                //System.out.println("cos_v2_v1 " + cos_v2_v1);
                //System.out.println("cos_v2_v2 " + cos_v2_v2);
                //System.out.println("cos_v3_v1 " + cos_v3_v1);
                //System.out.println("cos_v3_v2 " + cos_v3_v2);

                //formulas 7.288 - 7.289
                sin_v2_v1 =     (r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v2_v1*cos_v2_v1);
                //sin_v2_v2 = 0;//(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])/
                        //Math.abs(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])*
                        //Math.sqrt(1 - cos_v2_v2*cos_v2_v2);
                sin_v3_v1 =     (r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v3_v1*cos_v3_v1);
                sin_v3_v2 =     (r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])/
                        Math.abs(r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])*
                        Math.sqrt(1 - cos_v3_v2*cos_v3_v2);

                //System.out.println("sin_v2_v1 " + sin_v2_v1);
                //System.out.println("sin_v2_v2 " + sin_v2_v2);
                //System.out.println("sin_v3_v1 " + sin_v3_v1);
                //System.out.println("sin_v3_v2 " + sin_v3_v2);

                //System.out.println("v2 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1)));
                //System.out.println("v3 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1)));
                //System.out.println("v3 - v2 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v2, cos_v3_v2)));

                //condition before 7.288
                if(w_Vec.v[2] < 0){
                    //System.out.println("w_Vec.v[2] < 0");
                    //this is probably solved by method getAngleFromSinAndCos(sin, cos)
                    sin_v2_v1 = (-1) * sin_v2_v1;
                    //sin_v2_v2 = (-1) * sin_v2_v2;
                    sin_v3_v1 = (-1) * sin_v3_v1;
                    sin_v3_v2 = (-1) * sin_v3_v2;
                }

                //formula and condition 7.290
                //double v3_v1 = Math.acos(cos_v3_v1);
                //double v3_v1 = Math.acos(cos_v3_v1);
                double v3_v1 = getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1);
                //double v3_v1Sin = Math.asin(sin_v3_v1);
                //if(v3_v1Sin < 0) v3_v1 = 2*Math.PI - v3_v1;

                //System.out.println("v3_v1 " + Math.toDegrees(v3_v1) + " deg");
                if(v3_v1 > Math.PI){
                    //System.out.println("v3_v1 > PI !!!!");
                    //formulas 7.290
                    c_1 = (r_2/r_1)*(sin_v3_v2/sin_v3_v1);
                    c_3 = (r_2/r_3)*(sin_v2_v1/sin_v3_v1);
                    //formula 7.291
                    p = (c_1*r_1 + c_3*r_3 - r_2) / (c_1 + c_3 - 1);
                    //System.out.println("(c_1*r_1 + c_3*r_3 - r_2) " + (c_1*r_1 + c_3*r_3 - r_2));
                    //System.out.println("(c_1 + c_3 - 1) " + (c_1 + c_3 - 1));
                }
                //v3 - v1 <= Math.PI
                else{
                    //System.out.println("v3_v1 <= PI !!!!");
                    //formulas 7.292
                    cr_1 = (r_1/r_2)*(sin_v3_v1/sin_v3_v2);
                    cr_3 = (r_1/r_3)*(sin_v2_v1/sin_v3_v2);
                    //System.out.println("(sin_v3_v1/sin_v3_v2) " + (sin_v3_v1/sin_v3_v2));
                    //System.out.println("(sin_v2_v1/sin_v3_v2) " + (sin_v2_v1/sin_v3_v2));
                    //System.out.println("cr_1 " + cr_1 );
                    //System.out.println("cr_3 " + cr_3);
                    //formula 7.293
                    p = (r_1 + cr_3*r_3 - cr_1*r_2) / (1 + cr_3 - cr_1);
                    //System.out.println("(r_1 + cr_3*r_3 - cr_1*r_2) " + (r_1 + cr_3*r_3 - cr_1*r_2));
                    //System.out.println("(1 + cr_3 - cr_1) " + (1 + cr_3 - cr_1));
                }
                //System.out.println("----- Beware!!! p " + p);
                //formulas 7.294
                e_cosv_1 = p/r_1 - 1;
                e_cosv_2 = p/r_2 - 1;
                e_cosv_3 = p/r_3 - 1;

                //formula 7.295
                //double v2_v1 = Math.acos(cos_v2_v1);
                double v2_v1 = this.getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1);
                //condition before 7.295
                if(v2_v1 != Math.PI){
                    //formula 7.295
                    e_sinv_1 = (cos_v2_v1*e_cosv_1 - e_cosv_2)/sin_v2_v1;
                    e_sinv_2 = (-cos_v2_v1*e_cosv_2 + e_cosv_1)/sin_v2_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }
                else{
                    //formula 7.296
                    e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }
                if(v3_v1!= Math.PI){
                    //e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                    //e_sinv_3 = (-cos_v3_v2*e_cosv_3 + e_cosv_2)/sin_v3_v1;
                }

                //formula 7.297
                //e = Math.sqrt(e_cosv_1*e_cosv_1 + e_sinv_1*e_sinv_1);
                //System.out.println("e_1 " + e);
                e = Math.sqrt(e_cosv_2*e_cosv_2 + e_sinv_2*e_sinv_2);
                //System.out.println("e_2 " + e);
                //e = Math.sqrt(e_cosv_3*e_cosv_3 + e_sinv_3*e_sinv_3);
                //System.out.println("e_3 " + e);
                //e = 0.16419;
                //formula 7.298
                a = p/(1.0 - e*e);
                //System.out.println("a " + a);

                //e = 0.00001;
                //a = 63963481/Constants.R_Earth;

                //System.out.println("asin(e_sinv_2/e) " + Math.asin(e_sinv_2/e));
                //System.out.println("acos(e_cosv_2/e) " + Math.acos(e_cosv_2/e));
                //System.out.println("asin + acos " + (Math.asin(e_sinv_2/e)+ Math.acos(e_cosv_2/e)));
                //System.out.println("asin - acos " + (Math.asin(e_sinv_2/e)- Math.acos(e_cosv_2/e)));
                //System.out.println("pi          " + (Math.PI));
                //System.out.println("v_2 " + this.getAngleFromSinAndCos(e_sinv_2/e, e_cosv_2/e));

                //elliptic case
                if(e*e < 1.0){
                    //System.out.println("Elliptic case!!! ------!!!--------");
                    //formula 7.299
                    n = k*Math.sqrt(mi)*Math.pow(a,-1.5);
                    //System.out.println("n " + n);
                    //formula 7.300
                    s_e = r_2/p * Math.sqrt(1 - e*e) * e_sinv_2;
                    //System.out.println("s_e " + s_e);
                    //formula 7.301
                    c_e = r_2/p * (e*e + e_cosv_2);
                    //System.out.println("c_e " + c_e);
                    //formulas 7.302
                    sin_E3_E2 = r_3*sin_v3_v2/Math.sqrt(a*p) - (r_3/p)*(1 - cos_v3_v2)*s_e;
                    cos_E3_E2 = 1 - (r_3*r_2/(a*p))*(1 - cos_v3_v2);
                    //.out.println("sin_E3_E2 " + sin_E3_E2);
                    //System.out.println("cos_E3_E2 " + cos_E3_E2);
                    //formulas 7.303
                    sin_E2_E1 = r_1*sin_v2_v1/Math.sqrt(a*p) + (r_1/p)*(1 - cos_v2_v1)*s_e;
                    cos_E2_E1 = 1 - (r_2*r_1/(a*p))*(1 - cos_v2_v1);
                    //System.out.println("sin_E2_E1 " + sin_E2_E1);
                    //System.out.println("cos_E2_E1 " + cos_E2_E1);

                    e3_E2 = getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.PI*2 - getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.atan(sin_E3_E2 / cos_E3_E2);
                    e2_E1 = getAngleFromSinAndCos(sin_E2_E1, cos_E2_E1);
                    //e2_E1 = Math.atan(sin_E2_E1 / cos_E2_E1);

                    //System.out.println("e3_E2 " + e3_E2);
                    //System.out.println("e2_E1 " + e2_E1);

                    //formulas 7.304
                    m3_M2 = e3_E2 + 2*s_e*Math.pow(Math.sin(e3_E2/2),2) - c_e*sin_E3_E2;
                    m1_M2 = -e2_E1 + 2*s_e*Math.pow(Math.sin(e2_E1/2),2) + c_e*sin_E2_E1;

                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);
                    //formulas 7.305
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;

                    //System.out.println("k*(m1_M2/n) " + k*(m1_M2/n));
                    //System.out.println("tau_1 " + tau_1);
                    //System.out.println("k*(m3_M2/n) " + k*(m3_M2/n));
                    //System.out.println("tau_3 " + tau_3);
                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //hyperbolic case
                else{
                    //System.out.println("Hyperbolic case!!! ------!!!--------");
                    //System.out.println("a " + a);
                    //formula 7.306
                    n = k*Math.sqrt(mi)*Math.pow(-a, 1.5);
                    //System.out.println("n " + n);
                    //formulas 7.307
                    s_h = r_2/p * Math.sqrt(e*e - 1)*e_sinv_2;
                    c_h = r_2/p * (e*e + e_cosv_2);

                    //System.out.println("s_h " + s_h);
                    //System.out.println("c_h " + c_h);
                    //formulas 7.308
                    double sinh_F3_F2, sinh_F2_F1;
                    double f3_F2, f2_F1;
                    sinh_F3_F2 = r_3/Math.sqrt(-a*p)*sin_v3_v2 - r_3/p*(1 - cos_v3_v2)*s_h;
                    sinh_F2_F1 = r_1/Math.sqrt(-a*p)*sin_v2_v1 + r_1/p*(1 - cos_v2_v1)*s_h;

                    //System.out.println("sinh_F3_F2 " + sinh_F3_F2);
                    //System.out.println("sinh_F2_F1 " + sinh_F2_F1);

                    //formulas 7.309
                    f3_F2 = Math.log10(sinh_F3_F2 + Math.sqrt(sinh_F3_F2*sinh_F3_F2 + 1));
                    f2_F1 = Math.log10(sinh_F2_F1 + Math.sqrt(sinh_F2_F1*sinh_F2_F1 + 1));

                    //System.out.println("f3_F2 " + f3_F2);
                    //System.out.println("f2_F1 " + f2_F1);

                    //formulas 7.310
                    //m3_M2 = -(f3_F2) + 2*s_h*(Math.sqrt(1 + sinh_F3_F2*sinh_F3_F2) - 1)/2.0 +
                    //        c_h*sinh_F3_F2;
                    //m1_M2 = (f2_F1) + 2*s_h*(Math.sqrt(1 + sinh_F2_F1*sinh_F2_F1) - 1)/2.0 +
                    //        c_h*sinh_F2_F1;
                    m3_M2 = -(f3_F2) + 2*s_h*Math.pow(Math.sinh(f3_F2/2),2) + c_h*Math.sinh(f3_F2);
                    m1_M2 =  (f2_F1) + 2*s_h*Math.pow(Math.sinh(f2_F1/2),2) - c_h*Math.sinh(f2_F1);

                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);

                    //formulas 7.311
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;;

                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //1st/2nd/3rd case - l (low L)
                f_1Loop[l] = f_1;
                f_2Loop[l] = f_2;
            }

            for(int l = 0; l < 3; l++){
                 //System.out.println("f_1 f_1Loop[" + l + "]" + f_1Loop[l]);
                 //System.out.println("f_2 f_2Loop[" + l + "]" + f_2Loop[l]);
            }

            //System.out.println("-------------------- Loop 'for' OUT ------------------");

            //
            //double delta_r_PerCent = 0.04;

            //System.out.println("r_1g " + r_1g);
            //System.out.println("r_2g " + r_2g);

            //formula 7.312
            dF_1_DIV_dr_1 = (f_1Loop[1] - f_1Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.313
            dF_2_DIV_dr_1 = (f_2Loop[1] - f_2Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.314
            dF_1_DIV_dr_2 = (f_1Loop[2] - f_1Loop[0])/(delta_r_PerCent*r_2g);
            //formula 7.315
            dF_2_DIV_dr_2 = (f_2Loop[2] - f_2Loop[0])/(delta_r_PerCent*r_2g);

            //System.out.println("delta_r_PerCent*r_1g " + delta_r_PerCent*r_1g);
            //System.out.println("delta_r_PerCent*r_2g " + delta_r_PerCent*r_2g);

            //System.out.println("dF_1_DIV_dr_1 " + dF_1_DIV_dr_1);
            //System.out.println("dF_2_DIV_dr_1 " + dF_2_DIV_dr_1);
            //System.out.println("dF_1_DIV_dr_2 " + dF_1_DIV_dr_2);
            //System.out.println("dF_2_DIV_dr_2 " + dF_2_DIV_dr_2);

            //formula 7.316
            delta = dF_1_DIV_dr_1*dF_2_DIV_dr_2 - dF_2_DIV_dr_1*dF_1_DIV_dr_2;
            //formula 7.317
            delta_1 = dF_2_DIV_dr_2*f_1Loop[0] - dF_1_DIV_dr_2*f_2Loop[0];
            //System.out.println("dF_2_DIV_dr_2*f_1Loop[0] " + dF_2_DIV_dr_2*f_1Loop[0]);
            //System.out.println("dF_1_DIV_dr_2*f_2Loop[0] " + dF_1_DIV_dr_2*f_2Loop[0]);
            //formula 7.318
            delta_2 = dF_1_DIV_dr_1*f_2Loop[0] - dF_2_DIV_dr_1*f_1Loop[0];
            //System.out.println("dF_1_DIV_dr_1*f_2Loop[0] " + dF_1_DIV_dr_1*f_2Loop[0]);
            //System.out.println("dF_2_DIV_dr_1*f_1Loop[0] " + dF_2_DIV_dr_1*f_1Loop[0]);

            //formula 7.319
            deltar_1 = -delta_1/delta;
            //formula 7.320
            deltar_2 = -delta_2/delta;

            //System.out.println("delta " + delta);
            //System.out.println("delta_1 " + delta_1);
            //System.out.println("delta_2 " + delta_2);

            //System.out.println("deltar_1: " + deltar_1);
            //System.out.println("deltar_2: " + deltar_2);

            r_1 = r_1g + deltar_1;
            //System.out.println("r_1 " + r_1);
            r_2 = r_2g + deltar_2;
            //System.out.println("r_2 " + r_2);

            r_1g = r_1;
            r_2g = r_2;

            //System.out.println("iteration " + iteration);

            //loop 7.321 - 7.322
            //condition 7.321
            if((Math.abs(deltar_1) < epsilon)&&(Math.abs(deltar_2) < epsilon)){
                lookForValue = false;
                break;
            }
            //formulas 7.322
            else{
                //r_1 = r_1g + deltar_1;
                //System.out.println("r_2 " + r_2);
                //r_2 = r_2g + deltar_2;
                //System.out.println("r_2 " + r_2);

                //r_1g = r_1;
                //r_2g = r_2;
            }

            //System.out.println("r_1: " + r_1);
            //System.out.println("r_2: " + r_2);

            //System.out.println(zmazat);
            if(p < 0) {
                //System.out.println("---------------------------- p is negative!!!!");
                break;
            }
            if((r_1 < a_e)||(r_2 < a_e)) {
                //System.out.println("---------------------------- r is too small!!!!");
                break;
            }

            iteration++;
            //System.out.println("iteration " + iteration);
            if(iteration > iterationBoundery) lookForValue = false;
        }

        //formula 7.323
        fSerie = 1.0 - a/r_2*(1 - cos_E3_E2);
        //formula 7.324
        gSerie = tau_3 - Math.pow(a, 1.5)/Math.sqrt(mi)*(e3_E2 - sin_E3_E2);
        //System.out.println("FSerie 1 " + fSerie);
        //System.out.println("GSerie 1 " + gSerie);
        //System.out.println("Tau 1" + tau_3);

        //formula 7.325
        r_2_Dot = Vector.subtractVectors(r_Vec[2], Vector.multiplyVector(r_Vec[1], fSerie));
        r_2_Dot = Vector.multiplyVector(r_2_Dot, 1.0/gSerie);

        //State vector
        sv.r = r_Vec[1];
        sv.v = r_2_Dot;

        //System.out.println("Computed a = " + a);
        //System.out.println("Computed e = " + e);

        //vector array with vectors of 3 positions and 1 velocity
        Vector vectorArray[] = new Vector[4];
        vectorArray[0] = r_Vec[0];
        vectorArray[1] = r_Vec[1];
        vectorArray[2] = r_Vec[2];
        vectorArray[3] = r_2_Dot;

        //System.out.println("testLala2 Vel [x] " + vectorArray[3].v[0]);
        //System.out.println("testLala2 Vel [y] " + vectorArray[3].v[1]);
        //System.out.println("testLala2 Vel [z] " + vectorArray[3].v[2]);
        //System.out.println("vel_2 " + Vector.getSize(vectorArray[3]));

        return vectorArray;
    }

    /**
     * getVectors()
     *
     * Method to compute 3 position and 1 velocity vectors.
     *
     * INPUT:
     *  Observation observation[] - see constructor Observation
     *
     * OUTPUT:
     *  Vector vector[4] - Array of vectors:
     *               - 1st position vector for 1st time [m]
     *               - 1st position vector for 1st time [m]
     *               - 1st position vector for 1st time [m]
     *               - velocity Vector for 2nd time     [m]
     */
    public Vector[] getVectors(Observation observation[]){
        Vector vectors[] = new Vector[4];

        //Array R.A.
        double ra[] = new double[3];
        //Array declination
        double dec[] = new double[3];
        //Array time [Mjd]
        double timeJd[] = new double[3];
        //Array Geodetic
        Geodetic geodetic[] = new Geodetic[3];

        for(int i = 0; i<3; i++){
            ra[i] = observation[i].ra;
            dec[i] = observation[i].dec;
            timeJd[i] = new Time().getJdFromMjd(observation[i].timeMjd)*1440;

            geodetic[i] = new Geodetic();
            geodetic[i].lon = observation[i].lon;
            geodetic[i].lat = observation[i].lat;
            geodetic[i].altitude = observation[i].alt/Constants.R_Earth;
        }

        //constants
        double k2 = 0.07436574;  //[(e.r)^3/2 / min]
        double mi2 = 1.0;        //[e.m.] - Earth mass
        double a_e2 = 1.0;       //[e.r.] - Earth radius
        double dTheta_dTime2 = 4.3752695e-3; //[rad/min]
        double flattening2 = Constants.f_Earth;

        //CAREFULL !!!!
        boolean isMotionRetrograde = false;
        double r_1Appr = 20;
        double r_2Appr = 20.1;
        //0 revolutions
        int lambda = 0;

        //get the position vectors and velocity vector
        vectors = getStateVector_2(ra, dec, timeJd, geodetic, dTheta_dTime2, flattening, a_e2,
                mi2, k2, isMotionRetrograde, r_1Appr, r_2Appr, lambda);

        //from Earth radius to [m], and from Escobal units to [m/s], koeficient * 25936*0.3048
        //positions vectors
        for(int i=0; i < 3; i++){
            for(int j=0; j < 3; j++){
                vectors[i].v[j] = vectors[i].v[j]*Constants.R_Earth;
                //System.out.println("Geo vec   " + i  + " " + j + ": "+ vectors[i].v[j]);
            }
        }
        //velocity vector
        for(int i=0; i < 3; i++){
            vectors[3].v[i] = vectors[3].v[i]*25936*0.3048;
        }

        //geocentric vectors
        //return vectors;
        //double size = 10e14;
        //System.out.println("New vectors: ");
        //positions vectors helio
        //for(int i=0; i < 3; i++){
            //for(int j=0; j < 3; j++){
                //vectors[i].v[j] = vectors[i].v[j]*Constants.R_Earth;
                //System.out.println("Geo vec   " + i  + " " + j + ": "+ vectors[i].v[j]);
                //vectors[i] = Transformation.fromGeocentricToHeliocentric(vectors[i], new Time().getJdFromMjd(observation[i].timeMjd));
                //size = Vector.getSize(vectors[i]);
                //if(size<10e12) System.out.println("Helio vec " + i  + " " + j + ": "+ vectors[i].v[j]);
            //}
            //size = Vector.getSize(vectors[i]);
            //if(size<10e12) System.out.println("Size " + size/1000 + " km");
        //}

        //Vector eclipNormal = new Vector(3);
        //eclipNormal.v[0] = 0;
        //eclipNormal.v[1] = 0;
        //eclipNormal.v[2] = 1;
        //double inclination = Vector.getVectorsAngle(eclipNormal, Vector.getVectorProduct(vectors[0], vectors[2]));
        //if(size<10e12)System.out.println("Inclination " + Math.toDegrees(inclination));
        
        return vectors;
    }

     /**
     * getGreenwichSiderialTime()
     * Formula 1.27
     *
     * IN:
     *  double timeJd - jullian date [???]
     *
     * OUT:
     *  double gst [rad]
     */

    public double getGreenwichSiderialTime(double timeJd){
        /*timeJd = timeJd/1440;   //[min] -> [day]
        System.out.println("timeJd " + timeJd);
        double t_u = (timeJd - 2415020.0)/36525.0;
        //double t_u = (timeJd - 2415020.0)/36524.219879;
        System.out.println("t_u " + t_u);
        System.out.println("36000.7689*t_u " + 36000.7689*t_u);
        System.out.println("0.00038708*t_u*t_u " + 0.00038708*t_u*t_u);
        double gst = 99.6909833 + 36000.7689*t_u + 0.00038708*t_u*t_u;
        gst = Math.toRadians(gst);
        return gst;*/
        
        timeJd = timeJd/1440;   //[min] -> [day]
        //System.out.println("timeJd " + timeJd);
        double helpJd_1 = timeJd - (int)timeJd;
        double helpJd_2 = (int)timeJd;
        //System.out.println("helpJd_2 " + helpJd_2);
        if(helpJd_1 >= 0.5) {
            helpJd_1 = helpJd_1 - 0.5;
            helpJd_2 = helpJd_2 + 0.5;
        }
        else {
            helpJd_1 = helpJd_1 + 0.5;
            helpJd_2 = helpJd_2 - 0.5;
        }
        
        helpJd_1 = helpJd_1*1440;   //day -> min
        
        double t_u = (helpJd_2 - 2415020.0)/36525.0;
        //double t_u = (timeJd - 2415020.0)/36524.219879;
        //System.out.println("t_u " + t_u);
        //System.out.println("36000.7689*t_u " + 36000.7689*t_u);
        //System.out.println("0.00038708*t_u*t_u " + 0.00038708*t_u*t_u);
        double gst = 99.6909833 + 36000.7689*t_u + 0.00038708*t_u*t_u;
        gst = gst + + helpJd_1*0.25068447;
        gst = Math.toRadians(gst);
        return gst;
    }
    
    /**
     * getGreenwichSiderialTime()
     * Formula 1.27
     *
     * IN:
     *  double timeJd - jullian date [day] (days till midnight)
     *  double min - minutes since midnight
     *
     * OUT:
     *  double gst [rad]
     */

    public double getGreenwichSiderialTime(double timeJd, double min){
        //timeJd = timeJd/1440;   //[min] -> [day]
        //System.out.println("timeJd " + timeJd);
        double t_u = (timeJd - 2415020.0)/36525.0;
        //double t_u = (timeJd - 2415020.0)/36524.219879;
        //System.out.println("t_u " + t_u);
        //System.out.println("36000.7689*t_u " + 36000.7689*t_u);
        //System.out.println("0.00038708*t_u*t_u " + 0.00038708*t_u*t_u);
        double gst = 99.6909833 + 36000.7689*t_u + 0.00038708*t_u*t_u;
        gst = gst + min*0.25068447;
        gst = Math.toRadians(gst);
        return gst;
    }
    
    /**
     * getAngleFromSinAndCos()
     * 
     * Method to compute angle from cos nagle and sin angle, resulat is angle [rad]
     * in interval <0;2*PI>
     * 
     * INPUT:
     *  double sinAngle - sine of angle
     *  double cosAngle - cosine of angle
     * 
     * OUTPUT:
     *  double angle - [rad], <0;2*Pi>
     */
    
   public double getAngleFromSinAndCos(double sinAngle, double cosAngle){
       double angle = 0;
       if((sinAngle >= 0)&&(cosAngle >= 0)) angle = Math.atan(sinAngle/cosAngle);
       else if((sinAngle < 0)&&(cosAngle >= 0)) angle = 2*Math.PI + Math.atan(sinAngle/cosAngle);
       else if((sinAngle >= 0)&&(cosAngle < 0)) angle = Math.PI + Math.atan(sinAngle/cosAngle);
       else if((sinAngle < 0)&&(cosAngle < 0)) angle = Math.PI + Math.atan(sinAngle/cosAngle);
       
       return angle;
   }

   /**
    * getVelVecFrom3PosVec()
    *
    * Method to compute velocity vector for 2nd time from 3 position vectors and their associated times.
    * Using formulas in method getStateVector_2()
    *
    * INPUT:
    *  Vector r_Vec[] - array of 3 position vectors [in Earth radius]
    *  double jd[] - 3 associated times for given 3 position vectors, Julian date [day]
    *  double k - see Escobal
    *  double mi - see Escobal
    *  boolean isRetrograde - wheter or not is orbit retrograde
    *
    * OUPUT:
    *  Vector velVec2 - velocity vedtor for second time, in [Escobal units], to get SI multiply with  * 25936*0.3048 --> [m/s]
    */
   public Vector getVelVecFrom3PosVec(Vector r_Vec[], double jd[],
                            double k, double mi, boolean isRetrograde){
       Vector velVec2 = new Vector(3);

       /*
        System.out.println("testik x1 " + r_Vec[0].v[0]);
		System.out.println("y1 " + r_Vec[0].v[1]);
		System.out.println("z1 " + r_Vec[0].v[2]);
                System.out.println("x2 " + r_Vec[1].v[0]);
		System.out.println("y2 " + r_Vec[1].v[1]);
		System.out.println("z2 " + r_Vec[1].v[2]);
		System.out.println("x3 " + r_Vec[2].v[0]);
                System.out.println("y3 " + r_Vec[2].v[1]);
		System.out.println("z3 " + r_Vec[2].v[2]);
       */
       //times should be in [min]
       double jdMin[] = new double[3];
       for(int i =0; i <3; i++){
           jdMin[i] = jd[i]*1440;
       }
       // Formula 7.275
       double tau_3 = k * (jdMin[2] - jdMin[1]);
       //sizes of position vectors
       double rSize[] = new double[3];
       for(int i=0; i<3; i++){
           rSize[i] = r_Vec[i].getSize(r_Vec[i]);
       }
       r_1 = rSize[0];
       r_2 = rSize[1];
       r_3 = rSize[2];

       //vector W (see Escobal)
       w_Vec = Vector.getVectorProduct(r_Vec[0], r_Vec[1]);
       w_Vec = Vector.multiplyVector(w_Vec, 1.0/(r_1*r_2));
       int retrograde = 1;
       if(isRetrograde) retrograde = -1;
       //condition 20
       w_Vec = Vector.multiplyVector(w_Vec, retrograde);
       
       // formulas 7.287
       cos_v2_v1 = Vector.getScalarProduct(r_Vec[1], r_Vec[0])/(r_2*r_1);
       cos_v3_v1 = Vector.getScalarProduct(r_Vec[2], r_Vec[0])/(r_3*r_1);
       cos_v3_v2 = Vector.getScalarProduct(r_Vec[2], r_Vec[1])/(r_3*r_2);

       //formulas 7.288 - 7.289
       sin_v2_v1 =     (r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])/
                  Math.abs(r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])*
                  Math.sqrt(1 - cos_v2_v1*cos_v2_v1);
       sin_v3_v1 =     (r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])/
                  Math.abs(r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])*
                  Math.sqrt(1 - cos_v3_v1*cos_v3_v1);
       sin_v3_v2 =     (r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])/
                  Math.abs(r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])*
                  Math.sqrt(1 - cos_v3_v2*cos_v3_v2);

       //condition before 7.288
       if(w_Vec.v[2] < 0){
            //System.out.println("w_Vec.v[2] < 0");
            //this is probably solved by method getAngleFromSinAndCos(sin, cos)
            sin_v2_v1 = (-1) * sin_v2_v1;
            //sin_v2_v2 = (-1) * sin_v2_v2;
            sin_v3_v1 = (-1) * sin_v3_v1;
            sin_v3_v2 = (-1) * sin_v3_v2;
       }

       double v3_v1 = getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1);
       if(v3_v1 > Math.PI){
            //formulas 7.290
            c_1 = (r_2/r_1)*(sin_v3_v2/sin_v3_v1);
            c_3 = (r_2/r_3)*(sin_v2_v1/sin_v3_v1);
            //formula 7.291
            p = (c_1*r_1 + c_3*r_3 - r_2) / (c_1 + c_3 - 1);
           }
       else{
          //System.out.println("v3_v1 <= PI !!!!");
          //formulas 7.292
          cr_1 = (r_1/r_2)*(sin_v3_v1/sin_v3_v2);
          cr_3 = (r_1/r_3)*(sin_v2_v1/sin_v3_v2);
          p = (r_1 + cr_3*r_3 - cr_1*r_2) / (1 + cr_3 - cr_1);
       }
       //formulas 7.294
       e_cosv_1 = p/r_1 - 1;
       e_cosv_2 = p/r_2 - 1;
       e_cosv_3 = p/r_3 - 1;

       //formula 7.295
       double v2_v1 = this.getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1);
       //condition before 7.295
       if(v2_v1 != Math.PI){
            //formula 7.295
            e_sinv_1 = (cos_v2_v1*e_cosv_1 - e_cosv_2)/sin_v2_v1;
            e_sinv_2 = (-cos_v2_v1*e_cosv_2 + e_cosv_1)/sin_v2_v1;
       }
       else{
            //formula 7.296
            e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
       }

       //formula 7.297
       e = Math.sqrt(e_cosv_2*e_cosv_2 + e_sinv_2*e_sinv_2);
       //formula 7.298
       a = p/(1.0 - e*e);
       //System.out.println("a test " + a);

       sin_E3_E2 = r_3*sin_v3_v2/Math.sqrt(a*p) - (r_3/p)*(1 - cos_v3_v2)*s_e;
       cos_E3_E2 = 1 - (r_3*r_2/(a*p))*(1 - cos_v3_v2);
       e3_E2 = getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);

       //formula 7.323
       fSerie = 1.0 - a/r_2*(1 - cos_E3_E2);
       //formula 7.324
       gSerie = tau_3 - Math.pow(a, 1.5)/Math.sqrt(mi)*(e3_E2 - sin_E3_E2);

       //System.out.println("FSerie 2 " + fSerie);
       //System.out.println("GSerie 2 " + gSerie);
       //System.out.println("Tau 2 " + tau_3);

       //formula 7.325
       r_2_Dot = Vector.subtractVectors(r_Vec[2], Vector.multiplyVector(r_Vec[1], fSerie));
       r_2_Dot = Vector.multiplyVector(r_2_Dot, 1.0/gSerie);
       //System.out.println("vx " + r_2_Dot.v[0]);
       //System.out.println("vy " + r_2_Dot.v[1]);
       //System.out.println("vz " + r_2_Dot.v[2]);
       //from Escobal units to SI [m/s]
       //for(int i =0 ; i < 3; i++){
       //    velVec2.v[i] = r_2_Dot.v[i] * 25936*0.3048;
       //}
       velVec2 = r_2_Dot;

       //System.out.println("testLala Vel [x] " + velVec2.v[0]);
       //System.out.println("testLala Vel [y] " + velVec2.v[1]);
       //System.out.println("testLala Vel [z] " + velVec2.v[2]);
            
       //out Escobal units
       return velVec2;
   }

   /**
      * getPositionVector()
      *
      * Escobal's double r-iteration method to get all position vectors (geocentric) and 2nd time velocity
      * vectors of body from 3 angles (R.A., dec)
      *
      * INPUT:
      *     double ra[i] - right accesions of body for time[i], i = 0,1,2
      *     double dec[i] - declinations of body for time[i], i = 0,1,2
      *     double time[i] - times of observation of body, Julian date [min] i = 0,1,2
      *     Geodetic geodetic[i] - geodetic positions of observer, i = 0,1,2
      *     double dTheta_dTime - Earth rotation
      *     double f - flattening of Earth
      *     double a_e - planet radius (equator)
      *     double mi -
      *     double k  -
      *     boolean isMotionRetrograde - is motion retrograde
      *     double r_1Appr - 1st input approximation of distance r_1 [a_e]
      *     double r_2Appr - 1st input approximation of distance r_2 [a_e]
      *     int lambda - number of revolutions
      *
      * OUTPUT:
      *     Vector array - 3. Position vectors and 2nd velocity vector
      *     vectorArray[0] - 1st position
      *     vectorArray[1] - 2nd position
      *     vectorArray[2] - 3rd position
      *     vectorArray[3] - 2nd velocity
      */
    public Vector[] getPositionVector(double ra[], double dec[], double time[],
                Geodetic geodetic[], double dTheta_dTime, double flattening, double a_e,
                double mi, double k, boolean isMotionRetrograde, double r_1Appr, double r_2Appr){
        this.k = k;
        this.mi = mi;
        this.ra = ra;
        this.dec = dec;
        this.time = time;
        this.geodetic = geodetic;
        this.dTheta_dTime = dTheta_dTime;
        this.flattening = flattening;
        this.a_e = a_e;
        this.isMotionRetrograde = isMotionRetrograde;
        this.r_1Appr = r_1Appr;
        this.r_2Appr = r_2Appr;
        this.lambda = lambda;

        //condition 7.321
        //double epsilon = a_e*epsilonTolerance;//a_e*10e-6;
        //double epsilon = 10e-6;

        //initialization of main variable - info about position and velocity
        StateVector sv = new StateVector();

        // Formula 7.274
        //tau_1 = k * (time[0] - time[1]);

        // Formula 7.275
        //tau_3 = k * (time[2] - time[1]);

        //System.out.println("tau_1 " + tau_1);
        //System.out.println("tau_3 " + tau_3);

        //Greenwich siderial time - 1.27 [s]
        thetaG_2 = getGreenwichSiderialTime(time[1]);
        //System.out.println("thetaG_2 " + Math.toDegrees(thetaG_2) + " deg");

        for(int i = 0; i < 2; i++){
            // Formulas 7.276
            l_Vec[i] = new Vector(3);
            l_Vec[i].v[0] = Math.cos(dec[i])*Math.cos(ra[i]);
            l_Vec[i].v[1] = Math.cos(dec[i])*Math.sin(ra[i]);
            l_Vec[i].v[2] = Math.sin(dec[i]);

            //System.out.println("l_Vec[" + i + "] size " + Vector.getSize(l_Vec[i]));

            //Formulas 7.277
            g_1[i] = a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) +
                    geodetic[i].altitude;
            g_2[i] = (1 - flattening)*(1 - flattening)*a_e/(Math.sqrt(1 - (2*flattening - flattening*flattening)*Math.pow(Math.sin(geodetic[i].lat),2))) +
                    geodetic[i].altitude;

            //System.out.println("g_1[0] " + g_1[0]);
            //System.out.println("g_2[0] " + g_2[0]);

            //Formula 7.278
            theta[i] = thetaG_2 + dTheta_dTime*(time[i] - time[1]) + geodetic[i].lon;
            //System.out.println("theta[i] " + i + " " + theta[i]);

            //Formulas 7.279
            obsR_Vec[i] = new Vector(3);
            obsR_Vec[i].v[0] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.cos(theta[i]);
            obsR_Vec[i].v[1] = -g_1[i] * Math.cos(geodetic[i].lat) * Math.sin(theta[i]);
            obsR_Vec[i].v[2] = -g_2[i] * Math.sin(geodetic[i].lat);

            //System.out.println("observer x " + obsR_Vec[i].v[0] + " e.r.");
            //System.out.println("observer y " + obsR_Vec[i].v[1] + " e.r.");
            //System.out.println("observer z " + obsR_Vec[i].v[2] + " e.r.");

            //formula 7.280
            c_psi[i] = -2 * Vector.getScalarProduct(l_Vec[i], obsR_Vec[i]);
            //System.out.println("c_psi[" + i + "] " + c_psi[i]);
        }

        //before formula 7.281
        r_1g = r_1Appr;
        r_2g = r_2Appr;
        //r_1g = 2.0*a_e;//1.1*a_e;
        //r_2g = 2.1*a_e;//1.11*a_e;
        r_1 = r_1g;
        r_2 = r_2g;
        //System.out.println("r_1 " + r_1 + " e.r.");
        //System.out.println("r_2 " + r_2 + " e.r.");

        //observer position vector size obsR_Vec[i]
        double obsR_VecSize[] = new double[3];

        //loop 7.281 - 7.305/7.311
        double f_1Loop[] = new double[3];
        double f_2Loop[] = new double[3];

        boolean lookForValue = true;

        //before 7.312
        //double delta_r_PerCent = 0.0001;

        //condition 7.321
        //epsilon = a_e*10e-8;//a_e*10e-6;

        //number of iteration
        int iteration = 0;
        //int iterationBoundery = 50;

        //loop 2.281 - 7.322
        //while(lookForValue){
            //for(int l = 0; l < 3; l++){
                //System.out.println("------------- Case " + l + ". -----------------");
                //1st case - get F1(r_1,r_2) and F2(r_1,r_2)
                //if(l == 0){
                    r_1 = r_1g;
                    r_2 = r_2g;
                //}
                //2nd case - get F1{r_1+delta_r_1},r_2), F2(r_1+delta_r_1},r_2)
                //else if(l == 1) {
                //    r_1 = (1.0 + delta_r_PerCent)*r_1g;
                //    r_2 = r_2g;
                //}
                //3rd case - get F1{r_1,r_2+delta_r_1), F2(r_1,r_2+delta_r_1}
                //else if(l == 2) {
                //    r_1 = r_1g;
                //    r_2 = (1.0 + delta_r_PerCent)*r_2g;
                //}

                for(int i = 0; i < 2; i++){
                    obsR_VecSize[i] = Vector.getSize(obsR_Vec[i]);
                    //System.out.println("obsR_VecSize[" + i + "] " + obsR_VecSize[i] + " e.r.");
                    //formula 7.281
                    if(i == 0){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 *
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1))));//r_1 !!!!!!!
                        //System.out.println("r_1 " + r_1);
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 *
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_1*r_1)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    else if(i == 1){
                        ro[i] = 0.5 * (-c_psi[i] + Math.sqrt((c_psi[i]*c_psi[i] - 4 *
                                (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2))));//r_2 !!!!!!!
                        //System.out.println("Beginning sqrt[" + i + "] " + (c_psi[i]*c_psi[i] - 4 *
                        //        (obsR_VecSize[i]*obsR_VecSize[i] - r_2*r_2)) + " ");
                        //System.out.println("ro[" + i + "] " + ro[i] + " e.r.");
                    }
                    //formula 7.282
                    Vector l_ro_Vec = Vector.multiplyVector(l_Vec[i], ro[i]);
                    //r_Vec[i] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[i], ro[i]),
                    r_Vec[i] = Vector.subtractVectors(l_ro_Vec, obsR_Vec[i]);
                    //System.out.println("l_ro_Vec[" + i + "].v[1]    " + l_ro_Vec.v[1]);
                    //System.out.println("obsR_Vec[" + i + "].v[1] " + obsR_Vec[i].v[1]);
                    //System.out.println("r_Vec[" + i + "].v[1]    " + r_Vec[i].v[1]);
                }

                //formulas 7.283
                //is motion retrograde
                int retrograde = 1;
                if(isMotionRetrograde) retrograde = -1;

                r_1 = Vector.getSize(r_Vec[0]);
                //System.out.println("1. r_1 " + r_1);
                //System.out.println("1. r_1 " + Vector.getSize(r_Vec[0]));
                r_2 = Vector.getSize(r_Vec[1]);
                //System.out.println("1. r_2 " + r_2);
                //System.out.println("1. r_2 " + Vector.getSize(r_Vec[1]));

                //w_Vec.v[0] = (r_Vec[0].v[1]*r_Vec[1].v[2] - r_Vec[1].v[1]*r_Vec[0].v[2])/
                //                (r_1*r_2);
                //w_Vec.v[1] = (r_Vec[1].v[0]*r_Vec[0].v[2] - r_Vec[0].v[0]*r_Vec[1].v[2])/
                //                (r_1*r_2);
                //w_Vec.v[2] = (r_Vec[0].v[0]*r_Vec[1].v[1] - r_Vec[1].v[0]*r_Vec[0].v[1])/
                //                (r_1*r_2);

                /*
                w_Vec = Vector.getVectorProduct(r_Vec[0], r_Vec[1]);
                w_Vec = Vector.multiplyVector(w_Vec, 1.0/(r_1*r_2));

                //System.out.println("w_Vec.v[2] " + w_Vec.v[2]);

                //condition 20
                w_Vec = Vector.multiplyVector(w_Vec, retrograde);

                //formula 7.284
                double r3_w = Vector.getScalarProduct(obsR_Vec[2], w_Vec);
                double l3_w = Vector.getScalarProduct(l_Vec[2], w_Vec);
                ro[2] = r3_w/l3_w;

                //formula 7.285
                r_Vec[2] = Vector.subtractVectors(Vector.multiplyVector(l_Vec[2], ro[2]),
                           obsR_Vec[2]);

                //formula 7.286
                r_3 = Vector.getSize(r_Vec[2]);

                System.out.println("r_1 " + r_1 + " e.r.");
                System.out.println("r_2 " + r_2 + " e.r.");
                System.out.println("r_3 " + r_3 + " e.r.");

                // formulas 7.287
                cos_v2_v1 = Vector.getScalarProduct(r_Vec[1], r_Vec[0])/(r_2*r_1);
                //cos_v2_v2 = Vector.getScalarProduct(r_Vec[1], r_Vec[1])/(r_2*r_2);
                cos_v3_v1 = Vector.getScalarProduct(r_Vec[2], r_Vec[0])/(r_3*r_1);
                cos_v3_v2 = Vector.getScalarProduct(r_Vec[2], r_Vec[1])/(r_3*r_2);

                //System.out.println("cos_v2_v1 " + cos_v2_v1);
                //System.out.println("cos_v2_v2 " + cos_v2_v2);
                //System.out.println("cos_v3_v1 " + cos_v3_v1);
                //System.out.println("cos_v3_v2 " + cos_v3_v2);

                //formulas 7.288 - 7.289
                sin_v2_v1 =     (r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v2_v1*cos_v2_v1);
                //sin_v2_v2 = 0;//(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])/
                        //Math.abs(r_Vec[1].v[0] * r_Vec[1].v[1] - r_Vec[1].v[0] * r_Vec[1].v[1])*
                        //Math.sqrt(1 - cos_v2_v2*cos_v2_v2);
                sin_v3_v1 =     (r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])/
                        Math.abs(r_Vec[0].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[0].v[1])*
                        Math.sqrt(1 - cos_v3_v1*cos_v3_v1);
                sin_v3_v2 =     (r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])/
                        Math.abs(r_Vec[1].v[0] * r_Vec[2].v[1] - r_Vec[2].v[0] * r_Vec[1].v[1])*
                        Math.sqrt(1 - cos_v3_v2*cos_v3_v2);

                //System.out.println("sin_v2_v1 " + sin_v2_v1);
                //System.out.println("sin_v2_v2 " + sin_v2_v2);
                //System.out.println("sin_v3_v1 " + sin_v3_v1);
                //System.out.println("sin_v3_v2 " + sin_v3_v2);

                //System.out.println("v2 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1)));
                //System.out.println("v3 - v1 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1)));
                //System.out.println("v3 - v2 " + Math.toDegrees(getAngleFromSinAndCos(sin_v3_v2, cos_v3_v2)));

                //condition before 7.288
                if(w_Vec.v[2] < 0){
                    //System.out.println("w_Vec.v[2] < 0");
                    //this is probably solved by method getAngleFromSinAndCos(sin, cos)
                    sin_v2_v1 = (-1) * sin_v2_v1;
                    //sin_v2_v2 = (-1) * sin_v2_v2;
                    sin_v3_v1 = (-1) * sin_v3_v1;
                    sin_v3_v2 = (-1) * sin_v3_v2;
                }

                //formula and condition 7.290
                //double v3_v1 = Math.acos(cos_v3_v1);
                //double v3_v1 = Math.acos(cos_v3_v1);
                double v3_v1 = getAngleFromSinAndCos(sin_v3_v1, cos_v3_v1);
                //double v3_v1Sin = Math.asin(sin_v3_v1);
                //if(v3_v1Sin < 0) v3_v1 = 2*Math.PI - v3_v1;

                //System.out.println("v3_v1 " + Math.toDegrees(v3_v1) + " deg");
                if(v3_v1 > Math.PI){
                    //System.out.println("v3_v1 > PI !!!!");
                    //formulas 7.290
                    c_1 = (r_2/r_1)*(sin_v3_v2/sin_v3_v1);
                    c_3 = (r_2/r_3)*(sin_v2_v1/sin_v3_v1);
                    //formula 7.291
                    p = (c_1*r_1 + c_3*r_3 - r_2) / (c_1 + c_3 - 1);
                    //System.out.println("(c_1*r_1 + c_3*r_3 - r_2) " + (c_1*r_1 + c_3*r_3 - r_2));
                    //System.out.println("(c_1 + c_3 - 1) " + (c_1 + c_3 - 1));
                }
                //v3 - v1 <= Math.PI
                else{
                    //System.out.println("v3_v1 <= PI !!!!");
                    //formulas 7.292
                    cr_1 = (r_1/r_2)*(sin_v3_v1/sin_v3_v2);
                    cr_3 = (r_1/r_3)*(sin_v2_v1/sin_v3_v2);
                    //System.out.println("(sin_v3_v1/sin_v3_v2) " + (sin_v3_v1/sin_v3_v2));
                    //System.out.println("(sin_v2_v1/sin_v3_v2) " + (sin_v2_v1/sin_v3_v2));
                    //System.out.println("cr_1 " + cr_1 );
                    //System.out.println("cr_3 " + cr_3);
                    //formula 7.293
                    p = (r_1 + cr_3*r_3 - cr_1*r_2) / (1 + cr_3 - cr_1);
                    //System.out.println("(r_1 + cr_3*r_3 - cr_1*r_2) " + (r_1 + cr_3*r_3 - cr_1*r_2));
                    //System.out.println("(1 + cr_3 - cr_1) " + (1 + cr_3 - cr_1));
                }
                System.out.println("----- Beware!!! p " + p);
                //formulas 7.294
                e_cosv_1 = p/r_1 - 1;
                e_cosv_2 = p/r_2 - 1;
                e_cosv_3 = p/r_3 - 1;

                //formula 7.295
                //double v2_v1 = Math.acos(cos_v2_v1);
                double v2_v1 = this.getAngleFromSinAndCos(sin_v2_v1, cos_v2_v1);
                //condition before 7.295
                if(v2_v1 != Math.PI){
                    //formula 7.295
                    e_sinv_1 = (cos_v2_v1*e_cosv_1 - e_cosv_2)/sin_v2_v1;
                    e_sinv_2 = (-cos_v2_v1*e_cosv_2 + e_cosv_1)/sin_v2_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }
                else{
                    //formula 7.296
                    e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                }
                if(v3_v1!= Math.PI){
                    //e_sinv_2 = (cos_v3_v2*e_cosv_2 - e_cosv_3)/sin_v3_v1;
                    //System.out.println("e_sinv_2 " + e_sinv_2);
                    //e_sinv_3 = (-cos_v3_v2*e_cosv_3 + e_cosv_2)/sin_v3_v1;
                }

                //formula 7.297
                //e = Math.sqrt(e_cosv_1*e_cosv_1 + e_sinv_1*e_sinv_1);
                //System.out.println("e_1 " + e);
                e = Math.sqrt(e_cosv_2*e_cosv_2 + e_sinv_2*e_sinv_2);
                //System.out.println("e_2 " + e);
                //e = Math.sqrt(e_cosv_3*e_cosv_3 + e_sinv_3*e_sinv_3);
                //System.out.println("e_3 " + e);
                //e = 0.16419;
                //formula 7.298
                a = p/(1.0 - e*e);
                System.out.println("a " + a);

                //e = 0.00001;
                //a = 63963481/Constants.R_Earth;

                //System.out.println("asin(e_sinv_2/e) " + Math.asin(e_sinv_2/e));
                //System.out.println("acos(e_cosv_2/e) " + Math.acos(e_cosv_2/e));
                //System.out.println("asin + acos " + (Math.asin(e_sinv_2/e)+ Math.acos(e_cosv_2/e)));
                //System.out.println("asin - acos " + (Math.asin(e_sinv_2/e)- Math.acos(e_cosv_2/e)));
                //System.out.println("pi          " + (Math.PI));
                //System.out.println("v_2 " + this.getAngleFromSinAndCos(e_sinv_2/e, e_cosv_2/e));

                //elliptic case
                if(e*e < 1.0){
                    System.out.println("Elliptic case!!! ------!!!--------");
                    //formula 7.299
                    n = k*Math.sqrt(mi)*Math.pow(a,-1.5);
                    System.out.println("n " + n);
                    //formula 7.300
                    s_e = r_2/p * Math.sqrt(1 - e*e) * e_sinv_2;
                    //System.out.println("s_e " + s_e);
                    //formula 7.301
                    c_e = r_2/p * (e*e + e_cosv_2);
                    //System.out.println("c_e " + c_e);
                    //formulas 7.302
                    sin_E3_E2 = r_3*sin_v3_v2/Math.sqrt(a*p) - (r_3/p)*(1 - cos_v3_v2)*s_e;
                    cos_E3_E2 = 1 - (r_3*r_2/(a*p))*(1 - cos_v3_v2);
                    //.out.println("sin_E3_E2 " + sin_E3_E2);
                    //System.out.println("cos_E3_E2 " + cos_E3_E2);
                    //formulas 7.303
                    sin_E2_E1 = r_1*sin_v2_v1/Math.sqrt(a*p) + (r_1/p)*(1 - cos_v2_v1)*s_e;
                    cos_E2_E1 = 1 - (r_2*r_1/(a*p))*(1 - cos_v2_v1);
                    //System.out.println("sin_E2_E1 " + sin_E2_E1);
                    //System.out.println("cos_E2_E1 " + cos_E2_E1);

                    e3_E2 = getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.PI*2 - getAngleFromSinAndCos(sin_E3_E2, cos_E3_E2);
                    //e3_E2 = Math.atan(sin_E3_E2 / cos_E3_E2);
                    e2_E1 = getAngleFromSinAndCos(sin_E2_E1, cos_E2_E1);
                    //e2_E1 = Math.atan(sin_E2_E1 / cos_E2_E1);

                    //System.out.println("e3_E2 " + e3_E2);
                    //System.out.println("e2_E1 " + e2_E1);

                    //formulas 7.304
                    m3_M2 = e3_E2 + 2*s_e*Math.pow(Math.sin(e3_E2/2),2) - c_e*sin_E3_E2;
                    m1_M2 = -e2_E1 + 2*s_e*Math.pow(Math.sin(e2_E1/2),2) + c_e*sin_E2_E1;

                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);
                    //formulas 7.305
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;

                    //System.out.println("k*(m1_M2/n) " + k*(m1_M2/n));
                    //System.out.println("tau_1 " + tau_1);
                    //System.out.println("k*(m3_M2/n) " + k*(m3_M2/n));
                    //System.out.println("tau_3 " + tau_3);
                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //hyperbolic case
                else{
                    //System.out.println("Hyperbolic case!!! ------!!!--------");
                    System.out.println("a " + a);
                    //formula 7.306
                    n = k*Math.sqrt(mi)*Math.pow(-a, 1.5);
                    System.out.println("n " + n);
                    //formulas 7.307
                    s_h = r_2/p * Math.sqrt(e*e - 1)*e_sinv_2;
                    c_h = r_2/p * (e*e + e_cosv_2);

                    //System.out.println("s_h " + s_h);
                    //System.out.println("c_h " + c_h);
                    //formulas 7.308
                    double sinh_F3_F2, sinh_F2_F1;
                    double f3_F2, f2_F1;
                    sinh_F3_F2 = r_3/Math.sqrt(-a*p)*sin_v3_v2 - r_3/p*(1 - cos_v3_v2)*s_h;
                    sinh_F2_F1 = r_1/Math.sqrt(-a*p)*sin_v2_v1 + r_1/p*(1 - cos_v2_v1)*s_h;

                    //System.out.println("sinh_F3_F2 " + sinh_F3_F2);
                    //System.out.println("sinh_F2_F1 " + sinh_F2_F1);

                    //formulas 7.309
                    f3_F2 = Math.log10(sinh_F3_F2 + Math.sqrt(sinh_F3_F2*sinh_F3_F2 + 1));
                    f2_F1 = Math.log10(sinh_F2_F1 + Math.sqrt(sinh_F2_F1*sinh_F2_F1 + 1));

                    //System.out.println("f3_F2 " + f3_F2);
                    //System.out.println("f2_F1 " + f2_F1);

                    //formulas 7.310
                    //m3_M2 = -(f3_F2) + 2*s_h*(Math.sqrt(1 + sinh_F3_F2*sinh_F3_F2) - 1)/2.0 +
                    //        c_h*sinh_F3_F2;
                    //m1_M2 = (f2_F1) + 2*s_h*(Math.sqrt(1 + sinh_F2_F1*sinh_F2_F1) - 1)/2.0 +
                    //        c_h*sinh_F2_F1;
                    m3_M2 = -(f3_F2) + 2*s_h*Math.pow(Math.sinh(f3_F2/2),2) + c_h*Math.sinh(f3_F2);
                    m1_M2 =  (f2_F1) + 2*s_h*Math.pow(Math.sinh(f2_F1/2),2) - c_h*Math.sinh(f2_F1);

                    //System.out.println("m3_M2 " + m3_M2);
                    //System.out.println("m1_M2 " + m1_M2);

                    //formulas 7.311
                    f_1 = tau_1 - k*(m1_M2/n) + k*(2*Math.PI/n)*lambda;;
                    f_2 = tau_3 - k*(m3_M2/n) - k*(2*Math.PI/n)*lambda;;

                    //System.out.println("f_1 " + f_1);
                    //System.out.println("f_2 " + f_2);
                }

                //1st/2nd/3rd case - l (low L)
                f_1Loop[l] = f_1;
                f_2Loop[l] = f_2;
            }

            for(int l = 0; l < 3; l++){
                 //System.out.println("f_1 f_1Loop[" + l + "]" + f_1Loop[l]);
                 //System.out.println("f_2 f_2Loop[" + l + "]" + f_2Loop[l]);
            }

            System.out.println("-------------------- Loop 'for' OUT ------------------");

            //
            //double delta_r_PerCent = 0.04;

            System.out.println("r_1g " + r_1g);
            System.out.println("r_2g " + r_2g);

            //formula 7.312
            dF_1_DIV_dr_1 = (f_1Loop[1] - f_1Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.313
            dF_2_DIV_dr_1 = (f_2Loop[1] - f_2Loop[0])/(delta_r_PerCent*r_1g);
            //formula 7.314
            dF_1_DIV_dr_2 = (f_1Loop[2] - f_1Loop[0])/(delta_r_PerCent*r_2g);
            //formula 7.315
            dF_2_DIV_dr_2 = (f_2Loop[2] - f_2Loop[0])/(delta_r_PerCent*r_2g);

            System.out.println("delta_r_PerCent*r_1g " + delta_r_PerCent*r_1g);
            System.out.println("delta_r_PerCent*r_2g " + delta_r_PerCent*r_2g);

            //System.out.println("dF_1_DIV_dr_1 " + dF_1_DIV_dr_1);
            //System.out.println("dF_2_DIV_dr_1 " + dF_2_DIV_dr_1);
            //System.out.println("dF_1_DIV_dr_2 " + dF_1_DIV_dr_2);
            //System.out.println("dF_2_DIV_dr_2 " + dF_2_DIV_dr_2);

            //formula 7.316
            delta = dF_1_DIV_dr_1*dF_2_DIV_dr_2 - dF_2_DIV_dr_1*dF_1_DIV_dr_2;
            //formula 7.317
            delta_1 = dF_2_DIV_dr_2*f_1Loop[0] - dF_1_DIV_dr_2*f_2Loop[0];
            //System.out.println("dF_2_DIV_dr_2*f_1Loop[0] " + dF_2_DIV_dr_2*f_1Loop[0]);
            //System.out.println("dF_1_DIV_dr_2*f_2Loop[0] " + dF_1_DIV_dr_2*f_2Loop[0]);
            //formula 7.318
            delta_2 = dF_1_DIV_dr_1*f_2Loop[0] - dF_2_DIV_dr_1*f_1Loop[0];
            //System.out.println("dF_1_DIV_dr_1*f_2Loop[0] " + dF_1_DIV_dr_1*f_2Loop[0]);
            //System.out.println("dF_2_DIV_dr_1*f_1Loop[0] " + dF_2_DIV_dr_1*f_1Loop[0]);

            //formula 7.319
            deltar_1 = -delta_1/delta;
            //formula 7.320
            deltar_2 = -delta_2/delta;

            //System.out.println("delta " + delta);
            //System.out.println("delta_1 " + delta_1);
            //System.out.println("delta_2 " + delta_2);

            //System.out.println("deltar_1: " + deltar_1);
            //System.out.println("deltar_2: " + deltar_2);

            r_1 = r_1g + deltar_1;
            //System.out.println("r_1 " + r_1);
            r_2 = r_2g + deltar_2;
            //System.out.println("r_2 " + r_2);

            r_1g = r_1;
            r_2g = r_2;

            System.out.println("iteration " + iteration);

            //loop 7.321 - 7.322
            //condition 7.321
            if((Math.abs(deltar_1) < epsilon)&&(Math.abs(deltar_2) < epsilon)){
                lookForValue = false;
                break;
            }
            //formulas 7.322
            else{
                //r_1 = r_1g + deltar_1;
                //System.out.println("r_2 " + r_2);
                //r_2 = r_2g + deltar_2;
                //System.out.println("r_2 " + r_2);

                //r_1g = r_1;
                //r_2g = r_2;
            }

            System.out.println("r_1: " + r_1);
            System.out.println("r_2: " + r_2);

            //System.out.println(zmazat);
            if(p < 0) {
                System.out.println("---------------------------- p is negative!!!!");
                break;
            }
            if((r_1 < a_e)||(r_2 < a_e)) {
                System.out.println("---------------------------- r is too small!!!!");
                break;
            }

            iteration++;
            //System.out.println("iteration " + iteration);
            if(iteration > iterationBoundery) lookForValue = false;
        }

        //formula 7.323
        fSerie = 1.0 - a/r_2*(1 - cos_E3_E2);
        //formula 7.324
        gSerie = tau_3 - Math.pow(a, 1.5)/Math.sqrt(mi)*(e3_E2 - sin_E3_E2);
        //System.out.println("FSerie 1 " + fSerie);
        //System.out.println("GSerie 1 " + gSerie);
        //System.out.println("Tau 1" + tau_3);

        //formula 7.325
        r_2_Dot = Vector.subtractVectors(r_Vec[2], Vector.multiplyVector(r_Vec[1], fSerie));
        r_2_Dot = Vector.multiplyVector(r_2_Dot, 1.0/gSerie);

        //State vector
        sv.r = r_Vec[1];
        sv.v = r_2_Dot;

        System.out.println("Computed a = " + a);
        System.out.println("Computed e = " + e);
        */
        //vector array with vectors of 3 positions and 1 velocity
        Vector vectorArray[] = new Vector[2];
        vectorArray[0] = r_Vec[0];
        vectorArray[1] = r_Vec[1];
        //vectorArray[2] = r_Vec[2];
        //vectorArray[3] = r_2_Dot;

        //System.out.println("1st vector size " + vectorArray[0].getSize(vectorArray[0])/1000 + " km");
        //System.out.println("2nd vector size " + vectorArray[1].getSize(vectorArray[1])/1000 + " km");

        //System.out.println("testLala2 Vel [x] " + vectorArray[3].v[0]);
        //System.out.println("testLala2 Vel [y] " + vectorArray[3].v[1]);
        //System.out.println("testLala2 Vel [z] " + vectorArray[3].v[2]);
        //System.out.println("vel_2 " + Vector.getSize(vectorArray[3]));

        return vectorArray;
    }
    
    /**
     * main test method
     * 
     * Escobal excercise - 7.7 - reference orbit 7.7.2
     */
    public static void main(String args[]){
        //constants
        double k2 = 0.07436574;  //[(e.r)^3/2 / min]
        double mi2 = 1.0;        //[e.m.] - Earth mass
        double a_e2 = 1.0;       //[e.r.] - Earth radius
        
        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime2 = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening2 = Constants.f_Earth;
        
        /** Orbit 7.6.6 **/
        
        double time2[] = new double[3];  //[day]
        double ra2[] = new double[3];    //[deg]
        double dec2[] = new double[3];   //[deg]
        Geodetic geodetic2[] = new Geodetic[3];
        
        /*
        //initialization of variables
        time2[0] = (Time.getMjd(new Time(1959,3,2,18,4,57.67)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(1959,3,2,19,19,13.52)) + 2400000.5)*1440;    //[min]
        //System.out.println("time2[1] " + (Time.getMjd(new Time(1959,3,2,19,19,13.52)) + 2400000.5));        
        time2[2] = (Time.getMjd(new Time(1959,3,2,19,21,8.99)) + 2400000.5)*1440;    //[min]
        
        ra2[0] = Math.toRadians(223.12920833333333333333333333333);  //[rad]
        ra2[1] = Math.toRadians(86.880666666666666666666666666667);
        ra2[2] = Math.toRadians(98.893375);
        
        dec2[0] = Math.toRadians(23.999861111111111111111111111111);
        dec2[1] = Math.toRadians(2.1368888888888888888888888888889);
        dec2[2] = Math.toRadians(7.925);
        
        //Tokyo
        geodetic2[0] = new Geodetic(Math.toRadians(139.53525), Math.toRadians(35.67322), 8.153e-6);
        for(int i = 1; i <3; i++){
            //Olifantsfontein, altitude is in [e.r.]
            geodetic2[i] = new Geodetic(Math.toRadians(28.247528), Math.toRadians(-25.959639), 2.4e-4);
        }
        
        StateVector sv = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.1, a_e2*1.11, 0);
        
        //results
        //System.out.println("x =         " + sv.r.v[0]);
        //System.out.println("y =         " + sv.r.v[1]);
        //System.out.println("z =         " + sv.r.v[2]);
        
        //System.out.println("vx =         " + sv.v.v[0]);
        //System.out.println("vy =         " + sv.v.v[1]);
        //System.out.println("vz =         " + sv.v.v[2]);
        System.out.println("Escobal  a = " + 1.3023); 
        System.out.println("Escobal  e = " + 0.16419);
        
        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv.r,sv.v,time2[1], mi2);
        
        System.out.println("Escobal  a: 1.3023 e.r.");
        System.out.println("Computed a: " + kepler2.a + " e.r.");
        System.out.println("Escobal  e: 0.16419");
        System.out.println("Computed e: "+kepler2.e);
        System.out.println("Escobal  i: 32.878 °");
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+" °");
        System.out.println("Escobal  Omega: 136.53 °");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+" °");
        System.out.println("Escobal  omega: 203.95 °");
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+" °");
        System.out.println("Escobal  Ma: ??? °");
        System.out.println("Computed Ma: "+Math.toDegrees(kepler2.M)+" °\n");
        */
        
        /** Orbit no. VII. - 7.7.2 **/
        
        //double time2[] = new double[3];  //[day]
        //double ra2[] = new double[3];    //[deg]
        //double dec2[] = new double[3];   //[deg]
        //Geodetic geodetic2[] = new Geodetic[3];
        /*
        //initialization of variables
        time2[0] = 2438314.7916667*1440;    //[min]
        time2[1] = 2438314.8055556*1440;    //[min]
        time2[2] = 2438314.8194444*1440;    //[min]
        
        ra2[0] = Math.toRadians(153.6949);  //[rad]
        ra2[1] = Math.toRadians(186.7296);
        ra2[2] = Math.toRadians(202.9475);
        
        dec2[0] = Math.toRadians(36.2726);
        dec2[1] = Math.toRadians(24.8918);
        dec2[2] = Math.toRadians(9.8723);
        
        for(int i = 0; i <3; i++){
            //altitude is in [e.r.]
            geodetic2[i] = new Geodetic(Math.toRadians(250.0), Math.toRadians(40.0), 0.78393e-3);
        }
        
        
        StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*2.0, a_e2*2.1, 0);
        
        //results
        System.out.println("        x = " + sv2.r.v[0]);
        System.out.println("Escobal_x = -1.22192");
        System.out.println("        y = " + sv2.r.v[1]);
        System.out.println("Escobal_y = 0.0352894");
        System.out.println("        z = " + sv2.r.v[2]);
        System.out.println("Escobal_z = 1.54752");
        
        System.out.println("        vx = " + sv2.v.v[0]);
        System.out.println("Escobal_vx = -0.468449");
        System.out.println("        vy = " + sv2.v.v[1]);
        System.out.println("Escobal_vy = -0.513595");
        System.out.println("        vz = " + sv2.v.v[2]);
        System.out.println("Escobal_vz = -0.175849\n\n");
        */
        /*
        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv2.r,sv2.v,time2[1], mi2);
        
        System.out.println("a: "+kepler2.a + " e.r.");
        System.out.println("e: "+kepler2.e);
        System.out.println("i: "+Math.toDegrees(kepler2.incl)+" �");
        System.out.println("Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        System.out.println("omega: "+Math.toDegrees(kepler2.omega)+" �");
        System.out.println("Ma: "+Math.toDegrees(kepler2.M)+" �\n");
        */
        /** Reference orbit no. VIII. - 7.7.3 **/
        
        //double time2[] = new double[3];  //[day]
        //double ra2[] = new double[3];    //[deg]
        //double dec2[] = new double[3];   //[deg]
        //Geodetic geodetic2[] = new Geodetic[3];
        /*
        //initialization of variables
        time2[0] = 2438181.9583333*1440;    //[min]
        time2[1] = 2438182.1666667*1440;    //[min]
        time2[2] = 2438182.3750000*1440;    //[min]
        
        ra2[0] = Math.toRadians(96.7675);  //[rad]
        ra2[1] = Math.toRadians(182.5533);
        ra2[2] = Math.toRadians(215.0986);
        
        dec2[0] = Math.toRadians(34.1759);
        dec2[1] = Math.toRadians(35.2127);
        dec2[2] = Math.toRadians(8.4239);
        
        for(int i = 0; i <3; i++){
            //altitude is in [e.r.]
            geodetic2[i] = new Geodetic(Math.toRadians(353.7935), Math.toRadians(36.4594), 0.12386e-4);
        }
        
        
        StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        //sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*10, a_e2*10.1, 0);
        
        //results
        System.out.println("        x = " + sv2.r.v[0]);
        System.out.println("Escobal_x = -6.82944");
        System.out.println("        y = " + sv2.r.v[1]);
        System.out.println("Escobal_y = 0.388821");
        System.out.println("        z = " + sv2.r.v[2]);
        System.out.println("Escobal_z = 5.10408");
        
        System.out.println("        vx = " + sv2.v.v[0]);
        System.out.println("Escobal_vx = -0.194685");
        System.out.println("        vy = " + sv2.v.v[1]);
        System.out.println("Escobal_vy = -0.301251");
        System.out.println("        vz = " + sv2.v.v[2]);
        System.out.println("Escobal_vz = -0.0753532\n\n");
        */
        
        //test
       // double angle = Math.toRadians(320);
       // System.out.println(Math.toDegrees(new RIterationAnglesOnly().getAngleFromSinAndCos(Math.sin(angle), Math.cos(angle))));
        
        /**
         * Excercise 6 (s. 291)
         */
        
        //initialization of variables
        time2[0] = (Time.getMjd(new Time(1959,9,26,21,25,37.403)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(1959,9,26,21,26,37.862)) + 2400000.5)*1440;    //[min]
        //System.out.println("time2[1] " + time2[1]/1440);
        time2[2] = (Time.getMjd(new Time(1959,9,26,21,27,45.919)) + 2400000.5)*1440;    //[min]
        
        ra2[0] = Math.toRadians(254.673208333);  //[rad]
        ra2[1] = Math.toRadians(260.941125);
        ra2[2] = Math.toRadians(269.85375);
        
        dec2[0] = Math.toRadians(13.0927222);
        dec2[1] = Math.toRadians(13.28569444);
        dec2[2] = Math.toRadians(13.00425);
        
        //San Fernando, Spain
        for(int i = 0; i <3; i++){
            geodetic2[i] = new Geodetic(Math.toRadians(353.79486111), Math.toRadians(36.4638333), 3.76e-6);
        }
        
        StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        //sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.5, a_e2*1.51, 0);
        
        System.out.println("Escobal  a = " + 1.4400); 
        System.out.println("Escobal  e = " + 0.23167);
        
        System.out.println("        x = " + sv2.r.v[0]);
        System.out.println("        y = " + sv2.r.v[1]);
        System.out.println("        z = " + sv2.r.v[2]);
        
        System.out.println("        vx = " + sv2.v.v[0]);
        System.out.println("        vy = " + sv2.v.v[1]);
        System.out.println("        vz = " + sv2.v.v[2] + "\n");
        
        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv2.r,sv2.v,time2[1], mi2);
        
        System.out.println("Escobal  a: 1.44 e.r.");
        System.out.println("Computed a: "+kepler2.a + " e.r.");
        System.out.println("Computed a: "+kepler2.a*6378.15 + " km");
        System.out.println("Escobal  e: 0.23167");
        System.out.println("Computed e: "+kepler2.e);
        System.out.println("Escobal  i: 33.281 �");
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+" �");
        System.out.println("Escobal  Omega: 205.110 �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        System.out.println("Escobal  omega: 161.790 �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+" �");
        System.out.println("Escobal  Ma: ??? ");
        System.out.println("Computed Ma: "+Math.toDegrees(kepler2.M)+" �\n");

        Vector vectors[] = new Vector[4];
        vectors = new RIterationAnglesOnly().getStateVector_2(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.5, a_e2*1.51, 0);
        //TEST method getVelVecFrom3PosVec()
        System.out.println("v_x " + vectors[3].v[0]*25936*0.3048/1000 + " km/s");
        System.out.println("v_y " + vectors[3].v[1]*25936*0.3048/1000 + " km/s");
        System.out.println("v_z " + vectors[3].v[2]*25936*0.3048/1000 + " km/s");

        Vector testVelocity = new Vector(3);
        Vector vectors2[] = new Vector[3];
        for(int i = 0; i<3; i++){
            vectors2[i] = vectors[i];
            time2[i] = time2[i]/1440;
        }
        testVelocity = new RIterationAnglesOnly().getVelVecFrom3PosVec(vectors2, time2, k2, mi2, false);
        System.out.println("v_x " + testVelocity.v[0]*25936*0.3048/1000 + " km/s");
        System.out.println("v_y " + testVelocity.v[1]*25936*0.3048/1000 + " km/s");
        System.out.println("v_z " + testVelocity.v[2]*25936*0.3048/1000 + " km/s");
        /*
        //GST test
        //double jdGstTest = (Time.getMjd(new Time(1959,9,26,21,25,37.403)) + 2400000.5)*1440;    //[min]
        Time timeTest = new Time(2009,6,8,0,0,0.0);
        double jdGstTest = (Time.getMjd(new Time(2009,6,8,0,0,0.0)) + 2400000.5)*1440;    //[min]
        System.out.println("TimeMjd " + (Time.getMjd(timeTest)));
        jdGstTest = (Time.getMjd(2009,6,8,0,0,0.0) + 2400000.5)*1440;    //[min]
        System.out.println("TimeMjd " + (Time.getMjd(timeTest)));
        System.out.println("TimeJd " + (Time.getMjd(new Time(2009,6,8,0,0,0.0)) + 2400000.5));
        double gstTest = new RIterationAnglesOnly().getGreenwichSiderialTime(jdGstTest);
        while(gstTest > Math.PI*2) gstTest = gstTest - 2*Math.PI;
        System.out.println("GST " + Math.toDegrees(gstTest));
        */
        
        
        //double jdGstTest = (Time.getMjd(new Time(1962,10,12,10,15,30.0)) + 2400000.5)*1440;    //[min]
        //double gstTest = new RIterationAnglesOnly().getGreenwichSiderialTime(jdGstTest);
        //while(gstTest > Math.PI*2) gstTest = gstTest - 2*Math.PI;
        //System.out.println("jdGstTest " + jdGstTest/1440);
        //System.out.println("GST  " + Math.toDegrees(gstTest));
        //GST vs GMST
        //double mjdGmstTest = Time.getGMST(Time.getMjd(new Time(1962,10,12,10,15,30.0)));
        //System.out.println("GMST " + Math.toDegrees(mjdGmstTest));
        
        
        /*
        double jdGstTest = (Time.getMjd(new Time(1962,10,12,0,0,0.0)) + 2400000.5);    //[min]
        ///double gstTest = new RIterationAnglesOnly().getGreenwichSiderialTime(jdGstTest);
        double gstTest = new RIterationAnglesOnly().getGreenwichSiderialTime(jdGstTest, 615.5);
        while(gstTest > Math.PI*2) gstTest = gstTest - 2*Math.PI;
        //gstTest = gstTest + Math.toRadians(615.5*0.25068447);
        System.out.println("jdGstTest " + jdGstTest);
        System.out.println("GST " + Math.toDegrees(gstTest));
        */
        
        //Orbit VII.
        /*
        Vector positionVecTest = new Vector(3);
        positionVecTest.v[0] = -1.22192;
        positionVecTest.v[1] = 0.0352894;
        positionVecTest.v[2] = 1.54752;
                
        Vector roVecTest = Vector.addVectors(positionVecTest, obsR_Vec[1]);
        //roVecTest = Vector.multiplyVector(roVecTest, 1.0/Vector.getSize(roVecTest));
        System.out.println("obsR_Vec[1].v[1] " + obsR_Vec[1].v[1]);
        //double raTest, decTest;
        Vector eqaVecTest = Transformation.getEquatSphericalCoordinates2(roVecTest);
        while(thetaG_2 > Math.PI*2) thetaG_2 = thetaG_2 - 2*Math.PI;
        System.out.println("thetaG_2 " + Math.toDegrees(thetaG_2));

        System.out.println("ra_2In " + Math.toDegrees(ra2[1]));
        System.out.println("dec_2In " + Math.toDegrees(dec2[1]));
        System.out.println("ra_2 " + Math.toDegrees(eqaVecTest.v[0]));
        System.out.println("dec_2 " + Math.toDegrees(eqaVecTest.v[1]));
        System.out.println("ro_2 " + eqaVecTest.v[2]);
        */
        
        //Orbit VIII.
        /*
        Vector positionVecTest = new Vector(3);
        positionVecTest.v[0] = -6.82944;
        positionVecTest.v[1] = 0.388821;
        positionVecTest.v[2] = 5.10408;
                
        Vector roVecTest = Vector.addVectors(positionVecTest, obsR_Vec[1]);
        //roVecTest = Vector.multiplyVector(roVecTest, 1.0/Vector.getSize(roVecTest));
        System.out.println("obsR_Vec[1].v[1] " + obsR_Vec[1].v[1]);
        //double raTest, decTest;
        Vector eqaVecTest = Transformation.getEquatSphericalCoordinates2(roVecTest);
        while(thetaG_2 > Math.PI*2) thetaG_2 = thetaG_2 - 2*Math.PI;
        System.out.println("thetaG_2 " + Math.toDegrees(thetaG_2));

        System.out.println("ra_2In " + Math.toDegrees(ra2[1]));
        System.out.println("dec_2In " + Math.toDegrees(dec2[1]));
        System.out.println("ra_2 " + Math.toDegrees(eqaVecTest.v[0]));
        System.out.println("dec_2 " + Math.toDegrees(eqaVecTest.v[1]));
        System.out.println("ro_2 " + eqaVecTest.v[2]);
        */
        
        
        /**
         * Test from 09/06/2009 - Escobal vs Bucerius         
         */
        /*
        //initialization of variables
        time2[0] = (Time.getMjd(new Time(2008,3,29,18,40,0.0)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(2008,3,29,18,50,0.0)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(2008,3,29,18,55,0.0)) + 2400000.5)*1440;    //[min]
        
        ra2[0] = Math.toRadians(194.66032);  //[rad]
        ra2[1] = Math.toRadians(255.2985);
        ra2[1] = Math.toRadians(320.92274);
        
        dec2[0] = Math.toRadians(-7.40676);
        dec2[1] = Math.toRadians(62.32506);
        dec2[2] = Math.toRadians(50.84109);
        
        //AGO Modra, Slovakia
        for(int i = 0; i <3; i++){
            geodetic2[i] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733), 8.32e-5);
        }
        
        //StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.2, a_e2*1.21, 0);
        
        System.out.println("TLE      a = " + 1.2312); 
        System.out.println("TLE      e = " + 0.0023);
        
        System.out.println("        x = " + sv2.r.v[0]);
        System.out.println("        y = " + sv2.r.v[1]);
        System.out.println("        z = " + sv2.r.v[2]);
        
        System.out.println("        vx = " + sv2.v.v[0]);
        System.out.println("        vy = " + sv2.v.v[1]);
        System.out.println("        vz = " + sv2.v.v[2] + "\n");
        
        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv2.r,sv2.v,time2[1], mi2);
        
        System.out.println("TLE      a: 1.2312 e.r.");
        System.out.println("Computed a: "+kepler2.a + " e.r.");
        System.out.println("TLE      e: 0.0023");
        System.out.println("Computed e: "+kepler2.e);
        System.out.println("TLE      i: 101.88 �");
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+" �");
        System.out.println("TLE      Omega: 133.51 �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        System.out.println("TLE      omega: 165.34 �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+" �");
        System.out.println("TLE/Kpler Ma: 265.256 �");
        System.out.println("Computed  Ma: "+Math.toDegrees(kepler2.M)+" �\n");
        */
        
        /**
         * Identification of VFMO090316, observer S.Gajdos, AGO Modra
         */
        /*
        //initialization of variables
        time2[0] = (Time.getMjd(new Time(2009,3,16,22,10,55.73)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(2009,3,16,22,11,10.73)) + 2400000.5)*1440;    //[min]
        //time2[1] = (Time.getMjd(new Time(2009,3,16,22,11,47.58)) + 2400000.5)*1440;    //[min]
        time2[2] = (Time.getMjd(new Time(2009,3,16,22,13,46.33)) + 2400000.5)*1440;    //[min]
        
        ra2[0] = Math.toRadians(184.565125);  //[rad]
        ra2[1] = Math.toRadians(184.58958333333333333333333333333);
        //ra2[1] = Math.toRadians(184.658125);
        ra2[2] = Math.toRadians(184.86470833333333333333333333333);
        
        dec2[0] = Math.toRadians(17.120722222222222222222222222222);
        dec2[1] = Math.toRadians(17.127833333333333333333333333333);
        //dec2[1] = Math.toRadians(17.146222222222222222222222222222);
        dec2[2] = Math.toRadians(17.20125);
        
        //San Fernando, Spain
        for(int i = 0; i <3; i++){
            geodetic2[i] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733), 8.33e-5);
        }
        
        //StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.5, a_e2*1.51, 0);
        
        //System.out.println("Escobal  a = " + 1.4400); 
        //System.out.println("Escobal  e = " + 0.23167);
        
        System.out.println("        x = " + sv2.r.v[0]);
        System.out.println("        y = " + sv2.r.v[1]);
        System.out.println("        z = " + sv2.r.v[2]);
        
        System.out.println("        vx = " + sv2.v.v[0]);
        System.out.println("        vy = " + sv2.v.v[1]);
        System.out.println("        vz = " + sv2.v.v[2] + "\n");
        
        //Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv2.r,sv2.v,time2[1], mi2);
        
        //System.out.println("Escobal  a: 1.44 e.r.");
        System.out.println("Computed a: "+kepler2.a + " e.r.");
        System.out.println("Computed a: "+kepler2.a*6378.15 + " km");
        //System.out.println("Escobal  e: 0.23167");
        System.out.println("Computed e: "+kepler2.e);
        //System.out.println("Escobal  i: 33.281 �");
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+" �");
        //System.out.println("Escobal  Omega: 205.110 �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        //System.out.println("Escobal  omega: 161.790 �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+" �");
        //System.out.println("Escobal  Ma: ??? ");
        System.out.println("Computed Ma: "+Math.toDegrees(kepler2.M)+" �\n");
        */

        /*
        char degChar = 176;
        String stringDeg = degChar + "";
        
        // Tv camera, test Iridium 33- deb
        //initialization of variables
        //positions
        int amount = 5;
        Time time[]=new Time[amount];
        double meanAnom[] = new double[amount];

        ObservationEq obser[] = new ObservationEq[amount];

        Kepler kepler = new Kepler();

        //combinations of observations
        int a1,a2,a3;

        //time2[0] = (Time.getMjd(new Time(1959,3,2,18,4,57.67)) + 2400000.5)*1440;    //[min]
        //time2[1] = (Time.getMjd(new Time(1959,3,2,19,19,13.52)) + 2400000.5)*1440;    //[min]
        //time2[2] = (Time.getMjd(new Time(1959,3,2,19,21,8.99)) + 2400000.5)*1440;    //[min]

        //ra2[0] = Math.toRadians(223.12920833333333333333333333333);  //[rad]
        //ra2[1] = Math.toRadians(86.880666666666666666666666666667);
        //ra2[2] = Math.toRadians(98.893375);

        //dec2[0] = Math.toRadians(23.999861111111111111111111111111);
        //dec2[1] = Math.toRadians(2.1368888888888888888888888888889);
        //dec2[2] = Math.toRadians(7.925);

        Time t[] = new Time[amount];
        double ra[] = new double[amount];
        double dec[] = new double[amount];

        //1st position
        //t[0] = new Time(2009, 06, 17,21, 6, 0.72);
        //ra[0] = Math.toRadians(271.0440);  //[rad]
        //dec[0] = Math.toRadians(40.9084);
        //meanAnom[0] = 202.618;    //mean anomaly from SatEph

        //2nd position
        //t[1] = new Time(2009, 06, 17,21, 6, 17.48);
        //ra[1] = Math.toRadians(268.1360);  //[rad]
        //dec[1] = Math.toRadians(33.3225);
        //meanAnom[1] = 203.62;    //mean anomaly from SatEph

        //3rd position
        //t[2] = new Time(2009, 06, 17,21, 6, 25.86);
        //ra[2] = Math.toRadians(266.9550);  //[rad]
        //dec[2] = Math.toRadians(29.5086);
        //meanAnom[2] = 204.12;    //mean anomaly from SatEph

        //4th position
        //t[3] = new Time(2009, 06, 17,21, 6, 25.86);
        //ra[3] = Math.toRadians(265.9426);  //[rad]
        //dec[3] = Math.toRadians(25.7733);
        //meanAnom[3] = 204.623;    //mean anomaly from SatEph

        //5th position
        //t[4] = new Time(2009, 06, 17,21, 6, 59.5);
        //ra[4] = Math.toRadians(263.5278);  //[rad]
        //dec[4] = Math.toRadians(15.3045);
        //meanAnom[4] = 206.131;    //mean anomaly from SatEph

        //1st position
        t[0] = new Time(2009, 06, 17,21, 6, 0.76);
        ra[0] = Math.toRadians(271.0279);  //[rad]
        dec[0] = Math.toRadians(40.9783);
        meanAnom[0] = 202.62;    //mean anomaly from SatEph

        //2nd position
        t[1] = new Time(2009, 06, 17,21, 6, 25.91);
        ra[1] = Math.toRadians(266.9440);  //[rad]
        dec[1] = Math.toRadians(29.3831);
        meanAnom[1] = 204.123;    //mean anomaly from SatEph

        //3rd position
        t[2] = new Time(2009, 06, 17,21, 6, 59.5);
        ra[2] = Math.toRadians(263.5278);  //[rad]
        dec[2] = Math.toRadians(15.3045);
        meanAnom[2] = 206.131;    //mean anomaly from SatEph

        double ra3[] = new double[3];
        double dec3[] = new double[3];
        double timeJd[] = new double[3];

        a1 = 0;
        a2 = 1;
        a3 = 2;
                
        for(int i = 0; i < 3; i++){
            int a = 0;
            if(i == 0) a = a1;
            else if(i == 1) a = a2;
            else if(i == 2) a = a3;

            ra3[i] = ra[a];
            dec3[i] = dec[a];
            timeJd[i] = (Time.getMjd(t[a]) + 2400000.5)*1440;    //[min]
        }

        //Ago Modra
        for(int i = 0; i<3;i++){
            geodetic2[i] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733),531.1/6378150);
        }
               
        StateVector sv = new RIterationAnglesOnly().getStateVector(ra3, dec3, timeJd, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.1, a_e2*1.11, 0);

        //results
        //System.out.println("x =         " + sv.r.v[0]);
        //System.out.println("y =         " + sv.r.v[1]);
        //System.out.println("z =         " + sv.r.v[2]);

        //System.out.println("vx =         " + sv.v.v[0]);
        //System.out.println("vy =         " + sv.v.v[1]);
        //System.out.println("vz =         " + sv.v.v[2]);
        //System.out.println("Escobal  a = " + 1.3023);
        //System.out.println("Escobal  e = " + 0.16419);

        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv.r,sv.v,time2[1], mi2);

        System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");
        
        System.out.println("Computed a: " + kepler2.a + " e.r.");
        System.out.println("Computed a: " + kepler2.a*6378.15 + " km");
        System.out.println("Computed e: "+kepler2.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+ " " + stringDeg);
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+ " " + stringDeg);
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+ " " + stringDeg);
        System.out.println("Computed Ma   : "+Math.toDegrees(kepler2.M)+ " " + stringDeg + "\n");
        
        /*
        // Tv camera, orbit determination Norad 25723 - R/B
        //1 25723U 99022C   09256.64718178 -.00000716  00000-0 -20684-4 0  3832
        //2 25723 048.4407 324.3833 0021806 174.1355 185.9833 15.16115652572720

        //initialization of variables
        //positions
        int amount = 7;
        Time time[]=new Time[amount];
        double meanAnom[] = new double[amount];

        ObservationEq obser[] = new ObservationEq[amount];

        Kepler kepler = new Kepler();

        //combinations of observations
        int a1,a2,a3;

        Time t[] = new Time[amount];
        double ra[] = new double[amount];
        double dec[] = new double[amount];
        
        //1st position - AGO
        t[0] = new Time(2009,9,16,2,28,17.446);
        ra[0] = Math.toRadians(46.022);  //[rad]
        dec[0] = Math.toRadians(48.442);
        meanAnom[0] = 269.762;    //mean anomaly from SatEph

        //2nd position - AGO
        t[1] = new Time(2009,9,16,2,28,21.369);
        ra[1] = Math.toRadians(50.747);  //[rad]
        dec[1] = Math.toRadians(48.212);
        meanAnom[1] = 270.01;    //mean anomaly from SatEph

        //3rd position - AGO
        t[2] = new Time(2009,9,16,2,28,24.370);
        ra[2] = Math.toRadians(54.453);  //[rad]
        dec[2] = Math.toRadians(47.644);
        meanAnom[2] = 270.199;    //mean anomaly from SatEph

        //1st position - AGO - other version of positions
        //t[0] = new Time(2009,9,16,2,28,17.446);
        //ra[0] = Math.toRadians(46.015);  //[rad]
        //dec[0] = Math.toRadians(48.648);
        //meanAnom[0] = 269.762;    //mean anomaly from SatEph

        //2nd position - AGO - other version of positions
        //t[1] = new Time(2009,9,16,2,28,24.370);
        //ra[1] = Math.toRadians(54.479);  //[rad]
        //dec[1] = Math.toRadians(47.849);
        //meanAnom[1] = 270.199;    //mean anomaly from SatEph
        
        //1st position - Arbo
        t[3] = new Time(2009,9,16,2,28,17.058);
        ra[3] = Math.toRadians(32.747);  //[rad]
        dec[3] = Math.toRadians(47.905);
        meanAnom[3] = 269.737;    //mean anomaly from SatEph

        //2nd position Arbo
        t[4] = new Time(2009,9,16,2,28,20.902);
        ra[4] = Math.toRadians(37.330);  //[rad]
        dec[4] = Math.toRadians(48.043);
        meanAnom[4] = 269.98;    //mean anomaly from SatEph

        //3rd position Arbo
        t[5] = new Time(2009,9,16,2,28,23.699);
        ra[5] = Math.toRadians(40.485);  //[rad]
        dec[5] = Math.toRadians(48.036);
        meanAnom[5] = 270.157;    //mean anomaly from SatEph

        //1st position Ago - SatEph positions
        //t[4] = new Time(2009,9,16,2,28,17.446);
        //ra[4] = Math.toRadians(32.83895);  //[rad]
        //dec[4] = Math.toRadians(47.63012);
        //meanAnom[4] = 270.157;    //mean anomaly from SatEph

        //2nd position Ago - SatEph positions
        //t[5] = new Time(2009,9,16,2,28,24.370);
        //ra[5] = Math.toRadians(40.88983);  //[rad]
        //dec[5] = Math.toRadians(47.95598);
        //meanAnom[5] = 270.199;    //mean anomaly from SatEph

        //2nd position Arbo - SatEph positions
        //t[6] = new Time(2009,9,16,2,28,23.699);
        //ra[6] = Math.toRadians(40.0974);  //[rad]
        //dec[6] = Math.toRadians(47.949);
        //meanAnom[6] = 270.157;    //mean anomaly from SatEph
        
        double ra3[] = new double[3];
        double dec3[] = new double[3];
        double timeJd[] = new double[3];

        a1 = 1;
        a2 = 2;
        a3 = 5;

        for(int i = 0; i < 3; i++){
            int a = 0;
            if(i == 0) a = a1;
            else if(i == 1) a = a2;
            else if(i == 2) a = a3;

            ra3[i] = ra[a];
            dec3[i] = dec[a];
            timeJd[i] = (Time.getMjd(t[a]) + 2400000.5)*1440;    //[min]
        }

        //Ago Modra
        geodetic2[0] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733),531.1/6378150);
        geodetic2[1] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733),531.1/6378150);
        //geodetic2[2] = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733),531.1/6378150);
        //Arbo
        //geodetic2[0] = new Geodetic(Math.toRadians(18.3685), Math.toRadians(48.3235),185.0/6378150);
        //geodetic2[1] = new Geodetic(Math.toRadians(18.3685), Math.toRadians(48.3235),185.0/6378150);
        geodetic2[2] = new Geodetic(Math.toRadians(18.3685), Math.toRadians(48.3235),185.0/6378150);

        StateVector sv = new RIterationAnglesOnly().getStateVector(ra3, dec3, timeJd, geodetic2, dTheta_dTime2,
                        flattening2, a_e2, mi2, k2, false, a_e2*1.1, a_e2*1.11, 0);

        //results
        //System.out.println("x =         " + sv.r.v[0]);
        //System.out.println("y =         " + sv.r.v[1]);
        //System.out.println("z =         " + sv.r.v[2]);

        //System.out.println("vx =         " + sv.v.v[0]);
        //System.out.println("vy =         " + sv.v.v[1]);
        //System.out.println("vz =         " + sv.v.v[2]);
        //System.out.println("Escobal  a = " + 1.3023);
        //System.out.println("Escobal  e = " + 0.16419);

        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(sv.r,sv.v,time2[1], mi2);

        System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

        System.out.println("Computed a: " + kepler2.a + " e.r.");
        System.out.println("Computed a: " + kepler2.a*6378.15 + " km");
        System.out.println("Computed e: "+kepler2.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler2.incl)+ " " + stringDeg);
        System.out.println("Computed Omega: "+Math.toDegrees(kepler2.Omega)+ " " + stringDeg);
        System.out.println("Computed omega: "+Math.toDegrees(kepler2.omega)+ " " + stringDeg);
        System.out.println("Computed Ma   : "+Math.toDegrees(kepler2.M)+ " " + stringDeg + "\n");
        */
    }
}
 