/*
 * Class is for computing the geocentric position vectors of body with only
 * 3 topocentric positions on the celestial sphere - Az_i,h_i, or R.A._i, dec_i, for i = 1,2,3.
 * Laplace method is used. All formulas are from book P.Escobal - Methods of orbit
 * determination (1976). All formulas numbers in this code are identical with
 * formulas numbers used in Escobal's book.
 *
 * Author of Java source is Jiri Silha - 2009
 */

package OrbitDetermination;

import jahuwaldt.tools.math.*;
import org.netlib.math.complex.*;
import Compute.Silha.*;

/**
 *
 * @author Jiri Silha - 14/04/2009
 */
public class LaplaceAnglesOnly {
    
    //VARIABLES
    
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
    
    //Local variables
    
    /**
     * Formula 7.196
     */
    static double tau_1;
    
    /**
     * Formula 7.197
     */
    static double tau_3;
    
    /**
     * Formulas 7.198
     */
    static double s[] = new double[6];
    
    /**
     * Formulas 7.199
     */
    static Vector l_Vec[] = new Vector[3];
    
    /**
     * Formula 7.200
     */
    static Vector l2_VecDot = new Vector(3);
    
    /**
     * Formula 7.200
     */
    static Vector l2_VecDotDot = new Vector(3);
    
    /**
     * Formula 7.201
     */
    static double g_1[] = new double[3];
    
    /**
     * Formula 7.201
     */
    static double g_2[] = new double[3];
    
    /**
     * Formula 7.202
     */
    static double theta[] = new double[3];
    
    /**
     * Formula 7.202
     */
    static double thetaG0;
    
    /**
     * Formula 7.203
     */
    static Vector observerVec[] = new Vector[3];
    
    /**
     * Formula 7.204
     */
    static Vector observerVec_2_Dot = new Vector(3);
    
    /**
     * Formula 7.204
     */
    static Vector observerVec_2_DotDot = new Vector(3);
    
    /**
     * Formulas 7.207
     */
    static double delta, d_a, d_b, d_c, d_d;
    
    /**
     * Formulas 7.208
     */
    static double a_2_star, b_2_star, c_2_star, d_2_star;
    
    /**
     * Formula 7.209
     */
     static double c_psi;
     
     /**
      * Formulas 7.210
      */
     static double a,b,c;
     
     /**
      * Formula 7.210
      */
     static double observerVec_2_Size;
     
     /**
      * Formula 7.211
      */
     static double r_2[];
     
     /**
      * Formula 7.212
      */
     static double ro_2[];
     
     /**
      * Formula 7.213
      */
     static double ro_2_Dot[];
     
     /**
      * Formula 7.214
      */
     static Vector bodyVec_2 = new Vector(3);
     
     /**
      * Formula 7.215
      */
     static Vector bodyVec_2_Dot = new Vector(3);
     
     /**
      * State vector
      */
     //static StateVector stateVector[];
     
     /**
      * getPositionVectors()
      *  
      * Laplace method to get geocentric position vectors of body from 3 angles (R.A., dec)
      * 
      * INPUT:
      *     double ra[i] - right accesions of body for time[i], i = 0,1,2
      *     double dec[i] - declinations of body for time[i], i = 0,1,2
      *     double time[i] - times of observation of body, i = 0,1,2
      *     Geodetic geodetic[i] - geodetic positions of observer, i = 0,1,2
      *     double dTheta_dTime - Earth rotation
      *     double f - flattening of Earth
      *     double a_e - planet radius (equator)
      *     double mi - 
      *     double k  - 
      * 
      * OUTPUT:
      *     StateVector stateVector[j] - matrix of position vectors of body, 
      *                     j = 0,1,2 ... 8 (depends on no of real polynomial roots, see formula 7.211)
      */
     
     public StateVector[] getPositionVectors(double ra[], double dec[], double time[], 
                Geodetic geodetic[], double dTheta_dTime, double f, double a_e,
                double mi, double k){
         
         /**
          * State vector
          */
         StateVector stateVector[];
         
         /**
          * Array of body positions - matrix of position vectors of body, i = 0,1,2 (1st, 2nd, 3rd observation)
          * j = 0 - 8 (depends on no of real polynomial roots, see formula 7.211)
          */
         Vector bodyPositions[][] = new Vector[8][3];
         
         /**
          * Formula 7.196
          */
         tau_1 = k*(time[0] - time[1]);
         
         /**
          * Formula 7.197
          */
         tau_3 = k*(time[2] - time[1]);
         
         /**
          * Formulas 7.198
          */
         s[0] = -(tau_3/(tau_1*(tau_1 - tau_3)));
         s[1] = -(tau_3 + tau_1)/(tau_1*tau_3);
         s[2] = -(tau_1)/(tau_3*(tau_3 - tau_1));
         s[3] = 2/(tau_1*(tau_1 - tau_3));
         s[4] = 2/(tau_1*tau_3);
         s[5] = 2/(tau_3*(tau_3 - tau_1));
         
         /**
          * Formulas 7.199
          */
         for(int i = 0; i < 3; i++){
             l_Vec[i].v[0] = Math.cos(dec[i])*Math.cos(ra[i]);
             l_Vec[i].v[1] = Math.cos(dec[i])*Math.sin(ra[i]);
             l_Vec[i].v[2] = Math.sin(dec[i]);
             
         }
         
         /**
          * Formulas 7.200
          */
         Vector help_s1_L1 = new Vector(3);
         Vector help_s2_L2 = new Vector(3);
         Vector help_s3_L3 = new Vector(3);
         Vector help_s4_L1 = new Vector(3);
         Vector help_s5_L2 = new Vector(3);
         Vector help_s6_L3 = new Vector(3);
         
         help_s1_L1 = Vector.multiplyVector(l_Vec[0], s[0]);
         help_s2_L2 = Vector.multiplyVector(l_Vec[1], s[1]);
         help_s3_L3 = Vector.multiplyVector(l_Vec[2], s[2]);
         help_s4_L1 = Vector.multiplyVector(l_Vec[0], s[3]);
         help_s5_L2 = Vector.multiplyVector(l_Vec[1], s[4]);
         help_s6_L3 = Vector.multiplyVector(l_Vec[2], s[5]);
         
         l2_VecDot = Vector.addVectors(help_s1_L1, help_s2_L2);
         l2_VecDot = Vector.addVectors(l2_VecDot, help_s3_L3);
         
         l2_VecDotDot = Vector.addVectors(help_s4_L1, help_s5_L2);
         l2_VecDotDot = Vector.addVectors(l2_VecDotDot, help_s6_L3);
         
         /**
          * Formulas 7.201
          */
         for(int i = 0; i < 3; i++){
             g_1[i] = a_e/(Math.sqrt(1 - (2*f - f*f)*Math.sin(geodetic[i].lat)*
                     Math.sin(geodetic[i].lat))) + geodetic[i].altitude;
             g_2[i] = ((1-f)*(1-f)*a_e)/(Math.sqrt(1 - (2*f - f*f)*Math.sin(geodetic[i].lat)*
                     Math.sin(geodetic[i].lat))) + geodetic[i].altitude;
         }
         
         /**
          * Formula 7.202
          */
          //Greenwich siderial time - 1.27 [s]
          thetaG0 = getGreenwichSiderialTime(time[1]);

          for(int i = 0; i < 3; i++){
              theta[i] = thetaG0 + dTheta_dTime * (time[i] - time[1]) + geodetic[i].lon;
          }
          
          /**
           * Formulas 7.203
           */
          for(int i = 0; i < 3; i++){
              observerVec[i].v[0] = -g_1[i]*Math.cos(geodetic[i].lat)*Math.cos(theta[i]);
              observerVec[i].v[1] = -g_1[i]*Math.cos(geodetic[i].lat)*Math.sin(theta[i]);
              observerVec[i].v[2] = -g_2[i]*Math.sin(geodetic[i].lat);
          }
          
          /**
           * Condition, s.268, before formula 7.204
           */
          if((geodetic[0] == geodetic[1])||(geodetic[0]==geodetic[2])){
              /**
               * Formula 7.205
               */
              double dTheta_dT_div_k = dTheta_dTime/k;
              Vector helpVector = new Vector(3);
              helpVector.v[0] = -observerVec[1].v[1];
              helpVector.v[1] =  observerVec[1].v[0];
              helpVector.v[2] =  0;
              observerVec_2_Dot = Vector.multiplyVector(helpVector, dTheta_dT_div_k);
              helpVector.v[0] = -observerVec[1].v[0];
              helpVector.v[1] = -observerVec[1].v[1];
              helpVector.v[2] =  0;
              observerVec_2_DotDot = Vector.multiplyVector(helpVector, dTheta_dT_div_k*dTheta_dT_div_k);
          }
          
          else{
              /**
               * Formula 7.204
               */
             Vector help_s1_R1 = new Vector(3);
             Vector help_s2_R2 = new Vector(3);
             Vector help_s3_R3 = new Vector(3);
             Vector help_s4_R1 = new Vector(3);
             Vector help_s5_R2 = new Vector(3);
             Vector help_s6_R3 = new Vector(3);

             help_s1_R1 = Vector.multiplyVector(observerVec[0], s[0]);
             help_s2_R2 = Vector.multiplyVector(observerVec[1], s[1]);
             help_s3_R3 = Vector.multiplyVector(observerVec[2], s[2]);
             help_s4_R1 = Vector.multiplyVector(observerVec[0], s[3]);
             help_s5_R2 = Vector.multiplyVector(observerVec[1], s[4]);
             help_s6_R3 = Vector.multiplyVector(observerVec[2], s[5]);

             observerVec_2_Dot = Vector.addVectors(help_s1_R1, help_s2_R2);
             observerVec_2_Dot = Vector.addVectors(observerVec_2_Dot, help_s3_R3);

             observerVec_2_DotDot = Vector.addVectors(help_s4_R1, help_s5_R2);
             observerVec_2_DotDot = Vector.addVectors(observerVec_2_DotDot, help_s6_R3);
          }
          
          /**
           * Formulas 7.207 - determinants
           */
          Matrix delta_Matrix = new Matrix(3,3);
          Matrix d_a_Matrix = new Matrix(3,3);
          Matrix d_b_Matrix = new Matrix(3,3);
          Matrix d_c_Matrix = new Matrix(3,3);
          Matrix d_d_Matrix = new Matrix(3,3);
          
          /*delta_Matrix*/
          //1st column
          delta_Matrix.matrix[0][0] = l_Vec[1].v[0];
          delta_Matrix.matrix[0][1] = l_Vec[1].v[1];
          delta_Matrix.matrix[0][2] = l_Vec[1].v[2];
          //2nd column
          delta_Matrix.matrix[1][0] = l2_VecDot.v[0];
          delta_Matrix.matrix[1][1] = l2_VecDot.v[1];
          delta_Matrix.matrix[1][2] = l2_VecDot.v[2];
          //3rd column
          delta_Matrix.matrix[2][0] = l2_VecDotDot.v[0];
          delta_Matrix.matrix[2][1] = l2_VecDotDot.v[1];
          delta_Matrix.matrix[2][2] = l2_VecDotDot.v[2];
          
          /*d_a_Matrix*/
          //1st column
          d_a_Matrix.matrix[0][0] = l_Vec[1].v[0];
          d_a_Matrix.matrix[0][1] = l_Vec[1].v[1];
          d_a_Matrix.matrix[0][2] = l_Vec[1].v[2];
          //2nd column
          d_a_Matrix.matrix[1][0] = l2_VecDot.v[0];
          d_a_Matrix.matrix[1][1] = l2_VecDot.v[1];
          d_a_Matrix.matrix[1][2] = l2_VecDot.v[2];
          //3rd column
          d_a_Matrix.matrix[2][0] = observerVec_2_DotDot.v[0];
          d_a_Matrix.matrix[2][1] = observerVec_2_DotDot.v[1];
          d_a_Matrix.matrix[2][2] = observerVec_2_DotDot.v[2];
          
          /*d_b_Matrix*/
          //1st column
          d_b_Matrix.matrix[0][0] = l_Vec[1].v[0];
          d_b_Matrix.matrix[0][1] = l_Vec[1].v[1];
          d_b_Matrix.matrix[0][2] = l_Vec[1].v[2];
          //2nd column
          d_b_Matrix.matrix[1][0] = l2_VecDot.v[0];
          d_b_Matrix.matrix[1][1] = l2_VecDot.v[1];
          d_b_Matrix.matrix[1][2] = l2_VecDot.v[2];
          //3rd column
          d_b_Matrix.matrix[2][0] = observerVec[1].v[0];
          d_b_Matrix.matrix[2][1] = observerVec[1].v[1];
          d_b_Matrix.matrix[2][2] = observerVec[1].v[2];
          
          /*d_c_Matrix*/
          //1st column
          d_c_Matrix.matrix[0][0] = l_Vec[1].v[0];
          d_c_Matrix.matrix[0][1] = l_Vec[1].v[1];
          d_c_Matrix.matrix[0][2] = l_Vec[1].v[2];
          //2nd column
          d_c_Matrix.matrix[1][0] = observerVec_2_DotDot.v[0];;
          d_c_Matrix.matrix[1][1] = observerVec_2_DotDot.v[1];;
          d_c_Matrix.matrix[1][2] = observerVec_2_DotDot.v[2];;
          //3rd column
          d_c_Matrix.matrix[2][0] = l2_VecDotDot.v[0];
          d_c_Matrix.matrix[2][1] = l2_VecDotDot.v[1];
          d_c_Matrix.matrix[2][2] = l2_VecDotDot.v[2];
          
           /*d_d_Matrix*/
          //1st column
          d_d_Matrix.matrix[0][0] = l_Vec[1].v[0];
          d_d_Matrix.matrix[0][1] = l_Vec[1].v[1];
          d_d_Matrix.matrix[0][2] = l_Vec[1].v[2];
          //2nd column
          d_d_Matrix.matrix[1][0] = observerVec[1].v[0];
          d_d_Matrix.matrix[1][1] = observerVec[1].v[1];
          d_d_Matrix.matrix[1][2] = observerVec[1].v[2];
          //3rd column
          d_d_Matrix.matrix[2][0] = l2_VecDotDot.v[0];
          d_d_Matrix.matrix[2][1] = l2_VecDotDot.v[1];
          d_d_Matrix.matrix[2][2] = l2_VecDotDot.v[2];
          
          /**
           * Get the determinants of matrices, still formulas 7.207
           */
          delta = 2 * delta_Matrix.getDeterminant(delta_Matrix);
          d_a = d_a_Matrix.getDeterminant(d_a_Matrix);
          d_b = d_b_Matrix.getDeterminant(d_b_Matrix);
          d_c = d_c_Matrix.getDeterminant(d_c_Matrix);
          d_d = d_d_Matrix.getDeterminant(d_d_Matrix);
          
          /**
           * Formulas 7.208
           */
          a_2_star = 2*d_a/delta;
          b_2_star = 2*d_b/delta;
          c_2_star = d_c/delta;
          d_2_star = d_d/delta;
          
          /**
           * Formula 7.209
           */
          c_psi = -2 * (l_Vec[1].v[0]*observerVec[1].v[0] + 
                        l_Vec[1].v[1]*observerVec[1].v[1] +
                        l_Vec[1].v[2]*observerVec[1].v[2]);
          
          /**
           * Formulas 7.210
           */
          observerVec_2_Size = Vector.getSize(observerVec[1]);
          a = -(c_psi*a_2_star + a_2_star*a_2_star + observerVec_2_Size*observerVec_2_Size);
          b = -mi*(c_psi*b_2_star + 2*a_2_star*b_2_star);
          c = -mi*mi*b_2_star*b_2_star;
          
          /**
           * Formula 7.211 - solve the polynomial of 8th grade
           */
          //looking for polynomial roots
            double coef[] = new double[9];
            coef[0] = c;
            coef[1] = 0;
            coef[2] = 0;
            coef[3] = b;
            coef[4] = 0;
            coef[5] = 0;
            coef[6] = a;
            coef[7] = 0;
            coef[8] = 1;
            Polynomial poly = new Polynomial(coef);
            Complex complex[] = new Complex[8];
            try{
                complex = poly.zeros();
                double realRoot[] = new double[8]; 
                //looking for real roots
                int j = 0; 
                for(int i = 0; i < 8; i++){
                    if((complex[i].im() == 0)||(Math.abs(complex[i].im()) < 10e-10)){
                        realRoot[j] = complex[i].re();
                        //System.out.println("Real root no." + j + " :" + realRoot[j]);
                        j++;
                    }
                    else{
                        //System.out.println("Root is: " + complex[i]);
                    }
                }

                //looking for all r_2s
                if(j > 0) {
                    //System.out.println("260");
                    r_2 = new double[j];
                    bodyPositions = new Vector[j][3];
                    for(int l = 0; l<j; l++){
                        for(int i = 0; i < 3; i++){
                            bodyPositions[l][i] = new Vector(3);
                        }
                    }
                }
                else {
                    //System.out.println("No real roots of polynomial!!!");
                    return null;
                }
                
                /**
                * Formulas 7.212 - 7.215
                */
                stateVector = new StateVector[j];
                for(int l = 0; l < j; l++){
                    ro_2 = new double[j];
                    
                    /**
                     * Formula 7.212
                     */
                    ro_2[l] = a_2_star*a_2_star + ((mi*b_2_star)/(Math.pow(r_2[l],3)));
                    
                    /**
                     * Formula 7.213
                     */
                    ro_2_Dot[l] = c_2_star*c_2_star + ((mi*d_2_star)/(Math.pow(r_2[l],3)));
                    
                    /**
                     * Formula 7.214
                     */
                    stateVector[l].r = Vector.subtractVectors(Vector.multiplyVector(l_Vec[1], ro_2[l]), observerVec[1]);
                    
                    /**
                     * Formula 7.215
                     */
                    stateVector[l].r = Vector.addVectors(Vector.multiplyVector(l_Vec[1], ro_2_Dot[l]), 
                                                           Vector.multiplyVector(l2_VecDot, ro_2[l]));
                    stateVector[l].r =  Vector.subtractVectors( stateVector[l].r, 
                                            observerVec_2_Dot);
                }
         
            }
            catch(Exception e){
                System.out.println("Warning! " + e);
                return null;
            }
            
               
         return stateVector;
     }
    
    /**
     * getGreenwichSiderialTime()
     *
     * IN:
     *  double timeJd - jullian date [???]
     *
     * OUT:
     *  double gst [rad]
     */

    public double getGreenwichSiderialTime(double timeJd){
        double t_u = (timeJd - 2415020.0)/36525;
        double gst = 99.6909833 + 36000.7689*t_u + 0.00038708*t_u*t_u;
        gst = Math.toRadians(gst);
        return gst;
    }
}
