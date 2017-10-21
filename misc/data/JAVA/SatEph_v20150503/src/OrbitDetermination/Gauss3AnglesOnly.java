/**
 * Class is for computing the geocentric position vectors of body with only
 * 3 topocentric positions on the celestial sphere - Az_i,h_i, or R.A._i, dec_i, for i = 1,2,3.
 * Gauss method is used. All formulas are from book P.Escobal - Methods of orbit
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
 * @author Jiri Silha - 23.03.2009
 * /
 
/** Computing the geocentric position vectors of body*/ //-haluz
/** Computing the geocentric position vectors of body*/
public class Gauss3AnglesOnly {
    
    /**
     * Whether or not has the polynomial of 8th grade real root
     */
    public boolean hasPolyRealRoot = true;
    
    
    /**
     * getPositionVectors()
     *  
     * Method to get geocentric position vectors of body from 3 angles (R.A., dec)
     *
     * IN:
     *
     *  //body positions
     *  double dec[i] - declinations , i = 0,1,2, [rad]
     *  double ra[i]  - right ascension, i = 0,1,2 [rad]
     *  
     *  //observer position
     *  double lon[j] - geodetic longitude, j = 0,1,2 [rad]
     *  double lat[j] - geodetic latitude, j = 0,1,2 [rad]
     *  double alt[j] - geodetic latitude, j = 0,1,2 [m]
     * OR
     *  Geodetic geodetic[j] - geodetic coordinates, j = 0,1,2, see class Geodetic
     *
     *  //observation times
     *  dounle timeMjd[i] - i = 0,1,2, MJD - [sec]
     *
     *  //constants
     *  double dTheta_dTime - [rad/sec]
     *  double flattening - [-]
     *  double planetRadius - [m]
     *  double gm - gravitational coefficient, [m^3/s^2]
     *  double k
     *  double mi
     *  double a_e - radius of Earth [m]
     *
     * OUT:
     *
     *  Vector position[k] - position vectors of body, k = 0,1,2
     */

    public Vector[][] getPositionVectors(double dec[], double ra[], Geodetic geodetic[],
                double time[], double dTheta_dTime, double flattening, double a_e,
                double k, double mi){//, double gm, double a_e){

        //position vectors of body
        Vector position[] = new Vector[3];
        Vector position2[][];// = new Vector[][3];
        for(int l = 0; l < 3; l++){
            position[l] = new Vector(3);
        }

        //VARIABLES
        //Escobal formulas 7.104 - 7.110, - tau is from 1.9, or 1.32, t0 is second time
        double tau_1, tau_3, tau_13, a_1, b_1, a_3, b_3;
        //Escobal formulas 7.111 - 7.113, topocentric unit vector position of body
        Vector lVec[] = new Vector[3];
        for(int l = 0; l < 3; l++){
            lVec[l] = new Vector(3);
        }
        //1.27 - greenwich sidereal time
        double thetaG0;
        //7.114
        double theta[] = new double[3];
        //7.115
        double g_1[] = new double[3];
        //7.116
        double g_2[] = new double[3];
        //observer cartesian position - formula 7.117
        Vector observerPosition[] = new Vector[3];
        for(int i = 0; i < 3; i++){
            observerPosition[i] = new Vector(3);
        }
        //7.118
        double d;
        //7.119
        double a_11, a_12, a_13, a_21, a_22, a_23, a_31, a_32, a_33;
        //7.120 - 7.124
        Vector vecA = new Vector(3);
        Vector vecB = new Vector(3);
        Vector vecX = new Vector(3);
        Vector vecY = new Vector(3);
        Vector vecZ = new Vector(3);
        //7.125 - 7.130
        double a_2star, b_2star, c_psi, r_2square, a, b, c;
        //7.131, array, because r_2 is a root of 8th grade polynomial
        double r_2[];
        //7.132
        double u_2;
        //7.133 - 7.136
        double d_1, d_3, a_1star, b_1star, a_3star, b_3star;
        //7.137 - 7.139 - body distances from observer
        double ro[] = new double[3];
        //observer position - from 7.117
        Vector obsPos_1, obsPos_2, obsPos_3 = new Vector(3);

        //iteration condition [m]
        double epsilon = 10e6;
        
        //constants
        //gravitational constant 1.15
        //double k = Math.sqrt(gm);
        //double mi = gm;
        
        //time JD from MJD [s]
        //double timeJd[] = new double[3];
        //for(int i = 0; i < 3; i++){
         //   timeJd[i] = timeMjd[i] + 2400000.5*86400;
        //}
        
        //System.out.println("Beginning!");
        
        //Escobal's formulas
        
        //7.104 - 7.106
        tau_1 = k * (time[0] - time[1]);
        tau_3 = k * (time[2] - time[1]);
        tau_13 = tau_3 - tau_1;
        //7.107 - 7.110
        a_1 = tau_3/tau_13;
        b_1 = ((tau_13*tau_13 - tau_3*tau_3)*a_1)/6;
        
        //test
        //System.out.println("tau_1: " + tau_1);
        //System.out.println("tau_3: " + tau_3);
        //System.out.println("tau_13: " + tau_13);
        
        a_3 = (-1) * tau_1/tau_13;
        b_3 = ((tau_13*tau_13 - tau_1*tau_1)*a_3)/6;
        //7.111 - 7.113
        for(int i = 0; i < 3; i++){
            lVec[i].v[0] = Math.cos(dec[i]) * Math.cos(ra[i]);
            lVec[i].v[1] = Math.cos(dec[i]) * Math.sin(ra[i]);
            lVec[i].v[2] = Math.sin(dec[i]);
        }
        //Greenwich siderial time - 1.27 [s]
        thetaG0 = getGreenwichSiderialTime(time[1]);
        //7.114
        theta[0] = thetaG0 + dTheta_dTime * (time[0] - time[1]) + geodetic[0].lon;
        theta[1] = thetaG0 + dTheta_dTime * (time[1] - time[1]) + geodetic[1].lon;
        theta[2] = thetaG0 + dTheta_dTime * (time[2] - time[1]) + geodetic[2].lon;
        //7.155
        for(int i = 0; i < 3; i++){
            g_1[i] = a_e / Math.sqrt(1 - (2*flattening - flattening*flattening)
                        *Math.sin(geodetic[i].lat)*Math.sin(geodetic[i].lat))
                        + geodetic[i].altitude;
        }
        //7.116
        for(int i = 0; i < 3; i++){
            g_2[i] = (Math.pow(1 - flattening,2)*a_e) / Math.sqrt(1 - (2*flattening - flattening*flattening)
                        *Math.sin(geodetic[i].lat)*Math.sin(geodetic[i].lat))
                        + geodetic[i].altitude;
        }
        //7.117
        for(int i = 0; i < 3; i++){
            observerPosition[i].v[0] = (-1) * g_1[i] * Math.cos(geodetic[i].lat) * Math.cos(theta[i]);
            observerPosition[i].v[1] = (-1) * g_1[i] * Math.cos(geodetic[i].lat) * Math.sin(theta[i]);
            observerPosition[i].v[2] = (-1) * g_2[i] * Math.sin(geodetic[i].lat);
        }

        //7.118
        d =  lVec[0].v[0] * (lVec[1].v[1] * lVec[2].v[2] - lVec[1].v[2] * lVec[2].v[1])
          - lVec[1].v[0] * (lVec[0].v[1] * lVec[2].v[2] - lVec[0].v[2] * lVec[2].v[1])
          + lVec[2].v[0] * (lVec[0].v[1] * lVec[1].v[2] - lVec[0].v[2] * lVec[1].v[1]);

        //7.119
        a_11 =          (lVec[1].v[1] * lVec[2].v[2] - lVec[2].v[1] * lVec[1].v[2])/d;
        a_12 = (-1) *   (lVec[1].v[0] * lVec[2].v[2] - lVec[2].v[0] * lVec[1].v[2])/d;
        a_13 =          (lVec[1].v[0] * lVec[2].v[1] - lVec[2].v[0] * lVec[1].v[1])/d;
        a_21 = (-1) *   (lVec[0].v[1] * lVec[2].v[2] - lVec[2].v[1] * lVec[0].v[2])/d;
        a_22 =          (lVec[0].v[0] * lVec[2].v[2] - lVec[2].v[0] * lVec[0].v[2])/d;
        a_23 = (-1) *   (lVec[0].v[0] * lVec[2].v[1] - lVec[2].v[0] * lVec[0].v[1])/d;
        a_31 =          (lVec[0].v[1] * lVec[1].v[2] - lVec[1].v[1] * lVec[0].v[2])/d;
        a_32 = (-1) *   (lVec[0].v[0] * lVec[1].v[2] - lVec[1].v[0] * lVec[0].v[2])/d;
        a_33 =          (lVec[0].v[0] * lVec[1].v[1] - lVec[1].v[0] * lVec[0].v[1])/d;

        //7.120
        vecA.v[0] = a_1;
        vecA.v[1] = -1;
        vecA.v[2] = a_3;
        //7.121
        vecB.v[0] = b_1;
        vecB.v[1] = 0;
        vecB.v[2] = b_3;
        //7.122
        vecX.v[0] = observerPosition[0].v[0];
        vecX.v[1] = observerPosition[1].v[0];
        vecX.v[2] = observerPosition[2].v[0];
        //7.123
        vecY.v[0] = observerPosition[0].v[1];
        vecY.v[1] = observerPosition[1].v[1];
        vecY.v[2] = observerPosition[2].v[1];
        //7.124
        vecZ.v[0] = observerPosition[0].v[2];
        vecZ.v[1] = observerPosition[1].v[2];
        vecZ.v[2] = observerPosition[2].v[2];

        //7.125
        a_2star = (-1)*(  vecA.getScalarProduct(vecA.multiplyVector(vecA, a_21), vecX)
                    + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_22), vecY)
                    + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_23), vecZ)
                    );
        b_2star = (-1)*(  vecB.getScalarProduct(vecB.multiplyVector(vecB, a_21), vecX)
                    + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_22), vecY)
                    + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_23), vecZ)
                    );
        //7.126
        c_psi = (-2)*(vecX.v[1]*lVec[1].v[0] + vecY.v[1]*lVec[1].v[1] + vecZ.v[1]*lVec[1].v[2]);
        //7.127
        r_2square = vecX.v[1]*vecX.v[1] + vecY.v[1]*vecY.v[1] + vecZ.v[1]*vecZ.v[1];
        //7.128
        a = (-1)*(c_psi*a_2star + a_2star*a_2star + r_2square);
        //7.129
        b = (-1)*mi*(c_psi*b_2star + 2*a_2star*b_2star);
        //7.130
        c = (-1)*mi*mi*b_2star*b_2star;
        
        //System.out.println("Polynomial roots!");
        
        //7.131 - get r2, solution of polynomial - 8th grade
        //looking for polynomial roots
        double coef[] = new double[9];
        //a = 1;
        //b = 2;
        //c = 3;
        
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
                //if((complex[i].im() == 0)||(Math.abs(complex[i].im()) < 10e-10)){
                if(complex[i].im() == 0.0){
                    realRoot[j] = complex[i].re();
                    //System.out.println("Real root no." + j + " :" + realRoot[j]);
                    j++;
                }
                else{
                    //System.out.println("Root is: " + complex[i]);
                }
            }
            
            //test
            //j = 1;
           
            //looking for all r_2s
            if(j > 0) {
                //System.out.println("260");
                r_2 = new double[j];
                position2 = new Vector[j][3];
                for(int i = 0; i<j; i++){
                    //r_2[i] = realRoot[i]; 
                    r_2[i] = 12000000; //test
                    for(int l = 0; l < 3; l++){
                        //System.out.println("273");
                        position2[i][l] = new Vector(3);
                    }
                }
            }
            else {
                //System.out.println("No real roots of polynomial!!!");
                return null;
            }            
            
            for(int l = 0; l < j; l++){
                //System.out.println("280");
                //System.out.println("r_2 " + r_2[l]/1000 + " km");
                //7.132
                u_2 = mi/Math.pow(r_2[l],3);
                //7.133
                d_1 = a_1 + b_1*u_2;
                //7.134
                d_3 = a_3 + b_3*u_2;
                
                //test
                //System.out.println("u_2: " + u_2);
                //System.out.println("a_1: " + a_1);
                //System.out.println("a_3: " + a_3);
                //System.out.println("b_1: " + b_1);
                //System.out.println("b_3: " + b_3);
                //System.out.println("d_1: " + d_1);
                //System.out.println("d_3: " + d_3);
                                
                //7.135
                a_1star =  (  vecA.getScalarProduct(vecA.multiplyVector(vecA, a_11), vecX)
                            + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_12), vecY)
                            + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_13), vecZ)
                            );
                b_1star =  (  vecB.getScalarProduct(vecB.multiplyVector(vecB, a_11), vecX)
                            + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_12), vecY)
                            + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_13), vecZ)
                            );
                //7.136
                a_3star =  (  vecA.getScalarProduct(vecA.multiplyVector(vecA, a_31), vecX)
                            + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_32), vecY)
                            + vecA.getScalarProduct(vecA.multiplyVector(vecA, a_33), vecZ)
                            );
                b_3star =  (  vecB.getScalarProduct(vecB.multiplyVector(vecB, a_31), vecX)
                            + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_32), vecY)
                            + vecB.getScalarProduct(vecB.multiplyVector(vecB, a_33), vecZ)
                            );

                //7.137
                ro[0] = (a_1star + b_1star*u_2)/d_1;
                //7.138
                ro[1] = (a_2star + b_2star*u_2);
                //7.139
                ro[2] = (a_3star + b_3star*u_2)/d_3;

                for(int i = 0; i < 3; i++){
                    position[i] = Vector.subtractVectors(lVec[i].multiplyVector(lVec[i], ro[i]), observerPosition[i]);
                    //System.out.println("ros in loop " + ro[i]);
                }

                ro = getHerrickGibbsIteration(position, lVec, ro, tau_1, tau_3, tau_13, mi,
                      observerPosition, a_11, a_12, a_13, a_21, a_22, a_23, a_31, a_32, a_33,
                      epsilon);

                for(int i = 0; i < 3; i++){
                    //position[i] = lVec[i].subtructVectors(lVec[i].multiplyVector(lVec[i], ro[i]), observerPosition[i]);
                    position2[l][i] = Vector.subtractVectors(lVec[i].multiplyVector(lVec[i], ro[i]), observerPosition[i]);
                }
            }
        }
        catch(Exception e){
            //System.out.println("Exception in get position vectors method!");
            //System.out.println(e.getMessage());
            return null;
        }
        //System.out.println("Finishing!");
        return position2;
    }

    /**
     * getHerrickGibbsIteration()
     *
     * Method to compute ro distances betwen observer and body, iteration is
     * used
     * 
     * IN:
     *  //body geocentric position - 1st approximation from method getPositionVectors()
     *  Vector position[k] - position vectors of body, k = 0,1,2
     *  Vector lVec[k] - position vectors of body on celestial sphere, k = 0,1,2, 7.111 - 113
     *  double ro[] - ro distances, [m], k = 0,1,2
     *  double tau_1, tau_3, tau_13 - see method getPositionVectors()
     *  double mi - geocentric dimensionless ratios of masses formula 1.4, 1.13
     *  Vector observerPosition[] - geocentric observer position (vernal equinox coor system)
     *  double a_11, a_12, a_13, a_21, a_22, a_23, a_31, a_32, a_33 - formula 7.119
     *  double epsilon - tolerances, 7.158
     *
     * OUT:
     *  double ro_k - ro distances, [m], k = 0,1,2
     */

    public double[] getHerrickGibbsIteration(Vector position[], Vector lVec[], double ro[],
              double tau_1, double tau_3, double tau_13, double mi, Vector observerPosition[],
              double a_11, double a_12, double a_13, double a_21, double a_22,
              double a_23, double a_31, double a_32, double a_33, double epsilon){
              
        //System.out.println("Herrick Gibbs Iteration!");
        
        //help ro distances
        double ro_n[] = new double[3];
        double ro_n1[] = new double[3];

        //Escobal formulas 7.141 - 7.143
        double d_1, d_2, d_3;
        //7.144 - velocity vector
        Vector position_2Dot = new Vector(3);
        //7.145 - 7.146
        double position_2Size, position_2DotSize;
        //7.147
        double v_2;
        //7.148
        double a;
        //position vectors size
        double position_1Size, position_3Size;
        //7.149
        double f_1, f_3, g_1, g_3;
        //7.150
        double d_star;
        //7.151 - 7.153
        double c_1, c_2, c_3;
        //7.154
        Vector gVec = new Vector(3);

        position_1Size = position[0].getSize(position[0]);
        position_2Size = position[1].getSize(position[1]); //7.145
        position_3Size = position[2].getSize(position[2]);
        
        //test
        //System.out.println("position_1Size " + position_1Size/1000 + " km");
        //System.out.println("position_2Size " + position_2Size/1000 + " km");
        //System.out.println("position_3Size " + position_3Size/1000 + " km");
        //new Constants().G_const;//new Constants().GM_Earth;
        //7.141
        d_1 = tau_3*(mi/(12*Math.pow(position_1Size,3)) - (double)1.0/(tau_1*tau_13));
        //7.142
        d_2 = (tau_1 + tau_3)*(mi/(12*Math.pow(position_2Size,3)) - (double)1.0/(tau_1*tau_3));
        //7.143
        d_3 = (-1) * tau_1 * (mi/((12*Math.pow(position_3Size,3))) + (double)1.0/(tau_3*tau_13));
        
        //test
        /*
        d_1 = -(tau_3/(tau_1*(tau_1 - tau_3)));
        d_2 = -(tau_3 + tau_1)/(tau_1*tau_3);
        d_3 = -(tau_1)/(tau_3*(tau_3 - tau_1));
        position_2Dot = position_2Dot.addVectors(position_2Dot.addVectors(position_2Dot.multiplyVector(position[0],d_1),
                           position_2Dot.multiplyVector(position[1],d_2)),
                           position_2Dot.multiplyVector(position[2],d_3));
        */
        
        //7.144
        position_2Dot = position_2Dot.addVectors(position_2Dot.addVectors(position_2Dot.multiplyVector(position[0],(-1)*d_1),
                           position_2Dot.multiplyVector(position[1],d_2)),
                           position_2Dot.multiplyVector(position[2],d_3));
        //7.146
        position_2DotSize = position_2Dot.getScalarProduct(position_2Dot, position_2Dot.multiplyVector(position[1], (double)1.0/position_2Size));
        //7.147
        v_2 = position_2Dot.getSize(position_2Dot);
        //test
        //System.out.println("v_2 " + v_2/1000 + " km/s");
        
        //v_2 = 1000;
        //position_2Size = 11934876.150303153;
        
        //test - 7.148
        //pozor!!! gm ma byt mi, ale to je pred tym rovne mi = 1, what the f...?:-{
        //double semMajAxis = (position_2Size*new Constants().GM_Earth)/
        //                    (2*new Constants().GM_Earth - v_2*v_2*position_2Size);
        double semMajAxis = (position_2Size*mi)/
                            (2*mi - v_2*v_2*position_2Size);
        //System.out.println("r_2: " + position_2Size);
        //System.out.println("sma: " + semMajAxis);

        //7.149
        f_1 = getFSerie(v_2, position_2Size, position_2DotSize, tau_1);
        f_3 = getFSerie(v_2, position_2Size, position_2DotSize, tau_3);
        g_1 = getGSerie(v_2, position_2Size, position_2DotSize, tau_1);
        g_3 = getGSerie(v_2, position_2Size, position_2DotSize, tau_3);

        //7.150
        d_star = f_1*g_3 - f_3*g_1;
        //7.151
        c_1 = g_3/d_star;
        //7.152
        c_2 = -1.0;
        //7.153
        c_3 = -g_1/d_star;
        //7.154
        gVec = gVec.addVectors(gVec.addVectors(gVec.multiplyVector(observerPosition[0], c_1),
                gVec.multiplyVector(observerPosition[1], c_2)),
                gVec.multiplyVector(observerPosition[2], c_3));
        //7.155
        ro_n[0] = ((double)1.0/c_1)*(a_11*gVec.v[0] + a_12*gVec.v[1] + a_13*gVec.v[2]);
        //7.156
        ro_n[1] = (-1)*(a_21*gVec.v[0] + a_22*gVec.v[1] + a_23*gVec.v[2]);
        //1.157
        ro_n[2] = ((double)1.0/c_3)*(a_31*gVec.v[0] + a_32*gVec.v[1] + a_33*gVec.v[2]);
        
        ro_n1[0] = ro_n[0];
        ro_n1[1] = ro_n[1];
        ro_n1[2] = ro_n[2];

        ro_n[0] = ro[0];
        ro_n[1] = ro[1];
        ro_n[2] = ro[2];

        double condition = Math.abs(ro_n1[0] - ro_n[0]);

        //iteration from step 7.141 to 7.157
        //iteration condition
        while(condition > epsilon){
            //new positions vectors of body
            position[0] = position[0].subtractVectors(position[0].multiplyVector(lVec[0], ro_n1[0]),observerPosition[0]);
            position[1] = position[1].subtractVectors(position[1].multiplyVector(lVec[1], ro_n1[1]),observerPosition[1]);
            position[2] = position[2].subtractVectors(position[2].multiplyVector(lVec[2], ro_n1[2]),observerPosition[2]);

            //position vectors sizes
            position_1Size = position[0].getSize(position[0]);
            position_2Size = position[1].getSize(position[1]); //7.145
            position_3Size = position[2].getSize(position[2]);
        
            //7.141
            d_1 = tau_3*(mi/(12*Math.pow(position_1Size,3)) - (double)1.0/(tau_1*tau_13));
            //7.142
            d_2 = (tau_1 + tau_3)*(mi/(12*Math.pow(position_2Size,3)) - (double)1.0/(tau_1*tau_3));
            //7.143
            d_3 = (-1) * tau_1 * (mi/((12*Math.pow(position_3Size,3))) + (double)1.0/(tau_3*tau_13));
            
            //7.144
        position_2Dot = position_2Dot.addVectors(position_2Dot.addVectors(position_2Dot.multiplyVector(position[0],(-1)*d_1),
                           position_2Dot.multiplyVector(position[1],d_2)),
                           position_2Dot.multiplyVector(position[2],d_3));
        //7.146
        position_2DotSize = position_2Dot.getScalarProduct(position_2Dot, position_2Dot.multiplyVector(position[1], (double)1.0/position_2Size));
        //7.147
        v_2 = position_2Dot.getSize(position_2Dot);
        //test
        //System.out.println("v_2 " + v_2/1000 + " km/s");
        
        //v_2 = 1000;
        //position_2Size = 11934876.150303153;
        
        //test - 7.148
        //pozor!!! gm ma byt mi, ale to je pred tym rovne mi = 1, what the f...?:-{
        //double semMajAxis = (position_2Size*new Constants().GM_Earth)/
        //                    (2*new Constants().GM_Earth - v_2*v_2*position_2Size);
        semMajAxis = (position_2Size*mi)/
                            (2*mi - v_2*v_2*position_2Size);
        //System.out.println("r_2: " + position_2Size);
        //System.out.println("sma: " + semMajAxis);

        //7.149
        f_1 = getFSerie(v_2, position_2Size, position_2DotSize, tau_1);
        f_3 = getFSerie(v_2, position_2Size, position_2DotSize, tau_3);
        g_1 = getGSerie(v_2, position_2Size, position_2DotSize, tau_1);
        g_3 = getGSerie(v_2, position_2Size, position_2DotSize, tau_3);

        //7.150
        d_star = f_1*g_3 - f_3*g_1;
        //7.151
        c_1 = g_3/d_star;
        //7.152
        c_2 = -1.0;
        //7.153
        c_3 = -g_1/d_star;
        //7.154
        gVec = gVec.addVectors(gVec.addVectors(gVec.multiplyVector(observerPosition[0], c_1),
                gVec.multiplyVector(observerPosition[1], c_2)),
                gVec.multiplyVector(observerPosition[2], c_3));
            //7.155
            ro[0] = ((double)1.0/c_1)*(a_11*gVec.v[0] + a_12*gVec.v[1] + a_13*gVec.v[2]);
            //7.156
            ro[1] = (-1)*(a_21*gVec.v[0] + a_22*gVec.v[1] + a_23*gVec.v[2]);
            //1.157
            ro[2] = ((double)1.0/c_3)*(a_31*gVec.v[0] + a_32*gVec.v[1] + a_33*gVec.v[2]);          

            ro_n[0] = ro_n1[0];
            ro_n[1] = ro_n1[1];
            ro_n[2] = ro_n1[2];

            ro_n1[0] = ro[0];
            ro_n1[1] = ro[1];
            ro_n1[2] = ro[2];

            condition = Math.abs(ro_n1[0] - ro_n[0]);
            
            for(int i = 0; i < 3; i++){
                    //System.out.println("ros in loop 2 " + ro_n[i]);
                }
            
        }
       
        return ro_n1;
    }

    /**
     * getFSerie()
     *
     * Method to compute f serie - formula 3.228
     *
     * IN:
     *  double v - 7.147
     *  double r_2Size - 7.145
     *  double r_2DotSize - 7.146
     *  double tau - formula 1.32 - t0 is 2nd time 
     *  double mi - geocentric dimensionless ratios of masses formula 1.4, 1.13
     *
     *  OUT:
     *  double fSerie - [-]
     */

    public double getFSerie(double v, double r_2Size, double r_2DotSize,
                double tau){
        
        //System.out.println("F serie!");
                
        //f serie
        double fSerie = 0;
        //Escobal formulas 3.128 - 3.129
        double p,q;
        //3.212
        double mi = 1.0;
        double u = mi/(Math.pow(r_2Size,3));

        //3.218
        p = (r_2Size*r_2DotSize)/(r_2Size*r_2Size);
        q = (v*v - r_2Size*r_2Size*u)/(r_2Size*r_2Size);

        fSerie = 1 - 1/2*u*Math.pow(tau,2) + 1/2*u*p*Math.pow(tau,3) +
                1/24*(3*u*q - 15*u*p*p + u*u)*Math.pow(tau,4) +
                1/8*(7*u*p*p*p - 3*u*p*q - u*u*p)*Math.pow(tau,5) +
                1/720 * (630*u*p*p*q - 24*u*u*q - u*u*u - 45*u*q*q - 945*u*p*p*p*p +
                    210*u*u*p*p)*Math.pow(tau, 6) +
                1/5040*(882*u*u*p*q - 3150*u*u*p*p*p - 9450*u*p*p*p*q +
                    1575*u*p*q*q + 63*u*u*u*p + 10395*u*p*p*p*p*p) * Math.pow(tau,7) +
                1/40320 * (1107*u*u*q*q - 24570*u*u*p*p*q - 2205*u*u*u*p*p +
                    51975*u*u*p*p*p*p - 42535*u*p*p*q*q + 155925*u*p*p*p*p*q+
                    1575*u*q*q*q + 117*u*u*u*q - 135135*u*p*p*p*p*p*p + u*u*u*u) * Math.pow(tau, 8);


        return fSerie;
    }

    /**
     * getGSerie()
     *
     * Method to compute g serie - formula 3.229
     *
     * IN:
     *  double v - 7.147
     *  double r_2Size - 7.145
     *  double r_2DotSize - 7.146
     *  double tau - formula 1.32 - t0 is 2nd time
     *  double mi - geocentric dimensionless ratios of masses formula 1.4, 1.13
     *
     *  OUT:
     *  double fSerie - [-]
     */

    public double getGSerie(double v, double r_2Size, double r_2DotSize,
                double tau){
        
        //System.out.println("G serie!");
        
        //g serie
        double gSerie = 0;
        //Escobal formulas 3.128 - 3.129
        double p,q;
        //3.212
        double mi = 1.0;
        double u = mi/(Math.pow(r_2Size,3));

        //3.218
        p = (r_2Size*r_2DotSize)/(r_2Size*r_2Size);
        q = (v*v - r_2Size*r_2Size*u)/(r_2Size*r_2Size);

        gSerie = tau - 1/6*u*Math.pow(tau,3) + 1/4*u*p*Math.pow(tau,4) +
                1/120*(9*u*q - 45*u*p*p + u*u)*Math.pow(tau,5) +
                1/360*(210*u*p*p*p - 90*u*p*q - 15*u*u*p)*Math.pow(tau,6) +
                1/5040*(3150*u*p*p*q - 54*u*u*q - 225*u*q*q
                    - 4725*u*p*p*p*p + 630*u*u*p*p - u*u*u) * Math.pow(tau,7) +
                1/40320*(3024*u*u*p*q - 12600*u*u*p*p*p - 56700*u*p*p*p*q
                    + 9450*u*p*q*q + 62370*u*p*p*p*p*p + 126*u*u*u*p)* Math.pow(tau,8);
        
        return gSerie;
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
     * main test method
     */
    public static void main(String args[]){
        /*    
        //constants
            Constants c=new Constants();
		
            //observatory position - AGO Modra
            Geodetic geodObs=new Geodetic(Math.toRadians(17.2740),Math.toRadians(48.3733),531.1);
            Geodetic geodeticArray[] = new Geodetic[3];
            for(int i = 0; i < 3; i++){
                geodeticArray[i] = geodObs;
            }
            //Observation data - TIME, AZIMUTH, ELEVATION
            //horizontalne suradnice
            int amount = 4;
            double azi[]=new double[amount];
            double elev[]=new double[amount];
            Time time[]=new Time[amount];
            double timeMjd[] = new double[amount];
            
            double meanAnom[] = new double[amount];
            
            Observation obser[]=new Observation[amount];
            
            Kepler kepler = new Kepler();
            
            //combinations of obsarvations
            int a1,a2,a3;
            
            //Escobal 1.29
            double dTheta_dTime = 4.375269510e-3; //[rad/min]
            
            //Montenbruck flattening
            double flatt = new Constants().f_Earth;
            double e_Radius = new Constants().R_Earth;
            //geocentric dimensionless ratios of masses
            double mi = 1;
            double k = 0.07436574/60;   //e.r.^3/2/sec
            double gravCoef = new Constants().GM_Earth;
            
            //POSITION VECTORS FROM Gauss 3 angles method
            Vector bodyPosition[][];
             
            //  #30582 - SDP4
            //  1 30582U 07004C   08316.83333333 +.00000074 +00000-0 +00000-0 0 03227
            //  2 30582 011.5909 346.5875 8597791 341.0417 356.1707 00.50389163003859
            //
            
            time[0] = new Time(2008, 11, 19, 19, 26, 3.6);
            //azimuth/elevation
            //azi[0] = Math.toRadians(159.858);
            //elev[0] = Math.toRadians(15.929);
            //R.A./dec
            azi[0] = Math.toRadians(29.10276702);
            elev[0] = Math.toRadians(-23.2386643);
            meanAnom[0] = getAngle(1802.819761);

            time[1] = new Time(2008, 11, 19, 19, 41, 3.6);
            //azimuth/elevation
            //azi[1] = Math.toRadians(131.7429999);
            //elev[1] = Math.toRadians(18.94499999);
            //R.A./dec
            azi[1] = Math.toRadians(57.53671946);
            elev[1] = Math.toRadians(-10.1155372);
            meanAnom[1] = getAngle(1804.708985);

            time[2] = new Time(2008, 11, 19, 21, 11, 3.6);
            //azimuth/elevation
            //azi[2] = Math.toRadians(106.195);
            //elev[2] = Math.toRadians(18.549);
            //R.A./dec
            azi[2] = Math.toRadians(100.1123427);
            elev[2] = Math.toRadians(3.563228764);
            meanAnom[2] = getAngle(1816.044329);
            
            time[3] = new Time(2008, 11, 19, 21, 26, 3.6);
            //azimuth/elevation
            //azi[3] = Math.toRadians(106.778);
            //elev[3] = Math.toRadians(19.512);
            //R.A./dec
            azi[3] = Math.toRadians(102.8315162);
            elev[3] = Math.toRadians(3.951628306);
            meanAnom[3] = getAngle(1817.933553);
            
            
            for(int i = 0; i < amount; i++){
                obser[i] = new Observation(azi[i], elev[i], time[i]);
                timeMjd[i] = time[0].getMjd(time[i])*86400;   //from days to sec
            }

            a1 = 1;
            a2 = 2;
            a3 = 3;

            //kepler = (new OrbitDetermination()).getElementsMM(geodObs, obser[a1], obser[a2], obser[a3]);
            
            bodyPosition = (new Gauss3AnglesOnly()).getPositionVectors(elev, azi, geodeticArray, timeMjd, dTheta_dTime/60, 
                        flatt, e_Radius, gravCoef, k, mi, e_Radius);
            
            if(bodyPosition != null){
                System.out.println("length " + bodyPosition.length);
            
                System.out.println("Test 1");
            
                    for(int i = 0; i < bodyPosition.length; i++){
                    System.out.println("In loop");
                    kepler = kepler.getElements(gravCoef,timeMjd[0],timeMjd[1],bodyPosition[i][0],bodyPosition[i][2]);


                    System.out.print("Obs: " + a1 + " - " + a2 + " - " + a3 + "\n\n");

                    System.out.println("a: "+kepler.a/1000+" km");
                    System.out.println("e: "+kepler.e);
                    System.out.println("i: "+Math.toDegrees(kepler.incl)+" °");
                    System.out.println("Omega: "+Math.toDegrees(kepler.Omega)+" °");
                    System.out.println("omega: "+Math.toDegrees(kepler.omega)+" °");
                    System.out.println("Ma: "+Math.toDegrees(kepler.M)+" °\n");

                    System.out.println("SDP Ma: " + meanAnom[a1]);
                    System.out.println("SDP Ma: " + meanAnom[a2]);
                    System.out.println("SDP Ma: " + meanAnom[a3]+"\n");

                    System.out.println("delta Ma2-1: " + Math.abs(meanAnom[a2] - meanAnom[a1]));
                    System.out.println("delta Ma3-2: " + Math.abs(meanAnom[a3] - meanAnom[a2]));
                }
            }
         */
        
        //constants
        double k2 = 0.07436574;  //[(e.r)^3/2 / min]
        double mi2 = 1.0;        //[e.m.] - Earth mass
        double a_e2 = 1.0;       //[e.r.] - Earth radius
        
        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime2 = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening2 = Constants.f_Earth;
        
        /** Orbit no. VII. - 7.6.6 **/
        
        double time2[] = new double[3];  //[day]
        double ra2[] = new double[3];    //[deg]
        double dec2[] = new double[3];   //[deg]
        Geodetic geodetic2[] = new Geodetic[3];
        
        /*
        //initialization of variables
        time2[0] = (Time.getMjd(new Time(1959,3,2,18,4,57.67)) + 2400000.5)*1440;    //[min]
        time2[1] = (Time.getMjd(new Time(1959,3,2,19,19,13.52)) + 2400000.5)*1440;    //[min]
        System.out.println("time2[1] " + (Time.getMjd(new Time(1959,3,2,19,19,13.52)) + 2400000.5));        
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
                        flattening2, a_e2, mi2, k2, false, a_e2*1.1, a_e2*1.11);
        
        //results
        //System.out.println("x =         " + sv.r.v[0]);
        //System.out.println("y =         " + sv.r.v[1]);
        //System.out.println("z =         " + sv.r.v[2]);
        
        //System.out.println("vx =         " + sv.v.v[0]);
        //System.out.println("vy =         " + sv.v.v[1]);
        //System.out.println("vz =         " + sv.v.v[2]);
        System.out.println("Escobal  a = " + 1.3023); 
        System.out.println("Escobal  e = " + 0.16419); 
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
                        flattening2, a_e2, mi2, k2, false, a_e2*2.0, a_e2*2.1);
        
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
        
        /** Reference orbit no. VIII. - 7.7.3 **/
        
        //double time2[] = new double[3];  //[day]
        //double ra2[] = new double[3];    //[deg]
        //double dec2[] = new double[3];   //[deg]
        //Geodetic geodetic2[] = new Geodetic[3];
        
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
        
        
        //StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        //sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
                        //flattening2, a_e2, mi2, k2, false, a_e2*10, a_e2*10.1);
        
        //StateVector sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        //sv2 = new RIterationAnglesOnly().getStateVector(ra2, dec2, time2, geodetic2, dTheta_dTime2,
        //                flattening2, a_e2, mi2, k2, false, a_e2*10, a_e2*10.1);
        Vector sv2[][] = new Gauss3AnglesOnly().getPositionVectors(dec2, ra2, geodetic2, time2, dTheta_dTime2, flattening2, a_e2, k2, mi2);
        
        //results
        //System.out.println("        x = " + sv2[0][1].v[0]);
        //System.out.println("Escobal_x = -6.82944");
        //System.out.println("        y = " + sv2[0][1].v[1]);
        //System.out.println("Escobal_y = 0.388821");
        //System.out.println("        z = " + sv2[0][1].v[2]);
        //System.out.println("Escobal_z = 5.10408");
        
        //System.out.println("        vx = " + sv2[0].v.v[0]);
        //System.out.println("Escobal_vx = -0.194685");
        //System.out.println("        vy = " + sv2[0].v.v[1]);
        //System.out.println("Escobal_vy = -0.301251");
        //System.out.println("        vz = " + sv2[0].v.v[2]);
        //System.out.println("Escobal_vz = -0.0753532\n\n");
        /*
        //test polynomial
        double coef2[] = new double[9];
        coef2[0] = -0.320066;
        coef2[1] = 0;
        coef2[2] = 0;
        coef2[3] = -1.077159;
        coef2[4] = 0;
        coef2[5] = 0;
        coef2[6] = -1.000654;
        coef2[7] = 0;
        coef2[8] = 1;
        Polynomial poly2 = new Polynomial(coef2);
        Complex complex2[] = new Complex[8];
        
        try{
            complex2 = poly2.zeros();
            for(int i = 0; i < 8; i++){
                if(complex2[i].im() == 0) {
                    System.out.println("Root is: " + complex2[i]);
                }
            }
        }
        catch(Exception e){
            System.out.println(e);
        }
        System.out.println("Escobal: 1.269");
        */
        
    }
    
        /*
         *Method to get the angle 
         *IN: [deg]
         *
         *OUT: [deg]      
         */
	public static double getAngle(double angle){
            while(angle > 360) angle = angle - 360;
            
            return angle;
        }

}