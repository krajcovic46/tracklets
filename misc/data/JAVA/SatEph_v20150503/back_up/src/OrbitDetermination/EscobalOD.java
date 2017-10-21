/*
 * Class to compute orbital elements from position and velocity vectors of unknown body.
 * Chapters 3.7.2 - 3.8.1
 * All formulas are from book P.Escobal - Methods of orbit
 * determination (1976). All formulas numbers in this code are identical with
 * formulas numbers used in Escobal's book.
 *
 * Author of Java source is Jiri Silha - 2009
 */

package OrbitDetermination;

import Compute.Silha.*;
/**
 *
 * @author Jiri Silha - 02/06/2009
 */
public class EscobalOD {
     
    // main data keeper
     Kepler kepler;
    
    /*
     *incomming data
    */
    
    //position and velocity vectors
     Vector r, r_Dot;
    
    //corresponding universal time, mass
     double time, mi;
    
    /*
     * variables
     */
    //formuals 3.150, 3.151, 3.152, 3.155
    double r_Size, r_ScalarMult_rDot, r_Dot_Size, v_0;
    
    //formulas 3.154, 3.155, 3.159
    public double c_e, s_e, eccAnom_Dot, eccAnom;
        
    //formulas 3.162, 3.167
    public Vector uVec, vVec;
    
    //formula 3.166
    double p;
    
    //formulas 3.195, 3.196, 3.197, 3.198, 3.200, 3.204, 3.205
    double sin_i, cos_i, sin_l, cos_l, sin_u, cos_u, sin_v, cos_v;
    
    //formulas 3.199, 3.206
    double l, u, v;
    
    //formulas 3.201, 3.202
    double c_v, s_v;
    
    /**
     * getElementsFromPosAndVel()
     * 
     * INPUT:
     * Vector r - position vector of body for corresponding universal time
     * Vector r_Dot - velocity vector of body for corresponding universal time
     * double time - corresponding universal time
     * double mi
     * 
     * Output:
     * Kepler kepler - orbital elements of the body orbit (a, e, i omega, Omega, M)
     */
    public Kepler getElementsFromPosAndVel(Vector r, Vector r_Dot, double time,
                                            double mi){
        this.kepler = new Kepler();
        this.r = r;
        this.r_Dot = r_Dot;
        this.time = time;
        this.mi = mi;
        
        //formula 3.150
        r_Size = Vector.getSize(r);
        
        //formula 3.151
        r_ScalarMult_rDot = Vector.getScalarProduct(r, r_Dot);     
        r_Dot_Size = r_ScalarMult_rDot/r_Size;
        
        //formula 3.152
        v_0 = Vector.getSize(r_Dot);
        
        //formula 3.153 - modified
        kepler.a = r_Size*mi/(2*mi - v_0*v_0*r_Size);
        
        //formula 3.154
        c_e = 1 - r_Size/kepler.a;
        
        //formula 3.160
        s_e = r_ScalarMult_rDot/Math.sqrt(mi*kepler.a);
        
        //formula 3.161
        kepler.e = Math.sqrt(s_e*s_e + c_e*c_e);
        
        p = kepler.a*(1.0 - kepler.e*kepler.e);

        //formula 3.162
        uVec = Vector.multiplyVector(r, 1.0/r_Size);
        
        //formula 3.165
        vVec = Vector.multiplyVector(Vector.subtractVectors(Vector.multiplyVector(r_Dot, r_Size),Vector.multiplyVector(r, r_Dot_Size)),
                    1.0/Math.sqrt(mi*p));
        
        //formula 3.156
        eccAnom = new RIterationAnglesOnly().getAngleFromSinAndCos(s_e/kepler.e, c_e/kepler.e);
        kepler.M = eccAnom - s_e;
        
        //formula 3.195, 3.196
        sin_i = Math.sqrt(uVec.v[2]*uVec.v[2] + vVec.v[2]*vVec.v[2]);
        cos_i = Math.sqrt(Math.pow(uVec.v[0] + vVec.v[1],2) +
                          Math.pow(uVec.v[1] - vVec.v[0],2)) - 1;
        
        kepler.incl = new RIterationAnglesOnly().getAngleFromSinAndCos(sin_i, cos_i);
        
        //formulas 3.197, 3.198
        cos_l = (uVec.v[0] + vVec.v[1])/(1 + cos_i);
        sin_l = (uVec.v[1] - vVec.v[0])/(1 + cos_i);
        
        l = new RIterationAnglesOnly().getAngleFromSinAndCos(sin_l, cos_l);
        
        //formulas 3.200
        sin_u = uVec.v[2]/sin_i;
        cos_u = vVec.v[2]/sin_i;
        
        u = new RIterationAnglesOnly().getAngleFromSinAndCos(sin_u, cos_u);
        
        //formula 3.199
        kepler.Omega = l - u;
        while(kepler.Omega < 0) kepler.Omega = kepler.Omega + Math.PI*2;
        
        //formulas 3.201, 3.202
        c_v = p/r_Size - 1;
        s_v = r_Dot_Size*Math.sqrt(p/mi);
        
        //ecc 2
        //System.out.println("e1 " + kepler.e);
        //System.out.println("e2 " + Math.sqrt(c_v*c_v + s_v*s_v));
        
        //formulas 3.204, 3.205
        cos_v = 1/Math.sqrt(c_v*c_v + s_v*s_v) * (p/r_Size - 1);
        sin_v = r_Dot_Size/Math.sqrt(c_v*c_v + s_v*s_v)*Math.sqrt(p/mi);
        
        v = new RIterationAnglesOnly().getAngleFromSinAndCos(sin_v, cos_v);
        
        kepler.omega = u - v;
        while(kepler.omega < 0) kepler.omega = kepler.omega + Math.PI*2;
        
        return kepler;
    }
        
    /**
     * main method - test
     */
    public static void main(String args[]){
        //Reference orbit no.VII (7.7.2)
        //constants
        double mi = 1.0;
        //body position vector
        Vector r = new Vector(3);
        //body velocity vector
        Vector v = new Vector(3);
        //time of observation - JD [min]
        double timeJd = 2438314.8055556*1440;
        
        r.v[0] = -1.22192;
        r.v[1] = 0.0352894;
        r.v[2] = 1.54752;
        
        v.v[0] = -0.468449;
        v.v[1] = -0.513595;
        v.v[2] = -0.175849;
        
        //System.out.println(Vector.getSize(v) + " e.r./min");
        //System.out.println(Vector.getSize(v)*Constants.R_Earth/60/1000 + " km/s");
        //System.out.println(Vector.getSize(v)*25936 + " foot/sec");
        //System.out.println(Vector.getSize(v)*25936*0.3048/1000 + " km/s");
        
        Kepler kepler2 = new EscobalOD().getElementsFromPosAndVel(r,v,timeJd, mi);
        
        //System.out.println("a: "+kepler2.a + " e.r.");
        //System.out.println("e: "+kepler2.e);
        //System.out.println("i: "+Math.toDegrees(kepler2.incl)+" °");
        //System.out.println("Omega: "+Math.toDegrees(kepler2.Omega)+" °");
        //System.out.println("omega: "+Math.toDegrees(kepler2.omega)+" °");
        //System.out.println("Ma: "+Math.toDegrees(kepler2.M)+" °");
        double k2 = 0.07436574;  //[(e.r)^3/2 / min]
        double n = k2*Math.sqrt(mi)*Math.pow(kepler2.a,-1.5);
        //System.out.println("n " + n + "rev/min");
        //System.out.println("n " + n*1440/(2*Math.PI) + " rev/day\n");
        
       
        
        //Reference orbit no. VIII.
        //time of observation - JD [min]
        timeJd = 2438182.1666667*1440;
        
        r.v[0] = -6.82944;
        r.v[1] = 0.388821;
        r.v[2] = 5.10408;
        
        v.v[0] = -0.194685;
        v.v[1] = -0.301251;
        v.v[2] = -0.0753532;
        
        //System.out.println(Vector.getSize(v)*25936 + " foot/sec");
        //System.out.println(Vector.getSize(v)*25936*0.3048/1000 + " km/s");
        
        kepler2 = new EscobalOD().getElementsFromPosAndVel(r,v,timeJd, mi);
        
        //System.out.println("a: "+kepler2.a + " e.r.");
        //System.out.println("e: "+kepler2.e);
        //System.out.println("i: "+Math.toDegrees(kepler2.incl)+" �");
        //System.out.println("Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        //System.out.println("omega: "+Math.toDegrees(kepler2.omega)+" �");
        //System.out.println("Ma: "+Math.toDegrees(kepler2.M)+" �");
        k2 = 0.07436574;  //[(e.r)^3/2 / min]
        n = k2*Math.sqrt(mi)*Math.pow(kepler2.a,-1.5);
        //System.out.println("n " + n + "rev/min");
        //System.out.println("n " + n*1440/(2*Math.PI) + " rev/day\n\n");
        kepler2.a = 1.44;
        n = k2*Math.sqrt(mi)*Math.pow(kepler2.a,-1.5);
        //System.out.println("n " + n + "rev/min");
        //System.out.println("n " + n*1440/(2*Math.PI) + " rev/day\n\n");
        
        
        /**
         * Excerxise 4, chapter 4 
         */
        Vector posVecTest = new Vector(3);        
        Vector velVecTest = new Vector(3);
        posVecTest.v[0] = 1.5;
        posVecTest.v[1] = 2.0;
        posVecTest.v[2] = 0.0;
        velVecTest.v[0] = 0.2;
        velVecTest.v[1] = 0.3;
        velVecTest.v[2] = 0.1;
        
        EscobalOD escobal = new EscobalOD();
        
        kepler2 = escobal.getElementsFromPosAndVel(posVecTest,velVecTest,timeJd, mi);
        //System.out.println("Escobal  a: 1.51515151 e.r.");
        //System.out.println("Computed a: "+kepler2.a + " e.r.");
        //System.out.println("e: "+kepler2.e);
        //System.out.println("i: "+Math.toDegrees(kepler2.incl)+" �");
        //System.out.println("Omega: "+Math.toDegrees(kepler2.Omega)+" �");
        //System.out.println("omega: "+Math.toDegrees(kepler2.omega)+" �");
        //System.out.println("Ma: "+Math.toDegrees(kepler2.M)+" �");
        //System.out.println("Escobal  se: 0.73116346");
        //System.out.println("Computed se: " + escobal.s_e);
        //System.out.println("Escobal  ce: -0.6500000");
        //System.out.println("Computed ce: " + escobal.c_e);
        //System.out.println("Escobal  U_x: 0.6000000");
        //System.out.println("Computed U_x: " + escobal.uVec.v[0]);
        //System.out.println("Escobal  U_y: 0.8000000");
        //System.out.println("Computed U_y: " + escobal.uVec.v[1]);
        //System.out.println("Escobal  U_z: 0.0000000");
        //System.out.println("Computed U_z: " + escobal.uVec.v[2]);
        //System.out.println("Escobal  V_x: -0.15689293");
        //System.out.println("Computed V_x: " + escobal.vVec.v[0]);
        //System.out.println("Escobal  V_y: 0.11766969");
        //System.out.println("Computed V_y: " + escobal.vVec.v[1]);
        //System.out.println("Escobal  V_z: 0.98058081");
        //System.out.println("Computed V_z: " + escobal.vVec.v[2]);
    }

}
