/*
 * Test.java
 *
 * Created on March 21, 2008, 8:59 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package Compute.Silha;

import OrbitDetermination.*;

/**
 *
 * @author Jiri Silha
 */
public class Test {
    
    /** Creates a new instance of Test */
    public Test() {
    }
    
    public static void main (String [] args){
        /*
        //test triedy SunPosition
        Geodetic geodetic = new Geodetic(Math.toRadians(17.2740), Math.toRadians(48.3733), 100);
        double mjd = 54550.166666666664;
        Sun sp = new Sun();
        sp = sp.getSunLocalPositions(mjd, geodetic);
        System.out.println("Az: " + Math.toDegrees(sp.az));
        System.out.println("h: " + Math.toDegrees(sp.h)+"\n");

        System.out.println("ra: " + Math.toDegrees(sp.ra));
        System.out.println("dec: " + Math.toDegrees(sp.dec));
        */

        /*
        * Get the syntetic population of fragments after parent break up
        */
        //break up impulse kg.m/s, BEWARE!!! IMPULSE SHOULD BE 1.445, no 1.445/10
        double impulse = 15;
        //density of material kg/m^3
        double density = 2800;
        //diameter of object (object is sphere) [m], mass of object [kg]
        double diameter = 0.28;
        double mass;
        //correlation factor for mass, function of obcet's shape, object is never perfect filled shpere
        double corrMassFactor = 52;
        //pericentrum
        double qTest, q, Q;
        //TLE string
        String tle_String = "";
        //geoc. position vector of parent body [m], this case is Cosmos 2251
        Vector vec_PosParentBody = new Vector(3);

        //amount of given objects
        double amount = 469;

        /*
        //COSMOS 2251 data
        vec_PosParentBody.v[0] = -1468065.3559454195;
        vec_PosParentBody.v[1] = 1585916.6866229775;
        vec_PosParentBody.v[2] = 6812734.5420169365;
        //geoc. velocity vector of parent body [m], this case is Cosmos 2251
        Vector vec_VelParentBody = new Vector(3);
        vec_VelParentBody.v[0] = -6995.621549549612;
        vec_VelParentBody.v[1] = -2453.961487093762;
        vec_VelParentBody.v[2] = -934.0938610237473;
        */

        
        //Fengyun 1C data
        vec_PosParentBody.v[0] = -5935689.586472171;
        vec_PosParentBody.v[1] = -800896.1797053608;
        vec_PosParentBody.v[2] = 4056988.634341357;
        //geoc. velocity vector of parent body [m], this case is Fengyun 1C
        Vector vec_VelParentBody = new Vector(3);
        vec_VelParentBody.v[0] = -4236.852853661556;
        vec_VelParentBody.v[1] = 787.8079020037669;
        vec_VelParentBody.v[2] = -6042.272208332661;
        
        
        //relative velocity of object in respect to parent body [m/s]
        Vector vec_VelRelToParentBody = new Vector(3);
        //geoc. velocity of fragment [m/s]
        Vector vec_VelFragment = new Vector(3);
        //size of relative velocity of object in respect to parent body [m/s]
        double size_Vec_VelRelToParentBody;

        //geoc. orbital elements of object
        Kepler kepler = new Kepler();

        //spherical coordinates theta <0;2*PI> rad, phi <0;PI> rad
        double theta, phi;

        //amount of bed elements = hyperbol case, decayed orbit
        double wrong_count = 0;

        //coefficient B*Drag
        String bDragString = "37306-3"; //small particles; depends on atm. density
        //String bDragString = "51550-3";
        
        //get TLE of object
        int i = 1;  //23800, 10200, 4200, 1800 fragments
        //int i = 23801;  //23800, 10200, 4200, 1800 fragments
        //int i = 34001;  //23800, 10200, 4200, 1800 fragments
        //int i = 38201;  //23800, 10200, 4200, 1800 fragments
        //int i = 40001;  //23800, 10200, 4200, 1800 fragments
        //while (i <= 1297){
        //while (i <= 23800){
        //while (i <= 34000){
        //while (i <= 38200){
        //while (i <= 41297){
        while (i <= amount){
        //while (i <= 50000){
        //for(int i =0; i<10; i++){
            //get the mass of fragment
            mass = density*(double)4/(double)3*Math.PI*Math.pow(diameter/2,3)/corrMassFactor;
            //get the size of relative velocity
            size_Vec_VelRelToParentBody = impulse/mass;
            //size_Vec_VelRelToParentBody = 0;

            //System.out.println("rel vel " + size_Vec_VelRelToParentBody);

            //get theta, simple model
            //theta = Math.random()*Math.PI*2;
            //get phi, simple model
            //phi = Math.random()*Math.PI;
            //phi = Math.sin(phi)*Math.PI;

            //get relative velocity vector  - spherical coordinates solution, points are not rectangular distribution
            //vec_VelRelToParentBody.v[0] = size_Vec_VelRelToParentBody*Math.cos(theta)*Math.sin(phi);
            //vec_VelRelToParentBody.v[1] = size_Vec_VelRelToParentBody*Math.sin(theta)*Math.sin(phi);
            //vec_VelRelToParentBody.v[2] = size_Vec_VelRelToParentBody*Math.cos(phi);

            //get cartesian elements
            //get relative velocity vector - cartesian solution
            vec_VelRelToParentBody.v[0] = (Math.random()-0.5)*2; //from <-1;1>
            vec_VelRelToParentBody.v[1] = (Math.random()-0.5)*2; //from <-1;1>
            vec_VelRelToParentBody.v[2] = (Math.random()-0.5)*2; //from <-1;1>
            //get the size of obtained vector
            double vectorSize = vec_VelRelToParentBody.getSize(vec_VelRelToParentBody);
            //get unit vector
            vec_VelRelToParentBody = vec_VelRelToParentBody.multiplyVector(vec_VelRelToParentBody, (double)1/vectorSize);
            //System.out.println("Unit rel vel size: " + vec_VelRelToParentBody.getSize(vec_VelRelToParentBody));
            //get vector with size of relativ velocity
            vec_VelRelToParentBody = vec_VelRelToParentBody.multiplyVector(vec_VelRelToParentBody, size_Vec_VelRelToParentBody);
            //System.out.println("Rel vel size: " + vec_VelRelToParentBody.getSize(vec_VelRelToParentBody));

            //get geoc. velocity of fragment
            vec_VelFragment = vec_VelFragment.addVectors(vec_VelParentBody, vec_VelRelToParentBody);
            //System.out.println(vec_VelFragment.getSize(vec_VelFragment) + "\t" +vec_VelFragment.v[0] + "\t" + vec_VelFragment.v[1] + "\t" + vec_VelFragment.v[2]);
            //get the geocentric orbital elements
            kepler = kepler.getElementsFromPosAndVelVec(Constants.GM_Earth,vec_PosParentBody,vec_VelFragment);
            //System.out.println("a " + kepler.a/1000 + " km");
            if(kepler.a > 0){
                String nameFFragment = "Fragment 25730, " + (int)(diameter*100) + " cm, " + (int)(density/1000) + " g/cm^3, " + i;
                tle_String = kepler.getTLE(nameFFragment , i, "00000A",
                //tle_String = kepler.getTLE("Fragment 25730 (type 8) " + i , i, "00000A",
                            //Time.getMjd(2009, 2, 10, 16, 56, 0.0), Math.toDegrees(kepler.incl), //Cosmos time
                            Time.getMjd(2007, 1, 11, 22, 26, 10.0), Math.toDegrees(kepler.incl),  //Fengyun 1C time
                            Math.toDegrees(kepler.Omega), Math.toDegrees(kepler.omega),
                            Math.toDegrees(kepler.M), kepler.a, kepler.e, bDragString);
                //System.out.println(tle_String);

                //test if the TLE is correct!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                int lengthTleString = tle_String.length();
                String tle_line2 = tle_String.substring(lengthTleString - 69);
                //System.out.println("string lenght " + tle_String.length());
                //System.out.println("line2 " + tle_line2);
                Kepler keplerTest = new Kepler();
                keplerTest = keplerTest.fromTleToKepler(tle_line2);
                
                //if perigee is under 180 km altitude, object decayed
                qTest = keplerTest.a * (1 - keplerTest.e);
                q = kepler.a * (1 - kepler.e);
                Q = kepler.a * (1 + kepler.e);
                //System.out.println("q     " + q/1000 + " km" + "\t" + " a     " + kepler.a + "\t" + " e     " + kepler.e);
                //System.out.println("qTest " + qTest/1000 + " km" + "\t" + " aTest " + keplerTest.a + "\t" + " eTest " + keplerTest.e);
                if(qTest < (Constants.R_Earth+180000)) {
                    //System.out.println("q " + q/1000 + " km");
                    wrong_count++;
                }
                else {
                    //System.out.println(tle_String);
                    System.out.println(i + "\t" +kepler.a/1000 + "\t" + kepler.e + "\t" + Math.toDegrees(kepler.incl)
                                 //+ "\t" + sa.array[i].argument  + "\t" + sa.array[i].node + "\t" + sa.array[i].rcs + "\t" + sa.array[i].bDrag);
                                 + "\t" + Math.toDegrees(kepler.omega)  + "\t" + Math.toDegrees(kepler.Omega) + "\t" + (q-Constants.R_Earth)/1000 + "\t" + (Q-Constants.R_Earth)/1000 + "\t" + 2*Math.PI*Math.sqrt(Math.pow(kepler.a,3)/Constants.GM_Earth)/60);
                    i++;
                }
            }
            else wrong_count++;
        }
        System.out.println("wrong_count " + wrong_count);
    }
}
