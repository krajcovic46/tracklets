/*
 * Class to compute orbits from object positions on celestial sphere (geocentric, or heliocentric).
 */

package OrbitDetermination;

import Compute.Silha.*;

/**
 *
 * @author Jiri Silha - 22.10.2009
 */
public class OrbitDetermination {
    static int noGoodResults = 0;

    /**
     * getTheOdResults()
     *
     * Method to fill given TXT file with results of orbit determination process.
     * Method functions:
     *  1. Method will make every possible combination of 3 observations from given
     *  array (more than 2) of object's observations (see constructor Observation()).
     *  2. Method will compute 3 positions vectors and 2nd (for 2nd observation time)
     *  velocity vector for given observations with method that will be choosed
     * (Bucerius - 0, Escobal - 1).
     *  3. The vecors will be used to compute orbital elements of observed object.
     *  There are several methods to use (0 - Gauss (2 position vector),
     *  1 - Escobal (position and velocity vector for given time)).
     *  4. The results will be writen to given odresults/dirName/RESULTS.TXT file.
     *
     * INPUT:
     *  Observation observation[i] - i > 2, observations of object (see constructor Observation)
     *  int method_1 - which method should be use to compute vectors (Bucerius - 0, Escobal - 1)
     *  int method_2 - which method should be use to compute orbital elements (Gauss - 0, Escobal - 1)
     *  Kepler expKepler - expected orbital elements
     *  String dirName - name of directory, where should the results.txt file saved.
     *
     * OUTPUT:
     *  Data writen to odresults/dir_name/RESULTS.TXT file!!!
     */
    public static void getTheOdResults(Observation observation[], int method_1, int method_2,
                Kepler expKepler, String dirName){
        /*
        * 1. Making combinations
        */

        //make the arrray of integers for method getCombinations()
        int intArray[] = new int[observation.length];
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
        if(noGoodResults > 39999) noGoodResults = 0;
        for(int i = 0; i < vectorArray.length; i++){
            /*
            * 2. Computing vectors (3 positions, 1 velocity)
            */
            //selected observation, that will be combinated
            observation2[0] = observation[(int)vectorArray[i].v[0]];
            observation2[1] = observation[(int)vectorArray[i].v[1]];
            observation2[2] = observation[(int)vectorArray[i].v[2]];
            //3 postions and 1 velocity vector
            vectorsForOd = getVectors(observation2,method_1);

            /*
            System.out.println("method " + method_1);
            System.out.println("x1 " + vectorsForOd[0].v[0]);
            System.out.println("y1 " + vectorsForOd[0].v[1]);
            System.out.println("z1 " + vectorsForOd[0].v[2]);
            System.out.println("x2 " + vectorsForOd[1].v[0]);
            System.out.println("y2 " + vectorsForOd[1].v[1]);
            System.out.println("z2 " + vectorsForOd[1].v[2]);
            System.out.println("x3 " + vectorsForOd[2].v[0]);
            System.out.println("y3 " + vectorsForOd[2].v[1]);
            System.out.println("z3 " + vectorsForOd[2].v[2]);
            System.out.println("vx " + vectorsForOd[3].v[0]);
            System.out.println("vy " + vectorsForOd[3].v[1]);
            System.out.println("vz " + vectorsForOd[3].v[2]);

            System.out.println("Vel [x] " + vectorsForOd[3].v[0]);
            System.out.println("Vel [y] " + vectorsForOd[3].v[1]);
            System.out.println("Vel [z] " + vectorsForOd[3].v[2]);

            System.out.println("velocity " + method_1 + " " + Vector.getSize(vectorsForOd[3]));
            */
            
            /*
            * 3. Computing orbital elements
            */
            //0. case - Gauss method (2 position vectors)
            if(method_2 == 0){
                getKepler = getKepler.getElements(Constants.GM_Earth, observation2[0].timeMjd,
                        observation2[1].timeMjd, vectorsForOd[0], vectorsForOd[1]);
                getKepler.epochMJD = observation2[0].timeMjd;
            }
            //1. case - Escobal position and velocity vectors
            else if(method_2 == 1){
                //from SI units to Escobal units
                //position vectors from [m] to Earth radius
                for(int k=0; k < 3; k++){
                    for(int j=0; j < 3; j++){
                        vectorsForOd[k].v[j] = vectorsForOd[k].v[j]/Constants.R_Earth;
                    }
                }
                //velocity vector
                for(int k=0; k < 3; k++){
                    vectorsForOd[3].v[k] = vectorsForOd[3].v[k]/(25936*0.3048);
                }
                
                getKepler = new EscobalOD().getElementsFromPosAndVel(vectorsForOd[1], vectorsForOd[3],
                        Time.getJdFromMjd(observation2[1].timeMjd)*1440, 1.0);
                //from Earth radius to [m]
                getKepler.a = getKepler.a*Constants.R_Earth;
                getKepler.epochMJD = observation2[1].timeMjd;
            }

            /*
            * 4. Writing the OD results
            */
            //class to work with IO file commands
            WritingODData wod = new WritingODData();
            //create subdirectory and file whre should the results be writen
            wod.createFile(dirName);
            //checking wether or not is results file filled with expKepler - orbital elements that we expected.
            if(wod.isFileEmpty(dirName)) wod.fillResultsFileBasic(dirName, expKepler);
            wod.fillResultsFile(dirName, getKepler, vectorArray[i], method_1, method_2);
            //System.out.println("Data written in to ODResults/" + dirName + "/results.txt!");
            wod.fillResultsFileGood(dirName, getKepler, vectorArray[i], method_1, method_2);
            wod.fillResultsFileTLE(dirName, getKepler, vectorArray[i], method_1, method_2, noGoodResults);
            wod.fillResultsFileVectors(dirName, vectorArray[i], vectorsForOd, observation2);
            //set NORAD number for resutsTLE.txt file
            if(((getKepler.a/1000*(1 - getKepler.e))>3000)&&(getKepler.e<1)){
                noGoodResults++;
            }
            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGood.txt!");
            //get the heliocentric datas
            //get from Escobal units to SI
            if(method_2 == 1){
                for(int k=0; k < 3; k++){
                    for(int j=0; j < 3; j++){
                        vectorsForOd[k].v[j] = vectorsForOd[k].v[j]*Constants.R_Earth;
                    }
                }
            }
            //APRIL 2012
            //wod.fillResultsFileGoodHelio(dirName, vectorsForOd[0], vectorsForOd[2],
            //                         observation2[0].timeMjd, observation2[2].timeMjd, vectorArray[i], method_1, method_2);
            //wod.fillResultsFileGoodGeo(dirName, vectorsForOd[0], vectorsForOd[2],
            //                         observation2[0].timeMjd, observation2[2].timeMjd, vectorArray[i], method_1, method_2);
            //double elements[] = Transformation.getEclipticInclAndNodeFromGeocPosVectors(vectorsForOd[0],vectorsForOd[1],observation2[0].timeMjd,observation2[1].timeMjd);
            //System.out.println("TEST Ecl i    : " + Math.toDegrees(elements[0]));
            //System.out.println("TEST Ecl O    : " + Math.toDegrees(elements[1]));
        }
        System.out.println("Data written in to ODResults/" + dirName + "/results.txt!");
        System.out.println("Data written in to ODResults/" + dirName + "/resultsGood.txt!");
        System.out.println("Data written in to ODResults/" + dirName + "/resultsTLE.txt!");
        System.out.println("Data written in to ODResults/" + dirName + "/resultsVec.txt!");
    }

    /**
     * getVectors()
     *
     * Method to compute position and velocity vectors from more than 3 object's
     * observations (more than 1 different times of observations must be).
     * To compute the vectors the method can use 2 different methods (Bucerius - 0, Escobal - 1).
     *
     * Input:
     *  Observation observation[3], 3 observations of given object (see constructor Observation).
     *  int method - which method should be used to compute the vector.
     *                  (Bucerius - 0, Escobal - 1)
     *
     * Output:
     *  Vector vector[4] - 3 position vectors for 3 times, 1 velocity for 2nd time.
     */
    public static Vector[] getVectors(Observation observation[], int method){
        Vector vectors[] = new Vector[4];

        //0. case Bucerius
        if(method == 0){
            vectors = new MontenbruckMethod().getVectors(observation);
        }
        //1. case Escobal
        else if(method == 1){
            vectors = new RIterationAnglesOnly().getVectors(observation);
        }

        return vectors;
    }

    //testing methods
    public static void main (String args[]){
        //test of hole class
        
        //Observation observation[] = new Observation[6];
        
        /*
         * Molniya 3-47
         * 1 23642U 95042A   09295.27428251 -.00000258  00000-0  10000-3 0  4538
         * 2 23642  63.4104  21.2773 7386232 253.6609  20.2036  2.00605003104094
         */
        /*
        Observation observation[] = new Observation[6];
        //AGO Modra
        for(int i = 0; i<observation.length; i++){
            observation[i] = new Observation();
            observation[i].lon = Math.toRadians(17.2740);
            observation[i].lat = Math.toRadians(48.3733);
            observation[i].alt = 531.1;
        }

        //1. position
        observation[0].timeMjd = Time.getMjd(2009, 10, 23, 20, 32, 50);
        observation[0].ra = Math.toRadians(54.96743);
        observation[0].dec = Math.toRadians(38.24611);

        //2. postion
        observation[1].timeMjd = Time.getMjd(2009, 10, 23, 20, 35, 56);
        observation[1].ra = Math.toRadians(55.36907);
        observation[1].dec = Math.toRadians(38.82363);

        //3. postion
        observation[2].timeMjd = Time.getMjd(2009, 10, 23, 20, 37, 25);
        observation[2].ra = Math.toRadians(55.5611);
        observation[2].dec = Math.toRadians(39.09502);

        //4. postion
        observation[3].timeMjd = Time.getMjd(2009, 10, 23, 20, 38, 30);
        observation[3].ra = Math.toRadians(55.70131);
        observation[3].dec = Math.toRadians(39.29126);

        //5. postion
        observation[4].timeMjd = Time.getMjd(2009, 10, 23, 20, 42, 52);
        observation[4].ra = Math.toRadians(56.26615);
        observation[4].dec = Math.toRadians(40.06601);

        //6. postion
        observation[5].timeMjd = Time.getMjd(2009, 10, 23, 20, 45, 23);
        observation[5].ra = Math.toRadians(56.59158);
        observation[5].dec = Math.toRadians(40.50113);

        Kepler expKepler = new Kepler(Math.toRadians(253.6609), Math.toRadians(21.2773),
                            Math.toRadians(63.4104), 0.7386232,  26556.692, Math.toRadians(20.2));
        //getTheOdResults(observation, 0, 0, expKepler, "test00");
        //getTheOdResults(observation, 0, 1, expKepler, "test01");
        //getTheOdResults(observation, 1, 0, expKepler, "test10");
        //getTheOdResults(observation, 1, 1, expKepler, "test11");
        */

        /*
         *  SL-14 R/B
         *  1 24731U 97006G   09297.78808913  .00000026  00000-0  10000-3 0  4526
         *  2 24731 082.5996 197.2930 0041129 308.1645 051.5736 12.55221830581553
         */
        /*
        Observation observation[] = new Observation[6];
        //AGO Modra
        for(int i = 0; i<observation.length; i++){
            observation[i] = new Observation();
            observation[i].lon = Math.toRadians(17.2740);
            observation[i].lat = Math.toRadians(48.3733);
            observation[i].alt = 531.1;
        }

        //1. position
        observation[0].timeMjd = Time.getMjd(2009, 10, 26, 21, 18, 41);
        observation[0].ra = Math.toRadians(301.5174);
        observation[0].dec = Math.toRadians(70.89869);

        //2. postion
        observation[1].timeMjd = Time.getMjd(2009, 10, 26, 21, 18, 55);
        observation[1].ra = Math.toRadians(308.96548);
        observation[1].dec = Math.toRadians(69.54932);

        //3. postion
        observation[2].timeMjd = Time.getMjd(2009, 10, 26, 21, 19, 12);
        observation[2].ra = Math.toRadians(317.02926);
        observation[2].dec = Math.toRadians(67.42831);

        //4. postion
        observation[3].timeMjd = Time.getMjd(2009, 10, 26, 21, 19, 40);
        observation[3].ra = Math.toRadians(327.84956);
        observation[3].dec = Math.toRadians(62.93759);

        //5. postion
        observation[4].timeMjd = Time.getMjd(2009, 10, 26, 21, 20, 1);
        observation[4].ra = Math.toRadians(334.17972);
        observation[4].dec = Math.toRadians(58.90901);

        //6. postion
        observation[5].timeMjd = Time.getMjd(2009, 10, 26, 21, 20, 20);
        observation[5].ra = Math.toRadians(338.84924);
        observation[5].dec = Math.toRadians(54.88889);

        Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
                            Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation, 0, 0, expKepler, "test00");
        getTheOdResults(observation, 0, 1, expKepler, "test01");
        getTheOdResults(observation, 1, 0, expKepler, "test10");
        getTheOdResults(observation, 1, 1, expKepler, "test11");
        */

        /**
         * Excercise 6 (s. 291)
         */
        /*
        Observation observation[] = new Observation[3];
        //San Fernando, Spain
        for(int i = 0; i<observation.length; i++){
            observation[i] = new Observation();
            observation[i].lon = Math.toRadians(353.79486111);
            observation[i].lat = Math.toRadians(36.4638333);
            observation[i].alt = 3.76e-6*Constants.R_Earth;
        }

        //1. position
        observation[0].timeMjd = Time.getMjd(1959,9,26,21,25,37.403);
        observation[0].ra = Math.toRadians(254.673208333);
        observation[0].dec = Math.toRadians(13.0927222);

        //2. postion
        observation[1].timeMjd = Time.getMjd(1959,9,26,21,26,37.862);
        observation[1].ra = Math.toRadians(260.941125);
        observation[1].dec = Math.toRadians(13.28569444);

        //3. postion
        observation[2].timeMjd = Time.getMjd(1959,9,26,21,27,45.919);
        observation[2].ra = Math.toRadians(269.85375);
        observation[2].dec = Math.toRadians(13.00425);

        Kepler expKepler = new Kepler(Math.toRadians(161.790), Math.toRadians(205.11),
                            Math.toRadians(33.281), 0.23167,  8682.8, Math.toRadians(0));
        getTheOdResults(observation, 0, 0, expKepler, "test00");
        getTheOdResults(observation, 0, 1, expKepler, "test01");
        getTheOdResults(observation, 1, 0, expKepler, "test10");
        getTheOdResults(observation, 1, 1, expKepler, "test11");
        */

        /*
         * VFMO VFMO120091220
         */
        /*
        Observation observation2[] = new Observation[10];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }

        //1. position
        observation2[0].timeMjd = Time.getMjd(2009, 12, 20, 21, 36, 50.11);
        observation2[0].ra = Math.toRadians(63.11892);
        observation2[0].dec = Math.toRadians(10.61556);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2009, 12, 20, 21, 38, 20.11);
        observation2[1].ra = Math.toRadians(63.17183);
        observation2[1].dec = Math.toRadians(10.59528);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2009, 12, 20, 21, 39, 15.264);
        observation2[2].ra = Math.toRadians(63.20746);
        observation2[2].dec = Math.toRadians(10.58311);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2009, 12, 20, 21, 40, 45.12);
        observation2[3].ra = Math.toRadians(63.26075);
        observation2[3].dec = Math.toRadians(10.56256);

        //5. postion
        observation2[4].timeMjd = Time.getMjd(2009, 12, 20, 21, 41, 23.136);
        observation2[4].ra = Math.toRadians(63.28308);
        observation2[4].dec = Math.toRadians(10.55494);

        //6. postion
        observation2[5].timeMjd = Time.getMjd(2009, 12, 20, 21, 42, 52.992);
        observation2[5].ra = Math.toRadians(63.33737);
        observation2[5].dec = Math.toRadians(10.53467);

        //7. postion
        observation2[6].timeMjd = Time.getMjd(2009, 12, 20, 21, 43, 30.144);
        observation2[6].ra = Math.toRadians(63.36025);
        observation2[6].dec = Math.toRadians(10.52481);

        //8. postion
        observation2[7].timeMjd = Time.getMjd(2009, 12, 20, 21, 45, 0);
        observation2[7].ra = Math.toRadians(63.41383);
        observation2[7].dec = Math.toRadians(10.50608);

        //9. postion
        observation2[8].timeMjd = Time.getMjd(2009, 12, 20, 21, 45, 38.016);
        observation2[8].ra = Math.toRadians(63.4375);
        observation2[8].dec = Math.toRadians(10.49794);

        //10. postion
        observation2[9].timeMjd = Time.getMjd(2009, 12, 20, 21, 47, 7.872);
        observation2[9].ra = Math.toRadians(63.49088);
        observation2[9].dec = Math.toRadians(10.47661);

        Kepler expKepler2 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation2, 0, 0, expKepler2, "VFMO120091220_1");
        getTheOdResults(observation2, 0, 1, expKepler2, "VFMO120091220_2");
        getTheOdResults(observation2, 1, 0, expKepler2, "VFMO120091220_3");
        getTheOdResults(observation2, 1, 1, expKepler2, "VFMO120091220_4");
        */

        /*
         * VFMO VFMO090316
         */
         /*
        Observation observation2[] = new Observation[10];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }

        //1. position
        observation2[0].timeMjd = Time.getMjd(2009, 3, 16, 22, 10, 55.73);
        observation2[0].ra = Math.toRadians(184.56512);
        observation2[0].dec = Math.toRadians(17.12072);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2009, 3, 16, 22, 11, 10.73);
        observation2[1].ra = Math.toRadians(184.58958);
        observation2[1].dec = Math.toRadians(17.12783);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2009, 3, 16, 22, 11, 47.58);
        observation2[2].ra = Math.toRadians(184.65812);
        observation2[2].dec = Math.toRadians(17.14622);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2009, 3, 16, 22, 12, 2.58);
        observation2[3].ra = Math.toRadians(184.68076);
        observation2[3].dec = Math.toRadians(17.15228);

        //5. postion
        observation2[4].timeMjd = Time.getMjd(2009, 3, 16, 22, 12, 39.48);
        observation2[4].ra = Math.toRadians(184.74912);
        observation2[4].dec = Math.toRadians(17.17);

        //6. postion
        observation2[5].timeMjd = Time.getMjd(2009, 3, 16, 22, 12, 54.48);
        observation2[5].ra = Math.toRadians(184.77432);
        observation2[5].dec = Math.toRadians(17.17711);

        //7. postion
        observation2[6].timeMjd = Time.getMjd(2009, 3, 16, 22, 13, 31.33);
        observation2[6].ra = Math.toRadians(184.841);
        observation2[6].dec = Math.toRadians(17.19456);

        //8. postion
        observation2[7].timeMjd = Time.getMjd(2009, 3, 16, 22, 13, 46.33);
        observation2[7].ra = Math.toRadians(184.86472);
        observation2[7].dec = Math.toRadians(17.20125);

        Kepler expKepler3 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation2, 0, 0, expKepler3, "VFMO090316_1");
        getTheOdResults(observation2, 0, 1, expKepler3, "VFMO090316_2");
        getTheOdResults(observation2, 1, 0, expKepler3, "VFMO090316_3");
        getTheOdResults(observation2, 1, 1, expKepler3, "VFMO090316_4");
        */
        /*
         * 1 26609U 00072B   09353.89072269 -.00000314  00000-0  10000-3 0  3307
         * 2 26609 010.2682 350.3994 7976425 268.5232 009.9075 01.25588319 41956
         *
         * AMSAT OSCAR 40
         */
        
        Observation observation3[] = new Observation[10];
        //AGO Modra
        for(int i = 0; i<observation3.length; i++){
            observation3[i] = new Observation();
            observation3[i].lon = Math.toRadians(17.2740);
            observation3[i].lat = Math.toRadians(48.3733);
            observation3[i].alt = 531.1;
        }

        //1. position
        observation3[0].timeMjd = Time.getMjd(2009, 12, 20, 21, 36, 50.11);
        observation3[0].ra = Math.toRadians(63.32935);
        observation3[0].dec = Math.toRadians(5.55002);

        //2. postion
        observation3[1].timeMjd = Time.getMjd(2009, 12, 20, 21, 38, 20.11);
        observation3[1].ra = Math.toRadians(63.42324);
        observation3[1].dec = Math.toRadians(5.56413);

        //3. postion
        observation3[2].timeMjd = Time.getMjd(2009, 12, 20, 21, 39, 15.264);
        observation3[2].ra = Math.toRadians(63.48061);
        observation3[2].dec = Math.toRadians(5.5727);

        //4. postion
        observation3[3].timeMjd = Time.getMjd(2009, 12, 20, 21, 40, 45.12);
        observation3[3].ra = Math.toRadians(63.5738);
        observation3[3].dec = Math.toRadians(5.58653);

        //5. postion
        observation3[4].timeMjd = Time.getMjd(2009, 12, 20, 21, 41, 23.136);
        observation3[4].ra = Math.toRadians(63.61312);
        observation3[4].dec = Math.toRadians(5.59233);

        //6. postion
        observation3[5].timeMjd = Time.getMjd(2009, 12, 20, 21, 42, 52.992);
        observation3[5].ra = Math.toRadians(63.70584);
        observation3[5].dec = Math.toRadians(5.60593);

        //7. postion
        observation3[6].timeMjd = Time.getMjd(2009, 12, 20, 21, 43, 30.144);
        observation3[6].ra = Math.toRadians(63.74408);
        observation3[6].dec = Math.toRadians(5.6115);

        //8. postion
        observation3[7].timeMjd = Time.getMjd(2009, 12, 20, 21, 45, 0);
        observation3[7].ra = Math.toRadians(63.83634);
        observation3[7].dec = Math.toRadians(5.62487);

        //9. postion
        observation3[8].timeMjd = Time.getMjd(2009, 12, 20, 21, 45, 38.016);
        observation3[8].ra = Math.toRadians(63.87527);
        observation3[8].dec = Math.toRadians(5.63048);

        //10. postion
        observation3[9].timeMjd = Time.getMjd(2009, 12, 20, 21, 47, 7.872);
        observation3[9].ra = Math.toRadians(63.96708);
        observation3[9].dec = Math.toRadians(5.64362);

        Kepler expKepler3 = new Kepler();
        Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
                            Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation3, 0, 0, expKepler3, "2000_072B _1");
        getTheOdResults(observation3, 0, 1, expKepler3, "2000_072B _2");
        getTheOdResults(observation3, 1, 0, expKepler3, "2000_072B _3");
        getTheOdResults(observation3, 1, 1, expKepler3, "2000_072B _4");
        
        /*
         * Asteroid 2010 AL30
         */
        /*
        Observation observation4[] = new Observation[10];
        //AGO Modra
        for(int i = 0; i<observation4.length; i++){
            observation4[i] = new Observation();
            observation4[i].lon = Math.toRadians(17.2740);
            observation4[i].lat = Math.toRadians(48.3733);
            observation4[i].alt = 531.1;
        }

        //1. position
        observation4[0].timeMjd = Time.getMjd(2010, 1, 13, 13, 20, 0.0);
        observation4[0].ra = Math.toRadians(20.44937);
        observation4[0].dec = Math.toRadians(9.28192 );

        //2. postion
        observation4[1].timeMjd = Time.getMjd(2010, 1, 13, 13, 21, 0.0);
        observation4[1].ra = Math.toRadians(20.17801);
        observation4[1].dec = Math.toRadians(9.23617);

        //3. postion
        observation4[2].timeMjd = Time.getMjd(2010, 1, 13, 13, 22, 0.0);
        observation4[2].ra = Math.toRadians(19.90698);
        observation4[2].dec = Math.toRadians(9.19026);

        //4. postion
        observation4[3].timeMjd = Time.getMjd(2010, 1, 13, 13, 23, 0.0);
        observation4[3].ra = Math.toRadians(19.63629);
        observation4[3].dec = Math.toRadians(9.14420);

        //5. postion
        observation4[4].timeMjd = Time.getMjd(2010, 1, 13, 13, 24, 0.0);
        observation4[4].ra = Math.toRadians(19.36597);
        observation4[4].dec = Math.toRadians(9.09798);

        //6. postion
        observation4[5].timeMjd = Time.getMjd(2010, 1, 13, 13, 25, 0.0);
        observation4[5].ra = Math.toRadians(19.09601);
        observation4[5].dec = Math.toRadians(9.05162);

        //7. postion
        observation4[6].timeMjd = Time.getMjd(2010, 1, 13, 13, 26, 0.0);
        observation4[6].ra = Math.toRadians(18.82643);
        observation4[6].dec = Math.toRadians(9.00511);

        //8. postion
        observation4[7].timeMjd = Time.getMjd(2010, 1, 13, 13, 27, 0.0);
        observation4[7].ra = Math.toRadians(18.55725);
        observation4[7].dec = Math.toRadians(8.95846);

        //9. postion
        observation4[8].timeMjd = Time.getMjd(2010, 1, 13, 13, 28, 0.0);
        observation4[8].ra = Math.toRadians(18.28847);
        observation4[8].dec = Math.toRadians(8.91167);

        //10. postion
        observation4[9].timeMjd = Time.getMjd(2010, 1, 13, 13, 29, 0.0);
        observation4[9].ra = Math.toRadians(18.02010);
        observation4[9].dec = Math.toRadians(8.86475);

        Kepler expKepler4 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation4, 0, 0, expKepler2, "ast2010 AL30 _1");
        getTheOdResults(observation4, 0, 1, expKepler2, "ast2010 AL30 _2");
        getTheOdResults(observation4, 1, 0, expKepler2, "ast2010 AL30 _3");
        getTheOdResults(observation4, 1, 1, expKepler2, "ast2010 AL30 _4");
        */
        
        /**
          * Test bolid 03/04/2009
        */
        /*
        Observation observation4[] = new Observation[8];
        //AGO Modra
        for(int i = 0; i<observation4.length; i++){
            observation4[i] = new Observation();
            observation4[i].lon = Math.toRadians(17.2740);
            observation4[i].lat = Math.toRadians(48.3733);
            observation4[i].alt = 531.1;
        }

        Vector vectorHelp = new Vector(2);

        //1st position
        double timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 23.911));
        double azi = Math.toRadians(282.4642);
        double elev = Math.toRadians(63.0141);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[0].timeMjd = timeMjd;
        observation4[0].ra = vectorHelp.v[0];
        observation4[0].dec = vectorHelp.v[1];

        //2nd position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 23.951));
        azi = Math.toRadians(282.5062);
        elev = Math.toRadians(63.2516);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[1].timeMjd = timeMjd;
        observation4[1].ra = vectorHelp.v[0];
        observation4[1].dec = vectorHelp.v[1];

        //3rd position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 25.011));
        azi = Math.toRadians(284.3747);
        elev = Math.toRadians(70.9258);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[2].timeMjd = timeMjd;
        observation4[2].ra = vectorHelp.v[0];
        observation4[2].dec = vectorHelp.v[1];

        //4th position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 26.011));
        azi = Math.toRadians(288.4666);
        elev = Math.toRadians(79.2204);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[3].timeMjd = timeMjd;
        observation4[3].ra = vectorHelp.v[0];
        observation4[3].dec = vectorHelp.v[1];

        //5th position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 27,011));
        azi = Math.toRadians(324.8397);
        elev = Math.toRadians(87.7089);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[4].timeMjd = timeMjd;
        observation4[4].ra = vectorHelp.v[0];
        observation4[4].dec = vectorHelp.v[1];

        //6th position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 27.991));
        azi = Math.toRadians(87.3836);
        elev = Math.toRadians(82.4090);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[5].timeMjd = timeMjd;
        observation4[5].ra = vectorHelp.v[0];
        observation4[5].dec = vectorHelp.v[1];

        //7th position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 28.971));
        azi = Math.toRadians(94.0196);
        elev = Math.toRadians(73.3590);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[5].timeMjd = timeMjd;
        observation4[5].ra = vectorHelp.v[0];
        observation4[5].dec = vectorHelp.v[1];

        //8th position
        timeMjd = Time.getMjd(new Time(2009, 04, 03, 1, 23, 29.931));
        azi = Math.toRadians(95.8131);
        elev = Math.toRadians(65.782);
        vectorHelp = Transformation.fromHorizToEquatCoord_New(azi, elev, observation4[0].lon, observation4[0].lat, timeMjd);
        observation4[5].timeMjd = timeMjd;
        observation4[5].ra = vectorHelp.v[0];
        observation4[5].dec = vectorHelp.v[1];

        Kepler expKepler2 = new Kepler();
        getTheOdResults(observation4, 0, 0, expKepler2, "bolid090403_1");
        getTheOdResults(observation4, 0, 1, expKepler2, "bolid090403_2");
        getTheOdResults(observation4, 1, 0, expKepler2, "bolid090403_3");
        getTheOdResults(observation4, 1, 1, expKepler2, "bolid090403_4");
        */
        
        /*
         * VFMO100422
         */
         /*
        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }

        //1. position
        observation2[0].timeMjd = Time.getMjd(2010,04,22,19,55,57.38);
        observation2[0].ra = Math.toRadians(217.13628);
        observation2[0].dec = Math.toRadians(19.49119);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010,04,22,19,56,34.70);
        observation2[1].ra = Math.toRadians(217.21872);
        observation2[1].dec = Math.toRadians(19.51394);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010,04,22,19,58,04.70);
        observation2[2].ra = Math.toRadians(217.42942);
        observation2[2].dec = Math.toRadians(19.57269);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010,04,22,19,58,42.13);
        observation2[3].ra = Math.toRadians(217.51224);
        observation2[3].dec = Math.toRadians(19.59497);

        Kepler expKepler2 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation2, 0, 0, expKepler2, "VFMO100422_1");
        getTheOdResults(observation2, 0, 1, expKepler2, "VFMO100422_2");
        getTheOdResults(observation2, 1, 0, expKepler2, "VFMO100422_3");
        getTheOdResults(observation2, 1, 1, expKepler2, "VFMO100422_4");
        */

        /*
         *  Satellite INTELSAT 4A-F
         *  1 10778U 78035A   09115.32846199  .00000124  00000-0  10000-3 0  4338
         *  2 10778 013.8310 011.3753 0002637 172.4122 187.5454 01.00366373 59465
         */
         /*
        Observation observation2[] = new Observation[3];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }

        //1. position
        observation2[0].timeMjd = Time.getMjd(2009,04,27,20,54,51.32);
        observation2[0].ra = Math.toRadians(212.04662);
        observation2[0].dec = Math.toRadians(-11.81692);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2009,04,27,20,55,29.17);
        observation2[1].ra = Math.toRadians(212.19754);
        observation2[1].dec = Math.toRadians(-11.85361);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2009,04,27,20,55,48.17);
        observation2[2].ra = Math.toRadians(212.27916);
        observation2[2].dec = Math.toRadians(-11.87497);

        Kepler expKepler2 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation2, 0, 0, expKepler2, "VFMO090427_1");
        getTheOdResults(observation2, 0, 1, expKepler2, "VFMO090427_2");
        getTheOdResults(observation2, 1, 0, expKepler2, "VFMO090427_3");
        getTheOdResults(observation2, 1, 1, expKepler2, "VFMO090427_4");
        */

        /*
         *  CZ-3B R/B
         *  1 25405U 98044B   09089.84099089 -.00000138 +00000-0 +10000-3 0 08181
         *  2 25405 018.4429 107.8610 7154608 242.5817 029.8295 02.25512590088191
         */
        /*
        Observation observation2[] = new Observation[3];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }

        //1. position
        observation2[0].timeMjd = Time.getMjd(2009,03,30,22,59,47.76);
        observation2[0].ra = Math.toRadians(151.19554);
        observation2[0].dec = Math.toRadians(7.65944);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2009,03,30,23,0,47.76);
        observation2[1].ra = Math.toRadians(151.34846);
        observation2[1].dec = Math.toRadians(7.70511);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2009,03,30,23,01,49.04);
        observation2[2].ra = Math.toRadians(151.49992);
        observation2[2].dec = Math.toRadians(7.75222);

        Kepler expKepler2 = new Kepler();
        //Kepler expKepler = new Kepler(Math.toRadians(308.16), Math.toRadians(197.29),
        //                    Math.toRadians(82.6), 0.0041,  7820.866, Math.toRadians(51.57));
        getTheOdResults(observation2, 0, 0, expKepler2, "VFMO090330_1");
        getTheOdResults(observation2, 0, 1, expKepler2, "VFMO090330_2");
        getTheOdResults(observation2, 1, 0, expKepler2, "VFMO090330_3");
        getTheOdResults(observation2, 1, 1, expKepler2, "VFMO090330_4");
        
        */
    }

    /**
     * getAstroPosFromFile()
     *
     * Method to read the JPL ephemerides file and to fill arrray of Observation[] object.
     *
     * IN:
     *      String file - source file
     *
     * OUT:
     *      Observation observation[]
     */
    //public Observation[] getAstroPosFromFile(String file){
    //    Observation observation[];

    //    return observation;
    //}

    /**
     * getDataFromString()
     *
     * Method to get data from String - data --> Time, R.A., dec --> Observation.
     *
     * IN:
     * String dataString    - string with data - time, ra, dec
     * int intArray[] - arrray with positions of given datas in string
     *               0 - yearPos_1, 1 - yearPos_2, 2 - monthPos_1, 3 - monthPos_2,
     *               4 - dayPos_1,  5 - dayPos_2,  6 - hourPos_1,  7 -  hourPos_2,
     *               8 - minPos_1,  9 -  minPos_2, 10- secPos_1,  11 - secPos_2,
     *               12 - raPos_1, 13 -  raPos_2,  14 - decPos_1, 15 - decPos_2;
     *
     * OUT:
     *  Observation observation;
     */
    public static Observation getDataFromString(String dataString, int intArray[], boolean isRaInDeg){
        Observation observation = new Observation();
        Time time = new Time();

        try{
            time.year = (int)Double.parseDouble(dataString.substring(intArray[0], intArray[1]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong year input!");
        }
        try{
            time.month = (int)Double.parseDouble(dataString.substring(intArray[2], intArray[3]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong month input!");
        }
        try{
            time.day = (int)Double.parseDouble(dataString.substring(intArray[4], intArray[5]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong day input!");
        }
        try{
            time.hour = (int)Double.parseDouble(dataString.substring(intArray[6], intArray[7]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong hour input!");
        }
        try{
            time.min = (int)Double.parseDouble(dataString.substring(intArray[8], intArray[9]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong minute input!");
        }
        try{
            time.sec = Double.parseDouble(dataString.substring(intArray[10], intArray[11]));
        }
        catch(NumberFormatException nfe){
            time.sec = 0;
            System.out.println("Wrong seconds input!");
        }
        try{
            observation.ra = Double.parseDouble(dataString.substring(intArray[12], intArray[13]));
            if(isRaInDeg==false) observation.ra = observation.ra*15;
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong R.A. input!");
        }
        try{
            observation.dec = Double.parseDouble(dataString.substring(intArray[14], intArray[15]));
        }
        catch(NumberFormatException nfe){
            System.out.println("Wrong declination input!");
        }
        //System.out.println(time.year);
        //System.out.println(time.month);
        //System.out.println(time.day);
        //System.out.println(time.hour);
        //System.out.println(time.min);
        //System.out.println(time.sec);
        observation.timeMjd = time.getMjd(time);

        return observation;
    }

}
