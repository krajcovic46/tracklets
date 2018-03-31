/*
 * Class to comunicate between application and data files.
 */

package com.skrajcovic.orbitdetermination;


import com.skrajcovic.orbitdetermination.compute.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 *
 * @author Jiri Silha - 23.10.2009
 */
public class WritingODData {


    /**
     *
     * createFile()
     * Method to create non existing file RESULTS.TXT. All OD results will be in the directory
     * ODResults/dirName/results.txt
     *
     * INPUT:
     *  String dirName - directory name
     */
     public static void createFile(String dirName){
         //create directory
         new FileDownload().createDirectory("ODResults");
         //create subdirectory
         new FileDownload().createDirectory("ODResults\\" + dirName);
         //create file results.txt if does not exist
         File file = new File("ODResults\\" + dirName + "\\" + "results.txt");
         if(!file.exists()){
            try{
                file.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }
         //System.out.println(file.length());
         
         File fileGood = new File("ODResults\\" + dirName + "\\" + "resultsGood.txt");
         if(!fileGood.exists()){
            try{
                fileGood.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }

         File fileTLE = new File("ODResults\\" + dirName + "\\" + "resultsTLE.txt");
         if(!fileTLE.exists()){
            try{
                fileTLE.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }
         //System.out.println(fileGood.length());

         File fileGoodHelio = new File("ODResults\\" + dirName + "\\" + "resultsGoodHelio.txt");
         if(!fileGoodHelio.exists()){
            try{
                //APRIL 2012
                //fileGoodHelio.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }
         //System.out.println(fileGoodHelio.length());
         File fileGoodGeo = new File("ODResults\\" + dirName + "\\" + "resultsGoodGeo.txt");
         if(!fileGoodGeo.exists()){
            try{
                //APRIL 2012
                //fileGoodGeo.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }
         //System.out.println(fileGoodHelio.length());
         
         File fileVectors = new File("ODResults\\" + dirName + "\\" + "resultsVec.txt");
         if(!fileVectors.exists()){
            try{
                fileVectors.createNewFile();
            }
            catch(Exception e){
                System.out.println("Warning! " + e.getMessage());
            }
         }
     }

     /**
      * fillResultsFileBasic()
      *
      * Method to fill empty results.txt file file with basic informations -
      * column names, expected keplerian elements etc.
      *
      * INPUT:
      *     String dirName - directory name
      *     Kepler expKepler - expected orbittal elements
      * 
      */
     public static void fillResultsFileBasic(String dirName, Kepler expKepler) {
        File file = new File("ODResults\\" + dirName + "\\" + "results.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789\n";
            stringComments = stringComments + "Note  ---->  Method_1: 0-Bucerius, 1-Escobal, Method_2: 0-Gauss, 1-Escobal \n\n";
            stringComments = stringComments + getRightString("eph", 15);
            stringComments = stringComments + "  " + getRightString("a [km]", 11);
            stringComments = stringComments + "  " + getRightString("e",11);
            stringComments = stringComments + "  " + getRightString("i [deg]",11);
            stringComments = stringComments + "  " + getRightString("omega [deg]",11);
            stringComments = stringComments + "  " + getRightString("Omega [deg]",11);
            stringComments = stringComments + "  " + getRightString("Ma(t) [deg]",11);
            stringComments = stringComments + "  " + getRightString("Method_1",11);
            stringComments = stringComments + "  " + getRightString("Method_2",11);
            stringComments = stringComments + "  " + getRightString("q",11);
            stringComments = stringComments + "  " + getRightString("Epoch MJD [days]",16);
            stringComments = stringComments + "\n";
            
            //write expected kepler elements
            stringComments = stringComments + getRightString("Expected ", 15);
            stringComments = stringComments + "  " + getRightString(expKepler.a + "", 11);
            stringComments = stringComments + "  " + getRightString((float)((int)(expKepler.e*10e9))/10e9 + "",11);
            stringComments = stringComments + "  " + getRightString(Math.toDegrees(expKepler.incl) + "",11);
            stringComments = stringComments + "  " + getRightString(Math.toDegrees(expKepler.omega) + "",11);
            stringComments = stringComments + "  " + getRightString(Math.toDegrees(expKepler.Omega) + "",11);
            stringComments = stringComments + "  " + getRightString(Math.toDegrees(expKepler.M) + "",11);
            stringComments = stringComments + "  " + getRightString("-",11);
            stringComments = stringComments + "  " + getRightString("-",11);
            stringComments = stringComments + "  " + getRightString((expKepler.a*(1-expKepler.e)) + "",11);
            stringComments = stringComments + "\n";

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            raFile.writeBytes(stringComments);
            raFile.close();

            //System.out.println(stringComments);
            //out.write(stringComments + "\n");
            //out.close();
            //System.out.println("Data written to ODResults/" + dirName + "/results.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
            System.out.println("Data NOT written to ODResults/" + dirName + "/results.txt!");
        }
     }

     /**
      * fillResultsFile()
      *
      * Method to fill not empty results.txt file with orbital elements  -
      * ephemerides combinations, a, e, i , omega, Omega, M,etc.
      *
      * INPUT:
      *     String dirName - directory name
      *     Kepler kepler - computed orbital elements
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      *
      */
     public static void fillResultsFile(String dirName, Kepler kepler, Vector ephComb,
                                        int method_1, int method_2){
        File file = new File("ODResults\\" + dirName + "\\" + "results.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "\n";
            //write computeted orbital elements
            stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(kepler.a/1000 + "", 11);
            stringComments = stringComments + "\t" + getRightString((float)((int)(kepler.e*10e9))/10e9 + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.incl) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.Omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.M) + "",11);
            stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.a/1000*(1 - kepler.e)) + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.epochMJD + ""),16);
            stringComments = stringComments + "\n";

            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            raFile.writeBytes(stringComments);
            raFile.close();
            freader.close();

            //go to the end of the file and then write String
            //boolean doLoop = true;
            //while(doLoop){
            //    String line = breader.readLine();
            //    System.out.println(line);
                //if(line.equals("null")) doLoop = false;
            //}

            //out.write(stringComments + "\n");
            //out.close();
            //System.out.println("Data written in to ODResults/" + dirName + "/results.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
            System.out.println("Data NOT written to ODResults/" + dirName + "/results.txt!");
        }
     }

     /**
      * fillResultsFileGood()
      *
      * Method to fill not empty resultsGood.txt file with orbital elements  -
      * ephemerides combinations, a, e, i , omega, Omega, M,etc.
      *
      * INPUT:
      *     String dirName - directory name
      *     Kepler kepler - computed orbital elements
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      *
      */
     public static void fillResultsFileGood(String dirName, Kepler kepler, Vector ephComb,
                                        int method_1, int method_2){
        File file = new File("ODResults\\" + dirName + "\\" + "resultsGood.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "";
            //write computeted orbital elements
            stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(kepler.a/1000 + "", 11);
            stringComments = stringComments + "\t" + getRightString((float)((int)(kepler.e*10e9))/10e9 + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.incl) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.Omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.M) + "",11);
            stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.a/1000*(1 - kepler.e)) + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.epochMJD + ""),16);
            stringComments = stringComments + "\n";

            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            if(((kepler.a/1000*(1 - kepler.e))>3000)&&(kepler.e<1)){
                raFile.writeBytes(stringComments);
            }
            raFile.close();
            freader.close();
            //go to the end of the file and then write String
            //boolean doLoop = true;
            //while(doLoop){
            //    String line = breader.readLine();
            //    System.out.println(line);
                //if(line.equals("null")) doLoop = false;
            //}

            //out.write(stringComments + "\n");
            //out.close();
            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGood.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
            System.out.println("Data NOT written to ODResults/" + dirName + "/resultsGood.txt!");
        }
     }

     /**
      * fillResultsFileTLE()
      *
      * Method to fill empty resultsTLE.txt file with TLE data  -
      * ephemerides combinations, a, e, i , omega, Omega, M,etc.
      *
      * INPUT:
      *     String dirName - directory name
      *     Kepler kepler - computed orbital elements
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      *     int noGoodResults - number of good results for the NORAD number
      */
     public static void fillResultsFileTLE(String dirName, Kepler kepler, Vector ephComb,
                                        int method_1, int method_2, int noGoodResults){
        File file = new File("ODResults\\" + dirName + "\\" + "resultsTLE.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "";
            //write computeted orbital elements
            //stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            //stringComments = stringComments + "\t" + getRightString(kepler.a/1000 + "", 11);
            //stringComments = stringComments + "\t" + getRightString((float)((int)(kepler.e*10e9))/10e9 + "",11);
            //stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.incl) + "",11);
            //stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.omega) + "",11);
            //stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.Omega) + "",11);
            //stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.M) + "",11);
            //stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            //stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            //stringComments = stringComments + "\t" + getRightString((kepler.a/1000*(1 - kepler.e)) + "",11);
            //stringComments = stringComments + "\t" + getRightString((kepler.epochMJD + ""),16);
            //stringComments = stringComments + "\n";

            stringComments = kepler.getTLE(getRightString("OD-"+(int)ephComb.v[0] + "-" + (int)ephComb.v[1] + "-" + (int)ephComb.v[2] +(" -" +method_1+method_2),12),
                    60000+noGoodResults,"00000A", kepler.epochMJD, Math.toDegrees(kepler.incl), Math.toDegrees(kepler.Omega),
                    Math.toDegrees(kepler.omega), Math.toDegrees(kepler.M), kepler.a, kepler.e, 85848-3)+ "\n";
            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            if(((kepler.a/1000*(1 - kepler.e))>3000)&&(kepler.e<1)){
                raFile.writeBytes(stringComments);
            }
            raFile.close();
            freader.close();
            //go to the end of the file and then write String
            //boolean doLoop = true;
            //while(doLoop){
            //    String line = breader.readLine();
            //    System.out.println(line);
                //if(line.equals("null")) doLoop = false;
            //}

            //out.write(stringComments + "\n");
            //out.close();
            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGood.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
            System.out.println("Data NOT written to ODResults/" + dirName + "/resultsGood.txt!");
        }
     }
     
     /**
      * fillResultsFileVectors()
      *
      * Method to fill not empty resultsGood.txt file with orbital elements  -
      * ephemerides combinations, a, e, i , omega, Omega, M,etc.
      *
      * INPUT:
      *     String dirName - directory name
      *     Kepler kepler - computed orbital elements
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      *
      */
     public static void fillResultsFileVectors(String dirName, Vector ephComb, Vector vectorsForOd[], Observation observation2[]){
        File file = new File("ODResults\\" + dirName + "\\" + "resultsVec.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "";
            //write computeted orbital elements
            stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(observation2[0].timeMjd + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[0].v[0] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[0].v[1] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[0].v[2] + "", 15);
            stringComments = stringComments + "\t" + getRightString(observation2[1].timeMjd + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[1].v[0] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[1].v[1] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[1].v[2] + "", 15);
            stringComments = stringComments + "\t" + getRightString(observation2[2].timeMjd + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[2].v[0] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[2].v[1] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[2].v[2] + "", 15);
            stringComments = stringComments + "\t" + getRightString(observation2[1].timeMjd + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[3].v[0] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[3].v[1] + "", 15);
            stringComments = stringComments + "\t" + getRightString(vectorsForOd[3].v[2] + "", 15);
            stringComments = stringComments + "\n";
            /*stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(kepler.a/1000 + "", 11);
            stringComments = stringComments + "\t" + getRightString((float)((int)(kepler.e*10e9))/10e9 + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.incl) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.Omega) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(kepler.M) + "",11);
            stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.a/1000*(1 - kepler.e)) + "",11);
            stringComments = stringComments + "\t" + getRightString((kepler.epochMJD + ""),16);
            stringComments = stringComments + "\n";
            */
            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            //if(((kepler.a/1000*(1 - kepler.e))>3000)&&(kepler.e<1)){
            raFile.writeBytes(stringComments);
            //}
            raFile.close();
            freader.close();
            //go to the end of the file and then write String
            //boolean doLoop = true;
            //while(doLoop){
            //    String line = breader.readLine();
            //    System.out.println(line);
                //if(line.equals("null")) doLoop = false;
            //}

            //out.write(stringComments + "\n");
            //out.close();
            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGood.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
            System.out.println("Data NOT written to ODResults/" + dirName + "/resultsGood.txt!");
        }
     }
     
     /**
      * fillResultsFileGoodHelio()
      *
      * Method to fill not empty results.txt file with heliocentric orbital elements  -
      * ephemerides combinations, i, Omega, distance1, distance2, distance3
      *
      * INPUT:
      *     String dirName - directory name
      *     Vector geocPosVec_1 - 1st geocentric position vector [m]
      *     Vector geocPosVec_2 - 2nd geocentric position vector [m]
      *     double mjd_1        - get time for 1st position  [days]
      *     double mjd_2        - get time for 1st position  [days]
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      */
     public static void fillResultsFileGoodHelio(String dirName, Vector geocPosVec_1, Vector geocPosVec_2,
                                                double mjd_1, double mjd_2, Vector ephComb, int method_1, int method_2){
         //compute heliocentric elements
         double elements[] = Transformation.getEclipticInclAndNodeFromGeocPosVectors(geocPosVec_1,geocPosVec_2,mjd_1,mjd_2);
         if(elements[0]<0) elements[0] = elements[0]+Math.PI;
         if(elements[1]<0) elements[1] = elements[1]+Math.PI*2;
         //distance above Earth surface
         Geodetic geodPosVec_1 = Transformation.fromGeocToGeod(geocPosVec_1, Constants.R_Earth, Constants.f_Earth);
         Geodetic geodPosVec_2 = Transformation.fromGeocToGeod(geocPosVec_2, Constants.R_Earth, Constants.f_Earth);
         
        File file = new File("ODResults\\" + dirName + "\\" + "resultsGoodHelio.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "";
            //write computeted orbital elements
            stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(elements[0]) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(elements[1]) + "",11);
            stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            stringComments = stringComments + "\t" + getRightString((geodPosVec_1.altitude/1000) + "",11);
            stringComments = stringComments + "\t" + getRightString((geodPosVec_2.altitude/1000) + "",11);
            stringComments = stringComments + "\n";

            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            if((geodPosVec_1.altitude>0)&&(geodPosVec_2.altitude>0)&&(geodPosVec_1.altitude>geodPosVec_2.altitude)){
                raFile.writeBytes(stringComments);
            }
            raFile.close();
            freader.close();

            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGoodHelio.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
        }
     }

     /**
      * fillResultsFileGoodGeo()
      *
      * Method to fill not empty results.txt file with heliocentric orbital elements  -
      * ephemerides combinations, i, Omega, distance1, distance2, distance3
      *
      * INPUT:
      *     String dirName - directory name
      *     Vector geocPosVec_1 - 1st geocentric position vector [m]
      *     Vector geocPosVec_2 - 2nd geocentric position vector [m]
      *     double mjd_1        - get time for 1st position  [days]
      *     double mjd_2        - get time for 1st position  [days]
      *     Vector ephComb(3) - 3 ephemerides, that were used to compute elements
      *     int method_1 - which method was used to compute vectors
      *     int method_2 - which method was used to compute orbital elements from vectors
      */
     public static void fillResultsFileGoodGeo(String dirName, Vector geocPosVec_1, Vector geocPosVec_2,
                                                double mjd_1, double mjd_2, Vector ephComb, int method_1, int method_2){
         //compute heliocentric elements
         double elements[] = Transformation.getEquatorialInclAndNodeFromGeocPosVectors(geocPosVec_1,geocPosVec_2,mjd_1,mjd_2);
         if(elements[0]<0) elements[0] = elements[0]+Math.PI;
         if(elements[1]<0) elements[1] = elements[1]+Math.PI*2;
         //distance above Earth surface
         Geodetic geodPosVec_1 = Transformation.fromGeocToGeod(geocPosVec_1, Constants.R_Earth, Constants.f_Earth);
         Geodetic geodPosVec_2 = Transformation.fromGeocToGeod(geocPosVec_2, Constants.R_Earth, Constants.f_Earth);

        File file = new File("ODResults\\" + dirName + "\\" + "resultsGoodGeo.txt");
        //open file
        try{
            FileReader freader = new FileReader(file);
            //FileWriter fstream = new FileWriter(file);
            //BufferedWriter out = new BufferedWriter(fstream);
            //BufferedReader breader = new BufferedReader(freader);
            RandomAccessFile raFile = new RandomAccessFile(file, "rw");

            //write the basic data to results.txt
            String stringComments = "";
            stringComments = "";
            //write computeted orbital elements
            stringComments = stringComments + getRightString((int)ephComb.v[0] + " - " + (int)ephComb.v[1] + " - " + (int)ephComb.v[2], 15);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(elements[0]) + "",11);
            stringComments = stringComments + "\t" + getRightString(Math.toDegrees(elements[1]) + "",11);
            stringComments = stringComments + "\t" + getRightString(method_1 + "",11);
            stringComments = stringComments + "\t" + getRightString(method_2 + "",11);
            stringComments = stringComments + "\t" + getRightString((geodPosVec_1.altitude/1000) + "",11);
            stringComments = stringComments + "\t" + getRightString((geodPosVec_2.altitude/1000) + "",11);
            stringComments = stringComments + "\n";

            //System.out.println(stringComments);

            // Set write pointer to the end of the file
            raFile.seek(raFile.length());
            if((geodPosVec_1.altitude>0)&&(geodPosVec_2.altitude>0)&&(geodPosVec_1.altitude>geodPosVec_2.altitude)){
                raFile.writeBytes(stringComments);
            }
            raFile.close();
            freader.close();

            //System.out.println("Data written in to ODResults/" + dirName + "/resultsGoodGeo.txt!");
        }
        catch(IOException ioe){
            System.out.println("Warning! " + ioe.getMessage());
        }
     }
     
     /**
      * isFileEmpty()
      *
      * Method to check wether or not is file empty.
      *
      * INPUT:
      * String dirName
      *
      * OUTPUT:
      * boolean isFileEmpty
      */
     public static boolean isFileEmpty(String dirName){
         File file = new File("ODResults\\" + dirName + "\\" + "results.txt");
         boolean isFileEmpty = false;
         if(file.length() == 0) isFileEmpty = true;
         else isFileEmpty = false;
         return isFileEmpty;
     }

     /**
     * getRightString()
     *
     * Method to get the given type of string.
     *
     * IN:
     *  String stringToRepair
     *  int lengthOfString
     *
     * OUT:
     *  String rightString
     *
     */

    public static String getRightString(String stringToRepair, int lengthOfString){
        String rightString = "";
        int maxLength = 100;
        if(lengthOfString > maxLength) lengthOfString = maxLength;
        rightString = (stringToRepair + "                                                                      " +
                "                                        ").substring(0, lengthOfString);

        return rightString;
    }

     //test of methods
     public static void main(String args[]){
         String dirName = "test";
         Kepler kepler = new Kepler(Math.toRadians(10),Math.toRadians(100),Math.toRadians(90), 0.5, 12000, Math.toRadians(180));
         createFile(dirName);
         if(isFileEmpty(dirName)) fillResultsFileBasic(dirName, kepler);
         Vector vector = new Vector(3);
         vector.v[0] = 10;
         vector.v[1] = 13;
         vector.v[2] = 25;
         fillResultsFile(dirName, kepler, vector, 0, 0);
     }
}
