/*
 * ReadFile.java
 *
 * Created on February 28, 2008, 1:25 PM
 *
 * SK: Trieda sluzi na vycitanie udajov zo TLE suborov a satelitnej spravy
 * EN: Class for reading files. (TLE and satellite report)
 */

package InOut;

import java.io.*;
import java.util.*;
import SatEph.*;

/**
 *
 * @author Jiri Silha
 */
public class ReadFile {
    
    //Variables
    static BufferedReader file;    //read the file
    
    static StringTokenizer st;      //to separete the line for smaller parts

    //on which OS is application running
    String osName = System.getProperty("os.name");
    
    /**
     *Method: readFileTle()
     *
     *Method is to fill constructor SatelliteArray from TLE file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\file.txt")
     *      Satellite Array
     *      boolean isKelso    - which data base is source data base         
     *
     *OUT:  constructor filled SatelliteArray
     *
     */
    
    public SatelliteArray readFileTle(String filePos, SatelliteArray sa, boolean isKelso) {
        
        //if user is in Windows, read file like X:\\DIRECTORY\\file.txt
        if(osName.startsWith("Windows")){} //do nothing
        //UNIX, or Linux OS
        else{
            String helpString = "\\";
            while(filePos.indexOf("\\") != -1){
                int slashPos = filePos.indexOf(helpString);
                filePos = filePos.substring(0,slashPos) + "//" + filePos.substring(slashPos+1, filePos.length());
            }
        }        
        
        System.out.println("Reading TLE from " + filePos);
        //actual line
        String line = "";
        
        //auxiliary variable
        String line1 = "";
        String line2 = "";
        String noradString = "";
        int norad = 0;
        
        //number of lines
        int no = 0;
        
        //do cycle
        boolean doCycle = true;
        
        //array to fill        
        //FillConstructor methods
        FillConstructor fc = new FillConstructor();
       
        //try read file
        try{
            //read file from specify position
            file = new BufferedReader(new FileReader(filePos));
            //name of the object
            String objectName = "Unknown name";

            //read line
            while(doCycle){
                //beware of ":"
                line = file.readLine();
                line = line + "";
                //System.out.println(line);
                if(line.indexOf(":") != -1) {
                    String help = line.substring(line.indexOf(":") + 1);
                    line = line.substring(0, line.indexOf(":"));
                    line = line + " " + help;
                    //System.out.println("1: " + line);
                    //System.out.println("position " + line.indexOf(":"));
                    line.replaceAll(": ", "0 ");
                    //System.out.println("2: " + line);
                }
                //end of file
                if(line.equals("null")) doCycle = false;
                else{
                    //count no. of lines
                    no++;
                    
                    //condition
                    //thera are 2 types of TLE dates - with name, and without
                    if((line.substring(0,2).equals("1 "))||(line.substring(0,2).equals("2 "))){
                       //unpaired line (1st line of TLE datas)
                        if(no%2!=0) line1 = line;
                        //paired line (2nd line of TLE datas)
                        else {
                            line2 = line;
                            //now we can fill the constructor Satellite and SatelliteArray
                            //what object it is (NORAD cat No)
                            noradString = line1.substring(2,7);
                            norad = (int)Double.parseDouble(noradString);
                            //searching for lagrest catalog no.
                            if(norad > sa.largestNo) sa.largestNo = norad;
                            
                            //get name, when 3 line TLE format is using
                            if(sa.array[norad].name.equalsIgnoreCase("No name")) sa.array[norad].name = objectName;
                            //System.out.println("norad " + norad);
                            //System.out.println("name " + sa.array[norad].name);
                            //System.out.println(objectName + " test");
                            
                            //Fill the constructor
                            sa.array[norad] = fc.fillConstructorTle(line1,line2,sa.array[norad]);
                            sa.array[norad].catalogNo = norad;
                            sa.array[norad].line1 = line1;
                            sa.array[norad].line2 = line2;
                            
                            //if object has got TLE, is on orbit, till the Sat report says smth else
                            sa.array[norad].isOnOrbit = true;
                            //if(norad==96082) System.out.println("sa.array[norad].isOnOrbit " + sa.array[norad].isOnOrbit);
                            int selectedObj = norad;
                            //test - looking for D criterion for artifitial objects
                            //if(norad == selectedObj){
                                //System.out.println(selectedObj);
                                //System.out.println(sa.array[selectedObj].epoch - 9000);
                                //System.out.println(sa.array[selectedObj].ecc);
                                //System.out.println(sa.array[selectedObj].meanMot);
                                //semi major axis
                                //System.out.println(sa.array[selectedObj].smAxis/1000);
                                //mean altitude
                                //System.out.println(sa.array[selectedObj].smAxis/1000 - 6378.137);                         
                                //System.out.println(sa.array[selectedObj].incl);
                                //System.out.println(sa.array[selectedObj].argument);
                                //System.out.println(sa.array[selectedObj].node);
                                //System.out.println(sa.array[selectedObj].smAxis*(1 - sa.array[selectedObj].ecc)/1000);
                                //System.out.println(sa.array[selectedObj].smAxis*(1 + sa.array[selectedObj].ecc)/1000);
                            //}
                            
                        }
                    }
                    //if the line is with name 
                    else {
                        no = no - 1;
                        if(line.substring(0,2).equalsIgnoreCase("0 ")){
                            objectName = line.substring(2);
                        }
                        else {
                            objectName = line;
                        }
                    }
                }
            }
            file.close();
        }
        /*
        //if file can not be open, io exception is used
        catch(IOException e) {
            System.out.println("Error: can not read Tle file!");
            System.out.println(e);
        }
        catch (StringIndexOutOfBoundsException x)
        {
            System.out.println("Error: can not read Tle file!");
            System.out.println (x.getClass());
        }
        
        catch (NumberFormatException nfe){
            System.out.println("Error: can not read Tle file!");
            System.out.println (nfe);
        }
        */
        catch (Exception e){
            System.out.println("Error: can not read Tle file!");
            System.out.println (e);
        }
        
        //return filled constructor SatelliteArray
        //is from Kelso or DoD databaze
        if(isKelso)sa.tleKelsoAmount = no/2;
        else sa.tleAmount = no/2;

        int noOfTleObjects = 0;
        for(int i = 0; i <= sa.largestNo; i++){
            ///if(sa.array[i].epoch != 0) noOfTleObjects++;
            if(sa.array[i].epochMjd != 0) noOfTleObjects++;
        }
        sa.tleAmount = noOfTleObjects;
             
        //System.out.println("tleKelsoAmount " + sa.tleKelsoAmount);
        System.out.println("Done!");
        
        return sa;
    }
    
   
     /**
     *Method: readFileSatRep()
     *
     *Method is to fill constructor SatelliteArray from satellite report file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\file.txt")
     *      Satellite array    
     *
     *OUT:  constructor filled SatelliteArray
     *
     */
    
     public SatelliteArray readFileSatRep(String filePos, SatelliteArray sa) {
        
        System.out.println("Reading report from " + filePos);

         //if user is in Windows, read file like X:\\DIRECTORY\\file.txt
        if(osName.startsWith("Windows")){} //do nothing
        //UNIX, or Linux OS
        else{   //replace all backslashes
            String helpString = "\\";
            while(filePos.indexOf("\\") != -1){
                int slashPos = filePos.indexOf(helpString);
                filePos = filePos.substring(0,slashPos) + "//" + filePos.substring(slashPos+1, filePos.length());
            }
        }
        
        //actual line
        String line = "";
        
        //auxiliary variable
        String line1 = "";
        String line2 = "";
        String noradString = "";
        String lineDate = "";
        int norad = 0;
        
        //number of lines
        int no = 0;
        int no2 = 0;
        int noOfObjects = 0;
        
        //do cycle
        boolean doCycle = true;
        boolean doCycle2 = true;
        
        //array to fill
        //SatelliteArray sa = new SatelliteArray();
        
        //FillConstructor methods
        FillConstructor fc = new FillConstructor();

        //no of on Earth orbit objects
        int onOrbitObjects = 0;
        //try read file
	try{	
            
            //read file from specified position
            file = new BufferedReader(new FileReader(filePos));
            
            //Beginning of report
            //read line
            while(doCycle){
                line = file.readLine();
                line = line + "";
                no++;
                if(no == 15){
                    lineDate = line;
                    lineDate = lineDate.substring(26,lineDate.length());
                    sa.reportDate = lineDate;
                    //doCycle = false;
                }
                //if there are 5 this lines, report starts
                if(line.equals("-------------------------------------------------------------------------------")){
                    no2++;
                    if(no2 == 5) doCycle = false;
                }
                
                //if the number no is larger than for example 10 000, it is possible, that choosed file
                //is not Satellite report, and the while cycles wont stop, because the condition wont
                //be executed
                if(no > 10000) {
                    doCycle = false;
                    doCycle2 = false;
                }
            }
            
            //satellite datas
            //read line
            while(doCycle2){
                line = file.readLine();
                //System.out.println("line2 " + line);
                line = line + "";
                noOfObjects ++;
                
                //condition
                    //unpaired line (1st line of TLE datas)
                    if(noOfObjects % 2 != 0) line1 = line;
                    //paired line (2nd line of report datas)
                    else {
                        line2 = line;
                        //now we can fill the constructor Satellite and SatelliteArray
                        //what object it is (NORAD cat No)
                        if(line.equals("null")){}
                        
                        else {
                            noradString = line1.substring(16,23);
                            norad = (int)Double.parseDouble(noradString);
                            sa.array[norad] = fc.fillConstructorReport(line1,line2,sa.array[norad]);
                            //fill report lines
                            sa.array[norad].line1Rep = line1;
                            sa.array[norad].line2Rep = line2;

                            String helpString = (line2 + "                              ").substring(0, 90);
                            if((line2.indexOf("Decayed") != -1)&&(line1.indexOf("ELEMENTS") == -1)&&(line1.indexOf("ORBIT") == -1)
                                    &&(line1.indexOf("LUNA") == -1)&&(line1.indexOf("ISS") == -1)) //System.out.println(helpString + "                        " + line1);
                            //searching for lagrest catalog no.
                            if(norad > sa.largestNo) sa.largestNo = norad;
                            //no of objects on Earth orbit
                            if(sa.array[norad].isOnOrbit){
                                if((sa.array[norad].line1Rep.indexOf("NO ELEMENTS AVAILABLE")== -1))
                                        //||( sa.array[norad].apogee != 0))
                                {
                                     onOrbitObjects++;
                                     //System.out.println(sa.array[norad].line1Rep);
                                }
                            }
                        }
                    }
                //System.out.println(line);
                //end of file
                if(line.equals("null")) doCycle2 = false;
            }
            
            sa.allAmount = noOfObjects / 2 -1;
            sa.objectsOnOrbit = onOrbitObjects;
            
            file.close();
        }
        
         //if file can not be open, io exception is used
        catch(IOException e) {
            System.out.println("Error: problems read satellite report!");
        }
        
        System.out.println("Done!");
        
        //return filled constructor SatelliteArray
        return sa;
    }
     
     /**
     *Method: readFileSatRep()
     *
     *Method is to fill constructor SatelliteArray from satellite report file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\file.txt")
     *      Satellite array    
     *
     *OUT:  constructor filled SatelliteArray
     *
     */
    
     public SatelliteArray readFileSatRep_New(String filePos, SatelliteArray sa) {
        
        System.out.println("Reading report from " + filePos);

         //if user is in Windows, read file like X:\\DIRECTORY\\file.txt
        if(osName.startsWith("Windows")){} //do nothing
        //UNIX, or Linux OS
        else{   //replace all backslashes
            String helpString = "\\";
            while(filePos.indexOf("\\") != -1){
                int slashPos = filePos.indexOf(helpString);
                filePos = filePos.substring(0,slashPos) + "//" + filePos.substring(slashPos+1, filePos.length());                
            }
        }
        
        //actual line
        String line = "";
        
        //auxiliary variable
        String line1 = "";
        String line2 = "";
        String noradString = "";
        String lineDate = "";
        int norad = 0;
        
        //number of lines
        int no = 0;
        int no2 = 0;
        int noOfObjects = -1;
        
        //do cycle
        boolean doCycle = true;
        boolean doCycle2 = true;
        
        //array to fill
        //SatelliteArray sa = new SatelliteArray();
        
        //FillConstructor methods
        FillConstructor fc = new FillConstructor();

        //no of on Earth orbit objects
        int onOrbitObjects = 0;
        //try read file
	try{	
            
            //read file from specified position
            file = new BufferedReader(new FileReader(filePos));

            //satellite datas
            //read line
            while(doCycle2){
                line = file.readLine();
                //System.out.println("line2 " + line);
                line = line + "";
                if (noOfObjects < 0){
                    noOfObjects ++;
                }
                
                else{   
                    noOfObjects ++;                 
                
                    //now we can fill the constructor Satellite and SatelliteArray
                    //what object it is (NORAD cat No)
                    if(line.equals("null")){}
                        
                    else {
                        // Split the CSV string into sub-strings, get NORAD number
                        String array_csv[] = line.split(",");
                        array_csv[2] = array_csv[2].replaceAll("\"", "");                        
                        noradString = array_csv[2];
                        norad = (int)Double.parseDouble(noradString);
                        sa.array[norad] = fc.fillConstructorReportNew(line,sa.array[norad]);
                        //fill report lines
                        sa.array[norad].line1Rep = line;
                        sa.array[norad].line2Rep = "";
                    }
                }

                if(line.equals("null")) doCycle2 = false;
            }
            
            sa.allAmount = noOfObjects;
            sa.objectsOnOrbit = noOfObjects;
            
            file.close();
        }
        
         //if file can not be open, io exception is used
        catch(IOException e) {
            System.out.println("Error: problems read satellite report!");
        }
        
        System.out.println("Done!");
        
        //return filled constructor SatelliteArray
        return sa;
    }
     
     /**
     *Method: readStarCatalogueLandolt()
     *
     *Method is to fill constructor SatelliteArray from satellite report file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\data\\star_catalogue.txt")
     *      Satellite array    
     *
     *OUT:  constructor filled SatelliteArray
     *
     */
    
     public StarArray readStarCatalogueLandolt(String filePos, StarArray sta) {
        
         //Star name
         String star_name = "";
         //Right ascension hour
         double star_ra_hour = 0.0;
         //Right ascension min
         double star_ra_min = 0.0;
         //Right ascension sec
         double star_ra_sec = 0.0;
        //Declination sign
         double star_dec_sign = 0.0; 
        //Declination [deg]
         double star_dec_deg = 0.0;
         //Declination [deg]
         double star_dec_min = 0.0;
         //Declination [deg]
         double star_dec_sec = 0.0;
         //Magnitude
         double star_mag = 0.0;
         //String help
         String string_help = "";
         
        System.out.println("Reading star catalogue from " + filePos);

         //if user is in Windows, read file like X:\\DIRECTORY\\file.txt
        if(osName.startsWith("Windows")){} //do nothing
        //UNIX, or Linux OS
        else{   //replace all backslashes
            String helpString = "\\";
            while(filePos.indexOf("\\") != -1){
                int slashPos = filePos.indexOf(helpString);
                filePos = filePos.substring(0,slashPos) + "//" + filePos.substring(slashPos+1, filePos.length());
            }
        }
        
        //actual line
        String line = "";
               
        //number of lines
        int no = 0;
        int noOfObjects = 0;
        
        //do cycle
        boolean doCycle = true;
        boolean doCycle2 = true;
        
        //FillConstructor methods
        FillConstructor fc = new FillConstructor();

        //no of on Earth orbit objects
        int onOrbitObjects = 0;
        //try read file
	try{            
            //read file from specify position
            file = new BufferedReader(new FileReader(filePos));
            
            //Beginning of report
            //read line
            while(doCycle){
                line = file.readLine();
                //System.out.println(line);
                try{
                    if(line.equals("null")) doCycle = false;
                    else{                    
                        line = line + "";
                        star_name = line.substring(0,13);
                        star_ra_hour = Double.parseDouble(line.substring(13,15));
                        star_ra_min = Double.parseDouble(line.substring(16,18));
                        star_ra_sec = Double.parseDouble(line.substring(19,21));
                        string_help = line.substring(24,25);
                        if(string_help.contains("-")){
                            star_dec_sign = -1;
                        }
                        else{
                            star_dec_sign = 1;
                        }
                        star_dec_deg = Double.parseDouble(line.substring(25,27));
                        star_dec_min = Double.parseDouble(line.substring(28,30));
                        star_dec_sec = Double.parseDouble(line.substring(31,33));
                        star_mag = Double.parseDouble(line.substring(34,40));
                                                
                        //System.out.println(star_name);
                        //System.out.println(star_ra_hour);                
                        //System.out.println(star_ra_min);                
                        //System.out.println(star_ra_sec);
                        //System.out.println(star_dec_sign);
                        //System.out.println(star_dec_deg);
                        //System.out.println(star_dec_min);
                        //System.out.println(star_dec_sec);
                        //System.out.println(star_mag);
                        //Fill the star constructor
                        star_mag = (float)((int)(star_mag*10))/10;
                        string_help = star_name.substring(star_name.length()-1, star_name.length());
                        while(string_help.equalsIgnoreCase(" ")){
                            string_help = sta.array[no].name;
                            star_name = star_name.substring(0, star_name.length()-1);
                            string_help = star_name.substring(star_name.length()-1, star_name.length());
                        }
                        
                        sta.array[no].name = star_name + " (" + (float)(star_mag) + ")";
                        sta.array[no].ra = Math.toRadians((star_ra_hour + star_ra_min/60 + star_ra_sec/3600)*15);
                        sta.array[no].dec = Math.toRadians((star_dec_deg + star_dec_min/60 + star_dec_sec/3600)*star_dec_sign);
                        
                        //System.out.println(sta.array[no].name);
                        //System.out.println(Math.toDegrees(sta.array[no].ra));
                        //System.out.println(Math.toDegrees(sta.array[no].dec));
                        no++;
                    }
                }
                catch(Exception e){
                    doCycle = false;
                    //System.out.println("Reading file ERROR message: " + e.getMessage());                    
                }
            }
            sta.largestNo = no;
            file.close();
        }
        
         //if file can not be open, io exception is used
        catch(IOException e) {
            System.out.println("Error: problems reading star catalogue!");
        }
        
        System.out.println("Done!");
        
        //return filled constructor SatelliteArray
        return sta;
    }
     
     /**
     *Method: readFileTleKelso()
     *
     *Method is to fill constructor SatelliteArray from TLE file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\")
     *      Saellite Array
     *
     *OUT:  constructor filled SatelliteArray
     *
     */
     
      public SatelliteArray readFileTleKelso(String filePos, SatelliteArray sa) {
          //array of TLE file name
          String fileName[] = {"1999-025", "amateur", "cubesat", "dmc", "education",
                                "engineering", "galileo", "geo", "geodetic",
                                "globalstar", "glo-ops", "goes", "gorizont", 
                                "gps-ops", "intelsat", "iridium", "military", 
                                "molniya", "musson", "nnss", "noaa", "orbcomm", 
                                "other", "other-comm", "radar", "raduga", "resource",
                                "sarsat", "sbas", "science", "stations", "tdrss", 
                                "tle-new", "visual", "weather", "x-comm", "iridium-33-debris",
                                "cosmos-2251-debris"};
          String filePos2 = "";
          int amountHelp = 0;
          sa.tleKelsoAmount = 0;
          //filling array
          for(int i = 0; i < fileName.length; i++){
              filePos2 = filePos + fileName[i] + ".txt";
              sa = readFileTle(filePos2, sa, true);
              amountHelp = amountHelp + sa.tleKelsoAmount;
          }
          
          sa.tleKelsoAmount = amountHelp;
          return sa;
      }

      /**
     *Method: readFileTleKelsoSuppl()
     *
     *Method is to fill constructor SatelliteArray from TLE file
     *
     *IN:   String position of file ("X:\\DIRECTORY\\")
     *      Saellite Array
     *
     *OUT:  constructor filled SatelliteArray
     *
     */

      public SatelliteArray readFileTleKelsoSuppl(String filePos, SatelliteArray sa) {
          //array of TLE file name
          String fileName[] = {"gps", "glonass", "meteosat", "intelsat", "orbcomm"};
          String filePos2 = "";
          int amountHelp = 0;
          sa.tleKelsoAmount = 0;
          //filling array
          for(int i = 0; i < fileName.length; i++){
              //System.out.println("Before " + amountHelp + " - " + fileName[i]);
              filePos2 = filePos + fileName[i] + ".txt";
              sa = readFileTle(filePos2, sa, true);
              amountHelp = amountHelp + sa.tleKelsoAmount;
              //System.out.println("After " + amountHelp + " - " + fileName[i]);
          }

          sa.tleKelsoAmount = amountHelp;
          return sa;
      }

      /**
       * readFileSettings()
       *
       * Method to read settings files.
       *
       * IN:
       *    String filePos
       *
       * OUT:
       *    String double[11] - double array with every line of given file
       *                        mjd/exp/lon/lat/alt/RA/dec/size/width/height/type of coordinates - 0 - Hor, 1 - Eq
       */
      public double[] readFileSettings(String filePos){
            boolean doCycle = true;
            String array[] = new String[50];
            double arrayDouble[] = new double[11];
            try{
                //read file from specify position
                file = new BufferedReader(new FileReader(filePos));
                int i = 0;
                //read line
                while(doCycle){
                    String line = file.readLine();
                    if(line.equals("null")) doCycle = false;
                    else {
                        array[i] = line;
                        i++;
                    }
                }
                String array2[] = new String[i+1];
                for(int j = 0; j <= i; j++){
                    array2[j] = array[j];
                }
                //mjd time
                arrayDouble[0] = Double.parseDouble(array2[3].substring(12,array2[3].length()));
                //exposure time [s]
                arrayDouble[1] = Double.parseDouble(array2[4].substring(10,array2[4].length()));
                //obs position
                arrayDouble[2] = Double.parseDouble(array2[5].substring(10,18));
                arrayDouble[3] = Double.parseDouble(array2[5].substring(18,28));
                arrayDouble[4] = Double.parseDouble(array2[5].substring(28,array2[5].length()));
                //canvas position
                String typeOfCoor = array2[6].substring(17,21);
                if(typeOfCoor.equals("Hor")) arrayDouble[10] = 0;
                else arrayDouble[10] = 1;
                arrayDouble[5] = Double.parseDouble(array2[6].substring(22,32));
                arrayDouble[6] = Double.parseDouble(array2[6].substring(32,array2[6].length()));
                //canvas size
                arrayDouble[7] = Double.parseDouble(array2[7].substring(18,25));
                arrayDouble[8] = Double.parseDouble(array2[7].substring(25,31));
                arrayDouble[9] = Double.parseDouble(array2[7].substring(31,array2.length));

            }
            catch(Exception e){
                System.out.println(e);
            }

            
            return arrayDouble;
      }
}
