/*
 * FillConstructor.java
 *
 * Created on February 28, 2008, 3:11 PM
 *
 *  SK: Trieda obsahuje metody na vyplnovanie konstruktora Satellite
 *      a to aj z TLE dat aj zo satelitnej spravy.
 *  EN: Class with methods to fill constructor Saterllite. Sources are TLE database
 *      and satellite report.
 *
 * NOTES: Commentaries to the methods not necessary.
 */

package InOut;

import SatEph.*;
import Compute.Silha.*;

/**
 *
 * @author Jiri Silha
 */
public class FillConstructor {
    
    //TLE
    
    public double epoch = 0;                //Element set epoch
    
    public double der1mm = 0;               //1st derivato of mean motion
    
    public double der2mm = 0;               //2nd derivato of mean motion
    
    public double bDrag = 0;                //B*Drag term
    
    public int elementNo = 0;               //Element number
    
    public int revNoEp = 0;                 //Revolution Number at Epoch    
    
    public double incl = 0;                 //Orbit inclination (degrees)
    
    public double node = 0;                 //Right Ascension of Ascending Node (degrees)
    
    public double argument = 0;             //Argument of Perigee (degrees)
    
    public double meanAnom = 0;             //Mean Anomaly (degrees)
    
    public double ecc = 0;                  //Eccentricity
    
    public double meanMot = 0;              //Mean Motion (revolutions/day)
    
    
    //REPORT
   
    public String name = "NoObject";            //Name of object
    
    public String intDes = "NoObject";          //International designator
   
    public double rcs = 0;                      //Radar Cross Section
    
    public String source = "";                  //Organisation
    
    public String launched = "";                //Date of launch
        
    public String decayed = "";                 //Date of decayed
    
    public boolean isOnOrbit = true;            //Is the object on Earth orbit?
    
    public boolean isOrbitalDebris = false;     //Is the object orbital debris? (all debris object except R/B) 
                                                //Debris objects are considered to be any object with the string 'DEB', 
                                                //or 'COOLANT' or 'SHROUD' or 'WESTFORD NEEDLES' in the SatCat object common name.  
    
    public boolean isRB = false;                //Is the object rocket body (R/B)?
                                                //Rocket Bodies are considered to be any object which contains the strings 'R/B' or
                                                //'AKM' (Apogee Kick Motor) or 'PKM' (Perigee Kick Motor) but not 'DEB' in the SatCat object common name.    
    
    String string = "";
    
    Time actualTime = new Time();
    Kepler kepler = new Kepler();
    Constants constants = new Constants();
    
    /** Creates a new instance of FillConstructor */
    public FillConstructor() {
    }
    
     /**
     *Method: fillConstructorTle()
     *
     *Method is to fill constructor Satellite from TLE catalog
     *
     *IN:   line1
     *      line2
     *
     *OUT:  constructor Satellite
     *
     */
    
    public Satellite fillConstructorTle(String line1, String line2, Satellite sat){
        try{
            //Satellite sat = new Satellite();

            //filling constructor
            //line1
            //sat.epoch = getEpoch(line1);
            sat.epochMjd = getEpochMjd(line1);

            //System.out.println(line1);
            //System.out.println(sat.epoch);
            //System.out.println(sat.epochMjd);
            //sat.der1mm = (float)getDer1mm(line1);
            //sat.der2mm = (float)getDer2mm(line1);
            sat.bDrag = (float)getBDrag(line1);
            sat.elementNo = getElementNo(line1);
            if(sat.intDes.equalsIgnoreCase("No id")) sat.intDes = getIntDesTLE(line1);

            //line2
            sat.revNoEp = getRevNoEp(line2);
            sat.incl = (float) getIncl(line2);
            sat.node = (float) getNode(line2);
            sat.argument = (float) getArgument(line2);
            sat.meanAnom = (float) getMeanAnom(line2);
            sat.ecc = (float) getEcc(line2);
            sat.meanMot = (float) getMeanMot(line2);
            
            //0 line info about object name
            sat.isOrbitalDebris = getIsOrbitalDebris(sat.name);
            if(sat.isOrbitalDebris == false)sat.isRB = getIsRB(sat.name);

            //fill semi major axis
            sat.smAxis = (float) Kepler.getSMAxis2(constants.GM_Earth, constants.pi2*sat.meanMot/86400);
            //period
            sat.period = 1440/sat.meanMot;
        }

        catch(Exception e){
            System.out.println("Warning! Problem getting elements from TLE!" + e.getMessage());
        }
        
        return sat;
    }
    
    /**
    *   Small methods to get the each value from TLE lines
    *
    */
    
    //line1
    public double getEpoch(String string){
        string = string.substring(18,32);
        string = string.replaceAll(" ","0");
        double value = Double.parseDouble(string);
        return value;
    }
    
    //line1
    public double getEpochMjd(String string){
        double valueMjd = 0;
        try{
            //get epoch (YEAR DOY -> YYDDD.DDDDD...)
            string = string.substring(18,32);
            string = string.replaceAll(" ","0");
            //get YY
            String stringYear = string.substring(0,2);
            int valueYear = (int)(Double.parseDouble(stringYear));
            if(valueYear >= 57 ) valueYear = valueYear + 1900;
            else valueYear = valueYear + 2000; 
            //int valueYear = (int)(Double.parseDouble(stringYear))+2000;
            double yearMjd = actualTime.getMjd(valueYear,1,1,0,0,0.0);
            //get DDD.DDDD...
            string = string.substring(2,14);
            valueMjd = Double.parseDouble(string) + yearMjd-1;
        }
        catch(Exception e){
            System.out.println("Warning! " + e.getCause());
            System.out.println("Wrong string format! " + string);
            valueMjd = 51544.0;
        }
        return valueMjd;
    }
    
     public double getDer1mm(String string){
        string = string.substring(33,43);
        double value = 0;
        try{
            value = Double.parseDouble(string);
        }
        catch(Exception e){
            System.out.println(e);
            value = 0;
        }
        return value;
    }
     
     public double getDer2mm(String string){
        String string1, string2;
        string1 = string.substring(44,45);
        string2 = string.substring(45,52);
        //string = string.substring(44,52);
        string2 = string2.replaceAll("-","e-");
        string2 = string2.replace('+','e');
        string = string1+string2;
        double value = Double.parseDouble(string);
        return value;
    }
    
     public double getBDrag(String string){
        String string1, string2;
        double value;
        
        try{ 
            string1 = string.substring(53,54);
            string2 = string.substring(54,61);
            //System.out.println("1st BDrag:"+ string2);
            string2 = string2.replaceAll("-","e-");
            string2 = string2.replace('+','e');
            string2 = "0."+string2;
            string = string1+string2;
            //System.out.println("2nd BDrag:"+ string);
            value = Double.parseDouble(string);
         }
        catch(Exception e){
            value = 0.0;
            System.out.println("Warning! Problems getting B*Drag values!" + e.getMessage());
        }
        return value;
     }
     
    public int getElementNo(String string){
        int value = 0;
        try{
            string = string.substring(64,68);
            value = (int)Double.parseDouble(string);
        }
        catch(Exception e){
            value = 0;
            System.out.println("Warning! Problems getting element number value!" + e.getMessage());
        }
        return value;
    }
     
     //line2
     public int getRevNoEp(String string){
        string = string.substring(63,68);
        int value = (int)Double.parseDouble(string);
        return value;
    }
     
     public double getIncl(String string){
        string = string.substring(8,16);
        double value = Double.parseDouble(string);
        return value;
    }
     
     public double getNode(String string){
        string = string.substring(17,25);
        //System.out.println(string);
        double value = Double.parseDouble(string);
        return value;
    }
     
     public double getArgument(String string){
        string = string.substring(34,42);
        double value = Double.parseDouble(string);
        return value;
    }
     
      public double getMeanAnom(String string){
        string = string.substring(43,51);
        double value = Double.parseDouble(string);
        return value;
    }
      
       public double getEcc(String string){
        string = string.substring(26,33);
        double value = Double.parseDouble("0."+string);
        return value;
    }
       
        public double getMeanMot(String string){
        string = string.substring(52,63);
        double value = Double.parseDouble(string);
        return value;
    }
        
        public String getIntDesTLE(String line){
        String value = line.substring(9, 18);
        String yearString = value.substring(0, 2);
        value = value.substring(2);
        int year = 0;
        try{
            year = Integer.parseInt(yearString);
            if(year < 57) year = year + 2000;
            else year = year + 1900;
        }
        catch(Exception e){
            System.out.println(e);
        }
        
        value = year + "-" + value;
        
        return value;
    }
        
     /**
     *Method: fillConstructorReport()
     *
     *Method is to fill constructor Satellite from satellite report
     *
     *IN:   line1
     *      line2
     *
     *OUT:  constructor Satellite
     *
     */
    
    public Satellite fillConstructorReport(String line1, String line2, Satellite sat){
        //fill constructor
        
        //line1
        sat.intDes = getIntDes(line1);
        sat.source = getSource(line1);
        
        if((line1.indexOf("NO ELEMENTS AVAILABLE") == -1)&&
            (line1.indexOf("HELIOCENTRIC ORBIT (SUN)") == -1)&&
            (line1.indexOf("MARS LANDING") == -1)&&
            (line1.indexOf("MARS ORBIT") == -1)&&
            (line1.indexOf("MARS IMPACT") == -1)&&
            (line1.indexOf("DOCKED TO ISS") == -1)&&
            (line1.indexOf("DOCKED TO ETS 7") == -1)&&
            (line1.indexOf("NO INITIAL ELEMENTS") == -1)&&
            (line1.indexOf("LUNAR IMPACT") == -1)&&
            (line1.indexOf("LUNAR LANDING") == -1)&&
            (line1.indexOf("CIRCUMLUNAR") == -1)&&
            (line1.indexOf("BARYCENTRIC ORBIT (EARTH-MOON)") == -1)&&
            (line1.indexOf("JOVIAN ORBIT (JUPITER)") == -1)&&
            (line1.indexOf("VENUS ORBIT") == -1)&&
            (line1.indexOf("VENUS IMPACT") == -1)&&
            (line1.indexOf("VENUS LANDING") == -1)&&
            (line1.indexOf("ESCAPED SOLAR SYSTEM") == -1)&&
            (line1.indexOf("L2 LAGRANGIAN ORBIT") == -1)&&
            (line1.indexOf("SATURNIAN ORBIT (SATURN)") == -1)&&
            (line1.indexOf("COLLIDED WITH SATELLITE") == -1)&&
            (line1.indexOf("SELENOCENTRIC ORBIT (MOON)") == -1)){
            //there can be some chages in sat rep, which can cause the fall down of program
            try{
                //sat.apogee = getApogee(line1);
                //sat.perigee = getPerigee(line1);
            }
            catch(Exception ex){
                System.out.println(ex);
           }
        }
        
        else {
            //sat.apogee = 0;
            //sat.perigee = 0;
         }
        
        if(line1.indexOf("N/A") == -1){
            try{
                sat.rcs =  (float)getRcs(line1);
            }
            catch(Exception ex){
                 System.out.println(ex);
           }
        }
        else sat.rcs = 0;
        
        //line2
        sat.name = getName(line2);
        sat.launched = getLaunched(line2);
        sat.isOnOrbit = getIsOnOrbit(line2);
        
        if(sat.isOnOrbit == false){
            sat.decayed = getDecayed(line2);
        }
        else sat.decayed = "On orbit";
        
        sat.isOrbitalDebris = getIsOrbitalDebris(sat.name);
        if(sat.isOrbitalDebris == false)sat.isRB = getIsRB(sat.name);
        
        return sat;
    }
    
    /**
    *   Small methods to get the each value from sat. report lines
    *
    */
    
    //line1
    
    public String getIntDes(String line){
        String value = "NoObject";
        value = line.substring(0,15);
        return value;
    }
    
     public String getSource(String line){
        String value = "NoObject";
        value = line.substring(24,31);
        return value;
    }
    
    public double getApogee(String line){
        double value = 0;
        string = line.substring(48,56);
        //System.out.println("Test string " + string);
        value = Double.parseDouble(string);
        return value;
    }
    
     public double getPerigee(String line){
        double value = 0;
        string = line.substring(57,67);
        value = Double.parseDouble(string);
        return value;
    }
     
     //not all TLE nies have got the same structure, those are the exceptions
     public double getRcs(String line){
        double value = 0;
        if((line.indexOf("NO ELEMENTS AVAILABLE") != -1)||
            (line.indexOf("HELIOCENTRIC ORBIT (SUN)") != -1)||
            (line.indexOf("MARS LANDING") != -1)||
            (line.indexOf("MARS ORBIT") != -1)||
            (line.indexOf("MARS IMPACT") != -1)||
            (line.indexOf("DOCKED TO ISS") != -1)||
            (line.indexOf("DOCKED TO ETS 7") != -1)||
            (line.indexOf("NO INITIAL ELEMENTS") != -1)||
            (line.indexOf("LUNAR IMPACT") != -1)||
            (line.indexOf("LUNAR LANDING") != -1)||
            (line.indexOf("CIRCUMLUNAR") != -1)||
            (line.indexOf("BARYCENTRIC ORBIT (EARTH-MOON)") != -1)||
            (line.indexOf("JOVIAN ORBIT (JUPITER)") != -1)||
            (line.indexOf("VENUS ORBIT") != -1)||
            (line.indexOf("VENUS IMPACT") != -1)||
            (line.indexOf("VENUS LANDING") != -1)||
            (line.indexOf("ESCAPED SOLAR SYSTEM") != -1)||
            (line.indexOf("L2 LAGRANGIAN ORBIT") != -1)||
            (line.indexOf("SATURNIAN ORBIT (SATURN)") != -1)||
            (line.indexOf("COLLIDED WITH SATELLITE") != -1)||
            (line.indexOf("SELENOCENTRIC ORBIT (MOON)") != -1)){
            string = line.substring(65,72);
        }
        
        else  string = line.substring(68,79);
        
        value = Double.parseDouble(string);
            
        return value;
    }
     
     //line2
     
      public String getName(String line){
        String value = "";
        value = line.substring(0,38);
        return value;
    }
      
      public String getLaunched(String line){
        String value = "";
        value = line.substring(48,60);
        return value;
      }
      
      //is on orbit? that means all objectsm, which are note decayed, but they
      //can be on heliocenstric or lunar orbit!!!
      public boolean getIsOnOrbit(String line){
        boolean value = true;
        if(line.indexOf("Decayed") != -1){
            value = false;
        }
        return value;
      }
      
      public String getDecayed(String line){
        String value = "";
        value = line.substring(70,82);
        return value;
      }
    
      public boolean getIsOrbitalDebris(String name){
        boolean value = false;
        if((name.indexOf("DEB") != -1)||(name.indexOf("COOLANT") != -1)||
               (name.indexOf("SHROUD") != -1)||(name.indexOf("WESTFORD NEEDLES") != -1)){
               value = true;
        };
        return value;
      }
      
       public boolean getIsRB(String name){
        boolean value = false;
        if((name.indexOf("R/B") != -1)||(name.indexOf("AKM") != -1)||
               (name.indexOf("PKM") != -1)){
               value = true;
        };
        return value;
      }  
       
       /**
     *Method: fillConstructorReport()
     *
     *Method is to fill constructor Satellite from satellite report, new format in CSV
     *
     *IN:   line1
     *      line2
     *
     *OUT:  constructor Satellite
     *
     */
    
    public Satellite fillConstructorReportNew(String line, Satellite sat){
        //fill constructor
        
        // Split the CSV string into sub-strings
        String array_csv[] = line.split(",");
                            
        for(int i=0;i<array_csv.length;i++){                                
            array_csv[i] = array_csv[i].replaceAll("\"", "");
            //System.out.print(array_csv[i] + "; ");                            
        }
        
        //line1
        sat.intDes = array_csv[0];
        sat.source = array_csv[3];
        /*
        if((line1.indexOf("NO ELEMENTS AVAILABLE") == -1)&&
            (line1.indexOf("HELIOCENTRIC ORBIT (SUN)") == -1)&&
            (line1.indexOf("MARS LANDING") == -1)&&
            (line1.indexOf("MARS ORBIT") == -1)&&
            (line1.indexOf("MARS IMPACT") == -1)&&
            (line1.indexOf("DOCKED TO ISS") == -1)&&
            (line1.indexOf("DOCKED TO ETS 7") == -1)&&
            (line1.indexOf("NO INITIAL ELEMENTS") == -1)&&
            (line1.indexOf("LUNAR IMPACT") == -1)&&
            (line1.indexOf("LUNAR LANDING") == -1)&&
            (line1.indexOf("CIRCUMLUNAR") == -1)&&
            (line1.indexOf("BARYCENTRIC ORBIT (EARTH-MOON)") == -1)&&
            (line1.indexOf("JOVIAN ORBIT (JUPITER)") == -1)&&
            (line1.indexOf("VENUS ORBIT") == -1)&&
            (line1.indexOf("VENUS IMPACT") == -1)&&
            (line1.indexOf("VENUS LANDING") == -1)&&
            (line1.indexOf("ESCAPED SOLAR SYSTEM") == -1)&&
            (line1.indexOf("L2 LAGRANGIAN ORBIT") == -1)&&
            (line1.indexOf("SATURNIAN ORBIT (SATURN)") == -1)&&
            (line1.indexOf("COLLIDED WITH SATELLITE") == -1)&&
            (line1.indexOf("SELENOCENTRIC ORBIT (MOON)") == -1)){
            //there can be some chages in sat rep, which can cause the fall down of program
            try{
                //sat.apogee = getApogee(line1);
                //sat.perigee = getPerigee(line1);
            }
            catch(Exception ex){
                System.out.println(ex);
           }
        }
        
        else {
            //sat.apogee = 0;
            //sat.perigee = 0;
         }
        
        if(line1.indexOf("N/A") == -1){
            try{
                sat.rcs =  (float)getRcs(line1);
            }
            catch(Exception ex){
                 System.out.println(ex);
           }
        }
        else sat.rcs = 0;
        */        
        sat.name = array_csv[1];
        sat.launched = array_csv[10];
        sat.isOnOrbit = true;
        
        //RCS value, http://www.satobs.org/seesat/Aug-2014/0137.html
        // Small < 0.1m2, MEDIUM 0.1-1.0, LARGE >1.0        
        if(array_csv[8].equalsIgnoreCase("SMALL")){
            sat.rcs = 0.1f;
        }
        else if(array_csv[8].equalsIgnoreCase("MEDIUM")){
            sat.rcs = 0.5f;
        }
        else if(array_csv[8].equalsIgnoreCase("LARGE")){
            sat.rcs = 5f;
        }
        else{
            sat.rcs = 0;    //set to 0, now in new format SMALL, MEDIUM, LARGE only options
        }
        
        sat.isOrbitalDebris = getIsOrbitalDebris(sat.name);
        if(sat.isOrbitalDebris == false)sat.isRB = getIsRB(sat.name);
        
        return sat;
    }
}
