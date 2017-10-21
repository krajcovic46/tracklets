/*
 * SaveString.java
 *
 * Created on Nedela, 2008, marec 30, 3:12
 *
 * SK: trieda sluzi na ukladanie textu do konkretnych suborov
 * En: class to save String to files
 */

package InOut;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.OutputStreamWriter;
import java.io.FileOutputStream;
import java.awt.image.*;

import Compute.Silha.*;
import SatEph.*;

/**
 *
 * @author Jiri Silha
 */
public class SaveString {
    
    /** Creates a new instance of SaveString */
    public SaveString() {
    }
    
    FileWriter fw;
    
    /**
     * Method: save string to file
     *
     * IN: String file name
     * 
     * OUT: 
     *
     */
    
    public void saveFileString(String name, String text){
    	try{    
            //System.out.println("String to save "  + name);
		PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(name + ".dat"), "UTF-8"));
		pw.write(text);
		pw.flush();
		pw.close();
	}
	catch (Exception e) {
		System.out.println("Exception: " + e);
	}
    }
    
    /**
     * Method: save image and text files
     *
     * IN:  String file names
     *      String text
     *      RenderedImage image
     * 
     * OUT: 
     *
     */
    public void saveFiles(String name, String text, RenderedImage image){
        //saveFileString("save//"+name, text);
        (new SaveImage()).savePictureJpg("save//"+name, image);
    }
    
    /**
     * Method: get name of files (for image, info)
     *
     * IN:  int norad#
     *      Time time
     * 
     * OUT: String name of files
     *
     */
    public String getFilesName(int norad, Time time){
        String string = "#" + norad + "_" + getStringDateTime(time);
        return string;
    }
    
     /**
     * getInfoString()
     * Method: get test to be saved to file (info), this method get the string to save
      * and also it saves this string to file called fileName.dat
     *
     * IN:  String fileName - name of file, where the datas will be written
     *      Satellite sat
     *      int exp time
     *      Time time
     *      Geodetic observer
     *      SatEphData sad - array to hold the other positions of object
     * 
     * OUT: String info text to save
     *  
     */
    public String getInfoString(String fileName, Satellite sat, int exp, Time time, Geodetic observer,
            SatEphDataArray seda){
        String mainString = "";
        PrintWriter pw;
        try{
            //create file .dat
            pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream("save//" + fileName + ".dat"), "UTF-8"));
            //real exp time
            //exp = exp;

            /*
            float setEpoch = (float)time.getMjd(time);
            mainString = "----- NORAD# " + sat.catalogNo + " -----" + "\n\n";
            mainString = mainString + sat.line1 + "\n" + sat.line2 + "\n\n";
            mainString = mainString + "Start time:  " + getStringDateTime2(time) + "\n";
            mainString = mainString + "Exp time:    " + exp + " s" + "\n";
            pw.write(mainString + "\n");
            String string[] = new String[20];
            string[0] = "Az                  ";
            string[1] = "h                   ";
            string[2] = "R.A.                ";
            string[3] = "dec                  ";
            string[4] = "Ma(t)               ";
            string[5] = "v(t)               ";      //not yet
            string[6] = "Phase               ";
            //string[7] = "T since            ";
            string[8] = "Obs. r.              ";
            string[9] = "Center r.             ";
            //string[10] = "V1                        ";
            //string[11] = "H1                        ";
            //string[12] = "V2                        ";
            //string[13] = "H2                        ";
            string[14] = "Shadow                    ";

            string[19] = string[0].substring(0,12) + " " + string[1].substring(0,12) + " " + string[2].substring(0,12) + " " +
                        string[3].substring(0,12) + " " + string[4].substring(0,12) + " " + //string[5].substring(0,12) + " " +
                        string[6].substring(0,12) + " " + string[7].substring(0,10) + " " + string[8].substring(0,12) + " " +
                        string[9].substring(0,12) + " " + string[10].substring(0,6) + " " + string[11].substring(0,6) + " " +
                        string[12].substring(0,6) + " " + string[13].substring(0,6) + " " + string[14].substring(0,7);

            pw.write(string[19] + "\n");
            string[19] =("[deg]                     ").substring(0,12) + " " + ("[deg]                     ").substring(0,12) + " " + 
                        ("[deg]                     ").substring(0,12) + " " + ("[deg]                     ").substring(0,12) + " " + 
                        ("[deg]                     ").substring(0,12) + " " + //("[deg]                     ").substring(0,12) + " " +
                        ("[deg]                     ").substring(0,12) + " " + ("[day]                     ").substring(0,10) + " " + 
                        ("[km]                     ").substring(0,12) + " " +  ("[km]                     ").substring(0,12) + " " + 
                        ("                     ").substring(0,6) + " " + ("                   ").substring(0,6) + " " +
                        ("                     ").substring(0,6) + " " + ("                     ").substring(0,6) + " " +
                        ("                     ").substring(0,7);
            pw.write(string[19] + "\n\n");
            //mainString = mainString + "\n" + string[19] + "\n";
            //pw.write(mainString);

            for(int i = 0; i < exp; i++){
                string[0] = modifyNumber(Math.toDegrees(sat.array[i].azimuth),4) + "            ";
                string[1] = modifyNumber(Math.toDegrees(sat.array[i].altitude),4) + "            ";
                string[2] = modifyNumber(Math.toDegrees(sat.array[i].ra),4) + "            ";
                string[3] = modifyNumber(Math.toDegrees(sat.array[i].dec),4) + "            ";
                string[4] = modifyNumber(Math.toDegrees(sat.array[i].actMeanAnom),4) + "            ";
                string[5] = "                                    ";
                string[6] = modifyNumber(Math.toDegrees(sat.array[i].phaseAngle),4) + "            ";
                string[7] = getDeltaEpoch(sat.epochMjd, setEpoch) + "            ";
                string[8] = modifyNumber(sat.array[i].rangeObs/1000,3) + "            ";
                string[9] = modifyNumber(sat.array[i].rangeCenter/1000,3) + "            ";
                string[10] = modifyNumber(sat.array[i].visMag1,1) + "            ";
                string[11] = modifyNumber(sat.array[i].absMag1,1) + "            ";
                string[12] = modifyNumber(sat.array[i].visMag2,1) + "            ";
                string[13] = modifyNumber(sat.array[i].absMag2,1) + "            ";
                string[14] = sat.array[i].isInShadow + "                   ";

                string[19] = string[0].substring(0,12) + " " + string[1].substring(0,12) + " " + string[2].substring(0,12) + " " +
                        string[3].substring(0,12) + " " + string[4].substring(0,12) + " " + //string[5].substring(0,12) + " " +
                        string[6].substring(0,12) + " " + string[7].substring(0,10) + " " + string[8].substring(0,12) + " " +
                        string[9].substring(0,12) + " " + string[10].substring(0,6) + " " + string[11].substring(0,6) + " " +
                        string[12].substring(0,6) + " " + string[13].substring(0,6) + " " + string[14].substring(0,7);
                pw.write(string[19] + "\n");           
            }
             */


            float setEpoch = (float)time.getMjd(time);
            mainString = "123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890" + "\n";
            mainString = mainString + "----- NORAD# " + sat.catalogNo + " -----" + "\r\n\r\n";
            mainString = mainString + sat.name + "\r\n";
            mainString = mainString + sat.line1 + "\r\n" + sat.line2 + "\r\n\r\n";
            mainString = mainString + ("Object int ID:                ").substring(0,20) + sat.intDes + "\r\n\n";
            mainString = mainString + ("Start time:                ").substring(0,20) + getStringDateTime2(time) + " (UTC)\r\n";
            mainString = mainString + ("Start time:                ").substring(0,20) + Time.getMjd(time) + " days (MJD)\r\n";
            mainString = mainString + ("Observatory:               ").substring(0,20) + modifyNumber(Math.toDegrees(observer.lon),4) + " deg, " +
            modifyNumber(Math.toDegrees(observer.lat),4)+ " deg, " + observer.altitude + " m" +"\r\n";
            mainString = mainString + "Exp time:                  ".substring(0,20) + exp + " s" + "\r\n";
            mainString = mainString + "Visual mag:                ".substring(0,20) + modifyNumber(sat.array.movingVisMag1,1) + " ; " + modifyNumber(sat.array.movingVisMag1,1) + "\r\n";
            mainString = mainString + "Standard mag:              ".substring(0,20) + modifyNumber(sat.array.absMag1,1) + "\r\n";
            mainString = mainString + "Angular speed:             ".substring(0,20) + modifyNumber(Math.toDegrees(sat.array.angularSpeed)*3600,3) + " arcmin/min" + "\r\n";
            mainString = mainString + "Position angle:            ".substring(0,20) + modifyNumber(Math.toDegrees(sat.array.positionAngle),3) + " deg" + "\r\n";
            pw.write(mainString + "\r\n");
            String string[] = new String[20];
            string[0] = "Date                                            ";
            string[1] = "Az                  ";
            string[2] = "h                   ";
            string[3] = "R.A.                ";
            string[4] = "dec                  ";
            string[5] = "Ma(t)               ";
            //string[6] = "v(t)               ";      //not yet
            string[7] = "Phase               ";
            //string[8] = "T since            ";
            string[9] = "Obs. r.              ";
            string[10] = "Angular             ";
            //string[11] = "V1                        ";
            //string[12] = "H1                        ";
            //string[13] = "V2                        ";
            //string[14] = "H2                        ";
            string[15] = "Shadow                    ";

            string[19] = string[0].substring(0,25) + " " + string[1].substring(0,10) + " " + string[2].substring(0,10) + " " +
                        string[3].substring(0,10) + " " + string[4].substring(0,10) + " " + //string[5].substring(0,12) + " " +
                        //string[6].substring(0,12) + " " +
                        string[7].substring(0,6) + " " + //string[8].substring(0,12) + " " +
                        string[9].substring(0,10) + " " + string[10].substring(0,12) + " " + // string[11].substring(0,6) + " " +
                        //string[12].substring(0,6) + " " + string[13].substring(0,6) + " " + string[14].substring(0,7) + " " +
                        string[15].substring(0,6);

            pw.write(string[19] + "\r\n");
            string[19] =("(UTC)                                                 ").substring(0,25) + " " + ("[deg]                     ").substring(0,10) + " " +
                        ("[deg]                     ").substring(0,10) + " " + ("[deg]                     ").substring(0,10) + " " +
                        ("[deg]                     ").substring(0,10) + " " + //("[deg]                     ").substring(0,12) + " " +
                        ("[deg]                     ").substring(0,6) + " " + //("[day]                     ").substring(0,10) + " " +
                        ("[km]                     ").substring(0,10) + " " +  ("[arcm/min]                     ").substring(0,16) + " " +
                        //("                     ").substring(0,6) + " " + ("                   ").substring(0,6) + " " +
                        //("                     ").substring(0,6) + " " + ("                     ").substring(0,6) + " " +
                        ("                     ").substring(0,6);
            pw.write(string[19] + "\r\n\r\n");
            //mainString = mainString + "\n" + string[19] + "\n";
            //pw.write(mainString);

            //1st line
            //for(int i = 0; i <= exp; i++){
            
                string[0] = getStringDateTime2(time) + "                                           ";
                string[1] = modifyNumber(Math.toDegrees(sat.array.azimuth),5) + "              ";
                string[2] = modifyNumber(Math.toDegrees(sat.array.altitude),5) + "            ";
                string[3] = modifyNumber(Math.toDegrees(sat.array.ra),5) + "            ";
                string[4] = modifyNumber(Math.toDegrees(sat.array.dec),5) + "            ";
                string[5] = modifyNumber(Math.toDegrees(sat.array.actMeanAnom),4) + "            ";
                //string[6] = "                                    ";
                string[7] = modifyNumber(Math.toDegrees(sat.array.phaseAngle),4) + "            ";
                //string[8] = getDeltaEpoch(sat.epochMjd, setEpoch) + "            ";
                string[9] = modifyNumber(sat.array.rangeObs/1000,3) + "            ";
                //string[10] = modifyNumber(Vector.getSize(sat.array.positionVector)/1000,3) + "            ";
                string[10] = (modifyNumber(Math.toDegrees(sat.array.angularSpeed)*3600,3)) + "            ";
                //string[11] = modifyNumber(sat.array[i].visMag1,1) + "            ";
                //string[12] = modifyNumber(sat.array[i].absMag1,1) + "            ";
                //string[13] = modifyNumber(sat.array[i].visMag2,1) + "            ";
                //string[14] = modifyNumber(sat.array[i].absMag2,1) + "            ";
                string[15] = sat.array.isInShadow + "                   ";

                string[19] = string[0].substring(0,25) + " " + string[1].substring(0,10) + " " + string[2].substring(0,10) + " " +
                        string[3].substring(0,10) + " " + string[4].substring(0,10) + " " + //string[5].substring(0,12) + " " +
                        //string[6].substring(0,12) + " " +
                        string[7].substring(0,6) + " " + //string[8].substring(0,12) + " " +
                        string[9].substring(0,10) + " " + string[10].substring(0,12) + " " + //string[11].substring(0,6) + " " +
                        //string[12].substring(0,6) + " " + string[13].substring(0,6) + " " + string[14].substring(0,7) + " " +
                        string[15].substring(0,6);
                pw.write(string[19] + "\r\n");
            //}
            
            for(int i = 0; i <= seda.amount-1; i++){
                double time_2 = seda.sed[i].epochMjd;
                //System.out.println("angular " + Math.toDegrees(seda.sed[i].angularSpeed)*3600);
                //string[0] = (i+1) + "            ";
                //string[0] = getStringDateTime2(Time.getDateTime(Time.getMjd(time) + (double)(i+1)/(86400))) + "            ";
                string[0] = getStringDateTime2(Time.getDateTime(time_2)) + "            ";
                string[1] = modifyNumber(Math.toDegrees(seda.sed[i].azimuth),5) + "            ";
                string[2] = modifyNumber(Math.toDegrees(seda.sed[i].altitude),5) + "            ";
                string[3] = modifyNumber(Math.toDegrees(seda.sed[i].ra),5) + "            ";
                string[4] = modifyNumber(Math.toDegrees(seda.sed[i].dec),5) + "            ";
                string[5] = modifyNumber(Math.toDegrees(seda.sed[i].actMeanAnom),4) + "            ";
                //string[6] = "                                    ";
                string[7] = modifyNumber(Math.toDegrees(seda.sed[i].phaseAngle),4) + "            ";
                //string[8] = getDeltaEpoch(sat.epochMjd, setEpoch) + "            ";
                string[9] = modifyNumber(seda.sed[i].rangeObs/1000,3) + "            ";
                //string[10] = modifyNumber(Vector.getSize(seda.sed[i].positionVector)/1000,3) + "            ";
                string[10] = modifyNumber(Math.toDegrees(seda.sed[i].angularSpeed)*3600,3) + "            ";
                //string[11] = modifyNumber(sat.array[i].visMag1,1) + "            ";
                //string[12] = modifyNumber(sat.array[i].absMag1,1) + "            ";
                //string[13] = modifyNumber(sat.array[i].visMag2,1) + "            ";
                //string[14] = modifyNumber(sat.array[i].absMag2,1) + "            ";
                string[15] = seda.sed[i].isInShadow + "                   ";

                 string[19] = string[0].substring(0,25) + " " + string[1].substring(0,10) + " " + string[2].substring(0,10) + " " +
                        string[3].substring(0,10) + " " + string[4].substring(0,10) + " " + //string[5].substring(0,12) + " " +
                        //string[6].substring(0,12) + " " +
                        string[7].substring(0,6) + " " + //string[8].substring(0,12) + " " +
                        string[9].substring(0,10) + " " + string[10].substring(0,12) + " " + //string[11].substring(0,6) + " " +
                        //string[12].substring(0,6) + " " + string[13].substring(0,6) + " " + string[14].substring(0,7) + " " +
                        string[15].substring(0,6);
                pw.write(string[19] + "\r\n");
                //System.out.println(i);
            }
            
            pw.write("\r\n" + "SatEph");
            pw.flush();
            pw.close(); 
        }
        catch (Exception e){
            System.out.println(e.getMessage());
        }

        return mainString;
    }        
   
    /*
    * getStringDateTime
    *
    * Metoda sluzi na prevod datumu a casu na text. Vystupna hodnota bude String yyyy/mm/dd hh:mm:ss.s
    * Method to get the actual set time in String
    *  In: 	Time()
    *  Out:	String yymmddhhmmss
    */	
	public static String getStringDateTime(Time time){
		String stringTime = time.year+"";
                //System.out.println("00 " + time.year);
                stringTime = stringTime.substring(2,4);
                //System.out.println("01 " + stringTime);
		
		if(time.month>9) stringTime=stringTime+time.month+"";
			else stringTime=stringTime+"0"+time.month+"";
                 //System.out.println("02 " + stringTime);
		if(time.day>9) stringTime=stringTime +time.day+"_";
			else stringTime=stringTime +"0"+time.day+"_";
		if(time.hour>9) stringTime=stringTime+time.hour+"";
			else stringTime=stringTime +"0"+time.hour+"";
		if(time.min>9) stringTime=stringTime+time.min+"";
			else stringTime=stringTime +"0"+time.min+"";
		if(time.sec>9) stringTime=stringTime+Math.round(time.sec);
			else stringTime=stringTime+"0"+Math.round(time.sec);
		//System.out.println(stringTime);
		return stringTime;
	}
        
    /*
    * getStringDateTime2
    *
    * Metoda sluzi na prevod datumu a casu na text. Vystupna hodnota bude String yyyy/mm/dd hh:mm:ss.s
    * Method to get the actual set time in String
    *  In: 	Time()
    *  Out:	String yyyy-mm-dd hh:mm:ss.s
    */	
	public static String getStringDateTime2(Time time){
		String stringTime = time.year+"-";
		
		if(time.month>9) stringTime=stringTime+time.month+"-";
			else stringTime=stringTime+"0"+time.month+"-";
		if(time.day>9) stringTime=stringTime +time.day+" ";
			else stringTime=stringTime +"0"+time.day+" ";
		if(time.hour>9) stringTime=stringTime+time.hour+":";
			else stringTime=stringTime +"0"+time.hour+":";
		if(time.min>9) stringTime=stringTime+time.min+":";
			else stringTime=stringTime +"0"+time.min+":";
		if(time.sec>=10) stringTime=stringTime+(float)(time.sec);
			else stringTime=stringTime+"0"+(float)(time.sec);
		//System.out.println("2 " + stringTime);
		return stringTime;
	}
        
     /**
     *METHOD to getting days since set epoch
     *
     * IN:  double set epoch mjd [days]
     *      double set epoch    [days]
     *
     * OUT: double difference
     */
    
    public double getDeltaEpoch(double tleEpoch, double setEpoch){
        double deltaEpoch = setEpoch - tleEpoch;
        return deltaEpoch;
    }
    
    /**
     *METHOD modifiing the double number
     *
     * IN:  double number
     *      decimal position
     *
     * OUT: double number2
     */
    public double modifyNumber(double number, int decPos){
        double decPosition = Math.pow(10, decPos); 
        number = (float)(Math.round(number*decPosition))/decPosition;
        return number;
    }

    /**
     * getSettingsInfoString()
     *
     * IN:
     *  Time time   -   obs time
     *  double exp  -   exp time [s]
     *  Geodetic obsPos - observer position
     *  boolean isEqCoordinates - Hor/Eq coordinates
     *  Vector coordinates - hor or eq coordinates [deg]
     *  Vectror fieldSize - size of field, width, height [deg]
     *
     * OUT:
     *  String stringToSave - string with all the relevant data
     */
    public String getSettingsString(double mjdTime, double exp, Geodetic obsPos,
                boolean isEqCoordinates, Vector coordinates, Vector fieldSize){
        String stringToSave = "SatEph setting file";
        stringToSave = stringToSave + "\r\n\r\ntime (UTC):" + "\t\t" + getStringDateTime2((new Time()).getDateTime(mjdTime));
        stringToSave = stringToSave + "\r\ntime (MJD):" + "\t\t" + mjdTime;
        stringToSave = stringToSave + "\r\nexp [s]:" + "\t\t" + exp;
        stringToSave = stringToSave + "\r\nobs pos:" + "\t\t" + getRightString(Math.toDegrees(obsPos.lon) + "",8) + "\t" +
                    getRightString(Math.toDegrees(obsPos.lat) + "",8) + "\t" + getRightString(obsPos.altitude+ "", 8);
        String typeOfCoor = "Hor";
        if(isEqCoordinates) typeOfCoor = "Equ";
        stringToSave = stringToSave + "\r\nfield pos [deg]:" + "\t" + typeOfCoor + "\t\t" +
                    getRightString(Math.toDegrees(coordinates.v[0]) + "", 9) + "\t" +
                    getRightString(Math.toDegrees(coordinates.v[1]) + "", 9);
        stringToSave = stringToSave + "\r\nfield size [deg]:" + "\t" + fieldSize.v[0] + "\t\t" +
                    fieldSize.v[1] + "\t\t" + fieldSize.v[2];

        //System.out.println(stringToSave);
        return stringToSave;
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

    public String getRightString(String stringToRepair, int lengthOfString){
        String rightString = "";
        int maxLength = 100;
        if(lengthOfString > maxLength) lengthOfString = maxLength;
        rightString = (stringToRepair + "                                                                      " +
                "                                        ").substring(0, lengthOfString);

        return rightString;
    }
    
    /**
     * saveTLEFile()
     * 
     * Method to save TLEs from SatelliteArray to TLE file
     * 
     * IN:
     * SatelliteArray sa3 - array with satellites info
     * String fileName - name of file
     */
    
    public void saveTLEFile(SatelliteArray sa3, String fileName){
        String stringTLE = "";
        
        try{
            PrintWriter pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream("data//" + fileName + ".txt"), "UTF-8"));

            //fill String with TLEs
            for(int i = 0; i < sa3.largestNo; i++){
                //if(sa3.array[i].epoch != 0){
                if(sa3.array[i].epochMjd != 0){
                    stringTLE = sa3.array[i].line1 + "\n";
                    stringTLE = stringTLE + sa3.array[i].line2 + "\n";         
                    pw.write(stringTLE);
                }
            }            
            pw.flush();
            pw.close(); 
        }
        catch(Exception e){
            System.out.println(e.getMessage());
        }
        
        //saveFileString(stringTLE, fileName);
    }
   
    /**
     * getZMDFormat()
     * Method: method to modify the angle into the AIUB EPH input file formats
     *
     * IN:  double angle - angle (RA or DEC)
     *      
     * OUT: String string - new angle format in form of a String
     *  
     */
    String getZMDFormat(double angle){
        String string = "";
        double deg, min, sec;
        String deg_string, min_string, sec_string,sec_string2;
        
        //which sign the number has, + or -
        String sign = "+";
        if(angle < 0) {
            sign = "-";
            angle = Math.abs(angle);
        }
        
        deg = (int)angle;
        min = (int)((angle - deg)*60);
        sec = (((angle - deg)*60) - min)*60;        
        
        deg_string = deg + "";
        if(deg < 10){
            deg_string = deg_string.substring(0,1);
            deg_string = "0" + deg_string;
        }
        else if((deg < 100)&&(deg >= 10)){
            deg_string = deg_string.substring(0,2);
            //deg_string = "0" + deg_string;
        }
        else{
            deg_string = deg_string.substring(0,3);            
        }
        if(sign.equals("-")){
            deg_string = sign + deg_string;
        }
    
        min_string = min + "";
        if(min < 10){
            min_string = min_string.substring(0,1);
            min_string = "0" + min_string;
        }
        else if((min < 100)&&(min >= 10)){
            min_string = min_string.substring(0,2);
            //min_string = "0" + min_string;
        }
        else{
            min_string = min_string.substring(0,3);            
        }
        
        sec_string = sec + "";
        if(sec < 10){
            //min_string = min_string.substring(0,1);
            sec_string2 = sec_string.substring(2,sec_string.length());
            sec_string = sec_string.substring(0,1) + sec_string2;
            sec_string = "0" + sec_string;
        }
        else{ // if((sec < 100)&&(sec >= 10)){
            sec_string2 = sec_string.substring(3,sec_string.length());
            sec_string = sec_string.substring(0,2) + sec_string2;
            //sec_string = "0" + sec_string;
        }
        //else{
            //sec_string = sec_string.substring(0,3);            
        //}
        
        string = deg_string + "." + min_string + sec_string;
        
        return string;
    }
}
