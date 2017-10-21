//Jiri Silha - 09/08/08
// methods in class are supposed to convert string to no, or no to string

package Compute.Silha;

//import java.awt.*;
import javax.swing.*;

public class ReadCoordinates{
    
     /**
     * On which OS is application running
     */
    static String osName = System.getProperty("os.name");
        
    /**
     * Sign for deg
     */
    static String stringDeg = "";
    
	/**
	* Method: convert string to no.,
	* 
	* IN: String (deg min sec)||(h m s), or (deg.deg...)||(h.hhh...)
	*
	* OUT: double
	*/
	
	public static double getDegrees(String string, JFrame frame){
		double degrees = 0;
		boolean isNoPositive = true;
		double degrees2[] = new double[3];
		//there are 2 format types (deg min sec), or (deg.deg...), example: 141 52 63,8 (means 141 Deg 52 min 63,8 sec), or 141,88438...
		
		//modify the string
		//"," => "."
		string = string.replaceAll(",",".");
		//see the methods
		string = returnClearString(string);
		string = returnClearStringEnd(string);
		
		// have to find out wich format it is
		int help = string.indexOf(" ");
		
		//is number positive
		if(string.indexOf("-") != -1) isNoPositive = false;
		
		//condition
		//normal way
		if(help == -1) {
			try{
				degrees = Double.parseDouble(string);
			}
			catch( NumberFormatException nfe) {
                            //System.out.println(nfe);
                            JOptionPane.showMessageDialog(
                            frame,"" + nfe, "Invalid Number Format", JOptionPane.ERROR_MESSAGE
                        );
        }
		}
		//method will help
		else {
			degrees2 = getDegrees2(string);
			//to degrees
			if(isNoPositive){
				degrees = degrees2[0] + degrees2[1]/60 + degrees2[2]/3600;
			}
			else {
				degrees = degrees2[0] - degrees2[1]/60 - degrees2[2]/3600;
			}
		}
			
		return degrees;
	}
	
	/*
	* Method: getDegrees2 
	*
	* IN: String - format => "deg min sec.sec..."
	*
	* OUT : double
	*/
	public static double[] getDegrees2(String string){
		String deg, min, sec;
		double degrees[] = new double[3];
		
		//getting the values
		deg = string.substring(0, string.indexOf(" "));
		string = string.substring(string.indexOf(" ") + 1, string.length());
		min = string.substring(0, string.indexOf(" "));
		string = string.substring(string.indexOf(" "), string.length());
		sec = string;
		try{	
			degrees[0] = Math.round(Double.parseDouble(deg));
			degrees[1] = Math.round(Double.parseDouble(min));
			degrees[2] = Double.parseDouble(sec);
		}
		catch (Exception e) {
			System.out.println(e);
		}
		
		return degrees;
	}
	
	/*
	* Method: StringModify.java
	*
	*
	*/
	
	//metoda sluzi na vytvorenie stringu ktory nema na zaciatku medzery "   Jano"=>"Jano"
	public static String returnClearString(String string){
		while(string.substring(0,1).equals(" ")) string=string.substring(1,string.length());
		return string;
	}
	
	//metoda sluzi na vytvorenie stringu ktory nema na zaciatku medzery "Jano   "=>"Jano"
	public static String returnClearStringEnd(String string){
		while(string.substring((string.length()-1),string.length()).equals(" ")) string=string.substring(0,(string.length()-1));
		return string;
	}
	
	/*
	* Method: get from double the string Deg Min Sec.sec...
	*
	* IN: double degrees -> deg.deg....
	*
	* OUT: String Deg° Min' Sec.sec''
	*/
	public static String getDegreesString(double degrees){
		String string, deg, min, sec; 
                
                //which sign for deg should be used, there are 2 difference sign, one
                //for Windows, another for Linux - just using the ascii code, it is universal - 176
                char degChar = 176;
                stringDeg = degChar + "";
                
		boolean isNoPositive = true;
		if(degrees < 0) {
			isNoPositive = false;
			degrees = degrees * -1;
		}
		
		deg = (int)Math.floor(degrees) + "";
		min = (int)Math.floor((degrees - Math.floor(degrees)) * 60) + "";
		sec = (double)Math.round(((degrees - Math.floor(degrees)) * 60 - Math.floor((degrees - Math.floor(degrees)) * 60))*60*100)/(double)100  + "";
		string = deg + stringDeg + " " + min + "' " + sec + "''";
		
		if(isNoPositive == false) string = "-" + string;
			
		return string;
	}
	
	/*
	* Method: get from double the string Deg Min Sec.sec...
	*
	* IN: double degrees -> deg.deg....
	*
	* OUT: String Deg° Min' Sec.sec''
	*/
	public static String getHoursString(double degrees){
		String string, deg, min, sec;
                //degrees => hours
                degrees =  degrees / 15;
		boolean isNoPositive = true;
		if(degrees < 0) {
			isNoPositive = false;
			degrees = degrees * -1;
		}
		
		deg = (int)Math.floor(degrees) + "";
		min = (int)Math.floor((degrees - Math.floor(degrees)) * 60) + "";
		sec = (double)Math.round(((degrees - Math.floor(degrees)) * 60 - Math.floor((degrees - Math.floor(degrees)) * 60))*60*100)/(double)100  + "";
		string = deg + "h " + min + "m " + sec + "s";
		
		if(isNoPositive == false) string = "-" + string;
			
		return string;
	}
	
	public static void main (String args[]){
		String string = " -0 55 53.52";
		//String string = "  -141,88438 ";
		double number = getDegrees(string, new JFrame());
		System.out.println(number);
		
		System.out.println(getDegreesString(number));
	}
}