/*
 * Class to read files with ephemerides data - time, R.A., dec
 */

package OrbitDetermination;

import java.io.*;
import java.util.*;

/**
 *
 * @author Jiri Silha - 17.01.2010
 */
public class ReadFileOD {

     //Variables
    static BufferedReader file;    //read the file

    static StringTokenizer st;      //to separete the line for smaller parts

    //on which OS is application running
    String osName = System.getProperty("os.name");

    /**
     * getFileLines()
     *
     * Method to read how many lines the source file got.
     *
     * IN:
     * String filePos
     *
     * OUT:
     * int noOfLines
     */
    public int getFileLines(String filePos){
        int noOfLines = 0;

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

        //do cycle
        boolean doCycle = true;
        //string from given line of reading file
        String line = "";

        //try read file
        try{
            //read file from specify position
            file = new BufferedReader(new FileReader(filePos));
            //read line
            while(doCycle){
                try{
                    line = file.readLine();
                    if(line.equals("null")) break;
                }
                catch(Exception e){
                    file.close();
                    //System.out.println("End of file!");
                    break;
                }
                //if(line.equals("null")) break;
                noOfLines++;
            }
            file.close();
        }
        catch (Exception e){
            System.out.println("Warning! " + "End of file!");
            System.out.println (e);
        }
        return noOfLines;
    }

    /**
     * getTheStrings()
     *
     * Method to get all the lines from given file --> Strings
     *
     * IN:
     * String filePos - file position
     *
     * OUT:
     * String linesInStrings[] - array of Strings with given lines;
     */
    public String[] getTheStrings(String filePos){
        int noOfLines = getFileLines(filePos);
        String linesInStrings[] = new String[noOfLines];

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

        //do cycle
        boolean doCycle = true;
        //string from given line of reading file
        String line = "";

        noOfLines = 0;
        //try read file
        try{
            //read file from specify position
            file = new BufferedReader(new FileReader(filePos));
            //read line
            while(doCycle){
                try{
                    line = file.readLine();
                    if(line.equals("null")) break;
                    linesInStrings[noOfLines] = line;
                }
                catch(Exception e){
                    file.close();
                    //System.out.println("End of file!");
                    break;
                }
                //if(line.equals("null")) break;
                noOfLines++;
            }
            file.close();
        }
        catch (Exception e){
            System.out.println("Warning! " + "End of file!");
            System.out.println (e);
        }
        
        return linesInStrings;
    }
}
