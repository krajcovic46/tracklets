/**
 * Class to work with packed files - zip format.
 * Source - http://www.java-tips.org/java-se-tips/java.util.zip/how-to-extract-file-files-from-a-zip-file-3.html
 */

package InOut;

import java.io.*;
import java.util.zip.*;

public class  ZipFiles
{
    public static void main(String[] args) 
    {
        try
        {
            //String filename = args[0];
            //String filename = "ALL_TLE.ZIP";
            //String filename = "classfd.ZIP";
            String filename = "data/downloaded/tle/all_tle.zip";
            ZipFiles list = new ZipFiles( );
            list.getZipFiles(filename);
        }
        catch (Exception e)
        {
            //e.printStackTrace();
            System.out.println("ZipFiles " + e.getMessage());
        }
    }

    public void getZipFiles(String filename)
    {
        System.out.println("Reading " + filename);
        try
        {
            String destinationname = "";
            byte[] buf = new byte[1024];
            ZipInputStream zipinputstream = null;
            ZipEntry zipentry;
            zipinputstream = new ZipInputStream(
                new FileInputStream(filename));

            zipentry = zipinputstream.getNextEntry();
            File zipFile = new File(filename);
            String zipFilePos = zipFile.getParent();
            //System.out.println(zipFilePos);

            while (zipentry != null) 
            { 
                //for each entry to be extracted
                String entryName = zipentry.getName();
                System.out.println("Unpacking " + entryName + " from " + filename);
                int n;
                FileOutputStream fileoutputstream;
                File newFile = new File(zipFilePos + "/" + entryName);
                String directory = newFile.getParent();
                //System.out.println(directory);
                
                if(directory == null)
                {
                    if(newFile.isDirectory())
                        break;
                }
                
                //fileoutputstream = new FileOutputStream(destinationname+entryName);
                //change name of file ssrDATE.txt to ssr.txt
                if(entryName.indexOf("ssr") != -1) entryName = "ssr.txt";
                fileoutputstream = new FileOutputStream(directory + "/" + entryName);
                System.out.println("Saving " + entryName + " to " + directory + "/" + entryName);

                while ((n = zipinputstream.read(buf, 0, 1024)) > -1)
                    fileoutputstream.write(buf, 0, n);

                fileoutputstream.close(); 
                zipinputstream.closeEntry();
                zipentry = zipinputstream.getNextEntry();

            }//while

            zipinputstream.close();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }
}