/*
 * FileDownload.java
 *
 * Created on Utorok, 2008, april 15, 2:00
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 *
 * NOTES:  Source: http://schmidt.devlib.org/java/file-download.html 
 */

package com.skrajcovic.orbitdetermination;

import java.io.*;
import java.net.*;

/**
 *
 * @Author Jiri Silha 
 */
public class FileDownload {
    
    /** Source: http://schmidt.devlib.org/java/file-download.html */
    public FileDownload() {
    }
    
    /**
     * download()
     *
     * IN:  String adress           - url addres of file
     *      String localFileName    - local adres of saving place - file 
     *
     */
    
    public static void download(String address, String localFileName) {
        OutputStream out = null;
        URLConnection conn = null;
        InputStream  in = null;
        System.out.println("Downloading from " + address + " to " + localFileName);
	try {
            URL url = new URL(address);
            out = new BufferedOutputStream(
            new FileOutputStream(localFileName));
            conn = url.openConnection();
            in = conn.getInputStream();
            byte[] buffer = new byte[1024];
            int numRead;
            long numWritten = 0;
            while ((numRead = in.read(buffer)) != -1) {
                out.write(buffer, 0, numRead);
                numWritten += numRead;
            }
     } 
        catch (Exception exception) {
            System.out.println(exception);
            System.out.println("Can not download from " + address);
        } 
        finally {
            try {
		if (in != null) {
                    in.close();
                }
		if (out != null) {
                    out.close();
		}
            } 
            catch (IOException ioe) {}
        }
    }

    /**
     * download()
     *
     * IN:  String address           - url address of file
     *      String address2          - url of backup address of file
     *      String localFileName    - local adres of saving place - directory, and subdirectories
     *      String localFileName    - file name
     *
     */

    public static void download(String address, String address2, String fileLocation, String localFileName) {
        OutputStream out = null;
        URLConnection conn = null;
        InputStream  in = null;
        createDirectory(fileLocation);
        System.out.println("Downloading from " + address + " to " + fileLocation + localFileName);
	try {
            localFileName = fileLocation + localFileName;
            URL url = new URL(address);
            out = new BufferedOutputStream(
                new FileOutputStream(localFileName));
            conn = url.openConnection();
            in = conn.getInputStream();
            byte[] buffer = new byte[1024];
            int numRead;
            long numWritten = 0;
            while ((numRead = in.read(buffer)) != -1) {
                out.write(buffer, 0, numRead);
                numWritten += numRead;
                //System.out.println("numWritten 1" + numWritten);
            }

            //if object has very small size, try other TLEs sources
            if(numWritten>100000){
                System.out.println("Can not download from " + address);
                System.out.println("Trying to download from " + address2);
                try{
                    url = new URL(address2);
                    out = new BufferedOutputStream(
                    new FileOutputStream(localFileName));
                    conn = url.openConnection();
                    in = conn.getInputStream();
                    numWritten = 0;
                    while ((numRead = in.read(buffer)) != -1) {
                        out.write(buffer, 0, numRead);
                        numWritten += numRead;
                        //System.out.println("numWritten 2" + numWritten);
                    }
                }
                catch(Exception exception2) {
                    System.out.println(exception2);
                    System.out.println("Can not download from " + address);
                }
            }
        }
        catch (Exception exception) {
            System.out.println(exception);
            System.out.println("Can not download from " + address);
            System.out.println("Trying to download from " + address2);
            try{
                URL url = new URL(address2);
                out = new BufferedOutputStream(
                new FileOutputStream(localFileName));
                conn = url.openConnection();
                in = conn.getInputStream();
                byte[] buffer = new byte[1024];
                int numRead;
                long numWritten = 0;
                while ((numRead = in.read(buffer)) != -1) {
                    out.write(buffer, 0, numRead);
                    numWritten += numRead;
                    //System.out.println("numWritten 2" + numWritten);
                }
            }
            catch(Exception exception2) {
                System.out.println(exception2);
                System.out.println("Can not download from " + address);
            }
        }
        finally {
            try {
		if (in != null) {
                    in.close();
                }
		if (out != null) {
                    out.close();
		}
            }
            catch (IOException ioe) {}
        }
    }

    public static void download(String address) {
        int lastSlashIndex = address.lastIndexOf('/');
        if (lastSlashIndex >= 0 && lastSlashIndex < address.length() - 1) {
            download(address, address.substring(lastSlashIndex + 1));
        } 
        else {
            System.err.println("Could not figure out local file name for " + address);
        }
    }
        
     /**
     *Download TLE files from www.celestrak.com
     *
     */
     public static void dowloadKelsoTLE(){
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

     String fileUrlPos = "http://www.celestrak.com/NORAD/elements/";
     String filePos = "";
     String filePos2 = "data/downloaded/kelso/all/";
     createDirectory(filePos2);
     int amountHelp = 0;
     //filling array
     for(int i = 0; i < fileName.length; i++){
         filePos = fileUrlPos + fileName[i] + ".txt";
         download(filePos, filePos2 + fileName[i] + ".txt");
     }
   }

     /**
     *Download supplemental TLE files from www.celestrak.com
     *
     */
     public static void dowloadKelsoTLESuppl(){
     //array of TLE file name
     String fileName[] = {"gps", "glonass", "meteosat", "intelsat", "orbcomm"};

     String fileUrlPos = "http://www.celestrak.com/NORAD/elements/supplemental/";
     String filePos = "";
     String filePos2 = "data/downloaded/kelso/supplemental/";
     createDirectory(filePos2);
     int amountHelp = 0;
     //filling array
     for(int i = 0; i < fileName.length; i++){
         filePos = fileUrlPos + fileName[i] + ".txt";
         download(filePos, filePos2 + fileName[i] + ".txt");
     }
   }

     /**
     * createDirectory()
     *
     * IN: String directoryName - name of the directory
     *
     */

    public static void createDirectory(String dirName){
        try{
            boolean success = (new File(dirName)).mkdir();
            //System.out.println("Save directory");
        }
        catch(Exception e){
            System.out.println(e.getMessage());
        }
    }
}
