/*
 * SatelliteArray.java
 *
 * Created on February 28, 2008, 2:45 PM
 *
 * SK: Trieda sluzi na ukadanie pola premennych Satellite
 * EN: Class for Satellite variables array
 */

package SatEph;

/**
 *
 * @author Jiri Silha
 */
public class SatelliteArray {
    /**
     * array
     */
    public Satellite array[] = new Satellite[100000];
    
    /**
     * no of objects with TLE datas
     */
    public int tleAmount;
    
    /**
     * no of objects with Kelso TLE datas
     */
    public int tleKelsoAmount;
    
    /**
     * no of all objects
     */
    public int allAmount;
    
    /**
     * date of satellite report
     */
    public String reportDate;
    
    /**
     * largest no. in report
     */
    public int largestNo = 0;

    /**
     * no of objects on Earth orbit, included non TLE objects
     */
    public int objectsOnOrbit = 0;
    
    /** Creates a new instance of SatelliteArray */
    public SatelliteArray() {
        for(int i = 0; i<array.length;i++){
            array[i] = new Satellite();
            array[i].catalogNo = i;
        }
    }
    
}
