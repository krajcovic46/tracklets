/*
 * SatelliteArray.java
 *
 * Created on February 28, 2008, 2:45 PM
 *
 * EN: Class for Stars variables array
 */

package SatEph;

/**
 *
 * @author Jiri Silha
 */
public class StarArray {
    /**
     * Array
     */
    public Star array[] = new Star[1000];    
    /**
     * largest no. in report
     */
    public int largestNo = 0;
   
    /** Creates a new instance of SatelliteArray */
    public StarArray() {
        for(int i = 0; i<array.length;i++){
            array[i] = new Star();
        }
    }
    
}
