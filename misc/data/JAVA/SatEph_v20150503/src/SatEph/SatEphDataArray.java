/*
 * This class is supposed to hold information about satellite position data.
 * It is array.
 */

package SatEph;

/**
 *
 * @author Jiri
 */
public class SatEphDataArray {
    /**
     * Size of array
     */
    int sizeOfArray = 6000;

    /**
     * Array to hold the dynamic data of satellite
     */
    public SatEphData sed[] = new SatEphData[sizeOfArray];

    /**
     * The no. of data we want to use - exposition time, max value is 6000
     */
    public int amount;

    /**
     * getDelta()
     *
     * Method to compute delta time between 2 positions, it depends on exposition time
     * and array size - sizeOfArray
     *
     * INPUT:
     *  int expTime - exposition time in [s]
     *
     * OUTPUT:
     *  double delta time - delta [s] between 2 positions
     */
    public double getDelta(int expTime){
        double delta;
        if(expTime > sizeOfArray) {
            delta = (double)expTime/(double)sizeOfArray;
        }
        else {
            delta = 1;
        }
        //else
        return delta;
    }

}
