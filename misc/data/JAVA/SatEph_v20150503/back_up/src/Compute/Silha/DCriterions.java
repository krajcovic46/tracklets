/*
 * Dsh criterion - This class includes method to compute Soutworth - Hawkins criterion.
 * See f.e. G. B. Valsecchi,T. J. Jopek and Cl. Froeschlee - 
 * Meteoroid stream identi�cation: a new approach +- I. Theory (1999)
 * All formulas numbers are identical with formula numbers from this article Valsecchi(1999).
 */

package Compute.Silha;

/**
 *
 * @author Jiri Silha - 18-04-2009
 */
public class DCriterions {
    
    //Inside variables
    
    /**
     * Formula 1
     */
    double i_ab; 
    
     /**
     * Formula 1
     */
    double pi_ab;
    
    /**
     * Dsh criterion
     */
    //double d_sh;
    
    /**
     * getDshCriterion()
     * 
     * INPUT:
     * Orbital elements of 1st and 2nd objects
     *  double q_1, g_2 - perigee distance  [m] - [AU???]
     *  double e_1, e_2 - eccentricities    [-]
     *  double i_1, i_2 - inclinations      [rad]
     *  double p_1, p_2 - arguments of perigee [rad]
     *  double n_1, n_2 - longitude of ascending node [rad]
     */
    public double getDshCriterion(double q_1, double q_2, double e_1, double e_2,
                                  double i_1, double i_2, double p_1, double p_2,
                                  double n_1, double n_2){
        double d_sh;
        
        //Modified formula 2
        i_ab = 2 * Math.asin(0.5*Math.sqrt(
                Math.pow(2*Math.sin((i_1 - i_2)/2),2) +
                Math.sin(i_1)*Math.sin(i_2)*Math.pow(2*Math.sin((n_1 - n_2)/2),2)
                ));
                
        //Formula 3
        //sign +/-
        double k = 1;
        if(Math.abs(n_1 - n_2) > Math.PI) k = -1;
        else k = 1;
        
        pi_ab = p_1 - p_2 + k*2*Math.asin(Math.cos((i_1 + i_2)/2)*Math.sin((n_1 - n_2)/2)*
                1/Math.cos(i_ab/2));
                
        //S-H criterion
        d_sh = Math.sqrt((e_1 - e_2)*(e_1 - e_2) + (q_1 - q_2)*(q_1 - q_2)+
                Math.pow(2*Math.sin(i_ab/2),2) + 
                Math.pow((e_1 + e_2)/2*2*Math.sin(pi_ab/2),2));
        
        return d_sh;
    }
    
    /**
     * getDSilhaCriterion()
     * Modified Dsh criterion - omega, Omega not relevant anymore
     * INPUT:
     * Orbital elements of 1st and 2nd objects
     *  double q_1, g_2 - perigee distance  [m] - [AU???]
     *  double e_1, e_2 - eccentricities    [-]
     *  double i_1, i_2 - inclinations      [rad]
     *  double p_1, p_2 - arguments of perigee [rad]
     *  double n_1, n_2 - longitude of ascending node [rad]
     */
    public double getDSilhaCriterion(double q_1, double q_2, double e_1, double e_2,
                                  double i_1, double i_2, double p_1, double p_2,
                                  double n_1, double n_2){
        double d_sh;
        
        //Modified formula 2
        i_ab = 2 * Math.asin(0.5*Math.sqrt(
                Math.pow(2*Math.sin((i_1 - i_2)/2),2) +
                Math.sin(i_1)*Math.sin(i_2)*Math.pow(2*Math.sin((n_1 - n_2)/2),2)
                ));
                
        //Formula 3
        //sign +/-
        double k = 1;
        if(Math.abs(n_1 - n_2) > Math.PI) k = -1;
        else k = 1;
        
        pi_ab = 0;//p_1 - p_2 + k*2*Math.asin(Math.cos((i_1 + i_2)/2)*Math.sin((n_1 - n_2)/2)*
                //1/Math.cos(i_ab/2));
                
        //S-H criterion
        d_sh = Math.sqrt((e_1 - e_2)*(e_1 - e_2) + (q_1 - q_2)*(q_1 - q_2)+
                Math.pow(2*Math.sin(i_ab/2),2) + 
                Math.pow((e_1 + e_2)/2*2*Math.sin(pi_ab/2),2));
        
        return d_sh;
    }

    /**
     * getDhCriterion()
     *
     * INPUT:
     * Orbital elements of 1st and 2nd objects
     *  double q_1, g_2 - perigee distance  [m]/[AU]
     *  double e_1, e_2 - eccentricities    [-]
     *  double i_1, i_2 - inclinations      [rad]
     *  double p_1, p_2 - arguments of perigee [rad]
     *  double n_1, n_2 - longitude of ascending node [rad]
     */
    public double getDhCriterion(double q_1, double q_2, double e_1, double e_2,
                                  double i_1, double i_2, double p_1, double p_2,
                                  double n_1, double n_2){
        double d_h;

        //Modified formula 2
        i_ab = 2 * Math.asin(0.5*Math.sqrt(
                Math.pow(2*Math.sin((i_1 - i_2)/2),2) +
                Math.sin(i_1)*Math.sin(i_2)*Math.pow(2*Math.sin((n_1 - n_2)/2),2)
                ));

        //Formula 3
        //sign +/-
        double k = 1;
        if(Math.abs(n_1 - n_2) > Math.PI) k = -1;
        else k = 1;

        pi_ab = p_1 - p_2 + k*2*Math.asin(Math.cos((i_1 + i_2)/2)*Math.sin((n_1 - n_2)/2)*
                1/Math.cos(i_ab/2));

        //S-H criterion
        d_h = Math.sqrt((e_1 - e_2)*(e_1 - e_2) + 
                Math.pow((q_1 - q_2)*(q_1 - q_2)/((q_1 + q_2)*(q_1 + q_2)),2)+    //q in [m], or [AU] not relevant, q/q
                Math.pow(2*Math.sin(i_ab/2),2) +
                Math.pow((e_1 + e_2)/2*2*Math.sin(pi_ab/2),2));
        //System.out.println(Math.sqrt((e_1 - e_2)*(e_1 - e_2))+"");
        //System.out.println(Math.pow((q_1 - q_2)*(q_1 - q_2)/((q_1 + q_2)*(q_1 + q_2)),2)+"");
        //System.out.println(Math.pow(2*Math.sin(i_ab/2),2)+"");
        //System.out.println(Math.pow((e_1 + e_2)/2*2*Math.sin(pi_ab/2),2)+"");

        return d_h;
    }

    /**
     * getDCriterionSteel
     *
     * Method to compute D criterion between 2 orbits, omega and Omega are not importent, Steel 1991 and 1993
     *
     * IN:
     *      double q_1 - perigee of the 1st orbit [m]
     *      double q_2 - perigee of the 2nd orbit [m]
     *      double e_1 - eccentricity of the 1st orbit
     *      double e_2 - eccentricity of the 2nd orbit
     *      double i_1 - inclination of the 1st orbit [rad]
     *      double i_2 - inclination of the 2nd orbit [rad]
     *
     * OUT:
     *      double d_Steel
     */
    public double getDCriterionSteel(double q_1, double q_2, double e_1, double e_2,
                                    double i_1, double i_2){
        //factor for perigee, in original D criterion are q_1, q_2 in [AU]
        double g_factor = Constants.R_Earth*2;
        double d_Steel = Math.sqrt(Math.pow((q_1-q_2)/g_factor,2) + Math.pow((e_1-e_2),2) +
                         Math.pow((2*Math.sin((i_1-i_2)/2)),2));
        //System.out.println("Math.pow((q_1-q_2)/g_factor,2) " + Math.pow((q_1-q_2)/g_factor,2));
        //System.out.println("Math.pow((e_1-e_2),2) " + Math.pow((e_1-e_2),2));
        //System.out.println("Math.pow((2*Math.sin((i_1-i_2)/2)),2) " + Math.pow((2*Math.sin((i_1-i_2)/2)),2));
        return d_Steel;
    }
    
    /**
     * main method
     */
    public static void main(String args[]){
        /**
         * Orbital elements
         */
        //semi major axis
        double sma[] = new double[3];
        //eccentricity
        double ecc[] = new double[3];
        //perigee
        double per[] = new double[3];
        //inclination
        double inc[] = new double[3];
        //argument of perigee
        double arg[] = new double[3];
        //longitude of ascending node
        double nod[] = new double[3];
        //epoch - year must be the same [day]
        double epo[] = new double[3];
        
        /*
         * ISS TLE 1
         * 1 25544U 98067A   09104.08106641  .00010982  00000-0  86085-4 0  6684
         * 2 25544 051.6407 275.5284 0008674 236.5696 226.4625 15.72136261595864
         */
        sma[0] = 6730.930e3;
        ecc[0] = 0.0008674;
        per[0] = sma[0] * (1 + ecc[0]);
        inc[0] = Math.toRadians(051.6407);
        arg[0] = Math.toRadians(236.5696);
        nod[0] = Math.toRadians(275.5284);
        epo[0] = 104.08106641;
        
         /*
         * ISS TLE 2
         * 1 25544U 98067A   09091.68602787  .00012306  00000-0  96990-4 0  5909
         * 2 25544 051.6423 338.9854 0009269 184.1272 286.1737 15.71858385593916
         */
        sma[1] = 6731.723e3;
        ecc[1] = 0.0009269;
        per[1] = sma[1] * (1 + ecc[1]);
        inc[1] = Math.toRadians(051.6423);
        arg[1] = Math.toRadians(184.1272);
        nod[1] = Math.toRadians(338.9854);
        epo[1] = 91.68602787;
        
        /*
         * ISS TLE 3
         * 1 25544U 98067A   09106.49981230  .00010049  00000-0  79015-4 0  6843
         * 2 25544  51.6407 263.1420 0008692 245.9646 235.8437 15.72185824596245
         */
        sma[2] = 6730.788e3;
        ecc[2] = 0.0008692;
        per[2] = sma[2] * (1 + ecc[2]);
        inc[2] = Math.toRadians(51.6407);
        arg[2] = Math.toRadians(245.9646);
        nod[2] = Math.toRadians(263.1420);
        epo[2] = 106.49981230;
        
        //
        int i = 2;
        int j = 0;
        
        //D_sh criterion
        double d_shTest = new DCriterions().getDshCriterion(per[i]/new Constants().AU, per[j]/new Constants().AU, ecc[i], ecc[j], inc[i], inc[j],
                                        arg[i], arg[j], nod[i], nod[j]);
        
        System.out.println("dt: " + Math.abs(epo[i] - epo[j]) + " days");
        System.out.println("ISS D_sh: " + d_shTest);
    }
}
