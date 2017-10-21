/*
 * Class to compute circular orbit from two celestial positions.
 */

package OrbitDetermination;

import Compute.Silha.*;

/**
 *
 * @author Jiri Silha - 24.04.2010
 */
public class CircularOrbitDetermination {

    public static void main(String[] args){

        /**
         * VFMO100422
         */
        /*
        //assumption that observer - body distance is long, d > 20000 km
        double angSpeed = 8.05; //[arcmin/min]
        double meanMot = angSpeed*24/360;  //[rev/day]
        System.out.println("meanMot " + meanMot + " rev/day");
        double sma = Kepler.getSMAxis2(Constants.GM_Earth, meanMot*Math.PI*2/24/3600);
        System.out.println("sma " + sma/1000 + " km");
        double R_earth = Constants.R_Earth;
        double h_1 = Math.toRadians(40.35);
        double h_2 = Math.toRadians(40.57);
        double h_3 = Math.toRadians(40.72);
        double h_4 = Math.toRadians(40.78);
        //double sma = 63963481;  //[m]
        double r_1, r_2, r_3, r_4;
        r_1 = R_earth*(Math.cos(Math.toRadians(90) + h_1) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_1)));
        System.out.println("r_1 " + r_1/1000 + " km");
        System.out.println("r_1 " + r_1/R_earth + " R_earth");
        r_2 = R_earth*(Math.cos(Math.toRadians(90) + h_2) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_2)));
        System.out.println("r_2 " + r_2/1000 + " km");
        System.out.println("r_2 " + r_2/R_earth + " R_earth");
        r_3 = R_earth*(Math.cos(Math.toRadians(90) + h_3) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_3)));
        System.out.println("r_3 " + r_3/1000 + " km");
        System.out.println("r_3 " + r_3/R_earth + " R_earth");
        r_4 = R_earth*(Math.cos(Math.toRadians(90) + h_4) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_4)));
        System.out.println("r_4 " + r_4/1000 + " km");
        System.out.println("r_4 " + r_4/R_earth + " R_earth");

        Vector positionVec[] = new Vector[2];

        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }
        //1. position
        observation2[0].timeMjd = Time.getMjd(2010,04,22,19,55,57.38);
        observation2[0].ra = Math.toRadians(217.13628);
        observation2[0].dec = Math.toRadians(19.49119);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010,04,22,19,56,34.70);
        observation2[1].ra = Math.toRadians(217.21872);
        observation2[1].dec = Math.toRadians(19.51394);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010,04,22,19,58,04.70);
        observation2[2].ra = Math.toRadians(217.42942);
        observation2[2].dec = Math.toRadians(19.57269);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010,04,22,19,58,42.13);
        observation2[3].ra = Math.toRadians(217.51224);
        observation2[3].dec = Math.toRadians(19.59497);

        double ra[] = new double[2];
        double dec[] = new double[2];
        double time[] = new double[2];
        Geodetic geodetic[] = new Geodetic[2];

        //for(int i = 0; i< 2; i++){
            ra[0] = observation2[0].ra;
            dec[0] = observation2[0].dec;
            time[0] = Time.getJdFromMjd(observation2[0].timeMjd)*1440;
            geodetic[0] = new Geodetic(observation2[0].lon, observation2[0].lat, observation2[0].alt);
            ra[1] = observation2[2].ra;
            dec[1] = observation2[2].dec;
            time[1] = Time.getJdFromMjd(observation2[2].timeMjd)*1440;
            geodetic[1] = new Geodetic(observation2[2].lon, observation2[2].lat, observation2[2].alt);
        //}

        //constants
        double k = 0.07436574;  //[(e.r)^3/2 / min]
        double mi = 1.0;        //[e.m.] - Earth mass
        double a_e = 1.0;       //[e.r.] - Earth radius

        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening = Constants.f_Earth;
        
        positionVec = new RIterationAnglesOnly().getPositionVector(ra, dec, time, geodetic, dTheta_dTime,
                        flattening, a_e, mi, k, false, r_1, r_3);
        Kepler kepler = new Kepler().getElements(Constants.GM_Earth, observation2[0].timeMjd, observation2[1].timeMjd,
                positionVec[0], positionVec[1]);

        System.out.println("Computed a: "+kepler.a/1000+" km");
        System.out.println("Computed e: "+kepler.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" °");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" °");
        System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" °");
        System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" °\n");
        */
        
        /**
         * TEST 11007 / MOLNIYA 1-42 / 1978-080A
         *
           1 11007U 78080A   10112.89969172  .00000021  00000-0  10000-3 0   457
           2 11007 064.4211 047.5757 6891979 250.8522 027.3119 02.31251431246723
         */
        /*
        //assumption that observer - body distance is long, d > 20000 km
        double angSpeed = 85.431; //[arcmin/min]
        double meanMot = angSpeed*24/360;  //[rev/day]
        System.out.println("meanMot " + meanMot + " rev/day");
        double sma = Kepler.getSMAxis2(Constants.GM_Earth, meanMot*Math.PI*2/24/3600);
        System.out.println("sma " + sma/1000 + " km");
        double R_earth = Constants.R_Earth;
        double h_1 = Math.toRadians(44.77198);
        double h_2 = Math.toRadians(44.12071);
        double h_3 = Math.toRadians(42.48322);
        double h_4 = Math.toRadians(41.77346);
        //double sma = 63963481;  //[m]
        double r[] = new double[4];
        r[0] = R_earth*(Math.cos(Math.toRadians(90) + h_1) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_1)));
        System.out.println("r_1 " + r[0]/1000 + " km");
        System.out.println("r_1 " + r[0]/R_earth + " R_earth");
        r[1] = R_earth*(Math.cos(Math.toRadians(90) + h_2) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_2)));
        System.out.println("r_2 " + r[1]/1000 + " km");
        System.out.println("r_2 " + r[1]/R_earth + " R_earth");
        r[2] = R_earth*(Math.cos(Math.toRadians(90) + h_3) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_3)));
        System.out.println("r_3 " + r[2]/1000 + " km");
        System.out.println("r_3 " + r[2]/R_earth + " R_earth");
        r[3] = R_earth*(Math.cos(Math.toRadians(90) + h_4) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_4)));
        System.out.println("r_4 " + r[3]/1000 + " km");
        System.out.println("r_4 " + r[3]/R_earth + " R_earth");

        Vector positionVec[] = new Vector[2];

        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }
        //1. position
        observation2[0].timeMjd = Time.getMjd(2010,04,22,19,55,57.38);
        observation2[0].ra = Math.toRadians(220.17289);
        observation2[0].dec = Math.toRadians(28.34817);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010,04,22,19,56,34.70);
        observation2[1].ra = Math.toRadians(220.59492);
        observation2[1].dec = Math.toRadians(27.63287);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010,04,22,19,58,04.70);
        observation2[2].ra = Math.toRadians(221.62535);
        observation2[2].dec = Math.toRadians(25.84764);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010,04,22,19,58,42.13);
        observation2[3].ra = Math.toRadians(222.05952);
        observation2[3].dec = Math.toRadians(25.0793);

        double ra[] = new double[2];
        double dec[] = new double[2];
        double time[] = new double[2];
        Geodetic geodetic[] = new Geodetic[2];

        int posNo_1 = 0;
        int posNo_2 = 3;
        //for(int i = 0; i< 2; i++){
            ra[0] = observation2[posNo_1].ra;
            dec[0] = observation2[posNo_1].dec;
            time[0] = Time.getJdFromMjd(observation2[posNo_1].timeMjd)*1440;
            geodetic[0] = new Geodetic(observation2[posNo_1].lon, observation2[posNo_1].lat, observation2[posNo_1].alt);
            ra[1] = observation2[posNo_2].ra;
            dec[1] = observation2[posNo_2].dec;
            time[1] = Time.getJdFromMjd(observation2[posNo_2].timeMjd)*1440;
            geodetic[1] = new Geodetic(observation2[posNo_2].lon, observation2[posNo_2].lat, observation2[posNo_2].alt);
        //}

        //constants
        double k = 0.07436574;  //[(e.r)^3/2 / min]
        double mi = 1.0;        //[e.m.] - Earth mass
        double a_e = 1.0;       //[e.r.] - Earth radius

        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening = Constants.f_Earth;

        positionVec = new RIterationAnglesOnly().getPositionVector(ra, dec, time, geodetic, dTheta_dTime,
                        flattening, a_e, mi, k, false, r[posNo_1], r[posNo_2]);
        Kepler kepler = new Kepler().getElements(Constants.GM_Earth, observation2[0].timeMjd, observation2[1].timeMjd,
                positionVec[0], positionVec[1]);

        System.out.println("Computed a: "+kepler.a/1000+" km");
        System.out.println("Computed e: "+kepler.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" �");
        System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" �\n");
        */

       /**
         * TEST 27542 / SL-12 R/B(2)
         *
           1 27542U 02048C   10111.41666667 -.00000438 +00000-0 +10000-3 0 00171
           2 27542 089.5772 349.6806 8915048 249.1735 359.3654 00.36914236004948
         */
         /*
        //assumption that observer - body distance is long, d > 20000 km
        double angSpeed = 1.212; //[arcmin/min]
        double meanMot = angSpeed*24/360;  //[rev/day]
        System.out.println("meanMot " + meanMot + " rev/day");
        double sma = Kepler.getSMAxis2(Constants.GM_Earth, meanMot*Math.PI*2/24/3600);
        System.out.println("sma " + sma/1000 + " km");
        double R_earth = Constants.R_Earth;
        double h_1 = Math.toRadians(38.23611);
        double h_2 = Math.toRadians(38.3879);
        double h_3 = Math.toRadians(38.53801);
        double h_4 = Math.toRadians(38.62253);
        //double sma = 63963481;  //[m]
        double r[] = new double[4];
        r[0] = R_earth*(Math.cos(Math.toRadians(90) + h_1) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_1)));
        System.out.println("r_1 " + r[0]/1000 + " km");
        System.out.println("r_1 " + r[0]/R_earth + " R_earth");
        r[1] = R_earth*(Math.cos(Math.toRadians(90) + h_2) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_2)));
        System.out.println("r_2 " + r[1]/1000 + " km");
        System.out.println("r_2 " + r[1]/R_earth + " R_earth");
        r[2] = R_earth*(Math.cos(Math.toRadians(90) + h_3) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_3)));
        System.out.println("r_3 " + r[2]/1000 + " km");
        System.out.println("r_3 " + r[2]/R_earth + " R_earth");
        r[3] = R_earth*(Math.cos(Math.toRadians(90) + h_4) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_4)));
        System.out.println("r_4 " + r[3]/1000 + " km");
        System.out.println("r_4 " + r[3]/R_earth + " R_earth");

        Vector positionVec[] = new Vector[2];

        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }
        //1. position
        observation2[0].timeMjd = Time.getMjd(2010,04,25,1,50,59);
        observation2[0].ra = Math.toRadians(353.95469);
        observation2[0].dec = Math.toRadians(60.84025);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010,04,25,1,52,20);
        observation2[1].ra = Math.toRadians(353.9582);
        observation2[1].dec = Math.toRadians(60.86746);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010,04,25,1,53,40);
        observation2[2].ra = Math.toRadians(353.96155);
        observation2[2].dec = Math.toRadians(60.89432);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010,04,25,1,54,25);
        observation2[3].ra = Math.toRadians(353.96339);
        observation2[3].dec = Math.toRadians(60.90942);

        double ra[] = new double[2];
        double dec[] = new double[2];
        double time[] = new double[2];
        Geodetic geodetic[] = new Geodetic[2];

        int posNo_1 = 1;
        int posNo_2 = 2;
        //for(int i = 0; i< 2; i++){
            ra[0] = observation2[posNo_1].ra;
            dec[0] = observation2[posNo_1].dec;
            time[0] = Time.getJdFromMjd(observation2[posNo_1].timeMjd)*1440;
            geodetic[0] = new Geodetic(observation2[posNo_1].lon, observation2[posNo_1].lat, observation2[posNo_1].alt);
            ra[1] = observation2[posNo_2].ra;
            dec[1] = observation2[posNo_2].dec;
            time[1] = Time.getJdFromMjd(observation2[posNo_2].timeMjd)*1440;
            geodetic[1] = new Geodetic(observation2[posNo_2].lon, observation2[posNo_2].lat, observation2[posNo_2].alt);
        //}

        //constants
        double k = 0.07436574;  //[(e.r)^3/2 / min]
        double mi = 1.0;        //[e.m.] - Earth mass
        double a_e = 1.0;       //[e.r.] - Earth radius

        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening = Constants.f_Earth;

        positionVec = new RIterationAnglesOnly().getPositionVector(ra, dec, time, geodetic, dTheta_dTime,
                        flattening, a_e, mi, k, false, r[posNo_1], r[posNo_2]);
        Kepler kepler = new Kepler().getElements(Constants.GM_Earth, observation2[0].timeMjd, observation2[1].timeMjd,
                positionVec[0], positionVec[1]);

        System.out.println("Computed a: "+kepler.a/1000+" km");
        System.out.println("Computed e: "+kepler.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" �");
        System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" �\n");
        */
        /**
         * TEST 23723 / AMOS 5I (ASIASAT 2)
         *
           1 23723U 95064A   10111.75172963  .00000099  00000-0  10000-3 0  2125
           2 23723 000.0398 052.9447 0002405 338.9741 105.3766 01.00271407 52733
         */
         /*
        //assumption that observer - body distance is long, d > 20000 km
        double angSpeed = 14.919; //[arcmin/min]
        double meanMot = angSpeed*24/360;  //[rev/day]
        System.out.println("meanMot " + meanMot + " rev/day");
        double sma = Kepler.getSMAxis2(Constants.GM_Earth, meanMot*Math.PI*2/24/3600);
        System.out.println("sma " + sma/1000 + " km");
        double R_earth = Constants.R_Earth;
        double h_1 = Math.toRadians(34.46897);
        double h_2 = Math.toRadians(34.46894);
        double h_3 = Math.toRadians(34.46891);
        double h_4 = Math.toRadians(34.46889);
        //double sma = 63963481;  //[m]
        double r[] = new double[4];
        r[0] = R_earth*(Math.cos(Math.toRadians(90) + h_1) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_1)));
        System.out.println("r_1 " + r[0]/1000 + " km");
        System.out.println("r_1 " + r[0]/R_earth + " R_earth");
        r[1] = R_earth*(Math.cos(Math.toRadians(90) + h_2) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_2)));
        System.out.println("r_2 " + r[1]/1000 + " km");
        System.out.println("r_2 " + r[1]/R_earth + " R_earth");
        r[2] = R_earth*(Math.cos(Math.toRadians(90) + h_3) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_3)));
        System.out.println("r_3 " + r[2]/1000 + " km");
        System.out.println("r_3 " + r[2]/R_earth + " R_earth");
        r[3] = R_earth*(Math.cos(Math.toRadians(90) + h_4) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_4)));
        System.out.println("r_4 " + r[3]/1000 + " km");
        System.out.println("r_4 " + r[3]/R_earth + " R_earth");

        Vector positionVec[] = new Vector[2];

        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }
        //1. position
        observation2[0].timeMjd = Time.getMjd(2010,04,25,1,50,59);
        observation2[0].ra = Math.toRadians(257.50652);
        observation2[0].dec = Math.toRadians(-7.14458);

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010,04,25,1,52,20);
        observation2[1].ra = Math.toRadians(257.84481);
        observation2[1].dec = Math.toRadians(-7.14495);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010,04,25,1,53,40);
        observation2[2].ra = Math.toRadians(258.17892);
        observation2[2].dec = Math.toRadians(-7.1453);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010,04,25,1,54,25);
        observation2[3].ra = Math.toRadians(258.36686);
        observation2[3].dec = Math.toRadians(-7.1455);

        double ra[] = new double[2];
        double dec[] = new double[2];
        double time[] = new double[2];
        Geodetic geodetic[] = new Geodetic[2];

        int posNo_1 = 0;
        int posNo_2 = 3;
        //for(int i = 0; i< 2; i++){
            ra[0] = observation2[posNo_1].ra;
            dec[0] = observation2[posNo_1].dec;
            time[0] = Time.getJdFromMjd(observation2[posNo_1].timeMjd)*1440;
            geodetic[0] = new Geodetic(observation2[posNo_1].lon, observation2[posNo_1].lat, observation2[posNo_1].alt);
            ra[1] = observation2[posNo_2].ra;
            dec[1] = observation2[posNo_2].dec;
            time[1] = Time.getJdFromMjd(observation2[posNo_2].timeMjd)*1440;
            geodetic[1] = new Geodetic(observation2[posNo_2].lon, observation2[posNo_2].lat, observation2[posNo_2].alt);
        //}

        //constants
        double k = 0.07436574;  //[(e.r)^3/2 / min]
        double mi = 1.0;        //[e.m.] - Earth mass
        double a_e = 1.0;       //[e.r.] - Earth radius

        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening = Constants.f_Earth;

        positionVec = new RIterationAnglesOnly().getPositionVector(ra, dec, time, geodetic, dTheta_dTime,
                        flattening, a_e, mi, k, false, r[posNo_1], r[posNo_2]);
        Kepler kepler = new Kepler().getElements(Constants.GM_Earth, observation2[0].timeMjd, observation2[1].timeMjd,
                positionVec[0], positionVec[1]);

        System.out.println("Computed a: "+kepler.a/1000+" km");
        System.out.println("Computed e: "+kepler.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" �");
        System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" �\n");
        */

        /**
         * Test Asteroid 2010 AL30
         *  e 	.3070150178127762 	4.4089e-07
            a 	1.045143386619743 	6.9658e-07 	AU
            q 	.7242686711597774 	7.382e-07 	AU
            i 	3.812809187941254 	9.1348e-06 	deg
            node 	113.0042349280855 	8.9387e-07 	deg
            peri 	97.06208889462762 	0.00013755 	deg
            M 	114.0197780691696 	6.5698e-05 	deg
            tp 	2455276.893891800493
            (2010-Mar-21.39389181) 	5.4e-05 	JED
            period 	390.2673703226114
            1.07 	0.00039016
            1.068e-06 	d
            yr
            n 	.9224445274592362 	9.222e-07 	deg/d
        Q 	1.366018102079709 	9.1044e-07 	AU
         */

        //assumption that observer - body distance is long, d > 20000 km
        double angSpeed = 16.3; //[arcmin/min]
        double meanMot = angSpeed*24/360;  //[rev/day]
        System.out.println("meanMot " + meanMot + " rev/day");
        double sma = Kepler.getSMAxis2(Constants.GM_Earth, meanMot*Math.PI*2/24/3600);
        System.out.println("sma " + sma/1000 + " km");
        double R_earth = Constants.R_Earth;
        double h_1 = Math.toRadians(32.62);
        double h_2 = Math.toRadians(32.75);
        double h_3 = Math.toRadians(32.87);
        double h_4 = Math.toRadians(33);
        //double sma = 63963481;  //[m]
        double r[] = new double[4];
        r[0] = R_earth*(Math.cos(Math.toRadians(90) + h_1) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_1)));
        System.out.println("r_1 " + r[0]/1000 + " km");
        System.out.println("r_1 " + r[0]/R_earth + " R_earth");
        r[1] = R_earth*(Math.cos(Math.toRadians(90) + h_2) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_2)));
        System.out.println("r_2 " + r[1]/1000 + " km");
        System.out.println("r_2 " + r[1]/R_earth + " R_earth");
        r[2] = R_earth*(Math.cos(Math.toRadians(90) + h_3) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_3)));
        System.out.println("r_3 " + r[2]/1000 + " km");
        System.out.println("r_3 " + r[2]/R_earth + " R_earth");
        r[3] = R_earth*(Math.cos(Math.toRadians(90) + h_4) + Math.sqrt(sma*sma/(R_earth*R_earth) - Math.sin(Math.toRadians(90) + h_4)));
        System.out.println("r_4 " + r[3]/1000 + " km");
        System.out.println("r_4 " + r[3]/R_earth + " R_earth");

        Vector positionVec[] = new Vector[2];

        Observation observation2[] = new Observation[4];
        //AGO Modra
        for(int i = 0; i<observation2.length; i++){
            observation2[i] = new Observation();
            observation2[i].lon = Math.toRadians(17.2740);
            observation2[i].lat = Math.toRadians(48.3733);
            observation2[i].alt = 531.1;
        }
       //1. position
        observation2[0].timeMjd = Time.getMjd(2010, 1, 13, 13, 20, 0.0);
        observation2[0].ra = Math.toRadians(20.44937);
        observation2[0].dec = Math.toRadians(9.28192 );

        //2. postion
        observation2[1].timeMjd = Time.getMjd(2010, 1, 13, 13, 21, 0.0);
        observation2[1].ra = Math.toRadians(20.17801);
        observation2[1].dec = Math.toRadians(9.23617);

        //3. postion
        observation2[2].timeMjd = Time.getMjd(2010, 1, 13, 13, 22, 0.0);
        observation2[2].ra = Math.toRadians(19.90698);
        observation2[2].dec = Math.toRadians(9.19026);

        //4. postion
        observation2[3].timeMjd = Time.getMjd(2010, 1, 13, 13, 23, 0.0);
        observation2[3].ra = Math.toRadians(19.63629);
        observation2[3].dec = Math.toRadians(9.14420);

        double ra[] = new double[2];
        double dec[] = new double[2];
        double time[] = new double[2];
        Geodetic geodetic[] = new Geodetic[2];

        int posNo_1 = 0;
        int posNo_2 = 2;
        //for(int i = 0; i< 2; i++){
            ra[0] = observation2[posNo_1].ra;
            dec[0] = observation2[posNo_1].dec;
            time[0] = Time.getJdFromMjd(observation2[posNo_1].timeMjd)*1440;
            geodetic[0] = new Geodetic(observation2[posNo_1].lon, observation2[posNo_1].lat, observation2[posNo_1].alt);
            ra[1] = observation2[posNo_2].ra;
            dec[1] = observation2[posNo_2].dec;
            time[1] = Time.getJdFromMjd(observation2[posNo_2].timeMjd)*1440;
            geodetic[1] = new Geodetic(observation2[posNo_2].lon, observation2[posNo_2].lat, observation2[posNo_2].alt);
        //}

        //constants
        double k = 0.07436574;  //[(e.r)^3/2 / min]
        double mi = 1.0;        //[e.m.] - Earth mass
        double a_e = 1.0;       //[e.r.] - Earth radius

        //double dTheta_dTime2 = 1 + 1.0/365.24219879; //[rev/day]
        double dTheta_dTime = 4.3752695e-3; //[rad/min]
        //double dTheta_dTime2 = Math.toRadians(0.25068447); //[rad/min]
        double flattening = Constants.f_Earth;

        positionVec = new RIterationAnglesOnly().getPositionVector(ra, dec, time, geodetic, dTheta_dTime,
                        flattening, a_e, mi, k, false, r[posNo_1], r[posNo_2]);
        Kepler kepler = new Kepler().getElements(Constants.GM_Earth, observation2[0].timeMjd, observation2[1].timeMjd,
                positionVec[0], positionVec[1]);

        System.out.println("Computed a: "+kepler.a/1000+" km");
        System.out.println("Computed e: "+kepler.e);
        System.out.println("Computed i: "+Math.toDegrees(kepler.incl)+" �");
        System.out.println("Computed Omega: "+Math.toDegrees(kepler.Omega)+" �");
        System.out.println("Computed omega: "+Math.toDegrees(kepler.omega)+" �");
        System.out.println("Ma        : "+Math.toDegrees(kepler.M)+" �\n");

    }

}
