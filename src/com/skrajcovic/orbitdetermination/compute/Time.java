//------------------------------------------------------------------------------
//
//Time.java
// 
//Ucel (purpose) :
//
//	Vypocet casu a datumu (Time and date computation)
//
//Poznamka (notes):
//
//	Trieda je modifikacia triedy SAT_Time.h od O. Montenbruck, E. Gill (2005/04/14  OMO  Final version (2nd reprint))
//
// 2007/10/09 - Jiri Silha
//
//------------------------------------------------------------------------------

package com.skrajcovic.orbitdetermination.compute;

import java.util.GregorianCalendar;


//
// Deklaracia triedy (class declaration)
//

/** Time operations*/
public class Time{
	
	//
	//Deklaracia premennych pre konstruktor (variables declaration)
	//
	public double mjd;
	public int year,month,day, hour, min;
	public double sec;
	
	public Time(){}
	
	public Time(double mjd, int year, int month, int day, int hour, int min, double sec){
		if(sec >= 60){
            sec = sec - 60;
            min++;
        }
        if(min >= 60){
            min = min - 60;
            hour++;
        }
        if(hour >= 24){
            hour = hour -24;
            day++;
        }

        if(month > 12) {
            month = month - 12;
            year++;
        }

        this.mjd=mjd;
		this.year=year;
		this.month=month;
		this.day=day;
		this.hour=hour;
		this.min=min;
		this.sec=sec;
	}
	
	public Time(int year, int month, int day, int hour, int min, double sec){
		if(sec >= 60){
            sec = sec - 60;
            min++;
        }
        if(min >= 60){
            min = min - 60;
            hour++;
        }
        if(hour >= 24){
            hour = hour -24;
            day++;
        }

        if(month > 12) {
            month = month - 12;
            year++;
        }

        this.year=year;
		this.month=month;
		this.day=day;
		this.hour=hour;
		this.min=min;
		this.sec=sec;
	}
	
	/*
	*   Modifikovany juliansky datum z kalendarneho datumu a casu (Modified Julian Date from calendar date and time)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Rok(Year)      Calendar date components
	*   Mesiac(Month)
	*   Den(Day)
	*   Hodina(Hour)      Time components (optional)
	*   Min
	*   Sec
	*
	*   <navratna hodnota> (<return>) Modifikovany juliansky datum (Modified Julian Date) [days]
	*
	*/

	public static double getMjd ( int Year, int Month, int Day, int Hour, int Min, double Sec ){
		// Premenne (Variables)
		long    MjdMidnight;
		double  FracOfDay;
		int     b;


		if (Month<=2) { Month+=12; --Year;}
		
		if ( (10000L*Year+100L*Month+Day) <= 15821004L ){
			b = -2 + ((Year+4716)/4) - 1179;     // Julian calendar 
		}
		else{
			b = (Year/400)-(Year/100)+(Year/4);  // Gregorian calendar 
		}
		
		MjdMidnight = 365L*Year - 679004L + b + (int)(30.6001*(Month+1)) + Day;
		FracOfDay   = (Hour+Min/60.0+Sec/3600.0) / 24.0; 

		return MjdMidnight + FracOfDay;
	}
	
	/**
         * Modifikovany juliansky datum z kalendarneho datumu a casu (Modified Julian Date from calendar date and time)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Rok(Year)      Calendar date components
	*   Mesiac(Month)
	*   Den(Day)
	*   Hodina(Hour)      Time components (optional)
	*   Min
	*   Sec
	*
	*   <navratna hodnota> (<return>) Modifikovany juliansky datum (Modified Julian Date)
	*
	*/

	public static double getMjd ( Time time ){
		// Premenne (Variables)
		double  MjdMidnight;
		double  FracOfDay;
		int     b;
		
		int Year = time.year;
		int Month = time.month;
		int Day = time.day;
		int Hour = time.hour;
		int Min = time.min;
		double Sec = time.sec;
		
		if (Month<=2) { Month+=12; --Year;}
		
		if ( (10000L*Year+100L*Month+Day) <= 15821004L ){
			b = -2 + ((Year+4716)/4) - 1179;     // Julian calendar 
		}
		else{
			b = (Year/400)-(Year/100)+(Year/4);  // Gregorian calendar 
		}
		
		MjdMidnight = 365L*Year - 679004L + b + (int)(30.6001*(Month+1)) + Day;
		FracOfDay   = (Hour+Min/60.0+Sec/3600.0) / 24.0; 

		return MjdMidnight + FracOfDay;
	}

        /**
         * getJdFromMjd()
         *
         * Method to compute Julian date from Modified Julian Date
         *
         * INPUT:
         * double mjd - modified julian date in [day]
         *
         * OUTPUT:
         * double jd - julian date in [day]
         */
        public static double getJdFromMjd(double mjd){
            double jd = mjd + 2400000.5;
            return jd; 
        }

        /**
         * getMjdFromJd()
         *
         * Method to compute Modified Julian date from Julian Date
         *
         * INPUT:
         * double jd - julian date in [day]
         *
         * OUTPUT:
         * double mjd - modified julian date in [day]
         */
        public static double getMjdFromJd(double jd){
            double mjd = jd - 2400000.5;
            return mjd;
        }
	
	
	/**
	*
	* getDateTime
	*
	* Ucel(purpose):
	*
	*   Kalendarny datum a cas z modifikovaneho julianskeho datumu (calendar date and time from Modified Julian Date)
	*
	* Vstup/vystup (Input/output):
	*
	*   Mjd       Modified Julian Date
	*   
	*   <navratna hodnota> (<return>) Objekt typu Time()
	*   Rok(Year)      Calendar date components
	*   Mesiac(Month)
	*   Den(Day)
	*   Hodina(Hour)      Time components
	*   Min
	*   Sec
	*
	*/

	public static Time getDateTime ( double Mjd){
		// Variables
		long    a,b,c,d,e,f;
		double  Hours,x;
		int Year,Month, Day, Hour, Min;
		double Sec; 
		Time time;
		
		// Konvertovanie pocet julianskych dni do kalensarneho datumu (convert Julian day number to calendar date)
		a = (long)(Mjd+2400001.0);

		if ( a < 2299161 ) {  // Julian calendar
			b = 0;
			c = a + 1524;
		}
		else {                // Gregorian calendar
			b = (long)((a-1867216.25)/36524.25);
			c = a +  b - (b/4) + 1525;
		}

		d     = (long) ( (c-122.1)/365.25 );
		e     = 365*d + d/4;
		f     = (long) ( (c-e)/30.6001 );

		Day   = (int)(c - e - (int)(30.6001*f));
		Month =(int)( f - 1 - 12*(f/14));
		Year  = (int)(d - 4715 - ((7+Month)/10));

		Hours = 24.0*(Mjd-(int)(Mjd));

		Hour = (int)(Hours);
		x = (Hours-Hour)*60.0; 
		Min = (int)(x);  
		Sec = (x-Min)*60.0;
                Sec = (float)(Math.round((Sec*100)))/100;    
		//ukladame ziskane hodnoty do konstruktora
		time=new Time(Mjd,Year,Month, Day,Hour,Min,Sec);
		
		return time;
	}
	
	/**
         *  Ziskanie greenwichskeho stredneho siderickeho casu z jul. modif. datumu 
	*  ( Get greenwich mean siderial time from modiified julian date)
	*
	* Vstup/vystup (Input/output):
	*
	*   Mjd       Modified Julian Date  [year]
	*   
	*   <navratna hodnota> (<return>) GMST [rad]
	*
	*/
	
	public static double getGMST(double mjd_ut1){
		// Constants

		double secs =  86400.0;        // Seconds per day

		// Variables

		double mjd_0,ut1,t_0,t,gmst;
		Constants c = new Constants();
		
		// Mean Sidereal Time
                double yearLength = 36525.0;
		mjd_0 = Math.floor(mjd_ut1);//(int)(mjd_ut1);
		ut1   = secs*(mjd_ut1-mjd_0);          // [s]
		t_0   = (mjd_0   - c.MJD_J2000)/yearLength;
		t     = (mjd_ut1 - c.MJD_J2000)/yearLength;

        //test xsi + ecl angle
        //double xsi = Math.toRadians((-1) * (0.0048*Math.sin(Math.toRadians(125)) - 0.0004*Math.sin(Math.toRadians(201))));
        //double eps = Math.toRadians(23.44);

		gmst  = 24110.54841 + 8640184.812866*t_0 + 1.002737909350795*ut1
			+ (0.093104-6.2e-6*t)*t*t; // [s]

		gmst = c.pi2*((gmst/secs)-Math.floor(gmst/secs));// + xsi*Math.cos(eps);       // [rad], 0..2pi
        
        //test vallado
        /*
        while(gmst > c.pi2) gmst = gmst - c.pi2;
        System.out.println("gmst " + Math.toDegrees(gmst));
		//gmst = c.pi2*(0.799057273264 + 1.00273781191135448*((2400000.5 + mjd_ut1) - 2451545));       // [rad], 0..2pi
		gmst = c.pi2*(0.779057273264 + 1.00273781191135448*(mjd_ut1 - 51544.5));       // [rad], 0..2pi
        while(gmst > c.pi2) gmst = gmst - c.pi2;
        System.out.println("era " + Math.toDegrees(gmst));
        */
        
		return gmst;
	}
	
	/**
	*
	* getDateFromEpoch
	*
	* Ucel(purpose):
	*
	*  
	* Vstup/vystup (Input/output):
	*
	*   epoch
	*   
	*   <navratna hodnota> (<return>) 
	*
	*/
	
	public static Time getDateFromEpoch(double epoch){
		Time time=new Time();
		String epochString=epoch+"";
                try{
                    time.year=(int)Double.parseDouble(epochString.substring(0,1));
                    if(time.year>57) time.year=time.year+1900;
                    else time.year=time.year+2000;
                }
                catch(Exception e){
                    System.out.println("Warning! " + e.getMessage());
                }
		//int days		
		return time;
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
		if(time.sec>9) stringTime=stringTime+(float)(time.sec);
			else stringTime=stringTime+"0"+(float)(time.sec);
		//System.out.println("2 " + stringTime);
		return stringTime;
	}
	
	//
	//Test method
	//
	public static void main(String args[]){
            //System.out.println(getNoOfDaysInMonth(2100,1));
        }

        /**
         * getNoOfDaysInMonth()
         *
         * INPUT:
         * int year2
         * int month2 - from 0 - 11
         *
         * OUTPUT:
         * int noOfDays
         */
        public static int getNoOfDaysInMonth(int year2, int month2){
            GregorianCalendar cc = new GregorianCalendar();
            int daysInMonths[] = {31,28,31,30,31,30,31,31,30,31,30,31};
            if(cc.isLeapYear(year2)){
                daysInMonths[1]++;
            }
            int noOfDays = daysInMonths[month2];
            return noOfDays;
        }
	
}
