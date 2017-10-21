/*
 * Method to compute roots of cubic formula.
 */

package Compute.Silha;

// July 2007 - Jiri Silha

/*
*Trieda sluzi na vypocet korena kubicke (polynom 3.stupna) polynomu, vyuziva sa pri tom Cardenova metoda vypoctu kub. polynomu
*Vstup: a,b,c,d	-" a*x^3 + b*x^2 + c*x + d = 0 "
*Vystup: x
*/
public class CubicFormula{
	//premenne v zakladnej rovnici
	public double a,b,c,d;
	//vystupna hodnota
        double x;
	//pomocne premenne pri subtituciach
	double A,B;
		
	/**
         * getRoot()
         * 
         * Method to compute root of cubic formula
         * 
         * IN:
         * @param a
         * @param b
         * @param c
         * @param d
         * @return double root
         */
	public double getRoot(double a, double b, double c, double d){
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
            
            //System.out.println("a " + a);
            //System.out.println("b " + b);
            //System.out.println("c " + c);
            //System.out.println("d " + d);
		//premenne A a B si vyjadrime vdaka substitucii z Cardonovej metody
		A=c/a-b*b/(3*a*a);
		B=d/a+2*Math.pow((b/(3*a)),3)-b*c/(3*a*a);
		
		//pre koren x plati rovnica
		//pomocne premenne, koli zapornemu clenu pod odmocninou
		double K,L;
		
		K=-B/2+Math.sqrt((B*B/4+Math.pow(A/3,3)));
		//ak je K zaporne treba ju vynasobit -1, odmocnin a nasobit spat -1
		if(K<0) K=-1* Math.pow(-1*K,(double)1/3);
		else 	K=Math.pow(K,(double)1/3);
		
		//ak je L zaporne treba ju vynasobit -1, odmocnin a nasobit spat -1
		L=B/2+Math.sqrt((B*B/4+Math.pow(A/3,3)));
		if(L<0) L=-1*Math.pow(-1*L,(double)1/3);
		else 	L=Math.pow(L,(double)1/3);
		
		//vypocet korena
		x=K - L - b/(3*a);
		
		return x;
	}
	
	
	//test metody
	public static void main (String[] args){
		System.out.println(new CubicFormula().getRoot(1,-15,81,-175));
	}
}
