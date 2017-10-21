//------------------------------------------------------------------------------
//
//Vector.java
// 
//Ucel (purpose) :
//
//	Definicia vektora a metody na pracovanie s vektormi (vector definition and methods for vectors)
//
//
// 2007/10/10 - Jiri Silha
//
//------------------------------------------------------------------------------

package Compute.Silha;

//
// Deklaracia triedy (class declaration)
//

/** Mathematical operations with vectors*/
public class Vector{
	
	//
	//Deklaracia premennych pre konstruktor (variables declaration)
	//
	
	//double x,y,z;
	public double v[];
	
	//
	//Konstruktor (constructor)
	//
	
	public Vector(int x){
		v=new double[x]; 
	}
	
	/**
	*   Vypocet velkosti vektora (Computing the magnitude of vector)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector vector
	*
	*   <navratna hodnota> (<return>) Velkost vektora (size of vector)
	*
	*/
	
	public static double getSize (Vector vector){
		double size=Math.sqrt(vector.v[0]*vector.v[0]+vector.v[1]*vector.v[1]+vector.v[2]*vector.v[2]);
		return size;
	}

	/**
	*   Vypocet skalarneho sucinu 2 vektorov (Computing the scalar product of two vectors)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector vector1, Vector vector2 
	*
	*   <navratna hodnota> (<return>) Sucin (scalar product)
	*
	*/
	
	public static double getScalarProduct(Vector vector1, Vector vector2){
		double product = vector1.v[0]*vector2.v[0] + vector1.v[1]*vector2.v[1] + vector1.v[2]*vector2.v[2];
		return product;
	}
        public static double getScalarProduct_2(Vector vector1, Vector vector2){
		double product = Vector.getSize(vector1)*Vector.getSize(vector2)*Math.cos(Vector.getVectorsAngle(vector1, vector2));
		return product;
	}

	
	/**
	*
	* Vypocet vectoroveho sucinu 2 vektorov (Computing the vector product of two vectors)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector vector1, Vector vector2 
	*
	*   <navratna hodnota> (<return>) Vektorovy sucin (vector product)
	*
	*/
	
	public static Vector getVectorProduct(Vector vector1, Vector vector2 ){
		Vector vectorProduct=new Vector(3);
		vectorProduct.v[0]=vector1.v[1]*vector2.v[2]-vector1.v[2]*vector2.v[1];
		vectorProduct.v[1]=vector1.v[2]*vector2.v[0]-vector1.v[0]*vector2.v[2];
		vectorProduct.v[2]=vector1.v[0]*vector2.v[1]-vector1.v[1]*vector2.v[0];
		return vectorProduct;
	}
	
	/**
	*
	* Sucin vektora s konstantov (devide and multiply vector by scalar)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector  
	*  Skalar (scalar)
	*
	*   <navratna hodnota> (<return>) Vektor (vector)
	*
	*/
	
	public static Vector multiplyVector(Vector vector, double x ){
		Vector vector2=new Vector(3);
		vector2.v[0]=vector.v[0]*x;
		vector2.v[1]=vector.v[1]*x;
		vector2.v[2]=vector.v[2]*x;
		return vector2;
	}
	
	/**
	*
	* Scitanie 2 vektorov (add vector to vector)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector a 
	*  Vector b
	*
	*   <navratna hodnota> (<return>) Vektor (vector)
	*
	*/
	
	public static Vector addVectors(Vector vector1, Vector vector2){
		Vector vector3=new Vector(3);
		vector3.v[0]=vector1.v[0]+vector2.v[0];
		vector3.v[1]=vector1.v[1]+vector2.v[1];
		vector3.v[2]=vector1.v[2]+vector2.v[2];
		return vector3;
	}
	
	/**
	*  Odcitanie 2 vektorov (subtract 2 vectors)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector a 
	*  Vector b
	*
	*   <navratna hodnota> (<return>) Vektor (vector)
	*
	*/
	
	public static Vector subtractVectors(Vector vector1, Vector vector2){
		Vector vector3=new Vector(3);
		vector3.v[0] = vector1.v[0] - vector2.v[0];
		vector3.v[1] = vector1.v[1] - vector2.v[1];
		vector3.v[2] = vector1.v[2] - vector2.v[2];
		return vector3;
	}
	
	/**
	*  Prevod z polarnych do kartezskych suradnic (transformation from polar to cartesian coordinate system)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Polar Vector (azim [rad], elev[rad], r[m])
	*
	*   <navratna hodnota> (<return>) Vektor (vector)
	*
	*/
	
	public static Vector fromPolarToCartesian(double azim, double elev, double r){
		Vector vector=new Vector(3);
		vector.v[0]=r*Math.cos(azim)*Math.cos(elev);
		vector.v[1]=r*Math.sin(azim)*Math.cos(elev);
		vector.v[2]=r*Math.sin(elev);
		return vector;
	}
	
	/**
	*  Prevod do sferickych z kartezskych suradnic (transformation from cartesian to polar coordinate system)
	*
	* Vstup/Vystup (Input/output):
	*
	*  Vector vector
	*
	*   <navratna hodnota> (<return>) Polar Vector (azim [rad], elev[rad], r[m])
	*
	*/
	
	public static Vector fromCartesianToPolar(Vector r){
		Vector vector=new Vector(3);
		//azimut
		vector.v[0]=Math.atan2(r.v[1],r.v[0]);
		while((vector.v[0]<0)) vector.v[0]=vector.v[0]+Math.PI*2;
		while((vector.v[0]>Math.PI*2)) vector.v[0]=vector.v[0]-Math.PI*2;
		//vyska
		vector.v[1]=Math.atan2(r.v[2],(Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1])));
		//vzdialenost
		vector.v[2]=Math.sqrt(r.v[0]*r.v[0]+r.v[1]*r.v[1]+r.v[2]*r.v[2]);
		return vector;
	}
        
       /**
	*   Vypocet uhlu medzi dvoma vektormi (Computing angle between 2 vectors)
	* Vstup/Vystup (Input/output):
	*
	*  Vector vector1, Vector vector2
	*
	*   <navratna hodnota> (<return>) Uhol (angle between) [rad]
	*
	*/
	
        public static double getVectorsAngle(Vector vector1, Vector vector2){
		//System.out.println("!!!!!!!!!!!!!!!!BEWARE!!!!!!!!!!!!!!!!!!!");
                double a = getScalarProduct(vector1,vector2);
		double b = Math.abs(getSize(vector1)*getSize(vector2));
		double angle = Math.acos(a/b);
		return angle;
	}
	
	//
	//Test method
	//
	public static void main(String args[]){
		Vector vector1=new Vector(3);
		Vector vector2=new Vector(3);
		vector1.v[0]=-5;
		vector1.v[1]=5;
		vector1.v[2]=5;
		
		vector2.v[0]=10;
		vector2.v[1]=5;
		vector2.v[2]=15;

                double dotProduct = Vector.getScalarProduct(vector1, vector2);
                System.out.println("1st " + dotProduct);
                dotProduct = Vector.getScalarProduct_2(vector1, vector2);
                System.out.println("2nd " + dotProduct);
                double angle = Vector.getVectorsAngle(vector1, vector2);
                System.out.println("1st [rad]" + angle);
                System.out.println("2nd [deg]" + Math.toDegrees(angle));
		//Vector vector3=devideVector(vector1,5);
		//Vector vector3=subtructVectors(vector1,vector2);
		//System.out.println(vector.x+" "+vector.y+" "+vector.z);
		//System.out.println(vector3.v[0]+" "+vector3.v[1]+" "+vector3.v[2]);
		//System.out.println(35.1%7);
	}

        /**
      * getPositionVectorsDistance()
      *
      * Method to compute distance between 2 objects with given position vectors.
      *
      * INPUT:
      *     Vector posVec_1 - position vector of 1st object [vector units]
      *     Vector posVec_2 - position vector of 2nd object [vector units]
      *
      * OUTPUT:
      *     double distance - distance between two points [vector units]
      */
     public double getPositionVectorsDistance(Vector posVec_1, Vector posVec_2){
         double distance = 0;
         distance = Math.sqrt(Math.pow(posVec_1.v[0] - posVec_2.v[0],2) +
                              Math.pow(posVec_1.v[1] - posVec_2.v[1],2) +
                              Math.pow(posVec_1.v[2] - posVec_2.v[2],2));
         return distance;
     }
	
}