/*
*
*Matrix.java
* 
*Ucel (purpose) :
*
*	Definicia matice a metody na pracu s nimi (matrix definition and methods for matrix)
*
*
* 2007/10/11 - Jiri Silha
*
**/

package com.skrajcovic.orbitdetermination.compute;

/*
* Deklaracia triedy (class declaration)
*/
/** Mathematical operations with matrixes*/
public class Matrix{
	public double [][]matrix;
	
	public Matrix(int x, int y){
		matrix=new double[x][y];
	}
	
	/*
	METODY NA PRACOVANIE S MATICAMI
	*/
	
	/**
	*
	*   Otocenie matice okolo x - suradnice v  kladnom smere? (proti pohybu hod. ruciciek) ( rotate matrix around x - axis)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Uhol (angle)	  
	*
	*   <navratna hodnota> (<return>) Otocena matica 3x3 (matrix)
	*
	*/
	
	public static Matrix R_x(double Angle){
		double C = Math.cos(Angle);
		double S = Math.sin(Angle);
		Matrix U=new Matrix(3,3);
		U.matrix[0][0] = 1.0;  U.matrix[0][1] = 0.0;  U.matrix[0][2] = 0.0;
		U.matrix[1][0] = 0.0;  U.matrix[1][1] =  +C;  U.matrix[1][2] =  +S;
		U.matrix[2][0] = 0.0;  U.matrix[2][1] =  -S;  U.matrix[2][2] =  +C;
		return U;
	}

	/**
	*   Otocenie matice okolo y - suradnice v  kladnom smere? (proti pohybu hod. ruciciek) ( rotate matrix around y - axis)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Uhol (angle)	  
	*
	*   <navratna hodnota> (<return>) Otocena matica 3x3 (matrix)
	*
	*/
	
	public static Matrix R_y(double Angle){
		double C = Math.cos(Angle);
		double S = Math.sin(Angle);
		Matrix U=new Matrix(3,3);
		U.matrix[0][0] =  +C;  U.matrix[0][1] = 0.0;  U.matrix[0][2] =  -S;
		U.matrix[1][0] = 0.0;  U.matrix[1][1] = 1.0;  U.matrix[1][2] = 0.0;
		U.matrix[2][0] =  +S;  U.matrix[2][1] = 0.0;  U.matrix[2][2] =  +C;
		return U;
	}
	
	/**
	*   Otocenie matice okolo x - suradnice v  kladnom smere? (proti pohybu hod. ruciciek) ( rotate matrix around z - axis)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Uhol (angle)	  
	*
	*   <navratna hodnota> (<return>) Otocena matica 3x3 (matrix)
	*
	*/
	
	public static Matrix R_z(double Angle){
		double C = Math.cos(Angle);
		double S = Math.sin(Angle);
		Matrix U=new Matrix(3,3);
		U.matrix[0][0] =  +C;  U.matrix[0][1] =  +S;  U.matrix[0][2] = 0.0;
		U.matrix[1][0] =  -S;  U.matrix[1][1] =  +C;  U.matrix[1][2] = 0.0;
		U.matrix[2][0] = 0.0;  U.matrix[2][1] = 0.0;  U.matrix[2][2] = 1.0;
		return U;
	}

        /**
	*   Rotation matrix - incoming alpha, betam gamma
         * SOURCE:This program uses the divergence theorem to evaluate the volume of a solid
c  defined by triangular facets.
c
c Dr. Phillip D. Anz-Meador
c Sr. Science Specialist
c ESCG Jacobs
c m/c JE104
c 2224 Bay Area Blvd., Box 7
c Houston, TX  77058
c 1+281-483-4217 (voice)
c 1+281-244-5084 (fax)
c phillip.d.anz-meador@nasa.gov
c
c  Rev     Date                     Description
c -----  ---------   ----------------------------------------------
c BASIC  08/III/10   initial implementation & test vs. cylinder
c
	*
	* Vstup/Vystup (Input/output):
	*
	*   double aplha, beta, gamma   -   angles of rotation around x,y, z axis [rad]
	*
	*   <navratna hodnota> (<return>) Rotation matrix 3x3
	*
	*/

	public static Matrix rotationMatrix(double alpha, double beta, double gamma){
            //rotation matrix
            Matrix rotMat = new Matrix(3,3);
            //compute SINEs and COSINEs of angles [rad]
            double cr = Math.cos(alpha);
            double sr = Math.sin(alpha);
            double cp = Math.cos(beta);
            double sp = Math.sin(beta);
            double cy = Math.cos(gamma);
            double sy = Math.sin(gamma);
            //define rotation matrix T_rpy
           rotMat.matrix[0][0] = cy*cp;
           rotMat.matrix[0][1] = cy*sr*sp + sy*cr;
           rotMat.matrix[0][2] = -sp*cr*cy + sy*sr;
           rotMat.matrix[1][0] = -sy*cp;
           rotMat.matrix[1][1] = -sy*sp*sr + cy*cr;
           rotMat.matrix[1][2] = sy*sp*cr + cy*sr;
           rotMat.matrix[2][0] = sp;
           rotMat.matrix[2][1] = -cp*sr;
           rotMat.matrix[2][2] = cp*cr;
            return rotMat;
	}
	
	/**
	*
	*   Nasobenie 2 matic ( 2 matrixes multiply)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Matrix3x3 1
	*   Matrix3x3 2	
	*
	*   <navratna hodnota> (<return>) Matrix 3x3 
	*
	*/
	
	public Matrix getMatrixProduct(Matrix m1, Matrix m2){
		Matrix aux=new Matrix(3,3);
		double sum;
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				sum = 0.0;
				for (int k=0; k<3; k++) {
					sum += m1.matrix[i][k] * m2.matrix[k][j];
					aux.matrix[i][j] = sum;
				}
			}
		}
		return aux;	
	}
	
	/**
	*
	*  Otocenie matice (vymena riadkov a stlpcov) ( transpose matrix)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Matrix3x3
	*   
	*   <navratna hodnota> (<return>) Matrix 3x3 
	*
	*/
        
	public static Matrix transposeMatrix(Matrix m){
		 Matrix transMatrix=new Matrix(3,3);
			for ( int i=0; i<3; i++ ){
				for ( int j=0; j<3; j++ ){
					transMatrix.matrix[i][j] = m.matrix[j][i];
				}
			}
		return transMatrix;
	}
	
	/**
	*
	* Nasobenie matice a vektora (matrix and vector multiplication)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Matrix (3x3)
	*   Vector (3)
	*   
	*   <navratna hodnota> (<return>) Vector (3 )
	*
	*/
	
	public static Vector  matrixMultiplyVector(Matrix matrix, Vector vector){
		Vector aux=new Vector(3) ;
		double sum;
		for (int i=0; i<3; i++) {
			sum = 0.0;
			for (int j=0; j<3; j++) {
				sum += matrix.matrix[i][j] * vector.v[j];
				aux.v[i] = sum;
                                //System.out.println("VECTOR test " + aux.v[i]);
			}
		}
		return aux;
	}
	
	/**
	*
	* Nasobenie vektora maticou (vector and matrix multiplication)
	*
	* Vstup/Vystup (Input/output):
	*
	*   Matrix (3x3)
	*   Vector (3)
	*   
	*   <navratna hodnota> (<return>) Vector (3 )
	*
	*/
	
	public Vector vectorMultiplyMatrix(Matrix matrix, Vector vector){
		Vector aux=new Vector(3) ;
		double sum;
		for (int i=0; i<3; i++) {
			sum = 0.0;
			for (int j=0; j<3; j++) {
				sum += vector.v[j] * matrix.matrix[i][j];
				aux.v[i] = sum;
			}
		}
		return aux;
	}
        
        /**
         * getDeterminant()
         * 
         * Method to compute determinant of 3x3 matrix
         * 1st numbers are columns, and 2nd are lines 
         * 
         * IN:
         * Mathrix matrix
         * 
         * OUT:
         * double determinant
         * 
         */
        
        public double getDeterminant(Matrix matrix){
            double determinant = 0;
            determinant = matrix.matrix[0][0]*matrix.matrix[1][1]*matrix.matrix[2][2] - 
                          matrix.matrix[0][0]*matrix.matrix[1][2]*matrix.matrix[2][1] -
                          matrix.matrix[0][1]*matrix.matrix[1][0]*matrix.matrix[2][2] +
                          matrix.matrix[0][1]*matrix.matrix[1][2]*matrix.matrix[2][0] +
                          matrix.matrix[0][2]*matrix.matrix[1][0]*matrix.matrix[2][1] -
                          matrix.matrix[0][2]*matrix.matrix[1][1]*matrix.matrix[2][0];
            
            return determinant;
        }
        
        public static void main(String args[]){
            /*
            Matrix matrix = new Matrix(3,3);
            matrix.matrix[0][0] = -2;
            matrix.matrix[1][0] = 2;
            matrix.matrix[2][0] = -3;
            matrix.matrix[0][1] = -1;
            matrix.matrix[1][1] = 1;
            matrix.matrix[2][1] = 3;
            matrix.matrix[0][2] = 2;
            matrix.matrix[1][2] = 0;
            matrix.matrix[2][2] = -1;
            
            double determinant = matrix.getDeterminant(matrix);
            System.out.println("3x3 determinant " + determinant);
            */
            //reference vector, point of view
            Vector refVec = new Vector(3);
            //help vector has the smae direction as x axis (1,0,0)
            Vector helpVec = new Vector(3);
            //reference vector after applying rotation matrix
            Vector outVec = new Vector(3);
            //spherical coordinates of reference vector
            double theta, phi;
            Vector polarAngles = new Vector(3);
            //rotation angles for rotation matrix
            double alpha, beta, gamma;
            //rotation matrix
            Matrix rotMat = new Matrix(3,3);
            Matrix rotMat_2 = new Matrix(3,3);

            //initialization of variables
            double thetaFor = 285;
            double phiFor = 23;
            refVec.v[0] = -0.6078768187253742;
            refVec.v[1] = -0.7506653549675367;
            refVec.v[2] = 0.25881904510252074;
            refVec.v[0] = Math.cos(Math.toRadians(thetaFor))*Math.sin(Math.toRadians(phiFor));
            refVec.v[1] = Math.sin(Math.toRadians(thetaFor))*Math.sin(Math.toRadians(phiFor));
            refVec.v[2] = Math.cos(Math.toRadians(phiFor));

            helpVec.v[0] = 1;
            helpVec.v[1] = 0;
            helpVec.v[2] = 0;

            polarAngles = Vector.fromCartesianToPolar(refVec);
            //System.out.println("theta " + Math.toDegrees(polarAngles.v[0]));
            //System.out.println("phi " + Math.toDegrees(polarAngles.v[1]));
            theta = polarAngles.v[0];
            phi = polarAngles.v[1];
            System.out.println("theta " + Math.toDegrees(theta));
            System.out.println("phi " + Math.toDegrees(phi));
            alpha = 0;
            beta = -phi;
            gamma = theta;
            rotMat = Matrix.rotationMatrix(alpha,beta,gamma);
            //outVec = rotMat.vectorMultiplyMatrix(rotMat, refVec);

            System.out.println("refVec.v[0] " + refVec.v[0]);
            System.out.println("refVec.v[1] " + refVec.v[1]);
            System.out.println("refVec.v[2] " + refVec.v[2]+"\n");
            
            System.out.println("helpVec.v[0] " + helpVec.v[0]);
            System.out.println("helpVec.v[1] " + helpVec.v[1]);
            System.out.println("helpVec.v[2] " + helpVec.v[2]+"\n");

            //System.out.println("outVec.v[0] " + outVec.v[0]);
            //System.out.println("outVec.v[1] " + outVec.v[1]);
            //System.out.println("outVec.v[2] " + outVec.v[2]+"\n");

            //double angleVtoV = Vector.getVectorsAngle(helpVec, outVec);
            //System.out.println("angleVtoV " + Math.toDegrees(angleVtoV));

            outVec = rotMat.vectorMultiplyMatrix(Matrix.R_z(gamma), refVec);
            System.out.println("After 1st rotation around z");
            System.out.println("outVec.v[0] " + outVec.v[0]);
            System.out.println("outVec.v[1] " + outVec.v[1]);
            System.out.println("outVec.v[2] " + outVec.v[2]+"\n");
            outVec = rotMat.vectorMultiplyMatrix(Matrix.R_y(beta), outVec);
            System.out.println("After 2nd rotation around y");
            System.out.println("outVec.v[0] " + outVec.v[0]);
            System.out.println("outVec.v[1] " + outVec.v[1]);
            System.out.println("outVec.v[2] " + outVec.v[2]+"\n");
            double angleVtoV = Vector.getVectorsAngle(helpVec, outVec);
            System.out.println("angleVtoV " + Math.toDegrees(angleVtoV));

            //test 2
            refVec.v[0] = 1;
            refVec.v[1] = 0;
            refVec.v[2] = 0;
            double step = 10;
            int limit = (int)Math.round((360-step)/step);
            for(int i=0;i<limit;i++){
                for(int j=0;j<limit/2;j++){
                    //for(int k=0;k<limit;k++){
                        rotMat=rotMat.R_z(Math.toRadians(i*step));
                        helpVec = rotMat.matrixMultiplyVector(rotMat, refVec);
                        rotMat=rotMat.R_y(Math.toRadians(j*step));
                        helpVec = rotMat.matrixMultiplyVector(rotMat, helpVec);
                        //rotMat=rotMat.R_x(Math.toRadians(k*step));
                        //helpVec = rotMat.matrixMultiplyVector(rotMat, helpVec);
                        System.out.println(helpVec.v[0] + " " + helpVec.v[1] + " " + helpVec.v[2] + " " + Vector.getSize(helpVec));
                        //System.out.println(i + " " + j + " " + k);
                    //}
                }
            }
        }

}