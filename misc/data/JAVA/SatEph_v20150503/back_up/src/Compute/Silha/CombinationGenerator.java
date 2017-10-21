/*
 * CombinationGenerator.java
 *
 * Created on October 23, 2008, 12:21 PM
 *
 * Class to generate combinations of array of numbers.
 */

package Compute.Silha;

/**
 *
 * @author Jiri Silha
 */
public class CombinationGenerator {
    
    /** Creates a new instance of CombinationGenerator */
    public CombinationGenerator() {
    }
    
    /**
     * getCombitanions()
     *
     * Method to get matrix with combinations of given array of integers. 
     * Inthis method is the size of the tuples to combinate 3. So k = 3.
     *
     * IN: 
     *  int array[]         //array of integeres
     *  boolean repetition  //should or shouldn't the no. repeated
     * 
     * OUT:
     *  Vector[] - array of vectors, all combinations
     *
     * BEWARE!!! - max value of array length is appr. 170!!!
     */ 
    
    public static Vector[] getCombinations(int array[], boolean repetition){
        Vector vectors[];
        
        //in our case the size of the tuples to combinate is 3
        int k = 3;
        //with repetition
        if(repetition){
            //for now do nothing
            vectors = new Vector[3];
        }
        //without repetition
        else {
            int n = array.length;
            //System.out.println("array length " + n);
            
            //formula for no of combinations - without repitation the numbers
            //int noOfCombinations = (int) (factorial(n)/(factorial(n - k))/factorial(k));
            
            //correction of prevoius line, there is a max value appr. 170, because of factorial 170 is too big no.
            int helpComb = n;
            int helpComb2 = 1;

            while(helpComb2 < k){
                helpComb = helpComb * (n - helpComb2);
                helpComb2++;
            }
            
            int noOfCombinations = (int)((double)helpComb/factorial(k));
                      
            //System.out.println("noOfCombinations " + noOfCombinations);
            vectors = new Vector[noOfCombinations];
            //System.out.println("vectors length " + vectors.length);
            
            //filling matrix
            int amount = 0;
            for(int i = 0; i < (n-2); i++){
                for(int j = i + 1; j < (n-1); j++){
                    for(int l = j + 1; l < n; l++){
                        //
                        vectors[amount] = new Vector(3);
                        vectors[amount].v[0] = array[i]; 
                        vectors[amount].v[1] = array[j]; 
                        vectors[amount].v[2] = array[l];
                        
                        //System.out.println("vectors[i].v[0] " + vectors[i].v[0]);
                        //System.out.println("Array[i] " + array[i]);
                        //System.out.println("Amount " + amount);
                        
                        amount = amount + 1;
                        /*if(amount >= noOfCombinations) {
                            l = n;
                            j = n;
                            i = n;
                        }*/
                        //System.out.println("vectors[i].v[0] " + vectors[i].v[0]);
                        //System.out.println("Array[i] " + array[i]);
                        //System.out.println(amount);
                    }
                }
            }
            
        }
        
        return vectors;
    }
    
    
    /**
     * factorial{}
     *
     * Method to compute factorials.
     *
     * IN: 
     *  int no
     * OUT: 
     *  int factorial
     */
    
    public static double factorial(double no){
        double factorialNo = 0;
        if(no >= 0){
         if( no <= 1 )     // base case
             factorialNo = 1;
         else{
             factorialNo = no;
             while(no > 1){
                 factorialNo = factorialNo * (no - 1 );
                 no = no - 1;
                 //System.out.println(factorialNo);
             }
         }
        }
        
        else{
            //System.out.println("Negative number!!!");
        }
        return factorialNo;
    }
    
    
    //test method
    public static void main (String args[]){
        //System.out.println(factorial(70)/10e20 + "");
        //System.out.println(factorial(35) + "");
        int arrayLength = 6;
        int array[] = new int[arrayLength];
        for(int i = 0; i < array.length; i++){
            array[i] = i;
        }
        
        Vector vectors[] = getCombinations(array, false);
        //System.out.println(vectors.length);
        for(int i = 0; i < vectors.length; i++){
            System.out.println(vectors[i].v[0] + " " +  vectors[i].v[1] + " " + vectors[i].v[2]);
        }
    }
    
}
