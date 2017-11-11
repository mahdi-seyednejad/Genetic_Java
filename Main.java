/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package geneticimprove;

/**
 *
 * @author mehdi
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

//        FitFunctions myfunction = new FitFunctions(2, 4, -5, 5);
//        myfunction.test();




        /**
         * First you should set the parameters
         */
         int Iteration =500;
         int populationSize = 100;
         int numberOfDimensions = 4;
         int powerAccuracy = 8;
         int decimalAccuracy = (int) Math.pow(10, powerAccuracy);
//         int decimalAccuracy = 2;
         /**
         *
         * You may choose one of theses function and put their name as function:
         * 1- "Sphere"
         * 2- "Ackley"
         * 3- "Rastrigin"
         * 4- "Griewank"
         * 5- "Schwefel"
         * 6- "Rosenbrock"
         */
         String Function = "Sphere";
         GenAlg myGA = new GenAlg(populationSize, Iteration, Function, numberOfDimensions, decimalAccuracy);
         myGA.test();

         double survivalRate = 0.6;
         double mutationRate = 0.1;
         double elitismRate = 0.02;
         int poolSize =2;
         double tournamentProbability = 0.8;
         double desiredAccuracy;
//         desiredAccuracy = 0.05;
         desiredAccuracy = Math.pow(10, -7);
//         myGA.SGA(survivalRate, mutationRate, elitismRate, poolSize, tournamentProbability, desiredAccuracy);
//         myGA.DPGA(survivalRate, mutationRate, elitismRate, poolSize, tournamentProbability, desiredAccuracy);
//         myGA.ImprovedDPGA(survivalRate, mutationRate, elitismRate, poolSize, tournamentProbability, desiredAccuracy);


         int NumOfRun = 1;
         myGA.multipleRun(NumOfRun, survivalRate, mutationRate, elitismRate, poolSize, tournamentProbability, desiredAccuracy);

//         System.out.println(-Math.log10(Math.pow(10, -12))/Math.log10(650));
         


    }

}
