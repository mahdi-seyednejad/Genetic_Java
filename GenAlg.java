/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package geneticimprove;

import java.util.Random;

/**
 *
 * @author mehdi
 */
public class GenAlg {
    Plot plot;
    int iteration;
    Population basicPopulation;
    Population BasicChildren;
    Population SelectedBasicPop;
    Population MainPop;
    Population ReservePop;
    Population Mchildren;
    Population Rchildren;
    Population selectedMain;
    Population selectedReserve;
    Population complBasicPop;
    Population complMainPop;
    Population myMainPop;
    Population myReservePop;
    Population myMchildren;
    Population myRchildren;
    Population myselectedMain;
    Population myselectedReserve;

    FitFunctions Fitness;
    double [][] decimalPop;
    int [][] NewPop;
    int globalIndexMain;
    int globalIndexReserve;
    
//    int [][] selectedMain;
//    int [][] Newreserve;
    String functionName;
    int functionNum;
    double Xmin;
    double Xmax;
    int Dimension;
    int geneSize;
    int chromSize;
    int populationSize;
    int decimalAcc;

    public GenAlg() {
    }
    public GenAlg(int popSize , int Iteration , String Function, int numberOfDimensions, int decimalAccuracy)
    {
        iteration = Iteration;
        functionName = Function;
        Dimension = numberOfDimensions;
        if (Function.equalsIgnoreCase("Sphere"))
        {
            Xmin= -5.12;
            Xmax = 5.12;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }
        if (Function.equalsIgnoreCase("Ackley"))
        {
            Xmin= -32.768;
            Xmax = 32.768;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }
        if (Function.equalsIgnoreCase("Rastrigin"))
        {
            Xmin= -5.12;
            Xmax = 5.12;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }
        if (Function.equalsIgnoreCase("Griewank"))
        {
            Xmin= -600;
            Xmax = 600;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }
        if (Function.equalsIgnoreCase("Schwefel"))
        {
            Xmin= -500;
            Xmax = 500;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }
        if (Function.equalsIgnoreCase("Rosenbrock"))
        {
            Xmin= -5;
            Xmax = 10;
//            basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
//            ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        }

        basicPopulation = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        ReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        MainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        myReservePop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        myMainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        complBasicPop = new Population(popSize, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        complBasicPop.complementPopulation(basicPopulation);
        complMainPop = new Population(popSize/2, numberOfDimensions, Xmin, Xmax, decimalAccuracy);
        complMainPop.complementPopulation(MainPop);

        geneSize = basicPopulation.eachDimensionLength;
        chromSize = geneSize*Dimension;
        Fitness = new FitFunctions(Dimension, geneSize, Xmin, Xmax);
        globalIndexMain=0;
        globalIndexReserve=0;
        populationSize = popSize;
        decimalAcc = decimalAccuracy;
        SetFunction(Function);
    }// End of constructor
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
    public void SetFunction(String Fname)
    {
        if (Fname.equalsIgnoreCase("Sphere"))
            functionNum = 1;
        if (Fname.equalsIgnoreCase("Ackley"))
            functionNum = 2;
        if (Fname.equalsIgnoreCase("Rastrigin"))
            functionNum = 3;
        if (Fname.equalsIgnoreCase("Griewank"))
            functionNum = 4;
        if (Fname.equalsIgnoreCase("Schwefel"))
            functionNum = 5;
        if (Fname.equalsIgnoreCase("Rosenbrock"))
            functionNum = 6;
    }

    public void computeDecimalPopulation()
    {

    }// End of compute decimal


    /**
     * This method converts binary to decimal
     */
    public double bin2Decimal (int [] binNumber)
    {
        double realNum=0;
        int pow=0;
        for (int i=0; i<binNumber.length; i++)
        {
            pow = binNumber.length-1-i;
            if (binNumber[i]==1)
                realNum = realNum + Math.pow(2, pow);
        }
        return realNum;
    }// end of binari to real
    /**
     * This method performs selection to survive in standard GA
     */
    public void selectionStandardGA(double survivalRate, int poolSize, double tournProb)
    {
        int surviveSize=(int)(survivalRate*basicPopulation.popMatrix.length);
        int [] chosen = new int [basicPopulation.popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int countSurviver=0;
        int surviverIndex;
        while (countSurviver<surviveSize ) // Main basicPopulation selection
        {
            surviverIndex = basicPopulation.tournamentMinIndex(poolSize, tournProb);
            if (chosen[surviverIndex]==0)
            {
                SelectedBasicPop.insertIndividual(basicPopulation.popMatrix[surviverIndex]);
                chosen[surviverIndex]=1;
                countSurviver++;
            }
        }
    }// End of standard population selection
    public void crossoverAllStandard(double childRate, int poolSize, double tournProb)
    {
        int reproductionRate=(int)(childRate*basicPopulation.popMatrix.length);
        int [][] tempChild = new int [2] [basicPopulation.popMatrix[0].length];
        int p1=0;
        int p2=0;
        int childCount=0;
        do
        {
            do
            {
                p1 = basicPopulation.tournamentMinIndex(poolSize, tournProb);
                p2 = basicPopulation.tournamentMinIndex(poolSize, tournProb);
            } while (p1 == p2);
            tempChild = crossover2Point(basicPopulation.popMatrix[p1], basicPopulation.popMatrix[p2]);
            for (int i=0; i<tempChild.length; i++)
                BasicChildren.insertIndividual(tempChild[i]);

            childCount++;
        }while (childCount<reproductionRate);


    }
    /**
     * This method will calculate fitness function for its input basicPopulation
     */
    /**
     * This method performs selection to survive
     */
    public void selectionDualPop(double survivalRate, int poolSize, double tournProb)
    {
        int surviveSize=(int)(survivalRate*MainPop.popMatrix.length);
        int [] chosen = new int [MainPop.popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int countSurviver=0;
        int surviverIndex;
        while (countSurviver<surviveSize ) // Main basicPopulation selection
        {
            surviverIndex = MainPop.tournamentMinIndex(poolSize, tournProb);
            if (chosen[surviverIndex]==0)
            {
                selectedMain.insertIndividual(MainPop.popMatrix[surviverIndex]);
                chosen[surviverIndex]=1;
                countSurviver++;
            }
        }
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        countSurviver=0;
        while (countSurviver<surviveSize ) // Reserve basicPopulation selection
        {
            surviverIndex = ReservePop.tournamentMaxIndex(poolSize, tournProb);
            if (chosen[surviverIndex]==0)
            {
                selectedReserve.insertIndividual(ReservePop.popMatrix[surviverIndex]);
                chosen[surviverIndex]=1;
                countSurviver++;
            }
        }

    }// End of selectionDualPop
    /**
     * This method performs selection to survive
     */
    public void selectionDualPop(Population primaryPop, Population selectedMPop, Population reserve, Population selectedRpop, double survivalRate, int poolSize, double tournProb)
    {
        int surviveSize=(int)(survivalRate*primaryPop.popMatrix.length);
        int [] chosen = new int [primaryPop.popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int countSurviver=0;
        int surviverIndex;
        while (countSurviver<surviveSize ) // Main basicPopulation selection
        {
            surviverIndex = primaryPop.tournamentMinIndex(poolSize, tournProb);
            if (chosen[surviverIndex]==0)
            {
                selectedMPop.insertIndividual(primaryPop.popMatrix[surviverIndex]);
                chosen[surviverIndex]=1;
                countSurviver++;
            }
        }
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        countSurviver=0;
        while (countSurviver<surviveSize ) // Reserve Population selection
        {
            surviverIndex = reserve.tournamentMaxIndex(poolSize, tournProb);
            if (chosen[surviverIndex]==0)
            {
                selectedRpop.insertIndividual(reserve.popMatrix[surviverIndex]);
                chosen[surviverIndex]=1;
                countSurviver++;
            }
        }

    }// End of selectionDualPop
    public void fitnessForAllMain(Population pop)
    {
        for (int i=0; i<pop.popMatrix.length; i++)
            pop.fitness[i] = Fitness.appropriateFitness(pop.popMatrix[i], functionNum);

    }
    /**
     * This method will calculate fitness function for its input basicPopulation
     */
    public void fitnessForAllReserve()
    {
        for (int i=0; i<ReservePop.popMatrix.length; i++)
            ReservePop.fitness[i] = Fitness.diversityFrequencyOnGene(MainPop.popMatrix, ReservePop.popMatrix[i]);
    }// End of fitnessForAllReserve
    /**
     * This method will calculate fitness function for its input basicPopulation
     */
    public void fitnessForAllReserve(Population primaryPop, Population reserve)
    {
        for (int i=0; i<reserve.popMatrix.length; i++)
            reserve.fitness[i] = Fitness.diversityFrequencyOnGene(primaryPop.popMatrix, reserve.popMatrix[i]);
    }// End of fitnessForAllReserve


    /**
     * This method performs sufficient number of crossover on dual basicPopulation
     */
    public void crossoverAllDualPopulation(double crossRate, int PoolSize, double TourProb)
    {
        int crossTimes = (int) (crossRate*MainPop.popMatrix.length/2);
        for (int i=0; i<crossTimes; i++)
            crossoverDualPop(PoolSize, TourProb);
    }// End of crossoverAllDualPopulation
    /**
     * This method performs sufficient number of crossover on dual basicPopulation for Improved DPGA
     */
    public void crossoverAllDualPopImproved(double crossRate, Population PrimaryPop, Population ChildMain, Population Resrve, Population ChildReserve,  int PoolSize, double TourProb, double reserve2Chance)
    {
        int crossTimes = (int) (crossRate*MainPop.popMatrix.length/2);
        for (int i=0; i<crossTimes; i++)
            crossoverImprovesdDualPop(PrimaryPop, ChildMain, Resrve, ChildReserve, PoolSize, TourProb, reserve2Chance);
    }// End of crossoverAllDualPopulation

    /**
     * This method performs crossover on both main and reserve basicPopulation
     * It inserts the created childrenMain in the appropriate basicPopulation
     */
    public void crossoverDualPop(int PoolSize, double tourProb)
    {
        int [][] Mchild = new int [2][chromSize];
        double [] MchildFit = new double [2];
        int [][] Rchild = new int [2][chromSize];
        double [] RchildFit = new double [2];
        int [][] MRchild = new int [2][chromSize];
        double [] MRchildFit = new double [2];
        int [] parent1 = new int [chromSize];
        int [] parent2 = new int [chromSize];
        int p1;
        int p2;
        do
        {
            p1 = MainPop.tournamentMinIndex(PoolSize, tourProb);
            p2 = MainPop.tournamentMinIndex(PoolSize, tourProb);
        } while (p1==p2);
        Mchild = crossover2Point(MainPop.popMatrix[p1], MainPop.popMatrix[p1]);
        if (MainPop.fitness[p1]<MainPop.fitness[p2])
        {
            int temp=p1;
            p1=p2;
            p2=temp;
        }
        do
        {
            parent1 = tournamentCrossBreed(PoolSize, tourProb, MainPop.popMatrix[p2], ReservePop);
            parent2 = tournamentCrossBreed(PoolSize, tourProb, MainPop.popMatrix[p2], ReservePop);
        } while (Fitness.isEqualVector(parent1, parent2));
        Rchild = crossover2Point(parent1, parent2);

        if ( Fitness.HammingDistance(parent1, MainPop.popMatrix[p2]) < Fitness.HammingDistance(parent2, MainPop.popMatrix[p2]) )
            MRchild = crossover2Point(MainPop.popMatrix[p2], parent1);
        else
            MRchild = crossover2Point(MainPop.popMatrix[p2], parent2);
        //////////////////////////////////   Children are created  ////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        /////////////  Find better childrenMain for main basicPopulation
        for (int i=0; i<MchildFit.length; i++)
            switch (functionNum)
            {
                case 1: MchildFit[i] = Fitness.SphereFit(Mchild[i]);
                        break;
                case 2: MchildFit[i] = Fitness.AckleyFit(Mchild[i]);
                        break;
                case 3: MchildFit[i] = Fitness.RastriginFit(Mchild[i]);
                        break;
                case 4: MchildFit[i] = Fitness.GriewankFit(Mchild[i]);
                        break;
                case 5: MchildFit[i] = Fitness.SchwefelFit(Mchild[i]);
                        break;
                case 6: MchildFit[i] = Fitness.RosenbrockFit(Mchild[i]);
                        break;
            }

            for (int i=0; i<MRchildFit.length; i++)
            switch (functionNum)
            {
                case 1: MRchildFit[i] = Fitness.SphereFit(MRchild[i]);
                        break;
                case 2: MRchildFit[i] = Fitness.AckleyFit(MRchild[i]);
                        break;
                case 3: MRchildFit[i] = Fitness.RastriginFit(MRchild[i]);
                        break;
                case 4: MRchildFit[i] = Fitness.GriewankFit(MRchild[i]);
                        break;
                case 5: MRchildFit[i] = Fitness.SchwefelFit(MRchild[i]);
                        break;
                case 6: MRchildFit[i] = Fitness.RosenbrockFit(MRchild[i]);
                        break;
            }
            int [][] candidates = new int [4][chromSize];
            double [] candFit = new double [4];
            for (int i=0; i<2; i++)
            {
                candFit[i]=MchildFit[i];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = Mchild[i][j];
            }
            for (int i=2; i<4; i++)
            {
                candFit[i]=MRchildFit[i-2];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = MRchild[i-2][j];
            }
            int [][] selected = findBetterChildren(candidates, candFit, "min");
            for (int i=0; i<selected.length; i++ )
                Mchildren.insertIndividual(selected[i]);
            //////////////////////////////////  The best childrenMain are found and inserted in the main childrenMain basicPopulation
            ///////////////////////////////////////////////////////////////////////////
            for (int i=0; i<MRchildFit.length; i++)
                MRchildFit[i] = Fitness.diversityFrequencyOnGene(MainPop.popMatrix, MRchild[i]);
            for (int i=0; i<RchildFit.length; i++ )
                RchildFit[i] = Fitness.diversityFrequencyOnGene(MainPop.popMatrix, Rchild[i]);

            for (int i=0; i<2; i++)
            {
                candFit[i]=RchildFit[i];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = Rchild[i][j];
            }
            for (int i=2; i<4; i++)
            {
                candFit[i]=MRchildFit[i-2];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = MRchild[i-2][j];
            }
            selected = findBetterChildren(candidates, candFit, "max");
            for (int i=0; i<selected.length; i++ )
                Rchildren.insertIndividual(selected[i]);
            //////////////////////////////////  The best childrenMain are found and inserted in the reserve childrenMain basicPopulation

    }// End of crossover on dual basicPopulation
    /**
     * This method performs crossover on both main and reserve basicPopulation
     * It inserts the created childrenMain in the appropriate basicPopulation
     */
    public void crossoverImprovesdDualPop(Population primaryPop,Population childrenMain, Population reserve, Population childrenReserve, int PoolSize, double tourProb, double reserveChance)
    {
        int [][] Mchild = new int [2][chromSize];
        double [] MchildFit = new double [2];
        int [][] Rchild = new int [2][chromSize];
        double [] RchildFit = new double [2];
        int [][] MRchild = new int [2][chromSize];
        double [] MRchildFit = new double [2];
        int [] parent1 = new int [chromSize];
        int [] parent2 = new int [chromSize];
        int p1;
        int p2;
        do
        {
            p1 = primaryPop.tournamentMinIndex(PoolSize, tourProb);
            p2 = primaryPop.tournamentMinIndex(PoolSize, tourProb);
        } while (p1==p2);
        Mchild = crossover2Point(primaryPop.popMatrix[p1], primaryPop.popMatrix[p1]);
        if (primaryPop.fitness[p1]<primaryPop.fitness[p2])
        {
            int temp=p1;
            p1=p2;
            p2=temp;
        }
        do
        {
            parent1 = tournamentCrossBreed(PoolSize, tourProb, primaryPop.popMatrix[p2], reserve);
            parent2 = tournamentCrossBreed(PoolSize, tourProb, primaryPop.popMatrix[p2], reserve);
        } while (Fitness.isEqualVector(parent1, parent2)); // finding better parent from reserve
        Rchild = crossover2Point(parent1, parent2);

        if ( Fitness.HammingDistance(parent1, primaryPop.popMatrix[p2]) < Fitness.HammingDistance(parent2, primaryPop.popMatrix[p2]) )
            MRchild = crossover2Point(primaryPop.popMatrix[p2], parent1);
        else
            MRchild = crossover2Point(primaryPop.popMatrix[p2], parent2);
        //////////////////////////////////   Children are created  ////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        /////////////  Find better childrenMain for main basicPopulation
        for (int i=0; i<MchildFit.length; i++)
            switch (functionNum)
            {
                case 1: MchildFit[i] = Fitness.SphereFit(Mchild[i]);
                        break;
                case 2: MchildFit[i] = Fitness.AckleyFit(Mchild[i]);
                        break;
                case 3: MchildFit[i] = Fitness.RastriginFit(Mchild[i]);
                        break;
                case 4: MchildFit[i] = Fitness.GriewankFit(Mchild[i]);
                        break;
                case 5: MchildFit[i] = Fitness.SchwefelFit(Mchild[i]);
                        break;
                case 6: MchildFit[i] = Fitness.RosenbrockFit(Mchild[i]);
                        break;
            }

            for (int i=0; i<MRchildFit.length; i++)
            switch (functionNum)
            {
                case 1: MRchildFit[i] = Fitness.SphereFit(MRchild[i]);
                        break;
                case 2: MRchildFit[i] = Fitness.AckleyFit(MRchild[i]);
                        break;
                case 3: MRchildFit[i] = Fitness.RastriginFit(MRchild[i]);
                        break;
                case 4: MRchildFit[i] = Fitness.GriewankFit(MRchild[i]);
                        break;
                case 5: MRchildFit[i] = Fitness.SchwefelFit(MRchild[i]);
                        break;
                case 6: MRchildFit[i] = Fitness.RosenbrockFit(MRchild[i]);
                        break;
            }
            int [][] candidates = new int [4][chromSize];
            double [] candFit = new double [4];
            for (int i=0; i<2; i++)
            {
                candFit[i]=MchildFit[i];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = Mchild[i][j];
            }
            for (int i=2; i<4; i++)
            {
                candFit[i]=MRchildFit[i-2];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = MRchild[i-2][j];
            }
            int [][] selected = findBetterChildren(candidates, candFit, "min");// Choose candidates

            for (int i=0; i<selected.length; i++)// Give another chance to reserve childrenMain
                for (int k=0; k<MRchild.length; k++)
                    if (!isEqualVectors(selected[i], MRchild[k]))
                        if (myRandom()<reserveChance)
                        {
                            for (int j=0; j<selected[0].length; j++)
                                selected[i][j] = MRchild[k][j];
                            break;
                        }

            for (int i=0; i<selected.length; i++ )
                childrenMain.insertIndividual(selected[i]);
            //////////////////////////////////  The best childrenMain are found and inserted in the main childrenMain basicPopulation
            ///////////////////////////////////////////////////////////////////////////
            for (int i=0; i<MRchildFit.length; i++)
                MRchildFit[i] = Fitness.diversityFrequencyOnGene(primaryPop.popMatrix, MRchild[i]);
            for (int i=0; i<RchildFit.length; i++ )
                RchildFit[i] = Fitness.diversityFrequencyOnGene(primaryPop.popMatrix, Rchild[i]);

            for (int i=0; i<2; i++)
            {
                candFit[i]=RchildFit[i];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = Rchild[i][j];
            }
            for (int i=2; i<4; i++)
            {
                candFit[i]=MRchildFit[i-2];
                for (int j=0; j<chromSize; j++)
                    candidates[i][j] = MRchild[i-2][j];
            }
            selected = findBetterChildren(candidates, candFit, "max");
            for (int i=0; i<selected.length; i++ )
                childrenReserve.insertIndividual(selected[i]);
            //////////////////////////////////  The best childrenMain are found and inserted in the reserve childrenMain basicPopulation

    }// End of crossover on dual basicPopulation
    /**
     * This method returns true if vector1 and vector2 are equal
     */
    public boolean isEqualVectors(int [] vector1, int [] vector2)
    {
        for (int i=0; i<vector1.length; i++)
            if (vector1[i]!=vector2[i])
                return false;
        return true;
    }//  End of is equal
    /**
     * This method performs 2 points crossover
     */
    public int [][] crossover2Point(int []v1, int [] v2)
    {
        int [][] temp = new int [2][v1.length];
        int p1 = myRandomInt(v1.length);
        int p2 = myRandomInt(v2.length);
        int Ptemp;
        if (p2<p1)
        {
            Ptemp=p1;
            p1 = p2;
            p2 = Ptemp;
        }

        for (int i=0; i<=p1; i++)
            temp[0][i]= v1[i];
        for (int i=p1+1; i<=p2; i++)
            temp[0][i] = v2[i];
        for (int i=p2+1; i<v1.length; i++)
            temp[0][i] = v1[i];

        for (int i=0; i<=p1; i++)
            temp[1][i]= v2[i];
        for (int i=p1+1; i<=p2; i++)
            temp[1][i] = v1[i];
        for (int i=p2+1; i<v1.length; i++)
            temp[1][i] = v2[i];

        return temp;
    }// end of crossover 2 point
    /**
     * This method will find 2 best individual between its inputs and returns them
     */
    public int [][] findBetterChildren(int [][] Matrix, double [] fit, String MaxMin)
    {
        int rowSize = Matrix[0].length;
        int [][] mySelected = new int [2][rowSize];
        int [] tempVector = new int [rowSize];
        double tempFit =0;
        if (MaxMin.equalsIgnoreCase("min"))
            for (int i=0; i<fit.length; i++)
                for (int j=i; j<fit.length; j++)
                {
                    if (fit[i]>fit[j])
                    {
                        for (int k=0; k<rowSize; k++ )
                            tempVector[k]=Matrix[i][k];
                        for (int k=0; k<rowSize; k++ )
                            Matrix[i][k]=Matrix[j][k];
                        for (int k=0; k<rowSize; k++ )
                            Matrix[j][k]=tempVector[k];
                        tempFit = fit[i];
                        fit[i] = fit[j];
                        fit[j] = tempFit;
                    }
                }
        else
            for (int i=0; i<fit.length; i++)
                for (int j=i; j<fit.length; j++)
                {
                    if (fit[i]<fit[j])
                    {
                        for (int k=0; k<rowSize; k++ )
                            tempVector[k]=Matrix[i][k];
                        for (int k=0; k<rowSize; k++ )
                            Matrix[i][k]=Matrix[j][k];
                        for (int k=0; k<rowSize; k++ )
                            Matrix[j][k]=tempVector[k];
                        tempFit = fit[i];
                        fit[i] = fit[j];
                        fit[j] = tempFit;
                    }
                }//////////////////  Now, the matrix of candidates are sorted

        for (int i=0; i<mySelected.length; i++)
            for (int j=0; j<rowSize; j++)
                mySelected[i][j] = Matrix [i][j];
        return mySelected;
    }
    /**
     * This method generate a random integer number between 0 and n
     */
    public int myRandomInt(int n)
    {
        Random randomNumbers = new Random();
        int random1= randomNumbers.nextInt(n);
        return random1;

    }// end of rand int
    /**
     * This method generate a random number between 0 and 1
     */
    public double myRandom()
    {
        Random randomNumbers = new Random();
        double random1= randomNumbers.nextDouble();
        return random1;

    }// End of random (0,1)
    /**
     * This method pick an chIndex of an individual proportional to fitness value
     * base on Rollette Wheel
     */
    public int RolletteWheel (Population myPop)
    {
        double sum =0;
        double partSum=0;
        for (int i=0; i<myPop.fitness.length; i++)
            sum = sum + myPop.fitness[i];
        double rand = myRandom();
        rand = rand * sum;
        int dataIndex =-1;
        do
        {
            dataIndex ++;
            if (dataIndex == 60)
                System.out.println(dataIndex);
            partSum = partSum + myPop.fitness[dataIndex];
        }while (partSum<=rand);
//        dataIndex++;
        return dataIndex;
    } // End of Rollette Wheel

    /**
     * This method will perform mutation
     */
    public void mutation(Population pop, double mutRate)
    {
        int mutSize;
        if(mutRate<0 || mutRate>1) // it workd with the chromosome size
            mutSize = (int) ((1/chromSize) * pop.popMatrix.length);
        else
            mutSize = (int) (mutRate * pop.popMatrix.length);
        int mutTimes =0;
        int [] chosen = new int [pop.popMatrix.length];
        int geneIndex;
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int chIndex=0;
        while (mutTimes<mutSize)
        {
            chIndex = myRandomInt(pop.popMatrix.length);
            if (chosen[chIndex] == 0)
            {
                geneIndex = myRandomInt(pop.popMatrix[0].length);
                if (pop.popMatrix[chIndex][geneIndex] == 0)
                    pop.popMatrix[chIndex][geneIndex] = 1;
                else
                    pop.popMatrix[chIndex][geneIndex] = 0;
                chosen[chIndex]=1;
                mutTimes++;
            }
        }
    }// End of mutation with mutation ratio
    /**
     * This method will perform mutation base on chromozome size
     */
    public void mutation(Population pop)
    {
        int mutSize = (int) ((1/chromSize) * pop.popMatrix.length);
        int mutTimes =0;
        int [] chosen = new int [pop.popMatrix.length];
        int geneIndex;
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int chIndex=0;
        while (mutTimes<mutSize)
        {
            chIndex = myRandomInt(pop.popMatrix.length);
            if (chosen[chIndex] == 0)
            {
                geneIndex = myRandomInt(pop.popMatrix[0].length);
                if (pop.popMatrix[chIndex][geneIndex] == 0)
                    pop.popMatrix[chIndex][geneIndex] = 1;
                else
                    pop.popMatrix[chIndex][geneIndex] = 0;
                chosen[chIndex]=1;
                mutTimes++;
            }
        }
    }// End of mutation base on chromosome
    /**
     * This method prints results on each iteration
     */
    public void printResults(double [] results)
    {
        System.out.println();
        for (int i=0; i<results.length; i++)
            System.out.print("Iteratio#: "+i+" => "+ results[i]);
        System.out.println();
    }// End of print results
    /**
     * This method finds the closest vector (in basicPopulation "Rpop") to parent
     */
    public int [] findClosestHD(Population Rpop, int [] parent)
    {
        int [] hammingD = new int [Rpop.popMatrix.length];
        for (int i=0; i<hammingD.length; i++)
            hammingD[i]= Fitness.HammingDistance(Rpop.popMatrix[i], parent);
        int minHD = hammingD[0];
        int minHdIndex = 0;
        for (int i=1; i<hammingD.length; i++)
            if (hammingD[i]<minHD)
            {
                minHD=hammingD[i];
                minHdIndex = i ;
            }

        return Rpop.popMatrix[minHdIndex];
    } // end of find closest hamming distance

    public int [] tournamentCrossBreed(int poolSize, double prob, int [] parent, Population Rpop)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [Rpop.popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(Rpop.popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        double minHD = Fitness.HammingDistance(parent, Rpop.popMatrix[pool[0]]);
        int minHdIndex = pool[0];
        for (int i=1; i<pool.length; i++)
            if (Fitness.HammingDistance(parent, Rpop.popMatrix[pool[i]]) < minHD)
            {
                minHdIndex = pool[i];
                minHD = Fitness.HammingDistance(parent, Rpop.popMatrix[pool[i]]);
            }
        double rand2=0;
        if (poolSize == 2)
        {
            rand2 = myRandom();
            if (rand2>prob)
                if (minHdIndex == pool[0])
                    minHdIndex = pool[1];
                else
                    minHdIndex = pool[0];
        }


        return Rpop.popMatrix[minHdIndex];
    }// End of tournamentMin

    public void insertNewChildInNewPop(Population myNewPop, Population myChildPop)
    {
        for (int i=0; i<myChildPop.popMatrix.length; i++)
            myNewPop.insertIndividual(myChildPop.popMatrix[i]);
    }
    /**
     * This method updates the primary population with childrenMain and selected populations
     */
    public void updatePopulation(Population primary, Population child, Population selected)
    {
        int [][] matrix = new int [primary.popMatrix.length][primary.popMatrix[0].length];
        for (int i=0; i<child.popMatrix.length; i++)
            for (int j=0; j<child.popMatrix[0].length; j++)
                matrix[i][j] = child.popMatrix[i][j];
//        System.out.println();
        for (int i=0; i<selected.popMatrix.length; i++)
            for (int j=0; j<selected.popMatrix[0].length; j++)
                matrix[i+child.popMatrix.length][j] = selected.popMatrix[i][j];

        for (int i=0; i<matrix.length; i++)
            for (int j=0; j<matrix[0].length; j++)
                primary.popMatrix[i][j] = matrix[i][j];

        primary.setZeroInsertingIndex();
        child.setZeroInsertingIndex();
        selected.setZeroInsertingIndex();

    }// End of updatePopulation
    /**
     * This method applies elitism
     */
    public void keepElite(Population primaryPop, Population selectedPop, int numOfElites, String MaxMin)
    {
        int bestIndex;
        for (int i=0; i<numOfElites; i++)
        {
            if (MaxMin.equalsIgnoreCase("min"))
            bestIndex= primaryPop.minFitIndex();
            else
                bestIndex= primaryPop.maxFitIndex();
            selectedPop.insertIndividual(primaryPop.popMatrix[bestIndex]);
        }
        

    }// End of keepElite


    /**
     * This method will perform the standatd GA
     */
    public void SGA(double survivalRate, double mutationRate, double EliteRate, int poolSize, double tournamentProbability, double bestAccuracy)
    {
        int eliteSize;
        if (EliteRate>=1)
            eliteSize = (int) EliteRate;
        else
            eliteSize = (int) (EliteRate*basicPopulation.popMatrix.length);

        double reproductionRate = 1- survivalRate;
        int selectedPopSize = (int)(survivalRate*populationSize);

        selectedPopSize = selectedPopSize - eliteSize;

        SelectedBasicPop = new Population(selectedPopSize, Dimension, Xmin, Xmax, decimalAcc);
        int childSize= (int) (reproductionRate*populationSize);
        BasicChildren = new Population(childSize, Dimension, Xmin, Xmax, decimalAcc);
//        double [] Results = new double [iteration];
        double bestFitNow=100;
        int lastIteration = iteration;
        for (int iter=0; iter<iteration; iter++)
        {
            basicPopulation.setZeroInsertingIndex();
            SelectedBasicPop.setZeroInsertingIndex();
            fitnessForAllMain(basicPopulation);
            keepElite(basicPopulation, SelectedBasicPop, eliteSize, "min");
            selectionStandardGA(survivalRate, poolSize, tournamentProbability);
            crossoverAllStandard(reproductionRate, poolSize, tournamentProbability);
            mutation(BasicChildren, mutationRate);
            updatePopulation(basicPopulation, BasicChildren, SelectedBasicPop);
            bestFitNow = basicPopulation.minFit();
            basicPopulation.addBest(bestFitNow);
            basicPopulation.setLastGeneration(iter+1);
            if (bestFitNow<bestAccuracy)
            {
                lastIteration = iter+1;
                //break;
            }
                
        }
        double [] results = basicPopulation.bestResult();
        double [] myres= new double [results.length];
        double [] x = new double [results.length];
        for (int i=0; i<results.length; i++)
        {
            System.out.println(results[results.length-i-1]);
            myres[i]=results[results.length-i-1];
            x[i]=i+1;
        }

        plot = new Plot(functionName, x, myres);
//        printVector(results);
        System.out.println("SGA: Break in iteration number: "+ lastIteration + "  Best result: "+ bestFitNow);
//        double normal= -Math.log10(bestFitNow)/Math.log(lastIteration);
//        System.out.println("normalize evaluation: "+ normal);
//        normal= -Math.log10(bestFitNow)/Math.log10(lastIteration);
//        System.out.println("normalize evaluation (10): "+ normal);

    }// End of Standard GA
    /**
     * This method will perform the dual population GA
     */
    public void DPGA(double survivalRate, double mutationRate, double EliteRate, int poolSize, double tournamentProbability, double bestAccuracy)
    {
        int eliteSize;
        if (EliteRate>=1)
            eliteSize = (int) EliteRate;
        else
            eliteSize = (int) (EliteRate*MainPop.popMatrix.length);

        double reproductionRate = 1- survivalRate;
        int selectedPopSize = (int)(survivalRate*populationSize/2);

        selectedPopSize = selectedPopSize - eliteSize;

        selectedMain = new Population(selectedPopSize, Dimension, Xmin, Xmax, decimalAcc);
        int childSize= (int) (reproductionRate*populationSize/2);
        Mchildren = new Population(childSize, Dimension, Xmin, Xmax, decimalAcc);
        selectedReserve = new Population(selectedPopSize, Dimension, Xmin, Xmax, decimalAcc);
        Rchildren = new Population(childSize, Dimension, Xmin, Xmax, decimalAcc);
//        double [] Results = new double [iteration];
        double bestFitNow = 100;
        int lastIteration = iteration;
        for (int iter=0; iter<iteration; iter++)
        {
            MainPop.setZeroInsertingIndex();
            selectedMain.setZeroInsertingIndex();
            fitnessForAllMain(MainPop);
            keepElite(MainPop, selectedMain, eliteSize, "min");

            ReservePop.setZeroInsertingIndex();
            selectedReserve.setZeroInsertingIndex();
            fitnessForAllReserve();
            keepElite(ReservePop, selectedReserve, eliteSize, "max");

            selectionDualPop(survivalRate, poolSize, tournamentProbability);
            crossoverAllDualPopulation(reproductionRate, poolSize, tournamentProbability);

            mutation(Mchildren, mutationRate);
            mutation(Rchildren, mutationRate);

            updatePopulation(MainPop, Mchildren, selectedMain);
            updatePopulation(ReservePop, Rchildren, selectedReserve);

            bestFitNow = ReservePop.maxFit();
            ReservePop.addBest(bestFitNow);
            bestFitNow = MainPop.minFit();
            MainPop.addBest(bestFitNow);

            MainPop.setLastGeneration(iter+1);
            ReservePop.setLastGeneration(iter+1);
            if (bestFitNow<bestAccuracy)
            {
                lastIteration = iter+1;
                break;
            }

        }
        double [] results = MainPop.bestResult();
//        printVector(results);
        System.out.println("DPGA: Break in iteration number: "+ lastIteration + "  Best result: "+ bestFitNow);

    }// End of Standard DPGA
    /**
     * This method will perform the dual population GA
     */
    public void ImprovedDPGA(double survivalRate, double mutationRate, double EliteRate, int poolSize, double tournamentProbability, double bestAccuracy)
    {
        int eliteSize;
        if (EliteRate>=1)
            eliteSize = (int) EliteRate;
        else
            eliteSize = (int) (EliteRate*MainPop.popMatrix.length);

        double reproductionRate = 1- survivalRate;
        int selectedPopSize = (int)(survivalRate*populationSize/2);

        selectedPopSize = selectedPopSize - eliteSize;
        myMainPop.fillRandom();
        myReservePop.fillRandom();

        myselectedMain = new Population(selectedPopSize, Dimension, Xmin, Xmax, decimalAcc);
        int childSize= (int) (reproductionRate*populationSize/2);
        myMchildren = new Population(childSize, Dimension, Xmin, Xmax, decimalAcc);
        myselectedReserve = new Population(selectedPopSize, Dimension, Xmin, Xmax, decimalAcc);
        myRchildren = new Population(childSize, Dimension, Xmin, Xmax, decimalAcc);
//        double [] Results = new double [iteration];
        double reserveProb =0;
        double popDiversity=0;
        double bestFitNow = 100;
        int lastIteration = iteration;
        int compIndex;
        double myPopDiversity;
        double rand;
        int compTurn= (int) (iteration*0.1);
        int complNum=0;
        for (int iter=0; iter<iteration; iter++)
        {
            myMainPop.setZeroInsertingIndex();
            myselectedMain.setZeroInsertingIndex();
            complMainPop.complementPopulation(myMainPop);

            popDiversity =  myMainPop.populationDiversity();
            myMainPop.addDiversity();
            if (iter%compTurn == 0)
            {
                if (popDiversity < myMainPop.minPopDiversity())
                {
                    complNum = (int) (Math.log10(myMainPop.popMatrix[0].length)/Math.log10(myMainPop.popMatrix.length));
                    complNum++;
                    for (int count=0; count<complNum; count++)
                    {
                        compIndex= myRandomInt(myMainPop.popMatrix.length);
                        myselectedMain.insertIndividual(complMainPop.popMatrix[compIndex]);
                    }
                    myMainPop.emptyPopDivList();
                }
            }
            
            fitnessForAllMain(myMainPop);
//            myMainPop.addFitDiversity();
//            myMainPop.addDiversity();
            keepElite(myMainPop, myselectedMain, eliteSize, "min");



            myReservePop.setZeroInsertingIndex();
            myselectedReserve.setZeroInsertingIndex();
            fitnessForAllReserve(myMainPop, myReservePop);
            myReservePop.addDiversity();
            keepElite(myReservePop, myselectedReserve, eliteSize, "max");

            selectionDualPop(myMainPop, myselectedMain, myReservePop, myselectedReserve, survivalRate, poolSize, tournamentProbability);

            rand = myRandom()+1;// get a random number between 1 and 2
            reserveProb = rand*myMainPop.fitDiversity()/Math.abs(myMainPop.maxFit()-myMainPop.minFit()); // set the 2nd chance of reserve children to come in main

            crossoverAllDualPopImproved(reproductionRate, myMainPop, myMchildren, myReservePop, myRchildren, poolSize, tournamentProbability, reserveProb);
            

            mutation(myMchildren, mutationRate);
            mutation(myRchildren, mutationRate);

            updatePopulation(myMainPop, myMchildren, myselectedMain);
            complMainPop.complementPopulation(myMainPop);
            updatePopulation(myReservePop, myRchildren, myselectedReserve);

            bestFitNow = myReservePop.maxFit();
            myReservePop.addBest(bestFitNow);
            bestFitNow = myMainPop.minFit();
            myMainPop.addBest(bestFitNow);

            myMainPop.setLastGeneration(iter+1);
            myReservePop.setLastGeneration(iter+1);
            if (bestFitNow<bestAccuracy)
            {
                lastIteration = iter+1;
                break;
            }

        }
        double [] results = myMainPop.bestResult();
//        printVector(results);
        System.out.println("ImprovedDPGA: Break in iteration #: "+ lastIteration + "  Best result: "+ bestFitNow);

    }// End of Improved DPGA

    public void multipleRun(int run, double survivalRate, double mutationRate, double EliteRate, int poolSize, double tournamentProbability, double bestAccuracy)
    {
        double [] SGABestResult = new double [run];
        int [] SGABreak = new int [run];
        double [] DPGABestResult = new double [run];
        int [] DPGABreak = new int [run];
        double [] improvedDPGABestResult = new double [run];
        int [] improvedDPGABreak = new int [run];
        for (int r=0; r<run ; r++)
        {
            System.out.println(r+1);
            SGA(survivalRate, mutationRate, EliteRate, poolSize, tournamentProbability, bestAccuracy);
            SGABestResult[r] = basicPopulation.minFit();
            SGABreak[r]=basicPopulation.lastGeneration;
            basicPopulation.fillRandom();

            DPGA(survivalRate, mutationRate, EliteRate, poolSize, tournamentProbability, bestAccuracy);
            DPGABestResult[r] = MainPop.minFit();
            DPGABreak[r] = MainPop.lastGeneration;
            MainPop.fillRandom();
            ReservePop.fillRandom();

            ImprovedDPGA(survivalRate, mutationRate, EliteRate, poolSize, tournamentProbability, bestAccuracy);
            improvedDPGABestResult [r] = myMainPop.minFit();
            improvedDPGABreak[r] = myMainPop.lastGeneration;
            myMainPop.fillRandom();
            myReservePop.fillRandom();
        }
        double sgaFitAve = 0;
        int sgaBreakAve = 0;
        double dpgaFitAve = 0;
        int dpgaBreakAve = 0;
        double ImdpgaFitAve = 0;
        int ImdpgaBreakAve = 0;
        for (int r=0; r<run; r++)
        {
            sgaFitAve = sgaFitAve + SGABestResult[r];
            sgaBreakAve = sgaBreakAve + SGABreak[r];
            dpgaFitAve = dpgaFitAve + DPGABestResult[r];
            dpgaBreakAve = dpgaBreakAve + DPGABreak[r];
            ImdpgaFitAve = ImdpgaFitAve + improvedDPGABestResult[r];
            ImdpgaBreakAve = ImdpgaBreakAve + improvedDPGABreak[r];
        }
        sgaFitAve = sgaFitAve / run;
        sgaBreakAve = sgaBreakAve /run;
        dpgaFitAve = dpgaFitAve /run;
        dpgaBreakAve = dpgaBreakAve / run;
        ImdpgaFitAve = ImdpgaFitAve /run;
        ImdpgaBreakAve = ImdpgaBreakAve / run;
        System.out.println();
        System.out.println("Total average");
        double normal = - Math.log10(sgaFitAve)/Math.log10(sgaBreakAve);
        System.out.println("On average, SGA break on: " + sgaBreakAve+ " with the average best result: "+ sgaFitAve+ " NormalaizedRate: " + normal );
        normal = - Math.log10(dpgaFitAve)/Math.log10(dpgaBreakAve);
        System.out.println("On average, DPGA break on: " + dpgaBreakAve+ " with the average best result: "+ dpgaFitAve+ " NormalaizedRate: " + normal);
        normal = - Math.log10(ImdpgaFitAve)/Math.log10(ImdpgaBreakAve);
        System.out.println(" ImprovedDPGA break on: " + ImdpgaBreakAve+ " with the average best result: "+ ImdpgaFitAve+ " NormalaizedRate: " + normal);
    }






    public void printVector(double [] vector)
    {
        System.out.println();
        for (int i=0; i<vector.length; i++)
            System.out.println("  "+ vector[i]);
        System.out.println();
    }

    public void test()
    {
//        int [] temp = new int [] {1,0,0,1,0,0};
//        int [] temp2 = new int [] {1,1,0,1,1,1};
//        System.out.println(crossover2Point(temp, temp2));
//        basicPopulation.printPopulation();
    }



}
