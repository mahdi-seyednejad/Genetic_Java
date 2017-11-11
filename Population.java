/*
 *
 */

package geneticimprove;

import java.util.LinkedList;
import java.util.Random;

/**
 *
 */
public class Population {
    int [][] popMatrix;
    double [] fitness;
    int chromSize;
    int dimension;
    int eachDimensionLength;
    LinkedList<Double> linkOfBest;
    LinkedList<Double> listOfDiversity;
    LinkedList<Double> listOfFitDiversity;
    LinkedList<Double> listOfAverageFit;
    int selectedIndex;
    int insertingIndex;
    int lastGeneration;
    double fitFluctuation;


    /**
     * This is the main constructor.
     * It takes the following parameters as input and construct the population based on them
     * size of population, number of dimensions, Min, Max, accuracy of representation
     */
    public Population(int populationSize, int Dimension, double min, double max, int decimalAccuracy)
    {
        double bitNum = Math.abs(max-min)*decimalAccuracy;
        bitNum = log2(bitNum);
        eachDimensionLength = (int)bitNum +1;
        dimension = Dimension;
        chromSize = dimension*eachDimensionLength;
        popMatrix = new int [populationSize][chromSize];
        fitness = new double [populationSize];
        fillRandom();
        linkOfBest = new LinkedList<Double>();
        listOfDiversity = new LinkedList<Double>();
        listOfFitDiversity = new LinkedList<Double>();
        listOfAverageFit = new LinkedList<Double>();
        insertingIndex =0;
        lastGeneration = 0;


    } // End of the constructor


    /**
     * This method calculates Log2(x)
     */
    public double log2(double x)
    {
        return Math.log(x)/Math.log(2);
    } // End of Log2

    /**
     * This method will fill each individual in the population randomly
     */
    public void fillRandom()
    {
        for (int i=0; i<popMatrix.length; i++)
            for (int j=0; j<popMatrix[0].length; j++)
                popMatrix[i][j] = myRandomInt(2);
    } // End of fill randomly

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
     * This method will sort the population base on their fitness
     */
    public void sortPopulation_Fitness()
    {
        double [] tempFit = new double [fitness.length];
        for (int i=0; i<tempFit.length; i++)
            tempFit[i] = fitness[i];
        double temp;
        for (int i=0; i<fitness.length; i++)
            for (int k=i+1; k<fitness.length; k++)
                if (fitness[i]<fitness[k])
                {
                    exChangeOnPopulation(i, k);
                    temp = fitness[i];
                    fitness[i] = fitness[k];
                    fitness[k] = temp;
                }

    }// End of sort population

    /**
     * This method exchages two individuals which their indeces are determined as input
     */
    public void exChangeOnPopulation(int index1, int index2)
    {
        int [] temp= new int [popMatrix[0].length];
        for (int j=0; j<temp.length; j++)
            temp[j] = popMatrix[index1][j];
        for (int j=0; j<temp.length; j++)
            popMatrix[index1][j] = popMatrix[index2][j];
        for (int j=0; j<temp.length; j++)
            popMatrix[index2][j] = temp[j];
    }// End of exchange on population

    /**
     * This method will replace the population by a matrix in its input
     */
    public void replacePopulation( int [][] Matrix)
    {
        for (int i=1; i<Matrix.length; i++)
            for (int j=0; j<Matrix[0].length; j++)
                popMatrix[i][j] = Matrix[i][j];
    }// end of replacement

    /**
     * This method performs crossover and returns one parent as one vector.
     * It takes the size of the pool of tournamentMin. (if size=2 then it works as a binary tournamentMin)
     *
     */
    public int [] tournamentMin(int poolSize)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        double minFit=fitness[pool[0]];
        int minIndex=pool[0];
        for (int i=1; i<pool.length; i++)
            if (fitness[pool[i]]<minFit)
            {
                minIndex = pool[i];
                minFit = fitness[pool[i]];
            }
        selectedIndex = minIndex;

        return popMatrix[minIndex];
    }// End of tournamentMin
    /**
     * This method performs crossover and returns one parent as one vector.
     * It takes the size of the pool of tournamentMin. (if size=2 then it works as a binary tournamentMin)
     * if poolsize = 2 then it returns the better individual with the probability of "prob" and the other with 1-prob
     */
    public int [] tournamentMin(int poolSize , double prob)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        int temp;

        for (int i=0; i<poolSize; i++)
        {
            for (int j=i+1; j<poolSize; j++)
                if (fitness[pool[i]]>fitness[pool[j]])
                {
                    temp = pool[i];
                    pool[i] = pool[j];
                    pool[j] = temp;
                }

        }// now we have sorted pool
        double rand2=myRandom();
//        for (int i=0; i<pool.length; i++)


        double minFit=fitness[pool[0]];
        int minIndex=pool[0];
        for (int i=1; i<pool.length; i++)
            if (fitness[pool[i]]<minFit)
            {
                minIndex = pool[i];
                minFit = fitness[pool[i]];
            }
        
        if (poolSize == 2)
        {
            rand2 = myRandom();
            if (rand2>prob)
                if (minIndex == pool[0])
                    minIndex = pool[1];
                else
                    minIndex = pool[0];
        }
        selectedIndex = minIndex;
        return popMatrix[minIndex];
    }// End of tournamentMin
    /**
     * This method performs crossover and returns the index of parent
     * It takes the size of the pool of tournamentMin. (if size=2 then it works as a binary tournamentMin)
     * if poolsize = 2 then it returns the better individual with the probability of "prob" and the other with 1-prob
     */
    public int tournamentMinIndex(int poolSize , double prob)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        double minFit=fitness[pool[0]];
        int minIndex=pool[0];
        for (int i=1; i<pool.length; i++)
            if (fitness[pool[i]]<minFit)
            {
                minIndex = pool[i];
                minFit = fitness[pool[i]];
            }
        double rand2=0;
        if (poolSize == 2)
        {
            rand2 = myRandom();
            if (rand2>prob)
                if (minIndex == pool[0])
                    minIndex = pool[1];
                else
                    minIndex = pool[0];
        }


        return minIndex;
    }// End of tournamentMin
    /**
     * This method performs crossover and returns one parent as one vector.
     * It takes the size of the pool of tournamentMin. (if size=2 then it works as a binary tournamentMin)
     * if poolsize = 2 then it returns the better individual with the probability of "prob" and the other with 1-prob
     */
    public int [] tournamentMax(int poolSize , double prob)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        double maxFit=fitness[pool[0]];
        int maxIndex=pool[0];
        for (int i=1; i<pool.length; i++)
            if (fitness[pool[i]]>maxFit)
            {
                maxIndex = pool[i];
                maxFit = fitness[pool[i]];
            }
        double rand2=0;
        if (poolSize == 2)
        {
            rand2 = myRandom();
            if (rand2>prob)
                if (maxIndex == pool[0])
                    maxIndex = pool[1];
                else
                    maxIndex = pool[0];
        }
        selectedIndex = maxIndex;


        return popMatrix[maxIndex];
    }// End of tournamentMin
    /**
     * This method performs crossover and returns one parent as one vector.
     * It takes the size of the pool of tournamentMin. (if size=2 then it works as a binary tournamentMin)
     * if poolsize = 2 then it returns the better individual with the probability of "prob" and the other with 1-prob
     */
    public int tournamentMaxIndex(int poolSize , double prob)
    {
        int [] pool = new int [poolSize];
        int [] chosen = new int [popMatrix.length];
        for (int i=0; i<chosen.length; i++)
            chosen[i]=0;
        int rand;
        int countPool=0;
        while (countPool < poolSize)
        {
            rand = myRandomInt(popMatrix.length);
            if (chosen[rand] == 0)
            {
                pool[countPool] = rand;
                chosen [rand] = 1;
                countPool++;
            }
        }
        double maxFit=fitness[pool[0]];
        int maxIndex=pool[0];
        for (int i=1; i<pool.length; i++)
            if (fitness[pool[i]]>maxFit)
            {
                maxIndex = pool[i];
                maxFit = fitness[pool[i]];
            }
        double rand2=0;
        if (poolSize == 2)
        {
            rand2 = myRandom();
            if (rand2>prob)
                if (maxIndex == pool[0])
                    maxIndex = pool[1];
                else
                    maxIndex = pool[0];
        }

        return maxIndex;
    }// End of tournamentMin


    public int RolletteWheel ()
    {
        double sum =0;
        double partSum=0;
        for (int i=0; i<fitness.length; i++)
            sum = sum + fitness[i];
        double rand = myRandom();
        rand = rand * sum;
        int dataIndex =-1;
        do
        {
            dataIndex ++;
            if (dataIndex == 60)
                System.out.println(dataIndex);
            partSum = partSum + fitness[dataIndex];
        }while (partSum<=rand);
//        dataIndex++;
        return dataIndex;
    } // End of Rollette Wheel
    /**
     * This method will find the minimum fitness and returns its value
     */
    public double minFit()
    {
        double min;
        min = fitness[0];
        for (int i=1; i<fitness.length; i++)
            if (fitness[i]<min)
                min = fitness[i];
        return min;
    }// End minFit
    /**
     * This method will find the minimum fitness and returns its value
     */
    public int minFitIndex()
    {
        double min;
        min = fitness[0];
        int minIndex=0;
        for (int i=1; i<fitness.length; i++)
            if (fitness[i]<min)
            {
                min = fitness[i];
                minIndex = i;
            }
        return minIndex;
    }// End minFitIndex
    /**
     * This method will find the maximum fitness and returns its value
     */
    public double maxFit()
    {
        double max;
        max = fitness[0];
        for (int i=1; i<fitness.length; i++)
            if (fitness[i]>max)
                max = fitness[i];
        return max;
    }// End maxFit
    /**
     * This method will find the maximum fitness and returns its index
     */
    public int maxFitIndex()
    {
        double max;
        max = fitness[0];
        int maxIndex=0;
        for (int i=1; i<fitness.length; i++)
            if (fitness[i]> max)
            {
                max = fitness[i];
                maxIndex = i;
            }
        return maxIndex;
    }// End maxFitIndex

    /**
     * This method adds the best solution found as its input
     */
    public void addBest(double fit)
    {
        linkOfBest.addLast(fit);
    }// End of add best
    /**
     * This method will returns the best results which were inserted in this population
     */
    public double [] bestResult()
    {
        int listSize = linkOfBest.size();
        double [] tempResult = new double [listSize];
        double tempD=0;
        for (int i=0; i<listSize; i++)
        {
            tempD = linkOfBest.removeLast();
            tempResult[i] = tempD;
            linkOfBest.addFirst(tempD);
        }

        return tempResult;
    }// end of best results
    /**
     * This method inserts a vector in the population
     */
    public void insertIndividual(int [] vector)
    {
        if (insertingIndex<popMatrix.length)
        {
            for (int j=0; j<vector.length; j++)
                popMatrix[insertingIndex][j] = vector[j];
            insertingIndex++;
        }

    }// End of inserting

    public void setZeroInsertingIndex()
    {
        insertingIndex =0;
    }// End of set zero inserting index
    /**
     * This method will save the maximum number of generations that this population pass through
     */
    public void setLastGeneration(int generationNum)
    {
        lastGeneration = generationNum;
    }
    /**
     * This method will calculate this population diversity
     */
    public double populationDiversity()
    {
        double div=0;
        double [] freqVector = new double [popMatrix[0].length];
        for (int j=0; j<popMatrix[0].length; j++)
            for (int i=0; i<popMatrix.length; i++)
                freqVector[j] = freqVector[j] + popMatrix[i][j];
        for (int i=0; i<freqVector.length; i++)
            freqVector[i] = freqVector[i]/popMatrix.length;

        for (int i=0; i<freqVector.length; i++)
            div = div + (freqVector[i]*(1-freqVector[i]));
        div = div * Math.pow(popMatrix.length, 2)/popMatrix[0].length;

        return div;
    }// End of population diversity
    /**
     * This method adds new diversity to this population
     */
    public void addDiversity()
    {
        double diversity = populationDiversity();
        listOfDiversity.addFirst(diversity);
    }// End of add diversity to its list
    /**
     * This method will return this population diversity up to now
     */
    public double [] diversityVector()
    {
        int listSize = listOfDiversity.size();
        double [] tempResult = new double [listSize];
        double tempD=0;
        for (int i=0; i<listSize; i++)
        {
            tempD = listOfDiversity.removeLast();
            tempResult[i] = tempD;
            listOfDiversity.addFirst(tempD);
        }

        return tempResult;
    }// end of best results
    /**
     * This method will replace the complement of the inpup population
     */
    public void complementPopulation(Population primaryPop)
    {
        for (int i=0; i<popMatrix.length; i++)
            for (int j=0; j<popMatrix[0].length; j++)
                popMatrix[i][j] = Complement(primaryPop.popMatrix[i][j]);

    }// End of complementPopulation

    private int Complement(int a)
    {
        if (a==0)
            return 1;
        else
            return 0;
    }// End of complement
    /**
     * This method calculates the diversity of fitness in this population
     */
    public double fitDiversity()
    {
        double averageFit=0;
        double divFit=0;
        for (int i=0; i<fitness.length; i++)
            averageFit = averageFit + fitness[i];
        averageFit = averageFit / fitness.length;
        for (int i=0; i<fitness.length; i++)
            divFit = divFit + Math.pow(fitness[i] - averageFit, 2);

        divFit = Math.pow(divFit, 0.5)/fitness.length;

        return divFit;
    } // End of fitDiversity

    
    /**
     * This method adds current fitness diversity to its list
     */
    public void addFitDiversity()
    {
        listOfFitDiversity.addFirst(fitDiversity());
    }// addFitDiversity
    /**
     * This method will return this population diversity up to now
     */
    public double [] fitnessDiversityVector()
    {
        int listSize = listOfFitDiversity.size();
        double [] tempResult = new double [listSize];
        double tempD=0;
        for (int i=0; i<listSize; i++)
        {
            tempD = listOfFitDiversity.removeLast();
            tempResult[i] = tempD;
            listOfFitDiversity.addFirst(tempD);
        }

        return tempResult;
    }// end of best results
    /**
     * This method adds the curent average fitness to its list
     */
    public void setAverageFit()
    {
        double aveFit=0;
        for (int i=0; i<fitness.length; i++)
            aveFit = aveFit + fitness[i];
        aveFit = aveFit / fitness.length;
        listOfAverageFit.addFirst(aveFit);
    }// End of setAverageFit
    /**
     * This method computes the fitness fluctuation from generation G up to now
     */
    public void computeFitFluctuation(int startGeneration)
    {
        int listSize = listOfAverageFit.size();
        double [] tempResult = new double [listSize];
        double tempD=0;
        for (int i=0; i<listSize; i++)
        {
            tempD = listOfAverageFit.removeLast();
            tempResult[i] = tempD;
            listOfAverageFit.addFirst(tempD);
        }
        double aveFit=0;
        if (startGeneration<listSize)
        {
            for (int i=startGeneration; i<listSize; i++)
                aveFit = aveFit + tempResult[i];
            aveFit = aveFit/(Math.abs(startGeneration-listSize));
        }

        fitFluctuation = aveFit;
    }// End of compute Fitness Fluctuation

    public double minPopDiversity()
    {
        double [] tempDiv = diversityVector();
        return findMinValue(tempDiv);

    }

    public double findMaxValue(double [] vector)
    {
        double max = vector[0];
        for (int i=0; i<vector.length; i++)
            if (vector[i]>max)
                max = vector[i];
        return max;
    }
    public double findMinValue(double [] vector)
    {
        double min = vector[0];
        for (int i=0; i<vector.length; i++)
            if (vector[i]<min)
                min = vector[i];
        return min;
    }

    public void emptyPopDivList()
    {
        int size = listOfDiversity.size();
        for (int i=0; i<size ; i++)
            listOfDiversity.removeFirst();
    }






    public void printPopulation()
    {
        System.out.println();
        for (int i=0; i<popMatrix.length; i++)
        {
            System.out.print(i + " : " );
            for (int j=0;j<popMatrix[0].length; j++)
                System.out.print(" \t "+popMatrix[i][j]);
            System.out.println();
        }

    }








}
