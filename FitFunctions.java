/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package geneticimprove;

/**
 *
 * @author mehdi
 */
public class FitFunctions {
    int Dimension;
    int geneSize;
    double Xmin;
    double Xmax;

    public FitFunctions() {
    }
    public FitFunctions(int dimension, int eachGeneSize, double min, double max )
    {
        Dimension = dimension;
        geneSize = eachGeneSize;
        Xmin = min;
        Xmax = max;
    }


    /**
     * This method calculates Sphere function for an individual as input
     */
    public double SphereFit(int [] individual)
    {
        double sphereValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        sphereValue = Sphere(realTemp);

        return sphereValue;
    }// End of ackley fitness

    private double Sphere(double [] X)
    {
        double sum=0;
        for (int i=0; i<X.length; i++)
            sum = sum + Math.pow(X[i], 2);
//        if (sum<Math.pow(10, -18))
//            sum=0;
        return sum;
    }

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
    }// end of binari to Decimal

    /**
     * This method will convert the binary input to a real value in [min, max]
     */
    public double bin2Real(int [] binNumber)
    {
        double myX = Xmin + (bin2Decimal(binNumber)*(Xmax-Xmin)/(Math.pow(2, binNumber.length)-1));

        return myX;

    }// End of bin to real

    /**
     * This method calculates Ackley function for an individual as input
     */
    public double AckleyFit(int [] individual)
    {
        double achleyValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        achleyValue = Ackley(realTemp);

        return achleyValue;
    }// End of ackley fitness

    private double Ackley(double [] X)
    {
        double ackley=-1;
        int D = X.length;
        double sumCos = 0;
        double sumSqr=0;
        for (int i=0; i<X.length; i++)
            sumSqr = sumSqr + Math.pow(X[i], 2);
        for (int i=0; i<X.length; i++)
            sumCos = sumCos + Math.cos(2*Math.PI*X[i]);

//        ackley = 20 + Math.E - (20 * Math.exp(-0.2*Math.sqrt(sumSqr/D))) - Math.exp(sumCos/D);
        ackley = 20 + Math.pow(Math.E,1) - (20 * Math.pow(Math.E,(-0.2*Math.sqrt(sumSqr/D)))) - Math.pow(Math.E,(sumCos/D));
//        if (ackley == -4.440892098500626* Math.pow(10, -16))
//            ackley =0;
        if (ackley < Math.pow(10, -14))
            ackley =0;

        return ackley;
    } // End of Ackley on double

    /**
     * This method calculates Rastrigin function for an individual as input
     */
    public double RastriginFit(int [] individual)
    {
        double rastValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        rastValue = Rastrigin(realTemp);

        return rastValue;
    }// End of Rastrigin fitness

    /**
     * This method caculates the output of Rastrigin function on real values
     */
    private double Rastrigin(double [] X)
    {
        double rastrigin = 0;
        int D = X.length;
        double sum=0;
        double Cos=0;
        for (int i=0; i<D; i++)
        {
            Cos = Math.cos(2*Math.PI*X[i]);
            if (Math.abs(Cos)<Math.pow(10, -12))
                Cos = 0;
            sum = sum + (Math.pow(X[i], 2) - (10*Cos));
        }
            
        rastrigin = (10*D) + sum;

        return rastrigin;
    }// End of Rastrigin on double
    /**
     * This method calculates Griewank function for an individual as input
     */
    public double GriewankFit(int [] individual)
    {
        double GrieValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        GrieValue = Griewank(realTemp);

        return GrieValue;
    }// End of Griewank fitness

    /**
     * This method calculates Griewank function on real values
     */
    private double Griewank(double [] X)
    {
        double griewank = 0;
        double sum =0;
        double multiple=1;
        for (int i=0; i<X.length; i++)
            sum = sum + Math.pow(X[i], 2);
        sum = sum / 4000;
        for (int i=0; i<X.length; i++)
            multiple = multiple * Math.cos(       X[i]/(  Math.sqrt(i+1)  )    );

        griewank = 1 + sum - multiple;

        return griewank;
    }// End of Griewank on real values

    /**
     * This method calculates Schwefel function for an individual as input
     */
    public double SchwefelFit(int [] individual)
    {
        double SchewValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        SchewValue = Schwefel(realTemp);

        return SchewValue;
    }// End of Schwefel fitness

    /**
     * Schwefel on real values
     */
    private double Schwefel(double [] X)
    {
        double schwefel =0;
        double sum=0;
        for (int i=0; i<X.length; i++)
            sum = sum + (418.9829*(i+1)) - (X[i] * Math.sin(Math.sqrt(Math.abs(X[i]))) );
//            sum = sum  - (X[i] * Math.sin(Math.sqrt(Math.abs(X[i]))) );

        schwefel = sum;

        return Math.abs(schwefel) ;
    } // End of Schwefel on real values


    /**
     * This method calculates RosenBrock function for an individual as input
     */
    public double RosenbrockFit(int [] individual)
    {
        double RosenValue = -1;
        double [] realTemp = new double [Dimension];
        int realIndex=0;
        int [] binTemp = new int [geneSize];
        for (int d=0; d<Dimension; d++)
        {
            int start = d*geneSize;
            for (int i= 0; i<geneSize; i++)
                binTemp[i] = individual[i+start];
            realTemp[realIndex] = bin2Real(binTemp);
            realIndex++;
        }
        RosenValue = Rosenbrock(realTemp);

        return RosenValue;
    }// End of RosenBrock fitness
    /**
     * RosenBrock on real value
     */
    private double Rosenbrock (double [] X)
    {
        double rosenB =0;
        double temp;
        double sum1=0;
        double sum2=0;
        for (int i=0; i<X.length-1; i++)
        {
            temp = X[i+1] - Math.pow(X[i], 2);
            sum1 = sum1 + (100 * Math.pow(temp, 2));
        }
        for (int i=0; i<X.length; i++)
            sum2 = sum2 + Math.pow(1-X[i], 2);
        rosenB = sum1 + sum2;

        return rosenB;
    }// End of Rosenbrock on real values

    /**
     * This method calculates diversity
     * It is given a matrix and a vector as input and returns the normalaized hamming distance over
     * all rows of the matrix
     */
    public double diversityHamming(int [][] Matrix, int [] Vector)
    {
        int count =0;
        for (int i=0; i<Matrix.length; i++)
            count = count + HammingDistance(Matrix[i], Vector);

        double div = count;
        div = div/ Matrix.length;

        return div;
    }// End of diversity based on hamming distance
    /**
     * This method gives to vector as input and returns the hamming distance between them.
     */
    public int HammingDistance(int [] V1 , int [] V2)
    {
        int count =0;
        for (int i=0; i<V1.length; i++)
            if (V1[i] != V2[i])
                count++;
        return count;
    }// End of hamming distance

    /**
     * This method calculates diversity
     * It is given a matrix and a vector as input and returns the normalaized hamming distance over
     * all rows of the matrix
     */
    public double diversityByFrequency(int [][] Matrix, int [] Vector)
    {
        double count =0;
        for (int i=0; i<Matrix.length; i++)
            count = count + FrequencyOnEachIndividual(Matrix[i], Vector);

        double div = count;
        div = div/ Matrix.length;

        return div;
    }// End of diversity based on frequency


    /**
     * This metod will compute the frequency differences between two vectors
     */
    public double FrequencyOnEachIndividual(int [] V1, int [] V2)
    {
        double freq1=0;
        double freq2=0;
        for (int i=0; i<V1.length; i++)
        {
            if (V1[i]==1)
                freq1++;
            if (V2[i]==1)
                freq2++;
        }
        freq1 = freq1 / V1.length;
        freq2 = freq2 / V2.length;

        return Math.abs(freq1-freq2);
    }// End of Frequency Difference
    /**
     * This method will compute diversity base on frequency on each gene
     */
    public double diversityFrequencyOnGene(int [][] Matrix, int [] vector)
    {
        double sumFreq=0;
        double sumOnGene=0;
        double [] freqV = new double [vector.length];
        for (int j=0; j<Matrix[0].length; j++ )
        {
            sumOnGene=0;
            for (int i=0; i<Matrix.length; i++)
                sumOnGene = sumOnGene + Matrix[i][j];
            freqV[j] = sumOnGene;
        }
        for (int i=0; i<freqV.length; i++)
            freqV[i] = freqV[i]/Matrix.length;
        for (int i=0; i<freqV.length; i++)
            sumFreq = sumFreq + Math.abs(freqV[i]-vector[i]);

        return sumFreq;
    }// End of diversityFrequencyOnGene
    /**
     * This mathod returns true if the two vecrors in the input are similar
     */
    public boolean isEqualVector(int [] vector1, int [] vector2)
    {
        for (int i=0; i<vector1.length; i++)
            if (vector1[i] != vector2[i])
                return false;
        return true;
    }// end of is equal 2 vectors

    public double appropriateFitness(int [] vector, int functionNumber)
    {
        double result=100;
        switch (functionNumber)
            {
                case 1: result = SphereFit(vector);
                        break;
                case 2: result = AckleyFit(vector);
                        break;
                case 3: result = RastriginFit(vector);
                        break;
                case 4: result = GriewankFit(vector);
                        break;
                case 5: result = SchwefelFit(vector);
                        break;
                case 6: result = RosenbrockFit(vector);
                        break;
            }
        return result;
    }






    public void test()
    {
        int [] temp = new int [] {1,0,0,1 , 1,0,0,1};
        int [] temp2 = new int [] {1,0,1,1 , 1,1,0,1};
        int [] temp3 = new int [] {0,0,0,0 , 1,1,0,1};
        int [][] Mx = new int [2][temp.length];
        for (int i=0; i<Mx[0].length; i++)
            Mx[0][i] = temp[i];
        for (int i=0; i<Mx[0].length; i++)
            Mx[1][i] = temp2[i];
        double [] Xtemp = new double [] {1,1,1,1};//Math.pow(10, -10)};
        System.out.println(diversityFrequencyOnGene(Mx, temp3));

//        System.out.println(Math.cos(       1/(  Math.sqrt(4)  )    ));
//        System.out.println(Math.cos(0.5));
//        System.out.println("e: "+Math.E + "\texp(1): " + Math.exp(1)+ "\te^1: " + Math.pow(Math.E, 1));
    }




}

