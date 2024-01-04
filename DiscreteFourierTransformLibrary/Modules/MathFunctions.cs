using System;
using DiscreteFourierTransformLibrary.Models;
using System.Collections.Generic;

namespace DiscreteFourierTransformLibrary.Modules
{
    public class MathFunctions
    {
        /*
      * This module is created as a means to 
      * do discrete fourier Transforms on data
      *
      * The purpose of this library is to help with 
      * design of Fourier analysis for RoundTable's purposes
      *
      * Either to create a relatively quick correlation 
      * with the expectation of similar output for same input.
      * we can extract an optimized means for variable PID's 
      * to reach set points.
      *
      * Could also be used to find fault detection in devices
      * not behaving correctly
      *
      */

        public ModularElement GeneralRealMod(double input, double modNumber)
        {
            /*
            * SYNOPSIS: 
            * This functions performs
            * modular arithmetic with input
            * over modNumber
            * 
            *
            * DESCRIPTION:
            * 
            * input = modNumber*N + r
            * where r is in the range 0<=r<modNumber
            * and N is an integer.
            * 
            * The output of this function
            * contains all of the information in the
            * equation above
            * 
            * OUTPUTS:
            * ModularElement.OriginalInput = input
            * ModularElement.ModularDivisor = modNumber
            * ModularElement.IntegerMultiple = N
            * ModularElement.Remainder = r
            *
            *
            * DEPENDENCIES:
            *
            * 
            *
            * PARAMETER:
            *
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            ModularElement retVal = new ModularElement();
            retVal.OriginalInput = input;
            retVal.ModularDivisor = modNumber;
            double remainder = input;
            double intMult = 0;
            if (input >= 0)
            {
                

                while (!(remainder < modNumber))
                {
                    intMult += 1;

                    remainder = remainder - modNumber;


                }

                
            }
            else
            {
                while(remainder < 0)
                {
                    intMult += -1;

                    remainder = remainder + modNumber;

                }

                
            }
            retVal.Remainder = remainder;
            retVal.ModularDivisor = intMult;

            return retVal;


        }

        public List<int> PrimeNumbersLessthan(int MyNum)
        {
            /*
            * SYNOPSIS: 
            * This lists the prime numbers less
            * than the number given. 
            * 
            *
            * DESCRIPTION:
            * This finds all prime that are
            * less than or equal to MyNum. 
            * 
            * if the number is a composite then
            * the smallest possible factor is 2
            * while the largest is the number (j/2)
            *
            * This algorithm tests assumes the number
            * is prime until there is a factor 
            * between 1 and j/2 (i) that is a divisor.
            * the comparison used to determine if it 
            * is a divisor is when the number j % i == 0 
            * by definition of modular classes in the 
            * Domain Integral of integers.
            * 
            * this tests every number less than the input
            * MyNum to be prime.
            *
            * the order is O(MyNum^2)
            * 
            *
            *
            * DEPENDENCIES:
            *
            * OUTPUTS:
            * List<int>
            *
            * PARAMETER:
            *
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts
            var primes = new List<int>();

            for (int j = 2; j < MyNum + 1; j++)
            {

                bool isItPrime = true;
                for (int i = 2; i < (j / 2 + 1); i++)
                {
                    if (j % i == 0)
                    {
                        isItPrime = false;

                        break;
                    }
                }
                if (isItPrime)
                {
                    primes.Add(j);
                }



            }

            return primes;


        }
        public List<PrimeFactor> FindPrimeFactors(int x)
        {
            /*
            * SYNOPSIS: 
            * THis finds all prime factors
            * of a given number x
            *
            * DESCRIPTION:
            * if a number is a composite
            * then x = yz where y and z are not unitary
            * one number must be less than or equal to sqrt(x) "y" and 
            * the other is greater than or equal to sqrt(x) "z"
            * 
            * since z could be prime so we search among 
            * the integer that is less than 1 
            * 
            *
            * DEPENDENCIES:
            * PrimeNumbersLessthan(int MyNum)
            * 
            *
            * OUTPUTS:
            *
            * PARAMETER:
            * x; the number that you want to find the prime
            * factors of
            *
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            var doub = Convert.ToDouble(x);
            var primeMax = Math.Sqrt(doub);
            var primeMaxFloor = Math.Floor(primeMax);
            double primeMaxDoub;
            /* 
             * if the floor is equal to the sqrt of the integer
             * than the integer x = y^2 where y is an integer
             * if that is true then the primes of x are the same as
             * the primes of y. the difference between two numbders
             * is that the unique factorization of x has twice 
             * the power per prime as y.
             * 
             * The while loop below continually finds
             * the continually dwindles the number until
             * a point where it is not square rootable
             * 
             * so it finds the 1/2,1/4,1/8,1/16...
             * if possible. rootPower is the inverse
             * of that number
             * 
             * so x = doub^rootPower
             * thus we search for the primes of doub
             * to determine the primes of x.
             * 
             * doing so we initiate our coFactor = doub
             *
             */
            int rootPower = 1;
            while (primeMaxFloor == primeMax)
            {

                doub = Math.Sqrt(doub);
                primeMax = Math.Sqrt(doub);
                primeMaxFloor = Math.Floor(primeMax);
                rootPower = rootPower * 2;
            }

            primeMaxDoub = primeMaxFloor + 1.0;

            var primeMaxInt = Convert.ToInt32(primeMaxDoub);

            List<int> possiblePrimes = PrimeNumbersLessthan(primeMaxInt);


            /* 
             * What this loop does is that it finds
             * the prime factors of doub and its powers.
             * 
             * if doub = g*coFactor where both numbers are 
             * then the loop finds all prime factors 
             * of doub that are less than <= sqrt(doub)
             * and there maximum power.
             * 
             * as it iterates through, g absorbs all of the primes
             * and their powers causing the coFactor to shrink 
             * to a point where it is relatively prime to all
             * prime numbers <= sqrt(doub). 
             * 
             * This would imply if the coFactor is a composite
             * then it must be a composite of primes > sqrt(doub)
             * 
             * however since the coFactor <= doub then at least
             * one of its primes factors must be <=sqrt(doub) 
             * which is a contradiction. therefore the remaining 
             * coFactor is infact a prime number or 1 which means that 
             * prime number's maximal power is 1 for doub.
             * 
             * since we are interested in the maximal power of
             * a prime factor of x and that x = doub^rootPower
             * then the maximal power of x is rootPower times
             * the maximal power of doub.
             * 
             */
            List<PrimeFactor> primeFactors = new List<PrimeFactor>();
            int coFactor = Convert.ToInt32(doub);
            foreach (int prime in possiblePrimes)
            {
                if (coFactor % prime == 0)
                {
                    PrimeFactor primeFactor = new PrimeFactor();
                    primeFactor.Numb = prime;
                    int m = 1;
                    while (coFactor % Math.Pow(prime, m) == 0)
                    {
                        m++;
                    }
                    primeFactor.Pow = (m - 1) * rootPower;
                    primeFactors.Add(primeFactor);
                    double coFactDoub =
                    Convert.ToDouble(coFactor) / (Math.Pow(prime, m - 1));
                    coFactor = Convert.ToInt32(coFactDoub);

                }
            }


            if (coFactor != 1)
            {
                PrimeFactor primeFactor = new PrimeFactor();
                primeFactor.Numb = coFactor;
                primeFactor.Pow = rootPower;
                primeFactors.Add(primeFactor);

            }

            return primeFactors;
        }

        public List<PrimeFactorIndex> CompositeListGenerator(List<PrimeFactor> factors)
        {
            /*
            * SYNOPSIS: 
            * This generates an indexed list of prime factors 
            * for a composite number N
            * 
            * DESCRIPTION:
            * N = r0r1r2...rl = p1^m1p2^m2p3^m3...
            * so r0 to rm1-1 are all p1
            * rm1 to rm1_m2 are al p2 and so on
            * in the most reduced decomposition of N. 
            *
            * this method will list out each prime factor by its
            * where it lies. This way when we do our
            * Discrete Fourier Transform (DFT), we will organize
            * the prime numbers in a way that will allow us to 
            * define each subtransformation to calculate the 
            * DFT using the Cooley-Tukey
            * Fast Fourier Transform (FFT) methodology
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            Int32 index = 0;
            Int32 totalNumber = 1;
            var returnList = new List<PrimeFactorIndex>();
            foreach(PrimeFactor factor in factors)
            {
                    totalNumber =
                    totalNumber*Convert.ToInt32(Math.Pow(
                    Convert.ToDouble(factor.Numb),
                    Convert.ToDouble(factor.Pow)
                    ));
            }

            Int32 subComp = 1;
            foreach(PrimeFactor factor in factors)
            {
                for(Int32 subIndex=1; subIndex <= factor.Pow; subIndex++)
                {
                    subComp = subComp * factor.Numb;
                    var listFactor = new PrimeFactorIndex {
                        Prime = factor,
                        Index = index,
                        SubIndex = subIndex,
                        TotalComposite = totalNumber,
                        SubComposite = subComp,
                        LeftOverComposite = totalNumber/subComp
                    };
                    index += 1;

                    returnList.Add(listFactor);

                }
            }
            return returnList;

        }


        public double[] SortArray(double[] array, int leftIndex, int rightIndex)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This is quicksort for doubles
            *
            * Algorithm from 
            * https://code-maze.com/csharp-quicksort-algorithm/
            * 
            *
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            var i = leftIndex;
            var j = rightIndex;
            var pivot = array[leftIndex];

            while (i <= j)
            {
                while (array[i] < pivot)
                {
                    i++;
                }

                while (array[j] > pivot)
                {
                    j--;
                }

                if (i <= j)
                {
                    double temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                    i++;
                    j--;
                }
            }

            if (leftIndex < j)
                SortArray(array, leftIndex, j);

            if (i < rightIndex)
                SortArray(array, i, rightIndex);

            return array;
        }

        public int[] SortArrayInt(int[] array, int leftIndex, int rightIndex)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This is quicksort for ints
            *
            * Algorithm from 
            * https://code-maze.com/csharp-quicksort-algorithm/
            * 
            *
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            var i = leftIndex;
            var j = rightIndex;
            var pivot = array[leftIndex];

            while (i <= j)
            {
                while (array[i] < pivot)
                {
                    i++;
                }

                while (array[j] > pivot)
                {
                    j--;
                }

                if (i <= j)
                {
                    int temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                    i++;
                    j--;
                }
            }

            if (leftIndex < j)
                SortArrayInt(array, leftIndex, j);

            if (i < rightIndex)
                SortArrayInt(array, i, rightIndex);

            return array;
        }

        /*&&&&&&&&&&&&&&&&& FFT LINE &&&&&&&&&&&&&&&&&&
         * Below are the functions that do the standard 
         * Discrete Fourier Transform functions.
         */

        public System.Numerics.Complex[,] OneStepFFT(
            System.Numerics.Complex[,] preFFT,
            List<PrimeFactorIndex> primes,
            System.Numerics.Complex W,
            int stepi
            )
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This is a helper function
            * to perform FFT transformations
            * for generic N decomposed into the 
            * a List<PrimeFactorIndex>
            * from the method CompositeListGenerator
            * 
            * The Cooley-Tukey Method can be generailized
            * for N = r1r2...ri-1riri+1...rl
            * where {ri} are the primes composing N.
            *
            * This will be the most reduced form of N to
            * perform an FFT. The method that will be 
            * utilized by this program will be doing 
            * steps 3 through 7 iteratively to eventually
            * get to X. 
            *
            * In the original paper they only considered 
            * l = 2.
            * where A1 is calculated from the original
            * function A which we can call A0 in this calculation
            * They caculate A1 by breaking up k into two components
            * and summing over one component of k to get the 
            * A1 as a function of j and k
            * 0
            * The step that this function is doing is this calculation
            * below
            * Ai+1(k,ji,ji-1,...,j0)
            * = Sum{Ai(k'*N/(r1r2...ri)+k,ji-1,ji-2,...,j0)
            * *W^[k'[N/(r1r2...ri)](ji(r1r2..ri)+ji-1(r1r2...ri-1)
            * +...+j1r1+j0]
            * | where k'=0,1,...,(ri)-1}
            * 
            * k in this case is specifrically from 0 to [(N/(r1r2...ri))-1]
            * 0<=j0<r1
            * 0<=j1<r2
            * 0<=j2<r3
            * ...
            * 0<=ji-1<ri
            * where ji does not follow suit
            * due to being a remainder
            * 0<=ji<N/(r1r2...ri)
            * 
            * This being said our fourier transform will be
            * j=ji(r1r2...ri)+ji-1(r1r2...ri-1)+...
            * +j1r1+j0
            * where 
            * 0<=j0<r1
            * 0<=j1<r2
            * ...
            * 0<=ji-1<ri
            * with the exception of ji
            * where 
            * 0<=ji<N/(r1r2...ri)
            * 
            * x(j)=Sum{Ai(k,ji-1,ji-2,...,j1,j0)W^kj
            * |where k=0,1,...,N/(r1r2..ri)}
            *
            * %%%%%%%%%%%%%%%%%%%% Flattening %%%%%%%%%%%%%%%%%%%%
            * 
            * since j = ji(r1r2...ri)+ji-1(r1r2...ri-1)+... +j1r1+j0
            * There exsts a one to one and onto transformation between
            * j and the {ji}
            * 
            * we know that {ji} is onto since the sum
            * will cover every integer in [0,N)
            * (Euclidean division)
            * 
            * proof:
            * j0=j mod r1
            * j1r1=(j-j0) mod r2r1 
            * j2r1r2=(j-j0-j1r1) mod r3r2r1
            * ...
            * ji-1r1r2..ri-1=(j-j0-j1r1-...ji-2r1r2...ri-1) mod r1r2...ri
            * 
            * where ji is the only entry that does not follow pattern
            * jir1r2...ri=j-j0-j1r1-...ji-1r1r2...ri-1
            * 
            * thus we can find the unique {ji} for a given j
            * therefore we have proven that j and {ji} are one to one
            * and onto. 
            * 
            * In that case instead of Ai being a function
            * of Ai(k,ji-1,ji-2...,j1,j0)
            * we can define Ai as a function of just k and j
            * Ai(k,j)
            * 
            * In retrospect if we redefine {ji} recursively
            * like we have done above then we can define
            * Ai+1(k,j)=Sum{Ai(k+k'N/(r1r2...ri),j')W^kj|
            * k'=0,1,...,(ri)-1,j' = j mod r1r2...ri} 
            * 
            * %%%%%%%%%%%%%%%% Determining Size of Ai %%%%%%%%%%%%%%%%
            * the wavenumber has a count of N/r1r2...ri
            * the count of j will be r1r2...ri
            * 
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * System.Numerics.Complex[,] preFFT
            *  = This is the previously calculated Ai term 
            * to calculate Ai+1. This will always be in the
            * form Ai+1[k,j] where k has N/r1...ri terms
            * and j has r1...ri terms.
            * 
            * List<PrimeFactorIndex> primes
            *  = This is list of primes produced by the 
            *  method CompositeListGenerator for the number
            *  of data points N. 
            * 
            * System.Numerics.Complex W
            *  = This is the complex number 
            *  e^i2pi/N that you are using to 
            *  do the Fourier Transform.
            * 
            * int stepi
            *  = the step of the process that we
            *  are doing. 
            *  we will be using entry primes[stepi]
            *  to get the SubComposite and LeftOverComposite
            *  for the ranges. Since we use the "FirstStepFFT"
            *  to calculate A1 from our list of points. stepi will
            *  start at stepi=1 and will calculate Ai+1. 
            *  stepi will end at stepi=primes.Count - 2
            *  making kprime iterate over the last prime
            *  number making up N. Therefore will use 
            *  "LastStepFFT" at primes.Count - 1 effectively
            *  hence the for loop should be: 
            *  for(int stepi=1;stepi < primes.Count-1;stepi++)
            *  
            *  
            * 
            * ********Please Edit for consistency********
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements


            // the first step to stepi is usually zero
            // in computer convention so we might need
            // to add the line stepi = stepi-1
            // stepi = stepi -1;

            //This initiates the return of the function to be an array
            // with a known size
            System.Numerics.Complex[,] postFFT =
                new System.Numerics.Complex[primes[stepi].LeftOverComposite,
                primes[stepi].SubComposite];

            // This iterates through the jth index of our postFFT array
            for (Int32 j = 0; j < primes[stepi].SubComposite; j++)
            {
                // in are preFFT, j constitutes a smaller memory space
                Int32 jprime = j % primes[stepi].SubComposite;
                // this term is precomputed for ease of computation to calculate
                // the sum
                Int32 Omega = j * primes[stepi].LeftOverComposite;

                //This iterates through the k index that will be in postFFT

                for (Int32 kprime = 0;
                    kprime < primes[stepi].LeftOverComposite;
                    kprime++)
                {
                    postFFT[kprime, j] = System.Numerics.Complex.Zero;

                    // since the k index is a larger space in preFFT and
                    // will have to be Contracted, ktilda is the actual
                    //contracted part of our sum.
                    for (Int32 ktilda = 0;
                        ktilda < primes[stepi].Prime.Numb;
                        ktilda++)
                    {
                        

                        postFFT[kprime, j] = System.Numerics.Complex.Add(
                            postFFT[kprime, j],
                            System.Numerics.Complex.Multiply(
                            System.Numerics.Complex.Pow(W, ktilda * Omega)
                            , preFFT[ktilda * primes[stepi].LeftOverComposite
                            + kprime, jprime]
                            )
                        );
                    }


                }
            }

            return postFFT;

        }

        public System.Numerics.Complex[,] FirstStepFFT(
            System.Numerics.Complex[] preFFT,
            List<PrimeFactorIndex> primes,
            System.Numerics.Complex W)
        {
            /*
            * SYNOPSIS: 
            * This is to the first step of the 
            * DFT calculation. Could be an override
            * for OneStepFFT. 
            * 
            * DESCRIPTION:
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            System.Numerics.Complex[,] postFFT =
                new System.Numerics.Complex[primes[0].LeftOverComposite,
                primes[0].SubComposite];
            for (Int32 j = 0; j < primes[0].SubComposite; j++)
            {
                // this term is precomputed for ease of computation to calculate
                // the sum
                Int32 Omega = j * primes[0].LeftOverComposite;

                //This iterates through the k index that will be in postFFT

                for (Int32 kprime = 0;
                    kprime < primes[0].LeftOverComposite;
                    kprime++)
                {
                    postFFT[kprime, j] = System.Numerics.Complex.Zero;

                    // since the k index is a larger space in preFFT and
                    // will have to be Contracted, ktilda is the actual
                    //contracted part of our sum.
                    for (Int32 ktilda = 0;
                        ktilda < primes[0].Prime.Numb;
                        ktilda++)
                    {

                        postFFT[kprime, j] =System.Numerics.Complex.Add(
                            postFFT[kprime,j],
                            System.Numerics.Complex.Multiply(
                                System.Numerics.Complex.Pow(W, ktilda * Omega)
                                , preFFT[ktilda * primes[0].LeftOverComposite
                                + kprime]
                            )
                        );
                    }


                }


            }

            return postFFT;

        }



        public System.Numerics.Complex[] LastStepFFT(
            System.Numerics.Complex[,] preFFT,
            List<PrimeFactorIndex> primes,
            System.Numerics.Complex W
            )
        {
            /*
            * SYNOPSIS: 
            * This is to the first step of the 
            * DFT calculation. Could be an override
            * for OneStepFFT. 
            * 
            * DESCRIPTION:
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // the last index might have to be changed
            int lastInd = primes.Count - 1;
            System.Numerics.Complex[] postFFT =
                new System.Numerics.Complex[primes[lastInd].TotalComposite];
            for (Int32 j = 0; j < primes[lastInd].SubComposite; j++)
            {
                // this term is precomputed for ease of computation to calculate
                // the sum
                Int32 Omega = j * primes[lastInd].LeftOverComposite;

                //This iterates through the k index that will be in postFFT

                for (Int32 kprime = 0;
                    kprime < primes[lastInd].LeftOverComposite;
                    kprime++)
                {
                    preFFT[kprime, j] = System.Numerics.Complex.Zero;

                    // since the k index is a larger space in preFFT and
                    // will have to be Contracted, ktilda is the actual
                    //contracted part of our sum.
                    for (Int32 ktilda = 0;
                        ktilda < primes[lastInd].Prime.Numb;
                        ktilda++)
                    {

                        postFFT[j] = System.Numerics.Complex.Add(
                            postFFT[j],
                            System.Numerics.Complex.Multiply(
                                System.Numerics.Complex.Pow(W, ktilda * Omega)
                                , preFFT[ktilda * primes[lastInd].LeftOverComposite
                                + kprime,j]
                            )
                        );
                        
                    }


                }


            }

            return postFFT;
        }

        public System.Numerics.Complex[] FFTCalculation(
            System.Numerics.Complex[] InputArr
            )
        {
            /*
            * SYNOPSIS: 
            * This will be the function
            * that calculate the FFT for an
            * input array (InputArr)
            * given in complex numbers
            * 
            * if the InputArr = x(k)
            * then the output array is 
            * A(j)=Sum{x(k)*exp(-i2pikj/N)| where k=0,1,...,N-1}
            * where A(j) is an effective amplitude of the 
            * wavenumber j
            * 
            * 
            * DESCRIPTION:
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements


            // Step 1: Find the prime factorization of the
            // size of the InputArr to optimize FFT

            var primeFactorList = FindPrimeFactors(
                InputArr.Count());

            // Step 2: List the prime factors in a sequence
            // to to calculate the FFT components

            var primes = CompositeListGenerator(
                primeFactorList);

            // Step 3: Start doing the FFT by doing the first step
            System.Numerics.Complex W =
                System.Numerics.Complex.FromPolarCoordinates(
                    1, -2 * Math.PI / InputArr.Count()
                 );
            var FFT = FirstStepFFT(InputArr, primes, W);

            // Step 4: Iterate through the steps to transform
            // FFT continually to the last array

            for (Int32 i = 1; i < (primes.Count - 1); i++)
            {
                FFT = OneStepFFT(FFT, primes, W, i);
            }

            // Step 5: Once the last step is completed.
            // need to do the last iteration and ouput an array.
            var outputArr = LastStepFFT(FFT, primes, W);

            return outputArr;

        }

        public System.Numerics.Complex[] InverseFFT(
            System.Numerics.Complex[] InputArr
            )
        {
            /*
            * SYNOPSIS: 
            * This will be the function
            * that calculate the FFT for an
            * input array (InputArr)
            * given in complex numbers
            * 
            * if the InputArr = A(k)
            * then the output array is 
            * x(j)=Sum{A(k)*exp(i2pikj/N)| where k=0,1,...,N-1}
            * 
            * 
            * DESCRIPTION:
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements


            // Step 1: Find the prime factorization of the
            // size of the InputArr to optimize FFT

            var primeFactorList = FindPrimeFactors(
                InputArr.Count());

            // Step 2: List the prime factors in a sequence
            // to to calculate the FFT components

            var primes = CompositeListGenerator(
                primeFactorList);

            // Step 3: Start doing the FFT by doing the first step
            System.Numerics.Complex W =
                System.Numerics.Complex.FromPolarCoordinates(
                    1, 2 * Math.PI / InputArr.Count()
                 );
            var FFT = FirstStepFFT(InputArr, primes, W);

            // Step 4: Iterate through the steps to transform
            // FFT continually to the last array

            for (Int32 i = 1; i < (primes.Count - 1); i++)
            {
                FFT = OneStepFFT(FFT, primes, W, i);
            }

            // Step 5: Once the last step is completed.
            // need to do the last iteration and ouput an array.
            var outputArr = LastStepFFT(FFT, primes, W);

            // Step 6: Now we need to normalize the entries by dividing
            // All entries by N

            for(int i=0; i < outputArr.Count();i++)
            {
                outputArr[i] = System.Numerics.Complex.Divide(
                    outputArr[i],
                    outputArr.Count());
            }

            return outputArr;

        }


        public GraphDFT[] GraphingDFT(
            System.Numerics.Complex[] InputArr
            )
        {
            /*
            * SYNOPSIS: 
            * After a Discrete Fourier transform
            * creates the array (InputArr) with 
            * indicies from [0,N), it is analytically
            * more useful to recast the indicies to
            * (-N/2,N/2]. 
            *  
            * DESCRIPTION:
            * Ultimately we want to analyze the natural
            * frequencies of the data. If our data represents
            * an equally partitioned sample of a smooth curve 
            * with finite discontinuties, then the true curve 
            * can be modeled as 
            * x(t) = Sum{C(n)exp(i2pint/T)| n= -inf,...,-1,0,1,...inf}
            * where x(t) is our true underlying function over real t
            * and has a periodicity of T. 
            * 
            * If we perform a discrete fourier over t=0,...,N-1
            * then for k=0,...,N-1
            * A(k) = Sum{C(n)L(n/T-k/N;N)|n=-inf,...,-1,0,1,...,inf}
            * where
            * L(a;N):=(1-g(a)^N)/(1-g(a)) and
            * g(a):=exp[i2pia]
            * 
            * note: g(a+M)=g(a),L(a+M;N) for all integer M.
            * and |L(a;N)|^2=(1-cos(2piaN))/(1-cos(2pia))
            * 
            * A neat mathematical trick that comes from doing 
            * discrete fourier transforms is that when N->infinity
            * then 
            * |L(a;N)|^2 -> (N^2)D(a) where 
            * D(a) = 1 if a is an integer
            * and D(a)=0 if a is not an integer.
            * 
            * This is what they call a Dirac comb.
            * 
            * From a Data analytic perspective, 
            * the natural frequencies of the data will get extracted
            * proportionally to N at points where k/N-n/T-> an Integer,
            * while non-natural frequencies will diminish due to lack of
            * proportionality to N and D(a)=0 for a not an integer.
            * 
            * shifting k with the algorithm below does not change
            * the observed Dirac combs but will format the data
            * to a standard more in line with the full fourier series.
            * 
            * For the purposes of anomaly detection
            * |A(k)|^2 in (-N/2,N/2] will have even symmetry
            * Since C(n)=C(-n)* (complex conjugate). and will mathematically
            * make it better to determine the true period of the data
            * when you sort k by their amplitude. larger |A(k)| will
            * period T resonance.
            * 
            * Algo:
            * x(j) = Sum{A(k)exp(i2pikj/N)|k=0,1...,N-1}
            * Since exp(i2pikj)=exp(i2pi(k-N)j/N)
            * for all k, we can reassign the indicies
            * if 0<=k<=N/2 then k->k
            * if N/2<k<N then k->k-N
            * this makes k now index over (-N/2,N/2]
            * 
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            GraphDFT[] retVal = new GraphDFT[InputArr.Count()];

            for(int k=0; k < InputArr.Count();k++)
            {
                retVal[k].OriginalIndex = k;
                retVal[k].Amplitude = InputArr[k];
                retVal[k].Intensity =
                    System.Numerics.Complex.Multiply(
                    InputArr[k],
                    System.Numerics.Complex.Conjugate(InputArr[k])
                    ).Real;
                if (k <= InputArr.Count() / 2)
                {
                    retVal[k].AnalyticalIndex = k;
                    retVal[k].WaveNumber = k / InputArr.Count();
                }
                else
                {
                    int shiftk = k - InputArr.Count();
                    retVal[k].AnalyticalIndex = shiftk;
                    retVal[k].WaveNumber = shiftk / InputArr.Count();
                }

            }

            return retVal;

        }

        public GraphDFT[] GraphSort(GraphDFT[] array, int leftIndex, int rightIndex)
        {
            /*
            * SYNOPSIS: 
            * This will quickly sort from the least intense
            * wave to the most intense wave.
            * 
            * DESCRIPTION:
            * This will use quick sort to quickly sort 
            * the values from the Discrete Fourier transform
            * from the smallest intensity to the largest intensity
            * 
            * We will flip this array in a later function
            * to go from the largest intensity to the smallest intensity
            * to determine the common period of the data.
            * 
            * leftIndex initiates at 0
            * while rightIndex inititiates at array.Count()-1
            * 
            * 
            *
            * Algorithm from 
            * https://code-maze.com/csharp-quicksort-algorithm/
            * 
            *
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            var i = leftIndex;
            var j = rightIndex;
            var pivot = array[leftIndex];

            while (i <= j)
            {
                while (array[i].Intensity < pivot.Intensity)
                {
                    i++;
                }

                while (array[j].Intensity > pivot.Intensity)
                {
                    j--;
                }

                if (i <= j)
                {
                    GraphDFT temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                    i++;
                    j--;
                }
            }

            if (leftIndex < j)
                GraphSort(array, leftIndex, j);

            if (i < rightIndex)
                GraphSort(array, i, rightIndex);

            return array;
        }

        public NormGraphDFT[] NormalizeGraphDFT(GraphDFT[] array)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This takes the output of GraphSort 
            * flips the order so that
            * the array starts at the largest intensity
            * to the smallest intensity. 
            * It also normalizes the data so that 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            var retArr = new NormGraphDFT[array.Count()];

            double normSquared = 0;
            for(int k=1; k <= array.Count(); k++)
            {
                int newk = array.Count() - k;
                retArr[newk].OriginalGraph = array[k];
                normSquared += array[k].Intensity;


            }

            double norm = Math.Sqrt(normSquared);
            for(int k=0; k<retArr.Count();k++)
            {
                retArr[k].NormIntensity =
                    retArr[k].OriginalGraph.Intensity / normSquared;
                retArr[k].NormAmplitude =
                    retArr[k].OriginalGraph.Amplitude / norm;
                retArr[k].Norm = norm;
            }

            return retArr;

        }

        public GCF GetGCFOfTwoNumbers(int num1, int num2)
        {
           /*
           * SYNOPSIS: 
           * 
           * DESCRIPTION:
           * This get the Greatest common
           * factor of two numbers. So we can get the 
           * greatest common factor of a list.
           * This initial process will not be the fastest
           * for the first two numbers since
           * we will be picking up the prime numbers
           * that compose the greatest common factor.
           * However this will greatly reduce the time
           * when finding the greatest common factor of many 
           * numbers
           * ver getting the greatest common factor of a large
           * set of numbers. 
           * 
           * 
           * DEPENDENCIES:
           * 
           * OUTPUTS:
           *
           * PARAMETER:
           * 
           * EXAMPLE:
           */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            // Step 0: list factors
            int[] inputNumbers = new int[] { num1,num2};
            
            // Step 1: use the smallest number to determine GCF
            int c = 1;
            if(num1 <= num2)
            {
                c = num1;
            }
            else
            {
                c = num2;
            }

            // Step 2: list all prime numbers <= smallest of two numbers 
            var allprimes = PrimeNumbersLessthan(c);

            // defining our greatest common factor as retFactor
            int retFactor = 1;
            // defining the list of primeFactors defining our GCF
            List<PrimeFactor> primeFactors = new List<PrimeFactor>();

            // Step 3: Determine prime factor decomposition for GCF
            foreach (int q in allprimes)
            {
                // Checks to see if q is a prime factor of GCF
                if (num1 % q == 0 && num2 % q == 0)
                {
                    // Checks the power of the number

                    // N is our checker M is status quo and
                    // prime is the prime factor contribution

                    int N = q*q;
                    int M = q;
                    PrimeFactor prime = new PrimeFactor { Numb = q, Pow = 1 };

                    // The loop updates our checker, status quo and power
                    while (num1 % N == 0 && num2 % N == 0)
                    {
                        M = N;
                        N = q * N;
                        prime.Pow += 1;
                    }

                    // Generates the GCF and the primeFactor List
                    retFactor = retFactor*M;
                    primeFactors.Add(prime);
                }
            }

            GCF retVal = new GCF
            {
                InputNumbers = inputNumbers,
                Primes = primeFactors,
                Factor = retFactor
            };

            return retVal;

        }


        public GCF GreatestCommonFactor(int[] numbers)
        {
            /*
           * SYNOPSIS: 
           * 
           * DESCRIPTION:
           * This is an optimal way to calculate
           * the total number of the greatest common 
           * factor. 
           * First we sort the array from smallest number
           * to biggest to determine what must be the
           * largest possibility for the greatest common factor.
           * 
           * using the GetGCFOfTwoNumbers we get the greatest 
           * common factor of those two numbers and their prime 
           * factorization. 
           * since that the greatest common factor must have
           * a sub prime factorization then we use the prime
           * factorization of these two numbers to generate our
           * basis of comparison. this will greatly decrease
           * the list of numbers we try for divisibiliity.
           * 
           * DEPENDENCIES:
           * 
           * OUTPUTS:
           *
           * PARAMETER:
           * 
           * EXAMPLE:
           */
            numbers = SortArrayInt(numbers,0,numbers.Count()-1);
            

            GCF retValue = new GCF();
            retValue.InputNumbers = numbers;
            retValue.Factor = 1;

            GCF firstGCF = GetGCFOfTwoNumbers(numbers[0], numbers[1]);
            List<PrimeFactor> retPrimes =  firstGCF.Primes;
            for (int i=2;i<numbers.Count();i++)
            {

                List<PrimeFactor> tempPrimes = new List<PrimeFactor>();
                foreach (PrimeFactor prime in retPrimes)
                {
                    int N = prime.Numb;
                    int M = 1;
                    int power = 0;
                    for(int l=1;l<=prime.Pow;l++)
                    {
                        if (numbers[i]%N==0)
                        {
                            M = N;
                            N = N * prime.Numb;
                            power += 1;
                        }
                    }
                    if(power>0)
                    {
                        PrimeFactor newPrime = new PrimeFactor {
                            Numb = M,
                            Pow = power
                        };
                        tempPrimes.Add(newPrime);
                    }

                }
                retPrimes = tempPrimes;
            }

            if(retPrimes.Count>0)
            { 
                foreach(PrimeFactor prime in retPrimes)
                {
                    retValue.Primes.Add(prime);
                    for(int i = 0; i< prime.Pow;i++)
                    {
                        retValue.Factor =
                            retValue.Factor * prime.Numb;
                    }
                }
            }

            return retValue;

        }


        public int GCFEuclidAlgo(int[]numbers)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * if g is a factor of
            * a and b
            * then g is also a factor of
            * a-b
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            
            while (numbers.Count() > 1)
            { 
                // num1 <= num2
                int num1 = numbers[0];
                int num2 = numbers[1];
                // flips order if wrong
                
                if (numbers[0] > numbers[1])
                {
                    num2 = numbers[0];
                    num1 = numbers[1];
                }

                int rem = num2 % num1;
                int tempNum1 = num1;
                int tempNum2 = num2;
                if(num1 == num2)
                {
                    tempNum2 = num1;
                }
                else
                { 
                    while(rem!=0)
                    {
                        tempNum2 = rem;
                        rem = tempNum1 % rem;
                        tempNum1 = tempNum2;

                    }
                }

                int[] arr = new int[numbers.Count() - 1];
                arr[0] = tempNum2;
                for(int i=2;i<numbers.Count();i++)
                {
                    arr[i - 1] = numbers[i];
                }

                /*// This recursive step needs to be checked
                if (arr.Count() > 1)
                    GCFEuclidAlgo(arr);
                return arr;

                */

                numbers = arr;

            }

            return numbers[0];


        }


        public PeriodGraph[] CalculateCommonTime (NormGraphDFT[] graph)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This calculates all of the common
            * periods in order. where the 0th index
            * is the calculated period from the 0th index
            * while the ith index has the common period
            * from the 0th,1st,...,ith indicie.
            * 
            * The idea is that the graph was originally sorted
            * by the intensity of the where the most
            * intense entry is th 0th index and the 
            * least intense entry is at the graph.Count()-1
            * entry. This way the accumulative total
            * will be maximized for the common period calculated.
            * 
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            PeriodGraph[] _calcComTime = new PeriodGraph[graph.Count()];

            _calcComTime[0].AccumulatedSignal = graph[0].NormIntensity;
            _calcComTime[0].CommonPeriod =
                Math.ReciprocalEstimate(graph[0].OriginalGraph.WaveNumber);
            _calcComTime[0].FullNormGraph = graph;
            _calcComTime[0].ConsideredNormGraph = new NormGraphDFT[1]
            { graph[0] };
            _calcComTime[0].LogPeriod =
                Math.Log(_calcComTime[0].CommonPeriod);

            int gcfk = Math.Abs(graph[0].OriginalGraph.AnalyticalIndex);

            for(int i=1; i<graph.Count();i++)
            {
                int[] gcfTempArr = new int[]
                {
                    gcfk,
                    Math.Abs(graph[i].OriginalGraph.AnalyticalIndex)
                };

                gcfk = GCFEuclidAlgo(gcfTempArr);

                _calcComTime[i].CommonPeriod = graph.Count() / gcfk;

                _calcComTime[i].AccumulatedSignal = graph[i].NormIntensity
                    + _calcComTime[i - 1].AccumulatedSignal;
                _calcComTime[i].FullNormGraph = graph;
                _calcComTime[i].ConsideredNormGraph = graph[0..i];
                _calcComTime[i].LogPeriod =
                Math.Log(_calcComTime[i].CommonPeriod);


            }

            return _calcComTime;

        }

        public double[] SimpleDerivative (double[] array)
        {
            /*
            * SYNOPSIS: 
            * Estimate of the derivative 
            * in an array
            * 
            * DESCRIPTION:
            * This calculates the estimate of 
            * a derivative in an array with a 
            * step size of 1. 
            * 
            * If the underlying function that 
            * the array is sampling from is analytic
            * then the error is on the order of 
            * h^2 for all entries except the end points
            * where the error is on the order h 
            * where h is the step size.
            * 
            * using; 
            * f'(xi)= (f(xi+1)-f(xi-1))/(xi+1-xi-1)
            * for all points but end points which replace
            * the missing entry with the ith index
            * 
            * We are interested in relative comparisons
            * between the period and the cumulative total
            * thus the step size does not matter for this 
            * calculation and will be set to 1. 
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts
            double[] derivative = new double[array.Count()];

            derivative[0] = array[1] - array[0];
            derivative[array.Count() - 1] =
                array[array.Count() - 1] - array[array.Count() - 2];
            for(int i=1; i < array.Count()-1; i++)
            {
                derivative[i] = array[i + 1] - array[i - 1] / 2;
            }

            return derivative;
        }

        public double[] SimpleSecondDerivative(
            double[] originalData,
            double[] firstDerivative
            )
        {
            /*
            * SYNOPSIS: 
            * Estimate of the second derivative 
            * in an array
            * 
            * DESCRIPTION:
            * This calculates the estimate of 
            * the second derivative in an array with a 
            * step size of h=1. 
            * 
            * the firstDerivative is only to get 
            * an estimate of the second derivative at the
            * end points which may be inaccurate
            * this is calculated with 
            * f''(xi) = (f'(xi+1)-f'(xi))/h
            * 
            * 
            * the originalData will be used to generate most
            * of the points in the data. the step 
            * dependency is on order of h^2 for these calculations
            * these are calculated using the relationship
            * f''(xi) = (f(xi+1)+f(xi-1)-2f(xi))/(2h^2)
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            // Initiating second derivative array
            double[] secondDer = new double[originalData.Length];

            // Calculating endpoints
            secondDer[0] = firstDerivative[1] - firstDerivative[0];
            secondDer[originalData.Length - 1] =
                firstDerivative[originalData.Length - 1]
                - firstDerivative[originalData.Length - 2];

            for(int i=1; i < originalData.Length-1;i++)
            {
                secondDer[i] = originalData[i + 1] + originalData[i - 1]
                    - 2 * originalData[i];

            }

            return secondDer;



        }

        public Polynomial PolyAddition(
            Polynomial poly1,
            Polynomial poly2
            )
        {
            /*
            * SYNOPSIS: 
            * This calculates a polynomial addition
            * 
            * DESCRIPTION:
            * polynomials are arrays where
            * the indicies represent the power
            * and the array[index] gives the coefficient
            * 
            * this sums over shared indicies and shrinks 
            * the size of the polynomial if the highest 
            * degree cancels.
            * what indices can't be summed over will iterate
            * over remaining entries.
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts
            int deg = 1;
            int lastIndex = 0;
            bool p1biggerp2 = false;
            if (poly1.Degree > poly2.Degree)
            {
                deg = poly1.Degree;
                lastIndex = poly2.Degree;
                p1biggerp2 = true;
            }
            else
            {
                deg = poly2.Degree;
                lastIndex = poly1.Degree;
            }

            double[] sum = new double[deg + 1];

            for(int i=0;i<=lastIndex; i++)
            {
                sum[i] = poly1.Array[i] + poly2.Array[i];
            }

            if(poly1.Degree==poly2.Degree)
            {
                int trueDegree = lastIndex;
                for(int i=lastIndex; i>=0;i--)
                {
                    if (!(sum[i]==0))
                    {
                        trueDegree = i;
                        break;
                    }
                }
                Polynomial final = new Polynomial
                {
                    Array = sum[0..trueDegree],
                    Degree = trueDegree
                };

                return final;

            }
            else if(p1biggerp2)
            {
                for(int i=lastIndex+1;i<sum.Length;i++)
                {
                    sum[i] = poly1.Array[i];
                }
                Polynomial final = new Polynomial
                {
                    Array = sum,
                    Degree = sum.Length-1
                };
                return final;
            }
            else
            {
                for (int i = lastIndex + 1; i < sum.Length; i++)
                {
                    sum[i] = poly2.Array[i];
                }
                Polynomial final = new Polynomial
                {
                    Array = sum,
                    Degree = sum.Length - 1
                };
                return final;
            }


        }

        public Polynomial PolyMultiplication(
            Polynomial poly1,
            Polynomial poly2
            )
        {
            /*
            * SYNOPSIS: 
            * This multiplies polynomials
            * 
            * DESCRIPTION:
            * This multiplies the polynomials
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            // Step 1: figure out the degree 
            // thus the size of the array to construct the array

            
            int degree = poly1.Degree * poly2.Degree;
            bool multByConstant = false;
            bool constp1 = false;
            if(poly1.Degree ==0)
            {
                degree = poly2.Degree;
                multByConstant = true;
                constp1 = true; 
            }
            else if (poly2.Degree ==0)
            {
                degree = poly1.Degree;
                multByConstant = true;
                
            }
            Polynomial polynomial = new Polynomial
            {
                Degree = degree,
                Array = new double[degree + 1]
            };


            if(!multByConstant)
            { 
                for(int i=0; i <poly1.Array.Length;i++)
                {
                    for(int j = 0; j<poly2.Array.Length;i++)
                    {
                        polynomial.Array[i + j] +=
                            poly1.Array[i] * poly2.Array[j];
                    }
                }
            }
            else if(constp1)
            {
                for (int i = 0; i < poly2.Array.Length ;i++)
                {
                    polynomial.Array[i] = poly2.Array[i] * poly1.Array[0];
                }
            }
            else
            {
                for (int i = 0; i < poly1.Array.Length; i++)
                {
                    polynomial.Array[i] = poly1.Array[i] * poly2.Array[0];
                }
            }
            return polynomial;



        }

        public PolyDivisionOutput PolyDivision(
            Polynomial numerator,
            Polynomial divisor
            )
        {
            /*
            * SYNOPSIS: 
            * This method divdes one polynomial 
            * by another
            * 
            * DESCRIPTION:
            * This will output the results of polynomial division
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            PolyDivisionOutput _polyDivision = new PolyDivisionOutput();
            _polyDivision.Numerator = numerator;
            _polyDivision.Divisor = divisor;

            Polynomial remainder = numerator;
            Polynomial quotient = new Polynomial
            {
                Degree = 0,
                Array = new double[1]
            };
            while (remainder.Degree >= divisor.Degree)
            {
                int degMult = remainder.Degree = divisor.Degree;

                double coefMult =
                    remainder.Array[remainder.Degree]
                    / divisor.Array[divisor.Degree];

                Polynomial subQuotient = new Polynomial
                {
                    Degree = degMult,
                    Array = new double[degMult + 1]
                };
                quotient = PolyAddition(quotient, subQuotient);
                
                
                var subtractor = new Polynomial
                {
                    Degree = remainder.Degree,
                    Array = new double[remainder.Degree] // this assigns zeros
                };


                for(int i=0; i<divisor.Array.Length;i++)
                {
                    subtractor.Array[i + degMult] =
                        -1 *coefMult* divisor.Array[i];
                }

                remainder = PolyAddition(remainder, subtractor);

            }

            _polyDivision.Quotient = quotient;
            _polyDivision.Remainder = remainder;

            return _polyDivision;




        }

        public Polynomial ConvertArrayToPoly(double[] polynomial)
        {
            /*
            * SYNOPSIS: 
            * Converts an array to a polynomial
            * 
            * DESCRIPTION:
            * This is a quick way to convert an 
            * array to a polynomial with the assumption
            * that the array is ordereed where
            * polynomial[0] -> x^0 (const)
            * polymomial[1] -> x^1
            * polynomial[2] -> x^2 
            * ...
            * polynomial[N-1]->x^N-1
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts
            Polynomial _convert = new Polynomial
            {
                Array = polynomial,
                Degree = polynomial.Length-1
            };
            return _convert;
        }

        public double EvaluatePoly(Polynomial poly, double x)
        {
            /*
            * SYNOPSIS: 
            * This evaluates polynomials
            * 
            * DESCRIPTION:
            * This evaluates the polynomials
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            double y = poly.Array[0];

            for(int i=1;i<poly.Array.Length;i++)
            {
                y += poly.Array[i] * Math.Pow(x, i);
            }
            return y;
        }

        public Polynomial PolynomialInterpolation(double[] xArr, double[] yArr)
        {
            /*
            * SYNOPSIS: 
            * This interpolates a polynomial with the yaxis
            * representing yArr and the xaxis as xArr
            * 
            * DESCRIPTION:
            * This interpolates a polynomial with the yaxis
            * representing yArr and the xaxis as xArr
            * 
            * This uses a method where for N points 
            * there exists a unique N-1 polynomial that
            * contains all of the points. 
            * The equation is 
            * sum{yi*prod[(x-xj)/(xi-xj)|where j=0,..,N-1 and j!=i]
            * |where i=0,..,N-1} 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            //***** An error that could be added
            //is a case where there are two xArr entries that
            //are the same *****


            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts



            Polynomial finalPoly = new Polynomial
            {
                Degree = xArr.Length - 1,
                Array = new double[xArr.Length]
            };

            Polynomial[] constructor = new Polynomial[xArr.Length];
            for(int i = 0; i< xArr.Length;i++)
            {
                constructor[i].Degree = 1;
                constructor[i].Array = new double[2] { -xArr[i], 1 };
            };


            for(int i=0; i< xArr.Length;i++)
            {
                Polynomial subNomial = new Polynomial
                {
                    Degree = 0,
                    Array = new double[1] { 1 }
                };
                for(int j=0;j<i;j++)
                {
                    subNomial = PolyMultiplication(
                        subNomial,
                        constructor[j]
                        );
                }
                for(int j=i+1;j<xArr.Length;j++)
                {
                    subNomial = PolyMultiplication(
                        subNomial,
                        constructor[j]
                        );
                }

                double divisor = EvaluatePoly(subNomial, xArr[i]);

                for(int j=0;j<xArr.Length;j++)
                {
                    subNomial.Array[j] =
                        yArr[j]*subNomial.Array[j] / divisor;
                }

                finalPoly = PolyAddition(finalPoly, subNomial);

            }

            return finalPoly;


        }



        public Polynomial PolyDerivative(Polynomial poly)
        {
            /*
            * SYNOPSIS: 
            * This takes the derivative of a polynomial
            * 
            * DESCRIPTION:
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            Polynomial derivative = new Polynomial
            {
                Degree = poly.Degree - 1,
                Array = new double[poly.Array.Length - 1]
            };

            for(int i=1;i<poly.Array.Length;i++)
            {
                derivative.Array[i - 1] = poly.Array[i] * i;
            }

            return derivative;
        }


        public double[] PolyRoots(
            Polynomial poly,
            double convergence,
            double rangeMin,
            double rangeMax,
            int NumberOfAttempts
            )
        {
            /*
            * SYNOPSIS: 
            * This uses newton's method to
            * find roots.
            * 
            * DESCRIPTION:
            * This algorithm splits the range
            * into poly.Degree equal portions. each point
            * will start as a starting point for evaluating the 
            * root. 
            * 
            * DEPENDENCIES:
            * PolyDerivative
            * EvaluatePoly
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            // need to evaluate derivate at each point
            Polynomial derivative = PolyDerivative(poly);


            // The maximum number of roots possible is the degree of poly
            // hence split the region by the max number of roots.
            double range = rangeMax - rangeMin;
            double step = range / poly.Degree;
            double pos = 0.5 * step;
            List<double> roots = new List<double>();
            for(int i=0;i<poly.Degree;i++)
            {
                double root = pos;
                double polyRoot = EvaluatePoly(poly, root);
                double derRoot = EvaluatePoly(derivative, root);
                int Attempts = 0;
                while(
                    !(Math.Abs(polyRoot)<convergence)
                    && Attempts<NumberOfAttempts
                    )
                {
                    if(derRoot==0)
                    {
                        root = pos + convergence;
                        polyRoot = EvaluatePoly(poly, root);
                        derRoot = EvaluatePoly(derivative, root);
                    }
                    else
                    {
                        root = root - (polyRoot / derRoot);
                        polyRoot = EvaluatePoly(poly, root);
                        derRoot = EvaluatePoly(derivative, root);
                    }

                    Attempts += 1;
                }

                /*
                 * When we add a logging portion
                 * we should report logs for when
                 * the Attempts=NumberOfAttempts
                 * 
                 *
                 */

                bool sameRoot = false;
                foreach(double prevRoot in roots)
                {
                    // could change this to root==prevRoot
                    if(Math.Abs(root - prevRoot) < convergence)
                    {
                        sameRoot = true;
                    }
                }
                if(!sameRoot)
                {
                    roots.Add(root);
                }

            }

            double[] rootArr = new double[roots.Count()]; 

            for(int i=0;i<rootArr.Length;i++)
            {
                rootArr[i] = roots[i];
            }

            return rootArr;

        }

        public PolyCritical PolyCriticalPoints(
            Polynomial poly,
            double convergence,
            double rangeMin,
            double rangeMax,
            int NumberOfAttempts,
            bool getInflections
            )
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * 
            * DEPENDENCIES:
            * PolyDerivative
            * EvaluatePoly
            * PolyRoots
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts


            // Step 1: find derivatives
            Polynomial firstDer = PolyDerivative(poly);
            Polynomial secondDer = PolyDerivative(firstDer);

            double[] firstDerRoots = PolyRoots(
                firstDer,
                convergence,
                rangeMin,
                rangeMax,
                NumberOfAttempts
                );

            List<double> xminList = new List<double>();
            List<double> yminList = new List<double>();
            List<double> xmaxList = new List<double>();
            List<double> ymaxList = new List<double>();
            List<double> xsaddles = new List<double>();
            List<double> ysaddles = new List<double>();
            for(int i=0; i<firstDerRoots.Length;i++)
            {
                double ySecondDer = EvaluatePoly(secondDer, firstDerRoots[i]);
                double y = EvaluatePoly(poly, firstDerRoots[i]);
                if (ySecondDer>0)
                {
                    xminList.Add(firstDerRoots[i]);
                    yminList.Add(y);

                }
                else if (ySecondDer<0)
                {
                    xmaxList.Add(firstDerRoots[i]);
                    ymaxList.Add(y);
                }
                else
                {
                    double yless = EvaluatePoly(poly,
                        firstDerRoots[i] - convergence);
                    double ymore = EvaluatePoly(poly,
                        firstDerRoots[i] + convergence);
                    if(yless > y && ymore > y)
                    {
                        xminList.Add(firstDerRoots[i]);
                        yminList.Add(y);
                    }
                    else if (yless<y && ymore < y)
                    {
                        xmaxList.Add(firstDerRoots[i]);
                        yminList.Add(y);
                    }
                    else
                    {
                        xsaddles.Add(firstDerRoots[i]);
                        ysaddles.Add(y);
                    }
                }
            }
            int minCount = xminList.Count();
            int maxCount = xmaxList.Count();
            int saddleCount = xsaddles.Count();
            double[] xminArr = new double[minCount];
            double[] yminArr = new double[minCount];
            double[] xmaxArr = new double[maxCount];
            double[] ymaxArr = new double[maxCount];
            double[] xsadArr = new double[saddleCount];
            double[] ysadArr = new double[saddleCount];

            PolyCritical criticalPoints = new PolyCritical
            {
                Polynomial = poly,
                Convergence = convergence,
                RangeMin = rangeMin,
                RangeMax = rangeMax,

            };
            if (minCount>0)
            { 
                for(int i=0; i < minCount;i++)
                {
                    xminArr[i] = xminList[i];
                    yminArr[i] = yminList[i];
                }
                criticalPoints.MinXArr = xminArr;
                criticalPoints.MinYArr = yminArr;
            }

            if(maxCount>0)
            { 
                for (int i=0; i<maxCount;i++)
                {
                    xmaxArr[i] = xmaxList[i];
                    ymaxArr[i] = ymaxList[i];
                
                }
                criticalPoints.MaxXArr = xminArr;
                criticalPoints.MaxYArr = yminArr;
            }
            if(saddleCount>0)
            { 
                for (int i = 0; i < saddleCount; i++)
                {
                    xsadArr[i] = xsaddles[i];
                    ysadArr[i] = ysaddles[i];
                }
                criticalPoints.xSaddlePoints = xminArr;
                criticalPoints.ySaddlePoints = yminArr;
            }

            if(getInflections)
            {
                double[] inflectRoot = PolyRoots(
                    firstDer,
                    convergence,
                    rangeMin,
                    rangeMax,
                    NumberOfAttempts
                    );
                List<double> xinflections = new List<double>();
                List<double> yinflections = new List<double>();
                foreach(double root in inflectRoot)
                { 
                    double y=EvaluatePoly(secondDer, root);
                    double yless=EvaluatePoly(secondDer, root-convergence);
                    double ymore = EvaluatePoly(secondDer, root + convergence);

                    if(ymore>=0 && yless<=0)
                    {
                        xinflections.Add(root);
                        yinflections.Add(y);
                    }
                    else if(ymore<=0 && yless>=0)
                    {
                        xinflections.Add(root);
                        yinflections.Add(y);

                    }

                }
                int inflectCount = xinflections.Count();
                if (inflectCount>0)
                {
                    criticalPoints.xInflectionPoints = new double[inflectCount];
                    criticalPoints.yInflectionPoints = new double[inflectCount];

                    for (int i=0; i < xinflections.Count();i++)
                    {
                        criticalPoints.xInflectionPoints[i] =
                            xinflections[i];
                        criticalPoints.yInflectionPoints[i] = yinflections[i];
                    }


                }

            }



            return criticalPoints;


        }


        public DiscreteRoot[] FindDiscreetRoots(
            double[] data,
            double param,
            double rootConverge)
        {
            /*
            * SYNOPSIS: 
            * For discrete data of an array
            * this will find the an array of approximate
            * values of roots. 
            * 
            * DESCRIPTION:
            * Finds approximate roots of discrete data
            *
            * The idea behind this is that the data 
            * can at some point intersect with the x axis
            * what this will do is that it will list out the 
            * indicies that are the closest to the x axis within 
            * some parameter. 
            * 
            * (((This parameter should be determined from the 
            * varriance of the data to in some sense.)))
            * 
            * 
            * 
            * 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            List<List<int>> PossibleRoots = new List<List<int>>();
            var _discreetRoots = new List<DiscreteRoot>();
            bool tracking = false;
            int lastIndex = 0;
            // this condition lists the possible roots that
            // tthe function has
            // if multiple points in sequence are near zero
            // then we will fit a polynomial among all near zero points
            // to get a better estimate of root.
            for (int i = 0; i < data.Length; i++)
            {
                if (Math.Abs(data[i]) < param)
                {
                    if (tracking)
                    {
                        PossibleRoots[lastIndex].Add(i);
                    }
                    else
                    {
                        var subRoots = new List<int>();
                        subRoots.Add(i);
                        PossibleRoots.Add(subRoots);
                        tracking = true;
                        lastIndex = PossibleRoots.Count() - 1;
                    }
                }
                else
                {
                    tracking = false;
                }
            }


            // Next step is to fit polynomials

            foreach (List<int> rootSet in PossibleRoots)
            {
                // Initializing DiscreteRoot object
                DiscreteRoot discreteRoot = new DiscreteRoot();
                // sample near root array creation
                bool initialPoint = false;
                bool finalPoint = false;
                int arrSize = 0;

                if (rootSet.Contains(0) && rootSet.Contains(data.Length - 1))
                {
                    initialPoint = true;
                    finalPoint = true;
                    arrSize = rootSet.Count();
                }
                else if (rootSet.Contains(data.Length - 1))
                {
                    finalPoint = true;
                    arrSize = rootSet.Count() + 1;
                }
                else if (rootSet.Contains(0))
                {
                    initialPoint = true;
                    arrSize = rootSet.Count() + 1;
                }
                else
                {
                    arrSize = rootSet.Count() + 2;
                }

                double[] polySampX = new double[arrSize];
                double[] polySampY = new double[arrSize];
                discreteRoot.SampleIndices = new int[arrSize];
                int count = 0;
                if (initialPoint)
                {
                    foreach (int root in rootSet)
                    {
                        polySampX[count] = root;
                        polySampY[count] = data[root];
                        discreteRoot.SampleIndices[count] = root;
                        count++;
                    }
                }
                else
                {
                    polySampX[0] = rootSet[0] - 1;
                    polySampY[0] = data[rootSet[0] - 1];
                    discreteRoot.SampleIndices[count] = rootSet[0] - 1;
                    count++;
                    foreach (int root in rootSet)
                    {
                        polySampX[count] = root;
                        polySampY[count] = data[root];
                        discreteRoot.SampleIndices[count] = root;
                        count++;
                    }
                }

                if (!finalPoint)
                {
                    int xpos = rootSet[rootSet.Count() - 1] + 1;
                    polySampX[count] = xpos;
                    polySampY[count] = data[xpos];
                    discreteRoot.SampleIndices[count] = xpos;

                }
                discreteRoot.DataValues = polySampY;

                // Interpolate Polynomials

                Polynomial polynomial =
                    PolynomialInterpolation(polySampX, polySampY);

                discreteRoot.PolynomialFit = polynomial;
                double[] poliRoots = PolyRoots(
                    polynomial,
                    rootConverge,
                    polySampX[0],
                    polySampX[polySampX.Length - 1],
                    100
                    );
                
                var approxRootList = new List<double>();
                foreach(double rut in poliRoots)
                {
                    if (
                        polySampX[0] <= rut
                        &&
                        rut <= polySampX[polySampX.Length - 1]
                        ) 
                    {
                        approxRootList.Add(rut);
                    }
                }

                discreteRoot.RootIndex = new double[approxRootList.Count()];
                discreteRoot.PolyRootValue = new double[approxRootList.Count()];
                int count1 = 0;
                foreach(double root in approxRootList)
                {
                    discreteRoot.RootIndex[count1] = root;
                    discreteRoot.PolyRootValue[count1] = EvaluatePoly(
                        polynomial,
                        root
                        );
                    count1++;
                }

                _discreetRoots.Add(discreteRoot);

            }

            var discreetReturn = new DiscreteRoot[_discreetRoots.Count()];
            int count2 = 0;
            foreach(DiscreteRoot root in _discreetRoots)
            {
                discreetReturn[count2] = root;
                count2++;
            }
            return discreetReturn;


        }



    }
}
