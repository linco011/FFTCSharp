using System;
using DiscreteFourierTransformLibrary.Models;

namespace DiscreteFourierTransformLibrary.Modules
{
    public class FunctionalDraftSpace2
    {
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
                    // THIS IS PROBABLY NOT NEEDED
                    postFFT[kprime, j] = System.Numerics.Complex.Zero;

                    // since the k index is a larger space in preFFT and
                    // will have to be Contracted, ktilda is the actual
                    //contracted part of our sum.
                    for (Int32 ktilda = 0;
                        ktilda < primes[stepi].Prime.Numb;
                        ktilda++)
                    {

                        postFFT[kprime, j] += System.Numerics.Complex.Pow(W, ktilda * Omega)
                            * preFFT[ktilda * primes[stepi].LeftOverComposite + kprime, jprime];
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

                        postFFT[kprime, j] +=
                            System.Numerics.Complex.Pow(W, ktilda * Omega)
                            * preFFT[ktilda * primes[0].LeftOverComposite
                            + kprime];
                    }


                }


            }

            return postFFT;

        }



        public System.Numerics.Complex[] LastStepFFT(
            System.Numerics.Complex[,] preFFT,
            List<PrimeFactorIndex>primes,
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

            int lastInd = primes.Count - 1;
            System.Numerics.Complex[,] postFFT =
                new System.Numerics.Complex[primes[lastInd].LeftOverComposite,
                primes[lastInd].SubComposite];
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
                    postFFT[kprime, j] = System.Numerics.Complex.Zero;

                    // since the k index is a larger space in preFFT and
                    // will have to be Contracted, ktilda is the actual
                    //contracted part of our sum.
                    for (Int32 ktilda = 0;
                        ktilda < primes[lastInd].Prime.Numb;
                        ktilda++)
                    {

                        postFFT[j] += System.Numerics.Complex.Pow(W, ktilda * Omega)
                            * preFFT[ktilda * primes[lastInd].LeftOverComposite + kprime];
                    }


                }


            }

            return postFFT;
        }

        public System.Numerics.Complex[] FFTCalculation(
            System.Numerics.Complex[]InputArr
            )
        {
            /*
            * SYNOPSIS: 
            * This will be the function
            * that calculate the FFT for an
            * input array (InputArr)
            * given in complex numbers
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

            var primeFactorList = MathFunctions.FindPrimeFactors(
                InputArr.Count());

            // Step 2: List the prime factors in a sequence
            // to to calculate the FFT components

            var primes = MathFunctions.CompositeListGenerator(
                primeFactorList);

            // Step 3: Start doing the FFT by doing the first step
            System.Numerics.Complex W =
                System.Numerics.Complex.FromPolarCoordinates(
                    1, 2 * Math.PI / InputArr.Count()
                 );
            var FFT = FirstStepFFT(InputArr, primes, W);

            // Step 4: Iterate through the steps to transform
            // FFT continually to the last array

            for(Int32 i=1; i< (primes.Count - 1); i++)
            {
                FFT = OneStepFFT(FFT, primes, W, i);
            }

            // Step 5: Once the last step is completed.
            // need to do the last iteration and ouput an array.
            var outputArr = LastStepFFT(FFT, primes, W);

            return outputArr;

        }

    }

}