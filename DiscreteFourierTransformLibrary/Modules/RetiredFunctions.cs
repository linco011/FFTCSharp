using System;
using DiscreteFourierTransformLibrary.Models;

namespace DiscreteFourierTransformLibrary.Modules
{
	public class RetiredFunctions
	{
        public Dictionary<List<FFTModuleIndex>, System.Numerics.Complex> OneStepFFT(
            Dictionary<List<FFTModuleIndex>, System.Numerics.Complex> preFFT,
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
            * 
            * x(j)=Sum{Ai(k,ji-1,ji-2,...,j1,j0)W^kj
            * |where k=0,1,...,N/(r1r2..ri)}
            *
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * Dictionary<List<FFTModuleIndex>, System.Numerics.Complex>
            * preFFT = This is the previously calculated Ai term 
            * to calculate Ai+1. By convention 
            * List<FFTModuleIndex> will always start with k at 
            * 0. for simplicity the order will be flipped 
            * where Ai(k,j0,j1,...,ji-1) instead of Ai(k,ji-1,...,j0) like 
            * in the description above.
            * ********Please Edit for consistency********
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements
            var keys = new List<List<FFTModuleIndex>>();
            var listSubj = new List<Int32>();
            foreach (List<FFTModuleIndex> pos in preFFT.Keys)
            {
                Int32 subj = 0;
                Int32 mult = 1;
                var key = new List<FFTModuleIndex>();
                foreach (FFTModuleIndex i in pos)
                {
                    if (i.Prime is not null)
                    {
                        subj += i.Value * mult;
                        mult = mult * i.Bound;
                        key.Add(i);
                    }
                    else
                    {
                        FFTModuleIndex newK = new FFTModuleIndex
                        {
                            Value = i.Value %
                            primes[pos.Count].LeftOverComposite,
                            Bound = primes[pos.Count].LeftOverComposite,
                            Prime = null

                        };

                        key.Add(newK);


                    }

                }

                keys.Add(key);
                listSubj.Add(subj);


            }
            // RETHINK THIS
            /*
            for (Int32 jprime = 0;
                jprime < primes[keys.Count - 1].Prime.Numb;
                jprime++)
            {
                System.Numerics.Complex valueRet = 0;
                for (Int32 kqot = 0;
                    kqot < primes[keys.Count-1].Prime.Numb;
                    kqot++
                    )
                {
                    valueRet += preFFT[keys.Count - 1];
                    }
            }
            */
            Int32 eltNum = 0;
            foreach (List<FFTModuleIndex> key in keys)
            {
                var subj = listSubj[eltNum];
                Int32 Power = key[0].Value
                    * primes[key.Count - 1].LeftOverComposite * subj;
                double doubPow = Convert.ToDouble(Power);

                System.Numerics.Complex retValue = 0;
                // Think about this
                for (double kprime = 0; kprime <= key[key.Count - 1].Bound; kprime++)
                {
                    retValue += System.Numerics.Complex.Pow(W, doubPow * kprime) *
                        preFFT[preFFT.Keys[eltNum]]
                }
                eltNum++;


            }



        }


        public Dictionary<List<FFTModuleIndex>, System.Numerics.Complex> OneStepFFT(
            Dictionary<List<FFTModuleIndex>, System.Numerics.Complex> preFFT,
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
            * 
            * x(j)=Sum{Ai(k,ji-1,ji-2,...,j1,j0)W^kj
            * |where k=0,1,...,N/(r1r2..ri)}
            *
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * Dictionary<List<FFTModuleIndex>, System.Numerics.Complex>
            * preFFT = This is the previously calculated Ai term 
            * to calculate Ai+1. By convention 
            * List<FFTModuleIndex> will always start with k at 
            * 0. for simplicity the order will be flipped 
            * where Ai(k,j0,j1,...,ji-1) instead of Ai(k,ji-1,...,j0) like 
            * in the description above.
            * ********Please Edit for consistency********
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements
            var keys = new List<List<FFTModuleIndex>>();
            var listSubj = new List<Int32>();
            foreach (List<FFTModuleIndex> pos in preFFT.Keys)
            {
                Int32 subj = 0;
                Int32 mult = 1;
                var key = new List<FFTModuleIndex>();
                foreach (FFTModuleIndex i in pos)
                {
                    if (i.Prime is not null)
                    {
                        subj += i.Value * mult;
                        mult = mult * i.Bound;
                        key.Add(i);
                    }
                    else
                    {
                        FFTModuleIndex newK = new FFTModuleIndex
                        {
                            Value = i.Value %
                            primes[pos.Count].LeftOverComposite,
                            Bound = primes[pos.Count].LeftOverComposite,
                            Prime = null

                        };

                        key.Add(newK);


                    }

                }

                keys.Add(key);
                listSubj.Add(subj);


            }
            // RETHINK THIS
            /*
            for (Int32 jprime = 0;
                jprime < primes[keys.Count - 1].Prime.Numb;
                jprime++)
            {
                System.Numerics.Complex valueRet = 0;
                for (Int32 kqot = 0;
                    kqot < primes[keys.Count-1].Prime.Numb;
                    kqot++
                    )
                {
                    valueRet += preFFT[keys.Count - 1];
                    }
            }
            */
            Int32 eltNum = 0;
            foreach (List<FFTModuleIndex> key in keys)
            {
                var subj = listSubj[eltNum];
                Int32 Power = key[0].Value
                    * primes[key.Count - 1].LeftOverComposite * subj;
                double doubPow = Convert.ToDouble(Power);

                System.Numerics.Complex retValue = 0;
                // Think about this
                for (double kprime = 0; kprime <= key[key.Count - 1].Bound; kprime++)
                {
                    retValue += System.Numerics.Complex.Pow(W, doubPow * kprime) *
                        preFFT[preFFT.Keys[eltNum]]
                }
                eltNum++;


            }



        }


        public Dictionary<int, System.Numerics.Complex> DFT(
            Dictionary<int, System.Numerics.Complex>,
            List<PrimeFactor> primes)
        {
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * Most DFT methods depend on the
            * number of points N = 2^k where 
            * k is some integer.
            *
            * even though those number of points
            * have the most efficient DFT algorithm
            * it does not apply to all custom conditions
            * in the form of data aquisition.
            *
            * This algorithm is used to perform a DFT transform
            * in the spirit of CooleyTurkey's paper for a generalized
            * N = r1*r2 which can expanded to its prime factor 
            * decomposition.
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

            // N is the total Number of Points
            double N = 1;
            foreach (PrimeFactor prime in primes)
            {
                N = N * Math.Pow(prime.Numb, prime.Pow);
            }

            // W=exp(2*pi/N) which all fourier terms are
            // W^k where k is an 0<=k<N and is an integer
            System.Numerics.Complex W =
                System.Numerics.Complex.FromPolarCoordinates(
                    1, 2 * Math.PI / N
                 );


            foreach (PrimeFactor prime in primes)
            {
                for (int pow = 1; pow <= prime.Pow; pow++)
                {

                }
            }










        }

        public DiscreteRoot[] FindDiscreetRoots(double[] data, double param)
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

            List<int> PossibleRoots = new List<int>();
            var _discreetRoots = new List<DiscreteRoot>();
            for (int i = 1; i < data.Length - 1; i++)
            {
                if (Math.Abs(data[i]) < param)
                {
                    PossibleRoots.Add(i);
                }
            }
            var rootList = new List<DiscreteRoot>();

            if (PossibleRoots.Count() > 0)
            {
                foreach (int root in PossibleRoots)
                {
                    DiscreteRoot rootElt = new DiscreteRoot();
                    rootElt.ClosestIndex = root;
                    rootElt.LeftIndex = root - 1;
                    rootElt.RightIndex = root + 1;
                    double[] xArr = new double[3]
                    {
                        rootElt.LeftIndex,
                        rootElt.ClosestIndex,
                        rootElt.RightIndex
                    };
                    rootElt.LeftValue = data[rootElt.LeftIndex];
                    rootElt.ClosestValue = data[rootElt.ClosestIndex];
                    rootElt.RightValue = data[rootElt.RightIndex];
                    double[] yArr = new double[3]
                    {
                        rootElt.LeftValue,
                        rootElt.ClosestValue,
                        rootElt.RightValue
                    };

                    rootElt.PolynomialFit =
                    PolynomialInterpolation(xArr, yArr);
                    var rootsFirstGlance =
                        PolyRoots(
                            rootElt.PolynomialFit,
                            param,
                            rootElt.LeftIndex,
                            rootElt.RightIndex,
                            100
                            );

                    List<double> ruets = new List<double>();
                    foreach (double ruet in rootsFirstGlance)
                    {
                        if (ruet < rootElt.RightIndex &&
                            ruet > rootElt.LeftIndex)
                        {
                            ruets.Add(ruet);
                        }
                    }

                    if (ruets.Count() > 0)
                    {
                        double[] rudy = new double[ruets.Count()];
                        double[] polyVal = new double[rudy.Length];

                        int j = 0;
                        foreach (double ruet in ruets)
                        {
                            rudy[j] = ruet;
                            polyVal[j] =
                                EvaluatePoly(
                                    rootElt.PolynomialFit,
                                    ruet);
                            j++;
                        }
                        rootElt.RootIndex = rudy;
                        rootElt.PolyRootValue = polyVal;
                    };


                    rootList.Add(rootElt);


                }
            }

            int k = 0;
            DiscreteRoot[] rootArr = new DiscreteRoot[rootList.Count()];
            foreach (DiscreteRoot root in rootList)
            {
                rootArr[k] = root;
                k++;
            }









        }

    }
}

