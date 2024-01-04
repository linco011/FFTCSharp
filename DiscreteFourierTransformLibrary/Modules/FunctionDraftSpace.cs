using System;
using System.Collections.Generic;
using DiscreteFourierTransformLibrary.Models;

namespace DiscreteFourierTransformLibrary.Modules
{
	public class FunctionDraftSpace
	{
        public Dictionary<FFTModuleVector, System.Numerics.Complex> OneStepFFT(
            Dictionary<FFTModuleVector, System.Numerics.Complex> preFFT,
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
            List<FFTModuleVector> correlatedInput = new List<FFTModuleVector>();
            List<FFTModuleVector> keys = new List<FFTModuleVector> ();
            Dictionary<FFTModuleVector, System.Numerics.Complex> postFFT
                = new Dictionary<FFTModuleVector, System.Numerics.Complex>();
            // This iterates through the various positions to produce
            foreach (FFTModuleVector pos in preFFT.Keys)
            {
                

                FFTModuleVector newPos = new FFTModuleVector();

                // part of the vecotor j without the jith entry.
                foreach (FFTModuleIndex jcomp in pos.JVector)
                {
                    newPos.JVector.Add(jcomp);
                }
                // This sets the new wavenumber for Ai+1
                newPos.WaveNumber.Value =
                    newPos.WaveNumber.Value % primes[stepi].LeftOverComposite;
                newPos.WaveNumber.Bound = primes[stepi].LeftOverComposite;
                newPos.WaveNumber.Prime = null;

                newPos.Iteration = pos.Iteration + 1;


                // This creates the jith index for each key. 
                for (Int32 ji = 0; ji < primes[stepi].Prime.Numb; ji++)
                {
                    var addji = new FFTModuleVector();
                    addji = newPos;
                    var jiComp = new FFTModuleIndex
                    {
                        Value = ji,
                        Prime = primes[stepi],
                        Bound = primes[stepi].Prime.Numb
                    };
                    addji.JVector.Add(jiComp);

                    addji.SubJ = pos.SubJ + ji * primes[stepi].SubComposite;

                    keys.Add(addji);
                    // this way we can call the
                    //corresponding input for Ai-1 to get Ai with
                    // the same key
                    correlatedInput.Add(pos);
                }




            }

            // This needs to be checked if that the keys line up correctly
            // To be added in this fashion.
            for(Int32 keyInd=0; keyInd < keys.Count; keyInd++)
            {
                System.Numerics.Complex totalValue
                    = System.Numerics.Complex.Zero;
                for (Int32 k = 0; k < primes[stepi].Prime.Numb; k++)
                {
                    totalValue += System.Numerics.Complex.Pow(W,
                        k * primes[stepi].LeftOverComposite
                        * keys[keyInd].SubJ) * preFFT[correlatedInput[keyInd]];
                }

                postFFT.Add(keys[keyInd],totalValue);


            }





        }
    }
}

