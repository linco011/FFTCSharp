using System;
namespace DiscreteFourierTransformLibrary.Models
{
    /*
     * When we decompose the Discrete Fourier Transform
     * into N= r0r1r2...rl using the Cooley-Tukey method, 
     * we will be using the unique factorization of the integer
     * N as each r0, r1, r2...
     * 
     * This will mean that r0, r1... are all not necessarilly
     * distinct prime numbers. 
     * 
     * to distinguish which power of the set we are representing for a 
     * given PrimeFactor prime, we need an Int32 Index to represent
     * it to distinguish ri from ri+1 eventhough ri and ri+1 might be
     * the same prime
     * 
     * There is auxilary information here as well that might be useful
     * such as TotalComposite = N
     * 
     * SubComposite = r0r1r2...rk
     * 
     * LeftOverComposite = N/(r0r1r2...rk)= rk+1rk+2...rl
     * 
     * rk= PrimeFactorIndex.Prime.Numb
     * k = PrimeFactorIndex.Index
     * where k ranges from 0 to l.
     *
     * SubIndex = 1,2,...,PrimeFactorIndex.Prime.Pow
     * 
     * The value of rk=PrimeFactorIndex.Numb
     * 
     * 
     * 
     */
    public class PrimeFactorIndex
    {
        
        public PrimeFactor Prime { get; set; }

        public Int32 Index { get; set; }

        public Int32 SubIndex { get; set; }

        public Int32 TotalComposite {get;set;}

        public Int32 SubComposite { get; set; }

        public Int32 LeftOverComposite { get; set; }
    }
}

