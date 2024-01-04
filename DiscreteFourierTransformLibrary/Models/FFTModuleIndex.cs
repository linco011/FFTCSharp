using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class FFTModuleIndex
	{
		
        /*
		* For Tracking purposes
		* we need to have what the range
		* of our Index for our OneStepFFT
		* is.
		* for a given 
		* Ai(k,ji-1,ji-2,...,j1,j0)
		* 0<=j0<r1
		* 0<=j1<r2
		* ...
		* 0<=ji-1<ri
		* 0<=k<N/(r1r2...ri)
		* 
		* k and {ji} create
		* a list that can be a key
		* 
		* This class should be constructed
		* in a way that correlates
		* to PrimeFactorIndex
		* 
		* 
		* 
		*/

		// This is the Value used as the key
		public Int32 Value { get; set; }

		// Value can only be allowed to be < Bound
		// Bound=r1 or N/(r1r2...ri)
		public Int32 Bound { get; set; }

		// If this is k then this value is null
		// Otherwise this will list the associated PrimeFactorIndex
		public PrimeFactorIndex? Prime { get; set; }

    }
}

