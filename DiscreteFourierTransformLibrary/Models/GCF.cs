using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class GCF
	{
		/*
		 * This is the greatest 
		 * common factor class for a given set
		 * of numbers
		 */
		public int[] InputNumbers { get; set; }

		public List<PrimeFactor> Primes { get; set; }

		public int Factor { get; set; }
	}
}

