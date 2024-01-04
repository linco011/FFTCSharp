using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class Polynomial
	{
		/*
		 * This class signifies
		 * the object is a polynomial
		 * to make it easier for the user
		 * to understand what is going on. 
		 */

		// this is the degree of the polynomial (size of array - 1)
		public int Degree { get; set; }

        // this is the array where 0 is const or x^0, 1 is x^1,...
        public double[] Array { get; set; }
	}
}

