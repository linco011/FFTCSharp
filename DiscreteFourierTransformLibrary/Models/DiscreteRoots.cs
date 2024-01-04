using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class DiscreteRoot
	{
		/*
		 * This class contains
		 * the point to right and left
		 * of the root.
		 * 
		 * The polynomial fitted to find 
		 * the root 
		 * 
		 * And the root that it found
		 */

		// This list the indicies close to the root(s) 
		public int[] SampleIndices { get; set; }

		// The sample values of the array
		public double[] DataValues { get; set; }

		// the approximate index of the root
		public double[]? RootIndex { get; set; }

		// the polynomial used to find the root
		public Polynomial PolynomialFit { get; set; }

		// the value the polynomial gets at the root
		// SHOULD ALWAYS BE ZERO
		public double[]? PolyRootValue { get; set; }
	}
}

