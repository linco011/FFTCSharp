using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class ModularElement
	{
		/*
		 * Euclidean division can be 
		 * generalized to the entirity 
		 * of the real numbers
		 * 
		 * OriginalInput = ModularDivisor*IntegerMultiple + Remainder
		 * 
		 * where 0 <= Remainder < ModularDivisor
		 * 
		 * This is equivelant of 
		 * Remainder = OriginalInput % ModularDivisor
		 * 
		 * 
		 */

		public double OriginalInput { get; set; }

		public double ModularDivisor { get; set; }

		public double IntegerMultiple { get; set; }

		public double Remainder { get; set; }


	}
}

