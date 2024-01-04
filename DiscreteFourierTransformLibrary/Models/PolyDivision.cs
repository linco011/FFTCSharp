using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class PolyDivisionOutput
	{
		/*
		 * n(x) = q(x)d(x) + r(x)
		 * n=Numerator 
		 * q = quotient
		 * d = denominator/divisor
		 * r = remainder
		 * 
		 * deg r(x) < deg d(x)
		 * 
		 * 
		 */

		public Polynomial Quotient { get; set; }

		public Polynomial Numerator { get; set; }

		public Polynomial Divisor { get; set; }

		public Polynomial Remainder { get; set; }
	}
}

