using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class DerivativesPeriodGraph
	{
		PeriodGraph PeriodGraph { get; set; }

		double PerFirstDerivative { get; set; }

		double PerSecondDerivative { get; set; }

		double AccFirstDerivative { get; set; }

		double AccSecondDerivative { get; set; }

		double LogFirstDerivative { get; set; }

		double LogSecondDerivative { get; set; }


	}
}

