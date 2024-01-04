using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class PolyCritical
	{
		public Polynomial Polynomial { get; set; }

		public double RangeMin { get; set; }

		public double RangeMax { get; set; }

		public double Convergence { get; set; }

		public double[]? MinXArr { get; set; }

		public double[]? MinYArr { get; set; }

		public double[]? MaxXArr { get; set; }

		public double[]? MaxYArr { get; set; }

		public double[]? xSaddlePoints { get; set; }

		public double[]? ySaddlePoints { get; set; }

		public double[]? xInflectionPoints { get; set; }

		public double[]? yInflectionPoints { get; set; }


		
	}
}

