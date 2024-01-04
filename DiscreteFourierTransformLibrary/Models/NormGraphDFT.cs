using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class NormGraphDFT
	{
		public int Order {get; set;}

		public GraphDFT OriginalGraph { get; set; }

		public double Norm { get; set; }

		public double NormIntensity { get; set;}

		public System.Numerics.Complex NormAmplitude { get; set; }
	}
}

