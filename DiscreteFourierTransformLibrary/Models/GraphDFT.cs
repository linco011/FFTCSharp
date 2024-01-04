using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class GraphDFT
	{
		public int OriginalIndex { get; set; }

		public int AnalyticalIndex { get; set; }

		public double WaveNumber { get; set; }

		public System.Numerics.Complex Amplitude { get; set; }

		public double Intensity { get; set; }
	}
}

