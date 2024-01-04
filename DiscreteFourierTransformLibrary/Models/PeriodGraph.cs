using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class PeriodGraph
	{
		// This is the complete array of points 
		public NormGraphDFT[] FullNormGraph { get; set; }

		// This is the array of points that give the
		// the Accumulated signal and the CommonPeriod
		public NormGraphDFT[] ConsideredNormGraph { get; set; }

		//The accumulated Signal is the sum of the NormIntensities
		// of the graph in ConsideredNormGraph
		public double AccumulatedSignal { get; set; }

		// This is the common period produced by the
		// Analytical indecies in ConsideredNormGraph
		// T=N/GCF(k1,k2,k3...,ki)
		public double CommonPeriod { get; set; }

		// This graphs the natural log of the period
		// We care about relative change of the period
		// which can be calculated as the difference of log period
		public double LogPeriod { get; set; }
	}
}

