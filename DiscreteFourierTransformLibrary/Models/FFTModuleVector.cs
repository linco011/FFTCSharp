using System;
namespace DiscreteFourierTransformLibrary.Models
{
	public class FFTModuleVector
	{
		/*
		 * Instead of having a full list of FFTModuleIndex as the
		 * key. It would be better to have the wave number 
		 * seperate from the the j values.
		 * 
		 * This will allow us to calculate 
		 * the SubJ more efficiently without 
		 * having to calculate the entire thing from
		 * scratch for each JVector calculated
		 * 
		 * 
		 */


		public FFTModuleIndex WaveNumber { get; set; }

		public Int32 SubJ { get; set; }

		public Int32 Iteration { get; set; }

		public List<FFTModuleIndex> JVector { get; set:}


	}
}

