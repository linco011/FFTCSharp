using System;
using DiscreteFourierTransformLibrary.Models;
using System.Collections.Generic;
using DiscreteFourierTransformLibrary.Modules;
namespace DiscreteFourierTransformLibrary.Modules
{
	public class FindingPeriod
	{
		/*
		 * This is will find the period 
		 * for a given input of data. 
		 * 
		 * Algo: 
		 *  1.) FFT transform of data
		 *  2.) sort the plane wave contribution by
		 *  intensity
		 *  3.) Find the period for a given number of contributing
		 *  plane waves
		 *		3a.) find the gcf of the wavenumbers 
		 *		3b.) Number of data points/gcf= period
		 *	4.) Create an object for each index that is period (y axis)
		 *	versus accumulative total (x axis)
		 *	5.) graph should be concave down and will increase, plateau then
		 *	there will be a point where the graph will become concave up.
		 *  [If true, finding the inflection point will be the most important]
		 *  6.) Establish inflection point as the period
		 *  7.) get the number of points (N) nearest period (T)
		 *  8.) Since period is rational find the smallest x in integers
		 *  such that x,y are integers and xN=yT
		 *  9.) use x and N to define your comparison array and push that
		 *  as a return value
		 *  
		 *  Next Algo requirements:
		 *  1.) Need to calculate the variance and mean of comparison array
		 *  2.) define parameter of what is considered an anomally 
		 *  (value-mean>sqrt(variance))
		 *  3.) Define the number of sequential anomallies before reporting
		 *  an anomaly
		 *  4.) report anomaly after differences.
		 */

		// This is for real input (which is the case for most data other
		// than electrical)
		public PeriodGraph[] PeriodDistribution(
			double[]Input
			)
		{
            /*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * This takes in real data from
            * a source and converts the data into
            * its Fourier Transformed components
            * Organizes the graph from highest 
            * intensity to lowest intensity.
			* Finally will create a PeriodGraph[]
			* object that will calculate the common
			* period and cumulative total of signal
			* from the beginning of the distribution
			* to the end of the distribution
			* 
			* where the common period i will have
			* resonant freguencies from the zeroth entry
			* to the ith entry of interest.
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

            // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
            // Set of Info log statements
            // Set of Error Statements

            // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

            System.Numerics.Complex[] compInput =
				new System.Numerics.Complex[Input.Length];
			for(int i=0; i<Input.Length; i++)
			{
				compInput[i] = new System.Numerics.Complex
					(
						Input[i]
						,0.0
					);


			};

			System.Numerics.Complex[] fftInput=
				MathFunctions.FFTCalculation(compInput);

			GraphDFT[] unsorted =
				MathFunctions.GraphingDFT(fftInput);

			GraphDFT[] sorted =
				MathFunctions.GraphSort
				(
					unsorted,
					0,
					unsorted.Length - 1
				);
			NormGraphDFT[] normGraph =
				MathFunctions.NormalizeGraphDFT(sorted);

			PeriodGraph[] periodGraph =
				MathFunctions.CalculateCommonTime(normGraph);

			return periodGraph;

		}

		public DiscreteRoot[] FindPeriodInflections(
			PeriodGraph[]graph,
            double nearRoot,
            double rootConvergence
            )
        {
			/*
            * SYNOPSIS: 
            * 
            * DESCRIPTION:
            * The graph is expected to
            * from the period graph is expeccted
            * to have a period graph that will monotonically
            * increase. 
            *
			* If the system we are analyzing is expected to
			* have a variant steady state, then the steady state
			* should have a resonant period that contains most
			* of the information. This will imply there
			* will be a period will be plateued until 
			* the intensity gets small enough that we enter 
			* the noise range which would not have any useful
			* data. It is expected that the data will jump
			* where the period will be of the length of the whole
			* system.
			* 
			* In this case the period will have at least one point
			* where the concavity becomes concave down to concave up
			* Hence we take the second derivative of the data and 
			* find the roots of that data. The first root 
			* should be where the data has the most common period 
			* without noise which will be used for our calculations.
			* 
            * 
            * DEPENDENCIES:
            * 
            * OUTPUTS:
            *
            * PARAMETER:
            * 
            * EXAMPLE:
            */

			// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL Log & Errors
			// Set of Info log statements
			// Set of Error Statements

			// MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Method Starts

			// First we look at the system as a simple array.
			double[] logPeriodData = new double[graph.Length];
			for(int i=0; i<graph.Length;i++)
			{
				logPeriodData[i] = graph[i].LogPeriod;
			}

			// Find First Derivative dT/T
			double[] firstDer =
				MathFunctions.SimpleDerivative(logPeriodData);
			// Second Derivative (inflection should occur here)
			double[] secondDer =
				MathFunctions.SimpleSecondDerivative(
					logPeriodData,
					firstDer
					);

			DiscreteRoot[] roots =
				MathFunctions.FindDiscreetRoots(
					secondDer,
					nearRoot,
					rootConvergence
					);


			return roots;




        }




    }
}

