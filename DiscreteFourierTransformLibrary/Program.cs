using System;
using System.Collections.Generic;
using DiscreteFourierTransformLibrary.Models;
using DiscreteFourierTransformLibrary.Modules;

int x = 428490000; //3^4*2^4*5^4*23^2

MathFunctions myMath = new MathFunctions();

var listFactor = myMath.FindPrimeFactors(x);

foreach(PrimeFactor element in listFactor)
{
    Console.WriteLine("prime number is ");
    Console.WriteLine(element.Numb);
    Console.WriteLine("with the power bellow");
    Console.WriteLine(element.Pow);

}
