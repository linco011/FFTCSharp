# FFTCSharp
This is a libary that I have partially developed to perform Fast fourier transform on any data set. This is a C# draft of the project.
I plan on moving it to a lower level language so it could be implemented on embeded software devices. 
This differs from the standard fast fourier transform that has the underlying assumption that the number of data points is 2^n. Instead I 
follow cooley-tukey's orginal paper and expand it to any composite number of data points. 
This requires prime factor decomposition of N points of data to apply the fastest transformation. Prime factor decomposition is of O(N)
making it fast enough to justify not using the naive O(N^2) method of dft. 
Cooley-tukey present that you can break down a DFT as O(N(r1+r2+...+rl)) where {ri} are the set of non-distinct prime numbers composing N making it in general
faster then O(N^2) for all composite N. Only when N is prime would it be just as slow as standard DFT.  
