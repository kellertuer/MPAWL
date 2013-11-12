# The multivariate periodic anisotropic Wavelet Library

This `Mathematica` Library is an implementation of the periodic Wavelet Transform based on an integral regular matrix __M__ and its factorization into dilation matrices. Introducing the multivariate de la Vallée Poussin means, this Library provides many scaling functions for the levels of decomposition.

The underlying theory of patterns, its generating groups and the fast Fourier transform on these patterns is also implemented in this Library yielding a fast Wavelet Transform, when computing in Fourier coefficients. Further, several functions to work with the translates of a function with respect to the pattern

For the dyadic case, i.e. |det __J__<sub>k</sub>|=2 for all factor matrices, this Library also provides the construction of corresponding wavelets and an algorithm to decompose on several different factorizations at the same time.

Further, for the dyadic two-dimensional case, several visualization methods are given for the pattern, the wavelet and scaling functions and the obtained fractions of a function sampled on a pattern and decomposed with respect to the wavelets.

Several examples illustrate most of the implemented functions, each of which is equipped with an detailed `::usage` command.

## Also available on
This package is also available in the Wolfram Library Archive, see http://library.wolfram.com/infocenter/MathSource/8761/ 

## Software used in this package
The Smith normal form package `SmithFormV6.m` was written by Adriano Pascoletti and is used with his permission. The package can be obtained at [http://library.wolfram.com/infocenter/MathSource/7081/](http://library.wolfram.com/infocenter/MathSource/7081/).

The package to evaluate a Box spline is a transcription of MatLab implementation of `Efficient and stable evaluation of box-splines` written by Leif Kobbelt, which is available at [http://www.netlib.org/numeralgo/na11](http://www.netlib.org/numeralgo/na11).

## License
    MPAWL is free software : you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    MPAWL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with the MPAWL. If not, see <http://www.gnu.org/licenses/>.

