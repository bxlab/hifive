# include <cmath>
# include "_normal.hpp"
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

#define bool int
#define true 1
#define false 0

double stdnorm_cdf(double z)

//****************************************************************************80
//
// This method calculates the normal cumulative distribution function.
//
// It is based upon algorithm 5666 for the error function, from:
//
// Hart, J.F. et al, 'Computer Approximations', Wiley 1968
//
// The FORTRAN programmer was Alan Miller.  The documentation
// in the FORTRAN code claims that the function is "accurate
// to 1.e-15."
// Steve Verrill translated the FORTRAN code (the March 30, 1986 version)
// into Java.  This translation was performed on January 10, 2001.
//
// The method returns the value of the normal
// cumulative distribution function at z.
//
// @version .5 --- January 10, 2001
//
// Here is a copy of the documentation in the FORTRAN code:
//
// SUBROUTINE NORMP(Z, P, Q, PDF)
// Normal distribution probabilities accurate to 1.e-15.
// Z = no. of standard deviations from the mean.
// P, Q = probabilities to the left & right of Z.	P + Q = 1.
// PDF = the probability density.
//
// Based upon algorithm 5666 for the error function, from:
// Hart, J.F. et al, 'Computer Approximations', Wiley 1968
//
// Programmer: Alan Miller
// Latest revision - 30 March 1986
{
	double zabs;
	double p;
	double expntl,pdf;

	static double p0 = 220.2068679123761;
	static double p1 = 221.2135961699311;
	static double p2 = 112.0792914978709;
	static double p3 = 33.91286607838300;
	static double p4 = 6.373962203531650;
	static double p5 = .7003830644436881;
	static double p6 = .3526249659989109E-01;

	static double q0 = 440.4137358247522;
	static double q1 = 793.8265125199484;
	static double q2 = 637.3336333788311;
	static double q3 = 296.5642487796737;
	static double q4 = 86.78073220294608;
	static double q5 = 16.06417757920695;
	static double q6 = 1.755667163182642;
	static double q7 = .8838834764831844E-1;

	static double cutoff = 7.071;
	static double root2pi = 2.506628274631001;

	zabs = std::abs(z);

	//  |z| > 37

	if (z > 37.0) {
		p = 1.0;
		return p;
	}

	if (z < -37.0) {
		p = 0.0;
		return p;
	}

	//  |z| <= 37.

	expntl = exp(-.5*zabs*zabs);
	pdf = expntl/root2pi;

	//  |z| < cutoff = 10/sqrt(2).

	if (zabs < cutoff) {
		p = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs +
			 p2)*zabs + p1)*zabs + p0)/(((((((q7*zabs + q6)*zabs +
			 q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +
			 q0);
	} else {
		p = pdf/(zabs + 1.0/(zabs + 2.0/(zabs + 3.0/(zabs + 4.0/
			 (zabs + 0.65)))));
	}

	if (z < 0.0) {
		return p;
	} else {
		p = 1.0 - p;
		return p;
	}
}



//**************************************************

double stdnorm_pdf( double z )
{
  double value;
  value = 0.3989422804014327 * exp (-0.5 * pow(z,2.0));
  return value;
}