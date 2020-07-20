/* return fraction of line between two energies */

#include <math.h>


/*
  generated using:

#include <stdio.h>
#include <gsl/gsl_cdf.h>
  
int main()
{
  int i;
  for(i=0; i <= 60; ++i)
    {
      double x = i/60. * 3.;
      double fn = gsl_cdf_gaussian_Q(x, 1)*2;

      printf("%e %.10e\n ", x, fn);
    }

  printf("\n");
}
*/

#define GAUSS_NO_STEPS 60
#define GAUSS_MAX 3.
#define GAUSS_STEP (GAUSS_MAX / GAUSS_NO_STEPS)

static const double gauss_lookup[GAUSS_NO_STEPS+1] =
  {
    1.0000000000e+00, 9.6012238832e-01, 9.2034432545e-01, 8.8076461526e-01,
    8.4148058112e-01, 8.0258734863e-01, 7.6417715562e-01, 7.2633869765e-01,
    6.8915651678e-01, 6.5271044058e-01, 6.1707507745e-01, 5.8231937358e-01,
    5.4850623550e-01, 5.1569222161e-01, 4.8392730445e-01, 4.5325470475e-01,
    4.2371079717e-01, 3.9532508625e-01, 3.6812025069e-01, 3.4211225262e-01,
    3.1731050786e-01, 2.9371811275e-01, 2.7133212189e-01, 2.5014387127e-01,
    2.3013934044e-01, 2.1129954733e-01, 1.9360096917e-01, 1.7701598287e-01,
    1.6151331847e-01, 1.4705851922e-01, 1.3361440254e-01, 1.2114151600e-01,
    1.0959858340e-01, 9.8942936067e-02, 8.9130925517e-02, 8.0118313728e-02,
    7.1860638226e-02, 6.4313549591e-02, 5.7433119632e-02, 5.1176119043e-02,
    4.5500263896e-02, 4.0364430811e-02, 3.5728841126e-02, 3.1555214782e-02,
    2.7806895027e-02, 2.4448945310e-02, 2.1448220043e-02, 1.8773411070e-02,
    1.6395071849e-02, 1.4285621471e-02, 1.2419330652e-02, 1.0772291908e-02,
    9.3223760474e-03, 8.0491770855e-03, 6.9339476061e-03, 5.9595264701e-03,
    5.1102606609e-03, 4.3719229098e-03, 3.7316266008e-03, 3.1777392947e-03,
    2.6997960633e-03
  };

double line_frac(const double linecentre,
		 const double linewidth,
		 const double energy)
{
  const double sigma = fabs( (energy-linecentre)/linewidth );

  const unsigned minindex = (unsigned)( sigma * (1./GAUSS_STEP) );

  double remainder;
  double inter;

  /* if we go past the end of the array, return -1 */
  if( minindex >= GAUSS_NO_STEPS )
    return -1;

  remainder = (sigma - minindex*GAUSS_STEP) * (1./GAUSS_STEP);

  /* interpolate */
  inter = (1.-remainder)*gauss_lookup[minindex] +
    remainder*gauss_lookup[minindex+1];

  return 1.-inter;
}

#define EXP_NO_STEPS 60
#define EXP_MAX 3.
#define EXP_STEP (EXP_MAX / EXP_NO_STEPS)

/* generated using
from math import *
import sys

nosteps=60
max=3.
step=max/nosteps

for i in range(nosteps+1):
    x = i*step
    sys.stdout.write( "%e, " % ( exp( - x*x ), ) )
*/

static const double exp2_lookup[EXP_NO_STEPS+1] =
  {
    1.000000e+00, 9.975031e-01, 9.900498e-01, 9.777512e-01, 9.607894e-01,
    9.394131e-01, 9.139312e-01, 8.847059e-01, 8.521438e-01, 8.166865e-01,
    7.788008e-01, 7.389685e-01, 6.976763e-01, 6.554063e-01, 6.126264e-01,
    5.697828e-01, 5.272924e-01, 4.855369e-01, 4.448581e-01, 4.055545e-01,
    3.678794e-01, 3.320399e-01, 2.981973e-01, 2.664683e-01, 2.369278e-01,
    2.096114e-01, 1.845195e-01, 1.616212e-01, 1.408584e-01, 1.221507e-01,
    1.053992e-01, 9.049144e-02, 7.730474e-02, 6.571027e-02, 5.557621e-02,
    4.677062e-02, 3.916390e-02, 3.263076e-02, 2.705185e-02, 2.231491e-02,
    1.831564e-02, 1.495813e-02, 1.215518e-02, 9.828195e-03, 7.907054e-03,
    6.329715e-03, 5.041760e-03, 3.995846e-03, 3.151112e-03, 2.472563e-03,
    1.930454e-03, 1.499685e-03, 1.159229e-03, 8.915937e-04, 6.823281e-04,
    5.195747e-04, 3.936690e-04, 2.967858e-04, 2.226299e-04, 1.661699e-04,
    1.234098e-04
  };

/* interpolate the function exp(-x*x) */

double exp2_interpolate(const double x)
{
  const double fabsx = fabs(x);
  const unsigned minindex = (unsigned) ( fabsx * (1./EXP_STEP) );
  double remainder;

  if( minindex >= EXP_NO_STEPS )
    return -1;

  remainder = ( fabsx - minindex*EXP_STEP) * (1./EXP_STEP);

  /* interpolate */
  return (1.-remainder)*exp2_lookup[minindex] + 
    remainder*exp2_lookup[minindex+1];
}

/* #include <stdio.h> */
/* int main() */
/* { */
/*   double e; */
/*   double last = 0; */

/*   for( e=0; e<3.; e += 0.0001 ) */
/*     { */
/*       double f = line_frac(0, 1, e); */
/*       printf("%f %f\n", e, f-last); */
/*       last = f; */
/*     } */
/*   return 0; */
/* } */
