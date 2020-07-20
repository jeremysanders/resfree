#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "constants.h"
#include "linefrac.h"

/* xwrite - write output with a certain "chatter" level */
PROTOCCALLSFSUB2(XWRITE, xwrite, STRING, INT)

/* like printf, but does an xwrite instead */
void reson_xwrite(int level, const char* format, ...)
{
  char buffer[256];
    
  va_list arglist;
  va_start( arglist, format );
  vsnprintf( buffer, sizeof buffer, format, arglist );

  CCALLSFSUB2(XWRITE, xwrite, STRING, INT, (char*)buffer, level);

  va_end( arglist );
}

/* get filter for dataset in the form of a double */
PROTOCCALLSFFUN2(STRING, DGFILT, dgfilt, INT, INT)
double reson_get_dataset_filter( const unsigned dataset,
				 const unsigned no )
{
  const char* text = CCALLSFFUN2(DGFILT, dgfilt, INT, INT, no, dataset);
  double d;
  char* processed;

  if( text[0] == 0 )
    {
      reson_xwrite(5, "Error: reson: Filter %i not found in dataset %i",
		   no, dataset);
      return 1;
    }

  d = strtod( text, &processed );

  if( processed == text )
    {
      reson_xwrite(5, "Error: reson: Invalid number '%s' in filter %i,"
		   " dataset %i", text, no, dataset);
      return 1;
    }
  else
    {
      reson_xwrite(15, "Info: reson: Returned %f as filter %i from dataset %i",
		   d, no, dataset);
      return d;
    }
}


/* calculate luminosity distance (kpc) to object at redshift z */
PROTOCCALLSFFUN3(FLOAT, FZSQ, fzsq, FLOAT, FLOAT, FLOAT)
double calc_lumin_dist_kpc(const float z)
{
  const float q0 = get_q0();
  const float lambda0 = get_lambda0();
  const float H0 = get_H0();

  const double fz2 = CCALLSFFUN3(FZSQ, fzsq, FLOAT, FLOAT, FLOAT,
				 z, q0, lambda0);

  return sqrt(fz2) * (C_KM_S/H0) * 1000.;
}

/*
  Work out volume of shell (R1->R2) intersecting LOS (y1->y2)
  
  this is the integral:
   Int(y=y1,y2) Int(x=sqrt(R1^2-y^2),sqrt(R2^2-y^2)) 2*pi*y dx dy
   =
   Int(y=y1,y2) 2*pi*y*( sqrt(R2^2-y^2) - sqrt(R1^2-y^2) ) dy
   
   volume returned is half total volume (front only)
   
   routine below tested using monte-carlo integration
*/

double projection_vol(const double R1, const double R2,
		      const double y1, const double y2)
{
  const double y1_2 = y1*y1;
  const double y2_2 = y2*y2;
  const double R2_2 = R2*R2;
  const double R1_2 = R1*R1;

  const double p1 = trunc_sqrt( R1_2 - y2_2 );
  const double p2 = trunc_sqrt( R1_2 - y1_2 );
  const double p3 = trunc_sqrt( R2_2 - y2_2 );
  const double p4 = trunc_sqrt( R2_2 - y1_2 );

  return (2./3.) * M_PI * ( (p1*p1*p1 - p2*p2*p2) + (p4*p4*p4 - p3*p3*p3) );

/*   const double sqrt_R1_2_m_y1_2 = trunc_sqrt(R1_2 - y1_2); */
/*   const double sqrt_R1_2_m_y2_2 = trunc_sqrt(R1_2 - y2_2); */
/*   const double sqrt_R2_2_m_y1_2 = trunc_sqrt(R2_2 - y1_2); */
/*   const double sqrt_R2_2_m_y2_2 = trunc_sqrt(R2_2 - y2_2); */

/*   return( (2./3.) * M_PI * */
/*           ( R2_2 * (sqrt_R2_2_m_y1_2 - sqrt_R2_2_m_y2_2) + */
/*             R1_2 * (sqrt_R1_2_m_y2_2 - sqrt_R1_2_m_y1_2) + */
/*             y2_2 * (sqrt_R2_2_m_y2_2 - sqrt_R1_2_m_y2_2) + */
/*             y1_2 * (sqrt_R1_2_m_y1_2 - sqrt_R2_2_m_y1_2) ) */
/*           ); */

}

/* insert absorption line into spectrum */
void insert_absorb_line(float* const spec, const double coeff,
			const double energy_keV, const double width_keV)
{
  const double w1 = 1./width_keV;

  /* where the line is at the centre */
  const int centre_chan = ENERGY_TO_BIN( energy_keV );

  int chan;

  chan = centre_chan;
  /* positive part of the spectrum */
  while( chan < ENERGY_NOBINS )
    {
      const double x = ((chan*ENERGY_SPACING + ENERGY_SPACING*0.5 +
			 ENERGY_MIN) - energy_keV) * w1;

      const double e2 = exp2_interpolate( x );
      if( e2 < 0. )
	break;

      if( chan >= 0 )
	spec[chan] += (e2*coeff);
      ++chan;
    }

  chan = centre_chan - 1;
  /* negative part of the spectrum */
  while( chan >= 0 )
    {
      const double x = ( energy_keV -
			 ( chan*ENERGY_SPACING + ENERGY_SPACING*0.5 +
			   ENERGY_MIN )) * w1;

      const double e2 = exp2_interpolate( x );
      if( e2 < 0. )
	break;

      if( chan < ENERGY_NOBINS )
	spec[chan] += (e2*coeff);
      --chan;
    }
}

/* convert spectrum from one energy spacing to another, interpolating
   if necessary */
void regrade_spectrum(const unsigned ne_in,
		      const float* const ear_in,
		      const float* const spec_in,
		      const unsigned ne_out,
		      const float* const ear_out,
		      float* const spec_out)
{
  unsigned ne_in_step;
  unsigned ne_out_step;
  unsigned i;

  ne_in_step = 0;
  /* move to first applicable bin in the input corresponding to the output */
  while( ne_in_step < ne_in && ear_in[ne_in_step] < ear_out[0] )
    ++ne_in_step;

  /* blank output array */
  for(i = 0; i != ne_out; ++i)
    spec_out[i] = 0.;

  /* move to first applicable bin in the output corresponding to the input */
  /* FIXME: there's probably a problem on the first bin... */
  ne_out_step = 0;
  while( ne_out_step < ne_out && ear_out[ne_out_step] < ear_in[ne_in_step] )
    ++ne_out_step;

  /* step over each input bin and calculate contribution to output */
  while( ne_in_step < ne_in && ne_out_step < ne_out )
    {
      /* find the overlapping energy bounds */
      const double in1 = ear_in[ne_in_step];
      const double in2 = ear_in[ne_in_step+1];
      const double out1 = ear_out[ne_out_step];
      const double out2 = ear_out[ne_out_step+1];

      const double overlap1 = MAX( in1, out1 );
      const double overlap2 = MIN( in2, out2 );

      /* add fractional contribution of input */
      spec_out[ ne_out_step ] +=
	( spec_in[ne_in_step] * ( (overlap2-overlap1)/(in2-in1) ) );
      
      /* go to next appropriate bin */
      if( in2 <= out2 )
	++ne_in_step;

      if( in2 >= out2 )
	++ne_out_step;
    }
}
