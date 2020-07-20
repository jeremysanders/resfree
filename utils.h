#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <cfortran.h>
#include <assert.h>

#include "constants.h"

#ifdef __GNUC__
 /* define some functions inline if using gcc for speed */
# define INLINE inline
 /* check formatting for printf-style functions with gcc */
# define CHECKFORMAT(a,b,c) __attribute__ ((format (a, b, c)))
#else
# define INLINE
# define CHECKFORMAT(a,b,c)
#endif

/* the usual macros - are these defined somewhere else? */
#define MIN(x,y) ( ( (x) < (y) ) ? (x) : (y) )
#define MAX(x,y) ( ( (x) > (y) ) ? (x) : (y) )

/* obvious? */
#define SQUARE(x) ( (x)*(x) )
#define CUBE(x) ((x)*(x)*(x))

static INLINE double sqr(const double x)
{
  return x*x;
}

static INLINE double cube(const double x)
{
  return x*x*x;
}

/* these macros give warnings if the types are wrong */
#define ALLOCSAFE(type, no) (type*)( malloc((no)*sizeof(type)) )
#define REALLOCSAFE(ptr, type, no) (type*)( realloc(ptr, (no)*sizeof(type)) )

/* whether numbers are close */
static INLINE int near_double(double a, double b)
{
  if( fabs(a) < 1e-10 && fabs(b) < 1e-10 )
    return 1;

  return fabs((a-b)/a) < 1e-5;
}

/* a checking wrapper for realloc */
static INLINE void float_realloc(float** const ptr, const unsigned noitems)
{
  *ptr = REALLOCSAFE( *ptr, float, noitems ); assert( *ptr != NULL );
}

/* a checking wrapper for realloc */
static INLINE void double_realloc(double** const ptr, const unsigned noitems)
{
  *ptr = REALLOCSAFE( *ptr, double, noitems ); assert( *ptr != NULL );
}

/* call the apecrs model */
PROTOCCALLSFSUB6(APECRS, apecrs, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
static INLINE void F_apecrs(const float* ear, int ne, const float* param,
			    int ifl, float* photar, float* photer)
{
  CCALLSFSUB6(APECRS, apecrs, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV,
	      (float*)ear, ne, (float*)param, ifl, photar, photer);
}

/* call the zphabs model */
PROTOCCALLSFSUB6(XSZPHB, xszphb, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
static INLINE void F_zphabs(const float* ear, int ne, const float* param,
			    int ifl, float* photar, float* photer)
{
  CCALLSFSUB6(XSZPHB, xszphb, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV,
	      (float*)ear, ne, (float*)param, ifl, photar, photer);
}

/* get q0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGQ0, csmgq0)
static INLINE float get_q0(void)
{
  return CCALLSFFUN0(CSMGQ0, csmgq0);
}

/* get lambda0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGL0, csmgl0)
static INLINE float get_lambda0(void)
{
  return CCALLSFFUN0(CSMGL0, csmgl0);
}

/* get H0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGH0, csmgh0)
static INLINE float get_H0(void)
{
  return CCALLSFFUN0(CSMGH0, csmgh0);
}

/* calculate luminosity distance (kpc) to object at redshift z */
extern double calc_lumin_dist_kpc(const float z);

/* get the angular diameter distance (kpc) to object at redshift z */
static INLINE double calc_ang_diam_dist_kpc(const float z)
{
  const float lumin_dist = calc_lumin_dist_kpc(z);
  const float z1 = 1. + z;
  return lumin_dist / (z1*z1);
}

/* get number of datasets */
PROTOCCALLSFFUN0(INT, DGNDST, dgndst)
static INLINE unsigned get_no_datasets(void)
{
  return CCALLSFFUN0(DGNDST, dgndst);
}

/* check whether calculation is forced */
PROTOCCALLSFFUN0(LOGICAL,XGFCLC,xgfclc)
static INLINE int get_force_calc(void)
{
  return CCALLSFFUN0(XGFCLC, xgfclc);
}

/* force individual calculation of spectra in datasets */
PROTOCCALLSFSUB1(XPFCLC, xpfclc, LOGICAL)
static INLINE void set_force_calc(int OnorOff)
{
  CCALLSFSUB1(XPFCLC, xpfclc, LOGICAL, OnorOff);
}

/* call xspec xwrite but with printf style formatting */
extern void reson_xwrite(int level, const char* format, ...)
     CHECKFORMAT(printf, 2, 3);

/* get filter for dataset */
extern double reson_get_dataset_filter( const unsigned dataset,
					const unsigned no );

/* return sqrt of a if a > 0 else return 0 */
static INLINE double trunc_sqrt(const double a)
{
  return a > 0. ? sqrt(a) : 0.;
}

/* work out volume of intersection of shell radius R1->R2 to
   shell on sky of radius y1 and y2.

   volume returned is the front only */
extern double projection_vol(const double R1, const double R2,
			     const double y1, const double y2);

/* calculate the length between two shells of R1 and R2
   at projected radius r */
static INLINE double projection_len(const double R1, const double R2,
				    const double r)
{
  const double r_2 = r*r;
  return trunc_sqrt( R2*R2 - r_2 ) - trunc_sqrt( R1*R1 - r_2 );
}

/* insert absorption line into spectrum */
extern void insert_absorb_line(float* const spec,
			       const double coeff,
			       const double energy_keV,
			       const double width_keV);

/* convert spectrum from one energy spacing to another, interpolating
   if necessary */
extern void regrade_spectrum(const unsigned ne_in,
			     const float* const ear_in,
			     const float* const spec_in,
			     const unsigned ne_out,
			     const float* const ear_out,
			     float* const spec_out);


/* blank elements in spectrum */
static INLINE void blank_spectrum(float* const spec)
{
  unsigned i;
  for( i=0; i != ENERGY_NOBINS; ++i )
    spec[i] = 0;
}

/* add together elements in spectra */
static INLINE void add_spectrum( float* const spec1,
				 const float* const spec2 )
{
  unsigned i;
  for( i=0; i != ENERGY_NOBINS; ++i )
    spec1[i] += spec2[i];
}

/* add on a multiple of another spectrum */
static INLINE void add_mult_spectrum(float* const spec1,
				     const float* const spec2,
				     const float factor)
{
  unsigned i;
  for( i=0; i != ENERGY_NOBINS; ++i )
    spec1[i] += spec2[i]*factor;
}

#endif
