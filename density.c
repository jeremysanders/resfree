#include <math.h>
#include <string.h>

#include <cfortran.h>

/* test model to fit density rather than norm */

#ifdef __GNUC__
 /* define some functions inline if using gcc for speed */
# define INLINE inline
 /* check formatting for printf-style functions with gcc */
# define CHECKFORMAT(a,b,c) __attribute__ ((format (a, b, c)))
#else
# define INLINE
# define CHECKFORMAT(a,b,c)
#endif

/* call xspec xwrite but with printf style formatting */
extern void reson_xwrite(int level, const char* format, ...)
     CHECKFORMAT(printf, 2, 3);

#define C_KM_S 299792.458
#define KPC_CM 3.08568025e21
#define NH_NE (1./1.2)
#define RADIAN_ARCSEC ((M_PI/180.)/60./60.)

/* get q0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGQ0, csmgq0)
INLINE static float get_q0()
{
  return CCALLSFFUN0(CSMGQ0, csmgq0);
}

/* get lambda0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGL0, csmgl0)
INLINE static float get_lambda0()
{
  return CCALLSFFUN0(CSMGL0, csmgl0);
}

/* get H0 from xspec */
PROTOCCALLSFFUN0(FLOAT, CSMGH0, csmgh0)
INLINE static float get_H0()
{
  return CCALLSFFUN0(CSMGH0, csmgh0);
}

/* calculate luminosity distance (kpc) to object at redshift z */
PROTOCCALLSFFUN3(FLOAT, FZSQ, fzsq, FLOAT, FLOAT, FLOAT)
static double calc_lumin_dist_kpc(const float z)
{
  const float q0 = get_q0();
  const float lambda0 = get_lambda0();
  const float H0 = get_H0();

  const double fz2 = CCALLSFFUN3(FZSQ, fzsq, FLOAT, FLOAT, FLOAT,
				 z, q0, lambda0);

  return sqrt(fz2) * (C_KM_S/H0) * 1000.;
}

/* get the angular diameter distance (kpc) to object at redshift z */
INLINE static double calc_ang_diam_dist_kpc(const float z)
{
  const float lumin_dist = calc_lumin_dist_kpc(z);
  const float z1 = 1. + z;
  return lumin_dist / (z1*z1);
}

/* call the mekal model */
PROTOCCALLSFSUB6(XSMEKL, xsmekl, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
INLINE static void F_mekal(const float* ear, int ne, float T, float Z,
			   float z,
			   int ifl, float* photar, float* photer)
{
  float params[5];
  params[0] = T;
  params[1] = 1;
  params[2] = Z;
  params[3] = z;
  params[4] = 1;

  CCALLSFSUB6(XSMEKL, xsmekl, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV,
	      (float*)ear, ne, params, ifl, photar, photer);
}

static void density_model(const float* ear, int ne, const float* param,
			  int ifl, float* photar, float* photer);
FCALLSCSUB6(density_model, DENSITY, density, FLOATV, INT, FLOATV, INT, FLOATV,
	    FLOATV);

static void density_model(const float* ear, int ne,
			  const float* param,
			  int ifl,
			  float* photar,
			  float* photer)
{
  const float T = param[0];
  const float Z = param[1];
  const float z = param[2];
  const float radius_arc = param[3];
  const float angle = param[4];
  const float density = param[5];

  F_mekal(ear, ne, T, Z, z, ifl, photar, photer);

  {
    const double da_cm = calc_ang_diam_dist_kpc(z)*KPC_CM;
    const double dafact_cm = da_cm * (1.+z);

    const double radius_cm = radius_arc * RADIAN_ARCSEC * da_cm;

    const double vol_cm3 = (4./3.)*M_PI*pow(radius_cm, 3.) *
      (angle/360.);

    const double newnorm = 1e-14 / (dafact_cm*dafact_cm*4.*M_PI) *
      density * density * NH_NE * vol_cm3;

    int i;

    reson_xwrite(15, "da_cm = %e", da_cm);
    reson_xwrite(15, "vol = %e", vol_cm3);

    for( i=0; i != ne; ++i )
      {
	photar[i] *= newnorm;
      }
  }

}
