#include <math.h>
#include <string.h>

#include "utils.h"
#include "constants.h"
#include "elements.h"

const char * const element_list[NO_ELEMENTS] = {
  "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ",
  "F ", "Ne", "Na", "Mg", "Al", "Si", "P ", "S ",
  "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "V ",
  "Cr", "Mn", "Fe", "Co", "Ni"
};

const double atomic_weights[NO_ELEMENTS] = {
  /* from http://www.webelements.com/ */
  1.00794, /* H */
  4.002602, /* He */
  6.941, /* Li */
  9.012182, /* Be */
  10.811, /* B */
  12.0107, /* C */
  14.0067, /* N */
  15.9994, /* O */
  18.9984032, /* F */
  20.1797, /* Ne */
  22.989770, /* Na */
  24.3050, /* Mg */
  26.981538, /* Al */
  28.0855, /* Si */
  30.973761, /* P */
  32.065, /* S */
  35.453, /* Cl */
  39.948, /* Ar */
  39.0983, /* K */
  40.078, /* Ca */
  44.955910, /* Sc */
  47.867, /* Ti */
  50.9415, /* V */
  51.9961, /* Cr */
  54.938049, /* Mn */
  55.845, /* Fe */
  58.933200, /* Co */
  58.6934 /* Ni */
};

/* ANGR solar abundances
   Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197) */
const double solar_abundances[NO_ELEMENTS] = {
  1.00e+0, /* H */
  9.77e-2, /* He */
  1.45e-11, /* Li */
  1.41e-11, /* Be */
  3.98e-10, /* B */
  3.63e-4, /* C */
  1.12e-4, /* N */
  8.51e-4, /* O */
  3.63e-8, /* F */
  1.23e-4, /* Ne */
  2.14e-6, /* Na */
  3.80e-5, /* Mg */
  2.95e-6, /* Al */
  3.55e-5, /* Si */
  2.82e-7, /* P */
  1.62e-5, /* S */
  1.88e-7, /* Cl */
  3.63e-6, /* Ar */
  1.32e-7, /* K */
  2.29e-6, /* Ca */
  1.26e-9, /* Sc */
  9.77e-8, /* Ti */
  1.00e-8, /* V */
  4.68e-7, /* Cr */
  2.45e-7, /* Mn */
  4.68e-5, /* Fe */
  8.60e-8, /* Co */
  1.78e-6  /* Ni */
};

/* return element number corresponding to element */
/* Hydrogen is 1 */
unsigned identify_element(const char* name)
{
  unsigned i = 0;

  while( i != NO_ELEMENTS && strncmp(name, element_list[i], 2) != 0 )
    ++i;

  if( i != NO_ELEMENTS )
    return i+1;

  reson_xwrite(5, "Unidentified element found in ionic fraction file");

  return 0;
}

/* get fractional line width */
double get_line_width_frac(unsigned atomic_no,
			   double temperature_keV,
			   double turbulence_km_s)
{
  const double thermal = (2. * temperature_keV * KEV_ERG) /
    ( atomic_weights[atomic_no-1] * AMU_G );
  const double turb = turbulence_km_s * KM_CM;

  return sqrt( thermal + turb*turb ) * (1./C_CM_S);
}

