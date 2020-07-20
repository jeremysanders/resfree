#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "ionic_fractions.h"
#include "constants.h"
#include "elements.h"
#include "utils.h"
#include "resfree.h"

#define NO_IONS 29

#define MAZOTTA_FILENAME DATABASEDIR "/mazzotta_table2.dat"
#define TABLE_DELTA_LOGT 0.1
#define TABLE_START_LOGT 4.0

struct TemperatureList
{
  unsigned no_temps;
  double* temps;
};

static int init_tlist = 0;
static struct TemperatureList tlist[NO_ELEMENTS][NO_IONS];

/* read the mazotta list */
/* assumes increments of 0.1 in log T from 4.00 */
static void read_tlist()
{
  FILE *f;
  unsigned element, ion;

  for( element=0; element<NO_ELEMENTS; ++element )
    for( ion=0; ion<NO_IONS; ++ion )
      {
	tlist[element][ion].no_temps = 0;
	tlist[element][ion].temps = NULL;
      }

  f = fopen(MAZOTTA_FILENAME, "r");

  if( f == NULL )
    {
      reson_xwrite(5, "Error: reson: Unable to open ionic fraction file %s",
		   MAZOTTA_FILENAME);
      return;
    }

  while( ! feof(f) )
    {
      char buffer[256];
      struct TemperatureList* tl;
      unsigned element;
      unsigned linelen;
      unsigned ion;

      fgets( buffer, sizeof buffer, f );

      /* skip commented lines */
      if( buffer[0] == '#' || buffer[0] == 0 || buffer[0] == '\n' )
	continue;

      linelen = strlen(buffer);
      element = identify_element( buffer );

      for(ion=0; ion<NO_IONS; ++ion)
	{
	  char num[8];
	  double d;

	  const char* pos = buffer + ion*6 + 7;

	  if( ion*6+7 >= (linelen-1) )
	    {
	      /* insert effectively zero */
	      d = -999;

	    } else {
	      strncpy(num, pos, 6);
	      num[6] = 0; 

	      if( num[0] == ' ' )
		{
		  /* nothing, so assume -999 */
		  d = -999.;
		} else {
		  /* convert */
		  char* retn;
		  d = strtod( pos, &retn );
		  if( retn == num )
		    {
		      reson_xwrite(5, "Error: reson: Unable to convert number '%s'", num);
		      d = -999.;
		    }
		}

	    }

	  /* add new ionisation */
	  tl = & (tlist[element-1][ion]);
	  (tl -> no_temps)++;
	  tl->temps = REALLOCSAFE( tl->temps, double, tl -> no_temps );
	  tl->temps[tl->no_temps - 1] = d;
	}
    }

  fclose(f);
}

double get_ionisation_fraction(double T_keV,
			       unsigned element, unsigned ionisation)
{
  unsigned no, nextno;
  struct TemperatureList* tl;
  const double logT_K = log10( T_keV * KEV_K );
  double deltafrac;
  double frac;

  /* read in data if we haven't */
  if( ! init_tlist )
    {
      read_tlist();
      init_tlist = 1;
    }

  /* start of the table */
  if( logT_K < 4. )
    return 0.;

  /* work out where in the table we need to be */
  no = (unsigned)( ( logT_K - TABLE_START_LOGT ) / TABLE_DELTA_LOGT );

  tl = & (tlist[element-1][ionisation-1]);
  if( tl->no_temps == 0 )
    {
      reson_xwrite(5, "Error: reson: No data in ionisation state %i"
		   " of element %i", ionisation, element);
      return 0.;
    }

  /* final result is okay if we go past the end of the table */
  if( no >= tl->no_temps )
    no = tl->no_temps - 1;

  nextno = no+1;
  if( nextno >= tl->no_temps )
    nextno = tl->no_temps - 1;

  deltafrac = (logT_K - TABLE_START_LOGT - no*TABLE_DELTA_LOGT) /
    TABLE_DELTA_LOGT;

  /* weight temperatures */
  frac = tl->temps[no] * (1. - deltafrac) + tl->temps[nextno] * deltafrac;

  return pow(10., frac);
}

/* void reson_xwrite(int level, const char* format, ...) */
/* { */
/*   printf(format); */

/* } */


/* int main() */
/* { */
/*   int i; */
/*   for(i=1; i<80; ++i) */
/*     { */
/*       const double T_keV = i*0.1; */
/*       const double frac = get_ionisation_fraction( T_keV, 12, 11 ); */
/*       printf("%e %e\n", T_keV, frac); */
/*     } */


/*   return 0; */
/* } */
