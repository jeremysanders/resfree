#include <stdio.h>

#include "utils.h"
#include "aped_res_lines.h"
#include "resfree.h"

static const char* const input_filename =
 DATABASEDIR "/aped_resonance_lines.dat";

/* line data are stored here */
static ion_resonance_list ionlist = {0, NULL};

static void read_ions(void)
{
  unsigned num_lines = 0;
  ion_resonance_line_list* current_ion = NULL;

  FILE* f = fopen(input_filename, "r");
  if( f == NULL )
    {
      reson_xwrite(5, "Cannot open resonance line list file '%s'",
		   input_filename);
      return;
    }

  /* loop over lines of file */
  while( ! feof(f) )
    {
      char line[256];
      const char* getsret = fgets( line, sizeof line, f );
      if( getsret == NULL )
	break;

      /* skip commented-out or blank lines */
      if( line[0] == '#' || line[0] == '\n' || line[0] == 0 )
	continue;

      if( line[0] == '@' )
	{
	  /* a new element / ionisation */
	  unsigned element, ionisation;
	  const int ret = sscanf(line, "@ %u %u", &element, &ionisation);
	  if( ret != 2 )
	    {
	      reson_xwrite(5, "Illegal ion line in file '%s':\n%s",
			   input_filename, line);
	      break;
	    }

	  /* allocate memory for new ion */
	  ++ ionlist.no_ions;
	  ionlist.ions = REALLOCSAFE(ionlist.ions, ion_resonance_line_list,
				     ionlist.no_ions);
	  current_ion = & ionlist.ions[ionlist.no_ions-1];

	  current_ion->element = element;
	  current_ion->ion = ionisation;
	  current_ion->nolines = 0;
	  current_ion->energy = NULL;
	  current_ion->oscillator_strength = NULL;

	}
      else
	{
	  /* read in line (energy and  oscillator strength) */
	  double energy, wavelength, oscillator_strength;
	  const int ret = sscanf(line, "%lf %lf %lf", &energy, &wavelength,
				 &oscillator_strength);
	  if( ret != 3 || current_ion == NULL )
	    {
	      reson_xwrite(5, "Illegal line in file %s:\n%s",
			   input_filename, line);
	      break;
	    }

	  /* allocate more memory (yes, I know this is slow but who cares?) */
	  ++ current_ion->nolines;
	  current_ion->energy = REALLOCSAFE(current_ion->energy, double,
					    current_ion->nolines);	    
	  current_ion->oscillator_strength =
	    REALLOCSAFE(current_ion->oscillator_strength,
			double, current_ion->nolines);

	  current_ion->energy[current_ion->nolines-1] = energy;
	  current_ion->oscillator_strength[current_ion->nolines-1] =
	    oscillator_strength;

	  ++num_lines;
	}

    } /* loop over input lines in file */

  fclose(f);

  reson_xwrite(12, "resfree: read in %i resonance lines from database",
	       num_lines);
}


const ion_resonance_list* get_ion_resonance_list(void)
{
  if( ionlist.no_ions == 0 )
    read_ions();
    
  return &ionlist;
}
