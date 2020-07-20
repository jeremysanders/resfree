#include <assert.h>
#include <fitsio.h>

#include "linefrac.h"
#include "utils.h"
#include "constants.h"
#include "elements.h"
#include "resfree.h"

/* store the continuua for a particular element at a certain temperature */
typedef struct
{
  unsigned atomic_no;
  size_t no_energies;
  float* energies;
  float* continuua;
} continuum_record;

/* store the continuua at a particular temperature */
typedef struct
{
  size_t no_elements;
  continuum_record* continuum;
} cont_temp_record;

/* store set of continuua at different temperatures */
typedef struct
{
  size_t no_temp;
  float* temp;
  cont_temp_record* temp_rec;
} cont_dataset;

/* store lines at a particular temperature */
typedef struct
{
  size_t no_lines;

  unsigned* element;
  float* energy;
  float* epsilon;
} line_temp_record;

/* store sets of lines as a function of temperature */
typedef struct
{
  size_t no_temp;
  float* temp;
  line_temp_record* line_rec;
} line_dataset;

/* whether datasets loaded */
static int loaded_datasets = 0;

/* store the continuua and pseudo continuua */
static cont_dataset continuua;
static cont_dataset pseudo_continuua;

/* store the lines */
static line_dataset lines;

#define NO_APEC_ELEMENTS 14
/* atomic numbers of elements in the apec table */
const static unsigned apec_atomic_no_list[NO_APEC_ELEMENTS] =
  { 1, 2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28 };
static int atomic_reverse[NO_ELEMENTS];

/* this initialses a table so we can lookup elements by the atomic
   number in a table and get their position in the apec abundance table - 
   see apec_atomic_no_list */

static void init_reverse_atom_lookup(void)
{
  int i;
  for(i = 0; i != NO_ELEMENTS; ++i)
    atomic_reverse[i] = -1;

  for(i = 0; i != NO_APEC_ELEMENTS; ++i)
    atomic_reverse[ apec_atomic_no_list[i] - 1 ] = i;
}

/* add a new temperature record to a continuum dataset */
static cont_temp_record* new_cont_temp_rec(cont_dataset* ds,
					   float val, unsigned no_elements)
{
  /* make sure the file format doesn't change */
  assert( no_elements == NO_APEC_ELEMENTS );

  ds->no_temp++;
  ds->temp = REALLOCSAFE(ds->temp, float, ds->no_temp);
  assert( ds->temp != NULL );

  ds->temp_rec = REALLOCSAFE(ds->temp_rec, cont_temp_record, ds->no_temp);
  assert( ds->temp_rec != NULL );

  ds->temp[ds->no_temp-1] = val;
  ds->temp_rec[ds->no_temp-1].no_elements = no_elements;
  ds->temp_rec[ds->no_temp-1].continuum = ALLOCSAFE( continuum_record,
						     no_elements);
  assert( ds->temp_rec[ds->no_temp-1].continuum != NULL );

  return &ds->temp_rec[ds->no_temp-1];
}

static void read_continuum(const char* const filename,
			   const char* const no_continuum_name,
			   const char* const energy_name,
			   const char* const continuum_name,
			   cont_dataset* dataset)
{
  int status = 0;
  fitsfile* ff;

  /* open file */
  reson_xwrite(16, "fastapec::read_continuum: Opening '%s'", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  if( status )
    {
      reson_xwrite(5, "Could not open '%s'", filename);
      return;
    }

  /* loop over the HDUs */
  while( 1 )
    {
      int hdutype;
      char extname[64];

      reson_xwrite(16, "fastapec::read_continuum: Moving to next HDU");
      /* move to next HDU */
      fits_movrel_hdu(ff, 1, &hdutype, &status);
      /* there are no more remaining */
      if( status )
	break;

      fits_read_key(ff, TSTRING, "EXTNAME", extname, 0, &status);
      if( status )
	{
	  reson_xwrite(5, "Could not read EXTNAME header");
	  return;
	}

      /* we have an emissivity extension */
      if( strcmp(extname, "EMISSIVITY") == 0 )
	{
	  int no_col, cont_col, energy_col, atom_col;
	  int no_elements;
	  int element;
	  float temperature_K;
	  cont_temp_record* trec;

	  reson_xwrite(16, "fastapec::read_continuum: Found EMISSIVITY HDU");

	  /* get no of elements */
	  fits_read_key(ff, TINT, "NAXIS2", &no_elements, 0, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could not read"
			   " NAXIS2 keyword");
	      return;
	    }

	  /* get temperature */
	  fits_read_key(ff, TFLOAT, "TEMPERATURE", &temperature_K,
			0, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could not read"
			   " TEMPERATURE keyword");
	      return;
	    }
	  /* make new temperature record */
	  trec = new_cont_temp_rec(dataset,
				   temperature_K / KEV_K, no_elements);

	  /* lookup columns for reading */
	  fits_get_colnum(ff, CASEINSEN, (char*)no_continuum_name,
			  &no_col, &status);
	  fits_get_colnum(ff, CASEINSEN, (char*)energy_name,
			  &energy_col, &status);
	  fits_get_colnum(ff, CASEINSEN, (char*)continuum_name,
			  &cont_col, &status);
	  fits_get_colnum(ff, CASEINSEN, "Z", &atom_col, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could "
			   "not read %s, %s, %s or Z table headers",
			   no_continuum_name, energy_name, continuum_name);
	      return;
	    }

	  /* loop over rows */
	  for(element=0; element != no_elements; ++element)
	    {
	      continuum_record* cr = &trec->continuum[element];
	      int i;

	      /* get atomic no */
	      fits_read_col(ff, TINT, atom_col, element+1, 1, 1, 0,
			    &cr->atomic_no, 0, &status);

	      if( status )
		{
		  reson_xwrite(5, "fastapec::read_continuum: Could not read "
			       "Z column");
		  return;
		}

	      /* read no of continuua elements */
	      fits_read_col(ff, TINT, no_col, element+1, 1, 1, 0,
			    &i, 0, &status);
	      if( status )
		{
		  reson_xwrite(5, "fastapec::read_continuum: Could not read "
			       "%s column", no_continuum_name);
		  return;
		}
	      cr->no_energies = i;

	      /* do memory allocation */
	      cr->energies = ALLOCSAFE(float, cr->no_energies);
	      assert( cr->energies != NULL );
	      cr->continuua = ALLOCSAFE(float, cr->no_energies);
	      assert( cr->continuua != NULL );

	      /* read data */
	      fits_read_col(ff, TFLOAT, energy_col, element+1,
			    1, cr->no_energies, 0,
			    cr->energies, 0, &status);
	      if( status )
		{
		  reson_xwrite(5, "fastapec::read_continuum: Could not read "
			       "%s data", energy_name);
		  return;
		}
	      fits_read_col(ff, TFLOAT, cont_col, element+1,
			    1, cr->no_energies, 0,
			    cr->continuua, 0, &status);
	      if( status )
		{
		  reson_xwrite(5, "fastapec::read_continuum: Could not read "
			       "%s data", continuum_name);
		  return;
		}

/* 	      reson_xwrite(5, "Test: T=%e, el=%i, no=%i, 1=%g, 2=%g", */
/* 			   temperature_K, element, cr->no_energies, */
/* 			   cr->continuua[0], */
/* 			   cr->continuua[1]); */
	    }

	}
    }
  status = 0;
  fits_clear_errmsg();

  /* we're done */
  fits_close_file(ff, &status);
}

/* add a new line record to the line dataset */
static line_temp_record* new_line_temp_rec(float temp, unsigned no_lines)
{
  lines.no_temp++;
  lines.temp = REALLOCSAFE(lines.temp, float, lines.no_temp);
  assert( lines.temp != NULL );

  lines.line_rec = REALLOCSAFE(lines.line_rec, line_temp_record,
			       lines.no_temp);
  assert( lines.line_rec != NULL );

  lines.temp[lines.no_temp-1] = temp;
  lines.line_rec[lines.no_temp-1].no_lines = no_lines;
  lines.line_rec[lines.no_temp-1].element = ALLOCSAFE(unsigned, no_lines);
  assert( lines.line_rec[lines.no_temp-1].element != NULL );
  lines.line_rec[lines.no_temp-1].energy = ALLOCSAFE(float, no_lines);
  assert( lines.line_rec[lines.no_temp-1].energy != NULL );
  lines.line_rec[lines.no_temp-1].epsilon = ALLOCSAFE(float, no_lines);
  assert( lines.line_rec[lines.no_temp-1].epsilon != NULL );

  return &lines.line_rec[lines.no_temp-1];
}

static void read_lines(const char* const filename)
{
  int status = 0;
  fitsfile* ff;

  /* open file */
  reson_xwrite(16, "fastapec::read_continuum: Opening '%s'", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  if( status )
    {
      reson_xwrite(5, "Could not open '%s'", filename);
      return;
    }

  /* loop over the HDUs */
  while( 1 )
    {
      int hdutype;
      char extname[64];

      reson_xwrite(16, "fastapec::read_continuum: Moving to next HDU");
      /* move to next HDU */
      fits_movrel_hdu(ff, 1, &hdutype, &status);
      /* there are no more remaining */
      if( status )
	break;

      fits_read_key(ff, TSTRING, "EXTNAME", extname, 0, &status);
      if( status )
	{
	  reson_xwrite(5, "Could not read EXTNAME header");
	  return;
	}

      /* we have an emissivity extension */
      if( strcmp(extname, "EMISSIVITY") == 0 )
	{
	  int no_lines;	  
	  float temperature_K;
	  line_temp_record* trec;
	  int lambda_col, epsilon_col, element_col;
	  int i;

	  /* get number of lines */
	  fits_read_key(ff, TINT, "NAXIS2", &no_lines, 0, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_lines: Could not read"
			   " NAXIS2 keyword");
	      return;
	    }
	  
	  /* get temperature */
	  fits_read_key(ff, TFLOAT, "TEMPERATURE", &temperature_K,
			0, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could not read"
			   " TEMPERATURE keyword");
	      return;
	    }

/* 	  trec = new_line_temp_rec(temperature_K / KEV_K, 1); */
/* 	  trec->energy[0] = 2; */
/* 	  trec->epsilon[0] = 1e-15; */
/* 	  trec->element[0] = 26; */
/* 	  continue; */

	  /* make new temperature record */
	  trec = new_line_temp_rec(temperature_K / KEV_K, no_lines);

	  /* lookup columns for reading */
	  fits_get_colnum(ff, CASEINSEN, "Lambda", &lambda_col, &status);
	  fits_get_colnum(ff, CASEINSEN, "Epsilon", &epsilon_col, &status);
	  fits_get_colnum(ff, CASEINSEN, "Element", &element_col, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could "
			   "not read Lambda, Epsilon or "
			   "Element table headers");
	      return;
	    }

	  /* read lines */
	  fits_read_col(ff, TFLOAT, lambda_col, 1, 1, no_lines,
			0, trec->energy, 0, &status);
	  fits_read_col(ff, TFLOAT, epsilon_col, 1, 1, no_lines,
			0, trec->epsilon, 0, &status);
	  fits_read_col(ff, TINT, element_col, 1, 1, no_lines,
			0, trec->element, 0, &status);
	  if( status )
	    {
	      reson_xwrite(5, "fastapec::read_continuum: Could "
			   "not read Lambda, Epsilon or "
			   "Element table columns");
	      return;
	    }

/* 	  for(i=0; i != no_lines; ++i) */
/* 	    { */
/* 	      reson_xwrite(5, "%e %e %i", trec->energy[i], trec->epsilon[i], */
/* 			   trec->element[i]); */
/* 	    } */


	  /* convert lambdas (in Angs) to keV */
	  for(i=0; i != no_lines; ++i)
	    trec->energy[i] = KEV_ANGSTROM / trec->energy[i];
	}

    }

  /* we're done */
  fits_clear_errmsg();
  fits_close_file(ff, &status);
}

static void load_apec_datasets(void)
{
  continuua.no_temp = 0;
  continuua.temp = NULL;
  continuua.temp_rec = NULL;
  pseudo_continuua.no_temp = 0;
  pseudo_continuua.temp = NULL;
  pseudo_continuua.temp_rec = NULL;
  lines.no_temp = 0;
  lines.temp = NULL;
  lines.line_rec = NULL;

  read_continuum(XSPECDATADIR "/apec_v" APECVERSION "_coco.fits",
		 "N_Cont", "E_Cont", "Continuum", &continuua);
  read_continuum(XSPECDATADIR "/apec_v" APECVERSION "_coco.fits",
		 "N_Pseudo", "E_Pseudo", "Pseudo", &pseudo_continuua);
  read_lines(XSPECDATADIR "/apec_v" APECVERSION "_line.fits");

}

/***********************************************************************/

static void add_continuua(const continuum_record* const cr,
			  const double coeff,
			  const double z,
			  float* photar)
{
  const double z1 = 1 + z;
  const float* energies = cr->energies;
  const float* cont = cr->continuua;
  const int noenergies = cr->no_energies;
  const double minenergy = energies[0];
  const double maxenergy = energies[noenergies - 1];

  int input_index = 0;
  int i;

  for(i=0; i != ENERGY_NOBINS; ++i)
    {
      const double energy = ( ENERGY_MIN + ENERGY_SPACING*0.5 +
			      ENERGY_SPACING*i ) * z1;

      /* if energy is within model range */
      if( energy > minenergy && energy < maxenergy )
	{
	  /* loop until we're at the next available point to interpolate
	     from */
	  while( energies[input_index] < energy && input_index < noenergies )
	    input_index++;

	  {
	    const double f1 = (energies[input_index] - energy) *
	      cont[input_index-1];
	    const double f2 = (energy - energies[input_index-1]) *
	      cont[input_index];
	    const double dear = ENERGY_SPACING;
	    const double denergy = energies[input_index] -
	      energies[input_index-1];

	    const double delta = coeff * (f1 + f2) * dear * z1 / denergy;

	    photar[i] += delta;
	  }
	}
    }
}

static void add_lines( const line_temp_record* const lr,
		       const double T_keV,
		       const float* const abundances,
		       const double z,
		       const double turb_kmps,
		       const double mult,
		       float* photar )
{
  /* for easy access */
  const unsigned* const elements = lr->element;
  const float* const energies = lr->energy;
  const float* const epsilons = lr->epsilon;
  const size_t nolines = lr->no_lines;
  const double z_1_1 = 1. / (1. + z);
  size_t line;

  /* loop over lines */
  for(line=0; line != nolines; ++line)
    {
      /* restframe energy of the line */
      const double restframe_energy_keV = energies[line];
      /* redshifted energy of the line */
      const double line_energy_keV =  restframe_energy_keV * z_1_1;
      /* which bin the centre of line is in */      
      const int centre_bin = (int)( floor( (line_energy_keV - ENERGY_MIN) *
					   (1./ENERGY_SPACING) ) );

      const unsigned element = elements[line];
      const int apec_element = atomic_reverse[ element - 1 ];

      const double width_keV =  restframe_energy_keV * (1./C_CM_S) *
	sqrt( 2 * T_keV * KEV_ERG / ( atomic_weights[element-1]*AMU_G ) +
	      turb_kmps*turb_kmps*KM_CM*KM_CM );

      const double flux_2 = mult * epsilons[line] * abundances[apec_element] *
	0.5;

      int bin;
      double lastfrac;

      /* upper energy part of line */
      bin = centre_bin;
      lastfrac = 0;

      while( bin < ENERGY_NOBINS )
	{
	  const double energy = ENERGY_MIN + (bin+1)*ENERGY_SPACING;

	  /* we get -ve val if we go past the 'edge' of the line */
	  const double frac = line_frac( line_energy_keV, width_keV, energy );
	  const int lastbin = frac < 0;
	  const double thisfrac = lastbin ? 1 : frac;

	  /* add on flux to bin */
	  if( bin >= 0 )
	    photar[bin] += flux_2 * (thisfrac - lastfrac);

	  if( lastbin )
	    break;

	  lastfrac = frac;
	  ++bin;
	}

      /* lower energy part of line */
      bin = centre_bin;
      lastfrac = 0;

      while( bin >= 0 )
	{
	  const double energy = ENERGY_MIN + bin*ENERGY_SPACING;

	  /* we get -ve val if we go past the 'edge' of the line */
	  const double frac = line_frac( line_energy_keV, width_keV, energy );
	  const int lastbin = frac < 0;
	  const double thisfrac = lastbin ? 1 : frac;

	  /* add on flux to bin */
	  if( bin < ENERGY_NOBINS )
	    photar[bin] += flux_2 * (thisfrac - lastfrac);

	  if( lastbin )
	    break;

	  lastfrac = frac;
	  --bin;
	}

    } /* loop over lines */
}

/* search for temperature using binary search, should return lower bound */
static int binary_search_temperature(float T_keV)
{
  int high, low;

  T_keV = MAX( lines.temp[0], T_keV );
  T_keV = MIN( lines.temp[lines.no_temp-1], T_keV );

  high = lines.no_temp;
  low = 0;
  while( high-low > 1 )
    {
      const int index = (low+high)/2;

      if( T_keV < lines.temp[index] )
	high = index;
      else
	low = index;
    }

  return low;
}

/************************************************************************/
/* routines to make/use cached spectra */

static float* cached_continuum_H_He = NULL;
static float* cached_continuum_metals = NULL;

static void recalc_cached_continuua_at_abun(const float* const abundances,
					    const float redshift,
					    float* const outspectra)
{
   const unsigned notemps = continuua.no_temp;
   unsigned temp = 0;

   /* iterate over temperature components */
   for(temp=0; temp != notemps; ++temp)
     {
       float* const outspectrum = outspectra + ENERGY_NOBINS*temp;
       unsigned i;
       unsigned el;

       /* zero output */
       for(i=0; i != ENERGY_NOBINS; ++i)
	 outspectrum[i] = 0;

       /* iterate over elements */
       for(el=0; el != NO_APEC_ELEMENTS; ++el)
	 {
	   /* add continuum */
	   add_continuua( &continuua.temp_rec[temp].continuum[el],
			  abundances[el], redshift, outspectrum );

	   /* add pseudo-continuum */
	   add_continuua( &pseudo_continuua.temp_rec[temp].continuum[el],
			  abundances[el], redshift, outspectrum );
	 }
     }

}

static void recalc_cached_continuua(const float redshift)
{
  const unsigned notemps = continuua.no_temp;
  float abundances[NO_APEC_ELEMENTS];
  unsigned i;

  reson_xwrite(12, "fastapec:: recalculating cached spectra (z=%e)", redshift);

  /* allocate memory if necessary */
  if( cached_continuum_H_He == NULL )
    cached_continuum_H_He = ALLOCSAFE( float, notemps*ENERGY_NOBINS );

  if( cached_continuum_metals == NULL )
    cached_continuum_metals = ALLOCSAFE( float, notemps*ENERGY_NOBINS );

  /* calculate He and H continuua */
  for( i = 0; i != NO_APEC_ELEMENTS; i++ )
    abundances[i] = 0;

  abundances[0] = abundances[1] = 1;
  recalc_cached_continuua_at_abun( abundances, redshift,
				   cached_continuum_H_He );

  /* calculate metals continuua */
  for( i = 0; i != NO_APEC_ELEMENTS; i++ )
    abundances[i] = 1;

  abundances[0] = abundances[1] = 0;
  recalc_cached_continuua_at_abun( abundances, redshift,
				   cached_continuum_metals );

}

void fast_apec(const float T_keV,
	       const float z,
	       const float turb_kmps,
	       const float abundance,
	       float* photar)
{
  static float cached_redshift = -1;
  float abundances[NO_APEC_ELEMENTS];
  double fractions[2];
  int Tindex;
  int fraction_index;
  int i;

  /* load in data if not already loaded */
  if( ! loaded_datasets )
    {
      load_apec_datasets();
      init_reverse_atom_lookup();
      loaded_datasets = 1;
    }

  /* recalculate cached spectra if z has changed (or first time around) */
  if( cached_redshift != z )
    {
      recalc_cached_continuua( z );
      cached_redshift = z;
    }

  /* find nearest temperature using a binary search */
  /* Tindex = index from start of temperature array */
  Tindex = binary_search_temperature(T_keV);

  /* work out weighting factors */
  if( T_keV == lines.temp[Tindex] )
    {
      fractions[0] = 1;
      fractions[1] = 0;
    }
  else
    {
      const double d1 = T_keV - lines.temp[Tindex];
      const double d2 = lines.temp[Tindex+1] - T_keV;
      const double dt = d1+d2;

      fractions[0] = d2 / dt;
      fractions[1] = 1 - fractions[0];
    }
  reson_xwrite(15, "fastapec: Interpolating from %e to %e keV (fractions "
	       "%f and %f)", lines.temp[Tindex], lines.temp[Tindex+1],
	       fractions[0], fractions[1]);
  
  /* zero output */
  for(i=0; i != ENERGY_NOBINS; ++i)
    photar[i] = 0;

  /* set abundances */
  abundances[0] = abundances[1] = 1;
  for(i=2; i != NO_APEC_ELEMENTS; ++i)
    abundances[i] = abundance;

  /* iterate over two temperature components */
  /* loop between +0 and +1 on Tindex */
  for(fraction_index = 0; fraction_index <= 1; ++fraction_index)
    {
      const int index = Tindex + fraction_index;
      const double fraction = fractions[fraction_index];
      const double abun_fraction = abundance*fraction;

      if( index < lines.no_temp )
	{
	  { /* add hydrogen and helium */
	    const float* const continuum_He_H = cached_continuum_H_He +
	      index*ENERGY_NOBINS;

	    for(i=0; i != ENERGY_NOBINS; ++i)
	      photar[i] += continuum_He_H[i] * fraction;
	  }
	  { /* add metal continuum */
	    const float* const continuum_metals = cached_continuum_metals +
	      index*ENERGY_NOBINS;

	    for(i=0; i != ENERGY_NOBINS; ++i)
	      photar[i] += continuum_metals[i] * abun_fraction;
	  }
	  { /* add lines */
	    const line_temp_record* lr = &lines.line_rec[index];

	    add_lines( lr, lines.temp[index], abundances,
		       z, turb_kmps, fraction, photar );
	  }
	} /* if temperature within range */

    } /* interpolated temperatures */

  {
    /* correct for norm and time dilation */
    const double zf = 1./(1.+z) * 1e14;
    for(i=0; i != ENERGY_NOBINS; ++i)
      photar[i] *= zf;
  }
}


void make_fast_apec_model(const float T_keV,
			  const float z,
			  const float turb_kmps,
			  const float* const abundances,
			  float* photar)
{
  int i;
  int Tindex;
  double f[2];
  int findex;

  /* load in data if not already loaded */
  if( ! loaded_datasets )
    {
      load_apec_datasets();
      init_reverse_atom_lookup();
      loaded_datasets = 1;
    }

  /* zero output */
  for(i=0; i != ENERGY_NOBINS; ++i)
    photar[i] = 0;

  /* find nearest temperature using a binary search */
  Tindex = binary_search_temperature(T_keV);

  reson_xwrite(15, "fastapec: Tindex=%i", Tindex);

  /* work out weighting factors */
  if( T_keV == lines.temp[Tindex] )
    {
      f[0] = 1;
      f[1] = 0;
    } else {
      const double d1 = T_keV - lines.temp[Tindex];
      const double d2 = lines.temp[Tindex+1] - T_keV;
      const double dt = d1+d2;

      f[0] = d2 / dt;
      f[1] = 1 - f[0];

    }
  reson_xwrite(15, "fastapec: Interpolating from %e to %e keV (fractions "
	       "%f and %f)", lines.temp[Tindex], lines.temp[Tindex+1],
	       f[0], f[1]);

  /* loop between +0 and +1 on Tindex */
  for(findex = 0; findex != 2; ++findex)
    {
      if( findex + Tindex < lines.no_temp )
	{
	  const cont_temp_record* tr_cont = &continuua.temp_rec[Tindex+findex];
	  const cont_temp_record* tr_pseud = &pseudo_continuua.
	    temp_rec[Tindex+findex];
	  const line_temp_record* lr = &lines.line_rec[Tindex+findex];
	  unsigned el;

	  for(el=0; el != NO_APEC_ELEMENTS; ++el)
	    {
	      const double factor = f[findex]*abundances[el];

	      /* add continuua */
	      add_continuua( &tr_cont->continuum[el], factor, z, photar );
	      add_continuua( &tr_pseud->continuum[el], factor, z, photar );
	    }

	  /* add appropriate lines */
	  add_lines( lr, lines.temp[Tindex+findex], abundances,
		     z, turb_kmps, f[findex], photar );

	} /* if temperature within range */

    } /* interpolated temperatures */

  {
    /* correct for norm and time dilation */
    const double zf = 1./(1.+z) * 1e14;
    for(i=0; i != ENERGY_NOBINS; ++i)
      photar[i] *= zf;
  }
}
