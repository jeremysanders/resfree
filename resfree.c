#include <string.h>
#include <assert.h>

#include <cfortran.h>

#include "resfree.h"
#include "constants.h"
#include "utils.h"
#include "aped_res_lines.h"
#include "elements.h"
#include "ionic_fractions.h"
#include "fastapec.h"

static double angular_diameter_dist_cm = -1.;
static double angular_diameter_dist_kpc = -1.;

static unsigned total_subshells = (unsigned)(-1);
static unsigned total_shells = (unsigned)(-1);
static float* shell_spectra = NULL;
static float* shell_ear = NULL;
static float* shell_absorb_spectra = NULL;
static float* shell_scattered_spectra = NULL;

/* observed spectra for each subshell built up in here */
static float* shell_output_spectra = NULL;

/* ln absorption for unit NH here at redshift of object */
static float* log_zphabs_absorption = NULL;

/* these arrays are populated with the radii of the shells.
   The inner radius is item 0 */
static float* shell_radii_arcsec = NULL;
static double* shell_radii_cm = NULL;
static float* subshell_radii_arcsec = NULL;
static double* subshell_radii_cm = NULL;

/* physical parameters of the subshells */
static float* subshell_kT = NULL;
static float* subshell_Z = NULL;
static float* subshell_ne = NULL;
static float* subshell_nH = NULL;

/* get number of unique outer radii in loaded spectra */
static unsigned count_unique_shells(void)
{
  const unsigned nodatasets = get_no_datasets();
  double inner = -1;
  unsigned count = 0;
  unsigned i;

  for(i = 0; i < nodatasets; ++i)
    {
      const double r = reson_get_dataset_filter(i+1, 1);
      if( ! near_double(r, inner) )
	{
	  count++;
	  inner = r;
	}
    }

  return count;
}

/* reallocates the memory for all the temporary arrays */
static void realloc_shells(const float * const param)
{
  /*const unsigned nodatasets = get_no_datasets();*/
  const unsigned new_total_shells = count_unique_shells();
  const unsigned nosubshells = (unsigned)(param[PARAMS_SUBSHELLS]);
  const unsigned new_total_subshells = nosubshells * new_total_shells;

  if( new_total_shells > PARAMS_NOSHELLS )
    {
      reson_xwrite(5, "resfree: Too many shells. Maximum is %i",
		   PARAMS_NOSHELLS);
    }

  /* reallocate the memory */
  if( total_subshells != new_total_subshells || shell_spectra == NULL )
    {
      reson_xwrite(15, "resfree: reallocating memory");

      /* allocate memory for shells */
      float_realloc( &shell_spectra, new_total_subshells*ENERGY_NOBINS );
      /* allocate memory for output spectra */
      float_realloc( &shell_output_spectra, new_total_subshells*ENERGY_NOBINS );
      /* store absorption coefficients likewise */
      float_realloc( &shell_absorb_spectra, new_total_subshells*ENERGY_NOBINS );
      /* store scattered radiation */
      float_realloc( &shell_scattered_spectra, new_total_subshells*ENERGY_NOBINS );
      /* allocate space for shell and subshell radii */
      float_realloc( &shell_radii_arcsec, new_total_shells+1 );
      float_realloc( &subshell_radii_arcsec, new_total_subshells+1 );
      double_realloc( &shell_radii_cm, new_total_shells+1 );
      double_realloc( &subshell_radii_cm, new_total_subshells+1 );
      
      /* allocate space for precomputing the temperature, abundance and 
	 density */
      float_realloc( &subshell_kT, new_total_subshells );
      float_realloc( &subshell_Z, new_total_subshells );
      float_realloc( &subshell_ne, new_total_subshells );
      float_realloc( &subshell_nH, new_total_subshells );
      total_subshells = new_total_subshells;
      total_shells = new_total_shells;
    }

  /* allocate and fill shell ear array if haven't */
  if( shell_ear == NULL )
    {
      unsigned i;
      shell_ear = ALLOCSAFE( float, ENERGY_NOBINS+1 );
      assert( shell_ear != NULL );

      for( i = 0; i <= ENERGY_NOBINS; ++i ) 
	shell_ear[i] = ENERGY_MIN + ENERGY_SPACING*i;
    }

  /* allocate log absorption spectrum */
  if( log_zphabs_absorption == NULL )
    {
      log_zphabs_absorption = ALLOCSAFE( float, ABSORB_BINS );
    }
}

/* compute the radii of the subshells from the shells */
static void compute_subshell_radii(const float * const param)
{
  const unsigned nosubshells = (unsigned)(param[PARAMS_SUBSHELLS]);
  const unsigned nodatasets = get_no_datasets();
  unsigned shell, dataset, subshell;
  double lastradius = -1;

  /* now fill up radius array */
  shell_radii_arcsec[0] = param[PARAMS_INNER_ARCSEC];
  shell_radii_cm[0] = shell_radii_arcsec[0] *
    RADIAN_ARCSEC * angular_diameter_dist_cm;

  shell = 1;
  /* get radii from shells */
  for(dataset = 1; dataset <= nodatasets; ++dataset)
    {
      const float r = reson_get_dataset_filter(dataset, 1);
      /* skip datasets where radius doesn't change */
      if( ! near_double(r, lastradius) )
	{
	  shell_radii_arcsec[shell] = r;
	  shell_radii_cm[shell] = r * RADIAN_ARCSEC * angular_diameter_dist_cm;
	  shell++;
	  lastradius = r;
	}
    }

  /* work out subshell radii */
  for(shell=0; shell != total_shells; ++shell)
    {
      const double delta_arcsec = (shell_radii_arcsec[shell+1] -
				   shell_radii_arcsec[shell])/nosubshells;
      
      for(subshell=0; subshell != nosubshells; ++subshell)
	{
	  const unsigned index = subshell+shell*nosubshells;
	  
	  subshell_radii_arcsec[index] = shell_radii_arcsec[shell] +
	    delta_arcsec*subshell;
	  subshell_radii_cm[index] = subshell_radii_arcsec[index] *
	    RADIAN_ARCSEC * angular_diameter_dist_cm;
	}
    }
  
  subshell_radii_arcsec[total_subshells] = shell_radii_arcsec[total_shells];
  subshell_radii_cm[total_subshells] = shell_radii_cm[total_shells];
}

/* fill up absorption spectrum for unit absorption */
static void precompute_zphabs_absorption(const float * const param)
{
  float absorb_ear[ABSORB_BINS+1];
  float photer[ABSORB_BINS];
  float zphabs_params[2];
  unsigned i;

  zphabs_params[0] = 1.;
  zphabs_params[1] = param[PARAMS_REDSHIFT];

  /* absorption is calculated on a downsampled spectrum */
  /* this speeds up the routines which apply absorption significantly */
  for(i=0; i <= ABSORB_BINS; ++i)
    {
      absorb_ear[i] = ENERGY_MIN + i*ABSORB_DOWNSAMPLE*ENERGY_SPACING;
    }

  F_zphabs(absorb_ear, ABSORB_BINS,
	   zphabs_params, 1, log_zphabs_absorption, photer);
  for(i=0; i < ABSORB_BINS; ++i)
    {
      log_zphabs_absorption[i] = log(log_zphabs_absorption[i]);
    }
}

/* apply absorption of NH (10^22) to spectrum given */
static INLINE void apply_absorption(const float NH, float* const spectrum)
{
  if ( NH > 1e-6 )
    {
      /* absorption is downsampled */
      unsigned i, j;
      for( i=0; i < ABSORB_BINS; ++i )
	{
	  const float absorb = exp(NH * log_zphabs_absorption[i]);
	  const unsigned base = i * ABSORB_DOWNSAMPLE;
	  if(absorb < 0.99999)
	    for( j=0; j < ABSORB_DOWNSAMPLE; ++j )
	      spectrum[base+j] *= absorb;
	}
    }
}

/* fill up the structure shell_spectra with the spectra for the shells */
static void compute_emission_spectra(const float * const param)
{
  unsigned subshell;

  for(subshell=0; subshell != total_subshells; ++subshell)
    {
      float * const spectrum = shell_spectra + subshell*ENERGY_NOBINS;

      fast_apec( subshell_kT[subshell],
		 param[PARAMS_REDSHIFT],
		 param[PARAMS_TURB],
		 subshell_Z[subshell],
		 spectrum );
    }
}

/* compute the scattering per unit length as a fn of energy in a shell */
static void compute_scatter_in_shell(const float * const param,
				     const unsigned shell,
				     float * const absorb_out)
{
  const double redshift_factor = 1./( 1. + param[PARAMS_REDSHIFT] );

  /* get density, temperature and abundance */
  const double nH_cm3 = subshell_ne[shell]*NH_NE;
  const double kT_keV = subshell_kT[shell];
  const double Z_solar = subshell_Z[shell];
  const double turb_kmps = param[PARAMS_TURB];
 
  const ion_resonance_list* rl = get_ion_resonance_list();
  const unsigned no_ions = rl->no_ions;

  unsigned ion;
  /* iterate over each ion */
  for( ion=0; ion != no_ions; ++ion )
    {
      /* get list of lines */
      const ion_resonance_line_list* const list = &( rl->ions[ion] );
      const unsigned element = list->element;
      const unsigned ion = list->ion;
      const unsigned nolines = list->nolines;
      const double* const energies = list->energy;
      const double* const osc_strength = list->oscillator_strength;

      /* fractional width of lines of this element */
      const double frac_width = get_line_width_frac(element, kT_keV,
						    turb_kmps);
      /* element no density */
      const double element_cm3 = nH_cm3 * solar_abundances[element - 1] *
	Z_solar;

      /* ion no density */
      const double ion_fraction = get_ionisation_fraction(kT_keV, element,
							  ion);
      const double ion_cm3 = element_cm3 * ion_fraction;

      /* iterate over lines */
      unsigned i;
      for(i=0; i != nolines; ++i)
	{
	  const double energy_keV = energies[i] * redshift_factor;
	  const double width_keV = energies[i] * frac_width;
	  const double width_Hz = width_keV * KEV_HZ;

	  /* convert density to absorption coefficient */
	  const double absorp_const = M_PI * E_ESU * E_ESU /
	    ME_G / C_CM_S / PI_SQRT;
	  const double coeff = ion_cm3 * absorp_const * osc_strength[i] /
	    width_Hz;

	  if( coeff > 0. )
	    {
	      /* put the line in the spectrum */
	      insert_absorb_line( absorb_out, coeff, energy_keV, width_keV );
	    }

	} /* loop over lines in ion */

    } /* loop over ions */
}

/* compute scattering in all of the shells */
static void compute_scatter_spectra(const float * const param)
{
  const int scatsw = param[PARAMS_SCATSW] > 1e-10;
  unsigned shell;

  /* exit if no scattering */
  if( ! scatsw )
    return;

  /* iterate over all the subshells */
  for(shell=0; shell != total_subshells; ++shell)
    {
      float * const absorb_spectrum = shell_absorb_spectra +
	shell*ENERGY_NOBINS;
      
      /* blank output spectrum */
      blank_spectrum( absorb_spectrum );

      /* work out scattering if wanted */
      compute_scatter_in_shell(param, shell, absorb_spectrum);
    }
}

/* interpolate fit parameters to get physical parameters in the subshells */
static void precompute_subshell_vals(const float * const param)
{
  const unsigned nosubshells = (unsigned)(param[PARAMS_SUBSHELLS]);
  const int interpolate = param[PARAMS_INTERPOLATE] > 1e-10;

  int shell;
  for(shell = 0; shell != total_shells; ++shell)
    {
      int subshell;
      for(subshell = 0; subshell != nosubshells; ++subshell)
	{
	  /* number of the subshell */
	  const unsigned index = subshell + shell*nosubshells;

	  /* whether to interpolate temperatures or just set as constant */
	  if( interpolate )
	    {
	      /* calculate average radius of subshell and shell it is in */
	      const double delta = shell_radii_arcsec[shell+1] - 
		shell_radii_arcsec[shell];
	      const double radiussubshell = shell_radii_arcsec[shell] +
		delta*((subshell+0.5)/nosubshells);

	      /* find nearest two shells to the shell we're looking at */
	      /* shell1 < shell2 */
	      int shell1, shell2;

	      if(shell == 0)
		{
		  shell1 = 0; shell2 = 1;
		}
	      else if(shell == total_shells-1)
		{
		  shell1 = total_shells-2; shell2 = total_shells-1;
		}
	      else
		{
		  const double radiusshell = 0.5*(shell_radii_arcsec[shell] +
						  shell_radii_arcsec[shell+1]);

		  if(radiussubshell < radiusshell)
		    {
		      shell1 = shell-1; shell2 = shell;
		    }
		  else
		    {
		      shell1 = shell; shell2 = shell+1;
		    }
		}

	      {
		/* calculate logorithm of mean shell radii */
		const double logr1 = log( 0.5*(shell_radii_arcsec[shell1]+
					       shell_radii_arcsec[shell1+1]) );
		const double logr2 = log( 0.5*(shell_radii_arcsec[shell2]+
					       shell_radii_arcsec[shell2+1]) );
		const double deltalogr = logr1 - logr2;
		const double logr = log( radiussubshell );
		const double w1 = (logr-logr2)/deltalogr;
		const double w2 = (logr1-logr)/deltalogr;

		/* interpolate in log of radius between shells */
		subshell_ne[index] = exp(log(param[PARAMS_NE+shell1])*w1 +
					 log(param[PARAMS_NE+shell2])*w2);
		subshell_kT[index] = exp(log(param[PARAMS_KT+shell1])*w1 +
					 log(param[PARAMS_KT+shell2])*w2);

		subshell_Z [index] = ( param[PARAMS_Z +shell1]*w1 +
				       param[PARAMS_Z +shell2]*w2 );
		subshell_nH[index] = ( param[PARAMS_NH+shell1]*w1 +
				       param[PARAMS_NH+shell2]*w2 );
	      }
	    }
	  else
	    {
	      subshell_ne[index] = param[PARAMS_NE+shell];
	      subshell_kT[index] = param[PARAMS_KT+shell];
	      subshell_Z [index] = param[PARAMS_Z +shell];
	      subshell_nH[index] = param[PARAMS_NH+shell];
	    }
/* 	  const double delta = shell_radii_arcsec[shell+1] -  */
/* 	    shell_radii_arcsec[shell]; */
/* 	  const double radiussubshell = shell_radii_arcsec[shell] + */
/* 	    delta*((subshell+0.5)/nosubshells); */
/*  	  printf("%e %e %e %e %e\n", radiussubshell, subshell_ne[index], subshell_kT[index], */
/*  		 subshell_Z[index], subshell_nH[index]); */
	} /* loop subshells */
    } /* loop shells */
}

/* apply scattering to spectrum */
static void apply_scattering( unsigned shell,
			      float* const spectrum,
			      const float* const scat_spectrum,
			      const double deltay_cm )
{
  /* iterate over each spectral element */
  float* const scattered_spec = shell_scattered_spectra +
    shell*ENERGY_NOBINS;
  unsigned i;

  for( i=0; i != ENERGY_NOBINS; ++i )
    {
      const double absorb = deltay_cm * scat_spectrum[i];
      
      /* if there is significant absorption */
      if( absorb > 1e-6 )
	{
	  const double orig = spectrum[i];
	  /* optimisation for exp(-x) if x is small */
	  /* better hope absorb >= 0 */
	  const double afactor = exp(-absorb);
	  const double newval  = orig * afactor;
	  
	  /* remove the scattered radiation */
	  spectrum[i] = newval;
	  
	  /* keep track of the scattered radiation */
	  scattered_spec[i] += orig-newval;
	}
    }
}

/* add emission along a line-of-sight */
static void add_los_emission( const unsigned this_subshell,
			      const float * const param )
{
  /* whether scattering out of LoS is on */
  const int scatsw = param[PARAMS_SCATSW] > 1e-10;

  /* we need this later to convert norm of 1 into real norm */
  const double dist_factor_cm = angular_diameter_dist_cm *
    (1. + param[PARAMS_REDSHIFT]);
  const double norm_const = 1e-14 /
    (4.*M_PI * dist_factor_cm*dist_factor_cm);

  /* radii subshell is between */
  const double rho1_cm = subshell_radii_cm[this_subshell];
  const double rho2_cm = subshell_radii_cm[this_subshell+1];

  /* summed output spectrum as we iterate through */
  static float this_spec[ENERGY_NOBINS];

  int shell;
  int loop_zero_twice = 1; /* help us add two of the central bin */

  /* blank output */
  blank_spectrum(this_spec);

  for(shell = -(int)(total_subshells-1); shell < (int)(total_subshells); ++shell)
    {
      const unsigned absshell = abs(shell);

      const double vol_cm3 = projection_vol( subshell_radii_cm[absshell],
					     subshell_radii_cm[absshell+1],
					     rho1_cm, rho2_cm );

      /* alternatively: (less accurate for small nums of shells) */
      /*
      const double deltay_cm = projection_len( subshell_radii_cm[absshell],
					       subshell_radii_cm[absshell+1],
					       0.5*(rho1_cm+rho2_cm) );
      */

      /* whether there's actually volume to project */
      if( vol_cm3 > 1. )
	{
	  const double deltay_cm = vol_cm3 / (M_PI * ( rho2_cm * rho2_cm -
						       rho1_cm * rho1_cm ));

	  const double ne_cm3 = subshell_ne[absshell];
	  const double norm_factor = norm_const * ne_cm3 * ne_cm3 *
	    NH_NE * vol_cm3;
	  const double NH = subshell_nH[absshell] * (deltay_cm*(1./KPC_CM));

	  const unsigned delta = ENERGY_NOBINS*absshell;
	  const float* const input = shell_spectra + delta;

	  /* apply scattering out to previous sum */
	  if( scatsw )
	    {
	      apply_scattering( absshell, this_spec,
				shell_absorb_spectra + delta, deltay_cm );
	    }

	  /* apply absorption within the shell */
	  apply_absorption(NH, this_spec);

	  /* add on emission from this shell */
	  add_mult_spectrum(this_spec, input, (float)(norm_factor));
	}

      /* oh dear, we have to count the zeroth shell twice */
      if( shell == 0 && loop_zero_twice )
	{
	  shell -= 1;
	  loop_zero_twice = 0;
	}

    } /* all shells */

  /* now add on to output */
  add_spectrum(shell_output_spectra + this_subshell*ENERGY_NOBINS,
	       this_spec);
}

/* compute the output spectrum from each of the subshells */
static void compute_output_spectra(const float* const param)
{
  reson_xwrite(15, "resfree: recalculating output spectra");

  /* blank output spectra */
  {
    unsigned i;
    for( i = 0; i != total_subshells*ENERGY_NOBINS; ++i )
      shell_output_spectra[i] = 0;
    for( i = 0; i != total_subshells*ENERGY_NOBINS; ++i )
      shell_scattered_spectra[i] = 0;
  }

  /* calculate output spectra */
  {
    unsigned subshell;
    for(subshell = 0; subshell != total_subshells; ++subshell)
      {
	/* add emission along that line of sight */
	add_los_emission( subshell, param );
      }
  }
}

/* asin(x/y) that truncates at PI/2 if x/y exceeds 1 */
static INLINE double asin_ratio_trunc(const double x,
				      const double y)
{
  const double ratio = x/y;
  if( y < 1e-10 || ratio >= 1. )
    return M_PI * 0.5;
  else
    return asin(ratio);
}

/* Rayleigh scattering function */

/* radii of annulus on sky are R1->R2
   radii of emitting shell are B1->B2 */
static double calc_rayleigh_fraction(const double B1, const double B2,
                                     const double R1, const double R2)
{
  const double v =
    (4*sqr(B1)-7*sqr(R1))*trunc_sqrt(sqr(B1)-sqr(R1)) +
    (7*sqr(R1)-4*sqr(B2))*trunc_sqrt(sqr(B2)-sqr(R1)) +
    (7*sqr(R2)-4*sqr(B1))*trunc_sqrt(sqr(B1)-sqr(R2)) +
    (4*sqr(B2)-7*sqr(R2))*trunc_sqrt(sqr(B2)-sqr(R2)) +
    3*cube(R1)*( asin_ratio_trunc(R1,B2) - asin_ratio_trunc(R1,B1) ) +
    3*cube(R2)*( asin_ratio_trunc(R2,B1) - asin_ratio_trunc(R2,B2) );

  return v / (4*(cube(B1) - cube(B2)));
}

/* compute radiation scattered from shells into lines of sight */
/* also take account of absorption on scattered radiation */
static void compute_inscattered_radiation(const float* const param)
{
  unsigned subannulus;

  /* whether scattering out of LoS is on */
  if( param[PARAMS_SCATSW_IN] < 1e-10 )
    return;

  /* iterate over annuli on sky */
  for(subannulus = 0; subannulus < total_subshells; ++subannulus)
    {
      static float losspec[ENERGY_NOBINS];
      const double rho1_cm = subshell_radii_cm[subannulus];
      const double rho2_cm = subshell_radii_cm[subannulus+1];
      int first_zero = 1;
      int subshell;

      blank_spectrum(losspec);

      /* iterate over subshells from back of the cluster to the centre */
      for(subshell = -(int)(total_subshells)+1; subshell < (int)(total_subshells); ++subshell)
	{
	  /* subshells within a subannulus have no effect */
	  const unsigned absshell = abs(subshell);
	  if( absshell >= subannulus )
	    {
	      /* fraction of subshell scattered into subannulus */
	      const float rayleigh_factor =
		calc_rayleigh_fraction( subshell_radii_cm[absshell],
					subshell_radii_cm[absshell+1],
					rho1_cm, rho2_cm );

	      /* work out absorption */
	      const double vol_cm3 = projection_vol( subshell_radii_cm[absshell],
						     subshell_radii_cm[absshell+1],
						     rho1_cm, rho2_cm );
	      const double deltay_cm = vol_cm3 / (M_PI * ( rho2_cm * rho2_cm -
							   rho1_cm * rho1_cm ));
	      const float NH = subshell_nH[absshell] * (deltay_cm*(1./KPC_CM));

	      apply_absorption(NH, losspec);

	      /* add on scattered ratiation into this annulus */
	      /* 0.5 is to account for that we're only looking at one side of cluster */
	      add_mult_spectrum(losspec, shell_scattered_spectra + absshell*ENERGY_NOBINS,
				rayleigh_factor * 0.5);

	      /* we need to repeat the innermost annulus twice :-( */
	      if( subshell == 0 && first_zero )
		{
		  --subshell;
		  first_zero = 0;
		}
	    }

	} /* loop over subshells */

      /* add total computed spectrum to output */
      add_spectrum(shell_output_spectra + subannulus*ENERGY_NOBINS, losspec);
    } /* loop over subannuli */

}

/* reallocated / recalculate if the parameters to the model change */
static void handle_new_parameters(const float * const param)
{
  static float oldparams[PARAMS_NUMBER] = {-1.};

  /* check whether parameters have changed */
  unsigned i = 0;
  while( i != PARAMS_NUMBER && param[i] == oldparams[i] )
    ++i;

  /* if they have changed */
  if( i != PARAMS_NUMBER )
    {
      reson_xwrite(15, "resfree: calculating model");

      /* update */
      for( i=0; i != PARAMS_NUMBER; ++i )
	oldparams[i] = param[i];

      /* check whether reallocation is necessary */
      realloc_shells(param);

      /* work out cosmology */
      angular_diameter_dist_kpc = 
	calc_ang_diam_dist_kpc( param[PARAMS_REDSHIFT] );
      angular_diameter_dist_cm = KPC_CM * angular_diameter_dist_kpc;

      {
	/* update list of radii if cosmology/redshift has changed,
	   or if we have changed num of subshells */
	static double last_angular_diameter_dist_kpc = -1;
	static float last_inner_radius = -1;
	static int last_subshells = -1;
	if( ! near_double(last_angular_diameter_dist_kpc,
			  angular_diameter_dist_kpc) ||
	    ! near_double(last_inner_radius,
			  param[PARAMS_INNER_ARCSEC]) ||
	    last_subshells != (int)(param[PARAMS_SUBSHELLS]) )
	  {
	    last_angular_diameter_dist_kpc = angular_diameter_dist_kpc;
	    last_inner_radius = param[PARAMS_INNER_ARCSEC];
	    last_subshells = (int)(param[PARAMS_SUBSHELLS]);

	    compute_subshell_radii( param );
	    reson_xwrite(10, "resfree model, "
			 "written by Jeremy Sanders <jss@ast.cam.ac.uk>");
	    reson_xwrite(10, " version %s", RESVERSION);
	    reson_xwrite(10, "Angular diameter distance = %e kpc (%e cm)",
			 angular_diameter_dist_kpc, angular_diameter_dist_cm);
	    reson_xwrite(10, "Number of datasets = %i, number of shells = %i",
			 get_no_datasets(), total_shells);
	    reson_xwrite(15, "Number of subshells = %i, total subshells = %i",
			 last_subshells, total_subshells);

	    precompute_zphabs_absorption( param );
	  }
      }

      /* if calculation isn't forced in each dataset, force it */
      if( ! get_force_calc() ) 
	{
	  reson_xwrite(5, "resfree:: xset forcecalc not set -"
		       " automatically setting it to on");
	  set_force_calc(1);
	}

      /* work out physical parameters in subshells */
      precompute_subshell_vals(param);

      /* precompute shell spectra */
      compute_emission_spectra(param);

      /* precompute optical depths */
      compute_scatter_spectra(param);

      /* compute output spectra */
      compute_output_spectra(param);

      /* add on the contribution due to scattered radiation,
         into the LoS */
      compute_inscattered_radiation(param);
    }
}

static void resfree_model(const float * const ear,
			  const int ne,
			  const float * const param,
			  const int ifl,
			  float * const photar,
			  float * const photer)
{
  static float output_spectrum[ENERGY_NOBINS];
  const float frac_degrees = reson_get_dataset_filter(ifl, 2) / 360.;
  unsigned shell;

  /* recompute if parameters changed */
  handle_new_parameters(param);

  /* blank output spectrum */
  blank_spectrum( output_spectrum );

  {
    const float r = reson_get_dataset_filter(ifl, 1);

    shell = 0;
    while( shell < total_shells &&
	   ! near_double(r, shell_radii_arcsec[shell+1]) )
      shell++;

    if( shell == total_shells )
      {
	reson_xwrite(2, "resfree: Internal error in finding shell");
	return;
      }
  }

  {
    /* work out which subshells to add */
    const unsigned no_subshells = (unsigned)(param[PARAMS_SUBSHELLS]);
    const unsigned startshell = no_subshells * shell;
    const unsigned endshell = no_subshells * (shell+1);
    unsigned subshell;

    /* add up subshell spectra */
    for( subshell = startshell; subshell != endshell; ++subshell )
      {
	add_spectrum( output_spectrum,
		      shell_output_spectra + subshell*ENERGY_NOBINS );
      }
  }

  /* regrade output spectrum to what xspec wants */
  regrade_spectrum( ENERGY_NOBINS, shell_ear, output_spectrum,
		    ne, ear, photar );

  /* renormalise according to degrees */
  if( ! near_double(frac_degrees, 1.) )
    {
      unsigned i;
      for( i = 0; i != ne; ++i )
	photar[i] *= frac_degrees;
    }
}

/* allow resfree_model to be callable from fortran using name reson */
FCALLSCSUB6(resfree_model, RESFREE, resfree, FLOATV, INT, FLOATV, INT, FLOATV,
	    FLOATV)

