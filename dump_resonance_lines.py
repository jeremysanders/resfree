#!/usr/bin/env python

import string
import time

import pyfits
from numarray import *
from numarray.ieeespecial import *

atomdb = '/data/soft3/atomdb/v1.3.1'

# some defines
elementnames = [ "H", "He", "Li", "Be", "B", "C", "N", "O",
                 "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
                 "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
                 "Cr", "Mn", "Fe", "Co", "Ni" ]

romannumerals = [ 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
                  'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI',
                  'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII',
                  'XXIV', 'XXV', 'XXVI', 'XXVII', 'XXVIII' ]

# physical constants
A_keV = 12.39854
permittivity_vacuum_F_m_1 = 8.854187817e-12
electron_charge_C = 1.60217733e-19
electron_mass_kg = 9.1093897e-31
speed_light_m_s_1 = 2.99792458e8

# conversion factor from Einstein-A to oscillator strength
# Aji = -3*gamma*fji
# gamma = 1/(4pi epsilon0) * (2/3) * (e^2 w0^2 / m / c^3)
#       = gamma_factor * w0^2
# (w0 is angular frequency of photon)
gamma_factor = 1. / (4.*pi*permittivity_vacuum_F_m_1) * (2./3.) * \
               electron_charge_C **2 / electron_mass_kg / \
               speed_light_m_s_1**3

# where filenames for lines and levels are stored
linefile_map = {}
levelfile_map = {}

def read_filemap():
    """Get the list of LA_*_*.fits files used by APEC."""
    
    file = open( '%s/filemap' % atomdb, 'r' )
    for line in file:
        type, element, ion, filename = string.split(line)

        type = int(type)
        element = int(element)
        ion = int(ion)

        if type == 2:
            # we have an LV file (which contains the energy levels)
            # expand ATOMDB with real filename
            f = string.replace(filename, '$ATOMDB', atomdb)
            levelfile_map[ ( element, ion ) ] = f

        elif type == 3:
            # we have a LA file (which contains einstein A)
            f = string.replace(filename, '$ATOMDB', atomdb)
            linefile_map[ ( element, ion ) ] = f

def get_level_weights( element, ion ):
    """Return a list of statistical weights for the energy levels."""
    
    filename = levelfile_map[ (element, ion) ]
    print "Opening %s (levels) for element=%i, ion=%i" % \
          (filename, element, ion)

    fits = pyfits.open(filename, mode='readonly')
    tabhdu = fits[1]
    tabdat = tabhdu.data
    weights = tabdat.field('LEV_DEG')
    fits.close()

    return weights

def get_groundstate_transitions( outfile, element, ion ):
    """Collect a list of transitions from the groundstate."""

    filename = linefile_map[ (element, ion) ]
    print "Opening %s for element=%i, ion=%i" % (filename, element, ion)

    # write header to output file
    outfile.write( '\n# %s %s\n' %
                   ( elementnames[element-1], romannumerals[ ion-1 ] ) )
    outfile.write( '# %s\n' % filename )
    outfile.write( '@ %i %i\n' % (element, ion) )

    # open the file containing the list of lines
    fits = pyfits.open(filename, mode='readonly')

    tabhdu = fits[1]
    tabdat = tabhdu.data

    # get weights corresponding to the energy levels
    levelweights = get_level_weights( element, ion )

    # get upper and lower levels, compute statistical weights
    lowerlevels = tabdat.field('LOWER_LEV')
    upperlevels = tabdat.field('UPPER_LEV')
    weights = levelweights[ upperlevels - 1 ].astype(Float64) / \
              levelweights[ lowerlevels - 1 ]

    # get Einstein A, wavelength of line (observed & computed)
    einsteinAs = tabdat.field('EINSTEIN_A')
    wavelengths_observed_A = tabdat.field('WAVE_OBS')
    wavelengths_calc_A = tabdat.field('WAVELEN')
    fits.close()

    # replace unobserved lines (marked as nan) with calculated line energies
    wavelengths_A = wavelengths_observed_A.copy()
    nanlist = getnan( wavelengths_A )
    wavelengths_A[ nanlist ] = wavelengths_calc_A[ nanlist ]

    # convert to energy (keV)
    energies_keV = A_keV / wavelengths_A

    # angular freqencies of lines
    w0s = 2. * pi * speed_light_m_s_1 / (wavelengths_A * 1e-10)
    # gamma (to convert to oscillator strengths)
    gammas = gamma_factor * w0s**2
    # oscillator strengths
    oscillator_strengths = einsteinAs * weights / (3.*gammas)

    # select those which involve the ground state
    gs = (lowerlevels == 1)
    out = array( (energies_keV[gs], wavelengths_A[gs],
                  oscillator_strengths[gs]) )
    out = swapaxes(out, 0, 1)
    out = out.tolist()
    out.sort()

    # write out [energy (keV), wavelength (angstroms), oscillator strength]
    for l in out:
        outfile.write('%e %e %e\n' % (l[0], l[1], l[2]) )

read_filemap()

outfile = open('aped_resonance_lines.dat', 'w')
outfile.write('# APED resonance line list\n')
outfile.write('#  Input directory: %s\n' % atomdb)
outfile.write('#  Ran on %s\n' % time.strftime("%a, %d %b %Y %H:%M:%S +0000",
                                               time.gmtime()) )

ions = linefile_map.keys()
ions.sort()
for ion in ions:
    get_groundstate_transitions( outfile, *ion )
