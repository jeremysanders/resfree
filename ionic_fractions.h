#ifndef IONIC_FRACTIONS_H
#define IONIC_FRACTIONS_H

/* get the ionisation fraction of the element and ionisation.
   The element is numbered from 1 (H)
   The ionisation is numbered from 1 (Neutral) */
extern double get_ionisation_fraction(double T_keV,
				      unsigned element,
				      unsigned ionisation);

#endif
