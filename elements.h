#ifndef ELEMENTS_H
#define ELEMENTS_H

#define NO_ELEMENTS 28

/* list of elements (padded to two chars with spaces) (H=0) */
extern const char * const element_list[NO_ELEMENTS];

/* atomic weights (amu) (H=0) */
extern const double atomic_weights[NO_ELEMENTS];

/* relative solar abundances (ANGR) (H=0) */
extern const double solar_abundances[NO_ELEMENTS];

/* get number associated with name (H=1) */
extern unsigned identify_element(const char* name);

/* get fractional line width */
extern double get_line_width_frac(unsigned atomic_no,
				  double temperature_keV,
				  double turbulence_km_s);

#endif
