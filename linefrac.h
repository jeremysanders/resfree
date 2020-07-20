#ifndef LINEFRAC_H
#define LINEFRAC_H

extern double line_frac(const double linecentre,
			const double linewidth,
			const double energy);


/* interpolate the function exp(-x*x), return -1 if out of range */
extern double exp2_interpolate(const double x);

#endif
