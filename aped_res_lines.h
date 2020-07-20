#ifndef APED_RES_LINES_H
#define APED_RES_LINES_H

typedef struct
{
  unsigned element;   /* element lines are for */
  unsigned ion;       /* ionisation state (neutral == 1) */

  unsigned nolines;   /* number of lines */
  double* energy;     /* energy of lines */
  double* oscillator_strength; /* oscillator strengths of lines */

} ion_resonance_line_list;

typedef struct
{
  unsigned no_ions;               /* number of different element/ions */
  ion_resonance_line_list* ions;  /* array of line data */

} ion_resonance_list;

/* return ptr to list of ions, with their resonance lines */
const ion_resonance_list* get_ion_resonance_list(void);

#endif
