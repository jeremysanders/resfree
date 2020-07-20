#ifndef CONSTANTS_H
#define CONSTANTS_H

/* conversion factors */
#define KPC_CM 3.08568025e21
#define KM_CM 1e5

#define C_KM_S 299792.458
#define C_CM_S 29979245800.

#define ME_G 9.10938188e-28
#define AMU_G 1.6605402e-24
#define E_ESU 4.803e-10

#define NE_NH 1.2
#define NH_NE (1./NE_NH)

/* below wrong, but is what xspec uses in apec */
#define KEV_K 1.16059e7
/*#define KEV_K 11.6048e6*/
#define KEV_ERG 1.602177e-9
#define KEV_HZ 2.41798977e17
#define KEV_ANGSTROM 12.39854

#define RADIAN_ARCSEC ((M_PI/180.)/60./60.)
#define PI_SQRT 1.77245385091

/* in what range and resolution are the spectra stored */
#define ENERGY_MIN 0.3
#define ENERGY_MAX 9.
#define ENERGY_NOBINS 100000 /* should be 100000 */
#define ENERGY_SPACING ( (ENERGY_MAX-ENERGY_MIN)/ENERGY_NOBINS )

/* how to downsample absorption spectrum */
#define ABSORB_DOWNSAMPLE 100 /* should be 100 */
#define ABSORB_BINS (ENERGY_NOBINS/ABSORB_DOWNSAMPLE)

/* convert energy to bin */
#define ENERGY_TO_BIN(e) (int)( ((e)-ENERGY_MIN)/ENERGY_SPACING )
#define BIN_TO_ENERGY(e) ENERGY_MIN + ((e)*ENERGY_SPACING)

#endif
