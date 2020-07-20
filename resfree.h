#ifndef RESFREE_H
#define RESFREE_H

/* which version of the model this is */
#define RESVERSION "0.21"

/* user definable settings */
#define DATABASEDIR "/home/jss/code/xspec/resonnewscat"
#define XSPECDATADIR "/data/soft3/heasoft6.0/headas/spectral/xspec/manager"
#define APECVERSION "1.3.1"

/* these define where in the parameters we are */
#define PARAMS_NOSHELLS 16
#define PARAMS_NE 0
#define PARAMS_KT (PARAMS_NE + PARAMS_NOSHELLS)
#define PARAMS_Z  (PARAMS_KT + PARAMS_NOSHELLS)
#define PARAMS_NH (PARAMS_Z + PARAMS_NOSHELLS)
#define PARAMS_SCATSW (PARAMS_NH + PARAMS_NOSHELLS)
#define PARAMS_SCATSW_IN (PARAMS_SCATSW + 1)
#define PARAMS_SUBSHELLS (PARAMS_SCATSW_IN + 1)
#define PARAMS_REDSHIFT (PARAMS_SUBSHELLS + 1)
#define PARAMS_TURB (PARAMS_REDSHIFT + 1)
#define PARAMS_INNER_ARCSEC (PARAMS_TURB + 1)
#define PARAMS_INTERPOLATE (PARAMS_INNER_ARCSEC + 1)
#define PARAMS_NUMBER (PARAMS_INTERPOLATE + 1)

#endif
