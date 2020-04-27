#ifndef _TIDES_H_
#define _TIDES_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "rtklib_core.h"
#include "rtcm.h"
#include "gnss_math.h"
#include "model.h"

/* phase windup model --------------------------------------------------------*/
extern int model_phw(gtime_t time, int sat, const char *type, int opt, const double *rs, const double *rr, double *phw);

extern int model_phw_bnc(gtime_t time, int sat, const char *type, int opt, const double *rs, const double *rr, double *phw);

extern void tidedisp(gtime_t tutc, const double *rr, int opt, const double *odisp, double *dr);

#ifdef __cplusplus
}
#endif
#endif