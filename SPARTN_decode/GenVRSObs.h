#ifndef _GEN_VRS_OBS_H_
#define _GEN_VRS_OBS_H__

#include "rtcm.h"
#include "ephemeris.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


int gen_vobs_from_ssr(obs_t *obs_rov, sap_ssr_t *ssr, gad_ssr_t *gad, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev);


 /*--------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif
#endif