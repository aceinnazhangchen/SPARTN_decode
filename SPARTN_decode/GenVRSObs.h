#ifndef _GEN_VRS_OBS_H_
#define _GEN_VRS_OBS_H_

#include "rtcm.h"
#include "ephemeris.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


int gen_vobs_from_ssr(obs_t *obs_rov, sap_ssr_t *ssr, gad_ssr_t *gad, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev);

int gen_obs_from_ssr(gtime_t time, double* rcvpos, sap_ssr_t *ssr, gad_ssr_t* gad, vtec_t *vtec, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev, FILE *fLOG);

int gen_rtcm_vrsdata(obs_t * obs, rtcm_t * rtcm, unsigned char * buff);

int read_obs_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int stnID);
int sread_eph_rtcm(unsigned char * buffer, uint32_t len, gnss_rtcm_t * rtcm, uint32_t ns_gps, uint32_t ns_g);
int fread_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int ns_gps, int ns_g);
int sread_ssr_sapcorda(unsigned char * buffer, uint32_t len, raw_spartn_t * spartn, spartn_t * spartn_out, uint32_t * ssr_num);
int fread_ssr_sapcorda(FILE *fSSR, raw_spartn_t *raw_spartn, spartn_t *spartn, uint32_t *ssr_num);
int read_ssr_from_file(FILE *fRTCM, gnss_rtcm_t *rtcm);

 /*--------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif
#endif