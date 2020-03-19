#ifndef _EPHEMERIS_H_
#define _EPHEMERIS_H_

#include "rtcm.h"
#include "spartn.h"

#ifdef __cplusplus
extern "C" {
#endif

/* by Dr. Yudan Yi */
/*-----------------------------------------------------------*/
#define EPHOPT_BRDC    0                   /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC    1                   /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS    2                   /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC  3                   /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM  4                   /* ephemeris option: broadcast + SSR_COM */
#define EPHOPT_LEX     5                   /* ephemeris option: QZSS LEX ephemeris */
#define EPHOPT_SSRSAP  6                   /* ephemeris option: broadcast + SAPCORDA SSR */
    
typedef struct {
	int	   sat;        /*prn*/
    double rs[6];
    double dts[2];
    double var;
    int svh;    
	double azel[2];    /*azimuth,elevation*/
	double e[3];       /*partial deviation*/
	double tgd;        /* tgd*/
	double r;          /* vector */
	double rate;
	double tro;        /* tropospheric */
}vec_t;

/* compute satellit position */
void satposs(obs_t *obs, vec_t *vec, nav_t *nav, int ephopt);

/* compute satellit position using Sapcorda SSR*/
void satposs_sap(obs_t *obs, vec_t *vec, nav_t *nav, sap_ssr_t *ssr, int ephopt);

int compute_vector_data(obs_t* obs, vec_t* vec);

/*--------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif
#endif
