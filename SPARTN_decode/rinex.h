#ifndef _RINEX_H_
#define _RINEX_H_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rtcm.h"

#ifdef __cplusplus
extern "C" {
#endif



#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */

static const char syscodes[] = "GREJSCI"; /* satellite system codes */
static const char obscode[] = "CLDS";    /* obs type codes */
static const char frqcodes[] = "1256789"; /* frequency codes */

/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
	int n;                              /* number of index */
	int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
	int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
	unsigned char pri[MAXOBSTYPE];     /* signal priority (15-0) */
	unsigned char type[MAXOBSTYPE];     /* type (0:C,1:L,2:D,3:S) */
	unsigned char code[MAXOBSTYPE];     /* obs code (CODE_L??) */
	double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

int readrnxh(FILE * fp, double * ver, char * type, int * sys, int * tsys, char tobs[][MAXOBSTYPE][4], nav_t * nav, sta_t * sta, gtime_t *start_time, gtime_t *end_time);
int readrnxobsb(FILE * fp, double ver, int * tsys, char tobs[][MAXOBSTYPE][4], int * flag, obsd_t * data, gtime_t * time, sta_t * sta);
int readrnxobs(FILE * fp, double ver, int * tsys, char tobs[][MAXOBSTYPE][4], obs_t * obs, sta_t * sta);
int readrnxobs_one(FILE * fp, double ver, int * tsys, char tobs[][MAXOBSTYPE][4], obs_t * obs, sta_t * sta);
int readrnxnav(FILE * fp, double ver, int sys, nav_t * nav);

static int readrnxfp(FILE * fp, char * type, obs_t * obs, nav_t * nav, sta_t * sta);

static int readrnxfile(const char * file, char * type, obs_t * obs, nav_t * nav, sta_t * sta);

static void init_sta(sta_t * sta);

static int addobsdata(obs_t *obs, const obsd_t *data);
static void decode_navh(char *buff, nav_t *nav);
static void decode_gnavh(char *buff, nav_t *nav);
static void decode_hnavh(char *buff, nav_t *nav);
static void decode_obsh(FILE * fp, char * buff, double ver, int * tsys, char tobs[][MAXOBSTYPE][4], nav_t * nav, sta_t * sta);
static int decode_obsdata(FILE *fp, char *buff, double ver, sigind_t *index, obsd_t *obs);
static void set_index(double ver, int sys, char tobs[MAXOBSTYPE][4], sigind_t * ind);

static int decode_obsepoch(FILE * fp, char * buff, double ver, gtime_t * time, int * flag, int * sats);

static void setstr(char *dst, const char *src, int n);

extern double str2num(const char *s, int i, int n);
static void convcode(double ver, int sys, const char *str, char *type);
extern int satid2no(const char *id);
static int decode_geph(double ver, int sat, gtime_t toc, double *data, geph_t *geph);
static gtime_t adjday(gtime_t t, gtime_t t0);
static int decode_eph(double ver, int sat, gtime_t toc, const double *data, eph_t *eph);
extern int expath(const char *path, char *paths[], int nmax);
extern int sortobs(obs_t *obs);
extern void uniqnav(nav_t *nav);
static int uraindex(double value);

int readrnxt(const char * file, obs_t * obs, nav_t * nav, sta_t * sta);

int readrnxnavb(FILE * fp, double ver, int sys, int * type, eph_t * eph, geph_t * geph);


#ifdef __cplusplus
}
#endif

#endif