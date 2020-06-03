#ifndef _RTCM_H_
#define _RTCM_H_

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* by Dr. Yudan Yi */
#include "rtklib_core.h"

#ifndef WIN32
#define ARM_MCU
#endif

#define SECONDS_IN_WEEK (604800)

#ifdef ARM_MCU
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-const-variable="
#endif

/* disable the define for embeded */

/*-----------------------------------------------------------*/
/* from rtklib to decode RTCM3 */
#define RTCM2PREAMB 0x66 /* rtcm ver.2 frame preamble */
#define RTCM3PREAMB 0xD3 /* rtcm ver.3 frame preamble */
#define RANGE_MS (CLIGHT * 0.001) /* range in 1 ms */
#define ROUND_U(x) ((unsigned int)floor((x) + 0.5))

#define P2_10 0.0009765625          /* 2^-10 */
#define P2_34 5.820766091346740E-11 /* 2^-34 */
#define P2_46 1.421085471520200E-14 /* 2^-46 */
#define P2_59 1.734723475976810E-18 /* 2^-59 */
#define P2_66 1.355252715606880E-20 /* 2^-66 */

#define P2_5 0.03125                /* 2^-5 */
#define P2_6 0.015625               /* 2^-6 */
#define P2_11 4.882812500000000E-04 /* 2^-11 */
#define P2_15 3.051757812500000E-05 /* 2^-15 */
#define P2_17 7.629394531250000E-06 /* 2^-17 */
#define P2_19 1.907348632812500E-06 /* 2^-19 */
#define P2_20 9.536743164062500E-07 /* 2^-20 */
#define P2_21 4.768371582031250E-07 /* 2^-21 */
#define P2_23 1.192092895507810E-07 /* 2^-23 */
#define P2_24 5.960464477539063E-08 /* 2^-24 */
#define P2_27 7.450580596923828E-09 /* 2^-27 */
#define P2_29 1.862645149230957E-09 /* 2^-29 */
#define P2_30 9.313225746154785E-10 /* 2^-30 */
#define P2_31 4.656612873077393E-10 /* 2^-31 */
#define P2_32 2.328306436538696E-10 /* 2^-32 */
#define P2_33 1.164153218269348E-10 /* 2^-33 */
#define P2_35 2.910383045673370E-11 /* 2^-35 */
#define P2_38 3.637978807091710E-12 /* 2^-38 */
#define P2_39 1.818989403545856E-12 /* 2^-39 */
#define P2_40 9.094947017729280E-13 /* 2^-40 */
#define P2_43 1.136868377216160E-13 /* 2^-43 */
#define P2_48 3.552713678800501E-15 /* 2^-48 */
#define P2_50 8.881784197001252E-16 /* 2^-50 */
#define P2_55 2.775557561562891E-17 /* 2^-55 */

#ifndef NFREQ
#define NFREQ 2
#endif

#ifndef NEXOBS
#define NEXOBS 0
#endif

#define ENAGLO
#define ENAGAL
#define ENABDS

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   40                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   30                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO	    0
#define NSYSGLO     0
#endif

#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   40                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1

#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   199                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS SAIF */
#define MAXPRNQZS_S 189                 /* max satellite PRN number of QZSS SAIF */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif

#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   40                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1

#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNIRN   0	//add zc
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif

#ifdef ENASBSSSS
#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   142                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */
#else
#define MINPRNSBS   0                
#define MAXPRNSBS   0            
#define NSATSBS     0 
#endif

#define NSYS (NSYSGPS + NSYSGLO + NSYSGAL + NSYSCMP) /* only use GPS, GLO, GAL, BDS */

#define MAXSAT		(NSATGPS+NSATGLO+NSATGAL+NSATCMP)

/*----- modified by xuanxuan hu, 12/16/2019 -----*/
#define MAXFREQ 7 /* max NFREQ */

#define FREQ1 1.57542E9      /* L1/E1  frequency (Hz) */
#define FREQ2 1.22760E9      /* L2     frequency (Hz) */
#define FREQ5 1.17645E9      /* L5/E5a frequency (Hz) */
#define FREQ6 1.27875E9      /* E6/LEX frequency (Hz) */
#define FREQ7 1.20714E9      /* E5b    frequency (Hz) */
#define FREQ8 1.191795E9     /* E5a+b  frequency (Hz) */
#define FREQ9 2.492028E9     /* S      frequency (Hz) */
#define FREQ1_GLO 1.60200E9  /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO 0.56250E6  /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO 1.24600E9  /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO 0.43750E6  /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO 1.202025E9 /* GLONASS G3 frequency (Hz) */
#define FREQ1_CMP 1.561098E9 /* BeiDou B1 frequency (Hz) */
#define FREQ2_CMP 1.20714E9  /* BeiDou B2 frequency (Hz) */
#define FREQ3_CMP 1.26852E9  /* BeiDou B3 frequency (Hz) */

#define _SYS_NONE_    0x00                /* navigation system: none */
#define _SYS_GPS_     0x01                /* navigation system: GPS */
#define _SYS_SBS_     0x02                /* navigation system: SBAS */
#define _SYS_GLO_     0x04                /* navigation system: GLONASS */
#define _SYS_GAL_     0x08                /* navigation system: Galileo */
#define _SYS_QZS_     0x10                /* navigation system: QZSS */
#define _SYS_BDS_     0x20                /* navigation system: BeiDou */
#define _SYS_IRN_     0x40                /* navigation system: IRNSS */
#define _SYS_LEO_     0x80                /* navigation system: LEO */
#define _SYS_ALL_     0xFF                /* navigation system: all */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P    (GPS,GLO) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless (GPS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A        (GAL) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P) (GAL,QZS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1SAIF (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1I+Q (GPS,QZS,CMP) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5/E5aI    (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5/E5aQ    (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5/E5aI+Q/L5B+C (GPS,GAL,QZS,IRN,SBS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2I   (GAL,CMP) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2Q   (GAL,CMP) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2I+Q (GAL,CMP) */
#define CODE_L6A    30                  /* obs code: E6A        (GAL) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,CMP) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C    (GAL) */
#define CODE_L6S    35                  /* obs code: LEXS       (QZS) */
#define CODE_L6L    36                  /* obs code: LEXL       (QZS) */
#define CODE_L8I    37                  /* obs code: E5(a+b)I   (GAL) */
#define CODE_L8Q    38                  /* obs code: E5(a+b)Q   (GAL) */
#define CODE_L8X    39                  /* obs code: E5(a+b)I+Q (GAL) */
#define CODE_L2I    40                  /* obs code: B1I        (BDS) */
#define CODE_L2Q    41                  /* obs code: B1Q        (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define MAXCODE     55                  /* max number of obs code */

/* do not define different MAX for post-process verion vs. real-time, this will hide the performance issue, and potential bug */
/* increase MAXOBS, there are satellite in China now */
#define MAXOBS 24
#define GPS_ON
#define GLO_ON
//#define GAL_ON
//#define BDS_ON

//#define RTCM_SSR
//define _PC_
#ifdef _PC_
#define MAXEPH 100
#define MAXEPH_R 24
#else
#define MAXEPH   30
#define MAXEPH_R 24
#endif
#define MAXSSR 55

#define MAXANT 2

#define SQR(x)   ((x)*(x))

#define MAXEXFILE   1024                /* max number of expanded files */
#define NUMSYS      6                   /* number of systems */
#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */
#define FILEPATHSEP '\\'

#define BASE 1
#define ROVER 0
typedef struct {        /* station parameter type */
    char name[MAXANT]; /* marker name */
    char marker[MAXANT]; /* marker number */
    char antdes[MAXANT]; /* antenna descriptor */
    char antsno[MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver[MAXANT]; /* receiver firmware version */
    char recsno[MAXANT]; /* receiver serial number */
    int antsetup;       /* antenna setup id */
    int itrf;           /* ITRF realization year */
    int deltype;        /* antenna delta type (0:enu,1:xyz) */
    double pos[3];      /* station position (ecef) (m) */
    double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;         /* antenna height (m) */
} sta_t;

#ifdef _USE_PPP_
typedef struct
{
    double       Height;  /* m */
    unsigned int Degree;  /* 1-16 */
    unsigned int Order;   /* 1-16 */
    double       Sinus[MAXIONODEGREE][MAXIONOORDER];
    double       Cosinus[MAXIONODEGREE][MAXIONOORDER];
} IonoLayers;

typedef struct
{
    unsigned int EpochTime; /* GPS */
    unsigned int UpdateInterval;
    unsigned int SSRIOD;
    unsigned int NumLayers; /* 1-4 */
    double Quality;
    IonoLayers Layers[NUMIONOLAYERS];
} vtec_t;
#endif 


static const char *obscodes[] = {
    /* observation code strings */

    "", "1C", "1P", "1W", "1Y", "1M", "1N", "1S", "1L", "1E",   /*  0- 9 */
    "1A", "1B", "1X", "1Z", "2C", "2D", "2S", "2L", "2X", "2P", /* 10-19 */
    "2W", "2Y", "2M", "2N", "5I", "5Q", "5X", "7I", "7Q", "7X", /* 20-29 */
    "6A", "6B", "6C", "6X", "6Z", "6S", "6L", "8L", "8Q", "8X", /* 30-39 */
    "2I", "2Q", "6I", "6Q", "3I", "3Q", "3X", "1I", "1Q", "5A", /* 40-49 */
    "5B", "5C", "9A", "9B", "9C", "9X", "", "", "", ""          /* 50-59 */
};

typedef struct {        /* SSR correction type */
    unsigned char sat;
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];      /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data */
    int iodcrc;         /* issue of data crc for beidou/sbas */
    int ura;            /* URA indicator */
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    double deph [3];    /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk [3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;       /* high-rate clock corection (m) */
    double  cbias[NFREQ]; /* code biases (m) */
    double pbias[NFREQ]; /* phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    unsigned char update; /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* navigation data type */
    unsigned int n;         /* number of broadcast ephemeris */
    unsigned int ng;       /* number of glonass ephemeris */
	unsigned int n_gps;
	unsigned int n_gal;
	unsigned int n_bds;
	unsigned int n_qzs;
	unsigned int ns;
    unsigned int nsys[2];
    eph_t eph[MAXEPH];         /* GPS/QZS/GAL ephemeris */
    geph_t geph[MAXEPH_R];     /* GLONASS ephemeris */  
    unsigned char ephsat;
#ifdef RTCM_SSR
    ssr_t ssr[MAXSSR];        /* output of ssr corrections */
#endif
} nav_t;

typedef struct {             /* observation data */
    unsigned int n;          /* number of obervation data/allocated */
    obsd_t data[MAXOBS];     /* observation data records */
    gtime_t time;
    double pos[6];           /* station position (ecef) (m) */
	double refpos[6];        /* reference pos & vel for comparison purpose */
	unsigned char obsflag;   /* obs data complete flag (1:ok,0:not complete) */
    unsigned int staid;      /* station id */
} obs_t;

typedef struct {        /* RTCM control struct type */
    gtime_t time;       /* message time */
    char msmtype[6][128]; /* msm signal types */
    unsigned short lock[MAXSAT][NFREQ+NEXOBS]; /* lock time */
    unsigned int nbyte;          /* number of bytes in message buffer */ 
    unsigned int nbit;           /* number of bits in word buffer */ 
    unsigned int len;            /* message length (bytes) */
    unsigned int type; /* last rtcm type */
    unsigned char buff[1200]; /* message buffer */
	unsigned char key;

	double cp[MAXSAT][NFREQ + NEXOBS]; /* carrier-phase measurement, used in encode */
	gtime_t lltime[MAXSAT][NFREQ + NEXOBS]; /* last lock time */
	int seqno;          /* sequence number for rtcm 2 or iods msm */
} rtcm_t;


#define MAXSTN (2)

typedef struct {
    /* move the observation data struct out of rtcm definiton, to save more memory for PPP only mode */
    obs_t obs[MAXSTN];
    rtcm_t rcv[MAXSTN];
    nav_t  nav;
	double time;
} gnss_rtcm_t;

void trace   (int level, const char *format, ...);

int input_rtcm3_data(rtcm_t *rtcm, unsigned char data, obs_t *obs, nav_t *nav);

/* interface to GNSS db */
int input_rtcm3(unsigned char data, unsigned int stnID, gnss_rtcm_t *gnss);

unsigned int rtcm_getbitu(const unsigned char *buff, int pos, int len);
void setbitu(unsigned char *buff, int pos, int len, unsigned int data);
void setbits(unsigned char *buff, int pos, int len, int data);

/* glo frquent number function */
void set_glo_frq(unsigned char prn, int frq);
 int get_glo_frq(unsigned char prn);

 void set_week_number(int week);
 int  get_week_number();

/* time function */
gtime_t timeadd(gtime_t t, double sec);
double  timediff(gtime_t t1, gtime_t t2);
gtime_t epoch2time(const double *ep);
void    time2epoch(gtime_t t, double *ep);
gtime_t bdt2time(int week, double sec);
double  time2bdt(gtime_t t, int *week);
double  time2gpst(gtime_t t, int *week);
gtime_t utc2gpst(gtime_t t);
gtime_t gpst2utc(gtime_t t);
gtime_t gpst2time(int week, double sec);
gtime_t gpst2bdt(gtime_t t);
gtime_t bdt2gpst(gtime_t t);
void    time2str(gtime_t t, char *s, int n);
char   *time_str(gtime_t t, int n);
gtime_t timeget();
void    timeset(gtime_t t);

void adjweek(gtime_t *time, double tow);


int satno(int sys, int prn);

/* satellite function */
int  satsys(int sat, int *prn);
int  satidx(int sat, int *prn);
char satid (int sat, int *prn);
char sys2char(int sys);
double satwavelen(int sat, int frq);

unsigned char obs2code(int sys, const char * obs, int * freq);

int getcodepri(int sys, unsigned char code, const char * opt);

void ecef2pos(const double *r, double *pos);
void pos2ecef(const double *pos, double *r);

void set_approximate_time(int year, int doy, rtcm_t *rtcm);

int add_obs(obsd_t* obsd, obs_t* obs);
int add_eph(eph_t* eph, nav_t* nav);
int add_geph(geph_t* eph, nav_t* nav);

char *code2obs(int sys, unsigned char code, int *freq);
unsigned int rtk_crc24q(const unsigned char *buff, int len);

/*--------------------------------------------------------------------*/
int gen_rtcm3(rtcm_t* rtcm, obs_t *obs, int type, int sync);

#ifdef __cplusplus
}
#endif
#endif