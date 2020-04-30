/*------------------------------------------------------------------------------
* tides.c : tidal displacement corrections
*
*          Copyright (C) 2015-2017 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL use IERS tide model
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* version : $Revision:$ $Date:$
* history : 2015/05/10 1.0  separated from ppp.c
*           2015/06/11 1.1  fix bug on computing days in tide_oload() (#128)
*           2017/04/11 1.2  fix bug on calling geterp() in timdedisp()
*-----------------------------------------------------------------------------*/
#include "tides.h"
#define RHO_DEG     180.0 / PI
#define RHO_SEC     3600.0 * 180.0 / PI
#define MJD_J2000   51544.5
#define SQR(x)      ((x)*(x))
#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         398.6005e12     /* earth gravitational constant */
#define GMS         1.3271250e20    /* sun gravitational constant */
#define GMM         4.9027890e12    /* moon gravitational constant */

/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do { \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)



// Rectangular Coordinates -> North, East, Up Components
////////////////////////////////////////////////////////////////////////////
void xyz2neu(const double* Ell, const double* xyz, double* neu) 
{
    double sinPhi = sin(Ell[0]);
    double cosPhi = cos(Ell[0]);
    double sinLam = sin(Ell[1]);
    double cosLam = cos(Ell[1]);

    neu[0] = -sinPhi * cosLam * xyz[0]
        - sinPhi * sinLam * xyz[1]
        + cosPhi * xyz[2];

    neu[1] = -sinLam * xyz[0]
        + cosLam * xyz[1];

    neu[2] = +cosPhi * cosLam * xyz[0]
        + cosPhi * sinLam * xyz[1]
        + sinPhi * xyz[2];
}

/* time to day and sec -------------------------------------------------------*/
static double time2sec(gtime_t time, gtime_t *day)
{
    double ep[6], sec;
    time2epoch(time, ep);
    sec = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
    ep[3] = ep[4] = ep[5] = 0.0;
    *day = epoch2time(ep);
    return sec;
}

/* utc to gmst -----------------------------------------------------------------
* convert utc to gmst (Greenwich mean sidereal time)
* args   : gtime_t t        I   time expressed in utc
*          double ut1_utc   I   UT1-UTC (s)
* return : gmst (rad)
*-----------------------------------------------------------------------------*/
extern double utc2gmst(gtime_t t, double ut1_utc)
{
    const double ep2000[] = { 2000,1,1,12,0,0 };
    gtime_t tut, tut0;
    double ut, t1, t2, t3, gmst0, gmst;

    tut = timeadd(t, ut1_utc);
    ut = time2sec(tut, &tut0);
    t1 = timediff(tut0, epoch2time(ep2000)) / 86400.0 / 36525.0;
    t2 = t1 * t1; t3 = t2 * t1;
    gmst0 = 24110.54841 + 8640184.812866*t1 + 0.093104*t2 - 6.2E-6*t3;
    gmst = gmst0 + 1.002737909350795*ut;

    return fmod(gmst, 86400.0)*PI / 43200.0; /* 0 <= gmst <= 2*PI */
}

/* astronomical arguments: f={l,l',F,D,OMG} (rad) ----------------------------*/
static void ast_args(double t, double *f)
{
    static const double fc[][5] = { /* coefficients for iau 1980 nutation */
        { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
        { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
        {  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
        { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
        { 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939}
    };
    double tt[4];
    int i, j;

    for (tt[0] = t, i = 1; i < 4; i++) tt[i] = tt[i - 1] * t;
    for (i = 0; i < 5; i++) {
        f[i] = fc[i][0] * 3600.0;
        for (j = 0; j < 4; j++) f[i] += fc[i][j + 1] * tt[j];
        f[i] = fmod(f[i] * AS2R, 2.0*PI);
    }
}

/* iau 1980 nutation ---------------------------------------------------------*/
static void nut_iau1980(double t, const double *f, double *dpsi, double *deps)
{
    static const double nut[106][10] = {
        {   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
        {   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
        {   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
        {   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
        {   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
        {   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
        {   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
        {   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
        {   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
        {   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
        {  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
        {  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
        {   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
        {   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
        {  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
        {  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
        {   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
        {  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
        {   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
        {   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
        {   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
        {   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
        {   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
        {  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
        {   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
        {   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
        {  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
        {   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
        {   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
        {   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
        {   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
        {  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
        {   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
        {   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
        {   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
        {   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
        {  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
        {   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
        {   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
        {   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
        {   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
        {   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
        {   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
        {   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
        {   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
        {   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
        {   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
        {   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
        {   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
        {   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
        {   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
        {   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
        {  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
        {   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
        {   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
        {   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
        {   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
        {  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
        {   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
        {  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
        {   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
        {  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
        {   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
        {   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
        {   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
        {   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
        {   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
        {   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
        {   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
        {   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
        {   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
        {   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
        {  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
        {   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
        {   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
        {   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
        {   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
        {   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
        {  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
        {  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
        {   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
        {  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
        {   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
        {   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
        {   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
        {   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
        {   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
        {  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
        {  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
    };
    double ang;
    int i, j;

    *dpsi = *deps = 0.0;

    for (i = 0; i < 106; i++) {
        ang = 0.0;
        for (j = 0; j < 5; j++) ang += nut[i][j] * f[j];
        *dpsi += (nut[i][6] + nut[i][7] * t)*sin(ang);
        *deps += (nut[i][8] + nut[i][9] * t)*cos(ang);
    }
    *dpsi *= 1E-4*AS2R; /* 0.1 mas -> rad */
    *deps *= 1E-4*AS2R;
}

//static void NutMatrix(double Mjd_TT, double *dpsi, double *deps)
//{
//
//    const double ep[] = { 2000,1,1,12,0,0 };
//    double mjd;
//    mjd = 51544.5 + (timediff(gpst2utc(time), epoch2time(ep))) / 86400.0;
//
//    const double T = (mjd - MJD_J2000) / 36525.0;
//
//    double ls = 2.0*PI*Frac(0.993133 + 99.997306*T);
//    double D = 2.0*PI*Frac(0.827362 + 1236.853087*T);
//    double F = 2.0*PI*Frac(0.259089 + 1342.227826*T);
//    double N = 2.0*PI*Frac(0.347346 - 5.372447*T);
//
//    double dpsi = (-17.200*sin(N) - 1.319*sin(2 * (F - D + N)) - 0.227*sin(2 * (F + N))
//        + 0.206*sin(2 * N) + 0.143*sin(ls)) / RHO_SEC;
//    double deps = (+9.203*cos(N) + 0.574*cos(2 * (F - D + N)) + 0.098*cos(2 * (F + N))
//        - 0.090*cos(2 * N)) / RHO_SEC;
//
//    double eps = 0.4090928 - 2.2696E-4*T;
//
//    return;
//}


/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int geterp(const erp_t *erp, gtime_t time, double *erpv)
{
    const double ep[] = { 2000,1,1,12,0,0 };
    double mjd, day, a;
    int i, j, k;

    trace(4, "geterp:\n");

    if (erp->n <= 0) return 0;

    mjd = 51544.5 + (timediff(gpst2utc(time), epoch2time(ep))) / 86400.0;

    if (mjd <= erp->data[0].mjd) {
        day = mjd - erp->data[0].mjd;
        erpv[0] = erp->data[0].xp + erp->data[0].xpr*day;
        erpv[1] = erp->data[0].yp + erp->data[0].ypr*day;
        erpv[2] = erp->data[0].ut1_utc - erp->data[0].lod*day;
        erpv[3] = erp->data[0].lod;
        return 1;
    }
    if (mjd >= erp->data[erp->n - 1].mjd) {
        day = mjd - erp->data[erp->n - 1].mjd;
        erpv[0] = erp->data[erp->n - 1].xp + erp->data[erp->n - 1].xpr*day;
        erpv[1] = erp->data[erp->n - 1].yp + erp->data[erp->n - 1].ypr*day;
        erpv[2] = erp->data[erp->n - 1].ut1_utc - erp->data[erp->n - 1].lod*day;
        erpv[3] = erp->data[erp->n - 1].lod;
        return 1;
    }
    for (j = 0, k = erp->n - 1; j < k - 1;) {
        i = (j + k) / 2;
        if (mjd < erp->data[i].mjd) k = i; else j = i;
    }
    if (erp->data[j].mjd == erp->data[j + 1].mjd) {
        a = 0.5;
    }
    else {
        a = (mjd - erp->data[j].mjd) / (erp->data[j + 1].mjd - erp->data[j].mjd);
    }
    erpv[0] = (1.0 - a)*erp->data[j].xp + a * erp->data[j + 1].xp;
    erpv[1] = (1.0 - a)*erp->data[j].yp + a * erp->data[j + 1].yp;
    erpv[2] = (1.0 - a)*erp->data[j].ut1_utc + a * erp->data[j + 1].ut1_utc;
    erpv[3] = (1.0 - a)*erp->data[j].lod + a * erp->data[j + 1].lod;
    return 1;
}

/* eci to ecef transformation matrix -------------------------------------------
* compute eci to ecef transformation matrix
* args   : gtime_t tutc     I   time in utc
*          double *erpv     I   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *U        O   eci to ecef transformation matrix (3 x 3)
*          double *gmst     IO  greenwich mean sidereal time (rad)
*                               (NULL: no output)
* return : none
* note   : see ref [3] chap 5
*          not thread-safe
*-----------------------------------------------------------------------------*/
void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst)
{
    const double ep2000[] = { 2000,1,1,12,0,0 };
    static gtime_t tutc_;
    static double U_[9], gmst_;
    gtime_t tgps;
    double eps, ze, th, z, t, t2, t3, dpsi, deps, gast, f[5];
    double R1[9], R2[9], R3[9], R[9], W[9], N[9], P[9], NP[9];
    int i;

    if (fabs(timediff(tutc, tutc_)) < 0.01) { /* read cache */
        for (i = 0; i < 9; i++) U[i] = U_[i];
        if (gmst) *gmst = gmst_;
        return;
    }
    tutc_ = tutc;

    /* terrestrial time */
    tgps = utc2gpst(tutc_);
    t = (timediff(tgps, epoch2time(ep2000)) + 19.0 + 32.184) / 86400.0 / 36525.0;
    t2 = t * t; t3 = t2 * t;

    /* astronomical arguments */
    ast_args(t, f);

    /* iau 1976 precession */
    ze = (2306.2181*t + 0.30188*t2 + 0.017998*t3)*AS2R;
    th = (2004.3109*t - 0.42665*t2 - 0.041833*t3)*AS2R;
    z = (2306.2181*t + 1.09468*t2 + 0.018203*t3)*AS2R;
    eps = (84381.448 - 46.8150*t - 0.00059*t2 + 0.001813*t3)*AS2R;
    Rz(-z, R1); Ry(th, R2); Rz(-ze, R3);
    matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
    matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

    /* iau 1980 nutation */
    nut_iau1980(t, f, &dpsi, &deps);
    Rx(-eps - deps, R1); Rz(-dpsi, R2); Rx(eps, R3);
    matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
    matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

    /* greenwich aparent sidereal time (rad) */
    gmst_ = utc2gmst(tutc_, erpv[2]);
    gast = gmst_ + dpsi * cos(eps);
    gast += (0.00264*sin(f[4]) + 0.000063*sin(2.0*f[4]))*AS2R;

    /* eci to ecef transformation matrix */
    Ry(-erpv[0], R1); Rx(-erpv[1], R2); Rz(gast, R3);
    matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, W);
    matmul("NN", 3, 3, 3, 1.0, W, R3, 0.0, R); /* W=Ry(-xp)*Rx(-yp) */
    matmul("NN", 3, 3, 3, 1.0, N, P, 0.0, NP);
    matmul("NN", 3, 3, 3, 1.0, R, NP, 0.0, U_); /* U=W*Rz(gast)*N*P */

    for (i = 0; i < 9; i++) U[i] = U_[i];
    if (gmst) *gmst = gmst_;
}

/* sun and moon position in eci (ref [4] 5.1.1, 5.2.1) -----------------------*/
void sunmoonpos_eci(gtime_t tut, double *rsun, double *rmoon)
{
    const double ep2000[] = { 2000,1,1,12,0,0 };
    double t, f[5], eps, Ms, ls, rs, lm, pm, rm, sine, cose, sinp, cosp, sinl, cosl;

    t = timediff(tut, epoch2time(ep2000)) / 86400.0 / 36525.0;

    /* astronomical arguments */
    ast_args(t, f);

    /* obliquity of the ecliptic */
    eps = 23.43929111-0.0130042*t;
    sine = sin(eps*D2R); cose = cos(eps*D2R);

    /* sun position in eci */
    if (rsun)
    {
        Ms = 357.5277233 + 35999.05034*t;
        ls = 280.460 + 36000.770*t + 1.914666471*sin(Ms*D2R) + 0.019994643*sin(2.0*Ms*D2R);
        rs = AU * (1.000140612 - 0.016708617*cos(Ms*D2R) - 0.000139589*cos(2.0*Ms*D2R));
        sinl = sin(ls*D2R); cosl = cos(ls*D2R);
        rsun[0] = rs * cosl;
        rsun[1] = rs * cose*sinl;
        rsun[2] = rs * sine*sinl;
    }
    /* moon position in eci */
    if (rmoon)
    {
        lm = 218.32 + 481267.883*t + 6.29*sin(f[0]) - 1.27*sin(f[0] - 2.0*f[3]) +
            0.66*sin(2.0*f[3]) + 0.21*sin(2.0*f[0]) - 0.19*sin(f[1]) - 0.11*sin(2.0*f[2]);
        pm = 5.13*sin(f[2]) + 0.28*sin(f[0] + f[2]) - 0.28*sin(f[2] - f[0]) -
            0.17*sin(f[2] - 2.0*f[3]);
        rm = RE_WGS84 / sin((0.9508 + 0.0518*cos(f[0]) + 0.0095*cos(f[0] - 2.0*f[3]) +
            0.0078*cos(2.0*f[3]) + 0.0028*cos(2.0*f[0]))*D2R);
        sinl = sin(lm*D2R); cosl = cos(lm*D2R);
        sinp = sin(pm*D2R); cosp = cos(pm*D2R);
        rmoon[0] = rm * cosp*cosl;
        rmoon[1] = rm * (cose*cosp*sinl - sine * sinp);
        rmoon[2] = rm * (sine*cosp*sinl + cose * sinp);
    }
}

/* sun and moon position -------------------------------------------------------
* get sun and moon position in ecef
* args   : gtime_t tut      I   time in ut1
*          double *erpv     I   erp value {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *rsun     IO  sun position in ecef  (m) (NULL: not output)
*          double *rmoon    IO  moon position in ecef (m) (NULL: not output)
*          double *gmst     O   gmst (rad)
* return : none
*-----------------------------------------------------------------------------*/
void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,double *rmoon, double *gmst)
{
    gtime_t tut;
    double rs[3], rm[3], U[9], gmst_;

    // trace(4, "sunmoonpos: tutc=%s\n", time_str(tutc, 3));

    tut = timeadd(tutc, erpv[2]); /* utc -> ut1 */

    /* sun and moon position in eci */
    sunmoonpos_eci(tut, rsun ? rs : NULL, rmoon ? rm : NULL);

    /* eci to ecef transformation matrix */
    eci2ecef(tutc, erpv, U, &gmst_);

    /* sun and moon postion in ecef */
    if (rsun) matmul("NN", 3, 1, 3, 1.0, U, rs, 0.0, rsun);
    if (rmoon) matmul("NN", 3, 1, 3, 1.0, U, rm, 0.0, rmoon);
    if (gmst) *gmst = gmst_;
}


// Sun's position
///////////////////////////////////////////////////////////////////////////
//double Sun(double Mjd_TT) 
//{
//    double r_sun[3] = { 0.0 };
//    const double eps = 23.43929111 / RHO_DEG;
//    const double T = (Mjd_TT - MJD_J2000) / 36525.0;
//
//    double M = 2.0*PI * Frac(0.9931267 + 99.9973583*T);
//    double L = 2.0*PI * Frac(0.7859444 + M / 2.0 / PI + (6892.0*sin(M) + 72.0*sin(2.0*M)) / 1296.0e3);
//    double r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2 * M);
//
//    r_Sun << r * cos(L) << r * sin(L) << 0.0; r_Sun = rotX(-eps) * r_Sun;
//
//    return    rotZ(GMST(Mjd_TT))
//        * NutMatrix(Mjd_TT)
//        * PrecMatrix(MJD_J2000, Mjd_TT)
//        * r_Sun;
//}


/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
static void tide_pl(const double *eu, const double *rp, double GMp,
                    const double *pos, double *dr)
{
    const double H3=0.292,L3=0.015;
    double r,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,du,cosp,sinl,cosl;
    int i;
    
    trace(4,"tide_pl : pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    if ((r=norm(rp,3))<=0.0) return;
    
    for (i=0;i<3;i++) ep[i]=rp[i]/r;
    
    K2=GMp/GME*SQR(RE_WGS84)*SQR(RE_WGS84)/(r*r*r);
    K3=K2*RE_WGS84/r;
    latp=asin(ep[2]); lonp=atan2(ep[1],ep[0]);
    cosp=cos(latp); sinl=sin(pos[0]); cosl=cos(pos[0]);
    
    /* step1 in phase (degree 2) */
    p=(3.0*sinl*sinl-1.0)/2.0;
    H2=0.6078-0.0006*p;
    L2=0.0847+0.0002*p;
    a=dot(ep,eu,3);
    dp=K2*3.0*L2*a;
    du=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);
    
    /* step1 in phase (degree 3) */
    dp+=K3*L3*(7.5*a*a-1.5);
    du+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);
    
    /* step1 out-of-phase (only radial) */
    du+=3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1]-lonp);
    du+=3.0/4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1]-lonp));
    
    dr[0]=dp*ep[0]+du*eu[0];
    dr[1]=dp*ep[1]+du*eu[1];
    dr[2]=dp*ep[2]+du*eu[2];
    
    trace(5,"tide_pl : dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}
/* displacement by solid earth tide (ref [2] 7) ------------------------------*/
static void tide_solid(const double *rsun, const double *rmoon,
                       const double *pos, const double *E, double gmst, int opt,
                       double *dr)
{
    double dr1[3],dr2[3],eu[3],du,dn,sinl,sin2l;
    
    trace(3,"tide_solid: pos=%.3f %.3f opt=%d\n",pos[0]*R2D,pos[1]*R2D,opt);
    
    /* step1: time domain */
    eu[0]=E[2]; eu[1]=E[5]; eu[2]=E[8];
    tide_pl(eu,rsun, GMS,pos,dr1);
    tide_pl(eu,rmoon,GMM,pos,dr2);
    
    /* step2: frequency domain, only K1 radial */
    sin2l=sin(2.0*pos[0]);
    du=-0.012*sin2l*sin(gmst+pos[1]);
    
    dr[0]=dr1[0]+dr2[0]+du*E[2];
    dr[1]=dr1[1]+dr2[1]+du*E[5];
    dr[2]=dr1[2]+dr2[2]+du*E[8];
    
    /* eliminate permanent deformation */
    if (opt&8) {
        sinl=sin(pos[0]); 
        du=0.1196*(1.5*sinl*sinl-0.5);
        dn=0.0247*sin2l;
        dr[0]+=du*E[2]+dn*E[1];
        dr[1]+=du*E[5]+dn*E[4];
        dr[2]+=du*E[8]+dn*E[7];
    }
    trace(5,"tide_solid: dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}


/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
    const double args[][5]={
        {1.40519E-4, 2.0,-2.0, 0.0, 0.00},  /* M2 */
        {1.45444E-4, 0.0, 0.0, 0.0, 0.00},  /* S2 */
        {1.37880E-4, 2.0,-3.0, 1.0, 0.00},  /* N2 */
        {1.45842E-4, 2.0, 0.0, 0.0, 0.00},  /* K2 */
        {0.72921E-4, 1.0, 0.0, 0.0, 0.25},  /* K1 */
        {0.67598E-4, 1.0,-2.0, 0.0,-0.25},  /* O1 */
        {0.72523E-4,-1.0, 0.0, 0.0,-0.25},  /* P1 */
        {0.64959E-4, 1.0,-3.0, 1.0,-0.25},  /* Q1 */
        {0.53234E-5, 0.0, 2.0, 0.0, 0.00},  /* Mf */
        {0.26392E-5, 0.0, 1.0,-1.0, 0.00},  /* Mm */
        {0.03982E-5, 2.0, 0.0, 0.0, 0.00}   /* Ssa */
    };
    const double ep1975[]={1975,1,1,0,0,0};
    double ep[6],fday,days,t,t2,t3,a[5],ang,dp[3]={0};
    int i,j;
    
    trace(3,"tide_oload:\n");
    
    /* angular argument: see subroutine arg.f for reference [1] */
    time2epoch(tut,ep);
    fday=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    days=timediff(epoch2time(ep),epoch2time(ep1975))/86400.0+1.0;
    t=(27392.500528+1.000000035*days)/36525.0;
    t2=t*t; t3=t2*t;
    
    a[0]=fday;
    a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R; /* H0 */
    a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
    a[3]=(334.329653+4069.0340329577*t-0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
    a[4]=2.0*PI;
    
    /* displacements by 11 constituents */
    for (i=0;i<11;i++) {
        ang=0.0;
        for (j=0;j<5;j++) ang+=a[j]*args[i][j];
        for (j=0;j<3;j++) dp[j]+=odisp[j+i*6]*cos(ang-odisp[j+3+i*6]*D2R);
    }
    denu[0]=-dp[1];
    denu[1]=-dp[2];
    denu[2]= dp[0];
    
    trace(5,"tide_oload: denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
    const double ep2000[]={2000,1,1,0,0,0};
    double y,y2,y3;
    
    y=timediff(tut,epoch2time(ep2000))/86400.0/365.25;
    
    if (y<3653.0/365.25) { /* until 2010.0 */
        y2=y*y; y3=y2*y;
        *xp_bar= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
        *yp_bar=346.346+1.7896*y-0.10729*y2-0.000908*y3;
    }
    else { /* after 2010.0 */
        *xp_bar= 23.513+7.6141*y; /* (mas) */
        *yp_bar=358.891-0.6287*y;
    }
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
                      double *denu)
{
    double xp_bar,yp_bar,m1,m2,cosl,sinl;
    
    trace(3,"tide_pole: pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    /* iers mean pole (mas) */
    iers_mean_pole(tut,&xp_bar,&yp_bar);
    
    /* ref [7] eq.7.24 */
    m1= erpv[0]/AS2R-xp_bar*1E-3; /* (as) */
    m2=-erpv[1]/AS2R+yp_bar*1E-3;
    
    /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
    cosl=cos(pos[1]);
    sinl=sin(pos[1]);
    denu[0]=  9E-3*sin(pos[0])    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
    denu[1]= -9E-3*cos(2.0*pos[0])*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
    denu[2]=-33E-3*sin(2.0*pos[0])*(m1*cosl+m2*sinl); /* du= Sr      (m) */
    
    trace(5,"tide_pole : denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt,const double *odisp, double *dr)
{
    gtime_t tut;
    double pos[2],E[9],drt[3],denu[3],rs[3],rm[3],gmst,erpv[5]={0};
    int i;
#ifdef IERS_MODEL
    double ep[6],fhr;
    int year,mon,day;
#endif
    
    trace(3,"tidedisp: tutc=%s\n",time_str(tutc,0));
    
    //if (erp) {
    //    geterp(erp,utc2gpst(tutc),erpv);
    //}
    tut=timeadd(tutc,erpv[2]);
    
    dr[0]=dr[1]=dr[2]=0.0;
    
    if (norm(rr,3)<=0.0) return;
    
    pos[0]=asin(rr[2]/norm(rr,3));
    pos[1]=atan2(rr[1],rr[0]);
    xyz2enu(pos,E);
    
    if (opt&1) { /* solid earth tides */
        /* sun and moon position in ecef */
        sunmoonpos(tutc,erpv,rs,rm,&gmst);
        tide_solid(rs,rm,pos,E,gmst,opt,drt);
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    //if ((opt&2)&&odisp) { /* ocean tide loading */
    //    tide_oload(tut,odisp,denu);
    //    matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
    //    for (i=0;i<3;i++) dr[i]+=drt[i];
    //}
    //if ((opt&4)&&erp) { /* pole tide */
    //    tide_pole(tut,pos,erpv,denu);
    //    matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
    //    for (i=0;i<3;i++) dr[i]+=drt[i];
    //}
    //printf("tidedisp: dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}

/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta) < 1E-12&&fabs(mu) < 1E-12) return PI;
    return atan2(-tan(beta), sin(mu)) + PI;
}

/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu, double *yaw)
{
    *yaw = yaw_nominal(beta, mu);
    return 1;
}

/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
    const double *rs, double *exs, double *eys)
{
    double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
    double yaw, cosy, siny, erpv[5] = { 0 };
    int i;

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

    /* beta and orbit angle */
    matcpy(ri, rs, 6, 1);
    ri[3] -= OMGE * ri[1];
    ri[4] += OMGE * ri[0];
    cross3(ri, ri + 3, n);
    cross3(rsun, n, p);
    if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
        !normv3(p, ep)) return 0;
    beta = PI / 2.0 - acos(dot(esun, en, 3));
    E = acos(dot(es, ep, 3));
    mu = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
    if (mu < -PI / 2.0) mu += 2.0*PI;
    else if (mu >= PI / 2.0) mu -= 2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en, es, ex);
    cosy = cos(yaw);
    siny = sin(yaw);
    for (i = 0; i < 3; i++) 
    {
        exs[i] = -siny * en[i] + cosy * ex[i];
        eys[i] = -cosy * en[i] - siny * ex[i];
    }
    return 1;
}

/* phase windup model --------------------------------------------------------*/
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
    const double *rs, const double *rr, double *phw)
{
    double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
    int i;

    if (opt <= 0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;

    /* unit vector satellite to receiver */
    for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
    if (!normv3(r, ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);
    xyz2enu(pos, E);
    exr[0] =  E[1]; exr[1]  = E[4]; exr[2] =  E[7]; /* x = north */
    eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek, eys, eks);
    cross3(ek, eyr, ekr);
    for (i = 0; i < 3; i++) 
    {
        ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];
        dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];
    }
    cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
    if (cosp < -1.0)     cosp = -1.0;
    else if (cosp > 1.0) cosp =  1.0;

    if (fabs(fabs(cosp) - 1.0) < 1.0e-10) return 0;

    ph = acos(cosp) / 2.0 / PI;
    cross3(ds, dr, drs);
    if (dot(ek, drs, 3) < 0.0) ph = -ph;

    *phw = ph + floor(*phw - ph + 0.5); /* in cycle */

    if (*phw >  0.5) *phw -= 1.0;
    if (*phw < -0.5) *phw += 1.0;
    printf("phw: sat=%3i, phw=%.3f\n", sat, *phw);

    return 1;
}

/* phase windup model --------------------------------------------------------*/
extern int model_phw_bnc(gtime_t time, int sat, const char *type, int opt,
    const double *rs, const double *rr, double *phw)
{
    double exs[3], eys[3], sz[3],xsun[3], sx[3], sy[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph, erpv[5] = { 0.0 }, rsun[3];
    double neu[3], rx[3], ry[3];
    int i;

    /* unit vector satellite to receiver */
    for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
    if (!normv3(rs, sz)) return 0;

    for (i = 0; i < 3; i++) sz[i] = - sz[i];

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);
    if (!normv3(rsun, xsun)) return 0;

    cross3(sz, xsun, sy);
    cross3(sy, sz, sx);

    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);

    neu[0] = 1.0;
    neu[1] = 0.0;
    neu[2] = 0.0;
    xyz2neu(pos, neu, rx);

    neu[0] = 0.0;
    neu[1] = -1.0;
    neu[2] = 0.0;
    xyz2neu(pos, neu, ry);

    cross3(sz, sy, eks);
    cross3(sz, ry, ekr);

    /* unit vectors of receiver antenna */
    for (i = 0; i < 3; i++)
    {
        ds[i] = sx[i]  - sz[i] * dot(sz, sx, 3)  - eks[i];
        dr[i] = rx[i]  - sz[i] * dot(sz, rx, 3)  + ekr[i];
    }

    cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
    if (cosp < -1.0)     cosp = -1.0;
    else if (cosp > 1.0) cosp =  1.0;

    ph = acos(cosp) / 2.0 / PI;

    cross3(ds, dr, drs);
    if (dot(sz, drs, 3) < 0.0) ph = -ph;

    *phw = ph + floor(*phw - ph + 0.5); /* in cycle */

    //printf("phw: sat=%3i, phw=%.3f\n", sat, *phw);

    return 1;
}


/* phase windup model --------------------------------------------------------*/
extern int model_phw_sap(gtime_t time, int sat, const double *dSatPrecOrbitEcef_m, const double *dSatVelocity_mps, const double *rr, double *phw)
{
    double dr[3], ds[3], drs[3], r[3], dRcvrPosLLH[3], cosp, ph;
    double dMatNegRotY[9] = { 0, 0, 1, 1, 1, -1, 2, 2, 1};
    double dMatDirVec1[9], dMatDirVec2[9], R2[9], R3[9];
    int i;
    double dSatRcvrUnitVec_p[3];
    double dTempSatRcvrVec[3];

    // unit vector satellite to receiver 
    for (i = 0; i < 3; i++) r[i] = rr[i] - dSatPrecOrbitEcef_m[i];
    if (!normv3(r, dSatRcvrUnitVec_p)) return 0;

    // Compute local receiver coordinate system (Unit vectors)
    ecef2pos(rr, dRcvrPosLLH);

    double dRotValueR2_rad = (dRcvrPosLLH[0] - PI/2);
    double dRotValueR3_rad = (dRcvrPosLLH[1] - PI);

    Ry(dRotValueR2_rad, R2); 
    Rz(dRotValueR3_rad, R3);

    matmul("NN", 3, 3, 3, 1.0, dMatNegRotY, R2, 0.0, dMatDirVec1);
    // Matrix of unit vectors (NEU)
    matmul("NN", 3, 3, 3, 1.0, dMatDirVec1, R3, 0.0, dMatDirVec2);

    // Unit vector b (north) 
    double dTempRcvrNorthVec[3];
    double dRcvrNorthUnitVec_b[3];
    dTempRcvrNorthVec[0] = dMatDirVec2[0];
    dTempRcvrNorthVec[1] = dMatDirVec2[1];
    dTempRcvrNorthVec[2] = dMatDirVec2[2];
    if (!normv3(dTempRcvrNorthVec, dRcvrNorthUnitVec_b)) return 0;

    // Unit vector a (east) 
   double dTempRcvrEastVec[3];
   double dRcvrEastUnitVec_a[3];
   dTempRcvrEastVec[0] = dMatDirVec2[3];
   dTempRcvrEastVec[1] = dMatDirVec2[4];
   dTempRcvrEastVec[2] = dMatDirVec2[5];
   if (!normv3(dTempRcvrEastVec, dRcvrEastUnitVec_a)) return 0;
   //SAPALIB::C_Utilities::UnitVec(dTempRcvrEastVec, dRcvrEastUnitVec_a);

   // Compute local satellite coordinate system (Unit vectors)
   //-----------------------------------------------------------
  // Unit vector k (satellite mass center to Earth center)
   double dSatMcUnitVec_k[3];
   if (!normv3(dSatPrecOrbitEcef_m, dSatMcUnitVec_k)) return 0;

   dSatMcUnitVec_k[0] *= (-1.0);
   dSatMcUnitVec_k[1] *= (-1.0);
   dSatMcUnitVec_k[2] *= (-1.0);

   // Unit vector e (Velocity to Position Satellite) considering earth rotation effect.
// Note: This unit vector replaces the unit vector computed with the sun position.
   double dSatVelUnitVec_e[3];
   double dSatVel_eci[3];

   dSatVel_eci[0] = dSatVelocity_mps[0] - OMGE * dSatPrecOrbitEcef_m[1];
   dSatVel_eci[1] = dSatVelocity_mps[1] + OMGE * dSatPrecOrbitEcef_m[0];
   dSatVel_eci[2] = dSatVelocity_mps[2];
   if (!normv3(dSatVel_eci, dSatVelUnitVec_e)) return 0;

   // Unit vector j (k x e)
   double dXkeUnitVec_j[3];
   cross3(dSatMcUnitVec_k, dSatVelUnitVec_e, dXkeUnitVec_j);

   // Unit vector i (j x k)
   double dXjkUnitVec_i[3];
   cross3(dXkeUnitVec_j, dSatMcUnitVec_k, dXjkUnitVec_i);

   // Compute receiver dipole
   double dDotRcvr_p_a =dot(dSatRcvrUnitVec_p, dRcvrEastUnitVec_a,3);

   double dVecXRcvr_p_b[3];
   cross3(dSatRcvrUnitVec_p, dRcvrNorthUnitVec_b, dVecXRcvr_p_b);

   double dVecRcvrScale_p_dot[3];
   dVecRcvrScale_p_dot[0] = dSatRcvrUnitVec_p[0] * dDotRcvr_p_a;
   dVecRcvrScale_p_dot[1] = dSatRcvrUnitVec_p[1] * dDotRcvr_p_a;
   dVecRcvrScale_p_dot[2] = dSatRcvrUnitVec_p[2] * dDotRcvr_p_a;

   double dRcvrDipole[3];
   dRcvrDipole[0] = (dRcvrEastUnitVec_a[0] - dVecRcvrScale_p_dot[0]) + dVecXRcvr_p_b[0];
   dRcvrDipole[1] = (dRcvrEastUnitVec_a[1] - dVecRcvrScale_p_dot[1]) + dVecXRcvr_p_b[1];
   dRcvrDipole[2] = (dRcvrEastUnitVec_a[2] - dVecRcvrScale_p_dot[2]) + dVecXRcvr_p_b[2];

   // Compute satellite dipole
   double dDotSat_p_a = dot(dSatRcvrUnitVec_p, dXjkUnitVec_i,3);

   double dVecXSat_p_b[3];
   cross3(dSatRcvrUnitVec_p, dXkeUnitVec_j, dVecXSat_p_b);

   double dVecSatScale_p_dot[3];
   dVecSatScale_p_dot[0] = dSatRcvrUnitVec_p[0] * dDotSat_p_a;
   dVecSatScale_p_dot[1] = dSatRcvrUnitVec_p[1] * dDotSat_p_a;
   dVecSatScale_p_dot[2] = dSatRcvrUnitVec_p[2] * dDotSat_p_a;

   double dSatDipole[3];
   dSatDipole[0] = (dXjkUnitVec_i[0] - dVecSatScale_p_dot[0]) - dVecXSat_p_b[0];
   dSatDipole[1] = (dXjkUnitVec_i[1] - dVecSatScale_p_dot[1]) - dVecXSat_p_b[1];
   dSatDipole[2] = (dXjkUnitVec_i[2] - dVecSatScale_p_dot[2]) - dVecXSat_p_b[2];

   // Compute fractional part of cycle deltaPhi
   //------------------------------------------
   // Angle zeta
   double dVecXDipole[3];
   cross3(dSatDipole, dRcvrDipole, dVecXDipole);
   double dZeta = dot(dSatRcvrUnitVec_p, dVecXDipole,3);

   // deltaPhi
   double dSatDipoleNorm  = norm(dSatDipole,3);
   double dRcvrDipoleNorm = norm(dRcvrDipole,3);

   double dTheta = dot(dSatDipole, dRcvrDipole, 3) / (dSatDipoleNorm * dRcvrDipoleNorm);
   double dDeltaPhi_cyc = acos(dTheta) / (2*PI);

   if (dZeta < 0.0)
   {
       dDeltaPhi_cyc = dDeltaPhi_cyc * (-1.0);
   }

   // Compute phase wind-up correction
   //---------------------------------
   double dN_cyc = 0.0;
   double dPrevPhaseWindupCorr_cyc = 0.0;
   int   bIsPrevPhaseWindUpValid = 1;

   dPrevPhaseWindupCorr_cyc = *phw;

   if (bIsPrevPhaseWindUpValid)
   {
       // round to nearest integer
       dN_cyc = floor((dPrevPhaseWindupCorr_cyc - dDeltaPhi_cyc) + 0.5);
   }

   *phw = dDeltaPhi_cyc + dN_cyc;
    //*phw = ph + floor(*phw - ph + 0.5); /* in cycle */

    //printf("phw: sat=%3i, phw=%.3f\n", sat, *phw);
    return 1;
}