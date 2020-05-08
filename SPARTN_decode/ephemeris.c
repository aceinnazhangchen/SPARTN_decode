/*------------------------------------------------------------------------------
* ephemeris.c : satellite ephemeris and clock functions
*
*          Copyright (C) 2010-2018 by T.TAKASU, All rights reserved.
*
* references :
*     [1] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Edition 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performanc standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1.3, December, 2016
*     [8] Quasi-Zenith Satellite System Navigation Service Interface Control
*         Specification for QZSS (IS-QZSS) V1.1, Japan Aerospace Exploration
*         Agency, July 31, 2009
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 1.0), China Satellite
*         Navigation office, December 2012
*     [10] RTCM Standard 10403.1 - Amendment 5, Differential GNSS (Global
*         Navigation Satellite Systems) Services - version 3, July 1, 2011
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.1  moved from rtkcmn.c
*                           added api:
*                               eph2clk(),geph2clk(),seph2clk(),satantoff()
*                               satposs()
*                           changed api:
*                               eph2pos(),geph2pos(),satpos()
*                           deleted api:
*                               satposv(),satposiode()
*           2010/08/26 1.2  add ephemeris option EPHOPT_LEX
*           2010/09/09 1.3  fix problem when precise clock outage
*           2011/01/12 1.4  add api alm2pos()
*                           change api satpos(),satposs()
*                           enable valid unhealthy satellites and output status
*                           fix bug on exception by glonass ephem computation
*           2013/01/10 1.5  support beidou (compass)
*                           use newton's method to solve kepler eq.
*                           update ssr correction algorithm
*           2013/03/20 1.6  fix problem on ssr clock relativitic correction
*           2013/09/01 1.7  support negative pseudorange
*                           fix bug on variance in case of ura ssr = 63
*           2013/11/11 1.8  change constant MAXAGESSR 70.0 -> 90.0
*           2014/10/24 1.9  fix bug on return of var_uraeph() if ura<0||15<ura
*           2014/12/07 1.10 modify MAXDTOE for qzss,gal and bds
*                           test max number of iteration for Kepler
*           2015/08/26 1.11 update RTOL_ELPLER 1E-14 -> 1E-13
*                           set MAX_ITER_KEPLER for alm2pos()
*           2017/04/11 1.12 fix bug on max number of obs data in satposs()
*           2018/10/10 1.13 update reference [7]
*                           support ura value in var_uraeph() for galileo
*                           test eph->flag to recognize beidou geo
*                           add api satseleph() for ephemeris selection
*-----------------------------------------------------------------------------*/
#include "rtcm.h"
#include "spartn.h"
#include "ephemeris.h"
#include "gnss_math.h"
#include "model.h"
#include <math.h>

/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))
#define MAXDTOE     7200.0        /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0        /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0       /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0       /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0        /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_SBS 360.0         /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0       /* max time difference to ephem toe (s) for other */
#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */
#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */
#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */
#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */
#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 90.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR    90.0         /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */
#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,SBS */
//    0,0,1,0,0,0
	0,0,0,0,0,0
};

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
static double dot_(const double *a, const double *b, int n)
{
    double c = 0.0;

    while (--n >= 0) c += a[n] * b[n];
    return c;
}

/* variance by ura ephemeris -------------------------------------------------*/
static double var_uraeph(int sys, int ura)
{
    const double ura_value[]={   
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==_SYS_GAL_) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||15<ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* variance by ura ssr (ref [4]) ---------------------------------------------*/
static double var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}

/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
    double t;
    int i, prn, sys = satsys(eph->sat, &prn);
    
#ifdef _TRACE_
    trace(4,"eph2clk : time=%s sat=%c%02d\n",time_str(time,3),sys2char(sys),prn);
#endif    

    t=timediff(time,eph->toc);
    
    for (i=0;i<2;i++) {
        t-=eph->f0+eph->f1*t+eph->f2*t*t;
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}
/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
                    double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,prn,sys=satsys(eph->sat,&prn);
    
#ifdef _TRACE_
	trace(4,"eph2pos : time=%s sat=%c%02d\n",time_str(time,3),sys2char(sys),prn);
#endif
    
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);
    
    switch ((sys=satsys(eph->sat,&prn))) {
        case _SYS_GAL_: mu=MU_GAL; omge=OMGE_GAL; break;
        case _SYS_BDS_: mu=MU_CMP; omge=OMGE_CMP; break;
        default:        mu=MU_GPS; omge=OMGE;     break;
    }
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
#ifdef _TRACE_
        trace(2,"eph2pos: kepler iteration overflow sat=%c%02d\n",sys2char(sys),prn);
#endif
        return;
    }
    sinE=sin(E); cosE=cos(E);
    
#ifdef _TRACE_
    trace(4,"kepler: sat=%c%02d e=%8.5f n=%2d del=%10.3e\n",sys2char(sys),prn,eph->e,n,E-Ek);
#endif
    
    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);
    
    /* beidou geo satellite */
    if (sys==_SYS_BDS_&&(eph->flag==2||(eph->flag==0&&prn<=5))) {
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;
    
    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
    
    /* position and clock error variance */
    *var=var_uraeph(sys,eph->sva);
}
/* glonass orbit differential equations --------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,r2=dot_(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}
/* glonass position and velocity by numerical integration --------------------*/
static void glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;
    
    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
    double t;
    int i,prn,sys=satsys(geph->sat,&prn);

#ifdef _TRACE_    
    trace(4,"geph2clk: time=%s sat=%c%02d\n",time_str(time,3),sys2char(sys),prn);
#endif
    
    t=timediff(time,geph->toe);
    
    for (i=0;i<2;i++) {
        t-=-geph->taun+geph->gamn*t;
    }
    return -geph->taun+geph->gamn*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var)
{
    double t,tt,x[6];
    int i,prn,sys=satsys(geph->sat,&prn);
    
#ifdef _TRACE_
    trace(4,"geph2pos: time=%s sat=%c%02d\n",time_str(time,3),sys2char(sys),prn);
#endif
    
    t=timediff(time,geph->toe);
    
    *dts=-geph->taun+geph->gamn*t;
    
    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];
    
    *var=SQR(ERREPH_GLO);
}

/* select ephememeris --------------------------------------------------------*/
static const eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax,tmin;
    int i,j=-1,prn,sys=satsys(sat,&prn),sel=0;
    
#ifdef _TRACE_
    trace(4,"seleph  : time=%s sat=%c%02d iode=%d\n",time_str(time,3),sys2char(sys),prn,iode);
#endif
    
    sys=satsys(sat,NULL);
    switch (sys) {
        case _SYS_GPS_: tmax=MAXDTOE+1.0    ; sel=eph_sel[0]; break;
        case _SYS_GAL_: tmax=MAXDTOE_GAL    ; sel=eph_sel[2]; break;
        case _SYS_QZS_: tmax=MAXDTOE_QZS+1.0; sel=eph_sel[3]; break;
        case _SYS_BDS_: tmax=MAXDTOE_CMP+1.0; sel=eph_sel[4]; break;
        default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;
    
    for (i=0;i<nav->n;i++) 
    {
        if (nav->eph[i].sat!=sat) continue;
        if (iode>=0&&nav->eph[i].iode!=iode) continue;
        if (sys==_SYS_GAL_&&sel) 
        {
            if (sel==1&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==2&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
        if (iode>=0) 
			return nav->eph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) 
    {
#ifdef _TRACE_
        trace(3,"no broadcast ephemeris: %s sat=%c%02d iode=%3d\n",time_str(time,0),
               sys2char(sys),prn,iode);
#endif
        return NULL;
    }
    return nav->eph+j;
}
/* select glonass ephememeris ------------------------------------------------*/
static const geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1,prn,sys=satsys(sat,&prn);
    
#ifdef _TRACE_
    trace(4,"selgeph : time=%s sat=%c%02d iode=%2d\n",time_str(time,3),sys2char(sys),prn,iode);
#endif
    
    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
#ifdef _TRACE_
        trace(3,"no glonass ephemeris  : %s sat=%c%02d iode=%2d\n",time_str(time,0),
               sys2char(sys),prn,iode);
#endif
        return NULL;
    }
    return nav->geph+j;
}

/* satellite clock with broadcast ephemeris ----------------------------------*/
static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts)
{
    const eph_t  *eph;
    const geph_t *geph;
    int prn,sys=satsys(sat,&prn);
    
#ifdef _TRACE_
    trace(4,"ephclk  : time=%s sat=%c%02d\n",time_str(time,3),sys2char(sys),prn);
#endif
    
    sys=satsys(sat,NULL);
    
    if (sys==_SYS_GPS_||sys==_SYS_GAL_||sys==_SYS_QZS_||sys==_SYS_BDS_) 
    {
        if (!(eph=seleph(teph,sat,-1,nav))) return 0;
        *dts=eph2clk(time,eph);
    }
    else if (sys==_SYS_GLO_) 
    {
        if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,geph);
    }
    else return 0;
    
    return 1;
}
/* satellite position and clock by broadcast ephemeris -----------------------*/
static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  int iode, double *rs, double *dts, double *var, int *svh)
{
    const eph_t  *eph;
    const geph_t *geph;
    double rst[3],dtst[1],tt=1E-3;
    int i,prn,sys=satsys(sat,&prn);
    
#ifdef _TRACE_
	trace(4,"ephpos  : time=%s sat=%c%02d iode=%d\n",time_str(time,3),sys2char(sys),prn,iode);
#endif
    
    sys=satsys(sat,NULL);
    
    *svh=-1;
    
    if (sys==_SYS_GPS_||sys==_SYS_GAL_||sys==_SYS_QZS_||sys==_SYS_BDS_) 
    {
        if (!(eph=seleph(teph,sat,iode,nav))) return 0;
        eph2pos(time,eph,rs,dts,var);
        time=timeadd(time,tt);
        eph2pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==_SYS_GLO_) 
    {
        if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else return 0;
    
    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
    dts[1]=(dtst[0]-dts[0])/tt;
    
    return 1;
}

/* satellite position and clock with ssr correction --------------------------*/
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                      int opt, double *rs, double *dts, double *var, int *svh)
{
    ssr_t *ssr = NULL;
    eph_t *eph;
    double dorb[3], dclk0;
    double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk;
    int i,sys;

	if(nav->ns == 0) return 0;
   
    for (i = 0; i < nav->ns; i++)
    {
        ssr = nav->ssr + i;
        if (sat == ssr->sat) break;
    }
    
    if (!ssr->t0[0].time) {
   /*     trace(2,"no ssr orbit correction: %s sat=%3d\n",time_str(time,0),sat);  */
        return 0;
    }
    if (!ssr->t0[1].time) {
     /*   trace(2,"no ssr clock correction: %s sat=%3d\n",time_str(time,0),sat);*/
        return 0;
    }
    /* inconsistency between orbit and clock correction */
    if (ssr->iod[0]!=ssr->iod[1]) {
        trace(2,"inconsist ssr correction: %s sat=%3d iod=%d %d\n",
              time_str(time,0),sat,ssr->iod[0],ssr->iod[1]);
        *svh=-1;
        return 0;
    }
    t1=timediff(time,ssr->t0[0]);
    t2=timediff(time,ssr->t0[1]);
    t3=timediff(time,ssr->t0[2]);
    
    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1)>MAXAGESSR||fabs(t2)>MAXAGESSR) 
    {
        trace(2,"age of ssr error: %s sat=%3d t=%.0f %.0f\n",time_str(time,0),
              sat,t1,t2);
        *svh=-1;
        return 0;
    }
    if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
    if (ssr->udi[1]>=1.0) t2-=ssr->udi[0]/2.0;
    
    for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
    dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;
    
    /* ssr highrate clock correction (ref [4]) */
    if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) 
    {
        dclk+=ssr->hrclk;
    }
    if (norm(deph,3)>MAXECORSSR||fabs(dclk)>MAXCCORSSR) {
        trace(3,"invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
              time_str(time,0),norm(deph,3),dclk);
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time,teph,sat,nav,ssr->iode,rs,dts,var,svh)) return 0;
    
    /* satellite clock for gps, galileo and qzss */
    sys=satsys(sat,NULL);
    if (sys==_SYS_GPS_||sys==_SYS_GAL_||sys==_SYS_QZS_||sys==_SYS_BDS_) {
        if (!(eph=seleph(teph,sat,ssr->iode,nav))) return 0;
        
        /* satellite clock by clock parameters */
        tk=timediff(time,eph->toc);
        dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
        dts[1]=eph->f1+2.0*eph->f2*tk;
        
        /* relativity correction */
        dts[0]-=2.0*dot_(rs,rs+3,3)/CLIGHT/CLIGHT;
    }
    /* radial-along-cross directions in ecef */
    if (!normv3(rs+3,ea)) return 0;
    cross3(rs,rs+3,rc);
    if (!normv3(rc,ec)) {
        *svh=-1;
        return 0;
    }
    cross3(ea,ec,er);
    
    for (i=0;i<3;i++) {
        dorb[i] = rs[i];
        rs[i]+=-(er[i]*deph[0]+ea[i]*deph[1]+ec[i]*deph[2])+dant[i];
    }

    dclk0 = dts[0];
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    dts[0]+=dclk/CLIGHT;
    
    /* variance by ssr ura */
    *var=var_urassr(ssr->ura);
    
    //printf("satpos_ssr: %s sat=%3d orb=%14.3f %14.3f %14.3f %6.3f %6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%14.4f %6.3f %6.3f var=%6.3f\n",
    //    time_str(time, 2), sat, dorb[0], dorb[1], dorb[2], deph[0], deph[1], deph[2], t1, er[0], er[1], er[2], dclk0* 1E9, dclk,t2, *var);

    return 1;
}


extern void satposs_sap(obs_t *obs, vec_t *vec, nav_t *nav, sap_ssr_t *ssr, int ephopt)
{
    gtime_t teph = obs->time, time[MAXOBS] = { 0 };
    double dt, pr, e[3] = { 0 };
    int i, j, sys = 0, prn = 0, sap_sat, idx=-1;
    for (i = 0; i < obs->n; i++)
    {
        vec[i].sat = obs->data[i].sat;
        sys = satsys(obs->data[i].sat, &prn);
        /* search any pseudorange */
        for (j = 0, pr = 0.0; j < NFREQ; j++) if ((pr = obs->data[i].P[j]) != 0.0) break;
        if (j >= NFREQ)   continue;
        /* transmission time by satellite clock */
        time[i] = timeadd(obs->time, -pr / CLIGHT);
        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i], teph, obs->data[i].sat, nav, &dt))  continue;
        time[i] = timeadd(time[i], -dt);
        idx = -1;
        for (j = 0; j < nav->ns; j++)
        {
            sap_sat = (ssr[j].sys == 0 ? ssr[j].prn : MAXPRNGPS + ssr[j].prn);
            if (obs->data[i].sat == sap_sat)
            {
                idx = j;
                break;
            }
        }
        if (idx == -1) continue;
        /* satellite position and clock at transmission time */
        if (!satpos(time[i], teph, obs->data[i].sat, ephopt, nav, &ssr[idx], vec[i].rs + 0, vec[i].dts + 0, &vec[i].var, &vec[i].svh))
            continue;

        e[0] = vec[i].rs[0] - obs->pos[0];
        e[1] = vec[i].rs[1] - obs->pos[1];
        e[2] = vec[i].rs[2] - obs->pos[2];
        double r = norm(e, 3);
        /* transmission time by satellite clock */
        time[i] = timeadd(obs->time, -r / CLIGHT);

        /* satellite position and clock at transmission time */
        if (!satpos(time[i], teph, obs->data[i].sat, ephopt, nav, &ssr[idx], vec[i].rs + 0, vec[i].dts + 0, &vec[i].var, &vec[i].svh))
            continue;

        /* if no precise clock available, use broadcast clock instead */
        if (vec[i].dts[0] == 0.0)
        {
            if (!ephclk(time[i], teph, obs->data[i].sat, nav, vec[i].dts + 0)) continue;
            vec[i].var = SQR(STD_BRDCCLK);
        }

        //e[0] = vec[i].rs[0] - obs->pos[0];
        //e[1] = vec[i].rs[1] - obs->pos[1];
        //e[2] = vec[i].rs[2] - obs->pos[2];
        //r = norm(e, 3);
        //double dr = fabs(r-pr)/CLIGHT*3.0e3;
        //printf("satpos: %s, %3d, %14.3f, %14.3f, %14.3f,%10.3f, %10.3f, %10.3f, %14.3f\n", time_str(time[i], 6), vec[i].sat, vec[i].rs[0], vec[i].rs[1], vec[i].rs[2], vec[i].rs[3], vec[i].rs[4], vec[i].rs[5], vec[i].dts[0] * CLIGHT);
    }
}

static int match_nav_ssr(nav_t *nav, sap_ssr_t *ssr, int *inav, int *issr)
{
    int n=0, i, j;
    for (i = 0; i < SSR_NUM; ++i)
    {
        if (ssr[i].sys == 0)
        {
            for (j = 0; j < nav->n; j++)
            {
                if (ssr[i].sat == nav->eph[j].sat && ssr[i].iod[0]== nav->eph[j].iode)
                {
                    inav[n] = j;
                    issr[n] = i;
                    ++n;
                }
            }
        }
        else
        {
            for (j = 0; j < nav->ng; j++)
            {
                if (ssr[i].sat == nav->geph[j].sat  && ssr[i].iod[0] == nav->geph[j].iode)
                {
                    inav[n] = j+100;
                    issr[n] = i;
                    ++n;
                }
            }
        }
    }
    return n;
}

extern int satposs_sap_rcv(gtime_t teph, double *rcvpos, vec_t *vec, nav_t *nav, sap_ssr_t *ssr, int ephopt)
{
    gtime_t time[MAXOBS] = { 0 };
    int inav[MAXOBS] = { 0 }, issr[MAXOBS] = { 0 };
    double dt, pr, e[3] = { 0 },r=0.0;
    int i, j, sys = 0, prn = 0, nsat = 0;;
    double rho = 0.0;
    int nobs = match_nav_ssr(nav, ssr, inav, issr);

    for (i = 0; i < nobs; i++)
    {
        vec[i].sat = ssr[issr[i]].sat;
        sys = satsys(ssr[issr[i]].sat, &prn);
        double tt = 0.075, tmpt=0.0;
        while (1)
        {
        // Correction station position due to Earth Rotation
        // -------------------------------------------------
            double dPhi = OMGE * rho /  CLIGHT;
            double xRec = rcvpos[0] * cos(dPhi) - rcvpos[1] * sin(dPhi);
            double yRec = rcvpos[1] * cos(dPhi) + rcvpos[0] * sin(dPhi);
            double zRec = rcvpos[2];

            /* transmission time by satellite clock */
            time[i] = timeadd(teph, -tt);
            /* satellite clock bias by broadcast ephemeris */
            //if (!ephclk(time[i], teph, vec[i].sat, nav, &dt))  continue; // break;
            //time[i] = timeadd(time[i], -dt);
            /* satellite position and clock at transmission time */
            if (!satpos(time[i], teph, vec[i].sat, ephopt, nav, &ssr[issr[i]], vec[i].rs, vec[i].dts, &vec[i].var, &vec[i].svh))
                break;

            rho  = sqrt(SQR(vec[i].rs[0] - xRec) + SQR(vec[i].rs[1] - yRec) + SQR(vec[i].rs[2] - zRec));
            //tmpt = sqrt(SQR(vec[i].rs[0] - rcvpos[0]) + SQR(vec[i].rs[1] - rcvpos[1]) + SQR(vec[i].rs[2] - rcvpos[2])) ;
            tmpt = rho / CLIGHT;
            if (fabs(tmpt - tt) < 1.0e-8)
            {
                break;
            }
            else
                tt = tmpt;
        }
    
        /* if no precise clock available, use broadcast clock instead */
        if (vec[i].dts[0] == 0.0)
        {
            if (!ephclk(time[i], teph, vec[i].sat, nav, vec[i].dts + 0)) continue;
            vec[i].var = SQR(STD_BRDCCLK);
        }
        else
            nsat++;

        //printf("satpos: %s, %3d, %14.3f, %14.3f, %14.3f, %14.3f, %10.3f, %10.3f, %10.3f, %14.3f\n", time_str(time[i], 6), vec[i].sat, rho, vec[i].rs[0], vec[i].rs[1], vec[i].rs[2], vec[i].rs[3], vec[i].rs[4], vec[i].rs[5], vec[i].dts[0] * CLIGHT);

    }
    return nsat;
}


/* satellite position and clock with ssr correction --------------------------*/
static int satpos_sap_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav, const sap_ssr_t *ssr, double *rs, double *dts, double *var, int *svh)
{
    eph_t *eph;
    double t1, t2, t3, er[3], ea[3], ec[3], rc[3], dorb[3], deph[3], dclk0, dclk, dant[3] = { 0 }, tk, cur_time, ssr_time;
    int i, sys, week, sap_sat=0;
    if (ssr->t0[0]==0.0 || ssr->t0[1]==0.0)
    {
        /*     trace(2,"no ssr orbit correction: %s sat=%3d\n",time_str(time,0),sat);  */
        return 0;
    }


    t1 = fmod(time.time + time.sec, 43200.0) - fmod(ssr->t0[0], 43200.0);
    t2 = fmod(time.time + time.sec, 43200.0) - fmod(ssr->t0[1], 43200.0);
    t3 = fmod(time.time + time.sec, 43200.0) - fmod(ssr->t0[2], 43200.0);

    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1) > MAXAGESSR || fabs(t2) > MAXAGESSR)
    {
        trace(2, "age of ssr error: %s sat=%3d t=%.0f %.0f\n", time_str(time, 0),sat, t1, t2);
        *svh = -1;
        return 0;
    }

    for (i = 0; i < 3; i++) deph[i] = ssr->deph[i];
    dclk = ssr->dclk;

    if (norm(deph, 3) > MAXECORSSR || fabs(dclk) > MAXCCORSSR) {
        trace(3, "invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
            time_str(time, 0), norm(deph, 3), dclk);
        *svh = -1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time, teph, sat, nav, ssr->iod[0], rs, dts, var, svh)) return 0;

    /* satellite clock for gps, galileo and qzss */
    sys = satsys(sat, NULL);
    if (sys == _SYS_GPS_ || sys == _SYS_GAL_ || sys == _SYS_QZS_ || sys == _SYS_BDS_) 
    {
        if (!(eph = seleph(teph, sat, ssr->iod[0], nav))) return 0;

        /* satellite clock by clock parameters */
        tk = timediff(time, eph->toc);
        dts[0] = eph->f0 + eph->f1*tk + eph->f2*tk*tk;
        dts[1] = eph->f1 + 2.0*eph->f2*tk;

        /* relativity correction */
        dts[0] -= 2.0*dot_(rs, rs + 3, 3) / CLIGHT / CLIGHT;
    }
    /* radial-along-cross directions in ecef */
    if (!normv3(rs + 3, ea)) return 0;
    cross3(rs, rs + 3, rc);
    if (!normv3(rc, ec)) 
    {
        *svh = -1;
        return 0;
    }
    cross3(ea, ec, er);

    for (i = 0; i < 3; i++) 
    {
        dorb[i] = rs[i];
        rs[i] -= (er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2]);
    }

    dclk0 = dts[0];
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    dts[0] -= dclk / CLIGHT;

    /* variance by ssr ura */
    *var = var_urassr(ssr->ure);

    return 1;
}


/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock
* args   : gtime_t time     I   time (gpst)
*          gtime_t teph     I   time to select ephemeris (gpst)
*          int    sat       I   satellite number
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)
*                               {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts      O   sat clock {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                  const nav_t *nav, const sap_ssr_t *ssr, double *rs, double *dts, double *var,
                  int *svh)
{
    int prn,sys=satsys(sat,&prn);
#ifdef _TRACE_
    trace(4,"satpos  : time=%s sat=%c%02d ephopt=%d\n",time_str(time,3),sys2char(sys),prn,ephopt);
#endif
    
    *svh=0;
    switch (ephopt) 
    {
        case EPHOPT_BRDC  : return ephpos        (time,teph,sat,nav,-1,rs,dts,var,svh);
        case EPHOPT_SSRAPC: return satpos_ssr    (time,teph,sat,nav, 0,rs,dts,var,svh);
        case EPHOPT_SSRCOM: return satpos_ssr    (time,teph,sat,nav, 1,rs,dts,var,svh);
        case EPHOPT_SSRSAP: return satpos_sap_ssr(time, teph, sat, nav,ssr, rs, dts, var, svh);
    }
    *svh=-1;
    return 0;
}
/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : gtime_t teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void satposs(obs_t *obs, vec_t *vec, nav_t *nav, int ephopt)
{
    gtime_t teph = obs->time, time[MAXOBS] = { 0 };
    double dt, pr;
    int i, j, sys = 0, prn = 0;
#ifdef _TRACE_
	trace(3, "satposs : teph=%s n=%d ephopt=%d\n", time_str(teph, 3), obs->n, ephopt);
#endif

    for (i = 0; i < obs->n; i++) {
		vec[i].sat = obs->data[i].sat;
		sys = satsys(obs->data[i].sat, &prn);
        /* search any pseudorange */
        for (j = 0, pr = 0.0; j < NFREQ; j++) if ((pr = obs->data[i].P[j]) != 0.0) break;
        if (j >= NFREQ) {
#ifdef _TRACE_
            trace(2, "no pseudorange %s sat=%c%02d\n", time_str(obs->time, 3), sys2char(sys), prn);
#endif
            continue;
        }
        /* transmission time by satellite clock */
        time[i] = timeadd(obs->time, -pr / CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i], teph, obs->data[i].sat, nav, &dt)) {
#ifdef _TRACE_
            trace(3, "no broadcast clock %s sat=%c%02d\n", time_str(time[i], 3), sys2char(sys), prn);
#endif
            continue;
        }
        time[i] = timeadd(time[i], -dt);

        /* satellite position and clock at transmission time */
        if (!satpos(time[i], teph, obs->data[i].sat, ephopt, nav, NULL, vec[i].rs + 0, vec[i].dts + 0, &vec[i].var, &vec[i].svh))
        {
#ifdef _TRACE_
            trace(3, "no ephemeris %s sat=%3d\n", time_str(time[i], 3), sys2char(sys), prn);
#endif
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (vec[i].dts[0] == 0.0)
        {
            if (!ephclk(time[i], teph, obs->data[i].sat, nav, vec[i].rs + 6)) continue;
            vec[i].var = SQR(STD_BRDCCLK);
        }
    }
    for (i = 0; i < obs->n; i++)
    {
        if (vec[i].rs[0] == 0.0) continue;
		sys = satsys(obs->data[i].sat, &prn);
        printf("%s sat=%c%02d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f\n", time_str(time[i], 6), sys2char(sys), prn, vec[i].rs[0], vec[i].rs[1], vec[i].rs[2]
            , vec[i].dts[0] * 1E9, vec[i].var);
    }
}

/* select satellite ephemeris --------------------------------------------------
* select satellite ephemeris. call it before calling satpos(),satposs().
* args   : int    sys       I   satellite system (SYS_???)
*          int    sel       I   selection of ephemeris
*                                 _SYS_GAL_: 0:any,1:I/NAV,2:F/NAV
*                                 others : undefined
* return : none
* notes  : default ephemeris selection for galileo is any.
*-----------------------------------------------------------------------------*/
extern void satseleph(int sys, int sel)
{
    switch (sys) {
        case _SYS_GPS_: eph_sel[0]=sel; break;
        case _SYS_GLO_: eph_sel[1]=sel; break;
        case _SYS_GAL_: eph_sel[2]=sel; break;
        case _SYS_QZS_: eph_sel[3]=sel; break;
        case _SYS_BDS_: eph_sel[4]=sel; break;
        case _SYS_SBS_: eph_sel[5]=sel; break;
    }
}

extern int compute_vector_data(obs_t *obs, vec_t *vec)
{
    int i,sys, prn;
    int n = 0;
    vec_t  *vecd = NULL;
    double blh[3] = { 0.0 };

    ecef2pos(obs->pos, blh);
    //printf("rcvpos=%13.3f %13.3f %13.3f\n", obs->pos[0], obs->pos[1], obs->pos[2]);
    for (i = 0; i < obs->n; ++i)
    {
        sys = satsys(obs->data[i].sat, &prn);
        vecd = vec + i;
        if (norm(vecd->rs, 3) < 0.01) continue;
        /* compute geometric-range and azimuth/elevation angle */
        vecd->r = geodist(vecd->rs, obs->pos, vecd->e);
        //printf("%s sat=%c%2d rs=%13.3f dts=%12.3f\n", time_str(obs->time, 6), sys2char(sys), prn, vecd->r, vecd->dts[0] * CLIGHT);
		//vecd->r -= CLIGHT * vecd->dts[0];
        vecd->rate = geovel(vecd->rs, obs->pos, vecd->e);
        satazel(blh, vecd->e, vecd->azel);
        /* satellite clock-bias */
        vecd->sat = obs->data[i].sat;
        vecd->tro = tropmodel(blh, vecd->azel, 0.7);
        ++n;
    }
    return n;
}



