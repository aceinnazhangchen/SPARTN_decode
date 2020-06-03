#include "GenVRSObs.h"
#include <math.h>
#include <string.h>
#include "rtcm.h"
#include "gnss_math.h"
#include "rtklib_core.h"
#include "ephemeris.h"
#include "model.h"
#include "tides.h"

#ifdef ARM_MCU
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#define SC2RAD   3.1415926535898   /* semi-circle to radian (IS-GPS) */
#define AU       149597870691.0    /* 1 AU (m) */
#define AS2R     (D2R / 3600.0)    /* arc sec to radian */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_BDS   3.986004418E14   /* earth gravitational constant   ref [9] */

/*****************************************************************************
* From GLAB
* Name        : gravitationalDelayCorrection
* Description : Obtains the gravitational delay correction for the effect of
*               general relativity (red shift) to the GPS signal
* Parameters  :
* Name                           |Da|Unit|Description
* double  *receiverPosition       I  m    Position of the receiver
* double  *satellitePosition      I  m    Position of the satellite
* Returned value (double)         O  m    Gravitational delay correction
*****************************************************************************/
extern double ShapiroCorrection(const int sys, const double *rcvPos,const double *satPos)
{
    double	receiverModule;
    double	satelliteModule;
    double	distance;
    double  MU = MU_GPS;

    receiverModule  = sqrt(rcvPos[0] * rcvPos[0] + rcvPos[1] * rcvPos[1] + rcvPos[2] * rcvPos[2]);
    satelliteModule = sqrt(satPos[0] * satPos[0] + satPos[1] * satPos[1] + satPos[2] * satPos[2]);
    distance        = sqrt((satPos[0] - rcvPos[0])*(satPos[0] - rcvPos[0]) +
                           (satPos[1] - rcvPos[1])*(satPos[1] - rcvPos[1]) +
                           (satPos[2] - rcvPos[2])*(satPos[2] - rcvPos[2]));

    switch (sys)
    {
    case _SYS_GPS_:
        MU = MU_GPS;
        break;
    case _SYS_GLO_:
        MU = MU_GLO;
        break;
    case _SYS_GAL_:
        MU = MU_GAL;
        break;
    case _SYS_BDS_:
        MU = MU_BDS;
        break;
    default:
        MU = MU_GPS;
        break;
    }
    return 2.0*MU / (CLIGHT*CLIGHT)*log((satelliteModule + receiverModule + distance) / (satelliteModule + receiverModule - distance));
}


double tecu2meter(int sat, int frq)
{
    double tecu2m;
    double w = satwavelen(sat, frq);
    return 40.3*1.0e16 / SQR(CLIGHT / w);
}

void find_nearest_gridpoints_ionocoef(double *blh, int sat, gad_ssr_t *gad, int *areaId, int *gpt_idx)
{
    int i, j, k, l, lon_nc, lat_nc, dlon, dlat;
    double lon_sp, lat_sp;
    double blh_ref[3] = { 0 }, ned[3] = { 0 }, dist2D = 0.0, MinDist2D=1.0e8;
    int gad_areaId = -1;
    gpt_idx[0] = -1;
    gpt_idx[3] = -1;
    gpt_idx[4] = -1;

    for (j = 0; j < RAP_NUM; j++)
    {
        gad_areaId = gad[j].areaId;
        int aid = -1;
        for (i = 0; i < RAP_NUM; i++)
        {
            if (areaId[i] == 0)
            {
                break;
            }
            else if(areaId[i] == gad_areaId)
            {
                aid = i;
                break;
            }
        }
        if (aid == -1) continue;

        lat_nc = gad[j].nc_lat;
        lon_nc = gad[j].nc_lon;
        lat_sp = gad[j].spa_lat;
        lon_sp = gad[j].spa_lon;

        if (blh[0] * R2D > gad[j].rap_lat || blh[1] * R2D < gad[j].rap_lon)                                        continue;
        if (blh[0] * R2D < gad[j].rap_lat -lat_nc * lat_sp || blh[1] * R2D > gad[j].rap_lon + lon_nc * lon_sp)    continue;
        //printf("gad: sat=%3i,%3i,%7.2f,%7.2f,%7.2f,%7.2f\n", 
        //    sat,gad[j].areaId, gad[j].rap_lat, gad[j].rap_lon, gad[j].rap_lat- gad[j].nc_lat*gad[j].spa_lat, gad[j].rap_lon+ gad[j].nc_lon*gad[j].spa_lon);

        blh_ref[2] = blh[2];
        for (k = 0; k < lat_nc; k++)
        {
            for (l = 0; l < lon_nc; l++)
            {
                blh_ref[0] = (gad[j].rap_lat - k * lat_sp)* D2R;
                blh_ref[1] = (gad[j].rap_lon + l * lon_sp)* D2R;
                blhdiff(blh, blh_ref, ned);
                dist2D = sqrt(ned[0] * ned[0] + ned[1] * ned[1])/1.0e3;
                if (dist2D <= MinDist2D)
                {
                    gpt_idx[0] = j;
                    gpt_idx[1] = k;
                    gpt_idx[2] = l;
                    if (ned[0] > 0)  gpt_idx[3] = 1;
                    if (ned[1] > 0)  gpt_idx[4] = 1;
                    MinDist2D = dist2D;
                }
            }
        }
        break;
    }
}

void dist_inv_unit_weighting(double *blh, double *gpt_bl, double *wdi)
{
    int i;
    double blh_ref[3] = { 0.0 }, ned[3] = { 0.0 }, dist2D[4] = { 0.0 }, sumdis2D = 0.0;
    for (i = 0; i < 4; i++)
    {
        blh_ref[0] = gpt_bl[i * 2 + 0];
        blh_ref[1] = gpt_bl[i * 2 + 1];
        blh_ref[2] = blh[2];

        blhdiff(blh, blh_ref, ned);
        dist2D[i] = sqrt(ned[0] * ned[0] + ned[1] * ned[1]);
        sumdis2D += SQR(dist2D[i]);
    }

    for (i = 0; i < 4; i++) wdi[i] = SQR(dist2D[i]) / sumdis2D;
}


void high_prcision_slant_atm_polynomial(gtime_t time, double *blh, int sat, sap_ssr_t *ssr, gad_ssr_t *gad, double *azel, int *gpt_idx, double *stec, double *stro)
{
    int i, j, satidx=-1;
    int areaId;
    double icoef[12]  = { 0.0 };
    double tcoef[12]  = { 0.0 };
    double gpt_bl[8]  = { 0.0 };
    int    gpt_pos[8] = { 0 };
    double lon_sp = 0.0, lat_sp = 0.0, rap_lon=0.0, rap_lat=0.0,acp_lon=0.0, acp_lat=0.0;
    double Ip[4] = { 0.0 }, Tp[4] = { 0.0 };
    double ion = 0.0, trop = 0.0;
    double wdi[4] = { 0.0 };
    double m_h, m_w, Th = 0.0, Tw = 0.0;;
    /* find satid */
    for (i = 0; i < SSR_NUM; i++)
    {
        int sap_sat = (ssr[i].sys == 0 ? ssr[i].prn : MAXPRNGPS + ssr[i].prn);
        if (sap_sat == sat)
        {
            satidx = i;
            break;
        }
    }
    int aid = -1;
    areaId = gad[gpt_idx[0]].areaId;
    for (i = 0; i < RAP_NUM; i++)
    {
        if (ssr[satidx].areaId[i] == 0)
        {
            break;
        }
        else if (ssr[satidx].areaId[i] == areaId)
        {
            aid = i;
            break;
        }
    }
    if (aid == -1)  return;

    lat_sp  = gad[gpt_idx[0]].spa_lat;
    lon_sp  = gad[gpt_idx[0]].spa_lon;
    rap_lat = gad[gpt_idx[0]].rap_lat;
    rap_lon = gad[gpt_idx[0]].rap_lon;
    acp_lat = gad[gpt_idx[0]].rap_lat - gad[gpt_idx[0]].spa_lat*gad[gpt_idx[0]].nc_lat/2.0;
    acp_lon = gad[gpt_idx[0]].rap_lon + gad[gpt_idx[0]].spa_lon*gad[gpt_idx[0]].nc_lon/2.0;

    //if (gpt_idx[3] * gpt_idx[4] > 0)
    //{
    //    if (gpt_idx[3]<0) /*rap-> (1,1)    */
    //    {
    //        gpt_pos[0] = gpt_idx[1]-1;       gpt_pos[1] = gpt_idx[2]-1;
    //        gpt_pos[2] = gpt_idx[1]-1;       gpt_pos[3] = gpt_idx[2];
    //        gpt_pos[4] = gpt_idx[1];         gpt_pos[5] = gpt_idx[2]-1;
    //        gpt_pos[6] = gpt_idx[1];         gpt_pos[7] = gpt_idx[2];
    //    }
    //    else /*rap-> (0,0)    */
    //    {
    //        gpt_pos[0] = gpt_idx[1];         gpt_pos[1] = gpt_idx[2];
    //        gpt_pos[2] = gpt_idx[1];         gpt_pos[3] = gpt_idx[2]+1;
    //        gpt_pos[4] = gpt_idx[1]+1;       gpt_pos[5] = gpt_idx[2];
    //        gpt_pos[6] = gpt_idx[1]+1;       gpt_pos[7] = gpt_idx[2]+1;
    //    }
    //}
    //else 
    //{ 
    //    if (gpt_idx[3] < 0)  /*rap-> (0,1) */
    //    {
    //        gpt_pos[0] = gpt_idx[1];          gpt_pos[1] = gpt_idx[2] - 1;
    //        gpt_pos[2] = gpt_idx[1];          gpt_pos[3] = gpt_idx[2];
    //        gpt_pos[4] = gpt_idx[1] + 1;      gpt_pos[5] = gpt_idx[2] - 1;
    //        gpt_pos[6] = gpt_idx[1] + 1;      gpt_pos[7] = gpt_idx[2];
    //    }
    //    else   /*rap-> (1,0)    */
    //    {
    //        gpt_pos[0] = gpt_idx[1] - 1;      gpt_pos[1] = gpt_idx[2];
    //        gpt_pos[2] = gpt_idx[1] - 1;      gpt_pos[3] = gpt_idx[2] + 1;
    //        gpt_pos[4] = gpt_idx[1];          gpt_pos[5] = gpt_idx[2];
    //        gpt_pos[6] = gpt_idx[1];          gpt_pos[7] = gpt_idx[2] + 1;
    //    }
    //}

    //gpt_bl[0] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[0];
    //gpt_bl[1] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[1];
    //gpt_bl[2] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[2];
    //gpt_bl[3] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[3];
    //gpt_bl[4] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[4];
    //gpt_bl[5] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[5];
    //gpt_bl[6] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[6];
    //gpt_bl[7] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[7];

    //dist_inv_unit_weighting(blh, gpt_bl, wdi);

    icoef[0] = ssr[satidx].stec_coef[0 + aid * 3];
    icoef[1] = ssr[satidx].stec_coef[1 + aid * 3];
    icoef[2] = ssr[satidx].stec_coef[2 + aid * 3];
    tcoef[0] = ssr[satidx].tro_coef[0 + aid * 3];
    tcoef[1] = ssr[satidx].tro_coef[1 + aid * 3];
    tcoef[2] = ssr[satidx].tro_coef[2 + aid * 3];

    //printf("atmcor: sat=%3i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", sat, icoef[0], icoef[1], icoef[2], tcoef[0], tcoef[1], tcoef[2]);
    *stec = 0.0;
    *stro = 0.0;

    *stec = icoef[0] + icoef[1] * (blh[0] * R2D - acp_lat) + icoef[2] * (blh[1] * R2D - acp_lon);
    trop  = tcoef[0] + tcoef[1] * (blh[0] * R2D - acp_lat) + tcoef[2] * (blh[1] * R2D - acp_lon);
    //printf("atmcor: sat=%3i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", sat, Ip[i], Tp[i], wdi[i], gpt_bl[i * 2], blh[0] * R2D, gpt_bl[i * 2 + 1], blh[1] * R2D, acp_lat, acp_lon);

    //for (i = 0; i < 4; i++)
    //{
    //    wdi[i] = 0.25;
    //    Ip[i] = icoef[0] + icoef[1] * (gpt_bl[i * 2]- acp_lat) + icoef[2] * (gpt_bl[i * 2 + 1]- acp_lon);
    //    Tp[i] = tcoef[0] + tcoef[1] * (gpt_bl[i * 2]- acp_lat) + tcoef[2] * (gpt_bl[i * 2 + 1]- acp_lon);
    //    //Ip[i] = icoef[0] + icoef[1] * (gpt_bl[i * 2] - blh[0]*R2D) + icoef[2] * (gpt_bl[i * 2 + 1] - blh[1]*R2D);
    //    //Tp[i] = tcoef[0] + tcoef[1] * (gpt_bl[i * 2] - blh[0]*R2D) + tcoef[2] * (gpt_bl[i * 2 + 1] - blh[1]*R2D);
    //    *stec += wdi[i] * Ip[i];
    //    Tw    += wdi[i] * Tp[i];
    //    printf("atmcor: sat=%3i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", sat, Ip[i], Tp[i], wdi[i], gpt_bl[i * 2], blh[0] * R2D, gpt_bl[i * 2 + 1], blh[1] * R2D, acp_lat, acp_lon);
    //}
    Th = ssr[satidx].ave_htd;
    m_h=tropmapf(time, blh, azel, &m_w);
    //*stro = m_h*Th + m_w * Tw;
    *stro = m_h * Th + m_w * trop;
    //printf("atmcor: sat=%3i,%.3f,%.3f,%.3f,%.3f\n", sat, *stec, *stro, Th, Tw);
}

void compute_high_prcision_atm_corr(int sat, obs_t *obs_vrs, sap_ssr_t *ssr, gad_ssr_t *gad, double *azel, double maskElev, double *stec, double *stro)
{
    int i,j;
    double blh[3] = { 0 };
    int gpt_idx[5] = { 0 };
    ecef2pos(obs_vrs->pos, blh);
    int loc = -1;
    for (i = 0; i < SSR_NUM; i++)
    {
        ssr[i].sat = ssr[i].prn + ssr[i].sys * 40;
        if (ssr[i].sat == sat)
        {
            loc = i;
            break;
        }
    }
    if (loc == -1) return;

    find_nearest_gridpoints_ionocoef(blh, sat, gad, ssr[loc].areaId, gpt_idx);

    //printf("gad: sat=%3i,%3i,", sat, gad[gpt_idx[0]].areaId);
    int nlat = gpt_idx[1];
    int nlon = gpt_idx[2];
    double arp_lat = gad[gpt_idx[0]].rap_lat - nlat * gad[gpt_idx[0]].spa_lat;
    double arp_lon = gad[gpt_idx[0]].rap_lon + nlon * gad[gpt_idx[0]].spa_lon;
    //printf("%6.2f, %6.2f,", arp_lat,                                        arp_lon);
    //printf("%6.2f, %6.2f,", arp_lat + gpt_idx[3] * gad[gpt_idx[0]].spa_lat, arp_lon);
    //printf("%6.2f, %6.2f,", arp_lat,                                        arp_lon + gpt_idx[4] * gad[gpt_idx[0]].spa_lon);
    //printf("%6.2f, %6.2f\n",arp_lat + gpt_idx[3] * gad[gpt_idx[0]].spa_lat, arp_lon + gpt_idx[4] * gad[gpt_idx[0]].spa_lon);
    
    high_prcision_slant_atm_polynomial(obs_vrs->time, blh, sat, ssr, gad, azel, gpt_idx, stec, stro);
}


extern int gen_vobs_from_ssr(obs_t *obs_rov, sap_ssr_t *ssr, gad_ssr_t* gad, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev)
{
    int i,j,prn;
    double cbias[2] = { 0.0 }, pbias[2] = { 0.0 }, dr[3] = { 0.0 };
    double P[2] = { 0.0 }, L[2] = { 0.0 };
    double rs[6] = { 0.0 }, rr[6] = { 0.0 }, phw = 0.0, phw2 = 0.0;
    double stec = 0.0, strop = 0.0;
    double tecu2m1 = 0.0, tecu2m2 = 0.0;
    double w1 = 0.0, w2 = 0.0;
    double grav_delay = 0.0;
    double soltide = 0.0;

    /* earth tides correction */
    tidedisp(obs_rov->time, &obs_rov->pos, 1, NULL, dr);

    obs_vrs->time =  obs_rov->time;
    obs_vrs->n = obs_rov->n;
    for (i = 0; i < obs_rov->n; i++)
    {
       if (norm(vec_vrs[i].rs, 3) < 0.01) continue;

       int sys = satsys(obs_rov->data[i].sat, &prn);
       obs_vrs->data[i].sat = obs_rov->data[i].sat;
       obs_vrs->data[i].sys = sys;
       obs_vrs->data[i].prn = prn;

       soltide = vec_vrs[i].e[0] * dr[0]+ vec_vrs[i].e[1] * dr[1]+ vec_vrs[i].e[2] * dr[2];

       /* phase windup model */
       //model_phw(obs_rov->time, obs_rov->data[i].sat, NULL, 2, vec_vrs[i].rs, obs_vrs->pos, &phw);
       model_phw_bnc(obs_rov->time, obs_rov->data[i].sat, NULL, 2, vec_vrs[i].rs, obs_vrs->pos, &phw);

       /* gravitational delay correction */
       grav_delay = ShapiroCorrection(sys, obs_vrs->pos, vec_vrs[i].rs);

       /* slant tropospheric and ionospheric delay from HPAC*/
       compute_high_prcision_atm_corr(obs_vrs->data[i].sat, obs_vrs, ssr, gad, &vec_vrs[i].azel, maskElev, &stec, &strop);
       tecu2m1 = tecu2meter(obs_vrs->data[i].sat, 0);
       tecu2m2 = tecu2meter(obs_vrs->data[i].sat, 1);

       if (stec == 0.0 || strop == 0.0)  continue;

       int loc = -1;
       for (j = 0; j < SSR_NUM; j++)
       {
           ssr[j].sat = ssr[j].prn + ssr[j].sys * MAXPRNGPS;
           if (ssr[j].sat == obs_vrs->data[i].sat)
           {
               loc = j;
               break;
           }
       }
       if (loc == -1) continue;

       cbias[0] = ssr[loc].cbias[0];
       cbias[1] = ssr[loc].cbias[1];
       pbias[0] = ssr[loc].pbias[0];
       pbias[1] = ssr[loc].pbias[1];

       w1 = satwavelen(obs_vrs->data[i].sat, 0);
       w2 = satwavelen(obs_vrs->data[i].sat, 1);

       obs_vrs->data[i].LLI[0]  = 0;
       obs_vrs->data[i].LLI[1]  = 0;
       obs_vrs->data[i].code[0] = obs_rov->data[i].code[0];
       obs_vrs->data[i].code[1] = obs_rov->data[i].code[1];
       obs_vrs->data[i].SNR[0]  = 140;
       obs_vrs->data[i].SNR[1]  = 140;
       obs_vrs->data[i].P[0]    =  vec_vrs[i].r + strop + tecu2m1 * stec - grav_delay + cbias[0];
       obs_vrs->data[i].P[1]    =  vec_vrs[i].r + strop + tecu2m2 * stec - grav_delay + cbias[1];
       obs_vrs->data[i].L[0]    = (vec_vrs[i].r + strop - tecu2m1 * stec - grav_delay)/w1 + pbias[0] + phw;
       obs_vrs->data[i].L[1]    = (vec_vrs[i].r + strop - tecu2m2 * stec - grav_delay)/w2 + pbias[1] + phw;

       //printf("obs: %12i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%6.3f,%6.3f,%6.3f,%6.3f,%14.4f,%14.4f,%14.4f,%14.4f\n", 
       //obs_vrs->time.time, obs_vrs->data[i].sat, soltide, phw*w1, phw*w2, grav_delay, strop, cbias[0], cbias[1], pbias[0], pbias[1],
       //tecu2m1 * stec, tecu2m2 * stec, obs_vrs->data[i].P[0], obs_vrs->data[i].P[1], obs_vrs->data[i].L[0], obs_vrs->data[i].L[1]);
    }
    //printf("\n");
    return 1;
}

extern int gen_obs_from_ssr(gtime_t time, double* rcvpos, sap_ssr_t *ssr, gad_ssr_t* gad, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev, FILE *fLOG)
{
    obs_t obs_osr = { 0.0 };
    int i, j, prn;
    double cbias[2] = { 0.0 }, pbias[2] = { 0.0 }, dr[3] = { 0.0 }, dotl[3] = { 0.0 };
    double P[2] = { 0.0 }, L[2] = { 0.0 };
    double rs[6] = { 0.0 }, rr[6] = { 0.0 }, phw = 0.0, phw2 = 0.0;
    double stec = 0.0, strop = 0.0;
    double tecu2m1 = 0.0, tecu2m2 = 0.0;
    double w1 = 0.0, w2 = 0.0;
    double grav_delay = 0.0;
    double soltide = 0.0;
    double otload = 0.0;

    /* earth tides correction */
    tidedisp(time, rcvpos, 1, NULL, dr);

    //ComputeSolidEarthTideDisplacement(gpst2utc(time), rcvpos, dr);

    tide_oload_trm(time, rcvpos, dotl);

    obs_vrs->time = time;
    int obstime = obs_vrs->time.time;
    obstime = fmod(obstime, 43200);

    int nobs = 0;
    obs_vrs->obsflag = 1;

    for (i = 0; i < obs_vrs->n; i++)
    {
        if (norm(vec_vrs[i].rs, 3) < 0.01)     continue;

        if (vec_vrs[i].azel[1]*R2D < maskElev) continue;

        int sys = satsys(vec_vrs[i].sat, &prn);
        obs_vrs->data[i].sat = vec_vrs[i].sat;
        obs_vrs->data[i].sys = sys;
        obs_vrs->data[i].prn = prn;

        otload  = vec_vrs[i].e[0] * dotl[0] + vec_vrs[i].e[1] * dotl[1] + vec_vrs[i].e[2] * dotl[2];
        soltide = vec_vrs[i].e[0] * dr[0] + vec_vrs[i].e[1] * dr[1] + vec_vrs[i].e[2] * dr[2];
        vec_vrs[i].r -= soltide;

        /* phase windup model */
        //model_phw(obs_rov->time, obs_rov->data[i].sat, NULL, 2, vec_vrs[i].rs, obs_vrs->pos, &phw);
        //model_phw_bnc(time, vec_vrs[i].sat, NULL, 2, vec_vrs[i].rs, obs_vrs->pos, &phw);
        model_phw_sap(time, vec_vrs[i].sat, vec_vrs[i].rs, vec_vrs[i].rs + 3, obs_vrs->pos, &vec_vrs[i].phw);

        /* gravitational delay correction */
        grav_delay = ShapiroCorrection(sys, obs_vrs->pos, vec_vrs[i].rs);

        /* slant tropospheric and ionospheric delay from HPAC*/
        compute_high_prcision_atm_corr(vec_vrs[i].sat, obs_vrs, ssr, gad, &vec_vrs[i].azel, maskElev, &stec, &strop);
        tecu2m1 = tecu2meter(vec_vrs[i].sat, 0);
        tecu2m2 = tecu2meter(vec_vrs[i].sat, 1);
        if (stec == 0.0 || strop == 0.0)  continue;

        int loc = -1;
        for (j = 0; j < SSR_NUM; j++)
        {
            ssr[j].sat = ssr[j].prn + ssr[j].sys * MAXPRNGPS;
            if (ssr[j].sat == vec_vrs[i].sat)
            {
                loc = j;
                break;
            }
        }
        if (loc == -1) continue;

        cbias[0] = ssr[loc].cbias[0];
        cbias[1] = ssr[loc].cbias[1];
        pbias[0] = ssr[loc].pbias[0];
        pbias[1] = ssr[loc].pbias[1];

        w1 = satwavelen(obs_vrs->data[i].sat, 0);
        w2 = satwavelen(obs_vrs->data[i].sat, 1);

        obs_vrs->data[i].time    = time;
        obs_vrs->data[i].LLI[0]  = 0;
        obs_vrs->data[i].LLI[1]  = 0;
        obs_vrs->data[i].code[0] = CODE_L1C;
        obs_vrs->data[i].code[1] = CODE_L2C;
        obs_vrs->data[i].SNR[0]  = 140;
        obs_vrs->data[i].SNR[1]  = 140;
        obs_vrs->data[i].P[0]    =  vec_vrs[i].r - CLIGHT * vec_vrs[i].dts[0] + strop + tecu2m1 * stec  + cbias[0];          // -grav_delay
        obs_vrs->data[i].P[1]    =  vec_vrs[i].r - CLIGHT * vec_vrs[i].dts[0] + strop + tecu2m2 * stec  + cbias[1];
        obs_vrs->data[i].L[0]    = (vec_vrs[i].r - CLIGHT * vec_vrs[i].dts[0] + strop - tecu2m1 * stec + pbias[0]) / w1 + vec_vrs[i].phw;
        obs_vrs->data[i].L[1]    = (vec_vrs[i].r - CLIGHT * vec_vrs[i].dts[0] + strop - tecu2m2 * stec + pbias[1]) / w2 + vec_vrs[i].phw;
        nobs++;
        int dt1 = obstime - ssr[loc].t0[0];
        int dt2 = obstime - ssr[loc].t0[1];
        int dt3 = obstime - ssr[loc].t0[2];
        int dt4 = obstime - ssr[loc].t0[4];
        //printf("obs:%6i,%3i,%3i,%3i,%3i,%3i,%13.3f,%11.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n",
        //    obstime, dt1, dt2, dt3, dt4, obs_vrs->data[i].sat, vec_vrs[i].r, vec_vrs[i].dts[0]*CLIGHT, soltide, vec_vrs[i].phw*w1, vec_vrs[i].phw*w2, grav_delay, strop,
        //    tecu2m1 * stec, tecu2m2 * stec, cbias[0], cbias[1], pbias[0], pbias[1]);
    }

    double ep[6];
    int sys;
    time2epoch(obs_vrs->time, ep);
    if(fLOG)fprintf(fLOG, "%4.0f %2.0f %2.0f %2.0f %2.0f %4.1f %3i", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],nobs);

    for (i = 0; i < obs_vrs->n; i++)
    {
        if (obs_vrs->data[i].P[0] == 0 || obs_vrs->data[i].P[1] == 0 || obs_vrs->data[i].L[0] == 0 || obs_vrs->data[i].L[1] == 0)   continue;
        sys = satsys(obs_vrs->data[i].sat, &prn);
		if (fLOG)fprintf(fLOG, "%c%02d", sys2char(sys), prn);
    }
	if (fLOG)fprintf(fLOG, "\n");
    for (i = 0; i < obs_vrs->n; i++)
    {
        if (obs_vrs->data[i].P[0] == 0 || obs_vrs->data[i].P[1] == 0 || obs_vrs->data[i].L[0] == 0 || obs_vrs->data[i].L[1] == 0)   continue;
		if (fLOG)fprintf(fLOG, "%14.4f,%14.4f,%14.4f,%14.4f\n", obs_vrs->data[i].P[0], obs_vrs->data[i].P[1], obs_vrs->data[i].L[0], obs_vrs->data[i].L[1]);
    }
    return 1;
}


extern int get_match_epoch(obs_t *obs_ref, obs_t *obs_rov, vec_t *vec_ref, vec_t *vec_rov, int *iref, int *irov)
{
    /* get the matched satellite index from the base and rover recivers, and sort by the elevation */
    int n = 0, i = 0, j = 0;
    double elev[MAXOBS] = { 0.0 };
    obsd_t *pObsRov = NULL;
    obsd_t *pObsRef = NULL;
    vec_t  *pVecRov = NULL;
    vec_t  *pVecRef = NULL;
    for (i = 0; i < (int)obs_rov->n; ++i)
    {
        pObsRov = obs_rov->data + i;
        pVecRov = vec_rov + i;
        if (norm(pVecRov->rs, 3) < 0.1) continue;
        for (j = 0; j < (int)obs_ref->n; ++j)
        {
            pObsRef = obs_ref->data + j;
            pVecRef = vec_ref + j;
            if (norm(pVecRef->rs, 3) < 0.1) continue;
            if (pObsRov->sat != pObsRef->sat) continue;
            iref[n] = j;
            irov[n] = i;
            elev[n] = pVecRov->azel[1];
            ++n;
        }
    }
    for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (elev[i] < elev[j])
            {
                int temp1 = iref[i];
                int temp2 = irov[i];
                double v = elev[i];
                iref[i] = iref[j];
                irov[i] = irov[j];
                elev[i] = elev[j];
                iref[j] = temp1;
                irov[j] = temp2;
                elev[j] = v;
            }
        }
    }
    return n;
}


int gen_rtcm_vrsdata(obs_t* obs, rtcm_t* rtcm, unsigned char* buff)
{
	/* write generated VRS data to file, in realtime system, will write to the user */

	rtcm->time = obs->time;
	int satNUM[4] = { 0 }, ns = 0, type[4] = { 1074, 1084 , 1094, 1124 };
	int size_write = 0, j = 0;

	for (j = 0; j < (int)obs->n; ++j)
	{
		int prn = 0;
		int sys = satsys(obs->data[j].sat, &prn);
		if (sys == _SYS_GPS_) {
			++satNUM[0]; ++ns;
		}
		else if (sys == _SYS_GLO_) {
			++satNUM[1]; ++ns;
		}
		else if (sys == _SYS_GAL_) {
			++satNUM[2]; ++ns;
		}
		else if (sys == _SYS_BDS_) {
			++satNUM[3]; ++ns;
		}
	}
	if (ns > 0)
	{
		if (gen_rtcm3(rtcm, obs, 1005, 0) > 0)
		{
			if (buff)memcpy(buff + size_write, rtcm->buff, rtcm->nbyte);
			size_write += rtcm->nbyte;
		}

		int lastSYS = 0;
		for (j = 0; j < 2; ++j)
		{
			if (satNUM[j] > 0)
				lastSYS = j;
		}

		for (j = 0; j < 2; ++j)
		{
			int syn = lastSYS == j ? 0 : 1;
			if (gen_rtcm3(rtcm, obs, type[j], syn) > 0)
			{
				if (buff)memcpy(buff + size_write, rtcm->buff, rtcm->nbyte);
				size_write += rtcm->nbyte;
			}
		}
	}
	return size_write;
}

int read_obs_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int stnID)
{
	int ret = 0;
	int iter_ssr = 0;
	size_t readCount = 0;
	char buff = ' ';

	while (!feof(fRTCM))
	{
		memset(&buff, 0, sizeof(buff));
		readCount = fread(&buff, sizeof(char), 1, fRTCM);
		if (readCount < 1)
		{
			/* file error or eof of file */
			ret = -1;
			break;
		}

		ret = input_rtcm3(buff, stnID, rtcm);
		if (ret == 1)
		{
			break;
		}
	}
	return ret;
}

int sread_eph_rtcm(unsigned char* buffer, uint32_t len, gnss_rtcm_t *rtcm, uint32_t ns_gps, uint32_t ns_g)
{
	int ret = 0;
	int iter_ssr = 0;
	char c = ' ';
	uint32_t n;
	for (n = 0; n < len; ++n)
	{
		c = buffer[n];
		ret = input_rtcm3(c, 0, rtcm);
		if (ret == 2 && rtcm->nav.n_gps >= ns_gps && rtcm->nav.ng >= ns_g)
		{
			//break;
		}
	}
	return ret;
}

int fread_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int ns_gps, int ns_g)
{
	int ret = 0;
	int iter_ssr = 0;
	size_t readCount = 0;
	char buff = ' ';
	while (!feof(fRTCM))
	{
		memset(&buff, 0, sizeof(buff));
		readCount = fread(&buff, sizeof(char), 1, fRTCM);
		if (readCount < 1)
		{
			/* file error or eof of file */
			break;
		}
		ret = input_rtcm3(buff, 0, rtcm);
		if (ret == 2 && rtcm->nav.n_gps >= ns_gps && rtcm->nav.ng >= ns_g)
		{
			break;
		}
	}
	return ret;
}

int sread_ssr_sapcorda(unsigned char* buffer, uint32_t len, raw_spartn_t *spartn, spartn_t *spartn_out, uint32_t *ssr_num)
{
	int ret = 0;
	uint8_t c = 0;
	uint32_t n, i;

	for (n = 0; n < len; ++n)
	{
		c = buffer[n];
		int ret = input_spartn_data(spartn, spartn_out, c);
		int  areaId = spartn_out->ssr_gad[0].areaId;
		double t1 = spartn_out->ssr[0].t0[0];
		double t2 = spartn_out->ssr[0].t0[1];
		double t3 = spartn_out->ssr[0].t0[2];
		double t4 = spartn_out->ssr[0].t0[3];
		double t5 = spartn_out->ssr[0].t0[4];
		double t6 = spartn_out->ssr[0].t0[5];

		if (ret == 1 && t1*t2*t3*t4*t5*t6 > 0.0 && spartn_out->type == 0 && spartn_out->eos == 1)
		{
			//printf("\n");
			for (i = 0; i < SSR_NUM; i++)
			{
				if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys == 0)
					ssr_num[0]++;
				else if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys == 1)
					ssr_num[1]++;
			}
		}
	}
	return ret;
}

int fread_ssr_sapcorda(FILE *fSSR, raw_spartn_t *spartn, spartn_t *spartn_out, uint32_t *ssr_num)
{
	int ret = 0;
	char buff = ' ';
	size_t currentCount = 0;
	size_t frameSize = 0;
	size_t frameCount = 0;
	size_t readCount = 0;
	int i;

	while (!feof(fSSR))
	{
		memset(&buff, 0, sizeof(buff));
		readCount = fread(&buff, sizeof(char), 1, fSSR);
		if (readCount < 1)
		{
			/* file error or eof of file */
			return -1;
		}
		currentCount += readCount;
		frameSize += readCount;

		int ret = input_spartn_data(spartn, spartn_out, buff);
		int  areaId = spartn_out->ssr_gad[0].areaId;
		double t1 = spartn_out->ssr[0].t0[0];
		double t2 = spartn_out->ssr[0].t0[1];
		double t3 = spartn_out->ssr[0].t0[2];
		double t4 = spartn_out->ssr[0].t0[3];
		double t5 = spartn_out->ssr[0].t0[4];
		double t6 = spartn_out->ssr[0].t0[5];

		if (ret == 1 && t1*t2*t3*t4*t5*t6 > 0.0 && spartn_out->type == 0 && spartn_out->eos == 1)
		{
			for (i = 0; i < SSR_NUM; i++)
			{
				if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys == 0)
					ssr_num[0]++;
				else if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys == 1)
					ssr_num[1]++;
			}

			frameSize = 0;
			frameCount++;
			break;
		}
	}
	return 0;
}

int read_ssr_from_file(FILE *fRTCM, gnss_rtcm_t *rtcm)
{
	int ret = -1;
	int iter_ssr = 0;
	size_t currentCount = 0;
	size_t readCount = 0;
	char buff = ' ';
	//int *numofread = 0;
	while (!feof(fRTCM))
	{
		memset(&buff, 0, sizeof(buff));
		readCount = fread(&buff, sizeof(char), 1, fRTCM);
		if (readCount < 1)
		{
			/* file error or eof of file */
			break;
		}
		currentCount += readCount;
		ret = input_rtcm3(buff, 0, rtcm);

		if (ret == 10 && rtcm->nav.ns > 30)
		{
			break;
		}
	}
	return ret;
}