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
extern double gravitationalDelayCorrection(const int sys, const double *receiverPosition,
    const double *satellitePosition)
{
    double	receiverModule;
    double	satelliteModule;
    double	distance;
    double  MU = MU_GPS;

    receiverModule = sqrt(receiverPosition[0] * receiverPosition[0] + receiverPosition[1] * receiverPosition[1] +
        receiverPosition[2] * receiverPosition[2]);
    satelliteModule = sqrt(satellitePosition[0] * satellitePosition[0] + satellitePosition[1] * satellitePosition[1] +
        satellitePosition[2] * satellitePosition[2]);
    distance = sqrt((satellitePosition[0] - receiverPosition[0])*(satellitePosition[0] - receiverPosition[0]) +
        (satellitePosition[1] - receiverPosition[1])*(satellitePosition[1] - receiverPosition[1]) +
        (satellitePosition[2] - receiverPosition[2])*(satellitePosition[2] - receiverPosition[2]));

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

void find_nearest_gridpoints_ionocoef(double *blh, int sat, gad_ssr_t *gad, double *gpt_idx)
{
    int i, j, k, l, lon_nc, lat_nc, lon_sp, lat_sp, dlon, dlat;
    double blh_ref[3] = { 0 }, ned[3] = { 0 }, dist2D = 0.0, MinDist2D=1.0e8;

    gpt_idx[0] = -1;
    gpt_idx[3] = -1;
    gpt_idx[4] = -1;
    for (j = 0; j < RAP_NUM; j++)
    {
        blh_ref[1] = blh[2];
        lat_nc = gad[j].nc_lat;
        lon_nc = gad[j].nc_lon;
        for (k = 0; k < lat_nc; k++)
        {
            for (l = 0; l < lon_nc; l++)
            {
                blh_ref[0] = gad[j].rap_lat * D2R + k * gad[j].spa_lat;
                blh_ref[1] = gad[j].rap_lon * D2R + l * gad[j].spa_lon;
                blhdiff(blh, blh_ref, ned);
                dist2D = sqrt(ned[0] * ned[0] + ned[1] * ned[1]);
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
    }
}

void dist_inv_unit_weighting(double *blh, double *gpt_bl, double *wdi)
{
    int i;
    double blh_ref[3] = { 0.0 }, ned[3] = { 0.0 }, dist2D[4] = { 0.0 }, sumdis2D = 0.0;
    for (i = 0; i < 4; i++)
    {
        blh_ref[0] = gpt_bl[i * 2 + 0];
        blh_ref[2] = gpt_bl[i * 2 + 1];
        blh_ref[2] = blh[2];

        blhdiff(blh, blh_ref, ned);
        dist2D[i] = sqrt(ned[0] * ned[0] + ned[1] * ned[1]);
        sumdis2D += SQR(dist2D[i]);
    }

    for (i = 0; i < 4; i++) wdi[i] = SQR(dist2D[i]) / sumdis2D;
}


void high_prcision_slant_atm_polynomial(gtime_t time, double *blh, int sat, double *azel, sap_ssr_t *ssr, gad_ssr_t *gad, int *gpt_idx, double *stec, double *stro)
{
    int i, j, satidx=-1;
    double icoef[12]  = { 0.0 };
    double tcoef[12] = { 0.0 };
    double gpt_bl[8] = { 0.0 };
    int    gpt_pos[8] = { 0 };
    double lon_sp = 0.0, lat_sp = 0.0, rap_lon=0.0, rap_lat=0.0;
    double Ip[4] = { 0.0 }, Tp[4] = { 0.0 };
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
    lat_sp  = gad[gpt_idx[0]].spa_lat;
    lon_sp  = gad[gpt_idx[0]].spa_lon;
    rap_lat = gad[gpt_idx[0]].rap_lat;
    rap_lon = gad[gpt_idx[0]].rap_lon;

    if (gpt_idx[3] * gpt_idx[4] > 0)
    {
        if (gpt_idx[3]<0) /*rap-> (1,1)    */
        {
            gpt_pos[0] = gpt_idx[1]-1;       gpt_pos[1] = gpt_idx[2]-1;
            gpt_pos[2] = gpt_idx[1]-1;       gpt_pos[3] = gpt_idx[2];
            gpt_pos[4] = gpt_idx[1];         gpt_pos[5] = gpt_idx[2]-1;
            gpt_pos[6] = gpt_idx[1];         gpt_pos[7] = gpt_idx[2];
        }
        else /*rap-> (0,0)    */
        {
            gpt_pos[0] = gpt_idx[1];         gpt_pos[1] = gpt_idx[2];
            gpt_pos[2] = gpt_idx[1];         gpt_pos[3] = gpt_idx[2]+1;
            gpt_pos[4] = gpt_idx[1]+1;       gpt_pos[5] = gpt_idx[2];
            gpt_pos[6] = gpt_idx[1]+1;       gpt_pos[7] = gpt_idx[2]+1;
        }
    }
    else 
    { 
        if (gpt_idx[3] < 0)  /*rap-> (0,1) */
        {
            gpt_pos[0] = gpt_idx[1];          gpt_pos[1] = gpt_idx[2] - 1;
            gpt_pos[2] = gpt_idx[1];          gpt_pos[3] = gpt_idx[2];
            gpt_pos[4] = gpt_idx[1] + 1;      gpt_pos[5] = gpt_idx[2] - 1;
            gpt_pos[6] = gpt_idx[1] + 1;      gpt_pos[7] = gpt_idx[2];
        }
        else   /*rap-> (1,0)    */
        {
            gpt_pos[0] = gpt_idx[1] - 1;      gpt_pos[1] = gpt_idx[2];
            gpt_pos[2] = gpt_idx[1] - 1;      gpt_pos[3] = gpt_idx[2] + 1;
            gpt_pos[4] = gpt_idx[1];          gpt_pos[5] = gpt_idx[2];
            gpt_pos[6] = gpt_idx[1];          gpt_pos[7] = gpt_idx[2] + 1;
        }
    }

    gpt_bl[0] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[0];
    gpt_bl[1] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[1];
    gpt_bl[2] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[2];
    gpt_bl[3] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[3];
    gpt_bl[4] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[4];
    gpt_bl[5] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[5];
    gpt_bl[6] = gad[gpt_idx[0]].rap_lat + lat_sp * gpt_pos[6];
    gpt_bl[7] = gad[gpt_idx[0]].rap_lon + lon_sp * gpt_pos[7];

    dist_inv_unit_weighting(blh, gpt_bl, wdi);

    *stec = 0.0;
    *stro = 0.0;
    for (i = 0; i < 4; i++)
    {
        icoef[i * 3]     = ssr[satidx].stec_coef[0 + gpt_idx[0] * 3];
        icoef[i * 3 + 1] = ssr[satidx].stec_coef[1 + gpt_idx[0] * 3];
        icoef[i * 3 + 2] = ssr[satidx].stec_coef[2 + gpt_idx[0] * 3];
        tcoef[i * 3]     = ssr[satidx].tro_coef [0 + gpt_idx[0] * 3];
        tcoef[i * 3 + 1] = ssr[satidx].tro_coef [1 + gpt_idx[0] * 3];
        tcoef[i * 3 + 2] = ssr[satidx].tro_coef [2 + gpt_idx[0] * 3];

        Ip[i] = icoef[i * 3] + icoef[i * 3 + 1] * (gpt_bl[i * 2]- rap_lat) + icoef[i * 3 + 2] * (gpt_bl[i * 2 + 1]-rap_lon);
        Tp[i] = tcoef[i * 3] + tcoef[i * 3 + 1] * (gpt_bl[i * 2]- rap_lat) + tcoef[i * 3 + 2] * (gpt_bl[i * 2 + 1]-rap_lon);

        *stec += wdi[i] * Ip[i];
        Tw    += wdi[i] * Tp[i];
    }

    Th = ssr[satidx].ave_htd;
    m_h=tropmapf(time, blh, azel, &m_w);
    *stro = m_h*Th + m_w * Tw;
}

void compute_high_prcision_atm_corr(int sat, obs_t *obs_vrs, sap_ssr_t *ssr, gad_ssr_t *gad, double *azel, double maskElev, double *stec, double *stro)
{
    double blh[3] = { 0 };
    int gpt_idx[5] = { 0 };
    ecef2pos(obs_vrs, blh);
    find_nearest_gridpoints_ionocoef(blh, sat, gad, gpt_idx);
    high_prcision_slant_atm_polynomial(obs_vrs->time, blh, sat, ssr, gad, azel, gpt_idx, &stec, &stro);
}


extern int gen_vobs_from_ssr(obs_t *obs_rov, sap_ssr_t *ssr, gad_ssr_t* gad, obs_t *obs_vrs, vec_t *vec_vrs, double maskElev)
{
    int i;
    double P[2] = { 0.0 }, L[2] = { 0.0 };
    double rs[6] = { 0.0 }, rr[6] = { 0.0 }, phw = 0.0;
    double stec = 0.0, strop = 0.0;
    double dr[3] = { 0.0 };

    obs_vrs->n    = obs_rov->n;
    obs_vrs->time = obs_rov->time;
    memcpy(obs_vrs->pos, obs_rov->pos, sizeof(double) * 3);

    tidedisp(obs_rov->time, &obs_vrs->pos, 1, NULL, NULL, dr);

    for (i = 0; i < obs_rov->n; i++)
    {
       compute_high_prcision_atm_corr(obs_rov->data[i].sat, obs_vrs, ssr, gad, &vec_vrs[i].azel, maskElev, &stec, &strop);

       model_phw(obs_rov->time, obs_rov->data[i].sat,NULL,2, rs, rr, &phw);

      
    }
    return 1;
}

