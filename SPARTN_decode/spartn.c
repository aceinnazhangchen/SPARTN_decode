#include <stdio.h>
#include <memory.h>
#include "crc.h"
#include "bits.h"
#include "log.h"
#include "spartn.h"
#include "rtcm.h"
#include "GenVRSObs.h"
#include "ephemeris.h"
#include "gnss_math.h"
#include "tides.h"

#define Message_Type 2 //0:OCB 1:HPAC 2:GAD 3:LPAC
#define Leap_Sec 18.0
#define GLO_GPS_TD  10800
#define DAY_SECONDS 86400
#define SSR_SAP   

int read_obs_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int stnID);
int read_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int ssr_num);
int read_ssr_sapcorda(FILE *fSSR, raw_spartn_t *raw_spartn, spartn_t *spartn, int *ssr_num);
int read_ssr_from_file(FILE *fRTCM, gnss_rtcm_t *rtcm);

int decode_Dynamic_Key(raw_spartn_t* spartn) {
    return 1;
}

int decode_Group_Authentication(raw_spartn_t* spartn) {
    int i, bi, tab = 2;
    uint8_t* payload = spartn->buff + spartn->Payload_offset;
    int offset = 0;
    //Table 6.27 Message ID block 
    uint32_t SF089_Count_of_message_IDs = getbitu(payload, offset, 5);  offset += 5; log(LOG_DEBUG, tab, "SF089_Count_of_message_IDs = %d", SF089_Count_of_message_IDs);
    uint32_t SF090_Group_authentication_type = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF090_Group_authentication_type = %d", SF090_Group_authentication_type);
    uint32_t SF091_Computed_authentication_data_length = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF091_Computed_authentication_data_length = %d", SF091_Computed_authentication_data_length);
    //Table 6.28 
    uint32_t TF002_Message_type = getbitu(payload, offset, 7);  offset += 7; log(LOG_DEBUG, tab, "TF002_Message_type = %d", TF002_Message_type);
    uint32_t TF003_Message_sub_type = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "TF003_Message_sub_type = %d", TF003_Message_sub_type);
    uint32_t TF014_Encryption_sequence_number = getbitu(payload, offset, 6);  offset += 6; log(LOG_DEBUG, tab, "TF014_Encryption_sequence_number = %d", TF014_Encryption_sequence_number);
    //SF092 
    return 1;
}

int decode_spartn(raw_spartn_t* spartn) {
    switch (spartn->type)
    {
    case 0:
        return decode_OCB_message(spartn);
        break;
    case 1:
        return decode_HPAC_message(spartn);
        break;
    case 2:
        return decode_GAD_message(spartn);
        break;
    case 3:
        return decode_LPAC_message(spartn);
        break;
    case 4:
        //if (spartn->Subtype == 0) {
        //	decode_Dynamic_Key(spartn);
        //}
        //else if (spartn->Subtype == 1) {
        //	decode_Group_Authentication(spartn);
        //}
        break;
    }
    return 0;
}

int input_spartn_data(raw_spartn_t* spartn, spartn_t* spartn_out, uint8_t data) {
    int tab = 1;
    if (spartn->nbyte == 0) {
        if (data == SPARTN_PREAMB) {
            spartn->buff[spartn->nbyte++] = data;
        }
        return 0;
    }
    int CRC_Len = 0;
    if (spartn->nbyte > 4) {
        // type bit-len byte-len
        //	0		8		1
        //	1		16		2
        //	2		24		3
        //	3		32		4
        CRC_Len = spartn->CRC_type + 1;
    }
    int Time_tag_type_len = 0;
    if (spartn->nbyte > 6) {
        if (spartn->Time_tag_type == 0) {
            Time_tag_type_len = 2;
        }
        else if (spartn->Time_tag_type == 1) {
            Time_tag_type_len = 4;
        }
    }
    int EA_Len = 0;
    int EADL = 0;
    if (spartn->EAF == 1) {
        EA_Len = 2;
        if (spartn->AI > 1) {
            switch (spartn->EAL) {
            case 0:
                EADL = 8;//64
                break;
            case 1:
                EADL = 12;//96
                break;
            case 2:
                EADL = 16;//128
                break;
            case 3:
                EADL = 32;//256
                break;
            case 4:
                EADL = 64;//512
                break;
            }
        }
    }
    spartn->buff[spartn->nbyte++] = data;
    if (spartn->nbyte == 2) {//16
        spartn->type = getbitu(spartn->buff, 8, 7); log(LOG_DEBUG, tab, "type = %d", spartn->type);
    }
    else if (spartn->nbyte == 4) {//32
        spartn->len = getbitu(spartn->buff, 15, 10); log(LOG_DEBUG, tab, "len = %d BYTES", spartn->len);
        spartn->EAF = getbitu(spartn->buff, 25, 1); log(LOG_DEBUG, tab, "EAF = %d", spartn->EAF);
        spartn->CRC_type = getbitu(spartn->buff, 26, 2); log(LOG_DEBUG, tab, "CRC_type = %d", spartn->CRC_type);
        spartn->Frame_CRC = getbitu(spartn->buff, 28, 4);
        char Frame_CRC_Buffer[3] = { 0 };
        bitscopy(Frame_CRC_Buffer, 0, spartn->buff + 1, 0, 20);
        uint8_t Result_Frame_CRC = crc4_itu(Frame_CRC_Buffer, 3);
        log(LOG_DEBUG, tab, "Frame_CRC = %d : %d", spartn->Frame_CRC, Result_Frame_CRC);
        if (spartn->Frame_CRC != Result_Frame_CRC) {
            memset(spartn, 0, sizeof(raw_spartn_t));
            return -1;
        }
    }
    else if (spartn->nbyte == 5) {//40
        spartn->Subtype = getbitu(spartn->buff, 32, 4); log(LOG_DEBUG, tab, "Subtype = %d", spartn->Subtype);
        spartn->Time_tag_type = getbitu(spartn->buff, 36, 1); log(LOG_DEBUG, tab, "Time_tag_type = %d", spartn->Time_tag_type);
    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + Time_tag_type_len) {//48
        spartn->GNSS_time_type = getbitu(spartn->buff, 37, Time_tag_type_len * 8); log(LOG_DEBUG, tab, "GNSS_time_type = %d", spartn->GNSS_time_type);
        spartn->Solution_ID = getbitu(spartn->buff, 37 + Time_tag_type_len * 8, 7); log(LOG_DEBUG, tab, "Solution_ID = %d", spartn->Solution_ID);
        spartn->Solution_processor_ID = getbitu(spartn->buff, 44 + Time_tag_type_len * 8, 4); log(LOG_DEBUG, 1, "Solution_processor_ID = %d", spartn->Solution_processor_ID);
    }
    else if (spartn->EAF == 1 && Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len) {//64
        spartn->Encryption_ID = getbitu(spartn->buff, 48 + Time_tag_type_len * 8, 4); log(LOG_DEBUG, tab, "Encryption_ID = %d", spartn->Encryption_ID);
        spartn->ESN = getbitu(spartn->buff, 52 + Time_tag_type_len * 8, 6); log(LOG_DEBUG, tab, "ESN = %d", spartn->ESN);
        spartn->AI = getbitu(spartn->buff, 55 + Time_tag_type_len * 8, 3); log(LOG_DEBUG, tab, "AI = %d", spartn->AI);
        spartn->EAL = getbitu(spartn->buff, 58 + Time_tag_type_len * 8, 3); log(LOG_DEBUG, tab, "EAL = %d", spartn->EAL);
    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len) {
        spartn->Payload_offset = 6 + EA_Len + Time_tag_type_len;
        log(LOG_DEBUG, tab, "Payload_offset = %d", spartn->Payload_offset);
    }
    else if (spartn->EAF == 1 && spartn->AI > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL) {

    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL + CRC_Len) {
        spartn->Message_CRC = getbitu(spartn->buff, (spartn->nbyte - CRC_Len) * 8, CRC_Len * 8);
        log(LOG_DEBUG, tab, "nbyte = %d", spartn->nbyte);

        uint32_t Result_Message_CRC = crc24_radix(spartn->buff + 1, 5 + EA_Len + Time_tag_type_len + spartn->len + EADL);
        log(LOG_DEBUG, tab, "Message_CRC = %d : %d", spartn->Message_CRC, Result_Message_CRC);
        if (spartn->Message_CRC == Result_Message_CRC) {
           expanded_full_time(spartn);
            log(LOG_DEBUG, tab, "==========");
            spartn->spartn_out = spartn_out;
            spartn_out->type = spartn->type;
            spartn_out->Subtype = spartn->Subtype;
            spartn_out->len = spartn->len;
            decode_spartn(spartn);
            log(LOG_DEBUG, tab, "==========");
        }
        memset(spartn, 0, sizeof(raw_spartn_t));
        return 1;
    }
    return 0;
}

int check_nav(nav_t *nav, sap_ssr_t *sap_ssr)
{
    double time0 = sap_ssr[0].t0[1];
    int i,eph_n=0, geph_n = 0;
    for (i = 0; i < nav->n_gps; i++)
    {
        double dt = fabs(time0 - nav->eph[i].toe.time);
        if (dt < 7200.0) eph_n++;
    }

    for (i = 0; i < nav->ng; i++)
    {
        double dt = fabs(time0 - nav->geph[i].toe.time);
        if (dt < 1800.0) geph_n++;
    }

    if (eph_n == 0 || geph_n == 0)  return 1;
    return 0;
}

void print_dd_obs(obs_t *obs_ref, obs_t *obs_rov, vec_t *vec_ref, vec_t *vec_rov, int *iref, int *irov, int nsd, FILE *fLOG)
{
    obsd_t* pObsRef   = NULL;
    obsd_t* pObsRov   = NULL;
    obsd_t *pObsRef_0 = NULL;
    obsd_t *pObsRov_0 = NULL;
    int isd, f, s, sys, i, j, sat1, prn1, sat0, prn0;
    double sdP[2] = { 0.0 }, sdL[2] = { 0.0 };
    int nm = 0;
    int wn = 0;
    int num_valid_frq = 0;
    int num_sat = 0, refloc;
    double ws = time2gpst(obs_rov->time, &wn);
    int refsat[NSYS*NFREQ] = { 0 };
    for (s = 0; s < NSYS; ++s)
    {
        refloc = -1;
        num_sat = 0;
        for (isd = 0; isd < nsd; ++isd)
        {
            j = irov[isd];
            i = iref[isd];
            pObsRef_0 = obs_ref->data + i;
            pObsRov_0 = obs_rov->data + j;
            sat0 = pObsRov_0->sat;
            if (satidx(sat0, &prn0) != s)                           continue;
            if (pObsRov_0->L[0] == 0.0 || pObsRef_0->L[0] == 0.0)   continue;
            if (pObsRov_0->P[0] == 0.0 || pObsRef_0->P[0] == 0.0)   continue;
            if (pObsRov_0->L[1] == 0.0 || pObsRef_0->L[1] == 0.0)   continue;
            if (pObsRov_0->P[1] == 0.0 || pObsRef_0->P[1] == 0.0)   continue;
            num_sat++;
            refloc = isd;
            refsat[MI(s, 0, NFREQ)] = sat0;
            refsat[MI(s, 1, NFREQ)] = sat0;
            break;
        }
        if (refloc < 0) continue; /* cannot find the reference satellite */

        for (isd = 0; isd < nsd; ++isd)
        {
            if (isd == refloc) continue;
            j = irov[isd];
            i = iref[isd];
            pObsRef = obs_ref->data + i;
            pObsRov = obs_rov->data + j;
            sat1 = pObsRov->sat;
            sys = satsys(sat1, &prn1);
            if (satidx(sat1, &prn1) != s) continue;

            for (f = 0; f < NFREQ; ++f)
            {
                if (fabs(pObsRov->P[f]) < 0.001 || fabs(pObsRef->P[f]) < 0.001)   sdP[f] = 0.0;
                else                                                              sdP[f] = (pObsRov->P[f] - pObsRef->P[f]) - (pObsRov_0->P[f] - pObsRef_0->P[f]);

                if (fabs(pObsRov->L[f]) < 0.001 || fabs(pObsRef->L[f]) < 0.001)   sdL[f] = 0.0;
                else                                                              sdL[f] = (pObsRov->L[f] - pObsRef->L[f]) - (pObsRov_0->L[f] - pObsRef_0->L[f]);
                ++nm;
            }
 /*           fprintf(fLOG, "DD,%4i,%10.4f,%3i,%3i,%3i,%10.4f,%10.4f,%10.4f,%10.4f\n", wn, ws, sys, prn0, prn1, sdP[0], sdL[0], sdP[1], sdL[1]);*/
            printf("DD,%4i,%10.4f,%3i,%3i,%3i,%10.4f,%10.4f,%10.4f,%10.4f\n", wn, ws, sys, prn0, prn1, sdP[0], sdL[0], sdP[1], sdL[1]);
        }
    }
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

int gga_ssr2osr_main(FILE *fSSR, FILE *fEPH, double *ep, double *rovpos)  // int year, int doy, 
{
    gnss_rtcm_t rtcm      = { 0 };
    nav_t *nav            = &rtcm.nav;
    obs_t *rov            = rtcm.obs;
    obs_t obs_vrs         = { 0.0 };
    vec_t vec_vrs[MAXOBS] = { 0.0 };
    gtime_t time0         = epoch2time(ep);
    double cur_time       = (int)time0.time;
    int doy=time2doy(time0);
    int year = ep[0];
    set_approximate_time(year, doy, rtcm.rcv);
    if (fSSR == NULL)  return 0;
    if (fEPH == NULL)  return 0;
#ifdef TABLE_LOG
    open_ocb_table_file (NULL);
    open_hpac_table_file(NULL);
    open_gad_table_file (NULL);
    open_lpac_table_file(NULL);
#endif
    raw_spartn_t spartn;
    memset(&spartn, 0, sizeof(spartn));
    spartn_t spartn_out;
    memset(&spartn_out, 0, sizeof(spartn_t));
    sap_ssr_t *sap_ssr = &spartn_out.ssr;
    gad_ssr_t *sap_gad = &spartn_out.ssr_gad;
    OCB_t ocb = { 0 };
    HPAC_t hpac = { 0 };
    GAD_t gad = { 0 };
    LPAC_t lpac = { 0 };
    spartn_out.ocb  = &ocb;
    spartn_out.hpac = &hpac;
    spartn_out.gad  = &gad;
    spartn_out.lpac = &lpac;
    int i, j, nsat;
    int rov_ret, ret_nav, num_ssr = -1;
    double blh[3] = { 0.0 }, dr[3] = { 0.0 };
    gtime_t teph =epoch2time(ep);
    double obs_time = 0.0;
    int nc = 0;
    while (1)
    {
        nav->ns      = 0;
        nav->nsys[0] = 0;
        nav->nsys[1] = 0;
        read_ssr_sapcorda(fSSR, &spartn, &spartn_out, nav->nsys);
        if (feof(fSSR)) break;

        /* read broadcast eph data one byte */
        ret_nav = read_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
        if (ret_nav != 2)
        {
            /* can not find the complete epoch data */
            if (feof(fEPH)) break;
        }

        for (i = 0; i < spartn_out.ssr_offset; i++)
        {
            if (sap_ssr[i].t0[0] > 0.0) nav->ns++;
        }
        for (i = 0; i < nav->ns; i++)
        {
            int nav_iod = -1;
            int sys = sap_ssr[i].sys;
            if (sys == 0)
            {
                for (j = 0; j < nav->n; j++)
                {
                    if (sap_ssr[i].prn == nav->eph[j].sat)
                    {
                        nav_iod = nav->eph[j].iode;
                        break;
                    }
                }
            }
            else if (sys == 1)
            {
                for (j = 0; j < nav->ng; j++)
                {
                    if (sap_ssr[i].prn + 40 == nav->geph[j].sat)
                    {
                        nav_iod = nav->geph[j].iode;
                        break;
                    }
                }
            }
            if (nav_iod != sap_ssr[i].iod[0]) continue;
            printf("ocb:%6.0f,%6.0f,%3i,%3i,%2i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n", sap_ssr[i].t0[0], sap_ssr[i].t0[1], nav_iod, sap_ssr[i].iod[0], sys, sap_ssr[i].prn, sap_ssr[i].deph[0], sap_ssr[i].deph[1], sap_ssr[i].deph[2], sap_ssr[i].dclk, sap_ssr[i].cbias[0], sap_ssr[i].cbias[1], sap_ssr[i].cbias[2], sap_ssr[i].pbias[0], sap_ssr[i].pbias[1], sap_ssr[i].pbias[2]);
        }

        while (1)
        { 
            teph = timeadd(teph, 5.0);
            int time = teph.time;
            double time1 = fmod((double)time, 86400);
            if (fabs(time1 - sap_ssr[0].t0[1]) < 5.0)
                break;
        }
        
        nsat=satposs_sap_rcv(teph, rovpos, vec_vrs, nav, sap_ssr, EPHOPT_SSRSAP);

        obs_vrs.time = teph;
        obs_vrs.n    = nsat;
        memcpy(obs_vrs.pos, rovpos,3*sizeof(double));
        for (i=0;i<nsat;i++)
        {
            obs_vrs.data[i].sat = vec_vrs[i].sat;
        }
        nsat = compute_vector_data(&obs_vrs, vec_vrs);
        
        int vrs_ret = gen_obs_from_ssr(teph, rovpos, sap_ssr, sap_gad, &obs_vrs, vec_vrs, 0.0);

        nc++;
    }
}

int obs_ssr2osr_main(FILE *fSSR, FILE *fEPH, FILE *fROV, int year, int doy, double *ep, double *rovpos)
{
    gnss_rtcm_t rtcm = { 0 };
    nav_t *nav = &rtcm.nav;
    obs_t *rov = rtcm.obs;
    vec_t vec_vrs[MAXOBS] = { 0.0 };
    obs_t obs_vrs = { 0.0 };
    gtime_t time0 = epoch2time(ep);
    double cur_time = (int)time0.time;
    double cur_time0 = floor(cur_time / 86400.0) * 86400.0;
    if (fROV == NULL)  return 0;
    if (fSSR == NULL)  return 0;
    if (fEPH == NULL)  return 0;
#ifdef TABLE_LOG
    open_ocb_table_file(NULL);
    open_hpac_table_file(NULL);
    open_gad_table_file(NULL);
    open_lpac_table_file(NULL);
#endif

#ifdef SSR_SAP 
    raw_spartn_t spartn;
    memset(&spartn, 0, sizeof(spartn));
    spartn_t spartn_out;
    memset(&spartn_out, 0, sizeof(spartn_t));
    sap_ssr_t *sap_ssr = &spartn_out.ssr;
    gad_ssr_t *sap_gad = &spartn_out.ssr_gad;
    OCB_t  ocb  = { 0 };
    HPAC_t hpac = { 0 };
    GAD_t  gad  = { 0 };
    LPAC_t lpac = { 0 };
    spartn_out.ocb  = &ocb;
    spartn_out.hpac = &hpac;
    spartn_out.gad  = &gad;
    spartn_out.lpac = &lpac;
#endif
    set_approximate_time(year, doy, rtcm.rcv);
    int i, j, rov_ret, ret_nav, num_ssr = -1;
    double blh[3] = { 0.0 }, dr[3] = { 0.0 };
    while (1)
    {
        nav->ns = 0;
        nav->nsys[0] = 0;  nav->nsys[1] = 0;
#ifdef SSR_SAP 
        read_ssr_sapcorda(fSSR, &spartn, &spartn_out, nav->nsys);
#else
        num_ssr = read_ssr_from_file(fSSR, &rtcm);
#endif
        if (feof(fSSR)) break;
        /* read broadcast eph data one byte */
        ret_nav = read_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
        if (ret_nav != 2)
        {
            /* can not find the complete epoch data */
            if (feof(fEPH)) break;
        }
        
#ifdef SSR_SAP 
        for (i = 0; i < spartn_out.ssr_offset; i++)
        {
            if (sap_ssr[i].t0[0] > 0.0) nav->ns++;
        }
        for (i = 0; i < nav->ns; i++)
        {
            int nav_iod = -1;
            double eph_toe = 0.0;
            int sys = sap_ssr[i].sys;
            if (sys == 0)
            {
                for (j = 0; j < nav->n; j++)
                {
                    if (sap_ssr[i].prn == nav->eph[j].sat)
                    {
                        nav_iod = nav->eph[j].iode;
                        break;
                    }
                }
            }
            else if (sys == 1)
            {
                for (j = 0; j < nav->ng; j++)
                {
                    if (sap_ssr[i].prn + 40 == nav->geph[j].sat)
                    {
                        nav_iod = nav->geph[j].iode;
                        break;
                    }
                }
            }
            for (j = 0; j < 6; j++)
            {
                if (sap_ssr[i].t0[j] == 0)            continue;
                if (sap_ssr[i].t0[j] < DAY_SECONDS)   sap_ssr[i].t0[j] = sap_ssr[i].t0[j] + cur_time0;
            }
            if (nav_iod != sap_ssr[i].iod[0]) continue;
            printf("ocb:%7.0f,%7.0f,%2i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n", sap_ssr[i].t0[0]- cur_time0, sap_ssr[i].t0[1] - cur_time0, sys, sap_ssr[i].prn, sap_ssr[i].deph[0], sap_ssr[i].deph[1], sap_ssr[i].deph[2], sap_ssr[i].dclk, sap_ssr[i].cbias[0], sap_ssr[i].cbias[1], sap_ssr[i].cbias[2], sap_ssr[i].pbias[0], sap_ssr[i].pbias[1], sap_ssr[i].pbias[2]);
        }

        int idx = -1;
        while (1)
        {
            if (rov->time.time - sap_ssr[0].t0[1] < -0.1)
                rov_ret = read_obs_rtcm(fROV, &rtcm, 0);
            else
            {
                if (fabs(rov->time.time - sap_ssr[0].t0[1]) < 0.02)  idx = 1;
                break;
            }
        }
        if (idx == -1) continue;
#else
        int idx = -1;
        while (1)
        {
            if (rov->time.time - nav->ssr[0].t0[0].time < 0.1)
                rov_ret = read_obs_rtcm(fROV, &rtcm, 0);
            else
            {
                if (fabs(rov->time.time - nav->ssr[0].t0[0].time) < 5.02)  idx = 1;
                break;
            }
        }
        if (idx == -1) continue;
#endif
        memcpy(rov->pos, rovpos, 3 * sizeof(double));
        ecef2pos(rov->pos, blh);
        int time = rov->time.time;
        //printf("obs:%12.0f,%6.3f,%6.3f", fmod((double)time,86400), blh[0] * R2D, blh[1] * R2D);
        //for (i = 0; i < rov->n; i++)
        //{
        //    if (rov->data[i].sat > 70) continue;
        //    printf(",%3i", rov->data[i].sat);
        //}
        //printf("\n");
        memset(vec_vrs, 0, sizeof(vec_vrs));
        memset(&obs_vrs, 0, sizeof(obs_t));
        memcpy(obs_vrs.pos, rov->pos, 3 * sizeof(double));
#ifdef SSR_SAP 
        satposs_sap(rov, vec_vrs, nav, sap_ssr, EPHOPT_SSRSAP);
#else
        satposs(rov, vec_vrs, nav, EPHOPT_SSRAPC);
#endif
        compute_vector_data(rov, vec_vrs);
#ifdef SSR_SAP 
        int vrs_ret = gen_vobs_from_ssr(rov, sap_ssr, sap_gad, &obs_vrs, vec_vrs, 0.0);
#endif
        //int iref[MAXOBS] = { 0 }, irov[MAXOBS] = { 0 };
        //int nsd = get_match_epoch(rov, &obs_vrs, vec_vrs, vec_vrs, iref, irov);
        //if (nsd == 0) continue;
        //print_dd_obs(rov, &obs_vrs, vec_vrs, vec_vrs, iref, irov, nsd, NULL);
    }
#ifdef TABLE_LOG
    close_ocb_table_file();
    close_hpac_table_file();
    close_gad_table_file();
    close_lpac_table_file();
#endif
    return 1;
}


int main() 
{
    FILE *fSSR = { NULL };
    FILE *fEPH = { NULL };
    FILE *fROV = { NULL };
    //fROV = fopen("..\\20200317\\SF0320077g.dat", "rb");
    //fSSR = fopen("..\\20200317\\SPARTN20200317062037.raw","rb");
    //fEPH = fopen("..\\20200317\\Aux20200317062038.raw", "rb");
    //int year = 2020;
    //int doy = 77;
    //double ep[6] = { 2020,3,17,0,0,0 };

    //fROV = fopen("..\\20200316\\sf03076j.rtcm", "rb");
    //fSSR = fopen("..\\20200316\\SPARTN20200316091017.raw", "rb");
    //fEPH = fopen("..\\20200316\\Aux20200316091039.raw", "rb");
    //int year = 2020;
    //int doy = 76;
    //double ep[6] = { 2020,3,16,0,0,0 };

    //fROV = fopen("..\\20200405\\SF0320096h.dat", "rb");
    //fSSR = fopen("..\\20200405\\SSR20096h.dat", "rb");
    //fEPH = fopen("..\\20200405\\EPH20096h.dat", "rb");
    //int year = 2020;
    //int doy = 96;
    //double ep[6] = { 2020,4,5,0,0,0 };

    //fROV = fopen("..\\20200406\\SF0320097m.dat", "rb");
    //fSSR = fopen("..\\20200406\\SSR20097m.dat", "rb");
    //fEPH = fopen("..\\20200406\\EPH20097m.dat", "rb");
    //int year = 2020;
    //int doy = 97;
    //double ep[6] = { 2020,4,6,0,0,0 };

    //fROV = fopen("..\\20200407\\SF0320098l.dat", "rb");
    //fSSR = fopen("..\\20200407\\SSR20098l.dat", "rb");
    //fEPH = fopen("..\\20200407\\EPH20098l.dat", "rb");
    //int year = 2020;
    //int doy  = 98;
    //double ep[6] = { 2020,4,7,0,0,0 };

    fROV = fopen("..\\20200420\\SF0320111h.dat", "rb");
    fSSR = fopen("..\\20200420\\SPARTN20200420070130.raw", "rb");
    fEPH = fopen("..\\20200420\\Aux20200420070149.raw", "rb");
    int year = 2020;
    int doy = 111;
    //double ep[6] = { 2020,4,20,7,2,20 };
    double ep[6] = { 2020,4,20,7,1,50 };
    //07:47 : 00

    double rovpos[3] = { -2705297.408,-4283455.631,3861823.955 };

    //obs_ssr2osr_main(fSSR, fEPH, fROV, year, doy,ep, rovpos);

    gga_ssr2osr_main(fSSR, fEPH, ep, rovpos);

	return 0;
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


int read_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int ns_gps, int ns_g)
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
        if (ret == 2 && rtcm->nav.n_gps>= ns_gps && rtcm->nav.ng >= ns_g)
        {
            break;
        }
    }
    return ret;
}

int read_ssr_sapcorda(FILE *fSSR, raw_spartn_t *spartn, spartn_t *spartn_out, unsigned int *ssr_num)
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
        frameSize    += readCount;

        int ret = input_spartn_data(spartn, spartn_out, buff);
        int  areaId = spartn_out->ssr_gad[0].areaId;
        double t1 = spartn_out->ssr[0].t0[0];
        double t2 = spartn_out->ssr[0].t0[1];
        double t3 = spartn_out->ssr[0].t0[2];
        double t4 = spartn_out->ssr[0].t0[3];
        double t5 = spartn_out->ssr[0].t0[4];
        double t6 = spartn_out->ssr[0].t0[5];

        if (ret== 1 && t1*t2*t3*t4*t5*t6>0.0 && spartn_out->type ==0 && spartn_out->eos==1)
        {
            //printf("\n");
            for (i = 0; i < SSR_NUM; i++)
            {
                if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys==0)
                    ssr_num[0]++;
                else if (spartn_out->ssr[i].prn != 0 && spartn_out->ssr[i].sys == 1)
                    ssr_num[1]++;
            }
            frameSize = 0;
            frameCount++;
            break;
        }
    }
    return;
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