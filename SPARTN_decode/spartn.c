#include <stdio.h>
#include <memory.h>
#include "crc.h"
#include "bits.h"
#include "log.h"
#include "spartn.h"
#include "rtcm.h"
#include "GenVRSObs.h"
#include "gnss_math.h"
#include "tides.h"

#define Message_Type 2 //0:OCB 1:HPAC 2:GAD 3:LPAC
#define Leap_Sec 18.0
#define GLO_GPS_TD  10800

int decode_Dynamic_Key(raw_spartn_t* spartn) {
    return 1;
}

int decode_Group_Authentication(raw_spartn_t* spartn) {
    int tab = 2;
    uint8_t* payload = spartn->buff + spartn->Payload_offset;
    int offset = 0;
    //Table 6.27 Message ID block 
    uint32_t SF089_Count_of_message_IDs = getbitu(payload, offset, 5);  offset += 5; slog(LOG_DEBUG, tab, "SF089_Count_of_message_IDs = %d", SF089_Count_of_message_IDs);
    uint32_t SF090_Group_authentication_type = getbitu(payload, offset, 4);  offset += 4; slog(LOG_DEBUG, tab, "SF090_Group_authentication_type = %d", SF090_Group_authentication_type);
    uint32_t SF091_Computed_authentication_data_length = getbitu(payload, offset, 4);  offset += 4; slog(LOG_DEBUG, tab, "SF091_Computed_authentication_data_length = %d", SF091_Computed_authentication_data_length);
    //Table 6.28 
    uint32_t TF002_Message_type = getbitu(payload, offset, 7);  offset += 7; slog(LOG_DEBUG, tab, "TF002_Message_type = %d", TF002_Message_type);
    uint32_t TF003_Message_sub_type = getbitu(payload, offset, 4);  offset += 4; slog(LOG_DEBUG, tab, "TF003_Message_sub_type = %d", TF003_Message_sub_type);
    uint32_t TF014_Encryption_sequence_number = getbitu(payload, offset, 6);  offset += 6; slog(LOG_DEBUG, tab, "TF014_Encryption_sequence_number = %d", TF014_Encryption_sequence_number);
    //SF092 
    return 1;
}

int decode_spartn(raw_spartn_t* spartn, spartn_t* spartn_out) {
    switch (spartn->type)
    {
    case 0:
        return decode_OCB_message(spartn, spartn_out);
        break;
    case 1:
        return decode_HPAC_message(spartn, spartn_out);
        break;
    case 2:
        return decode_GAD_message(spartn, spartn_out);
        break;
    case 3:
        return decode_LPAC_message(spartn, spartn_out);
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
        spartn->type = getbitu(spartn->buff, 8, 7); slog(LOG_DEBUG, tab, "type = %d", spartn->type);
    }
    else if (spartn->nbyte == 4) {//32
        spartn->len = getbitu(spartn->buff, 15, 10); slog(LOG_DEBUG, tab, "len = %d BYTES", spartn->len);
        spartn->EAF = getbitu(spartn->buff, 25, 1); slog(LOG_DEBUG, tab, "EAF = %d", spartn->EAF);
        spartn->CRC_type = getbitu(spartn->buff, 26, 2); slog(LOG_DEBUG, tab, "CRC_type = %d", spartn->CRC_type);
        spartn->Frame_CRC = getbitu(spartn->buff, 28, 4);
        char Frame_CRC_Buffer[3] = { 0 };
        bitscopy(Frame_CRC_Buffer, 0, spartn->buff + 1, 0, 20);
        uint8_t Result_Frame_CRC = crc4_itu(Frame_CRC_Buffer, 3);
        slog(LOG_DEBUG, tab, "Frame_CRC = %d : %d", spartn->Frame_CRC, Result_Frame_CRC);
        if (spartn->Frame_CRC != Result_Frame_CRC) {
            memset(spartn, 0, sizeof(raw_spartn_t));
            return -1;
        }
    }
    else if (spartn->nbyte == 5) {//40
        spartn->Subtype = getbitu(spartn->buff, 32, 4); slog(LOG_DEBUG, tab, "Subtype = %d", spartn->Subtype);
        spartn->Time_tag_type = getbitu(spartn->buff, 36, 1); slog(LOG_DEBUG, tab, "Time_tag_type = %d", spartn->Time_tag_type);
    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + Time_tag_type_len) {//48
        spartn->GNSS_time_type = getbitu(spartn->buff, 37, Time_tag_type_len * 8); slog(LOG_DEBUG, tab, "GNSS_time_type = %d", spartn->GNSS_time_type);
        spartn->Solution_ID = getbitu(spartn->buff, 37 + Time_tag_type_len * 8, 7); slog(LOG_DEBUG, tab, "Solution_ID = %d", spartn->Solution_ID);
        spartn->Solution_processor_ID = getbitu(spartn->buff, 44 + Time_tag_type_len * 8, 4); slog(LOG_DEBUG, 1, "Solution_processor_ID = %d", spartn->Solution_processor_ID);
    }
    else if (spartn->EAF == 1 && Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len) {//64
        spartn->Encryption_ID = getbitu(spartn->buff, 48 + Time_tag_type_len * 8, 4); slog(LOG_DEBUG, tab, "Encryption_ID = %d", spartn->Encryption_ID);
        spartn->ESN = getbitu(spartn->buff, 52 + Time_tag_type_len * 8, 6); slog(LOG_DEBUG, tab, "ESN = %d", spartn->ESN);
        spartn->AI = getbitu(spartn->buff, 55 + Time_tag_type_len * 8, 3); slog(LOG_DEBUG, tab, "AI = %d", spartn->AI);
        spartn->EAL = getbitu(spartn->buff, 58 + Time_tag_type_len * 8, 3); slog(LOG_DEBUG, tab, "EAL = %d", spartn->EAL);
    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len) {
        spartn->Payload_offset = 6 + EA_Len + Time_tag_type_len;
        slog(LOG_DEBUG, tab, "Payload_offset = %d", spartn->Payload_offset);
    }
    else if (spartn->EAF == 1 && spartn->AI > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL) {

    }
    else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL + CRC_Len) {
        spartn->Message_CRC = getbitu(spartn->buff, (spartn->nbyte - CRC_Len) * 8, CRC_Len * 8);
        slog(LOG_DEBUG, tab, "nbyte = %d", spartn->nbyte);

        uint32_t Result_Message_CRC = crc24_radix(spartn->buff + 1, 5 + EA_Len + Time_tag_type_len + spartn->len + EADL);
        slog(LOG_DEBUG, tab, "Message_CRC = %d : %d", spartn->Message_CRC, Result_Message_CRC);
        if (spartn->Message_CRC == Result_Message_CRC) {
           expanded_full_time(spartn);
            slog(LOG_DEBUG, tab, "==========");
            spartn_out->type = spartn->type;
            spartn_out->Subtype = spartn->Subtype;
			spartn_out->time = spartn->GNSS_time_type;
            spartn_out->len = spartn->len;
            decode_spartn(spartn, spartn_out);
            slog(LOG_DEBUG, tab, "==========");
        }
        memset(spartn, 0, sizeof(raw_spartn_t));
        return 1;
    }
    return 0;
}

int check_nav(nav_t *nav, sap_ssr_t *sap_ssr)
{
    double time0 = sap_ssr[0].t0[1];
    unsigned int i,eph_n=0, geph_n = 0;
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
