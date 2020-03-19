#include <stdio.h>
#include <memory.h>
#include "crc.h"
#include "bits.h"
#include "log.h"
#include "spartn.h"
#include "rtcm.h"
#include "GenVRSObs.h"
#include "ephemeris.h"

#define Message_Type 2 //0:OCB 1:HPAC 2:GAD 3:LPAC

int read_obs_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm, int stnID);
int read_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm);
//int read_ssr_sapcorda(FILE *fRTCM, gnss_rtcm_t *rtcm);

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

int main() 
{
    FILE *fSSR = { NULL };
    FILE *fEPH = { NULL };
    FILE *fROV = { NULL };
    fROV = fopen("..\\20200317\\SF0320077g.dat", "rb");
    fSSR = fopen("..\\20200317\\SPARTN20200317062037.raw","rb");
    fEPH = fopen("..\\20200317\\Aux20200317062038.raw", "rb");
    if (fROV == NULL)  return 0;
	if (fSSR == NULL)  return 0;
    if (fEPH == NULL)  return 0;
#ifdef TABLE_LOG
	open_ocb_table_file (NULL);
	open_hpac_table_file(NULL);
	open_gad_table_file (NULL);
	open_lpac_table_file(NULL);
#endif

    gnss_rtcm_t rtcm = { 0 };
    raw_spartn_t raw_spartn = { 0 };
    spartn_t spartn = { 0 };
    nav_t *nav = &rtcm.nav;
    obs_t *rov = rtcm.obs;
    sap_ssr_t *ssr = &spartn.ssr;
    gad_ssr_t *gad = &spartn.gad;
    vec_t vec_vrs[MAXOBS];
    obs_t *obs_vrs = NULL;

    while (1)
    {
      int rov_ret = read_obs_rtcm(fROV, &rtcm, 0);
        /* read broadcast eph data one byte */
      int ret_nav = read_eph_rtcm(fEPH, &rtcm);

      if (ret_nav != 2) 
      {
          /* can not find the complete epoch data */
          if (feof(fEPH)) break;
          continue;
      }

      /* no brdc */
      if (nav->n == 0) continue;

      int ret_ssr = read_ssr_sapcorda(fSSR, ssr, gad);

      satposs_sap(rov, vec_vrs, nav, ssr, EPHOPT_SSRSAP);
      compute_vector_data(&rov, vec_vrs);

     memset(&obs_vrs,0, sizeof(obs_t));
     int vrs_ret = gen_vobs_from_ssr(rov, ssr, gad, obs_vrs, vec_vrs, 10.0);


    }

#ifdef TABLE_LOG
	close_ocb_table_file();
	close_hpac_table_file();
	close_gad_table_file();
	close_lpac_table_file();
#endif
	fclose(fSSR);
	system("pause");
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
       // *numofread += (int)readCount;
        ret = input_rtcm3(buff, stnID, rtcm);
        if (ret == 1)
        {
            break;
        }
    }
    return ret;
}

int read_eph_rtcm(FILE *fRTCM, gnss_rtcm_t *rtcm)
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
        if (ret == 2 && rtcm->nav.n_gps>0 && rtcm->nav.ng > 0)
        {
            break;
        }
    }
    return ret;
}

int read_ssr_sapcorda(FILE *fSSR, raw_spartn_t *raw_spartn, spartn_t *spartn)
{
    //raw_spartn_t raw_spartn = { 0 };
    //spartn_t spartn = { 0 };
    //memset(&raw_spartn.spartn_out, 0, sizeof(spartn_t));
    //memset(&spartn    , 0, sizeof(spartn    ));
    //ssr = spartn->ssr;
    //gad = spartn->gad;

    int ret = 0;
    char buff = ' ';
    size_t currentCount = 0;
    size_t frameSize    = 0;
    size_t frameCount   = 0;
    size_t readCount    = 0;
    while (!feof(fSSR))
    {
        memset(&buff, 0, sizeof(buff));
        readCount = fread(&buff, sizeof(char), 1, fSSR);
        if (readCount < 1)
        {
            /* file error or eof of file */
            break;
        }
        currentCount += readCount;
        frameSize += readCount;
        int ret = input_spartn_data(raw_spartn, spartn, buff);

        if (ret == 1 && spartn->ssr[0].prn!=0)
        {
            frameSize = 0;
            frameCount++;
            break;
        }
    }

    return ret;
}