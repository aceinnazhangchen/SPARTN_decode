#include "spartn.h"
#include "bits.h"
#include "log.h"
#define Leap_Sec 18.0
#define GLO_GPS_TD  10800

void decode_GPS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF011_Type = getbitu(data, offset, 2);  offset += 2; slog(LOG_DEBUG, tab, "SF011_Type = %d", SF011_Type);
	uint8_t SF011_Len = 0;
	switch (SF011_Type) {
	case 0:SF011_Len = 32; break;
	case 1:SF011_Len = 44; break;
	case 2:SF011_Len = 56; break;
	case 3:SF011_Len = 64; break;
	}
	*satellite_mask_len = SF011_Len;
	slog(LOG_DEBUG, tab, "SF011_Len = %d", SF011_Len);
	uint8_t SF011[8] = { 0 };
	bitscopy(SF011, 0, data, offset, SF011_Len); offset += SF011_Len;
	bits_to_bytes_array(SF011, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[65] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	slog(LOG_DEBUG, tab, "SF011 = %s", str_satellite_mask);
	*pos = offset;
}

void decode_GLONASS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF012_Type = getbitu(data, offset, 2);  offset += 2; slog(LOG_DEBUG, tab, "SF012_Type = %d", SF012_Type);
	uint8_t SF012_Len = 0;
	switch (SF012_Type) {
	case 0:SF012_Len = 24; break;
	case 1:SF012_Len = 36; break;
	case 2:SF012_Len = 48; break;
	case 3:SF012_Len = 63; break;
	}
	*satellite_mask_len = SF012_Len;
	slog(LOG_DEBUG, tab, "SF012_Len = %d", SF012_Len);
	uint8_t SF012[8] = { 0 };
	bitscopy(SF012, 0, data, offset, SF012_Len); offset += SF012_Len;
	bits_to_bytes_array(SF012, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[64] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	slog(LOG_DEBUG, tab, "SF012 = %s", str_satellite_mask);
	*pos = offset;
}

sap_ssr_t* suitable_ssr(spartn_t* spartn, int prn, int sys) {
    int j = 0, n = 0;
    sap_ssr_t* ssr = NULL;
    for (j = 0; j < SSR_NUM; j++) {
        ssr = &spartn->ssr[j];
        if (ssr->prn == 0) {
            spartn->ssr_offset++;
            break;
        }
        if (ssr->prn == prn && ssr->sys == sys) {
            break;
        }
        ssr = NULL;
    }
    double early_time = 0.0;
    if (ssr == NULL) {
        for (j = 0; j < SSR_NUM; j++) {
            ssr = &spartn->ssr[j];
            if (early_time == 0.0) {
                early_time = ssr->t0[1];
            }
            else {
                if (early_time > ssr->t0[1]) {
                    early_time = ssr->t0[1];
                    n = j;
                }
            }
        }
        ssr = &spartn->ssr[n];
        memset(ssr, 0, sizeof(sap_ssr_t));
    }
    return ssr;
}


void transform_spartn_ssr(raw_spartn_t* raw_spartn)
{
    int i = 0, j = 0, n = 0, m=0;
    if (!raw_spartn->spartn_out) return;
    spartn_t* spartn = raw_spartn->spartn_out;
    if (!spartn->ocb) return;
    if (!spartn->hpac) return;
    if (!spartn->gad) return;
    if (!spartn->lpac) return;
    OCB_t* ocb = spartn->ocb;
    HPAC_t* hpac = spartn->hpac;
    GAD_t* gad = spartn->gad;
    LPAC_t* lpac = spartn->lpac;
    sap_ssr_t* ssr = NULL;
    OCB_Satellite_t* sat_obc = NULL;
    gad_ssr_t* ssr_gad = NULL;
    //to full time
    //if (raw_spartn->Time_tag_type) {
    //    spartn->day = raw_spartn->GNSS_time_type / DAY_SECOND;
    //}
    //else {
    //    raw_spartn->GNSS_time_type = spartn->day*DAY_SECOND + raw_spartn->GNSS_time_type;
    //}

    if (spartn->type == 0) {
        //if (spartn->eos == 1) {
        //    spartn->ssr_offset = 0;
        //}
        for (i = 0; i < ocb->satellite_num; ++i) {
            //if (spartn->ssr_offset >= SSR_NUM)break;
            //ssr = &spartn->ssr[spartn->ssr_offset];
            sat_obc = &ocb->satellite[i];
            ssr = suitable_ssr(spartn, sat_obc->PRN_ID, raw_spartn->Subtype);
            
            ssr->prn = sat_obc->PRN_ID;
            ssr->sys = raw_spartn->Subtype;
            ssr->sat = ssr->sys*40+ ssr->prn;


            if (sat_obc->orbit.SF018_SF019_IODE != 0) {
                ssr->t0[0] = (double)raw_spartn->GNSS_time_type + (ssr->sys?(-GLO_GPS_TD + Leap_Sec):0);
                ssr->iod[0] = sat_obc->orbit.SF018_SF019_IODE;
                ssr->deph[0] = sat_obc->orbit.SF020_radial;
                ssr->deph[1] = sat_obc->orbit.SF020_along;
                ssr->deph[2] = sat_obc->orbit.SF020_cross;
                ssr->yaw_ang = sat_obc->orbit.SF021_Satellite_yaw;
            }

            ssr->t0[1] = (double)raw_spartn->GNSS_time_type + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
            ssr->iod[1] = sat_obc->clock.SF022_IODE_continuity;
            ssr->ure = sat_obc->clock.SF024_User_range_error;
            ssr->dclk = sat_obc->clock.SF020_Clock_correction;

            int update_num = 0;
            for (j = 0; j < Bias_Effective_Len; ++j) {
                if (sat_obc->GPS_bias.SF027_code_bias[j] == 1) {
                    ssr->cbias[j] = sat_obc->GPS_bias.SF029_Code_bias_correction[j];
                    update_num++;
                }
            }
            if (update_num > 0) {
                ssr->t0[2] = (double)raw_spartn->GNSS_time_type + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
            }
            update_num = 0;
            for (j = 0; j < Bias_Effective_Len; ++j) {
                if (sat_obc->GPS_bias.SF025_phase_bias[j] == 1) {
                    ssr->iod[j + 2]  = sat_obc->GPS_bias.Phase_bias[j].SF015_Continuity_indicator;
                    ssr->fix_flag[j] = sat_obc->GPS_bias.Phase_bias[j].SF023_Fix_flag;
                    ssr->pbias[j]    = sat_obc->GPS_bias.Phase_bias[j].SF020_Phase_bias_correction;
                    update_num++;
                }
            }
            if (update_num > 0) {
                ssr->t0[3] = (double)raw_spartn->GNSS_time_type + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
            }

            //spartn->ssr_offset++;
        }

        spartn->eos = ocb->header.SF010_EOS;
    }
    else if (spartn->type == 1) {
        for (i = 0; i < hpac->header.SF030_Area_count; ++i) {
            for (j = 0; j < SSR_NUM; ++j) {
                ssr = &spartn->ssr[j];
                ssr->rap_num = 0;
                if (ssr->prn == 0) break;
                if (ssr->sys != spartn->Subtype) continue;
                ssr->t0[4] = (double)raw_spartn->GNSS_time_type + (ssr->sys?(-GLO_GPS_TD + Leap_Sec):0);
                ssr->t0[5] = (double)raw_spartn->GNSS_time_type + (ssr->sys?(-GLO_GPS_TD + Leap_Sec):0);

                for (m = 0; m < RAP_NUM; ++m) {
                    if (ssr->areaId[m] == hpac->atmosphere[i].area.SF031_Area_ID) break;
                    if (ssr->areaId[m] == 0) break;
                }

                ssr->areaId[m] = hpac->atmosphere[i].area.SF031_Area_ID;
                ssr->ave_htd   = hpac->atmosphere[i].troposphere.SF043_Area_average_vertical_hydrostatic_delay;
                ssr->tro_coef[(m + ssr->rap_num) * 3 + 0] = hpac->atmosphere[i].troposphere.small_coefficient.SF045_T00;
                ssr->tro_coef[(m + ssr->rap_num) * 3 + 1] = hpac->atmosphere[i].troposphere.small_coefficient.SF046_T01;
                ssr->tro_coef[(m + ssr->rap_num) * 3 + 2] = hpac->atmosphere[i].troposphere.small_coefficient.SF046_T10;
                for (n = 0; n < hpac->atmosphere[i].ionosphere.ionosphere_satellite_num; ++n) {
                    if (ssr->prn == hpac->atmosphere[i].ionosphere.ionosphere_satellite[n].PRN_ID) {
                        ssr->stec_coef[(m + ssr->rap_num) * 3 + 0] = hpac->atmosphere[i].ionosphere.ionosphere_satellite[n].small_coefficient.SF057_C00;
                        ssr->stec_coef[(m + ssr->rap_num) * 3 + 1] = hpac->atmosphere[i].ionosphere.ionosphere_satellite[n].small_coefficient.SF058_C01;
                        ssr->stec_coef[(m + ssr->rap_num) * 3 + 2] = hpac->atmosphere[i].ionosphere.ionosphere_satellite[n].small_coefficient.SF058_C10;
                        break;
                    }
                }
            }         
        }
        ssr->rap_num += hpac->header.SF030_Area_count;
    }
    else if (spartn->type == 2) {
        for (i = 0; i < gad->header.SF030_Area_count; ++i) {
            for (j = 0; j < RAP_NUM; ++j) {
                ssr_gad = &spartn->ssr_gad[j];
                if (ssr_gad->areaId == gad->areas[i].SF031_Area_ID) {
                    break;
                }
                if (ssr_gad->areaId == 0) {
                    break;
                }
            }
            ssr_gad->areaId = gad->areas[i].SF031_Area_ID;
            ssr_gad->rap_lon = gad->areas[i].SF033_Area_reference_longitude;
            ssr_gad->rap_lat = gad->areas[i].SF032_Area_reference_latitude;
            ssr_gad->nc_lon = gad->areas[i].SF035_Area_longitude_grid_node_count;
            ssr_gad->nc_lat = gad->areas[i].SF034_Area_latitude_grid_node_count;
            ssr_gad->spa_lon = gad->areas[i].SF037_Area_longitude_grid_node_spacing;
            ssr_gad->spa_lat = gad->areas[i].SF036_Area_latitude_grid_node_spacing;
        }
    }
}


void expanded_full_time(raw_spartn_t* raw_spartn) {
    //to full time
    if (raw_spartn->Subtype == 0) {
    	if (raw_spartn->Time_tag_type) {
            raw_spartn->GNSS_time_type = raw_spartn->GNSS_time_type % DAY_SECOND;
    	}
    	else {
    		//if (GPS_DAYS > 0) {
    		//	raw_spartn->GNSS_time_type = GPS_DAYS * DAY_SECOND + raw_spartn->GNSS_time_type;
    		//	raw_spartn->Time_tag_type = 1;//modify Time_tag_type
    		//}
    	}
    }
    else if (raw_spartn->Subtype == 1) {//
    	if (raw_spartn->Time_tag_type) {
            raw_spartn->GNSS_time_type = raw_spartn->GNSS_time_type % DAY_SECOND;
    	}
    	else {
    		//if (GLONASS_DAYS > 0) {
    		//	raw_spartn->GNSS_time_type = GLONASS_DAYS * DAY_SECOND + raw_spartn->GNSS_time_type;
    		//	raw_spartn->Time_tag_type = 1;//modify Time_tag_type
    		//}
    	}
    }
    //raw_spartn->GNSS_time_type += 1262304000;
}