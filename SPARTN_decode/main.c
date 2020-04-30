#include <stdio.h>
#include "spartn.h"
#include "log.h"
#include "rinex.h"
#include "ephemeris.h"

#define SPARTN_2_RTCM
//#define READ_RTCM
//#define READ_RINEX

//int main() {
//
//	FILE* fp;
//	fp = fopen("..\\raw_data\\SPARTN20200317062037.raw", "rb");
//	if (fp == NULL) {
//		return 0;
//	}
//#ifdef TABLE_LOG
//	open_ocb_table_file(NULL);
//	open_hpac_table_file(NULL);
//	open_gad_table_file(NULL);
//	open_lpac_table_file(NULL);
//#endif
//	char buff = ' ';
//	size_t currentCount = 0;
//	size_t frameSize = 0;
//	size_t frameCount = 0;
//	size_t readCount = 0;
//	raw_spartn_t spartn;
//	memset(&spartn, 0, sizeof(spartn));
//	spartn_t spartn_out;
//	memset(&spartn_out, 0, sizeof(spartn_t));
//	OCB_t ocb = { 0 };
//	HPAC_t hpac = { 0 };
//	GAD_t gad = { 0 };
//	LPAC_t lpac = { 0 };
//	spartn_out.ocb = &ocb;
//	spartn_out.hpac = &hpac;
//	spartn_out.gad = &gad;
//	spartn_out.lpac = &lpac;
//	log(LOG_INFO, 0, "%5s,%3s,%4s,%4s,%3s,%5s,%3s,%4s,%11s,%3s,%4s,%3s,%3s,%2s,%3s,%4s", "frame", "crc", "type", "len", "EAF", "crc_t","sub","ti_t","time","sid","spid","EId","ESN","AI","EAL", "EADL");
//	int i;
//	while (!feof(fp)) {
//		memset(&buff, 0, sizeof(buff));
//		readCount = fread(&buff, sizeof(char), 1, fp);
//		if (readCount < 1)
//		{
//			/* file error or eof of file */
//			break;
//		}
//		currentCount += readCount;
//		frameSize += readCount;
//		if (buff == SPARTN_PREAMB) {
//			if (spartn.nbyte == 0) {
//				log(LOG_DEBUG, 0, "frameSize = %d", frameSize);
//			}
//		}
//		int ret = input_spartn_data(&spartn, &spartn_out, buff);
//		if (ret != 0) {
//			frameCount++;
//			if (ret == 1) {
//				char IVS[32] = { 0 };
//				for (i = 0; i < 12; ++i) {
//					sprintf(IVS+i*2,"%02X", spartn.IV[i]);
//					//Bp(spartn.IV[i]);
//				}
//				//printf("\n");
//				log(LOG_INFO, 0, "%5d,%3d,%4d,%4d,%3d,%5d,%3d,%4d,%11d,%3d,%4d,%3d,%3d,%2d,%3d,%4d,%s", frameCount, ret, spartn.type, spartn.len, spartn.EAF, spartn.CRC_type, spartn.Subtype, spartn.Time_tag_type,spartn.GNSS_time_type
//				,spartn.Solution_ID,spartn.Solution_processor_ID, spartn.Encryption_ID, spartn.ESN, spartn.AI, spartn.EAL, spartn.EADL*8, IVS);
//			}
//			else{
//				log(LOG_INFO, 0, "%5d,%3d,%4d,%4d,%3d", frameCount, ret, spartn.type, spartn.len, spartn.EAF);
//			}
//			frameSize = 0;
//			log(LOG_DEBUG, 0, "frame = %d", frameCount);
//			log(LOG_DEBUG, 0, "=======================");
//			memset(&spartn, 0, sizeof(raw_spartn_t));
//		}
//	}
//#ifdef TABLE_LOG
//	close_ocb_table_file();
//	close_hpac_table_file();
//	close_gad_table_file();
//	close_lpac_table_file();
//#endif
//	fclose(fp);
////	system("pause");
//	return 0;
//}

int gga_ssr2osr_main(FILE *fSSR, FILE *fEPH, FILE *fRTCM, double *ep, double *rovpos)  // int year, int doy, 
{
	gnss_rtcm_t rtcm = { 0 };
	nav_t *nav = &rtcm.nav;
	obs_t *rov = rtcm.obs;
	obs_t obs_vrs = { 0.0 };
	vec_t vec_vrs[MAXOBS] = { 0.0 };
	gtime_t time0 = epoch2time(ep);
	double cur_time = (int)time0.time;
	int doy = time2doy(time0);
	int year = ep[0];
	set_approximate_time(year, doy, rtcm.rcv);
	if (fSSR == NULL)  return 0;
	if (fEPH == NULL)  return 0;
#ifdef TABLE_LOG
	open_ocb_table_file(NULL);
	open_hpac_table_file(NULL);
	open_gad_table_file(NULL);
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
	spartn_out.ocb = &ocb;
	spartn_out.hpac = &hpac;
	spartn_out.gad = &gad;
	spartn_out.lpac = &lpac;
	int i, j, nsat;
	int rov_ret, ret_nav, num_ssr = -1;
	double blh[3] = { 0.0 }, dr[3] = { 0.0 };
	gtime_t teph = epoch2time(ep);
	double obs_time = 0.0;
	int nc = 0;
	while (1)
	{
		nav->ns = 0;
		nav->nsys[0] = 0;
		nav->nsys[1] = 0;
		fread_ssr_sapcorda(fSSR, &spartn, &spartn_out, nav->nsys);
		if (feof(fSSR)) break;

		/* read broadcast eph data one byte */
		ret_nav = fread_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
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
			//printf("ocb:%6.0f,%6.0f,%3i,%3i,%2i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n", sap_ssr[i].t0[0], sap_ssr[i].t0[1], nav_iod, sap_ssr[i].iod[0], sys, sap_ssr[i].prn, sap_ssr[i].deph[0], sap_ssr[i].deph[1], sap_ssr[i].deph[2], sap_ssr[i].dclk, sap_ssr[i].cbias[0], sap_ssr[i].cbias[1], sap_ssr[i].cbias[2], sap_ssr[i].pbias[0], sap_ssr[i].pbias[1], sap_ssr[i].pbias[2]);
		}

		while (1)
		{
			teph = timeadd(teph, 5.0);
			int time = teph.time;
			double time1 = fmod((double)time, 86400);
			if (fabs(time1 - sap_ssr[0].t0[1]) < 5.0)
				break;
		}

		nsat = satposs_sap_rcv(teph, rovpos, vec_vrs, nav, sap_ssr, EPHOPT_SSRSAP);

		obs_vrs.time = teph;
		obs_vrs.n = nsat;
		memcpy(obs_vrs.pos, rovpos, 3 * sizeof(double));
		for (i = 0; i < nsat; i++)
		{
			obs_vrs.data[i].sat = vec_vrs[i].sat;
		}
		nsat = compute_vector_data(&obs_vrs, vec_vrs);

		int vrs_ret = gen_obs_from_ssr(teph, rovpos, sap_ssr, sap_gad, &obs_vrs, vec_vrs, 0.0);
		for (i = 0; i < obs_vrs.n; ++i) {
			printf("obs: %12i,%3i,%14.4f,%14.4f,%14.4f,%14.4f\n",
				obs_vrs.time.time, obs_vrs.data[i].sat, obs_vrs.data[i].P[0], obs_vrs.data[i].P[1], obs_vrs.data[i].L[0], obs_vrs.data[i].L[1]);
		}
		printf("\n");

		rtcm_t out_rtcm = { 0 };
		unsigned char buffer[1200] = { 0 };
		int len = gen_rtcm_vrsdata(&obs_vrs, &out_rtcm, buffer);

		//obs_t obs_test = { 0 };
		//decode_rtcm3(&out_rtcm, &obs_test, NULL);

		fwrite(buffer, 1, len, fRTCM);
		nc++;
	}
	return 0;
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
	OCB_t  ocb = { 0 };
	HPAC_t hpac = { 0 };
	GAD_t  gad = { 0 };
	LPAC_t lpac = { 0 };
	spartn_out.ocb = &ocb;
	spartn_out.hpac = &hpac;
	spartn_out.gad = &gad;
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
		fread_ssr_sapcorda(fSSR, &spartn, &spartn_out, nav->nsys);
#else
		num_ssr = read_ssr_from_file(fSSR, &rtcm);
#endif
		if (feof(fSSR)) break;
		/* read broadcast eph data one byte */
		ret_nav = fread_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
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
			//printf("ocb:%7.0f,%7.0f,%2i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n", sap_ssr[i].t0[0]- cur_time0, sap_ssr[i].t0[1] - cur_time0, sys, sap_ssr[i].prn, sap_ssr[i].deph[0], sap_ssr[i].deph[1], sap_ssr[i].deph[2], sap_ssr[i].dclk, sap_ssr[i].cbias[0], sap_ssr[i].cbias[1], sap_ssr[i].cbias[2], sap_ssr[i].pbias[0], sap_ssr[i].pbias[1], sap_ssr[i].pbias[2]);
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

typedef struct
{
	double ver_;
	int sys_;
	int tsys_;
	char tobs[NUMSYS][MAXOBSTYPE][4];
	char type_;
	gtime_t start_time;
	gtime_t end_time;
}rinex_header;

int main()
{
#ifdef SPARTN_2_RTCM
	FILE *fSSR = { NULL };
	FILE *fEPH = { NULL };
	FILE *fROV = { NULL };
	FILE *fRTCM = NULL;

	fROV = fopen("..\\20200420\\SF0320111h.dat", "rb");
	fSSR = fopen("..\\20200420\\SPARTN20200420070130.raw", "rb");
	fEPH = fopen("..\\20200420\\Aux20200420070149.raw", "rb");
	fRTCM = fopen("..\\20200420\\SPARTN20200420070130.rtcm", "wb");
	int year = 2020;
	int doy = 111;
	//double ep[6] = { 2020,4,20,7,2,20 };
	double ep[6] = { 2020,4,20,7,1,50 };

	double rovpos[3] = { -2705297.408,-4283455.631,3861823.955 };

	//obs_ssr2osr_main(fSSR, fEPH, fROV, year, doy,ep, rovpos);

	gga_ssr2osr_main(fSSR, fEPH, fRTCM, ep, rovpos);

	if (fSSR) fclose(fSSR);
	if (fEPH) fclose(fEPH);
	if (fROV) fclose(fROV);
	if (fRTCM) fclose(fRTCM);
#endif
#ifdef  READ_RTCM
	FILE * fRTCM_IN = fopen("..\\20200420\\SLIB_1110_20200420070225.rtcm", "rb");
	FILE * fRTCM_OUT = fopen("..\\20200420\\SLIB_1110_20200420070225.rtcm_1", "wb");
	//FILE * fRTCM = fopen("..\\20200420\\SF0320111h.dat_1", "rb");
	//FILE * fRTCM_OUT = fopen("..\\20200420\\SF0320111h.dat_2", "wb");
	gnss_rtcm_t rtcm = { 0 };
	set_approximate_time(2020, 111, rtcm.rcv);
	obs_t* obs_vrs = &rtcm.obs[0];
	int ret = 0;
	int i = 0;
	while (ret != -1) {
		ret = read_obs_rtcm(fRTCM_IN, &rtcm, 0);
		if (ret == 1) {
			for (i = 0; i < obs_vrs->n; i++)
			{
				printf("obs: %12i,%3i,%14.4f,%14.4f,%14.4f,%14.4f\n",
					obs_vrs->time.time, obs_vrs->data[i].sat, obs_vrs->data[i].P[0], obs_vrs->data[i].P[1], obs_vrs->data[i].L[0], obs_vrs->data[i].L[1]);
			}
			printf("\n");

			rtcm_t out_rtcm = { 0 };
			unsigned char buffer[1200] = { 0 };
			int len = gen_rtcm_vrsdata(obs_vrs, &out_rtcm, buffer);

			obs_t obs_test = { 0 };
			decode_rtcm3(&out_rtcm, &obs_test, NULL);
			fwrite(buffer, 1, len, fRTCM_OUT);
		}
	}
	if (fRTCM_IN) fclose(fRTCM_IN);
	if (fRTCM_OUT) fclose(fRTCM_OUT);
#endif //  READ_RTCM
#ifdef READ_RINEX
	FILE * fRINEX_IN = fopen("..\\20200420\\SLIB_1110_20200420070225.rxo", "r");
	FILE * fRTCM_OUT = fopen("..\\20200420\\SLIB_1110_20200420070225.rtcm", "wb");
	rinex_header header = { 0 };
	nav_t navs = { 0 };
	sta_t sta = { 0 };
	obs_t obs = { 0 };
	int i = 0;
	int ret = 0;
	if (!readrnxh(fRINEX_IN, &header.ver_, &header.type_, &header.sys_, &header.tsys_, header.tobs, &navs, &sta, &(header.start_time), &(header.end_time)))
	{
		return 1;
	}
	while (ret != -1) {
		ret = readrnxobs_one(fRINEX_IN, header.ver_, &header.tsys_, header.tobs, &obs, &sta);
		if (ret > 0) {
			for (i = 0; i < obs.n; i++)
			{
				printf("obs: %12i,%3i,%14.4f,%14.4f,%14.4f,%14.4f\n",
					obs.time.time, obs.data[i].sat, obs.data[i].P[0], obs.data[i].P[1], obs.data[i].L[0], obs.data[i].L[1]);
			}
			printf("\n");

			rtcm_t out_rtcm = { 0 };
			unsigned char buffer[1200] = { 0 };
			int len = gen_rtcm_vrsdata(&obs, &out_rtcm, buffer);
			fwrite(buffer, 1, len, fRTCM_OUT);
		}
	}
	if (fRINEX_IN) fclose(fRINEX_IN);
	if (fRTCM_OUT) fclose(fRTCM_OUT);
#endif // READ_RINEX
	return 0;
}