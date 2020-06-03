#include <stdio.h>
#include "spartn.h"
#include "log.h"
#include "rinex.h"
#include "ephemeris.h"
#include "GenVRSObs.h"
#include "model.h"

#define SPARTN_2_RTCM
//#define READ_RTCM
//#define READ_RINEX

int match_ssr_nav_iode(sap_ssr_t *sap_ssr, nav_t *nav)
{
	int i, j, ng = 0;
	for (i = 0; i < nav->ns; i++)
	{
		int nav_iod = -1;
		int sys = sap_ssr[i].sys;
		if (sys == 0)
		{
			for (j = 0; j < nav->n; j++)
			{
				if (sap_ssr[i].prn == nav->eph[j].sat && sap_ssr[i].iod[0] == nav->eph[j].iode)
				{
					ng++;
					break;
				}
			}
		}
		else if (sys == 1)
		{
			for (j = 0; j < nav->ng; j++)
			{
				if (sap_ssr[i].prn + 40 == nav->geph[j].sat && sap_ssr[i].iod[0] == nav->geph[j].iode)
				{
					ng++;
					break;
				}
			}
		}
	}

	double ratio = (double)ng / nav->ns;
	if (ng >= 10 && ratio > 0.8)
		return 1;
	else
		return 0;
}

int gga_ssr2osr_main(FILE *fSSR, FILE *fEPH, FILE *fRTCM, FILE *fLOG, double *ep, double *rovpos)
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

	printf("spartn_t = %zd\n", sizeof(spartn_t));

	//printf("OCB_t = %zd\n",sizeof(OCB_t));
	//printf("HPAC_t = %zd\n", sizeof(HPAC_t));
	//printf("GAD_t = %zd\n", sizeof(GAD_t));
	//printf("LPAC_t = %zd\n", sizeof(LPAC_t));

	//printf("OCB_Satellite_t = %zd\n",sizeof(OCB_Satellite_t));
	//printf("HPAC_atmosphere_t = %zd\n", sizeof(HPAC_atmosphere_t));
	//printf("GAD_area_t = %zd\n", sizeof(GAD_area_t));
	//printf("LPAC_area_t = %zd\n", sizeof(LPAC_area_t));

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
		nav->ns = nav->nsys[0] + nav->nsys[1];
		if (feof(fSSR)) break;
		int epffEPH = 0;
		while (1)
		{
			/* read broadcast eph data one byte */
			ret_nav = fread_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
			if (ret_nav != 2)
			{
				/* can not find the complete epoch data */
				if (feof(fEPH))
				{
					epffEPH = 1;
					break;
				}
			}
			if (match_ssr_nav_iode(sap_ssr, nav))
			{
				break;
			}
		}
		if (epffEPH == 1) break;

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
			double nav_toe = (sys == 0) ? fmod(nav->eph[j].toe.time, 86400) : fmod(nav->geph[j].toe.time, 86400);
			fprintf(fLOG, "ocb:%6.0f,%6.0f,%6.0f,%6.0f,%6.0f,%3i,%3i,%2i,%3i,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n",
				sap_ssr[i].t0[0], sap_ssr[i].t0[1], sap_ssr[i].t0[2], sap_ssr[i].t0[4], nav_toe, nav_iod, sap_ssr[i].iod[0], sys, sap_ssr[i].prn,
				sap_ssr[i].deph[0], sap_ssr[i].deph[1], sap_ssr[i].deph[2], sap_ssr[i].dclk,
				sap_ssr[i].cbias[0], sap_ssr[i].cbias[1], sap_ssr[i].cbias[2], sap_ssr[i].pbias[0], sap_ssr[i].pbias[1], sap_ssr[i].pbias[2]);
		}

		while (1)
		{
			int time = teph.time;
			double time1 = fmod((double)time, DAY_SECONDS);
			if (time1 - sap_ssr[0].t0[1] < 20.0 && time1>sap_ssr[0].t0[1])
				break;

			teph = timeadd(teph, 1.0);
		}
		memset(&obs_vrs, 0, sizeof(obs_vrs));

		nsat = satposs_sap_rcv(teph, rovpos, vec_vrs, nav, sap_ssr, EPHOPT_SSRSAP);

		obs_vrs.time = teph;
		obs_vrs.n = nsat;

		memcpy(obs_vrs.pos, rovpos, 3 * sizeof(double));
		for (i = 0; i < nsat; i++)
		{
			obs_vrs.data[i].sat = vec_vrs[i].sat;
		}
		nsat = compute_vector_data(&obs_vrs, vec_vrs);

		if (nsat == 0)  continue;

		int vrs_ret = gen_obs_from_ssr(teph, rovpos, sap_ssr, sap_gad, &obs_vrs, vec_vrs, 0.0, fLOG);

		rtcm_t out_rtcm = { 0 };
		unsigned char buffer[1200] = { 0 };
		int len = gen_rtcm_vrsdata(&obs_vrs, &out_rtcm, buffer);

		fwrite(buffer, 1, len, fRTCM);

		//obs_t obs_test = { 0 };
		//decode_rtcm3(&out_rtcm, &obs_test, NULL);
		nc++;
	}
	return 0;
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

void correction_diff(FILE *fCOR, FILE *fLOG, FILE *fDIF)
{
    sap_cor_dif cor_dif = { 0 };
    char cor1[1024] = { 0 }, cor2[1024] = { 0 };
    char buffer[512] = { 0 };
    if (fCOR != NULL) fgets(cor1, sizeof(cor1), fCOR);
    if (fLOG != NULL) fgets(cor2, sizeof(cor2), fLOG);

}

int process(const char *fname, const char *root_dir)
{
    FILE *fSSR = { NULL };
    FILE *fEPH = { NULL };
    FILE *fROV = { NULL };
    FILE *fLOG = { NULL };
    FILE *fINI = fopen(fname, "r");
    char buffer[255];
    char fname1[255] = { 0 };
    char fname2[255] = { 0 };
    char fname3[255] = { 0 };
    char fname4[255] = { 0 };
    char inp_dir[255] = { 0 };
    char out_dir[255] = { 0 };
    int year = 0, doy = 0, line = 0, type = 0, isPrint = 0, isObs = 0, num = 0;
    double ep[6] = { 0 };
    double rovpos[3] = { 0.0 };
    double refpos[3] = { 0.0 };
    if (root_dir) {
        strncpy(inp_dir, root_dir, strlen(root_dir));
    }

    while (fINI != NULL && !feof(fINI))
    {
        memset(buffer, 0, sizeof(buffer));
        if (fgets(buffer, sizeof(buffer), fINI) == NULL) break;
        if (strlen(buffer) < 2) continue;
        if (buffer[0] == '#' || buffer[0] == ';') continue;
        memset(fname1, 0, sizeof(fname1));
        memset(fname2, 0, sizeof(fname2));
        memset(fname3, 0, sizeof(fname3));
        memset(fname4, 0, sizeof(fname4));
        year = 0;
        doy = 0;
        type = 0;
        num = sscanf(buffer, "%i", &type);
        num = sscanf(buffer, "%i", &type);
        switch (type) {
        case 0: /* RTK data process */
        {
            strncpy(fname1, inp_dir, strlen(inp_dir));
            strncpy(fname2, inp_dir, strlen(inp_dir));
            strncpy(fname3, inp_dir, strlen(inp_dir));
            strncpy(fname4, inp_dir, strlen(inp_dir));
            num = sscanf(buffer, "%i,%[^\,],%[^\,],%[^\,],%[^\,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &type, fname1 + strlen(inp_dir), fname2 + strlen(inp_dir), fname3 + strlen(inp_dir), fname4 + strlen(inp_dir),
                &refpos[0], &refpos[1], &refpos[2], &ep[0], &ep[1], &ep[2], &ep[3], &ep[4], &ep[5]);
            fSSR = fopen(fname1, "rb");
            fEPH = fopen(fname2, "rb");
            FILE * fRTCM_OUT = fopen(fname3, "wb");
            FILE * fLOG = fopen(fname4, "wb");
            gga_ssr2osr_main(fSSR, fEPH, fRTCM_OUT, fLOG, ep, refpos);

        }break;
        case 1: /* RTK data process */
        {
            gnss_rtcm_t rtcm = { 0 };
            nav_t *nav = &rtcm.nav;
            strncpy(fname1, inp_dir, strlen(inp_dir));
            num = sscanf(buffer, "%i,%[^\,],%lf,%lf,%lf,%lf,%lf,%lf", &type, fname1 + strlen(inp_dir),
                &ep[0], &ep[1], &ep[2], &ep[3], &ep[4], &ep[5]);
            fEPH = fopen(fname1, "rb");
            while (1)
            {
                nav->ns = 0;
                nav->nsys[0] = 0;
                nav->nsys[1] = 0;
                /* read broadcast eph data one byte */
                int ret_nav = fread_eph_rtcm(fEPH, &rtcm, nav->nsys[0], nav->nsys[1]);
                if (ret_nav != 2)
                {
                    /* can not find the complete epoch data */
                    if (feof(fEPH)) break;
                }
            }

        }break;
        case 4: /* input directory, effective after this command */
        {
            if (root_dir) {
                num = sscanf(buffer, "%i,%[^\,]", &type, inp_dir + strlen(root_dir));
            }
            else {
                num = sscanf(buffer, "%i,%[^\,]", &type, inp_dir);
            }

            char* temp = strchr(inp_dir, '\n');
            if (temp != NULL) temp[0] = '\0';
            if (strlen(inp_dir) > 0)
            {
                if (inp_dir[strlen(inp_dir) - 1] != '\\')
                {
                    inp_dir[strlen(inp_dir)] = '\\';
                }
            }
        } break;
        default:
            break;
        }
        ++line;
    }
    if (fINI != NULL) fclose(fINI);
    return 0;
}

int main()
{
    process("..\\data.ini", NULL);
    return 0;
}
