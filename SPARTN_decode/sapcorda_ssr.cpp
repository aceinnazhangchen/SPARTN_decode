#include "log.h"
#include "sapcorda_ssr.h"
#include "ephemeris.h"
#include "GenVRSObs.h"
#include "export_ssr.h"
#include <string>
#include "stringex.h"

static double dmm2deg(double dmm)
{
	return floor(dmm / 100.0) + fmod(dmm, 100.0) / 60.0;
}

sapcorda_ssr * sapcorda_ssr::m_instance = NULL;
std::once_flag      sapcorda_ssr::m_flag;

sapcorda_ssr * sapcorda_ssr::getInstance()
{
	if (m_instance == NULL)
	{
		try
		{
			std::call_once(m_flag, createInstance);
		}
		catch (...)
		{
			printf("CreateInstance error\n");
		}
	}
	return m_instance;
}

void sapcorda_ssr::createInstance()
{
	m_instance = new(std::nothrow) sapcorda_ssr();
	if (NULL == m_instance)
	{
		throw std::exception();
	}
}

sapcorda_ssr::sapcorda_ssr()
{
	memset(&m_spartn, 0, sizeof(m_spartn));
	memset(&m_spartn_out, 0, sizeof(m_spartn_out));
	memset(&m_obs_vrs, 0, sizeof(m_obs_vrs));
	memset(&ocb, 0, sizeof(ocb));
	memset(&hpac, 0, sizeof(hpac));
	memset(&gad, 0, sizeof(gad));
	memset(&lpac, 0, sizeof(lpac));
	m_spartn_out.ocb = &ocb;
	m_spartn_out.hpac = &hpac;
	m_spartn_out.gad = &gad;
	m_spartn_out.lpac = &lpac;
}

sapcorda_ssr::~sapcorda_ssr()
{
}

void sapcorda_ssr::input_ssr_stream(unsigned char * buffer, uint32_t len)
{
	nav_t *nav = &m_rtcm.nav;
	sread_ssr_sapcorda(buffer, len, &m_spartn,&m_spartn_out, nav->nsys);
}

void sapcorda_ssr::input_eph_stream(unsigned char * buffer, uint32_t len)
{
	nav_t *nav = &m_rtcm.nav;
	sread_eph_rtcm(buffer, len, &m_rtcm, nav->nsys[0], nav->nsys[1]);
}

void input_ssr(unsigned char * buffer, uint32_t len)
{
	sapcorda_ssr::getInstance()->input_ssr_stream(buffer, len);
}

void input_eph(unsigned char * buffer, uint32_t len)
{
	sapcorda_ssr::getInstance()->input_eph_stream(buffer, len);
}

void input_gga(char * buffer, unsigned char*out_buffer, uint32_t *len)
{
	std::string gga = buffer;
	//double rovpos[3] = { -2705297.408,-4283455.631,3861823.955 };

	double pos[3] = { 0 };
	double xyz[3] = { 0 };

	std::vector<std::string> gga_split = split(gga, ",");
	if (gga_split.size() != 15)
	{
		return;
	}
	char  ns = 'N', ew = 'E';
	double lat = atof(gga_split[2].c_str());
	ns = gga_split[3][0];
	double lon = atof(gga_split[4].c_str());
	ew = gga_split[5][0];
	double alt = atof(gga_split[9].c_str());
	double msl = atof(gga_split[11].c_str());

	pos[0] = (ns == 'N' ? 1.0 : -1.0)*dmm2deg(lat)*D2R;
	pos[1] = (ew == 'E' ? 1.0 : -1.0)*dmm2deg(lon)*D2R;
	pos[2] = alt + msl;

	pos2ecef(pos, xyz);

	merge_ssr_to_obs(xyz,out_buffer,len);
}

unsigned char* merge_ssr_to_obs(double* rovpos, unsigned char*out_buffer, uint32_t *len)
{
	vec_t vec_vrs[MAXOBS] = { 0.0 };
	gtime_t teph = timeget();
	sap_ssr_t *sap_ssr = sapcorda_ssr::getInstance()->m_spartn_out.ssr;
	gad_ssr_t *sap_gad = sapcorda_ssr::getInstance()->m_spartn_out.ssr_gad;
	nav_t *nav = &sapcorda_ssr::getInstance()->m_rtcm.nav;
	int i, j, nsat;
	uint8_t ssr_offset = sapcorda_ssr::getInstance()->m_spartn_out.ssr_offset;
	for (i = 0; i < ssr_offset; i++)
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

	nsat = satposs_sap_rcv(teph, rovpos, vec_vrs, nav, sap_ssr, EPHOPT_SSRSAP);

	obs_t* obs_vrs = &sapcorda_ssr::getInstance()->m_obs_vrs;
	obs_vrs->time = teph;
	obs_vrs->n = nsat;
	memcpy(obs_vrs->pos, rovpos, 3 * sizeof(double));
	for (i = 0; i < nsat; i++)
	{
		obs_vrs->data[i].sat = vec_vrs[i].sat;
	}
	nsat = compute_vector_data(obs_vrs, vec_vrs);

	int vrs_ret = gen_obs_from_ssr(teph, rovpos, sap_ssr, sap_gad, obs_vrs, vec_vrs, 0.0);
	for (i = 0; i < obs_vrs->n; ++i) {
		printf("obs: %12i,%3i,%14.4f,%14.4f,%14.4f,%14.4f\n",
			obs_vrs->time.time, obs_vrs->data[i].sat, obs_vrs->data[i].P[0], obs_vrs->data[i].P[1], obs_vrs->data[i].L[0], obs_vrs->data[i].L[1]);
	}

	rtcm_t out_rtcm = { 0 };
	*len = gen_rtcm_vrsdata(obs_vrs, &out_rtcm, out_buffer);
	return out_buffer;
}
