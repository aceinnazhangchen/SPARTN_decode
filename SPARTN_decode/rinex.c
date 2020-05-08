#include "rinex.h"
#ifdef WIN32
#include<windows.h>
#endif
#include <assert.h>

#define MAXPOSHEAD  1024                /* max head line position */

/* adjust time considering week handover -------------------------------------*/
static gtime_t adjweek_(gtime_t t, gtime_t t0)
{
	double tt = timediff(t, t0);
	if (tt < -302400.0) return timeadd(t, 604800.0);
	if (tt > 302400.0) return timeadd(t, -604800.0);
	return t;
}

/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int str2time(const char *s, int i, int n, gtime_t *t)
{
	double ep[6];
	char str[256], *p = str;

	if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < i) return -1;
	for (s += i; *s&&--n >= 0;) *p++ = *s++;
	*p = '\0';
	if (sscanf(str, "%lf %lf %lf %lf %lf %lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5) < 6)
		return -1;
	if (ep[0] < 100.0) ep[0] += ep[0] < 80.0 ? 2000.0 : 1900.0;
	*t = epoch2time(ep);
	return 0;
}

/* read rinex header ---------------------------------------------------------*/
int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
	char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta, gtime_t *start_time, gtime_t *end_time)
{
//	double bias;
	char buff[MAXRNXLEN], *label = buff + 60;
	int i = 0, block = 0;

	trace(3, "readrnxh:\n");

	*ver = 2.10; *type = ' '; *sys = _SYS_GPS_;

	while (fgets(buff, MAXRNXLEN, fp))
	{

		if (strlen(buff) <= 60) continue;

		else if (strstr(label, "RINEX VERSION / TYPE"))
		{
			*ver = str2num(buff, 0, 9);
			*type = *(buff + 20);

			/* satellite system */
			switch (*(buff + 40))
			{
			case ' ':
			case 'G': *sys = _SYS_GPS_;  *tsys = TSYS_GPS; break;
			case 'R': *sys = _SYS_GLO_;  *tsys = TSYS_UTC; break;
			case 'E': *sys = _SYS_GAL_;  *tsys = TSYS_GAL; break; /* v.2.12 */
			case 'S': *sys = _SYS_SBS_;  *tsys = TSYS_GPS; break;
			case 'J': *sys = _SYS_QZS_;  *tsys = TSYS_QZS; break; /* v.3.02 */
			case 'C': *sys = _SYS_BDS_;  *tsys = TSYS_CMP; break; /* v.2.12 */
			case 'I': *sys = _SYS_IRN_;  *tsys = TSYS_IRN; break; /* v.3.03 */
			case 'M': *sys = _SYS_NONE_; *tsys = TSYS_GPS; break; /* mixed */
			default:
				trace(2, "not supported satellite system: %c\n", *(buff + 40));
				break;
			}
			continue;
		}
		else if (strstr(label, "PGM / RUN BY / DATE")) continue;
		else if (strstr(label, "COMMENT"))
		{ /* opt */
			/* read cnes wl satellite fractional bias */
			if (strstr(buff, "WIDELANE SATELLITE FRACTIONAL BIASES") ||
				strstr(buff, "WIDELANE SATELLITE FRACTIONNAL BIASES"))
			{
				block = 1;
			}
			//else if (block)
			//{
			//	/* cnes/cls grg clock */
			//	if (!strncmp(buff, "WL", 2) && (sat = satid2no(buff + 3)) &&
			//		sscanf(buff + 40, "%lf", &bias) == 1)
			//	{
			//		nav->wlbias[sat - 1] = bias;
			//		nav->upd_from = 1;
			//	}
			//	/* cnes ppp-wizard clock */
			//	else if ((sat = satid2no(buff + 1)) && sscanf(buff + 6, "%lf", &bias) == 1) {
			//		nav->wlbias[sat - 1] = bias;
			//		nav->upd_from = 1;
			//	}
			//}
			continue;
		}
		else if (strstr(label, "TIME OF FIRST OBS")) {
			if (start_time)
			{
				str2time(buff, 2, 44, start_time);
			}
		}
		else if (strstr(label, "TIME OF LAST OBS")){
			if (end_time)
			{
				str2time(buff, 2, 44, end_time);
			}
		}
		/* file type */
		switch (*type)
		{
		case 'O': decode_obsh(fp, buff, *ver, tsys, tobs, nav, sta); break;
		case 'N': decode_navh(buff, nav); break;
		case 'G': decode_gnavh(buff, nav); break;
		case 'H': decode_hnavh(buff, nav); break;
		case 'J': decode_navh(buff, nav); break; /* extension */
		case 'L': decode_navh(buff, nav); break; /* extension */
		}
		if (strstr(label, "END OF HEADER")) return 1;

		if (++i >= MAXPOSHEAD && *type == ' ') break; /* no rinex file */
	}
	return 0;
}

/* read rinex obs data body --------------------------------------------------*/
int readrnxobsb(FILE *fp, double ver, int *tsys, char tobs[][MAXOBSTYPE][4], int *flag, obsd_t *data, gtime_t *time,
	sta_t *sta)
{
	//gtime_t time = { 0 };
	sigind_t index[7] = { {0} };
	char buff[MAXRNXLEN];
	int i = 0, n = 0, nsat = 0, sats[MAXOBS] = { 0 };

	/* set signal index */
	set_index(ver, _SYS_GPS_, tobs[0], index);
	set_index(ver, _SYS_GLO_, tobs[1], index + 1);
	set_index(ver, _SYS_GAL_, tobs[2], index + 2);
	set_index(ver, _SYS_QZS_, tobs[3], index + 3);
	set_index(ver, _SYS_BDS_, tobs[5], index + 5);

	/* read record */
	while (fgets(buff, MAXRNXLEN, fp))
	{

		/* decode obs epoch */
		if (i == 0) {
			nsat = decode_obsepoch(fp, buff, ver, time, flag, sats);
			trace(2, buff);
			if (nsat <= 0) {
				continue;
			}
		}
		else if (*flag <= 2 || *flag == 6)
		{

			data[n].sat = (unsigned char)sats[i - 1];

			/* decode obs data */
			if (decode_obsdata(fp, buff, ver, index, data + n) && n < MAXOBS) n++;
		}
		else if (*flag == 3 || *flag == 4) { /* new site or header info follows */

			/* decode obs header */
			decode_obsh(fp, buff, ver, tsys, tobs, NULL, sta);
		}
		if (++i > nsat) return n;
	}
	return -1;
}

/* read rinex obs ------------------------------------------------------------*/
int readrnxobs(FILE *fp, double ver, int *tsys, char tobs[][MAXOBSTYPE][4], obs_t *obs, sta_t *sta)
{
	obsd_t *data;
	int i, n, flag = 0, stat = 0;

	//trace(4, "readrnxobs: rcv=%d ver=%.2f tsys=%d\n", rcv, ver, tsys);

	if (!obs||fp==NULL) return 0;

	if (!(data = (obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) return 0;

	/* read rinex obs data body */
	while ((n = readrnxobsb(fp, ver, tsys, tobs, &flag, data, &obs->time, sta)) >= 0 && stat >= 0)
	{
		obs->n = 0;
		for (i = 0; i < n; i++)
		{
			/* save obs data */
			if ((stat = addobsdata(obs, data + i)) < 0) break;
		}
	}
	trace(4, "readrnxobs: nobs=%d stat=%d\n", obs->n, stat);

	free(data);

	return stat;
}
/* read rinex obs ------------------------------------------------------------*/
int readrnxobs_one(FILE *fp, double ver, int *tsys, char tobs[][MAXOBSTYPE][4], obs_t *obs, sta_t *sta)
{
	obsd_t data[MAXOBS];
	int i, n, flag = 0;

	if (!obs) return 0;

	/* read rinex obs data body */
	if ((n = readrnxobsb(fp, ver, tsys, tobs, &flag, data, &obs->time, sta)) >= 0)
	{
		obs->n = 0;
		for (i = 0; i < n; i++)
		{
			/* save obs data */
			addobsdata(obs, data + i);
		}
	}
	trace(4, "readrnxobs: nobs=%d stat=%d\n", obs->n);

	return n;
}

/* decode nav header ---------------------------------------------------------*/
static void decode_navh(char *buff, nav_t *nav)
{
//	int i, j;
	char *label = buff + 60;

	trace(4, "decode_navh:\n");

}

/* decode gnav header --------------------------------------------------------*/
static void decode_gnavh(char *buff, nav_t *nav)
{
	//char *label = buff + 60;

	//trace(4, "decode_gnavh:\n");

	//if (strstr(label, "CORR TO SYTEM TIME")); /* opt */
	//else if (strstr(label, "LEAP SECONDS")) { /* opt */
	//	if (nav) nav->leaps = (int)str2num(buff, 0, 6);
	//}
}
/* decode geo nav header -----------------------------------------------------*/
static void decode_hnavh(char *buff, nav_t *nav)
{
	//char *label = buff + 60;

	//trace(4, "decode_hnavh:\n");

	//if (strstr(label, "CORR TO SYTEM TIME")); /* opt */
	//else if (strstr(label, "D-UTC A0,A1,T,W,S,U")); /* opt */
	//else if (strstr(label, "LEAP SECONDS")) { /* opt */
	//	if (nav) nav->leaps = (int)str2num(buff, 0, 6);
	//}
}

/* add obs data --------------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
	trace(2, "obsd: sat=%d,C=%.3f,L=%.3f,D=%.3f,S=%d,C=%.3f,L=%.3f,D=%.3f,S=%d\n"
		, data->sat, data->P[0], data->L[0], data->D[0], data->SNR[0], data->P[1], data->L[1], data->D[1], data->SNR[1]);
	obs->data[obs->n++] = *data;
	return 1;
}


/* decode obs header ---------------------------------------------------------*/
static void decode_obsh(FILE *fp, char *buff, double ver, int *tsys,
	char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta)
{
	/* default codes for unknown code */
	const char *defcodes[] = {
		"CWX    ",  /* GPS: L125____ */
		"CC     ",  /* GLO: L12_____ */
		"X XXXX ",  /* GAL: L1_5678_ */
		"CXXX   ",  /* QZS: L1256___ */
		"C X    ",  /* SBS: L1_5____ */
		"X  XX  ",  /* BDS: L1__67__ */
		"  A   A"   /* IRN: L__5___9 */
	};
	double del[3];
	int i, j, k, n, nt;
	const char *p;
	char *label = buff + 60, str[4];

	if (strstr(label, "MARKER NAME")) {
		if (sta) setstr(sta->name, buff, 60);
	}
	else if (strstr(label, "MARKER NUMBER")) { /* opt */
		if (sta) setstr(sta->marker, buff, 20);
	}
	else if (strstr(label, "MARKER TYPE")); /* ver.3 */
	else if (strstr(label, "OBSERVER / AGENCY"));
	else if (strstr(label, "REC # / TYPE / VERS")) {
		if (sta) {
			setstr(sta->recsno, buff, 20);
			setstr(sta->rectype, buff + 20, 20);
			setstr(sta->recver, buff + 40, 20);
		}
	}
	else if (strstr(label, "ANT # / TYPE"))
	{
		if (sta)
		{
			setstr(sta->antsno, buff, 20);
			setstr(sta->antdes, buff + 20, 20);
		}

	}
	else if (strstr(label, "APPROX POSITION XYZ")) {
		if (sta) {
			for (i = 0, j = 0; i < 3; i++, j += 14) sta->pos[i] = str2num(buff, j, 14);
		}
	}
	else if (strstr(label, "ANTENNA: DELTA H/E/N")) {
		if (sta) {
			for (i = 0, j = 0; i < 3; i++, j += 14) del[i] = str2num(buff, j, 14);
			sta->del[2] = del[0]; /* h */
			sta->del[0] = del[1]; /* e */
			sta->del[1] = del[2]; /* n */
		}

	}
	else if (strstr(label, "ANTENNA: DELTA X/Y/Z")); /* opt ver.3 */
	else if (strstr(label, "ANTENNA: PHASECENTER")); /* opt ver.3 */
	else if (strstr(label, "ANTENNA: B.SIGHT XYZ")); /* opt ver.3 */
	else if (strstr(label, "ANTENNA: ZERODIR AZI")); /* opt ver.3 */
	else if (strstr(label, "ANTENNA: ZERODIR XYZ")); /* opt ver.3 */
	else if (strstr(label, "CENTER OF MASS: XYZ")); /* opt ver.3 */
	else if (strstr(label, "SYS / # / OBS TYPES")) { /* ver.3 */
		if (!(p = strchr(syscodes, buff[0]))) {
			trace(2, "invalid system code: sys=%c\n", buff[0]);
			return;
		}
		i = (int)(p - syscodes);
		n = (int)str2num(buff, 3, 3);
		for (j = nt = 0, k = 7; j < n; j++, k += 4) {
			if (k > 58) {
				if (!fgets(buff, MAXRNXLEN, fp)) break;
				k = 7;
			}
			if (nt < MAXOBSTYPE - 1) setstr(tobs[i][nt++], buff + k, 3);
		}
		*tobs[i][nt] = '\0';

		/* change beidou B1 code: 3.02 */
		if (i == 5 && fabs(ver - 3.02) < 1e-3) {
			for (j = 0; j < nt; j++) if (tobs[i][j][1] == '2') tobs[i][j][1] = '1';
		}
		/* if unknown code in ver.3, set default code */
		for (j = 0; j < nt; j++) {
			if (tobs[i][j][2]) continue;
			if (!(p = strchr(frqcodes, tobs[i][j][1]))) continue;
			tobs[i][j][2] = defcodes[i][(int)(p - frqcodes)];
			trace(2, "set default for unknown code: sys=%c code=%s\n", buff[0],
				tobs[i][j]);
		}
	}
	else if (strstr(label, "WAVELENGTH FACT L1/2")); /* opt ver.2 */
	else if (strstr(label, "# / TYPES OF OBSERV")) { /* ver.2 */
		n = (int)str2num(buff, 0, 6);
		for (i = nt = 0, j = 10; i < n; i++, j += 6) {
			if (j > 58) {
				if (!fgets(buff, MAXRNXLEN, fp)) break;
				j = 10;
			}
			if (nt >= MAXOBSTYPE - 1) continue;
			if (ver <= 2.99) {
				setstr(str, buff + j, 2);
				convcode(ver, _SYS_GPS_, str, tobs[0][nt]);
				convcode(ver, _SYS_GLO_, str, tobs[1][nt]);
				convcode(ver, _SYS_GAL_, str, tobs[2][nt]);
				convcode(ver, _SYS_QZS_, str, tobs[3][nt]);
				convcode(ver, _SYS_BDS_, str, tobs[5][nt]);
			}
			nt++;
		}
		*tobs[0][nt] = '\0';
	}
	else if (strstr(label, "SIGNAL STRENGTH UNIT")); /* opt ver.3 */
	else if (strstr(label, "INTERVAL")); /* opt */
	else if (strstr(label, "TIME OF FIRST OBS")) {
		if (!strncmp(buff + 48, "GPS", 3)) *tsys = TSYS_GPS;
		else if (!strncmp(buff + 48, "GLO", 3)) *tsys = TSYS_UTC;
		else if (!strncmp(buff + 48, "GAL", 3)) *tsys = TSYS_GAL;
		else if (!strncmp(buff + 48, "QZS", 3)) *tsys = TSYS_QZS; /* ver.3.02 */
		else if (!strncmp(buff + 48, "BDT", 3)) *tsys = TSYS_CMP; /* ver.3.02 */
	}

}

/* set signal index ----------------------------------------------------------*/
static void set_index(double ver, int sys, char tobs[MAXOBSTYPE][4], sigind_t *ind)
{
	const char *p;
//  char str[8];
    char *optstr = "";
//	double shift;
	int i, j, k, n;

	for (i = n = 0; *tobs[i]; i++, n++) {
		ind->code[i] = obs2code(sys, tobs[i] + 1, ind->frq + i);
		ind->type[i] = (p = strchr(obscode, tobs[i][0])) ? (int)(p - obscode) : 0;
		ind->pri[i] = getcodepri(sys, ind->code[i], NULL);
		ind->pos[i] = -1;

		/* frequency index for beidou */
		if (sys == _SYS_BDS_)
		{
			if (ind->frq[i] == 5) ind->frq[i] = 2; /* B2 */
			else if (ind->frq[i] == 4) ind->frq[i] = 3; /* B3 */
		}
	}
	/* parse phase shift options */
	switch (sys) {
	case _SYS_GPS_: optstr = "-GL%2s=%lf"; break;
	case _SYS_GLO_: optstr = "-RL%2s=%lf"; break;
	case _SYS_GAL_: optstr = "-EL%2s=%lf"; break;
	case _SYS_QZS_: optstr = "-JL%2s=%lf"; break;
	case _SYS_BDS_: optstr = "-CL%2s=%lf"; break;

	}
	//for (p = opt; p && (p = strchr(p, '-')); p++) {
	//    if (sscanf(p, optstr, str, &shift) < 2) continue;
	//    for (i = 0; i < n; i++) {
	//        if (strcmp(code2obs(sys,ind->code[i], NULL), str)) continue;
	//        ind->shift[i] = shift;
	//        trace(2, "phase shift: sys=%2d tobs=%s shift=%.3f\n", sys,
	//            tobs[i], shift);
	//    }
	//}
	/* assign index for highest priority code */
	for (i = 0; i < NFREQ; i++) {
		for (j = 0, k = -1; j < n; j++) {
			if (ind->frq[j] == i + 1 && ind->pri[j] && (k<0 || ind->pri[j]>ind->pri[k])) {
				k = j;
			}
		}
		if (k < 0) continue;

		for (j = 0; j < n; j++) {
			if (ind->code[j] == ind->code[k]) ind->pos[j] = i;
		}
	}
	/* assign index of extended obs data */
	for (i = 0; i < NEXOBS; i++) {
		for (j = 0; j < n; j++) {
			if (ind->code[j] && ind->pri[j] && ind->pos[j] < 0) break;
		}
		if (j >= n) break;

		for (k = 0; k < n; k++) {
			if (ind->code[k] == ind->code[j]) ind->pos[k] = NFREQ + i;
		}
	}
	for (i = 0; i < n; i++) {
		if (!ind->code[i] || !ind->pri[i] || ind->pos[i] >= 0) continue;
		trace(4, "reject obs type: sys=%2d, obs=%s\n", sys, tobs[i]);
	}
	ind->n = n;

#if 0 /* for debug */
	for (i = 0; i < n; i++) {
		trace(2, "set_index: sys=%2d,tobs=%s code=%2d pri=%2d frq=%d pos=%d shift=%5.2f\n",
			sys, tobs[i], ind->code[i], ind->pri[i], ind->frq[i], ind->pos[i],
			ind->shift[i]);
	}
#endif
}

/* decode obs epoch ----------------------------------------------------------*/
static int decode_obsepoch(FILE *fp, char *buff, double ver, gtime_t *time,
	int *flag, int *sats)
{
	int i, j, n;
	char satid[8] = "";

	trace(4, "decode_obsepoch: ver=%.2f\n", ver);

	if (ver <= 2.99) { /* ver.2 */
		if ((n = (int)str2num(buff, 29, 3)) <= 0) return 0;

		/* epoch flag: 3:new site,4:header info,5:external event */
		*flag = (int)str2num(buff, 28, 1);

		if (3 <= *flag&&*flag <= 5) return n;

		if (str2time(buff, 0, 26, time)) {
			trace(2, "rinex obs invalid epoch: epoch=%26.26s\n", buff);
			return 0;
		}
		for (i = 0, j = 32; i < n; i++, j += 3) {
			if (j >= 68) {
				if (!fgets(buff, MAXRNXLEN, fp)) break;
				j = 32;
			}
			if (i < MAXOBS) {
				strncpy(satid, buff + j, 3);
				sats[i] = satid2no(satid);
			}
		}
	}
	else { /* ver.3 */
		if ((n = (int)str2num(buff, 32, 3)) <= 0) return 0;

		*flag = (int)str2num(buff, 31, 1);

		if (3 <= *flag&&*flag <= 5) return n;

		if (buff[0] != '>' || str2time(buff, 1, 28, time))
		{
			trace(2, "rinex obs invalid epoch: epoch=%29.29s\n", buff);
			return 0;
		}
	}
	trace(4, "decode_obsepoch: time=%s flag=%d\n", time_str(*time, 3), *flag);
	return n;
}

/* read rinex nav/gnav/geo nav -----------------------------------------------*/
int readrnxnav(FILE *fp, double ver, int sys, nav_t *nav)
{
	eph_t eph;
	geph_t geph;
	int stat, type;
	trace(3, "readrnxnav: ver=%.2f sys=%d\n", ver, sys);

	if (!nav) return 0;

	/* read rinex navigation data body */
	while ((stat = readrnxnavb(fp, ver, sys, &type, &eph, &geph)) >= 0) {
		/* add ephemeris to navigation data */
		if (stat) {
			switch (type)
			{
			case 1: stat = add_geph(&geph, nav); break;
			default: stat = add_eph(&eph, nav); break;
			}
			if (!stat) return 0;
		}
	}
	return nav->n > 0 || nav->ng > 0 || nav->ns > 0;
}

/* read rinex file -----------------------------------------------------------*/
static int readrnxfp(FILE *fp, char *type, obs_t *obs, nav_t *nav, sta_t *sta)
{
	double ver;
	int sys, tsys = TSYS_GPS;
	char tobs[NUMSYS][MAXOBSTYPE][4] = { {""} };

	//trace(3, "readrnxfp: flag=%d index=%d\n", flag, index);

	/* read rinex header */
	if (!readrnxh(fp, &ver, type, &sys, &tsys, tobs, nav, sta,NULL,NULL)) return 0;

	/* read rinex body */
	switch (*type)
	{
	case 'O': return readrnxobs(fp, ver, &tsys, tobs, obs, sta);
	case 'N': return readrnxnav(fp, ver, sys, nav);
	case 'G': return readrnxnav(fp, ver, _SYS_GLO_, nav);
	}
	trace(2, "unsupported rinex type ver=%.2f type=%c\n", ver, *type);
	return 0;
}

/* uncompress and read rinex file --------------------------------------------*/
static int readrnxfile(const char *file, char *type, obs_t *obs, nav_t *nav, sta_t *sta)
{
	FILE *fp;
	int stat;
	//char tmpfile[1024];

	trace(3, "readrnxfile: file=%s \n", file);

	if (sta) init_sta(sta);

	if (!(fp = fopen(file, "r"))) {
		return 0;
	}
	/* read rinex file */
	stat = readrnxfp(fp, type, obs, nav, sta);

	fclose(fp);

	return stat;
}

/* initialize station parameter ----------------------------------------------*/
static void init_sta(sta_t *sta)
{
	int i;
	*sta->name = '\0';
	*sta->marker = '\0';
	*sta->antdes = '\0';
	*sta->antsno = '\0';
	*sta->rectype = '\0';
	*sta->recver = '\0';
	*sta->recsno = '\0';
	sta->antsetup = sta->itrf = sta->deltype = 0;
	for (i = 0; i < 3; i++) sta->pos[i] = 0.0;
	for (i = 0; i < 3; i++) sta->del[i] = 0.0;
	sta->hgt = 0.0;
}

/* set string without tail space ---------------------------------------------*/
static void setstr(char *dst, const char *src, int n)
{
	char *p = dst;
	const char *q = src;
	while (*q&&q < src + n) *p++ = *q++;
	*p-- = '\0';
	while (p >= dst && *p == ' ') *p-- = '\0';
}

/* read rinex format obs and nav data  -----------------------------------------------------*/
extern int readobsnav(const char **infile, int n, obs_t *obs, nav_t *nav, sta_t *sta)
{
	int i, ind = 0, nobs = 0;
	int nepoch;

	for (i = 0; i < n; i++)
	{
		/* read rinex obs and nav file */
		if (readrnxt(infile[i], obs, nav, sta) < 0)
		{
			trace(1, "insufficient memory\n");
			return 0;
		}
	}

	if (obs->n <= 0) {
		trace(1, "\n");
		return 0;
	}
	if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
		trace(1, "\n");
		return 0;
	}
	/* sort observation data */
	nepoch = sortobs(obs);

	/* delete duplicated ephemeris */
	uniqnav(nav);

	return 1;
}

/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
*-----------------------------------------------------------------------------*/
extern double str2num(const char *s, int i, int n)
{
	double value;
	char str[256], *p = str;

	if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < n) return 0.0;
	for (s += i; *s&&--n >= 0; s++) *p++ = *s == 'd' || *s == 'D' ? 'E' : *s;
	*p = '\0';
	return sscanf(str, "%lf", &value) == 1 ? value : 0.0;
}
/* convert rinex obs type ver.2 -> ver.3 -------------------------------------*/
static void convcode(double ver, int sys, const char *str, char *type)
{
	strcpy(type, "   ");

	if (!strcmp(str, "P1")) { /* ver.2.11 GPS L1PY,GLO L2P */
		if (sys == _SYS_GPS_) sprintf(type, "%c1W", 'C');
		else if (sys == _SYS_GLO_) sprintf(type, "%c1P", 'C');
	}
	else if (!strcmp(str, "P2")) { /* ver.2.11 GPS L2PY,GLO L2P */
		if (sys == _SYS_GPS_) sprintf(type, "%c2W", 'C');
		else if (sys == _SYS_GLO_) sprintf(type, "%c2P", 'C');
	}
	else if (!strcmp(str, "C1")) { /* ver.2.11 GPS L1C,GLO L1C/A */
		if (ver >= 2.12); /* reject C1 for 2.12 */
		else if (sys == _SYS_GPS_) sprintf(type, "%c1C", 'C');
		else if (sys == _SYS_GLO_) sprintf(type, "%c1C", 'C');
		else if (sys == _SYS_GAL_) sprintf(type, "%c1X", 'C'); /* ver.2.12 */
		else if (sys == _SYS_QZS_) sprintf(type, "%c1C", 'C');
		else if (sys == _SYS_SBS_) sprintf(type, "%c1C", 'C');
	}
	else if (!strcmp(str, "C2")) {
		if (sys == _SYS_GPS_) {
			if (ver >= 2.12) sprintf(type, "%c2W", 'C'); /* L2P(Y) */
			else           sprintf(type, "%c2X", 'C'); /* L2C */
		}
		else if (sys == _SYS_GLO_) sprintf(type, "%c2C", 'C');
		else if (sys == _SYS_QZS_) sprintf(type, "%c2X", 'C');
		else if (sys == _SYS_BDS_) sprintf(type, "%c1X", 'C'); /* ver.2.12 B1 */
	}
	else if (ver >= 2.12&&str[1] == 'A') { /* ver.2.12 L1C/A */
		if (sys == _SYS_GPS_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_GLO_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_SBS_) sprintf(type, "%c1C", str[0]);
	}
	else if (ver >= 2.12&&str[1] == 'B') { /* ver.2.12 GPS L1C */
		if (sys == _SYS_GPS_) sprintf(type, "%c1X", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c1X", str[0]);
	}
	else if (ver >= 2.12&&str[1] == 'C') { /* ver.2.12 GPS L2C */
		if (sys == _SYS_GPS_) sprintf(type, "%c2X", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c2X", str[0]);
	}
	else if (ver >= 2.12&&str[1] == 'D') { /* ver.2.12 GLO L2C/A */
		if (sys == _SYS_GLO_) sprintf(type, "%c2C", str[0]);
	}
	else if (ver >= 2.12&&str[1] == '1') { /* ver.2.12 GPS L1PY,GLO L1P */
		if (sys == _SYS_GPS_) sprintf(type, "%c1W", str[0]);
		else if (sys == _SYS_GLO_) sprintf(type, "%c1P", str[0]);
		else if (sys == _SYS_GAL_) sprintf(type, "%c1X", str[0]); /* tentative */
		else if (sys == _SYS_BDS_) sprintf(type, "%c1X", str[0]); /* extension */
	}
	else if (ver < 2.12&&str[1] == '1') {
		if (sys == _SYS_GPS_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_GLO_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_GAL_) sprintf(type, "%c1X", str[0]); /* tentative */
		else if (sys == _SYS_QZS_) sprintf(type, "%c1C", str[0]);
		else if (sys == _SYS_SBS_) sprintf(type, "%c1C", str[0]);
	}
	else if (str[1] == '2') {
		if (sys == _SYS_GPS_) sprintf(type, "%c2W", str[0]);
		else if (sys == _SYS_GLO_) sprintf(type, "%c2P", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c2X", str[0]);
		else if (sys == _SYS_BDS_) sprintf(type, "%c1X", str[0]); /* ver.2.12 B1 */
	}
	else if (str[1] == '5') {
		if (sys == _SYS_GPS_) sprintf(type, "%c5X", str[0]);
		else if (sys == _SYS_GAL_) sprintf(type, "%c5X", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c5X", str[0]);
		else if (sys == _SYS_SBS_) sprintf(type, "%c5X", str[0]);
	}
	else if (str[1] == '6') {
		if (sys == _SYS_GAL_) sprintf(type, "%c6X", str[0]);
		else if (sys == _SYS_QZS_) sprintf(type, "%c6X", str[0]);
		else if (sys == _SYS_BDS_) sprintf(type, "%c6X", str[0]); /* ver.2.12 B3 */
	}
	else if (str[1] == '7') {
		if (sys == _SYS_GAL_) sprintf(type, "%c7X", str[0]);
		else if (sys == _SYS_BDS_) sprintf(type, "%c7X", str[0]); /* ver.2.12 B2 */
	}
	else if (str[1] == '8') {
		if (sys == _SYS_GAL_) sprintf(type, "%c8X", str[0]);
	}
	trace(3, "convcode: ver=%.2f sys=%2d type= %s -> %s\n", ver, sys, str, type);
}

/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number
* args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn,Inn or Snn)
* return : satellite number (0: error)
* notes  : 120-142 and 193-199 are also recognized as sbas and qzss
*-----------------------------------------------------------------------------*/
extern int satid2no(const char *id)
{
	int sys, prn;
	char code;

	if (sscanf(id, "%d", &prn) == 1) {
		if (MINPRNGPS <= prn && prn <= MAXPRNGPS) sys = _SYS_GPS_;
		else if (MINPRNSBS <= prn && prn <= MAXPRNSBS) sys = _SYS_SBS_;
		else if (MINPRNQZS <= prn && prn <= MAXPRNQZS) sys = _SYS_QZS_;
		else return 0;
		return satno(sys, prn);
	}
	if (sscanf(id, "%c%d", &code, &prn) < 2) return 0;

	switch (code) {
	case 'G': sys = _SYS_GPS_; prn += MINPRNGPS - 1; break;
	case 'R': sys = _SYS_GLO_; prn += MINPRNGLO - 1; break;
	case 'E': sys = _SYS_GAL_; prn += MINPRNGAL - 1; break;
	case 'J': sys = _SYS_QZS_; prn += MINPRNQZS - 1; break;
	case 'C': sys = _SYS_BDS_; prn += MINPRNCMP - 1; break;
	case 'I': sys = _SYS_IRN_; prn += MINPRNIRN - 1; break;
	case 'L': sys = _SYS_LEO_; prn += MINPRNLEO - 1; break;
	case 'S': sys = _SYS_SBS_; prn += 100; break;
	default: return 0;
	}
	return satno(sys, prn);
}
/* decode obs data -----------------------------------------------------------*/
static int decode_obsdata(FILE *fp, char *buff, double ver, sigind_t *index, obsd_t *obs)
{
	sigind_t *ind;
	double val[MAXOBSTYPE] = { 0 };
	unsigned char lli[MAXOBSTYPE] = { 0 };
	char satid[8] = "";
	int i, j, n, m, stat = 1, p[MAXOBSTYPE], k[16], l[16];

	//trace(4, "decode_obsdata: ver=%.2f\n", ver);

	if (ver > 2.99) { /* ver.3 */
		strncpy(satid, buff, 3);
		obs->sat = (unsigned char)satid2no(satid);
	}
	if (!obs->sat) {
		trace(4, "decode_obsdata: unsupported sat sat=%s\n", satid);
		stat = 0;
	}
	else if (!(satsys(obs->sat, NULL))) {
		stat = 0;
	}
	/* read obs data fields */
	switch (satsys(obs->sat, NULL)) {
	case _SYS_GLO_: ind = index + 1; break;
	case _SYS_GAL_: ind = index + 2; break;
	case _SYS_QZS_: ind = index + 3; break;
	case _SYS_SBS_: ind = index + 4; break;
	case _SYS_BDS_: ind = index + 5; break;
	default:      ind = index; break;
	}
	for (i = 0, j = ver <= 2.99 ? 0 : 3; i < ind->n; i++, j += 16) {

		if (ver <= 2.99&&j >= 80) { /* ver.2 */
			if (!fgets(buff, MAXRNXLEN, fp)) break;
			j = 0;
		}
		if (stat) {
			val[i] = str2num(buff, j, 14) + ind->shift[i];
			lli[i] = (unsigned char)str2num(buff, j + 14, 1) & 3;
		}
	}
	if (!stat) return 0;

	for (i = 0; i < NFREQ + NEXOBS; i++) {
		obs->P[i] = obs->L[i] = 0.0; obs->D[i] = 0.0f;
		obs->SNR[i] = obs->LLI[i] = obs->code[i] = 0;
	}
	/* assign position in obs data */
	for (i = n = m = 0; i < ind->n; i++) {

		p[i] = ver <= 2.11 ? ind->frq[i] - 1 : ind->pos[i];

		if (ind->type[i] == 0 && p[i] == 0) k[n++] = i; /* C1? index */
		if (ind->type[i] == 0 && p[i] == 1) l[m++] = i; /* C2? index */
	}
	if (ver <= 2.11) {

		/* if multiple codes (C1/P1,C2/P2), select higher priority */
		if (n >= 2) {
			if (val[k[0]] == 0.0&&val[k[1]] == 0.0) {
				p[k[0]] = -1; p[k[1]] = -1;
			}
			else if (val[k[0]] != 0.0&&val[k[1]] == 0.0) {
				p[k[0]] = 0; p[k[1]] = -1;
			}
			else if (val[k[0]] == 0.0&&val[k[1]] != 0.0) {
				p[k[0]] = -1; p[k[1]] = 0;
			}
			else if (ind->pri[k[1]] > ind->pri[k[0]]) {
				p[k[1]] = 0; p[k[0]] = NEXOBS < 1 ? -1 : NFREQ;
			}
			else {
				p[k[0]] = 0; p[k[1]] = NEXOBS < 1 ? -1 : NFREQ;
			}
		}
		if (m >= 2) {
			if (val[l[0]] == 0.0&&val[l[1]] == 0.0) {
				p[l[0]] = -1; p[l[1]] = -1;
			}
			else if (val[l[0]] != 0.0&&val[l[1]] == 0.0) {
				p[l[0]] = 1; p[l[1]] = -1;
			}
			else if (val[l[0]] == 0.0&&val[l[1]] != 0.0) {
				p[l[0]] = -1; p[l[1]] = 1;
			}
			else if (ind->pri[l[1]] > ind->pri[l[0]]) {
				p[l[1]] = 1; p[l[0]] = NEXOBS < 2 ? -1 : NFREQ + 1;
			}
			else {
				p[l[0]] = 1; p[l[1]] = NEXOBS < 2 ? -1 : NFREQ + 1;
			}
		}
	}
	/* save obs data */
	for (i = 0; i < ind->n; i++) {
		if (p[i] < 0 || val[i] == 0.0) continue;
		switch (ind->type[i]) {
		case 0: obs->P[p[i]] = val[i]; obs->code[p[i]] = ind->code[i]; break;
		case 1: obs->L[p[i]] = val[i]; obs->LLI[p[i]] = lli[i];       break;
		case 2: obs->D[p[i]] = (float)val[i];                        break;
		case 3: obs->SNR[p[i]] = (unsigned char)(val[i] * 4.0 + 0.5);    break;
		}
	}
	//trace(4, "decode_obsdata: time=%s sat=%2d\n", time_str(obs->time, 0), obs->sat);
	return 1;
}

/* decode glonass ephemeris --------------------------------------------------*/
static int decode_geph(double ver, int sat, gtime_t toc, double *data,
	geph_t *geph)
{
	geph_t geph0 = { 0 };
	gtime_t tof;
	double tow, tod;
	int week, dow;

	trace(4, "decode_geph: ver=%.2f sat=%2d\n", ver, sat);

	if (satsys(sat, NULL) != _SYS_GLO_) {
		trace(3, "glonass ephemeris error: invalid satellite sat=%2d\n", sat);
		return 0;
	}
	*geph = geph0;

	geph->sat = sat;

	/* toc rounded by 15 min in utc */
	tow = time2gpst(toc, &week);
	toc = gpst2time(week, floor((tow + 450.0) / 900.0) * 900);
	dow = (int)floor(tow / 86400.0);

	/* time of frame in utc */
	tod = ver <= 2.99 ? data[2] : fmod(data[2], 86400.0); /* tod (v.2), tow (v.3) in utc */
	tof = gpst2time(week, tod + dow * 86400.0);
	tof = adjday(tof, toc);

	geph->toe = utc2gpst(toc);   /* toc (gpst) */
	geph->tof = utc2gpst(tof);   /* tof (gpst) */

	/* iode = tb (7bit), tb =index of UTC+3H within current day */
	geph->iode = (int)(fmod(tow + 10800.0, 86400.0) / 900.0 + 0.5);

	geph->taun = -data[0];       /* -taun */
	geph->gamn = data[1];       /* +gamman */

	geph->pos[0] = data[3] * 1E3; geph->pos[1] = data[7] * 1E3; geph->pos[2] = data[11] * 1E3;
	geph->vel[0] = data[4] * 1E3; geph->vel[1] = data[8] * 1E3; geph->vel[2] = data[12] * 1E3;
	geph->acc[0] = data[5] * 1E3; geph->acc[1] = data[9] * 1E3; geph->acc[2] = data[13] * 1E3;

	geph->svh = (int)data[6];
	geph->frq = (int)data[10];
	geph->age = (int)data[14];

	/* some receiver output >128 for minus frequency number */
	if (geph->frq > 128) geph->frq -= 256;

	//if (geph->frq < MINFREQ_GLO || MAXFREQ_GLO < geph->frq) {
	//	trace(2, "rinex gnav invalid freq: sat=%2d fn=%d\n", sat, geph->frq);
	//}
	return 1;
}

/* adjust time considering week handover -------------------------------------*/
static gtime_t adjday(gtime_t t, gtime_t t0)
{
	double tt = timediff(t, t0);
	if (tt < -43200.0) return timeadd(t, 86400.0);
	if (tt > 43200.0) return timeadd(t, -86400.0);
	return t;
}

/* decode ephemeris ----------------------------------------------------------*/
static int decode_eph(double ver, int sat, gtime_t toc, const double *data, eph_t *eph)
{
	eph_t eph0 = { 0 };
	int sys;

	trace(4, "decode_eph: ver=%.2f sat=%2d\n", ver, sat);

	sys = satsys(sat, NULL);

	if (!(sys&(_SYS_GPS_ | _SYS_GAL_ | _SYS_QZS_ | _SYS_BDS_ | _SYS_IRN_))) {
		trace(3, "ephemeris error: invalid satellite sat=%2d\n", sat);
		return 0;
	}
	*eph = eph0;

	eph->sat = sat;
	eph->toc = toc;

	eph->f0 = data[0];
	eph->f1 = data[1];
	eph->f2 = data[2];

	eph->A = SQR(data[10]); eph->e = data[8]; eph->i0 = data[15]; eph->OMG0 = data[13];
	eph->omg = data[17]; eph->M0 = data[6]; eph->deln = data[5]; eph->OMGd = data[18];
	eph->idot = data[19]; eph->crc = data[16]; eph->crs = data[4]; eph->cuc = data[7];
	eph->cus = data[9]; eph->cic = data[12]; eph->cis = data[14];

	if (sys == _SYS_GPS_ || sys == _SYS_QZS_) {
		eph->iode = (int)data[3];      /* IODE */
		eph->iodc = (int)data[26];      /* IODC */
		eph->toes = data[11];      /* toe (s) in gps week */
		eph->week = (int)data[21];      /* gps week */
		eph->toe = adjweek_(gpst2time(eph->week, data[11]), toc);
		eph->ttr = adjweek_(gpst2time(eph->week, data[27]), toc);

		eph->code = (int)data[20];      /* GPS: codes on L2 ch */
		eph->svh = (int)data[24];      /* sv health */
		eph->sva = uraindex(data[23]);  /* ura (m->index) */
		eph->flag = (int)data[22];      /* GPS: L2 P data flag */

		eph->tgd[0] = data[25];      /* TGD */
		if (sys == _SYS_GPS_) {
			eph->fit = data[28];        /* fit interval (h) */
		}
		else {
			eph->fit = data[28] == 0.0 ? 1.0 : 2.0; /* fit interval (0:1h,1:>2h) */
		}
	}
	else if (sys == _SYS_GAL_) { /* GAL ver.3 */
		eph->iode = (int)data[3];      /* IODnav */
		eph->toes = data[11];      /* toe (s) in galileo week */
		eph->week = (int)data[21];      /* gal week = gps week */
		eph->toe = adjweek_(gpst2time(eph->week, data[11]), toc);
		eph->ttr = adjweek_(gpst2time(eph->week, data[27]), toc);

		eph->code = (int)data[20];      /* data sources */
									  /* bit 0 set: I/NAV E1-B */
									  /* bit 1 set: F/NAV E5a-I */
									  /* bit 2 set: F/NAV E5b-I */
									  /* bit 8 set: af0-af2 toc are for E5a.E1 */
									  /* bit 9 set: af0-af2 toc are for E5b.E1 */
		eph->svh = (int)data[24];      /* sv health */
									  /* bit     0: E1B DVS */
									  /* bit   1-2: E1B HS */
									  /* bit     3: E5a DVS */
									  /* bit   4-5: E5a HS */
									  /* bit     6: E5b DVS */
									  /* bit   7-8: E5b HS */
		eph->sva = uraindex(data[23]); /* ura (m->index) */

		eph->tgd[0] = data[25];      /* BGD E5a/E1 */
		eph->tgd[1] = data[26];      /* BGD E5b/E1 */
	}
	else if (sys == _SYS_BDS_) { /* BeiDou v.3.02 */
		eph->toc = bdt2gpst(eph->toc);  /* bdt -> gpst */
		eph->iode = (int)data[3];      /* AODE */
		eph->iodc = (int)data[28];      /* AODC */
		eph->toes = data[11];      /* toe (s) in bdt week */
		eph->week = (int)data[21];      /* bdt week */
		eph->toe = bdt2gpst(bdt2time(eph->week, data[11])); /* bdt -> gpst */
		eph->ttr = bdt2gpst(bdt2time(eph->week, data[27])); /* bdt -> gpst */
		eph->toe = adjweek_(eph->toe, toc);
		eph->ttr = adjweek_(eph->ttr, toc);

		eph->svh = (int)data[24];      /* satH1 */
		eph->sva = uraindex(data[23]);  /* ura (m->index) */

		eph->tgd[0] = data[25];      /* TGD1 B1/B3 */
		eph->tgd[1] = data[26];      /* TGD2 B2/B3 */
	}
	else if (sys == _SYS_IRN_) { /* IRNSS v.3.03 */
		eph->iode = (int)data[3];      /* IODEC */
		eph->toes = data[11];      /* toe (s) in irnss week */
		eph->week = (int)data[21];      /* irnss week */
		eph->toe = adjweek_(gpst2time(eph->week, data[11]), toc);
		eph->ttr = adjweek_(gpst2time(eph->week, data[27]), toc);
		eph->svh = (int)data[24];      /* sv health */
		eph->sva = uraindex(data[23]);  /* ura (m->index) */
		eph->tgd[0] = data[25];      /* TGD */
	}
	if (eph->iode < 0 || 1023 < eph->iode) {
		trace(2, "rinex nav invalid: sat=%2d iode=%d\n", sat, eph->iode);
	}
	if (eph->iodc < 0 || 1023 < eph->iodc) {
		trace(2, "rinex nav invalid: sat=%2d iodc=%d\n", sat, eph->iodc);
	}
	return 1;
}

/* expand file path ------------------------------------------------------------
* expand file path with wild-card (*) in file
* args   : char   *path     I   file path to expand (captal insensitive)
*          char   *paths    O   expanded file paths
*          int    nmax      I   max number of expanded file paths
* return : number of expanded file paths
* notes  : the order of expanded files is alphabetical order
*-----------------------------------------------------------------------------*/
extern int expath(const char *path, char *paths[], int nmax)
{
	int i, j, n = 0;
	char tmp[1024];
#ifdef WIN32
	WIN32_FIND_DATA file;
	HANDLE h;
	char dir[1024] = "", *p;

	trace(3, "expath  : path=%s nmax=%d\n", path, nmax);

	if ((p = strrchr(path, '\\'))) {
		strncpy(dir, path, p - path + 1); dir[p - path + 1] = '\0';
	}
	if ((h = FindFirstFile((LPCTSTR)path, &file)) == INVALID_HANDLE_VALUE) {
		strcpy(paths[0], path);
		return 1;
	}
	sprintf(paths[n++], "%s%s", dir, file.cFileName);
	while (FindNextFile(h, &file) && n < nmax) {
		if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
		sprintf(paths[n++], "%s%s", dir, file.cFileName);
	}
	FindClose(h);
#else
	struct dirent *d;
	DIR *dp;
	const char *file = path;
	char dir[1024] = "", s1[1024], s2[1024], *p, *q, *r;

	trace(3, "expath  : path=%s nmax=%d\n", path, nmax);

	if ((p = strrchr(path, '/')) || (p = strrchr(path, '\\'))) {
		file = p + 1; strncpy(dir, path, p - path + 1); dir[p - path + 1] = '\0';
	}
	if (!(dp = opendir(*dir ? dir : "."))) return 0;
	while ((d = readdir(dp))) {
		if (*(d->d_name) == '.') continue;
		sprintf(s1, "^%s$", d->d_name);
		sprintf(s2, "^%s$", file);
		for (p = s1; *p; p++) *p = (char)tolower((int)*p);
		for (p = s2; *p; p++) *p = (char)tolower((int)*p);

		for (p = s1, q = strtok_r(s2, "*", &r); q; q = strtok_r(NULL, "*", &r)) {
			if ((p = strstr(p, q))) p += strlen(q); else break;
		}
		if (p&&n < nmax) sprintf(paths[n++], "%s%s", dir, d->d_name);
	}
	closedir(dp);
#endif
	/* sort paths in alphabetical order */
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			if (strcmp(paths[i], paths[j]) > 0) {
				strcpy(tmp, paths[i]);
				strcpy(paths[i], paths[j]);
				strcpy(paths[j], tmp);
			}
		}
	}
	for (i = 0; i < n; i++) trace(3, "expath  : file=%s\n", paths[i]);

	return n;
}

/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t *obs    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
	//int i, j, n;

	//trace(3, "sortobs: nobs=%d\n", obs->n);

	//if (obs->n <= 0) return 0;

	//qsort(obs->data, obs->n, sizeof(obsd_t), cmpobs);

	///* delete duplicated data */
	//for (i = j = 0; i < obs->n; i++) {
	//	if (obs->data[i].sat != obs->data[j].sat ||
	//		obs->data[i].rcv != obs->data[j].rcv ||
	//		timediff(obs->data[i].time, obs->data[j].time) != 0.0) {
	//		obs->data[++j] = obs->data[i];
	//	}
	//}
	//obs->n = j + 1;

	//for (i = n = 0; i < obs->n; i = j, n++) {
	//	for (j = i + 1; j < obs->n; j++) {
	//		if (timediff(obs->data[j].time, obs->data[i].time) > DTTOL) break;
	//	}
	//}
	//return n;
	return 0;
}

/* unique ephemerides ----------------------------------------------------------
* unique ephemerides in navigation data and update carrier wave length
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
	trace(3, "uniqnav: neph=%d ngeph=%d nseph=%d\n", nav->n, nav->ng, nav->ns);

	/* unique ephemeris */
	//uniqeph(nav);
	//uniqgeph(nav);
	//uniqseph(nav);

	/* update carrier wave length */
	//for (i = 0; i < MAXSAT; i++) for (j = 0; j < NFREQ; j++) {
	//	nav->lam[i][j] = satwavelen(i + 1, j, nav);
	//}
}
static const double ura_eph[] = {         /* ura values (ref [3] 20.3.3.3.1.1) */
	2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
	3072.0,6144.0,0.0
};
/* ura value (m) to ura index ------------------------------------------------*/
static int uraindex(double value)
{
	int i;
	for (i = 0; i < 15; i++) if (ura_eph[i] >= value) break;
	return i;
}
extern int readrnxt(const char *file, obs_t *obs, nav_t *nav, sta_t *sta)
{
	int i, n, stat = 0;
//	const char *p;
	char type = ' ', *files[MAXEXFILE] = { 0 };

	if (!*file)
	{
		return readrnxfp(stdin, &type, obs, nav, sta);
	}
	for (i = 0; i < MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return -1;
		}
	}
	// expand wild-card 
	if ((n = expath(file, files, MAXEXFILE)) <= 0) {
		for (i = 0; i < MAXEXFILE; i++) free(files[i]);
		return 0;
	}
	/* read rinex files */
	for (i = 0; i < n&&stat >= 0; i++) {
		stat = readrnxfile(files[i], &type, obs, nav, sta);
	}

	for (i = 0; i < MAXEXFILE; i++) free(files[i]);

	return stat;
}


/* read rinex navigation data body -------------------------------------------*/
int readrnxnavb(FILE *fp, double ver, int sys,
	int *type, eph_t *eph, geph_t *geph)
{
	gtime_t toc;
	double data[64];
	int i = 0, j, prn, sat = 0, sp = 3;
	char buff[MAXRNXLEN] = {0}, id[8] = "", *p;

	trace(4, "readrnxnavb: ver=%.2f sys=%d\n", ver, sys);


	while (fgets(buff, MAXRNXLEN, fp)) {

		if (i == 0) {

			/* decode satellite field */
			if (ver >= 3.0 || sys == _SYS_GAL_ || sys == _SYS_QZS_) { /* ver.3 or GAL/QZS */
				strncpy(id, buff, 3);
				sat = satid2no(id);
				sp = 4;
				if (ver >= 3.0) sys = satsys(sat, NULL);
			}
			else {
				prn = (int)str2num(buff, 0, 2);

				if (sys == _SYS_GLO_)
				{
					sat = satno(_SYS_GLO_, prn);
				}
				else if (93 <= prn && prn <= 97) { /* extension */
					sat = satno(_SYS_QZS_, prn + 100);
				}
				else sat = satno(_SYS_GPS_, prn);
			}
			/* decode toc field */
			if (str2time(buff + sp, 0, 19, &toc))
			{
				trace(2, "rinex nav toc error: %23.23s\n", buff);
				return 0;
			}
			/* decode data fields */
			for (j = 0, p = buff + sp + 19; j < 3; j++, p += 19) {
				data[i++] = str2num(p, 0, 19);
			}
		}
		else {
			/* decode data fields */
			for (j = 0, p = buff + sp; j < 4; j++, p += 19) {
				data[i++] = str2num(p, 0, 19);
			}
			/* decode ephemeris */
			if (sys == _SYS_GLO_ && i >= 15) {
				if (!(sys)) return 0;
				*type = 1;
				return decode_geph(ver, sat, toc, data, geph);
			}
			else if (i >= 31) {
				if (!(sys)) return 0;
				*type = 0;
				return decode_eph(ver, sat, toc, data, eph);
			}
		}
	}
	return -1;
}

