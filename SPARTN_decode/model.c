
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "model.h"
#include "rtcm.h"
#include "gnss_math.h"

#ifndef RAD2DEG
#define RAD2DEG (180.0/PI)
#endif


void float_to_char(double data, int32_t totalNum, uint8_t decimalsNum, char *buf);
char *itoa_user(int num, char *str, int radix);
void RealToArray(double real, char *a, unsigned char id, unsigned char dd); 

/* time to day of year ---------------------------------------------------------
* convert time to day of year
* args   : gtime_t t        I   gtime_t struct
* return : day of year (days)
*-----------------------------------------------------------------------------*/
extern double time2doy(gtime_t t)
{
    double ep[6];

    time2epoch(t, ep);
    ep[1] = ep[2] = 1.0; ep[3] = ep[4] = ep[5] = 0.0;
    return timediff(t, epoch2time(ep)) / 86400.0 + 1.0;
}

extern void blh2C_en(const double *blh, double C_en[3][3])
{
	/* blh => C_en */
	double lat = blh[0], lon = blh[1]; /*, ht = blh[2];*/
	C_en[0][0] = -sin(lat) * cos(lon);
	C_en[1][0] = -sin(lat) * sin(lon);
	C_en[2][0] = cos(lat);
	C_en[0][1] = -sin(lon);
	C_en[1][1] = cos(lon);
	C_en[2][1] = 0.0;
	C_en[0][2] = -cos(lat) * cos(lon);
	C_en[1][2] = -cos(lat) * sin(lon);
	C_en[2][2] = -sin(lat);
	return;
}

extern void xyz2ned(double C_en[3][3], double *xyz, double *covXYZ, double *ned, double *covNED)
{
	double temp[3][3] = {0};
	int i, j, k;
	ned[0] = C_en[0][0] * xyz[0] + C_en[1][0] * xyz[1] + C_en[2][0] * xyz[2];
	ned[1] = C_en[0][1] * xyz[0] + C_en[1][1] * xyz[1] + C_en[2][1] * xyz[2];
	ned[2] = C_en[0][2] * xyz[0] + C_en[1][2] * xyz[1] + C_en[2][2] * xyz[2];
	if (covXYZ != NULL && covNED != NULL)
	{
		/* covNED = C_en'*covXYZ*C_en */
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
			{
				temp[i][j] = 0.0;
				for (k = 0; k < 3; ++k)
				{
					temp[i][j] += C_en[k][i] * covXYZ[SMI(k, j)];
				}
			}
		}

		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 3; ++j)
			{
				covNED[SMI(i, j)] = 0.0;
				for (k = 0; k < 3; ++k)
				{
					covNED[SMI(i, j)] += temp[i][k] * C_en[k][j];
				}
			}
		}
	}
	return;
}

extern void blhdiff(double *blh, double *blh_ref, double *ned)
{
	double C_en[3][3] = {0};
	double xyz[3] = {0};
	double xyz_ref[3] = {0};
	double dxyz[3] = {0};
	blh2C_en(blh_ref, C_en);
	pos2ecef(blh_ref, xyz_ref);
	pos2ecef(blh, xyz);
	dxyz[0] = xyz[0] - xyz_ref[0];
	dxyz[1] = xyz[1] - xyz_ref[1];
	dxyz[2] = xyz[2] - xyz_ref[2];
	ned[0] = C_en[0][0] * dxyz[0] + C_en[1][0] * dxyz[1] + C_en[2][0] * dxyz[2];
	ned[1] = C_en[0][1] * dxyz[0] + C_en[1][1] * dxyz[1] + C_en[2][1] * dxyz[2];
	ned[2] = C_en[0][2] * dxyz[0] + C_en[1][2] * dxyz[1] + C_en[2][2] * dxyz[2];
	return;
}

/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void xyz2enu_(const double *pos, double *E)
{
	double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);

	E[0] = -sinl;
	E[3] = cosl;
	E[6] = 0.0;
	E[1] = -sinp * cosl;
	E[4] = -sinp * sinl;
	E[7] = cosp;
	E[2] = cosp * cosl;
	E[5] = cosp * sinl;
	E[8] = sinp;
}
/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2enu(const double *pos, const double *r, double *e)
{
	double E[9];

	xyz2enu_(pos, E);
	matmul("NN", 3, 1, 3, 1.0, E, r, 0.0, e);
}

/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void xyz2enu(const double *pos, double *E)
{
	double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);

	E[0] = -sinl;
	E[3] = cosl;
	E[6] = 0.0;
	E[1] = -sinp * cosl;
	E[4] = -sinp * sinl;
	E[7] = cosp;
	E[2] = cosp * cosl;
	E[5] = cosp * sinl;
	E[8] = sinp;
}

/* transform local enu coordinate covariance to xyz-ecef -----------------------
* transform local enu covariance to xyz-ecef coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *Q        I   covariance in local enu coordinate
*          double *P        O   covariance in xyz-ecef coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covecef(const double *pos, const double *Q, double *P)
{
	double E[9], EQ[9];

	xyz2enu(pos, E);
	matmul("TN", 3, 3, 3, 1.0, E, Q, 0.0, EQ);
	matmul("NN", 3, 3, 3, 1.0, EQ, E, 0.0, P);
}

extern double satazel(const double *pos, const double *e, double *azel)
{
	double az = 0.0, el = PI / 2.0, enu[3];

	ecef2enu(pos, e, enu);
	az = dot(enu, enu, 2) < 1E-12 ? 0.0 : atan2(enu[0], enu[1]);
	if (az < 0.0)
		az += 2 * PI;
	el = asin(enu[2]);

	azel[0] = az;
	azel[1] = el;

	return el;
}

extern double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	for (i = 0; i < 3; i++)
		e[i] = rs[i] - rr[i];
	r = norm(e, 3);
	for (i = 0; i < 3; i++)
		e[i] /= r;
	return r + OMGE * (rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT;
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position and vel(ecef at transmission) (m)
*          double *rr       I   receiver velocity (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geovel(const double *rs, const double *rr, double *e)
{
	double rate = 0.0;
	int i;
	double vs[3] = {0.0};

	for (i = 0; i < 3; i++)
	{
		vs[i] = rs[3 + i] - rr[3 + i];
	}
	rate = dot(vs, e, 3);
	rate += OMGE / CLIGHT * (rs[4] * rr[0] + rs[1] * rr[3] - rs[3] * rr[1] - rs[0] * rr[4]);
	return rate;
}

extern double tropmodel(const double *blh, const double *azel, double humi)
{
	const double temp0 = 15.0; /* temparature at sea level */
	double hgt, pres, temp, e, z, trph, trpw;

	if (blh[2] < -100.0 || 1E4 < blh[2] || azel[1] <= 0)
		return 0.0;

	/* standard atmosphere */
	hgt = blh[2] < 0.0 ? 0.0 : blh[2];

	pres = 1013.25 * pow(1.0 - 2.2557E-5 * hgt, 5.2568);
	temp = temp0 - 6.5E-3 * hgt + 273.16;
	e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));

	/* saastamoninen model */
	z = PI / 2.0 - azel[1];
	trph = 0.0022768 * pres / (1.0 - 0.00266 * cos(2.0 * blh[0]) - 0.00028 * hgt / 1E3) / cos(z);
	trpw = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);
	return trph + trpw;
}

double interpc(const double coef[], double lat)
{
    int i = (int)(lat / 15.0);
    if (i < 1) return coef[0]; else if (i > 4) return coef[4];
    return coef[i - 1] * (1.0 - lat / 15.0 + i) + coef[i] * (lat / 15.0 - i);
}

double mapf(double el, double a, double b, double c)
{
    double sinel = sin(el);
    return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
}

double nmf(gtime_t time, const double pos[], const double azel[],double *mapfw)
{
    /* ref [5] table 3 */
    /* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
    const double coef[][5] = {
        { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
        { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
        { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

        { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
        { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
        { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

        { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
        { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
        { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
    };
    const double aht[] = { 2.53E-5, 5.49E-3, 1.14E-3 }; /* height correction */

    double y, cosy, ah[3], aw[3], dm, el = azel[1], lat = pos[0] * R2D, hgt = pos[2];
    int i;

    if (el <= 0.0) {
        if (mapfw) *mapfw = 0.0;
        return 0.0;
    }
    /* year from doy 28, added half a year for southern latitudes */
    y = (time2doy(time) - 28.0) / 365.25 + (lat < 0.0 ? 0.5 : 0.0);

    cosy = cos(2.0*PI*y);
    lat = fabs(lat);

    for (i = 0; i < 3; i++) {
        ah[i] = interpc(coef[i], lat) - interpc(coef[i + 3], lat)*cosy;
        aw[i] = interpc(coef[i + 6], lat);
    }
    /* ellipsoidal height is used instead of height above sea level */
    dm = (1.0 / sin(el) - mapf(el, aht[0], aht[1], aht[2]))*hgt / 1E3;

    if (mapfw) *mapfw = mapf(el, aw[0], aw[1], aw[2]);

    return mapf(el, ah[0], ah[1], ah[2]) + dm;
}

/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref [5] (NMF) and [9] (GMF)
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
double tropmapf(gtime_t time, const double pos[], const double azel[],
    double *mapfw)
{
#ifdef IERS_MODEL
    const double ep[] = { 2000,1,1,12,0,0 };
    double mjd, lat, lon, hgt, zd, gmfh, gmfw;
#endif

    if (pos[2]<-1000.0 || pos[2]>20000.0) {
        if (mapfw) *mapfw = 0.0;
        return 0.0;
    }
#ifdef IERS_MODEL
    mjd = 51544.5 + (timediff(time, epoch2time(ep))) / 86400.0;
    lat = pos[0];
    lon = pos[1];
    hgt = pos[2] - geoidh(pos); /* height in m (mean sea level) */
    zd = PI / 2.0 - azel[1];

    /* call GMF */
    gmf_(&mjd, &lat, &lon, &hgt, &zd, &gmfh, &gmfw);

    if (mapfw) *mapfw = gmfw;
    return gmfh;
#else
    return nmf(time, pos, azel, mapfw); /* NMF */
#endif
}

/*----------------------------------------------------------------------------*
/* output GGA given xyz */
/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args   : double deg       I   degree
*          double *dms      O   degree-minute-second {deg,min,sec}
*          int    ndec      I   number of decimals of second
* return : none
*-----------------------------------------------------------------------------*/
void deg2dms(double deg, double *dms, int ndec)
{
	double sign = deg < 0.0 ? -1.0 : 1.0, a = fabs(deg);
	double unit = pow(0.1, ndec);
	dms[0] = floor(a);
	a = (a - dms[0]) * 60.0;
	dms[1] = floor(a);
	a = (a - dms[1]) * 60.0;
	dms[2] = floor(a / unit + 0.5) * unit;
	if (dms[2] >= 60.0)
	{
		dms[2] = 0.0;
		dms[1] += 1.0;
		if (dms[1] >= 60.0)
		{
			dms[1] = 0.0;
			dms[0] += 1.0;
		}
	}
	dms[0] *= sign;
}

extern int print_nmea_gga(double *ep, double *xyz, int nsat, int type, double dop, double age, char *buff)
{
	double h, pos[3], dms1[3], dms2[3];
	char *p = (char *)buff, *q, sum;
	char buf[20] = {0};

	if ((xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]) < 1.0)
	{
		strcpy(buff,"$GPGGA,,,,,,,,,,,,,,");
		for (q = (char *)buff + 1, sum = 0; *q; q++)
			sum ^= *q;
		strcat(buff, "*");
		memset(buf, 0, 20);
		itoa_user(sum, buf, 16);
		strcat(buff, buf);
		strcat(buff, "\r\n");

	}
	else
	{
		ecef2pos(xyz, pos);
		h = 0.0; 
		deg2dms(fabs(pos[0]) * RAD2DEG, dms1, 7);
		deg2dms(fabs(pos[1]) * RAD2DEG, dms2, 7);

		strcpy(buff,"$GPGGA,");

		RealToArray(ep[3] * 10000 + ep[4] * 100 + ep[5] + 0.001, buf, 6, 2);
		strcat(buff, buf);
		strcat(buff, ",");

		memset(buf, 0, 20);
		RealToArray(dms1[0] * 100 + dms1[1] + dms1[2] / 60.0, buf, 4, 7);
		strcat(buff, buf);
		strcat(buff, pos[0] >= 0 ? ",N," : ",S,");

		memset(buf, 0, 20);
		RealToArray(dms2[0] * 100 + dms2[1] + dms2[2] / 60.0, buf, 5, 7);
		strcat(buff, buf);
		strcat(buff, pos[1] >= 0 ? ",E," : ",W,");

		memset(buf, 0, 20);
		// RealToArray(type,buf,1,0);
		itoa_user(type, buf, 10);
		strcat(buff, buf);
		strcat(buff, ",");

		memset(buf, 0, 20);
		RealToArray(nsat, buf, 2, 0);
		buf[2] = 0;
		strcat(buff, buf);
		strcat(buff, ",");

		memset(buf, 0, 20);
		RealToArray(dop, buf, 0, 1);
		strcat(buff, buf);
		strcat(buff, ",");

		memset(buf, 0, 20);
		RealToArray(pos[2] - h, buf, 0, 3);
		strcat(buff, buf);
		strcat(buff, ",M,");

		memset(buf, 0, 20);
		RealToArray(h, buf, 0, 3);
		strcat(buff, buf);
		strcat(buff, ",M,");

		memset(buf, 0, 20);
		RealToArray(age, buf, 0, 1);
		strcat(buff, buf);
		strcat(buff, ",");

		for (q = (char *)buff + 1, sum = 0; *q; q++)
			sum ^= *q; /* check-sum */

		strcat(buff, "*");
		memset(buf, 0, 20);
		itoa_user(sum, buf, 16);
		strcat(buff, buf);

		strcat(buff, "\r\n");
	}
	return strlen(buff);
}

void RealToArray(double real, char *a, unsigned char id, unsigned char dd) //id:integer digits，dd:decimal digits
{
	uint64_t temp1;
	uint64_t temp2;
	char i;
	if (id == 0)
	{
		if (abs(real) < 10)
		{
			id = 1;
		}
		else if (abs(real) < 100)
		{
			id = 2;
		}
		else if (abs(real) < 1000)
		{
			id = 3;
		}
		else if (abs(real) < 10000)
		{
			id = 4;
		}
		else
		{
			id = 8;
		}
	}

	if (real >= 0)
	{
		temp1 = real * pow(10, dd);
		temp2 = ((int)(real)) * pow(10, dd);
		temp2 = temp1 - temp2; 
		temp1 = (int)real;	 
		for (i = id; i > 0; i--)
		{
			*(a + i - 1) = (int)temp1 % 10 + 0x30;
			temp1 = (int)(temp1 / 10);
		}
		*(a + id) = '.'; 
		for (i = id + dd; i > id; i--)
		{
			*(a + i) = temp2 % 10 + 0x30;
			temp2 = temp2 / 10;
		}
	}
	else
	{
		real = 0 - real;
		temp1 = real * pow(10, dd);
		temp2 = ((int)(real)) * pow(10, dd);
		temp2 = temp1 - temp2; 
		temp1 = (int)real;	 
		*a = 0x2D;			   
		for (i = id; i > 0; i--)
		{
			*(a + i) = (int)temp1 % 10 + 0x30;
			temp1 = (int)(temp1 / 10);
		}
		*(a + id + 1) = '.'; 
		for (i = id + dd + 1; i > id + 1; i--)
		{
			*(a + i) = temp2 % 10 + 0x30;
			temp2 = temp2 / 10;
		}
	}
}

char *itoa_user(int num, char *str, int radix)
{ 
	char index[] = "0123456789ABCDEF";
	unsigned unum; 
	int i = 0, j, k;
	
	if (radix == 10 && num < 0) 
	{
		unum = (unsigned)-num;
		str[i++] = '-';
	}
	else
		unum = (unsigned)num; 
	do
	{
		str[i++] = index[unum % (unsigned)radix];
		unum /= radix;
	} while (unum);
	str[i] = '\0';

	if (str[0] == '-')
		k = 1; 
	else
		k = 0;

	for (j = k; j <= (i - 1) / 2; j++)
	{
		char temp;
		temp = str[j];
		str[j] = str[i - 1 + k - j];
		str[i - 1 + k - j] = temp;
	}
	return str;
}


int print_nmea_gst(double *ep, float *var_llh, char *buff)
{
	char *p = (char *)buff, *q, sum;

	p += sprintf(p, "$GNGST,%02.0f%02.0f%05.2f,%2.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f",
				 ep[3], ep[4], ep[5], 0.0, 0.0, 0.0, var_llh[1], var_llh[0], var_llh[2]);
	for (q = (char *)buff + 1, sum = 0; *q; q++)
		sum ^= *q; /* check-sum */

	p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);

	return p - (char *)buff;
}