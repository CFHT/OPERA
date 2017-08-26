/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaHelio
 Version: 1.0
 Description: Heliocentric-related library routines.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <math.h>

#include "globaldefines.h"
#include "operaError.h"

/*!
 * \brief Heliocentric library.
 * \file operaHelio.c/
 * \ingroup libraries
 */

#include "libraries/operaLibCommon.h"
#include "libraries/operaHelio.h"

/* From
 * SKY CALCULATOR PROGRAM 
 * John Thorstensen, Dartmouth College.  
 * Credits:  
 * The julian date and sidereal time routines were 
 * originally coded in PL/I by  Steve Maker of Dartmouth College.  
 * They were based on routines in the old American Ephemeris.
 * Many of the routines were coded from Jean Meeus' "Astronomical
 * Formulae for Calculators", published by Willman-Bell.  This is
 * an extraordinarily helpful little book!
 */

/* 
 * operaHelio
 * \author John Thorstensen
 * \brief Heliocentric library.
 * \ingroup libraries
 */

/* converts a sexigesimal structure into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_to_dec(struct time_coord *sexigesimal)
{
	double x;
	x = sexigesimal->sign * (sexigesimal->hh + sexigesimal->mm / 60. + sexigesimal->ss / 3600.);
	return(x);
}

/* function for converting decimal to sexigesimalylonian hh mm ss.ss */
void dec_to_sexigesimal (double deci, struct time_coord *sexigesimal)
{
	int hr_int, min_int;
	
	if (deci >= 0.) sexigesimal->sign = 1; 
	else {
		sexigesimal->sign = -1;
		deci = -1. * deci;
	}
	hr_int = deci;   /* use conversion conventions to truncate */
	sexigesimal->hh = hr_int;
	min_int = 60. * (deci - sexigesimal->hh);
	sexigesimal->mm = min_int;
	sexigesimal->ss = 3600. * (deci - sexigesimal->hh - sexigesimal->mm / 60.);
}


/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_to_dec(const char *sexigesimal_str)
{
	double x;
	struct time_coord sexigesimal;
	sexigesimal.sign = (sexigesimal_str[0]=='-'?-1:1);
	if (sexigesimal.sign == -1) {
		sexigesimal_str++;
	}
	char *copy = malloc(sizeof(char)*strlen(sexigesimal_str));
	memcpy(copy, sexigesimal_str, strlen(sexigesimal_str));
	char *p = strchr(copy, ':');
	*p = '\0';
	sexigesimal.hh = atof(copy);
	*p = ':';
	p++;
	char *e = strchr(p, ':');
	*e = '\0';
	sexigesimal.mm = atof(p);
	*e = ':';
	p = e++;
	sexigesimal.ss = atof(p);
	x = sexigesimal.sign * (sexigesimal.hh + sexigesimal.mm / 60. + sexigesimal.ss / 3600.);
	free(copy);
	return x;
}

/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_hour(const char *sexigesimal_str)
{
	double x;
	char *copy = malloc(sizeof(char)*strlen(sexigesimal_str));
	memcpy(copy, sexigesimal_str, strlen(sexigesimal_str));
	char *p = strchr(copy, ':');
	*p = '\0';
	x = atof(copy);
	free(copy);
	return x;
}

/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_minute(const char *sexigesimal_str)
{
	double x;
	char *copy = malloc(sizeof(char)*strlen(sexigesimal_str));
	memcpy(copy, sexigesimal_str, strlen(sexigesimal_str));
	char *p = strchr(copy, ':');
	p++;
	char *e = strchr(p, ':');
	*e = '\0';
	x = atof(p);
	free(copy);
	return x;
}

/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_second(const char *sexigesimal_str)
{
	double x;
	char *copy = malloc(sizeof(char)*strlen(sexigesimal_str));
	memcpy(copy, sexigesimal_str, strlen(sexigesimal_str));
	char *p = strchr(copy, ':');
	p++;
	p = strchr(p, ':');
	p++;
	x = atof(p);
	free(copy);
	return x;
}

/* 
 computes the geocentric coordinates from the geodetic 
 (standard map-type) longitude, latitude, and height. 
 These are assumed to be in decimal hours, decimal degrees, and
 meters respectively.  Notation generally follows 1992 Astr Almanac, 
 p. K11
 */
void geocentric(double geolong, double geolat, double height, double *x_geo, double *y_geo, double *z_geo)
{
	
	double denom, C_geo, S_geo;
	
	geolat = geolat / DEG_IN_RADIAN;
	geolong = geolong / HRS_IN_RADIAN;      
	denom = (1. - FLATTEN) * sin(geolat);
	denom = cos(geolat) * cos(geolat) + denom*denom;
	C_geo = 1. / sqrt(denom);
	S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
	C_geo = C_geo + height / EQUAT_RAD;  /* deviation from almanac
										  notation -- include height here. */
	S_geo = S_geo + height / EQUAT_RAD;
	*x_geo = C_geo * cos(geolat) * cos(geolong);
	*y_geo = C_geo * cos(geolat) * sin(geolong);
	*z_geo = S_geo * sin(geolat);
}

/* Given a julian date in 1900-2100, returns the correction
 delta t which is:
 TDT - UT (after 1983 and before 1998)
 ET - UT (before 1983)
 an extrapolated guess  (after 2001). 
 
 For dates in the past (<= 2001 and after 1900) the value is linearly
 interpolated on 5-year intervals; for dates after the present,
 an extrapolation is used, because the true value of delta t
 cannot be predicted precisely.  Note that TDT is essentially the
 modern version of ephemeris time with a slightly cleaner 
 definition.  
 
 Where the algorithm shifts there will be a small (< 0.1 sec)
 discontinuity.  Also, the 5-year linear interpolation scheme can 
 lead to errors as large as 0.5 seconds in some cases, though
 usually rather smaller.   One seldom has actual UT to work with anyway,
 since the commonly-used UTC is tied to TAI within an integer number
 of seconds.  */

double etcorrection(double jd)
{
	
	double jd1900 = 2415019.5;
	double dates[22];
	double delts[22];  /* can't initialize this look-up table
						with stupid old sun compiler .... */
	double year, delt = 0.0;
	int i;
	
	/* this stupid patch for primitive sun C compilers .... 
	 do not allow automatic initialization of arrays! */
	
	for(i = 0; i <= 20; i++) dates[i] = 1900 + (double) i * 5.;
	dates[21] = 2001.;  /* the last accurately tabulated one in the
						 2003 Almanac ... */
	
	delts[0] = -2.72;  delts[1] = 3.86; delts[2] = 10.46;
	delts[3] = 17.20;  delts[4] = 21.16; delts[5] = 23.62;
	delts[6] = 24.02;  delts[7] = 23.93; delts[8] = 24.33;
	delts[9] = 26.77;  delts[10] = 29.15; delts[11] = 31.07;
	delts[12] = 33.15;  delts[13] = 35.73; delts[14] = 40.18;
	delts[15] = 45.48;  delts[16] = 50.54; delts[17] = 54.34;
	delts[18] = 56.86;  delts[19] = 60.78; delts[20] = 63.83;
	delts[21] = 64.09;
	
	year = 1900. + (jd - jd1900) / 365.25;
	
	if(year < 2001. && year >= 1900.) {
		i = (year - 1900) / 5;
		delt = delts[i] + 
		((delts[i+1] - delts[i])/(dates[i+1] - dates[i])) * (year - dates[i]);
	}
	
	else if (year >= 2001. && year < 2100.)
		delt = 31.69 + (2.164e-3) * (jd - 2436935.4);  /* rough extrapolation */
	/* the 31.69 is adjusted to give 64.09 sec at the start of 2001. */
	else if (year < 1900) {
		//printf("etcorrection ... no ephemeris time data for < 1900.\n");
		delt = 0.;
	}
	
	else if (year >= 2100.) {
		//printf("etcorrection .. very long extrapolation in delta T - inaccurate.\n");
		delt = 180.; /* who knows? */
	} 
	
	return(delt);
}

double jd_el;     /* ************** */
element_t elements[MAX_ELEMENTS];

/*
 This routine takes the position
 x,y,z and velocity xdot,ydot,zdot, assumed heliocentric,
 and corrects them to the solar system barycenter taking into
 account the nine major planets.  Routine evolved by inserting
 planetary data (given above) into an earlier, very crude
 barycentric correction. 
 */

void barymetric_correction(double jd, double *x, double *y, double *z, double *xdot, double *ydot, double *zdot)
{
	
	int p;
	double xp, yp, zp, xvp, yvp, zvp;
	
	double xc=0.,yc=0.,zc=0.,xvc=0.,yvc=0.,zvc=0.;
	
	compute_elements(jd);
	
	for(p=1;p<=9;p++) { /* sum contributions of the planets */
		planet_xyz(p,jd,&xp,&yp,&zp);
		xc = xc + elements[p].mass * xp;  /* mass is fraction of solar mass */
		yc = yc + elements[p].mass * yp;
		zc = zc + elements[p].mass * zc;
		planet_velocity(p,jd,&xvp,&yvp,&zvp);
		xvc = xvc + elements[p].mass * xvp;
		yvc = yvc + elements[p].mass * yvp;
		zvc = zvc + elements[p].mass * zvc;
		// diagnostic ... nice place to check planets if needed
#ifdef PRINT_DEBUG
		double xo, yo, zo;  /* for diagn */
		printf("%d :", p);
		xo = xp;
		yo = yp;
		zo = zp;
		ecliptic_rotation(jd,&xo,&yo,&zo);
		xyz2000(jd,xo,yo,zo);
		printf("    ");
		xo = xvp;
		yo = yvp;
		zo = zvp;
		ecliptic_rotation(jd,&xo,&yo,&zo);
		xyz2000(jd,xo,yo,zo);
#endif
	}
	/* normalize properly and rotate corrections to equatorial coords */
	xc = xc / SS_MASS;
	yc = yc / SS_MASS;
	zc = zc / SS_MASS;     /* might as well do it right ... */
	
	ecliptic_rotation(jd, &xc, &yc, &zc);
	
#ifdef PRINT_DEBUG
	printf("posn corrn:");   diagnostic commented out
	xyz2000(jd,xc,yc,zc);
	printf(" vel corrn:");
	xyz2000(jd,(1.0e9 * xvc),(1.0e9 * yvc),(1.0e9 * zvc));
#endif	
	xvc = xvc * KMS_AUDAY / SS_MASS;  
	yvc = yvc * KMS_AUDAY / SS_MASS;
	zvc = zvc * KMS_AUDAY / SS_MASS;
	ecliptic_rotation(jd, &xvc, &yvc, &zvc);
	
	/* add them in */
	*x = *x - xc;  /* these are in AU -- */
	*y = *y - yc;
	*z = *z - zc;
	/*      xyz2000(jd,*x,*y,*z);   */
	*xdot = *xdot - xvc;
	*ydot = *ydot - yvc;
	*zdot = *zdot - zvc;
	
#ifdef PRINT_DEBUG
	// diagnostic -- trash variables  -- scale for direct comparison with almanac
	xp = 1.0e9 * *xdot / KMS_AUDAY;    
	yp = 1.0e9 * *ydot / KMS_AUDAY;
	zp = 1.0e9 * *zdot / KMS_AUDAY;
	xyz2000(jd,xp,yp,zp);
#endif
}

/* Planetary part, added 1992 August.  The intention of this is
 to compute low-precision planetary positions for general info
 and to inform user if observation might be interfered with by
 a planet -- a rarity, but it happens.  Also designed to make
 handy low-precision planet positions available for casual
 planning purposes.  Do not try to point blindly right at the
 middle of a planetary disk with these routines!  */

void compute_elements(double jd) 
{
	
	double T, Tsq, Tcb, d;
	double ups, P, Q, S, V, W, G, H, zeta/*, psi*/; /* Meeus p. 110 ff. */
	double sinQ,sinZeta,cosQ,cosZeta,sinV,cosV,sin2Zeta,cos2Zeta;
	
	jd_el = jd;   /* true, but not necessarily; set explicitly */
	d = jd - 2415020.;   
	T = d / 36525.;
	Tsq = T * T;
	Tcb = Tsq * T;
	
	/* computes and loads mean elements for planets.  */
	
	/* Mercury, Venus, and Mars from Explanatory Suppl., p. 113 */
	
	strcpy(elements[1].name,"Mercury");
	elements[1].incl = 7.002880 + 1.8608e-3 * T - 1.83e-5 * Tsq;
	elements[1].Omega = 47.14594 + 1.185208 * T + 1.74e-4 * Tsq;
	elements[1].omega = 75.899697 + 1.55549 * T + 2.95e-4 * Tsq;
	elements[1].a = .3870986;
	elements[1].daily = 4.0923388;
	elements[1].ecc = 0.20561421 + 0.00002046 * T;
	elements[1].L_0 = 178.179078 + 4.0923770233 * d  +
	0.0000226 * pow((3.6525 * T),2.);
	
	strcpy(elements[2].name,"Venus  ");
	elements[2].incl = 3.39363 + 1.00583e-03 * T - 9.722e-7 * Tsq;
	elements[2].Omega = 75.7796472 + 0.89985 * T + 4.1e-4 * Tsq;
	elements[2].omega = 130.16383 + 1.4080 * T + 9.764e-4 * Tsq;
	elements[2].a = .723325;
	elements[2].daily = 1.60213049;
	elements[2].ecc = 0.00682069 - 0.00004774 * T;
	elements[2].L_0 = 342.767053 + 1.6021687039 * 36525 * T +
	0.000023212 * pow((3.6525 * T),2.);
	
	/* Earth from old Nautical Almanac .... */
	
	strcpy(elements[3].name,"Earth  ");
	elements[3].ecc = 0.01675104 - 0.00004180*T + 0.000000126*Tsq;
	elements[3].incl = 0.0;
	elements[3].Omega = 0.0;
	elements[3].omega = 101.22083 + 0.0000470684*d + 0.000453*Tsq + 0.000003*Tcb;
	elements[3].a = 1.0000007;;
	elements[3].daily = 0.985599;
	elements[3].L_0 = 358.47583 + 0.9856002670*d - 0.000150*Tsq - 0.000003*Tcb +
	elements[3].omega;
	
	strcpy(elements[4].name,"Mars   ");
	elements[4].incl = 1.85033 - 6.75e-04 * T - 1.833e-5 * Tsq;
	elements[4].Omega = 48.786442 + .770992 * T + 1.39e-6 * Tsq;
	elements[4].omega = 334.218203 + 1.840758 * T + 1.299e-4 * Tsq;
	elements[4].a = 1.5236915;
	elements[4].daily = 0.5240329502 + 1.285e-9 * T;
	elements[4].ecc = 0.09331290 - 0.000092064 * T - 0.000000077 * Tsq;
	elements[4].L_0 = 293.747628 + 0.5240711638 * d  +
	0.000023287 * pow((3.6525 * T),2.);
	
	/* Outer planets from Jean Meeus, Astronomical Formulae for
	 Calculators, 3rd edition, Willman-Bell; p. 100. */
	
	strcpy(elements[5].name,"Jupiter");
	elements[5].incl = 1.308736 - 0.0056961 * T + 0.0000039 * Tsq;
	elements[5].Omega = 99.443414 + 1.0105300 * T + 0.0003522 * Tsq 
	- 0.00000851 * Tcb;
	elements[5].omega = 12.720972 + 1.6099617 * T + 1.05627e-3 * Tsq
	- 3.43e-6 * Tcb;
	elements[5].a = 5.202561;
	elements[5].daily = 0.08312941782;
	elements[5].ecc = .04833475  + 1.64180e-4 * T - 4.676e-7*Tsq -
	1.7e-9 * Tcb;
	elements[5].L_0 = 238.049257 + 3036.301986 * T + 0.0003347 * Tsq -
	1.65e-6 * Tcb;
	
	/* The outer planets have such large mutual interactions that
	 even fair accuracy requires lots of perturbations --- here
	 are some of the larger ones, from Meeus' book. */
	
	ups = 0.2*T + 0.1;
	P = (237.47555 + 3034.9061 * T) / DEG_IN_RADIAN;
	Q = (265.91650 + 1222.1139 * T) / DEG_IN_RADIAN;
	S = (243.51721 + 428.4677 * T) / DEG_IN_RADIAN;
	V = 5*Q - 2*P;
	W = 2*P - 6*Q + 3*S;
	zeta = Q - P;
	//psi = S - Q;
	sinQ = sin(Q);
	cosQ = cos(Q);
	sinV = sin(V);
	cosV = cos(V);
	sinZeta = sin(zeta);
	cosZeta = cos(zeta);
	sin2Zeta = sin(2*zeta);
	cos2Zeta = cos(2*zeta);
	
	elements[5].L_0 = elements[5].L_0 
	+ (0.331364 - 0.010281*ups - 0.004692*ups*ups)*sinV
	+ (0.003228 - 0.064436*ups + 0.002075*ups*ups)*cosV
	- (0.003083 + 0.000275*ups - 0.000489*ups*ups)*sin(2*V)
	+ 0.002472 * sin(W) + 0.013619 * sinZeta + 0.018472 * sin2Zeta 
	+ 0.006717 * sin(3*zeta) 
	+ (0.007275  - 0.001253*ups) * sinZeta * sinQ
	+ 0.006417 * sin2Zeta * sinQ  
	- (0.033839 + 0.001253 * ups) * cosZeta * sinQ 
	- (0.035681 + 0.001208 * ups) * sinZeta * sinQ;
	/* only part of the terms, the ones first on the list and
	 selected larger-amplitude terms from farther down. */
	
	elements[5].ecc = elements[5].ecc + 1e-7 * (
												(3606 + 130 * ups - 43 * ups*ups) * sinV 
												+ (1289 - 580 * ups) * cosV - 6764 * sinZeta * sinQ 
												- 1110 * sin2Zeta * sin(Q) 
												+ (1284 + 116 * ups) * cosZeta * sinQ 
												+ (1460 + 130 * ups) * sinZeta * cosQ 
												+ 6074 * cosZeta * cosQ);
	
	elements[5].omega = elements[5].omega 
	+ (0.007192 - 0.003147 * ups) * sinV
	+ ( 0.000197*ups*ups - 0.00675*ups - 0.020428) * cosV
	+ 0.034036 * cosZeta * sinQ + 0.037761 * sinZeta * cosQ;
	
	elements[5].a = elements[5].a + 1.0e-6 * (
											  205 * cosZeta - 263 * cosV + 693 * cos2Zeta + 312 * sin(3*zeta)
											  + 147 * cos(4*zeta) + 299 * sinZeta * sinQ 
											  + 181 * cos2Zeta * sinQ + 181 * cos2Zeta * sinQ
											  + 204 * sin2Zeta * cosQ + 111 * sin(3*zeta) * cosQ 
											  - 337 * cosZeta * cosQ - 111 * cos2Zeta * cosQ
											  );
	
	strcpy(elements[6].name,"Saturn ");
	elements[6].incl = 2.492519 - 0.00034550*T - 7.28e-7*Tsq;
	elements[6].Omega = 112.790414 + 0.8731951*T - 0.00015218*Tsq - 5.31e-6*Tcb ;
	elements[6].omega = 91.098214 + 1.9584158*T + 8.2636e-4*Tsq;
	elements[6].a = 9.554747;
	elements[6].daily = 0.0334978749897;
	elements[6].ecc = 0.05589232 - 3.4550e-4 * T - 7.28e-7*Tsq;
	elements[6].L_0 = 266.564377 + 1223.509884*T + 0.0003245*Tsq - 5.8e-6*Tcb 
	+ (0.018150*ups - 0.814181 + 0.016714 * ups*ups) * sinV 
	+ (0.160906*ups - 0.010497 - 0.004100 * ups*ups) * cosV
	+ 0.007581 * sin(2*V) - 0.007986 * sin(W) 
	- 0.148811 * sinZeta - 0.040786*sin2Zeta 
	- 0.015208 * sin(3*zeta) - 0.006339 * sin(4*zeta) 
	- 0.006244 * sinQ
	+ (0.008931 + 0.002728 * ups) * sinZeta * sinQ 
	- 0.016500 * sin2Zeta * sinQ
	- 0.005775 * sin(3*zeta) * sinQ 
	+ (0.081344 + 0.003206 * ups) * cosZeta * sinQ 
	+ 0.015019 * cos2Zeta * sinQ 
	+ (0.085581 + 0.002494 * ups) * sinZeta * cosQ
	+ (0.025328 - 0.003117 * ups) * cosZeta * cosQ
	+ 0.014394 * cos2Zeta * cosQ;   /* truncated here -- no
									 terms larger than 0.01 degrees, but errors may
									 accumulate beyond this.... */
	elements[6].ecc = elements[6].ecc + 1.0e-7 * (
												  (2458 * ups - 7927.) * sinV + (13381. + 1226. * ups) * cosV
												  + 12415. * sinQ + 26599. * cosZeta * sinQ 
												  - 4687. * cos2Zeta * sinQ - 12696. * sinZeta * cosQ 
												  - 4200. * sin2Zeta * cosQ +(2211. - 286*ups) * sinZeta*sin(2*Q)
												  - 2208. * sin2Zeta * sin(2*Q) 
												  - 2780. * cosZeta * sin(2*Q) + 2022. * cos2Zeta*sin(2*Q) 
												  - 2842. * sinZeta * cos(2*Q) - 1594. * cosZeta * cos(2*Q)
												  + 2162. * cos2Zeta*cos(2*Q) );  /* terms with amplitudes
																				   > 2000e-7;  some secular variation ignored. */
	elements[6].omega = elements[6].omega 
	+ (0.077108 + 0.007186 * ups - 0.001533 * ups*ups) * sinV
	+ (0.045803 - 0.014766 * ups - 0.000536 * ups*ups) * cosV
	- 0.075825 * sinZeta * sinQ - 0.024839 * sin2Zeta*sinQ
	- 0.072582 * cosQ - 0.150383 * cosZeta * cosQ +
	0.026897 * cos2Zeta * cosQ;  /* all terms with amplitudes 
								  greater than 0.02 degrees -- lots of others! */
	elements[6].a = elements[6].a + 1.0e-6 * (
											  2933. * cosV + 33629. * cosZeta - 3081. * cos2Zeta 
											  - 1423. * cos(3*zeta) + 1098. * sinQ - 2812. * sinZeta * sinQ 
											  + 2138. * cosZeta * sinQ  + 2206. * sinZeta * cosQ 
											  - 1590. * sin2Zeta*cosQ + 2885. * cosZeta * cosQ 
											  + 2172. * cos2Zeta * cosQ);  /* terms with amplitudes greater
																			than 1000 x 1e-6 */
	
	strcpy(elements[7].name,"Uranus ");
	elements[7].incl = 0.772464 + 0.0006253*T + 0.0000395*Tsq;
	elements[7].Omega = 73.477111 + 0.4986678*T + 0.0013117*Tsq;
	elements[7].omega = 171.548692 + 1.4844328*T + 2.37e-4*Tsq - 6.1e-7*Tcb;
	elements[7].a = 19.21814;
	elements[7].daily = 1.1769022484e-2;
	elements[7].ecc = 0.0463444 - 2.658e-5 * T;
	elements[7].L_0 = 244.197470 + 429.863546*T + 0.000316*Tsq - 6e-7*Tcb;
	/* stick in a little bit of perturbation -- this one really gets
	 yanked around.... after Meeus p. 116*/
	G = (83.76922 + 218.4901 * T)/DEG_IN_RADIAN;
	H = 2*G - S;
	elements[7].L_0 = elements[7].L_0 + (0.864319 - 0.001583 * ups) * sin(H)
	+ (0.082222 - 0.006833 * ups) * cos(H)
	+ 0.036017 * sin(2*H);
	elements[7].omega = elements[7].omega + 0.120303 * sin(H) 
	+ (0.019472 - 0.000947 * ups) * cos(H)
	+ 0.006197 * sin(2*H);
	elements[7].ecc = elements[7].ecc + 1.0e-7 * (
												  20981. * cos(H) - 3349. * sin(H) + 1311. * cos(2*H));
	elements[7].a = elements[7].a - 0.003825 * cos(H);
	
	/* other corrections to "true longitude" are ignored. */
	
	strcpy(elements[8].name,"Neptune");
	elements[8].incl = 1.779242 - 9.5436e-3 * T - 9.1e-6*Tsq;
	elements[8].Omega = 130.681389 + 1.0989350 * T + 2.4987e-4*Tsq - 4.718e-6*Tcb;
	elements[8].omega = 46.727364 + 1.4245744*T + 3.9082e-3*Tsq - 6.05e-7*Tcb;
	elements[8].a = 30.10957;
	elements[8].daily = 6.020148227e-3;
	elements[8].ecc = 0.00899704 + 6.33e-6 * T;
	elements[8].L_0 = 84.457994 + 219.885914*T + 0.0003205*Tsq - 6e-7*Tcb;
	elements[8].L_0 = elements[8].L_0 
	- (0.589833 - 0.001089 * ups) * sin(H) 
	- (0.056094 - 0.004658 * ups) * cos(H)
	- 0.024286 * sin(2*H);
	elements[8].omega = elements[8].omega + 0.024039 * sin(H) 
	- 0.025303 * cos(H);
	elements[8].ecc = elements[8].ecc + 1.0e-7 * (
												  4389. * sin(H) + 1129. * sin(2.*H)
												  + 4262. * cos(H) + 1089. * cos(2.*H));
	elements[8].a = elements[8].a + 8.189e-3 * cos(H); 
	
	/* crummy -- osculating elements a la Sept 15 1992 */
	
	d = jd - 2448880.5;  /* 1992 Sep 15 */       
	T = d / 36525.;
	strcpy(elements[9].name,"Pluto  ");
	elements[9].incl = 17.1426;
	elements[9].Omega = 110.180;
	elements[9].omega = 223.782;
	elements[9].a = 39.7465;
	elements[9].daily = 0.00393329;
	elements[9].ecc = 0.253834;
	elements[9].L_0 = 228.1027 + 0.00393329 * d;
#ifdef PRINT_DEBUG
	/printf("inc Om om : %f %f %f\n",elements[9].incl,elements[9].Omega,elements[9].omega);
	printf("a  dail ecc: %f %f %f\n",elements[9].a,elements[9].daily,elements[9].ecc);
	printf("L_0 %f\n",elements[9].L_0);
#endif
	elements[1].mass = 1.660137e-7;  /* in units of sun's mass, IAU 1976 */
	elements[2].mass = 2.447840e-6;  /* from 1992 *Astron Almanac*, p. K7 */
	elements[3].mass = 3.040433e-6;  /* earth + moon */
	elements[4].mass = 3.227149e-7;
	elements[5].mass = 9.547907e-4;
	elements[6].mass = 2.858776e-4; 
	elements[7].mass = 4.355401e-5;  
	elements[8].mass = 5.177591e-5;
	elements[9].mass = 7.69e-9;  /* Pluto+Charon -- ? */
	
}

/*
 * produces ecliptic x,y,z coordinates for planet number 'p'
 * at date jd.
 */

void planet_xyz(int p, double jd, double *x, double *y, double *z)
{
	double M, omnotil, nu, r;
	double e, LL, Om, om, ii;
	
	/* see 1992 Astronomical Almanac, p. E 4 for these formulae. */
	
	ii = elements[p].incl/DEG_IN_RADIAN;
	e = elements[p].ecc;
	
	LL = (elements[p].daily * (jd - jd_el) + elements[p].L_0) / DEG_IN_RADIAN;
	Om = elements[p].Omega / DEG_IN_RADIAN;
	om = elements[p].omega / DEG_IN_RADIAN;
	
	M = LL - om;
	omnotil = om - Om;
	nu = M + (2.*e - 0.25 * pow(e,3.)) * sin(M) +
	1.25 * e * e * sin(2 * M) +
	1.08333333 * pow(e,3.) * sin(3 * M);
	r = elements[p].a * (1. - e*e) / (1 + e * cos(nu));
	
	*x = r * 
	(cos(nu + omnotil) * cos(Om) - sin(nu +  omnotil) * 
	 cos(ii) * sin(Om));
	*y = r * 
	(cos(nu +  omnotil) * sin(Om) + sin(nu +  omnotil) * 
	 cos(ii) * cos(Om));
	*z = r * sin(nu +  omnotil) * sin(ii);
}


void planet_velocity(int p, double jd, double *vx, double *vy, double *vz)
{
	/* numerically evaluates planet velocity by brute-force
	 numerical differentiation. Very unsophisticated algorithm. */
	
	double dt; /* timestep */
	double x1,y1,z1,x2,y2,z2;
	
	dt = 0.1 / elements[p].daily; /* time for mean motion of 0.1 degree */
	planet_xyz(p, (jd - dt), &x1, &y1, &z1);
	planet_xyz(p, (jd + dt), &x2, &y2, &z2);
	*vx = 0.5 * (x2 - x1) / dt;    
	*vy = 0.5 * (y2 - y1) / dt;
	*vz = 0.5 * (z2 - z1) / dt;   
	/* answer should be in ecliptic coordinates, in AU per day.*/
}


/*
 finds heliocentric correction for given jd, ra, dec, ha, and lat.
 tcor is time correction in seconds, vcor velocity in km/s, to 
 be added to the observed values.  
 Input ra and dec assumed to be at current epoch
 */
void heliocentric_correction(double jd, double ra, double dec, double ha, double lat, double elevsea, double *tcor, double *vcor)
{
	double x, y, z, xdot, ydot, zdot;
	double xobj,yobj,zobj;
	double ras, decs, dists, jd1, jd2, x1, x2, y1, y2, z1, z2;
	double topora, topodec;
	double x_geo, y_geo, z_geo;		/* geocentric coords of observatory */
	double a=499.0047837;			/* light travel time for 1 AU, sec  */
	
	dec=dec/DEG_IN_RADIAN; /* pass by value! */
	ra = ra/HRS_IN_RADIAN;
	ha = ha/HRS_IN_RADIAN;
	
	xobj = cos(ra) * cos(dec);
	yobj = sin(ra) * cos(dec);
	zobj = sin(dec);
	
#ifdef PRINT_DEBUG
	//diagnostic -- temporarily trash jd1 
	 jd1 = jd1 + etcorrection(jd1) / SEC_IN_DAY;
	 printf("TDT: %20f\n",jd1);
#endif	
	jd1 = jd - EARTH_DIFF;
	jd2 = jd + EARTH_DIFF;
	
	accurate_solar_ephemeris(jd1,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x1,&y1,&z1);
	accurate_solar_ephemeris(jd2,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x2,&y2,&z2);
	accurate_solar_ephemeris(jd,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x,&y,&z);
	
#ifdef PRINT_DEBUG
	printf("ra dec distance:");  diagnostic -- commented out
	put_coords(ras,3,0);
	printf(" ");
	put_coords(decs,2,1);
	printf(" %f\n",dists);
#endif	
	xdot = KMS_AUDAY*(x2 - x1)/(2.*EARTH_DIFF);  /* numerical differentiation */
	ydot = KMS_AUDAY*(y2 - y1)/(2.*EARTH_DIFF);  /* crude but accurate */
	zdot = KMS_AUDAY*(z2 - z1)/(2.*EARTH_DIFF);
	
#ifdef PRINT_DEBUG
	printf("Helio earth:");  diagnostic -- commmented out
	xyz2000(jd,x,y,z);
	xyz2000(jd,xdot,ydot,zdot);
#endif	
	barymetric_correction(jd,&x,&y,&z,&xdot,&ydot,&zdot);
	*tcor = a * (x*xobj + y*yobj + z*zobj);
	*vcor = xdot * xobj + ydot * yobj + zdot * zobj;
	/* correct diurnal rotation for elliptical earth including obs. elevation */
	geocentric(0., lat, elevsea, &x_geo, &y_geo, &z_geo);
	/* longitude set to zero arbitrarily so that x_geo = perp. distance to axis */
	*vcor = *vcor - 0.4651011 * x_geo * sin(ha) * cos(dec);
	/* 0.4651011 = 6378.137 km radius * 2 pi / (86164.1 sec per sidereal day) */
	/* could add time-of-flight across earth's radius here -- but rest of
	 theory is not good to 0.02 seconds anyway. */
}

/*
 implemenataion of Jean Meeus' more accurate solar
 ephemeris.  For ultimate use in helio correction! From
 Astronomical Formulae for Calculators, pp. 79 ff.  This
 gives sun's position wrt *mean* equinox of date, not
 *apparent*.  Accuracy is << 1 arcmin.  Positions given are
 geocentric ... parallax due to observer's position on earth is 
 ignored. This is up to 8 arcsec; routine is usually a little 
 better than that. 
 // -- topocentric correction *is* included now. -- //
 Light travel time is apparently taken into
 account for the ra and dec, but I don't know if aberration is
 and I don't know if distance is simlarly antedated. 
 
 x, y, and z are heliocentric equatorial coordinates of the
 EARTH, referred to mean equator and equinox of date.
 */

void accurate_solar_ephemeris(double jd, double lst, double geolat, double *ra, double *dec, double *dist, double *topora, double *topodec, double *x, double *y, double *z)  
{
	double L, T, Tsq, Tcb;
	double M, e, Cent, nu, sunlong;
	double /*Lrad,*/ Mrad, nurad, R;
	double A, B, C, D, E, H;
	double xtop, ytop, ztop, topodist, l, m, n, xgeo, ygeo, zgeo;
	
	jd = jd + etcorrection(jd)/SEC_IN_DAY;  /* might as well do it right .... */
	T = (jd - 2415020.) / 36525.;  /* 1900 --- this is an oldish theory*/
	Tsq = T*T;
	Tcb = T*Tsq;
	L = 279.69668 + 36000.76892*T + 0.0003025*Tsq;
	M = 358.47583 + 35999.04975*T - 0.000150*Tsq - 0.0000033*Tcb;
	e = 0.01675104 - 0.0000418*T - 0.000000126*Tsq;
	
	L = circulo(L);
	M = circulo(M);
	/*      printf("raw L, M: %15.8f, %15.8f\n",L,M); */
	
	A = 153.23 + 22518.7541 * T;  /* A, B due to Venus */
	B = 216.57 + 45037.5082 * T;
	C = 312.69 + 32964.3577 * T;  /* C due to Jupiter */
	/* D -- rough correction from earth-moon 
	 barycenter to center of earth. */
	D = 350.74 + 445267.1142*T - 0.00144*Tsq;  
	E = 231.19 + 20.20*T;    /* "inequality of long period .. */
	H = 353.40 + 65928.7155*T;  /* Jupiter. */
	
	A = circulo(A) / DEG_IN_RADIAN;
	B = circulo(B) / DEG_IN_RADIAN;
	C = circulo(C) / DEG_IN_RADIAN;
	D = circulo(D) / DEG_IN_RADIAN;
	E = circulo(E) / DEG_IN_RADIAN;
	H = circulo(H) / DEG_IN_RADIAN;
	
	L = L + 0.00134 * cos(A) 
	+ 0.00154 * cos(B)
	+ 0.00200 * cos(C)
	+ 0.00179 * sin(D)
	+ 0.00178 * sin(E);
	
	//Lrad = L/DEG_IN_RADIAN;
	Mrad = M/DEG_IN_RADIAN;
	
	Cent = (1.919460 - 0.004789*T -0.000014*Tsq)*sin(Mrad)
		+ (0.020094 - 0.000100*T) * sin(2.0*Mrad)
		+ 0.000293 * sin(3.0*Mrad);
		sunlong = L + Cent;
	
	
	nu = M + Cent;
	nurad = nu / DEG_IN_RADIAN;
	
	R = (1.0000002 * (1 - e*e)) / (1. + e * cos(nurad));
	R = R + 0.00000543 * sin(A)
	+ 0.00001575 * sin(B)
	+ 0.00001627 * sin(C)
	+ 0.00003076 * cos(D)
	+ 0.00000927 * sin(H);
#ifdef PRINT_DEBUG
	printf("solar longitude: %10.5f  Radius vector %10.7f\n",sunlong,R);
	printf("eccentricity %10.7f  eqn of center %10.5f\n",e,Cent);
#endif	
	sunlong = sunlong/DEG_IN_RADIAN;
	
	*dist = R;
	*x = cos(sunlong);  /* geocentric */
	*y = sin(sunlong);
	*z = 0.;
	ecliptic_rotation(jd, x, y, z);
	
	/*      --- code to include topocentric correction for sun .... */
	
	geocentric(lst,geolat,0.,&xgeo,&ygeo,&zgeo);
	
	xtop = *x - xgeo*EQUAT_RAD/ASTRO_UNIT;
	ytop = *y - ygeo*EQUAT_RAD/ASTRO_UNIT;
	ztop = *z - zgeo*EQUAT_RAD/ASTRO_UNIT;
	
	topodist = sqrt(xtop*xtop + ytop*ytop + ztop*ztop);
	
	l = xtop / (topodist);
	m = ytop / (topodist);
	n = ztop / (topodist);
	
	*topora = atan_circ(l,m) * HRS_IN_RADIAN;
	*topodec = asin(n) * DEG_IN_RADIAN; 
	
	*ra = atan_circ(*x,*y) * HRS_IN_RADIAN;
	*dec = asin(*z) * DEG_IN_RADIAN; 
	
	*x = *x * R * -1;  /* heliocentric */
	*y = *y * R * -1;
	*z = *z * R * -1;
	
}

/*
 rotates ecliptic rectangular coords x, y, z to
 equatorial (all assumed of date.)
 */

void ecliptic_rotation(double jd, double *x, double *y, double *z)
{
	double incl;
	double ypr,zpr;
	double T;
	
	T = (jd - J2000) / 36525;  /* centuries since J2000 */
	
	incl = (23.439291 + T * (-0.0130042 - 0.00000016 * T))/DEG_IN_RADIAN; 
	/* 1992 Astron Almanac, p. B18, dropping the 
	 cubic term, which is 2 milli-arcsec! */
	ypr = cos(incl) * *y - sin(incl) * *z;
	zpr = sin(incl) * *y + cos(incl) * *z;
	*y = ypr;
	*z = zpr;
	/* x remains the same. */       
}

/*
 assuming x is an angle in degrees, returns 
 modulo 360 degrees.
 */

double circulo(double x)
{
	int n;
	
	n = (int)(x / 360.);
	return(x - 360. * n);
}       

/* 
 * returns radian angle 0 to 2pi for coords x, y --
 * get that quadrant right !! 
 */

double atan_circ(double x, double y)
{
	
	double theta;
	
	if((x == 0.) && (y == 0.)) return(0.);  /* guard ... */
	
	theta = atan2(y,x);  /* turns out there is such a thing in math.h */
	while(theta < 0.) 
		theta += TWOPI;
	return(theta);
}


