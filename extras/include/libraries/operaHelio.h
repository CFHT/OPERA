#ifndef OPERAHELIO_H
#define OPERAHELIO_H

/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMath
 Version: 1.0
 Description: ThisC library implements mathematics routines..
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: eder@cfht.hawaii.edu
 
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

/*! 
 * operaHelio
 * \author John Thorstensen
 * \brief Heliocentric library.
 * \file operaHelio.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
#define  KMS_AUDAY         1731.45683633	/* km per sec in 1 AU/day */
#define  SPEED_OF_LIGHT    299792.458		/* in km per sec ... exact. */
#define  SS_MASS           1.00134198		/* solar system mass in solar units */
#define  J2000             2451545.			/* Julian date at standard epoch */
#define  SEC_IN_DAY        86400.
#define  FLATTEN           0.003352813		/* flattening of earth, 1/298.257 */
#define  EQUAT_RAD         6378137.			/* equatorial radius of earth, meters */
#define  ASTRO_UNIT        1.4959787066e11	/* 1 AU in meters */
#define  RSUN              6.96000e8		/* IAU 1976 recom. solar radius, meters */
#define  RMOON             1.738e6			/* IAU 1976 recom. lunar radius, meters */
#define  KZEN              0.172			/* zenith extinction, mag, for use in lunar sky brightness calculations. */
#define FIRSTJD            2415387.			/* 1901 Jan 1 -- calendrical limit */
#define LASTJD             2488070.			/* 2099 Dec 31 */
#define  EARTH_DIFF        0.05            /* used in numerical
											differentiation to find earth velocity -- this value gives
											about 8 digits of numerical accuracy on the VAX, but is 
											about 3 orders of magnitude larger than the value where roundoff
											errors become apparent. */
#define MAX_ELEMENTS	10

struct time_coord
{
	short sign;  /* carry sign explicitly since -0 not neg. */
	double hh;
	double mm;
	double ss;
};

/* converts a sexigesimal structure into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_to_dec(struct time_coord *sexigesimal);

/* function for converting decimal to sexigesimalylonian hh mm ss.ss */
void dec_to_sexigesimal (double deci, struct time_coord *sexigesimal);

/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_to_dec(const char *sexigesimal);	
/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_hour(const char *sexigesimal);	
/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_minute(const char *sexigesimal);	
/* converts a sexigesimal string into 
 * double-precision floating point ("decimal") number. 
 */
double sexigesimal_str_second(const char *sexigesimal);
	/* 
 computes the geocentric coordinates from the geodetic 
 (standard map-type) longitude, latitude, and height. 
 These are assumed to be in decimal hours, decimal degrees, and
 meters respectively.  Notation generally follows 1992 Astr Almanac, 
 p. K11
 */
void geocentric(double geolong, double geolat, double height, double *x_geo, double *y_geo, double *z_geo);

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

double etcorrection(double jd);

/*
 This routine takes the position
 x,y,z and velocity xdot,ydot,zdot, assumed heliocentric,
 and corrects them to the solar system barycenter taking into
 account the nine major planets.  Routine evolved by inserting
 planetary data (given above) into an earlier, very crude
 barycentric correction. 
 */

void barymetric_correction(double jd, double *x, double *y, double *z, double *xdot, double *ydot, double *zdot);

/* Planetary part, added 1992 August.  The intention of this is
 to compute low-precision planetary positions for general info
 and to inform user if observation might be interfered with by
 a planet -- a rarity, but it happens.  Also designed to make
 handy low-precision planet positions available for casual
 planning purposes.  Do not try to point blindly right at the
 middle of a planetary disk with these routines!  */

/* elements of planetary orbits */
typedef struct element {
	char name[9];
	double incl;
	double Omega;
	double omega;
	double a;
	double daily;
	double ecc;
	double L_0;
	double mass;
} element_t;

void compute_elements(double jd);

/*
 * produces ecliptic x,y,z coordinates for planet number 'p'
 * at date jd.
 */

void planet_xyz(int p, double jd, double *x, double *y, double *z);
void planet_velocity(int p, double jd, double *vx, double *vy, double *vz);
/*
 finds heliocentric correction for given jd, ra, dec, ha, and lat.
 tcor is time correction in seconds, vcor velocity in km/s, to 
 be added to the observed values.  
 Input ra and dec assumed to be at current epoch
 */
void heliocentric_correction(double jd, double ra, double dec, double ha, double lat, double elevsea, double *tcor, double *vcor);

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

void accurate_solar_ephemeris(double jd, double lst, double geolat, double *ra, double *dec, double *dist, double *topora, double *topodec, double *x, double *y, double *z);

/*
 rotates ecliptic rectangular coords x, y, z to
 equatorial (all assumed of date.)
 */

void ecliptic_rotation(double jd, double *x, double *y, double *z);
/*
 assuming x is an angle in degrees, returns 
 modulo 360 degrees.
 */

double circulo(double x);
/* 
 * returns radian angle 0 to 2pi for coords x, y --
 * get that quadrant right !! 
 */

double atan_circ(double x, double y);
	
#ifdef __cplusplus
}
#endif
		
		
#endif
