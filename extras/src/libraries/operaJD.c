/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaJD
 Version: 1.0
 Description: Julian Date-related library routines.
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
 * \brief Julian Date library.
 * \file operaJD.c
 * \ingroup libraries
 */

#include "libraries/operaLibCommon.h"
#include "libraries/operaJD.h"

/*  \fn long getJulianDay(int y, int m, int d)
 *  \brief Calculate Julian day
 *  \param y: year
 *  \param m: month
 *  \param d: day
 *  \return the Julian day
 *
 */
long getJulianDay(int y, int m, int d) {
	y += 8000;
	if (m < 3) {
		y--;
		m += 12;
	}
	return (y*365)+(y/4)-(y/100)+(y/400)-1200820+(m*153+3)/5-92+d-1;
}

/*
 * Julian Date library
 *	- calculate the JD, MJD, HJD dates
 *
 * (c) 2002 Richard Ogley.
 *
 * Version 0.1   06/02/2002  RNO.
 *
 * Does not rely on the Starlink heliocentric and barycentric position
 * and velocity routines.  This is system independent.
 */

/* 
 * operaJD
 * \author Richard Ogley
 * \brief Julian Date library.
 * \ingroup libraries
 */

/*
 * Julian_date
 *           - procedure to calculate the Julian date
 *
 * Inputs
 *    greg_time - a time structure of the form tm in <time.h>
 *
 * Returns
 *   Julian date
 *
 */

double Julian_date(struct tm *greg_time)
{
	
	double date;  /* A running variable calculating the Julian date */
	int month;    /* The month in the year */
	int day;      /* The day in the month */
	int year;     /* The years from 1900 */
	
	/* Get the day, month and year from the greg_time structure and form
     them into 1/1/2000 rather than 1/0/100 from greg_time. */
	
	day = greg_time->tm_mday;
	month = greg_time->tm_mon + 1;
	year = greg_time->tm_year + 1900;
	
	/* Calculate the Julian days */
	
	date = day - 32076 + 
    1461*(year + 4800 + (month - 14)/12)/4 +
    367*(month - 2 - (month - 14)/12*12)/12 - 
    3*((year + 4900 + (month - 14)/12)/100)/4;
	
	/* Add the fractional hours, mins and seconds */
	
	date += (greg_time->tm_hour + 12.0)/24.0;
	date += (greg_time->tm_min)/1440.0;
	date += (greg_time->tm_sec)/86400.0;
	
	return date;
}

/*
 * modified_Julian_date
 *           - procedure to calculate the modified Julian date
 *
 * Inputs
 *   jd  - the Julian date
 *
 * Returns
 *   the modified Julian date
 *
 */

double modified_Julian_date(double jd)
{
	
	return jd - 2400000.5;
}

/*
 * heliocentric_modified_Julian_date
 *           - procedure to calculate the heliocentric modified Julian
 *             date from the MJD value.  Requests the co-ordinates of
 *             the source in the procedure.
 *
 * Inputs
 *   mjd  - the modified Julian date
 *
 * Returns
 *   the modified heliocentric Julian date
 *
 */

double heliocentric_modified_Julian_date(double mjd, double ra, double dec)
{
	
	const double ausec = 499.01265;   /* Time it takes light to travel 1 AU */
	const double degtorad = 0.01745329251; /* Conversion factor from degrees to radians */
	
	struct coordinates {
		double ra;   /* The right ascension in degrees */
		double dec;  /* The declination in degrees */
		double x;    /* The X position in AU */
		double y;    /* The Y position in AU */
		double z;    /* The Z position in AU */
	};
	
	struct coordinates earth;  /* Earth co-ordinates */
	struct coordinates source; /* Source co-ordinates */
	
	double correction_secs; /* Correction factor in seconds from MJD to HMJD */
	double cel;   /* intermediate calculation from spherical to cartesian co-ordinates */
	double hmjd;  /* The modified heliocentric Julian date */
	
	/* Defaults */
	
	earth.ra = earth.dec = 0.0;
	
	/* Mauna Kea */
	
	source.ra = ra;		// 10.36478;
	source.dec = dec;	// 19.8267;
	
	/* Attempt to find the RA and Dec of the earth using the astronomical
     calculator book. */
	
	heliocentric_ra_dec(mjd, &earth.ra, &earth.dec);
	
	/* Calculate the heliocentric co-ordinates as X, Y and Z terms */
	cel = cos(earth.dec * degtorad);
	earth.x = cos(earth.ra * degtorad) * cel;
	earth.y = sin(earth.ra * degtorad) * cel;
	earth.z = sin(earth.dec * degtorad);
	
	/* Calculate the X,Y,Z co-ordinates of your source */
	cel = cos(source.dec * degtorad);
	source.x = cos(source.ra * degtorad) * cel;
	source.y = sin(source.ra * degtorad) * cel;
	source.z = sin(source.dec * degtorad);
	
	/* Calculate the correction in seconds for the light travel time
     between the source and the earth vectors in a heliocentric
     reference frame. */
	
	correction_secs = ausec * 
    (earth.x * source.x + earth.y * source.y + earth.z * source.z);
	
	/* Modify the MJD in a heliocentric reference frame */
	
	hmjd = mjd + (correction_secs / (24.0 * 3600));
	
	return hmjd;
}

/*
 * heliocentric_ra_dec
 *           - procedure to calculate heliocentric right ascension and
 *             declination of the earth at a given date.
 *
 * Inputs
 *   mjd  - the modified Julian date
 *
 * Returns
 *   *ra  - the right ascension of the earth
 *   *dec - the declination of the earth
 *
 */

void heliocentric_ra_dec(double mjd, double *ra, double *dec)
{
	
	const double eccentricity = 0.016718;		/* Eccentricity of the Earth's orbit*/
	const double ecliptic_long = 278.833540;	/* The longitude of the ecliptic at 1 Jan 1980 0:00 UT */
	const double perigee_long = 282.596403;		/* The longitude of perigee at 1 Jan 1980 00:00 UT */
	const double deg_to_rad = PI / 180.0;		/* A degrees to radians conversion */
	const double tropical_year = 365.24219572;	/* The length of the tropical year in days */
	const double obliquity = 23.441884;			/* The obliquity of the orbit */
	const double mjd_1980 = 44238.0;			/* The MJD on 1 Jan 1980 00:00 UT */
	
	double mean_anomoly;						/* The mean anomoly of the sun */
	double days_from_1980;						/* The number of days from 1 Jan 1980 */
	double solar_longitude;						/* The longitude of the sun */
	double number_of_deg;						/* The number of degrees in longitude the sun has travelled */
	double equation_of_centres;					/* The value for the equation of centres */
	double x;									/* An X position */
	double y;									/* A  Y position */
	double beta;								/* The ecliptic longitude */ 
	double number_of_rotations;					/* An integer number of solar orbits */
	
	/* Calculate the number of days from 1 Jan 1980 */
	
	days_from_1980 = (mjd - mjd_1980);
	
	/* Calculate the number of degrees around in the orbit travelled in this time */
	
	number_of_deg = (360.0 / tropical_year) * days_from_1980;
	
	/* Adjust so the number of degrees is between 0 and 360 */
	
	if ( (number_of_deg < 0.0) || (number_of_deg > 360.0)) {
		number_of_rotations = number_of_deg / 360.0;
		number_of_rotations = floor(number_of_rotations);
		number_of_deg -= number_of_rotations * 360.0;
	}
	
	/* Calculate the mean anomoly */
	
	mean_anomoly = number_of_deg - perigee_long + ecliptic_long;
	
	/* Since the orbit is elliptical and not circular, calculate the equation of centres */
	
	equation_of_centres = (360.0 / PI) * eccentricity *
    sin(mean_anomoly * deg_to_rad);
	
	/* Calculate the solar longitude */
	
	solar_longitude = number_of_deg + equation_of_centres + ecliptic_long;
	if (solar_longitude > 360.0)
		solar_longitude -= 360.0;
	
	/* The ecliptic latitude is zero for the Sun. */
	
	beta = 0.0;
	
	/* Calculate the RA and Dec of the sun */
	
	*dec = asin( (sin(beta * deg_to_rad) * cos(obliquity * deg_to_rad)) +
				(cos(beta * deg_to_rad) * sin(obliquity * deg_to_rad) *
				 sin(solar_longitude * deg_to_rad)) );
	
	*dec /= deg_to_rad;
	x = cos(solar_longitude * deg_to_rad);
	y = (sin(solar_longitude * deg_to_rad) * 
		 cos(obliquity * deg_to_rad)) -
		(tan(beta * deg_to_rad) * sin(obliquity * deg_to_rad));
	*ra = atan(y/x);
	*ra /= deg_to_rad;
	
	if (*ra < 0.0)
		*ra += 360.0;
	
	*ra /= 15.0;
	
	/* Convert from geocentric to heliocentric co-ordinates for the Earth*/
	
	*dec *= -1.0;
	*ra  -= 12.0;
	if (*ra < 0.0) {
		*ra += 24.0;
	}
	
}


