#ifndef OPERAJD_H
#define OPERAJD_H

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
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

#ifdef __cplusplus
extern "C" {
#endif
	
#include <time.h>
#include <math.h>

struct time_coord
{
	short sign;  /* carry sign explicitly since -0 not neg. */
	double hh;
	double mm;
	double ss;
};

/*! \brief Julian Date library. */
/*! \file operaJD.c */
/*! \package operaJD */

#include "libraries/operaLibCommon.h"
#include "libraries/operaJD.h"

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

/*! 
 * operaJD
 * \author Richard Ogley
 * \brief Julian Date library.
 * \file operaJD.h
 * \ingroup libraries
 */

/*! \fn long getJulianDay(int y, int m, int d)
 *  \brief Calculate Julian day
 *  \param y: year
 *  \param m: month
 *  \param d: day
 *  \return the Julian day
 *
 */
long getJulianDay(int y, int m, int d);
	
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

double Julian_date(struct tm *greg_time);
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

double modified_Julian_date(double jd);

/*
 * heliocentric_modified_Julian_date
 *           - procedure to calculate the heliocentric modified Julian
 *             date from the MJD value.  Requests the co-ordinates of
 *             the source in the procedure.
 *
 * Inputs
 *   mjd  - the modified Julian date
 *   ra = right ascension in decimal format
 *   dec = declination in decimal format
 *
 * Returns
 *   the modified heliocentric Julian date
 *
 */

double heliocentric_modified_Julian_date(double mjd, double ra, double dec);
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

void heliocentric_ra_dec(double mjd, double *ra, double *dec);

#ifdef __cplusplus
}
#endif
		
		
#endif
