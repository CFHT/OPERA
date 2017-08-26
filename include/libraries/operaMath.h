#ifndef LIBOPERAMATH_H
#define LIBOPERAMATH_H
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

/*! 
 * \brief math library.
 * \file operaMatch.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
	// for Linux...
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define TENTOMINUSSIX 1.0e-6
#define JMAX 20

double DiffPolynomialFunction(double x, const double *p, int n_par);
double LengthofPolynomial(double a, double b, const double *p, int n_par);
double PolyLengthFunc(double x, const double *p, int n_par);
double PolyLengthTrapezoid(const double *p, int n_par, double a, double b, int n, double initials);
double PlanckLaw(double Temperature, double Wavelength);
	
#ifdef __cplusplus
}
#endif
		
		
#endif
