/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMath
 Version: 1.0
 Description: This C library implements mathematics routines..
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

#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaMath.h"
#include "libraries/operaFit.h"

/*!
 * operaMath
 * \author Eder Martioli
 * \brief operaMath library.
 * \details {This library contains the basic routines for math tools.}
 * \file operaMath.c
 * \ingroup libraries
 */


/* 
 * double DiffPolynomialFunction(double x, const double *p, int n_par)
 * \brief This function returns the value of the first derivative of a polynomial function.
 * \param x is a double input value for which the given polynomial is evaluated  
 * \param p is a const double pointer that contains the coefficients of the given polynomial
 * \param n_par is an int that defines the order of the polynomial
 * \return double value for the first derivative of the polynomial 
 */
double DiffPolynomialFunction(double x, const double *p, int n_par)
{
	double fpoly = p[1];
	
	for(unsigned i=2;i<n_par;i++) {
		fpoly += i*p[i]*pow(x,(double)i-1); 
	}
	return fpoly;
}
/* 
 * double LengthofPolynomial(double x0, double xf, const double *p, int n_par)
 * \brief This function returns the value of the length of a polynomial function in the interval x0:xf.
 * \param x0 is a double input initial value for which the polynomial is evaluated
 * \param xf is a double input final value for which the polynomial is evaluated
 * \param p is a const double pointer that contains the coefficients of the given polynomial
 * \param n_par is an int that defines the order of the polynomial
 * \return double value for the length of the polynomial
 */
double LengthofPolynomial(double a, double b, const double *p, int n_par)
{
	int j;
	double s,st,ost=0.0,os=0.0, initials = 0.0;
	
	st = PolyLengthTrapezoid(p,n_par,a,b,1,initials);	// set to 1 to initialize
	for (j=1; j<=JMAX; j++) {
		st = PolyLengthTrapezoid(p,n_par,a,b,j,st); // now iterate and refine initial value
		s = (4.0*st-ost)/3.0;
		if (j > 5) {
			if (fabs(s-os) < TENTOMINUSSIX*fabs(os) || (s == 0.0 && os == 0.0)) {
				return s;
			}
		}
		os = s;
		ost = st;
	}
#ifdef PRINT_DEBUG
	fprintf(stderr, "operaMath::LengthofPolynomial: Too many steps (%d > %d) in routine LengthofPolynomial\n", j, JMAX);
#endif
	return NAN;
	
}

double PolyLengthFunc(double x, const double *p, int n_par)
{	
	return sqrt(1 + DiffPolynomialFunction(x,p,n_par)*DiffPolynomialFunction(x,p,n_par));
}

/*
 * DT May 14 2013 - Note the nasty side-effect here. the variable "s" MUST be static and
 * this function MUST be called with n == 1 first so that s gets initialized and then
 * is usable on line 133. Not good coding...
 */
double PolyLengthTrapezoid(const double *p, int n_par, double a, double b, int n, double initials)
{
	double x,tnm,sum,del,s;
	int it,j;
	
	if (n == 1) {
		s = 0.5*(b-a)*(PolyLengthFunc(a,p,n_par)+PolyLengthFunc(b,p,n_par));
	} else {
		for (it=1,j=1; j<n-1; j++) {
			it <<= 1;
		}
		tnm = it;
		del = (b-a)/tnm;
		x = a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
			sum += PolyLengthFunc(x,p,n_par);
		}
		s = 0.5*(initials+(b-a)*sum/tnm);	// NAB the s is used unitialized but it is static!
	}
	return s;
}

/* 
 * double PlanckLaw(double Temperature, double Wavelength)
 * \brief This function returns the spectral radiance of a black body in W/(m^2 dlambda).
 * \param Temperature is a float input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a float input that represents the wavelength at which the black body is observed in nanometers
 * \return float value for the spectral radiance of the black body in W/(m^2 dlambda)
 */
double PlanckLaw(double Temperature, double Wavelength)
{
    double SpectralRadiance;
    
    SpectralRadiance = (2.0 * M_PI * PLANCK_CONSTANT * pow(SPEED_OF_LIGHT_M,2)) / (pow(Wavelength * 1.0e-9,5) * (exp(PLANCK_CONSTANT * SPEED_OF_LIGHT_M / (Wavelength * 1.0e-9 * BOLTZMANN_CONSTANT * Temperature)) - 1.0));
    
    return SpectralRadiance;
}
