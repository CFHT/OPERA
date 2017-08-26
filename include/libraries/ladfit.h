#ifndef LADFIT_H
#define LADFIT_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: ladfit
 Version: 1.0
 Description: Least Avergae Deviation FIT of a first order polynomial.
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

/*! 
 * \brief least average deviation fit,  in C.
 * \file ladfit.h
 * \ingroup libraries
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "libraries/operaLibCommon.h"
#include "libraries/operaStats.h"

	/*! 
	 * static inline float MedianFunction(float *d, float *dm, double b, float x[], float y[], double *a, float *absdev, int size, float eps)
	 */
	static inline float MedianFunction(float *d, float *dm, double b, float x[], float y[], double *a, float *absdev, int size, float eps)
	{
		float sum = 0.0, adev = 0.0;
		int i, j = 0;
		
		for (i=0; i<size; i++) {
			dm[i] = y[i] - ((float)b * x[i]);									// y - b*x
		}
		*a = operaArrayMedian(size, dm);									// a = MEDIAN(y - b*x, /EVEN)
		for (i=0; i<size; i++) {
			d[i] = y[i] - (((float)b * x[i]) + (float)*a);						// d = y - (b * x + a)
			adev += fabs(d[i]);
			if ( y[i] != 0.0 )
				d[i] /= fabs(y[i]);
			if (fabs(d[i]) > eps) {
				sum += (d[i] >= 0.0 ? x[i] : -x[i]);
				j++;
			}
		}
		
		*absdev = adev;
		return sum;
	}
	
	/*! 
	 * static inline float MedianFunction_d(double *d, double *dm, double b, double x[], double y[], double *a, double *absdev, int size, double eps)
	 */
	static inline float MedianFunction_d(double *d, double *dm, double b, const double x[], const double y[], double *a, double *absdev, int size, double eps)
	{
		double sum = 0.0, adev = 0.0;
		int i, j = 0;
		
		for (i=0; i<size; i++) {
			dm[i] = y[i] - (b * x[i]);											// y - b*x
		}
		*a = operaArrayMedian_d(size, dm);									// a = MEDIAN(y - b*x, /EVEN)
		for (i=0; i<size; i++) {
			d[i] = y[i] - ((b * x[i]) + *a);									// d = y - (b * x + a)
			adev += fabs(d[i]);
			if ( y[i] != 0.0 )
				d[i] /= fabs(y[i]);
			if (fabs(d[i]) > eps) {
				sum += (d[i] >= 0.0 ? x[i] : -x[i]);
				j++;
			}
		}
		
		*absdev = adev;
		return sum;
	}
	
	/*! 
 * void ladfit(float x[], float y[], int nX, float *absdev, float *slope, float *xoffset)
 * \brief This function fits the paired data {X(i), Y(i)} to the linear model,
 * \brief y = A + Bx, using a "robust" least absolute deviation method. The 
 * \brief result is a two-element vector containing the model parameters, slope
 * \brief and xoffset.
 * \param x An n-element vector of type float.
 * \param y An n-element vector of typef loat.
 * \param absdev the mean absolute deviation for each data-point in the y-direction.
 * \param slope the B in Bx
 * \param offset the A in equation above
 * \verbatim
 *       Define two n-element vectors of paired data.
 *         x = [-3.20, 4.49, -1.66, 0.64, -2.43, -0.89, -0.12, 1.41, $
 *               2.95, 2.18,  3.72, 5.26]
 *         y = [-7.14, -1.30, -4.26, -1.90, -6.19, -3.98, -2.87, -1.66, $
 *              -0.78, -2.61,  0.31,  1.74]
 *       Compute the model parameters, A and B.
 *         result = ladfit(x, y, absdev = absdev)
 *       The result should be the two-element vector:
 *         [-3.15301, 0.930440]
 *       The keyword parameter should be returned as:
 *         absdev = 0.636851
 * \endverbatim
 * \return void
 */

	void ladfitWithError(float x[], float y[], int nX, float *a, float *aError, float *b, float *bError, float *absdev);
	
	void ladfit(float x[], float y[], int nX, float *a, float *b, float *absdev);
	
	void ladfitWithError_d(const double x[], const double y[], int nX, double *a, double *aError, double *b, double *bError, double *absdev);
	
	void ladfit_d(const double x[], const double y[], int nX, double *a, double *b, double *absdev);
	
#ifdef __cplusplus
}
#endif
	
#endif
