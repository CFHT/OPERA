/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
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
 * ladfit
 * \author Doug Teeple
 * \brief This class implements least average deviation fit.
 * \file ladfit.c
 * \ingroup libraries
 */

#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/ladfit.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaStats.h"

/* 
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

void ladfitWithError(float x[], float y[], int nX, float *a, float *aError, float *b, float *bError, float *absdev)
{
	double aa, bb, f, f1, f2, b1, b2, del, delb, chisqr = 0.0, sigb, sigy;
	double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
 	double sumsqrX = 0.0, sumX = 0.0, delta;
	float *dev = (float *)malloc(nX*sizeof(float));
	float *d = (float *)malloc(nX*sizeof(float));
	float *dm = (float *)malloc(nX*sizeof(float));
	float eps = EPS;
	
	for (int i=0; i<nX; i++) {
		sx += (double)x[i];
		sy += (double)y[i];
		sxy += (double)x[i]*(double)y[i];
		sxx += (double)x[i]*(double)x[i];
	}
	
	del = (nX * sxx) -( sx * sx);
	
	if (del == 0.0) {					// All X's are the same
		*b = operaArrayMedian(nX, y);
		*a = 0.0; 						// Bisect the range w/ a flat line
		return;
	}
	
	aa = (sxx * sy - sx * sxy) / del;	// Least squares solution y = x * aa + bb
	bb = (nX * sxy - sx * sy) / del;
	
	for (int i=0; i<nX; i++) {
		register double t = (double)y[i] - (aa + (bb * (double)x[i]));
		chisqr += t * t;
	}
	sigb = sqrt(chisqr / del);			// Standard deviation
	
	b1 = bb;
	f1 = MedianFunction(d, dm, b1, x, y, &aa, absdev, nX, eps);
	
	//  Quick return. The initial least squares gradient is the LAD solution.
	if (f1 == 0.0) {
		goto done;
	}
	
	delb = ((f1 >= 0) ? 3.0 : -3.0) * sigb;
	
	b2 = b1 + delb;
	f2 = MedianFunction(d, dm, b2, x, y, &aa, absdev, nX, eps);
	while (f1*f2 > 0.0) {     // Bracket the zero of the function
		b1 = b2;
		f1 = f2;
		b2 = b1 + delb;
		f2 = MedianFunction(d, dm, b2, x, y, &aa, absdev, nX, eps);
	}
	
	//  In case we finish early.
	bb = b2;
	f = f2;
	
	// Narrow tolerance to refine 0 of fcn.
	sigb = 0.01 * sigb;
	while (fabs(b2-b1) > sigb && f != 0.0) { // bisection of interval b1,b2.
		bb = 0.5 * (b1 + b2);
		if (bb == b1 || bb == b2)
			break;
		f = MedianFunction(d, dm, bb, x, y, &aa, absdev, nX, eps);
		if (f*f1 >= 0.0) {
			f1 = f;
			b1 = bb;
		} else {
			f2 = f;
			b2 = bb;
		}
	}
	
done:
	
	for (int i=0; i<nX; i++) {
		dev[i] = fabs(y[i] - ((float)aa + ((float)bb * x[i])));
		sumsqrX += x[i] * x[i];
		sumX += x[i];
	}
	sigy = operaArrayMedianQuick(nX, dev) /  0.674433; // quick is OK since dev is a copy
	delta = (nX * sumsqrX) - (sumX * sumX);
	*absdev = *absdev / nX;
	*aError = sigy * sqrt(sumsqrX / delta);
	*bError = sigy * sqrt(nX / delta);
	*b = bb;
	*a = aa;
	free(dev);
}

void ladfit(float x[], float y[], int nX, float *a, float *b, float *absdev)
{
	double aa, bb, f, f1, f2, b1, b2, del, delb, chisqr = 0.0, sigb;
	double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
	float *d = (float *)malloc(nX*sizeof(float));
	float *dm = (float *)malloc(nX*sizeof(float));
	float eps = EPS;
	
	for (int i=0; i<nX; i++) {
		sx += (double)x[i];
		sy += (double)y[i];
		sxy += (double)x[i]*(double)y[i];
		sxx += (double)x[i]*(double)x[i];
	}
	
	del = (nX * sxx) -( sx * sx);
	
	if (del == 0.0) {					// All X's are the same
		*b = operaArrayMedian(nX, y);
		*a = 0.0; 						// Bisect the range w/ a flat line
		return;
	}
	
	aa = (sxx * sy - sx * sxy) / del;	// Least squares solution y = x * aa + bb
	bb = (nX * sxy - sx * sy) / del;
	
	for (int i=0; i<nX; i++) {
		register double t = (double)y[i] - (aa + (bb * (double)x[i]));
		chisqr += t * t;
	}
	sigb = sqrt(chisqr / del);			// Standard deviation
	
	b1 = bb;
	f1 = MedianFunction(d, dm, b1, x, y, &aa, absdev, nX, eps);
	
	//  Quick return. The initial least squares gradient is the LAD solution.
	if (f1 == 0.0) {
		goto done;
	}
	
	delb = ((f1 >= 0) ? 3.0 : -3.0) * sigb;
	
	b2 = b1 + delb;
	f2 = MedianFunction(d, dm, b2, x, y, &aa, absdev, nX, eps);
	while (f1*f2 > 0.0) {     // Bracket the zero of the function
		b1 = b2;
		f1 = f2;
		b2 = b1 + delb;
		f2 = MedianFunction(d, dm, b2, x, y, &aa, absdev, nX, eps);
	}
	
	//  In case we finish early.
	bb = b2;
	f = f2;
	
	// Narrow tolerance to refine 0 of fcn.
	sigb = 0.01 * sigb;
	while (fabs(b2-b1) > sigb && f != 0.0) { // bisection of interval b1,b2.
		bb = 0.5 * (b1 + b2);
		if (bb == b1 || bb == b2)
			break;
		f = MedianFunction(d, dm, bb, x, y, &aa, absdev, nX, eps);
		if (f*f1 >= 0.0) {
			f1 = f;
			b1 = bb;
		} else {
			f2 = f;
			b2 = bb;
		}
	}
done:   	
	*absdev = *absdev / nX;
	*b = bb;
	*a = aa;
}

/* 
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

void ladfitWithError_d(const double x[], const double y[], int nX, double *a, double *aError, double *b, double *bError, double *absdev)
{
	double aa, bb, f, f1, f2, b1, b2, del, delb, chisqr = 0.0, sigb, sigy;
	double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
 	double sumsqrX = 0.0, sumX = 0.0, delta;
	double *dev = (double *)malloc(nX*sizeof(double));
	double *d = (double *)malloc(nX*sizeof(double));
	double *dm = (double *)malloc(nX*sizeof(double));
	float eps = EPS;
	
	for (int i=0; i<nX; i++) {
		sx += x[i];
		sy += y[i];
		sxy += x[i]*y[i];
		sxx += x[i]*x[i];
	}
	
	del = (nX * sxx) -( sx * sx);
	
	if (del == 0.0) {					// All X's are the same
		*b = operaArrayMedian_d(nX, y);
		*a = 0.0; 						// Bisect the range w/ a flat line
		return;
	}
	
	aa = (sxx * sy - sx * sxy) / del;	// Least squares solution y = x * aa + bb
	bb = (nX * sxy - sx * sy) / del;
	
	for (int i=0; i<nX; i++) {
		register double t = y[i] - (aa + (bb * x[i]));
		chisqr += t * t;
	}
	sigb = sqrt(chisqr / del);			// Standard deviation
	
	b1 = bb;
	f1 = MedianFunction_d(d, dm, b1, x, y, &aa, absdev, nX, eps);
	
	//  Quick return. The initial least squares gradient is the LAD solution.
	if (f1 == 0.0) {
		goto done;
	}
	
	delb = ((f1 >= 0) ? 3.0 : -3.0) * sigb;
	
	b2 = b1 + delb;
	f2 = MedianFunction_d(d, dm, b2, x, y, &aa, absdev, nX, eps);
	while (f1*f2 > 0.0) {     // Bracket the zero of the function
		b1 = b2;
		f1 = f2;
		b2 = b1 + delb;
		f2 = MedianFunction_d(d, dm, b2, x, y, &aa, absdev, nX, eps);
	}
	
	//  In case we finish early.
	bb = b2;
	f = f2;
	
	// Narrow tolerance to refine 0 of fcn.
	sigb = 0.01 * sigb;
	while (fabs(b2-b1) > sigb && f != 0.0) { // bisection of interval b1,b2.
		bb = 0.5 * (b1 + b2);
		if (bb == b1 || bb == b2)
			break;
		f = MedianFunction_d(d, dm, bb, x, y, &aa, absdev, nX, eps);
		if (f*f1 >= 0.0) {
			f1 = f;
			b1 = bb;
		} else {
			f2 = f;
			b2 = bb;
		}
	}
	
done:
	
	for (int i=0; i<nX; i++) {
		dev[i] = fabs(y[i] - (aa + (bb * x[i])));
		sumsqrX += x[i] * x[i];
		sumX += x[i];
	}
	sigy = operaArrayMedianQuick_d(nX, dev) /  0.674433; // quick is OK since dev is a copy
	delta = (nX * sumsqrX) - (sumX * sumX);
	*absdev = *absdev / nX;
	*aError = sigy * sqrt(sumsqrX / delta);
	*bError = sigy * sqrt(nX / delta);
	*b = bb;
	*a = aa;
	free(dev);
}

void ladfit_d(const double x[], const double y[], int nX, double *a, double *b, double *absdev)
{
	double aa, bb, f, f1, f2, b1, b2, del, delb, chisqr = 0.0, sigb;
	double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;
	double *d = (double *)malloc(nX*sizeof(double));
	double *dm = (double *)malloc(nX*sizeof(double));
	float eps = EPS;
	
	for (int i=0; i<nX; i++) {
		sx += x[i];
		sy += y[i];
		sxy += x[i]*y[i];
		sxx += x[i]*x[i];
	}
	
	del = (nX * sxx) -( sx * sx);
	
	if (del == 0.0) {					// All X's are the same
		*b = operaArrayMedian_d(nX, y);
		*a = 0.0; 						// Bisect the range w/ a flat line
		return;
	}
	
	aa = (sxx * sy - sx * sxy) / del;	// Least squares solution y = x * aa + bb
	bb = (nX * sxy - sx * sy) / del;
	
	for (int i=0; i<nX; i++) {
		register double t = y[i] - (aa + (bb * x[i]));
		chisqr += t * t;
	}
	sigb = sqrt(chisqr / del);			// Standard deviation
	
	b1 = bb;
	f1 = MedianFunction_d(d, dm, b1, x, y, &aa, absdev, nX, eps);
	
	//  Quick return. The initial least squares gradient is the LAD solution.
	if (f1 == 0.0) {
		goto done;
	}
	
	delb = ((f1 >= 0) ? 3.0 : -3.0) * sigb;
	
	b2 = b1 + delb;
	f2 = MedianFunction_d(d, dm, b2, x, y, &aa, absdev, nX, eps);
	while (f1*f2 > 0.0) {     // Bracket the zero of the function
		b1 = b2;
		f1 = f2;
		b2 = b1 + delb;
		f2 = MedianFunction_d(d, dm, b2, x, y, &aa, absdev, nX, eps);
	}
	
	//  In case we finish early.
	bb = b2;
	f = f2;
	
	// Narrow tolerance to refine 0 of fcn.
	sigb = 0.01 * sigb;
	while (fabs(b2-b1) > sigb && f != 0.0) { // bisection of interval b1,b2.
		bb = 0.5 * (b1 + b2);
		if (bb == b1 || bb == b2)
			break;
		f = MedianFunction_d(d, dm, bb, x, y, &aa, absdev, nX, eps);
		if (f*f1 >= 0.0) {
			f1 = f;
			b1 = bb;
		} else {
			f2 = f;
			b2 = bb;
		}
	}
done:   	
	*absdev = *absdev / nX;
	*b = bb;
	*a = aa;
}

