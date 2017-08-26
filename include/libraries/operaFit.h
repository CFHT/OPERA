#ifndef OPERAFIT_H
#define OPERAFIT_H

/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFit
 Version: 1.0
 Description: This C library implements LM curve fitting for commonly 
 used functions such as polynomials, gaussians, etc.  
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
 * operaFit 
 * \author Eder Martioli
 * \brief Curve fitting routines in C.
 * \file operaFit.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
#include "libraries/operaLibCommon.h"	// definition of a CMatrix
	
	/*
	 * The functions below uses operaLMFit library (LMFIT)
	 */
	void operaLMFitPolynomial(unsigned m_dat, const double *x, const double *y, int n_par, double par[], double *chi2);
	void operaLMFitLaurentPolynomial(unsigned m_dat, double *x, double *y, int minn_par, int maxn_par, double *par, double *chi2);
	void operaLMFitGaussian(unsigned m_dat, double *x, double *y, double *a, double *x0, double *sig, double *chi2);
	void operaLMFitMultipleGaussian(unsigned m_dat, double *x, double *y, unsigned ngauss, double a[], double x0[], double sig[], double *chi2);
	void operaLMFit2DPolynomial(unsigned m_dat, double *x, double *y, double *fxy, int n_par, double par[], double *chi2);
	void operaLMFit2DGaussian(unsigned m_dat, double *x, double *y, double *fxy, double *a, double *x0, double *y0, double *sigx, double *sigy, double *chi2);

	/*
	 * The functions below uses mpfit library (MPFIT)
	 */
	int MPPolyFunc(int m, int n, double *p, double *dy, double **dvec, void *vars);
	int operaMPFitPolynomial(unsigned m_dat, const double *x, const double *y, const double *ey, int n_par, double *par, double *epar, double *chi2);
	int MPGaussFunc(int m, int n, double *p, double *dy, double **dvec, void *vars);
	int	operaMPFitGaussian(unsigned m_dat, double *x, double *y, double *ey, double *a, double *ea, double *x0, double *ex0, double *sig, double *esig, double *chi2);
    int	operaMPFitMultipleGaussian(unsigned m_dat, double *x, double *y, double *ey, unsigned ngauss, double a[], double ea[], double x0[], double ex0[], double sig[], double esig[], double *chi2);
    int MPGaussFuncWithBaseline(int m, int n, double *p, double *dy, double **dvec, void *vars);	
	/*
	 * The functions below uses nr library, although all necessary functions to make them work are in operaFit.c
	 */
	// spline and 2D-spline interpolation and related functions
	void operaFitSpline(unsigned nin, const float *xin, const float *yin, unsigned nout, const float *xout, float *yout);
	void operaFitSplineDouble(unsigned nin, const double *xin, const double *yin, unsigned nout, const double *xout, double *yout);
	int cubicspline(const float *x, const float *y, unsigned n, float yp1, float ypn, float *y2);
	int cubicsplineDouble(const double *x, const double *y, unsigned n, double yp1, double ypn, double *y2);
	int splineinterpolate(const float *xa, const float *ya, const float *y2a, int n, float x, float *y);
	int splineinterpolateDouble(const double *xa, const double *ya, const double *y2a, int n, double x, double *y);
	void operaFit2DSpline(unsigned nxin, float *xin, unsigned nyin, float *yin, float *fxyin, unsigned nxout, float *xout, unsigned nyout, float *yout, float *fxyout);
	void cubicspline2D(float *x1a, float *x2a, CMatrix ya, int m, int n, CMatrix y2a);
	void splineinterpolate2D(float *x1a, float *x2a, CMatrix ya, CMatrix y2a, int m, int n, float x1, float x2, float *y);
	
	// Least Average Deviation Fit
	void ladfit(float x[], float y[], int nX, float *a, float *b, float *absdev);
	void ladfitWithError(float x[], float y[], int nX, float *a, float *aError, float *b, float *bError, float *absdev);
	
	// wrappers
	double PolynomialFunction(double x, const double *p, int n_par);
    double LaurentPolynomialFunction(double x, const double *p, int n_par);
   // double LaurentPolynomialFunction(double x, const double *p, int minn_par, int maxn_par); // this format is not possible due to limitation of LMFit library
	double GaussianFunction(double x, const double *p, int n_totalpar);
	double Polynomial2DFunction(double x, double y, const double *p, int n_par);
	double Gaussian2DFunction(double x, double y, const double *p, int n_par);
    double GaussianWithBaseline(double x, const double *p, int n_totalpar);
#ifdef __cplusplus
}
#endif


#endif
