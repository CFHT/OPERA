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

#include <math.h>

/*!
 * \brief Fitting library.
 * \file operaFit.c
 */

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaLMFit.h"
#include "libraries/operaFit.h"
#include "libraries/mpfit.h"

/*!
 * operaFit
 * \author Eder Martioli
 * \brief Fitting library.
 * \details {This library contains the routines for curve fitting using the Levenberg-Marquardt 
 * Least Squares method for usual functions such as polynomials, gaussians, etc.}
 * \ingroup libraries
 */

#undef PRINT_DEBUG

/* 
 * double PolynomialFunction(double x, const double *p, int n_par)
 * \brief This function returns the value of a given polynomial function.
 * \param x is a double input value for which the given polynomial is evaluated  
 * \param p is a const double pointer that contains the coefficients of the given polynomial
 * \param n_par is an int that defines the order of the polynomial
 * \return double value of the polynomial 
 */
double PolynomialFunction(double x, const double *p, int n_par)
{
	double total = 0;
	for(unsigned i = n_par; i > 0; i--) total = x*total + p[i-1];
	return total;
}

/* 
 * void operaLMFitPolynomial(unsigned m_dat, double *x, double *y, int n_par, double *par, double *chi2)
 * \brief This function performs Polynomial Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param n_par is an int for the order of the polynomial to be fit
 * \param par is a double pointer for the coefficients of the polynomial
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
void operaLMFitPolynomial(unsigned m_dat, const double *x, const double *y, int n_par, double *par, double *chi2) 
{
	int i;  
	/* auxiliary parameters */
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
	
	/* perform the fit */
	lmcurve_fit( n_par, par, m_dat, x, y, PolynomialFunction, &control, &status);
#ifdef PRINT_DEBUG
	printf( "\n *** Results:\n" );
	printf( "Status after %d function evaluations:\n  %s\n",
			status.nfev, lm_infmsg[status.info] );
#endif	
	double WSSR,fcalc;
	WSSR = 0;
	
	for(i=0;i<m_dat;i++)
	{
		fcalc =  PolynomialFunction(x[i],par,n_par);
		WSSR = WSSR + pow((fcalc - y[i]),2.0);
	}   
	
	*chi2 = WSSR/(double)(m_dat-n_par);
#ifdef PRINT_DEBUG	
	double rms = sqrt(WSSR/(double)(m_dat - n_par));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-n_par));    
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-n_par));	 
	printf("\n *** Obtained parameters:\n");
	for ( i = 0; i < n_par; ++i)
		printf("  par[%i] = %12g\n", i, par[i]);
	printf("\n *** Obtained norm:\n  %12g\n", status.fnorm );
	printf("\n *** Fit polynomial: \n");    	
	printf("f(x) =");
	for(i=0;i<n_par;i++) {
		if(i==0)
			printf(" %lf",par[i]);
		else if (i==1)
			printf(" + %lf*x",par[i]);
		else
			printf(" + %lf*x**(%d)",par[i],i);		
	}    	
	printf("\n");
	printf("\n *** Data:\n");	
	printf("x\ty\tf(x)\ty-f(x)\n");
	for (i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i], PolynomialFunction(x[i],par,n_par), y[i] - PolynomialFunction(x[i],par,n_par) );
#endif
}

/*
 * double LaurentPolynomialFunction(double x, const double *p, int minn_par, int maxn_par)
 * \brief This function returns the value of a given Laurent polynomial function.
 * \param x is a double input value for which the given polynomial is evaluated
 * \param p is a const double pointer that contains the coefficients of the given polynomial
 * \param minn_par is an int that defines the lowest order of the polynomial
 * \param maxn_par is an int that defines the highest order of the polynomial
 * \return double value of the polynomial
 */
// The format below would be preferable, however it is not possible due to a limitation of LMFit lib
//double LaurentPolynomialFunction(double x, const double *p, int minn_par, int maxn_par)
double LaurentPolynomialFunction(double x, const double *p, int n_par)
{
    // Below is a workaround in order to calculate the Laurent Polynomial using only one
    // parameter n_par. The problem is that we have to assume the polynomial only has
    // terms with negative power.
    int maxn_par = 0;
    // Since  n_par = (maxn_par - minn_par + 1)
    int minn_par = 1 - n_par;
    
    double fpoly = 0;
    
	for (int i=minn_par; i<=maxn_par; i++) {
        if(i==0) {
            fpoly += p[(unsigned)(i-minn_par)];
        } else {
            fpoly += p[(unsigned)(i-minn_par)]*pow(x, (double)i);
        }
	}
	return fpoly;
}

/*
 * void operaLMFitLaurentPolynomial(unsigned m_dat, double *x, double *y, int minn_par, int maxn_par, double *par, double *chi2)
 * \brief This function performs Polynomial Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param minn_par is an int that defines the lowest order of the polynomial
 * \param maxn_par is an int that defines the highest order of the polynomial
 * \param par is a double pointer for the coefficients of the polynomial
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
void operaLMFitLaurentPolynomial(unsigned m_dat, double *x, double *y, int minn_par, int maxn_par, double *par, double *chi2)
{
	int i;
	/* auxiliary parameters */
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
	
    int n_par = (maxn_par - minn_par) + 1;
    
    if(maxn_par!= 0 || minn_par >0){
        printf("operaLMFitLaurentPolynomial: Error : maxn_par must be 0. This is due to a limitation of LMFit library.");
        return;
    }
    
	/* perform the fit */
	lmcurve_fit(n_par, par, m_dat, x, y, LaurentPolynomialFunction, &control, &status);
    
#ifdef PRINT_DEBUG
	printf( "\n *** Results:\n" );
	printf( "Status after %d function evaluations:\n  %s\n",
           status.nfev, lm_infmsg[status.info] );
#endif
	double WSSR,fcalc;
	WSSR = 0;
	
	for(i=0;i<m_dat;i++)
	{
		fcalc =  LaurentPolynomialFunction(x[i],par,n_par);
		WSSR = WSSR + pow((fcalc - y[i]),2.0);
	}
	
	*chi2 = WSSR/(double)(m_dat-n_par);
#ifdef PRINT_DEBUG
	double rms = sqrt(WSSR/(double)(m_dat - n_par));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-n_par));
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-n_par));
	printf("\n *** Obtained parameters:\n");
	for ( i = 0; i < n_par; ++i)
		printf("  par[%i] = %12g\n", i, par[i]);
	printf("\n *** Obtained norm:\n  %12g\n", status.fnorm );
	printf("\n *** Fit Laurent polynomial: \n");
	printf("f(x) =");
	for(int i=minn_par;i<maxn_par;i++) {
		if(i==0)
			printf(" %lf",par[i]);
		else if (i==1)
			printf(" + %lf*x",par[i]);
		else
			printf(" + %lf*x**(%d)",par[i],i);
	}
	printf("\n");
	printf("\n *** Data:\n");
	printf("x\ty\tf(x)\ty-f(x)\n");
	for (i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i], LaurentPolynomialFunction(x[i],par,n_par), y[i] - LaurentPolynomialFunction(x[i],par,n_par) );
#endif
}


/* This is the private data structure which contains the data points
 and their uncertainties */
struct vars_struct {
  double *x;
  double *y;
  double *ey;
};

/* 
 * void operaMPFitPolynomial(unsigned m_dat, double *x, double *y, double *ey, int n_par, double *par, double *epar, double *chi2) 
 * \brief This function performs Polynomial Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param ey is a double pointer with y-error data values 
 * \param n_par is an int for the order of the polynomial to be fit
 * \param par is a double pointer for the coefficients of the polynomial
 * \param par is a double pointer for the errors on the coefficients of the polynomial
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
int operaMPFitPolynomial(unsigned m_dat, const double *x, const double *y, const double *ey, int n_par, double *par, double *epar, double *chi2) 
{
  struct vars_struct v;
  int status;
  
  mp_result* resultptr = 0;
  if(epar) {
    mp_result result;
    memset(&result,0,sizeof(result));       /* Zero results structure */
    result.xerror = epar;
    resultptr = &result;
  }
  
  mp_par pars[MAXPOLYNOMIAL];
  memset(pars, 0, sizeof(pars));       /* Initialize constraint structure */
  
  v.x = x;
  v.y = y;
  v.ey = ey;
  
  /* Call fitting function for 10 data points and 2 parameters */
  status = mpfit(MPPolyFunc, m_dat, n_par, par, pars, 0, (void *) &v, resultptr);
  
  if(chi2 && resultptr) *chi2 = resultptr->bestnorm/(double)(m_dat-n_par);
  
  return status;
}

int MPPolyFunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey, f;
	
  x = v->x;
  y = v->y;
  ey = v->ey;
	
  for (i=0; i<m; i++) {
    f = PolynomialFunction(x[i],p,n);     /* Linear fit function; note f = a - b*x */
    dy[i] = (y[i] - f)/ey[i];
  }
	
  return 0;
}

/* 
 * void operaMPFitGaussian(unsigned m_dat, double *x, double *y, double *ey, double *a, double *ea, double *x0, double *ex0, double *sig, double *esig, double *chi2)
 * \brief This function performs Gaussian Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param ey is a double pointer with y-error data values 
 * \param a is a double pointer for the amplitude of the Gaussian
 * \param ea is a double pointer for the error on the amplitude of the Gaussian 
 * \param x0 is a double pointer for the center of the Gaussian
 * \param ex0 is a double pointer for the error on the center of the Gaussian 
 * \param sig is a double pointer for the spread of the Gaussian
 * \param esig is a double pointer for the error on the spread of the Gaussian 
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
int	operaMPFitGaussian(unsigned m_dat, double *x, double *y, double *ey, double *a, double *ea, double *x0, double *ex0, double *sig, double *esig, double *chi2)
{
	int n_par = 3;
	double *par, *epar;
	par = (double*) malloc(n_par * sizeof(double));	
	epar = (double*) malloc(n_par * sizeof(double));	
	
	par[0] = *a;
	par[1] = *x0;
	par[2] = *sig;

	epar[0] = *ea;
	epar[1] = *ex0;
	epar[2] = *esig;
	
    mp_par pars[3];
    memset(pars, 0, sizeof(pars));       /* Initialize constraint structure */
    
  struct vars_struct v;
  int status;
  mp_result result;
	
  memset(&result,0,sizeof(result));       /* Zero results structure */
	
  result.xerror = epar;
	
  v.x = x;
  v.y = y;
  v.ey = ey;
	
  /* Call fitting function for 10 data points and 2 parameters */
  status = mpfit(MPGaussFunc, m_dat, n_par, par, pars, 0, (void *) &v, &result);
	
	*a = par[0];
	*x0 = par[1];
	*sig = par[2];
	
	*ea = result.xerror[0];
	*ex0 = result.xerror[1];
	*esig =	result.xerror[2];
	
	*chi2 = result.bestnorm/(double)(m_dat-n_par);	
	
	free(par);
	free(epar);
	return status;
}

int MPGaussFunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey, f;
	
  x = v->x;
  y = v->y;
  ey = v->ey;
	
  for (i=0; i<m; i++) {
    f = GaussianFunction(x[i],p,n);     /* Linear fit function; note f = a - b*x */
    dy[i] = (y[i] - f)/ey[i];
  }
	
  return 0;
}

int MPGaussFuncWithBaseline(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
    int i;
    struct vars_struct *v = (struct vars_struct *) vars;
    double *x, *y, *ey, f;
	
    x = v->x;
    y = v->y;
    ey = v->ey;
	
    for (i=0; i<m; i++) {
        f = GaussianWithBaseline(x[i],p,n);     /* Linear fit function; note f = a - b*x */
        dy[i] = (y[i] - f)/ey[i];
    }
	
    return 0;
}

/* 
 * int operaMPFitMultipleGaussian(unsigned m_dat, double *x, double *y, double *ey, unsigned ngauss, double a[], double ea[], double x0[], double ex0[], double sig[], double esig[], double *chi2)
 * \brief This function performs Multiple Gaussian Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param ey is a double pointer with y-error data values 
 * \param naguss is an unsigned for the number of gaussians to fit the data 
 * \param a is a double pointer for the amplitude(s) of the Gaussian(s)
 * \param ea is a double pointer for the error(s) on the amplitude(s) of the Gaussian(s) 
 * \param x0 is a double pointer for the center(s) of the Gaussian(s)
 * \param ex0 is a double pointer for the error(s) on the center(s) of the Gaussian(s) 
 * \param sig is a double pointer for the spread(s) of the Gaussian(s)
 * \param esig is a double pointer for the error(s) on the spread(s) of the Gaussian(s) 
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
int	operaMPFitMultipleGaussian(unsigned m_dat, double *x, double *y, double *ey, unsigned ngauss, double a[], double ea[], double x0[], double ex0[], double sig[], double esig[], double *chi2)
{
	int n_par = 3;
    int totalnpar = n_par*(int)ngauss;
	double *par, *epar;
	par = (double*) malloc(totalnpar * sizeof(double));	
    epar = (double*) malloc(totalnpar * sizeof(double));	
    
    mp_par pars[MAXGAUSSCOEFFS];
    memset(pars, 0, sizeof(pars));       /* Initialize constraint structure */
    	
    for(unsigned i=0;i<ngauss;i++) {
        par[0+i*(unsigned)n_par] = a[i];
        par[1+i*(unsigned)n_par] = x0[i];
        par[2+i*(unsigned)n_par] = sig[i];
    
        epar[0+i*(unsigned)n_par] = ea[i];
        epar[1+i*(unsigned)n_par] = ex0[i];
        epar[2+i*(unsigned)n_par] = esig[i];
	}

    struct vars_struct v;
    int status;
    mp_result result;
	
    memset(&result,0,sizeof(result));       /* Zero results structure */
	
    result.xerror = epar;
	
    v.x = x;
    v.y = y;
    v.ey = ey;
    
    /* Call fitting function for 10 data points and 2 parameters */
    status = mpfit(MPGaussFunc, m_dat, totalnpar, par, pars, 0, (void *) &v, &result);
	
    for(unsigned i=0;i<ngauss;i++) {    
        a[i] = par[0+i*(unsigned)n_par];
        x0[i] = par[1+i*(unsigned)n_par];
        sig[i] = par[2+i*(unsigned)n_par];
	
        ea[i] = result.xerror[0+i*(unsigned)n_par];
        ex0[i] = result.xerror[1+i*(unsigned)n_par];
        esig[i] = result.xerror[2+i*(unsigned)n_par];
//        printf("%lf\t%lf\t%lf\t",par[0+i*(unsigned)n_par],par[1+i*(unsigned)n_par],par[2+i*(unsigned)n_par]);
	}
	*chi2 = result.bestnorm/(double)(m_dat-totalnpar);	
	
	free(par);
	free(epar);
	return status;
}

/* 
 * double GaussianFunction(double x, const double *p, int n_totalpar)
 * \brief This function returns the value for a sum of Gaussians: 
 * \brief f(x) = A1*exp(-((x-x01)/sqrt(2)*sig1)^2) + A2*exp(-((x-x02)/sqrt(2)*sig2)^2) + ...
 * \param x is a double input value for which the given Gaussian is evaluated  
 * \param p is a const double pointer for the coefficients of the Gaussian
 * \param n_par is an int for the number of coefficients (must be multiple of 3!)
 * \return double value of Gaussian
 */
double GaussianFunction(double x, const double *p, int n_totalpar)
{	
    int n_par = 3;
    
    int NumberOfGaussians = n_totalpar/n_par;
    
    double gaussfunc = 0;
    
    for(int i=0; i<NumberOfGaussians; i++){
        double amp = p[0+i*n_par];
        double x0 = p[1+i*n_par];
        double sig = p[2+i*n_par];
        
        gaussfunc += amp*exp(-(x-x0)*(x-x0)/(2.0*sig*sig));
    }
    
    return gaussfunc;
}

/* 
 * double GaussianWithBaseline(double x, const double *p, int n_totalpar)
 * \brief This function returns the value for a sum of Gaussians plus a linear trend: 
 * \brief f(x) = A1*exp(-((x-x01)/sqrt(2)*sig1)^2) + A2*exp(-((x-x02)/sqrt(2)*sig2)^2) + ... + (a + b*x)
 * \param x is a double input value for which the given Gaussian is evaluated  
 * \param p is a const double pointer for the coefficients of the Gaussian and the slope
 * \param n_totalpar is an int for the number of coefficients (n_totalpar-2 must be multiple of 3!)
 * \return double value of Gaussian
 */
double GaussianWithBaseline(double x, const double *p, int n_totalpar)
{	
    int n_par = 3;
    
    int NumberOfGaussians = (n_totalpar - 2)/n_par;
    
    double gaussfunc = 0;
    
    for(int i=0; i<NumberOfGaussians; i++){
        double amp = p[0+i*n_par];
        double x0 = p[1+i*n_par];
        double sig = p[2+i*n_par];
        
        gaussfunc += amp*exp(-(x-x0)*(x-x0)/(2.0*sig*sig));
    }
    
    double baselineIntercept = p[n_totalpar - 2];
    double baselineSlope =  p[n_totalpar - 1]; 
    
    gaussfunc += baselineIntercept + baselineSlope*x;
    
    return gaussfunc;
}

/* 
 * void operaLMFitMultipleGaussian(unsigned m_dat, double *x, double *y, unsigned ngauss, double a[], double x0[], double sig[], double *chi2)
 * \brief This function performs Gaussian Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param a is a double pointer for the amplitude of the Gaussian
 * \param x0 is a double pointer for the center of the Gaussian
 * \param sig is a double pointer for the spread of the Gaussian
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
void operaLMFitMultipleGaussian(unsigned m_dat, double *x, double *y, unsigned ngauss, double a[], double x0[], double sig[], double *chi2)
{
	int n_par = 3;
    int totalnpar = n_par*(int)ngauss;
	double *par;
	par = (double*) malloc(totalnpar * sizeof(double));		
	
    for(unsigned i=0;i<ngauss;i++) {
        par[0+i*(unsigned)n_par] = a[i];
        par[1+i*(unsigned)n_par] = x0[i];
        par[2+i*(unsigned)n_par] = sig[i];
//        printf("%lf\t%lf\t%lf\t",par[0+i*(unsigned)n_par],par[1+i*(unsigned)n_par],par[2+i*(unsigned)n_par]);
        
	}
 	
	/* auxiliary parameters */
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
	
	/* perform the fit */
	lmcurve_fit(totalnpar, par, m_dat, x, y, GaussianFunction, &control, &status);
#ifdef PRINT_DEBUG	
	printf( "\n *** Results:\n" );
	printf( "Status after %d function evaluations:\n  %s\n",
           status.nfev, lm_infmsg[status.info] );
#endif	
	double WSSR,fcalc;
	WSSR = 0;
	
	for(unsigned i=0;i<m_dat;i++)
	{
		fcalc =  GaussianFunction(x[i],par,totalnpar);
		WSSR = WSSR + pow((fcalc - y[i]),2.0);
	}   
	
	*chi2 = WSSR/(double)(m_dat - totalnpar);
#ifdef PRINT_DEBUG	
	double rms = sqrt(WSSR/(double)(m_dat - totalnpar));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-totalnpar));    
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-totalnpar));	 
	printf("\n *** Obtained parameters:\n");
	for (unsigned i = 0; i < totalnpar; ++i) {
		printf("  par[%i] = %12g\n", i, par[i]);
    }
	printf("\n *** Obtained norm:\n  %12g\n", status.fnorm );
	printf("\n *** Fit gaussian: \n");    
    printf("f(x) =	");
    for(unsigned i=0;i<ngauss;i++) {
        printf("+ %lf*exp(-(x-%lf)*(x-%lf)/(2*%lf*%lf))",par[0+i*(unsigned)n_par],par[1+i*(unsigned)n_par],par[1+i*(unsigned)n_par],par[2+i*(unsigned)n_par],par[2+i*(unsigned)n_par]); 	
    }
    printf("\n");           
#endif
    for(unsigned i=0;i<ngauss;i++) {           
        a[i] = par[0+i*(unsigned)n_par];
        x0[i] = par[1+i*(unsigned)n_par];
        sig[i] = par[2+i*(unsigned)n_par];	
//        printf("%lf\t%lf\t%lf\t",par[0+i*(unsigned)n_par],par[1+i*(unsigned)n_par],par[2+i*(unsigned)n_par]);        
    }    
#ifdef PRINT_DEBUG	
	printf("\n *** Data:\n");	
	printf("x\ty\tf(x)\ty-f(x)\n");
	for (unsigned i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i], GaussianFunction(x[i],par,totalnpar), y[i] - GaussianFunction(x[i],par,totalnpar) );
#endif
}

/* 
 * void operaLMFitGaussian(unsigned m_dat, double *x, double *y, double *a, double *x0, double *sig, double *chi2)
 * \brief This function performs Gaussian Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param a is a double pointer for the amplitude of the Gaussian
 * \param x0 is a double pointer for the center of the Gaussian
 * \param sig is a double pointer for the spread of the Gaussian
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
void operaLMFitGaussian(unsigned m_dat, double *x, double *y, double *a, double *x0, double *sig, double *chi2) 
{
	int i;  
	int n_par = 3;
	double *par;
	par = (double*) malloc(n_par * sizeof(double));	
	
	par[0] = *a;
	par[1] = *x0;
	par[2] = *sig;
	
	/* auxiliary parameters */
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
	
	/* perform the fit */
	lmcurve_fit( n_par, par, m_dat, x, y, GaussianFunction, &control, &status);
#ifdef PRINT_DEBUG	
	printf( "\n *** Results:\n" );
	printf( "Status after %d function evaluations:\n  %s\n",
			status.nfev, lm_infmsg[status.info] );
#endif	
	double WSSR,fcalc;
	WSSR = 0;
	
	for(i=0;i<m_dat;i++)
	{
		fcalc =  GaussianFunction(x[i],par,n_par);
		WSSR = WSSR + pow((fcalc - y[i]),2.0);
	}   
	
	*chi2 = WSSR/(double)(m_dat-n_par);
#ifdef PRINT_DEBUG	
	double rms = sqrt(WSSR/(double)(m_dat - n_par));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-n_par));    
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-n_par));	 
	printf("\n *** Obtained parameters:\n");
	for ( i = 0; i < n_par; ++i)
		printf("  par[%i] = %12g\n", i, par[i]);
	printf("\n *** Obtained norm:\n  %12g\n", status.fnorm );
	printf("\n *** Fit gaussian: \n");    	
	printf("f(x) =%lf*exp(-(x-%lf)*(x-%lf)/(2*%lf*%lf))\n",par[0],par[1],par[1],par[2],par[2]); 	
#endif
	*a = par[0];
	*x0 = par[1];
	*sig = par[2];	
#ifdef PRINT_DEBUG	
	printf("\n *** Data:\n");	
	printf("x\ty\tf(x)\ty-f(x)\n");
	for (i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i], GaussianFunction(x[i],par,n_par), y[i] - GaussianFunction(x[i],par,n_par) );
#endif
}

/* 
 * double Polynomial2DFunction(double x, double y, const double *p, int n_par)
 * \brief This function returns the value of a given 2D polynomial surface function.
 * \param x is a double input value for which the given polynomial is evaluated  
 * \param y is a double input value for which the given polynomial is evaluated
 * \param p is a const double pointer that contains the coefficients of the given 2D polynomial
 * \param n_par is an int that defines the order of the 2D polynomial
 * \return double value of the 2D polynomial 
 */
double Polynomial2DFunction(double x, double y, const double *p, int n_par)
{
	double fxy;
	int i,j,k;
	int m,n;
	
	n = n_par/2;
	m = n_par - n;
	
	k=0;
	fxy=0;
	for(j=0;j<m;j++){
		for(i=0;i<n;i++){	
			fxy += p[k++]*pow(x,(double)i)*pow(y,(double)j);
		}
	}
	return fxy;
}

/* \brief poly2Ddata_struct - a data structure to transmit model data to function evaluation */
typedef struct {
	double *x;
	double *y;
	double *fxy;
	double (*Polynomial2DFunction)(double x, double y, const double *par, int n_par);
} poly2Ddata_struct;

/* 
 * void evaluate_surface(int n_par, const double *par, int m_dat, const void *data, double *fvec, int *info )
 * \brief This function performs 2D polynomial function evaluation, determination of residues.
 * \param n_par is an int that ..,
 * \param par is a double pointer that ..,
 * \param m_dat is an int that ..,
 * \param data is a double pointer that ..,
 * \param fvex is a double pointer that ..,
 * \param info is a double pointer that ..,
 * \return void
 */
/* function evaluation, determination of residues */
void evaluate_surface(int n_par, const double *par, int m_dat, const void *data, double *fvec, int *info )
{
	/* for readability, perform explicit type conversion */
	poly2Ddata_struct *mydata;
	mydata = (poly2Ddata_struct*)data;
	
	int i;
	for ( i = 0; i < m_dat; i++ )
		fvec[i] = mydata->fxy[i] - mydata->Polynomial2DFunction(mydata->x[i], mydata->y[i], par, n_par);
}

/* 
 * void operaLMFit2DPolynomial(unsigned m_dat, double *x, double *y, double *fxy, int n_par, double *par, double *chi2)
 * \brief This function performs 2D Polynomial Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param fxy is a double pointer with f(x,y) data values 
 * \param n_par is an int for the order of the 2D polynomial to be fit
 * \param par is a double pointer for the coefficients of the 2D polynomial
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
void operaLMFit2DPolynomial(unsigned m_dat, double *x, double *y, double *fxy, int n_par, double *par, double *chi2) 
{
	int i;
	
	poly2Ddata_struct data = {x, y, fxy, Polynomial2DFunction};
	
	/* auxiliary parameters */
	
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
	
	/* perform the fit */
#ifdef PRINT_DEBUG	
	lmmin( n_par, par, m_dat, (const void*) &data,
		  evaluate_surface, &control, &status, lm_printout_std );
	printf( "\n *** Results:\n" );
	printf( "Status after %d function evaluations:\n  %s\n",
		   status.nfev, lm_infmsg[status.info] );
#else
	lmmin( n_par, par, m_dat, (const void*) &data,
		  evaluate_surface, &control, &status, NULL);
#endif	
	double WSSR,fcalc;
	WSSR = 0;
	
	for(i=0;i<m_dat;i++)
	{
		fcalc =  Polynomial2DFunction(x[i],y[i],par,n_par);
		WSSR += pow((fcalc - fxy[i]),2.0);
	}   
	
	*chi2 = WSSR/(double)(m_dat - n_par);
#ifdef PRINT_DEBUG	
	double rms = sqrt(WSSR/(double)(m_dat - n_par));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-n_par));    
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-n_par));	 
	printf("\n *** Data:\n");	
	printf("x\ty\tz\tf(x,y)\tz-f(x,y)\n");
	for (i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i],fxy[i],Polynomial2DFunction(x[i],y[i],par,n_par), fxy[i] - Polynomial2DFunction(x[i],y[i],par,n_par) );
#endif	
}

/* 
 * double Gaussian2DFunction(double x, double y, const double *p, int n_par)
 * \brief This function returns the value for a 2D Gaussian: f(x,y) = A*exp(-(((x-x0)/sqrt(2)*sigx)^2 + ((y-y0)/sqrt(2)*sigy)^2))
 * \param x is a double input x-value for which the given 2D Gaussian is evaluated  
 * \param y is a double input y-value for which the given 2D Gaussian is evaluated
 * \param p is a const double pointer for the coefficients of the 2D Gaussian
 * \param n_par is an int for the number of coefficients (must be 5!)
 * \return double value of 2D Gaussian
 */
double Gaussian2DFunction(double x, double y, const double *p, int n_par)
{
	double fxy;
	fxy = p[0]*exp(-((x-p[1])*(x-p[1])/(2*p[3]*p[3])+(y-p[2])*(y-p[2])/(2*p[4]*p[4])));
	
	return fxy;
}

/* \brief gauss2Ddata_struct - a data structure to transmit model data to function evalution */
typedef struct {
	double *x;
	double *y;	
	double *fxy;
	double (*Gaussian2DFunction)(double x, double y, const double *p, int n_par);
} gauss2Ddata_struct;

/* 
 * void evaluate_2DGaussian(int n_par, const double *par, int m_dat, const void *data, double *fvec, int *info )
 * \brief This function performs 2D Gaussian function evaluation, determination of residues.
 * \param n_par is an unsigned that ..,
 * \param par is a double pointer that ..,
 * \param m_dat is an int that ..,
 * \param data is a double pointer that ..,
 * \param fvec is a double pointer that ..,
 * \param info is an int pointer that ..,
 * \return void
 */
void evaluate_2DGaussian(int n_par, const double *par, int m_dat, const void *data, double *fvec, int *info )
{
	/* for readability, perform explicit type conversion */
	gauss2Ddata_struct *mydata;
	mydata = (gauss2Ddata_struct*)data;
	
	for (int i = 0; i < m_dat; i++ )
		fvec[i] = mydata->fxy[i] - mydata->Gaussian2DFunction(mydata->x[i], mydata->y[i], par, n_par);
}

/* 
 * void operaLMFit2DGaussian(unsigned m_dat, double *x, double *y, double *fxy, double *a, double *x0, double *y0, double *sigx, double *sigy, double *chi2)
 * \brief This function performs 2D Gaussian Fitting using the Levenberg-Marquardt method.
 * \param m_dat is an unsigned for the number of data points
 * \param x is a double pointer with x data values
 * \param y is a double pointer with y data values
 * \param a is a double pointer for the amplitude of the Gaussian
 * \param x0 is a double pointer for the x-center of the Gaussian
 * \param y0 is a double pointer for the y-center of the Gaussian 
 * \param sigx is a double pointer for the spread of the Gaussian in the x-direction
 * \param sigy is a double pointer for the spread of the Gaussian in the y-direction
 * \param chi2 is a double pointer that returns the reduced chi-square of the fit
 * \return void
 */
/*** Fit 2D Gaussian ***/
void operaLMFit2DGaussian(unsigned m_dat, double *x, double *y, double *fxy, double *a, double *x0, double *y0, double *sigx, double *sigy, double *chi2) 
{
	int n_par = 5;
	double *par;
	par = (double*) malloc(n_par * sizeof(double));
	
	par[0] = *a;
	par[1] = *x0;
	par[2] = *y0;
	par[3] = *sigx;
	par[4] = *sigy;
	
	gauss2Ddata_struct data = { x, y, fxy, Gaussian2DFunction};
	
	/* auxiliary parameters */
	
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	control.printflags = 0; // monitor status (+1) and parameters (+2)
#ifdef PRINT_DEBUG
	/* perform the fit */
	lmmin( n_par, par, m_dat, (const void*) &data,
		  evaluate_2DGaussian, &control, &status, lm_printout_std );
#else
	/* perform the fit */
	lmmin( n_par, par, m_dat, (const void*) &data,
		  evaluate_2DGaussian, &control, &status, NULL );
#endif	
	int i;
	double WSSR,fcalc;
	WSSR = 0;
	
	for(i=0;i<m_dat;i++)
	{
		fcalc =  Gaussian2DFunction(x[i],y[i],par,n_par);
		WSSR = WSSR + pow((fcalc - fxy[i]),2.0);
	}   
	
	*chi2 = WSSR/(double)(m_dat-n_par);
	*a = par[0];
	*x0 = par[1];
	*y0 = par[2];
	*sigx = par[3];
	*sigy = par[4];	
#ifdef PRINT_DEBUG	
	double rms = sqrt(WSSR/(double)(m_dat - n_par));
	printf("\n *** Statistics:\n");
	printf("Chi-Square = %lf, DOF = %d\n",WSSR,(m_dat-n_par));    
	printf("RMS of residuals (stdfit) = sqrt(chi2/ndf): %lf\n",rms);
	printf("variance of residuals (reduced chi-square) = chi2/ndf: %lf\n",WSSR/(double)(m_dat-n_par));	 
	printf("\n *** Data:\n");	
	printf("x\ty\tz\tf(x,y)\tz-f(x,y)\n");
	for (i = 0; i < m_dat; ++i)
		printf( "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",x[i], y[i],fxy[i],Gaussian2DFunction(x[i],y[i],par,n_par), fxy[i] - Gaussian2DFunction(x[i],y[i],par,n_par) );
}
#endif
}

/*** Spline interpolation ***/

/* 
 * void operaFitSpline(unsigned nin, float *xin,float *yin, unsigned nout, float *xout, float *yout)
 * \brief This function performs interpolation of an input array using a sufficiently smooth piecewise-polynomial function (cubic spline)
 * \param nin is an unsigned for the number of input data points
 * \param xin is a float array for the input x data points 
 * \param yin is a float array for the input y data points  
 * \param nout is an unsigned for the number of output data points
 * \param xout is a float array for the output x data points, which should be provided. 
 * \param yout is a float array for the output interpolated y data points, which is calculated
 * \return void
 */

void operaFitSpline(unsigned nin, const float *xin, const float *yin, unsigned nout, const float *xout, float *yout)
{
	float yp1 = (yin[1] - yin[0])/(xin[1] - xin[0]);
	float ypn = (yin[nin-1] - yin[nin-2])/(xin[nin-1] - xin[nin-2]);
	float *y2 = (float *)malloc(nin*sizeof(float));
	
	// Call cubicspline to get second derivatives
	cubicspline(xin, yin, nin, yp1, ypn, y2);
	
	// Call splineinterpolate for interpolations
	for (unsigned i=0; i<nout; i++) {
		splineinterpolate(xin, yin, y2, nin, xout[i] , &yout[i]);
	}
	free(y2);
}

void operaFitSplineDouble(unsigned nin, const double *xin, const double *yin, unsigned nout, const double *xout, double *yout)
{
	double yp1 = (yin[1] - yin[0])/(xin[1] - xin[0]);
	double ypn = (yin[nin-1] - yin[nin-2])/(xin[nin-1] - xin[nin-2]);
	double *y2 = (double *)malloc(nin*sizeof(double));
	
	// Call cubicspline to get second derivatives
	cubicsplineDouble(xin, yin, nin, yp1, ypn, y2);
	
	// Call splineinterpolate for interpolations
	for (unsigned i=0; i<nout; i++) {
		splineinterpolateDouble(xin, yin, y2, nin, xout[i] , &yout[i]);
	}
	free(y2);
}

/* 
 * void operaFit2DSpline(unsigned nxin, float *xin, unsigned nyin, float *yin,double *fxyin, unsigned nxout, float *xout, unsigned nyout, float *yout, float *fxyout)
 * \brief This function performs interpolation of an input 2D-array using a sufficiently smooth piecewise-polynomial function (cubic spline)
 * \param nxin is an unsigned for the number of x input data points
 * \param xin is a float array for the input x data points 
 * \param nyin is an unsigned for the number of y input data points 
 * \param yin is a float array for the input y data points  
 * \param fxyin is a float array for the input fxy data points  
 * \param nxout is an unsigned for the number of x output data points.
 * \param xout is a float array for the output x data points, which should be provided.
 * \param nyout is an unsigned for the number of y output data points.
 * \param yout is a float array for the output y data points, which should be provided. 
 * \param fxyout is a float array for the output interpolated fxy data points, which is calculated
 * \return void
 */

void operaFit2DSpline(unsigned nxin, float *xin, unsigned nyin, float *yin, float *fxyin, unsigned nxout, float *xout, unsigned nyout, float *yout, float *fxyout)
{
	unsigned ii, ipix;
	float x, y; 
	
	
	CMatrix z  = newCMatrix(nxin, nyin); // CMatrix of data
	CMatrix z2 = newCMatrix(nxin, nyin); // CMatrix containing 2nd derivatives

	// get a pointer to the base of the data and copy in the vector
	float *zbase = z[0];
	for (unsigned j=0; j<nxin*nyin; j++) {
		*zbase++ = *fxyin++;
	}  
	
	// Call cubicspline2D to get second derivatives 
	cubicspline2D(xin, yin, z, nxin, nyin, z2);  
	
	// Call splineinterpolate2D for interpolations 
	// ??? fixme
	for (unsigned j=1; j<=nyout; j++) {
		ii = (j-1)*nxout;  
		y = yout[j-1];
		for (unsigned i=1; i<=nxout;i ++) {
			x = xout[i-1];     
			ipix = ii + i - 1;    
			splineinterpolate2D(xin, yin, z, z2, nxin, nyin, x, y, &fxyout[ipix]);
		}
	}  
	
	deleteCMatrix(z);
	deleteCMatrix(z2);
}

/*
 cubicspline constructs a cubic spline given a set of x and y values, through
 these values.
 */
int cubicspline(const float *x, const float *y, unsigned n, float yp1, float ypn, float *y2)
{
	float p, qn, sigma, un;
	
	float *tmp = (float *)malloc(n*sizeof(float));
	if (yp1 > 0.99e30) {
		y2[0] = 0.0;
		tmp[0] = 0.0;
	} else {
		y2[0] = -0.5;
		tmp[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (unsigned i=1; i<=n-2; i++) {
		sigma = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sigma*y2[i-1]+2.0;
		y2[i] = (sigma-1.0)/p;
		tmp[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		tmp[i] = (6.0*tmp[i]/(x[i+1]-x[i-1])-sigma*tmp[i-1])/p;
	}
	if (ypn > 0.99e30) {
		qn = 0.0;
		un = 0.0;
	} else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*tmp[n-2])/(qn*y2[n-2]+1.0);
	for (int k=n-2; k>=0; k--) {
		y2[k] = y2[k]*y2[k+1]+tmp[k];
	}
	free(tmp);
	return operaErrorCodeOK;
}

int cubicsplineDouble(const double *x, const double *y, unsigned n, double yp1, double ypn, double *y2)
{
	double p, qn, sigma, un;
	
	double *tmp = (double *)malloc(n*sizeof(double));
	if (yp1 > 0.99e30) {
		y2[0] = 0.0;
		tmp[0] = 0.0;
	} else {
		y2[0] = -0.5;
		tmp[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (unsigned i=1; i<=n-2; i++) {
		sigma = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sigma*y2[i-1]+2.0;
		y2[i] = (sigma-1.0)/p;
		tmp[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		tmp[i] = (6.0*tmp[i]/(x[i+1]-x[i-1])-sigma*tmp[i-1])/p;
	}
	if (ypn > 0.99e30) {
		qn = 0.0;
		un = 0.0;
	} else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*tmp[n-2])/(qn*y2[n-2]+1.0);
	for (int k=n-2; k>=0; k--) {
		y2[k] = y2[k]*y2[k+1]+tmp[k];
	}
	free(tmp);
	return operaErrorCodeOK;
}


/*
 splineinterpolate uses the cubic spline generated with spline to interpolate values
 in the XY table.
 */
int splineinterpolate(const float *xa, const float *ya, const float *y2a, int n, float x, float *y)
{
	int low = -1;
	int high = -1;
	int k;
	float h, b, a;
	
	if ( low < 0 ){
		low = 0;
		high = n-1;
	} else {
		if (x < xa[low])
			low = 0;
		if (x > xa[high])
			high = n-1;
	}
	while (high-low > 1) {
		k = (high+low) >> 1;
		if (xa[k] > x)
			high = k;
		else
			low = k;
	}
	h = xa[high]-xa[low];
	if (h == 0.0) {
		return operaSplineInterpolationAxesEqual;
	} else {
		a = (xa[high]-x)/h;
		b = (x-xa[low])/h;
		*y = a*ya[low]+b*ya[high]+((a*a*a-a)*y2a[low]+(b*b*b-b)*y2a[high]) * (h*h) / 6.0;
		return operaErrorCodeOK;
	}
}

int splineinterpolateDouble(const double *xa, const double *ya, const double *y2a, int n, double x, double *y)
{
	int low = -1;
	int high = -1;
	int k;
	double h, b, a;
	
	if ( low < 0 ){
		low = 0;
		high = n-1;
	} else {
		if (x < xa[low])
			low = 0;
		if (x > xa[high])
			high = n-1;
	}
	while (high-low > 1) {
		k = (high+low) >> 1;
		if (xa[k] > x)
			high = k;
		else
			low = k;
	}
	h = xa[high]-xa[low];
	if (h == 0.0) {
		return operaSplineInterpolationAxesEqual;
	} else {
		a = (xa[high]-x)/h;
		b = (x-xa[low])/h;
		*y = a*ya[low]+b*ya[high]+((a*a*a-a)*y2a[low]+(b*b*b-b)*y2a[high]) * (h*h) / 6.0;
		return operaErrorCodeOK;
	}
}

void splineinterpolate2D(float *x1a, float *x2a, CMatrix ya, CMatrix y2a, int nx, int ny, float x1, float x2, float *y)
{
	float *ytmp  = (float *)malloc(nx*sizeof(float));
	float *yytmp = (float *)malloc(nx*sizeof(float));
	
	for (unsigned j=1; j<=ny; j++) {
		splineinterpolate(x2a, ya[j], y2a[j], nx, x2, &yytmp[j]);
	}
	cubicspline(x1a, yytmp, nx, 1.0e30, 1.0e30, ytmp);
	splineinterpolate(x1a, yytmp, ytmp, nx, x1, y);
	
	free(yytmp);
	free(ytmp);
}

void cubicspline2D(float *x1a, float *x2a, CMatrix ya, int nx, int ny, CMatrix y2a)
{
	for (unsigned j=1; j<=ny; j++)
		cubicspline(x2a, ya[j], ny, 1.0e30, 1.0e30, y2a[j]);
}
