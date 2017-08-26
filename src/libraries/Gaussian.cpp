/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: Polynomial
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

/*!
 * Gaussian
 * \author Eder Martioli
 * \brief This class encapsulates the Gaussian.
 * \file Gaussian.cpp
 * \ingroup libraries
 */

#include <math.h>

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"
#include "libraries/Gaussian.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"
#include "libraries/mpfit.h"
#include "libraries/operaLib.h"     // for itos

/*
 * Constructors / Destructors
 */

Gaussian::Gaussian() :
numberOfPeaks(0),
maxnumberOfPeaks(0),
amplitudeVector(NULL),
amplitudeErrors(NULL),
sigmaVector(NULL),
sigmaErrors(NULL),
centerVector(NULL),
centerErrors(NULL),
gausschisqr(0.0)
{
	throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
}

Gaussian::Gaussian(unsigned NumberOfPeaks) :
numberOfPeaks(0),
maxnumberOfPeaks(0),
amplitudeVector(NULL),
amplitudeErrors(NULL),
sigmaVector(NULL),
sigmaErrors(NULL),
centerVector(NULL),
centerErrors(NULL),
gausschisqr(0.0)
{
	if (NumberOfPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    numberOfPeaks = maxnumberOfPeaks = NumberOfPeaks;
    
    amplitudeVector = new double[numberOfPeaks];    
    amplitudeErrors = new double[numberOfPeaks];                
    
    centerVector = new double[numberOfPeaks]; 
    centerErrors = new double[numberOfPeaks];  
    
    sigmaVector = new double[numberOfPeaks]; 
    sigmaErrors = new double[numberOfPeaks];      

      
}

Gaussian::Gaussian(unsigned NumberOfPeaks, const double* AmplitudeVector, const double* SigmaVector, const double* CenterVector) :
numberOfPeaks(0),
maxnumberOfPeaks(0),
amplitudeVector(NULL),
amplitudeErrors(NULL),
sigmaVector(NULL),
sigmaErrors(NULL),
centerVector(NULL),
centerErrors(NULL),
gausschisqr(0.0)
{
	if (NumberOfPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    numberOfPeaks = maxnumberOfPeaks = NumberOfPeaks;
    amplitudeVector = new double[numberOfPeaks];             
    centerVector = new double[numberOfPeaks]; 
    sigmaVector = new double[numberOfPeaks]; 
    
    for(unsigned i=0;i<numberOfPeaks;i++) {
        amplitudeVector[i]=AmplitudeVector[i];
        centerVector[i]=CenterVector[i]; 
        sigmaVector[i]=SigmaVector[i];  
    }
    
    amplitudeErrors = new double[numberOfPeaks];             
    centerErrors = new double[numberOfPeaks]; 
    sigmaErrors = new double[numberOfPeaks];     
}

Gaussian::Gaussian(unsigned NumberOfPeaks, 
                   const double* AmplitudeVector, const double* SigmaVector, const double* CenterVector,
                   const double* AmplitudeErrors, const double* SigmaErrors, const double* CenterErrors) :
numberOfPeaks(0),
maxnumberOfPeaks(0),
amplitudeVector(NULL),
amplitudeErrors(NULL),
sigmaVector(NULL),
sigmaErrors(NULL),
centerVector(NULL),
centerErrors(NULL),
gausschisqr(0.0)
{
	if (NumberOfPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    numberOfPeaks = maxnumberOfPeaks = NumberOfPeaks;
    
    amplitudeVector = new double[numberOfPeaks];             
    centerVector = new double[numberOfPeaks]; 
    sigmaVector = new double[numberOfPeaks]; 
    
    amplitudeErrors = new double[numberOfPeaks];             
    centerErrors = new double[numberOfPeaks]; 
    sigmaErrors = new double[numberOfPeaks]; 
    
    for(unsigned i=0;i<numberOfPeaks;i++) {
        amplitudeVector[i]=AmplitudeVector[i];
        centerVector[i]=CenterVector[i]; 
        sigmaVector[i]=SigmaVector[i];  
        
        amplitudeErrors[i]=AmplitudeErrors[i];
        centerErrors[i]=CenterErrors[i]; 
        sigmaErrors[i]=SigmaErrors[i];  
    }    
}

Gaussian::~Gaussian(){
	
	if (amplitudeVector) {
		delete[] amplitudeVector;
        amplitudeVector = NULL;
    }
	
	if (amplitudeErrors) {
		delete[] amplitudeErrors;
        amplitudeErrors = NULL;
    }
    
	if (centerVector) {
		delete[] centerVector;
        centerVector = NULL;
    }
    
	if (centerErrors) {
		delete[] centerErrors;
        centerErrors = NULL;
    }
    
	if (sigmaVector) {
		delete[] sigmaVector;
        sigmaVector = NULL;
    }
    
	if (sigmaErrors) {
		delete[] sigmaErrors;
        sigmaErrors = NULL;
    }
    numberOfPeaks = maxnumberOfPeaks = 0;
}

/*
 * double getAmplitude(unsigned index)
 * \brief This function gets the amplitude value at index.
 * \brief usage: double a3 = getAmplitude(3);
 * \param index is a unsigned index of the gaussian amplitude to get
 * \return void
 */
double Gaussian::getAmplitude(unsigned index) const{
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return amplitudeVector[index];
}

/*
 * void setAmplitude(double amplitude, unsigned index)
 * \brief This function sets the amplitude value at index.
 * \brief usage: setAmplitude((double)a, 3);
 * \param amplitude is a double value of the amplitude for the gaussian at index
 * \return void
 */
void Gaussian::setAmplitude(double amplitude, unsigned index){
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    amplitudeVector[index] = amplitude;
}

/*
 * double getSigma(unsigned index)
 * \brief This function gets the sigma value at index.
 * \brief usage: double sig3 = getSigma(3);
 * \param index is a unsigned index of the gaussian sigma to get
 * \return void
 */
double Gaussian::getSigma(unsigned index) const{
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return sigmaVector[index];
}

/*
 * void setSigma(double sigma, unsigned index)
 * \brief This function sets the sigma value at index.
 * \brief usage: setSigma((double)sig, 3);
 * \param sigma is a double value of the sigma for the gaussian at index
 * \return void
 */
void Gaussian::setSigma(double sigma, unsigned index){
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    sigmaVector[index] = sigma;
}

/*
 * double getCenter(unsigned index)
 * \brief This function gets the x center value at index.
 * \brief usage: double x03 = getCenter(3);
 * \param index is a unsigned index of the gaussian center to get
 * \return void
 */
double Gaussian::getCenter(unsigned index) const{
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return centerVector[index];
}

/*
 * void setCenter(double center, unsigned index)
 * \brief This function sets the center value at index.
 * \brief usage: setCenter((double)x0, 3);
 * \param center is a double value of the center for the gaussian at index
 * \return void
 */
void Gaussian::setCenter(double center, unsigned index){
	if (index > numberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	centerVector[index] = center; 
}

/*
 * double EvaluateGaussian(double x)
 * \brief This function returns the value of a given gaussian function at x.
 * \brief usage: double value = EvaluateGaussian((double)x);
 * \param x is a double input value for which the given gaussian is evaluated
 * \return double value of evaluation of the gaussian
 */
double Gaussian::EvaluateGaussian(double x) const{
	double gaussfunc = 0;
	if (numberOfPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned index=0;index<numberOfPeaks;index++){
        gaussfunc += amplitudeVector[index]*exp(-(x-centerVector[index])*(x-centerVector[index])/(2*sigmaVector[index]*sigmaVector[index]));
    }
    return gaussfunc;
}

/*
 * double* getAmplitudeVector(void);
 * \brief This function returns the double *amplitudeVector.
 * \brief usage: double *vec = getAmplitudeVector();
 * \return double * - the vector of gaussian amplitude coefficients
 */
const double* Gaussian::getAmplitudeVector(void) const {
    return amplitudeVector;
}
double* Gaussian::getAmplitudeVector(void){
    return amplitudeVector;
}

/*
 * double* getAmplitudeErrorVector(void);
 * \brief This function returns the double *amplitudeErrorVector.
 * \brief usage: double *errs = getAmplitudeErrorVector();
 * \return double * - the vector of gaussian amplitude coefficient errors
 */
const double* Gaussian::getAmplitudeErrorVector(void) const {
    return amplitudeErrors;	
}
double* Gaussian::getAmplitudeErrorVector(void){
    return amplitudeErrors;	
}

/*
 * double* getSigmaVector(void);
 * \brief This function returns the double *sigmaVector.
 * \brief usage: double *vec = getSigmaVector();
 * \return double * - the vector of gaussian sigma coefficients
 */
const double* Gaussian::getSigmaVector(void) const {
      return sigmaVector;  
}
double* Gaussian::getSigmaVector(void){
      return sigmaVector;  
}
/*
 * double* getSigmaErrorVector(void);
 * \brief This function returns the double *sigmaErrorVector.
 * \brief usage: double *errs = getSigmaErrorVector();
 * \return double * - the vector of gaussian sigma coefficient errors
 */
const double* Gaussian::getSigmaErrorVector(void) const {
	return sigmaErrors;
}
double* Gaussian::getSigmaErrorVector(void){
	return sigmaErrors;
}

/*
 * double* getCenterVector(void);
 * \brief This function returns the double *centerVector.
 * \brief usage: double *vec = getCenterVector();
 * \return double * - the vector of gaussian center coefficients
 */
const double* Gaussian::getCenterVector(void) const {
    return centerVector;
}
double* Gaussian::getCenterVector(void){
    return centerVector;
}

/*
 * double* getCenterErrorVector(void);
 * \brief This function returns the double *centerErrorVector.
 * \brief usage: double *errs = getCenterErrorVector();
 * \return double * - the vector of gaussian center coefficient errors
 */
const double* Gaussian::getCenterErrorVector(void) const {
    return centerErrors;
}
double* Gaussian::getCenterErrorVector(void){
    return centerErrors;
}

/*
 * unsigned getNumberOfPeaks();
 * \brief This function returns the unsigned number of peaks (=number of Gaussians).
 * \brief usage: unsigned npar = getNumberOfPeaks<float>();
 * \brief usage: unsigned npar = getNumberOfPeaks();
 * \return unsigned - the number of peaks
 */
unsigned Gaussian::getNumberOfPeaks() const{
	return numberOfPeaks;
}

/*
 * void setNumberOfPeaks(unsigned Order);
 * \brief This function sets the unsigned number of peaks (=number of Gaussians).
 * \brief usage: setNumberOfPeaks(3);
 * \return void
 */
void Gaussian::setNumberOfPeaks(unsigned NumberOfPeaks){
	if (NumberOfPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NumberOfPeaks > maxnumberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    numberOfPeaks = NumberOfPeaks;
}

/*
 * void setGaussianChisqr(double Chisqr)
 * \brief This function sets chisqr.
 * \brief usage: setGaussianChisqr(0.98);
 * \return void
 */

void Gaussian::setGaussianChisqr(double Chisqr){
    gausschisqr = Chisqr;
}

/*
 * double getGaussianChisqr(void)
 * \brief This function gets chisqr.
 * \brief usage: double c = getGaussianChisqr();
 * \return double
 */
double Gaussian::getGaussianChisqr(void) const{
	return gausschisqr;
}

void Gaussian::setBaselineSlope(double BaselineSlope) {
    baselineSlope = BaselineSlope;
}
void Gaussian::setBaselineSlopeError(double BaselineSlopeError) {
    baselineSlopeError = BaselineSlopeError;
}
void Gaussian::setBaselineIntercept(double BaselineIntercept) {
    baselineIntercept = BaselineIntercept;
}
void Gaussian::setBaselineInterceptError(double BaselineInterceptError) {
    baselineInterceptError = BaselineInterceptError;
}
void Gaussian::setUseBaseline(bool UseBaseline) {
    useBaseline = UseBaseline;
}
double Gaussian::getBaselineSlope(void) const {
    return baselineSlope;
}
double Gaussian::getBaselineSlopeError(void) const {
    return baselineSlopeError;
}
double Gaussian::getBaselineIntercept(void) const {
    return baselineIntercept;    
}
double Gaussian::getBaselineInterceptError(void) const {
    return baselineInterceptError;    
}
bool Gaussian::getUseBaseline(void) const {
    return useBaseline;   
}

/*
 * void MPFitModeltoData(unsigned nPeaks)
 * \brief This function performs a least-squares MP fit of the model to the data.
 * \brief usage: MPFitModeltoData(np);
 * \return void
 */
void Gaussian::MPFitModeltoData(unsigned nPeaks, unsigned NumberOfDataPoints, double *Xdata, double *Ydata, double *Yerrors) {  
	if (nPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (nPeaks > maxnumberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    setNumberOfPeaks(nPeaks);

    unsigned npar = 3*getNumberOfPeaks();
    
    if(npar >= NumberOfDataPoints) {
		throw operaException("Gaussian: invalid DOF (npars="+itos(npar)+") >= (nDataPoints="+itos(NumberOfDataPoints)+").",operaErrorCodeNOTIMPLEMENTED, __FILE__, __FUNCTION__, __LINE__);
    }
    
	double *par = (double*) malloc(npar * sizeof(double));
    double *epar = (double*) malloc(npar * sizeof(double));	
    
    mp_par pars[MAXGAUSSCOEFFS];    // create constraint structure
    memset(pars, 0, sizeof(pars));  // initialize constraint structure

    for(unsigned i=0;i<getNumberOfPeaks();i++) {
        par[0+i*3] = amplitudeVector[i];
        par[1+i*3] = centerVector[i];
        par[2+i*3] = sigmaVector[i];

        epar[0+i*3] = amplitudeErrors[i];
        epar[1+i*3] = centerErrors[i];
        epar[2+i*3] = sigmaErrors[i];
        
        // keep center fixed for the 1st fit
        pars[1+i*3].fixed = 1; 
        
        // force positive sigma: 
        pars[2+i*3].limited[0] = 1; 
        pars[2+i*3].limits[0] = 0.0;        
	}    

    struct vars_struct {
        double *x;
        double *y;
        double *ey;
    } v;    
    
    //int status;
    mp_result result;                   // create result structure
    memset(&result,0,sizeof(result));   // initialize result structure
    
    result.xerror = epar;
	
    v.x = Xdata;
    v.y = Ydata;
    v.ey = Yerrors;

    // Call fitting function (1st pass)
    /*status = */mpfit(MPGaussFunc, NumberOfDataPoints, npar, par, pars, 0, (void *) &v, &result);

    for(unsigned i=0;i<getNumberOfPeaks();i++) {    
        pars[1+i*3].fixed = 0; // free center for 2nd fit
	}

    // Call fitting function (2nd pass)
    /*status = */mpfit(MPGaussFunc, NumberOfDataPoints, npar, par, pars, 0, (void *) &v, &result);

    for(unsigned i=0;i<getNumberOfPeaks();i++) {            
        amplitudeVector[i] = par[0+i*3];
        centerVector[i] = par[1+i*3];
        sigmaVector[i] = par[2+i*3];
        
        amplitudeErrors[i] = result.xerror[0+i*3];
        centerErrors[i] = result.xerror[1+i*3];
        sigmaErrors[i] = result.xerror[2+i*3];
	}

	gausschisqr = result.bestnorm/(double)(NumberOfDataPoints-npar);
	
	free(par);
    free(epar);
}

/*
 * void MPFitModeltoData(unsigned nPeaks, double BaselineIntercept, double BaselineSlope)
 * \brief This function performs a least-squares MP fit of the model to the data.
 * \brief usage: MPFitModeltoData(np);
 * \return void
 */
void Gaussian::MPFitModeltoData(unsigned nPeaks, double BaselineIntercept, double BaselineSlope, unsigned NumberOfDataPoints, double *Xdata, double *Ydata, double *Yerrors) {  
	if (nPeaks == 0) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (nPeaks > maxnumberOfPeaks) {
		throw operaException("Gaussian: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    setNumberOfPeaks(nPeaks);
    setUseBaseline(true);
    
    mp_config config;
    memset(&config,0,sizeof(config));   // initialize config structure
    
    /* maxfev = Maximum number of function evaluations, or 0 for no limit
     Default: 0 (no limit) */
    config.maxfev = 200;
    
    /*maxiter = Maximum number of iterations.*/    
    /*ftol = Relative chi-square convergence criterium Default: 1e-10 */
    /*xtol = Relative parameter convergence criterium  Default: 1e-10 */
    /*gtol = Orthogonality convergence criterium       Default: 1e-10 */    
    
    unsigned npar = 3*getNumberOfPeaks() + 2;
    
    if(npar >= NumberOfDataPoints) {
		throw operaException("Gaussian: invalid DOF (npars="+itos(npar)+") >= (nDataPoints="+itos(NumberOfDataPoints)+").",operaErrorCodeNOTIMPLEMENTED, __FILE__, __FUNCTION__, __LINE__);
    }
    
	double *par = (double*) malloc(npar * sizeof(double));
    double *epar = (double*) malloc(npar * sizeof(double));	
    
    mp_par pars[MAXGAUSSCOEFFS];    // create constraint structure
    memset(pars, 0, sizeof(pars));  // initialize constraint structure
    
    for(unsigned i=0;i<getNumberOfPeaks();i++) {
        par[0+i*3] = amplitudeVector[i];
        par[1+i*3] = centerVector[i];
        par[2+i*3] = sigmaVector[i];
        
        epar[0+i*3] = amplitudeErrors[i];
        epar[1+i*3] = centerErrors[i];
        epar[2+i*3] = sigmaErrors[i];
        
        // keep centers, sigma and amplitudes fixed for the 1st fit
        pars[0+i*3].fixed = 1;        
        pars[1+i*3].fixed = 1;
        pars[2+i*3].fixed = 1;  
	}         
    par[npar-2] = BaselineIntercept;
    par[npar-1] = BaselineSlope;

    epar[npar-2] = BaselineIntercept*0.2;
    epar[npar-1] = BaselineSlope*0.2;    
    
    struct vars_struct {
        double *x;
        double *y;
        double *ey;
    } v;    
    
    //int status;
    mp_result result;                   // create result structure
    memset(&result,0,sizeof(result));   // initialize result structure
    
    result.xerror = epar;
	
    v.x = Xdata;
    v.y = Ydata;
    v.ey = Yerrors;
    
    // Call fitting function (1st pass)
    /*status = */mpfit(MPGaussFuncWithBaseline, NumberOfDataPoints, npar, par, pars, &config, (void *) &v, &result);    
    for(unsigned i=0;i<getNumberOfPeaks();i++) {   
        pars[0+i*3].fixed = 0; // free amplitudes for 2nd pass        
        pars[1+i*3].fixed = 0; // free centers for 2nd pass
        pars[2+i*3].fixed = 0; // free sigma for 2nd pass   
        
        // set sigma constraints: 
        pars[2+i*3].limited[0] = 1; 
        pars[2+i*3].limited[1] = 1;        
        pars[2+i*3].limits[0] = 0.0;   
        pars[2+i*3].limits[1] = par[2+i*3]*3;    
        
        // set amplitude constraints: 
        pars[0+i*3].limited[0] = 1;  
        pars[0+i*3].limits[0] = 0.0;   
	}
    // fix background for the 2nd pass
//    pars[npar-2].fixed = 1; 
//    pars[npar-1].fixed = 1;    
    
    // Call fitting function (2nd pass)
    /*status = */mpfit(MPGaussFuncWithBaseline, NumberOfDataPoints, npar, par, pars, &config, (void *) &v, &result);
 
    // free background for the 3rd pass
//    pars[npar-2].fixed = 0; 
//    pars[npar-1].fixed = 0;    
    
    // Call fitting function (3rd pass)
//    status = mpfit(MPGaussFuncWithBaseline, nDataPoints, npar, par, pars, &config, (void *) &v, &result);    
    
    for(unsigned i=0;i<getNumberOfPeaks();i++) {            
        amplitudeVector[i] = par[0+i*3];
        centerVector[i] = par[1+i*3];
        sigmaVector[i] = par[2+i*3];
        
        amplitudeErrors[i] = result.xerror[0+i*3];
        centerErrors[i] = result.xerror[1+i*3];
        sigmaErrors[i] = result.xerror[2+i*3];
	}
    baselineIntercept = par[npar-2];
    baselineSlope = par[npar-1];
    
    baselineInterceptError = result.xerror[npar-2];
    baselineSlopeError = result.xerror[npar-1];    
    
	gausschisqr = result.bestnorm/(double)(NumberOfDataPoints-npar);
	
	free(par);
    free(epar);
}
