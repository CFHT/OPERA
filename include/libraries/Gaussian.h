#ifndef GAUSSIAN_H
#define GAUSSIAN_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: Gaussian
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
 
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

#include "libraries/operaLibCommon.h"

/*! 
 * \sa class Gaussian
 * \brief Encapsulation of the gaussian.
 * \file Gaussian.h
 * \ingroup libraries
 */

class Gaussian {
	
private:
	unsigned numberOfPeaks;
	unsigned maxnumberOfPeaks;
	
    double *amplitudeVector;
	double *amplitudeErrors;
    
	double *sigmaVector;
	double *sigmaErrors;
	
    double *centerVector;
	double *centerErrors;
	
	double baselineSlope;
	double baselineSlopeError;
	double baselineIntercept;
	double baselineInterceptError;    
    bool useBaseline;
    
    double gausschisqr;
	
public:
	/*
	 * Constructors / Destructors
	 */
	Gaussian();	
	Gaussian(unsigned NumberOfPeaks);	
	Gaussian(unsigned NumberOfPeaks, const double* AmplitudeVector, const double* SigmaVector, const double* CenterVector);	
	Gaussian(unsigned NumberOfPeaks, 
             const double* AmplitudeVector, const double* SigmaVector, const double* CenterVector,
             const double* AmplitudeErrors, const double* SigmaErrors, const double* CenterErrors);	
    
	~Gaussian();
	
	/*!
	 * double getAmplitude(unsigned index)
	 * \brief This function gets the amplitude value at index.
	 * \brief usage: double a3 = getAmplitude(3);
	 * \param index is a unsigned index of the gaussian amplitude to get
	 * \return void
	 */
	double getAmplitude(unsigned index) const;
    
	/*!
	 * void setAmplitude(double amplitude, unsigned index)
	 * \brief This function sets the amplitude value at index.
	 * \brief usage: setAmplitude((double)a, 3);
	 * \param amplitude is a double value of the amplitude for the gaussian at index
	 * \return void
	 */
	void setAmplitude(double amplitude, unsigned index);
    
	/*!
	 * double getSigma(unsigned index)
	 * \brief This function gets the sigma value at index.
	 * \brief usage: double sig3 = getSigma(3);
	 * \param index is a unsigned index of the gaussian sigma to get
	 * \return void
	 */
	double getSigma(unsigned index) const;
    
	/*!
	 * void setSigma(double sigma, unsigned index)
	 * \brief This function sets the sigma value at index.
	 * \brief usage: setSigma((double)sig, 3);
	 * \param sigma is a double value of the sigma for the gaussian at index
	 * \return void
	 */
	void setSigma(double sigma, unsigned index);
    
 	/*!
	 * double getCenter(unsigned index)
	 * \brief This function gets the x center value at index.
	 * \brief usage: double x03 = getCenter(3);
	 * \param index is a unsigned index of the gaussian center to get
	 * \return void
	 */
	double getCenter(unsigned index) const;
    
	/*!
	 * void setCenter(double center, unsigned index)
	 * \brief This function sets the center value at index.
	 * \brief usage: setCenter((double)x0, 3);
	 * \param center is a double value of the center for the gaussian at index
	 * \return void
	 */
	void setCenter(double center, unsigned index);	   

	/*!
	 * double EvaluateGaussian(double x)
	 * \brief This function returns the value of a given gaussian function at x.
	 * \brief usage: double value = EvaluateGaussian((double)x);
	 * \param x is a double input value for which the given gaussian is evaluated
	 * \return double value of evaluation of the gaussian
	 */
	double EvaluateGaussian(double x) const;
    
	/*!
	 * double* getAmplitudeVector(void);
	 * \brief This function returns the double *amplitudeVector.
	 * \brief usage: double *vec = getAmplitudeVector();
	 * \return double * - the vector of gaussian amplitude coefficients
	 */
	const double* getAmplitudeVector(void) const;
	double* getAmplitudeVector(void);
	/*!
	 * double* getAmplitudeErrorVector(void);
	 * \brief This function returns the double *amplitudeErrorVector.
	 * \brief usage: double *errs = getAmplitudeErrorVector();
	 * \return double * - the vector of gaussian amplitude coefficient errors
	 */
	const double* getAmplitudeErrorVector(void) const;
	double* getAmplitudeErrorVector(void);
    
	/*!
	 * double* getSigmaVector(void);
	 * \brief This function returns the double *sigmaVector.
	 * \brief usage: double *vec = getSigmaVector();
	 * \return double * - the vector of gaussian sigma coefficients
	 */
	const double* getSigmaVector(void) const;
	double* getSigmaVector(void);
	/*!
	 * double* getSigmaErrorVector(void);
	 * \brief This function returns the double *sigmaErrorVector.
	 * \brief usage: double *errs = getSigmaErrorVector();
	 * \return double * - the vector of gaussian sigma coefficient errors
	 */
	const double* getSigmaErrorVector(void) const;
	double* getSigmaErrorVector(void);
    
	/*!
	 * double* getCenterVector(void);
	 * \brief This function returns the double *centerVector.
	 * \brief usage: double *vec = getCenterVector();
	 * \return double * - the vector of gaussian center coefficients
	 */
	const double* getCenterVector(void) const;
	double* getCenterVector(void);
    
	/*!
	 * double* getCenterErrorVector(void);
	 * \brief This function returns the double *centerErrorVector.
	 * \brief usage: double *errs = getCenterErrorVector();
	 * \return double * - the vector of gaussian center coefficient errors
	 */
	const double* getCenterErrorVector(void) const;
	double* getCenterErrorVector(void);
    
	/*!
	 * unsigned getNumberOfPeaks();
	 * \brief This function returns the unsigned number of peaks (=number of Gaussians).
	 * \brief usage: unsigned npar = getNumberOfPeaks<float>();
	 * \brief usage: unsigned npar = getNumberOfPeaks();
	 * \return unsigned - the number of peaks
	 */
	unsigned getNumberOfPeaks() const;
	/*!
	 * void setNumberOfPeaks(unsigned Order);
	 * \brief This function sets the unsigned number of peaks (=number of Gaussians).
	 * \brief usage: setNumberOfPeaks(3);
	 * \return void
	 */
	void setNumberOfPeaks(unsigned NumberOfPeaks);		

	/*!
	 * void setGaussianChisqr(double Chisqr)
	 * \brief This function sets chisqr.
	 * \brief usage: setGaussianChisqr(0.98);
	 * \return void
	*/
	void setGaussianChisqr(double Chisqr);
	/*!
	 * double getGaussianChisqr(void)
	 * \brief This function gets chisqr.
	 * \brief usage: double c = getGaussianChisqr();
	 * \return double
	 */
	double getGaussianChisqr(void) const;
    
   /*!
     * void MPFitModeltoData(void)
     * \brief This function performs a least-squares MP fit of the model to the data.
     * \brief usage: MPFitModeltoData(np);
     * \return void
     */
    void MPFitModeltoData(unsigned nPeaks, unsigned NumberOfDataPoints, double *Xdata, double *Ydata, double *Yerrors);
    void MPFitModeltoData(unsigned nPeaks, double BaselineIntercept, double BaselineSlope, unsigned NumberOfDataPoints, double *Xdata, double *Ydata, double *Yerrors);
    
    void setBaselineSlope(double BaselineSlope);
    void setBaselineSlopeError(double BaselineSlopeError);
    void setBaselineIntercept(double BaselineIntercept);
    void setBaselineInterceptError(double BaselineInterceptError);
    void setUseBaseline(bool UseBaseline);
    double getBaselineSlope(void) const;
    double getBaselineSlopeError(void) const;
    double getBaselineIntercept(void) const;
    double getBaselineInterceptError(void) const;
    bool getUseBaseline(void) const;
    
};
#endif
