
#ifndef OPERAFFT_H
#define OPERAFFT_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: liboperaFFT
 a wrapper for the open source library FFTW3
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

/*! 
 * operaFFT
 * \author Doug Teeple / Eder Martioli
 * \brief Common C language Fast Fourier Transform functions.
 * \file operaFFT.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
	/* prototypes */
	
	void operaFFTForward(unsigned np,double *y_Re, double *y_Im,double *fft_y_Re,double *fft_y_Im);
	void operaFFTBackward(unsigned np,double *y_Re, double *y_Im,double *fft_y_Re,double *fft_y_Im);
	void operaFFTLowPass(unsigned np,float *xin, float *xout, float cutfreq);
	void operaFFTHighPass(unsigned np,float *xin, float *xout, float cutfreq);
	void operaFFTBandPass(unsigned np,float *xin, float *xout, float lowfreq, float highfreq);
	void operaFFTLowPassDouble(unsigned np, double *xin, double *xout, float cutfreq);
	void operaFFTHighPassDouble(unsigned np, double *xin, double *xout, float cutfreq);
	void operaFFTBandPassDouble(unsigned np, double *xin, double *xout, float lowfreq, float highfreq);
	void operaFFTPowSpc(unsigned np,float *xin,float *yin, float *freq, float *fftpow);
	void operaFFTPowSpcDouble(unsigned np, double *xin, double *yin, double *freq, double *fftpow);
	
	void operaZeroPadSymmetricFunc(unsigned npin, float *xin, float *yin, unsigned npout, float *xout, float *yout);
	float operaConvolve(unsigned np1,float *x1,float *y1,unsigned np2,float *x2,float *y2,float x);
	float operaXCorrelation(unsigned np1,float *x1,float *y1,unsigned np2,float *x2,float *y2,float x);
	double operaXCorrelationDouble(unsigned np1,double *x1,double *y1,unsigned np2,double *x2,double *y2,double x);
	
	/*void operaFFTConvolution
	 void operaFFTCrossCorrel
	 void operaFFTAutoCorrel
	 void operaFFTWienerFilter
	 void operaFFTWienerDeconvol
	 */
#ifdef __cplusplus
}
#endif


#endif
