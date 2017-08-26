/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFFT
 Version: 1.0
 Description: FFT-related library routines.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Oct/2011
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


#include "fftw3.h"

#include "globaldefines.h"
#include "operaError.h"

/*
 * \brief FFT library.
 * \file operaFFT.c
 * \ingroup libraries
 */

#include "libraries/operaFFT.h"
#include "libraries/operaFit.h"
#include "libraries/operaStats.h"

/*** 1D FFT Forward ***/
void operaFFTForward(unsigned np,double *y_Re, double *y_Im,double *fft_y_Re,double *fft_y_Im)
{
	unsigned i;
	
	fftw_complex *in, *out;  
	fftw_plan p;
	
	in = (fftw_complex*) malloc(sizeof(fftw_complex) * np);
	out = (fftw_complex*) malloc(sizeof(fftw_complex) * np);
	
	p = fftw_plan_dft_1d(np,in,out,FFTW_FORWARD,FFTW_MEASURE);
	
	for(i=0;i<np;i++)
	{
		in[i][0] = y_Re[i];
		in[i][1] = y_Im[i];	
	}
	
	fftw_execute(p);
	
	for(i=0;i<np;i++)
	{
		fft_y_Re[i] = out[i][0];
		fft_y_Im[i] = out[i][1];  
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}

/*** 1D FFT Backward ***/
void operaFFTBackward(unsigned np,double *y_Re, double *y_Im,double *fft_y_Re,double *fft_y_Im)
{
	unsigned i;
	
	fftw_complex *in, *out;  
	fftw_plan p;
	
	in = (fftw_complex*) malloc(sizeof(fftw_complex) * np);
	out = (fftw_complex*) malloc(sizeof(fftw_complex) * np);
	
	p = fftw_plan_dft_1d(np,in,out,FFTW_BACKWARD,FFTW_MEASURE);
	
	for(i=0;i<np;i++)
	{
		in[i][0] = y_Re[i];
		in[i][1] = y_Im[i];	
	}
	
	fftw_execute(p);
	
	for(i=0;i<np;i++)
	{
		fft_y_Re[i] = out[i][0];
		fft_y_Im[i] = out[i][1];  
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}

/*** Low pass FFT filter ***/

void operaFFTLowPass(unsigned np,float *xin, float *xout, float cutfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq > cutfreq && freq < 1 - cutfreq) {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		} else {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;			
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

/*** High pass FFT filter ***/

void operaFFTHighPass(unsigned np,float *xin, float *xout, float cutfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq < cutfreq || freq > 1 - cutfreq) {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		} else {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;			
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

/*** Band pass FFT filter ***/

void operaFFTBandPass(unsigned np,float *xin, float *xout, float lowfreq, float highfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq > lowfreq && freq < highfreq) {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;
		} else if	(freq > 1-highfreq && freq < 1-lowfreq) {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;
		} else {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

void operaFFTLowPassDouble(unsigned np, double *xin, double *xout, float cutfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq > cutfreq && freq < 1 - cutfreq) {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		} else {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;			
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

/*** High pass FFT filter ***/

void operaFFTHighPassDouble(unsigned np, double *xin, double *xout, float cutfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq < cutfreq || freq > 1 - cutfreq) {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		} else {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;			
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

/*** Band pass FFT filter ***/

void operaFFTBandPassDouble(unsigned np, double *xin, double *xout, float lowfreq, float highfreq) {
	
	double *x_Re = (double *) malloc (np * sizeof(double));
	double *x_Im = (double *) malloc (np * sizeof(double));		
	double *fft_x_Re = (double *) malloc (np * sizeof(double));
	double *fft_x_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		x_Re[i] = (double)xin[i];
		x_Im[i] = 0;
	}
	
	operaFFTForward(np,x_Re, x_Im,fft_x_Re,fft_x_Im);
	
	float freq;
	
	for(unsigned i=0;i<np;i++) {
		freq = (float)i/(float)np;
		if(freq > lowfreq && freq < highfreq) {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;
		} else if	(freq > 1-highfreq && freq < 1-lowfreq) {
			fft_x_Re[i] /= (double)np;
			fft_x_Im[i] /= (double)np;
		} else {
			fft_x_Re[i] = 0;
			fft_x_Im[i] = 0;
		}
	}	
	
	double *xout_Re = (double *) malloc (np * sizeof(double));
	double *xout_Im = (double *) malloc (np * sizeof(double));	
	
	operaFFTBackward(np,fft_x_Re,fft_x_Im,xout_Re,xout_Im);
	for(unsigned i=0;i<np;i++) {
		xout[i] = xout_Re[i];		
#ifdef PRINT_DEBUG
		freq = (float)i/(float)np;
		printf("%d\t%f\t%f\t%lf\t%lf\t%lf\t%.15lf\n",i,freq,xin[i], fft_x_Re[i], fft_x_Im[i], xout_Re[i], xout_Im[i]);
#endif
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(x_Re);
	free(x_Im);		
	free(fft_x_Re);
	free(fft_x_Im);	
	free(xout_Re);
	free(xout_Im);	
}

/*** FFT Power Spectrum ***/

void operaFFTPowSpc(unsigned np,float *xin,float *yin, float *freq, float *fftpow) {
	
	double *y_Re = (double *) malloc (np * sizeof(double));
	double *y_Im = (double *) malloc (np * sizeof(double));		
	double *fft_y_Re = (double *) malloc (np * sizeof(double));
	double *fft_y_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		freq[i] = ((float)i/(float)np)*((xin[np-1] - xin[0])/(float)np);
		y_Re[i] = (double)yin[i];
		y_Im[i] = 0;
	}
	
	operaFFTForward(np,y_Re, y_Im,fft_y_Re,fft_y_Im);
	
	for(unsigned i=0;i<np;i++) {
		fft_y_Re[i] /= (double)np;
		fft_y_Im[i] /= (double)np;		
		fftpow[i] = (1./(float)np)*(float)(fft_y_Re[i]*fft_y_Re[i] + fft_y_Im[i]*fft_y_Im[i]);
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(y_Re);
	free(y_Im);		
	free(fft_y_Re);
	free(fft_y_Im);	
}

void operaFFTPowSpcDouble(unsigned np, double *xin, double *yin, double *freq, double *fftpow) {
	
	double *y_Re = (double *) malloc (np * sizeof(double));
	double *y_Im = (double *) malloc (np * sizeof(double));		
	double *fft_y_Re = (double *) malloc (np * sizeof(double));
	double *fft_y_Im = (double *) malloc (np * sizeof(double));	
	
	for(unsigned i=0;i<np;i++) {
		freq[i] = ((float)i/(float)np)*((xin[np-1] - xin[0])/(float)np);
		y_Re[i] = (double)yin[i];
		y_Im[i] = 0;
	}
	
	operaFFTForward(np,y_Re, y_Im,fft_y_Re,fft_y_Im);
	
	for(unsigned i=0;i<np;i++) {
		fft_y_Re[i] /= (double)np;
		fft_y_Im[i] /= (double)np;		
		fftpow[i] = (1./(float)np)*(float)(fft_y_Re[i]*fft_y_Re[i] + fft_y_Im[i]*fft_y_Im[i]);
	}
	// DT Sept 6 2013 -- fix memory leaks
	free(y_Re);
	free(y_Im);		
	free(fft_y_Re);
	free(fft_y_Im);	
}


/*** ZeroPadding ***/
void	operaZeroPadSymmetricFunc(unsigned npin, float *xin, float *yin, unsigned npout, float *xout, float *yout) 
{
#ifdef IMPLEMENTED
	int odd_in = 0;
	int odd_out = 0;	
	
#ifdef PRINT_DEBUG
	printf("%d\t%d\n",npin, npout);
#endif	
	// test parity
	if (npin - 2*ceil((float)npin/2))
		odd_in = 1;
	
	if (npout - 2*ceil((float)npout/2))
		odd_out = 1;
	
#ifdef PRINT_DEBUG
	printf("%d\t%d\n",odd_in, odd_out);
#endif
#endif
}

/*** Convolution between 2 functions (f*g)(x) where f= y1(x1) and g=y2(x2) ***/
float operaConvolve(unsigned np1, float *x1, float *y1, unsigned np2, float *x2, float *y2, float x)
{
	
	// is x in the range?
	/*	if(x < x1[0] || x > x1[np1-1])
	 exit(EXIT_FAILURE);
	 */
	unsigned i;
	
	// create x,y for f1 with the same sampling as f2		
	float *x1out = (float *) malloc( np2 * sizeof(float));
	float *y1out = (float *) malloc( np2 * sizeof(float));
	
	//	printf("%u\t%f\t%f\n",np1,x1[0],x1[np1-1]);
	
	unsigned np=0;
	
	for(i=0; i<np2; i++) {
		x1out[np2-i-1] = x - x2[i];		
		
		if(x1out[np2-i-1] >= x1[0] && x1out[np2-i-1] <= x1[np1-1])
			np++;
	}	
	
	float *xfit = (float *) malloc( np * sizeof(float));
	float *yfit = (float *) malloc( np * sizeof(float));	
	
	np=0;
	for(i=0; i<np2; i++) {
		//		 printf("%d\t%f",i,x1out[i]);	
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			xfit[np] = x1out[i];
			//			printf("\t%u",np);
			np++;	
		}	
		//		printf("\n");
	}	
	
	operaFitSpline(np1, x1, y1, np, xfit, yfit);		
	
	np=0;
	for(i=0; i<np2; i++) {
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			y1out[i] = yfit[np];
			np++;	
		}	else {
			y1out[i] = 0;
		}
	}	
	
	/*	for(i=0; i<np2; i++) {
	 printf("%f\t%f\t%f\t%f\n",x2[np2-i-1], y2[np2-i-1], x1out[i], y1out[i]);
	 }		
	 */
	float conv = 0;
	
	for(i=0; i<np2; i++) {
		conv += y1out[i]*y2[np2-i-1];
	}
	free(x1out);
	free(y1out);
	free(xfit);
	free(yfit);
	return conv;
}

/*** Cross correlation between 2 functions (f*g)(x) where f= y1(x1) and g=y2(x2) ***/
float operaXCorrelation(unsigned np1,float *x1,float *y1,unsigned np2,float *x2,float *y2,float x)
{
	// is x in the range?
	/*	if(x < x1[0] || x > x1[np1-1])
	 exit(EXIT_FAILURE);
	 */
	unsigned i;
	
	// create x,y for f1 with the same sampling as f2		
	float *x1out = (float *) malloc( np2 * sizeof(float));
	float *y1out = (float *) malloc( np2 * sizeof(float));
	
	//printf("%u\t%f\t%f\n",np1,x1[0],x1[np1-1]);
	
	unsigned np=0;
	
	for(i=0; i<np2; i++) {
		x1out[i] = x + x2[i];		
		
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1])
			np++;
	}	
	
	float *xfit = (float *) malloc( np2 * sizeof(float));
	float *yfit = (float *) malloc( np2 * sizeof(float));	
	
	np=0;
	for(i=0; i<np2; i++) {
		//		 printf("%d\t%f",i,x1out[i]);	
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			xfit[np] = x1out[i];
			//			printf("\t%u",np);
			np++;	
		}	
		//		printf("\n");
	}	
	
	operaFitSpline(np1, x1, y1, np, xfit, yfit);		
	
	np=0;
	for(i=0; i<np2; i++) {
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			y1out[i] = yfit[np];
			np++;	
		} else {
			y1out[i] = 0;
		}
	}	
	
	/*	for(i=0; i<np2; i++) {
	 printf("%f\t%f\t%f\t%f\n",x2[np2-i-1], y2[np2-i-1], x1out[i], y1out[i]);
	 }		
	 */
	float xcorr = 0;
	
	for(i=0; i<np2; i++) {
		xcorr += y1out[i]*y2[i];
	}
	free(x1out);
	free(y1out);
	free(xfit);
	free(yfit);
	return xcorr;
}

double operaXCorrelationDouble(unsigned np1,double *x1,double *y1,unsigned np2,double *x2,double *y2,double x)
{
	// is x in the range?
	/*	if(x < x1[0] || x > x1[np1-1])
	 exit(EXIT_FAILURE);
	 */
	unsigned i;
	
	// create x,y for f1 with the same sampling as f2		
	double *x1out = (double *) malloc( np2 * sizeof(double));
	double *y1out = (double *) malloc( np2 * sizeof(double));
	
	//printf("%u\t%f\t%f\n",np1,x1[0],x1[np1-1]);
	
	unsigned np=0;
	
	for(i=0; i<np2; i++) {
		x1out[i] = x + x2[i];		
		
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1])
			np++;
	}	
	
	double *xfit = (double *) malloc( np2 * sizeof(double));
	double *yfit = (double *) malloc( np2 * sizeof(double));	
	
	np=0;
	for(i=0; i<np2; i++) {
		//		 printf("%d\t%f",i,x1out[i]);	
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			xfit[np] = x1out[i];
			//			printf("\t%u",np);
			np++;	
		}	
		//		printf("\n");
	}	
	
	operaFitSplineDouble(np1, x1, y1, np, xfit, yfit);		
	
	np=0;
	for(i=0; i<np2; i++) {
		if(x1out[i] >= x1[0] && x1out[i] <= x1[np1-1]){		
			y1out[i] = yfit[np];
			np++;	
		} else {
			y1out[i] = 0;
		}
	}	
	
	/*	for(i=0; i<np2; i++) {
	 printf("%f\t%f\t%f\t%f\n",x2[np2-i-1], y2[np2-i-1], x1out[i], y1out[i]);
	 }		
	 */
	float xcorr = 0;
	
	for(i=0; i<np2; i++) {
		xcorr += y1out[i]*y2[i];
	}
	free(x1out);
	free(y1out);
	free(xfit);
	free(yfit);
	return xcorr;
}

