#ifndef OPERACCD_H
#define OPERACCD_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: liboperaFITSSubImage
 Class: operaCCD
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
 * operaCCD 
 * \author Eder Martioli
 * \brief CCD manipulation routines in C.
 * \file operaCCD.h
 * \ingroup libraries
 */

#define MAXORDERS 200
#define MAXIPPOINTS 100

#define DETECTTHRESHOLD 0.4
#define NPINSAMPLE 20
#define NSAMPLES 3

#define NUMBEROFSAMPLESTORECENTERORDER 200
#define FRACTIONOFIPTOSCANFORRECENTERORDER 2
#define MINIMUMCROSSCORRELATIONTORECENTERORDER 0.1

//#ifdef __cplusplus
//extern "C" {
//#endif
	
	/* prototypes for gain calculations*/
    void operaMaskPixbyCountRange(unsigned npixels, float *array, float *badpixmask, float *newmask, float minvalue, float maxvalue, unsigned char invert);
    void operaMarkPixbyCountRange(unsigned npixels, float *array, float *previousmarks, float *newmarks, float minvalue, float maxvalue, float index2mark);
    void operaCCDGainNoise(unsigned npixels, unsigned nbias, float *biasdata[], unsigned nflat, float *flatdata[], float *badpixdata, float lowcount, float highcount, float maxbins, float minnpixelsperbin, float *gain, float *gainerror, float *bias, float *noise);

    /* prototypes for geometry/findorders calculations*/
    unsigned operaCCDDetectPeaksWithErrorsUsingIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr);
    unsigned operaCCDDetectPeaksWithErrorsUsingGaussian(unsigned np, float *x,float *y,float sigma, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr);
    unsigned operaCCDDetectPeaksWithErrorsUsingGaussianDouble(unsigned np, double *x,double *y,double sigma, double noise, double gain, double threshold,double *xmean, double *ymean, double *xmeanerr);
    unsigned operaCCDDetectPeaksWithErrorsUsingTopHat(unsigned np, float *x,float *y,unsigned width, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr);

    unsigned operaCCDDetectPeaksWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean);
    unsigned operaCCDDetectPeaksWithGaussian(unsigned np, float *x,float *y,float sigma, float noise, float gain, float threshold,float *xmean, float *ymean);
    unsigned operaCCDDetectPeaksWithTopHat(unsigned np, float *x,float *y,unsigned width, float noise, float gain, float threshold,float *xmean, float *ymean);

    int operaCCDDetectMissingOrders(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,int AbsPrevOrdNumber,float x0prev);
    
    int operaCCDRecenterOrderUsingXCorrWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float *ipx, float noise, float gain, float *xmean, float *ymean, float *xmeanerr);
    
    int operaCCDRecenterOrderWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float *ipx, float noise, float gain, float *xmean, float *ymean, float *xmeanerr);

    void operaCCDFitIP(unsigned np, float *x,float *y, unsigned nords, float *xmean, float *ymean, float *ipfunc, float *ipx,float *iperr, unsigned slit);

    void operaMedianWidthFromSetOfLines(unsigned np, float *mx, float *my, float *myerr, unsigned nlines, float *xlines, float *ylines, float *medianWidth);
    
    unsigned operaCCDDetectPeaksByXCorrWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean);
    
    int operaCCDDetectOrderMapBasedOnSpacingPolynomial(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,unsigned minordertouse, unsigned maxorders);

    unsigned matchMeasuredOrdersWithMap(unsigned np,float *fx,float *fy,float slit, unsigned npar,double *par,unsigned nords,float *xmean,float *ymean, float *orderMap, float *ordSepPred, float *ordXPosPred,float *ordYValue, unsigned *ordIndex, unsigned j, unsigned nAccHops,unsigned minordertouse, unsigned maxorders, float *xiSqr);

    float operaXCorrelationWithRefOrders(unsigned np,float *fx,float *fy,unsigned nrefs, float *xref,float *yref,unsigned nords, float *xmean,float *ymean, float slit, float xrange, float xstep);

    int operaCCDDetectMissingOrdersUsingRefMap(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,unsigned nrefs,float *xref,float *yref,int *AbsRefOrdNumber,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,float xrange,float xstep);

    int operaCCDDetectMissingOrdersUsingNearMap(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,unsigned nrefs,float *xref,float *yref,int *AbsRefOrdNumber,float *xord,float *yord, float *xerrord, int *AbsOrdNumber);


//#ifdef __cplusplus
//}
//#endif
		
		
#endif
