#ifndef LIBOPERAIMAGE_H
#define LIBOPERAIMAGE_H
/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaImage
 Version: 1.0
 Description: ThisC library implements low level image routines..
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
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
 * operaImage
 * \author Eder Martioli
 * \brief Image manipulation routines in C.
 * \file operaImage.h
 * \ingroup libraries
 */
#ifdef __cplusplus
extern "C" {
#endif
	
#include <math.h>

/*
 * The static inline functions below operates img by a constant  (they overrride original img values)
 */
static inline void operaSumImbyConstant(long npixels, float *img, float number){while (npixels--) *img++ += number;}
static inline void operaSubtractImbyConstant(long npixels, float *img, float number){while (npixels--) *img++ -= number;}
static inline void operaMultiplyImbyConstant(long npixels, float *img, float number){while (npixels--) *img++ *= number;}
static inline void operaDivideImbyConstant(long npixels, float *img, float number){while (npixels--) *img++ /= number;}
static inline void operaImSubtractIm(long npixels, float *img1, float *img2) {while (npixels--) *img1++ -= *img2++;}
float *medianCombineFloat(unsigned depth, long npixels, float *master, float *arrays[]);
unsigned short *operaArrayMedianCombineUSHORT(unsigned depth, long npixels, unsigned short *master, unsigned short *arrays[]);
void operaImMean(unsigned depth, long npixels, float *master, float *arrays[]);	
void operaImWeightedMean(unsigned depth, long npixels, float *master, float *arrays[], float *weights[]);
void operaImSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master); 
void operaImWeightedSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *weights[], float *master);
void operaImAvgSigClip(unsigned depth, long npixels, float *master, float *arrays[], unsigned nsig);
float operaCCDVarDiff(unsigned depth, long npixels, float *arrays[], float *weight);
void operaImVarDiff(unsigned depth, long npixels, float *arrays[], float *diffvarimg);

/*! 
 * static inline void operaImMeanQuick(unsigned depth, long npixels, float *master, float *arrays[])
 * \brief Mean combine a series of arrays (images) into the master.
 * \param depth is an unsigned that ...
 * \param master is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \return void
 */
static inline void operaImMeanQuick(unsigned depth, long npixels, float *master, float *arrays[]) 
{
	long np;
	for (unsigned d=0; d<depth; d++) {
		float *p = master;		
		float *arr = arrays[d];
		np = npixels;
		
		while (np--)
			*p++ += *arr++;
	}	
	np = npixels;
	float *p = master;	
	while (np--)
		*p++ /= (float)depth;	
}	

/*! 
 * static inline void operaImWeightedMeanQuick(unsigned depth, long npixels, float *master, float *arrays[], float *weights[])
 * \brief Mean combine a series of arrays (images) weighted by weights into the master.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param master is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \param weights is a float pointer pointer that ...
 * \return void
 */
static inline void operaImWeightedMeanQuick(unsigned depth, long npixels, float *master, float *arrays[], float *weights[]) 
{
	float *weightsum = (float *)malloc(sizeof(float)*npixels);
	memset(weightsum, sizeof(float)*npixels, 0);
	
	for (unsigned d=0; d<depth; d++) {
		long np = npixels;
		
		float *p = master;
		float *arr = arrays[d];
		float *pwi = weights[d];
		float *ws = weightsum;		
		
		while (np--) {
			(*p++) = (*arr++) * (*pwi);
			(*ws++) += *pwi++;
		}
	}
	
	long np = npixels;
	float *p = master;
	float *ws = weightsum;	
	while (np--) {
		*p++ /= *ws++;
	}
	free(weightsum);
}	

/*! 
 * static inline void operaImSigQuick(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master)
 * \brief Calculate the standard deviation image from a series of images with respect to a master image.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param sigarray is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \param master is a float pointer that ...
 * \return void
 */
static inline void operaImSigQuick(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master) 
{
	long np;
	for (unsigned d=0; d<depth; d++) {
		float *p = master;		
		float *arr = arrays[d];
		float *sig = sigarray;
		np = npixels;
		
		while (np--){
			*sig++ += (*arr - *p)*(*arr - *p);
			arr++;
			p++;
		}	
	}	
	np = npixels;
	float *sig = sigarray;	
	while (np--) {
		*sig = sqrt(*sig/(float)depth);	
		sig++;
	}	
}	

/*! 
 * static inline void operaImWeightedSigQuick(unsigned depth, long npixels, float *sigarray, float *arrays[], float *weights[], float *master)
 * \brief Calculate the weighted standard deviation image from a series of images with respect to a master image.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param sigarray is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \param master is a float pointer that ...
 * \return void
 */
static inline void operaImWeightedSigQuick(unsigned depth, long npixels, float *sigarray, float *arrays[], float *weights[], float *master) 
{
	float *weightsum = (float *)malloc(sizeof(float)*npixels);
	memset(weightsum, sizeof(float)*npixels, 0);
	
	long np;
	for (unsigned d=0; d<depth; d++) {
		float *p = master;		
		float *arr = arrays[d];
		float *sig = sigarray;
		float *pwi = weights[d];
		float *ws = weightsum;	
		
		np = npixels;
		
		while (np--){
			*sig++ += (*pwi)*(*arr - *p)*(*arr - *p);
			*ws++ += *pwi++;			
			arr++;
			p++;
		}	
	}	
	
	np = npixels;
	float *sig = sigarray;	
	float *ws = weightsum;	
	
	while (np--) {
		*sig = sqrt( *sig / (*ws++) );	
		sig++;
	}	
	free(weightsum);	
}	

/*! 
 * static inline void operaImAvgSigClipQuick(unsigned depth, long npixels, float *master, float *arrays[], unsigned nsig)
 * \brief Average sigma clip combine a series of arrays (images) into the master.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param sigarray is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \param master is a float pointer that ...
 * \return void
 */
static inline void operaImAvgSigClipQuick(unsigned depth, long npixels, float *master, float *arrays[], unsigned nsig) 
{
	long np;
	unsigned d;
	
	float *avg = (float *)malloc(sizeof(float)*npixels);
	memset(avg, sizeof(float)*npixels, 0);	
	
	float *sigarray = (float *)malloc(sizeof(float)*npixels);
	memset(sigarray, sizeof(float)*npixels, 0);	
	
	for (d=0; d<depth; d++) {
		float *p = avg;		
		float *arr = arrays[d];
		np = npixels;
		
		while (np--)
			*p++ += *arr++;
	}	
	
	np = npixels;
	float *p = avg;	
	while (np--) {
		*p++ /= (float)depth;	
	}	
	
	for (d=0; d<depth; d++) {
		float *p = avg;		
		float *arr = arrays[d];
		float *sig = sigarray;
		np = npixels;
		
		while (np--){
			*sig++ += (*arr - *p)*(*arr - *p);
			arr++;
			p++;
		}	
	}	
	np = npixels;
	float *sig = sigarray;	
	while (np--) {
		*sig = sqrt(*sig/(float)depth);	
		sig++;
	}		
	
	float *weightsum = (float *)malloc(sizeof(float)*npixels);
	memset(weightsum, sizeof(float)*npixels, 0);	
	
	for (d=0; d<depth; d++) {
		float *arr = arrays[d];
		float *sig = sigarray;
		float *p = avg;
		float *upp = master;	
		float *w = weightsum;		
		np = npixels;		
		
		while (np--) {	
			if(*arr > (*p - *sig*(float)nsig) && *arr < (*p + *sig*(float)nsig)) {
				*upp += *arr;
				*w += 1;
			}
			arr++;
			p++;
			sig++;
			upp++;
			w++;
		}		
	}
	
	np = npixels;
	float *upp = master;	
	float *w = weightsum;		
	float *pold = avg;
	
	while (np--) {
		if(*upp == 0) {
			*upp++ = *pold; // if all points are clipped then take regular mean
		}	else {
			*upp++ /= *w;	
		}
		pold++;
		w++;
	}		
	
	free(weightsum);
	free(avg);
	free(sigarray);
}
	
#ifdef __cplusplus
}
#endif
		
#endif

