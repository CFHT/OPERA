#ifndef LIBOPERASTATS_H
#define LIBOPERASTATS_H
/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaStats
 Version: 1.0
 Description: ThisC library implements basic statistics routines..
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
 * operaStats 
 * \author Eder Martioli
 * \author Doug Teeple
 * \brief Image statistics routines in C.
 * \file operaStats.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
#include "operaError.h"			// for error codes
#include <stdlib.h>				// for random
#include <math.h>
#include <string.h>				// for memset

// for Linux...
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define ELEM_SWAP(a,b) {register float t=(a);(a)=(b);(b)=t;}
#define ELEM_SWAP_d(a,b) {register double t=(a);(a)=(b);(b)=t;}
#define ELEM_USHORTSWAP(a,b) {register unsigned short t=(a);(a)=(b);(b)=t;}
	
#define OneSigma 1.0
#define TwoSigma 2.0
#define ThreeSigma 3.0
#define FourSigma 4.0
#define FiveSigma 5.0

float operaArrayMean(unsigned np, const float *array);
double operaArrayMean_d(unsigned np, const double *array);
float operaArrayWeightedMean(unsigned np, const float *array, const float *weigh);
float operaArrayAvgSigmaClip(unsigned np, const float *array, unsigned nsig);
float operaArraySigma(unsigned np, const float *array);
float operaArrayWeightedSigma(unsigned np, const float *array, const float *weigh);
float operaArrayMedian(unsigned np, const float *array);
double operaArrayMedian_d(unsigned np, const double *array);
float operaArrayMedianSigma(unsigned np, const float *array, float median);
float operaArrayChisqr(unsigned np, const float *array, float central, unsigned dof);    
float operaArrayMaxValue(unsigned np, const float *array);
float operaArrayMinValue(unsigned np, const float *array);
void operaArrayHeapSort(unsigned np, float *arr);
void operaArrayIndexSort(int n, float x[], int sindex[]);
void operaArrayIndexSort_d(int n, double x[], int index[]);
float operaUniformRand(float xmin, float xmax);
float operaGaussRand(float xcen, float sig);
unsigned operaCountPixels(unsigned np, const float *array, float minvalue, float maxvalue);
double operaCrossCorrelation(unsigned np, const double *array1, const double *array2);
float operaCrossCorrelation_f(unsigned np, const float *array1, const float *array2);
operaErrorCode operaHistogram(unsigned int *outarray, unsigned int nbins, const float *inarray, unsigned int incount, unsigned int binsize, float minvalue, float maxvalue);
void weight_mean_error(unsigned nx, float *x, float *sigmax, float *xmean, float *xsigma);
float robust_sigma(unsigned n, float *inarray, unsigned reference);
float biweight_mean(unsigned n, float *vector, float *weights);
void robust_mean(unsigned n, float *vector, float *sig, float *robmean, float *robsigma);
/*!
 * float operaArrayIndexedMeanQuick(unsigned np, const float *array, const float *weigh)
 * \brief This function ....
 * \param np is an unsigned that ..,
 * \param array is a const float pointer that ..,
 * \param weigh is a const float pointer that ..,
 * \return float
 */
static inline void operaArrayIndexedMeanQuick(unsigned np, const float *array, const float *indexmask, unsigned nb, float *meanbin) {
	
	float *sumbin = (float *)malloc(sizeof(float)*nb);
	
	memset(sumbin, 0, sizeof(float)*nb);
	memset(meanbin, 0, sizeof(float)*nb);
	
	while (np--){
		if(*indexmask) {
			meanbin[(unsigned)*indexmask-1] += *array;
			sumbin[(unsigned)*indexmask-1]++;
		}
		indexmask++;
		array++;
	}
	while(nb--) {
		*meanbin++ /= *sumbin++;
	}	
}

/*! 
 * float operaArrayIndexedSigQuick(unsigned np, const float *array, const float *weigh)
 * \brief This function ....
 * \param np is an unsigned that ..,
 * \param array is a const float pointer that ..,
 * \param weigh is a const float pointer that ..,
 * \return float
 */
static inline void operaArrayIndexedSigmaQuick(unsigned np, const float *array, const float *indexmask, unsigned nb, float *sigbin) {
	
	float *sumbin = (float *)malloc(sizeof(float)*nb);
	float *meanbin = (float *)malloc(sizeof(float)*nb);	
	unsigned i=np;
	unsigned j=0;
	
	memset(sumbin, 0, sizeof(float)*nb);
	memset(meanbin, 0, sizeof(float)*nb);
	memset(sigbin, 0, sizeof(float)*nb);
	
	if(np){	
		
		while (i--){
			if(*indexmask) {
				j = (unsigned)*indexmask-1;
				meanbin[j] += *array;
				sumbin[j]++;
			}
			indexmask++;
			array++;
		}
		
		j=nb;
		while(j--) {
			meanbin[j] /= sumbin[j];
		}
		
		i=np;
		while(i--) {
			array--;	
			indexmask--;
			if(*indexmask) {
				j = (unsigned)*indexmask-1;
				sigbin[j] += (*array - meanbin[j])*(*array - meanbin[j]);
			}	
		}
		j=nb;
		while(j--) {
			sigbin[j] = sqrt(sigbin[j]/sumbin[j]);
		}		
	}
}

/*! 
 * static inline float operaArrayMedianQuick(unsigned np, float *arr)
 * \brief Optimized inline function returning a median -- WARNING: DESTRUCTIVE -- the input is modified.
 * \param arr - the array to get median from
 * \param np - number of pixels
 * \return median
 */
static inline float operaArrayMedianQuick(unsigned np, /*MODIFIED!!!*/float *arr) {
	unsigned low = 0, high = np-1; 
	unsigned median = (low + high) / 2; 
	unsigned odd = (np & 1);
	unsigned middle, ll, hh; 
	
	if (np == 0) {
		return FP_NAN;
	}
	if (np == 1) {
		return arr[0];
	}
	if (np == 2) {
		return (arr[0]+arr[1])/2.0;
	}
	while (1) { 
		
		if (high <= low) {/* One element only */ 
			return arr[median] ; 
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]) ; 
			if (odd)
				return arr[median] ; 
			else
				return (arr[low] + arr[high]) / 2.0;
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_SWAP(arr[middle], arr[high]) 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]) 
				if (arr[middle] > arr[low]) 
					ELEM_SWAP(arr[middle], arr[low]) 
					
				/* Swap low item (now in position middle) into position (low+1) */ 
					ELEM_SWAP(arr[middle], arr[low+1]) ; 
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]) ; 
			do hh--; while (arr[hh] > arr[low]) ; 
			
			if (hh < ll)
				break;
			ELEM_SWAP(arr[ll], arr[hh])
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_SWAP(arr[low], arr[hh]) 
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	} 	
}

/*! 
 * static inline float operaArrayMedianQuick(unsigned np, float *arr)
 * \brief Optimized inline function returning a median -- WARNING: DESTRUCTIVE -- the input is modified.
 * \param arr - the array to get median from
 * \param np - number of pixels
 * \return median
 */
static inline float operaArrayMedianQuick_d(unsigned np, /*MODIFIED!!!*/double *arr) {
	unsigned low = 0, high = np-1; 
	unsigned median = (low + high) / 2; 
	unsigned odd = (np & 1);
	unsigned middle, ll, hh; 
	
	if (np == 0) {
		return NAN;
	}
	if (np == 1) {
		return arr[0];
	}
	if (np == 2) {
		return (arr[0]+arr[1])/2.0;
	}
	while (1) { 
		
		if (high <= low) {/* One element only */ 
			return arr[median] ; 
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]) ; 
			if (odd)
				return arr[median] ; 
			else
				return (arr[low] + arr[high]) / 2.0;
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_SWAP_d(arr[middle], arr[high]) 
			if (arr[low] > arr[high]) 
				ELEM_SWAP_d(arr[low], arr[high]) 
				if (arr[middle] > arr[low]) 
					ELEM_SWAP_d(arr[middle], arr[low]) 
					
				/* Swap low item (now in position middle) into position (low+1) */ 
					ELEM_SWAP_d(arr[middle], arr[low+1]) ; 
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]) ; 
			do hh--; while (arr[hh] > arr[low]) ; 
			
			if (hh < ll)
				break;
			ELEM_SWAP_d(arr[ll], arr[hh])
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_SWAP_d(arr[low], arr[hh]) 
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	} 	
}

/*! 
 * static inline unsigned short operaArrayMedianQuickUSHORT(unsigned np, float *arr)
 * \brief Optimized inline function returning a median -- WARNING: DESTRUCTIVE -- the input is modified.
 * \param arr - the array to get median from
 * \param np - number of pixels
 * \return unsigned short median
 */
static inline unsigned short operaArrayMedianQuickUSHORT(unsigned np, /*MODIFIED!!!*/unsigned short *arr) {
	unsigned low = 0, high = np-1; 
	unsigned median = (low + high) / 2; 
	unsigned odd = (np & 1);
	unsigned middle, ll, hh; 
	
	if (np == 0) {
		return 0;
	}
	if (np == 1) {
		return arr[0];
	}
	if (np == 2) {
		return (arr[0]+arr[1])/2;
	}
	while (1) { 
		
		if (high <= low) {/* One element only */ 
			return arr[median] ; 
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_USHORTSWAP(arr[low], arr[high]) ; 
			if (odd)
				return arr[median] ; 
			else
				return (arr[low] + arr[high]) / 2;
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_USHORTSWAP(arr[middle], arr[high]) 
			if (arr[low] > arr[high]) 
				ELEM_USHORTSWAP(arr[low], arr[high]) 
				if (arr[middle] > arr[low]) 
					ELEM_USHORTSWAP(arr[middle], arr[low]) 
					
				/* Swap low item (now in position middle) into position (low+1) */ 
					ELEM_USHORTSWAP(arr[middle], arr[low+1]) ; 
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]) ; 
			do hh--; while (arr[hh] > arr[low]) ; 
			
			if (hh < ll)
				break;
			ELEM_USHORTSWAP(arr[ll], arr[hh])
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_USHORTSWAP(arr[low], arr[hh]) 
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	} 	
}

	/*! 
	 * float operaArrayMedianSigmaQuick(unsigned np, const float *array, float median)
	 * \brief destructive Calculate median deviation of float array
	 * \param np is an unsigned for the number of elements in array
	 * \param array is a float pointer with data
	 * \param median is a float input for the median of array
	 * \return float median deviation of array
	 */
static inline float operaArrayMedianSigmaQuick(unsigned np, float *array, float median) {
		
		/* Calculate the median deviation at each location in the array */
		for (unsigned i = 0; i < np; i++) {
			array[i] = (float)fabs(array[i] - median);
		}	
		/*
		 * Return the median deviation
		 *
		 * 0.674433 is magic number such that this deviation is the same as a
		 * classic standard deviation assuming a normal distribution function
		 * (gaussian)
		 */
		return operaArrayMedianQuick(np, array) / 0.674433;
	}
	
#ifdef __cplusplus
}
#endif
		
#endif

