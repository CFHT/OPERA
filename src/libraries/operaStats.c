/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaStats
 Version: 1.0
 Description: This C library implements basic statistics routines.
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
 along with this program.  if (not, see:
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
 * \brief Statistics library.
 * \file operaStats.c
 * \package operaStats 
 */

#include <stdlib.h>		// random, sometimes...
#include <math.h>

#ifndef random
#define random rand		// for Linux
#endif

#include "globaldefines.h"
#include "libraries/operaStats.h"

/*! 
 * operaStats
 * \author Eder Martioli
 * \brief Statistics library.
 * \details {This library contains the basic routines for statistics of arrays.}
 * \file operaStats.cpp
 * \ingroup libraries
 */

/*
 * Statistics routines for Array.
 */

/*
 * float operaArrayMean(unsigned np, const float *array)
 * \brief Calculate arithmetic mean of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \return float mean of array
 */

float operaArrayMean(unsigned np, const float *array) {
	float xmean = 0.0;
	unsigned i = 0;
	while (i < np) { xmean += *array++; i++; }
	if (!np)
		return 0.0;
	else
		return xmean/(float)np;
}

/*
 * float operaArrayTotal(unsigned np, const float *array)
 * \brief Calculate total of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \return float total of array
 */

float operaArrayTotal(unsigned np, const float *array) {
	float total = 0.0;
	while (np--) { total += *array++;}
    return total;
}

/*
 * double operaArrayMean_d(unsigned np, const double *array)
 * \brief Calculate arithmetic mean of double array
 * \param np is an unsigned for the number of elements in array
 * \param array is a double pointer with data
 * \return double mean of array
 */
double operaArrayMean_d(unsigned np, const double *array) {
	double xmean = 0.0;
	unsigned i=np;
	while (i--) xmean += *array++;
	if (!np)
		return 0.0;
	else
		return xmean/(double)np;
}

/* 
 * float operaArrayWeightedMean(unsigned np, const float *array, const float *weigh)
 * \brief Calculate weighted mean of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \param weigh is a float pointer with weights  
 * \return float mean of array
 */
float operaArrayWeightedMean(unsigned np, const float *array, const float *weights) {
	float xmean = 0.0;
	float npf = 0.0;	
	
	while (np--){
		xmean += (*array++)*(*weights);
		npf += *weights++;
	}
	if (npf == 0.0)
		return 0.0;
	else
		return xmean/npf;
}

/* 
 * float operaArrayAvgSigmaClip(unsigned np, const float *array, unsigned nsig)
 * \brief Calculate average of array clipped by nsig x sdtdev of original array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \param nsig is an unsigned for size for which data is clipped in sigma units
 * \return float avgsigclip of array
 */
float operaArrayAvgSigmaClip(unsigned np, const float *array, unsigned nsig) 
{
	float xsig = 0.0;	
	float xmean = 0.0;
	float clipedmean = 0.0;	
	const float *p = array;
	unsigned clipednp = 0;
	unsigned i=np;
	
	if (np) {
		while (i--) xmean += *p++;
		xmean/=(float)np;
		i=np;
		p=array;
		while (i--) {
			xsig += (*p - xmean)*(*p - xmean);
			p++;
		}	
		xsig = sqrt(xsig/(float)np);		
		
		i=np;
		p=array;
		while (i--) {
			if((*p > xmean - (float)nsig*xsig) && (*p < xmean + (float)nsig*xsig)) {
				clipedmean += *p;
				clipednp++;
			}
			p++;			
		}			
		return(clipedmean/(float)clipednp);	
	}
	return 0.0;
}

/* 
 * float operaArraySigma(unsigned np, const float *array)
 * \brief Calculate standard deviation around mean of array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \return float sigma
 */
float operaArraySigma(unsigned np, const float *array) {
	
	float xsig = 0.0;	
	float xmean = 0.0;
	const float *p = array;
	unsigned i=np;
	
	if (np) {
		while (i--) {
			xmean += *p++;
		}
		xmean /= (float)np;
		i=np;
		p = array;
		while (i--) {
			xsig += (*p - xmean)*(*p - xmean);
			p++;			
		}	
		xsig /= (float)np;
	}
	return sqrt(xsig);	
}

/* 
 * float operaArrayWeightedSigma(unsigned np, const float *array, const float *weigh)
 * \brief Calculate weighted standard deviation of weighted mean of array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \param weigh is a float pointer with weights  
 * \return float weighted sigma
 */
float operaArrayWeightedSigma(unsigned np, const float *array, const float *weights) {
	float xmean = 0.0;
	float xsig = 0.0;
	float npf = 0.0;	
	unsigned i=np;
	
	if (np) {
		while (i--){
			xmean += (*array)*(*weights);		
			npf += *weights;
			array++;
			weights++;			
		}	
		xmean /= npf;				// NOTE: This may throw a floating point exception if (weigh is all zeroes
		i=np;
		while (i--){
			array--;
			weights--;			
			xsig += (*weights)*(*array - xmean)*(*array - xmean);
		}	
		xsig /= npf;
	}
	return sqrt(xsig);	
}

/* 
 * float operaArrayMedian(unsigned np, const float *array)
 * \brief This function non-destructively finds the median of an array.
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \return float median of array
 */
float operaArrayMedian(unsigned np, const float *array) {
	unsigned low = 0, high = np-1; 
	unsigned median = (low + high) / 2; 
	unsigned odd = (np & 1);
	unsigned middle, ll, hh; 
	
	if (np == 0) {
		return FP_NAN;
	}
	if (np == 1) {
		return array[0];
	}
	if (np == 2) {
		return (array[0]+array[1])/2.0;
	}
	float *arr = (float *)malloc(np * sizeof(float)); 	
	if (arr) {
		memcpy (arr, array, np * sizeof(float));
	while (1) { 
		
		if (high <= low) { /* One element only */ 
			float ret = arr[median];  
			free(arr);
			return ret;
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]);  
			if (odd) {
				float ret = arr[median];  
				free(arr);
				return ret;
			} else {
				float ret = (arr[low] + arr[high]) / 2.0;
				free(arr);
				return ret;
			}
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_SWAP(arr[middle], arr[high]); 
		if (arr[low] > arr[high]) 
			ELEM_SWAP(arr[low], arr[high]);
		if (arr[middle] > arr[low]) 
			ELEM_SWAP(arr[middle], arr[low]);
					
		/* Swap low item (now in position middle) into position (low+1) */ 
		ELEM_SWAP(arr[middle], arr[low+1]);  
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]);  
			do hh--; while (arr[hh] > arr[low]);  
			
			if (hh < ll)
				break;
			ELEM_SWAP(arr[ll], arr[hh]);
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_SWAP(arr[low], arr[hh]);
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	}
	}
	return FP_NAN;
}

/* 
 * double operaArrayMedian(unsigned np, const double *array)
 * \brief This function non-destructively finds the median of an array.
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \return float median of array
 */
double operaArrayMedian_d(unsigned np, const double *array) {
	unsigned low = 0, high = np-1; 
	unsigned median = (low + high) / 2; 
	unsigned odd = (np & 1);
	unsigned middle, ll, hh; 
	
	if (np == 0) {
		return FP_NAN;
	}
	if (np == 1) {
		return array[0];
	}
	if (np == 2) {
		return (array[0]+array[1])/2.0;
	}
	double *arr = (double *)malloc(np * sizeof(double)); 	
	if (arr) {
		memcpy (arr, array, np * sizeof(double));
	while (1) { 
		
		if (high <= low) {/* One element only */ 
			double ret = arr[median];  
			free(arr);
			return ret;
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_SWAP_d(arr[low], arr[high]);  
			if (odd) {
				double ret = arr[median];  
				free(arr);
				return ret;
			} else {
				double ret = (arr[low] + arr[high]) / 2.0;
				free(arr);
				return ret;
			}
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_SWAP_d(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])
				ELEM_SWAP_d(arr[low], arr[high]) ;
		if (arr[middle] > arr[low]) 
			ELEM_SWAP_d(arr[middle], arr[low]) ;
					
		/* Swap low item (now in position middle) into position (low+1) */ 
		ELEM_SWAP_d(arr[middle], arr[low+1]);  
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]);  
			do hh--; while (arr[hh] > arr[low]);  
			
			if (hh < ll)
				break;
			ELEM_SWAP_d(arr[ll], arr[hh]);
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_SWAP_d(arr[low], arr[hh]) 
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	} 
	}
	return FP_NAN;
}

/* 
 * float operaArrayMedianSigma(unsigned np, const float *array, float median)
 * \brief Calculate median deviation of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \param median is a float input for the median of array
 * \return float median deviation of array
 */
float operaArrayMedianSigma(unsigned np, const float *array, float median) {
	unsigned i;
	if (np == 0) {
		return FP_NAN;
	}
	if (np == 1) {
		return (array[0] - median )/ 0.674433;
	}
	if (np == 2) {
		return ((array[0] - median) + (array[1] - median))/2.0 / 0.674433;
	}
	float *arr = (float *)malloc(np * sizeof(float)); 
	if (arr) {
		memcpy (arr, array, np * sizeof(float));
		
		/* Calculate the median deviation at each location in the array */
		for (i = 0; i < np; i++) {
			arr[i] = (float)fabs(arr[i] - median);
		}	
		/*
		 * Return the median deviation
		 *
		 * 0.674433 is magic number such that this deviation is the same as a
		 * classic standard deviation assuming a normal distribution function
		 * (gaussian)
		 */
		float medsig =  operaArrayMedian(np, arr) / 0.674433;
		free(arr);
		return medsig;		
	}
	return FP_NAN;
}


/* 
 * float operaArrayChisqr(unsigned np, const float *central, float median, unsigned dof)
 * \brief Calculate the chisqr of float array with respect to a central value
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data
 * \param central is a float input for the central values of array
 * \param np is an unsigned for the number of degrees of freedom 
 * \return float chisqr of array
 */
float operaArrayChisqr(unsigned np, const float *array, float central, unsigned dof) {
	float chisqr = 0.0;
	unsigned i=np;
	while (i--) {
        chisqr += (*array - central)*(*array - central); 
        array++;
    }
    
	if (!np)
		return 0.0;
	else
        return chisqr/(float)dof;
}

/* 
 * float operaArrayMaxValue(unsigned np, const float *array)
 * \brief Find maximum value of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \return float maximum value of array
 */
float operaArrayMaxValue(unsigned np, const float *array) {
	float xmax = -3.4e+38;
	
	while (np--) {
		if(*array > xmax)
			xmax = *array;	
		array++;
	}	
	return xmax;	
}

/* 
 * float operaArrayMinValue(unsigned np, const float *array)
 * \brief Find minimum value of float array
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \return float minimum value of array
 */
float operaArrayMinValue(unsigned np, const float *array) {
	float xmin = 3.4e+38;
	while (np--){
		if(*array < xmin)
			xmin = *array;	
		array++;
	}	
	return xmin;
}

/*
 * Statistical tools
 */

/* 
 * void operaArrayHeapSort(unsigned np, float *arr)
 * \brief This function destructively heap sorts an array.
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \return void
 */
void operaArrayHeapSort(unsigned np, float *arr) { // based on heapsort
    unsigned int n = np, i = n/2, parent, child;
    int t;	 
	
    for (;;) { /* Loops until arr is sorted */
		if (i > 0) { /* First stage - Sorting the heap */
			i--;           /* Save its index to i */
			t = arr[i];    /* Save parent value to t */
		} else {     /* Second stage - Extracting elements in-place */
			n--;           /* Make the new heap smaller */
			if (n == 0) return; /* When the heap is empty, we are done */
			t = arr[n];    /* Save last value (it will be overwritten) */
			arr[n] = arr[0]; /* Save largest value at the end of arr */
		}
		
		parent = i; /* We will start pushing down t from parent */
		child = i*2 + 1; /* parent's left child */
		
		/* Sift operation - pushing the value of t down the heap */
		while (child < n) {
			if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
				child++; /* Choose the largest child */
			}
			if (arr[child] > t) { /* if (any child is bigger than the parent */
				arr[parent] = arr[child]; /* Move the largest child up */
				parent = child; /* Move parent pointer to this child */
				child = parent*2 + 1; /* Find the next child */
			} else {
				break; /* t's place is found */
			}
		}
		arr[parent] = t; /* We save t in the heap */
    }
}

/* 
 * void operaArrayIndexSort(unsigned n, const float *x, unsigned *sindex) 
 * \brief This function sorts a float array in increasing order using index association (doesn't change original array)
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \param sindex is an unsigned pointer for the index association of the sorted array
 * \return void
 */
//void spline_dsrt(int n, double x[], int index[])
void operaArrayIndexSort(int n, float x[], int index[])
{
    for (unsigned i=0; i < n; i++) index[i]=i;
    int sifting, parent, current, child;
    parent=n/2;
    for (int heapsize=n; heapsize > 1;) {
		// Still "heapifying" the array - starting from last parent node, ending at root
        if (parent > 0) {
            sifting=index[--parent]; // Shift parent back by one then mark it as the node to sift down
        }
        
        // Sort the heap by swapping the root to the end of the remaining heap, then repair the heap
        else {
			heapsize--;
            sifting=index[heapsize]; // Mark previous last node as value to be sifted down
            index[heapsize]=index[0]; // Take root (largest remaining value) and move to previous end
        }
        
        // Sift marked value down
        current = parent;
        child = 2*current+1; // Left child of current
        while (child < heapsize) {
            if (child+1 < heapsize && x[index[child]] < x[index[child+1]]) child++; // Select the larger child
            if (x[sifting] >= x[index[child]]) break; // If the value being sifted down is larger than the children, we are done sifting
			index[current]=index[child]; // Move the value of the larger child up, then go down and repeat
			current = child;
			child = 2*current+1;
        }
        index[current]=sifting; // Take the value being sifted and store it in the selected location
    }
}

/* 
 * void operaArrayIndexSort_d(unsigned n, const double *x, unsigned *sindex) 
 * \brief This function sorts a double array in increasing order using index association (doesn't change original array)
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \param sindex is an unsigned pointer for the index association of the sorted array
 * \return void
 */
//void spline_dsrt(int n, double x[], int index[])
void operaArrayIndexSort_d(int n, double x[], int index[])
{
    for (unsigned i=0; i < n; i++) index[i]=i;
    int sifting, parent, current, child;
    parent=n/2;
    for (int heapsize=n; heapsize > 1;) {
		// Still "heapifying" the array - starting from last parent node, ending at root
        if (parent > 0) {
            sifting=index[--parent]; // Shift parent back by one then mark it as the node to sift down
        }
        
        // Sort the heap by swapping the root to the end of the remaining heap, then repair the heap
        else {
			heapsize--;
            sifting=index[heapsize]; // Mark previous last node as value to be sifted down
            index[heapsize]=index[0]; // Take root (largest remaining value) and move to previous end
        }
        
        // Sift marked value down
        current = parent;
        child = 2*current+1; // Left child of current
        while (child < heapsize) {
            if (child+1 < heapsize && x[index[child]] < x[index[child+1]]) child++; // Select the larger child
            if (x[sifting] >= x[index[child]]) break; // If the value being sifted down is larger than the children, we are done sifting
			index[current]=index[child]; // Move the value of the larger child up, then go down and repeat
			current = child;
			child = 2*current+1;
        }
        index[current]=sifting; // Take the value being sifted and store it in the selected location
    }
}

/* 
 * float operaUniformRand(float xmax, float xmin) 
 * \brief This function produces a random number with uniform distribution.
 * \note For better randomness initialize srandom(time(NULL)) prior to call this function
 * \param xmin is a float to define the lower limit of the distribution (rand > xmin)
 * \param xmax is a float to define the higher limit of the distribution (rand < xmax)
 * \return float random number
 */
float operaUniformRand(float xmin, float xmax) {
    float randBetweenZeroAndOne = (float)random()/(float)RAND_MAX;
    float randBetweenXmiFP_NANdXmax = xmin + (xmax - xmin)*randBetweenZeroAndOne;
	return randBetweenXmiFP_NANdXmax;
}

/* 
 * float operaGaussRand(float xcen, float sig)
 * \brief This function produces a random number with normal distribution.
 * \note For better randomness initialize srandom(time(NULL)) prior to call this function
 * \param xcen is a float to define the center of the distribution
 * \param sig is a float to define the spread of the distribution (symmetric)
 * \return float random number
 */
float operaGaussRand(float xcen, float sig) {
	float y1,y2,x1,x2;
	
	x1 = ((float)random()/(float)RAND_MAX);
	x2 = ((float)random()/(float)RAND_MAX);
	
	y1 = sqrt( - 2. * log(x1) ) * cos( 2.*M_PI*x2 );
	y2 = sqrt( - 2. * log(x1) ) * sin( 2.*M_PI*x2 );
	
	return(xcen + sig*(y1+y2)/2.);
}

/* 
 * operaCountPixels(unsigned np, const float *array, float minvalue, float maxvalue)
 * \brief Count number of elements in array lying between minvalue and maxvalue
 * \param np is an unsigned for the number of elements in array
 * \param array is a float pointer with data 
 * \param minvalue is a float for the minimum value to count an element
 * \param maxvalue is a float for the maximum value to count an element
 * \return unsigned number of elements
 */

unsigned operaCountPixels(unsigned np, const float *array, float minvalue, float maxvalue) {
	unsigned count = 0;
	while (np--) {
		if(*array > minvalue && *array < maxvalue) {
			count++;
		}
		array++;
	}	
	return count;
}

/* 
 * double operaCrossCorrelation(unsigned np, const double *array1, const double *array2)
 * \brief Calculate cross-correlation between two arrays
 * \param np is an unsigned for the number of elements in array
 * \param array1 is a double pointer with data 
 * \param array2 is a double pointer with data  
 * \return double value for cross-correlation 
 */

double operaCrossCorrelation(unsigned np, const double *array1, const double *array2) {
    if (np) {
        double mean1 = 0.0;	
        double mean2 = 0.0;		
        double numerator = 0;
        double denominator1 = 0;    
        double denominator2 = 0;
        
        unsigned i=np;
        
        while (i--){
            mean1 += *array1++;
            mean2 += *array2++;        
        }
        mean1 /= (double)np;
        mean2 /= (double)np;        
        i=np;
        while (i--){
            array1--;
            array2--;       
            numerator += (*array1 - mean1) * (*array2 - mean2);
            denominator1 += (*array1 - mean1) * (*array1 - mean1);
            denominator2 += (*array2 - mean2) * (*array2 - mean2);
        }        
        
        double crosscorrelation = numerator/(sqrt(denominator1)*sqrt(denominator2));
        return crosscorrelation;
   } else {
        return EXIT_FAILURE;
    }
}

/*
 * double operaCrossCorrelation_f(unsigned np, const float *array1, const float *array2)
 * \brief Calculate cross-correlation between two arrays
 * \param np is an unsigned for the number of elements in array
 * \param array1 is a float pointer with data
 * \param array2 is a float pointer with data
 * \return float value for cross-correlation
 */

float operaCrossCorrelation_f(unsigned np, const float *array1, const float *array2) {
    if (np) {
        float mean1 = 0.0;
        float mean2 = 0.0;
        float numerator = 0;
        float denominator1 = 0;
        float denominator2 = 0;
        
        unsigned i=np;
        
        while (i--){
            mean1 += *array1++;
            mean2 += *array2++;
        }
        mean1 /= (float)np;
        mean2 /= (float)np;
        i=np;
        while (i--){
            array1--;
            array2--;
            numerator += (*array1 - mean1) * (*array2 - mean2);
            denominator1 += (*array1 - mean1) * (*array1 - mean1);
            denominator2 += (*array2 - mean2) * (*array2 - mean2);
        }
        
        float crosscorrelation = numerator/(sqrt(denominator1)*sqrt(denominator2));
        return crosscorrelation;
    } else {
        return EXIT_FAILURE;
    }
}

/*
 * void operaHistogram(const float *inarray, unsigned int *outarray, float binsize, unsigned int nbins, float minvalue, float maxvalue)
 * \brief computes the density function of an array. The array values must be positive,
 *  minvalue and maxvalue must also be positive.
 * \param outarray is the output array, must be an array of unsigned ints
 * \param nbins is the number of elements in outarray (bins), must be an unsigned int
 * \param inarray is the input array
 * \param incount is the number of elements in inarray
 * \param binsize is the size of the bins
 * \param minvalue is the minimum value to be included in the histogram
 * \param maxvalue is the maximum value to be included in the histogram
 * \return array equal to the density function of inarray
*/
operaErrorCode operaHistogram(unsigned int *outarray, unsigned int nbins, const float *inarray, unsigned int incount, unsigned int binsize, float minvalue, float maxvalue) {
	unsigned int bin;
	if (maxvalue < 0.0 || minvalue < 0 || minvalue >= maxvalue) {
		return operaErrorInputValuesDisagree; // Check that binsize, nbins, minvalue, maxvalue work together. if (not, return error
	}
	if ((maxvalue-minvalue)/binsize != nbins) {
		return operaErrorInputValuesDisagree; // Check that binsize, nbins, minvalue, maxvalue work together. if (not, return error
	}
	memset(outarray, 0, sizeof(unsigned int)*nbins); // Set all values in outarray to 0
	while (incount) { // Sort each element of inarray in to bins
		if (inarray[incount] >= minvalue && inarray[incount] <= maxvalue) {	// check if (value is in within minvalue and maxvalue
			bin = inarray[incount]/binsize; // rounds down, taking in to account that arrays start at 0
			outarray[bin]++; // increment the corresponding bin (outarray element)
		}
		incount--;
	}
	return 0; // No errors
}

//	Calculate the center and dispersion (like mean and sigma) of a 
//	distribution using bisquare weighting.
//
//	Sigma = An outlier-resistant measure of the dispersion about the 
//	      center, analogous to the standard deviation. 
//
//	Weights = The weights applied to the data in the last iteration, 
//                 floating point vector
//
// NOTES:
//       Since a sample mean  scaled by sigma/sqrt(N), has a Student's T 
//       distribution, the half-width of the  95% confidence interval for 
//       the sample mean  can be determined as follows: 
//          ABS( T_CVF( .975, .7*(N-1) )*SIGMA/SQRT(N) ) 
//       where N = number of  points, &&  0.975 = 1 - (1 - 0.95)/2. 
float biweight_mean(unsigned n, float *vector, float *weights) {
	const unsigned maxiterations = 20; // Allow 20 iterations, this should nearly always be sufficient
	float eps = 1.0e-24;
	
	float close_enough = 0.03*sqrt(0.5/(n-1)); // compare to fractional change in width
	
	float diff = 1.0e30;
	unsigned iteration = 0;
	
	// As an initial estimate of the center, use the median:
	float y0 = operaArrayMedian(n, vector);
	
	// Calculate the weights:
	float *deviation = malloc(n*sizeof(float));
	for (unsigned i=0; i<n; i++) {
		deviation[i] = vector[i]-y0;
	}
	float sigma = robust_sigma(n, deviation, 0);
	
	if (sigma < eps) {
        float limit = 3.0*sigma;
        for (unsigned i=0; i<n; i++) {
            weights[i] = (abs(deviation[i]) <= limit?1.0:0.0);
        }
		diff = 0.0; // (skip rest of routine)
	}
	
	// Repeat:
	float *tmparray = malloc(n*sizeof(float));
    float prev_sigma = 0.0;
	while ( (diff > close_enough) && (iteration++ < maxiterations) ) {
        float totalweight = 0.0;
        for (unsigned i=0; i<n; i++) {
            tmparray[i] = pow((vector[i]-y0)/(6.0*sigma), 2);
            tmparray[i] = (tmparray[i] < 1.0 ? 1.0 : 0.0);
            weights[i] = pow(1.0-tmparray[i], 2);
            totalweight += weights[i];
		}
        for (unsigned i=0; i<n; i++) {
            weights[i] = weights[i]/totalweight;
        }
        y0 = 0.0;
        for (unsigned i=0; i<n; i++) {
            y0 += weights[i] * vector[i];
        }
        for (unsigned i=0; i<n; i++) {
            deviation[i] = vector[i] - y0;
        }
        prev_sigma = sigma;
        sigma = robust_sigma(n, deviation, 1);
        if (sigma > eps) {
            diff = fabs(prev_sigma-sigma)/prev_sigma;
        } else {
            diff = 0.0;
        }
	}
    free(tmparray);
    free(deviation);
	return y0;
}
//
// This calculates a robust, outlier resistant MEAN with
// optional weighting.
//
// INPUTS:
//	vector    The vector of values for which to compute the mean of.
//	=sig      Optional uncertainties in the vector values.  These
//	          will be used for weighting.
//
// OUTPUTS:
//	robmean   The robust mean of the vector values.
//	robsigma   The robust standard deviation of vector.
//
void robust_mean(unsigned n, float *vector, float *sigma_in, float *robustmean, float *robustsigma) {
	unsigned maxiteration = 20;
	float eps = 1.0e-24;
	
	// Do we have uncertainties?
	// if not then unweighted
    float *sigma = sigma_in;
    if (sigma == NULL) {
        sigma = malloc(sizeof(float)*n);
        for (unsigned i=0; i<n; i++) {
            sigma[i] = 1.0;
        }
		*robustmean = operaArrayMedian(n , vector);
		*robustsigma = robust_sigma(n, vector, 0);
	} else {	// No weights, use median/robust_sigma
		weight_mean_error(n, vector, sigma_in, robustmean, robustsigma);
	}
	
	// Iterate until it converged
	char done = 0;
    float *inliers = malloc(sizeof(float)*n);
    float *inliersigma = malloc(sizeof(float)*n);
    unsigned count = 0;
    unsigned old_count = 0;
    float new_mean, new_sigma;
	while (!done) {
		// Remove outliers from the whole array
		// Make sure we have a decent SIGMA
		if (*robustsigma > 0.0) {
            for (unsigned i=0; i<n; i++) {
                if (vector[i] <= 2.5 * *robustsigma) {
                    inliers[count] = vector[i];
                    inliersigma[count] = sigma[i];
                    count++;
                }
            }
		}
		// No good points left, loosen sigma limit
		if (count == 0) {
            for (unsigned i=0; i<n; i++) {
                if (vector[i] <= 3.0 * *robustsigma) {
                    inliers[count] = vector[i];
                    inliersigma[count] = sigma[i];;
                    count++;
                }
            }
		}
		// No good points left, loosen sigma limit again
		if (count == 0) {
            for (unsigned i=0; i<n; i++) {
                if (vector[i] <= 4.0 * *robustsigma) {
                    inliers[count] = vector[i];
                    inliersigma[count] = sigma[i];;
                    count++;
                }
            }
		}
		// No good points left, bail
		if (count == 0) {
            *robustsigma = 0.0;
            *robustmean = 0.0;
            free(inliers);
            free(inliersigma);
            return;
		}
		//weight_mean_error(unsigned nx, float *x, float *sigmax, float *xmean, float *xsigma)
		// Use WMEANERR to get weighted mean
		weight_mean_error(count, inliers, inliersigma, &new_mean, &new_sigma);
		
		// What's the change
		float dmean = fabs(*robustmean-new_mean);
		float dsigma = fabs(*robustsigma-new_sigma);
		unsigned dcount = abs(count-old_count);
		
		// Are we done?
		if ((dmean < eps && dsigma < eps && dcount == 0 && count > 0) || (count > maxiteration)) {
			done = 1;
		}
		
		// New best estimate
		old_count = count;
		*robustmean = new_mean;
		*robustsigma = new_sigma;
	}
    free(inliers);
    free(inliersigma);
}
// 	Calculate a resistant estimate of the dispersion of a distribution.
//  EXPLANATION:
// 	For an uncontaminated distribution, this is identical to the standard
// 	deviation.
// 
//  CALLING SEQUENCE:
// 	result = robust_sigma( Y, [ /ZERO, GOODVEC = ] )
// 
//  INPUT: 
// 	Y = Vector of quantity for which the dispersion is to be calculated
// 
//  OPTIONAL INPUT KEYWORD:
// 	/ZERO - if (set, the dispersion is calculated w.r.t. 0.0 rather than the
// 		central value of the vector. if (Y is a vector of residuals, this
// 		should be set.
// 
//  OUTPUT:
// 	robust_sigma returns the dispersion. In case of failure, returns 
// 	value of -1.0
// 
//  PROCEDURE:
// 	Use the median absolute deviation as the initial estimate,) {weight 
// 	points using Tukey's Biweight. See, for example, "Understanding Robust
// 	and Exploratory Data Analysis," by Hoaglin, Mosteller and Tukey, John
// 	Wiley & Sons, 1983, or equation 9 in Beers et al. (1990, AJ, 100, 32)
// 

float robust_sigma(unsigned n, float *inarray, unsigned reference) {
	float eps = 1.0e-20;
	float sigma = 0.0;
	float *tmparray = malloc(sizeof(float)*n);
	
	if (!reference) {
		float y0  = operaArrayMedian(n, inarray);
		for (unsigned i=0; i<n; i++) {
			inarray[i] -= y0;
		}
	}
	
	//  First, the median absolute deviation median_absolute_deviation about the median:
	float median_absolute_deviation = fabs(operaArrayMedian(n, inarray) / 0.6745);
	
	//  if (the median_absolute_deviation=0, try the MEAN absolute deviation:
	if (median_absolute_deviation < eps) {
		median_absolute_deviation = fabs(operaArrayMean(n, inarray) / 0.80);
	}
	if (median_absolute_deviation < eps) {
		free(tmparray);
		return sigma;
	}
	
	//  Now the biweighted value:
	
	unsigned count = 0;
	for (unsigned i=0; i<n; i++) {
		tmparray[i] = inarray[i] / 6.0 * median_absolute_deviation;
		tmparray[i] *= tmparray[i] ;
		if (tmparray[i] <= 1.0) {
			count++;
		}
	}
	if (count < 3) {
		free(tmparray);
		return -1.0;
	}
	
	float numerator_total = 0.0;
	float denominator_total = 0.0;
	for (unsigned i=0; i<n; i++) {
		if (tmparray[i] <= 1.0) {
			numerator_total += pow(inarray[i], 2) * pow(1.0-tmparray[i], 4);
			denominator_total += (1.0 - tmparray[i]) * (1.0 - 5.0 * tmparray[i]);
		}
	}
	sigma = n * numerator_total / (denominator_total*(denominator_total-1.0));
	
	free(tmparray);
	
	return (sigma > 0.0?sqrt(sigma):0.0);
}
//	Calculate the mean and estimated error for a set of weighted data points
// DESCRIPTION:
//	This routine is adapted from Program 5-1, XFIT, from "Data Reduction
//	and Error Analysis for the Physical Sciences", p. 76, by Philip R.
//	Bevington, McGraw Hill.  This routine computes the weighted mean using
//	Instrumental weights (w=1/sigma^2).  The final uncertainty is
//	insensitive to a multiplicative constant on the weights.
// INPUTS:
//	x      - Array of data points
//	sigmax - Array of errors in x
// OUTPUTS:
//	xmean   - weighted mean
//	xsigma  - uncertainty of xmean.  Weighted St.Dev.
//            NOT Stdev. of the Mean.
void weight_mean_error(unsigned nx, float *x, float *sigmax, float *xmean, float *xsigma) {
	if (nx == 0) {
		*xmean  = 0.0;
		*xsigma = 0.0;
    } else if (nx == 1) {
        *xmean  = x[0];
        *xsigma = sigmax[0];
	} else {
        float *weight = malloc(sizeof(float)*nx);
        for (unsigned i=0; i<nx; i++) {
            weight[i] = 1.0 / pow(sigmax[i],2);
        }
		float sum = operaArrayTotal(nx, weight);
		if (sum == 0.0) {
			sum = 1.0;
            for (unsigned i=0; i<nx; i++) {
                weight[i] = 1.0;
            }
		}
        float sumx = 0.0;
        for (unsigned i=0; i<nx; i++) {
            sumx += weight[i] * x[i];
        }
		*xmean  = sumx / sum;
        float total = 0.0;
        for (unsigned i=0; i<nx; i++) {
            total += weight[i] * pow(x[i]-*xmean, 2);
        }
		*xsigma = sqrt(total * nx / (( nx-1.0) * sum));
        free(weight);
	}
	
}
