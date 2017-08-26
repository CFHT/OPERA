/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaIMage
 Version: 1.0
 Description: This C library implements low level image routines..
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

/* \brief  Image routines */

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"

/*!
 * \brief image library.
 * \file operaImage.c
 */

/*!
 * operaImage
 * \author Doug Teeple & Eder Martioli
 * \brief image library.
 * \brief Image routines
 * \ingroup libraries
 */

/* 
 * void operaImMean(unsigned depth, long npixels, float *master, float *arrays[])
 * \brief Mean combine a series of arrays (images) into the master
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param master is a float pointer that returns the resulting image
 * \param arrays is an array of float pointers that contains the input images
 * \return void
 */
void operaImMean(unsigned depth, long npixels, float *master, float *arrays[]) 
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

/* 
 * void operaImWeightedMean(unsigned depth, long npixels, float *master, float *arrays[], float *weights[])
 * \brief Mean combine a series of arrays (images) weighted by weights into the master.
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param master is a float pointer that returns the resulting image
 * \param arrays is an array of float pointers that contains the input images 
 * \param weights is an array of float pointers that contains the input weight images
 * \return void
 */
void operaImWeightedMean(unsigned depth, long npixels, float *master, float *arrays[], float *weights[]) 
{
	float *weightsum = (float *)malloc(sizeof(float)*npixels);
	if (weightsum) {
		memset(weightsum, 0, sizeof(float)*npixels);
		
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
}	

/* 
 * void operaImSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master)
 * \brief Calculate the standard deviation image from a series of images with respect to a master image.
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param sigarray is a float pointer that returns the resulting standard deviation image 
 * \param arrays is an array of float pointers that contains the input images 
 * \param master is a float pointer for the input master image
 * \return void
 */
void operaImSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master) 
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

/* 
 * void operaImWeightedSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *master)
 * \brief Calculate the weighted standard deviation image from a series of images with respect to a master image.
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param sigarray is a float pointer that returns the resulting standard deviation image 
 * \param arrays is an array of float pointers that contains the input images 
 * \param arrays is an array of float pointers that contains the input weight images
 * \param master is a float pointer for the input master image
 * \return void
 */
void operaImWeightedSig(unsigned depth, long npixels, float *sigarray, float *arrays[], float *weights[], float *master) 
{
	float *weightsum = (float *)malloc(sizeof(float)*npixels);
	if (weightsum) {
		memset(weightsum, 0, sizeof(float)*npixels);
		
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
}	

/* 
 * void operaImAvgSigClip(unsigned depth, long npixels, float *master, float *arrays[], unsigned nsig)
 * \brief Average sigma clip combine a series of arrays (images) into the master.
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param master is a float pointer that returns the resulting image
 * \param arrays is an array of float pointers that contains the input images 
 * \param nsig is an unsigned for the size to clip data in sigma units 
 * \return void
 */
void operaImAvgSigClip(unsigned depth, long npixels, float *master, float *arrays[], unsigned nsig) 
{
	long np;
	unsigned d;
	
	float *avg = (float *)malloc(sizeof(float)*npixels);
	if (avg) {
		memset(avg, 0, sizeof(float)*npixels);	
		
		float *sigarray = (float *)malloc(sizeof(float)*npixels);
		if (sigarray) {
			memset(sigarray, 0, sizeof(float)*npixels);	
			
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
			if (weightsum) {
				memset(weightsum, 0, sizeof(float)*npixels);	
				
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
		}
	}
}
/* 
 * void operaImVarDiff(unsigned depth, long npixels, float *arrays[], float *diffvarimg)
 * \brief Calculate the variance of the difference between consecutive images.
 * \param depth is an unsigned for the number of input images
 * \param npixels is a long for the number of elements in array 
 * \param arrays is an array of float pointers that contains the input images 
 * \param diffvarimg is a float pointer that returns the resulting image 
 * \return void
 */
void operaImVarDiff(unsigned depth, long npixels, float *arrays[], float *diffvarimg) {
	
	long np;
	
	for (unsigned d=0; d<depth-1; d++) {
		
		float *arr2 = arrays[d+1];		
		float *arr1 = arrays[d];
		float *diff = diffvarimg;	
		
		np = npixels;
		
		while (np--){
			*diff++ += (*arr2 - *arr1)*(*arr2 - *arr1);
			arr1++;
			arr2++;			
		}		
	}
	np = npixels;
	float *diff = diffvarimg;	
	
	while (np--) {
		*diff++ /= (float)(depth - 1);	
	}	
}

/* 
 * float operaCCDVarDiff(unsigned depth, long npixels, float *arrays[], float *weight)
 * \brief ....
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param arrays is a float pointer pointer that ...
 * \param weight is a float pointer that ...
 * \return float
 */
float operaCCDVarDiff(unsigned depth, long npixels, float *arrays[], float *weight) 
{	
	float *diffarr = (float *)malloc(sizeof(float)*npixels);
	if (diffarr) {
		memset(diffarr, 0, sizeof(float)*npixels);
		
		long np;
		
		for (unsigned d=0; d<depth-1; d++) {
			
			float *arr2 = arrays[d+1];		
			float *arr1 = arrays[d];
			float *diff = diffarr;	
			np = npixels;
			
			while (np--){
				*diff++ += (*arr2++ - *arr1++);
			}		
		}
		np = npixels;
		float *diff = diffarr;	
		float vardiff = 0;
		float *w = weight;	
		float weightsum = 0;	
		
		while (np--) {
			vardiff += (*w)*(*diff/(float)(depth - 1))*(*diff/(float)(depth - 1));	
			weightsum += *w++;
			diff++;
		}	
		free(diffarr);
		return vardiff/weightsum;
		
	}
	return FP_NAN;
}	

/* 
 * float *medianCombineFloat(unsigned depth, long npixels, float *output, float *images[])
 * \brief median combine a series of images.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param master is a float pointer that ...
 * \param arrays is a float pointer pointer that ...
 * \return float*
 */
float *medianCombineFloat(unsigned depth, long npixels, float *master, float *arrays[]) {
	float *pfi = malloc(depth*sizeof(float));
	if (pfi) {
		long np = 0;
		float *p = master;
		while (np < npixels) {
			for (unsigned d=0; d<depth; d++) {
				float *pfd = arrays[d];
				pfi[d] = pfd[np];
			}
			*p++ = operaArrayMedianQuick(depth, pfi);
			np++;
		}
		free(pfi);
		return master;		
	}
	return NULL;
}

/* 
 * unsigned short *operaArrayMedianCombineUSHORT(unsigned depth, long npixels, unsigned short *master, unsigned short *arrays[])
 * \brief median combine a series of images.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param master is a  unsigned short pointer that ...
 * \param arrays is a  unsigned short pointer pointer that ...
 * \return unsigned short *
 */
unsigned short *operaArrayMedianCombineUSHORT(unsigned depth, long npixels, unsigned short *master, unsigned short *arrays[]) {
	unsigned short *pfi = malloc(depth*sizeof(unsigned short));
	if (pfi) {
		long np = 0;
		unsigned short *p = master;
		while (np < npixels) {
			for (unsigned d=0; d<depth; d++) {
				unsigned short *pud = arrays[d];
				pfi[d] = pud[np];
			}
			*p++ = operaArrayMedianQuickUSHORT(depth, pfi);
			np++;
		}
		free(pfi);
		return master;
	}
	return NULL;
}

/* 
 * unsigned short *operaArrayMedianCombineParUSHORT(unsigned depth, long npixels, unsigned short *output, unsigned short *images[])
 * \brief parallel median combine a series of images.
 * \note creates 8,000,000 threads.
 * \note not used.
 * \param depth is an unsigned that ...
 * \param npixels is a long that ...
 * \param output is a  unsigned short pointer that ...
 * \param images is a  unsigned short pointer pointer that ...
 * \return unsigned short *
 */
#include <pthread.h>

typedef struct thread_args {
	unsigned short median;		// return value
	unsigned depth;			// depth
	unsigned short *parr;		// stack of pixels at this x,y
} thread_args_t;

void *operaParArrayMedian(void *argument) {
	thread_args_t *thread_args_s = (thread_args_t *)argument;
	thread_args_s->median = operaArrayMedianQuickUSHORT((unsigned)thread_args_s->depth, (unsigned short*)thread_args_s->parr);
	return NULL;
}
unsigned short *operaArrayMedianCombineParUSHORT(unsigned depth, long npixels, unsigned short *output, unsigned short *images[]) {
	long np = npixels;
	unsigned short *p = output;
	pthread_t *threads = (pthread_t *)malloc(np*sizeof(pthread_t*));
	thread_args_t *thread_args = (thread_args_t *)malloc(np*sizeof(thread_args_t));
	int rc;
	
	/* create all threads */
	for (unsigned i=0; i<np; ++i) {
		unsigned short *pfi = (unsigned short *)malloc(depth*sizeof(unsigned short));
		for (unsigned d=0; d<depth; d++) {
			unsigned short *pud = images[d];
			pfi[d] = pud[np];
		}
		thread_args[i].median = 0;	// the median
		thread_args[i].depth = depth;	// the depth
		thread_args[i].parr = pfi;	// the array of pixels in the stack
		if ((rc = pthread_create(&threads[i], NULL, operaParArrayMedian, (void *) &thread_args[i])) != 0)
			return NULL;
	}
	
	/* wait for all threads to complete */
	for (unsigned i=0; i<np; ++i) {
		if ((rc = pthread_join(threads[i], NULL) != 0))
			return NULL;
		*p++ = thread_args[i].median;
		free(thread_args[i].parr);
	}
	return output;
}



