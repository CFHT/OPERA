/*******************************************************************
****                  MODULE FOR OPERA v1.0                     ****
********************************************************************
Module name: operaStatsLibTest
Version: 1.0
Description: Test the imageStats library functions.
Author(s): CFHT OPERA team
Affiliation: Canada France Hawaii Telescope 
Location: Hawaii USA
Date: Jan/2011
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaStats.h"


/*! \file operaStatsLibTest.c */

/*! 
 * operaStatsLibTest
 * \author Doug Teeple
 * \brief Test the imageStats library functions.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
		unsigned i;
		float xmean, xmedian, xsig, xmedsig, xmax, xmin;
		float x[8] = {10,5,8,7,10,11,3,9};
		float w[8] = {0.8,1.0,0.4,0.9,1.0,0.75,0.9,0.95};
	
		printf("Array:"); 		
		for(i=0;i<8;i++)
			printf(" %.0f[%.2f]",x[i],w[i]);
		printf("\n"); 	
		
		xmean = operaArrayMean(8,x);
		printf("Mean: %f\n", xmean);
		xmean = operaArrayWeightedMean(8,x,w);
		printf("Weighted mean: %f\n", xmean);	
		xmedian = operaArrayMedian(8,x);
		printf("Median: %f\n", xmedian);		
		xsig = operaArraySigma(8,x);
		printf("Sig: %f\n", xsig);	
		xsig = operaArrayWeightedSigma(8, x, w);
		printf("Weighted sig: %f\n", xsig);	
		xmedsig = operaArrayMedianSigma(8, x, xmedian);	
		printf("MedSig: %f\n", xmedsig);	
		xmax = operaArrayMaxValue(8, x);
		printf("Max: %f\n", xmax);
		xmin = operaArrayMinValue(8, x);
		printf("Min: %f\n", xmin);	
	
		srand(time(NULL));
	
		const float inarray[] = {0,0,1,2,1,4,5,4,5,3,3,3,4,5,6,7,8,10,-1,-1};
		unsigned int incount = 20;
		unsigned int outarray[3];
		unsigned int nbins = 3;
		unsigned int binsize = 3;
		float minvalue = 0;
		float maxvalue = 9;
		operaErrorCode err = operaHistogram(outarray, nbins, inarray, incount, binsize, minvalue, maxvalue);
		// print the results
		for (unsigned n = 0; n < 3; n++)
			printf("%d %d\n", outarray[n] , err);
		// should get 5, 9, 3	

	/*
		for(i=0;i<1000;i++) {
			printf("%u\t%f\n",i,operaUniformRand(10, 5, 8));
		}			
		
		for(i=0;i<1000;i++) {
			printf("%u\t%f\n",i,operaGaussRand(10, 3));		
		}		
*/	
  return EXIT_SUCCESS;
}  
