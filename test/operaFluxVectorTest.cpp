/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFluxvectorTest
 Version: 1.0
 Description: Perform various tests on the operaFITSImage class.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

#include <stdio.h>		// for printf

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFluxVector.h"

/*! \file operaFluxvectorTest.cpp */

using namespace std;

/*! 
 * operaFluxvectorTest
 * \author Doug Teeple
 * \ingroup test
 */
int main()
{
	const unsigned length = 5;
	
	double i1EFluxes[] = {2000.0, 1450.0, 123.0, 4324.9, 1235.0};
	double i1EVariances[] = {0.25, 11.2, 0.2, 9.8, 12.2};
	double i1AFluxes[] = {2500.0, 1050.0, 1230.0, 434.9, 235.0};
	double i1AVariances[] = {0.25, 1.2, 0.2, 1.8, 2.2};
	
	double i2EFluxes[] = {2400.0, 1650.0, 103.0, 4424.9, 1035.0};
	double i2EVariances[] = {0.25, 11.2, 0.2, 9.8, 12.2};
	double i2AFluxes[] = {2200.0, 1250.0, 1030.0, 494.9, 225.0};
	double i2AVariances[] = {0.25, 1.2, 0.2, 1.8, 2.2};
	
	double i3EFluxes[] = {2200.0, 1430.0, 163.0, 4524.9, 1295.0};
	double i3EVariances[] = {0.25, 11.2, 0.2, 9.8, 12.2};
	double i3AFluxes[] = {2560.0, 1050.0, 1230.0, 434.9, 235.0};
	double i3AVariances[] = {0.25, 1.2, 0.2, 1.8, 2.2};
	
	double i4EFluxes[] = {2300.0, 1420.0, 173.0, 4624.9, 1265.0};
	double i4EVariances[] = {0.25, 11.2, 0.2, 9.8, 12.2};
	double i4AFluxes[] = {2520.0, 1550.0, 1238.0, 404.9, 205.0};
	double i4AVariances[] = {0.25, 1.2, 0.2, 1.8, 2.2};
	
	/*
	 * Populate the 
	 vectors with the E/A data
	 */
	operaFluxVector i1E(i1EFluxes, i1EVariances, length, ToZero);
	operaFluxVector i1A(i1AFluxes, i1AVariances, length, ToZero);
	operaFluxVector i2E(i2EFluxes, i2EVariances, length, ToZero);
	operaFluxVector i2A(i2AFluxes, i2AVariances, length, ToZero);
	operaFluxVector i3E(i3EFluxes, i3EVariances, length, ToZero);
	operaFluxVector i3A(i3AFluxes, i3AVariances, length, ToZero);
	operaFluxVector i4E(i4EFluxes, i4EVariances, length, ToZero);
	operaFluxVector i4A(i4AFluxes, i4AVariances, length, ToZero);
	
	/*
	 * Show the difference between content copy and address copy
	 */
	operaFluxVector* ptrTest1 = new operaFluxVector(length);
	operaFluxVector* ptrTest2 = new operaFluxVector(length);
	
	*ptrTest1 = i1E;		// copy contents
	*ptrTest2 = *ptrTest1;	// copy contents
	
	ptrTest2 = ptrTest1;	// copy address

	/*
	 * Populate the rn vectors with the E data
	 */
	operaFluxVector r1(length, ToOne);
	operaFluxVector r2(length, ToOne);
	operaFluxVector r3(length, ToOne);
	operaFluxVector r4(length, ToOne);
	
	/* 
	 * STEP 1 - calculate ratio of beams for each exposure
	 * r1 = i1E / i1A
	 * r2 = i2E / i2A
	 * r3 = i3E / i3A
	 * r4 = i4E / i4A
	 */
	r1 = i1E / i1A;
	r2 = i2E / i2A;
	r3 = i3E / i3A;
	r4 = i4E / i4A;
	/*	
	 * STEP 2 - calculate the quantity R (Eq #2 on page 663 of paper)
	 * 	R^4 = (r1 * r4) / (r2 * r3)
	 * 	R = the 4th root of R^4
	 */	
	operaFluxVector R(length, ToZero);
	R = Sqrt(Sqrt((r1 * r4) / (r2 * r3)));
	/*	
	 * STEP 3 - calculate the Stokes parameter (Eq #1 on page 663 of paper)
	 * 	aka "the polarization" P
	 * 	P/I = (R-1) / (R+1)
	 */
	operaFluxVector PoverI(length);
	PoverI = (R-1.0) / (R+1.0);
	/*
	 * STEP 4 - calculate the first NULL spectrum (Eq #3 on page 663 of paper)
	 * 	R2^4 = (r1 * r2) / (r3 * r4)
	 * 	R2 = the 4th root of R^4
	 * 	and repeat Step 3 to get "P/I" which is then the NULL spectrum.
	 */
	operaFluxVector R2(length, ToZero);
	R2 = Sqrt(Sqrt((r1 * r2) / (r3 * r4)));
	
	operaFluxVector PoverI2(length);
	PoverI2 = (R2-1.0) / (R2+1.0);
    
    /*
	 * CollapseInfinitiesDivisionToOne test
	 */
    
    operaFluxVector Infinity1(length, ToOne);
    operaFluxVector Infinity2(length, ToOne);
    operaFluxVector CollapseInfinitiesDivisionToOneResult(length);
    
    Infinity1 = (double)INFINITY;
    Infinity2 = (double)INFINITY;
    
    // deprecated CollapseInfinitiesDivisionToOneResult = collapseInfinitiesDivisionToOne(&Infinity1, &Infinity2);
	
	CollapseInfinitiesDivisionToOneResult = Infinity1 / Infinity2;
	printf("r1=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", r1[i]->first, r1[i]->second);
	}
	printf("\n");
	
	printf("r2=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", r2[i]->first, r2[i]->second);
	}
	printf("\n");
	
	printf("r3=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", r3[i]->first, r3[i]->second);
	}
	printf("\n");
	
	printf("r4=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", r4[i]->first, r4[i]->second);
	}
	printf("\n");
	
	printf("R=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", R[i]->first, R[i]->second);
	}
	printf("\n");
	
	printf("PoverI=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", PoverI[i]->first, PoverI[i]->second);
	}
	printf("\n");
	
	printf("PoverI2=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", PoverI2[i]->first, PoverI2[i]->second);
	}
	printf("\n");
	
	printf("R errors=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f, ", R.geterror(i));
	}
	printf("\n");
    
    printf("CollapseInfinitiesDivisionToOneResult=");
	for (unsigned i=0; i<length; i++) {
		printf("%.2f %.2f, ", CollapseInfinitiesDivisionToOneResult[i]->first, CollapseInfinitiesDivisionToOneResult[i]->second);
	}
	printf("\n");

	return 0;
}
