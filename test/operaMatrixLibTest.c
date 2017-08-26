/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaMatrixLibTest
 Version: 1.0
 Description: Test the operaMatrix library functions.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2012
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
#include "libraries/operaMatrix.h"
#include "libraries/operaLibCommon.h"

/*! \file operaMatrixLibTest.c */

/*! 
 * operaMatrixLibTest
 * \author Eder Martioli
 * \brief Test the operaMatrix library functions.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
	
	// Input 4X4 matrix to test operaMatrix routines.
	// Results have been checked with the oline tool:
	// http://www.bluebit.gr/matrix-calculator/
	
	unsigned NX = 4;
	unsigned NY = 4;
	
	CMatrix matrix = newCMatrix(NX,NY);
	CMatrix outputmatrix = newCMatrix(NX,NY);
	
	matrix[0][0] = 3;
	matrix[0][1] = 12;
	matrix[0][2] = 1;
	matrix[0][3] = 13;
	
	matrix[1][0] = 5;
	matrix[1][1] = 7;
	matrix[1][2] = -13;
	matrix[1][3] = 3;
	
	matrix[2][0] = -21;
	matrix[2][1] = 19;
	matrix[2][2] = 11;	
	matrix[2][3] = 5;		
	
	matrix[3][0] = 22;
	matrix[3][1] = 4;
	matrix[3][2] = 2;	
	matrix[3][3] = 0;		
	/*
	 3.000  12.000   1.000  13.000
	 5.000   7.000 -13.000   3.000
	 -21.000  19.000  11.000   5.000
	 22.000   4.000   2.000   0.000
	 */
	printf("The input matrix A is:\n\n");	
	printMatrix(matrix);
	
	printf("\nDeterminant of A: %f\n",MatrixDeterminant(matrix));
	
	CMatrix matrix_t = MatrixTranspose(matrix, outputmatrix);
	printf("\nTranspose of A:\n\n");
	printMatrix(matrix_t);
	
	printf("\nTrace of A: %f\n",MatrixTrace(matrix));	
	
	CMatrix matrix2 = newCMatrix(1,NX);
	CMatrix outputmatrix2 = newCMatrix(1,NY);
	matrix2[0][0] = 1;
	matrix2[1][0] = 2;
	matrix2[2][0] = -1;
	matrix2[3][0] = 0;
	
	printf("\nThe input matrix B is:\n\n");
	printMatrix(matrix2);	
	
	CMatrix productmatrix = MatrixMultiplication(matrix,matrix2,outputmatrix2);	
	
	printf("\nThe product (A x B) is:\n\n");
	printMatrix(productmatrix);	
	
	CMatrix inversematrix = MatrixInverse(matrix,outputmatrix);
	printf("\nInverse of A:\n\n");
	printMatrix(inversematrix);
	
	/*
	 * Now test the polynomial matrix
	 */
	printf("Testing PolynomialMatrix...\n");
	
	PolynomialMatrix polym = newPolynomialMatrix(20, 40);
	
	for (unsigned y=0; y<40; y++) {
		for (unsigned x=0; x<20; x++) {
			PolynomialCoeffs_t *pp = polym[y][x];
			pp->orderofPolynomial = 3;
			pp->p[0] = 1.0;
			pp->p[1] = 0.1;
			pp->p[2] = 0.0;
		}
	}
	
	for (unsigned y=0; y<40; y++) {
		for (unsigned x=0; x<20; x++) {
			PolynomialCoeffs_t *pp = polym[y][x];
			printf("x=%d y=%u order=%d p[0]=%f p[1]=%f p[2]=%f \n", x, y, pp->orderofPolynomial, pp->p[0], pp->p[1], pp->p[2]);
		}
	}
	
	/*
	 * Now test the Cube
	 */
	CCube cCube = newCCube(10, 10, 8);	// x, y, z
	
	printf("Testing Ccube...\n");
	
	for (unsigned z=0; z<8; z++) {
		for (unsigned y=0; y<10; y++) {
			for (unsigned x=0; x<10; x++) {
				cCube[z][y][x] = z * 1000.0 + y * 100.00 + x;
			}
		}
	}
	
	for (unsigned z=0; z<8; z++) {
		for (unsigned y=0; y<10; y++) {
			printf("z=%u y=%u ", z, y);
			for (unsigned x=0; x<10; x++) {
				printf("x=%u %.0f ", x, cCube[z][y][x]);
			}
			printf("\n");
		}
	}
	
	PolynomialMatrix polymatrix = newPolynomialMatrix(10, 10);
	// tests here...
	polymatrix = polymatrix;
	
	return EXIT_SUCCESS;
}  






