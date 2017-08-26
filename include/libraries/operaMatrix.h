#ifndef LIBOPERAMATRIX_H
#define LIBOPERAMATRIX_H
/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMatrix
 Version: 1.0
 Description: This C library implements matrix routines..
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

#include "libraries/operaLibCommon.h"
#include "libraries/operaVector.h"
#include <cmath>

/*! \brief Opera Matrix library. */
/*! \file operaMatrix.cpp */

/*! 
 * operaMatrix
 * \author Eder Martioli
 * \brief operaMatrix library.
 * \brief This library contains the basic routines for matrix tools.
 * \brief The main use of the CMatrix is to have matrices available in C
 * \brief whose bounds are known only at runtime.
 * \ingroup libraries
 */

template <typename T>
class Matrix {
private:
	std::vector<std::vector<T> > data;
	unsigned r;
	unsigned c;
public:
	Matrix(unsigned rows, unsigned cols) : data(rows, std::vector<T>(cols)), r(rows), c(cols) { }
	std::vector<T>& operator[](unsigned i) { return data[i]; }
	const std::vector<T>& operator[](unsigned i) const { return data[i]; }
	unsigned rows() const { return r; }
	unsigned cols() const { return c; }
};

template <typename T>
class Cube {
private:
	std::vector<Matrix<T> > data;
	unsigned s;
	unsigned r;
	unsigned c;
public:
	Cube(unsigned slices, unsigned rows, unsigned cols) : data(slices, Matrix<T>(rows, cols)), s(slices), r(rows), c(cols) { }
	Matrix<T>& operator[](unsigned i) { return data[i]; }
	const Matrix<T>& operator[](unsigned i) const { return data[i]; }
	unsigned slices() const { return s; }
	unsigned rows() const { return r; }
	unsigned cols() const { return c; }
	void clear() { data.clear(); }
};

typedef Matrix<double> DMatrix;
typedef Cube<double> DCube;

/**********************************************************************************/
/**********************************************************************************/
/**** NOTE WELL:                                                               ****/
/**** Theses functions return pointers to CMatrices created inside the         ****/
/**** functions. The caller is responsible for calling deleteCMatrix() when    ****/
/**** finished with the reutrned pointer.                                      ****/
/**********************************************************************************/
/**********************************************************************************/
	
	
/*! 
 * float MatrixDeterminant(CMatrix inputmatrix)
 * \brief  This function calculates the value of determinant of a square matrix 
 * \brief  it applies the recursive definition of determinate using expansion by minors.
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the determinant  
 * \note Source: http://paulbourke.net/miscellaneous/determinant/
 */
float MatrixDeterminant(CMatrix inputmatrix);
/*! 
 * CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the transpose of a given matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for transpose matrix
 */ 
CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix);
/*! 
 * float MatrixTrace(CMatrix inputmatrix)
 * \brief  This function calculates the trace of a square matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the trace  
 * \note Input matrix must be square 
 */
float MatrixTrace(CMatrix inputmatrix);
/*! 
 * CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the cofactor matrix
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for cofactor matrix
 * \note Input matrix must be square  
 */ 
CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixMultiplication(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue)
 * \brief  This function multiply a matrix by a constant 
 * \param  inputmatrix is a CMatrix type 
 * \param  float constantValue  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2)
 * \brief  This function calculates the sum between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2)
 * \brief  This function calculates the subtraction between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix);
/*! 
 * CMatrix MatrixInverse(CMatrix inputmatrix)
 * \brief  This function calculates the inverse matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for output inverse matrix
 */ 
CMatrix MatrixInverse(CMatrix inputmatrix, CMatrix outputmatrix);
/*! 
 * CMatrix RotationMatrix2D(float angleInDegrees) {
 * \brief  This function produces an operator matrix to rotate 
 * \brief  points in the xy-Cartesian plane counterclockwise through 
 * \brief  a given angle about the origin of the Cartesian coordinate system. 
 * \param  float angleInDegrees is the rotation angle in degrees. 
 * \return CMatrix (2X2) for rotation matrix
 */ 
CMatrix RotationMatrix2D(float angleInDegrees, CMatrix outputmatrix);
/*! 
 * void printMatrix(CMatrix inputmatrix)
 * \brief  This function prints the values of a matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return void
 */
void printMatrix(CMatrix matrix);

CMatrix MatrixIdentity(CMatrix InputsquareMatrix);
float* MatrixGetDiagonal(CMatrix inputmatrix, float *v);
CMatrix MatrixEigenValue(CMatrix inputmatrix, CMatrix eigenVectorMatrix, float *eigenValues);

void jacobi_eigenvalue(int n, double a[], int it_max, double v[], double d[], int *it_num, int *rot_num);
void r8mat_diag_get_vector ( int n, double a[], double v[] );
void r8mat_identity ( int n, double a[] );
double r8mat_is_eigen_right ( int n, int k, double a[], double x[], double lambda[] );
double r8mat_norm_fro ( int m, int n, double a[] );
     
#endif
