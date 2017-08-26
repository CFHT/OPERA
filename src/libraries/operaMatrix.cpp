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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaMatrix.h"
#include <numeric>
#include <functional>

/* 
 * operaMatrix
 * \author Eder Martioli
 * \brief operaMatrix library.
 * \details {This library contains the basic routines for matrix tools.}
 * \file operaMatrix.cpp
 * \ingroup libraries
 */

/**********************************************************************************/
/**********************************************************************************/
/**** NOTE WELL:                                                               ****/
/**** Theses functions return pointers to CMatrices created inside the         ****/
/**** functions. The caller is responsible for calling deleteCMatrix() when    ****/
/**** finished with the returned CMtarix.                                      ****/
/**********************************************************************************/
/**********************************************************************************/

/* 
 * void printMatrix(CMatrix inputmatrix)
 * \brief  This function prints the values of a matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return void
 */
void printMatrix(CMatrix matrix){
	for(unsigned j=0; j<getCMatrixRows(matrix); j++) {		
		for(unsigned i=0; i<getCMatrixCols(matrix); i++) {
			printf("%6.7f\t",matrix[j][i]); 
		}
		printf("\n");
	}	
}

/* 
 * float MatrixDeterminant(CMatrix inputmatrix)
 * \brief  This function calculates the value of determinant of a square matrix 
 * \brief  it applies the recursive definition of determinate using expansion by minors.
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the determinant  
 * \notes Source: http://paulbourke.net/miscellaneous/determinant/
 */
float MatrixDeterminant(CMatrix inputmatrix)
{
	
	float det;
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);
	
	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixDeterminant ", MatrixNotSquare);
		return FP_NAN;
	}
	
	unsigned n = NXPoints;

	unsigned i,j,j1,j2;
	
	det = 0;
	
	if (n == 1) {
		
		det = inputmatrix[0][0];
		
	} else if (n == 2)  {
		
		det = (inputmatrix[0][0] * inputmatrix[1][1]) - (inputmatrix[0][1] * inputmatrix[1][0]);		

	} else {
		
		det = 0;		
		
		for (j1 = 0 ; j1 < n ; j1++) {
			
			CMatrix minormatrix = newCMatrix((n-1),(n-1));
			
			for (i = 1 ; i < n ; i++) {
				j2 = 0 ;              
				for (j = 0 ; j < n ; j++) {
					if (j == j1) continue;
					minormatrix[j2][i-1] = inputmatrix[j][i];
					j2++;
				}
			}
			
			det += pow(-1.0,1.0 + j1 + 1.0) * inputmatrix[j1][0] * MatrixDeterminant(minormatrix);
			
			deleteCMatrix(minormatrix);
		}
	}
	
	return(det);
}

 /* 
 * CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the transpose of a given matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for transpose matrix
 */ 
CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	if(getCMatrixCols(inputmatrix) != getCMatrixRows(outputmatrix)) {
		operaPError("operaMatrix:MatrixTranspose ", MatrixInvalidDimensions);
		return NULL;
	}
	if(getCMatrixCols(outputmatrix) != getCMatrixRows(inputmatrix)) {
		operaPError("operaMatrix:MatrixTranspose ", MatrixInvalidDimensions);
		return NULL;
	}
	unsigned NXPoints = getCMatrixCols(outputmatrix);	
	unsigned NYPoints = getCMatrixRows(outputmatrix);

	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			outputmatrix[j][i] = inputmatrix[i][j];
		}
	}		
	
	return outputmatrix;
}

/* 
 * float MatrixTrace(CMatrix inputmatrix)
 * \brief  This function calculates the trace of a square matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the trace  
 * \notes Input matrix must be square 
 */
float MatrixTrace(CMatrix inputmatrix)
{
	float tracevalue;
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);
	
	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixTrace ", MatrixNotSquare);
		return FP_NAN;
	}
	
	tracevalue=0;
	
	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			if(i==j) {
				tracevalue += inputmatrix[j][i];
			}
		}
	}
	
	return(tracevalue);
}


/* 
 * CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the cofactor matrix
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for cofactor matrix
 * \notes Input matrix must be square  
 */ 

CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);	

	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixCofactor ", MatrixNotSquare);
		return NULL;
	}
	
	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			
			CMatrix minormatrix = newCMatrix(NXPoints-1,NYPoints-1);
			if (minormatrix) {
				/* Form the adjoint a_ij */
				unsigned i1 = 0;
				for (unsigned ii=0;ii<NXPoints;ii++) {
					if (ii == i)
						continue;
					unsigned j1 = 0;
					for (unsigned jj=0;jj<NYPoints;jj++) {
						if (jj == j)
							continue;
						minormatrix[j1][i1] = inputmatrix[jj][ii];
						j1++;
					}
					i1++;
				}
				
				/* Calculate the determinate */
				float det = MatrixDeterminant(minormatrix);
				
				/* Fill in the elements of the cofactor */
				outputmatrix[j][i] = pow(-1.0,i+j+2.0) * det;
				
				deleteCMatrix(minormatrix);					
			}
		}
	}
	return outputmatrix;
}


/* 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix) {
	return MatrixTranspose(MatrixCofactor(inputmatrix, outputmatrix), outputmatrix);
}

/* 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixMultiplication(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	

	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	unsigned NCOLS3 = getCMatrixCols(outputmatrix);	
	unsigned NROWS3 = getCMatrixRows(outputmatrix);		
	
	if(NCOLS1 != NROWS2) {
		operaPError("operaMatrix:MatrixMultiplication ", MatrixInvalidDimensions);
		return NULL;
	}	
	unsigned NCOMMON = NCOLS1;
	
	unsigned NCOLSRES = NCOLS2;	
	unsigned NROWSRES = NROWS1;	
	
	if(NCOLSRES != NCOLS3 || NROWSRES != NROWS3) {
		operaPError("operaMatrix:MatrixMultiplication ", MatrixInvalidDimensions);
		return NULL;
	}
    
    float sum = 0;
                                    
	for (unsigned j=0;j<NROWSRES;j++) {
		for (unsigned i=0;i<NCOLSRES;i++) {
			sum = 0;
			for (unsigned k=0;k<NCOMMON;k++) {
                sum += inputmatrix1[j][k]*inputmatrix2[k][i];
			}
            outputmatrix[j][i] = sum;
		}
	}
    
	return outputmatrix;
}

/* 
 * CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue, CMatrix outputmatrix)
 * \brief  This function multiply a matrix by a constant 
 * \param  inputmatrix is a CMatrix type 
 * \param  float constantValue  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue, CMatrix outputmatrix) {
	
	unsigned NCOLS = getCMatrixCols(inputmatrix);	
	unsigned NROWS = getCMatrixRows(inputmatrix);	
	
	for (unsigned j=0;j<NROWS;j++) {
		for (unsigned i=0;i<NCOLS;i++) {
			outputmatrix[j][i] = inputmatrix[j][i]*constantValue;
		}
	}				
	return outputmatrix;
}

/* 
 * CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix)
 * \brief  This function calculates the sum between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	
	
	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	if(NCOLS1 != NCOLS2 || NROWS1 != NROWS2) {
		operaPError("operaMatrix:MatrixAddition ", MatrixInvalidDimensions);
		return NULL;
	}	
	
	unsigned NCOLSRES = NCOLS1;	
	unsigned NROWSRES = NROWS1;	
	
	for (unsigned j=0;j<NROWSRES;j++) {
		for (unsigned i=0;i<NCOLSRES;i++) {
			outputmatrix[j][i] = inputmatrix1[j][i] + inputmatrix2[j][i];
		}
	}	
	return outputmatrix;
}

/* 
 * CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2)
 * \brief  This function calculates the subtraction between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	
	
	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	unsigned NCOLS3 = getCMatrixCols(outputmatrix);	
	unsigned NROWS3 = getCMatrixRows(outputmatrix);		
	
	if(NCOLS1 != NCOLS2 || NROWS1 != NROWS2 || NROWS1 != NROWS3 || NCOLS1 != NCOLS3) {
		operaPError("operaMatrix:MatrixSubtraction ", MatrixInvalidDimensions);
		return NULL;
	}	
	
	for (unsigned j=0;j<NROWS3;j++) {
		for (unsigned i=0;i<NCOLS3;i++) {
			outputmatrix[j][i] = inputmatrix1[j][i] - inputmatrix2[j][i];
		}
	}		
	return outputmatrix;
}

/* 
 * CMatrix MatrixInverse(CMatrix inputmatrix)
 * \brief  This function calculates the inverse matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for output inverse matrix
 */ 
CMatrix MatrixInverse(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	float det = MatrixDeterminant(inputmatrix);
	
	if(det == 0) {
		operaPError("operaMatrix:MatrixInverse ", MatrixZeroDeterminant);
		return NULL;
	}	
	
	return MatrixMultiplicationbyConstant(MatrixAdjoint(inputmatrix, outputmatrix), (1.0/det), outputmatrix);
}


/* 
 * CMatrix RotationMatrix2D(float angleInDegrees, CMatrix outputmatrix) {
 * \brief  This function produces an operator matrix to rotate 
 * \brief  points in the xy-Cartesian plane counterclockwise through 
 * \brief  a given angle about the origin of the Cartesian coordinate system. 
 * \param  float angleInDegrees is the rotation angle in degrees. 
 * \return CMatrix (2X2) for rotation matrix
 */ 
CMatrix RotationMatrix2D(float angleInDegrees, CMatrix outputmatrix) {
    
    unsigned NCOLS = 2;
    unsigned NROWS = 2;
    
	unsigned NCOLS2 = getCMatrixCols(outputmatrix);	
	unsigned NROWS2 = getCMatrixRows(outputmatrix);		

	if(NCOLS != NCOLS2 || NROWS != NROWS2) {
		operaPError("operaMatrix:MatrixSubtraction ", MatrixInvalidDimensions);
		return NULL;
	}	
    float angleInRadians = angleInDegrees*M_PI/180.0;
    
    outputmatrix[0][0] = cos(angleInRadians);
    outputmatrix[0][1] = -sin(angleInRadians);
    outputmatrix[1][0] = sin(angleInRadians);
    outputmatrix[1][1] = cos(angleInRadians); 
	
    return outputmatrix;
}

/*
 * CMatrix MatrixIdentity(CMatrix InputsquareMatrix) 
 * \brief  This function produces an identity matrix 
 */
CMatrix MatrixIdentity(CMatrix InputsquareMatrix) {
	unsigned NCOLS = getCMatrixCols(InputsquareMatrix);
	unsigned NROWS = getCMatrixRows(InputsquareMatrix);
    
    if(NCOLS != NROWS) {
		operaPError("operaMatrix:MatrixIdentity", MatrixInvalidDimensions);
		return NULL;
    }
    
    for (unsigned i=0;i<NROWS;i++) {
        for (unsigned j=0;j<NCOLS;j++) {
            if (i == j) {
                InputsquareMatrix[j][i] = 1.0;
            } else {
                InputsquareMatrix[j][i] = 0.0;
            }
		}
	}
    return InputsquareMatrix;
}

float *MatrixGetDiagonal (CMatrix inputmatrix, float *v) {
   
	unsigned NCOLS = getCMatrixCols(inputmatrix);
	unsigned NROWS = getCMatrixRows(inputmatrix);
    if(NCOLS != NROWS) {
		operaPError("operaMatrix:MatrixGetDiagonal", MatrixInvalidDimensions);
		return NULL;
    }
    unsigned np = NROWS;
    v = (float *) malloc(np*sizeof(float));
    
    unsigned count = 0;
    for (unsigned i=0;i<NROWS;i++) {
        for (unsigned j=0;j<NCOLS;j++) {
            if (i == j) {
                v[count++] = inputmatrix[j][i];
            }
		}
	}
    return v;
}


CMatrix MatrixEigenValue (CMatrix inputmatrix, CMatrix eigenVectorMatrix, float *eigenValues) {
    
	unsigned NCOLS = getCMatrixCols(inputmatrix);
	unsigned NROWS = getCMatrixRows(inputmatrix);
    if(NCOLS != NROWS || getCMatrixCols(eigenVectorMatrix) != NCOLS || getCMatrixRows(eigenVectorMatrix) != NROWS) {
		operaPError("operaMatrix:MatrixGetDiagonal", MatrixInvalidDimensions);
		return NULL;
    }
    unsigned np = NROWS;
    
    double *a = (double *) malloc(np*np*sizeof(double));
    
    for (unsigned i=0;i<NROWS;i++) {
        for (unsigned j=0;j<NCOLS;j++) {
            a[i*NROWS + j] = (double)inputmatrix[j][i];
		}
	}
    
    double *d = (double *) malloc(np*sizeof(double));
    double *v = (double *) malloc(np*np*sizeof(double));
    
    // DT: not used! double error_frobenius;
    int it_max = 100;
    int it_num;
    int n = np;
    int rot_num;
    
    jacobi_eigenvalue (n, a, it_max, v, d, &it_num, &rot_num);

    for (unsigned i=0;i<NROWS;i++) {
        for (unsigned j=0;j<NCOLS;j++) {
            eigenVectorMatrix[j][i] = (float)v[i*NROWS + j];
		}
        eigenValues[i] = (float)d[i];
	}

    free(a);
    free(v);
    free(d);
    
    return eigenVectorMatrix;
}


/*
printMatrix
MatrixDeterminant
MatrixTranspose 	
MatrixTrace
MatrixCofactor
MatrixAdjoint
MatrixMultiplication
MatrixMultiplicationbyConstant
MatrixAddition
MatrixSubtraction
MatrixInverse 	
Eigenvalues
Eigenvectors
 
CholeskyDecomposition
Moore-PenroseInverse
QRDecomposition

LUDecomposition 	
SingularValueDecomposition 	
SystemsofLinearEquations
MatrixRank 	
LQDecomposition
*/


// Source: http://people.sc.fsu.edu/~jburkardt/cpp_src/jacobi_eigenvalue/jacobi_eigenvalue.html

//****************************************************************************80

void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], double d[], int *it_num, int *rot_num )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
//
//  Discussion:
//
//    This function computes the eigenvalues and eigenvectors of a
//    real symmetric matrix, using Rutishauser's modfications of the classical
//    Jacobi rotation method with threshold pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2013
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix, which must be square, real,
//    and symmetric.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, double V[N*N], the matrix of eigenvectors.
//
//    Output, double D[N], the eigenvalues, in descending order.
//
//    Output, int &IT_NUM, the total number of iterations.
//
//    Output, int &ROT_NUM, the total number of rotations.
//
{
    double *bw;
    double c;
    double g;
    double gapq;
    double h;
    int i;
    int j;
    int k;
    int l;
    int m;
    int p;
    int q;
    double s;
    double t;
    double tau;
    double term;
    double termp;
    double termq;
    double theta;
    double thresh;
    double w;
    double *zw;
    
    r8mat_identity ( n, v );
    
    r8mat_diag_get_vector ( n, a, d );
    
    bw = (double *) malloc(n*sizeof(double));
    zw = (double *) malloc(n*sizeof(double));

    for ( i = 0; i < n; i++ )
    {
        bw[i] = d[i];
        zw[i] = 0.0;
    }
    it_num = 0;
    rot_num = 0;
    
    while ( *it_num < it_max )	// DT May 8 2014 was comparing a pointer to an int
    {
        it_num = it_num + 1;
        //
        //  The convergence threshold is based on the size of the elements in
        //  the strict upper triangle of the matrix.
        //
        thresh = 0.0;
        for ( j = 0; j < n; j++ )
        {
            for ( i = 0; i < j; i++ )
            {
                thresh = thresh + a[i+j*n] * a[i+j*n];
            }
        }
        
        thresh = sqrt ( thresh ) / ( double ) ( 4 * n );
        
        if ( thresh == 0.0 )
        {
            break;
        }
        
        for ( p = 0; p < n; p++ )
        {
            for ( q = p + 1; q < n; q++ )
            {
                gapq = 10.0 * fabs ( a[p+q*n] );
                termp = gapq + fabs ( d[p] );
                termq = gapq + fabs ( d[q] );
                //
                //  Annihilate tiny offdiagonal elements.
                //
                if ( 4 < *it_num &&	// DT May 8 2014 was comparing a pointer to an int
                    termp == fabs ( d[p] ) &&
                    termq == fabs ( d[q] ) )
                {
                    a[p+q*n] = 0.0;
                }
                //
                //  Otherwise, apply a rotation.
                //
                else if ( thresh <= fabs ( a[p+q*n] ) )
                {
                    h = d[q] - d[p];
                    term = fabs ( h ) + gapq;
                    
                    if ( term == fabs ( h ) )
                    {
                        t = a[p+q*n] / h;
                    }
                    else
                    {
                        theta = 0.5 * h / a[p+q*n];
                        t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
                        if ( theta < 0.0 )
                        {
                            t = - t;
                        }
                    }
                    c = 1.0 / sqrt ( 1.0 + t * t );
                    s = t * c;
                    tau = s / ( 1.0 + c );
                    h = t * a[p+q*n];
                    //
                    //  Accumulate corrections to diagonal elements.
                    //
                    zw[p] = zw[p] - h;
                    zw[q] = zw[q] + h;
                    d[p] = d[p] - h;
                    d[q] = d[q] + h;
                    
                    a[p+q*n] = 0.0;
                    //
                    //  Rotate, using information from the upper triangle of A only.
                    //
                    for ( j = 0; j < p; j++ )
                    {
                        g = a[j+p*n];
                        h = a[j+q*n];
                        a[j+p*n] = g - s * ( h + g * tau );
                        a[j+q*n] = h + s * ( g - h * tau );
                    }
                    
                    for ( j = p + 1; j < q; j++ )
                    {
                        g = a[p+j*n];
                        h = a[j+q*n];
                        a[p+j*n] = g - s * ( h + g * tau );
                        a[j+q*n] = h + s * ( g - h * tau );
                    }
                    
                    for ( j = q + 1; j < n; j++ )
                    {
                        g = a[p+j*n];
                        h = a[q+j*n];
                        a[p+j*n] = g - s * ( h + g * tau );
                        a[q+j*n] = h + s * ( g - h * tau );
                    }
                    //
                    //  Accumulate information in the eigenvector matrix.
                    //
                    for ( j = 0; j < n; j++ )
                    {
                        g = v[j+p*n];
                        h = v[j+q*n];
                        v[j+p*n] = g - s * ( h + g * tau );
                        v[j+q*n] = h + s * ( g - h * tau );
                    }
                    rot_num = rot_num + 1;
                }
            }
        }
        
        for ( i = 0; i < n; i++ )
        {
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }
    //
    //  Restore upper triangle of input matrix.
    //
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < j; i++ )
        {
            a[i+j*n] = a[j+i*n];
        }
    }
    //
    //  Ascending sort the eigenvalues and eigenvectors.
    //
    for ( k = 0; k < n - 1; k++ )
    {
        m = k;
        for ( l = k + 1; l < n; l++ )
        {
            if ( d[l] < d[m] )
            {
                m = l;
            }
        }
        
        if ( m != k )
        {
            t    = d[m];
            d[m] = d[k];
            d[k] = t;
            for ( i = 0; i < n; i++ )
            {
                w        = v[i+m*n];
                v[i+m*n] = v[i+k*n];
                v[i+k*n] = w;
            }
        }
    }
    
    free(bw);
    free(zw);
    
    return;
}
//****************************************************************************80

void r8mat_diag_get_vector ( int n, double a[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N*N], the N by N matrix.
//
//    Output, double V[N], the diagonal entries
//    of the matrix.
//
{
    int i;
    
    for ( i = 0; i < n; i++ )
    {
        v[i] = a[i+i*n];
    }
    
    return;
}
//****************************************************************************80

void r8mat_identity ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IDENTITY sets the square matrix A to the identity.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double A[N*N], the N by N identity matrix.
//
{
    int i;
    int j;
    int k;
    
    k = 0;
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            if ( i == j )
            {
                a[k] = 1.0;
            }
            else
            {
                a[k] = 0.0;
            }
            k = k + 1;
        }
    }
    
    return;
}

//****************************************************************************80

double r8mat_is_eigen_right ( int n, int k, double a[], double x[],
                             double lambda[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.
//
//  Discussion:
//
//    An R8MAT is a matrix of doubles.
//
//    This routine computes the Frobenius norm of
//
//      A * X - X * LAMBDA
//
//    where
//
//      A is an N by N matrix,
//      X is an N by K matrix (each of K columns is an eigenvector)
//      LAMBDA is a K by K diagonal matrix of eigenvalues.
//
//    This routine assumes that A, X and LAMBDA are all real.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int K, the number of eigenvectors.
//    K is usually 1 or N.
//
//    Input, double A[N*N], the matrix.
//
//    Input, double X[N*K], the K eigenvectors.
//
//    Input, double LAMBDA[K], the K eigenvalues.
//
//    Output, double R8MAT_IS_EIGEN_RIGHT, the Frobenius norm
//    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
//    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
//
{
    double *c;
    double error_frobenius;
    int i;
    int j;
    int l;
    
    c = (double *) malloc(n*k*sizeof(double));
    
    for ( j = 0; j < k; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            c[i+j*n] = 0.0;
            for ( l = 0; l < n; l++ )
            {
                c[i+j*n] = c[i+j*n] + a[i+l*n] * x[l+j*n];
            }
        }
    }
    
    for ( j = 0; j < k; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            c[i+j*n] = c[i+j*n] - lambda[j] * x[i+j*n];
        }
    }
    
    error_frobenius = r8mat_norm_fro ( n, k, c );
    
    free(c);
    
    return error_frobenius;
}
//****************************************************************************80

double r8mat_norm_fro ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
//
{
    int i;
    int j;
    double value;
    
    value = 0.0;
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            value = value + pow ( a[i+j*m], 2 );
        }
    }
    value = sqrt ( value );
    
    return value;
}
//****************************************************************************80
