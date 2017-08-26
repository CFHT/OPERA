/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaMuellerMatrix
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Tele]scope
 
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

#ifndef OPERAMUELLERMATRIX_H
#define OPERAMUELLERMATRIX_H

#include <stdio.h>      // for printMuellerMatrix, printVarianceMatrix

#include "operaError.h"
#include "libraries/operaFluxVector.h"			// for FluxVector
#include "libraries/operaStokesVector.h"		// for StokesVector
#include "libraries/operaLibCommon.h"			// for doubleValue

/*!
 * \file operaMuellerMatrix.h
 * \brief This file holds the declaration of the class operaMuellerMatrix.
 * \ingroup libraries
 */

using namespace std;

/*!
 * \author Andre Venne
 * \brief This class encapsulates the Mueller matrix.
 * \sa class operaStokesVector, class operaPolarimetry
 * 
 * This class holds all the parameters and their variance to create a Mueller matrix. It also holds the Mueller matrix itself and its variance matrix.
 * The variances are propagated as followed :
 * 
 * F = F(a,b)
 * DF = Pow(dF/da,2) * Da + Pow(dF/db,2) *Db
 * 
 * where DF is the resulting variance, Da and Db are the variance of the variables a and b, dF/da and dF/db are the partial derivatives of F.
 * The variables are supposed uncorrelated.
 */
class operaMuellerMatrix {
	
private:
    bool istemp;
    
	double p;           // amplitude attenuation coefficient 0 <= p <= 1   |   px^2 + py^2 = p^2
    double alpha;       // trigonometric angle 0 <= alpha <= 90 deg  |   px = p cos(alpha), py = p sin(alpha)
    double phi;         // phase shift between orthogonal components
    double theta;       // rotation angle of components
    
    double pVariance;
    double alphaVariance;
    double phiVariance;
    double thetaVariance;

    double muellerMatrix[4][4];
    double varianceMatrix[4][4];
	
public:
	/*
	 * Constructors / Destructors
	 */
    
    /*!
     * \brief Basic operaMuellerMatrix constructor.
     * \param Istemp An optional bool defaults to false
     * \return void
     */
    operaMuellerMatrix(bool Istemp=false);
    
    /*!
     * \brief Basic operaMuellerMatrix destructor.
     * \return void
     */
    ~operaMuellerMatrix();
    
    /*
     * Creators
     */
    
    /*!
     * \brief Creates a rotated polarizer.
     * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotated polarizer.
     * \param P A double value, 0 <= P <= 1
     * \param PVariance A double value
     * \param Alpha A double value, 0 <= Alpha <= 90 deg
     * \param AlphaVariance A double value
     * \param Theta A double value
     * \param ThetaVariance A double value
     * \return void
     */
	void createRotatedPolarizer(double P, double PVariance, double Alpha, double AlphaVariance, double Theta, double ThetaVariance);                // for a Rotated Polarizer
    
    /*!
     * \brief Creates a rotated retarder.
     * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotated retarder.
     * \param Phi A double value
     * \param PhiVariance A double value
     * \param Theta A double value
     * \param ThetaVariance A double value
     * \return void
     */
    void createRotatedRetarder(double Phi, double PhiVariance, double Theta, double ThetaVariance);                                                 // for a Rotated Retarder
    
    /*!
     * \brief Creates a rotator.
     * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotator.
     * \param Theta A double value
     * \param ThetaVariance A double value
     * \return void
     */
    void createRotator(double Theta, double ThetaVariance);                                                                                         // for a Rotator
    
    /*!
     * \brief Creates an attenuator.
     * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents an attenuator.
     * \param P A double value, 0 <= P <= 1
     * \param PVariance A double value
     * \return void
     */
    void createAttenuator(double P, double PVariance);                                                                                              // for an Attenuator
	
    /*
	 * Getters/Setters
	 */
    
    /*!
     * \brief Sets the amplitude of attenuation.
     * \details A function that sets the amplitude of attenuation, which is defined as: px^2 + py^2 = p^2.
     * \param P A double value, 0 <= P <= 1
     * \return void
     */
    void setP(double P);
    
    /*!
     * \brief Gets the amplitude of attenuation.
     * \details A function that gets the amplitude of attenuation, which is defined as: px^2 + py^2 = p^2.
     * \return A double value
     */
    double getP(void);
    
    /*!
     * \brief Sets the variance of the amplitude of attenuation.
     * \details A function that sets the variance of the amplitude of attenuation.
     * \param PVariance A double value
     * \return void
     */
    void setPVariance(double PVariance);
    
    /*!
     * \brief Gets the variance of the amplitude of attenuation.
     * \details A function that gets the variance of the amplitude of attenuation.
     * \return A double value
     */
    double getPVariance(void);
    
    /*!
     * \brief Sets the trigonometric angle.
     * \details A function that sets the trigonometric angle, which is defined as: px = p cos(alpha), py = p sin(alpha).
     * \param Alpha A double value, 0 <= Alpha <= 90 deg
     * \return void
     */
    void setAlpha(double Alpha);
    
    /*!
     * \brief Gets the trigonometric angle.
     * \details A function that gets the trigonometric angle, which is defined as: px = p cos(alpha), py = p sin(alpha).
     * \return A double value
     */
    double getAlpha(void);
    
    /*!
     * \brief Sets the variance of the trigonometric angle.
     * \details A function that sets the variance of the trigonometric angle.
     * \param AlphaVariance A double value
     * \return void
     */
    void setAlphaVariance(double AlphaVariance);
    
    /*!
     * \brief Gets the variance of the trigonometric angle.
     * \details A function that gets the variance of the trigonometric angle.
     * \return A double value
     */
    double getAlphaVariance(void);
    
    /*!
     * \brief Sets the phase shift between orthogonal components.
     * \details A function that sets the phase shift between orthogonal components px and py.
     * \param Phi A double value
     * \return void
     */
    void setPhi(double Phi);
    
    /*!
     * \brief Gets the phase shift between orthogonal components.
     * \details A function that gets the phase shift between orthogonal components px and py.
     * \return A double value
     */
    double getPhi(void);
    
    /*!
     * \brief Sets the variance of the phase shift between orthogonal components.
     * \details A function that sets the variance of the phase shift between orthogonal components px and py.
     * \param PhiVariance A double value
     * \return void
     */
    void setPhiVariance(double PhiVariance);
    
    /*!
     * \brief Gets the variance of the phase shift between orthogonal components.
     * \details A function that gets the variance of the phase shift between orthogonal components px and py.
     * \return A double value
     */
    double getPhiVariance(void);
    
    /*!
     * \brief Sets the rotation angle of components.
     * \details A function that sets the rotation angle of components px and py.
     * \param Theta A double value
     * \return void
     */
    void setTheta(double Theta);
    
    /*!
     * \brief Gets the rotation angle of components.
     * \details A function that gets the rotation angle of components px and py.
     * \return A double value
     */
    double getTheta(void);
    
    /*!
     * \brief Sets the variance of the rotation angle of components.
     * \details A function that sets the variance of the rotation angle of components px and py.
     * \param ThetaVariance A double value
     * \return void
     */
    void setThetaVariance(double ThetaVariance);
    
    /*!
     * \brief Gets the variance of the rotation angle of components.
     * \details A function that gets the variance of the rotation angle of components px and py.
     * \return A double value
     */
    double getThetaVariance(void);
    
    /*!
     * \brief Sets the Mueller matrix.
     * \details A function that sets the Mueller matrix from a 4x4 matrix.
     * \param MuellerMatrix A 4x4 matrix of double values
     * \return void
     */
    void setMuellerMatrix(double MuellerMatrix[4][4]);
    
    /*!
     * \brief Gets a Mueller matrix element.
     * \details A function that gets a Mueller matrix element at a specified row and column.
     * \param indexX An unsigned index to the row position
     * \param indexY An unsigned index to the column position
     * \return A double value
     */
    double getMuellerMatrixElement(unsigned indexX, unsigned indexY);
    
    /*!
     * \brief Sets the variance matrix.
     * \details A function that sets the variance matrix from a 4x4 matrix.
     * \param VarianceMatrix A 4x4 matrix of double values
     * \return void
     */
    void setVarianceMatrix(double VarianceMatrix[4][4]);
    
    /*!
     * \brief Gets a variance matrix element.
     * \details A function that gets a variance matrix element at a specified row and column.
     * \param indexX An unsigned index to the row position
     * \param indexY An unsigned index to the column position
     * \return A double value
     */
    double getVarianceMatrixElement(unsigned indexX, unsigned indexY);
    
    /*
	 * Matrix Operations
	 */
    
    /*!
     * \brief Prints the Mueller matrix.
     * \details A function that prints the Mueller matrix to the terminal in a 4x4 matrix representation.
     * \return void
     */
    void printMuellerMatrix(void);
    
    /*!
     * \brief Prints the variance matrix.
     * \details A function that prints the variance matrix to the terminal in a 4x4 matrix representation.
     * \return void
     */
    void printVarianceMatrix(void);
    
    /*!
     * \brief Calculates the determinant of the Mueller matrix.
     * \details A function that calculates the determinant of the Mueller matrix. It also calculates the variance of the determinant using the variance matrix.
     * \return A doubleValue structure.
     */
    doubleValue matrixDeterminant(void);
    
    /*!
     * \brief Calculates the cofactor matrix of the Mueller matrix.
     * \details A function that calculates the cofactor matrix of the Mueller matrix. It also calculates the variance of the cofactor matrix elements using the variance matrix.
     * \param Istemp Optional bool defaults to false
     * \return An operaMuellerMatrix address.
     */
    operaMuellerMatrix& matrixCofactor(bool Istemp=false);
    
    /*!
     * \brief Calculates the transpose matrix of the Mueller matrix.
     * \details A function that calculates the transpose matrix of the Mueller matrix. It also calculates the variance of the transpose matrix elements using the variance matrix.
     * \param Istemp Optional bool defaults to false
     * \return An operaMuellerMatrix address.
     */
    operaMuellerMatrix& matrixTranspose(bool Istemp=false);
    
    /*!
     * \brief Calculates the adjoint matrix of the Mueller matrix.
     * \details A function that calculates the adjoint matrix of the Mueller matrix. It also calculates the variance of the adjoint matrix elements using the variance matrix.
     * \details It calls matrixCofactor and matrixTranspose.
     * \param Istemp Optional bool defaults to false
     * \return An operaMuellerMatrix address.
     */
    operaMuellerMatrix& matrixAdjoint(bool Istemp=false);
    
    /*!
     * \brief Calculates the inverse matrix of the Mueller matrix.
     * \details A function that calculates the inverse matrix of the Mueller matrix. It also calculates the variance of the inverse matrix elements using the variance matrix.
     * \details It calls matrixDeterminant and matrixAdjoint.
     * \param Istemp Optional bool defaults to false
     * \return An operaMuellerMatrix address.
     */
    operaMuellerMatrix& matrixInverse(bool Istemp=false);
    
    /*
	 * Operators
	 */
    
    /*!
	 * \brief Assignment operator.
     * \details The operator copies the amplitude attenuation coefficient p, the trigonometric angle alpha, the phase shift phi, the rotation angle theta, the Mueller matrix and their variance from the right side of the operator to the left side.
	 * \param b An operaMuellerMatrix pointer
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator=(operaMuellerMatrix* b) {
		this->p = b->p;
        this->alpha = b->alpha;
        this->phi = b->phi;
        this->theta = b->theta;
        this->pVariance = b->pVariance;
        this->alphaVariance = b->alphaVariance;
        this->phiVariance = b->phiVariance;
        this->thetaVariance = b->thetaVariance;
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				this->muellerMatrix[row][column] = b->muellerMatrix[row][column];
                this->varianceMatrix[row][column] = b->varianceMatrix[row][column];
			}
		}
        if (b->istemp)
            delete b;
		return *this;
	};
    
	/*!
	 * \brief Assignment operator.
     * \details The operator copies the amplitude attenuation coefficient p, the trigonometric angle alpha, the phase shift phi, the rotation angle theta, the Mueller matrix and their variance from the right side of the operator to the left side.
	 * \param b An operaMuellerMatrix address
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator=(operaMuellerMatrix& b) {
		this->p = b.p;
        this->alpha = b.alpha;
        this->phi = b.phi;
        this->theta = b.theta;
        this->pVariance = b.pVariance;
        this->alphaVariance = b.alphaVariance;
        this->phiVariance = b.phiVariance;
        this->thetaVariance = b.thetaVariance;
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				this->muellerMatrix[row][column] = b.muellerMatrix[row][column];
                this->varianceMatrix[row][column] = b.varianceMatrix[row][column];
			}
		}
        if (b.istemp)
            delete &b;
		return *this;
	};
    
    /*!
	 * \brief Addition operator.
     * \details The operator adds the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaMuellerMatrix pointer
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c + operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator+(operaMuellerMatrix* b) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				t->muellerMatrix[row][column] = this->muellerMatrix[row][column] + b->muellerMatrix[row][column];
                t->muellerMatrix[row][column] = this->varianceMatrix[row][column] + b->varianceMatrix[row][column];
			}
		}
        if (b->istemp)
            delete b;
		return *t;
	};
    
    /*!
	 * \brief Addition operator.
     * \details The operator adds the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaMuellerMatrix address
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c + operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator+(operaMuellerMatrix& b) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				t->muellerMatrix[row][column] = this->muellerMatrix[row][column] + b.muellerMatrix[row][column];
                t->muellerMatrix[row][column] = this->varianceMatrix[row][column] + b.varianceMatrix[row][column];
			}
		}
        if (b.istemp)
            delete &b;
		return *t;
	};
    
    /*!
	 * \brief Subtraction operator.
     * \details The operator subtracts the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaMuellerMatrix pointer
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c - operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator-(operaMuellerMatrix* b) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				t->muellerMatrix[row][column] = this->muellerMatrix[row][column] - b->muellerMatrix[row][column];
                t->muellerMatrix[row][column] = this->varianceMatrix[row][column] + b->varianceMatrix[row][column];
			}
		}
        if (b->istemp)
            delete b;
		return *t;
	};
    
    /*!
	 * \brief Subtraction operator.
     * \details The operator subtracts the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaMuellerMatrix address
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c - operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator-(operaMuellerMatrix& b) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
				t->muellerMatrix[row][column] = this->muellerMatrix[row][column] - b.muellerMatrix[row][column];
                t->muellerMatrix[row][column] = this->varianceMatrix[row][column] + b.varianceMatrix[row][column];
			}
		}
        if (b.istemp)
            delete &b;
		return *t;
	};
    
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the matrix on the left side of the operator with the matrix on the right side.
	 * \param b An operaMuellerMatrix pointer
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c * operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator*(operaMuellerMatrix* b) {
        operaMuellerMatrix *t = new operaMuellerMatrix(true);
        for (unsigned row = 0 ; row < 4 ; row++) {
            for (unsigned tcolumn = 0 ; tcolumn < 4 ; tcolumn++) {
                for (unsigned column = 0 ; column < 4 ; column++) {
                    t->muellerMatrix[row][tcolumn] += this->muellerMatrix[row][column] * b->muellerMatrix[column][tcolumn];
                    t->varianceMatrix[row][tcolumn] += pow(this->muellerMatrix[row][column],2) * b->varianceMatrix[column][tcolumn] + pow(b->muellerMatrix[column][tcolumn],2) * this->varianceMatrix[row][column];
                }
            }
		}
        if (b->istemp)
            delete b;
        if (this->istemp) {
            *this = *t;
            return *this;
        }
        else
            return *t;
	};
    
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the matrix on the left side of the operator with the matrix on the right side.
	 * \param b An operaMuellerMatrix address
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c * operaMuellerMatrix b;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator*(operaMuellerMatrix& b) {
        operaMuellerMatrix *t = new operaMuellerMatrix(true);
        for (unsigned row = 0 ; row < 4 ; row++) {
            for (unsigned tcolumn = 0 ; tcolumn < 4 ; tcolumn++) {
                for (unsigned column = 0 ; column < 4 ; column++) {
                    t->muellerMatrix[row][tcolumn] += this->muellerMatrix[row][column] * b.muellerMatrix[column][tcolumn];
                    t->varianceMatrix[row][tcolumn] += pow(this->muellerMatrix[row][column],2) * b.varianceMatrix[column][tcolumn] + pow(b.muellerMatrix[column][tcolumn],2) * this->varianceMatrix[row][column];
                }
            }
		}
        if (b.istemp)
            delete &b;
        if (this->istemp) {
            *this = *t;
            return *this;
        }
        else
            return *t;
	};
    
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the matrix on the left side of the operator with the Stokes vector on the right side.
	 * \param b An operaStokesVector pointer
	 * \note Usage: operaStokesVector a = operaMuellerMatrix c * operaStokesVector b;
	 * \return An operaStokesVector address
	 */
	operaStokesVector operator*(operaStokesVector* b) {
        operaStokesVector t(b->length);
		for (unsigned row = 0 ; row < 4 ; row++) {
			t.stokesVector[row] = 0.0;	
		}
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
                double *tstokes = t.stokesVector[row].getfluxpointer();
                double *tvariances = t.stokesVector[row].getvariancepointer();
                double *stokes = b->stokesVector[column].getfluxpointer();
                double *variances = b->stokesVector[column].getvariancepointer();
                unsigned n = b->length;
                while (n--) {
                    *tstokes += this->muellerMatrix[row][column] * *stokes;
                    *tvariances += pow(*stokes,2) * this->varianceMatrix[row][column] + pow(this->muellerMatrix[row][column],2) * *variances;
                    tstokes++; tvariances++;
                    stokes++; variances++;
                }
			}
		}
        return t;
	};
    
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the matrix on the left side of the operator with the Stokes vector on the right side.
	 * \param b An operaStokesVector address
	 * \note Usage: operaStokesVector a = operaMuellerMatrix c * operaStokesVector b;
	 * \return An operaStokesVector address
	 */
	operaStokesVector operator*(operaStokesVector& b) {
        operaStokesVector t(b.length);
		for (unsigned row = 0 ; row < 4 ; row++) {
			t.stokesVector[row] = 0.0;	
		}
		for (unsigned row = 0 ; row < 4 ; row++) {
			for (unsigned column = 0 ; column < 4 ; column++) {
                double *tstokes = t.stokesVector[row].getfluxpointer();
                double *tvariances = t.stokesVector[row].getvariancepointer();
                double *stokes = b.stokesVector[column].getfluxpointer();
                double *variances = b.stokesVector[column].getvariancepointer();
                unsigned n = b.length;
                while (n--) {
                    *tstokes += this->muellerMatrix[row][column] * *stokes;
                    *tvariances += pow(*stokes,2) * this->varianceMatrix[row][column] + pow(this->muellerMatrix[row][column],2) * *variances;
                    tstokes++; tvariances++;
                    stokes++; variances++;
                }
			}
		}
        return t;
	};
    
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the matrix on the left side of the operator with the doubleValue structure on the right side.
     * \details Each element is multiplied by the double value in the structure. The variances follow the normal error propagation rule.
	 * \param d A doubleValue structure
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c * doubleValue d;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator*(doubleValue d) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
        for (unsigned row = 0 ; row < 4 ; row++) {
            for (unsigned column = 0 ; column < 4 ; column++) {
                t->muellerMatrix[row][column] = this->muellerMatrix[row][column] * d.value;
                t->varianceMatrix[row][column] = pow(d.value,2) * this->varianceMatrix[row][column] + pow(this->muellerMatrix[row][column],2) * d.error;
            }
		}
		return *t;
	};
    
    /*! 
	 * \brief Division operator.
     * \details The operator divides the matrix on the left side of the operator by the doubleValue structure on the right side.
     * \details Each element is divided by the double value in the structure. The variances follow the normal error propagation rule.
	 * \param d A doubleValue structure
	 * \note Usage: operaMuellerMatrix a = operaMuellerMatrix c / doubleValue d;
	 * \return An operaMuellerMatrix address
	 */
	operaMuellerMatrix& operator/(doubleValue d) {
		operaMuellerMatrix *t = NULL;
        if (this->istemp)
			t = this; 
        else
            t = new operaMuellerMatrix(true);
        for (unsigned row = 0 ; row < 4 ; row++) {
            for (unsigned column = 0 ; column < 4 ; column++) {
                t->muellerMatrix[row][column] = this->muellerMatrix[row][column] / d.value;
                t->varianceMatrix[row][column] = pow(d.value,-2) * this->varianceMatrix[row][column] + pow(this->muellerMatrix[row][column]/pow(d.value,2),2) * d.error;
            }
		}
		return *t;
	};
};
#endif
