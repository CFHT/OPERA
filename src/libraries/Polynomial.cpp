/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: Polynomial
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

/*!
 * Polynomial
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the polynomial object.
 * \file Polynomial.cpp
 * \ingroup libraries
 */

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/Polynomial.h"
#include "libraries/operaLibCommon.h"

Polynomial::Polynomial() : polychisqr(0.0) { }

Polynomial::Polynomial(const operaVector& coeffs) : polynomialVector(coeffs), polynomialErrors(coeffs.size()), polychisqr(0.0) { }

Polynomial::Polynomial(unsigned OrderOfPolynomial) : polynomialVector(OrderOfPolynomial), polynomialErrors(OrderOfPolynomial), polychisqr(0.0) { }

Polynomial::Polynomial(unsigned OrderOfPolynomial, const double* CoefficientVector) : polynomialVector(CoefficientVector, OrderOfPolynomial), polynomialErrors(OrderOfPolynomial), polychisqr(0.0) { }

Polynomial::Polynomial(unsigned OrderOfPolynomial, const double* CoefficientVector, const double* CoefficientErrorVector) : polynomialVector(CoefficientVector, OrderOfPolynomial), polynomialErrors(CoefficientErrorVector, OrderOfPolynomial), polychisqr(0.0) { }

Polynomial::Polynomial(const PolynomialCoeffs_t* Coefficients) : polynomialVector(Coefficients->orderofPolynomial), polynomialErrors(Coefficients->orderofPolynomial), polychisqr(Coefficients->polychisqr) {
	for (unsigned i=0; i<polynomialVector.size(); i++) {
		polynomialVector[i] = Coefficients->p[i];
		polynomialErrors[i] = Coefficients->e[i];
	}
}

Polynomial::Polynomial(const doublePolynomialCoeffs_t* Coefficients) : polynomialVector(Coefficients->orderofPolynomial), polynomialErrors(Coefficients->orderofPolynomial), polychisqr(Coefficients->polychisqr) {
	for (unsigned i=0; i<polynomialVector.size(); i++) {
		polynomialVector[i] = Coefficients->p[i];
		polynomialErrors[i] = Coefficients->e[i];
	}
}

const operaVector& Polynomial::getCoefficients() const {
	return polynomialVector;
}

const operaVector& Polynomial::getErrors() const {
	return polynomialErrors;
}

void Polynomial::setCoefficients(const operaVector& coeffs, const operaVector& errors) {
	polynomialVector = coeffs;
	polynomialErrors = errors;
}

void Polynomial::setCoefficients(const operaVector& coeffs) {
	polynomialVector = coeffs;
}

double Polynomial::Get(const unsigned index) const {
	if (index > polynomialVector.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialVector[index];
}

void Polynomial::Set(const double x, const unsigned index) {
	if (index > polynomialVector.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialVector[index] = x;
}

double Polynomial::Evaluate(const double x) const {
	double total = 0;
	for(unsigned i = polynomialVector.size(); i > 0; i--) total = x*total + polynomialVector[i-1];
	return total;
}

double* Polynomial::getVector() {
	return polynomialVector.datapointer();
}
const double* Polynomial::getVector() const {
	return polynomialVector.datapointer();
}

double* Polynomial::getErrorVector() {
	return polynomialErrors.datapointer();
}
const double* Polynomial::getErrorVector() const {
	return polynomialErrors.datapointer();
}

unsigned Polynomial::getOrderOfPolynomial() const {
	return polynomialVector.size();
}

void Polynomial::resize(unsigned Order) {
	polynomialVector.resize(Order);
	polynomialErrors.resize(Order);
}

double Polynomial::getCoefficient(unsigned index) const {
	if (index > polynomialVector.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialVector[index];
}

void Polynomial::setCoefficient(unsigned index, double value) {
	if (index > polynomialVector.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialVector[index] = value;
}

double Polynomial::getCoefficientError(unsigned index) const {
	if (index > polynomialErrors.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialErrors[index];
}

void Polynomial::setCoefficientError(unsigned index, double value) {
	if (index > polynomialErrors.size()) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialErrors[index] = value;
}

void Polynomial::setPolynomialCoeffs(const PolynomialCoeffs_t* pcoefficients) {
	polynomialVector.resize(pcoefficients->orderofPolynomial);
	polynomialErrors.resize(pcoefficients->orderofPolynomial);
	polychisqr = pcoefficients->polychisqr;
	for (unsigned i=0; i<polynomialVector.size(); i++) {
		polynomialVector[i] = pcoefficients->p[i];
		polynomialErrors[i] = pcoefficients->e[i];
	}
}

void Polynomial::setChisqr(double Chisqr) {
	polychisqr = Chisqr;
}

double Polynomial::getChisqr(void) const {
	return polychisqr;
}

void Polynomial::printEquation(ostream *pout) const {
    if (pout != NULL) {
        *pout << "f(x) =";
        for(unsigned i=0;i<polynomialVector.size();i++) {
            if(i==0) {
                *pout << " " << polynomialVector[i];
            } else if (i==1) {
                *pout << " + " << polynomialVector[i] << "*x";
            } else {
                *pout << " + " << polynomialVector[i] << "*x**" << i;
            }
        }
        *pout << endl;
    }
}

double Polynomial::operator()(double x) const {
	return Evaluate(x);
}
