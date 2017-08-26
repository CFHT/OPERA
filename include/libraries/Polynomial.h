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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "libraries/operaLibCommon.h"
#include "libraries/operaVector.h"

/*! 
 * \sa class Polynomial
 * \brief Encapsulates a polynomial. 
 * \file Polynomial.h
 * \ingroup libraries
 */
class Polynomial {
	
private:
	operaVector polynomialVector;
	operaVector polynomialErrors;
	double polychisqr;
	
public:
	/*
	 * Constructors
	 */
	Polynomial();
	Polynomial(const operaVector& coeffs);
	Polynomial(unsigned OrderOfPolynomial);	
	Polynomial(unsigned OrderOfPolynomial, const double* CoefficientVector);
	Polynomial(unsigned OrderOfPolynomial, const double* CoefficientVector, const double* CoefficientErrorVector);
	Polynomial(const PolynomialCoeffs_t* Coefficients);	
	Polynomial(const doublePolynomialCoeffs_t* Coefficients);	
	
	/*!
	 * const operaVector& getCoefficients()
	 * \brief This function returns an immutable reference to the vector of coefficients.
	 * \return const operaVector& - the immutable vector of polynomial coefficients
	 */
	const operaVector& getCoefficients() const;
	/*!
	 * const operaVector& getErrors()
	 * \brief This function returns an immutable reference to the vector of coefficient errors.
	 * \return const operaVector& - the immutable vector of polynomial coefficient errors
	 */
	const operaVector& getErrors() const;
	/*!
	 * void setCoefficients(const operaVector& coeffs, const operaVector& errors)
	 * \brief This function sets the coefficent values and errors.
	 * \param coeffs are the values to use for the coefficients
	 * \param errors are the values to use for the coefficient errors
	 * \return void
	 */
	void setCoefficients(const operaVector& coeffs, const operaVector& errors);
	/*!
	 * void setCoefficients(const operaVector& coeffs)
	 * \brief This function sets the coefficent values.
	 * \param coeffs are the values to use for the coefficients
	 * \return void
	 */
	void setCoefficients(const operaVector& coeffs);
	/*!
	 * double Get(const unsigned index)
	 * \brief This function gets the polynomial value at index.
	 * \note usage: float p3 = Get(3);
	 * \param index is a const unsigned index of the coefficient to get
	 * \return void
	 */
	double Get(const unsigned index) const;
	/*!
	 * void Set(const double x, const unsigned index)
	 * \brief This function sets the polynomial value at index.
	 * \note usage: Set((double)x, 3);
	 * \param x is a double value of the polynomial coefficient at index
	 * \return void
	 */
	void Set(const double x, const unsigned index);	
	/*!
	 * double Evaluate(double x)
	 * \brief This function returns the value of a given polynomial function.
	 * \note usage: double value = Evaluate((double)x);
	 * \param x is a double input value for which the given polynomial is evaluated
	 * \return double value of evaluation of the polynomial
	 */
	double Evaluate(const double x) const;
	/*!
	 * double double* getVector();
	 * \brief This function returns the double *polynomialVector.
	 * \note usage: double *vec = getvector();
	 * \return double * - the vector of polynomial order coefficients
	 */
	double* getVector();
	const double* getVector() const;
	/*!
	 * double double* getErrorVector();
	 * \brief This function returns the double *polynomialErrorVector.
	 * \note usage: double *errs = getvector();
	 * \return double * - the vector of polynomial coefficient errors
	 */
	double* getErrorVector();	
	const double* getErrorVector() const;	
	/*!
	 * unsigned getOrderOfPolynomial();
	 * \brief This function returns the unisgned polynomial Order.
	 * \note usage: unsigned npar = getOrderOfPolynomial<float>();
	 * \note usage: unsigned npar = getOrderOfPolynomial();
	 * \return double * - the vector of polynomial order coefficients
	 */
	unsigned getOrderOfPolynomial() const;
	/*!
	 * void resize(unsigned Order);
	 * \brief This function sets the unsigned polynomial Order.
	 * \note usage: resize(3);
	 * \return void
	 */
	void resize(unsigned Order);	
	/*!
	 * void double getCoefficient(unsigned index);
	 * \brief This function gets a coefficent at the index.
	 * \note usage: double coeff = getCoefficient(3);
	 * \return double
	 */
	double getCoefficient(unsigned index) const;
	/*!
	 * void setCoefficient(unsigned index, double value);
	 * \brief This function sets a coefficent at the index.
	 * \return void
	 */
	void setCoefficient(unsigned index, double value);	
	/*!
	 * void double getCoefficientError(unsigned index);
	 * \brief This function gets a coefficent error at the index.
	 * \note usage: double err = getCoefficientError(3);
	 * \return double
	 */
	double getCoefficientError(unsigned index) const;
	/*!
	 * void setCoefficientError(unsigned index, double value);
	 * \brief This function sets a coefficent error at the index.
	 * \note usage: setCoefficientError(3, 0.001);
	 * \return void
	 */
	void setCoefficientError(unsigned index, double value);	
	/*!
	 * void setPolynomialCoeffs(PolynomialCoeffs_t* coeffs);
	 * \brief This function sets a PolynomialCoeffs_t struct.
	 * \return void
	 */
	void setPolynomialCoeffs(const PolynomialCoeffs_t* coeffs);
	/*!
	 * void setChisqr(double Chisqr)
	 * \brief This function sets chisqr.
	 * \note usage: setChisqr(0.98);
	 * \return void
	*/
	void setChisqr(double Chisqr);	
	/*!
	 * double getChisqr(void)
	 * \brief This function gets chisqr.
	 * \note usage: double c = getChisqr();
	 * \return double
	 */
	double getChisqr(void) const;
	/*!
	 * void printEquation(ostream *pout)
	 * \brief This function prints the polynomial equation in Gnuplot format.
	 * \note usage: printEquation(&cout);
	 * \return void
	 */
	void printEquation(ostream *pout) const;
	/*!
	 * double operator()(double x)
	 * \brief Evaluates the polynomial function at x.
	 * \note usage: double value = polynomial((double)x);
	 * \param x is a double input value for which the given polynomial is evaluated
	 * \return double value of evaluation of the polynomial
	 */
	double operator()(double x) const;
};

#endif
