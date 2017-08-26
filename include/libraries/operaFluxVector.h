#ifndef OPERAFLUXVECTOR_H
#define OPERAFLUXVECTOR_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaFluxVector
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Mar/2012
 
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

#include "libraries/operaVector.h"
#include "libraries/operaVectorOperations.h"
#include <utility>	// for pair

/*!
 * \file operaFluxVector.h
 */

/*!
 * \brief Definition of the value of each "tends towards" optional field.
 */
typedef enum {ToDefault, ToINF, ToNAN, ToZero, ToOne} TendsTowards_t;

/*!
 * \author Doug Teeple
 * \author Andre Venne
 * \brief This class encapsulates the flux vector.
 * \ingroup libraries
 * \details
 * 
 * The flux vector stores fluxes and variances. The operators are defined to
 * operate on operaFluxVectors and propagate the variances.
 * 
 * The variances are calculated as followed :
 * 
 * F = F(a,b)
 * DF = Pow(dF/da,2) * Da + Pow(dF/db,2) * Db
 * 
 * where DF is the resulting variance, Da and Db are the variance of the fluxes a and b, dF/da and dF/db are the partial derivatives of F.
 * The fluxes are supposed uncorrelated.
 *
 * The flux vector also has an optional "tends towards" field which allows control of INF / INF situations,
 * where the results can tend towards INF, NaN, 0.0, 1.0 or default result.
 */
class operaFluxVector {
	
private:
	TendsTowards_t towards;
    operaVector flux;
    operaVector variance;
	
public:
	/*!
     * \brief operaFluxVector constructor.
     * \details This constructor creates an operaFluxVector.
     * \param tendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(TendsTowards_t towards=ToDefault);
	
	/*!
     * \brief operaFluxVector constructor.
     * \details This constructor creates an operaFluxVector.
     * \param length Number of elements in the operaFluxVector
     * \param tendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(unsigned length, TendsTowards_t towards=ToDefault);
    
    /*!
     * \brief operaFluxVector constructor from flux and variance vectors.
     * \details This constructor creates an operaFluxVector by copying the contents of two equal sized operaVectors holding the flux and variance.
     * \param fluxes An operaVector to copy into the flux
     * \param variances An operaVector to copy into the variance
     * \param towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(const operaVector& fluxes, const operaVector& variances, TendsTowards_t towards=ToDefault);
    
    /*!
     * \brief operaFluxVector constructor from flux and variance arrays.
     * \details This constructor creates an operaFluxVector by copying the contents of two equal sized arrays holding the flux and variance.
     * \param fluxes An array with of the given length
     * \param variances An array with of the given length
     * \param length Number of elements in the operaFluxVector
     * \param towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(double *fluxes, double *variances, unsigned length, TendsTowards_t towards=ToDefault);
    
    /*!
     * \brief operaFluxVector constructor from an operaFluxVector.
     * \details This constructor creates an operaFluxVector by copying the contents of an existing operaFluxVector.
     * \param b An operaFluxVector to be copied
     * \param towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(const operaFluxVector &b, TendsTowards_t towards=ToDefault);
	
	/*!
     * \brief operaFluxVector constructor from an operaFluxVector.
     * \details This constructor creates an operaFluxVector by copying the flux from of an operaVector, leaving the variance at 0.
     * \param b An operaVector to be copied
     * \param towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(const operaVector &b, TendsTowards_t towards=ToDefault);
	
    /*!
     * \brief Removes all elements from the operaFluxVector.
     * \details The flux and variance vectors will be empty after this.
     */
    void clear();
	
	/*!
     * \brief Resizes the operaFluxVector.
     * \details Resizes the flux and variance vectors, and sets length to the new size. Elements are preserved up to min(length, newlength).
	 * \param newlength The new length to resize to
     * \return void
     */
	void resize(unsigned newlength);
	
	/*!
     * \brief Removes all elments outside of the specified range from the operaFluxVector.
     * \param range The index range of elements to keep.
     * \return void
     */
	void trim(operaIndexRange range);
	
	/*!
     * \brief Inserts new values into the operaFluxVector.
     * \details Inserts new flux and variance values at the end of the operaFluxVector. Length increases by one.
	 * \param newflux The flux value to insert.
	 * \param newvariance The variance value to insert.
     * \return void
     */
	void insert(double newflux, double newvariance);
	
	/*!
     * \brief Reverses the operaFluxVector.
     * \details Reverses the order of the elements of the flux and variance vectors.
     * \return void
     */
	void reverse();
	
	/*!
	 * \brief Re-arranges the order of elements in the operaFluxVector.
	 * \details Changes the order of elements in the flux and variance vectors so that the element at position indexmap[n] is moved to position n
	 * \param indexmap A map from the desired new index to the current index
	 * \return void
	 */
	void reorder(const operaIndexMap& indexmap);
	
    /*!
     * \brief Gets the flux vector.
     * \return An immutable reference to the flux vector
     */
	const operaVector &getflux() const;
	
    /*!
     * \brief Gets the variance vector.
     * \return An immuatable reference to the variance vector.
     */
	const operaVector &getvariance() const;
	
	/*!
     * \brief Sets the flux vector.
     * \param newflux The new flux vector of the same size.
     * \return void
     */
	void setflux(const operaVector& newflux);
    
    /*!
     * \brief Sets the variance vector.
     * \param newvariance The new variance vector of the same size.
     * \return void
     */
	void setvariance(const operaVector& newvariance);
	
	/*!
     * \brief Sets the flux vector to a constant value.
     * \param Flux The new flux value
     * \return void
     */
	void setflux(double newflux);
    
    /*!
     * \brief Sets the variance vector to a constant value.
     * \param Variance The new variance value
     * \return void
     */
	void setvariance(double newvariance);
    
    /*!
     * \brief Gets the flux vector as an array.
     * \return A pointer to the flux vector
     */
	double* getfluxpointer();
	const double* getfluxpointer() const;
	
    /*!
     * \brief Gets the variance vector as an array.
     * \return A pointer to the variance vector
     */
	double* getvariancepointer();
	const double* getvariancepointer() const;
	
    /*!
     * \brief Gets a flux vector element.
     * \param index The index to get
     * \return The flux at index
     */
	double getflux(unsigned index) const;
	
    /*!
     * \brief Gets a variance vector element.
     * \param index The index to get
     * \return The variance at index
     */
	double getvariance(unsigned index) const;
	
    /*!
     * \brief Sets a flux vector element.
     * \param Flux The new flux value
     * \param index The index to set
     * \return void
     */
	void setflux(double newflux, unsigned index);
    
    /*!
     * \brief Sets a variance vector element.
     * \param Variance The new variance value
     * \param index The index to set
     * \return void
     */
	void setvariance(double newvariance, unsigned index);
    
    /*!
     * \brief Gets the number of elements in the operaFluxVector.
     * \return The length
     */
	unsigned getlength() const;
    
    /*!
     * \brief Gets the TendsTowards_t value of operaFluxVector.
     * \return The TendsTowards_t value to which the operaFluxVector will tend to in case of infinity operations.
     */
    TendsTowards_t gettowards(void) const;
	
    /*!
     * \brief Gets the error of a flux element.
     * \details Gets the square root of a variance vector element, which is the error of the flux element.
     * \param index The index to get
     * \return The error at index
     */
	double geterror(unsigned index) const;
	
	/*! 
	 * \brief Indexing operator.
	 * \param i The index to get
	 * \note Usage: pair<double,double> p = FluxVector[i];
     * \note To print : cout << p.first << p.second;
	 * \return A pair of doubles containing the flux and variance at the given index
	 */
	std::pair<double,double> operator[](unsigned index) const;
	
	/*!
     * \brief Assignment operator.
     * \details This constructor sets the flux vector to the given operaVector and every element of the variance vector to 0.
     * \param b The operaVector to be copied
     */
	operaFluxVector& operator=(const operaVector& b);
	
    /*! 
	 * \brief Assignment operator.
     * \details The operator sets every element of the flux vector to the given double and every element of the variance vector to 0.
     * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator=(double d);
	
    /*!
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The resulting variances will be given by var(a+b) = var(a) + var(b).
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */	
	operaFluxVector& operator+=(const operaFluxVector& b);
	
	/*! 
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator+=(double d);
	
	/*!
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The resulting variances will be given by var(a-b) = var(a) + var(b).
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator-=(const operaFluxVector& b);
	
    /*! 
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator-=(double d);
	
    /*!
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The resulting variances will be given by var(a*b) = var(a) * b^2 + var(b) * a^2.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator*=(const operaFluxVector& b);
	
    /*! 
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the double from the right side of the operator to every value of the flux vector on the left side.
     * \details The resulting variances will be given by var(a*d) = var(a) * d^2.
     * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator*=(double d);
	
    /*!
	 * \brief Division/assignment operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side and copies them to the left side.
     * \details The resulting variances will be given by var(a/b) = (var(a) * b^2 + var(b) * a^2) / b^4.
     * \details In the case of infinite values, the resulting flux and variance values will be determined by the value of towards.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator/=(const operaFluxVector& b);
	
    /*! 
	 * \brief Division/assignment operator.
     * \details The operator divides every value of the flux vector on the left side of the operator by the double on the right side and copies them to the left side.
     * \details The resulting variances will be given by var(a/d) = var(a) / d^2.
     * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator/=(double d);

	/*! 
	 * \brief Subtraction operator.
	 * \details Subtracts each flux value of an operaFluxVector from a constant flux value.
	 * \details The resulting variances will be the same as the initial variances.
	 * \param d The flux value
	 * \param a An operaFluxVector
	 * \return The resulting difference operaFluxVector
	 */
	friend operaFluxVector operator-(double d, operaFluxVector a);
	
	/*! 
	 * \brief Divison operator.
	 * \details Divides a flux value by an operaFluxVector without modifying it.
	 * \details The resulting variances will be given by var(d/a) = var(a) * d^2 / a^4.
	 * \param d The flux value
	 * \param a An operaFluxVector
	 * \return The resulting difference operaFluxVector
	 */
	friend operaFluxVector operator/(double d, operaFluxVector a);
	
	/*!
	 * \brief Square root function.
	 * \details The function applies a square root to every element in the flux vector.
	 * \details The resulting variances will be given by var(sqrt(b)) = var(b) / 4b.
	 * \param b An operaFluxVector
	 * \return The resulting square root operaFluxVector
	 */
	friend operaFluxVector Sqrt(operaFluxVector b);

	/*!
	 * \brief Power function.
	 * \details The function raises every element in the flux vector to the specified power.
	 * \details The resulting variances will be given by var(b^d) = var(b) * b^(2(d-1)) * d^2.
	 * \param b An operaFluxVector
	 * \param d The power
	 * \return The resulting power operaFluxVector
	 */
	friend operaFluxVector Pow(operaFluxVector b, double d);

	/*!
	 * \brief Sum function.
	 * \details Calculuates the sum of a flux vector.
	 * \param b An operaFluxVector
	 * \return A pair of doubles containing the sums for the fluxes and variances respspectively
	 */
	friend std::pair<double,double> Sum(const operaFluxVector& b);

	/*!
	 * \brief Mean function.
	 * \details Calculuates the mean of a flux vector.
	 * \param b An operaFluxVector
	 * \return A pair of doubles containing the means of the fluxes and variances respspectively
	 */
	friend std::pair<double,double> Mean(const operaFluxVector& b);
};

/*!
 * \brief Addition operator.
 * \details Adds two operaFluxVectors without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting sum operaFluxVector
 */
inline operaFluxVector operator+(operaFluxVector a, const operaFluxVector& b) { return a += b; }

/*! 
 * \brief Addition operator.
 * \details Adds a flux value to an operaFluxVector without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting sum operaFluxVector
 */
inline operaFluxVector operator+(operaFluxVector a, double d) { return a += d; }
inline operaFluxVector operator+(double d, operaFluxVector a) { return a += d; }

 /*!
 * \brief Subtraction operator.
 * \details Subtracts one operaFluxVector from another without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting difference operaFluxVector
 */
inline operaFluxVector operator-(operaFluxVector a, const operaFluxVector& b) { return a -= b; }

/*! 
 * \brief Subtraction operator.
 * \details Subtracts a flux value from an operaFluxVector without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting difference operaFluxVector
 */
inline operaFluxVector operator-(operaFluxVector a, double d) { return a -= d; }

/*!
 * \brief Multiplication operator.
 * \details Multiplies two operaFluxVectors without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting product operaFluxVector
 */
inline operaFluxVector operator*(operaFluxVector a, const operaFluxVector& b) { return a *= b; }

/*! 
 * \brief Multiplication operator.
 * \details Multilpies an operaFluxVector by a flux value without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting product operaFluxVector
 */
inline operaFluxVector operator*(operaFluxVector a, double d) { return a*=d; }
inline operaFluxVector operator*(double d, operaFluxVector a) { return a*=d; }

/*!
 * \brief Division operator.
 * \details Divides one operaFluxVector by another without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting quotient operaFluxVector
 */
inline operaFluxVector operator/(operaFluxVector a, const operaFluxVector& b) { return a /= b; }

/*! 
 * \brief Division operator.
 * \details Divides an operaFluxVector by a flux value without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting quotient operaFluxVector
 */
inline operaFluxVector operator/(operaFluxVector a, double d) { return a/=d; }

#endif
