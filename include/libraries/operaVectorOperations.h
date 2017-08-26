#ifndef OPERAVECTOROPERATIONS_H
#define OPERAVECTOROPERATIONS_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaVectorOperations
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
#include <functional>

/*!
 * \brief Performs a unary operation on all elements of an operaVector, storing the results in that vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \return A reference to the vector post-operation.
 */
template <typename T, typename V> V& ApplyOperation(T funct, V& a);

/*!
 * \brief Performs a binary operation on all elements of two operaVectors, storing the results in the first vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \param b The second operand for the operation.
 * \return A reference to the vector post-operation.
 */
template <typename T, typename V> V& ApplyOperation(T funct, V& a, const V& b);

/*!
 * \brief Performs a binary operation with a constant value on all elements of an operaVector, storing the results in that vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \param d The second operand for the operation.
 * \return A reference to the vector post-operation.
 */
template <typename T, typename V> V& ApplyOperation(T funct, V& a, typename V::value_type d);

/*!
 * \brief Performs an operation that returns a single value from all elements of an operaVector.
 * \param funct The binary operation to apply between each element and the current result value.
 * \param a The vector of elements to use.
 * \param initalValue The value to start with for the result.
 * \return The value that results from applying funct with each value of the vector.
 * \details For example, the product of elements of a vector could be calculated using muliplication as funct, and an initial value of 1.
 */
template <typename T, typename V> typename V::value_type ScalarOperation(T funct, const V& a, const typename V::value_type initalValue);

/*!
 * \brief Performs an operation that returns a single value from all elements two operaVectors. Equivalent to applying a scalar operation to the result of a binary operation.
 * \param funct The operation to use for the scalar operation.
 * \param pairwiseOp The operation to use for the binary operation.
 * \param a The first vector of elements to use.
 * \param a The second vector of elements to use.
 * \param initalValue The initial value for the scalar operation.
 * \return The value that results from applying funct with each value of the vector that results from applying pairwiseOp on the two input vectors.
 * \details For example, the dot product of two vectors could be calculated using additon as funct, multiplication as pairwiseOp, and an initial value of 0.
 */
template <typename T, typename U, typename V> typename V::value_type ScalarOperation(T funct, U pairwiseOp, const V& a, const V& b, const typename V::value_type initalValue);

/*!
 * \brief Removes the type ambiguity from a function pointer for a unary operation by converting it into a functor.
 * \param funct The function pointer.
 * \return A function wrapper object.
 */
template <typename T> std::pointer_to_unary_function<T, T> ToFunctor(T (*funct)(T));

/*!
 * \brief Removes the type ambiguity from a function pointer for a binary operation by converting it into a functor.
 * \param funct The function pointer.
 * \return A function wrapper object.
 */	
template <typename T> std::pointer_to_binary_function<T, T, T> ToFunctor(T (*funct)(T, T));

/*!
 * \brief Performs an element-wise addition of one operaVector to another and assigns the result to the vector.
 * \param a The vector to be added to.
 * \param b The vector to be added.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator+=(Vector<T>& a, const Vector<T>& b);

/*!
 * \brief Adds a value to each element of a vector and assigns the result to the vector.
 * \param a The vector to be added to.
 * \param d The value to be added.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator+=(Vector<T>& a, T d);

/*!
 * \brief Performs an element-wise subtraction of one operaVector from another and assigns the result to the vector.
 * \param a The vector to be subtracted from.
 * \param b The vector to be subtracted.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator-=(Vector<T>& a, const Vector<T>& b);

/*!
 * \brief Subtracts a value from each element of the vector and assigns the result to the vector.
 * \param a The vector to be subtracted from.
 * \param d The value to be subtracted.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator-=(Vector<T>& a, T d);

/*!
 * \brief Performs an element-wise multiplication of one operaVector by another and assigns the result to the vector.
 * \param a The vector to be multiplied.
 * \param b The vector to be multiplied by.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator*=(Vector<T>& a, const Vector<T>& b);

/*!
 * \brief Mulitiplies each element of the vector by a value and assigns the result to the vector.
 * \param a The vector to be multiplied.
 * \param d The value to be multiplied by.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator*=(Vector<T>& a, T d);

/*!
 * \brief Performs an element-wise division of one operaVector by another and assigns the result to the vector.
 * \param a The vector to be divided.
 * \param b The vector to be divided by.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator/=(Vector<T>& a, const Vector<T>& b);

/*!
 * \brief Divides each element of the vector by a value and assigns the result to the vector.
 * \param a The vector to be divided.
 * \param d The value to be divided by.
 * \return A reference to the vector post-operation.
 */
template <typename T> Vector<T>& operator/=(Vector<T>& a, T d);

/*!
 * \brief Subtracts each element of the vector from a constant value.
 * \param d The value to be subtracted from.
 * \param a The vector to be subtracted.
 * \return A The resulting vector.
 */
template <typename T> Vector<T> operator-(T d, Vector<T> a);

/*!
 * \brief Divides a value by each element of a vector.
 * \param d The value to be divided.
 * \param a The vector to divide by.
 * \return A The resulting vector.
 */
template <typename T> Vector<T> operator/(T d, Vector<T> a);

/*!
 * \brief Gets the index of the minimum element of the vector.
 * \param a The vector.
 * \return The index of the minimum.
 */
template <typename V> unsigned MinIndex(const V& a);

/*!
 * \brief Gets the index of the maximum element of the vector.
 * \param a The vector.
 * \return The index of the maximum.
 */
template <typename V> unsigned MaxIndex(const V& a);

/*!
 * \brief Calculates the minimum element of the vector.
 * \param a The vector.
 * \return The minimum.
 */
template <typename V> typename V::value_type Min(const V& a);

/*!
 * \brief Calculates the minimum element of the vector.
 * \param a The vector.
 * \return The maximum.
 */
template <typename V> typename V::value_type Max(const V& a);

/*!
 * \brief Calculates the median of the elements in the vector.
 * \param a The vector.
 * \return The median.
 */
template <typename V> typename V::value_type Median(V a);

/*!
 * \brief Destructively calculates the median of the elements in the vector.
 * \param a The vector (warning: order of elements will be modified).
 * \return The median.
 */
template <typename V> typename V::value_type MedianQuick(V& a);

/*!
 * \brief Gets the absolute value of each element of an operaVector.
 * \param a The vector to take the absolute value of.
 * \return The resultant vector.
 */
template <typename V> V Abs(V a);

/*!
 * \brief Gets the square root of each element of an operaVector.
 * \param a The vector to take the square root of.
 * \return The resultant vector.
 */
template <typename V> V Sqrt(V a);

/*!
 * \brief Raises each element of an operaVector to a power.
 * \param a The vector to raise to a power.
 * \param d The power to raise each element to.
 * \return The resultant vector.
 */
template <typename V> V Pow(V a, typename V::value_type d);

/*!
 * \brief Calculates the inner-product (dot product) of one vector with another vector.
 * \param b An operaVector.
 * \param b An operaVector of the same size.
 * \return The calculated inner-product.
 */
template <typename V> typename V::value_type InnerProduct(const V& a, const V& b);

/*!
 * \brief Calculates the magnitude/norm of the vector.
 * \param a The vector.
 * \return The magnitude.
 */
template <typename V> typename V::value_type Magnitude(const V& a);

/*!
 * \brief Calculates the sum of the elements in the vector.
 * \param a The vector.
 * \return The sum.
 */
template <typename V> typename V::value_type Sum(const V& a);

/*!
 * \brief Calculates the mean of the elements in the vector.
 * \param a The vector.
 * \return The mean.
 */
template <typename V> typename V::value_type Mean(const V& a);

/*!
 * \brief Calculates the mean of the squares of elements in the vector.
 * \param a The vector.
 * \return The mean square.
 */
template <typename V> typename V::value_type MeanSquare(const V& a);

/*!
 * \brief Calculates the square root of the mean square of the vector.
 * \param a The vector.
 * \return The root mean square.
 */
template <typename V> typename V::value_type RMS(const V& a);

/*!
 * \brief Calculates the variance of the elements in the vector.
 * \param a The vector.
 * \return The variance.
 */
template <typename V> typename V::value_type Variance(const V& a);

/*!
 * \brief Calculates the variance of the elements in the vector from a given mean.
 * \param a The vector.
 * \param mean The already calculated mean.
 * \return The variance.
 */
template <typename V> typename V::value_type Variance(V a, typename V::value_type mean);

/*!
 * \brief Calculates the standard deviation of the elements in the vector.
 * \param a The vector.
 * \return The standard deviation.
 */
template <typename V> typename V::value_type StdDev(const V& a);

/*!
 * \brief Calculates the standard deviation of the elements in the vector from a given mean.
 * \param a The vector.
 * \param mean The already calculated mean.
 * \return The standard deviation.
 */
template <typename V> typename V::value_type StdDev(const V& a, typename V::value_type mean);

/*!
 * \brief Calculates the median absolute deviation of the elements in the vector.
 * \param a The vector.
 * \return The median absolute deviation.
 */
template <typename V> typename V::value_type MedianAbsDev(const V& a);

/*!
 * \brief Calculates the median absolute deviation of the elements in the vector from a given median.
 * \param a The vector.
 * \param median The already calculated median.
 * \return The median absolute deviation.
 */
template <typename V> typename V::value_type MedianAbsDev(V a, typename V::value_type median);

/*!
 * \brief Calculates the median standard deviation of the elements in the vector.
 * \param a The vector.
 * \details Note: divides by 0.674433 to convert from absolute deviation to standard deviation.
 * \return The median standard deviation.
 */
template <typename V> typename V::value_type MedianStdDev(const V& a);

/*!
 * \brief Calculates the median standard deviation of the elements in the vector.
 * \param a The vector.
 * \param median The already calculated median.
 * \details Note: divides by 0.674433 to convert from absolute deviation to standard deviation.
 * \return The median standard deviation.
 */
template <typename V> typename V::value_type MedianStdDev(const V& a, typename V::value_type median);

/*!
 * \brief Calculates the chi-square of the elements in the vector.
 * \param v The vector.
 * \param central The center to subtract from the vector.
 * \param degsfreedom The degrees of freedom.
 * \return The chi-square p-value.
 */
template <typename V> typename V::value_type ChiSqr(V v, typename V::value_type central, unsigned degsfreedom);

/*!
 * \brief Filter out elements matching specified condition from a vector.
 * \param a The vector to filter.
 * \return The resultant vector.
 */
template <typename T, typename V> V FilterElements(const V& a, T condition);

/*!
 * \brief Filter out NaNs from a vector.
 * \param a The vector to filter.
 * \return The resultant vector.
 */
template <typename V> V FilterNans(const V& a);

/*!
 * \brief Performs a unary operation on all elements of an operaVector, returning the result.
 * \param funct The operation to apply.
 * \param a The operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T, typename V> V Operation(T funct, V a) { return ApplyOperation(funct, a); }

/*!
 * \brief Performs a binary operation on all elements of two operaVectors, returning the result.
 * \param funct The operation to apply.
 * \param a The first operand for the operation.
 * \param b The second operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T, typename V> V Operation(T funct, V a, const V& b) { return ApplyOperation(funct, a, b); }

/*!
 * \brief Performs a binary operation with a constant value on all elements of an operaVector, returning the result.
 * \param funct The operation to apply.
 * \param a The first operand for the operation.
 * \param d The second operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T, typename V> V Operation(T funct, V a, typename V::value_type d) { return ApplyOperation(funct, a, d); }

/*!
 * \brief Performs a binary operation with a constant value on all elements of an operaVector, returning the result.
 * \param funct The operation to apply.
 * \param d The first operand for the operation.
 * \param a The second operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T, typename V> V Operation(T funct, typename V::value_type d, V a) { return ApplyOperation(funct, d, a); }

/*!
 * \brief Performs an element-wise addition of two operaVectors without modifying either operaVector.
 * \param a The first vector to be added.
 * \param b The second vector to be added.
 * \return The result of the addition.
 */
template <typename T> Vector<T> operator+(Vector<T> a, const Vector<T>& b) { return a += b; }

/*!
 * \brief Performs an element-wise subtraction of two operaVectors without modifying either operaVector.
 * \param a The vector of elements to subtract from.
 * \param b The vector of elements to subtract.
 * \return The result of the subtraction.
 */
template <typename T> Vector<T> operator-(Vector<T> a, const Vector<T>& b) { return a -= b; }

/*!
 * \brief Performs an element-wise multiplication of two operaVectors without modifying either operaVector.
 * \param a The first vector to be multiplied.
 * \param b The second vector to be multiplied.
 * \return The result of the multiplication.
 */
template <typename T> Vector<T> operator*(Vector<T> a, const Vector<T>& b) { return a *= b; }

/*!
 * \brief Performs an element-wise division of two operaVectors without modifying either operaVector.
 * \param a The vector of elements to be divided.
 * \param b The vector of elements to divide by.
 * \return The result of the division.
 */
template <typename T> Vector<T> operator/(Vector<T> a, const Vector<T>& b) { return a /= b; }

/*!
 * \brief Adds a value to each element of an operaVectorwithout modifying the original operaVector.
 * \param a The vector to be added.
 * \param d The value to be added.
 * \return The result of the addition.
 */
template <typename T> Vector<T> operator+(Vector<T> a, T d) { return a += d; }

/*!
 * \brief Adds a value to each element of an operaVectorwithout modifying the original operaVector.
 * \param d The value to be added.
 * \param a The vector to be added.
 * \return The result of the addition.
 */
template <typename T> Vector<T> operator+(T d, Vector<T> a) { return a += d; }

/*!
 * \brief Subtracts a value from each element of an operaVector without modifying the original operaVector.
 * \param a The vector to be subtracted from.
 * \param d The value to be subtracted.
 * \return The result of the subtraction.
 */
template <typename T> Vector<T> operator-(Vector<T> a, T d) { return a -= d; }

/*!
 * \brief Multiplies each element of an operaVector by a value without modifying the original operaVector.
 * \param a The vector to be multiplied.
 * \param d The value to be multiplied by.
 * \return The result of the multiplication.
 */
template <typename T> Vector<T> operator*(Vector<T> a, T d) { return a *= d; }

/*!
 * \brief Multiplies each element of an operaVector by a value without modifying the original operaVector.
 * \param d The value to be multiplied by.
 * \param a The vector to be multiplied.
 * \return The result of the multiplication.
 */
template <typename T> Vector<T> operator*(T d, Vector<T> a) { return a *= d; }

/*!
 * \brief Divides each element of an operaVector by a value without modifying the original operaVector.
 * \param a The vector to be divided.
 * \param d The value to be divided by.
 * \return The result of the division.
 */
template <typename T> Vector<T> operator/(Vector<T> a, T d) { return a /= d; }

#include "libraries/operaVectorOperations.tpp"

#endif
