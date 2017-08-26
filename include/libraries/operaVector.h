#ifndef OPERAVECTOR_H
#define OPERAVECTOR_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaVector
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

#include <vector>

/*!
 * \brief This class provides a mapping from one set of indexes to another.
 * \ingroup libraries
 * \details
 * This container is designed to allow the re-indexing of a Vector.
 * In particular, it stores a one-to-one mapping between an old and a new set of indexes.
 * Note that the class has no public interface outside of its default methods.
 */
class operaIndexMap {
private:
	std::vector <unsigned> index;
	operaIndexMap(unsigned size);
public:
	template <typename T> friend class Vector;
};

/*!
 * \brief This class encapsulates the indexes for a range of elements.
 * \ingroup libraries
 * \details
 * This container is designed to allow the trimming of a Vector.
 * In particular, it stores the indexes for a subrange of Vector elements.
 * Note that the class has no public constructors aside from the copy constructor.
 */
class operaIndexRange {
private:
	unsigned start;
	unsigned end;
	operaIndexRange(unsigned startindex, unsigned endindex);
public:
	/*!
     * \brief Gets the size of the index range.
     * \return The number of elements in the range.
     */
	unsigned size();
	
	template <typename T> friend class Vector;
};

/*!
 * \brief This class encapsulates a vector of data.
 * \ingroup libraries
 * \details
 */
template <typename T>
class Vector {
private:
	std::vector <T> data;
	
	/*!
     * \brief Private constructor to prevent implicit conversions from non-integers values.
     * \details Normally, if a double is passed as an operaVector, it is implictly converted to int.
     */
	template <class U> Vector(U);
public:
	typedef T value_type;
	typedef typename std::vector<value_type>::iterator iterator;
	typedef typename std::vector<value_type>::const_iterator const_iterator;
	typedef const T& const_reference;
	
	/*!
     * \brief Creates an empty Vector.
     */
	Vector();
	
	/*!
     * \brief Creates a Vector by copying in an existing array.
     * \param dataarray The array to be copied.
     * \param length The length of the array.
     */
	Vector(const value_type* dataarray, unsigned length);
	
	/*!
     * \brief Creates a Vector with the specified size.
     * \param length The length of the vector.
     */
	Vector(unsigned length);
	Vector(int length);
	
	/*!
     * \brief Sets each element of the vector to the specified value.
     * \param value The value to set the vector to.
     */
	Vector& operator=(value_type value);
	
	/*!
     * \brief Gets the element at the given index.
     * \param i The index to access.
     * \return A reference to the element.
     */
	value_type& operator[](unsigned i);
	
	/*!
     * \brief Gets the element at the given index.
     * \param i The index to access.
     * \return A const reference to the element.
     */
	const value_type& operator[](unsigned i) const;
	
	/*!
     * \brief Gets an STL vector iterator for the beginning of the vector.
     * \return The begin iterator.
     */
	iterator begin();
	
	/*!
     * \brief Gets an STL vector constant iterator for the beginning of the vector.
     * \return The begin const_iterator.
     */
	const_iterator begin() const;
	
	/*!
     * \brief Gets an STL vector iterator for the end of the vector.
     * \return The end iterator.
     */
	iterator end();
	
	/*!
     * \brief Gets an STL vector constant iterator for the end of the vector.
     * \return The end const_iterator.
     */
	const_iterator end() const;
	
	/*!
     * \brief Gets the element at the first index in the vector.
     * \return The value of the first element.
     */
	value_type first() const;
	
	/*!
     * \brief Gets the element at the last index in the vector.
     * \return The value of the last element.
     */
	value_type last() const;
	
	/*!
     * \brief Gets the size of the vector.
     * \return The number of elements in the vector.
     */
	unsigned size() const;
	
	/*!
     * \brief Checks whether or not the vector is empty.
     * \return True if the vector is empty.
     */
	bool empty() const;
	
	/*!
     * \brief Removes all elements from the vector.
     * \details The vector will be empty after this.
     */
	void clear();
	
	/*!
     * \brief Resizes the vector to the specified size.
     * \param newsize The new length of the vector.
     * \details Existing elements in the vector are preserved up to the newsize.
     */
	void resize(unsigned newsize);
	
	/*!
     * \brief Gets the index range of elements between the specified values.
     * \param min The value that all elements in the subrange must be greater than or equal to.
     * \param max The value that all elements in the subrange must be less than or equal to.
     * \details The vector must already be sorted, and min must be less than or equal to max.
     * \return The index range of elements between min and max.
     */
	operaIndexRange subrange(value_type min, value_type max) const;
	
	/*!
     * \brief Removes all elments outside of the specified range from the vector.
     * \param range The index range of elements to keep.
     */
	void trim(operaIndexRange range);
	
	/*!
     * \brief Sets each element of the vector to the specified value.
     * \param value The value to fill the vector with.
     */
	void fill(value_type value);
	
	/*!
     * \brief Inserts a new value into the vector.
     * \param newdata The value to be inserted.
     * \details The value is inserted at the end of the vector and the size is increased by one.
     */
	void insert(value_type newdata);
    
    /*!
     * \brief An alias for the insert method.
     * \param newdata The value to be inserted.
     * \details Useful for standard library templates such as back_inserter.
     */
	void push_back(value_type newdata);
	
	/*!
     * \brief Inserts all the values from another vector into this vector.
     * \param newdata The vetor to be appended.
     * \details The values are inserted in the same order at the end of the vector.
     */
	void append(const Vector<T>& newdata);
	
	/*!
     * \brief Reverses the order of the elements in the vector.
     */
	void reverse();
	
	/*!
     * \brief Sorts the elements in the vector.
     * \details Elements are sorted in ascending order.
     */
	void sort();
	
	/*!
     * \brief Creates an operaIndexMap mapping from the current element order to their sorted order.
     * \details If the index map given by this method is used to reorder the vector, it will be in sorted order.
     * \return The index map indicating the sorted order of elements in the vector.
     */
	operaIndexMap indexsort() const;
	
	/*!
     * \brief Rearranges the order of elements in the vector according to an index map.
     * \param indexmap The index map indicating how the current ordering relates to the new ordering.
     * \details Multiple vectors of the same size can be reordered using the same index map.
     * The resulting elements will maintain their position relative to the other reordered vectors.
     * This method can be used with indexsort to emulate the sorting of tuples on a key value.
     */
	void reorder(const operaIndexMap &indexmap);
	
	/*!
     * \brief Copies the values from an array into the vector.
     * \param dataarray The array to be copied.
     * \details Assumes the array is equal in size to the current vector.
     */
	void copyfrom(value_type* dataarray);
	
	/*!
     * \brief Gets a pointer to the data stored in the vector.
     * \details The pointer can be treated as if the vector is a typical array.
     * The pointer may become invalidated after calling another non-const method of the vector.
     * \return A pointer to the the array which can be modified.
     */
	value_type* datapointer();
	
	/*!
     * \brief Gets a pointer to the data stored in the vector.
     * \details The pointer can be treated as if the vector is a typical array.
     * The pointer may become invalidated after calling another non-const method of the vector.
     * \return A pointer the array which cannot be modified.
     */
	const value_type* datapointer() const;
};

typedef Vector<double> operaVector;

#endif
