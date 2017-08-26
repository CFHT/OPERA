#include "globaldefines.h"
#include "libraries/operaVector.h"
#include "libraries/operaException.h"
#include <algorithm>

template <typename T>
class operaIndexComparator {
private:
    const std::vector<T> &data;
public:
    operaIndexComparator(const std::vector<T> &data) : data(data) {}
    bool operator()(unsigned i, unsigned j) { return data[i] < data[j]; }
};

operaIndexMap::operaIndexMap(unsigned size) : index(size) {
	for(unsigned i = 0; i < size; i++) index[i] = i;
}

operaIndexRange::operaIndexRange(unsigned startindex, unsigned endindex) : start(startindex), end(endindex) { }

unsigned operaIndexRange::size() {
	return end - start;
}

template <typename T>
Vector<T>::Vector() { }

template <typename T>
Vector<T>::Vector(unsigned length) : data(length) { }

template <typename T>
Vector<T>::Vector(int length) : data(length) { }

template <typename T>
Vector<T>::Vector(const T* dataarray, unsigned length) : data(dataarray, dataarray + length) { }

template <typename T>
Vector<T>& Vector<T>::operator=(T value) {
	fill(value);
	return *this;
}

template <typename T>
T& Vector<T>::operator[](unsigned i) {
#ifdef RANGE_CHECK
    if (i >= data.size()) {
		throw operaException("Vector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data[i];
}

template <typename T>
const T& Vector<T>::operator[](unsigned i) const {
#ifdef RANGE_CHECK
    if (i >= data.size()) {
		throw operaException("Vector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data[i];
}

template <typename T>
typename Vector<T>::iterator Vector<T>::begin() {
	return data.begin();
}

template <typename T>
typename Vector<T>::const_iterator Vector<T>::begin() const {
	return data.begin();
}
	
template <typename T>
typename Vector<T>::iterator Vector<T>::end() {
	return data.end();
}
	
template <typename T>
typename Vector<T>::const_iterator Vector<T>::end() const {
	return data.end();
}

template <typename T>
T Vector<T>::first() const {
#ifdef RANGE_CHECK
    if (data.empty()) {
		throw operaException("Vector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data.front();
}

template <typename T>
T Vector<T>::last() const {
#ifdef RANGE_CHECK
    if (data.empty()) {
		throw operaException("Vector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data.back();
}

template <typename T>
unsigned Vector<T>::size() const {
	return data.size();
}

template <typename T>
bool Vector<T>::empty() const {
	return data.empty();
}

template <typename T>
void Vector<T>::clear() {
	data.clear();
}

template <typename T>
void Vector<T>::resize(unsigned newsize) {
	data.resize(newsize);
}

template <typename T>
operaIndexRange Vector<T>::subrange(T min, T max) const {
	return operaIndexRange(std::lower_bound(data.begin(), data.end(), min) - data.begin(), std::upper_bound(data.begin(), data.end(), max) - data.begin());
}

template <typename T>
void Vector<T>::trim(operaIndexRange range) {
#ifdef RANGE_CHECK
    if (range.start >= data.size() || range.end > data.size() || range.start > range.end) {
		throw operaException("Vector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	std::vector<value_type> newdata(data.begin() + range.start, data.begin() + range.end);
	data.swap(newdata);
}

template <typename T>
void Vector<T>::fill(T value) {
	std::fill(data.begin(), data.end(), value);
}

template <typename T>
void Vector<T>::insert(T newdata) {
	data.push_back(newdata);
}

template <typename T>
void Vector<T>::push_back(T newdata) {
	data.push_back(newdata);
}

template <typename T>
void Vector<T>::append(const Vector& newdata) {
	data.insert(data.end(), newdata.data.begin(), newdata.data.end());
}

template <typename T>
void Vector<T>::reverse() {
	std::reverse(data.begin(), data.end());
}

template <typename T>
void Vector<T>::sort() {
	std::stable_sort(data.begin(), data.end());
}

template <typename T>
operaIndexMap Vector<T>::indexsort() const {
	operaIndexMap indexmap(data.size());
	std::sort(indexmap.index.begin(), indexmap.index.end(), operaIndexComparator<value_type>(data));
	return indexmap;
}

template <typename T>
void Vector<T>::reorder(const operaIndexMap &indexmap) {
	if(indexmap.index.size() != data.size()) throw operaException("Vector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	std::vector<value_type> newdata(data.size());
	for(unsigned i = 0; i < newdata.size(); i++) newdata[i] = data[indexmap.index[i]];
	data.swap(newdata);
}

template <typename T>
void Vector<T>::copyfrom(T* dataarray) {
	std::copy(dataarray, dataarray + data.size(), data.begin());
}

template <typename T>
T* Vector<T>::datapointer() {
	return &data[0];
}

template <typename T>
const T* Vector<T>::datapointer() const {
	return &data[0];
}

template class Vector<double>;
template class Vector<float>;
template class Vector<unsigned>;
