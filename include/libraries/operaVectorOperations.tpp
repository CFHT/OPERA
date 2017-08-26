#include "libraries/operaException.h"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>

template <typename T, typename V>
V& ApplyOperation(T funct, V& a) {
	std::transform(a.begin(), a.end(), a.begin(), funct);
	return a;
}

template <typename T, typename V>
V& ApplyOperation(T funct, V& a, const V& b) {
	if(a.size() != b.size()) throw operaException("operaVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	std::transform(a.begin(), a.end(), b.begin(), a.begin(), funct);
	return a;
}

template <typename T, typename V>
V& ApplyOperation(T funct, V& a, typename V::value_type d) {
	std::transform(a.begin(), a.end(), a.begin(), std::bind2nd(funct, d));
	return a;
}

template <typename T, typename V>
V& ApplyOperation(T funct, typename V::value_type d, V& a) {
	std::transform(a.begin(), a.end(), a.begin(), std::bind1st(funct, d));
	return a;
}

template <typename T, typename V>
typename V::value_type ScalarOperation(T funct, const V& a, typename V::value_type initalValue) {
	return std::accumulate(a.begin(), a.end(), initalValue, funct);
}

template <typename T, typename U, typename V>
typename V::value_type ScalarOperation(T funct, U pairwiseOp, const V& a, const V& b, typename V::value_type initalValue) {
	if(a.size() != b.size()) throw operaException("operaVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	return std::inner_product(a.begin(), a.end(), b.begin(), initalValue, funct, pairwiseOp);
}

template <typename T>
Vector<T>& operator+=(Vector<T>& a, const Vector<T>& b) {
	return ApplyOperation(std::plus<T>(), a, b);
}

template <typename T>
Vector<T>& operator+=(Vector<T>& a, T d) {
	return ApplyOperation(std::plus<T>(), a, d);
}

template <typename T>
Vector<T>& operator-=(Vector<T>& a, const Vector<T>& b) {
	return ApplyOperation(std::minus<T>(), a, b);
}

template <typename T>
Vector<T>& operator-=(Vector<T>& a, T d) {
	return ApplyOperation(std::minus<T>(), a, d);
}

template <typename T>
Vector<T>& operator*=(Vector<T>& a, const Vector<T>& b) {
	return ApplyOperation(std::multiplies<T>(), a, b);
}

template <typename T>
Vector<T>& operator*=(Vector<T>& a, T d) {
	return ApplyOperation(std::multiplies<T>(), a, d);
}

template <typename T>
Vector<T>& operator/=(Vector<T>& a, const Vector<T>& b) {
	return ApplyOperation(std::divides<T>(), a, b);
}

template <typename T>
Vector<T>& operator/=(Vector<T>& a, T d) {
	return ApplyOperation(std::divides<T>(), a, d);
}

template <typename T>
Vector<T> operator-(T d, Vector<T> a) {
	return ApplyOperation(std::minus<T>(), d, a);
}

template <typename T>
Vector<T> operator/(T d, Vector<T> a) {
	return ApplyOperation(std::divides<T>(), d, a);
}

template <typename V>
V Abs(V a) {
	return ApplyOperation(ToFunctor<typename V::value_type>(std::abs), a);
}

template <typename V>
V Sqrt(V a) {
	return ApplyOperation(ToFunctor<typename V::value_type>(std::sqrt), a);
}

template <typename V>
V Pow(V a, typename V::value_type d) {
	return ApplyOperation(ToFunctor<typename V::value_type>(std::pow), a, d);
}

template <typename V>
typename V::value_type InnerProduct(const V& a, const V& b) {
	return ScalarOperation(std::plus<typename V::value_type>(), std::multiplies<typename V::value_type>(), a, b, 0);
}

template <typename V>
typename V::value_type Magnitude(const V& a) {
	return std::sqrt(InnerProduct(a, a));
}

template <typename V>
typename V::value_type Sum(const V& a) {
	return ScalarOperation(std::plus<typename V::value_type>(), a, 0);
}

template <typename V>
typename V::value_type Mean(const V& a) {
	return Sum(a) / a.size();
}

template <typename V>
typename V::value_type MeanSquare(const V& a) {
	return InnerProduct(a, a) / a.size();
}

template <typename V>
typename V::value_type RMS(const V& a) {
	return std::sqrt(MeanSquare(a));
}

template <typename V>
typename V::value_type Variance(const V& a) {
	return Variance(a, Mean(a));
}

template <typename V>
typename V::value_type Variance(V a, typename V::value_type mean) {
	a -= mean;
	return MeanSquare(a);
}

template <typename V>
typename V::value_type StdDev(const V& a) {
	return StdDev(a, Mean(a));
}

template <typename V>
typename V::value_type StdDev(const V& a, typename V::value_type mean) {
	return std::sqrt(Variance(a, mean));
}

template <typename V>
unsigned MinIndex(const V& a) {
	return std::min_element(a.begin(), a.end()) - a.begin();
}

template <typename V>
unsigned MaxIndex(const V& a) {
	return std::max_element(a.begin(), a.end()) - a.begin();
}

template <typename V>
typename V::value_type Min(const V& a) {
	return *(std::min_element(a.begin(), a.end()));
}

template <typename V>
typename V::value_type Max(const V& a) {
	return *(std::max_element(a.begin(), a.end()));
}

template <typename V>
typename V::value_type Median(V a) {
	return MedianQuick(a); //passes a copy
}

template <typename V>
typename V::value_type MedianQuick(V& a) {
	unsigned middle_index = a.size() / 2; //Middle index if size is odd, index after the middle if size is even
	typename V::iterator middle = a.begin() + middle_index;
	std::nth_element(a.begin(), middle, a.end()); //Partial sort putting element at middle in correct position, all elements less before it, all elements greater after it
	if(a.size() % 2 == 1) return *middle;
	typename V::value_type prev = *(std::max_element(a.begin(), middle)); //Find the largest element that is less than middle
	return (prev + *middle) / 2.0;
}

template <typename V>
typename V::value_type MedianAbsDev(const V& a) {
	return MedianAbsDev(a, Median(a));
}

template <typename V>
typename V::value_type MedianAbsDev(V a, typename V::value_type median) {
	a -= median;
	ApplyOperation(ToFunctor<typename V::value_type>(std::fabs), a);
	return MedianQuick(a);
}

template <typename V>
typename V::value_type MedianStdDev(const V& a) {
	return MedianStdDev(a, Median(a));
}

template <typename V>
typename V::value_type MedianStdDev(const V& a, typename V::value_type median) {
	return MedianAbsDev(a, median) / 0.674433; //magic number to convert from absolute deviation to standard deviation, assuming a normal distribution
}

template <typename V>
typename V::value_type ChiSqr(V v, typename V::value_type central, unsigned degsfreedom) {
	v -= central;
	return InnerProduct(v,v)/degsfreedom;
}

template <typename T, typename V>
V FilterElements(const V& a, T condition) {
	V temp;
	std::remove_copy_if(a.begin(), a.end(), std::back_inserter(temp), condition);
	return temp;
}

struct IsNan { // Better to make a functor since isnan isn't standard before C++11, and can be a macro
	bool operator()(double d) { return isnan(d); }
};

template <typename V>
V FilterNans(const V& a) {
	return FilterElements(a, IsNan());
}

template <typename T>
std::pointer_to_unary_function<bool, T> ToFunctorC(bool (*funct)(T)) {
	return std::ptr_fun(funct);
}

template <typename T>
std::pointer_to_unary_function<T, T> ToFunctor(T (*funct)(T)) {
	return std::ptr_fun(funct);
}

template <typename T>
std::pointer_to_binary_function<T, T, T> ToFunctor(T (*funct)(T, T)) {
	return std::ptr_fun(funct);
}
