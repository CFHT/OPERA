#ifndef VLARRAY
#define VLARRAY

/* A variable-length array that can be treated as a 1D or 2D array */
template <typename T>
class VLArray {
public:
	VLArray(const int length) : rows(1), cols(length) { arr = new T[cols]; }
	VLArray(const int rows, const int cols) : rows(rows), cols(cols) { arr = new T[rows * cols]; }
	~VLArray() { delete[] arr; }
	const T& operator() (const int i, const int j) const { return arr[i * cols + j]; }
	const T& operator() (const int i) const { return arr[i]; }
	T& operator() (int i, int j) { return arr[i * cols + j]; }
	T& operator() (int i) { return arr[i]; }
	const int Length() { return rows * cols; }
	const int& Rows() { return rows; }
	const int& Cols() { return cols; }
private:
	int rows;
	int cols;
	T* arr;
	VLArray(const VLArray& rhs);
};

#endif
