#ifndef OPERAIOFORMATS_H
#define OPERAIOFORMATS_H

#include "libraries/operaSpectralOrderVector.h"

class FormatData {
private:
	stringstream ss;
	bool newline;
public:
	FormatData();
	string tostring() const;
	void skip(unsigned count);
	template <typename T> void insert(const T& value);
	template <typename T> T extract();
};

template <typename T> FormatData& operator<<(FormatData& data, const T& value);
template <typename T> FormatData& operator>>(FormatData& data, T& value);
FormatData& operator<<(FormatData& data, ostream& (*manip)(ostream&));

enum HeaderValue { newline = '\n', dots = '.' };

class FormatHeader {
private:
	string fullname;
	std::vector<string> headernames;
	std::vector<string> notes;
public:
	FormatHeader();
	FormatHeader(string fullname);
	void note(string note);
	string tostring() const;
	friend FormatHeader& operator<<(FormatHeader& header, string name);
	friend FormatHeader& operator<<(FormatHeader& header, HeaderValue value);
	friend FormatHeader& operator<<(FormatHeader& header, ostream& (*manip)(ostream&));
};

namespace operaIOFormats {
	void WriteCustomFormat(string formatname, const FormatHeader& formatheader, const FormatData& formatdata, string filename);
	
	void ReadCustomFormat(string formatname, FormatData& formatdata, string filename);
	
	/*! 
	 * \sa method void WriteSpectralOrder(operaSpectralOrderVector& orders, string filename, operaSpectralOrder_t format);
	 * \details Writes the specified format information from a spectral order vector to a file
	 * \details the optional order argument permits incremental addition to the output file, where zero means write all.
	 * \return none.
	 */
	void WriteFromSpectralOrders(const operaSpectralOrderVector& orders, string filename, operaSpectralOrder_t format);
	
	/*! 
	 * \sa method void ReadIntoSpectralOrders(operaSpectralOrderVector& orders, string filename);
	 * \brief augment an existing spectral order vector with information from a file
	 * \param filename - string.
	 * \return none.
	 */
	void ReadIntoSpectralOrders(operaSpectralOrderVector& orders, string filename);
}

template <typename T> void FormatData::insert(const T& value) {
	ostringstream temp;
	temp << value;
	if(temp.str() == "\n") newline = true;
	else {
		if(!newline) ss << ' ';
		newline = false;
	}
	ss << value;
}

template <typename T> T FormatData::extract() {
	T value;
	if(ss >> value) return value;
	throw operaException("FormatData: couldn't extract value to specified type ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
}

template <typename T> FormatData& operator<<(FormatData& data, const T& value) {
	data.insert(value);
	return data;
}

template <typename T> FormatData& operator>>(FormatData& data, T& value) {
	value = data.extract<T>();
	return data;
}

#endif
