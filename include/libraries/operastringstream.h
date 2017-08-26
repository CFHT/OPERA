#ifndef OPERASTRINGSTREAM_H
#define OPERASTRINGSTREAM_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 ******************************************************************
 Library name: operaLib
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <sstream>
#include <string>

/*! 
 * \brief Float - reads nans and infs safely in a stringstream.
 * \file operastringstream.h
 * \ingroup libraries
 */

using namespace std;

class Float {
public:
	Float() {f = 0.0;};
	Float(float F) { f = F; };
	
	friend istream& operator >>(istream &is, Float &trouble);
	friend ostream& operator <<(ostream &os, const Float &trouble);
	
	float f;
	
	Float& operator=(float& F) {
		f = F;
        return *this;
	};
};

/*
 * Read from a stringstream, but before reading a floating point
 * directly, first put it into a local stringstream, so that if
 * it fails the original stream did not fail. If the local stream
 * read failed (i.e. nan or inf), then assign the internal floating
 * value to the result. Note that the input stream is advanced
 * in the process.
 */
inline istream& operator >>(istream &is, Float &result)
{
	result.f = 0.0;
	string str;
	is >> str;
	stringstream ss(str);
	ss >> result.f;
	if (ss.fail()) {
		if (str == "nan") {
			result.f = numeric_limits<float>::quiet_NaN();
		} else if (str == "inf" || str == "-inf") {
			result.f = numeric_limits<float>::infinity();
		}
	}
	return is;
}

inline ostream& operator <<(ostream &os, const Float &trouble)
{
	os << trouble.f;
	return os;
}

/*! 
 * \brief Double - reads nans and infs safely in a stringstream.
 * \file operastringstream.h
 * \ingroup libraries
 */

class Double {
public:
	Double() {d = 0.0;};
	Double(double D) { d = D; };
	
	friend istream& operator >>(istream &is, Double &trouble);
	friend ostream& operator <<(ostream &os, const Double &trouble);
	
	double d;
	
	Double& operator=(double& D) {
		d = D;
        return *this;
	};
};

inline istream& operator >>(istream &is, Double &result)
{
	result.d = 0.0;
	string str;
	is >> str;
	stringstream ss(str);
	ss >> result.d;
	if (ss.fail()) {
		if (str == "nan") {
			result.d = numeric_limits<double>::quiet_NaN();
		} else if (str == "inf" || str == "-inf") {
			result.d = numeric_limits<double>::infinity();
		}
	}
	return is;
}

inline ostream& operator <<(ostream &os, const Double &trouble)
{
	os << trouble.d;
	return os;
}

#endif

