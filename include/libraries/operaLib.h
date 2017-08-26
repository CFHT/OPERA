#ifndef OPERALIB_H
#define OPERALIB_H

/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
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

#include <stdarg.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <locale>
#include <limits>

/*! 
 * \brief general C++ library routines.
 * \file operaLib.h
 * \ingroup libraries
 */

using namespace std;

/* 
 * trimFITSKeyword(const char *in)
 * \brief Trim the single quotes and trailing spaces from a FITS keyword value.
 * \param in
 * \return string trimmed value
 */
string trimFITSKeyword(const char *in);

/*
 * systemf(const char *form, ...);
 * like printf except it does a system call
 */
void systemf(const char *form, ...);

/*
 * itos(int i)
 * integer to C++ string
 */
inline string itos(int i){
	stringstream out;
	out << i;
	return out.str();
}

/*
 * ftos(float f)
 * float to C++ string
 */
inline string ftos(float f){
	stringstream out;
	out << f;
	return out.str();
}

/*
 * dtos(double d)
 * double to C++ string
 */
inline string dtos(double d){
	stringstream out;
	out << d;
	return out.str();
}

/* 
 * stot (string to type): 
 * The third parameter of stot() should be one of std::hex, std::dec or std::oct
 * 
 * if (stot<int>(i, string("ff"), hex)) {
 *	cout << i << endl;
 * }
 * if (stot<float>(f, string("123.456"), dec)) {
 *	cout << f << endl;
 * }
 */ 
template <class T>
bool stot(T& t, const string& s, ios_base& (*f)(ios_base&) = dec)
{
	istringstream iss(s);
	return !(iss >> f >> t).fail();
}

/* 
 * bool parsefloating(T& result, stringstream &ssin): 
 * Parse a stringstream for a floating point number,
 * handles nan and inf in stringstreams, returns
 * true on encountering a nan or inf, else false.
 * Has side-effect of advancing the stream properly.
 *
 */ 
template <class T>
bool parsefloating(T& result, stringstream& ssin)
{
	result = 0.0;
	string str;
	ssin >> str;
	stringstream ss(str);
	ss >> result;
	if (ss.fail()) {
		if (str == "nan") {
			result = numeric_limits<double>::quiet_NaN();
			return true;
		} else if (str == "inf" || str == "-inf") {
			result = numeric_limits<double>::infinity();
			return true;
 		}
	}
	return false;
}

/*
 * string open_temp(std::string path, std::ofstream& f)
 * create a unique temporary file
 * usage: 
 *    ofstream tempfile;
 *    open_temp("/tmp/XXXXXX", tempfile);
 *    if(tempfile.is_open()) {
 *		tempfile << ...;
 *    }
 */
string open_temp(std::string path, std::ofstream& f);

/*
 * bool fileexists(string filename );
 * check to see if a file exists without opening it
 */
bool fileexists(string filename);

/*
*
* directoryexists(string directoryname)
* Check if a directory exists without opening it.
*/
bool directoryexists(string directoryname);

/*
 * bool getRealFileName(string &filename, string &actualfilename)
 * Follow links to get real filename.
 */
bool getRealFileName(string &filename, string &actualfilename);

/*
 * string lowerCase(string &in)
 * Convert in to lower case.
 */
string lowerCase(string &in);

/*
 * string upperCase(string &in)
 * Convert in to upper case.
 */
string upperCase(string &in);

/*
 * int mkpath(const char *path, mode_t mode)
 * make a full path.
 */
int mkpath(const char *path, mode_t mode);

#endif
