/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaLib - common C++ library functions
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

/*!
 * operaLib
 * \author Doug Teeple
 * \brief This class encapsulates C++ helper functions.
 * \file operaLib.cpp
 * \ingroup libraries
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <vector>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLib.h"

using namespace std;

/*
 * systemf(const char *format, ...);
 * like printf except it does a system call
 */
void systemf(const char *format, ...)
{
	static char line[32768];
	va_list argp;
	
	va_start(argp, format);
	vsnprintf(line, 32768, format, argp);
	system(line);
	va_end(argp);
}

/* 
 * trimFITSKeyword(const char *in)
 * \brief Trim the single quotes and trailing spaces from a FITS keyword value.
 * \param in
 * \return string trimmed value
 */
string trimFITSKeyword(const char *in) {
	string *str = new string(in);
	if (in[0] == '\'') {
		str->erase(str->end()-1);
		str->erase(str->begin());
		const char *end = in + str->length();
		while (*end-- == ' ')
			str->erase(str->end()-1);		
	}
	return *str;
}

/*
 * string open_temp(string path, ofstream& f)
 * create and open a unique temporary file
 * usage: 
 *    ofstream tempfile;
 *    open_temp("/tmp/XXXXXX", tempfile);
 *    if(tempfile.is_open()) {
 *		tempfile << ...;
 *    }
 */
string open_temp(string path, ofstream& f) {
	
    vector<char> destinationPath(path.begin(), path.end());
    destinationPath.push_back('\0');
	
    int fd = mkstemp(&destinationPath[0]);
    if (fd != -1) {
        path.assign(destinationPath.begin(), destinationPath.end() - 1);
        f.open(path.c_str(), std::ios_base::trunc | std::ios_base::out);
        close(fd);
    }
    return path;
}

/*
 * bool fileexists(string filename )
 * Check if a file exists without opening it.
 */
bool fileexists(string filename) {
	struct stat buffer;
	if ( stat(filename.c_str(), &buffer) ) 
		return false;
	return true;
}

/*
 * directoryexists(string directoryname)
 * Check if a directory exists without opening it.
 */
bool directoryexists(string directoryname) {
	struct stat buffer;
	if (stat(directoryname.c_str(), &buffer) == 0) {
		if ((buffer.st_mode & S_IFDIR) != 0) {
			return true;
		}
	}
	return false;
}

/*
 * string lowerCase(string &in)
 * Convert in to lower case.
 */
string lowerCase(string &in) {
	locale loc;
	for (size_t i=0; i < in.length(); i++ ) {
		in[i] = tolower(in[i], loc);
	}
	return in;
}

/*
 * bool getRealFileName(string &filename, string &actualfilename)
 * Follow links to get real filename.
 */
bool getRealFileName(string &filename, string &actualfilename) {
	char buff[2048];
	int count = readlink(filename.c_str(), buff, sizeof(buff));
	if (count > 0) {
		buff[count] = '\0';
		actualfilename = string(buff);
		return true;
	}
	return false;
}

/*
 * string upperCase(string &in)
 * Convert in to upper case.
 */
string upperCase(string &in) {
	locale loc;
	for (size_t i=0; i < in.length(); i++ ) {
		in[i] = toupper(in[i], loc);
	}
	return in;
}

typedef struct stat Stat;

static int do_mkdir(const char *path, mode_t mode)
{
    Stat            st;
    int             status = 0;
	
    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }
	
    return(status);
}

int mkpath(const char *path, mode_t mode)
{
    char	*pp;
    char	*sp;
    int		status;
    char	*copypath = strdup(path);
	
    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
	
    if (status == 0)
        status = do_mkdir(path, mode);
	
    free(copypath);
    return (status);
}


