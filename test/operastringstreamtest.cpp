/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operastringstreamtest.cpp
 Version: 1.0
 Description: Can we read Nans and infs?
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA 
 Date: Mar/2013
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2013 Opera Pipeline team, Canada France Hawaii Telescope
 
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

#include <iostream>
#include <limits>
#include <libraries/operastringstream.h>

using namespace std;

int main()
{
	double inf = std::numeric_limits<double>::infinity();
	double qnan = std::numeric_limits<double>::quiet_NaN();

	cout << "This test check whether inf and nan can be read and written to/from a stream" << endl;
	cout << inf << ' ' << qnan <<  endl;
	
	Double dnan(0.0);
	Double dinf(0.0);
	double one = 0.0;
	double two = 0.0;
	
	stringstream qnanss(stringstream::in | stringstream::out);
	stringstream infss(stringstream::in | stringstream::out);
	
	qnanss << "nan 1.0";
	infss << "inf 1.0";
	
	qnanss >> dnan >> one;
	infss >> dinf >> two;
	
	cout << dinf << ' ' << one << ' ' << dnan <<  ' ' << two <<  endl;
	
	return 0;
}
