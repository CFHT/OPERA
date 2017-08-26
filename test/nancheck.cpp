/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: nancheck.cpp
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

#include <cstdio>
#include <iostream>
#include <limits>

using namespace std;

int main()
{
	double myd = std::numeric_limits<double>::infinity();
	double myd2 = std::numeric_limits<double>::quiet_NaN();
	double myd3 = std::numeric_limits<double>::signaling_NaN();

	cout << "This test check whether inf and nan can be read from a stream and scanf" << endl;
	cout << "It prints inf, nan, nan and waits for input; type: inf nan nan" << endl;
	cout << "and it will print out what it gets two ways, first istream and then scanf" << endl;
	cout << myd << ' ' << myd2 <<  ' ' << myd3 << endl;
	cin >> myd >> myd2 >> myd3;
	cout << "inf? " << myd << endl;
	cout << "nan? " << myd2 << endl;
	cout << "nan? " << myd3 << endl;
	scanf("%lf %lf %lf", &myd, &myd2, &myd3);
	cout << "inf? " << myd << endl;
	cout << "nan? " << myd2 << endl;
	cout << "nan? " << myd3 << endl;
	return 0;
}
