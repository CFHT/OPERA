/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMJD
 Version: 1.0
 Description: Calculate the modified Julian date
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
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

#include "libraries/operaArgumentHandler.h"
#include "libraries/operaDateTime.h"
#include "libraries/operaException.h"

/*! \file operaMJD.cpp */

/*! 
 * operaMJD
 * \brief Calculate the modified Julian date
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

using namespace std;

int main(int argc, char *argv[]) {
    operaArgumentHandler args;
    
    double jd = 0.0;
    string datetime;
    string time;
    string date;
    args.AddOptionalArgument("julian", jd, 0.0, "input Julian date");
    args.AddOptionalArgument("time", time, string(), "input time (use with date, overrides JD) hh:mm:ss.sss");
    args.AddOptionalArgument("date", date, string(), "input date (use with time, overrides JD) yyyy-mm-dd");
    args.AddOptionalArgument("datetime", datetime, string(), "input cominbed date and time (overrides all) yyyy-mm-ddThh:mm:ss.sss");
    
	try {
		args.Parse(argc, argv);
		
		if (!datetime.empty()) {
			DateTime dt;
			if(!dt.SetFromString(datetime)) return EXIT_FAILURE;
			jd = dt.ToJulianDate();
		}
		else if(!time.empty() && !date.empty()) {
			DateTime dt;
			if(!dt.SetFromStrings(date, time)) return EXIT_FAILURE;
			jd = dt.ToJulianDate();
		}
		else if(jd == 0.0) {
			return EXIT_FAILURE;
		}
		cout.precision(7); //max precision to represent all unique values of double
		cout << fixed << JDtoMJD(jd) << endl;
	}
	catch (operaException e) {
		cerr << "operaMJD: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMJD: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
