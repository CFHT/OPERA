/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operasaturated
 Version: 1.0
 Description: Determine the level of saturation of an image.
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

#include <getopt.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaImage.h"
#include "libraries/operaImageVector.h"

/*! \file operasaturated.cpp */

using namespace std;

/*! 
 * operasaturated
 * \author Doug Teeple
 * \brief Determine the level of saturation of an image.
 * \arg argc
 * \arg argv
 * \note --saturation=<float> - what is consideredd saturation (35000.0 ADU if not given) 
 * \note --count - number of pixels required to falg as saturated (default 500)
 * \return EXIT_STATUS ==1 if number of pixels abice saturation ADU exceeds count
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout << " Usage: operasaturated  [--count=<int> ] [--saturated=<float>] [<file name>]+ -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	int success = EXIT_SUCCESS;
	float saturation = 35000.0;
	unsigned count = 500;
	unsigned pixelcount = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"saturated",		1, NULL, 's'},	// saturation in ADU
		{"count",			1, NULL, 'c'},	// count to be considered saturated

		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "s:c:pvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 's':		// saturation
					saturation = atof(optarg);
					break;
				case 'c':		// count
					count  = atoi(optarg);
					break;
					
				case 'v':
					verbose = true;
					break;
				case 'd':
					debug = true;
					break;
				case 't':
					trace = true; 
					break;         
				case 'h':
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
				default:
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
			}	// switch
		}	// while
		
		if (optind == argc) {
			throw operaException("filename not given ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		string filename = argv[optind];
		operaFITSImage *fitsimage = new operaFITSImage(filename, tfloat, READONLY);
		//operaImageVector vector = where(*fitsimage >= saturation, &pixelcount);
		if (pixelcount > count) {
			cout << pixelcount << endl;
			success = EXIT_FAILURE;
		}
		delete fitsimage;
	}
	catch (operaException e) {
		cerr << "operasaturated: " << e.getFormattedMessage() << endl;
		success = EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operasaturated: " << operaStrError(errno) << endl;
		success = EXIT_FAILURE;
	}
	return success;
}

