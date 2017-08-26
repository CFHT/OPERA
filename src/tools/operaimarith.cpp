/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaimarith
 Version: 1.0
 Description: do simple image arithmetic.
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
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"

#include "libraries/operaImage.h"

/*! \file operaimarith.cpp */

using namespace std;

/*! 
 * operaimarith
 * \author Doug Teeple
 * \brief do simple image arithmetic.
 * \arg argc
 * \arg argv
 * \note arg0 = keyword 
 * \note arg1, arg2, ... FITS files
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout << " Usage: operaimarith [--add|--subtract|--multiply|--divide] <FITSImage>+ -[asmdvth]\n";
	cout << " Usage: performs simple image arithmetic.\n";
	cout << " Usage: Use --printfilename to identify which header the keyword is from.\n";
	cout << " Usage: The extension is optional ([n]) and can be the extension name.\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	enum op_t {
		add, subtract, multiply, divide
	} op = add;
	operaFITSImage *image, resultimage;
	int imagecount = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"add",				0, NULL, 'a'},
		{"multiply",		0, NULL, 'm'},
		{"subtract",		0, NULL, 's'},
		{"divide",			0, NULL, 'o'},
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "amsovdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'a':
					op = add;
					break;
				case 'm':
					op = multiply;
					break;
				case 's':
					op = subtract;
					break;
				case 'o':	
					op = divide;
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
		
		while (optind < argc) {
			image = new operaFITSImage(argv[optind], READONLY);
			if (imagecount == 0) {
				resultimage = new operaFITSImage(*image);
			} else {
				switch (op) {
					case add:
						resultimage += *image;
						break;
					case multiply:
						resultimage *= *image;
						break;
					case subtract:
						resultimage -= *image;
						break;
					case divide:
						resultimage /= *image;
						break;
					default:
						break;
				}
			}
			imagecount++;
			delete image;
		}
		resultimage.operaFITSImageSave();
		resultimage.operaFITSImageClose();
		optind++;
	}
	catch (operaException e) {
		cerr << "operaimarith: " << e.getFormattedMessage() << '\n';
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaimarith: " << operaStrError(errno) << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

