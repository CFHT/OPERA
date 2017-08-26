/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFITSProductTest
 Version: 1.0
 Description: Perform various tests on the operaFITSProductTest class.
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

#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaEspadonsImage.h"
#include "libraries/operaFITSProduct.h"

#include "libraries/operaImage.h"

/*! \file operaFITSProductTest.cpp */

using namespace std;

/*! 
 * operaFITSProductTest
 * \author Doug Teeple
 * \brief Perform various tests on the operaFITSProduct class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaFITSProductTest --product=<filename> --basedon=<filename> -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string productfilename, basedonfilename;
	string basedon;
	operaFITSProduct *product;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"product",			1, NULL, 'p'},	// output product
		{"basedon",			1, NULL, 'b'},	// based on image
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "p:b:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'p':		// output product
					productfilename = optarg;
					break;
				case 'b':		// based on image
					basedonfilename = optarg;
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
		
		if (productfilename.empty()) {
			throw operaException("operaFITSProductTest: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}


		// open the input file
		bool viewingOnly = true;
		if (!basedonfilename.empty()) {
			unsigned rows = 2000;
			unsigned cols = 3;
			product = new operaFITSProduct(productfilename, basedonfilename, MODE_UNKNOWN, cols, rows);
			viewingOnly = false;
			cout << "operaFITSProductTest: creating a new product " << productfilename << " based on  " << basedonfilename << endl;
		} else {
			product = new operaFITSProduct(productfilename, READONLY);
			cout << "operaFITSProductTest: reading an existing product " << productfilename << endl;
		}
		unsigned width = product->getcols();
		unsigned lines = product->getrows();
		cout << productfilename << ": " << width << "x" << lines << " " << "(" << product->getnaxis2() << "x" << product->getnaxis1() << ") ";
		switch (product->getmode()) {
			case MODE_STAR_ONLY:
				cout << "STAR_ONLY";
				break;
			case MODE_STAR_PLUS_SKY:
				cout << "STAR_PLUS_SKY";
				break;
			case MODE_POLAR:
				cout << "POLAR";
				break;
			default:
				break;
		}
		cout << "\n";
		for (unsigned w=0; w<width; w++) {
			cout << setw(14) << product->getcolname(w) << ' ';
		}
		cout << endl;
		if (viewingOnly) {
			for (unsigned w=0; w<width; w++) {
				for (unsigned l=0; l<20; l++) {
					cout << setw(14) << setprecision(4) << (float)*product[w][l] << ' ';
				}
				cout << endl;
			}
			
		}
		delete product;
	}
	catch (operaException e) {
		cerr << "operaFITSProductTest: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaFITSProductTest: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

