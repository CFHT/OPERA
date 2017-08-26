
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: catz.cpp
 Version: 1.0
 Description: cat for zipped files.
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
#include "libraries/operaLib.h"				// for tempfile
#include "libraries/operaException.h"
#include "libraries/gzstream.h"

/*! \file catz.cpp */

using namespace std;

/*! 
 * catz
 * \author Doug Teeple
 * \brief cat for zipped files.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

void printUsageSyntax() {
	
	cout << " Usage: catz <zipped filename> -- uncompresses a zipped file to stdout\n";
	cout << " Usage: catz --compress <unzipped filename> ><zipped filename> -- compresses a regular file to stdout\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	unsigned compress = false;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"compress",		optional_argument, NULL, 'c'},
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "c:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'c':
					compress = true;
					break;
					
				case 'v':
					verbose = 1;
					break;
				case 'p':
					plot = 1;
					break;
				case 'd':
					debug = 1;
					break;
				case 't':
					trace = 1; 
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
		
		if (optind >= argc || !strcmp(argv[optind], "-")) {
			// stdin
			ofstream fout;
			string tmp = open_temp("/tmp/XXXXXX", fout);
			while (!cin.eof() && cin.good()) {
				string dataline;
				getline(cin, dataline);
				fout << dataline << endl;
			}
			fout.close();
			igzstream fin(tmp.c_str());
			if (fin.is_open()) {
				string dataline;
				while (fin.good()) {
					getline(fin, dataline);
					if (!fin.eof()) {
						cout << dataline << endl;
					}
				}
				fin.close();
			}
			remove(tmp.c_str());
		} else {
			while (optind < argc) {
				igzstream fin(argv[optind++]);
				if (fin.is_open()) {
					string dataline;
					while (fin.good()) {
						getline(fin, dataline);
						if (!fin.eof()) {
							cout << dataline << endl;
						}
					}
					fin.close();
				}
			}
		}
	}
	catch (operaException e) {
		cerr << "catz: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "catz: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

