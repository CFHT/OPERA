/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaConfigurationAccessTest
 Version: 1.0
 Description: Perform various tests on the oepraConfigurationAccess functions.
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
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "operaConfigurationAccesstest.h"
#include "libraries/operaException.h"

#include "libraries/operaParameterAccess.h"
#include "libraries/operaConfigurationAccess.h"

/*! \file operaConfigurationAccessTest.cpp */

using namespace std;

/*! 
 * operaConfigurationAccessTest
 * \author Doug Teeple
 * \brief Test the configuration access library.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 */

/**
 * \defgroup test Test Modules
 */

/*
 * main entry point for operaConfigurationAccessTest module.
 */
int main(int argc, char *argv[])
{
	int opt;
	string keywords[1000];
	char *value = (char *)malloc(1024);
	unsigned keys = 0;
	operaErrorCode e;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"keyword",			1, NULL, 'k'},	// input keyword
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "k:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'k':		// keywords
					keywords[keys++] = optarg;
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
		

		for (unsigned k=0; k<keys; k++) {
			if ((e=operaConfigurationAccessGet(keywords[k].c_str(), &value))) {
				throw operaException("operaConfigurationAccessTest: ", e, __FILE__, __FUNCTION__, __LINE__);
			} else if (value == NULL || strlen(value) == 0) {
				cout << keywords[k] << " not found\n";
			} else {
				cout << keywords[k] << " = " << value << '\n';
			}
		}
	}
	catch (operaException e) {
		cerr << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaConfigurationAccessTest: default error...\n";
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaConfigurationAccessTest [--keyword=<filename>]+ -[dvth]\n";
}	
