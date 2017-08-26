/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateInstrumentEnvironmentSetup
 Version: 1.0
 Description: Finds the location of orders.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Feb/2013
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

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaParameterAccess.h"
#include "libraries/operaConfigurationAccess.h"

#include "libraries/operaInstrumentEnvironmentSetup.h"
#include "operaCreateInstrumentEnvironmentSetup.h"

/*! \file operaCreateInstrumentEnvironmentSetup.cpp */

using namespace std;

/*! 
 * operaCreateInstrumentEnvironmentSetup
 * \author Eder Martioli
 * \brief Module to create instrument and environment setup classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
	int opt;
    
    string outputObservingConditonsFile;
    string outputObjectInTheSkyFile;
    string outputSpectrographFile;
    string outputTelescopeFile;
    
	int debug=0, verbose=0, trace=0;
    
	struct option longopts[] = {
		{"outputObservingConditonsFile",1, NULL, 'O'},
		{"outputObjectInTheSkyFile",1, NULL, 'J'},
		{"outputSpectrographFile",1, NULL, 'S'},
		{"outputTelescopeFile",1, NULL, 'T'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "O:J:S:T:v::d::t::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
			case 'O':		// outputObservingConditonsFile
				outputObservingConditonsFile = optarg;
				break;
			case 'J':		// outputObjectInTheSkyFile
				outputObjectInTheSkyFile = optarg;
				break;
			case 'S':		// outputSpectrographFile
				outputSpectrographFile = optarg;
				break;
            case 'T':		// outputTelescopeFile
				outputTelescopeFile = optarg;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 't':
				trace = 1;
				break;         
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	try {
		
        if (verbose) {
			cout << "operaCreateInstrumentEnvironmentSetup: outputObservingConditonsFile = " << outputObservingConditonsFile << endl;
			cout << "operaCreateInstrumentEnvironmentSetup: outputObjectInTheSkyFile = " << outputObjectInTheSkyFile << endl;
			cout << "operaCreateInstrumentEnvironmentSetup: outputSpectrographFile = " << outputSpectrographFile << endl;
			cout << "operaCreateInstrumentEnvironmentSetup: outputTelescopeFile = " << outputTelescopeFile << endl;
        }

        operaInstrumentEnvironmentSetup *instrumentEnvironment = new operaInstrumentEnvironmentSetup();
#if 0
		... get rid of warnings...
        operaObservingConditions *ObservingConditions = instrumentEnvironment->getObservingConditions();
        operaObjectInTheSky *ObjectInTheSky = instrumentEnvironment->getObjectInTheSky();
        operaSpectrograph *Spectrograph = instrumentEnvironment->getSpectrograph();
        operaTelescope *Telescope = instrumentEnvironment->getTelescope();
        
        
        
        
        
        // .....
        
        
        
#endif        
        /*
		 * write out the geometry info
		 */
		if (!outputObservingConditonsFile.empty()) {
            instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputObservingConditonsFile, ObsCond);
		}

		if (!outputObjectInTheSkyFile.empty()) {
            instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputObjectInTheSkyFile, SkyObj);
		}

		if (!outputSpectrographFile.empty()) {
            instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputSpectrographFile, spectrograph);
		}

		if (!outputTelescopeFile.empty()) {
            instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputTelescopeFile, telescope);
		}

        delete instrumentEnvironment;

	}
	catch (operaException e) {
		cerr << "operaCreateInstrumentEnvironmentSetup: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaCreateInstrumentEnvironmentSetup: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateInstrumentEnvironmentSetup: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --outputObservingConditonsFile=<OBSCOND_FILE>"
    " --outputObjectInTheSkyFile=<SKYOBJ_FILE>"
    " --outputSpectrographFile=<SPECTROGRAPH_FILE>"
    " --outputTelescopeFile=<TELESCOPE_FILE>\n\n"
	" Example: "+string(modulename)+"  \n\n"
	"  -O, --outputObservingConditonsFile,  Output observing conditions file\n"
	"  -J, --outputObjectInTheSkyFile,  Output sky object file\n"
	"  -S, --outputSpectrographFile,  Output spectrograph file\n"
	"  -T, --outputTelescopeFile,  Output telescope file\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -o, --outputGeomFile=<GEOM_FILE>, Geometry output calibration file name\n\n";
}
