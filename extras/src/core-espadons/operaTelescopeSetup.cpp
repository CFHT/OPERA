/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaTelescopeSetup
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
#include "libraries/operaTelescope.h"

#include "core-espadons/operaTelescopeSetup.h"

/*! \file operaTelescopeSetup.cpp */

using namespace std;

/*! 
 * operaTelescopeSetup
 * \author Eder Martioli
 * \brief Module to create telescope setup class.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
    
    string outputTelescopeFile;
    
    double latitude=19.8266666;              // in deg, (+) North, (-) South
    double longitude=-155.4716666;           // in deg, +East
    double elevation=4207;                   // in meters

    double CollectingArea=8.414677;              // in m^2
    double Aperture=3.6;                    // in meters
    double FocalRatio=8.0;                  // f-number:  f/focalratio
    
    telescopeMount_t TelescopeMount = equatorial;    // supported types: equatorial, altazimuth
    opticalCoating_t OpticalCoating = aluminium;    // supported types: aluminium (Al), silver (Ag),
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"outputTelescopeFile",1, NULL, 'T'},
		{"latitude",1, NULL, 'l'},
		{"longitude",1, NULL, 'g'},
		{"elevation",1, NULL, 'e'},
		{"CollectingArea",1, NULL, 'A'},
		{"Aperture",1, NULL, 'D'},
		{"FocalRatio",1, NULL, 'F'},
		{"TelescopeMount",1, NULL, 'M'},
		{"OpticalCoating",1, NULL, 'C'},
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "T:l:g:e:A:D:F:M:C:p::v::d::t::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
            case 'T':		// outputTelescopeFile
				outputTelescopeFile = optarg;
				break;
                
			case 'l':
				latitude = atof(optarg);
				break;
			case 'g':
				longitude = atof(optarg);
				break;
			case 'e':
				elevation = atof(optarg);
				break;
			case 'A':
				CollectingArea = atof(optarg);
				break;
			case 'D':
				Aperture = atof(optarg);
				break;
			case 'F':
				FocalRatio = atof(optarg);
				break;
			case 'M':
				TelescopeMount = (telescopeMount_t)atoi(optarg);
				break;
    		case 'C':
				OpticalCoating = (opticalCoating_t)atoi(optarg);
				break;
			case 'p':
				plot = 1;
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
		// we need an outputSpectrographFile
		if (outputTelescopeFile.empty()) {
			throw operaException("operaTelescopeSetup: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
        if (verbose) {
			cout << "operaTelescopeSetup: outputTelescopeFile = " << outputTelescopeFile << endl;
			cout << "operaTelescopeSetup: latitude = " << latitude << endl;
			cout << "operaTelescopeSetup: longitude = " << longitude << endl;
			cout << "operaTelescopeSetup: elevation = " << elevation << endl;
			cout << "operaTelescopeSetup: CollectingArea = " << CollectingArea << endl;
			cout << "operaTelescopeSetup: Aperture = " << Aperture << endl;
			cout << "operaTelescopeSetup: FocalRatio = " << FocalRatio << endl;
			cout << "operaTelescopeSetup: TelescopeMount = " << TelescopeMount << endl;
			cout << "operaTelescopeSetup: OpticalCoating = " << OpticalCoating << endl;
        }

        operaInstrumentEnvironmentSetup *instrumentEnvironment = new operaInstrumentEnvironmentSetup();
        operaTelescope *Telescope = instrumentEnvironment->getTelescope();

        Telescope->setlatitude(latitude);
        Telescope->setlongitude(longitude);
        Telescope->setelevation(elevation);
        Telescope->setCollectingArea(CollectingArea);
        Telescope->setAperture(Aperture);
        Telescope->setFocalRatio(FocalRatio);
        Telescope->setTelescopeMount(TelescopeMount);
        Telescope->setOpticalCoating(OpticalCoating);
        
        /*
		 * write out the telescope info
		 */
        instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputTelescopeFile, telescope);

        delete instrumentEnvironment;

	}
	catch (operaException e) {
		cerr << "operaTelescopeSetup: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaTelescopeSetup: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaTelescopeSetup: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --outputTelescopeFile=<TELESCOPE_FILE>"
    " --latitude=<DBL_VALUE>"
    " --longitude=<DBL_VALUE>"
    " --elevation=<DBL_VALUE>"
    " --CollectingArea=<DBL_VALUE>"
    " --Aperture=<DBL_VALUE>"
    " --FocalRatio=<DBL_VALUE>"
    " --TelescopeMount=<telescopeMount_t>"
    " --OpticalCoating=<opticalCoating_t>\n\n"

	" Example: "+string(modulename)+" --outputTelescopeFile=cfht.telescope.gz --latitude=19.8267 --longitude=-155.472 --elevation=4207 --CollectingArea=8.41468 --Aperture=3.6 --FocalRatio=8 --TelescopeMount=1 --OpticalCoating=1 -v \n\n"
	"  -T, --outputTelescopeFile,  Output telescope file\n"
	"  -l, --latitude,  Latitude in deg, (+) North, (-) South\n"
	"  -g, --longitude,  Longitude in deg, +East\n"
	"  -e, --elevation, Elevation in meters\n"
	"  -A, --CollectingArea,  Primary collecting area in m^2\n"
	"  -D, --Aperture,  Telescope primary aperture (diameter) in meters\n"
	"  -F, --FocalRatio,  f-number:  f/focalratio\n"
	"  -M, --TelescopeMount,  supported types: 1=equatorial, 2=altazimuth\n"
	"  -C, --OpticalCoating,  supported types: 1=aluminium (Al), 2=silver (Ag)\n"
	"  -h, --help  display help message\n"
	"  -p, --plot,  Turn on plotting\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n\n";
}
