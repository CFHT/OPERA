/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaObjectInTheSkySetup
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
#include "libraries/operaObjectInTheSky.h"

#include "core-espadons/operaObjectInTheSkySetup.h"

/*! \file operaObjectInTheSkySetup.cpp */

using namespace std;

/*! 
 * operaObjectInTheSkySetup
 * \author Eder Martioli
 * \brief  Create observing conditions and objectInTheSky setup classes for photometric standards.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
    
    string outputObjectInTheSkyFile;
    // operaObjectInTheSky:
    string sourceID("Vega (Default)");                                       // source identifier: name, catalog number, etc
    double RA=18.61564889;                                          // in degrees
    double Dec=+38.78368889;                                        // in degrees
    struct ProperMotion {                                           // in mas/yr
		double muRA, muDec;
	} ProperMotion = {200.94, 286.23};
    double Parallax = 130.23;                                       // in mas
    double V_magnitude = VEGA_VBAND_MAGNITUDE;                      // in Vega magnitude in V-band (Johnson)
    double EffectiveTemperature = VEGA_EFF_TEMPERATURE_K;           // in Kelvin
    double RadialVelocity = -13.9;                                  // in km/s
    operaSpectralType_t SpectralType=A_type;                        // main sequence (V): O,B,A,F,G,K,M,L,T,Y
    
    int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"outputObjectInTheSkyFile",1, NULL, 'B'},
		{"sourceID",1, NULL, 'I'},
		{"RA",1, NULL, 'R'},
		{"Dec",1, NULL, 'D'},
		{"ProperMotion",1, NULL, 'm'},
 		{"Parallax",1, NULL, 'P'},
		{"V_magnitude",1, NULL, 'V'},
		{"EffectiveTemperature",1, NULL, 'T'},
		{"RadialVelocity",1, NULL, 'r'},
		{"SpectralType",1, NULL, 'S'},
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "B:I:R:D:m:P:V:T:r:S:p::v::d::t::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
			case 'B':		// outputObjectInTheSkyFile
				outputObjectInTheSkyFile = optarg;
				break;
			case 'I':
				sourceID = optarg;
				break;
			case 'R':
				RA = atof(optarg);
				break;
			case 'D':
				Dec = atof(optarg);
				break;
			case 'm':
				if (strlen(optarg))
					sscanf(optarg, "%lf %lf", &ProperMotion.muRA, &ProperMotion.muDec);
				break;
			case 'P':
				Parallax = atof(optarg);
				break;
			case 'V':
				V_magnitude = atof(optarg);
				break;
			case 'T':
				EffectiveTemperature = atof(optarg);
				break;
			case 'r':
				RadialVelocity = atof(optarg);
				break;
			case 'S':
				SpectralType = (operaSpectralType_t)atoi(optarg);
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
        
		// we need an outputObjectInTheSkyFile
		if (outputObjectInTheSkyFile.empty()) {
			throw operaException("operaObjectInTheSkySetup: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (verbose) {
			cout << "operaObjectInTheSkySetup: outputObjectInTheSkyFile = " << outputObjectInTheSkyFile << endl;
			cout << "operaObjectInTheSkySetup: sourceID = " << sourceID << endl;
			cout << "operaObjectInTheSkySetup: RA = " << RA << endl;
			cout << "operaObjectInTheSkySetup: Dec = " << Dec << endl;
			cout << "operaObjectInTheSkySetup: ProperMotion \"muRA muDec\" = \"" << ProperMotion.muRA << " " << ProperMotion.muDec << "\"" << endl;
            cout << "operaObjectInTheSkySetup: Parallax = " << Parallax << endl;
			cout << "operaObjectInTheSkySetup: V_magnitude = " << V_magnitude << endl;
			cout << "operaObjectInTheSkySetup: EffectiveTemperature = " << EffectiveTemperature << endl;
			cout << "operaObjectInTheSkySetup: RadialVelocity = " << RadialVelocity << endl;
			cout << "operaObjectInTheSkySetup: SpectralType = " << SpectralType << endl;
        }

        operaInstrumentEnvironmentSetup *instrumentEnvironment = new operaInstrumentEnvironmentSetup();
        operaObjectInTheSky *ObjectInTheSky = instrumentEnvironment->getObjectInTheSky();

        ObjectInTheSky->setsourceID(string(sourceID.c_str()));
        ObjectInTheSky->setRA(RA);
        ObjectInTheSky->setDec(Dec);
        ObjectInTheSky->setProperMotionRA(ProperMotion.muRA);
        ObjectInTheSky->setProperMotionDec(ProperMotion.muDec);
        ObjectInTheSky->setParallax(Parallax);
        ObjectInTheSky->setV_magnitude(V_magnitude);
        ObjectInTheSky->setEffectiveTemperature(EffectiveTemperature);
        ObjectInTheSky->setRadialVelocity(RadialVelocity);
        ObjectInTheSky->setSpectralType(SpectralType);
        
        /*
		 * write out ObjectInTheSky info
		 */
        instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputObjectInTheSkyFile, SkyObj);
        
        delete instrumentEnvironment;

	}
	catch (operaException e) {
		cerr << "operaObjectInTheSkySetup: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaObjectInTheSkySetup: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaObjectInTheSkySetup: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --outputObjectInTheSkyFile=<SKYOBJ_FILE>"
    " --sourceID=<STRING>"
    " --RA=<DBL_VALUE>"
    " --Dec=<DBL_VALUE>"
    " --ProperMotion=<DBL_VALUE DBL_VALUE>"
    " --Parallax=<DBL_VALUE>"
    " --V_magnitude=<DBL_VALUE>"
    " --EffectiveTemperature=<DBL_VALUE>"
    " --RadialVelocity=<DBL_VALUE>"
    " --SpectralType=<operaSpectralType_t>\n\n"
	" Example: "+string(modulename)+" -v \n\n"
	"  -B, --outputObjectInTheSkyFile,  Output sky object file\n"
	"  -I, --sourceID,  Source identifier: name, catalog number, etc\n"
	"  -R, --RA,  Right ascension in degrees\n"
	"  -D, --Dec,  Declination in degrees\n"
	"  -m, --ProperMotion,  Proper motion (RA and Dec) in mas/yr\n"
	"  -p, --Parallax,  Parallax in mas \n"
	"  -V, --V_magnitude,  V-band magnitude in Vega mag (Johnson)\n"
	"  -T, --EffectiveTemperature,  Effective black body temperature in Kelvin\n"
	"  -r, --RadialVelocity,  Radial velocity in km/s\n"
	"  -S, --SpectralType,  Spectral type. Support the following\n"  
	"                       types of main sequence stars:\n"
	"                                      O_type=1,\n"
	"                                      B_type=2,\n"
	"                                      A_type=3,\n"
	"                                      F_type=4,\n"
	"                                      G_type=5,\n"
	"                                      K_type=6,\n"
	"                                      M_type=7,\n"
	"                                      L_type=8,\n"
	"                                      T_type=9,\n"
	"                                      Y_type=10,\n"
	"  -h, --help  display help message\n"
	"  -p, --plot,  Turn on plotting\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n\n";
}
