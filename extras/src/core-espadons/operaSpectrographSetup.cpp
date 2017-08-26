/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaSpectrographSetup
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
#include "libraries/operaSpectrograph.h"

#include "core-espadons/operaSpectrographSetup.h"

/*! \file operaSpectrographSetup.cpp */

using namespace std;

/*! 
 * operaSpectrographSetup
 * \author Eder Martioli
 * \brief Module to create Spectrograph setup class.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
    
    string outputSpectrographFile;
    
    double InjectionHoleDiameter = 1.6;               // in arcsec
    opticalFiber_t  OpticalFiber = STUPolymicro;               // Supported types: FBPPolymicro=1, STUPolymicro=2
    double fiberLength = 33;                         // in meters
    double fiberCoreDiameter = 100;                   // in microns
    unsigned numberOfInputFibers = 2;               // 1 for star-only, 2 for polar or star+sky
    unsigned numberOfSlices = 3;                    // Espadons currently only uses 3 or 6 slices
    double spectralResolution = 65000;                  // R=lambda/dlambda
    spectrographCCD_t SpectrographCCD = Olapa;                  // Supported types: Olapa, EEV1
    EspadonsCCDReadoutSpeed_t EspadonsCCDReadoutSpeed = normalmode;  // Supported modes: slowmode, normalmode, fastmode
    EspadonsInstrumentMode_t EspadonsInstrumentMode = polarimetric;    // polarimetric, staronly, starplussky, GRACES_staronly, GRACES_starplussky
   	struct pixelsize {
		double x, y;
	} pixelsize = {13, 13};
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"outputSpectrographFile",      1, NULL, 'S'},
        {"InjectionHoleDiameter",       1, NULL, 'I'},
        {"OpticalFiber",                1, NULL, 'F'},
        {"fiberLength",                 1, NULL, 'L'},
        {"fiberCoreDiameter",           1, NULL, 'D'},
        {"numberOfInputFibers",         1, NULL, 'n'},
        {"numberOfSlices",              1, NULL, 's'},
        {"spectralResolution",          1, NULL, 'R'},
        {"SpectrographCCD",             1, NULL, 'C'},
        {"EspadonsCCDReadoutSpeed",     1, NULL, 'r'},
        {"EspadonsInstrumentMode",      1, NULL, 'm'},
		{"pixelsize",                   1, NULL, 'P'},
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "S:I:F:L:D:n:s:R:C:r:m:P:p::v::d::t::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
			case 'S':		// outputSpectrographFile
				outputSpectrographFile = optarg;
				break;
			case 'I':
				InjectionHoleDiameter = atof(optarg);
				break;
			case 'F':
				OpticalFiber = (opticalFiber_t)atoi(optarg);
				break;
			case 'L':
				fiberLength = atof(optarg);
				break;
			case 'D':
				fiberCoreDiameter = atof(optarg);
				break;
			case 'n':
				numberOfInputFibers = atoi(optarg);
				break;
			case 's':
				numberOfSlices = atoi(optarg);
				break;
			case 'R':
				spectralResolution = atof(optarg);
				break;
			case 'C':
				SpectrographCCD = (spectrographCCD_t)atoi(optarg);
				break;
			case 'r':
				EspadonsCCDReadoutSpeed = (EspadonsCCDReadoutSpeed_t)atoi(optarg);
				break;
			case 'm':
				EspadonsInstrumentMode = (EspadonsInstrumentMode_t)atoi(optarg);
				break;
			case 'P':		// pixel size in microns
				if (strlen(optarg))
					sscanf(optarg, "%lf %lf", &pixelsize.x, &pixelsize.y);
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
		if (outputSpectrographFile.empty()) {
			throw operaException("operaSpectrographSetup: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
        if (verbose) {
			cout << "operaSpectrographSetup: outputSpectrographFile = " << outputSpectrographFile << endl;
			cout << "operaSpectrographSetup: InjectionHoleDiameter = " << InjectionHoleDiameter << endl;
			cout << "operaSpectrographSetup: OpticalFiber = " << OpticalFiber << endl;
			cout << "operaSpectrographSetup: fiberLength = " << fiberLength << endl;
			cout << "operaSpectrographSetup: fiberCoreDiameter = " << fiberCoreDiameter << endl;
			cout << "operaSpectrographSetup: numberOfInputFibers = " << numberOfInputFibers << endl;
			cout << "operaSpectrographSetup: numberOfSlices = " << numberOfSlices << endl;
			cout << "operaSpectrographSetup: spectralResolution = " << spectralResolution << endl;
			cout << "operaSpectrographSetup: SpectrographCCD = " << SpectrographCCD << endl;
			cout << "operaSpectrographSetup: EspadonsCCDReadoutSpeed = " << EspadonsCCDReadoutSpeed << endl;
			cout << "operaSpectrographSetup: EspadonsInstrumentMode = " << EspadonsInstrumentMode << endl;
			cout << "operaSpectrographSetup: pixelsize {x,y} = {" << pixelsize.x << "," << pixelsize.y << "}" << endl;
        }

        operaInstrumentEnvironmentSetup *instrumentEnvironment = new operaInstrumentEnvironmentSetup();

        operaSpectrograph *Spectrograph = instrumentEnvironment->getSpectrograph();
                
        
        Spectrograph->setInjectionHoleDiameter(InjectionHoleDiameter);

        Spectrograph->setOpticalFiber(OpticalFiber);

        Spectrograph->setfiberLength(fiberLength);

        Spectrograph->setfiberCoreDiameter(fiberCoreDiameter);

        Spectrograph->setnumberOfInputFibers(numberOfInputFibers);

        Spectrograph->setnumberOfSlices(numberOfSlices);

        Spectrograph->setspectralResolution(spectralResolution);

        Spectrograph->setSpectrographCCD(SpectrographCCD);

        Spectrograph->setEspadonsCCDReadoutSpeed(EspadonsCCDReadoutSpeed);

        Spectrograph->setEspadonsInstrumentMode(EspadonsInstrumentMode);

        Spectrograph->setx_pixelsize(pixelsize.x);
        
        Spectrograph->sety_pixelsize(pixelsize.y);
        
        /*
		 * write out Spectrograph info
		 */

        instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputSpectrographFile, spectrograph);

        delete instrumentEnvironment;

	}
	catch (operaException e) {
		cerr << "operaSpectrographSetup: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaSpectrographSetup: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaSpectrographSetup: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --outputSpectrographFile=<SPECTROGRAPH_FILE>"
    " --InjectionHoleDiameter=<DBL_VALUE>"
    " --OpticalFiber=<opticalFiber_t>"
    " --fiberLength=<DBL_VALUE>"
    " --fiberCoreDiameter=<DBL_VALUE>"
    " --numberOfInputFibers=<UNS_VALUE>"
    " --numberOfSlices=<UNS_VALUE>"
    " --spectralResolution=<DBL_VALUE>"
    " --SpectrographCCD=<spectrographCCD_t>"
    " --EspadonsCCDReadoutSpeed=<EspadonsCCDReadoutSpeed_t>"
    " --EspadonsInstrumentMode=<EspadonsInstrumentMode_t>"
    " --pixelsize=<DBL_VALUE,DBL_VALUE>\n\n"
	" Example: "+string(modulename)+" --outputSpectrographFile=Espadons_pol_normal.spectrograph.gz --InjectionHoleDiameter=1.6 --OpticalFiber=2 --fiberLength=33 --fiberCoreDiameter=100 --numberOfInputFibers=2 --numberOfSlices=3 --spectralResolution=65000 --SpectrographCCD=1 --EspadonsCCDReadoutSpeed=2 --EspadonsInstrumentMode=1 --pixelsize=\"13 13\" -v  \n\n"
	"  -S, --outputSpectrographFile,  Output spectrograph file\n"
	"  -I, --InjectionHoleDiameter, Diameter of injection hole in arcsec\n"
	"  -F, --OpticalFiber, Optical Fiber type. Supported types: 1=FBPPolymicro, 2=STUPolymicro\n"
	"  -L, --fiberLength, Fiber length in meters\n"
	"  -D, --fiberCoreDiameter, Fiber core diameter in microns\n"
	"  -n, --numberOfInputFibers, Number of input fibers: 1 for star-only, 2 for polar or star+sky\n"
	"  -s, --numberOfSlices, Number of slices. Espadons currently only uses 3 or 6 slices. \n"
	"  -R, --spectralResolution, Nominal spectral resolution R=lambda/dlambda\n"
	"  -C, --SpectrographCCD, CCD detector. Supported types: 1=Olapa, 2=EEV1\n"
	"  -r, --EspadonsCCDReadoutSpeed, CCD readout speed. Supported modes: 1=slowmode, 2=normalmode, 3=fastmode\n"
	"  -m, --EspadonsInstrumentMode, Supported modes: 1=polarimetric, 2=staronly, 3=starplussky, 4=GRACES_staronly, 5=GRACES_starplussky\n"
	"  -P, --pixelsize, \"<xsize> <ysize>\" in microns\n"
	"  -h, --help  display help message\n"
	"  -p, --plot,  Turn on plotting\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n\n";
}
