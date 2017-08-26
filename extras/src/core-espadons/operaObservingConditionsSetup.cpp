/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaObservingConditionsSetup
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

#include "libraries/operaObservingConditions.h"

#include "libraries/operaInstrumentEnvironmentSetup.h"
#include "core-espadons/operaObservingConditionsSetup.h"

/*! \file operaObservingConditionsSetup.cpp */

using namespace std;

/*! 
 * operaObservingConditionsSetup
 * \author Eder Martioli
 * \brief Module to create observing conditions setup class.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
    
    string outputObservingConditionsFile;
	string observercomments;
	string qccomments;
    
    // operaObservingConditions:
    double JDTime = 0.0;                    // Time in Julian days at start of exposure
    double exposureTime = 0.0;              // Exposure time in seconds
    double imageQuality = 0.8;              // seeing in arcsec
    double airmass = 1.3;                   // airmass@zenith = 1
    bool photometric = false;               // based on extinction (true or false)
    double moonphase = 0.5;                 // phase from 0 (newmoon) to 1 (fullmoon)
    double moonalt = 30;                    // in degrees
    double angularDistFromMoon = 60;        // in degrees
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"outputObservingConditionsFile",1, NULL, 'o'},
		{"JDTime",1, NULL, 'J'},
		{"exposureTime",1, NULL, 'e'},
		{"imageQuality",1, NULL, 'Q'},
		{"airmass",1, NULL, 'a'},
		{"photometric",1, NULL, 'P'},
		{"moonphase",1, NULL, 'M'},
		{"moonalt",1, NULL, 'z'},
		{"angularDistFromMoon",1, NULL, 'D'},
		{"observercomments",1, NULL, 'O'},
		{"qccomments",1, NULL, 'q'},
        {"plot",		optional_argument, NULL, 'p'},
        {"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:J:e:Q:a:P:M:z:D:O:q:p::v::d::t::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
			case 'o':		// outputObservingConditionsFile
				outputObservingConditionsFile = optarg;
				break;
			case 'J':
				JDTime = atof(optarg);
				break;
			case 'e':
				exposureTime = atof(optarg);
				break;
			case 'Q':
				imageQuality = atof(optarg);
				break;
			case 'a':
				airmass = atof(optarg);
				break;
			case 'P':
				photometric = (atoi(optarg)?true:false);
				break;
			case 'M':
				moonphase = atof(optarg);
				break;
			case 'z':
				moonalt = atof(optarg);
				break;
			case 'D':
				angularDistFromMoon = atof(optarg);
				break;
			case 'O':
				observercomments = optarg;
				break;
			case 'q':
				qccomments = optarg;
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
		// we need an outputObservingConditionsFile
		if (outputObservingConditionsFile.empty()) {
			throw operaException("operaObservingConditionsSetup: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

        moonphase_t moonphase_obscond = quartermoon;    // 5 phases: newmoon=0, crescentmoon=1, quartermoon=2, gibbousmoon=3, fullmoon=4
        
        if(moonphase >= 0 && moonphase < 0.2) {
            moonphase_obscond = newmoon;
        } else if(moonphase >= 0.2 && moonphase < 0.4) {
            moonphase_obscond = crescentmoon;
        }  else if(moonphase >= 0.4 && moonphase < 0.6) {
            moonphase_obscond = quartermoon;
        }  else if(moonphase >= 0.6 && moonphase < 0.8) {
            moonphase_obscond = gibbousmoon;
        }  else if(moonphase >= 0.8 && moonphase <= 1.0) {
            moonphase_obscond = fullmoon;
        }
        
        double zenithalDistofMoon = 90 - moonalt;         // in degrees
        
        if (verbose) {
			cout << "operaObservingConditionsSetup: outputObservingConditionsFile = " << outputObservingConditionsFile << endl;
			cout << "operaObservingConditionsSetup: JDTime = " << JDTime << endl;
			cout << "operaObservingConditionsSetup: exposureTime = " << exposureTime << endl;
			cout << "operaObservingConditionsSetup: imageQuality = " << imageQuality << endl;
			cout << "operaObservingConditionsSetup: airmass = " << airmass << endl;
			cout << "operaObservingConditionsSetup: photometric = " << photometric << endl;
			cout << "operaObservingConditionsSetup: moonphase = " << moonphase << endl;
			cout << "operaObservingConditionsSetup: zenithalDistofMoon = " << zenithalDistofMoon << endl;
            cout << "operaObservingConditionsSetup: moonalt = " << moonalt << endl;       
			cout << "operaObservingConditionsSetup: angularDistFromMoon = " << angularDistFromMoon << endl;
			cout << "operaObservingConditionsSetup: observercomments = " << observercomments << endl;
			cout << "operaObservingConditionsSetup: qccomments = " << qccomments << endl;
        }
     
        operaInstrumentEnvironmentSetup *instrumentEnvironment = new operaInstrumentEnvironmentSetup();

        operaObservingConditions *ObservingConditions = instrumentEnvironment->getObservingConditions();

        ObservingConditions->setJDTime(JDTime);
        ObservingConditions->setexposureTime(exposureTime);
        ObservingConditions->setimageQuality(imageQuality);
        ObservingConditions->setairmass(airmass);
        ObservingConditions->setphotometric(photometric);
        ObservingConditions->setmoonphase(moonphase_obscond);
        ObservingConditions->setzenithalDistofMoon(zenithalDistofMoon);
        ObservingConditions->setangularDistFromMoon(angularDistFromMoon);
        ObservingConditions->setobservercomments(observercomments);
        ObservingConditions->setqccomments(qccomments);
        
        /*
		 * write out ObservingConditons info
		 */
        instrumentEnvironment->WriteInstrumentEnvironmentSetup(outputObservingConditionsFile, ObsCond);
        
        delete instrumentEnvironment;

	}
	catch (operaException e) {
		cerr << "operaObservingConditionsSetup: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaObservingConditionsSetup: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaObservingConditionsSetup: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --outputObservingConditionsFile=<OBSCOND_FILE>"
    " --JDTime=<DBL_VALUE>"
    " --exposureTime=<DBL_VALUE>"
    " --imageQuality=<DBL_VALUE>"
    " --airmass=<DBL_VALUE>"
    " --photometric=<BOOL_VALUE>"
    " --moonphase=<moonphase>"
    " --moonalt=<DBL_VALUE>"
    " --angularDistFromMoon=<DBL_VALUE>"
    " --observercomments=\"...\"\n\n"
	" Example: "+string(modulename)+" --JDTime=56344.666228 --exposureTime=390 --airmass=1.639 --moonphase=0.88 --moonalt=-29.00 --angularDistFromMoon=116.80 --outputObservingConditionsFile=/data/uhane5/opera//spectra/13AQ02-Feb20/1609811.obscond.gz -v -t \n\n"
	"  -O, --outputObservingConditionsFile,  Output observing conditions file\n"
	"  -J, --JDTime,  Time in Julian days at start of exposure\n"
	"  -e, --exposureTime,  Exposure time in seconds\n"
	"  -Q, --imageQuality,  Seeing in arcsec\n"
	"  -a, --airmass,  Airmass@zenith == 1\n"
	"  -P, --photometric,  Photometric flag based on extinction < 0.5 mag (true=1 or false=0)\n"
	"  -M, --moonphase, Phase defined only between 0 and 1, where newmoon=0 and fullmoon=1\n"
	"  -z, --moonalt, Altitude of Moon in degrees\n"
	"  -D, --angularDistFromMoon,  Angular distance from Moon to target in degrees\n"
	"  -h, --help  display help message\n"
	"  -p, --plot,  Turn on plotting\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n\n";
}
