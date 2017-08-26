/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirPickTwilightFlats
 Version: 1.0
 Description: Module to Pick twilight flats from a natural sequence, creating a twilight flat list
 to start up with an OPERA module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Jul/2013
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2013  Opera Pipeline team, Canada France Hawaii Telescope
 
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
#include <vector>

#include "operaError.h"
#include "core-wircam/wirPickTwilightFlats.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"
#include "libraries/operaWIRCamImage.h"
#include "libraries/operaHelio.h"

/* \file wirPickTwilightFlats.cpp */
/* \package core_wircam */

using namespace std;

/*
 * wirPickTwilightFlats
 * \author Doug Teeple
 * \brief Module to Pick twilight flats from a natural sequence, creating a twilight flat list.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core_wircam
 */
/*
 */

#define MAXTWILIGHTFLATS 50

int main(int argc, char *argv[])
{
	int opt;
	string name_twilightflats[MAXTWILIGHTFLATS];
	string name_flatlist;
	string filter;
	unsigned flat_count = 0;
	string iiwiversion = "3.0";
	string procdate;
	unsigned int maxtwilightflats = MAXTWILIGHTFLATS;
	
	unsigned int minnum = 0;
	float minadu = 0;
	float maxadu = 0;
	float maxtime = 0;
	float minflux = 0;
	float minslope = 0;
	
	vector<string> flatlistvector;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
	
	struct option longopts[] = {
		{"name_twilightflats",	1, NULL, 'i'},
		{"name_flatlist",		1, NULL, 'o'},
		{"maxtwilightflats",	1, NULL, 'm'},
		// flat parameters
		{"filter",				1, NULL, 'r'},
		{"minadu",				1, NULL, 'a'},
		{"maxadu",				1, NULL, 'x'},
		{"maxtime",				1, NULL, 'e'},
		{"minnum",				1, NULL, 'l'},
		{"minflux",				1, NULL, 'f'},
		{"minslope",			1, NULL, 's'},
		
		{"plotfilename",		1, NULL, 'P'},
		{"datafilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:m:r:a:x:e:l:f:s:P:F:S:v::d::t::p::h",
							 longopts, NULL))  != -1) {
		switch(opt) {
             case 'i':
                name_twilightflats[flat_count++] = optarg;
                if (flat_count > MAXTWILIGHTFLATS) {
                    throw operaException("wirPickTwilightFlats: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
            break;
            case 'o':
                name_flatlist = optarg;
            break;
            case 'm':
                maxtwilightflats = atoi(optarg);
                if (maxtwilightflats > MAXTWILIGHTFLATS) {
                    throw operaException("wirPickTwilightFlats: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
				break;
            case 'r':
                filter = optarg;
 				break;
            case 's':
                minslope = atof(optarg);
 				break;
            case 'a':
                minadu = atof(optarg);
 				break;
            case 'x':
                maxadu = atof(optarg);
 				break;
            case 'e':
                maxtime = atof(optarg);
 				break;
            case 'l':
                minnum = atoi(optarg);
 				break;
            case 'f':
                minflux = atof(optarg);
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
		if (name_flatlist.empty()) {
			throw operaException("wirPickTwilightFlats: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat_count == 0) {
			throw operaException("wirPickTwilightFlats: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat_count < MAXTWILIGHTFLATS) {
			maxtwilightflats = flat_count;
		}
		
		if (verbose) {
			for (unsigned i=0; i<flat_count; i++) {
				cout << "wirPickTwilightFlats: sky image[" << itos(i) << "] = " << name_twilightflats[i] << endl;
			}
			cout << "wirPickTwilightFlats: name_flatlist = " << name_flatlist << endl;
			for (unsigned i=0; i<flat_count; i++) {
				cout << "wirPickTwilightFlats: name_twilightflats[" << i << "] = " << name_twilightflats[i] << endl;
			}
			cout << "wirPickTwilightFlats: maxtwilightflats = " << maxtwilightflats << endl;
			cout << "wirPickTwilightFlats: minslope = " << minslope << endl;
			cout << "wirPickTwilightFlats: minadu = " << minadu << endl;
			cout << "wirPickTwilightFlats: maxadu = " << maxadu << endl;
			cout << "wirPickTwilightFlats: maxtime = " << maxtime << endl;
			cout << "wirPickTwilightFlats: minnum = " << minnum << endl;
			cout << "wirPickTwilightFlats: minflux = " << minflux << endl;
            if (plot) {
                cout << "wirPickTwilightFlats: plotfilename = " << plotfilename << endl;
                cout << "wirPickTwilightFlats: datafilename = " << datafilename << endl;
                cout << "wirPickTwilightFlats: scriptfilename = " << scriptfilename << endl;
            }
		}
		
        ofstream *fdata = NULL;
        ofstream flatlist(name_flatlist.c_str());
		operaWIRCamImage sky(name_twilightflats[0], tfloat, READONLY);
		string last_filter = sky.operaFITSGetHeaderValue("FILTER");
		float last_juliandate = sky.operaFITSGetFloatHeaderValue("MJDATE") * 24.0 * 60.0;
		float last_skylevel = sky.getAverageSkyLevel();
		float totaladu = last_skylevel;
		if (last_skylevel > minadu && filter == last_filter) {
			flatlistvector.push_back(name_twilightflats[0]);
			if (verbose) {
				cout << "wirPickTwilightFlats: adding " << name_twilightflats[0] << " filter " << last_filter << " JD " << last_juliandate << " SKY LEVEL " << last_skylevel << " total ADU " << totaladu << endl;
			}
		}
		for (unsigned i = 1; i < flat_count; i++) {
			operaWIRCamImage sky(name_twilightflats[i], tfloat, READONLY);
			string newfilter = sky.operaFITSGetHeaderValue("FILTER");
			float juliandate = sky.operaFITSGetFloatHeaderValue("MJDATE") * 24.0 * 60.0;
			float skylevel = sky.getAverageSkyLevel();
			float timedifference = fabs(juliandate - last_juliandate);
			float slope = fabs(last_skylevel - skylevel) / timedifference;
			if (skylevel > minadu && timedifference < maxtime && slope > minslope && filter == newfilter) {
				flatlistvector.push_back(name_twilightflats[i]);
				if (verbose) {
					cout << "wirPickTwilightFlats: adding " << name_twilightflats[i] << " filter " << last_filter << " JD " << last_juliandate << " SKY LEVEL " << last_skylevel << " total ADU " << totaladu << endl;
				}
			}
			last_skylevel = skylevel;
			last_juliandate = juliandate;
			totaladu += last_skylevel;
			sky.operaFITSImageClose();
 		}
		if (flatlistvector.size() < minnum || totaladu < minflux) {
            flatlistvector.empty();                      // remove what is there
			if (verbose) {
				cout << "wirPickTwilightFlats: requirements not satisfied, removing sequence... size = " << flatlistvector.size() << " minnum " << minnum << " totaladu " << totaladu << " minflux " << minflux << endl;
			}
		}
		// now select the images closes in time up to the MAXTWILIGHTFLATS limit
		// and write out to the flatlist
		for (unsigned i=0; i < flatlistvector.size(); i++) {
			flatlist << flatlistvector[i] << endl;
		}
        flatlist.close();
		
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
        }
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateTwilightFlatPickPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }
	}
	catch (operaException e) {
		cerr << "wirPickTwilightFlats: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirPickTwilightFlats: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";		
}

void GenerateTwilightFlatPickPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName,"w");
    
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"set view 0,0\n");
    fprintf(fgnu,"set iso 100\n");
    fprintf(fgnu,"set samples 100\n");
    fprintf(fgnu,"set pm3d at s\n");
    fprintf(fgnu,"set ticslevel 0\n");   
    
    fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
    fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
	
    if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
		fprintf(fgnu,"replot\n");       
		fclose(fgnu);   
		systemf("gnuplot -persist %s",gnuScriptFileName);
    } else {
		fclose(fgnu);  
	}
}
