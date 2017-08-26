/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirPickSkies
 Version: 1.0
 Description: Module to Pick skies from a natural sequence, creating a sky list
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
#include "core-wircam/wirPickSkies.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"
#include "libraries/operaWIRCamImage.h"
#include "libraries/operaHelio.h"

/* \file wirPickSkies.cpp */
/* \package core_wircam */

using namespace std;

/*
 * wirPickSkies
 * \author Doug Teeple
 * \brief Module to Pick skies from a natural sequence, creating a sky list.
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
 * observingstrategy
 * A string, either DP, NDP, WDP, ST or STDP that applies to all the odometers in
 * list_naturalsequence.
 * maxskies
 * The maximum number of images that can be associated to an odometer for sky construction.
 * rejectcurrentimage
 * If set, the current image is rejected from the association. Set to 1 by default.
 * rejectcurrentdp
 * If set, all images at the current DP position are rejected from the
 * association. Set to 1 by default. Example: if DP6 and rejectcurrentdp=1 then the pool
 * DPs to use is 5.
 * mindp
 * The hard minimum number of different dither positions to flag the association as sucessfull.
 * For example, a DP3 has 3 different dither positions (even if it's repeated 100 times).
 * Setting mindp to 4 would make the association fail, at 2 it will work but
 * sky construction will have only 2 sets of positions for medianing. Usually, at least 3 is
 * required, better 5 or more.
 * dpradius
 * Defines the minimum sky angular distance for two images to be considered as separate dither
 * positions. In arcseconds with default=15. Anything below about 5 arcseconds becomes
 * suspicious because wings of a PSF may still cross and leave some sky pixels undefined in
 * the sky construction. Use dpradius<5 for pathetic cases (bad observing strategy).
 * equalsampleperdp
 * If set (default) then the same number of images must be used for each acceptable DP.
 * For example, if it's a DP7 repeated twice, the parameter imposes that you either include
 * a single loop of DP7 or the two of them, but you could not include the first loop plus
 * a couple of the 2nd loop. The parameter does not impose that all of the DP7 be used but
 * rather that those DP included all be included with equal amount so no bias in medianing
 * the images be introduced.
 * maxtime
 * In minutes, set the maximum time difference between the current image and an image to
 * associate. A good default number is 15 minutes. A float.
 */

#define MAXSKIES 50

int main(int argc, char *argv[])
{
	int opt;
	string name_skies[MAXSKIES];
	unsigned int sky_grades[MAXSKIES];
	string name_odometer;
	string name_skylist;
	string name_defaultsky;
	string name_badpix;
	string observingstrategy;
	unsigned sky_count = 0;
	unsigned grade_count = 0;
	string iiwiversion = "3.0";
	string procdate;
	bool rejectcurrentimage = true;
	bool rejectcurrentdp = true;
	bool equalsampleperdp = true;
	float dpradius = 15.0;
	float maxtime = 15.0;
	unsigned int grade = 1;
	unsigned int mindp = 2;
	unsigned int maxskies = MAXSKIES;
	
	vector<string> skylistvector;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
	
	struct option longopts[] = {
		{"odometer",			1, NULL, '1'},
		{"name_skies",			1, NULL, 'i'},
		{"sky_grades",			1, NULL, 'k'},
		{"name_badpix",			1, NULL, 'B'},
		{"name_skylist",		1, NULL, 'o'},
		{"name_defaultsky",		1, NULL, 'u'},
		{"observingstrategy",	1, NULL, 's'},
		{"maxskies",			1, NULL, 'm'},
		{"rejectcurrentimage",	1, NULL, 'c'},
		{"rejectcurrentdp",		1, NULL, 'j'},
		{"equalsampleperdp",	1, NULL, 'e'},
		{"dpradius",			1, NULL, 'r'},
		{"maxtime",				1, NULL, 'x'},
		{"mindp",				1, NULL, 'y'},
		{"grade",				1, NULL, 'g'},
		
		{"plotfilename",		1, NULL, 'P'},
		{"datafilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "1:i:o:k:u:s:m:c:j:e:r:x:y:B:g:v::d::t::p::h",
							 longopts, NULL))  != -1) {
		switch(opt) {
            case '1':
                name_odometer = optarg;
            break;
            case 'i':
                name_skies[sky_count++] = optarg;
                if (sky_count > MAXSKIES) {
                    throw operaException("wirPickSkies: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
            break;
            case 'k':
                sky_grades[grade_count++] = atoi(optarg);
                if (grade_count > MAXSKIES) {
                    throw operaException("wirPickSkies: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
            break;
            case 'o':
                name_skylist = optarg;
            break;
            case 'u':
                name_defaultsky = optarg;
            break;
            case 'B':
                name_badpix = optarg;
            break;
            case 's':
                observingstrategy = optarg;
            break;
            case 'm':
                maxskies = atoi(optarg);
                if (maxskies > MAXSKIES) {
                    throw operaException("wirPickSkies: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
            break;
            case 'c':
                rejectcurrentimage = atoi(optarg) == 1;
            break;
            case 'j':
                rejectcurrentdp = atoi(optarg) == 1;
            break;
            case 'e':
                equalsampleperdp = atoi(optarg) == 1;
            break;
            case 'r':
				dpradius = atof(optarg);
            break;
            case 'x':
                maxtime = atof(optarg);
            break;
            case 'y':
                mindp = atoi(optarg);
            break;
            case 'g':
                grade = atoi(optarg);
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
		if (name_odometer.empty()) {
			throw operaException("wirPickSkies: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (name_skylist.empty()) {
			throw operaException("wirPickSkies: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (sky_count == 0) {
			throw operaException("wirPickSkies: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (sky_count < maxskies) {
			maxskies = sky_count;
		}
		
		if (verbose) {
			for (unsigned i=0; i<sky_count; i++) {
				cout << "wirPickSkies: sky image[" << itos(i) << "] = " << name_skies[i] << endl;
			}
			cout << "wirPickSkies: name_odometer = " << name_odometer << endl;
			cout << "wirPickSkies: name_skylist = " << name_skylist << endl;
			cout << "wirPickSkies: name_defaultsky = " << name_defaultsky << endl;
			cout << "wirPickSkies: badpix = " << name_badpix << endl;
			cout << "wirPickSkies: observingstrategy = " << observingstrategy << endl;
			cout << "wirPickSkies: maxskies = " << maxskies << endl;
			cout << "wirPickSkies: rejectcurrentimage = " << rejectcurrentimage << endl;
			cout << "wirPickSkies: rejectcurrentdp = " << rejectcurrentdp << endl;
			cout << "wirPickSkies: equalsampleperdp = " << equalsampleperdp << endl;
			cout << "wirPickSkies: dpradius = " << dpradius << endl;
			cout << "wirPickSkies: maxtime = " << maxtime << endl;
			cout << "wirPickSkies: sky_count = " << sky_count << endl;
			cout << "wirPickSkies: grade = " << grade << endl;
			for (unsigned i=0; i<sky_count; i++) {
				cout << "wirPickSkies: name_skies[" << i << "] = " << name_skies[i] << endl;
				cout << "wirPickSkies: sky_grades[" << i << "] = " << sky_grades[i] << endl;
			}
            if (plot) {
                cerr << "wirPickSkies: plotfilename = " << plotfilename << endl;
                cerr << "wirPickSkies: datafilename = " << datafilename << endl;
                cerr << "wirPickSkies: scriptfilename = " << scriptfilename << endl;
            }
		}
		
        ofstream *fdata = NULL;
        ofstream skylist(name_skylist.c_str());
		operaWIRCamImage odometer(name_odometer, tfloat, READONLY);
		float odometer_juliandate = odometer.operaFITSGetFloatHeaderValue("MJDATE") * 24.0 * 60.0;
		string Object_RA = odometer.operaFITSGetHeaderValue("OBJRA");		//'17:40:57.97'
		string Object_DEC = odometer.operaFITSGetHeaderValue("OBJDEC");
		double odometer_ra = 15.0 * sexigesimal_str_to_dec(Object_RA.c_str());
		double odometer_dec = sexigesimal_str_to_dec(Object_DEC.c_str());
		double WIRCamFieldOfView = 20.0 * 60.0;
		int k = -1, j = -1;
		
		// locate the anchor point, which is this odometer
		for (unsigned i=0; i<sky_count; i++) {
			if (name_odometer == name_skies[i]) {
				k = j = i;
			}
		}
		if (k < 0) {
			throw operaException("wirPickSkies: ", operaErrorNoAnchor, __FILE__, __FUNCTION__, __LINE__);
		}
		// index outward from the anchor point
		unsigned i = 0;
		for (unsigned n = 0; n < sky_count; n++) {
			if (i%2)
				i = k--;
			else
				i = j++;
			if (i >= 0 && i < sky_count) {
                if (rejectcurrentimage && name_skies[i] == name_odometer) {
                    // exclude self, i.e. do nothing
                } else { // include self
                    operaWIRCamImage sky(name_skies[i], tfloat, READONLY);
                    float sky_juliandate = sky.operaFITSGetFloatHeaderValue("MJDATE") * 24.0 * 60.0;
                    string Object_RA = sky.operaFITSGetHeaderValue("OBJRA");		//'17:40:57.97'
                    string Object_DEC = sky.operaFITSGetHeaderValue("OBJDEC");
					// The target type can be:
					// NDP:			"TARGET" or "SKY"
					// DP/WDP/UDP:	"TARGET+SKY"		i.e. use surrounding images
					// ST:			"TARGET" or "SKY"
                    string TargetType = sky.operaFITSGetHeaderValue("TRGTYPE");
                    double sky_ra = 15.0 * sexigesimal_str_to_dec(Object_RA.c_str());
                    double sky_dec = sexigesimal_str_to_dec(Object_DEC.c_str());
                    double timedifference = fabs(sky_juliandate - odometer_juliandate);
                    double calculated_dpradius = 3600.0 * sqrt((odometer_ra - sky_ra)*(odometer_ra - sky_ra)+(odometer_dec - sky_dec)*(odometer_dec - sky_dec));
                    if (observingstrategy == "DP") {
                        // 1. the time is less the given time
                        // 2. Grade is better or equal
                        // 3. DP is larger than the dpradius
                        // 4. DP shoule be less than a WIRCamFieldOfView
                        if (timedifference <= maxtime &&
                            sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius &&
                            calculated_dpradius <= WIRCamFieldOfView) {
							skylistvector.push_back(name_skies[i]);
                        }
                    } else if (observingstrategy == "WDP") {
                        // 1. the time is less the given time
                        // 2. Grade is better or equal
                        // 3. DP is larger than the dpradius
                        if (timedifference <= maxtime &&
                            sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius) {
                            skylistvector.push_back(name_skies[i]);
                        }
                    } else if (observingstrategy == "NDP") {
                        // 1. the time is less the given time
                        // 2. Grade is better or equal
                        // 3. DP radius is greater than 20
                        if (timedifference <= maxtime &&
                            sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius &&
							TargetType == "SKY") {
							skylistvector.push_back(name_skies[i]);
                        }
                    } else if (observingstrategy == "UDP") {
                        // 1. the time is less the given time
                        // 2. Grade is better or equal
                        // 3. DP is larger than the dpradius
                        if (timedifference <= maxtime &&
                            sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius) {
							skylistvector.push_back(name_skies[i]);
                        }
                    } else if (observingstrategy == "ST") {
                        // 1. Grade is better or equal
                        // 2. DP is larger than the dpradius
                        if (sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius && 
							TargetType == "SKY") {
                            skylistvector.push_back(name_skies[i]);
                        }
                    } else if (observingstrategy == "P22") {
                        // 1. the time is less the given time
                        // 2. Grade is better or equal
                        // 3. DP is larger than the dpradius
                        if (timedifference <= maxtime &&
                            sky_grades[i] <= grade &&
                            calculated_dpradius > dpradius) {
							skylistvector.push_back(name_skies[i]);
                        }
                    }
                    sky.operaFITSImageClose();
                }
			}
		}
		if (skylistvector.size() < mindp) {
            skylistvector.empty();                      // remove what is there
			skylistvector.push_back(name_defaultsky);   // and use the default sky
		}
		// now select the images closes in time up to the maxskies limit
		// and write out to the skylist
		for (unsigned i=0; i < skylistvector.size() && i < maxskies; i++) {
			skylist << skylistvector[i] << endl;
		}
        skylist.close();
		odometer.operaFITSImageClose();
		
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
        }
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateSkyPickPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }
	}
	catch (operaException e) {
		cerr << "wirPickSkies: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirPickSkies: " << operaStrError(errno) << endl;
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

void GenerateSkyPickPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
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
