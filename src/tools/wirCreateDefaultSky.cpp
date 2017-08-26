/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: wirCreateDefaultSky
 Version: 1.0
 Description: Create a default sky image for WIRCam composed of all 0's.
 Author(s): Megan Tannock
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: May/2013
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

#include <string.h>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaWIRCamImage.h"

/*! \file wirCreateDefaultSky.cpp */

using namespace std;

/*! 
 * wirCreateDefaultSky
 * \author Megan Tannock
 * \brief Create a default sky image for WIRCam composed of all 0's.
 * \return EXIT_STATUS
 */

/* Print out the proper program usage syntax */
void printUsageSyntax(char *whoami) {
	
	cerr << " Usage: " << whoami << " -[dvth] \n";
}

int main(int argc, char *argv[])
{        
	int opt;
	
	int plot=0, verbose=0, debug=0, trace=0;
	
	string defaultskydefaultskyname = "defaultsky.fits";
    
	struct option longopts[] = {
		
		{"plot",		optional_argument, NULL, 'p'},
        {"display",		optional_argument, NULL, 'D'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}
    };
	
	try {
        while((opt = getopt_long(argc, argv, "p::v::d::t::h", longopts, NULL))  != -1)
        {
            switch(opt)
            {
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
        
        if (verbose)
			cout << "wirCreateDefaultSky: Creating default sky image " << defaultskydefaultskyname << endl;

        operaWIRCamImage defaultsky(defaultskydefaultskyname, WIRCAM_NAXIS1, WIRCAM_NAXIS2, 1, WIRCAM_EXTENSIONS, tfloat, cNone, false);
        
		// Normalize to 1 (set all pixels to 1)
		defaultsky = defaultsky + 1.0;
		
        // Save new image
        if (verbose)
			cout << "wirCreateDefaultSky: Saving outImage.fits..." << endl;
 
		defaultsky.operaWIRCamImageSave();
        
        // Close files
        if (verbose)
			cout << "wirCreateDefaultSky: Closing files..." << endl;

        defaultsky.operaFITSImageClose();
        
    }
    catch (operaException e) {
        cerr << "wirCreateDefaultSky: " << e.getFormattedMessage() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "wirCreateDefaultSky: " << operaStrError(errno) << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
