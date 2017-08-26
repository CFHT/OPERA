/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateBadPixelMask
 Version: 1.0
 Description: Create a FITS bad pixel mask from a bad pixel map.
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
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaCreateBadpixMask.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaEspadonsImage.h"

/*! \brief command line interface to create a bad pixel mask FITS file from a map. */
/*! \file operaCreateBadpixMask.cpp */

/*!
 * \defgroup tools Tools
 */


using namespace std;

/*! 
 * operaCreateBadpixMask
 * \author Doug Teeple
 * \brief command line interface to create a bad pixel mask FITS file from a map.
 * \arg argc
 * \arg argv
 * \note --mask=...
 * \note --map=...
 * \return EXIT_STATUS
 */

int main(int argc, char *argv[])
{
	int opt;
	unsigned nx = ESPADONS_DEFAULT_NAXIS1;
	unsigned ny = ESPADONS_DEFAULT_NAXIS2;
	string badpixelmask;
	string badpixelmap;
	eCompression compression = cNone;
	
	int debug=0, verbose=0, trace=0;
	
	struct option longopts[] = {
		{"mask",1, NULL, 'm'},
		{"map",1, NULL, 'p'},
		{"naxis1",1, NULL, '1'},
		{"naxis2",1, NULL, '2'},
		{"compressiontype", 1, NULL, 'C'},			
		
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "m:p:1:2:C:vdth", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'p':		// badpixelmap
				badpixelmap = optarg;
				break;    
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;            
			case '1':		// naxis1
				nx = atoi(optarg);
				break;            
			case '2':		// naxis2
				ny = atoi(optarg);
				break;            
			case 'C':
				compression = (eCompression)atoi(optarg);	
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
				printUsageSyntax();
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax();
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	try {
		
		if (badpixelmap.empty()) {
			throw operaException("operaCreateBadPixelMask: badpixelmap not given ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (badpixelmask.empty()) {
			badpixelmask = badpixelmap + ".fits";
			cout << "Mask name not given, using " << badpixelmask << endl;
		}
		
		unsigned bad = 0;
		unsigned short *xl = NULL;
		unsigned short *xu = NULL;
		unsigned short *yl = NULL;
		unsigned short *yu = NULL;
		unsigned short count = 0;
		
		remove(badpixelmask.c_str());
		operaFITSImage *badpiximage = new operaFITSImage(badpixelmask, nx, ny, tushort, compression);
			
		ifstream flist(badpixelmap.c_str());		
		if (flist.is_open())
		{
			string line;
			while (flist.good()) {
				getline (flist, line);
				if (line.size() == 0 ||line[0] == '#')
					continue;
				if (bad == 0) {
					bad = atoi(line.c_str());
					xl = new unsigned short[bad];
					xu = new unsigned short[bad];
					yl = new unsigned short[bad];
					yu = new unsigned short[bad];
					if (verbose) {
						cout << "Found " << bad << " entries in map " << badpixelmap << endl;
					}
				} else {
					if (sscanf(line.c_str(), "%hu %hu %hu %hu", &xl[count], &xu[count], &yl[count], &yu[count])) {
						if (xl[count] > nx) xl[count] = nx;
						if (xl[count] > 0) xl[count]--;
						if (xu[count] > nx) xu[count] = nx;
						if (xu[count] > 0) xu[count]--;
						if (yl[count] > ny) yl[count] = ny;
						if (yl[count] > 0) yl[count]--;
						if (yu[count] > ny) yu[count] = ny;
						if (yu[count] > 0) yu[count]--;
						count++;
					}
				} 
			}	
			flist.close();
		} else
			throw operaException("operaCreateBadPixelMask: could not find bad pixel map "+badpixelmap+" ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		
		if (verbose) {
			cout << "count = " << count << " bad = " << bad << " nx = " << nx << " ny = " << ny << '\n';
			for (unsigned i = 0; i < count; i++)
				cout << i << ": " << xl[i] << ' ' << xu[i] << ' ' << yl[i] << ' ' << yu[i] << '\n';
		}
		if (!count) {
			throw operaException("operaCreateBadPixelMask: could not read data from bad pixel map "+badpixelmap+" ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		// set everything to one
		for (unsigned i = 0; i < nx; i++) {
			for (unsigned j = 0; j < ny; j++) {
				badpiximage->setpixel((unsigned short)1,i,j);
			}
		}
		
		// reset the bad pixels to zero
		for (unsigned i = 0; i < count; i++) {
			for (unsigned j = xl[i]; j <= xu[i]; j++) {
				for (unsigned k = yl[i]; k <= yu[i]; k++) {
					badpiximage->setpixel((unsigned short)0,j,k);
				}
			}
		}

		delete[] xl;
		delete[] xu;
		delete[] yl;
		delete[] yu;
		badpiximage->operaFITSImageSave();
		badpiximage->operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "operaCreateBadPixelMask: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateBadPixelMask: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaCreateBadPixelMask  --mask=<filenpath> --map=<filepath> [naxis1=<unsigned> --naxis2=<unsigned>] --compressiontype=<cfitsio compression type>\n";
}	
