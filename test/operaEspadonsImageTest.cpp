/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaEspadonsImageTest
 Version: 1.0
 Description: Perform various tests on the operaEspadonsImage class.
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
#include <iomanip>

#include "globaldefines.h"
#include "operaError.h"
#include "operaEspadonsImageTest.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"

#include "libraries/operaImage.h"

/*! \file operaEspadonsImageTest.cpp */

using namespace std;

/*! 
 * operaEspadonsImageTest
 * \author Doug Teeple
 * \brief Test the espadonsimage class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
	int opt;
	string image;
	string output;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"image",			1, NULL, 'i'},	// input image
		{"output",			1, NULL, 'o'},	//  output fits file
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// image in
					image = optarg;
					break;
				case 'o':		// output
					output = optarg;
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
		
		if (image.empty()) {
			throw (operaErrorNoInput);
		}
		if (output.empty()) {
			throw (operaErrorNoOutput);
		}
		
		DATASEC_t datasec = ESPADONS_DEFAULT_DATASEC;
		// open the input file
		cerr << "Opening " << image << '\n';
		operaEspadonsImage In(image, tfloat, READONLY);

		cerr << "DETECTOR= " << In.operaFITSGetHeaderValue("DETECTOR") << '\n';
		cerr << "AMPLIST= " << In.operaFITSGetHeaderValue("AMPLIST") << '\n';
		cerr << "INSTMODE= " << In.operaFITSGetHeaderValue("INSTMODE") << '\n';
		cerr << "OBSTYPE= " << In.operaFITSGetHeaderValue("OBSTYPE") << '\n';
		cerr << "EREADSPD= " << In.operaFITSGetHeaderValue("EREADSPD") << '\n';
		cerr << "CMMTSEQ= " << In.operaFITSGetHeaderValue("CMMTSEQ") << '\n';

		cerr << "IMTYPE= " << In.getimtype() << ' ' << In.getimtypestring() << '\n';
		cerr << "DETECTOR= " << In.getdetector() << ' ' << In.getdetectorstring() << '\n';
		cerr << "AMPLIFIER= " << In.getamplifier() << ' ' << In.getamplifierstring() << '\n';
		cerr << "MODE= " << In.getmode() << ' ' << In.getmodestring() << '\n';
		cerr << "SPEED= " << In.getspeed() << ' ' << In.getspeedstring() << '\n';
		cerr << "STOKES= " << In.getstokes() << ' ' << In.getstokesstring() << '\n';
		cerr << "POLAR QUAD= " << In.getpolarquad() << ' ' << In.getpolarquadstring() << '\n';
		
 		// clone to the output file
		cerr << "Creating " << output << ' ' << In.getnaxis1() << 'x' << In.getnaxis2() << '\n';
		operaEspadonsImage Out(output, In.getnaxis1(), In.getnaxis2(), datasec, tfloat, 0);
		cerr << "Copying header from " << image << " to "<< output << '\n';
		Out.operaEspadonsImageCopyHeader(&In);
		operaFITSSubImage *subImage = In.getDatasecSubImage();

		cerr << "Checking operators\n";
		*subImage += 100.0;		// add bias
		*subImage *= *subImage < 35000.0;	// remove saturated pixels
		
		// print some pixel values from the center of the image
		cerr << "print some pixel values from the center of the image\n";
		for (unsigned y=2048; y<2056; y++) {
			for (unsigned x=1024; x<1032; x++) {
				cout << setw(8) << setprecision(4) << In[y][x] << ' ';	// output the pixel value
				Out[y][x] = In[y][x];
			}
			cout << endl;
		}
		cout << endl;
		// now get the datasec subImage
		cerr << "now get the datasec subImage\n";
		cerr << "print some pixel values from the center of the datasecsubimage\n";
		for (unsigned y=2048; y<2056; y++) {
			for (unsigned x=1024; x<1032; x++) {
				cout << setw(8) << setprecision(4) << (float)(*subImage)[y][x] << ' ';			// output the pixel value
			}
			cout << endl;
		}
		cout << endl;
		// now, if we ever change the subImage, which is a copy of the image
		// we need to set it back in to the espadons image
		Out.operaFITSImageSetData(*subImage);
		// Note that we use the operaFITSImage classes here
		cerr << "Saving " << output << '\n';
		Out.operaFITSImageSave();
		cerr << "Closing " << output << '\n';
		Out.operaFITSImageClose();
		cerr << "Done " << '\n';
	}
	catch (operaErrorCode errorcode) {
		operaPError("operaEspadonsImageTest", errorcode);
	}
	catch (string extra) {
		cerr << "operaEspadonsImageTest: " << extra << '\n';
	}
	catch (operaException e) {
		cerr << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaEspadonsImageTest: default error...\n";
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaTestFITSImage [--image=<filename>]+ --output=<file name> -[dvth]\n";
}	
