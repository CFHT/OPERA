/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operabiasinjector
 Version: 1.0
 Description: add a bias to amp B.
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
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"

/* \file operabiasinjector.cpp */

using namespace std;

/*!
 * operabiasinjector
 * \author Doug Teeple
 * \brief add a bias to amp B.
 * \arg argc
 * \arg argv
 * \note arg0 = keyword 
 * \note arg1, arg2, ... FITS files
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout << " Usage: operabiasinjector --input=<FITSImage>  --bias=<float> -[vdth]\n";
	cout << " Usage: adds a bias to amp B.\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string inputfilename;
	string outputfilename;
	int bias = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"bias",			1, NULL, 'b'},
		{"input",			1, NULL, 'i'},
		{"output",			1, NULL, 'o'},
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "b:i:o:dth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'b':
					bias = atoi(optarg);
					break;
				case 'i':
					inputfilename = optarg;
					break;
				case 'o':
					outputfilename = optarg;
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
		
		// we need an input...
		if (inputfilename.empty()) {
			throw operaException("operabiasinjector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputfilename.empty()) {
			throw operaException("operabiasinjector: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		remove(outputfilename.c_str());
		operaFITSImage *input = new operaFITSImage(inputfilename, READONLY);
		unsigned naxis1 = input->getnaxis1();
		unsigned naxis2 = input->getnaxis2();
		operaFITSImage *output = new operaFITSImage(outputfilename, naxis1, naxis2, tushort, cNone);
		output->operaFITSImageCopyHeader(input);
		/*
		 * copy amp A
		 */
		for (unsigned y=0; y<naxis2; y++) {
			for (unsigned x=0; x<naxis1/2; x++) {
				output->setpixel((unsigned short)input->getpixelUSHORT(x,y), x, y);
			}
		}
		/*
		 * bias amp B
		 */
		for (unsigned y=0; y<naxis2; y++) {
			for (unsigned x=naxis1/2; x<naxis1; x++) {
				output->setpixel((unsigned short)(input->getpixelUSHORT(x,y) + bias), x, y);
			}
		}
		
		output->operaFITSImageSave();
		output->operaFITSImageClose();
		input->operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "operabiasinjector: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operabiasinjector: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

