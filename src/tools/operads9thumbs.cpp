/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operads9thumbs
 Version: 1.0
 Description: Makes a 4 up display of the last 4 images.
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaStats.h"			// median

/*! \file operads9thumbs.cpp */

using namespace std;

/*!
 * operads9thumbs
 * \author Doug Teeple
 * \brief Makes a 4 up display of the last 4 images.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operads9thumbs --input[1234]=<FITS image filename>] --output=<filename>.png-[dvth]\n";
}

int main(int argc, char *argv[])
{
	int opt;
	string inputfilename1;
	string inputfilename2;
	string inputfilename3;
	string inputfilename4;
	string outputfilename;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"input1",			1, NULL, '1'},			// input
		{"input2",			1, NULL, '2'},			// input
		{"input3",			1, NULL, '3'},			// input
		{"input4",			1, NULL, '4'},			// input
		{"output",			1, NULL, 'o'},			// output
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "1:2:3:4:o:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case '1':
					inputfilename1 = optarg;
					break;
				case '2':
					inputfilename2 = optarg;
					break;
				case '3':
					inputfilename3 = optarg;
					break;
				case '4':
					inputfilename4 = optarg;
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
		
		if (inputfilename4.empty() || outputfilename.empty()) {
			throw operaException("operads9thumbs: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		string basefilename = outputfilename;
		string directory = "./";
		if (outputfilename.find_last_of("/") != string::npos) {
			basefilename = outputfilename.substr(outputfilename.find_last_of("/")+1);
			directory = outputfilename.substr(0, outputfilename.find_last_of("/"));
		}
		string tempoutfilename = directory + '.' + basefilename;
		if (verbose) {
			cout << "operads9thumbs: inputfilename1: " << inputfilename1 << endl;
			cout << "operads9thumbs: inputfilename2: " << inputfilename2 << endl;
			cout << "operads9thumbs: inputfilename3: " << inputfilename3 << endl;
			cout << "operads9thumbs: inputfilename4: " << inputfilename4 << endl;
			cout << "operads9thumbs: outputfilename: " << outputfilename << endl;
			cout << "operads9thumbs: tempoutfilename: " << tempoutfilename << endl;
		}
		operaFITSImage out(tempoutfilename, 2080, 4640, tfloat, cNone);
		if (!inputfilename1.empty() && !inputfilename2.empty() && !inputfilename3.empty() && !inputfilename4.empty()) {
			if (inputfilename1.length() && inputfilename2.length() && inputfilename3.empty() && inputfilename4.length()) {
				operaFITSImage in1(inputfilename1, tfloat, cNone);		
				operaFITSImage in2(inputfilename2, tfloat, cNone);		
				operaFITSImage in3(inputfilename3, tfloat, cNone);		
				operaFITSImage in4(inputfilename4, tfloat, cNone);		
				out = operaArrayMedian(in4.getnpixels(), (float *)in4.getpixels());
				unsigned xmax = in4.getnaxis1();
				unsigned ymax = in4.getnaxis2();
				for (unsigned y=0; y<ymax; y++) {
					for (unsigned x=0; x<xmax; x++) {
						out[y][x/4] = in1[y][x/4];
						out[y][(xmax/4)+x/4] = in2[y][x/4];
						out[y][(xmax/4)*2+x/4] = in3[y][x/4];
						out[y][(xmax/4)*3+x/4] = in4[y][x/4];
					}
				}
				out.operaFITSImageSave();
				out.operaFITSImageClose();
				rename(tempoutfilename.c_str(), outputfilename.c_str());
				return EXIT_SUCCESS;
			}
		}
		if (!inputfilename2.empty() && !inputfilename3.empty() && !inputfilename4.empty()) {
			if (inputfilename2.length() && inputfilename3.length() && inputfilename4.length()) {
				operaFITSImage in2(inputfilename2, tfloat, cNone);		
				operaFITSImage in3(inputfilename3, tfloat, cNone);		
				operaFITSImage in4(inputfilename4, tfloat, cNone);		
				out = operaArrayMedian(in4.getnpixels(), (float *)in4.getpixels());
				unsigned xmax = in4.getnaxis1();
				unsigned ymax = in4.getnaxis2();
				for (unsigned y=0; y<ymax; y++) {
					for (unsigned x=0; x<xmax; x++) {
						out[y][(xmax/4)+x/4] = in2[y][x/4];
						out[y][(xmax/4)*2+x/4] = in3[y][x/4];
						out[y][(xmax/4)*3+x/4] = in4[y][x/4];
					}
				}
				out.operaFITSImageSave();
				out.operaFITSImageClose();
				rename(tempoutfilename.c_str(), outputfilename.c_str());
				return EXIT_SUCCESS;
			}
		}
		if (!inputfilename3.empty() && !inputfilename4.empty()) {
			if (inputfilename3.length() && inputfilename4.length()) {
				operaFITSImage in3(inputfilename3, tfloat, cNone);		
				operaFITSImage in4(inputfilename4, tfloat, cNone);		
				unsigned xmax = in4.getnaxis1();
				unsigned ymax = in4.getnaxis2();
				out = operaArrayMedian(in4.getnpixels(), (float *)in4.getpixels());
				for (unsigned y=0; y<ymax; y++) {
					for (unsigned x=0; x<xmax; x++) {
						out[y][(xmax/4)*2+x/4] = in3[y][x/4];
						out[y][(xmax/4)*3+x/4] = in4[y][x/4];
					}
				}
				out.operaFITSImageSave();
				out.operaFITSImageClose();
				rename(tempoutfilename.c_str(), outputfilename.c_str());
				return EXIT_SUCCESS;
			}
		}
		if (!inputfilename4.empty()) {
			if (inputfilename4.length()) {
				operaFITSImage in4(inputfilename4, tfloat, cNone);		
				out = operaArrayMedian(in4.getnpixels(), (float *)in4.getpixels());
				unsigned xmax = in4.getnaxis1();
				unsigned ymax = in4.getnaxis2();
				for (unsigned y=0; y<ymax; y++) {
					for (unsigned x=0; x<xmax; x++) {
						out[y][(xmax/4)*3+x/4] = in4[y][x/4];
					}
				}
				out.operaFITSImageSave();
				out.operaFITSImageClose();
				rename(tempoutfilename.c_str(), outputfilename.c_str());
				return EXIT_SUCCESS;
			}
		}
	}
	catch (operaException e) {
		cerr << "operads9thumbs: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operads9thumbs: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

