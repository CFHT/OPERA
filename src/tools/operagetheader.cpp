/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operagetheader
 Version: 1.0
 Description: Get the values of a hedaer keyowrd from a list of files.
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
#include "libraries/operaMultiExtensionFITSImage.h"

#include "libraries/operaImage.h"

/*! \file operagetheader.cpp */

using namespace std;

/*! 
 * operagetheader
 * \author Doug Teeple
 * \brief  Get the values of a hedaer keyowrd from a list of files.
 * \arg argc
 * \arg argv
 * \note arg0 = keyword 
 * \note arg1, arg2, ... FITS files
 * \return EXIT_STATUS
 * \ingroup tools
 */

#define MAXKEYWORDS 1000

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout << " Usage: operagetheader [--printfilename|-p ] [--keyword|-k=<keyword>]+ [--info|-i] [--pixels|-p] [<file name>[n]]+ -[dvth]\n";
	cout << " Usage: prints list if keywords from list of filenames, or whole header if no filename given.\n";
	cout << " Usage: Use --printfilename to identify which header the keyword is from.\n";
	cout << " Usage: The extension is optional ([n]) and can be the extension name.\n";
	cout << " Usage: --pixels dumps the pixels\n";
	cout << " Usage: --info dumps the basic header info\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string keywords[MAXKEYWORDS];
	unsigned keyIndex = 0;
	bool printfilename  = false;
	bool info  = false;
	bool pixels  = false;
	int x = -1, y = -1;
	unsigned extension = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"keyword",			1, NULL, 'k'},	// input keyword
		{"printfilename",	0, NULL, 'p'},	// also print the filename
		{"info",			0, NULL, 'i'},	// dump basic info
		{"pixels",			0, NULL, 's'},	// print the pixels
		{"x",				1, NULL, 'x'},	// print the pixel at x
		{"y",				1, NULL, 'y'},	// print the pixel at y
		{"extension",		1, NULL, 'e'},	// a particular extension
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "k:x:y:e:sipvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'k':		// keyword
					keywords[keyIndex++] = optarg;
					break;
				case 'p':		// printfilename
					printfilename  = true;
					break;
				case 'i':		// brief info only
					info  = true;
					break;
				case 's':		// pixels
					pixels  = true;
					break;
				case 'x':		// pixel at x
					pixels  = true;
					x  = atoi(optarg);
					break;
				case 'y':		// pixel at y
					pixels  = true;
					y  = atoi(optarg);
					break;
				case 'e':		// a single extension
					extension = atoi(optarg);
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
		
		while (optind < argc) {
			ifstream ifile(argv[optind]);
			if (ifile.good()) {
				if (keyIndex == 0) {		// then dump whole header
					if (pixels) {
						operaFITSImage *fitsfile = new operaFITSImage(argv[optind], tfloat, READONLY, cNone, true);
						if (x >= 0 && y >= 0) {
							cout << *fitsfile[y][x] << endl;
						} else if (x >= 0) {
							for (unsigned col=0; col<fitsfile->getnaxis1(); col++) {
								cout << *fitsfile[col][x] << ' ';
							}
							cout << endl;
						} else if (y >= 0) {
							for (unsigned row=0; row<fitsfile->getnaxis2(); row++) {
								cout << *fitsfile[y][row] << ' ';
							}
							cout << endl;
						} else {
							for (unsigned row=0; row<fitsfile->getnaxis2(); row++) {
								for (unsigned col=0; col<fitsfile->getnaxis1(); col++) {
									cout << *fitsfile[col][row] << ' ';
								}
								cout << endl;
							}
						}
						fitsfile->operaFITSImageClose();
						delete fitsfile;
					} else {
						fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
						char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
						int status = 0;			/* CFITSIO status value MUST be initialized to zero! */
						bool single = false;
						int hdupos, numHDU, hduType, nkeys;
						
						if (!fits_open_file(&fptr, argv[optind], READONLY, &status)) {
							
							fits_get_num_hdus(fptr, &numHDU, &status);  /* Get the number of HDUs */
							
							if (numHDU == 1)
								single = true;
							
							fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
							
							/* List only a single header if a specific extension was given */ 
							if (hdupos != 1 || strchr(argv[optind], '[')) 
								single = true;
							
							if (info) {
								
								char zvalue[FLEN_VALUE], comment[FLEN_COMMENT];
								
								cout << argv[optind] << ':' << endl;
								cout << numHDU << " extension" << (numHDU>1?"s":"") << endl;
								
								fits_get_hdu_type(fptr, &hduType, &status);
								switch(hduType) {
									case IMAGE_HDU : cout << "HDU Type: Image" << endl;        break;
									case ASCII_TBL : cout << "HDU Type: ASCII Table" << endl;  break;
									case BINARY_TBL: cout << "HDU Type: Binary Table" << endl; break;
									default        : cout << "HDU Type: UNKNOWN" << endl;      break;
								}
								
								fits_read_keyword(fptr, "NAXIS1", zvalue, comment, &status);
								cout << "Naxis 1 = " << string(zvalue) << endl;
								fits_read_keyword(fptr, "NAXIS2", zvalue, comment, &status);
								cout << "Naxis 2 = " << string(zvalue) << endl;
								fits_read_keyword(fptr, "BITPIX", zvalue, comment, &status);
								cout << "BITPIX = " << string(zvalue) << ' ' << string(comment) << endl;
								
							} else {
								for (; !status; hdupos++)  /* Main loop through each extension */ {
									fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */
									
									for (int i = 1; i <= nkeys; i++) { /* Read and print each keywords */
										
										if (fits_read_record(fptr, i, card, &status))
											break;
										if (!single) {
											cout << "[" << hdupos << "] ";
										}
										cout << card << endl;
									}
									
									cout << "END" << endl << endl;  /* terminate listing with END */
									
									if (single) 
										break;  /* quit if only listing a single header */
									
									fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
								}
							}
							
							if (status == END_OF_FILE)  
								status = 0; /* Reset after normal error */
							
							fits_close_file(fptr, &status);
						}
					}
				} else {
					unsigned XDimension, YDimension, ZDimension, Extensions;
					edatatype Datatype;
					long Npixels;
					bool isMef = false;
					
					getFITSImageInformation(argv[optind], &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
					isMef = Extensions > 1;
					operaFITSImage *fitsfile = new operaFITSImage(argv[optind], tushort, READONLY, cNone, true);
					for (unsigned i=0; i<keyIndex; i++) {
						if (printfilename)
							cout << argv[optind] << ' ';
						// the header may not exist, so print an error but keep going
						try {
							if (isMef) {
								cout << ((operaMultiExtensionFITSImage*)fitsfile)->operaFITSGetHeaderValue(keywords[i], extension) << ' ';
							} else {
								cout << fitsfile->operaFITSGetHeaderValue(keywords[i]) << ' ';
							}
						}
						catch (operaException e) {
							if (verbose)
								cout << endl << e.getFormattedMessage() << endl;
							// note issue, but keep going....
						}
					}
					cout << endl;
					fitsfile->operaFITSImageClose();
					delete fitsfile;
				}
			} else {
				cerr << "operagetheader: "+string(argv[optind])+" does not exist." << endl;	
			}
			
			optind++;
		}
	}
	catch (operaException e) {
		cerr << "operagetheader: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operagetheader: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

