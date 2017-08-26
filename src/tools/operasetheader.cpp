/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operasetheader
 Version: 1.0
 Description: Sets the values of a header keyowrd (or deletes) in a list of files.
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
#include <iostream>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaImage.h"

/*! \file operasetheader.cpp */

using namespace std;

/*! 
 * operasetheader
 * \author Doug Teeple
 * \brief Sets the values of a header keyowrd (or deletes) in a list of files.
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
	
	cout << " Usage: operasetheader [--keyword=<keyword>]+ [--value=<value>]+ [--comment=<comment>]+ [--delete] [--extension=<n>] [<file name>[n]]+ -[dvth]\n";
	cout << " Usage: sets keyword to value with optional comment.\n";
	cout << " Usage: if comment is given without keyword the it is added as a COMMENT.\n";
	cout << " Usage: --delete deletes the keyword (not used in conjunction wuth set).\n";
	cout << " Usage: --extension sets the keyword in an extension (default is all).\n";
	cout << " Usage: The extension is optional ([n]) and can be the extension name.\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string keywords[MAXKEYWORDS];
	string values[MAXKEYWORDS];
	string comments[MAXKEYWORDS];
	bool   deletes[MAXKEYWORDS];
	unsigned keyIndex = 0;
	unsigned valueIndex = 0;
	unsigned commentIndex = 0;
	unsigned deleteIndex = 0;
	unsigned extension = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"keyword",			1, NULL, 'k'},	// input keyword
		{"value",			1, NULL, 'a'},	// value
		{"comment",			1, NULL, 'c'},	// add a comment
		{"delete",			0, NULL, 'x'},	// delete the keyword
		{"extension",		1, NULL, 'e'},	// a particular extension
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		for (unsigned i=0; i<MAXKEYWORDS; i++) {
			deletes[i] = false;
		}
		while ((opt = getopt_long(argc, argv, "k:a:c:e:xvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'k':		// keyword
					keywords[keyIndex++] = optarg;
					break;
				case 'a':		// value
					values[valueIndex++] = optarg;
					break;
				case 'c':		//add a comment
					comments[commentIndex++] = optarg;
					break;
				case 'x':		// delete the keyword
					deletes[deleteIndex++] = true;
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
			unsigned XDimension, YDimension, ZDimension, Extensions;
			edatatype Datatype;
			long Npixels;
			bool isMef = false;
			
			ifstream ifile(argv[optind]);
			if (ifile.good()) {
				getFITSImageInformation(argv[optind], &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
				isMef = Extensions > 1;
				operaFITSImage *fitsfile = new operaFITSImage(argv[optind], tfloat, cNone, true);
				for (unsigned i=0; i<keyIndex; i++) {
					if (deletes[i]) {
						try {
							if (isMef) {
								((operaMultiExtensionFITSImage*)fitsfile)->operaFITSDeleteHeaderKey(keywords[i], extension);
							} else {
								fitsfile->operaFITSDeleteHeaderKey(keywords[i]);
							}
							if (verbose) {
								cout << "Deleted " << keywords[i] << " from " << argv[optind] << endl;
							}
						}
						catch (operaException e) {
							// it is OK if the header is missing...
						}
					} else if (values[i].empty()) {
						if (comments[i].empty()) {
							cerr << "operasetheader: " << keywords[i] << " needs an associated comment or value or deletion flag -- skipped..." << endl;
						} else {
							if (isMef) {
								((operaMultiExtensionFITSImage*)fitsfile)->operaFITSAddComment(comments[i], extension);
							} else {
								fitsfile->operaFITSAddComment(comments[i]);
							}
							if (verbose) {
								cout << "Added comment '" << comments[i] << "' to " << argv[optind] << endl;
							}
						}
					} else {
						try {
							if (isdigit(values[i].c_str()[0]) || values[i].c_str()[0] == '-' || values[i].c_str()[0]  == '+' ) {
								if (values[i].find(".") != string::npos) {
									if (isMef) {
										((operaMultiExtensionFITSImage*)fitsfile)->operaFITSSetHeaderValue(keywords[i], (float)atof(values[i].c_str()), (comments[i].empty()?"":comments[i]), extension);
									} else {
										fitsfile->operaFITSSetHeaderValue(keywords[i], (float)atof(values[i].c_str()), (comments[i].empty()?"":comments[i]));
									}
								} else {
									if (values[i].c_str()[0] == '-') {
										if (isMef) {
											((operaMultiExtensionFITSImage*)fitsfile)->operaFITSSetHeaderValue(keywords[i], (short)atoi(values[i].c_str()), (comments[i].empty()?"":comments[i]), extension);
										} else {
											fitsfile->operaFITSSetHeaderValue(keywords[i], (short)atoi(values[i].c_str()), (comments[i].empty()?"":comments[i]));
										}
									} else {
										if (isMef) {
											((operaMultiExtensionFITSImage*)fitsfile)->operaFITSSetHeaderValue(keywords[i], (unsigned short)atoi(values[i].c_str()), (comments[i].empty()?"":comments[i]), extension);
										} else {
											fitsfile->operaFITSSetHeaderValue(keywords[i], (unsigned short)atoi(values[i].c_str()), (comments[i].empty()?"":comments[i]));
										}
									}
								}
							} else {
								if (isMef) {
									((operaMultiExtensionFITSImage*)fitsfile)->operaFITSSetHeaderValue(keywords[i], values[i].c_str(), (comments[i].empty()?"":comments[i]), extension);
								} else {
									fitsfile->operaFITSSetHeaderValue(keywords[i], values[i].c_str(), (comments[i].empty()?"":comments[i]));
								}
							}
							
							if (verbose) {
								cout << "Set " << keywords[i] << " to "  << values[i] << (comments[i].empty()?"":comments[i]) << " in " << argv[optind] << (isMef?(" extension "+itos(extension)):"") << endl;
							}
						}
						catch (operaException e) {
							// it is OK if the header is missing...
						}
					}
					fitsfile->operaFITSImageClose();
					delete fitsfile;
				}
			} else {
				cerr << "operasetheader: "+string(argv[optind])+" does not exist." << endl;	
			}
			optind++;
		}
	}
	catch (operaException e) {
		cerr << "operasetheader: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operasetheader: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

