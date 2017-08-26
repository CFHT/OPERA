/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operacompress
 Version: 1.0
 Description: Fix the acidently MEF'd headers.
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"

/*! \file operacompress.cpp */

using namespace std;

/*! 
 * operacompress
 * \author Doug Teeple
 * \brief Compress a FITS Image.
 * \arg argc
 * \arg argv
 * \note arg0 = keyword 
 * \note arg1, arg2, ... FITS files
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout << " Usage: operacompress [<file name>]+ -[grspn] -[dvth]\n";
	cout << " Usage: operacompress [<file name>]+ --gzip --extension=<extension>\n";
	cout << " Usage: operacompress [<file name>]+ --rice --extension=<extension>\n";
	cout << " Usage: operacompress [<file name>]+ --hcompress --extension=<extension>\n";
	cout << " Usage: operacompress [<file name>]+ --plio --extension=<extension>\n";
	cout << " Usage: operacompress [<file name>]+ --none --extension=<extension>\n";
	cout << " Usage: operacompress [<file name>]+ --uncompress\n";
	cout << " Usage: Compress a FITS Image (RICE compression is default, extension .fz).\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	eCompression compression = cRICE;
	string extension = ".fz";
	string typeofcompressionstring = "No";
	bool uncompress = false;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		
		{"gzip",			0, NULL, 'g'},
		{"rice",			0, NULL, 'r'},
		{"hcompress",		0, NULL, 's'},
		{"plio",			0, NULL, 'p'},
		{"none",			0, NULL, 'n'},
		{"extension",		1, NULL, 'e'},
		{"compressiontype",	1, NULL, 'C'},
		{"uncompress",		0, NULL, 'u'},
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "grspne:uvdth", longopts, NULL))  != -1) {
			switch (opt) {
					
				case 'g':
					compression = cGZIP;
					extension = ".gz";
					typeofcompressionstring = "gzip";
					break;
				case 'r':
					compression = cRICE;
					extension = ".fz";
					typeofcompressionstring = "RICE";
					break;
				case 's':
					compression = cHCOMPRESS;
					extension = ".hz";
					typeofcompressionstring = "HCOMPRESS";
					break;
				case 'p':
					compression = cPLIO;
					extension = ".pz";
					typeofcompressionstring = "PLIO";
					break;
				case 'n':
					compression = cNone;
					extension = "";
					typeofcompressionstring = "No";
					break;
				case 'e':
					extension = optarg;
					break;
				case 'C':
					compression = (eCompression)atoi(optarg);
					switch (compression) {
						case cNone:
							extension = "";
							break;
						case cRICE:
							extension = ".fz";
							typeofcompressionstring = "RICE";
							break;
						case cGZIP:
							extension = ".gz";
							typeofcompressionstring = "gzip";
							break;
						case cHCOMPRESS:
							extension = ".hz";
							typeofcompressionstring = "HCOMPRESS";
							break;
						case cPLIO:
							extension = ".pz";
							typeofcompressionstring = "PLIO";
							break;
						default:
							break;
					}
					break;    
				case 'u':
					uncompress = true;
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
			string infilename = argv[optind++];
			string outfilename;
			if (uncompress) {
				outfilename = infilename;
				size_t found = outfilename.rfind(".");
				if (found != string::npos) {
					outfilename = outfilename.substr(0, found);
				}
				if (verbose) {
					cout << "operacompress: Uncompressing " << infilename << " to " << outfilename << endl;
				}
			} else {
				outfilename = infilename+extension;
				if (verbose) {
					cout << "operacompress: Compressing " << infilename << " to " << outfilename << " using " << typeofcompressionstring << " compression." << endl;
				}
			}
#ifdef I_WISH_THIS_WOULD_WORK
			operaFITSImage *fitsfile = new operaFITSImage(infilename, tushort, READONLY, compression);
			remove(outfilename.c_str());
			operaFITSImage *outfile = new operaFITSImage(outfilename, fitsfile->getnaxis1(), fitsfile->getnaxis2(), tushort, compression);
			outfile->operaFITSImageCopyHeader(fitsfile);
			outfile = fitsfile;
			outfile->operaFITSImageSave();
			outfile->operaFITSImageClose();
			fitsfile->operaFITSImageClose();
			//delete outfile; //strange, causes a segfault
			delete fitsfile;
#else
			fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
			int status = 0, ii = 1, iteration = 0, single = 0, hdupos;
			int hdutype, bitpix, bytepix, naxis = 0, nkeys, datatype = 0, anynul;
			long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
			long first, totpix = 0, npix;
			double *array, bscale = 1.0, bzero = 0.0, nulval = 0.;
			char card[81];
			
			remove(outfilename.c_str());
			/* Open the input file and create output file */
			fits_open_file(&infptr, infilename.c_str(), READONLY, &status);
			fits_create_file(&outfptr, outfilename.c_str(), &status);
			if (status != 0) {    
				fits_report_error(stderr, status);
				return(status);
			}
			fits_set_compression_type(outfptr, compression, &status);
			if (status != 0) {    
				fits_report_error(stderr, status);
				return(status);
			}
			
			fits_get_hdu_num(infptr, &hdupos);  /* Get the current HDU position */
			
			/* Copy only a single HDU if a specific extension was given */ 
			if (hdupos != 1 || strchr(argv[1], '[')) single = 1;
			
			for (; !status; hdupos++)  /* Main loop through each extension */
			{
				
				fits_get_hdu_type(infptr, &hdutype, &status);
				
				if (hdutype == IMAGE_HDU) {
					
					/* get image dimensions and total number of pixels in image */
					for (ii = 0; ii < 9; ii++)
						naxes[ii] = 1;
					
					fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
					
					totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
				}
				
				if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0) { 
					
					/* just copy tables and null images */
					fits_copy_hdu(infptr, outfptr, 0, &status);
					
				} else {
					
					/* Explicitly create new image, to support compression */
					fits_create_img(outfptr, bitpix, naxis, naxes, &status);
					
					/* copy all the user keywords (not the structural keywords) */
					fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
					
					for (ii = 1; ii <= nkeys; ii++) {
						fits_read_record(infptr, ii, card, &status);
						if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
							fits_write_record(outfptr, card, &status);
					}
					
					switch(bitpix) {
						case BYTE_IMG:
							datatype = TBYTE;
							break;
						case SHORT_IMG:
							datatype = TSHORT;
							break;
						case LONG_IMG:
							datatype = TLONG;
							break;
						case FLOAT_IMG:
							datatype = TFLOAT;
							break;
						case DOUBLE_IMG:
							datatype = TDOUBLE;
							break;
					}
					
					bytepix = abs(bitpix) / 8;
					
					npix = totpix;
					iteration = 0;
					
					/* try to allocate memory for the entire image */
					/* use double type to force memory alignment */
					array = (double *) calloc(npix, bytepix);
					
					/* if allocation failed, divide size by 2 and try again */
					while (!array && iteration < 10)  {
						iteration++;
						npix = npix / 2;
						array = (double *) calloc(npix, bytepix);
					}
					
					if (!array)  {
						printf("operacompress: Memory allocation error\n");
						return(EXIT_FAILURE);
					}
					
					/* turn off any scaling so that we copy the raw pixel values */
					fits_set_bscale(infptr,  bscale, bzero, &status);
					fits_set_bscale(outfptr, bscale, bzero, &status);
					
					first = 1;
					while (totpix > 0 && !status) {
						/* read all or part of image then write it back to the output file */
						fits_read_img(infptr, datatype, first, npix, &nulval, array, &anynul, &status);
						fits_write_img(outfptr, datatype, first, npix, array, &status);
						totpix = totpix - npix;
						first  = first  + npix;
					}
					free(array);
				}
				
				if (single) 
					break;  /* quit if only copying a single HDU */
				fits_movrel_hdu(infptr, 1, NULL, &status);  /* try to move to next HDU */
			}
			
			if (status == END_OF_FILE)
				status = 0; /* Reset after normal error */
			
			fits_close_file(outfptr,  &status);
			fits_close_file(infptr, &status);
			
			/* if error occurred, print out error message */
			if (status)
				fits_report_error(stderr, status);
#endif
		}
	}
	catch (operaException e) {
		cerr << "operacompress: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operacompress: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

