/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaExtractProducts.cpp
 Version: 1.0
 Description: Extract calibration and spectal-polarimetric data from FITS packages.
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

/*! \file operaExtractProducts.cpp */

using namespace std;

/*! 
 * operaExtractProducts
 * \author Doug Teeple
 * \brief Extract calibration and spectal-polarimetric data from FITS packages.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 * \verbatim
 *
 * Extract an opera/upena reduced fits Table 
 *  C++ Version
 *
 *
 *	   headers for column names
 *	 
 *    (1) Spectroscopy Star only mode
 *        First column = wavelength in nanometres
 *        Second column = intensity
 *        Third column = error bar
 *    
 *    (2) Polarimetry
 *        1st col = wavelength in nanometres
 *        2d  col = intensity
 *        3rd col = polarisation (Q or U or V or W)
 *        4th col = Check Spectra #1
 *        5th col = Check Spectra #2
 *        6th col = error bar
 *    
 *    (3) Spectroscopy Star + Sky
 *        1st col = wavelength
 *        2d  col = star spectra (sky subtracted)
 *        3rd col = star + sky spectra
 *        4th col = sky spectra
 *        5, 6, 7 = error bars for each column 2, 3, 4
 *
 * \endverbatim
 */

#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>

#include <fstream>
#include <iomanip>
#include <string>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/gzstream.h"

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] <product filenames> [--zip=.gz] [--overwrite] [--calibrationdir=...] [--mode=] [--speed=] [--detector=] [--amplifier=] [--oset=] [--libreesprit] \n\n"    
	"  Extracts calibration and spectral-polarimetric data from FITS packages i.fits p.fits m.fits \n"
	"  -h, --help,	   Display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,    Turn on debug messages\n"
	"  -t, --trace,    Turn on trace messages\n"
	"  -p, --plot,     Plot output \n"
	"\n";
}

int main(int argc, char *argv[])
{
	int opt;
	int debug=0, verbose=0, trace=0, plot=0;
	string productname;
	string zipped = ".gz";
	string calibrationdir;
	string directory;
	string mode;
	string speed;
	string detector;
	string amplifier;
	string oset = "";
	bool LibreEspritFormat = false;			// for backwards compatibility
	bool overwrite = false;
	
	struct option longopts[] = {
		{"zip",				1, NULL, 'z'},	// intermediate products are zipped, i.e. ".gz"	
		{"libreesprit",		0, NULL, 'l'},	// Libre-Esprit format	
		{"calibrationdir",	1, NULL, 'c'},	// for calibration files
		{"spectradir",		1, NULL, 't'},	// for spectra files
		{"mode",			1, NULL, 'm'},	
		{"speed",			1, NULL, 's'},	
		{"detector",		1, NULL, 'e'},	
		{"amplifier",		1, NULL, 'a'},	
		{"oset",			1, NULL, 'o'},	
		{"overwrite",		0, NULL, 'w'},	// overwrite if calibrtion already exists
		
		{"plot",		optional_argument,	NULL, 'p'},       
		{"verbose",		optional_argument,	NULL, 'v'},
		{"debug",		optional_argument,	NULL, 'd'},
		{"trace",		optional_argument,	NULL, 't'},
		{"help",		no_argument,		NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "lzwc:m:s:e:o:a:r:p::v::d::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'z':
				zipped = optarg;
				break;						
			case 'l':
				LibreEspritFormat = true;
				zipped = "";
				break;						
			case 'w':
				overwrite = true;
				break;						
			case 'c':
				calibrationdir = optarg;
				break;						
			case 'r':
				directory = optarg;
				break;						
			case 'm':
				mode = optarg;
				break;						
			case 's':
				speed = optarg;
				break;						
			case 'e':
				detector = optarg;
				break;						
			case 'o':
				oset = optarg;
				break;						
			case 'a':
				amplifier = optarg;
				break;						
				
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
	
	try {
		if (verbose) {
			cout << "operaExtractProducts: zipped=" << zipped << endl;
			cout << "operaExtractProducts: LibreEspritFormat=" << LibreEspritFormat << endl;
		}
		if (optind >= argc) {
			printUsageSyntax(argv[0]);
		} else {
			while (optind < argc) {
				
				productname = argv[optind++];
				string basefilename = productname;
				if (basefilename.find_last_of("/") != string::npos) {
					basefilename = basefilename.substr(basefilename.find_last_of("/")+1);
				}
				if (verbose) {
					cout << "operaExtractProducts: basefilename=" << basefilename << endl;
				}
				if (directory.empty()) {
					directory = "./";
				}
				if (productname.find("m.fits") != string::npos) {
					if (calibrationdir.empty()) {
						calibrationdir = "./";
					}
					if (verbose) {
						cout << "operaExtractProducts: calibrationdir=" << calibrationdir << endl;
					}
					operaMEFFITSProduct in(productname, READONLY);
					basefilename = basefilename.substr(0, basefilename.find("m.fits"));
					if (mode.empty()) {
						string rawmode = in.operaFITSGetHeaderValue("INSTMODE", 0);
						if (rawmode.find("Polarimetry,") != string::npos) {
							mode = "pol";
						} else if (rawmode.find("star+sky,") != string::npos) {
							mode = "sp1";
						} else {
							mode = "sp2";
						}
					}
					if (speed.empty()) {
						speed = in.operaFITSGetHeaderValue("EREADSPD", 0);
						speed = speed.substr(0, speed.find(":"));
					}
					if (detector.empty()) {
						detector = in.operaFITSGetHeaderValue("DETECTOR", 0);
					}
					if (amplifier.empty()) {
						try {	// optional keyword, doesn't always exist
							amplifier = in.operaFITSGetHeaderValue("AMPLIST", 0);
							if (amplifier.find(",") != string::npos) {
								amplifier.erase(amplifier.find(","), 1);
							}
						} catch (...) {
							amplifier = "";
						}
					}
					in.operaFITSImageClose();
					
					string calfilename = detector + amplifier + "_" + mode + oset + "_" + speed;
					string outfilename;

					if (verbose) {
						cout << "operaExtractProducts: modes: " << calfilename << endl;
					}
					operaSpectralOrderVector spectralOrders;
					
					outfilename = directory + basefilename + "i.e" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, OptimalBeamSpectrum)) {
						spectralOrders.WriteSpectralOrders(outfilename, OptimalBeamSpectrum);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = directory + basefilename + ".p" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Polarimetry)) {
						spectralOrders.WriteSpectralOrders(outfilename, Polarimetry);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = directory + basefilename + "ie.sn" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, SNR)) {
						spectralOrders.WriteSpectralOrders(outfilename, SNR);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".geom" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Geom)) {
						spectralOrders.WriteSpectralOrders(outfilename, Geom);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".ordp" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Orderspacing)) {
						spectralOrders.WriteSpectralOrders(outfilename, Orderspacing);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".disp" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Disp)) {
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					spectralOrders.WriteSpectralOrders(outfilename, Disp);
					
					outfilename = directory + basefilename + "i.rvel" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, RVel)) {
						spectralOrders.WriteSpectralOrders(outfilename, RVel);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = directory + basefilename + "i.tell" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Tell)) {
						spectralOrders.WriteSpectralOrders(outfilename, Wave);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = directory + basefilename + "p.rvel" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, PRVel)) {
						spectralOrders.WriteSpectralOrders(outfilename, PRVel);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = directory + basefilename + "p.tell" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, PTell)) {
						spectralOrders.WriteSpectralOrders(outfilename, Wave);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".wcal" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Wave)) {
						spectralOrders.WriteSpectralOrders(outfilename, Wave);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".aper" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Aperture)) {
						spectralOrders.WriteSpectralOrders(outfilename, Aperture);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".prof" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, Prof)) {
						spectralOrders.WriteSpectralOrders(outfilename, Prof);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					if (spectralOrders.ReadSpectralOrders(productname, Fcal)) {
						unsigned sequence = spectralOrders.getSequence();
						if (sequence > 0) {
							outfilename = calibrationdir + "masterfluxcalibration_" + calfilename + itos(sequence) + ".fcal" + zipped;
						} else {
							outfilename = calibrationdir + "masterfluxcalibration_" + calfilename + ".fcal" + zipped;
						}
						spectralOrders.WriteSpectralOrders(outfilename, Fcal);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
					
					outfilename = calibrationdir + calfilename + ".gain" + zipped;
					if (spectralOrders.ReadSpectralOrders(productname, GainNoise)) {
							spectralOrders.WriteSpectralOrders(outfilename, GainNoise);
						if (verbose) {
							cout << "operaExtractProducts: " << outfilename << endl;
						}
					}
				
				} else if (productname.find("p.fits") != string::npos) {
					if (verbose) {
						cout << "operaExtractProducts: dir=" << directory << endl;
					}
					operaFITSProduct in(productname, READONLY);
					basefilename = basefilename.substr(0, basefilename.find("p.fits"));
					string outfilenamebase = directory + basefilename;
					instrumentmode_t instrumentmode;
					string object;
					string mode = in.operaFITSGetHeaderValue("INSTMODE");
					if (in.getnaxis1() > in.getnaxis2()) {
						in.rotate90();
					}
					unsigned rows = (unsigned)in.getYDimension();
					unsigned columns = 5;
					if (mode.find("Polarimetry") != string::npos) {
						instrumentmode = MODE_POLAR;
					} else {
						throw operaException("operaExtractProduct: "+mode+' ', operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);	
					}
					/*
					 * pu
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "pu.s";
						fout.open(outfilename.c_str());
						object = in.operaFITSGetHeaderValue("OBJECT");
						if (verbose) {
							cout << "operaExtractProducts: mode=" << mode << endl;
							cout << "operaExtractProducts: object='" << object << "'"<< endl;
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						for (unsigned row=0; row<rows; row++) {
							fout << fixed << setprecision(4) << in[row][18] << ' ';
							fout << scientific << in[row][19] << ' ';
							fout << scientific << in[row][20] << ' ';
							fout << scientific << in[row][21] << ' ';
							fout << scientific << in[row][22] << ' ';
							fout << scientific << in[row][23] << ' ';
							fout << endl;
						}
						fout.close();
					}
					/*
					 * pn
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "pn.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						for (unsigned row=0; row<rows; row++) {
							fout << fixed << setprecision(4) << in[row][12] << ' ';
							fout << scientific << in[row][13] << ' ';
							fout << scientific << in[row][14] << ' ';
							fout << scientific << in[row][15] << ' ';
							fout << scientific << in[row][16] << ' ';
							fout << scientific << in[row][17] << ' ';
							fout << endl;
						}
						fout.close();
					}
					/*
					 * puw
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "puw.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						for (unsigned row=0; row<rows; row++) {
							fout << fixed << setprecision(4) << in[row][6] << ' ';
							fout << scientific << in[row][7] << ' ';
							fout << scientific << in[row][8] << ' ';
							fout << scientific << in[row][9] << ' ';
							fout << scientific << in[row][10] << ' ';
							fout << scientific << in[row][11] << ' ';
							fout << endl;
						}
						fout.close();
					}
					/*
					 * pnw
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "pnw.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						for (unsigned row=0; row<rows; row++) {
							fout << fixed << setprecision(4) << in[row][0] << ' ';
							fout << scientific << in[row][1] << ' ';
							fout << scientific << in[row][2] << ' ';
							fout << scientific << in[row][3] << ' ';
							fout << scientific << in[row][4] << ' ';
							fout << scientific << in[row][5] << ' ';
							fout << endl;
						}
						fout.close();
					}
					in.operaFITSImageClose();
				} else if (productname.find("i.fits") != string::npos) {
					if (verbose) {
						cout << "operaExtractProducts: dir=" << directory << endl;
					}
					operaFITSProduct in(productname, READONLY);
					basefilename = basefilename.substr(0, basefilename.find("i.fits"));
					string outfilenamebase = directory + basefilename;
					if (in.getnaxis1() > in.getnaxis2()) {
						in.rotate90();
					}
					unsigned rows = (unsigned)in.getYDimension();
					instrumentmode_t instrumentmode;
					string mode = in.operaFITSGetHeaderValue("INSTMODE");
					unsigned columns = 0;
					if (mode.find("Polarimetry") != string::npos) {
						instrumentmode = MODE_POLAR;
						columns = 2;
					} else if (mode.find("Spectroscopy, star+sky") != string::npos) {
						instrumentmode = MODE_STAR_PLUS_SKY;
						columns = 6;
					} else if (mode.find("Spectroscopy, star only") != string::npos) {
						instrumentmode = MODE_STAR_ONLY;
						columns = 2;
					} else {
						throw operaException("operaExtractProduct: "+mode+' ', operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);	
					}
					string object = in.operaFITSGetHeaderValue("OBJECT");
					/*
					 * iu
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "iu.s";
						if (verbose) {
							cout << "operaExtractProducts: mode=" << mode << endl;
							cout << "operaExtractProducts: object='" << object << "'"<< endl;
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						switch (instrumentmode) {
							case MODE_POLAR:
							case MODE_STAR_ONLY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][9] << ' ';
									fout << scientific << in[row][10] << ' ';
									fout << scientific << in[row][11] << ' ';
									fout << endl;
								}
								break;
							case MODE_STAR_PLUS_SKY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][21] << ' ';
									fout << scientific << in[row][22] << ' ';
									fout << scientific << in[row][23] << ' ';
									fout << scientific << in[row][24] << ' ';
									fout << scientific << in[row][25] << ' ';
									fout << scientific << in[row][26] << ' ';
									fout << scientific << in[row][27] << ' ';
									fout << endl;
								}
								break;
							default:
								break;
						}
						fout.close();
					}
					/*
					 * in
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "in.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						switch (instrumentmode) {
							case MODE_POLAR:
							case MODE_STAR_ONLY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][6] << ' ';
									fout << scientific << in[row][7] << ' ';
									fout << scientific << in[row][8] << ' ';
									fout << endl;
								}
								break;
							case MODE_STAR_PLUS_SKY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][14] << ' ';
									fout << scientific << in[row][15] << ' ';
									fout << scientific << in[row][16] << ' ';
									fout << scientific << in[row][17] << ' ';
									fout << scientific << in[row][18] << ' ';
									fout << scientific << in[row][19] << ' ';
									fout << scientific << in[row][20] << ' ';
									fout << endl;
								}
								break;
							default:
								break;
						}
						fout.close();
					}
					/*
					 * iuw
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "iuw.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						switch (instrumentmode) {
							case MODE_POLAR:
							case MODE_STAR_ONLY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][3] << ' ';
									fout << scientific << in[row][4] << ' ';
									fout << in[row][5] << ' ';
									fout << endl;
								}
								break;
							case MODE_STAR_PLUS_SKY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][7] << ' ';
									fout << scientific << in[row][8] << ' ';
									fout << scientific << in[row][9] << ' ';
									fout << scientific << in[row][10] << ' ';
									fout << scientific << in[row][11] << ' ';
									fout << scientific << in[row][12] << ' ';
									fout << endl;
								}
								break;
							default:
								break;
						}
						fout.close();
					}
					/*
					 * inw
					 */
					{
						ofstream fout;
						string outfilename = outfilenamebase + "inw.s";
						if (verbose) {
							cout << "operaExtractProducts: outfilename=" << outfilename << endl;
						}
						fout.open(outfilename.c_str());
						fout << "***Reduced spectrum of '" << object << "'" << endl;
						fout << rows << ' ' << columns << endl;
						switch (instrumentmode) {
							case MODE_POLAR:
							case MODE_STAR_ONLY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][0] << ' ';
									fout << scientific << in[row][1] << ' ';
									fout << scientific << in[row][2] << ' ';
									fout << endl;
								}
								break;
							case MODE_STAR_PLUS_SKY:
								for (unsigned row=0; row<rows; row++) {
									fout << fixed << setprecision(4) << in[row][0] << ' ';
									fout << scientific << in[row][1] << ' ';
									fout << scientific << in[row][2] << ' ';
									fout << scientific << in[row][3] << ' ';
									fout << scientific << in[row][4] << ' ';
									fout << scientific << in[row][5] << ' ';
									fout << scientific << in[row][6] << ' ';
									fout << endl;
								}
								break;
							default:
								break;
						}
						fout.close();
					}
					in.operaFITSImageClose();
					if (verbose) {
						cout << "operaExtractProducts: complete." << endl;
					}
				}
			}
		}
	}
	catch (operaException e) {
		cerr << "operaExtractProducts: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtractProducts: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	exit(EXIT_SUCCESS);
}

