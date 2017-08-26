/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: espqlh
 Version: 1.0
 Description: A quicklook tool.
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

#include <sys/stat.h>
#include <getopt.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <algorithm>			// for max

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaIOFormats.h"

/*! \file espqlh.cpp */

using namespace std;

/*! 
 * espqlh
 * \author Doug Teeple
 * \brief A quicklook tool for use at CFHT while observing.
 * \details Calculates SNR on the fly while observing. Also
 * \details automatically performs the necessary calibrations.
 * \details Runs as a daemon in director, and is started by
 * \details the espqldaemon script which is restarted every night.
 * \arg argc
 * \arg argv
 * \ingroup tools
 * \return EXIT_STATUS
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: espqlh --night=<basename of images directory> [--director] [--cd|--directory=<directory>] [--upenadir=<directory>] --zip=.gz [--ordernumber=n] --saturation=f(65535.0) --maxsaturated=n(500) --updatelogbook -[dvth]\n";
}	

static int oldtime = 0;

bool fileChanged(string filename) {
	struct stat attrib;
	
	stat(filename.c_str(), &attrib);
	if (oldtime != attrib.st_mtime) {
		oldtime = attrib.st_mtime;
		return true;
	}
	return false;
}

int main(int argc, char *argv[])
{
	int opt;
	string upenadir = "/data/niele/espadons/opera/";	//  "/data/uhane5/opera/"; 
	string directory = "/data/niele/espadons/";
	string operainstalldir = "~/opera-1.0/";
	string optargs = "";
	string night;
	string zipped = ".gz";
	bool director = false;
	bool updatelogbook = false;
	float saturation = 65535.0;
	unsigned maxSaturatedCount = 25;
	unsigned ordernumber = 0;
	string orderstring = "";
	const string status = "status: ";
	const string info4 = "info4: ";
	const string logonly = "logonly: ";
	int upperlowerbounds = 1000;
	bool waitabit = true;
	unsigned waittime = 1;
	bool isPolar = false;
	float polarAccumulatedSNR = 0.0;
	//const float apertureheight = 0.8;
	float ccdbin = 1.0;	
	// SNR (e-/ccd pixel) = SNR (e-/spec elem) x [ 1.0 / sqrt(Np . As) ],  where
	// Np = extractionAperture->getSubpixels()->getNPixels()
	// As = extractionAperture->getSubpixels()->getSubpixelArea()
	
	
	int debug=false, verbose=false, trace=false, plot=false;
	
	struct option longopts[] = {
		{"cd",				1, NULL, 'c'},	// where to look for inputs
		{"directory",		1, NULL, 'y'},	// where to look for inputs
		{"night",			1, NULL, 'n'},	// night directory
		{"upenadir",		1, NULL, 'u'},	// upena output directory
		{"director",		0, NULL, 'r'},	// are we running in director?		
		{"zip",				1, NULL, 'z'},	// intermediate products are zipped, i.e. ".gz"	
		{"ordernumber",		1, NULL, 'o'},	// show only specific order		
		{"saturation",		1, NULL, 's'},	// set saturation limit	value in ADU	
		{"maxsaturated",	1, NULL, 'm'},	// set saturation count		
		{"wait",			1, NULL, 'w'},	// how long to sleep between images		
		{"updatelogbook",	0, NULL, 'l'},	// update the CFHT logbook with peak SNR		
		{"upperlowerbounds",1, NULL, 'b'},	// upper lower bounds for smoothing		
		{"ccdbin",			0, NULL, 'i'},	// ccd bin stats		
		
		{"plot",			0, NULL, 'p'},
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	while ((opt = getopt_long(argc, argv, "c:u:ro:s:y:n:m:w:lbi:z:pvdth", longopts, NULL))  != -1) {
		switch (opt) {
			case 'c':
			case 'y':
				directory = optarg;
				break;
			case 'n':
				//directory = directory + "/" + optarg;
				night = optarg;
				break;
			case 'u':		// upena directory used when executing in director (a local summit path)
				upenadir = optarg;
				break;	
			case 'r':
				director = true;
				break;						
			case 'o':
				ordernumber = atoi(optarg);
				orderstring = string("--ordernumber=")+optarg;
				break;						
			case 's':
				saturation = atof(optarg);
				break;						
			case 'm':
				maxSaturatedCount = atoi(optarg);
				break;						
			case 'z':
				zipped = optarg;
				break;						
			case 'l':
				updatelogbook = true;
				break;						
			case 'w':
				waittime = atoi(optarg);
				break;						
			case 'b':
				upperlowerbounds = atoi(optarg);
				break;						
			case 'i':
				ccdbin = sqrt(2.6/1.8);
				break;						
				
			case 'p':
				plot = true;
				optargs += " -p ";
				break;
			case 'v':
				if (!optarg) {
					verbose = true;
				} else {
					verbose = atoi(optarg);
				}
				optargs += " -v ";
				break;
			case 'd':
				debug = true;
				optargs += " -d ";
				break;
			case 't':
				trace = true; 
				optargs += " -t ";
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
	try {
		if (night.empty()) {
			throw operaException("espqlh: Please specify --night=<directory> ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		string filename = directory+"/current.fits";
		char *prefix = NULL;
		
		prefix = getenv("opera");
		if (prefix != NULL) {
			operainstalldir = string(prefix);
		}
		directory += night + "/";
		if (!fileexists(directory)) {
			throw operaException("espqlh: Directory does not exist: "+directory, operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (verbose) {
			cout << (director?logonly:"") << "Starting quicklook (espqlh) in directory: '" << directory << "' output to: '" << upenadir << "' OPERA home: '" << operainstalldir << "'" << endl;
		}
		/*
		 * loop forever
		 */
		while (true) {
			try {
				if (fileexists(filename.c_str())) {
					if (fileChanged(filename)) {
						string actualfilename;
						if (getRealFileName(filename, actualfilename)) {
							if (fileexists(actualfilename.c_str())) {
								string basefilename = actualfilename;
								if (actualfilename.find_last_of("/") != string::npos) {
									basefilename = actualfilename.substr(actualfilename.find_last_of("/")+1);
								}
								operaFITSImage image(actualfilename, tfloat, READONLY, cNone);
								string etype = image.operaFITSGetHeaderValue("EXPTYPE");
								string newfilename = directory + "/" + basefilename; // might not be in this dir...
								if (fileexists(newfilename.c_str())) {
									string odometer;
									waitabit = false;
									if (etype == "FLAT" || etype == "ALIGN" || etype == "OBJECT") {
										if (etype == "OBJECT") {
											odometer = basefilename.substr(0, basefilename.find("o.fits"));
										} else if (etype == "FLAT") {
											odometer = basefilename.substr(0, basefilename.find("f.fits"));
										} else if (etype == "ALIGN") {
											odometer = basefilename.substr(0, basefilename.find("a.fits"));
										}
										if (verbose) {
											cout << (director?logonly:"") << "quicklook processing image " << basefilename << endl;
										}
										unsigned saturatedCount = 0;
										float peakPixelValue = 0.0;
										unsigned maxx = image.getnaxis1();
										unsigned maxy = image.getnaxis2();
										for (unsigned j=0; j<maxy; j++) {
											for (unsigned i=0; i<maxx; i++) {
												float fluxValue = image[j][i];
												if (fluxValue >= saturation) {
													saturatedCount++;
												}
												if (fluxValue >= peakPixelValue) {
													peakPixelValue = fluxValue;
												}
											}
										}
										if (etype == "FLAT" || etype == "ALIGN") {
											cout << (director?info4:"") << basefilename << ": " << (saturatedCount==0?"zero":itos(saturatedCount)) << " saturated pixels, peak saturation: " << (int)(peakPixelValue/saturation*100.0) << "%" << endl;
										}
										if (saturatedCount > maxSaturatedCount) {
											cout << (director?info4:"") << "warning: " << basefilename << ": " << saturatedCount << " pixels above saturation limit of " << saturation << " ADU." << endl;
										}
										unsigned polarSequence = false;
										if (etype == "OBJECT") {
											string mode = image.operaFITSGetHeaderValue("INSTMODE");		// Polarimetry, star+sky, star-only,
											isPolar = mode.find("Polarimetry") != string::npos;
											if (isPolar) {
												string polarsequencestring = image.operaFITSGetHeaderValue("CMMTSEQ"); // QUIV exposure i, sequence m of n (we are looking for i, should be 1,2,3,4 but could be 1,2)
												polarSequence = isPolar?atoi(polarsequencestring.substr(11,1).c_str()):0;
											}
											if (verbose) {
												string speed = image.operaFITSGetHeaderValue("EREADSPD");		// Fast: 4.70e noise, 1.80e/ADU, 32s
												string detector = image.operaFITSGetHeaderValue("DETECTOR");	// OLAPA
												string pi_name = image.operaFITSGetHeaderValue("PI_NAME");
												string runid = image.operaFITSGetHeaderValue("RUNID");
												string object = image.operaFITSGetHeaderValue("OBJECT");
												cout << basefilename << ": " << mode << ' ' << speed << ' ' << runid << ' ' << pi_name << ' ' << object << (isPolar?" polar sequence " + itos(polarSequence):"") << endl;
												cout << (director?logonly:"") << "Starting analysis of " << basefilename << endl;
												printf("opera DATADIR=%s upenadir=%s FILE=%s spectrum %s %s\n", directory.c_str(), upenadir.c_str(), basefilename.c_str(), orderstring.c_str(), optargs.c_str());
											} 
											systemf("export opera=%s ; %s/bin/opera DATADIR=%s upenadir=%s FILE=%s spectrum %s %s >>/tmp/quicklook.log", operainstalldir.c_str(), operainstalldir.c_str(), directory.c_str(), upenadir.c_str(), basefilename.c_str(), orderstring.c_str(), optargs.c_str());
											if (verbose) {
												cout << (director?logonly:"") << "Reduction of " << basefilename << " complete." << endl;
											}
											string spectrumpath = upenadir + "/spectra/" + night + "/" + odometer + ".q" + zipped;
											string snrpath = upenadir + "/spectra/" + night + "/" + odometer + "q.sn" + zipped;
											float peak = 0.0;
											operaSpectralOrderVector spectralOrderVector;
											operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, spectrumpath);
											operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, snrpath);
											unsigned minorder = spectralOrderVector.getMinorder();
											unsigned maxorder = spectralOrderVector.getMaxorder();
											if (ordernumber != 0) {
												minorder = ordernumber;
												maxorder = ordernumber;
											}
											unsigned orderofmax = 0;
											for (unsigned order=minorder; order<=maxorder; order++) {
												operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
												operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
												if (spectralElements->getHasFluxSNR()) {
													float smoothed = spectralOrder->getCentralSmoothedSNR(upperlowerbounds);
													if (peak < smoothed) {
														peak = smoothed;
														orderofmax = order;
													}
													if (verbose) {
														cout << (director?logonly:"") << "Order " << orderofmax << " : Peak SNR for  " << basefilename << " order " << order << " : " << smoothed << " / " << ftos(smoothed*ccdbin) << endl;
													}
												}
											}
											peak *= ccdbin;	// per spectral bin / per CCD bin
											if (isPolar) {
												polarAccumulatedSNR += peak*ccdbin;
											} else {
												polarAccumulatedSNR = 0.0;
											}
											if (updatelogbook) {
												try {
													systemf("%s/bin/wiropdb \"update op..xexp set snr=%d where _obsid=%d\"", operainstalldir.c_str(), (int)(peak+0.5), atoi(odometer.c_str()));
													if (verbose)
														cout << (director?info4:"") << basefilename << " logbook updated.";
												}
												catch (...) {
													if (verbose)
														cout << (director?info4:"") << basefilename << " logbook not updated.";
												}
											}
											cout << (director?info4:"") << basefilename << ": peak saturation: " << (int)(peakPixelValue/saturation*100.0) << "%, " << "SNR: " << ((int)peak) << " / " << itos((int)(peak*ccdbin+0.5)) << (isPolar?(" ACC "+itos(polarSequence)+"/4: "+itos((int)polarAccumulatedSNR)):"") << endl;
										}
										if (isPolar && polarSequence == 4) {
											polarAccumulatedSNR = 0.0;
										}
									} else {
										cout << (director?info4:"") << basefilename << endl;
									}

								}
							}
						}
					}
				} else {
					oldtime = 0;
				}

			}
			catch (operaException e) {
				//cerr << "espqlh: " << e.getFormattedMessage() << endl;
				// keep going....
			}
			catch (...) {
				//cerr << "espqlh: " << operaStrError(errno) << endl;
				// keep going....
			}
			if (waitabit) {
				sleep(waittime);
			} else {
				waitabit = true;
			}
		}
	}
	catch (operaException e) {
		cerr << "espqlh: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "espqlh: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


