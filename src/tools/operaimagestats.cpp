/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaimagestats
 Version: 1.0
 Description: Print various image statistics.
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
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaStats.h"
#include "libraries/Polynomial.h"

/*! \file operaimagestats.cpp */

using namespace std;

/*! 
 * operaimagestats
 * \author Doug Teeple
 * \brief Print various image statistics.
 * \arg argc
 * \arg argv
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operaimagestats [--median ] [--peak] [--saturated] [--saturation=<ADU> (65535)] [--maxsnr] [--snr] [--snrpoint] [--fluxpoint] [--order=<n>] [--spectrum] [--count=<n>] [<file name>]+ -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	string imagefilename;
	string geomfilename;
	string wavefilename;
	string specfilename;
	string aperfilename;
	
	bool median = false;
	bool peak = false;
	bool saturated = false;
	bool peakunsaturated = false;
	bool smoothedpeaksnr = false;
	bool maxsnr = false;
	bool snr = false;
	bool spectrum = false;
	bool peakunsaturatedorder = false;
	bool normalize = false;
	bool snrpoint = false;
	bool fluxpoint = false;
	bool noprintorder = false;
	bool getorder = false;
	bool orders = false;
	bool tilt = false;
	bool isbin = false;						// for fluxpoint
	unsigned long saturation = 65535;
	unsigned order = 0;
	unsigned count = 0;
	float x = 0.0;
	float y = 0.0;
	float wl = 0.0;
	float aperture = 26.0;
	unsigned long naxis1 = 2080;
	unsigned long naxis2 = 4096;
	int upperlowerbounds = 10;
	
	bool debug=false, verbose=false, trace=false, plot=false;
	
	struct option longopts[] = {
		{"median",					0, NULL, 'm'},
		{"peak",					0, NULL, 'k'},
		{"saturated",				0, NULL, 's'},
		{"peakunsaturated",			0, NULL, 'a'},
		{"saturation",				1, NULL, 'l'},
		{"peakunsaturatedorder",	0, NULL, 'e'},
		{"smoothedpeaksnr",			0, NULL, 'S'},
		{"maxsnr",					0, NULL, 'r'},
		{"order",					1, NULL, 'o'},
		{"snr",						0, NULL, 'n'},
		{"spectrum",				0, NULL, 'u'},
		{"normalize",				0, NULL, 'z'},
		{"snrpoint",				0, NULL, 'i'},
		{"fluxpoint",				0, NULL, 'f'},
		{"x",						1, NULL, 'x'},
		{"y",						1, NULL, 'y'},
		{"wl",						1, NULL, 'W'},
		{"naxis1",					1, NULL, '1'},
		{"naxis2",					1, NULL, '2'},
		{"count",					1, NULL, 'c'},
		{"imagefilename",			1, NULL, 'G'},
		{"geomfilename",			1, NULL, 'g'},
		{"specfilename",			1, NULL, 'U'},
		{"wavefilename",			1, NULL, 'w'},
		{"noprintorder",			0, NULL, 'N'},
		{"getorder",				0, NULL, 'O'},
		{"orders",					0, NULL, 'R'},
		{"aperture",				1, NULL, 'A'},
		{"aperturetilt",			1, NULL, 'T'},
		{"aperfilename",			1, NULL, 'P'},
		{"upperlowerbounds",		1, NULL, 'b'},		

		{"plot",					0, NULL, 'p'},
		{"verbose",					0, NULL, 'v'},
		{"debug",					0, NULL, 'd'},
		{"trace",					0, NULL, 't'},
		{"help",					0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "mksSfil:o:c:x:y:W:1:2:G:g:U:w:A:TP:b:NORrpnueazpvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'm':
					median = true;
					break;
				case 'k':
					peak = true;
					break;
				case 's':
					saturated = true;
					break;
				case 'r':
					maxsnr = true;
					break;
				case 'n':
					snr = true;
					break;
				case 'u':
					spectrum = true;
					break;
				case 'e':
					peakunsaturatedorder = true;
					break;
				case 'a':
					peakunsaturated = true;
					break;
				case 'S':
					smoothedpeaksnr = true;
					break;
				case 'z':
					normalize = true;
					break;
				case 'i':
					snrpoint = true;
					break;
				case 'f':
					fluxpoint = true;
					break;
				case 'N':
					noprintorder = true;
					break;
				case 'O':
					getorder = true;
					break;
				case 'R':
					orders = true;
					break;
				case 'T':
					tilt = true;
					break;
				case 'l':
					saturation = atol(optarg);
					break;
				case 'G':
					imagefilename = optarg;
					break;
				case 'g':
					geomfilename = optarg;
					break;
				case 'w':
					wavefilename = optarg;
					break;
				case 'U':
					specfilename = optarg;
					break;
				case 'P':
					aperfilename = optarg;
					break;
				case 'c':
					count = atoi(optarg);
					break;
				case 'o':
					order = atoi(optarg);
					break;
				case 'W':
					wl = atof(optarg);
					break;
				case 'A':
					aperture = atof(optarg);
					break;
				case 'x':
					x = atof(optarg);
					break;
				case 'y':
					y = atof(optarg);
					break;
				case '1':
					naxis1 = atoi(optarg);
					break;
				case '2':
					naxis2 = atoi(optarg);
					break;
				case 'b':
					upperlowerbounds = atoi(optarg);
					break;
					
				case 'v':
					verbose = true;
					break;
				case 'p':
					plot = true;
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
		
		// given a wl, get the order it is in
		if (getorder && !wavefilename.empty()) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, wavefilename);
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned maxorder = spectralOrderVector.getMaxorder();
			bool found = false;
			if (verbose) {
				cout << "operaImageStats: getorder: wl= " << wl << endl;
			}
			for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
				double length = 6000.0;
				if (!geomfilename.empty()) {
					if (spectralOrder->gethasGeometry()) {
						operaGeometry *geometry = spectralOrder->getGeometry();
						length = geometry->CalculateAndSetOrderLength();
					}
				}
				if (spectralOrder->gethasWavelength()) {
					operaWavelength *wavelength = spectralOrder->getWavelength();
					wavelength->setDmin(0.0);
					wavelength->setDmax(length);
					double wlmin = wavelength->getinitialWavelength();
					double wlmax = wavelength->getfinalWavelength();
					if (verbose) {
						cout << "operaImageStats: order= " << myorder << " wlmin=" << wlmin << " wlmax=" << wlmax << " wl=" << wl << endl;
					}
					if (wl >= wlmin && wl <= wlmax) {
						cout << myorder << endl;
						found = true;
					}								
				}
			}
			exit(found?EXIT_SUCCESS:EXIT_FAILURE);
		}
		
		if (optind == argc) {
			throw operaException("filename not given ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		while (optind < argc) {
			string filename = argv[optind++];
			isbin = filename.find(".bin") != string::npos;
			if (saturated) {
				operaFITSImage fitsimage(filename, tfloat, READONLY, cNone);
				float *pixels = (float *)fitsimage.getpixels();
				unsigned n = fitsimage.getnpixels();
				unsigned count = 0;
				while (n--) {
					if (*pixels++ >= saturation) {
						count++;
					}
				}
				cout << count << ' ' << (int)(count*100/fitsimage.getnpixels()) << ' ';
				fitsimage.operaFITSImageClose();
			}
			if (median) {
				operaFITSImage fitsimage(filename, tfloat, READONLY, cNone);
				float med = operaArrayMedianQuick(fitsimage.getnpixels(), (float *)fitsimage.getpixels());
				cout << med << ' ';
				fitsimage.operaFITSImageClose();
			}
			if (peak) {
				operaFITSImage fitsimage(filename, tfloat, READONLY, cNone);
				float *pixels = (float *)fitsimage.getpixels();
				unsigned n = fitsimage.getnpixels();
				float max = 0.0;
				while (n--) {
					if (max < *pixels) {
						max = *pixels;
					}
					pixels++;
				}
				cout << max << ' ';
				fitsimage.operaFITSImageClose();
			}
			if (peakunsaturated) {
				operaFITSImage fitsimage(filename, tfloat, READONLY, cNone);
				float *pixels = (float *)fitsimage.getpixels();
				unsigned n = fitsimage.getnpixels();
				float max = 0.0;
				while (n--) {
					if (max < *pixels && *pixels < saturation) {
						max = *pixels;
					}
					pixels++;
				}
				cout << max << ' ' << (int)(max*100/saturation) << ' ';
				fitsimage.operaFITSImageClose();
			}
			if (peakunsaturatedorder) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				float maxflux = 0.0;
				float maxorderflux = 0.0;
				unsigned mymaxorder = 0;
				for (unsigned order=minorder; order <= maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							if (spectralElements->getFlux(k) > maxflux) {
								maxflux = spectralElements->getFlux(k);
							}
						}
						if (maxflux > maxorderflux) {
							maxorderflux = maxflux;
							mymaxorder = order;
						}
					}
					maxflux = 0.0;
				}
				cout << mymaxorder << ' ';
			}
			if (maxsnr) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				float maxsnr = -999.0;
				for (unsigned order=minorder; order <= maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
					if (spectralOrder->gethasSNR()) {
						if (verbose) {
							cout << order << ' ' << spectralOrder->getCenterSNR();
						}
						if (spectralOrder->getCenterSNR() > maxsnr) {
							maxsnr = spectralOrder->getCenterSNR();
						}
					}
				}
				cout << (int)maxsnr << ' ';
			}
			if (smoothedpeaksnr) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				float peak = -999.0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
					if (spectralOrder->gethasSNR()) {
						float smoothed = spectralOrder->getCentralSmoothedSNR(upperlowerbounds);
						peak = max(peak, smoothed);
						if (verbose) {
							cout << "Smoothed (" << upperlowerbounds << ") Smoothed SNR for  " << filename << " order " << order << " : " << smoothed << " peak " << peak << endl;
						}
					}
				}
				cout << (int)peak << ' ';
			}
			if (snr) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				if (!specfilename.empty()) {
					operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, specfilename);
				}
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasSNR()) {
						if (spectralOrder->gethasSpectralElements() && (order == 0 || myorder == order)) {
							operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
							// get the peak if we are normalizing
							unsigned skip = 0;
							float maxvalue = 1.0;
							if (count > 0) {
								skip = spectralElements->getnSpectralElements() / count;
							}
							if (normalize) {
								for (unsigned k=0; k<spectralElements->getnSpectralElements(); k+=skip) {
									if (spectralElements->getFluxSNR(k) > maxvalue)  {
										maxvalue = spectralElements->getFluxSNR(k);
									}
								}
								
							}
							for (unsigned k=0; k<spectralElements->getnSpectralElements(); k+=skip) {
								cout << (noprintorder?"":itos(myorder)) << ' ' << spectralElements->getwavelength(k) << ' ' << (normalize?spectralElements->getFluxSNR(k)/maxvalue:(int)spectralElements->getFluxSNR(k)) << endl;
							}
						}
					}
				}
			}
			// given an x,y image coordinate, find the SNR at that point
			if (snrpoint && !geomfilename.empty()) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, geomfilename);
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasSNR() && spectralOrder->gethasGeometry()) {
						if (spectralOrder->gethasSpectralElements()) {
							operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
							operaGeometry *geometry = spectralOrder->getGeometry();
							if (y > geometry->getYmin() && y < geometry->getYmax()) {
								Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
								float x0 = geometryPolynomial->Evaluate(y);
								if (fabs(x-x0) < aperture) {
									unsigned k = (unsigned)y; // an approximation, should be distance
									float wl = 0.0;
									if (!wavefilename.empty()) {
										operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, wavefilename);		// wavelength
										if (spectralOrder->gethasWavelength()) {
											operaWavelength *wavelength = spectralOrder->getWavelength();
											Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
											wl = wavelengthPolynomial->Evaluate(y);	
										}
									}
									cout << (noprintorder?"":itos(myorder)) << ' ' << spectralElements->getwavelength(k) << ' ' << (int)spectralElements->getFluxSNR(k) << ' ' << (int)wl << ' ';
									break;
								}								
							}
						}
					}
				}
			}
			// given an x,y image coordinate, find the flux at that point
			if (fluxpoint && !filename.empty()) {
				if (isbin) {
					unsigned short buffer;
					streampos pos = (streampos)((unsigned long)x * naxis1 + (unsigned long)y);
					ifstream inbin(filename.c_str(), ios::in | ios::binary);
					inbin.seekg(pos);
					inbin.read((char *)&buffer, sizeof(unsigned short));
					cout << buffer;
					inbin.close();
				} else {
					operaFITSImage inimage(filename, tfloat, READONLY);
					cout << (float)inimage[(unsigned)y][(unsigned)x];
					inimage.operaFITSImageClose();
				}
			}
			// given an x,y image coordinate, find the order at that point
			if (getorder && !geomfilename.empty()) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, geomfilename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasGeometry()) {
						operaGeometry *geometry = spectralOrder->getGeometry();
						if (y > geometry->getYmin() && y < geometry->getYmax()) {
							Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
							float x0 = geometryPolynomial->Evaluate(y);
							if (fabs(x-x0) < aperture) {
								cout << myorder << ' ';
								break;
							}								
						}
					}
				}
			}
			// give the orders xmin, ymin, wlmin, xmax, ymax wlmax values
			if (orders && !geomfilename.empty() && !wavefilename.empty()) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, geomfilename);	// geometry
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, wavefilename);		// wavelength
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasGeometry() && spectralOrder->gethasWavelength()) {
						operaGeometry *geometry = spectralOrder->getGeometry();
						operaWavelength *wavelength = spectralOrder->getWavelength();
						Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
						float xmin = geometryPolynomial->Evaluate(geometry->getYmin());
						float xmax = geometryPolynomial->Evaluate(geometry->getYmax());
						Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
						float wlmin = wavelengthPolynomial->Evaluate(geometry->getYmin());
						float wlmax = wavelengthPolynomial->Evaluate(geometry->getYmax());
						cout << myorder << ' ' << xmin << ' ' << (int)geometry->getYmin() << ' ' << (int)wlmin << ' ' << (int)xmax << ' ' << (int)geometry->getYmax() << ' ' << (int)wlmax << endl;
					}
				}
			}
			// get the aperture tilt
			if (tilt && !aperfilename.empty()) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, aperfilename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasExtractionApertures()) {
						cout << (double)spectralOrder->getTiltInDegrees().value << endl;
						break;
					}
				}
			}
			if (spectrum) {
				operaSpectralOrderVector spectralOrderVector;
				operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, filename);
				unsigned minorder = spectralOrderVector.getMinorder();
				unsigned maxorder = spectralOrderVector.getMaxorder();
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasSpectralElements() && (order == 0 || myorder == order)) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						// get the peak if we are normalizing
						unsigned skip = 0;
						float maxvalue = 0.0;
						if (count > 0) {
							skip = spectralElements->getnSpectralElements() / count;
						}
						if (normalize) {
							for (unsigned k=0; k<spectralElements->getnSpectralElements(); k+=skip) {
								if (spectralElements->getFlux(k) > maxvalue)  {
									maxvalue = spectralElements->getFlux(k);
								}
							}
						}
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k+=skip) {
							if (spectralElements->getFlux(k) < saturation) {
								cout << spectralElements->getwavelength(k) << ' ' << (normalize?spectralElements->getFlux(k)/maxvalue:spectralElements->getFlux(k)) << endl;
							}
						}
					}						
				}
			}
			cout << endl;
		}
	}
	catch (operaException e) {
		cerr << "operaimagestats: " << e.getFormattedMessage() << endl;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "operaimagestats: " << operaStrError(errno) << endl;
		exit(EXIT_FAILURE);
	}
	return EXIT_SUCCESS;
}

