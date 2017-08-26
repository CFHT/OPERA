/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFITSDisplayImage
 Version: 1.0
 Description: Creates a FITS image of the SNR or spectrum values.
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
#include <algorithm>	// for min

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaEspadonsImage.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaStats.h"
#include "libraries/Polynomial.h"

/*! \file operaFITSDisplayImage.cpp */

using namespace std;

/*!
 * operaFITSDisplayImage
 * \author Doug Teeple
 * \brief Creates a FITS image of the SNR or spectrum values.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operaFITSDisplayImage [--snr | --spectrum ] --input=<filename> --output=<filename> --order=${spectralorder} --mode=${mode} --gapbetweenbeams=${gapbetweenbeams} --geomfilename=<filename> --wavefilename=<filename> -[dpvth]\n";
}

static inline int sgn(int i) {
	return (i<0?-1:1);
}

static inline int min(int x, int y) {
	if (x<y) {
		return x;
	}
	return y;
}
/*
 * draw a line - Bresenham's line algorithm
 */
static void drawline(operaFITSImage &image, unsigned short intensity, unsigned short x1, unsigned short y1, unsigned short x2, unsigned short y2) { 
	int n, deltax, deltay, sgndeltax, sgndeltay, deltaxabs, deltayabs, x, y, drawx, drawy; 
	
	deltax = x2 - x1; 
	deltay = y2 - y1; 
	deltaxabs = abs(deltax); 
	deltayabs = abs(deltay); 
	sgndeltax = sgn(deltax); 
	sgndeltay = sgn(deltay); 
	x = deltayabs >> 1; 
	y = deltaxabs >> 1; 
	drawx = x1; 
	drawy = y1; 
	
	image.setpixel(intensity, drawx, drawy); 
	
	if (deltaxabs >= deltayabs) { 
		for(n = 0; n < deltaxabs; n++){ 
			y += deltayabs; 
			if (y >= deltaxabs) { 
				y -= deltaxabs; 
				drawy += sgndeltay; 
			} 
			drawx += sgndeltax; 
			image.setpixel(intensity, drawx, drawy); 
		} 
	} else { 
		for(n = 0; n < deltayabs; n++) { 
			x += deltaxabs; 
			if (x >= deltayabs) { 
				x -= deltayabs; 
				drawx += sgndeltax; 
			} 
			drawy += sgndeltay; 
			image.setpixel(intensity, drawx, drawy); 
		} 
	} 
}

int main(int argc, char *argv[])
{
	int opt;
	
	string inputfilename;
	string outputfilename;
	string geomfilename;
	string wavefilename;
	
	string mode = "sp2";
	
	float aperture = 26.0;			// pixels
	float gapbetweenbeams = 4.0;	// pixels
	float tiltangle = 0.0;
	
	unsigned order = 0;
	unsigned count = 2080;
	bool spectrum = false;
    bool snr = false;
	
	bool debug=false, verbose=false, trace=false, plot = false;
	
	struct option longopts[] = {
		{"input",			1, NULL, 'i'},			// input epsadons image
		{"output",			1, NULL, 'o'},			// output
		{"spectrum",		0, NULL, 's'},			// bool
		{"snr",				0, NULL, 'r'},			// bool
		{"order",			1, NULL, 'e'},			// a particular order
		{"geomfilename",	1, NULL, 'g'},
		{"wavefilename",	1, NULL, 'w'},
		{"mode",			1, NULL, 'm'},
		{"count",			1, NULL, 'c'},
		{"aperture",		1, NULL, 'A'},
		{"gapbetweenbeams",	1, NULL, 'G'},
		{"tilt",			1, NULL, 'T'},
		
		{"plot",			0, NULL, 'p'},
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "o:i:s:r:e:g:w:m:A:G:T:c:vdpth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// output filename
					inputfilename = optarg;
					break;
				case 'o':		// output filename
					outputfilename = optarg;
					break;
				case 's':		// output filename
					spectrum = true;
					break;
				case 'r':		// snr filename
					snr = true;
					break;
				case 'e':
					order = (unsigned)atoi(optarg);
					break;
				case 'g':
					geomfilename = optarg;
					break;
				case 'w':
					wavefilename = optarg;
					break;
				case 'c':
					count = atoi(optarg);
					break;
				case 'T':
					tiltangle = atof(optarg);
					break;
				case 'm':
					mode = optarg;
					break;
				case 'A':
					break;
					aperture = atof(optarg);
				case 'G':
					gapbetweenbeams = atof(optarg);
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
		
		if (inputfilename.empty() || outputfilename.empty()) {
			throw operaException("operaFITSDisplayImage: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		/*
		 * put the output in a temp so that quicklook doesn't see an empty file
		 */
		string basefilename = outputfilename;
		string directory = "./";
		if (outputfilename.find_last_of("/") != string::npos) {
			basefilename = outputfilename.substr(outputfilename.find_last_of("/")+1);
			directory = outputfilename.substr(0, outputfilename.find_last_of("/"));
		}
		string tempoutfilename = directory + '.' + basefilename;
		if (mode == "sp2") {
			gapbetweenbeams = 0.0;
		}
		operaFITSImage out(tempoutfilename, ESPADONS_DEFAULT_NAXIS1, ESPADONS_DEFAULT_NAXIS2, tushort, cNone);
		if (spectrum) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, inputfilename);
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, geomfilename);		// geometry
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, wavefilename);		// wavelength
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned minactualorder = minorder;
			unsigned maxorder = spectralOrderVector.getMaxorder();
			unsigned maxactualorder = maxorder;
			unsigned lastx = 0;
			unsigned lasty = 0;
			unsigned x = 0;
			unsigned totalx = 0;
			unsigned short maxvalue = 0;
			unsigned short minvalue = 65535;
			unsigned short color = 0;
			for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
				if (spectralOrder->gethasSpectralElements()) {
					if (minactualorder == minorder) {
						minactualorder = myorder;
					}
					maxactualorder = myorder;
					if (order == 0 || order == myorder) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						// get the total number of elements in x direction
						totalx += spectralElements->getnSpectralElements();
						// get the peak for scaling
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							if ((unsigned short)spectralElements->getFlux(k) > maxvalue)  {
								maxvalue = (unsigned short)spectralElements->getFlux(k);
							}
							if ((unsigned short)spectralElements->getFlux(k) < minvalue)  {
								minvalue = (unsigned short)spectralElements->getFlux(k);
							}
						}
					}
				}
			}
			unsigned orders = maxactualorder-minactualorder+1;
			if (verbose) {
				cout << "operaFITSDisplayImage: order="<< order << " count=" << count << " totalx=" << totalx << " maxvalue=" << maxvalue << " color=" << color << endl;
			}
			if (order == 0) {	// all orders
				unsigned lasty = 0;
				float scaley = (float)(ESPADONS_DEFAULT_NAXIS2-1) / (float)(maxvalue-minvalue);
				float scalex = (float)(ESPADONS_DEFAULT_NAXIS1-1) / (float)totalx;
				unsigned skip = 1;
				if (count > 0) {
					skip = totalx / count;
				}
				float fskip = (float)skip;
				fskip *= scalex;
				skip = (unsigned)fskip;
				unsigned xstop = totalx / ESPADONS_DEFAULT_NAXIS1 / skip;
				if (verbose) {
					cout << "operaFITSDisplayImage: skip=" << skip << " scalex=" << scalex << " scaley=" << scaley << " xstop=" << xstop << endl;
				}
				for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
					operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						if (verbose) {
							cout << "operaFITSDisplayImage: order=" << myorder << " lastx=" << lastx << " color=" << color << " xstop+lastx=" << (xstop+lastx) << " skip=" << skip << endl;
						}
						color = (ESPADONS_DEFAULT_NAXIS1/orders)*(orders-myorder);
						unsigned orderx = 0;
						for (x=lastx; x+skip<xstop+lastx; x+=skip) {
							float fy = lasty;
							for (unsigned kk=orderx; kk<orderx+xstop; kk++) {
								fy += spectralElements->getFlux(kk) - minvalue;
							}
							fy /= xstop;	// take the mean
							unsigned y = min(ESPADONS_DEFAULT_NAXIS2-1, (unsigned)(fy * scaley));
							if (debug) {
								cout << "operaFITSDisplayImage: x=" << x << " orderx=" << orderx << " lastx=" << lastx << " lasty=" << lasty << " xstop=" << xstop << " y=" << y <<" lasty=" << lasty << endl;
							}
							drawline(out, color, x, lasty, x+skip, y);
							lasty = y;
							color -= skip;
							orderx += xstop * skip;
						}
						lastx = x;
					}
				}
			} else {	// a particular order
				float scaley = (float)(ESPADONS_DEFAULT_NAXIS2-1) / (float)(maxvalue-minvalue);
				float scalex = (float)(ESPADONS_DEFAULT_NAXIS1-1) / (float)totalx;
				unsigned xstop = ESPADONS_DEFAULT_NAXIS1;
				unsigned skip = 1;
				if (count > 0) {
					skip = totalx / count;
				}
				float fskip = (float)skip;
				fskip *= scalex;
				skip = (unsigned)fskip;
				if (verbose) {
					cout << "operaFITSDisplayImage: skip="<< skip << " scalex=" << scalex << " scaley=" << scaley << " xstop=" << xstop << endl;
				}
				color = (ESPADONS_DEFAULT_NAXIS1/orders)*(orders-order);
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
					if (verbose) {
						cout << "operaFITSDisplayImage: order="<< order << " lastx=" << lastx << " color=" << color << " xstop+lastx=" << (xstop+lastx) << endl;
					}
					for (x=0; x+skip<xstop; x+=skip) {
						float fy = 0.0;
						for (unsigned kk=x; kk<x+skip; kk++) {
							fy += spectralElements->getFlux(kk) - minvalue;
						}
						fy /= skip;	// take the mean
						unsigned y = min(ESPADONS_DEFAULT_NAXIS2-1, (unsigned)(fy * scaley));
						if (debug) {
							cout << "operaFITSDisplayImage: x=" << x << " lastx=" << lastx << " x+skip=" << (x+skip) << " y=" << y <<" lasty=" << lasty << endl;
						}
						drawline(out, color, min(x,ESPADONS_DEFAULT_NAXIS1-1), lasty, min(x+skip,ESPADONS_DEFAULT_NAXIS1-1), y);
						lasty = y;
						color -= skip;
					}
				}
			}
		}
		if (snr) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, inputfilename);
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, geomfilename);		// geometry
			operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, wavefilename);		// wavelength
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned maxorder = spectralOrderVector.getMaxorder();
			for (unsigned myorder=minorder; myorder <= maxorder; myorder++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(myorder);
				if (spectralOrder->gethasSNR()) {
					operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
					if (spectralOrder->gethasGeometry() && spectralOrder->gethasWavelength()) {
						operaGeometry *geometry = spectralOrder->getGeometry();
						Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
						unsigned count = 0;
						for (unsigned y=(unsigned)geometry->getYmin(); y<(unsigned)geometry->getYmax(); y++) {
							if (count >= spectralElements->getnSpectralElements()) {
								break;
							}
							unsigned x = (unsigned)geometryPolynomial->Evaluate(count);
                            unsigned width = (unsigned)((aperture-gapbetweenbeams)/2);
							unsigned start = (unsigned)(x-aperture/2);
							unsigned stop = (unsigned)(start + width);
							drawline(out, (unsigned short)spectralElements->getFluxSNR(count), start, y, stop, y);
							start = (unsigned)(stop + gapbetweenbeams);
							stop = start + width;
							drawline(out, (unsigned short)spectralElements->getFluxSNR(count), start, y, stop, y);
							count++;
						}
					}
				}
			}
		}
		out.operaFITSImageSave();
		out.operaFITSImageClose();
		rename(tempoutfilename.c_str(), outputfilename.c_str());
	}
	catch (operaException e) {
		cerr << "operaFITSDisplayImage: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaFITSDisplayImage: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

