
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPlotOut.cpp
 Version: 1.0
 Description: extract geometry and wave information suitable for gnuplot.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA 
 Date: Dec/2012
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2012  Opera Pipeline team, Canada France Hawaii Telescope
 
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
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"
#include "libraries/operaWavelength.h"
#include "libraries/operaCCD.h"       // for MAXORDERS

/*! \file operaPlotOut.cpp */

using namespace std;

/*! 
 * operaPlotOut
 * \author Doug Teeple
 * \brief extract geometry, aperture and wave information suitable for gnuplot.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

#define NOTPROVIDED -999

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaPlotOut -[dvth] --geom=... [--wave=...] [--aper=...] [--ordernumber=<n>[,<n2>]*]*\n";
	cerr << " Generates an extended data file suitable for gnuplot.\n";
	cerr << " format: order x y wl\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
    int ordernumber = NOTPROVIDED;	

	string geomfilename;
	string wavefilename;
	string aperfilename;
	string outputfilename;
	unsigned orders[MAXORDERS];
	unsigned orderindex = 0;
	
	struct option longopts[] = {
		
		{"geom",			1, NULL, 'g'},
		{"wave",			1, NULL, 'w'},
		{"aper",			1, NULL, 'a'},
		{"output",			1, NULL, 'o'},
        {"ordernumber",		1, NULL, 'O'},
		
		{"plot",			optional_argument, NULL, 'p'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "g:w:a:o:O:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'g':
					geomfilename = optarg;
					break;
				case 'o':
					outputfilename = optarg;
					break;
				case 'w':
					wavefilename = optarg;
					break;
				case 'a':
					aperfilename = optarg;
					break;
				case 'O':
					if (strchr(optarg, ',') != NULL) {	// comma separated list of orders
						char ostr[1024];
						strncpy(ostr, optarg, sizeof(ostr));
						char *pch, *pstart = ostr;
						pch = strchr(ostr, ',');
						while (pch != NULL) {
							*pch = '\0';
							ordernumber = atoi(pstart);
							orders[orderindex++] = ordernumber;
							pstart = pch + 1;
							pch = strchr(pstart, ',');
						}
						// last entry...
						ordernumber = atoi(pstart);
						orders[orderindex++] = ordernumber;
					} else {	// a single order
						ordernumber = atoi(optarg);
						orders[orderindex++] = ordernumber;
					}
					break;
					
				case 'v':
					verbose = 1;
					break;
				case 'p':
					plot = 1;
					break;
				case 'd':
					debug = 1;
					break;
				case 't':
					trace = 1; 
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
		
		if (geomfilename.empty()) {
			throw operaException("operaPlotOut: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (outputfilename.empty()) {
			outputfilename = geomfilename.substr(0, geomfilename.find_last_of(".")) + ".bin";
		}
		
		if (verbose) {
			cout << "operaPlotOut: outputfilename = " << outputfilename << endl;
			cout << "operaPlotOut: geomfilename = " << geomfilename << endl;
			cout << "operaPlotOut: wavefilename = " << wavefilename << endl;
			cout << "operaPlotOut: aperfilename = " << aperfilename << endl;
		}
		
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, geomfilename);	// get the geometry
		
		if (!wavefilename.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, wavefilename); // read wavelength calibration reference first guess
		}
		if (!aperfilename.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, aperfilename); // read aperture
		}
		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();
		
		if (ordernumber == NOTPROVIDED) {
			for (unsigned i=0; i<(maxorder-minorder+1); i++) {
				orders[i] = minorder+i;
				orderindex = i;
			}
		}

		if (verbose) {
			cout << "operaPlotOut: orders = ";
			for (unsigned i = 0; i < orderindex; i++) {
				cout << orders[i] << ' ';
			}
			cout << endl;
		}
		// output format: order x y wl beam 0:(aper x aper y aper width aper height) beam 1:(aper x aper y aper width aper height)
		if (!outputfilename.empty()) {
			ofstream fout;
			fout.open(outputfilename.c_str());
			for (unsigned i = 0; i < orderindex; i++) {
				unsigned order = orders[i];
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasGeometry()) {
					operaGeometry *geometry = spectralOrder->getGeometry();
					operaWavelength *wavelength = spectralOrder->getWavelength();
					Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
					float miny = geometry->getYmin();
					float maxy = geometry->getYmax();
					for (float y = miny; y <= maxy; y+=1.0) {
						float x = geometryPolynomial->Evaluate(y);
						fout << order << ' ' << x << ' ' << y;
						if (spectralOrder->gethasWavelength()) {
							fout << ' ' << wavelength->evaluateWavelength(y);
						} else {
							fout << ' ' <<  0.0;
						}
						if (spectralOrder->gethasExtractionApertures()) {
							for(unsigned beam=0;beam<spectralOrder->getnumberOfBeams(); beam++) {
								operaExtractionAperture<Line> *beamAperture = spectralOrder->getExtractionApertures(beam);
								const Line *beamLineAperture = beamAperture->getApertureShape();
								fout << beamLineAperture->getMidPoint().getXcoord() << ' '
								<< beamLineAperture->getMidPoint().getYcoord() << ' '            
								<< beamLineAperture->getWidth() << ' '
								<< beamLineAperture->getLength() << ' ';
							}     
							fout << endl;
						} else {
							for(unsigned beam=0;beam<spectralOrder->getnumberOfBeams(); beam++) {
								fout << ' ' <<  0.0 << ' ' <<  0.0 << ' ' <<  0.0 << ' ' <<  0.0;
							}
						}
						
						fout << endl;
					}
					fout << endl;	// empty line between orders
				}
			} 
			fout.close();
		}
	}
	catch (operaException e) {
		cerr << "operaPlotOut: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaPlotOut: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

