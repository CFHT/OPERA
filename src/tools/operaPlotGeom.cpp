/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPlotGeom
 Version: 1.0
 Description: Plots a geometry output file.
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
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaWavelength.h"
#include "libraries/operaGeometry.h"
#include "libraries/Polynomial.h"
#include "libraries/operaPNG.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFit.h"
#include "libraries/operaLibCommon.h"

/*! \file operaPlotGeom.cpp */

using namespace std;

/*! 
 * operaPlotGeom
 * \author Doug Teeple
 * \brief Plots a geometry output file..
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operaPlotGeom [--wave=<wavelength filename>] --geom=<geometry filename> --output=<filename>.png -[dvth]\n";
}

int main(int argc, char *argv[])
{
	int opt;
	string fontpath = "/usr/X11/lib/X11/fonts/TTF/luxisr.ttf";
	string geometryfilename, wavelengthfilename, uncalibratedfilename, masterflatfilename, outputfilename;
	int firstorder = 18;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"geom",			1, NULL, 'g'},			// geometry polynomial coefficients 
		{"wave",			1, NULL, 'w'},			// wavelengh calibration polynomial coefficients 
		{"uncalibrated",	1, NULL, 'u'},			// uncalibrated...s an uncalibrated spectrum
		{"masterflat",		1, NULL, 'm'},			// masterflat...fits from which polynomials were derived
		{"output",			1, NULL, 'o'},			// output
		{"firstorder",		1, NULL, 'f'},			// firstorder
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"plot",			0, NULL, 'p'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "g:f:w:o:u:m:vdtph", longopts, NULL))  != -1) {
			switch (opt) {
				case 'g':		// input geometry polynomial with coefficients
					geometryfilename = optarg;
					break;
				case 'w':		// input wavelengh calibration polynomial with coefficients
					wavelengthfilename = optarg;
					break;
				case 'u':		// input spectrum
					uncalibratedfilename = optarg;
					break;
				case 'o':		// output
					outputfilename = optarg;
					break;
				case 'm':		// masterflat filename
					masterflatfilename = optarg;
					break;
				case 'f':		// firstorder
					firstorder = atoi(optarg);
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
		
		if (geometryfilename.empty()) {
			throw operaException("operaPlotGeom: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (outputfilename.empty()) {
			throw operaException("operaPlotGeom: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		const int x_boundary_width = 20;
		const int y_boundary_width = 20;
		const int x_margin = 400;
		const int y_margin = 400;
		const int plotwidth = MAXESPADONSX+x_margin*2+x_boundary_width*2;
		const int plotheight = MAXESPADONSY+y_margin*2+y_boundary_width*2;
		const int x_min = x_margin;
		const int x_max = plotwidth - x_margin;
		const int y_min = y_margin;
		const int y_max = plotheight - y_margin;
		const double backgroundcolour =  1.0;
		double r =  1.0;
		double g =  0.0;
		double b =  0.0;
		
		operaPNG image(outputfilename.c_str(), plotwidth, plotheight, backgroundcolour);
		// draw axes in black
		r =  0.0;
		g =  0.0;
		b =  0.0;
		image.filledsquare(x_min-x_boundary_width*2, y_min-y_boundary_width*2, x_max+x_boundary_width, y_min-y_boundary_width, r, g, b);		// bottom
		image.filledsquare(x_min-x_boundary_width*2, y_min-y_boundary_width*2, x_min-x_boundary_width, y_max, r, g, b);		// left
		image.filledsquare(x_max, y_min-y_boundary_width*2, x_max+x_boundary_width, y_max+y_boundary_width, r, g, b);		// right
		image.filledsquare(x_min-x_boundary_width*2, y_max, x_max, y_max+y_boundary_width, r, g, b);		// top
		int fontsize = 80;
		string title = "OPERA Geometry Plot ";
		
		if (!masterflatfilename.empty()) {
			operaFITSImage *input = new operaFITSImage(masterflatfilename, tfloat, READONLY);
			unsigned npixels = input->getnaxis1() * input->getnaxis2();
            float median2 = operaArrayMedian(npixels, (float *)input->getpixels()) * 2.0;
            float maxpixel = operaArrayMaxValue(npixels, (float *)input->getpixels());
            float minpixel = operaArrayMinValue(npixels, (float *)input->getpixels());
			if (verbose) {
				cout << title << endl;
                cout << "naxis1= " << input->getnaxis1() << " naxis2= " << input->getnaxis2() << '\n';
                cout << "minpixel= " << minpixel << " maxpixel= " << maxpixel << " median x 2.0= " << median2 << '\n';
			}
            unsigned x1 = 0;
            unsigned y1 = 0;
            unsigned x2 = input->getnaxis1();
            unsigned y2 = input->getnaxis2();
			for (unsigned x=x1; x<x2; x++) {
				for (unsigned y=y1; y<y2; y++) {
					float pixelvalue = *input[y][x] - minpixel;
					r = g = b = pixelvalue / median2;
					image.line(x+y_margin+x_boundary_width, y+y_margin+y_boundary_width, x+y_margin+x_boundary_width, y+y_margin+y_boundary_width, r, g, b);
				}
			}
			delete input;
		}
		
		if (!uncalibratedfilename.empty()) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, uncalibratedfilename);
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned maxorder = spectralOrderVector.getMaxorder();
			r =  0.0;
			g =  1.0;
			b =  0.0;
			for (unsigned order=minorder; order<=maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
 				for (unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
					double x = spectralElements->getphotoCenterX(i);
					double y = spectralElements->getphotoCenterY(i);
					image.filledsquare((int)(x+x_margin+x_boundary_width), (int)(y+y_margin+y_boundary_width), (int)(x+x_margin+x_boundary_width+2), (int)(y+y_margin+y_boundary_width+2), r, g, b);
				}
			}
		}
		if (!geometryfilename.empty()) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, geometryfilename);
			r =  0.0;
			g =  0.0;
			b =  1.0;
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned maxorder = spectralOrderVector.getMaxorder();
			for (unsigned order=minorder; order<=maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
				operaGeometry *geometry = spectralOrder->getGeometry();
				Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
				float x0 = geometryPolynomial->Evaluate(geometry->getYmin());
				float xp = x0;
				for (float y=geometry->getYmin(); y<=geometry->getYmax(); y+=1.0) {
					xp = geometryPolynomial->Evaluate(y);
					image.line((int)(xp+x_margin+x_boundary_width), (int)(y+y_margin+y_boundary_width), (int)(xp+x_margin+x_boundary_width), (int)(y+1+y_margin+y_boundary_width), r, g, b);
				}
				char w_char[128];
				sprintf(w_char, "%d", order);
				int textwidth = image.gettextwidth((char *)fontpath.c_str(), 15, (char *)w_char);
				image.plottext((char *)fontpath.c_str(), 15, (int)(x0-textwidth/2+x_margin+x_boundary_width), (int)(y_max-40), 0.0, (char *)w_char, 0.0, 0.0, 1.0);
				image.plottext((char *)fontpath.c_str(), 15, (int)(xp-textwidth/2+x_margin+x_boundary_width), (int)(y_margin+y_boundary_width-30), 0.0, (char *)w_char, 0.0, 0.0, 1.0);
			}
		}
		
		if (!geometryfilename.empty() && !wavelengthfilename.empty()) {
			operaSpectralOrderVector spectralOrderVector;
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, geometryfilename);
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, wavelengthfilename);
			r =  0.0;
			g =  1.0;
			b =  1.0;
			unsigned minorder = spectralOrderVector.getMinorder();
			unsigned maxorder = spectralOrderVector.getMaxorder();
			for (unsigned order=minorder; order<=maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
				spectralOrder->CalculateWavelengthSolution();
				operaWavelength *wavelength = spectralOrder->getWavelength();
				operaGeometry *geometry = spectralOrder->getGeometry();
				//Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
				float wl0 = wavelength->getinitialWavelength();
				float wlf = wavelength->getfinalWavelength();
				Polynomial *geometryPolynomial = geometry->getCenterPolynomial();
				float x0 = geometryPolynomial->Evaluate(geometry->getYmin());
				float xf = geometryPolynomial->Evaluate(geometry->getYmax());
				//for (float d=wavelength->getDmin(); d<=wavelength->getDmax(); d+=1.0) {
				//	wlf = wavelengthPolynomial->Evaluate(d);
				//image.line(xp+x_margin+x_boundary_width, y+y_margin+y_boundary_width, xp+x_margin+x_boundary_width, y+1+y_margin+y_boundary_width, r, g, b);
				//}
				char w_char[128];
				sprintf(w_char, "%.0f", wl0);
				int textwidth = image.gettextwidth((char *)fontpath.c_str(), 10, (char *)w_char);
				image.plottext((char *)fontpath.c_str(), 10, (int)(x0-textwidth/2+x_margin+x_boundary_width), (int)(y_margin+y_boundary_width/2), 0.0, (char *)w_char, r, g, b);
				sprintf(w_char, "%.0f", wlf);
				textwidth = image.gettextwidth((char *)fontpath.c_str(), 10, (char *)w_char);
				image.plottext((char *)fontpath.c_str(), 10, (int)(xf-textwidth/2+x_margin+x_boundary_width), (int)(y_max-y_boundary_width/2), 0.0, (char *)w_char, r, g, b);
			}
		}
		r = 0.0;
		g = 0.0;
		b = 0.0;
		int textwidth = image.gettextwidth((char *)fontpath.c_str(), fontsize, (char *)title.c_str());
		image.plottext((char *)fontpath.c_str(), fontsize, plotwidth-y_margin-textwidth, plotheight-y_margin/2, 0.0, (char *)title.c_str(), r, g, b);
#if 0
		// label the axes
		char w_char[128];
		x1 = 0;
		fontsize = 100;
		for (int i=x_min; i<x_max; i+=(x_max-x_min)/x_ticks) {
			sprintf(w_char, "%.0f", *input[wavelength][x1]);
			image.plottext((char *)fontpath.c_str(), fontsize, i, y_min-y_min/2, 0.0, (char *)w_char, r, g, b);
			x1 += skip*(x_max-x_min)/x_ticks;
		}
		sprintf(w_char, "nm");
		image.plot_text((char *)fontpath.c_str(), fontsize, (x_max-x_min)/2, y_min-y_min*3/4, 0.0, (char *)w_char, r, g, b);
		float theintensity = minintensity;
		float intensityincr = intensityrange/y_ticks;
		for (int i=y_min; i<y_max; i+=(y_max-y_min)/y_ticks) {
			sprintf(w_char, "%.2f", theintensity);
			image.plottext((char *)fontpath.c_str(), fontsize, x_min-x_min*3/4, i, 0.0, (char *)w_char, r, g, b);
			theintensity += intensityincr;
		}
#endif
		//  image.resize(plotwidth, plotheight);
		image.close();
	}
	catch (operaException e) {
		cerr << "operaPlotGeom: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPlotGeom: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

