/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPlot
 Version: 1.0
 Description: Plots a fits product into a PNG files.
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
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaEspadonsImage.h"
#include "libraries/operaPNG.h"
#include "libraries/operaLibCommon.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

/*! \file operaPlot.cpp */

using namespace std;

/*!
 * operaPlot
 * \author Doug Teeple
 * \brief Plots a fits product into a PNG file.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operaPlot [--product=<FITS Product filename> | --input=<FITS image filename>] | --spectrum=<spectrum filename> --output=<filename>.png --notwavelengthcorrected --unnormalized -[dvth]\n";
}

int main(int argc, char *argv[])
{
	int opt;
	string fontpath = "/usr/X11/lib/X11/fonts/TTF/luxisr.ttf";
	string productfilename, inputfilename, outputfilename, spectrumfilename, snrfilename;
	string basedon;
	operaFITSProduct *product = NULL;
	operaEspadonsImage *input = NULL;
	int wavelengthcorrected = 1;
	int normalized = 1;
	bool polarimetry = false;
    bool orthometric = false;
	
	bool debug=false, verbose=false, trace=false, plot=true;
	
	struct option longopts[] = {
		{"product",			1, NULL, 'P'},			// input product
		{"input",			1, NULL, 'i'},			// input epsadons image
		{"spectrum",		1, NULL, 's'},			// .s spectrum
		{"snr",				1, NULL, 'r'},			// snr filename
		{"output",			1, NULL, 'o'},			// output
		{"unnormalized",	0, NULL, 'n'},			// normalized
		{"notwavelengthcorrected",0, NULL, 'w'},	// wavelengthcorrected
		{"orthometric",     0, NULL, 'm'},          // orthometric projection
		
		{"verbose",			0, NULL, 'v'},
		{"plot",			0, NULL, 'p'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "P:o:i:s:r:nwmpvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'P':		// input product
					productfilename = optarg;
					break;
				case 'i':		// output filename
					inputfilename = optarg;
					break;
				case 'o':		// output filename
					outputfilename = optarg;
					break;
				case 's':		// output filename
					spectrumfilename = optarg;
					break;
				case 'r':		// snr filename
					snrfilename = optarg;
					break;
				case 'n':
					normalized = 0;
					break;
				case 'w':
					wavelengthcorrected = 0;
					break;
				case 'm':
					orthometric = true;
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
		
		if (productfilename.empty() && inputfilename.empty() && spectrumfilename.empty()) {
			throw operaException("operaPlot: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        /*
         * PART I - PLOT PRODUCT SPECTRA
         */
		if (!spectrumfilename.empty() || !productfilename.empty()) {
			if (!spectrumfilename.empty()) {
				if (spectrumfilename.find(".s") == string::npos)
					throw operaException("operaPlot: spectrum must end in .s ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
				
				product = new operaFITSProduct();
				product->readText(spectrumfilename);
				if (outputfilename.empty()) {
					outputfilename = spectrumfilename.substr(0,spectrumfilename.find(".s")) + ".png";
				}
			}
			if (!productfilename.empty()) {
				if (productfilename.find("p.fits") != string::npos) {
					polarimetry = true;
					if (outputfilename.empty()) {
						outputfilename = spectrumfilename.substr(0,spectrumfilename.find(".fits")) + ".png";
					}
				} else if (productfilename.find("i.fits") == string::npos) {
					if (outputfilename.empty()) {
						if (outputfilename.empty()) {
							outputfilename = spectrumfilename.substr(0,spectrumfilename.find(".fits")) + ".png";
						}
					}
				}
				product = new operaFITSProduct(productfilename, READONLY);
			}
			int lines = product->getrows();
			string mode;
			switch (product->getmode()) {
				case MODE_STAR_ONLY:
					mode = "STAR_ONLY";
					break;
				case MODE_STAR_PLUS_SKY:
					mode = "STAR_PLUS_SKY";
					break;
				case MODE_POLAR:
					mode = "POLAR";
					break;
				default:
					break;
			}
			/*
			 * Color table we are trying to reproduce...
			 *
			 *    Wavelength (nm)
			 *    Red 		622 - 780
			 *    Orange 	597 - 622
			 *    Yellow 	577 - 597
			 *    Green 	492 - 577
			 *    Blue 		455 - 492
			 *    Violet 	390 - 455
			 */
			string title;
			if (!productfilename.empty()) 
				try {
					string object = product->operaFITSGetHeaderValue("OBJECT");
					title = object;
					title += " " + mode + " ";
					if (normalized)
						title += " normalized ";
					if (wavelengthcorrected)
						title += " wavelength corrected ";
					title += product->operaFITSGetHeaderValue("HSTTIME") + " ";
				} catch (...) {
					// do nothing
				}
			int wavelength = 0;
			int intensity = 1;
			int errorbar = 2;
			// for polarimetry
			int stokes = 0;
			int checkN1 = 0;
			int checkN2 = 0;
			if (wavelengthcorrected) {
				if (normalized) {
					wavelength = 0;
					intensity = 1;
					errorbar = 2;
					if (polarimetry) {
						stokes = 3;
						checkN1 = 4;
						checkN2 = 5;
						errorbar = 6;
					}
				} else {
					wavelength = 3;
					intensity = 4;
					errorbar = 5;
					if (polarimetry) {
						wavelength = 7;
						intensity = 8;
						stokes = 9;
						checkN1 = 10;
						checkN2 = 11;
						errorbar = 12;
					}
				}
			} else {
				if (normalized) {
					wavelength = 6;
					intensity = 7;
					errorbar = 8;
					if (polarimetry) {
						wavelength = 13;
						intensity = 14;
						stokes = 15;
						checkN1 = 16;
						checkN2 = 17;
						errorbar = 18;
					}
				} else {
					wavelength = 9;
					intensity = 10;
					errorbar = 11;
					if (polarimetry) {
						wavelength = 19;
						intensity = 20;
						stokes = 21;
						checkN1 = 22;
						checkN2 = 23;
						errorbar = 24;
					}
				}
			}
			const int plotwidth = 10000;
			const int plotheight = 6000;
			const int x_margin = (polarimetry?600:500);
			const int y_margin = 500;
			const int x_min = x_margin;
			const int x_max = plotwidth - x_margin;
			const int y_min = y_margin;
			const int y_max = plotheight - y_margin;
			const int x_ticks = 20;
			const int y_ticks = 10;
			const double backgroundcolour =  1.0;
			double r =  1.0;
			double g =  0.0;
			double b =  0.0;
			operaPNG image(outputfilename.c_str(), plotwidth, plotheight, backgroundcolour);
			// draw axes in black
			const int x_boundary_width = 20;
			const int y_boundary_width = 20;
			r =  0.0;
			g =  0.0;
			b =  0.0;
			image.filledsquare(x_min, y_min, x_max, y_min+y_boundary_width, r, g, b);		// bottom
			image.filledsquare(x_min, y_min, x_min+x_boundary_width, y_max, r, g, b);		// left
			image.filledsquare(x_max, y_min, x_max+x_boundary_width, y_max, r, g, b);		// right
			image.filledsquare(x_min, y_max, x_max, y_max+y_boundary_width, r, g, b);		// top
			for (int i=x_min; i<x_max; i+=(x_max-x_min)/x_ticks) {
				image.filledsquare(i, y_min, i+x_boundary_width, y_min+(y_max-y_min)/60, r, g, b);			// bottom
			}
			for (int i=y_min; i<y_max; i+=(y_max-y_min)/y_ticks) {
				image.filledsquare(x_min, i, x_min+(x_max-x_min)/120, i+y_boundary_width, r, g, b);	// left
			}
			int fontsize = 120;
			// find the min and max wavelength / intensity
			float minwavelength = 1.0e20;
			float maxwavelength = 0.0;
			float minintensity = 1.0e20;
			float maxintensity = 0.0;
			for (int l=0; l<lines; l++) {
				if ((float)*product[wavelength][l] < minwavelength) {
					minwavelength = *product[wavelength][l];
				};
				if ((float)*product[wavelength][l] > maxwavelength) {
					maxwavelength = *product[wavelength][l];
				};
				if ((float)*product[intensity][l] < minintensity) {
					minintensity = *product[intensity][l];
				};
				if ((float)*product[intensity][l] > maxintensity) {
					maxintensity = *product[intensity][l];
				};
			}
			int skip = lines / (x_max - x_min);
			float sincostick = (PI / (float)(x_max - x_min)*1.3);
			float angle = 0.0;
			float intensityrange = maxintensity - minintensity;
			int y1 = y_min;
			int y2 = y_min;
			int ebar = 0;
			int stoke = 0;
			//int checkn1 = 0;
			//int checkn2 = 0;
			int x1 = 0;
			r =  0.0;
			g =  0.0;
			b =  0.0;
			float er = 1.0;
			float eg = 0.0;
			float eb = 0.0;
			if (verbose) {
				cout << title << endl;
				if (polarimetry) {
					cout << "polar input data" << endl;
				}
				cout << "min intensity = " << minintensity << " max intensity = " << maxintensity << endl;
				cout << "min wavelength = " << minwavelength << " max wavelength = " << maxwavelength << endl;
				cout << "intensity range = " << intensityrange << endl;
				cout << "lines = " << lines << endl;
			}
			for (int x=x_min; x<x_max; x++) {
				y2 = (int)((*product[intensity][x1]+fabs(minintensity)) / intensityrange) * (y_max - y_min) + y_min;
				ebar = (int)((*product[errorbar][x1]) / intensityrange) * (y_max - y_min) + y_min;
				image.line(x, y1, x+1, y2, r, g, b);
				if (polarimetry) {
					stoke = (int)((*product[stokes][x1]+fabs(minintensity)) / intensityrange) * (y_max - y_min) + y_min;
					image.filledsquare(x, stoke, x, stoke+3, 0.0, 0.0, 1.0);
					//checkn1 = ((*product[checkN1][x1]+fabs(minintensity)) / intensityrange) * (y_max - y_min) + y_min;
					//image.filledsquare(x, checkn1, x, checkn1+3, 0.0, 1.0, 0.0);
					//checkn2 = ((*product[checkN2][x1]+fabs(minintensity)) / intensityrange) * (y_max - y_min) + y_min;
					//image.filledsquare(x, checkn2, x, checkn2+3, 1.0, 0.0, 0.0);
				}
				image.filledsquare(x, y2+ebar/2, x, y2+ebar/2+3, er, eg, eb);
				y1 = y2;
				x1 += skip;
				if (x < (x_max - x_min) / 2) {
					b = cos(angle);
					g = sin(angle);
				} else {
					g = sin(angle);
					r = cos(angle+PI);
				}
				angle += sincostick;
			}
			r = 0.0;
			g = 0.0;
			b = 0.0;
			int textwidth = image.gettextwidth((char *)fontpath.c_str(), fontsize, (char *)title.c_str());
			image.plottext((char *)fontpath.c_str(), fontsize, (plotwidth-y_margin*2-textwidth)/2, plotheight-y_margin/2, 0.0, (char *)title.c_str(), r, g, b);
			// label the axes
			char w_char[128];
			x1 = 0;
			fontsize = 100;
			for (int i=x_min; i<x_max; i+=(x_max-x_min)/x_ticks) {
				sprintf(w_char, "%.0f", *product[wavelength][x1]);
				image.plottext((char *)fontpath.c_str(), fontsize, i, y_min-y_min/2, 0.0, (char *)w_char, r, g, b);
				x1 += skip*(x_max-x_min)/x_ticks;
			}
			sprintf(w_char, "nm");
			image.plottext((char *)fontpath.c_str(), fontsize, (x_max-x_min)/2, y_min-y_min*3/4, 0.0, (char *)w_char, r, g, b);
			float theintensity = minintensity;
			float intensityincr = intensityrange/y_ticks;
			for (int i=y_min; i<y_max; i+=(y_max-y_min)/y_ticks) {
				sprintf(w_char, "%.2f", theintensity);
				image.plottext((char *)fontpath.c_str(), fontsize, x_min-x_min*3/4, i, 0.0, (char *)w_char, r, g, b);
				theintensity += intensityincr;
			}
			//  image.resize(plotwidth, plotheight);
			image.close();
			delete product;
		}
		/*
		 * PART II - PLOT SNR table
		 ***
		 40 1
		 372.000000 6.000000
		 377.000000 8.000000
		 */
		if (!snrfilename.empty()) {
			if (snrfilename.find(".sn") == string::npos)
				throw operaException("operaPlot: snr filename must end in .sn", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
			if (outputfilename.empty()) {
				outputfilename = spectrumfilename.substr(0,spectrumfilename.find(".sn")) + ".png";
			}
		}
		/*
		 * PART III - PLOT A REGULAR IMAGE
		 */
		if (!inputfilename.empty()) {
			input = new operaEspadonsImage(inputfilename, tfloat, READONLY);
			if (outputfilename.empty()) {
				outputfilename = spectrumfilename.substr(0,spectrumfilename.find(".fits")) + ".png";
			}
			const int x_boundary_width = 20;
			const int y_boundary_width = 20;
			const int x_margin = 400;
			const int y_margin = 400;
			const int plotwidth = (orthometric?2080*6+x_margin*2+x_boundary_width*2:2080+x_margin*2+x_boundary_width*2);
			const int plotheight = (orthometric?4640*2+y_margin*2+y_boundary_width*2:4640+y_margin*2+y_boundary_width*2);
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
			image.filledsquare(x_min, y_min, x_max, y_min+y_boundary_width, r, g, b);		// bottom
			image.filledsquare(x_min, y_min, x_min+x_boundary_width, y_max, r, g, b);		// left
			image.filledsquare(x_max, y_min, x_max+x_boundary_width, y_max+y_boundary_width, r, g, b);		// right
			image.filledsquare(x_min, y_max, x_max, y_max+y_boundary_width, r, g, b);		// top
			int fontsize = 120;
			string title;
			string object = input->operaFITSGetHeaderValue("OBJECT");
			title = object;
			//title += product->operaFITSGetHeaderValue("HSTTIME") + " ";
			operaFITSSubImage *datasec = input->getDatasecSubImage();
			float median2 = operaArrayMedian(input->getndpixels(), (float *)datasec->getpixels()) * 2.0;
			float maxpixel = operaArrayMaxValue(input->getndpixels(), (float *)datasec->getpixels());
			float minpixel = operaArrayMinValue(input->getndpixels(), (float *)datasec->getpixels());
			if (verbose) {
				cout << title << endl;
				cout << "nx= " << input->getnx() << " ny= " << input->getny() << '\n';
				cout << "naxis1= " << input->getnaxis1() << " naxis2= " << input->getnaxis2() << '\n';
				cout << "maxpixel= " << maxpixel << " median x 2.0= " << median2 << '\n';
			}
			unsigned x1 = input->getnaxis1() - input->getnx();  // skip overscan
			unsigned y1 = input->getnaxis2() - input->getny();  // skip overscan
			unsigned x2 = input->getnx();
			unsigned y2 = input->getny();
			if (orthometric) {
				unsigned xp = 0;
				for (unsigned y=y1; y<y2; y+=20) {
					int yval = 0;
					int yl = 0;
					for (unsigned x=x1; x<x2; x++) {
						float pixelvalue = *datasec[y][x] - minpixel;
						if (pixelvalue > median2) {
							b = 0.0;
							r = g = 1.0;
						} else {
							b = pixelvalue / median2;
							r = g = 1.0 - b;
						}
						pixelvalue = (*datasec[y][x] - minpixel) / median2 * 800.0;
						yval = (int)pixelvalue;
						if (debug) {
							cout << " " << pixelvalue << " r=" << r << " g=" << g << " b=" << b;
						}
						image.line(x+xp+x_margin+x_boundary_width, y+yl+y_margin+y_boundary_width, x+xp+y_margin+x_boundary_width, y+yval+y_margin+y_boundary_width, r, g, b);
						yl = yval;
					}
					xp+=36;
					if (debug) {
						cout << endl;
					}
				}
			} else {    // straight image
				for (unsigned x=x1; x<x2; x++) {
					for (unsigned y=y1; y<y2; y++) {
						float pixelvalue = *datasec[y][x] - minpixel;
						if (pixelvalue > median2) {
							b = 0.0;
							r = g = 1.0;
						} else {
							b = pixelvalue / median2;
							r = g = 1.0 - b;
						}
						if (debug) {
							cout << " " << pixelvalue << " r=" << r << " g=" << g << " b=" << b;
						}
						image.line(x+y_margin+x_boundary_width, y+y_margin+y_boundary_width, x+y_margin+x_boundary_width, y+y_margin+y_boundary_width, r, g, b);
					}
					if (debug) {
						cout << endl;
					}
				}
			}
			r = 0.0;
			g = 0.0;
			b = 0.0;
			int textwidth = image.gettextwidth((char *)fontpath.c_str(), fontsize, (char *)title.c_str());
			image.plottext((char *)fontpath.c_str(), fontsize, (plotwidth-y_margin*2-textwidth)/2, plotheight-y_margin/2, 0.0, (char *)title.c_str(), r, g, b);
#if 0
			// label the axes
			char w_char[128];
			x1 = 0;
			fontsize = 100;
			for (int i=x_min; i<x_max; i+=(x_max-x_min)/x_ticks) {
				sprintf(w_char, "%.0f", *product[wavelength][x1]);
				image.plottext((char *)fontpath.c_str(), fontsize, i, y_min-y_min/2, 0.0, (char *)w_char, r, g, b);
				x1 += skip*(x_max-x_min)/x_ticks;
			}
			sprintf(w_char, "nm");
			image.plottext((char *)fontpath.c_str(), fontsize, (x_max-x_min)/2, y_min-y_min*3/4, 0.0, (char *)w_char, r, g, b);
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
			delete input;
		}
	}
	catch (operaException e) {
		cerr << "operaPlot: " << e.getFormattedMessage() << '\n';
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPlot: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

