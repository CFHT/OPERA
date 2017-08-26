/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPlotWave
 Version: 1.0
 Description: Plots a wavelength calibration output file.
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
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaFITSProduct.h"
#include "core-espadons/operaWavelengthCalibration.h"	// for MAXREFWAVELENGTHS
#include "libraries/operaPNG.h"
#include "libraries/operaLibCommon.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaFit.h"
#include "libraries/operaParameterAccess.h"

/*! \file operaPlotWave.cpp */

using namespace std;

/*! 
 * operaPlotWave
 * \author Doug Teeple
 * \brief Plots a wavelength calibration output file..
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout << " Usage: operaPlotWave --wave=<wavelength filename> [--tharref=...] --output=<filename>.png -[dvth]\n";
}

int main(int argc, char *argv[])
{
	int opt;
	string fontpath = "/usr/X11/lib/X11/fonts/TTF/luxisr.ttf";	// works on Linux and Mac
	string waveplotfilename, outputfilename, tharreffilename;
	int firstorder = 18;
	bool linesonly = false;
	
	bool debug=false, verbose=false, trace=false, plot=false;
	
	struct option longopts[] = {
		{"wave",			1, NULL, 'w'},			// wavelength datapoints 
		{"tharref",			1, NULL, 'r'},			// thorium argon reference spectrum 
		{"output",			1, NULL, 'o'},			// output
		{"firstorder",		1, NULL, 'f'},			// firstorder
		{"linesonly",		0, NULL, 'l'},			// just do the lines, not the spectra
		
		{"verbose",			0, NULL, 'v'},
		{"plot",			0, NULL, 'p'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "w:o:f:r:pvdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'w':		// input geometry with coefficients
					waveplotfilename = optarg;
					break;
				case 'o':		// output
					outputfilename = optarg;
					break;
				case 'f':		// firstorder
					firstorder = atoi(optarg);
					break;
				case 'r':		// tharreffilename
					tharreffilename = optarg;
					break;
				case 'l':
					linesonly = true;
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
		
		if (waveplotfilename.empty()) {
			throw operaException("operaPlotWave: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (outputfilename.empty()) {
			throw operaException("operaPlotWave: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		const int x_boundary_width = 20;
		const int y_boundary_width = 20;
		const int x_margin = 400;
		const int y_margin = 400;
		const int plotwidth = MAXESPADONSY+x_margin*2+x_boundary_width*2;
		const int plotheight = 5000+y_margin*2+y_boundary_width*2;
		const int x_min = x_margin;
		const int x_max = plotwidth - x_margin;
		const int y_min = y_margin;
		const int y_max = plotheight - y_margin;
		float minx, maxx;
		float wl0 = 0.0;
		float wlf = 0.0;
		const double backgroundcolour =  1.0;
		double r =  1.0;
		double g =  0.0;
		double b =  0.0;
		string dataline;
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
		string title = "OPERA WCAL Plot ";
		/*
		 * The wave file contains four sets of data
		 * First is the raw spectrum (distance and flux)
		 * Second is the detected lines (unscaled distance and flux)
		 * Third is the detected lines (scaled distance and flux)
		 * Fourth is the atlas lines (distance and flux)
		 * These should correlate
		 */
		if (!waveplotfilename.empty()) {
			ifstream flist(waveplotfilename.c_str());		
			if (flist.is_open()) {
				int order = 0;
				int lines = 0;
				int totalLines = 0;
				int cols = 0;
				float distance = 0.0;
				float flux = 0.0;
				float distance_last = 0.0;
				float flux_last = 0.0;
				int loops = 0;
				while (flist.good()) {
					getline (flist, dataline);
					if (strlen(dataline.c_str()))  {
						if (dataline.c_str()[0] == '#') {
							// skip comments
						} else {
							if (lines == 0) {
								sscanf(dataline.c_str(), "*** %d %d %d %f %f %f %f", &totalLines, &cols, &order, &minx, &maxx, &wl0, &wlf);
								loops++;
								if (loops == 1) {
									r =  0.0;
									g =  0.0;
									b =  1.0;
									const int legendfontsize = 60;
									const int legendleftmargin = 650;
									char buff[24];
									sprintf(buff, "Order %d, First Order %d", (int)order, firstorder);
									image.plottext((char *)fontpath.c_str(), fontsize, x_margin+x_boundary_width, plotheight-(y_margin+y_boundary_width)/2, 0.0, (char *)buff, r, g, b);
									
									image.plottext((char *)fontpath.c_str(), legendfontsize, x_max-legendleftmargin, y_max-120, 0.0, (char *)"Legend", 0.0,0.0,0.0);
									
									image.plottext((char *)fontpath.c_str(), legendfontsize, x_max-legendleftmargin, y_max-240, 0.0, (char *)"uncal flux", 0.0,0.0,0.0);
									image.line(x_max-legendleftmargin-100, y_max-200, x_max-legendleftmargin-200, y_max-200, 1.0,0.0,0.0);
									
									image.plottext((char *)fontpath.c_str(), legendfontsize, x_max-legendleftmargin, y_max-320, 0.0, (char *)"uncal lines", 0.0,0.0,0.0);
									image.line(x_max-legendleftmargin-100, y_max-280, x_max-legendleftmargin-200, y_max-280, 0.0,1.0,0.0);
									
									image.plottext((char *)fontpath.c_str(), legendfontsize, x_max-legendleftmargin, y_max-400, 0.0, (char *)"scaled lines", 0.0,0.0,0.0);
									image.line(x_max-legendleftmargin-100, y_max-360, x_max-legendleftmargin-200, y_max-360, 0.0,0.0,1.0);
									
									image.plottext((char *)fontpath.c_str(), legendfontsize, x_max-legendleftmargin, y_max-480, 0.0, (char *)"atlas", 0.0,0.0,0.0);
									image.line(x_max-legendleftmargin-100, y_max-440, x_max-legendleftmargin-200, y_max-440, 1.0,1.0,0.0);
								}
							} else {
								sscanf(dataline.c_str(), "%f %f", &distance, &flux);
								distance = (float)(x_max-x_min) * (distance-minx)/(maxx-minx);		// scale for plotting purposes
								flux *= (y_max-y_min);		// scale flux for plotting purposes
								if (loops == 1) {			// orange thar reference lines
									r =  255.0/255.0;
									g =  165.0/255.0;
									b =  0.0/255.0;
									image.line((int)(distance+x_margin+x_boundary_width), (int)(flux+y_margin+y_boundary_width), (int)(distance+x_margin+x_boundary_width), (int)(y_margin+y_boundary_width), r, g, b);
									if (lines == 1) {
										// label the axes
										char w_char[128];
										const int x_ticks = 10;
										int x1 = 0;
										float value = wl0;
										fontsize = 60;
										for (int i=x_min; i<x_max; i+=(x_max-x_min)/x_ticks) {
											sprintf(w_char, "%.0f", value);
											image.plottext((char *)fontpath.c_str(), fontsize, i, y_min-y_margin/2, 0.0, (char *)w_char, 0.0, 0.0, 0.0);
											x1 += (x_max-x_min)/x_ticks;
											value += (wlf-wl0)/x_ticks;
										}
									}
								} else if (loops == 2 && linesonly == false) {	// do red lines - raw spectrum
									if (lines == 1) {
										distance_last = distance;
										flux_last = flux;
									}
									r =  220.0/255.0;
									g =  20.0/255.0;
									b =  60.0/255.0;
									if (distance_last)
										image.line((int)(distance+x_margin+x_boundary_width), (int)(flux+y_margin+y_boundary_width), (int)(distance_last+x_margin+x_boundary_width), (int)(flux_last+y_margin+y_boundary_width), r, g, b);
								} else if (loops == 3) {	// do light blue vertical lines  - the raw detected lines in the spectrum
									r =  127.0/255.0;
									g =  255.0/255.0;
									b =  212.0/255.0;
									image.line((int)(distance+x_margin+x_boundary_width), (int)(flux+y_margin+y_boundary_width), (int)(distance+x_margin+x_boundary_width), (int)(y_margin+y_boundary_width), r, g, b);
								} else if (loops == 4) {	// do blue vertical lines - the scaled detected lines in the spectrum
									r =  70.0/255.0;
									g =  130.0/255.0;
									b =  180.0/255.0;
									image.line((int)(distance+x_margin+x_boundary_width), (int)(flux+y_margin+y_boundary_width), (int)(distance+x_margin+x_boundary_width), (int)(y_margin+y_boundary_width), r, g, b);
								}
							}
						}
						lines++;
						if (lines == totalLines+1) {
							lines = 0;
							distance_last = 0.0;
							flux_last = 0.0;
							if (loops == 4)
								break;
						} else {
							distance_last = distance;
							flux_last = flux;
						}
					}									
				}	
				flist.close();
			}
		}
		//
		// do the thorium argon reference if given in green
		//
		float minthrefflux = BIG;
		float maxthrefflux = 0.0;
		if (!tharreffilename.empty() && linesonly == false) {
			ifstream flist(tharreffilename.c_str());		
			if (flist.is_open()) {
				float wavelength = 0.0;
				float flux = 0.0;
				while (flist.good()) {
					getline (flist, dataline);
					if (strlen(dataline.c_str())) {
						if (dataline.c_str()[0] == '#') {
							// skip comments
						} else {
							sscanf(dataline.c_str(), "%f %f", &wavelength, &flux);
							wavelength *= 0.1;
							if (flux < minthrefflux) {
								minthrefflux = flux;
							}
							if (flux > maxthrefflux) {
								maxthrefflux = flux;
							}
						}
					}
				}
				flist.close();
			}
		}
		if (!tharreffilename.empty() && linesonly == false) {
			ifstream flist(tharreffilename.c_str());		
			if (flist.is_open()) {
				int lines = 0;
				float wavelength = 0.0;
				float flux = 0.0;
				float wavelength_last = 0.0;
				float flux_last = 0.0;
				while (flist.good()) {
					getline (flist, dataline);
					if (strlen(dataline.c_str())) {
						if (dataline.c_str()[0] == '#') {
							// skip comments
						} else {
							sscanf(dataline.c_str(), "%f %f", &wavelength, &flux);
							wavelength *= 0.1;
							if (wavelength >= wl0 && wavelength <= wlf) {
								flux = (float)(y_max-y_min) * flux / (maxthrefflux-minthrefflux);		// scale for plotting purposes
								wavelength = (float)(x_max-x_min) * (wavelength-wl0)/(wlf-wl0);		// scale for plotting purposes
								if (lines == 1) {
									wavelength_last = wavelength;
									flux_last = flux;
								}
								r =  20.0/255.0;
								g =  220.0/255.0;
								b =  60.0/255.0;
								if (wavelength_last)
									image.line((int)(wavelength+x_margin+x_boundary_width), (int)(flux+y_margin+y_boundary_width), (int)(wavelength_last+x_margin+x_boundary_width), (int)(flux_last+y_margin+y_boundary_width), r, g, b);
								lines++;
								wavelength_last = wavelength;
								flux_last = flux;							
							}
						}
					}
				}
				flist.close();
			}
		}
		//
		// now plot the title
		//
		r = 0.0;
		g = 0.0;
		b = 0.0;
		fontsize = 80;
		int textwidth = image.gettextwidth((char *)fontpath.c_str(), fontsize, (char *)title.c_str());
		image.plottext((char *)fontpath.c_str(), fontsize, (plotwidth-x_margin-textwidth)/2, plotheight-y_margin/2, 0.0, (char *)title.c_str(), r, g, b);
		// and label the axes
		// x
		image.plottext((char *)fontpath.c_str(), fontsize, (plotwidth-x_margin)/2, y_min-y_margin*3/4, 0.0, (char *)"Distance/WL(nm)", r, g, b);
		// y
		image.plottext((char *)fontpath.c_str(), fontsize, x_margin/4, (plotheight-y_margin)/2, 0.0, (char *)"Flux", r, g, b);
#if 0
		float theintensity = miny;
		float intensityincr = (maxy-miny)/y_ticks;
		const int y_ticks = 10;
		for (int i=y_min; i<y_max; i+=(y_max-y_min)/y_ticks) {
			sprintf(w_char, "%.2f", theintensity);
			image.plottext((char *)fontpath.c_str(), fontsize, x_min-x_min*3/4, i, 0.0, (char *)w_char, r, g, b);
			theintensity += intensityincr;
		}
#endif
		image.close();
	}
	catch (operaException e) {
		cerr << "operaPlotWave: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPlotWave: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

