/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaNormalize
 Version: 1.0
 Description: Normalize a spectrum 
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

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"

#include "core-espadons/operaNormalize.h"

#define NOTPROVIDED -999

/*! \file operaNormalize.cpp */
/*
 
 How the continuum is normalized:
 
 (i) a fit as a function of lambda (wavelength) over the whole spectrum
 
 (ii) a 2d fit of the shape of the order, in case there are residuals not
 compatible with the fit done in (i) - this would be the case at the
 edge of all the orders, where we have a drop in flux
 this is done with a low order polynomial fit, in both directions
 (along and perpendicular to the dispersion), therefore, only the
 low frequency variations are found and fit
 
 Donati says he chose to do it that way because of stability and robustness
 reasons.
 
 Donati also notes that this method will work in general, but that if one
 is interested in perfect normalization, then one would have to do a higher
 order polynomial fit on a reference star, and use that to normalize.
 
 Donati also mentions something that we already know, which is that for
 emission line spectra, it's better to NOT normalize. It's something I've
 told some PIs over the years.
 
 Something to use as a starting point, and as a point of comparison.
 
*/

using namespace std;

/*! 
 * operaNormalize
 * \author Doug Teeple
 * \brief Normalize a spectrum.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */
void GenerateNormalization3DBeamSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, unsigned numberOfBeams, bool plotContinuum, bool display);
void GenerateNormalization3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool plotContinuum, bool display);
void GenerateNormalizationSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool display);

int main(int argc, char *argv[])
{
	int opt;
	
	string inputSpectraFile;
	string outputSpectraFile;
	string spectrumtypename;
	string inputWaveFile;    
	
	operaSpectralOrder_t spectralOrderType = OptimalBeamSpectrum;
	
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
    
    unsigned binsize = 100;    
    bool usePolynomial = FALSE;
    unsigned orderOfPolynomial = 5;
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string continuumDataFilename;	
	string scriptfilename;	
    
    bool generate3DPlot = TRUE;
    bool generateBeamPlot = FALSE;
    bool plotContinuum = FALSE;
    
    
	struct option longopts[] = {
		{"inputSpectraFile",	1, NULL, 'i'},
		{"inputWaveFile",       1, NULL, 'w'},        
		{"outputSpectraFile",	1, NULL, 'o'},
		{"spectrumtype",		1, NULL, 'T'},	
		{"spectrumtypename",	1, NULL, 'N'},	
		{"ordernumber",			1, NULL, 'O'},	
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},               
		{"binsize",				1, NULL, 'b'},    
		{"usePolynomial",		1, NULL, 'l'},
		{"orderOfPolynomial",	1, NULL, 'r'},
		{"plotfilename",		1, NULL, 'P'},
        {"generate3DPlot",      1, NULL, 'E'},
        {"generateBeamPlot",    1, NULL, 'B'},
        {"plotContinuum",       1, NULL, 'c'},
		{"spectrumDataFilename",		1, NULL, 'F'},
		{"continuumDataFilename",		1, NULL, 'C'},        
		{"scriptfilename",		1, NULL, 'S'},  
		{"interactive",			0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:w:o:T:N:O:M:X:b:l:r:P:E:B:c:F:C:S:I:v::d::t::p::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputSpectraFile = optarg;	
				break;
			case 'w':
				inputWaveFile = optarg;
				break;
			case 'o':		// output
				outputSpectraFile = optarg;
				break;
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'N':		// spectrum type name for verbose mode
				spectrumtypename = optarg;
				break;
			case 'O':
				ordernumber = atoi(optarg);
				break;				
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;                                
			case 'b':		// binsize
				binsize = atoi(optarg);
				break;
			case 'l':		
				usePolynomial = (atoi(optarg)?true:false);
				break;
			case 'r':
				orderOfPolynomial = atoi(optarg);
				break;                
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break;
			case 'E':
				generate3DPlot = (atoi(optarg)?true:false);
				break;
			case 'B':
				generateBeamPlot = (atoi(optarg)?true:false);
				break;
			case 'c':
				plotContinuum = (atoi(optarg)?true:false);
				break;
			case 'F':
				spectrumDataFilename = optarg;
				break; 	
			case 'C':
				continuumDataFilename = optarg;
				break;                 
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
				interactive = (atoi(optarg)?true:false);
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
		// we need an input...
		if (inputSpectraFile.empty()) {
			throw operaException("operaNormalize: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputSpectraFile.empty()) {
			throw operaException("operaNormalize: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cerr << "operaNormalize: inputSpectraFile = " << inputSpectraFile << endl;
			cerr << "operaNormalize: inputWaveFile = " << inputWaveFile << endl;
			cerr << "operaNormalize: outputSpectraFile = " << outputSpectraFile << endl; 
			cerr << "operaNormalize: spectrumtype = " << spectralOrderType << endl; 
			cerr << "operaNormalize: spectrumtypename = " << spectrumtypename << endl; 
			cerr << "operaNormalize: ordernumber = " << ordernumber << endl; 
			cerr << "operaNormalize: binsize = " << binsize << endl;
			cerr << "operaNormalize: usePolynomial = " << usePolynomial << endl; 
			cerr << "operaNormalize: orderOfPolynomial = " << orderOfPolynomial << endl;
			cerr << "operaNormalize: generate3DPlot = " << generate3DPlot << endl;
			cerr << "operaNormalize: generateBeamPlot = " << generateBeamPlot << endl;
			cerr << "operaNormalize: plotContinuum = " << generateBeamPlot << endl;
            cerr << "operaNormalize: plotfilename = " << plotfilename << endl;
			cerr << "operaNormalize: spectrumDataFilename = " << spectrumDataFilename << endl; 
			cerr << "operaNormalize: continuumDataFilename = " << continuumDataFilename << endl; 
			cerr << "operaNormalize: scriptfilename = " << scriptfilename << endl; 
			cerr << "operaNormalize: interactive = " << interactive << endl;
		}
        ofstream *fspecdata = NULL;
        ofstream *fcontinuumdata = NULL;
        
        if (!spectrumDataFilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(spectrumDataFilename.c_str());  
        }    
        
        if (!continuumDataFilename.empty()) {
            fcontinuumdata = new ofstream();
            fcontinuumdata->open(continuumDataFilename.c_str());  
        }          

		operaSpectralOrderVector spectralOrders(inputSpectraFile);
        if (inputWaveFile.empty()) {
            spectralOrders.ReadIntoSpectralOrders(inputWaveFile);
        }
        if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        if(!maxorderprovided) {
            maxorder = spectralOrders.getMaxorder();            
        }        
        
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
        
		if (verbose)
			cerr << "operaNormalize: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        unsigned numberOfBeams = 0; // for plotting
        unsigned numberOfOrders = 0; // for plotting
        unsigned numberOfprintouts; // for plotting
        
        if(generate3DPlot) {
            numberOfprintouts= 2; // for 3D 
        } else {
            numberOfprintouts= 1; // for 2D 
        }
        
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if(order == minorder) {
                numberOfBeams = spectralOrder->getnumberOfBeams(); // for plotting
            }
            
            if (verbose)
                cerr << "operaNormalize: applying normalization to order " << order << endl;
            
            spectralOrder->applyNormalization(binsize,orderOfPolynomial,usePolynomial,fspecdata,fcontinuumdata,TRUE, numberOfprintouts);
            if (fspecdata != NULL) {
                *fspecdata << endl;
            }
            numberOfOrders++;
        }
        
		/* 
		 * and write it out
		 */
		spectralOrders.WriteSpectralOrders(outputSpectraFile, spectralOrderType);

        if(generate3DPlot) {
            if (fspecdata != NULL) {
                fspecdata->close();
                if (!scriptfilename.empty()) {
                    if(generateBeamPlot) {
                        GenerateNormalization3DBeamSpecPlot(scriptfilename,plotfilename,spectrumDataFilename, numberOfBeams, plotContinuum, interactive);
                    } else {
                        GenerateNormalization3DSpecPlot(scriptfilename,plotfilename,spectrumDataFilename, plotContinuum, interactive);
                    }
                }
            }
        } else {
            if (fspecdata != NULL) {
                fspecdata->close();
                if (!scriptfilename.empty()) {
                    GenerateNormalizationSpecPlot(scriptfilename,plotfilename,spectrumDataFilename, interactive);
                }
            }
        }
	}
	catch (operaException e) {
		cerr << "operaNormalize: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaNormalize: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputSpectraFile=<BEAMSPEC_FILE>"
	" --inputWaveFile=<WAVE_FILE>"    
	" --outputSpectraFile=<BEAMSPEC_FILE>"
	" --spectrumtype=<UNS_VALUE>"
	" --spectrumtypename=<STRING>"
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --binsize=<UNS_VALUE>"
	" --usePolynomial=<BOOL>"
	" --orderOfPolynomial=<UNS_VALUE>"
	" --generate3DPlot=<BOOL>"
	" --generateBeamPlot=<BOOL>"
	" --plotContinuum=<BOOL>"
	" --plotfilename=<EPS_FILE>"
	" --spectrumDataFilename=<DATA_FILE>"
	" --continuumDataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputSpectraFile=1599047.e.gz --outputSpectraFile=1599047.en.gz --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --ordernumber=-999 --usePolynomial=0 --spectrumDataFilename=normspec.dat --continuumDataFilename=normcontinuum.dat --scriptfilename=norm.gnu --spectrumDataFilename=norm.dat --plotfilename=norm.eps --binsize=110 --generate3DPlot=1 --generateBeamPlot=0 --plotContinuum=0 -v\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputSpectraFile=<BEAMSPEC_FILE>, Input beam spectrum file to normalize \n"
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file\n"    
	"  -o, --outputSpectraFile=<BEAMSPEC_FILE>, Output beam spectrum file \n"
	"  -T, --spectrumtype=<UNS_VALUE>, Option for spectrum type \n"
	"  -N, --spectrumtypename=<STRING>, Option for spectrum type \n"
	"  -O, --ordernumber=<UNS_VALUE>, Absolute order number to extract (default=all)\n"
	"  -M, --minorder=<UNS_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<UNS_VALUE>, Define maximum order number\n"
	"  -b, --binsize=<UNS_VALUE>, Number of points to bin for continuum estimate \n"
	"  -l, --usePolynomial=<BOOL>, Use polynomial function to fit continuum \n"
	"  -r, --orderOfPolynomial=<UNS_VALUE>, Order of polynomial to fit continuum \n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -E, --generate3DPlot=<BOOL>, Switch to generate 3D or 2D plot spectra\n"
	"  -B, --generateBeamPlot=<BOOL>, Switch to generate plot of beams or full slit spectra\n"
	"  -c, --plotContinuum=<BOOL>, Switch to generate plot of continuum or normalized line spectra\n"
	"  -F, --spectrumDataFilename=<DATA_FILE>\n"
	"  -C, --continuumDataFilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateNormalization3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool plotContinuum, bool display) {
    
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "set view 0,0" << endl;
    
    *fgnu << "set palette gray" << endl;
    *fgnu << "set palette gamma 2.0" << endl;
    *fgnu << "set pm3d map" << endl;
    *fgnu << "unset ztics" << endl;
    
    *fgnu << "set xrange[-200:*]" << endl;
    
    *fgnu << "\nset xlabel \"distance (pixels)\"" << endl;
    *fgnu << "set ylabel \"order number\"" << endl;
    
    unsigned fluxColumn;
    if(plotContinuum) {
        *fgnu << "set cblabel \"continuum flux\"" << endl;        
        fluxColumn = 9;
    } else {
        *fgnu << "set cblabel \"normalized flux\"" << endl;
        *fgnu << "set zrange[0:1]" << endl;
        fluxColumn = 10;
    }
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.8*$1 - 0.4):" << fluxColumn <<" w pm3d" << endl;
        *fgnu << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.8*$1 - 0.4):" << fluxColumn <<" w pm3d" << endl;
        
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}

void GenerateNormalization3DBeamSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, unsigned numberOfBeams, bool plotContinuum, bool display) {
    
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "set view 0,0" << endl;
    
    *fgnu << "set palette gray" << endl;
    *fgnu << "set palette gamma 2.0" << endl;
    *fgnu << "set pm3d map" << endl;
    *fgnu << "unset ztics" << endl;
    
    *fgnu << "set xrange[-200:*]" << endl;
    
    *fgnu << "\nset xlabel \"distance (pixels)\"" << endl;
    *fgnu << "set ylabel \"order number\"" << endl;
    
    unsigned fluxColumnForFirstBeam;
    if(plotContinuum) {
        *fgnu << "set cblabel \"continuum flux\"" << endl;
        fluxColumnForFirstBeam = 14;
    } else {
        *fgnu << "set cblabel \"normalized flux\"" << endl;

        *fgnu << "set zrange[0:1]" << endl;
        fluxColumnForFirstBeam = 13;
    }
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            *fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        *fgnu << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            *fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}


/*
 * Generate 2D plot for normalized spectra
 */
void GenerateNormalizationSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool display) {
    
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;

    *fgnu << "set xrange[-200:*]" << endl;
    *fgnu << "\nset xlabel \"distance (pixels)\"" << endl;

    *fgnu << "set ylabel \"order number - (1.0 - norm flux)\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "plot \"" << datafilename << "\" u 6:($2+$10-1.0) w l lt 3" << endl;
       
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << endl;

        *fgnu << "plot \"" << datafilename << "\" u 6:($2+$10-1.0) w l lt 3" << endl;

        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }    
}

