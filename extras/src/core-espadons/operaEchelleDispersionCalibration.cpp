/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaEchelleDispersionCalibration
 Version: 1.0
 Description: Calculate the Echelle Dispersion Polynomial.
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
#include <iostream>
#include <string>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaWavelength.h" // for MAXORDEROFWAVELENGTHPOLYNOMIAL
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "core-espadons/operaEchelleDispersionCalibration.h"

/*! \file operaEchelleDispersionCalibration.cpp */

using namespace std;

#define NOTPROVIDED -999

/*! 
 * operaEchelleDispersionCalibration
 * \author Eder Martioli
 * \brief Calculate the Echelle Dispersion Polynomial.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
	string inputWaveFile;
	string outputDispFile;
    string outputWaveFile;
    
    int minorderOfLaurentPolynomial = -3;
    int maxorderOfLaurentPolynomial = 0;
	
    unsigned binsizeToRemoveOutliers = 7;
    float thresholdToRemoveOutliers = 2;
    
    unsigned maxorderofpolynomial = 4;      // maximum degree of polynomial for wavelength solution

    bool updateOrderOfPolynomial = false;
    
    unsigned minoutputorder = 22;
    bool minorderprovided = false;
    unsigned maxoutputorder = 61;
    bool maxorderprovided = false;
    int ordernumber = NOTPROVIDED;	    
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0;
    
    int plot=0;
    string plotfilename;	
	string datafilename;	
	string scriptfilename;	

	struct option longopts[] = {
		{"inputWaveFile",1, NULL, 'w'},
		{"outputWaveFile",1, NULL, 'u'},
		{"outputDispFile",1, NULL, 'o'},
		{"minorderOfLaurentPolynomial",1, NULL, 'M'},
		{"maxorderOfLaurentPolynomial",1, NULL, 'L'},
		{"binsizeToRemoveOutliers",1, NULL, 'B'},
		{"thresholdToRemoveOutliers",1, NULL, 'T'},
		{"maxorderofpolynomial",1, NULL, 'Y'},          // maximum degree of polynomial for wavelength solution
        {"ordernumber",1, NULL, 'O'},
		{"minoutputorder",1, NULL, 'N'},
		{"maxoutputorder",1, NULL, 'X'},
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},  
		{"interactive",0, NULL, 'I'},
        {"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "w:u:o:M:L:B:T:Y:O:N:X:P:F:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'w':
				inputWaveFile = optarg;
				break;
			case 'u':
				outputWaveFile = optarg;
				break;
			case 'o':		// output
				outputDispFile = optarg;
				break;
			case 'M':
				minorderOfLaurentPolynomial = atoi(optarg);
				break;
			case 'L':
				maxorderOfLaurentPolynomial = atoi(optarg);
				break;
            case 'B':
				binsizeToRemoveOutliers = atoi(optarg);
				break;
			case 'T':
				thresholdToRemoveOutliers = atof(optarg);
				break;
			case 'Y':
				maxorderofpolynomial = atoi(optarg);
				break;
			case 'O':
				ordernumber = atoi(optarg);
				break;
			case 'N':
				minoutputorder = atoi(optarg);
                minorderprovided = true;
				break;
			case 'X':
				maxoutputorder = atoi(optarg);
                maxorderprovided = true;
				break;
			case 'P':
				plotfilename = optarg;
				break; 		                
			case 'F':
				datafilename = optarg;
				break; 	
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
				interactive = true;
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
		// we need an input...
		if (inputWaveFile.empty()) {
			throw operaException("operaEchelleDispersionCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputDispFile.empty() && outputWaveFile.empty()) {
			throw operaException("operaEchelleDispersionCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}

		if (verbose) {
			cout << "operaEchelleDispersionCalibration: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaEchelleDispersionCalibration: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaEchelleDispersionCalibration: outputDispFile = " << outputDispFile << endl;
			cout << "operaEchelleDispersionCalibration: minorderOfLaurentPolynomial = " << minorderOfLaurentPolynomial << endl;
			cout << "operaEchelleDispersionCalibration: maxorderOfLaurentPolynomial = " << maxorderOfLaurentPolynomial << endl;
			cout << "operaEchelleDispersionCalibration: binsizeToRemoveOutliers = " << binsizeToRemoveOutliers << endl;
			cout << "operaEchelleDispersionCalibration: thresholdToRemoveOutliers = " << thresholdToRemoveOutliers << endl;
			cout << "operaEchelleDispersionCalibration: maxorderofpolynomial = " << maxorderofpolynomial << endl;
            if(ordernumber != NOTPROVIDED) {
                cerr << "operaWavelengthCalibration: ordernumber = " << ordernumber << endl;
            }
            if(plot) {
                cout << "operaEchelleDispersionCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaEchelleDispersionCalibration: datafilename = " << datafilename << endl;
                cout << "operaEchelleDispersionCalibration: scriptfilename = " << scriptfilename << endl;                
            }
		}
        
        ofstream *fdata = NULL;
        ofstream *forigdata = NULL;
        
		operaSpectralOrderVector spectralOrdervector(inputWaveFile);
		//unsigned NumberofBeams = 0;
		
		unsigned minorder = spectralOrdervector.getMinorder();
		unsigned maxorder = spectralOrdervector.getMaxorder();

        if(maxorder <= minorder) {
			throw operaException("operaEchelleDispersionCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        LaurentPolynomial *dispersionPolynomial[MAXORDEROFWAVELENGTHPOLYNOMIAL];
        unsigned numberOfDataPoints = (unsigned)(maxorder - minorder + 1);
        string coeffOriginalDatafilenames[MAXORDEROFWAVELENGTHPOLYNOMIAL]; // for plot
        string coeffDatafilenames[MAXORDEROFWAVELENGTHPOLYNOMIAL]; // for plot

        for(unsigned wlcoeffIndex=0; wlcoeffIndex < maxorderofpolynomial; wlcoeffIndex++) {

            dispersionPolynomial[wlcoeffIndex] = new LaurentPolynomial(minorderOfLaurentPolynomial,maxorderOfLaurentPolynomial);

            dispersionPolynomial[wlcoeffIndex]->createDataVectors(numberOfDataPoints);
        
            unsigned index = 0;
            double dispersion = 0;
            double dispersionError = 0;
            
            if (!datafilename.empty()) {
                fdata = new ofstream();
                forigdata = new ofstream();
                
				string basefilename = datafilename;
				string directory =  "./";
				if (datafilename.find_last_of("/") != string::npos) {
					basefilename = datafilename.substr(datafilename.find_last_of("/")+1);
					directory = datafilename.substr(0, datafilename.find_last_of("/"));
				}
                coeffDatafilenames[wlcoeffIndex] = directory + "/" + itos(wlcoeffIndex) + "_" + basefilename;
                coeffOriginalDatafilenames[wlcoeffIndex] = directory + "/" + itos(wlcoeffIndex) + "_origData_" + basefilename;

                forigdata->open(coeffOriginalDatafilenames[wlcoeffIndex].c_str());
                fdata->open(coeffDatafilenames[wlcoeffIndex].c_str());

            }
            
            for (unsigned order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrdervector.GetSpectralOrder(order);
				if (spectralOrder->gethasWavelength()) {
					Polynomial *wavelengthPolynomial =  spectralOrder->getWavelength()->getWavelengthPolynomial();

					double *par = (double *)wavelengthPolynomial->getVector();
					double *errpar = (double *)wavelengthPolynomial->getErrorVector();
					unsigned npar = wavelengthPolynomial->getOrderOfPolynomial();
					
					if(wlcoeffIndex < npar){
						dispersion = par[wlcoeffIndex];
						dispersionError = (errpar[wlcoeffIndex]*errpar[wlcoeffIndex]);
						
						if(forigdata != NULL) {
								*forigdata << wlcoeffIndex << " " << order << " " << dispersion << " " << dispersionError << endl;
						}
						dispersionPolynomial[wlcoeffIndex]->setDataValues(index,(double)order, dispersion, sqrt(dispersionError));
						index++;
					}
                }
            }
            
            dispersionPolynomial[wlcoeffIndex]->setnDataPoints(index);
            
            for(unsigned i=0;i<dispersionPolynomial[wlcoeffIndex]->getNumberOfCoefficients();i++) {
                dispersionPolynomial[wlcoeffIndex]->setCoefficient(i,dispersion);
            }
            
            dispersionPolynomial[wlcoeffIndex]->removeOutLiersFromDataSet(binsizeToRemoveOutliers,thresholdToRemoveOutliers);
            
            if(fdata != NULL) {
                for (unsigned i=0; i<dispersionPolynomial[wlcoeffIndex]->getnDataPoints(); i++) {
                    *fdata << wlcoeffIndex << " "
                    << dispersionPolynomial[wlcoeffIndex]->getXdataValue(i) << " "
                    << dispersionPolynomial[wlcoeffIndex]->getYdataValue(i) << " "
                    << dispersionPolynomial[wlcoeffIndex]->getYerrorValue(i) << endl;
                }
            }
            
            dispersionPolynomial[wlcoeffIndex]->FitModeltoData();
            
            if(verbose){
                cout << "operaEchelleDispersionCalibration: dispersion solution => ";
                dispersionPolynomial[wlcoeffIndex]->printEquation(&cout);
            }
            
            double rms = dispersionPolynomial[wlcoeffIndex]->calculateRMSofResiduals();
            if(verbose)
                cout << "operaEchelleDispersionCalibration: RMS=" <<  rms << endl;
            
            PolynomialCoeffs_t *pc = dispersionPolynomial[wlcoeffIndex]->getLaurentPolynomialCoeffs();
            spectralOrdervector.setDispersionPolynomial(wlcoeffIndex,minorderOfLaurentPolynomial,maxorderOfLaurentPolynomial,pc);
            if(fdata != NULL) {
               fdata->close();
               forigdata->close();
            }
        }
        
        if(!datafilename.empty()) {
            if (!scriptfilename.empty()) {
                GenerateEchelleDispersionPlot(scriptfilename,plotfilename,coeffDatafilenames,coeffOriginalDatafilenames,maxorderofpolynomial, dispersionPolynomial, interactive);
            }
        }
        
        if(!minorderprovided) {
            minoutputorder = spectralOrdervector.getMinorder();
        }
        if(!maxorderprovided) {
            maxoutputorder = spectralOrdervector.getMaxorder();
        }
		
		if(ordernumber != NOTPROVIDED) {
			minoutputorder = ordernumber;
			maxoutputorder = ordernumber;
		}
		
		if (verbose) {
			cerr << "operaWavelengthCalibration: minorder = " << minoutputorder << " maxorder = " << maxoutputorder << endl;
		}
        
        /*
         * Below it uses the echelle dispersion calibration to update the wavelength solution
         * Note that this will only be called if there is an associated output calibration file
         */
        unsigned neworderofpolynomial = maxorderofpolynomial;
        
        if (!outputWaveFile.empty()) {
			spectralOrdervector.setMinorder(minoutputorder);
			spectralOrdervector.setMaxorder(maxoutputorder);
            for (unsigned order=minoutputorder; order<=maxoutputorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrdervector.GetSpectralOrder(order);
				if (!spectralOrder->gethasWavelength()) {
					spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
				}
				spectralOrder->sethasWavelength(true);
                operaWavelength *wavelength =  spectralOrder->getWavelength();
                Polynomial *wavelengthPolynomial =  wavelength->getWavelengthPolynomial();
                
                if(updateOrderOfPolynomial) {
                    wavelengthPolynomial->setOrderOfPolynomial(neworderofpolynomial); // set new order of output polynomial
                } else {
                    neworderofpolynomial = wavelengthPolynomial->getOrderOfPolynomial();  // don't update order of polynomial
                }
                
                for(unsigned wlcoeffIndex=0; wlcoeffIndex < maxorderofpolynomial; wlcoeffIndex++) {
                    double wavelengthPolynomialCoeff = dispersionPolynomial[wlcoeffIndex]->Evaluate(double(order));
                    double wavelengthPolynomialError = dispersionPolynomial[wlcoeffIndex]->calculateRMSofResiduals();
                    
                    wavelengthPolynomial->setCoefficient(wlcoeffIndex, wavelengthPolynomialCoeff);
                    wavelengthPolynomial->setCoefficientError(wlcoeffIndex, wavelengthPolynomialError);
                }
                if (verbose) {
					wavelengthPolynomial->printEquation(&cout);
				}
            }
        }
        
        for(unsigned wlcoeffIndex=0; wlcoeffIndex < maxorderofpolynomial; wlcoeffIndex++) {
            delete dispersionPolynomial[wlcoeffIndex];
        }
        spectralOrdervector.setnumberOfDispersionPolynomials(maxorderofpolynomial);
        
        if (!outputDispFile.empty()) {
			spectralOrdervector.WriteSpectralOrders(outputDispFile, Disp);
		}
        if (!outputWaveFile.empty()) {
            spectralOrdervector.WriteSpectralOrders(outputWaveFile, Wave);
		}
        
	}
	catch (operaException e) {
		cerr << "operaEchelleDispersionCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaEchelleDispersionCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdthp]" +
	" --inputWaveFile=<WAVE_FILE>"
	" --outputWaveFile=<WAVE_FILE>"
	" --outputDispFile=<DISP_FILE>"
	" --minorderOfLaurentPolynomial=<INT_VALUE>"
	" --maxorderOfLaurentPolynomial=<INT_VALUE>"
	" --binsizeToRemoveOutliers=<UNS_VALUE>"
	" --thresholdToRemoveOutliers=<FLT_VALUE>"
	" --maxorderofpolynomial=<UNS_VALUE>"
    " --ordernumber=<INT_VALUE"
	" --minoutputorder=<INT_VALUE>"
	" --maxoutputorder=<INT_VALUE>"
    " --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputWaveFile=/Users/edermartioli/opera/calibrations/PolarData/OLAPAa_pol_Normal.wcal.gz --minorderOfLaurentPolynomial=-3 --maxorderOfLaurentPolynomial=0 --outputDispFile=/Users/edermartioli/opera/calibrations/PolarData/OLAPAa_pol_Normal.disp.gz --datafilename=disp.dat --scriptfilename=disp.gnu\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file\n"
	"  -u, --outputWaveFile=<WAVE_FILE>, Output wavelength calibration file with updated solution for all orders \n"
	"  -o, --outputDispFile=<DISP_FILE>, Output echelle dispersion calibration file to store final solution\n"
	"  -M, --minorderOfLaurentPolynomial=<INT_VALUE>, Minimum power of Laurent Polynomial solution. WARNING: Current version only supports fitting for min power negative. \n"
	"  -L, --maxorderOfLaurentPolynomial=<INT_VALUE>, Maximum power of Laurent Polynomial solution. WARNING: Current version only supports fitting for max power = 0. \n"
	"  -B, --binsizeToRemoveOutliers=<UNS_FILE>, Bin size used to estimate dispersion and remove outliers\n"
	"  -T, --thresholdToRemoveOutliers=<FLT_FILE>, Threshold above which outliers are excluded in units of absolute deviation\n"
    "  -Y, --maxorderofpolynomial=<UNS_VALUE>, Maximum degree of polynomial for input wavelength solution\n"    
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minoutputorder=<INT_VALUE>, Define minimum output order number\n"
	"  -X, --maxoutputorder=<INT_VALUE>, Define maximum output order number\n";
}

void GenerateEchelleDispersionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileNames[], string origdataFileNames[], unsigned npolynomials, LaurentPolynomial *polynomials[], bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key\n" << endl;
    *fgnu << "set ylabel \"y (nm, nm/pixel, nm/pixel^2,..)\"" << endl;
    *fgnu << "set xlabel \"order number\"" << endl;
    
    *fgnu << "set pointsize 1.0" << endl;
    
    for(unsigned k=0;k<npolynomials;k++) {
        *fgnu << "poly" << k;
        polynomials[k]->printEquation(fgnu);
    }
    
    *fgnu << "NX=2; NY=2" << endl;
    *fgnu << "DX=0.1; DY=0.1; SX=0.4; SY=0.4" << endl;
    *fgnu << "set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
    *fgnu << "set size SX*NX+DX*2,SY*NY+DY*2" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "set multiplot" << endl;
        *fgnu << "set size 0.9*SX,0.9*SY" << endl;
        
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2tics; set ylabel \"{/Symbol l}(x=0) (nm)\"" << endl;
        *fgnu << "set origin DX,DY+SY" << endl;
        *fgnu << "plot \"" << origdataFileNames[0] << "\" u 2:3 w p pt 6, \"" << dataFileNames[0] << "\" u 2:3 w p pt 7, poly0f(x)" << endl;
        *fgnu << endl;

        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"d{/Symbol l}/dx (x=0) (nm/pixel)\"" << endl;
        *fgnu << "set origin DX+SX,DY+SY" << endl;
        *fgnu << "plot \"" << origdataFileNames[1] << "\" u 2:3 w p pt 6, \"" << dataFileNames[1] << "\" u 2:3 w p pt 7, poly1f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2label; set ylabel \"(1/2)*d^2{/Symbol l}/dx^2 (x=0) (nm/pixel^2)\"" << endl;
        *fgnu << "set origin DX,DY" << endl;
        *fgnu << "plot \"" << origdataFileNames[2] << "\" u 2:3 w p pt 6, \"" << dataFileNames[2] << "\" u 2:3 w p pt 7, poly2f(x)" << endl;
        *fgnu << endl;

        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"(1/6)*d^3{/Symbol l}/dx^3 (x=0)  (nm/pixel^3)\"" << endl;
        *fgnu << "set origin DX+SX,DY" << endl;
        *fgnu << "plot \"" << origdataFileNames[3] << "\" u 2:3 w p pt 6, \"" << dataFileNames[3] << "\" u 2:3 w p pt 7, poly3f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset multiplot" << endl;
        
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
        *fgnu << "set multiplot" << endl;
        *fgnu << "set size 0.9*SX,0.9*SY" << endl;
        
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2tics; set ylabel \"{/Symbol l}(x=0) (nm)\"" << endl;
        *fgnu << "set origin DX,DY+SY" << endl;
        *fgnu << "plot \"" << origdataFileNames[0] << "\" u 2:3 w p pt 6, \"" << dataFileNames[0] << "\" u 2:3 w p pt 7, poly0f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"d{/Symbol l}/dx (x=0) (nm/pixel)\"" << endl;
        *fgnu << "set origin DX+SX,DY+SY" << endl;
        *fgnu << "plot \"" << origdataFileNames[1] << "\" u 2:3 w p pt 6, \"" << dataFileNames[1] << "\" u 2:3 w p pt 7, poly1f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2label; set ylabel \"(1/2)*d^2{/Symbol l}/dx^2 (x=0) (nm/pixel^2)\"" << endl;
        *fgnu << "set origin DX,DY" << endl;
        *fgnu << "plot \"" << origdataFileNames[2] << "\" u 2:3 w p pt 6, \"" << dataFileNames[2] << "\" u 2:3 w p pt 7, poly2f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"(1/6)*d^3{/Symbol l}/dx^3 (x=0)  (nm/pixel^3)\"" << endl;
        *fgnu << "set origin DX+SX,DY" << endl;
        *fgnu << "plot \"" << origdataFileNames[3] << "\" u 2:3 w p pt 6, \"" << dataFileNames[3] << "\" u 2:3 w p pt 7, poly3f(x)" << endl;
        *fgnu << endl;
        
        *fgnu << "unset multiplot" << endl;
        
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
