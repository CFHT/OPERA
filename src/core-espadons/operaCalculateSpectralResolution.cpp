/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCalculateSpectralResolution
 Version: 1.0
 Description: Calculate spectral resolution.
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

#include "libraries/operaIOFormats.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include <fstream>

/*! \file operaCalculateSpectralResolution.cpp */

using namespace std;

void GenerateSpectralResolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileNames[], string origdataFileNames[], unsigned npolynomials, LaurentPolynomial *polynomials[], bool display);

/*! 
 * operaCalculateSpectralResolution
 * \author Eder Martioli
 * \brief Calculate the spectral resolution.
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
	operaArgumentHandler args;
	
	string inputWaveFile;
	string outputWaveFile;
	string outputResolutionFile;
    int minorderOfLaurentPolynomial = -3;
    int maxorderOfLaurentPolynomial = 0;
    unsigned maxorderofpolynomial = 4; // maximum degree of polynomial for wavelength solution
    unsigned binsizeToRemoveOutliers = 7;
    double thresholdToRemoveOutliers = 2;
    int minoutputorder = NOTPROVIDED;
    int maxoutputorder = NOTPROVIDED;
    int ordernumber = NOTPROVIDED;
    string plotfilename;
	string datafilename;
	string scriptfilename;
	bool interactive = false;
	
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddOptionalArgument("outputWaveFile", outputWaveFile, "", "Output wavelength calibration file with updated solution for all orders");
	args.AddOptionalArgument("outputResolutionFile", outputResolutionFile, "", "Output echelle dispersion calibration file to store final solution");
	args.AddRequiredArgument("minorderOfLaurentPolynomial", minorderOfLaurentPolynomial, "Minimum power of Laurent Polynomial solution. WARNING: Current version only supports fitting for min power negative");
	args.AddRequiredArgument("maxorderOfLaurentPolynomial", maxorderOfLaurentPolynomial, "Maximum power of Laurent Polynomial solution. WARNING: Current version only supports fitting for max power = 0");
	args.AddRequiredArgument("maxorderofpolynomial", maxorderofpolynomial, "Maximum degree of polynomial for input wavelength solution");
	args.AddOptionalArgument("binsizeToRemoveOutliers", binsizeToRemoveOutliers, 7, "?");
	args.AddOptionalArgument("thresholdToRemoveOutliers", thresholdToRemoveOutliers, 4, "?");
	args.AddOrderLimitArguments(ordernumber, minoutputorder, maxoutputorder, NOTPROVIDED);
	args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
	
	try {
		args.Parse(argc, argv);
		// we need an input...
		if (inputWaveFile.empty()) {
			throw operaException("operaCalculateSpectralResolution: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputResolutionFile.empty() && outputWaveFile.empty()) {
			throw operaException("operaCalculateSpectralResolution: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}

		if (args.verbose) {
			cout << "operaCalculateSpectralResolution: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaCalculateSpectralResolution: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaCalculateSpectralResolution: outputResolutionFile = " << outputResolutionFile << endl;
			cout << "operaCalculateSpectralResolution: minorderOfLaurentPolynomial = " << minorderOfLaurentPolynomial << endl;
			cout << "operaCalculateSpectralResolution: maxorderOfLaurentPolynomial = " << maxorderOfLaurentPolynomial << endl;
			cout << "operaCalculateSpectralResolution: maxorderofpolynomial = " << maxorderofpolynomial << endl;
			cout << "operaCalculateSpectralResolution: binsizeToRemoveOutliers = " << binsizeToRemoveOutliers << endl;
			cout << "operaCalculateSpectralResolution: thresholdToRemoveOutliers = " << thresholdToRemoveOutliers << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaWavelengthCalibration: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaCalculateSpectralResolution: plotfilename = " << plotfilename << endl;
                cout << "operaCalculateSpectralResolution: datafilename = " << datafilename << endl;
                cout << "operaCalculateSpectralResolution: scriptfilename = " << scriptfilename << endl;                
            }
		}
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);
		
		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();

        if(maxorder <= minorder) {
			throw operaException("operaCalculateSpectralResolution: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
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
            
            ofstream fdata;
            ofstream forigdata;
            if (!datafilename.empty()) {
                string basefilename = datafilename;
				string directory =  "./";
				if (datafilename.find_last_of("/") != string::npos) {
					basefilename = datafilename.substr(datafilename.find_last_of("/")+1);
					directory = datafilename.substr(0, datafilename.find_last_of("/"));
				}
                coeffDatafilenames[wlcoeffIndex] = directory + "/" + itos(wlcoeffIndex) + "_" + basefilename;
                coeffOriginalDatafilenames[wlcoeffIndex] = directory + "/" + itos(wlcoeffIndex) + "_origData_" + basefilename;
				
                forigdata.open(coeffOriginalDatafilenames[wlcoeffIndex].c_str());
                fdata.open(coeffDatafilenames[wlcoeffIndex].c_str());
				
            }
            
            for (unsigned order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasWavelength()) {
					Polynomial *wavelengthPolynomial =  spectralOrder->getWavelength()->getWavelengthPolynomial();

					double *par = (double *)wavelengthPolynomial->getVector();
					double *errpar = (double *)wavelengthPolynomial->getErrorVector();
					unsigned npar = wavelengthPolynomial->getOrderOfPolynomial();
					
					if(wlcoeffIndex < npar){
						dispersion = par[wlcoeffIndex];
						dispersionError = (errpar[wlcoeffIndex]*errpar[wlcoeffIndex]);
						if(forigdata.is_open()) {
							forigdata << wlcoeffIndex << " " << order << " " << dispersion << " " << dispersionError << endl;
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

            if(fdata.is_open()) {
                for (unsigned i=0; i<dispersionPolynomial[wlcoeffIndex]->getnDataPoints(); i++) {
                    fdata << wlcoeffIndex << " "
                    << dispersionPolynomial[wlcoeffIndex]->getXdataValue(i) << " "
                    << dispersionPolynomial[wlcoeffIndex]->getYdataValue(i) << " "
                    << dispersionPolynomial[wlcoeffIndex]->getYerrorValue(i) << endl;
                }
            }
            
            dispersionPolynomial[wlcoeffIndex]->FitModeltoData();
            
            if(args.verbose){
                cout << "operaCalculateSpectralResolution: dispersion solution => ";
                dispersionPolynomial[wlcoeffIndex]->printEquation(&cout);
            }
            
            double rms = dispersionPolynomial[wlcoeffIndex]->calculateRMSofResiduals();
            if(args.verbose)
                cout << "operaCalculateSpectralResolution: RMS=" <<  rms << endl;
            
            PolynomialCoeffs_t *pc = dispersionPolynomial[wlcoeffIndex]->getLaurentPolynomialCoeffs();
            spectralOrders.setDispersionPolynomial(wlcoeffIndex,minorderOfLaurentPolynomial,maxorderOfLaurentPolynomial,pc);
            if(fdata.is_open()) {
				fdata.close(); 
				forigdata.close();
           }
        }
        
        if(!datafilename.empty()) {
            if (!scriptfilename.empty()) {
                GenerateSpectralResolutionPlot(scriptfilename,plotfilename,coeffDatafilenames,coeffOriginalDatafilenames,maxorderofpolynomial, dispersionPolynomial, interactive);
            }
        }
        
        UpdateOrderLimits(ordernumber, minoutputorder, maxoutputorder, spectralOrders);
		if (args.verbose) cerr << "operaWavelengthCalibration: minorder = " << minoutputorder << " maxorder = " << maxoutputorder << endl;
        
        /*
         * Below it uses the echelle dispersion calibration to update the wavelength solution
         * Note that this will only be called if there is an associated output calibration file
         */
        if (!outputWaveFile.empty()) {
            for (unsigned order=minoutputorder; order<=maxoutputorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (!spectralOrder->getWavelength()) {
					spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
				}
				spectralOrder->sethasWavelength(true);
                operaWavelength *wavelength =  spectralOrder->getWavelength();
                Polynomial *wavelengthPolynomial =  wavelength->getWavelengthPolynomial();
                
                wavelengthPolynomial->resize(maxorderofpolynomial);
                
                for(unsigned wlcoeffIndex=0; wlcoeffIndex < maxorderofpolynomial; wlcoeffIndex++) {
                    double wavelengthPolynomialCoeff = dispersionPolynomial[wlcoeffIndex]->Evaluate(double(order));
                    double wavelengthPolynomialError = dispersionPolynomial[wlcoeffIndex]->calculateRMSofResiduals();
                    
                    wavelengthPolynomial->setCoefficient(wlcoeffIndex, wavelengthPolynomialCoeff);
                    wavelengthPolynomial->setCoefficientError(wlcoeffIndex, wavelengthPolynomialError);
                }
                if (args.verbose) {
					wavelengthPolynomial->printEquation(&cout);
				}
            }
        }
        
        for(unsigned wlcoeffIndex=0; wlcoeffIndex < maxorderofpolynomial; wlcoeffIndex++) {
            delete dispersionPolynomial[wlcoeffIndex];
        }
        spectralOrders.setnumberOfDispersionPolynomials(maxorderofpolynomial);
        
        if (!outputResolutionFile.empty()) {
			operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputResolutionFile, Disp);
		}
        if (!outputWaveFile.empty()) {
            operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputWaveFile, Wave);
		}
        
	}
	catch (operaException e) {
		cerr << "operaCalculateSpectralResolution: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCalculateSpectralResolution: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

void GenerateSpectralResolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileNames[], string origdataFileNames[], unsigned npolynomials, LaurentPolynomial *polynomials[], bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key\n" << endl;
    fgnu << "set ylabel \"y (nm, nm/pixel, nm/pixel^2,..)\"" << endl;
    fgnu << "set xlabel \"order number\"" << endl;
    
    fgnu << "set pointsize 1.0" << endl;
    
    for(unsigned k=0;k<npolynomials;k++) {
        fgnu << "poly" << k;
        polynomials[k]->printEquation(&fgnu);
    }
    
    fgnu << "NX=2; NY=2" << endl;
    fgnu << "DX=0.1; DY=0.1; SX=0.4; SY=0.4" << endl;
    fgnu << "set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
    fgnu << "set size SX*NX+DX*2,SY*NY+DY*2" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "set multiplot" << endl;
        fgnu << "set size 0.9*SX,0.9*SY" << endl;
        
        fgnu << "unset xlabel" << endl;
        fgnu << "unset y2tics; set ytics" << endl;
        fgnu << "unset y2tics; set ylabel \"{/Symbol l}(x=0) (nm)\"" << endl;
        fgnu << "set origin DX,DY+SY" << endl;
        fgnu << "plot \"" << origdataFileNames[0] << "\" u 2:3 w p pt 6, \"" << dataFileNames[0] << "\" u 2:3 w p pt 7, poly0f(x)" << endl;
        fgnu << endl;
		
        fgnu << "unset xlabel" << endl;
        fgnu << "unset ytics; set y2tics" << endl;
        fgnu << "unset ylabel; set y2label \"d{/Symbol l}/dx (x=0) (nm/pixel)\"" << endl;
        fgnu << "set origin DX+SX,DY+SY" << endl;
        fgnu << "plot \"" << origdataFileNames[1] << "\" u 2:3 w p pt 6, \"" << dataFileNames[1] << "\" u 2:3 w p pt 7, poly1f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        fgnu << "unset y2tics; set ytics" << endl;
        fgnu << "unset y2label; set ylabel \"(1/2)*d^2{/Symbol l}/dx^2 (x=0) (nm/pixel^2)\"" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "plot \"" << origdataFileNames[2] << "\" u 2:3 w p pt 6, \"" << dataFileNames[2] << "\" u 2:3 w p pt 7, poly2f(x)" << endl;
        fgnu << endl;
		
        fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        fgnu << "unset ytics; set y2tics" << endl;
        fgnu << "unset ylabel; set y2label \"(1/6)*d^3{/Symbol l}/dx^3 (x=0)  (nm/pixel^3)\"" << endl;
        fgnu << "set origin DX+SX,DY" << endl;
        fgnu << "plot \"" << origdataFileNames[3] << "\" u 2:3 w p pt 6, \"" << dataFileNames[3] << "\" u 2:3 w p pt 7, poly3f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset multiplot" << endl;
        
        if (display) {
            fgnu << "\nset terminal x11" << endl;
            fgnu << "set output" << endl;
            fgnu << "replot" << endl;
        } else {
            fgnu << "\n#set terminal x11" << endl;
            fgnu << "#set output" << endl;
            fgnu << "#replot" << endl;
        }
    } else {
        fgnu << "set multiplot" << endl;
        fgnu << "set size 0.9*SX,0.9*SY" << endl;
        
        fgnu << "unset xlabel" << endl;
        fgnu << "unset y2tics; set ytics" << endl;
        fgnu << "unset y2tics; set ylabel \"{/Symbol l}(x=0) (nm)\"" << endl;
        fgnu << "set origin DX,DY+SY" << endl;
        fgnu << "plot \"" << origdataFileNames[0] << "\" u 2:3 w p pt 6, \"" << dataFileNames[0] << "\" u 2:3 w p pt 7, poly0f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset xlabel" << endl;
        fgnu << "unset ytics; set y2tics" << endl;
        fgnu << "unset ylabel; set y2label \"d{/Symbol l}/dx (x=0) (nm/pixel)\"" << endl;
        fgnu << "set origin DX+SX,DY+SY" << endl;
        fgnu << "plot \"" << origdataFileNames[1] << "\" u 2:3 w p pt 6, \"" << dataFileNames[1] << "\" u 2:3 w p pt 7, poly1f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        fgnu << "unset y2tics; set ytics" << endl;
        fgnu << "unset y2label; set ylabel \"(1/2)*d^2{/Symbol l}/dx^2 (x=0) (nm/pixel^2)\"" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "plot \"" << origdataFileNames[2] << "\" u 2:3 w p pt 6, \"" << dataFileNames[2] << "\" u 2:3 w p pt 7, poly2f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        fgnu << "unset ytics; set y2tics" << endl;
        fgnu << "unset ylabel; set y2label \"(1/6)*d^3{/Symbol l}/dx^3 (x=0)  (nm/pixel^3)\"" << endl;
        fgnu << "set origin DX+SX,DY" << endl;
        fgnu << "plot \"" << origdataFileNames[3] << "\" u 2:3 w p pt 6, \"" << dataFileNames[3] << "\" u 2:3 w p pt 7, poly3f(x)" << endl;
        fgnu << endl;
        
        fgnu << "unset multiplot" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}
