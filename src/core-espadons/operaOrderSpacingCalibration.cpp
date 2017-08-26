/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaOderSpacingCalibration
 Version: 1.0
 Description: this module performs the measurement of spacing between orders
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

#include <fstream>
#include "libraries/operaIOFormats.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaOrderSpacingCalibration.cpp */

using namespace std;

/*! 
 * operaOrderSpacingCalibration
 * \author Doug Teeple
 * \brief Order spacing measurements.
 * \arg argc
 * \arg argv
 * \note --orderspacingoutput=...
 * \note --masterbias=...
 * \note --masterflat=...
 * \note --badpixelmask=...
 * \note --inputGainFile=...
 * \note --subformat=...
 * \note --detectionMethod=...
 * \note --FFTfilter=...
 * \note --aperture=...
 * \note --numberOfsamples=...
 * \note --plotfilename=...
 * \note --datafilename=...
 * \note --scriptfilename=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

void GenerateOrderSpacingPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display);

int main(int argc, char *argv[])
{
	operaArgumentHandler args;

	string orderspacingoutput;
	string masterbias; 
	string masterflat; 
	string badpixelmask;
	string inputGainFile;
	double gain = 1;
    double noise = 5;
	string subformat_str;
	int detectionMethod = 1; // 1. Gaussian, 2. IP, 3. Top-hat
	bool FFTfilter = false;
	int aperture = 20;
    unsigned nsamples = 30;
    unsigned sampleCenterPosition = 1;  // This position is with respect to the dispersion direction (rows for Espadons)
    unsigned referenceOrderNumber = 55;          // Number of reference order for order number identification
    double referenceOrderSeparation =  67.0;         // Order separation in pixels for order number identification
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
    bool interactive = false;
    
    args.AddRequiredArgument("orderspacingoutput", orderspacingoutput, "Order spacing output file name");
    args.AddOptionalArgument("masterbias", masterbias, "", "FITS image with masterbias");
    args.AddRequiredArgument("masterflat", masterflat, "FITS image with masterflat");
    args.AddRequiredArgument("badpixelmask", badpixelmask, "FITS image with badpixel mask");
    args.AddOptionalArgument("inputGainFile", inputGainFile, "", "Input noise/gain/bias file");
    args.AddOptionalArgument("defaultgain", gain, 1, "Gain value if inputGainFile not provided");
    args.AddOptionalArgument("defaultnoise", noise, 5, "Noise value if inputGainFile not provided");
    args.AddRequiredArgument("subformat", subformat_str, "Image subformat to be inspected \"x1 x2 y1 y2\"");
    args.AddOptionalArgument("detectionMethod", detectionMethod, 1, "Method for detecting orders: 1 = Gaussian, 2 = IP, 3 = Top-hat");
    args.AddSwitch("FFTfilter", FFTfilter, "Activate Fourier smoothing filter");
    args.AddRequiredArgument("aperture", aperture, "Aperture size in pixel units");
    args.AddRequiredArgument("numberOfsamples", nsamples, "Number of row samples for detecting orders");
    args.AddRequiredArgument("sampleCenterPosition", sampleCenterPosition, "Detector position to center samples, along the dispersion direction (rows for Espadons)");
    args.AddRequiredArgument("referenceOrderNumber", referenceOrderNumber, "Number of reference order for order number identification");
    args.AddRequiredArgument("referenceOrderSeparation", referenceOrderSeparation, "Order separation in pixels for order number identification");
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
	
	try {
		args.Parse(argc, argv);
		
		struct subformat {
			unsigned x1, x2;
			unsigned y1, y2;
		} subformat = {8, 2040, 3, 4600};
		SplitStringIntoVals(subformat_str, subformat.x1, subformat.x2, subformat.y1, subformat.y2);

		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaOrderSpacingCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need some place to put data...
		if (orderspacingoutput.empty()) {
			throw operaException("operaOrderSpacingCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (args.verbose) {
			cout << "operaOrderSpacingCalibration: orderspacingoutput = " << orderspacingoutput << endl;
			cout << "operaOrderSpacingCalibration: masterbias = " << masterbias << endl;            
			cout << "operaOrderSpacingCalibration: masterflat = " << masterflat << endl; 	
			cout << "operaOrderSpacingCalibration: badpixelmask = " << badpixelmask << endl;
			cout << "operaOrderSpacingCalibration: inputGainFile = " << inputGainFile << endl;            
			cout << "operaOrderSpacingCalibration: subformat = " << subformat.x1 << " " << subformat.x2 << " "<< subformat.y1  << " "<< subformat.y2 << "\n";
			cout << "operaOrderSpacingCalibration: detectionMethod = " << detectionMethod << endl;
			cout << "operaOrderSpacingCalibration: FFTfilter = " << FFTfilter << endl;
			cout << "operaOrderSpacingCalibration: aperture = " << aperture << endl;
			cout << "operaOrderSpacingCalibration: nsamples = " << nsamples << endl;
			cout << "operaOrderSpacingCalibration: referenceOrderNumber = " << referenceOrderNumber << endl;
			cout << "operaOrderSpacingCalibration: referenceOrderSeparation = " << referenceOrderSeparation << endl;
            if(args.plot) {
                cout << "operaOrderSpacingCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaOrderSpacingCalibration: datafilename = " << datafilename << endl;
                cout << "operaOrderSpacingCalibration: scriptfilename = " << scriptfilename << endl;
                cout << "operaOrderSpacingCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
            }            
		}
        
        ofstream fdata;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        
		// Open input images and load data into an in-memory image
		unsigned x1 = subformat.x1;
		unsigned y1 = subformat.y1;
		unsigned nx = subformat.x2 - subformat.x1;
		unsigned ny = subformat.y2 - subformat.y1;
		if (args.verbose) {
			cout << "operaOrderSpacingCalibration: x1,y1,nx,ny = " << x1 << ' ' << y1 << ' ' << nx  << ' ' << ny << '\n'; 
		}
                
        operaFITSImage flat(masterflat, tfloat, READONLY);

		if (!masterbias.empty()){
			operaFITSImage bias(masterbias, tfloat, READONLY);
			flat -= bias; // remove bias from masterflat
		}
        
        operaFITSImage *badpix = NULL;
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *badpix = 1.0;
        }
		
        float slit = (float)aperture;
        
		operaSpectralOrderVector spectralOrders(MAXORDERS, ny, ny, 0);        
		
		// Read gain and noise from input file
        if (!inputGainFile.empty()) {
			unsigned amp = 0;
            operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputGainFile);
            gain = spectralOrders.getGainBiasNoise()->getGain(amp);
            noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		}
		if (args.verbose) cout << "operaOrderSpacingCalibration: gain="<< gain << " noise=" << noise << endl;
        
        // The values below are used to apply a sigma clipping to clean outliers.
        unsigned cleanbinsize = 7;  // number of points in the bin
        float nsigcut = 2.5;        // clipping region in units of sigma
        
        spectralOrders.fitOrderSpacingPolynomial(flat, *badpix, slit, nsamples, sampleCenterPosition, referenceOrderNumber, (float)referenceOrderSeparation, detectionMethod, FFTfilter, (float)gain, (float)noise, subformat.x1, subformat.x2, subformat.y1, subformat.y2, cleanbinsize, nsigcut, &fdata);
        // Note: spacing polynomial must be changed to order number versus order separation
        // then one can project all separations as a map and identify orders later in
        // in geometry. The way it works now it is relying on a single point to identify
        // all orders and this has shown to be unreliable.
        
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateOrderSpacingPlot(scriptfilename,plotfilename,datafilename, interactive);
            }
        }
        
		// flush out the order spacing output
		if (!orderspacingoutput.empty()) operaIOFormats::WriteFromSpectralOrders(spectralOrders, orderspacingoutput, Orderspacing);
		
        if(badpix) delete badpix;
		flat.operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "operaOrderSpacingCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaOrderSpacingCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaOrderSpacingCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

void GenerateOrderSpacingPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName.c_str(),"w");
    
    fprintf(fgnu,"reset\n");
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"\nset xlabel \"order number \"\n");
    fprintf(fgnu,"set ylabel \"order separation (pixels)\"\n");
    
    if(!outputPlotEPSFileName.empty()) {
        fprintf(fgnu,"\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
        fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName.c_str());
        
        fprintf(fgnu,"\nplot \"%s\" u 1:2:3 with yerr,\"%s\" u 1:4 w l\n",dataFileName.c_str(),dataFileName.c_str());
        
        if (display) {
            fprintf(fgnu,"\nset terminal x11\n");
            fprintf(fgnu,"set output\n");
            fprintf(fgnu,"replot\n");
        } else {
            fprintf(fgnu,"\n#set terminal x11\n");
            fprintf(fgnu,"#set output\n");
            fprintf(fgnu,"#replot\n"); 
        }
    } else {
        fprintf(fgnu,"\nplot \"%s\" u 1:2:3 with yerr,\"%s\" u 1:4 w l\n",dataFileName.c_str(),dataFileName.c_str());
        
        fprintf(fgnu,"\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
        fprintf(fgnu,"#set output \"outputPlotEPSFileName.eps\"\n");
        fprintf(fgnu,"#replot\n");
        fprintf(fgnu,"#set terminal x11\n");
        fprintf(fgnu,"#set output\n");
    }
    
    fclose(fgnu);
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}
