/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaGain
 Version: 1.0
 Description: Calculate gain and noise.
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
#include "libraries/operaCCD.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaGain.cpp */

using namespace std;

/*!
 * operaGain
 * \author Doug Teeple & Eder Martioli
 * \brief Output gain and noise based on at least 2 flat images.
 * \arg argc
 * \arg argv
 * \arg --biasimgs=...
 * \arg --flatimgs=...
 * \arg --badpixelmask=...
 * \arg --output=...
 * \arg --defaultgain=...
 * \arg --defaultnoise=...
 * \note --subwindow="x1 nx y1 ny"
 * \note --defaultgain=...
 * \note --defaultnoise=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 */

/**
 * \defgroup core Core Modules
 */

int main(int argc, char *argv[])
{
	const unsigned MAXIMAGES=1000;
	const unsigned MAXNAMPS=2;
	
	operaArgumentHandler args;
    
	string listofbiasimgs;
	string listofflatimgs;
	string biaslistfile;
	string flatlistfile;
	string output; 
	string badpixelmask;
	unsigned namps = 1; // 1 for EEV1 and OLAPAa, 2 for OLAPAab	
	double defaultgain = 1.0;
	double defaultnoise = 0.0;
	double minimumAllowedGain = 0.0;
	unsigned gainMinPixPerBin = 100;
	unsigned gainMaxNBins = 1000;
	unsigned gainLowestCount = 1000;
	unsigned gainHighestCount = 25000;
	string subwindow_str;
	string datasec_str;
	string dseca_str;
	string dsecb_str;
	//unsigned maximages = 12;    // otherwise we run out of memory on 32 bit systems...
	
	args.AddOptionalArgument("biasimgs", listofbiasimgs, "", "List of bias fits files, seperated by spaces");
	args.AddOptionalArgument("flatimgs", listofflatimgs, "", "List of flat fits files, seperated by spaces");
	args.AddOptionalArgument("biaslistfile", biaslistfile, "", "File containing list of bias fits files");
	args.AddOptionalArgument("flatlistfile", flatlistfile, "", "File containing list of flat fits files");
	args.AddRequiredArgument("output", output, "Output file");
	args.AddRequiredArgument("badpixelmask", badpixelmask, "Badpixel mask fits file");	
	args.AddRequiredArgument("numberofamplifiers", namps, "Number of amplifiers (1 or 2)");
	args.AddRequiredArgument("defaultgain", defaultgain, "Default gain value");
	args.AddRequiredArgument("defaultnoise", defaultnoise, "Default noise value");
	args.AddRequiredArgument("gainMinPixPerBin", gainMinPixPerBin, "Minimum number of pixels allowed per bin");
	args.AddRequiredArgument("gainMaxNBins", gainMaxNBins, "Maximum number of bins to use");
	args.AddRequiredArgument("gainLowestCount", gainLowestCount, "Lowest pixel count value to be considered");
	args.AddRequiredArgument("gainHighestCount", gainHighestCount, "Highest pixel count value to be considered");
	args.AddOptionalArgument("minimumAllowedGain", minimumAllowedGain, 0.0, "Minimum allowed gain");
	args.AddRequiredArgument("subwindow", subwindow_str, "Subwindow where data will used for the gain/noise calculation \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DATASEC", datasec_str, "", "Valid data region on the amplifier (1 amp mode) \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DSECA", dseca_str, "", "Valid data region on amplifier A (2 amp mode) \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DSECB", dsecb_str, "", "Valid data region on amplifier B (2 amp mode) \"x1 x2 y1 y2\"");
	
	try {
		args.Parse(argc, argv);

		struct subwindow {
			unsigned x0, nx, y0, ny;
		} subwindow = {0,0,0,0};
		DATASEC_t datasec = {1,2068,1,4608};
		DATASEC_t dseca = {21,1044,1,4608};
		DATASEC_t dsecb = {1045,2068,1,4608};
		
		SplitStringIntoVals(subwindow_str, subwindow.x0, subwindow.nx, subwindow.y0, subwindow.ny);
		SplitStringIntoVals(datasec_str, datasec.x1, datasec.x2, datasec.y1, datasec.y2);
		SplitStringIntoVals(dseca_str, dseca.x1, dseca.x2, dseca.y1, dseca.y2);
		SplitStringIntoVals(dsecb_str, dsecb.x1, dsecb.x2, dsecb.y1, dsecb.y2);

		string biasimgs[MAXIMAGES];
		string flatimgs[MAXIMAGES];
		unsigned biasimgIndex = 0;
		unsigned flatimgIndex = 0;
		SplitStringIntoArray(listofbiasimgs, biasimgs, biasimgIndex, MAXIMAGES); // Split list of bias images into array
		SplitStringIntoArray(listofflatimgs, flatimgs, flatimgIndex, MAXIMAGES); // Split list of flat images into array
		ReadStringsFromFileIntoArray(biaslistfile, biasimgs, biasimgIndex, MAXIMAGES); // Read list of images from file
		ReadStringsFromFileIntoArray(flatlistfile, flatimgs, flatimgIndex, MAXIMAGES); // Read list of images from file
		
		// we need at least 1 bias...
		if (biasimgIndex < 1) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need at least 2 flats to calc the gain...
		if (flatimgIndex < 2) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// only accept number of amplifiers = 1 or 2...
		if (namps != 1 && namps != 2) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		//need DATASEC for 1 amplifier mode, but no DSECA and DSECB
		if (namps == 1 && (datasec_str.empty() || !dseca_str.empty() || !dsecb_str.empty())) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		//need DSECA and DSECB for 2 amplifier mode, but no DATASEC
		if (namps == 2 && (dseca_str.empty() || dsecb_str.empty() || !datasec_str.empty())) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output...
		if (output.empty()) {
			throw operaException("operaGain: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// this input wasn't required before, but we are being stricter about arguments now
		if(badpixelmask.empty()) {
			throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (args.debug) {
			for(unsigned i=0; i<biasimgIndex; i++) cout << "operaGain: input bias: " << i << ' ' << biasimgs[i] << endl;
			for(unsigned i=0; i<flatimgIndex; i++) cout << "operaGain: input flat: " << i << ' ' << flatimgs[i] << endl;
		}
		
		if (args.verbose) {
			cout << "operaGain: subwindow = " << subwindow.x0 << " " << subwindow.nx << " " << subwindow.y0  << " " << subwindow.ny << endl; 
			cout << "operaGain: gainMinPixPerBin = " << gainMinPixPerBin << " gainMaxNBins = " << gainMaxNBins << " gainLowestCount = " << gainLowestCount << " gainHighestCount = " << gainHighestCount << " minimumAllowedGain = " << minimumAllowedGain << endl; 
		}
		
		operaSpectralOrderVector spectralOrders;

		for (unsigned amp=0; amp<namps; amp++) { // loop over all possible amplifiers
			float gain = defaultgain;
			float noise = defaultnoise;
			float gainError = 0.0;
			float bias = 0.0;
			
			if (namps == 2) {
				if (amp == 0) datasec = dseca;
				else datasec = dsecb;
			}
			if (args.verbose) cout << "operaGain: amp " << amp+1 << " datasec = " << datasec.x1 << " : " << datasec.x2 << " , " << datasec.y1  << " : " << datasec.y2 << endl;
			
			int amp_x0 = datasec.x1 + subwindow.x0;
			int amp_y0 = datasec.y1 + subwindow.y0;
			int amp_xf = amp_x0 + subwindow.nx;
			int amp_yf = amp_y0 + subwindow.ny;
			long npixels = subwindow.ny * subwindow.nx;
			
			if(amp_xf > datasec.x2 || amp_yf > datasec.y2) {
				throw "operaGain: error: subwindow exceeds area of datasec\n";
			}
			
			if (args.verbose) {
				cout << "operaGain: amp = " << amp+1 << ": xw1 = " << amp_x0 << " xw2 = " << amp_xf << " yw1 = " << amp_y0  << " yw2 = " << amp_yf << " npixels = " << npixels << endl; 
			}
			
			// open badpixelmask and load data into a vector
			operaFITSImage badpix(badpixelmask, tfloat, READONLY);
			float *badpixdata = (float *)badpix.operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
			
			// open bias images and load data into vectors
			float *biasdata[MAXIMAGES];
			for (unsigned i=0; i<biasimgIndex; i++) {	
                operaFITSImage biasIn(biasimgs[i], tfloat, READONLY);					
				biasdata[i] = (float *)biasIn.operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
			}
			
			// open flat images and load data into vectors
			float *flatdata[MAXIMAGES];
			for (unsigned i=0; i<flatimgIndex; i++) {	
				operaFITSImage flatIn(flatimgs[i], tfloat, READONLY);			
				flatdata[i] = (float *)flatIn.operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
			}
			
			operaCCDGainNoise(npixels, biasimgIndex, biasdata, flatimgIndex, flatdata, badpixdata, gainLowestCount, gainHighestCount, gainMaxNBins, gainMinPixPerBin, &gain, &gainError, &bias, &noise);
			
			free(badpixdata);
			for (unsigned i=0; i<biasimgIndex; i++) free(biasdata[i]);
			for (unsigned i=0; i<flatimgIndex; i++) free(flatdata[i]);
			
			if (isnan(gain)) throw operaException("operaGain: gain: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gain)) throw operaException("operaGain: gain: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(noise)) throw operaException("operaGain: noise: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(noise)) throw operaException("operaGain: noise: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(gainError)) throw operaException("operaGain: gainError: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gainError)) throw operaException("operaGain: gainError: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (gain < minimumAllowedGain) throw operaException("operaGain: gain under minimumAllowedGain", operaErrorNegativeValues, __FILE__, __FUNCTION__, __LINE__);	
			
			spectralOrders.getGainBiasNoise()->setGain(amp, gain);
			spectralOrders.getGainBiasNoise()->setGainError(amp, gainError);
			spectralOrders.getGainBiasNoise()->setNoise(amp, noise);
			spectralOrders.getGainBiasNoise()->setBias(amp, bias);
			spectralOrders.getGainBiasNoise()->setDatasec(amp, datasec);
			if (args.verbose) {
				cout << "operaGain: amp " << amp+1 << " gain " << gain << " gainError " << gainError << " Noise " << noise  << " Bias " << bias << endl; 
			}
		}
		spectralOrders.getGainBiasNoise()->setAmps(namps);
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, output, GainNoise);
	}
	catch (operaException e) {
		cout << "operaGain: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cout << "operaGain: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
