/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaBias
 Version: 1.0
 Description: Calculate the chip bias.
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

/*! \file operaBias.cpp */

using namespace std;

/*!
 * operaBias
 * \author Doug Teeple / Eder Martioli
 * \brief Calculate bias levels of each amp and update the GainNoiseBias stats file.
 * \arg argc
 * \arg argv
 * \note --bias=...
 * \note 
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOutput
 * \ingroup tools
 * \return EXIT_STATUS
 */

#define MAXIMAGES 1000

void CalculateMedianBiases(string filename, const bool dooverscan, unsigned short& medianBiasA, unsigned short& medianBiasB) {
	const unsigned overscan = 20;
	const unsigned AmpBStartCol = 1045;
	const unsigned AmpBEndCol = 2068;
	const unsigned CCDX = 2080;
	const unsigned CCDY = 4608;
	operaFITSImage bias(filename, READONLY);
	if (dooverscan) {
		// ampA
		unsigned short *pixels = new unsigned short[CCDY * overscan];
		unsigned short *p = pixels;
		for (unsigned y=0; y<CCDY; y++) {
			for (unsigned x=0; x<overscan; x++) {
				*p++ = bias.getpixelUSHORT(x, y);
			}
		}
		medianBiasA = operaArrayMedianQuickUSHORT(CCDY, pixels);	// scrambles the pixels
		// ampB
		p = pixels;
		for (unsigned y=0; y<CCDY; y++) {
			for (unsigned x=AmpBEndCol; x<CCDX; x++) {
				*p++ = bias.getpixelUSHORT(x, y);
			}
		}
		medianBiasB = operaArrayMedianQuickUSHORT(CCDY, pixels);	// scrambles the pixels
		delete[] pixels;
	} else {
		// ampA
		unsigned short *pixels = new unsigned short[(AmpBStartCol-overscan) * CCDY];
		unsigned short *p = pixels;
		for (unsigned y=0; y<CCDY; y++) {
			for (unsigned x=overscan; x<AmpBStartCol; x++) {
				*p++ = bias.getpixelUSHORT(x, y);
			}
		}
		medianBiasA = operaArrayMedianQuickUSHORT((AmpBStartCol-overscan) * CCDY, pixels);	// scrambles the pixels
		// ampB
		p = pixels;
		for (unsigned y=0; y<CCDY; y++) {
			for (unsigned x=AmpBStartCol; x<AmpBEndCol; x++) {
				*p++ = bias.getpixelUSHORT(x, y);
			}
		}
		medianBiasB = operaArrayMedianQuickUSHORT((AmpBEndCol-AmpBStartCol) * CCDY, pixels);	// scrambles the pixels
		delete[] pixels;
	}
	bias.operaFITSImageClose();
}

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	string biaslist;
	string gainNoiseBias;
	string output;
	bool dooverscan; // whether to use the median of the whole array or just the last column of ampA and the first column of ampB (reduce troughing)
	unsigned nampsexpected = 1;
	
	args.AddRequiredArgument("bias", biaslist, "list of biases, separated by spaces");
	args.AddRequiredArgument("gain", gainNoiseBias, "gain / noise / bias file");
	args.AddRequiredArgument("output", output, "output .bias file");
	args.AddSwitch("overscan", dooverscan, "do overscan only");
	args.AddRequiredArgument("numberofamplifiers", nampsexpected, "how many amps are expected");
	
	try  {
		args.Parse(argc, argv);
		vector<string> biasimgs;
		if (!biaslist.empty()) { //split biaslist on spaces into vector
			stringstream ss(biaslist);
			string temp;
			while (getline(ss, temp, ' ')) biasimgs.push_back(temp);
		}
		
		if (output.empty()) {
			throw operaException("operaBias: please specify an output ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (biasimgs.empty() || biasimgs[0].empty()) {
			throw operaException("operaBias: please specify a bias ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		operaSpectralOrderVector spectralOrders;
		unsigned namps = 0;
		if (!gainNoiseBias.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, gainNoiseBias);
			namps = spectralOrders.getGainBiasNoise()->getAmps();
		}
		if (args.verbose) {
			cout << "operaBias: gain " << gainNoiseBias << endl;
			cout << "operaBias: namps " << namps << endl;
			cout << "operaBias: nampsexpected " << nampsexpected << endl;
			cout << "operaBias: do overscan only " << dooverscan << endl;
		}
		
		for (unsigned i=0; i<biasimgs.size(); i++) {
			if (args.verbose) {
				cout << "operaBias: bias " << biasimgs[i] << endl;
			}
			// get the two median bias levels
			unsigned short medianBiasA;
			unsigned short medianBiasB;
			CalculateMedianBiases(biasimgs[i], dooverscan, medianBiasA, medianBiasB);
			
			if (isnan(medianBiasA))
				throw operaException("operaBias: bias 0 : ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(medianBiasA))
				throw operaException("operaBias: bias 0: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(medianBiasB))
				throw operaException("operaBias: bias 1: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(medianBiasB))
				throw operaException("operaBias: bias 1: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (args.verbose) {
				cout << "operaBias: Changing Bias ampA from " << spectralOrders.getGainBiasNoise()->getBias(0) << " to " << medianBiasA  << endl;					
				if (nampsexpected > 1) {
					cout << "operaBias: Changing Bias ampB from " << spectralOrders.getGainBiasNoise()->getBias(1) << " to " << medianBiasB  << endl;					
				}
			}
			// set the gain/noise/bias
			spectralOrders.getGainBiasNoise()->setBias(0, medianBiasA);
			spectralOrders.getGainBiasNoise()->setBias(1, medianBiasB);
		}
		if (namps < 2) spectralOrders.getGainBiasNoise()->setAmps(nampsexpected);
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, output, GainNoise);
	}
	catch (operaException e) {
		cerr << "operaBias: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaBias: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
