/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaSNR
 Version: 1.0
 Description: Calculate SNR stats.
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

/*! \file operaSNR.cpp */

using namespace std;

/*! 
 * operaSNR
 * \author Doug Teeple
 * \brief Signal to Noise calculation.
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

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string inputfilename;
	string outputfilename;
	string wcalfilename; 
	string object;
	unsigned spectralOrderType_val = SNR;
	bool centralsnr = false;
	double spectralbinsize;
	
	args.AddRequiredArgument("input", inputfilename, "Input spectrum file");
	args.AddRequiredArgument("output", outputfilename, "Output SNR file");
	args.AddRequiredArgument("wavelengthCalibration", wcalfilename, "Wavelength calibration file");
	args.AddRequiredArgument("object", object, "Object name, needed for Libre-Esprit output");
	args.AddRequiredArgument("spectrumtype", spectralOrderType_val, "Spectrum type");
	args.AddRequiredArgument("centralsnr", centralsnr, "Use central SNR");
	args.AddRequiredArgument("spectralbinsize", spectralbinsize, "Spectral bin size in pixels");
	
	try {
		args.Parse(argc, argv);
		
		operaSpectralOrder_t spectralOrderType = operaSpectralOrder_t(spectralOrderType_val);
		
		if (inputfilename.empty()) {
			throw operaException("operaSNR: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (outputfilename.empty()) {
			throw operaException("operaSNR: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (wcalfilename.empty()) {
			throw operaException("operaSNR: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		if (args.verbose) {
			cout << "operaSNR: input = " << inputfilename << endl;
			cout << "operaSNR: object = " << object << endl;
			cout << "operaSNR: output = " << outputfilename << endl;
			cout << "operaSNR: spectrum type = " << spectralOrderType << endl;
			cout << "operaSNR: centralsnr = " << centralsnr << endl;
			cout << "operaSNR: wavelength calibration file = " << wcalfilename << endl;
			cout << "operaSNR: spectralbinsize = " << spectralbinsize << endl;
		}
		
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputfilename);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, wcalfilename);
		spectralOrders.setObject(object);
		
		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();
		for (unsigned order=minorder; order <= maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
				if (spectralElements->getnSpectralElements() > 0) {
					spectralOrder->setWavelengthsFromCalibration();
					spectralOrder->sethasCenterSNROnly(centralsnr);
					spectralOrder->calculateSNR();	// has side effect of retaining center SNR
					spectralOrder->setsnrSpectralBinSize(spectralbinsize);
				}
			}
		}
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputfilename, spectralOrderType);
	}
	catch (operaException e) {
		cerr << "operaSNR: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaSNR: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
