/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaGenerateLEFormats
 Version: 1.0
 Description:  Module to generate LE compatible formats.
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
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaGenerateLEFormats.cpp */

using namespace std;

/*! 
 * operaGenerateLEFormats
 * \author Eder Martioli
 * \brief  Module to generate LE compatible formats
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
    
	string inputOperaSpectrum; 
	string outputLEfilename;
    string LEorderwavelength;
    string object = "Nowhere";
    unsigned LibreEspritSpectrumType_val = LibreEspritsp2Spectrum;
    unsigned fluxType_val = RawFluxInElectronsPerElement;
    unsigned wavelengthType_val = ThArCalibratedInNM;
    bool removePolarContinuum;
	int ordernumber = NOTPROVIDED;
    int minorder = 0;
    int maxorder = 0;
	
	args.AddRequiredArgument("inputOperaSpectrum", inputOperaSpectrum, "Extended opera spectrum (.spc)");
	args.AddRequiredArgument("outputLEfilename", outputLEfilename, "Libre-Esprit spectrum (.s)");
	args.AddOptionalArgument("LEorderwavelength", LEorderwavelength, "", "Table with LE order wavelength ranges");
	args.AddRequiredArgument("object", object, "Object name");
	args.AddRequiredArgument("LibreEspritSpectrumType", LibreEspritSpectrumType_val, "Spectrum type");
	args.AddRequiredArgument("fluxType", fluxType_val, "Flux type");
	args.AddRequiredArgument("wavelengthType", wavelengthType_val, "Wavelength type");
	args.AddOptionalArgument("removePolarContinuum", removePolarContinuum, false, "Use continuum polarization removal");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	
	try {
		args.Parse(argc, argv);
		
		operaSpectralOrder_t LibreEspritSpectrumType = operaSpectralOrder_t(LibreEspritSpectrumType_val);
		/*  Available LibreEspritSpectrumType options for LE formats are:
			LibreEspritpolarimetry
			LibreEspritpolSpectrum
			LibreEspritsp1Spectrum
			LibreEspritsp2Spectrum
		*/
		operaFluxType_t fluxType = operaFluxType_t(fluxType_val);
		/*  Available operaFluxType_t options are:
			1 = RawFluxInElectronsPerElement
			2 = NormalizedFluxToContinuum
			3 = CalibratedFluxNormalizedToRefWavelength
		*/
		operaWavelengthType_t wavelengthType = operaWavelengthType_t(wavelengthType_val);
		/*  Available operaWavelengthType_t options are:
			1 = ThArCalibratedInNM
			2 = TelluricCorrectedWavelengthInNM
			3 = RVCorrectedWavelengthInNM
			4 = RVAndTelluricCorrectedWavelengthInNM
		*/
		
		// we need an input .e spectrum...
		if (inputOperaSpectrum.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputLEfilename.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (args.verbose) {
			cout << "operaGenerateLEFormats: inputOperaSpectrum = " << inputOperaSpectrum << endl; 
            cout << "operaGenerateLEFormats: outputLEfilename = " << outputLEfilename << endl;
            cout << "operaGenerateLEFormats: LEorderwavelength = " << LEorderwavelength << endl;
			cout << "operaGenerateLEFormats: object = " << object << endl;
			cout << "operaGenerateLEFormats: LibreEspritSpectrumType = " << LibreEspritSpectrumType << endl;
			cout << "operaGenerateLEFormats: fluxType = " << fluxType << endl;
			cout << "operaGenerateLEFormats: wavelengthType = " << wavelengthType << endl;
		}
        
		/*
		 * Down to business, read in all the source and calibration data.
		 */        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputOperaSpectrum);
		
		UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaGenerateLEFormats: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        // Trim orders
        if(!LEorderwavelength.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, LEorderwavelength);
			spectralOrders.trimOrdersByWavelengthRanges(minorder, maxorder);
		}

		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
			if (spectralOrder->gethasSpectralElements()) {
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
				
                if (spectralElements->getHasExtendedBeamFlux()){
                    switch (fluxType) {
                        case RawFluxInElectronsPerElement:
							spectralOrder->CopyRawFluxIntoFluxVector();
                            break;
                        case NormalizedFluxToContinuum:
							spectralOrder->CopyNormalizedFluxIntoFluxVector();
                            break;
                        case CalibratedFluxNormalizedToRefWavelength:
							spectralOrder->CopyFcalFluxIntoFluxVector();
                            break;
                        default:
                            break;
                    }
                    switch (wavelengthType) {
                        case ThArCalibratedInNM:
                            break;
                        case TelluricCorrectedWavelengthInNM:
                            spectralElements->copyFROMtell();
                            break;
                        case RVCorrectedWavelengthInNM:
                            spectralOrder->applyWavelengthCorrectionFromExtendedRvel();
                            break;
                        case RVAndTelluricCorrectedWavelengthInNM:
                            spectralElements->copyFROMtell();
                            spectralOrder->applyWavelengthCorrectionFromExtendedRvel();
                            break;
                        default:
                            break;
                    }
                }
			}
			if (spectralOrder->gethasPolarimetry() && removePolarContinuum) {
				operaPolarimetry *polarimetry = spectralOrder->getPolarimetry();
				if (polarimetry->getHasContinuumRemoved()) {
					polarimetry->copyFROMcontinuumremoved();
				}
			}
		}        
 		// output wavelength/flux calibrated spectrum...
		spectralOrders.setObject(object);
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputLEfilename, LibreEspritSpectrumType);
	}
	catch (operaException e) {
		cerr << "operaGenerateLEFormats: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaGenerateLEFormats: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
