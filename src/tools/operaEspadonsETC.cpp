/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaEspadonsETC.cpp
 Version: 1.0
 Description: ESPaDOnS Exposure Time Calculator.
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
#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaLibCommon.h"

#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaIOFormats.h"
#include "libraries/GainBiasNoise.h"                // for operaGainBiasNoise

#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaSpectralTools.h"           // for calculateBlackBodyVFlux and IntegrateSpectralElementOfBlackBody

#include "libraries/operaObservingConditions.h"		// for operaObservingConditions
#include "libraries/operaObjectInTheSky.h"          // for operaObjectInTheSky
#include "libraries/operaSpectrograph.h"            // for operaSpectrograph
#include "libraries/operaTelescope.h"               // for operaTelescope

#include "libraries/operaInstrumentEnvironmentSetup.h"

#include "tools/operaEspadonsETC.h"

#define NOTPROVIDED -999

#define MAXIMUMEXPOSURETIME_S 2400

/*! \file operaEspadonsETC.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*!
 * operaEspadonsETC
 * \author Eder Martioli
 * \brief Espadons Exposure Time Calculator.
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
	int opt;

    string inputSpectrum;
    string inputWaveFile;
    string inputApertureFile;
    string inputGainFile;
    
    string outputEspadonsETCFile;
    string outputHTMLEETCFile;

    string inputSpectrographFile;
    string inputTelescopeFile;
    
    string inputObservingConditionsFile;
    string inputObjectInTheSkyFile;

    string inputAtmosphericTransmissionFile;
    
    
    string inputETCCalibrationFile;
    
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
        
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;
	string dataFilename;
	string scriptfilename;
	
    // Calculator options:
    spectrographCCD_t CCDoption = undefinedCCD; // {undefinedCCD=0, Olapa=1, EEV1=2}
    
    double ExposureTime = 300;
    bool calculateExposureTime = true;
    
    double SignalToNoise = 100;
    bool providedSNR = false;
    
    double wavelength_nm = 372;  // in nm
    bool providedWavelength = false;
    
    //    Source model
    double Vmagnitude = 6.0; // in Vega magnitudes
    double effectiveTemperature = 20000; // in deg K
    
    //   Observing conditions
    double imageQuality = 0.8;
    double airmass = 1.3;
    
    //    Instrument configuration
    EspadonsInstrumentMode_t observingMode = undefinedInstrumentMode; // {undefinedInstrumentMode=0, polarimetric=1, staronly=2, starplussky=3, GRACES_staronly=4, GRACES_starplussky=5}
    EspadonsCCDReadoutSpeed_t CCDReadoutMode = undefinedReadoutMode; // {undefinedReadoutMode=0, slowmode=1, normalmode=2, fastmode=3}
    
    //    Sky brightness
    moonphase_t moonPhase = quartermoon;   //{newmoon=0, crescentmoon=1, quartermoon=2, gibbousmoon=3, fullmoon=4}
    double angularDistFromMoon = 0.0;  // deg
    double zenithalDistOfMoon = 0.0;   // deg
    
	struct option longopts[] = {
		{"inputSpectrum",                   1, NULL, 's'},
		{"inputWaveFile",                   1, NULL, 'W'},
		{"inputApertureFile",               1, NULL, 'R'},
		{"inputGainFile",                   1, NULL, 'G'},
		{"outputEspadonsETCFile",           1, NULL, 'o'},
		{"inputETCCalibrationFile",         1, NULL, 'i'},        
		{"outputHTMLEETCFile",              1, NULL, 'H'},
		{"inputSpectrographFile",           1, NULL, 'E'},
		{"inputTelescopeFile",              1, NULL, 'T'},
		{"inputObservingConditionsFile",	1, NULL, 'C'},
		{"inputAtmosphericTransmissionFile",1, NULL, 'A'},
		{"inputObjectInTheSkyFile",         1, NULL, 'J'},
		{"CCDoption",                       1, NULL, 'c'},
		{"ExposureTime",                    1, NULL, 'e'},
        {"calculateExposureTime",           1, NULL, 'x'},
        {"SignalToNoise",                   1, NULL, 'n'},
        {"wavelength_nm",                   1, NULL, 'w'},
        {"Vmagnitude",                      1, NULL, 'V'},
        {"effectiveTemperature",            1, NULL, 'f'},
        {"imageQuality",                    1, NULL, 'Q'},
        {"airmass",                         1, NULL, 'm'},
        {"observingMode",                   1, NULL, 'B'},
        {"CCDReadoutMode",                  1, NULL, 'r'},
        {"moonPhase",                       1, NULL, 'p'},
        {"angularDistFromMoon",             1, NULL, 'a'},
        {"zenithalDistOfMoon",              1, NULL, 'z'},
        
		{"ordernumber",			1, NULL, 'O'},
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},
		{"plotfilename",		1, NULL, 'P'},
		{"dataFilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},
		{"interactive",			0, NULL, 'I'},
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "s:W:R:G:o:H:i:E:T:C:J:A:c:e:x:n:w:V:f:Q:m:B:r:p:a:z:O:M:X:P:F:C:S:I::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt)
		{
			case 's':
				inputSpectrum = optarg;
				break;
			case 'W':
				inputWaveFile = optarg;
				break;
			case 'R':
				inputApertureFile = optarg;
				break;
			case 'G':
				inputGainFile = optarg;
				break;
			case 'o':
				outputEspadonsETCFile = optarg;
				break;
			case 'H':
				outputHTMLEETCFile = optarg;
				break;
			case 'i':
				inputETCCalibrationFile = optarg;
				break;
			case 'E':
				inputSpectrographFile = optarg;
				break;
			case 'T':
				inputTelescopeFile = optarg;
				break;
			case 'C':
				inputObservingConditionsFile = optarg;
				break;
			case 'J':
				inputObjectInTheSkyFile = optarg;
				break;
			case 'A':
				inputAtmosphericTransmissionFile = optarg;
				break;                
			case 'c':
				CCDoption = (spectrographCCD_t)atoi(optarg);
				break;
			case 'e':
				ExposureTime = atof(optarg);
				break;
			case 'x':
				calculateExposureTime = (atoi(optarg)?true:false);
				break;
			case 'n':
				SignalToNoise = atof(optarg);
                providedSNR = true;
				break;
			case 'w':
				wavelength_nm = atof(optarg);
                providedWavelength = true;
				break;
			case 'V':
				Vmagnitude = atof(optarg);
				break;
			case 'f':
				effectiveTemperature = atof(optarg);
				break;
			case 'Q':
				imageQuality = atof(optarg);
				break;
			case 'm':
				airmass = atof(optarg);
				break;

            case 'B':
				observingMode = (EspadonsInstrumentMode_t)atoi(optarg);
				break;
			case 'r':
				CCDReadoutMode = (EspadonsCCDReadoutSpeed_t)atoi(optarg);
				break;
			case 'p':
				moonPhase = (moonphase_t)atoi(optarg);
				break;
			case 'a':
				angularDistFromMoon = atof(optarg);
				break;
			case 'z':
				zenithalDistOfMoon = atof(optarg);
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
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break;
			case 'F':
				dataFilename = optarg;
				break;
			case 'S':
				scriptfilename = optarg;
				break;
			case 'I':		// for interactive plots
				interactive = true;
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
        // we need an input ETC calibration file...
		if (inputETCCalibrationFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need a Spectrograph file...
		if (inputSpectrographFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a Telescope file...
		if (inputTelescopeFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
       
		// we need an atmospheric transmission file...
		if (inputAtmosphericTransmissionFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

        // we need a reference spectrum (*.e)...
		if (inputSpectrum.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need an aperture file...
		if (inputApertureFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need a gain/noise/bias file...
		if (inputGainFile.empty()) {
			throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

		if (verbose) {
            cout << "operaEspadonsETC: inputSpectrum = " << inputSpectrum << endl;
            cout << "operaEspadonsETC: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaEspadonsETC: inputApertureFile = " << inputApertureFile << endl;
            cout << "operaEspadonsETC: inputGainFile = " << inputGainFile << endl;
            cout << "operaEspadonsETC: outputEspadonsETCFile = " << outputEspadonsETCFile << endl;
            cout << "operaEspadonsETC: outputHTMLEETCFile = " << outputHTMLEETCFile << endl;
            cout << "operaEspadonsETC: inputETCCalibrationFile = " << inputETCCalibrationFile << endl;
            cout << "operaEspadonsETC: inputSpectrographFile = " << inputSpectrographFile << endl;
            cout << "operaEspadonsETC: inputTelescopeFile = " << inputTelescopeFile << endl;
            cout << "operaEspadonsETC: inputObservingConditionsFile = " << inputObservingConditionsFile << endl;
            cout << "operaEspadonsETC: inputObjectInTheSkyFile = " << inputObjectInTheSkyFile << endl;
            cout << "operaEspadonsETC: inputAtmosphericTransmissionFile = " << inputAtmosphericTransmissionFile << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaEspadonsETC: ordernumber = " << ordernumber << endl;
            }
            if(plot) {
                cout << "operaEspadonsETC: plotfilename = " << plotfilename << endl;
                cout << "operaEspadonsETC: dataFilename = " << dataFilename << endl;
                cout << "operaEspadonsETC: scriptfilename = " << scriptfilename << endl;
                if(interactive) {
                    cout << "operaEspadonsETC: interactive = YES" << endl;
                } else {
                    cout << "operaEspadonsETC: interactive = NO" << endl;
                }
            }
            
		}
        ofstream *fdata = NULL;
        
        if (!dataFilename.empty()) {
            fdata = new ofstream();
            fdata->open(dataFilename.c_str());
        }
        
        FILE *htmlfile = NULL;
        
        if (!outputHTMLEETCFile.empty()) {
            fopen(outputHTMLEETCFile.c_str(),"w");
        }

        operaInstrumentEnvironmentSetup instrumentEnvironment(inputSpectrographFile);
        instrumentEnvironment.ReadInstrumentEnvironmentSetup(inputTelescopeFile);
        operaSpectrograph *Spectrograph = instrumentEnvironment.getSpectrograph();
        operaTelescope *Telescope = instrumentEnvironment.getTelescope();

        operaObservingConditions *ObservingConditions = NULL;
        operaObjectInTheSky *ObjectInTheSky = NULL;
        
		if (!inputObservingConditionsFile.empty()) {
            instrumentEnvironment.ReadInstrumentEnvironmentSetup(inputObservingConditionsFile);
            
            ObservingConditions = instrumentEnvironment.getObservingConditions();
           
            ExposureTime = ObservingConditions->getexposureTime();
            imageQuality = ObservingConditions->getimageQuality();
            airmass = ObservingConditions->getairmass();
            moonPhase = ObservingConditions->getmoonphase();
            angularDistFromMoon = ObservingConditions->getangularDistFromMoon();
            zenithalDistOfMoon = ObservingConditions->getzenithalDistofMoon();
        } else {
            
            ObservingConditions = instrumentEnvironment.getObservingConditions();
            
            ObservingConditions->setexposureTime(ExposureTime);
            ObservingConditions->setimageQuality(imageQuality);
            ObservingConditions->setairmass(airmass);
            ObservingConditions->setmoonphase(moonPhase);
            ObservingConditions->setangularDistFromMoon(angularDistFromMoon);
            ObservingConditions->setzenithalDistofMoon(zenithalDistOfMoon);
        }
        
		if (!inputObjectInTheSkyFile.empty()) {
            instrumentEnvironment.ReadInstrumentEnvironmentSetup(inputObjectInTheSkyFile);
            
            ObjectInTheSky = instrumentEnvironment.getObjectInTheSky();
            
            Vmagnitude = ObjectInTheSky->getV_magnitude();
            effectiveTemperature = ObjectInTheSky->getEffectiveTemperature();
        } else {
            
            ObjectInTheSky = instrumentEnvironment.getObjectInTheSky();
            
            ObjectInTheSky->setV_magnitude(Vmagnitude);
            ObjectInTheSky->setEffectiveTemperature(effectiveTemperature);
        }
        
        if (CCDoption == undefinedCCD) { // no input CCDoption
            CCDoption = Spectrograph->getSpectrographCCD();
        }
        if (observingMode == undefinedInstrumentMode) { // no input observingMode
            observingMode = Spectrograph->getEspadonsInstrumentMode();
        }
        if (CCDReadoutMode == undefinedReadoutMode) { // no input CCDReadoutMode
            CCDReadoutMode = Spectrograph->getEspadonsCCDReadoutSpeed();
        }
        
   		if (verbose) {
            cout << "operaEspadonsETC: CCDoption = " << CCDoption << endl;
            cout << "operaEspadonsETC: ExposureTime = " << ObservingConditions->getexposureTime() << endl;
            cout << "operaEspadonsETC: calculateExposureTime = " << calculateExposureTime << endl;
            cout << "operaEspadonsETC: SignalToNoise = " << SignalToNoise << endl;
            cout << "operaEspadonsETC: wavelength_nm = " << wavelength_nm << endl;
            cout << "operaEspadonsETC: Vmagnitude = " << ObjectInTheSky->getV_magnitude() << endl;
            cout << "operaEspadonsETC: effectiveTemperature = " << ObjectInTheSky->getEffectiveTemperature() << endl;
            cout << "operaEspadonsETC: imageQuality = " << ObservingConditions->getimageQuality() << endl;
            cout << "operaEspadonsETC: airmass = " << ObservingConditions->getairmass() << endl;
            cout << "operaEspadonsETC: observingMode = " << observingMode << endl;
            cout << "operaEspadonsETC: CCDReadoutMode = " << CCDReadoutMode << endl;
            cout << "operaEspadonsETC: moonPhase = " << ObservingConditions->getmoonphase() << endl;
            cout << "operaEspadonsETC: angularDistFromMoon = " << ObservingConditions->getangularDistFromMoon() << endl;
            cout << "operaEspadonsETC: zenithalDistOfMoon = " << ObservingConditions->getzenithalDistofMoon() << endl;
		}
     
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputGainFile);        
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputApertureFile);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);
        
        //GainBiasNoise *gainBiasNoise = spectralOrders.getGainBiasNoise();
        
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
			cout << "operaEspadonsETC: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        
        /*
         * Calculate Exposure Time for given S/N at given wavelength_nm
         */
        if(calculateExposureTime && providedSNR && providedWavelength) {
            double wavelength_m = wavelength_nm*1e-9;
            double dwl_m = wavelength_m/Spectrograph->getspectralResolution();
            double wl0 = wavelength_m - dwl_m/2;
            double wlf = wavelength_m + dwl_m/2;
           
            double klamb = ObservingConditions->getAtmosphericTransmissionAtMaunaKea(inputAtmosphericTransmissionFile,wavelength_nm);
            double ipie = Spectrograph->CalculateIPIE(ObservingConditions->getimageQuality());
            double ccdqe = Spectrograph->CCDQuantumEfficiency(wavelength_nm);
            double fiberThru = Spectrograph->FiberThroughput(wavelength_nm);
            double opticsThru = Spectrograph->OpticsThroughput(wavelength_nm);
            
            double TotTrans = klamb*ipie*ccdqe*fiberThru*opticsThru;
            
            if(Spectrograph->getEspadonsInstrumentMode() == GRACES_staronly) {
                TotTrans *= GracesToCFHTRelativeTransmission(wavelength_nm);
            }
            
            double SrcEmit = ObjectInTheSky->getSpectralBinFlux(wl0,wlf);

            double SkyFlux = ObservingConditions->calculateSkyFlux(wavelength_nm,dwl_m,Spectrograph->getInjectionHoleDiameter());

            double ExpTime = (SignalToNoise*SignalToNoise)*((SrcEmit + SkyFlux)/pow(SrcEmit,2))/(Telescope->getCollectingArea()*TotTrans);
            
            ObservingConditions->setexposureTime(ExpTime);
        }
            
        if(htmlfile != NULL) {
            fprintf(htmlfile,"%s%c%c\n","Content-Type:text/html;charset=iso-8859-1",13,10);
            fprintf(htmlfile,"<html><head>\n");
            fprintf(htmlfile,"<meta http-equiv=\"content-type\" content=\"text/html; charset=ISO-8859-1\">\n");
            fprintf(htmlfile,"<title>ESPADONS CFHT Exposure Time Calculator </title>\n");
            fprintf(htmlfile,"</head><body background=\"etcform_files/bg.jpg\" bgcolor=\"Silver\">\n");
            fprintf(htmlfile,"<b><h2><font color=\"blue\">Parameters :</font></h2></b>\n");
            
            fprintf(htmlfile,"<font color=\"green\">Instrument mode</font>: %s <br>\n",Spectrograph->getFullModeName().c_str());
            fprintf(htmlfile,"<font color=\"green\">V magnitude</font>:  %.2f <br>\n", ObjectInTheSky->getV_magnitude());
            fprintf(htmlfile,"<font color=\"green\">Effective temperature</font>:  %.0f K<br>\n",ObjectInTheSky->getEffectiveTemperature());
            fprintf(htmlfile,"<font color=\"green\">Image Quality (IQ)</font>:  %.2f arcsec<br>\n", ObservingConditions->getimageQuality());
            fprintf(htmlfile,"<font color=\"green\">Airmass</font>:  %.2f <br>\n",ObservingConditions->getairmass());
            fprintf(htmlfile,"<font color=\"green\">Zenith dist of source</font>:  %.1f deg <br>\n",ObservingConditions->getZenithalDistanceOfSource());
            fprintf(htmlfile,"<font color=\"green\">CCD readout mode</font>: %s (noise=%.1fe, gain=%.2fe/adu, readout time=%.0fs) <br>\n",(Spectrograph->getCCDName()).c_str(), Spectrograph->getNominalNoise(), Spectrograph->getNominalGain(), Spectrograph->getNominalReadoutTime());
            fprintf(htmlfile,"<font color=\"green\">Moon</font>: %d <br>\n",ObservingConditions->getmoonphase());
            
            if( ObservingConditions->getangularDistFromMoon() > 0.25 &&  ObservingConditions->getangularDistFromMoon() <= 10 ) {
                fprintf(htmlfile,"<font color=\"green\">Angular distance from Moon</font>: %.0f deg<br>\n", ObservingConditions->getangularDistFromMoon());
                fprintf(htmlfile,"<font color=\"red\" size=\"large\"> WARNING! Object is less than 10 deg from the Moon.</font><br>\n");
            } else if( ObservingConditions->getangularDistFromMoon() >= 0 &&  ObservingConditions->getangularDistFromMoon() <= 0.25) {
                fprintf(htmlfile,"<font color=\"green\">Angular distance from Moon</font>: %.0f deg<br>\n", ObservingConditions->getangularDistFromMoon());
                fprintf(htmlfile,"<font color=\"red\" size=\"large\"> DISTANCE NOT ACCEPTED! Object is behind the Moon.</font><br>\n");
                fclose(htmlfile);
                throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            } else if( ObservingConditions->getangularDistFromMoon() < 0 ||  ObservingConditions->getangularDistFromMoon() > 180 ) {
                fprintf(htmlfile,"<font color=\"green\">Angular distance from Moon</font>: %.0f deg<br>\n", ObservingConditions->getangularDistFromMoon());
                fprintf(htmlfile,"<font color=\"red\" size=\"large\">ERROR: Angular distance from Moon must be between 0 and 180 deg.</font><br>\n");
                fclose(htmlfile);
                throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            } else {
                fprintf(htmlfile,"<font color=\"green\">Angular distance from Moon</font>: %.0f deg<br>\n", ObservingConditions->getangularDistFromMoon());
            }
            
            if(ObservingConditions->getzenithalDistofMoon() >= 0 && ObservingConditions->getzenithalDistofMoon() <= 90)
            {
                fprintf(htmlfile,"<font color=\"green\">Zenithal distance of Moon</font>:  %.0f deg<br>\n",ObservingConditions->getzenithalDistofMoon());
            } else {
                fprintf(htmlfile,"<font color=\"green\">Zenithal distance of Moon</font>:  %.0f deg<br>\n",ObservingConditions->getzenithalDistofMoon());
                fprintf(htmlfile,"<font color=\"red\" size=\"large\">ERROR: zenithal dist. of Moon must be between 0 and 90 deg.</font><br>\n");
                fclose(htmlfile);
                throw operaException("operaEspadonsETC: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            }
            fprintf(htmlfile,"<b><h2><font color=\"blue\">Results :</font></h2></b>\n");
            fprintf(htmlfile,"<font color=\"red\">CAUTION!!</font>  Average seeing over long exposures is probably not better than 0.8 arcsec. <br>\n");
            fprintf(htmlfile,"<font color=\"green\">Instrument pinhole injection efficiency</font>: %.6f <br>\n",Spectrograph->CalculateIPIE(ObservingConditions->getimageQuality()));
            fprintf(htmlfile,"<font color=\"green\">Readout noise per spectral bin</font>:  %.1f e<br><br>\n",(Spectrograph->getNominalNoise())*sqrt(1.*Spectrograph->getNominalAperture()*(double)(Spectrograph->getNumberOfExposures())));
            fprintf(htmlfile,"<font color=\"green\"> Sky V magnitude</font>:  %.1f mag/asec^2<br><br>\n",ObservingConditions->getSkyZeroVMagnitude(wavelength_nm));
            if(calculateExposureTime) {
                if(ObservingConditions->getexposureTime() > MAXIMUMEXPOSURETIME_S) {
                    ObservingConditions->setexposureTime(MAXIMUMEXPOSURETIME_S);
                    fprintf(htmlfile,"<font color=\"red\">WARNING: </font> Exposure time has been truncated to the maximum value.<br>\n");
                    fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                    if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                        fprintf(htmlfile,"<B><font color=\"green\">Exposure time per frame: </font> %.0f s <B><br>\n",ObservingConditions->getexposureTime());
                        fprintf(htmlfile,"<B><font color=\"red\">Total exposure time for 4 polarimetric frames: </font> %.0f s <B><br>\n",4*ObservingConditions->getexposureTime());
                    } else {
                        fprintf(htmlfile,"<B><font color=\"red\">Exposure time: </font> %.0f s <B><br>\n",ObservingConditions->getexposureTime());
                    }
                    fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
                } else if (ObservingConditions->getexposureTime() < 1) {
                    fprintf(htmlfile,"<font color=\"red\">WARNING:  Calculated exposure time is below 1 s (minimum accepted is 0.1 s).</font><br>\n");
                    fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                    
                    if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                        fprintf(htmlfile,"<B><font color=\"green\">Exposure time per frame: </font> %.2f s <B><br>\n",ObservingConditions->getexposureTime());
                        fprintf(htmlfile,"<B><font color=\"red\">Total exposure time for 4 polarimetric frames: </font> %.2f s <B><br>\n",4*ObservingConditions->getexposureTime());
                    } else {
                        fprintf(htmlfile,"<B><font color=\"red\">Exposure time: </font> %.2f s <B><br>\n",ObservingConditions->getexposureTime());
                    }
                    fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
                } else {
                    fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                    if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                        fprintf(htmlfile,"<B><font color=\"green\">Calculated exposure time per frame: </font> %.1f s <B><br>\n",ObservingConditions->getexposureTime());
                        fprintf(htmlfile,"<B><font color=\"red\">Total exposure time for 4 polarimetric frames: </font> %.1f s <B><br>\n",4*ObservingConditions->getexposureTime());
                    } else {
                        fprintf(htmlfile,"<B><font color=\"red\">Calculated exposure time: </font> %.1f s <B><br>\n",ObservingConditions->getexposureTime());
                    }
                    fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
                }
            } else {
                fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                    fprintf(htmlfile,"<B><font color=\"green\"> Exposure time per frame: </font> %.1f s <B><br>\n",ObservingConditions->getexposureTime());
                    fprintf(htmlfile,"<B><font color=\"red\"> Total exposure time for 4 polarimetric frames: </font> %.1f s <B><br>\n",4*ObservingConditions->getexposureTime());
                } else {
                    fprintf(htmlfile,"<B><font color=\"red\">  Exposure time: </font> %.1f s <B> <br>\n",ObservingConditions->getexposureTime());
                }
                fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
            }
            
            if (ObjectInTheSky->getV_magnitude() >= 0.0 && ObjectInTheSky->getV_magnitude() < 8) {
                
                double satexptime = Spectrograph->getSaturationExptime(ObjectInTheSky->getV_magnitude());
                
                if(ObservingConditions->getexposureTime() >= satexptime) {
                    fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                    fprintf(htmlfile,"<B><font color=\"red\"> WARNING: exposure will likely saturate! </font> <B><br>\n");
                    fprintf(htmlfile,"<B><font color=\"black\"> Source with V = %.2f mag saturates at EXPTIME = %.0f s<B><br>\n",ObjectInTheSky->getV_magnitude(),satexptime);
                    fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
                } else	{
                    fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                    fprintf(htmlfile,"<B><font color=\"black\"> Source with V = %.2f mag saturates at EXPTIME = %.0f s<B><br>\n",ObjectInTheSky->getV_magnitude(),satexptime);
                    fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
                }
            } else if (ObjectInTheSky->getV_magnitude() < 0.0) {
                fprintf(htmlfile,"<table border=\"1.5\" cellpadding=\"0\"> <tbody><tr> <td>");
                fprintf(htmlfile,"<B><font color=\"red\"> WARNING: exposure will likely saturate! </font> <B><br>\n");
                fprintf(htmlfile,"<B><font color=\"black\"> Source with V = 0.0 mag saturates at EXPTIME = 2 s<B><br>\n");
                fprintf(htmlfile,"</td></tr></tbody></table><br>\n");
            }
            
            fprintf(htmlfile,"<font color=\"blue\"><b>Spectral dependence of S/N:</b></font> <br>\n");
            fprintf(htmlfile,"<br>\n");
            fprintf(htmlfile,"<font color=\"green\">CCD pixel bin size (CCD pxl bin)</font>:  %.2f km/s<br>\n",2.6);
            fprintf(htmlfile,"<font color=\"green\">Spectral bin size (bin)</font>:  %.2f km/s<br>\n",1.8);
            fprintf(htmlfile,"<p align=\"CENTER\">\n");
            fprintf(htmlfile,"<table border=\"1\">\n");
            
            if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                fprintf(htmlfile,"<tbody><tr><td></td><td></td><td colspan=\"2\"><font color=\"blue\"><b>emitted photons</b></font></td><td colspan=\"2\"><font color=\"blue\"><b>transmission</b></font></td><td colspan=\"2\"><font color=\"blue\"><b>collected photons</b></font></td><td colspan=\"5\"><font color=\"blue\"><b>Results</b></font></td></tr>\n");
                fprintf(htmlfile,"<tr><td><b>order</b></td><td><b>wave</b></td><td><b>object</b></td><td><b>sky</b></td><td><b>atmosphere</b> </td><td><b>total</b></td><td><b>object</b></td><td><b>sky</b></td><td><b>S/N w/o readout</b></td><td colspan=\"2\"><b>S/N (per frame)</b></td><td colspan=\"2\"><b>S/N (4x polar)</b></td></tr>\n");
                fprintf(htmlfile,"<tr><td><br></td><td>(nm)</td><td>(ph/s/m2/bin)</td><td>(ph/s/m2/bin)</td><td><br></td><td><br></td><td>(ph/bin)</td><td>(ph/bin)</td><td>(/bin)</td><td>(/bin)</td><td>(/CCD pxl bin)</td><td>(/bin)</td><td>(/CCD pxl bin)</td></tr>\n");
            } else {
                fprintf(htmlfile,"<tbody><tr><td></td><td></td><td colspan=\"2\"><font color=\"blue\"><b>emitted photons</b></font></td><td colspan=\"2\"><font color=\"blue\"><b>transmission</b></font></td><td colspan=\"2\"><font color=\"blue\"><b>collected photons</b></font></td><td colspan=\"3\"><font color=\"blue\"><b>Results</b></font></td></tr>\n");
                fprintf(htmlfile,"<tr><td><b>order</b></td><td><b>wave</b></td><td><b>object</b></td><td><b>sky</b></td><td><b>atmosphere</b> </td><td><b>total</b></td><td><b>object</b></td><td><b>sky</b></td><td><b>S/N w/o readout</b></td><td colspan=\"2\"><b>S/N</b></td></tr>\n");
                fprintf(htmlfile,"<tr><td><br></td><td>(nm)</td><td>(ph/s/m2/bin)</td><td>(ph/s/m2/bin)</td><td><br></td><td><br></td><td>(ph/bin)</td><td>(ph/bin)</td><td>(/bin)</td><td>(/bin)</td><td>(/CCD pxl bin)</td></tr>\n");
            }            
        }
        
        /*
         * Calculate S/N for central wavelength of all orders
         */
		for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);

            if (spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasExtractionApertures() &&
                spectralOrder->gethasWavelength()) {
                
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                operaWavelength *Wavelength = spectralOrder->getWavelength();
                SpectralElements->setwavelengthsFromCalibration(Wavelength);
                
                double apertureAreaInPixelUnits = 0;
                for(unsigned beam=0;beam<spectralOrder->getnumberOfBeams();beam++) {
                    operaExtractionAperture<Line> *aperture = spectralOrder->getExtractionApertures(beam);
                    const PixelSet *aperturePixels = aperture->getSubpixels();
                    apertureAreaInPixelUnits += (double)aperturePixels->getNPixels()*(double)aperturePixels->getSubpixelArea();
                }
                
                unsigned centralIndex = SpectralElements->getnSpectralElements()/2;
                //unsigned xcenter = (unsigned)floor(SpectralElements->getphotoCenterX(centralIndex));
                //unsigned ycenter = (unsigned)floor(SpectralElements->getphotoCenterY(centralIndex));
                
                double centralWavelength_nm = SpectralElements->getwavelength(centralIndex);
                double centralWavelength_m = centralWavelength_nm*1e-9;
                double dwl_m = (1e-9)*(SpectralElements->getwavelength(centralIndex+1) + SpectralElements->getwavelength(centralIndex-1))/2.0;
                double wl0 = centralWavelength_m - dwl_m/2;
                double wlf = centralWavelength_m + dwl_m/2;

                double klamb = ObservingConditions->getAtmosphericTransmissionAtMaunaKea(inputAtmosphericTransmissionFile,centralWavelength_nm);
                double ipie = Spectrograph->CalculateIPIE(ObservingConditions->getimageQuality());                
                double ccdqe = Spectrograph->CCDQuantumEfficiency(centralWavelength_nm);
                double fiberThru = Spectrograph->FiberThroughput(centralWavelength_nm);
                
                double opticsThru = Spectrograph->OpticsThroughput(centralWavelength_nm);
                
                double TotTrans = klamb*ipie*ccdqe*fiberThru*opticsThru;
                
                if(Spectrograph->getEspadonsInstrumentMode() == GRACES_staronly) {
                    TotTrans *= GracesToCFHTRelativeTransmission(centralWavelength_nm);
                }
                
                double SrcEmit = ObjectInTheSky->getSpectralBinFlux(wl0,wlf);
                double SkyFlux = ObservingConditions->calculateSkyFlux(centralWavelength_nm,dwl_m,Spectrograph->getInjectionHoleDiameter());
                
                double SrcCole = SrcEmit*Telescope->getCollectingArea()*ObservingConditions->getexposureTime()*TotTrans;
                double SkyColFlux = SkyFlux*Telescope->getCollectingArea()*ObservingConditions->getexposureTime()*TotTrans;
                
                // noise per spectral bin (npsb) - just compute as it was 1 pixel x Spectrograph->getNominalAperture()
//                double npsb = (Spectrograph->getNominalNoise())*sqrt(1.*Spectrograph->getNominalAperture()*(double)(Spectrograph->getNumberOfExposures()));
//                double noise = gainBiasNoise->getNoise(xcenter,ycenter);
                double npsb = /*???DT Mar 2013 uninitialized npsb**/sqrt(apertureAreaInPixelUnits*(double)(Spectrograph->getNumberOfExposures()));
                double SNRnoread = SrcCole/sqrt(SrcCole + SkyColFlux);
                double SNR = SrcCole/sqrt(SrcCole + SkyColFlux + npsb*npsb);
                
                if(htmlfile != NULL) {
                    if (Spectrograph->getEspadonsInstrumentMode() == polarimetric) {
                        fprintf(htmlfile,"<tr><td> #%d </td><td> %.2f </td><td> %.3e </td><td> %.3e </td><td> %.3f </td><td> %.3f </td><td> %.3e </td><td> %.3e </td><td> %.0f </td><td><font color=\"red\"> %.0f </b></td><td><font color=\"red\"> %.0f </b></td><td><font color=\"red\"> %.0f </b></td><td><font color=\"red\"> %.0f </b></td></tr>\n",order,centralWavelength_nm,SrcEmit,SkyFlux,klamb,TotTrans,SrcCole,SkyColFlux,SNRnoread,SNR,SNR*1.2, 2*SNR, 2*SNR*1.2);
                    } else {
                        fprintf(htmlfile,"<tr><td> #%d </td><td> %.2f </td><td> %.3e </td><td> %.3e </td><td> %.3f </td><td> %.3f </td><td> %.3e </td><td> %.3e </td><td> %.0f </td><td><font color=\"red\"> %.0f </b></td><td><font color=\"red\"> %.0f </b></td></tr>\n",order,centralWavelength_nm,SrcEmit,SkyFlux,klamb,TotTrans,SrcCole,SkyColFlux,SNRnoread,SNR,SNR*1.2);
                    }
                }
            }
        }
        
        if(htmlfile != NULL) {
            fprintf(htmlfile,"</tbody></table>\n");
            fprintf(htmlfile,"</p>  </span></body></html>\n");
        }
        
		/*
		 * and write it out
		 */
		//spectralOrders.WriteSpectralOrders(outputEspadonsETCFile, eetc);
		
        
        if (fdata != NULL) {
			fdata->close();
            
            if (!scriptfilename.empty()) {
                //GenerateCreateFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),dataFilename.c_str(),continuumDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
        
        if (htmlfile != NULL) {
			fclose(htmlfile);
        }
        
	}
	catch (operaException e) {
		cerr << "operaEspadonsETC: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaEspadonsETC: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --inputSpectrum=<SPEC_FILE>"
    " --inputWaveFile=<WCAL_FILE>"
    " --inputApertureFile=<APER_FILE>"
    " --inputGainFile=<GAIN_FILE>"
    " --outputEspadonsETCFile=<EETC_FILE>"
    " --outputHTMLEETCFile=<HTML_FILE>"
    " --inputETCCalibrationFile=<Fcal_FILE>"
    " --inputSpectrographFile=<SPECTROGRAPH_FILE>"
    " --inputTelescopeFile=<TELESCOPE_FILE>"
    " --inputObservingConditionsFile=<OBSCOND_FILE>"
    " --inputObjectInTheSkyFile=<SKYOBJ_FILE>"
    " --inputAtmosphericTransmissionFile=<ATMTRANS_FILE>"    
    " --CCDoption=<spectrographCCD_t>"
    " --ExposureTime=<DBL_VALUE>"
    " --calculateExposureTime=<BOOL>"
    " --SignalToNoise=<DBL_VALUE>"
    " --wavelength_nm=<DBL_VALUE>"
    " --Vmagnitude=<DBL_VALUE>"
    " --effectiveTemperature=<DBL_VALUE>"
    " --imageQuality=<DBL_VALUE>"
    " --airmass=<DBL_VALUE>"
    " --observingMode=<EspadonsInstrumentMode_t>"
    " --CCDReadoutMode=<EspadonsCCDReadoutSpeed_t>"
    " --moonPhase=<moonphase_t>"
    " --angularDistFromMoon=<DBL_VALUE>"
    " --zenithalDistOfMoon=<DBL_VALUE>"    
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --plotfilename=<EPS_FILE>"
	" --dataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" -v \n\n"
	"  -h, --help,  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -s, --inputSpectrum=<SPEC_FILE>, Reference spectrum \n"
	"  -W, --inputWaveFile=<WCAL_FILE>, Input wavelength calibration file\n"
	"  -R, --inputApertureFile=<APER_FILE>, Input aperture calibration file\n"
	"  -G, --inputGainFile=<GAIN_FILE>, Input gain/noise/bias calibration file\n"
    "  -o, --outputEspadonsETCFile, Output Espadons ETC data file\n"
    "  -H, --outputHTMLEETCFile, Output Espadons ETC html file\n"
    "  -i, --inputETCCalibrationFile, Input ETC calibration file\n"
    "  -E, --inputSpectrographFile, Input spectrograph file\n"
    "  -T, --inputTelescopeFile, Input telescope file\n"
    "  -C, --inputObservingConditionsFile, Input observing conditions file\n"
    "  -J, --inputObjectInTheSkyFile, Input sky object file\n"
    "  -A, --inputAtmosphericTransmissionFile, Atmospheric transmission spectrum data file (2 cols: lambda(A), trans)\n"
    "  -c, --CCDoption, Supported options: Olapa, EEV1\n"
    "  -e, --ExposureTime, Given exposure time (in seconds) to calcule S/N table \n"
    "  -x, --calculateExposureTime, Boolean to calculate exposure time for a given S/N \n"
    "  -n, --SignalToNoise, Desired signal-to-noise\n"
    "  -w, --wavelength_nm, Central wavelength of order for desired S/N \n"
    "  -V, --Vmagnitude, V-band (Johnson) magnitude\n"
    "  -f, --effectiveTemperature, Black body temperature in K\n"
    "  -Q, --imageQuality, Seeing (in arcsec)\n"
    "  -m, --airmass, \n"
    "  -B, --observingMode, Supported options: polarimetric, staronly, starplussky, GRACES_staronly, GRACES_starplussky\n"
    "  -r, --CCDReadoutMode, Supported options: slowmode, normalmode, fastmode\n"
    "  -p, --moonPhase, Supported options: newmoon, crescentmoon, quartermoon, gibbousmoon, fullmoon\n"
    "  -a, --angularDistFromMoon, Angular distance between target and the Moon (in degrees)\n"
    "  -z, --zenithalDistOfMoon, Zenithal distance of Moon (in degrees)\n"
	"  -O, --ordernumber=<UNS_VALUE>, Absolute order number to extract (default=all)\n"
	"  -M, --minorder=<UNS_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<UNS_VALUE>, Define maximum order number\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --dataFilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateCreateFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    *fgnu << "set ylabel \"flux\"" << endl;
    
    *fgnu << "set pointsize 1.0" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << dataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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
        *fgnu << "\nplot \"" << dataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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

double GracesToCFHTRelativeTransmission(double wavelength_nm) {
    
    double wl[40],graces2cfht[40];
    
    wl[0]=372; graces2cfht[0]=1.04535;
    wl[1]=378; graces2cfht[1]=1.05484;
    wl[2]=384; graces2cfht[2]=1.07185;
    wl[3]=391; graces2cfht[3]=1.09912;
    wl[4]=398; graces2cfht[4]=1.12494;
    wl[5]=405; graces2cfht[5]=1.15065;
    wl[6]=412; graces2cfht[6]=1.17407;
    wl[7]=420; graces2cfht[7]=1.19989;
    wl[8]=428; graces2cfht[8]=1.23031;
    wl[9]=436; graces2cfht[9]=1.25723;
    wl[10]=444; graces2cfht[10]=1.28354;
    wl[11]=453; graces2cfht[11]=1.31406;
    wl[12]=463; graces2cfht[12]=1.34038;
    wl[13]=472; graces2cfht[13]=1.37431;
    wl[14]=482; graces2cfht[14]=1.40745;
    wl[15]=493; graces2cfht[15]=1.44485;
    wl[16]=504; graces2cfht[16]=1.47930;
    wl[17]=515; graces2cfht[17]=1.51443;
    wl[18]=527; graces2cfht[18]=1.55196;
    wl[19]=540; graces2cfht[19]=1.59336;
    wl[20]=553; graces2cfht[20]=1.62862;
    wl[21]=567; graces2cfht[21]=1.66214;
    wl[22]=581; graces2cfht[22]=1.70047;
    wl[23]=597; graces2cfht[23]=1.73253;
    wl[24]=613; graces2cfht[24]=1.76547;
    wl[25]=630; graces2cfht[25]=1.80951;
    wl[26]=648; graces2cfht[26]=1.83806;
    wl[27]=667; graces2cfht[27]=1.87160;
    wl[28]=687; graces2cfht[28]=1.91638;
    wl[29]=709; graces2cfht[29]=1.95915;
    wl[30]=732; graces2cfht[30]=2.00700;
    wl[31]=756; graces2cfht[31]=2.08058;
    wl[32]=782; graces2cfht[32]=2.13913;
    wl[33]=810; graces2cfht[33]=2.19853;
    wl[34]=840; graces2cfht[34]=2.21699;
    wl[35]=872; graces2cfht[35]=2.20332;
    wl[36]=907; graces2cfht[36]=2.10251;
    wl[37]=945; graces2cfht[37]=1.97123;
    wl[38]=986; graces2cfht[38]=2.04329;
    wl[39]=1031; graces2cfht[39]=2.10000;

    unsigned nin = 40;
    
    double yp1 = (graces2cfht[1] - graces2cfht[0])/(wl[1] - wl[0]);
	double ypn = (graces2cfht[nin-1] - graces2cfht[nin-2])/(wl[nin-1] - wl[nin-2]);
	double *y2 = (double *)malloc(nin*sizeof(double));
	
	// Call cubicspline to get second derivatives
	cubicsplineDouble(wl, graces2cfht, nin, yp1, ypn, y2);
	
    double outputgraces2cfht;
    
	// Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, graces2cfht, y2, nin, wavelength_nm , &outputgraces2cfht);
    
    free(y2);
    
    return outputgraces2cfht;
}


