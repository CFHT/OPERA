/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaApplyTransform.cpp
 Version: 1.0
 Description: Apply a transform to derive an output from an m.fits product.
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

/*! \file operaApplyTransform.cpp */

using namespace std;

/*!
 * operaApplyTransform
 * \author Doug Teeple
 * \brief Apply a transform to derive an output from an m.fits product.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup tools
 */

#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>

#include <fstream>
#include <string>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/gzstream.h"

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] <list of m.fits, i.fits, p.fits images>\n\n"
	"  Applies a transform to intensity or polar data stored in the fits products,\n"
	"  or extracts a calibration.\n"
	"  -y, --spectrumtype=...             Type of output: (default raw spectrum)\n"
    "    RawSpectrum = 1                - uncalibrated / standard vectors of x,y coordinates and flux\n"
    "    CalibratedSpectrum = 9         - wavelength calibrated / standard vectors of wl and flux (default)\n"
    "    LibreEspritSpectrum = 17       - wavelength calibrated / LE vectors of wl and flux\n"
    "    LibreEspritpolarimetry = 21    - wavelength calibrated / LE vectors of wl and polar\n"
    "    Geom = 23                      - geometry vectors / polynomial\n"
    "    Wave = 24                      - wavelength vectors / polynomial\n"
    "    Spec = 25                      - spectral element vectors / polynomial\n"
    "    Prof = 26                      - instrument profile vectors / polynomial\n"
    "    Disp = 27                      - Dispersion polynomial\n"
    "    SNR = 28                       - SNR for each order\n"
    "    Polarimetry = 29               - polarimetry\n"
    "    Orderspacing = 30              - the geometry order spacing polynomial\n"
    "    GainNoise = 32                 - gain noise calibration data\n"
    "    Aperture = 33                  - aperture data\n"
    "    Fcal = 34                      - flux calibration data\n"
    "    RVel = 35                      - Radial Velocity Correction\n"
    "    Tell = 37                      - Telluric wave Correction\n"
    "\n"
	"  -o, --object=<object name>, for Libre-Esprit output (optional can be gotten from headers)\n"
	"  -w, --wavelengthCalibration=1|0, wavelength calibration polynomials\n"
	"  -V, --ApplyRadialVelocityCorrection=1|0, Barycentric Radial Velocity Correction\n"
	"  -N, --normalize=1|0, apply flux normalization\n"
	"  -b, --normalizationBinsize=<float>,  binsize for normalization\n"
	"  -B, --orderBin=<int>, number or orders to bin for continuum evaluation\n"
	"  -A, --AbsoluteCalibration=<bool>, perform absolute flux calibration\n"    
	"  -l, --usePolynomial=1|0, option to use polynomial for normalization\n"
	"  -r, --orderOfPolynomial=<unsigned>, option to set degree of polynomial for normalization\n"
    "\n"
	"  -h, --help,	   Display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,    Turn on debug messages\n"
	"  -t, --trace,    Turn on trace messages\n"
	"  -p, --plot,     Plot output \n"
	"\n";
}

int main(int argc, char *argv[])
{
	int opt;
	int debug=0, verbose=0, trace=0, plot=0;
    
	string productname;
	string basefilename;
	string object;
	string directory = "./";
	string zipped = ".gz";
	string mode;
	string speed;
	string detector;
	string amplifier;
	string oset = "";
	operaSpectralOrder_t spectralOrderType = RawBeamSpectrum;
    /*
     * Parameters for normalization
     */
    bool normalize = false;
    unsigned normalizationBinsize = 100;
    bool usePolynomial = FALSE;
    bool ApplyTelluricCorrection = FALSE;
    unsigned orderOfPolynomial = 5;
    
    /*
     * Barycentric radial velocity correction
     */
	bool ApplyRadialVelocityCorrection = false;
	
    /*
     * Parameters for flux calibration
     */
    bool ApplyFluxCorrection = false;
	float exposureTime = 0.0;
    bool AbsoluteCalibration = false;
    int orderBin = 2;
    
	bool ApplyWavelengthCorrection = false;
    
	struct option longopts[] = {
		
		{"spectrumtype",				1, NULL, 'y'},	// spectrum type
		{"ApplyWavelengthCorrection",		1, NULL, 'w'},	// wavelength calibration file (.wcal)
 		{"ApplyRadialVelocityCorrection",	1, NULL, 'V'},  // Barycentric wavelength correction file (.rvel)
 		{"telluriccorrection",			1, NULL, 'T'},  // wavelength calibration file (.tell)
		
		{"normalize",					1, NULL, 'N'},	// apply flux normalization
 		{"normalizationBinsize",		1, NULL, 'b'},	// binsize for normalization
		{"usePolynomial",				1, NULL, 'l'},	// option to use polynomial for normalization
		{"orderOfPolynomial",			1, NULL, 'r'},	// option to set degree of polynomial for normalization
		
        {"fluxCalibration",				1, NULL, 'C'},	// apply flux calibration; file (.fcal)
        {"etime",						1, NULL, 'E'},	// needed for flux calibration
		{"orderBin",                    1, NULL, 'B'},  // needed for flux calibration
		{"AbsoluteCalibration",         1, NULL, 'A'},  // absolute or relative flux calibration

		{"object",						1, NULL, 'o'},	// needed for Libre-Esprit output
		
		{"plot",		optional_argument,	NULL, 'p'},
		{"verbose",		optional_argument,	NULL, 'v'},
		{"debug",		optional_argument,	NULL, 'd'},
		{"trace",		optional_argument,	NULL, 't'},
		{"help",		no_argument,		NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "y:w:V:T:N:b:l:r:C:E:B:A:o:p::v::d::t::h",
							 longopts, NULL))  != -1)
        {
		switch(opt)
            {
                case 'y':		// spectrum type
                    spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
                    switch (spectralOrderType) {
                        case RawSpectrum:
                        case CalibratedRawSpectrum:
                        case LibreEspritSpectrum:
                        case LibreEspritpolarimetry:
                        case Prof:
                        case Geom:
                        case Wave:
                        case Spec:
                        case Disp:
                        case SNR:
                        case Polarimetry:
                        case GainNoise:
                        case Aperture:
                        case Fcal:
                        case RVel:
                        case Tell:
                            break;
                            
                        default:
                            printUsageSyntax(argv[0]);
                            exit(EXIT_SUCCESS);
                            break;
                    }
				break;
                case 'w':
                    ApplyWavelengthCorrection = true;
				break;
                case 'V':		// for telluric wl correction
                    ApplyRadialVelocityCorrection = true;
				break;
                case 'T':		// for telluric wl correction
                    ApplyTelluricCorrection = true;
				break;
                
                case 'N':		// for normalization
                    normalize = true;
				break;
                case 'b':		// normalization binsize
                    normalizationBinsize = atoi(optarg);
				break;
                case 'l':
                    usePolynomial = true;
				break;
                case 'r':
                    orderOfPolynomial = atoi(optarg);
				break;
                
                case 'C':       // for flux calibration
                    ApplyFluxCorrection = true;
				break;
                case 'E':
                    exposureTime = atof(optarg);
				break;
                case 'B':
                    orderBin = atoi(optarg);
                    break;
                case 'A':
                    AbsoluteCalibration = true;
                    break;
                case 'o':
                object = optarg;
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
		while (optind < argc) {
			basefilename = productname = string(argv[optind++]);
            if (verbose) {
                cout << "operaApplyTransform: object = " << object << endl;
                cout << "operaApplyTransform: spectrum type = " << spectralOrderType << endl;
                cout << "operaApplyTransform: wavelength cdorrection = " << ApplyWavelengthCorrection << endl;
                cout << "operaApplyTransform: binsize for normalization = " << normalizationBinsize << endl;
                cout << "operaApplyTransform: apply flux calibration file = " << ApplyFluxCorrection << endl;
                cout << "operaApplyTransform: exposure time = " << exposureTime << endl;
                cout << "operaApplyTransform: order bin = " << orderBin << endl;
                cout << "operaApplyTransform: absolute calibration = " << AbsoluteCalibration << endl;
                cout << "operaApplyTransform: ApplyRadialVelocityCorrection = " << ApplyRadialVelocityCorrection << endl;
            }
            if (basefilename.find_last_of("/") != string::npos) {
                basefilename = basefilename.substr(basefilename.find_last_of("/")+1);
                if (directory.empty()) {
                    directory = productname.substr(0, productname.find_last_of("/")+1);
                }
            }
            if (verbose) {
                cout << "operaApplyTransform: productname=" << productname << endl;
                cout << "operaApplyTransform: basefilename=" << basefilename << endl;
                cout << "operaApplyTransform: directory=" << directory << endl;
            }
            if (productname.find("m.fits") != string::npos) {
                operaMEFFITSProduct in(productname, READONLY);
                basefilename = basefilename.substr(0, basefilename.find("m.fits"));
                if (mode.empty()) {
                    string rawmode = in.operaFITSGetHeaderValue("INSTMODE", 0);
                    if (rawmode.find("Polarimetry,") != string::npos) {
                        mode = "pol";
                    } else if (rawmode.find("star+sky,") != string::npos) {
                        mode = "sp1";
                    } else {
                        mode = "sp2";
                    }
                }
                if (speed.empty()) {
                    speed = in.operaFITSGetHeaderValue("EREADSPD", 0);
                    speed = speed.substr(0, speed.find(":"));
                }
                if (detector.empty()) {
                    detector = in.operaFITSGetHeaderValue("DETECTOR", 0);
                }
                if (amplifier.empty()) {
                    try {	// optional keyword, doesn't always exist
                        amplifier = in.operaFITSGetHeaderValue("AMPLIST", 0);
                        if (amplifier.find(",") != string::npos) {
                            amplifier.erase(amplifier.find(","));
                        }
                    } catch (...) {
                        amplifier = "";
                    }
                }
                if (object.empty()) {
                    object = in.operaFITSGetHeaderValue("OBJECT", 0);
                }
                if (verbose) {
                    cout << "operaApplyTransform: mode=" << mode << endl;
                    cout << "operaApplyTransform: speed=" << speed << endl;
                    cout << "operaApplyTransform: detector=" << detector << endl;
                    cout << "operaApplyTransform: amplifier=" << amplifier << endl;
                    cout << "operaApplyTransform: object=" << object << endl;
                }
                in.operaFITSImageClose();
                
                string calfilename = detector + amplifier + "_" + mode + oset + "_" + speed;
                string outfilename;
                
                operaSpectralOrderVector spectralOrders;
                
                switch (spectralOrderType) {
                        /*
                         * These cases should extract the beam spectrum or polarimetry
                         * and apply the requested calibrations.
                         */
                    case RawSpectrum:
                    case StandardSpectrum:
                    case OperaOptimalSpectrum:
                    case RawBeamSpectrum:
                    case StandardBeamSpectrum:
                    case OptimalBeamSpectrum:
                    case OperaOptimalBeamSpectrum:
                    {
                        outfilename = directory + basefilename + "i.e" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, RawBeamSpectrum))
                            spectralOrders.WriteSpectralOrders(outfilename, spectralOrderType);
                    }
                        break;
                    case CalibratedRawSpectrum:
                    case CalibratedStandardSpectrum:
                    case CalibratedOperaOptimalSpectrum:
                    case CalibratedRawBeamSpectrum:
                    case CalibratedStandardBeamSpectrum:
                    case CalibratedOptimalBeamSpectrum:
                    case CalibratedOperaOptimalBeamSpectrum:
                    {
                        outfilename = directory + basefilename + "i.e" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        

                        
                        if (spectralOrders.ReadSpectralOrders(productname, RawBeamSpectrum)) {
                            int minorder = spectralOrders.getMinorder();
                            int maxorder = spectralOrders.getMaxorder();
                            if (ApplyRadialVelocityCorrection) {
                                ApplyRadialVelocityCorrection = spectralOrders.ReadSpectralOrders(productname, RVel);
                            }
                            if (ApplyWavelengthCorrection) {
                                ApplyWavelengthCorrection = spectralOrders.ReadSpectralOrders(productname, Wave);
                            }
                            if (ApplyTelluricCorrection) {
                                ApplyTelluricCorrection = spectralOrders.ReadSpectralOrders(productname, Tell);
                            }
                            
                            /*
                             * Flux Calibration Stuff ...
                             */
                            unsigned NumberofBeams = 0;
                            double uncalibratedContinuumFluxForNormalization = 0;
                            double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS];
                            double spectralBinConstant = 0;
                            double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
                            
                            if (ApplyFluxCorrection) {
                                ApplyFluxCorrection = spectralOrders.ReadSpectralOrders(productname, Fcal);
                                for (int order=minorder; order<=maxorder; order++) {
                                    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                                    if (spectralOrder->gethasSpectralElements()) {
                                        NumberofBeams = spectralOrder->getnumberOfBeams();
                                        break;
                                    }
                                }
                                unsigned nsigcut = 3;
                                spectralOrders.getContinuumFluxesForNormalization(&uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization,normalizationBinsize, orderBin, nsigcut);
                                spectralBinConstant = exposureTime;
                                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                                    BeamSpectralBinConstant[beam] = exposureTime;
                                }
                            }
                            for (int order=minorder; order<=maxorder; order++) {
                                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                                if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
                                    
                                    if (mode == "sp1") {
                                        double SkyOverStarFiberAreaRatio = (2.2*2.2)/(1.6*1.6);
                                        for(unsigned beam = 0; beam < NumberofBeams; beam++) {
                                            if(beam < floor(float(NumberofBeams)/2)) {
                                                BeamSpectralBinConstant[beam] = exposureTime;
                                            } else if(float(beam) >= float(NumberofBeams)/2) { // Sky Fiber
                                                BeamSpectralBinConstant[beam] = exposureTime/SkyOverStarFiberAreaRatio;
                                            }
                                        }
                                        operaFluxVector skyFlux(spectralOrder->getSpectralElements()->getnSpectralElements());
                                        operaFluxVector starPlusSkyFlux(spectralOrder->getSpectralElements()->getnSpectralElements());
                                        // deprecated spectralOrder->calculateStarAndSkyElements(NULL);
                                    }
                                    operaWavelength *wavelength = spectralOrder->getWavelength();
                                    Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
                                    operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                                    
                                    if(normalize) {
                                        spectralOrder->applyNormalization(normalizationBinsize, orderOfPolynomial, usePolynomial, NULL, NULL, TRUE, 0);
                                    } else if (!normalize && ApplyFluxCorrection && spectralOrder->gethasSpectralEnergyDistribution()) {
                                        //breaking things in this module for the greater good. spectralOrder->applyFluxCalibration(spectralBinConstant, BeamSpectralBinConstant, uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization, AbsoluteCalibration, NULL);
                                    }
                                    unsigned elements = (unsigned)ceil((float)spectralElements->getnSpectralElements());
                                    while (elements--) {
                                        spectralElements->setwavelength(wavelengthPolynomial->Evaluate(spectralElements->getdistd(elements)), elements);
                                    }
                                    spectralElements->setHasWavelength(true);
                                    spectralElements->setHasDistance(false);
                                    if (ApplyRadialVelocityCorrection) {
                                        spectralOrder->applyRvelWavelengthCorrection(spectralOrders.getRadialVelocityCorrection());
                                    }
                                }
                            }
                            // output a wavelength calibrated spectrum...
                            spectralOrders.setObject(object);
                            spectralOrders.WriteSpectralOrders(outfilename, spectralOrderType);
                        }
                    }
                        break;
                    case LibreEspritSpectrum:
                    case LibreEspritsp1Spectrum:
                    case LibreEspritsp2Spectrum:
                    case LibreEspritpolSpectrum:
                    {
                        outfilename = directory + basefilename + "i.s";
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        spectralOrders.setObject(object);
                        if (spectralOrders.ReadSpectralOrders(productname, RawBeamSpectrum))
                            spectralOrders.WriteSpectralOrders(outfilename, spectralOrderType);
                    }
                        break;
                    case LibreEspritpolarimetry:
                    {
                        outfilename = directory + basefilename + "p.s";
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        spectralOrders.setObject(object);
                        if (spectralOrders.ReadSpectralOrders(productname, Polarimetry))
                            spectralOrders.WriteSpectralOrders(outfilename, Polarimetry);
                    }
                        break;
                    case Polarimetry:
                    {
                        outfilename = directory + basefilename + "p.e" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Polarimetry))
                            spectralOrders.WriteSpectralOrders(outfilename, Polarimetry);
                    }
                        break;
                        /*
                         * These cases just extract the calibration
                         */
                    case Geom:
                    {
                        outfilename = directory + calfilename + ".geom" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Geom))
                            spectralOrders.WriteSpectralOrders(outfilename, Geom);
                    }
                        break;
                    case Wave:
                    {
                        outfilename = directory + calfilename + ".wcal" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Wave))
                            spectralOrders.WriteSpectralOrders(outfilename, Wave);
                    }
                        break;
                    case Prof:
                    {
                        outfilename = directory + calfilename + ".prof" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Prof))
                            spectralOrders.WriteSpectralOrders(outfilename, Prof);
                    }
                        break;
                    case Disp:
                    {
                        outfilename = directory + calfilename + ".disp" + zipped;
                        if (spectralOrders.ReadSpectralOrders(productname, Disp))
                            spectralOrders.WriteSpectralOrders(outfilename, Disp);
                    }
                        break;
                    case SNR:
                    {
                        outfilename = directory + calfilename + ".sn" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, SNR))
                            spectralOrders.WriteSpectralOrders(outfilename, SNR);
                    }
                        break;
                    case Orderspacing:
                    {
                        outfilename = directory + calfilename + ".ordp" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Orderspacing))
                            spectralOrders.WriteSpectralOrders(outfilename, Orderspacing);
                    }
                        break;
                    case GainNoise:
                    {
                        outfilename = directory + calfilename + ".gain" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, GainNoise))
                            spectralOrders.WriteSpectralOrders(outfilename, GainNoise);
                    }
                        break;
                    case Aperture:
                    {
                        outfilename = directory + calfilename + ".aper" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Aperture))
                            spectralOrders.WriteSpectralOrders(outfilename, Aperture);
                    }
                        break;
                    case Fcal: 
                    {
                        outfilename = directory + calfilename + ".fcal" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Fcal))
                            spectralOrders.WriteSpectralOrders(outfilename, Fcal);
                    }
                        break;
                    case RVel: 
                    {
                        outfilename = directory + calfilename + ".rvel" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, RVel))
                            spectralOrders.WriteSpectralOrders(outfilename, RVel);
                    }
                        break;
                    case PRVel: 
                    {
                        outfilename = directory + calfilename + ".rvel" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, RVel))
                            spectralOrders.WriteSpectralOrders(outfilename, RVel);
                    }
                        break;
                    case Tell: 
                    {
                        outfilename = directory + calfilename + ".tell" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Tell))
                            spectralOrders.WriteSpectralOrders(outfilename, Tell);
                    }
                        break;
                    case PTell: 
                    {
                        outfilename = directory + calfilename + ".tell" + zipped;
                        if (verbose) {
                            cout << "operaApplyTransform: " << outfilename << endl;
                        }
                        if (spectralOrders.ReadSpectralOrders(productname, Tell))
                            spectralOrders.WriteSpectralOrders(outfilename, Tell);
                    }
                        break;
                    default:
                        printUsageSyntax(argv[0]);
                        exit(EXIT_SUCCESS);
                        break;
                }
            } else if (productname.find("p.fits") != string::npos) {
                if (verbose) {
                    cout << "operaApplyTransform: dir=" << directory << endl;
                }
                operaFITSProduct in(productname, READONLY);
                basefilename = basefilename.substr(0, basefilename.find("p.fits"));
                string outfilenamebase = directory + basefilename;
                instrumentmode_t instrumentmode;
                string object;
                string mode = in.operaFITSGetHeaderValue("INSTMODE");
                if (in.getnaxis1() > in.getnaxis2()) {
                    in.rotate90();
                }
                unsigned rows = (unsigned)in.getYDimension();
                unsigned columns = 5;
                if (mode.find("Polarimetry") != string::npos) {
                    instrumentmode = MODE_POLAR;
                } else {
                    throw operaException("operaApplyTransform: "+mode+' ', operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
                }
                /*
                 * pu
                 */
                if (!normalize && !ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "pu.s";
                    fout.open(outfilename.c_str());
                    object = in.operaFITSGetHeaderValue("OBJECT");
                    if (verbose) {
                        cout << "operaApplyTransform: mode=" << mode << endl;
                        cout << "operaApplyTransform: object='" << object << "'"<< endl;
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    for (unsigned row=0; row<rows; row++) {
                        fout << in[row][18] << ' ';
                        fout << in[row][19] << ' ';
                        fout << in[row][20] << ' ';
                        fout << in[row][21] << ' ';
                        fout << in[row][22] << ' ';
                        fout << in[row][23] << ' ';
                        fout << endl;
                    }
                    fout.close();
                }
                /*
                 * pn
                 */
                if (normalize && !ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "pn.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    for (unsigned row=0; row<rows; row++) {
                        fout << in[row][12] << ' ';
                        fout << in[row][13] << ' ';
                        fout << in[row][14] << ' ';
                        fout << in[row][15] << ' ';
                        fout << in[row][16] << ' ';
                        fout << in[row][17] << ' ';
                        fout << endl;
                    }
                    fout.close();
                }
                /*
                 * puw
                 */
                if (!normalize && ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "puw.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    for (unsigned row=0; row<rows; row++) {
                        fout << in[row][6] << ' ';
                        fout << in[row][7] << ' ';
                        fout << in[row][8] << ' ';
                        fout << in[row][9] << ' ';
                        fout << in[row][10] << ' ';
                        fout << in[row][11] << ' ';
                        fout << endl;
                    }
                    fout.close();
                }
                /*
                 * pnw
                 */
                if (normalize && ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "pnw.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    for (unsigned row=0; row<rows; row++) {
                        fout << in[row][0] << ' ';
                        fout << in[row][1] << ' ';
                        fout << in[row][2] << ' ';
                        fout << in[row][3] << ' ';
                        fout << in[row][4] << ' ';
                        fout << in[row][5] << ' ';
                        fout << endl;
                    }
                    fout.close();
                }
                in.operaFITSImageClose();
            } else if (productname.find("i.fits") != string::npos) {
                if (verbose) {
                    cout << "operaApplyTransform: dir=" << directory << endl;
                }
                operaFITSProduct in(productname, READONLY);
                basefilename = basefilename.substr(0, basefilename.find("i.fits"));
                string outfilenamebase = directory + basefilename;
                if (in.getnaxis1() > in.getnaxis2()) {
                    in.rotate90();
                }
                unsigned rows = (unsigned)in.getYDimension();
                instrumentmode_t instrumentmode;
                string mode = in.operaFITSGetHeaderValue("INSTMODE");
                unsigned columns = 0;
                if (mode.find("Polarimetry") != string::npos) {
                    instrumentmode = MODE_POLAR;
                    columns = 6;
                } else if (mode.find("Spectroscopy, star+sky") != string::npos) {
                    instrumentmode = MODE_STAR_PLUS_SKY;
                    columns = 6;
                } else if (mode.find("Spectroscopy, star only") != string::npos) {
                    instrumentmode = MODE_STAR_ONLY;
                    columns = 2;
                } else {
                    throw operaException("operaApplyTransform: "+mode+' ', operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
                }
                string object = in.operaFITSGetHeaderValue("OBJECT");
                /*
                 * iu
                 */
                if (!normalize && !ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "iu.s";
                    if (verbose) {
                        cout << "operaApplyTransform: mode=" << mode << endl;
                        cout << "operaApplyTransform: object='" << object << "'"<< endl;
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    switch (instrumentmode) {
                        case MODE_STAR_ONLY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][9] << ' ';
                                fout << in[row][10] << ' ';
                                fout << in[row][11] << ' ';
                                fout << endl;
                            }
                            break;
                        case MODE_POLAR:
                        case MODE_STAR_PLUS_SKY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][21] << ' ';
                                fout << in[row][22] << ' ';
                                fout << in[row][23] << ' ';
                                fout << in[row][24] << ' ';
                                fout << in[row][25] << ' ';
                                fout << in[row][26] << ' ';
                                fout << in[row][27] << ' ';
                                fout << endl;
                            }
                            break;
                        default:
                            break;
                    }
                    fout.close();
                }
                /*
                 * in
                 */
                if (normalize && !ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "in.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    switch (instrumentmode) {
                        case MODE_STAR_ONLY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][6] << ' ';
                                fout << in[row][7] << ' ';
                                fout << in[row][8] << ' ';
                                fout << endl;
                            }
                            break;
                        case MODE_POLAR:
                        case MODE_STAR_PLUS_SKY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][14] << ' ';
                                fout << in[row][15] << ' ';
                                fout << in[row][16] << ' ';
                                fout << in[row][17] << ' ';
                                fout << in[row][18] << ' ';
                                fout << in[row][19] << ' ';
                                fout << in[row][20] << ' ';
                                fout << endl;
                            }
                            break;
                        default:
                            break;
                    }
                    fout.close();
                }
                /*
                 * iuw
                 */
                if (!normalize && ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "iuw.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    switch (instrumentmode) {
                        case MODE_STAR_ONLY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][3] << ' ';
                                fout << in[row][4] << ' ';
                                fout << in[row][5] << ' ';
                                fout << endl;
                            }
                            break;
                        case MODE_POLAR:
                        case MODE_STAR_PLUS_SKY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][7] << ' ';
                                fout << in[row][8] << ' ';
                                fout << in[row][9] << ' ';
                                fout << in[row][10] << ' ';
                                fout << in[row][11] << ' ';
                                fout << in[row][12] << ' ';
                                fout << endl;
                            }
                            break;
                        default:
                            break;
                    }
                    fout.close();
                }
                /*
                 * inw
                 */
                if (normalize && ApplyWavelengthCorrection) {
                    ofstream fout;
                    string outfilename = outfilenamebase + "inw.s";
                    if (verbose) {
                        cout << "operaApplyTransform: outfilename=" << outfilename << endl;
                    }
                    fout.open(outfilename.c_str());
                    fout << "***Reduced spectrum of '" << object << "'" << endl;
                    fout << rows << ' ' << columns << endl;
                    switch (instrumentmode) {
                        case MODE_STAR_ONLY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][0] << ' ';
                                fout << in[row][1] << ' ';
                                fout << in[row][2] << ' ';
                                fout << endl;
                            }
                            break;
                        case MODE_POLAR:
                        case MODE_STAR_PLUS_SKY:
                            for (unsigned row=0; row<rows; row++) {
                                fout << in[row][0] << ' ';
                                fout << in[row][1] << ' ';
                                fout << in[row][2] << ' ';
                                fout << in[row][3] << ' ';
                                fout << in[row][4] << ' ';
                                fout << in[row][5] << ' ';
                                fout << in[row][6] << ' ';
                                fout << endl;
                            }
                            break;
                        default:
                            break;
                    }
                    fout.close();
                }
            } else {
                throw operaException("operaApplyTransform: "+productname+" ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
    }
	catch (operaException e) {
		cerr << "operaApplyTransform: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaApplyTransform: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	exit(EXIT_SUCCESS);
}

