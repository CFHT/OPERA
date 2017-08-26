/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaPolarTest
 Version: 1.0
 Description: Perform various tests on the operaPolar module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jun/2012
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
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <fitsio.h>
#include <getopt.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaMath.h"   // for PlanckLaw
#include "libraries/operaException.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaFluxVector.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaMuellerMatrix.h"
#include "libraries/operaPolarimetry.h"
#include "libraries/operaStats.h"
#include "libraries/operaSpectrumSimulation.h"

#include "core-espadons/operaPolar.h"

#define NOTPROVIDED -999

/*!
 * \file operaPolarTest.cpp
 * \brief Simulate polarimetry.
 * \details This file holds the implementation for the simulation of polarimetry.
 */

using namespace std;

void GenerateExtractionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display);

/*!
 * \author Andre Venne
 * \brief Simulate polarimetry.
 * \details This module simulates spectrum to test that the polarimetry analysis works for simulated data.
 * \param argc
 * \param argv
 * \note --output=...
 * \note --stokesparameter=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \sa operaPolar, operaSpectrumSimulation
 */
int main(int argc, char *argv[])
{
    srandom(time(NULL));
    
	int opt;
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	string outputfilename, plotoutfilename; 
    
    typedef enum { Difference=1, Ratio, DifferenceWithBeamSwapped} method_t;
    method_t method = Ratio;
    
    stokes_parameter_t StokesParameter = StokesI;
    
    simulation_t SimulationType = BlackBody;
    
    SpectrumVariables_t SpectrumVariables;
    
    SpectrumVariables.StokesParameter = StokesParameter;
    
    SpectrumVariables.ErrorType = GaussianError;
    SpectrumVariables.ErrorPercentage = 0.05;
    
    SpectrumVariables.Temperature = 5000.0; // K
    SpectrumVariables.MinWavelength = 502.0; // nm
    SpectrumVariables.MaxWavelength = 504.0; // nm
    SpectrumVariables.WavelengthIncrement = 0.0005; // nm
    
    SpectrumVariables.StokesParameterFractionPolarizationOfBlackBody = 0.1;
    
    double SunSpectrumConstantValue = 1.0;
    SpectrumVariables.ReflectedJupiterFractionOfSunSpectrum = 0.01;
    
    SpectrumVariables.NumberOfPolarizationPeaks = 1;
    double *AmplitudeOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
    double *SigmaOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
    double *CenterOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
    AmplitudeOfPolarizationVector[0] = 0.03;
    SigmaOfPolarizationVector[0] = 0.1;
    CenterOfPolarizationVector[0] = 503;
    unsigned AmplitudeOfPolarizationIndex = 0;
    unsigned SigmaOfPolarizationIndex = 0;
    unsigned CenterOfPolarizationIndex = 0;
    
    SpectrumVariables.ZeemanSplit = 2.0;  // zeeman split in units of sigma
    SpectrumVariables.linedepth = 0.03;          // line depth with respect to continuum
    SpectrumVariables.linewidth = 0.1;          // line width in nm
    SpectrumVariables.linecenter = 503;          // line center in nm
    
	unsigned namps = 1;
	unsigned NumberOfExposures = 4;
	int ordernumber = NOTPROVIDED;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
    
    bool display = false;
	
    struct option longopts[] = {
		{"output",1, NULL, 'o'},
        {"stokesparameter",1, NULL, 's'},
		{"method",1, NULL, 'm'},
		{"numberofexposures",1, NULL, 'c'},
		{"ordernumber",1, NULL, 'r'},
		{"numberofamplifiers",	1, NULL, 'N'},
		{"simulationtype",	1, NULL, 'a'},
		{"errortype",	1, NULL, 'e'},
		{"errorpercentage",	1, NULL, 'f'},
		{"temperature",	1, NULL, 'T'},
		{"minwavelength",	1, NULL, 'n'},
		{"maxwavelength",	1, NULL, 'x'},
		{"wavelengthincrement",	1, NULL, 'w'},
		{"blackbodyfractionpolarization",	1, NULL, 'b'},
		{"sunspectrum",	1, NULL, 'g'},
		{"jupiterfractionreflection",	1, NULL, 'j'},
		{"numberofpolarizationpeaks",	1, NULL, 'k'},
		{"amplitudeofpolarization",	1, NULL, 'l'},
		{"sigmaofpolarization",	1, NULL, 'z'},
		{"centerofpolarization",	1, NULL, 'q'},
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
        {"scriptfilename",1, NULL, 'S'},
		
		{"plot",		optional_argument, NULL, 'p'},
        {"display",		optional_argument, NULL, 'D'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}
    };
	
	while((opt = getopt_long(argc, argv, "o:s:m:c:r:N:a:e:f:T:n:x:w:b:g:j:k:l:z:q:P:F:S:p::D::v::d::t::h", longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'o':
				outputfilename = optarg;
				break;
            case 's':
				StokesParameter = (stokes_parameter_t)atoi(optarg);
                SpectrumVariables.StokesParameter = StokesParameter;
				break;
			case 'm':
				method = (method_t)atoi(optarg);
				break;
			case 'c':
				NumberOfExposures = atoi(optarg);
				break;  
			case 'r':
				ordernumber = atoi(optarg);
				break;
            case 'N':
				namps = atoi(optarg);
				break;
            case 'a':
				SimulationType = (simulation_t)atoi(optarg);
				break;
            case 'e':
                SpectrumVariables.ErrorType = (error_type_t)atoi(optarg);
				break;
            case 'f':
				SpectrumVariables.ErrorPercentage = atof(optarg);
				break;
            case 'T':
				SpectrumVariables.Temperature = atof(optarg);
				break;
            case 'n':
				SpectrumVariables.MinWavelength = atof(optarg);
				break;
            case 'x':
				SpectrumVariables.MaxWavelength = atof(optarg);
				break;
            case 'w':
				SpectrumVariables.WavelengthIncrement = atof(optarg);
				break;
            case 'b':
				SpectrumVariables.StokesParameterFractionPolarizationOfBlackBody = atof(optarg);
				break;
            case 'g':
				SunSpectrumConstantValue = atof(optarg);
				break;
            case 'j':
				SpectrumVariables.ReflectedJupiterFractionOfSunSpectrum = atof(optarg);
				break;
            case 'k': {
				SpectrumVariables.NumberOfPolarizationPeaks = atoi(optarg);
                
                if (AmplitudeOfPolarizationVector)
                    delete AmplitudeOfPolarizationVector;
                AmplitudeOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
                
                if (SigmaOfPolarizationVector)
                    delete SigmaOfPolarizationVector;
                SigmaOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
                
                if (CenterOfPolarizationVector)
                    delete CenterOfPolarizationVector;
                CenterOfPolarizationVector = new double[SpectrumVariables.NumberOfPolarizationPeaks];
            }
				break;
            case 'l': {
                if (AmplitudeOfPolarizationIndex < SpectrumVariables.NumberOfPolarizationPeaks) {
                    AmplitudeOfPolarizationVector[AmplitudeOfPolarizationIndex] = atof(optarg);
                    AmplitudeOfPolarizationIndex++;
                }
            }
				break;
            case 'z': {
                if (SigmaOfPolarizationIndex < SpectrumVariables.NumberOfPolarizationPeaks) {
                    SigmaOfPolarizationVector[SigmaOfPolarizationIndex] = atof(optarg);
                    SigmaOfPolarizationIndex++;
                }
            }
				break;
            case 'q': {
                if (CenterOfPolarizationIndex < SpectrumVariables.NumberOfPolarizationPeaks) {
                    CenterOfPolarizationVector[CenterOfPolarizationIndex] = atof(optarg);
                    CenterOfPolarizationIndex++;
                }
            }
				break;
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break;
            case 'F':
				datafilename = optarg;
				break;
            case 'S':
				scriptfilename = optarg;
				break;
			case 'p':
				plot = 1;
				break;
            case 'D':
				display = true;
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
        /* Stokes parameter check */
		if (StokesParameter != StokesQ && StokesParameter != StokesU && StokesParameter != StokesV) {
			throw operaException("operaPolar: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        /* Method check */
        if (method != 1 && method != 2 && method != 3) {
            throw operaException("operaPolar: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);    
        }
        /* Number of exposures check */
        if (NumberOfExposures != 2 && NumberOfExposures != 4) {
            throw operaException("operaPolar: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);    
        }
        /* Output file name check */
		if (outputfilename.empty()) {
			throw operaException("operaPolar: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        operaSpectralOrderVector outputorderVector;
        
        ofstream *fdata = NULL ;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }
        
        unsigned length = (unsigned)( (SpectrumVariables.MaxWavelength - SpectrumVariables.MinWavelength) / SpectrumVariables.WavelengthIncrement );
        SpectrumVariables.length = length;
        
        operaFluxVector SunSpectrum(length);
        SunSpectrum = SunSpectrumConstantValue;
        SpectrumVariables.SunSpectrum = &SunSpectrum;
        
        SpectrumVariables.AmplitudeOfPolarizationVector = AmplitudeOfPolarizationVector;
        SpectrumVariables.SigmaOfPolarizationVector = SigmaOfPolarizationVector;
        SpectrumVariables.CenterOfPolarizationVector = CenterOfPolarizationVector;
        
        unsigned minorder = 21;
        unsigned maxorder = 62;
        
        if(ordernumber != NOTPROVIDED) {
            minorder = ordernumber;
            maxorder = ordernumber;
        }
        
        // angleP[9] = {NOTPROVIDED, P1, P2, P3, P4, P5, P6, P7, P8}
        double angleP[9] = {NOTPROVIDED, 0.0, M_PI/8, M_PI/4, 3*M_PI/8, M_PI/2, 5*M_PI/8, 3*M_PI/4, 7*M_PI/8};
        
        if (NumberOfExposures == 2) {
            for (unsigned order=minorder; order<=maxorder; order++) {
                if (verbose)
                    cerr << "operaPolar: Processing order number: " << order << endl;
                
                /*
                 * Populate vectors with the E/A data
                 */
                
                operaFluxVector ESPaDOnSBeamE(length);
                operaFluxVector ESPaDOnSBeamA(length);
                
                operaFluxVector i1E(length);
                operaFluxVector i1A(length);
                
                operaFluxVector i2E(length);
                operaFluxVector i2A(length);
                
                if (StokesParameter == StokesQ) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                }
                else if (StokesParameter == StokesU) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[2];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[2];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[6];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[6];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                }
                else if (StokesParameter == StokesV) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                }
                
                /*
                 * Declare the number of pair of exposures
                 * Create the vectors
                 * The argument ToOne in the declaration of R is there so that in case of R->INFINITY, the result for P/I is 1.0
                 */
                
                double PairOfExposures = (double)NumberOfExposures / 2.0;
                
                operaFluxVector Intensity(length);
                
                operaFluxVector PoverI(length);
                
                operaFluxVector N1(length);
                operaFluxVector N2(length);
                
                operaFluxVector r1(length);
                operaFluxVector r2(length);
                
                operaFluxVector G1(length);         // Difference method
                operaFluxVector G2(length);         // Difference method
                
                operaFluxVector D1(length);         // Difference method
                
                operaFluxVector R1(length);         // Ratio method
                
                operaFluxVector R(length,ToOne);    // Ratio method
                
                /* 
                 * Calculate the Intensity (Stokes I)
                 *  Intensity = (i1E + i1A + i2E + i2A) / NumberOfExposures
                 */
                
                Intensity = (i1E + i1A + i2E + i2A) / (double)NumberOfExposures;
                
                /*
                 * Difference method
                 */
                
                if (method == Difference) {
                    
                    /* 
                     * STEP 1 - calculate the quantity Gn (Eq #12-14 on page 997 of Bagnulo et al. 2009 paper)
                     *  n being the pair of exposures
                     *  G1 = (i1E - i1A)/(i1E + i1A)
                     *  G2 = (i2E - i2A)/(i2E + i2A)
                     */
                    
                    G1 = (i1E - i1A)/(i1E + i1A);
                    G2 = (i2E - i2A)/(i2E + i2A);
                    
                    /*	
                     * STEP 2 - calculate the quantity Dm (Eq #18 on page 997 of Bagnulo et al. 2009 paper)
                     *  m being the pair of exposures
                     * 	D1 = G1 - G2
                     */
                    
                    D1 = G1 - G2;
                    
                    /*	
                     * STEP 3 - calculate the degree of Stokes parameter (Eq #19 on page 997 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = D1 / (2.0 * PairOfExposures)
                     */
                    
                    PoverI = D1 / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 4 - cannot calculate the NULL spectra with only 1 pair of exposures
                     */
                    
                }
                
                /*
                 * Ratio method
                 */
                
                if (method == Ratio) {
                    
                    /* 
                     * STEP 1 - calculate ratio of beams for each exposure (Eq #12 on page 997 of Bagnulo et al. 2009 paper)
                     *  r1 = i1E / i1A
                     *  r2 = i2E / i2A
                     */
                    
                    r1 = i1E / i1A;
                    r2 = i2E / i2A;
                    
                    /*	
                     * STEP 2 - calculate the quantity Rm (Eq #23 on page 998 of Bagnulo et al. 2009 paper)
                     *  m being the pair of exposures
                     * 	R1 = r1 / r2
                     */
                    
                    R1 = r1 / r2;
                    
                    /*	
                     * STEP 3 - calculate the quantity R (Part of Eq #24 on page 998 of Bagnulo et al. 2009 paper)
                     * 	R = Pow( R1 , 1.0/(2.0 * PairOfExposures) )
                     */
                    
                    R = Pow( R1 , 1.0/(2.0 * PairOfExposures) );
                    
                    /*
                     * STEP 4 - calculate the degree of Stokes parameter (Simplification with STEP 2 of Eq #24 on page 998 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = (R - 1.0) / (R + 1.0)
                     */
                    
                    PoverI = (R - 1.0) / (R + 1.0);
                    
                    /*
                     * STEP 5 - cannot calculate the NULL spectra with only 1 pair of exposures
                     */
                    
                }
                
                /*
                 * Difference method with beam swapped
                 */
                
                if (method == DifferenceWithBeamSwapped) {
                    
                    /* 
                     * STEP 1 - calculate the quantity Gn (Eq #12-14 on page 997 of Bagnulo et al. 2009 paper) with beam i1A and i2E, and i3A and i4E swapped
                     *  n being the pair of exposures
                     *  G1 = (i1E - i2E)/(i1E + i2E)
                     *  G2 = (i1A - i2A)/(i1A + i2A)
                     */
                    
                    G1 = (i1E - i2E)/(i1E + i2E);
                    G2 = (i1A - i2A)/(i1A + i2A);
                    
                    /*	
                     * STEP 2 - calculate the quantity Dm (Eq #18 on page 997 of Bagnulo et al. 2009 paper)
                     *  m being the pair of exposures
                     * 	D1 = G1 - G2
                     */
                    
                    D1 = G1 - G2;
                    
                    /*	
                     * STEP 3 - calculate the degree of Stokes parameter (Eq #19 on page 997 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = D1 / (2.0 * PairOfExposures)
                     */
                    
                    PoverI = D1 / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 4 - cannot calculate the NULL spectra with only 1 pair of exposures
                     */
                    
                }
                
                /* Writting to the data file */
                if (fdata != NULL) {
                    fdata->precision(6);
                    *fdata << fixed;
                    *fdata << "# operaPolar: <index> <Degree of Polarization> <Intensity> <i1E> <i1A> <i2E> <i2A> <r1 = i1E / i1A> <r2 = i2E / i2A> <R1 = r1 / r2>\n";
					for (unsigned index = 0 ; index < length ; index++) {
						*fdata << "operaPolar:\t"
                        << index << "\t"
                        << SpectrumVariables.MinWavelength + (SpectrumVariables.WavelengthIncrement)*(double)index << "\t"
						<< PoverI.getflux(index) << "\t"
						<< Intensity.getflux(index) << "\t"
						<< i1E.getflux(index) << "\t"
						<< i1A.getflux(index) << "\t"
						<< i2E.getflux(index) << "\t"
						<< i2A.getflux(index) << "\t"
						<< r1.getflux(index) << "\t"
						<< r2.getflux(index) << "\t"
						<< R1.getflux(index) << "\t"
						<< endl;
					}
				}
                
                /* Writting to the output file */
                outputorderVector.GetSpectralOrder(order)->deletePolarimetry();
                outputorderVector.GetSpectralOrder(order)->createPolarimetry(length);
                outputorderVector.GetSpectralOrder(order)->sethasPolarimetry(true);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setStokesParameter(StokesI, &Intensity);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setDegreeOfPolarization(StokesParameter, &PoverI);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->calculatePolarization();
            }
        }
        else if (NumberOfExposures == 4) {
            for (unsigned order=minorder; order<=maxorder; order++) {
                if (verbose)
                    cerr << "operaPolar: Processing order number: " << order << endl;
                
                /*
                 * Populate vectors with the E/A data
                 * The inputs 3 and 4 are swapped to fit with the Bagnulo algorithm
                 */
                
                operaFluxVector ESPaDOnSBeamE(length);
                operaFluxVector ESPaDOnSBeamA(length);
                
                operaFluxVector i1E(length);
                operaFluxVector i1A(length);
                
                operaFluxVector i2E(length);
                operaFluxVector i2A(length);
                
                operaFluxVector i3E(length);
                operaFluxVector i3A(length);
                
                operaFluxVector i4E(length);
                operaFluxVector i4A(length);
                
                if (StokesParameter == StokesQ) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[3];
                    SpectrumVariables.Rhomb2Angle = angleP[5];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[3];
                    SpectrumVariables.Rhomb2Angle = angleP[5];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i4E.setVector(ESPaDOnSBeamE);
                    i4A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[7];
                    SpectrumVariables.Rhomb2Angle = angleP[7];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[7];
                    SpectrumVariables.Rhomb2Angle = angleP[7];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i3E.setVector(ESPaDOnSBeamE);
                    i3A.setVector(ESPaDOnSBeamA);
                }
                else if (StokesParameter == StokesU) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[2];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[2];
                    SpectrumVariables.Rhomb2Angle = angleP[1];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[6];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[6];
                    SpectrumVariables.Rhomb2Angle = angleP[3];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[4];
                    SpectrumVariables.Rhomb2Angle = angleP[5];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[4];
                    SpectrumVariables.Rhomb2Angle = angleP[5];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i4E.setVector(ESPaDOnSBeamE);
                    i4A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[8];
                    SpectrumVariables.Rhomb2Angle = angleP[7];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[8];
                    SpectrumVariables.Rhomb2Angle = angleP[7];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i3E.setVector(ESPaDOnSBeamE);
                    i3A.setVector(ESPaDOnSBeamA);
                }
                else if (StokesParameter == StokesV) {
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[1];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i1E.setVector(ESPaDOnSBeamE);
                    i1A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[5];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i2E.setVector(ESPaDOnSBeamE);
                    i2A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[3];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[3];
                    SpectrumVariables.Rhomb2Angle = angleP[2];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i4E.setVector(ESPaDOnSBeamE);
                    i4A.setVector(ESPaDOnSBeamA);
                    
                    SpectrumVariables.Beam = 0;
                    SpectrumVariables.Rhomb1Angle = angleP[7];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamE.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamE));
                    
                    SpectrumVariables.Beam = 1;
                    SpectrumVariables.Rhomb1Angle = angleP[7];
                    SpectrumVariables.Rhomb2Angle = angleP[4];
                    
                    ESPaDOnSBeamA.setVector(SimulatedSpectrum(SimulationType, SpectrumVariables, ESPaDOnSBeamA));
                    
                    i3E.setVector(ESPaDOnSBeamE);
                    i3A.setVector(ESPaDOnSBeamA);
                }
                
                /*
                 * Declare the number of pair of exposures
                 * Create the vectors
                 * The argument ToOne in the declaration of R is there so that in case of R->INFINITY, the result for P/I is 1.0
                 */
                
                double PairOfExposures = (double)NumberOfExposures / 2.0;
                
                operaFluxVector Intensity(length);
                
                operaFluxVector PoverI(length);
                
                operaFluxVector N1(length);
                operaFluxVector N2(length);
                
                operaFluxVector r1(length);
                operaFluxVector r2(length);
                operaFluxVector r3(length);
                operaFluxVector r4(length);
                
                operaFluxVector G1(length);         // Difference method
                operaFluxVector G2(length);         // Difference method
                operaFluxVector G3(length);         // Difference method
                operaFluxVector G4(length);         // Difference method
                
                operaFluxVector D1(length);         // Difference method
                operaFluxVector D2(length);         // Difference method
                
                operaFluxVector D1s(length);        // Difference method
                operaFluxVector D2s(length);        // Difference method
                
                operaFluxVector R1(length);         // Ratio method
                operaFluxVector R2(length);         // Ratio method
                
                operaFluxVector R1s(length);        // Ratio method
                operaFluxVector R2s(length);        // Ratio method
                
                operaFluxVector R(length,ToOne);    // Ratio method
                
                operaFluxVector RN1(length);        // Ratio method
                operaFluxVector RN2(length);        // Ratio method
                
                /* 
                 * Calculate the Intensity (Stokes I)
                 *  Intensity = (i1E + i1A + i2E + i2A + i3E + i3A + i4E + i4A) / NumberOfExposures
                 */
                
                Intensity = (i1E + i1A + i2E + i2A + i3E + i3A + i4E + i4A) / (double)NumberOfExposures;
                
                /*
                 * Difference method
                 */
                
                if (method == Difference) {
                    
                    /* 
                     * STEP 1 - calculate the quantity Gn (Eq #12-14 on page 997 of Bagnulo et al. 2009 paper)
                     *  n being the pair of exposures
                     *  G1 = (i1E - i1A)/(i1E + i1A)
                     *  G2 = (i2E - i2A)/(i2E + i2A)
                     *  G3 = (i3E - i3A)/(i3E + i3A)
                     *  G4 = (i4E - i4A)/(i4E + i4A)
                     */
                    
                    G1 = (i1E - i1A)/(i1E + i1A);
                    G2 = (i2E - i2A)/(i2E + i2A);
                    G3 = (i3E - i3A)/(i3E + i3A);
                    G4 = (i4E - i4A)/(i4E + i4A);
                    
                    /*	
                     * STEP 2 - calculate the quantity Dm (Eq #18 on page 997 of Bagnulo et al. 2009 paper) and the quantity Dms with exposure 2 and 4 swapped
                     *  m being the pair of exposures
                     * 	D1 = G1 - G2
                     *  D2 = G3 - G4
                     *  
                     * 	D1s = G1 - G4
                     *  D2s = G3 - G2
                     */
                    
                    D1 = G1 - G2;
                    D2 = G3 - G4;
                    
                    D1s = G1 - G4;
                    D2s = G3 - G2;
                    
                    /*	
                     * STEP 3 - calculate the degree of Stokes parameter (Eq #19 on page 997 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = (D1 + D2) / (2.0 * PairOfExposures)
                     */
                    
                    PoverI = (D1 + D2) / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 4 - calculate the first NULL spectrum (Eq #20 on page 997 of Bagnulo et al. 2009 paper)
                     * 	N1 = (D1 - D2) / (2.0 * PairOfExposures)
                     */
                    
                    N1 = (D1 - D2) / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 5 - calculate the second NULL spectrum (Eq #20 on page 997 of Bagnulo et al. 2009 paper) with exposure 2 and 4 swapped
                     * 	N2 = (D1s - D2s) / (2.0 * PairOfExposures)
                     */
                    
                    N2 = (D1s - D2s) / (2.0 * PairOfExposures);
                }
                
                /*
                 * Ratio method
                 */
                
                if (method == Ratio) {
                    
                    /* 
                     * STEP 1 - calculate ratio of beams for each exposure (Eq #12 on page 997 of Bagnulo et al. 2009 paper)
                     *  r1 = i1E / i1A
                     *  r2 = i2E / i2A
                     *  r3 = i3E / i3A
                     *  r4 = i4E / i4A
                     */
                    
                    r1 = i1E / i1A;
                    r2 = i2E / i2A;
                    r3 = i3E / i3A;
                    r4 = i4E / i4A;
                    
                    /*	
                     * STEP 2 - calculate the quantity Rm (Eq #23 on page 998 of Bagnulo et al. 2009 paper) and the quantity Rms with exposure 2 and 4 swapped
                     *  m being the pair of exposures
                     * 	R1 = r1 / r2
                     * 	R2 = r3 / r4
                     *  
                     * 	R1s = r1 / r4
                     * 	R2s = r3 / r2
                     */
                    
                    R1 = r1 / r2;
                    R2 = r3 / r4;
                    
                    R1s = r1 / r4;
                    R2s = r3 / r2;
                    
                    /*	
                     * STEP 3 - calculate the quantity R (Part of Eq #24 on page 998 of Bagnulo et al. 2009 paper)
                     * 	R = Pow( R1 * R2 , 1.0/(2.0 * PairOfExposures) )
                     */
                    
                    R = Pow( R1 * R2 , 1.0/(2.0 * PairOfExposures) );
                    
                    /*
                     * STEP 4 - calculate the degree of Stokes parameter (Simplification with STEP 2 of Eq #24 on page 998 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = (R - 1.0) / (R + 1.0)
                     */
                    
                    PoverI = (R - 1.0) / (R + 1.0);
                    
                    /*	
                     * STEP 5 - calculate the quantity RN1 (Part of Eq #25-26 on page 998 of Bagnulo et al. 2009 paper)
                     * 	RN1 = Pow( R1 / R2 , 1.0/(2.0 * PairOfExposures) )
                     */
                    
                    RN1 = Pow( R1 / R2 , 1.0/(2.0 * PairOfExposures) );
                    
                    /*
                     * STEP 6 - calculate the first NULL spectrum (Simplification with STEP 4 of Eq #25-26 on page 998 of Bagnulo et al. 2009 paper)
                     * 	N1 = (RN1 - 1.0) / (RN1 + 1.0)
                     */
                    
                    N1 = (RN1 - 1.0) / (RN1 + 1.0);
                    
                    /*	
                     * STEP 7 - calculate the quantity RN2 (Part of Eq #25-26 on page 998 of Bagnulo et al. 2009 paper) with exposure 2 and 4 swapped
                     * 	RN2 = Pow( R1s / R2s , 1.0/(2.0 * PairOfExposures) )
                     */
                    
                    RN2 = Pow( R1s / R2s , 1.0/(2.0 * PairOfExposures) );
                    
                    /*
                     * STEP 8 - calculate the second NULL spectrum (Simplification with STEP 4 of Eq #25-26 on page 998 of Bagnulo et al. 2009 paper) with exposure 2 and 4 swapped
                     * 	N2 = (RN2 - 1.0) / (RN2 + 1.0)
                     */
                    
                    N2 = (RN2 - 1.0) / (RN2 + 1.0);
                }
                
                /*
                 * Difference method with beam swapped
                 */
                
                if (method == DifferenceWithBeamSwapped) {
                    
                    /* 
                     * STEP 1 - calculate the quantity Gn (Eq #12-14 on page 997 of Bagnulo et al. 2009 paper) with beam i1A and i2E, and i3A and i4E swapped
                     *  n being the pair of exposures
                     *  G1 = (i1E - i2E)/(i1E + i2E)
                     *  G2 = (i1A - i2A)/(i1A + i2A)
                     *  G3 = (i3E - i4E)/(i3E + i4E)
                     *  G4 = (i3A - i4A)/(i3A + i4A)
                     */
                    
                    G1 = (i1E - i2E)/(i1E + i2E);
                    G2 = (i1A - i2A)/(i1A + i2A);
                    G3 = (i3E - i4E)/(i3E + i4E);
                    G4 = (i3A - i4A)/(i3A + i4A);
                    
                    /*	
                     * STEP 2 - calculate the quantity Dm (Eq #18 on page 997 of Bagnulo et al. 2009 paper) and the quantity Dms with exposure 2 and 4 swapped
                     *  m being the pair of exposures
                     * 	D1 = G1 - G2
                     *  D2 = G3 - G4
                     *  
                     * 	D1s = G1 - G4
                     *  D2s = G3 - G2
                     */
                    
                    D1 = G1 - G2;
                    D2 = G3 - G4;
                    
                    D1s = G1 - G4;
                    D2s = G3 - G2;
                    
                    /*	
                     * STEP 3 - calculate the degree of Stokes parameter (Eq #19 on page 997 of Bagnulo et al. 2009 paper)
                     * 	aka "the degree of polarization" P
                     * 	P/I = (D1 + D2) / (2.0 * PairOfExposures)
                     */
                    
                    PoverI = (D1 + D2) / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 4 - calculate the first NULL spectrum (Eq #20 on page 997 of Bagnulo et al. 2009 paper)
                     * 	N1 = (D1 - D2) / (2.0 * PairOfExposures)
                     */
                    
                    N1 = (D1 - D2) / (2.0 * PairOfExposures);
                    
                    /*
                     * STEP 5 - calculate the second NULL spectrum (Eq #20 on page 997 of Bagnulo et al. 2009 paper) with exposure 2 and 4 swapped
                     * 	N2 = (D1s - D2s) / (2.0 * PairOfExposures)
                     */
                    
                    N2 = (D1s - D2s) / (2.0 * PairOfExposures);
                }
                
                /* Writting to the data file */
                if (fdata != NULL) {
                    fdata->precision(6);
                    *fdata << fixed;
                    *fdata << "# operaPolar: <index> <Degree of Polarization> <Intensity> <i1E> <i1A> <i2E> <i2A> <i3E> <i3A> <i4E> <i4A> <r1 = i1E / i1A> <r2 = i2E / i2A> <r3 = i3E / i3A> <r4 = i4E / i4A> <R1 = r1 / r2> <R2 = r3 / r4> <First Null Flux> <Second Null Flux>\n";
					for (unsigned index = 0 ; index < length ; index++) {
						*fdata << "operaPolar:\t"
                        << index << "\t"
                        << SpectrumVariables.MinWavelength + (SpectrumVariables.WavelengthIncrement)*(double)index << "\t"                        
						<< PoverI.getflux(index) << "\t"
						<< Intensity.getflux(index) << "\t"
						<< i1E.getflux(index) << "\t"
						<< i1A.getflux(index) << "\t"
						<< i2E.getflux(index) << "\t"
						<< i2A.getflux(index) << "\t"
						<< i3E.getflux(index) << "\t"
						<< i3A.getflux(index) << "\t"
						<< i4E.getflux(index) << "\t"
						<< i4A.getflux(index) << "\t"
						<< r1.getflux(index) << "\t"
						<< r2.getflux(index) << "\t"
						<< r3.getflux(index) << "\t"
						<< r4.getflux(index) << "\t"
						<< R1.getflux(index) << "\t"
						<< R2.getflux(index) << "\t"
						<< N1.getflux(index) << "\t"
						<< N2.getflux(index) << "\t"
						<< endl;
					}
				}
                
                /* Writting to the output file */
                outputorderVector.GetSpectralOrder(order)->deletePolarimetry();
                outputorderVector.GetSpectralOrder(order)->createPolarimetry(length);
                outputorderVector.GetSpectralOrder(order)->sethasPolarimetry(true);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setStokesParameter(StokesI, &Intensity);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setDegreeOfPolarization(StokesParameter, &PoverI);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->calculatePolarization();
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setFirstNullPolarization(StokesParameter, &N1);
                outputorderVector.GetSpectralOrder(order)->getPolarimetry()->setSecondNullPolarization(StokesParameter, &N2);
            }
        }

		outputorderVector.WriteSpectralOrders(outputfilename, Polarimetry);
        
        if (AmplitudeOfPolarizationVector)
            delete AmplitudeOfPolarizationVector;
        if (SigmaOfPolarizationVector)
            delete SigmaOfPolarizationVector;
        if (CenterOfPolarizationVector)
            delete CenterOfPolarizationVector;
	}
	catch (operaException e) {
		cerr << "operaPolar: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaPolar: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*! 
 * \brief Print out the proper program usage syntax.
 * \param modulename Name of the module, which is operaPolar
 */
static void printUsageSyntax(char * modulename) {
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-pDvdth]"
    " --output=<FILE_NAME>"
    " --stokesparameter=<STOKES_PARAMETER>"
    " --method=<UNS_VALUE>"
    " --numberofexposures=<UNS_VALUE>"
    " --ordernumber=<INT_VALUE>"
    " --numberofamplifiers=<UNS_VALUE>"
    " --plotfilename=<FILE_NAME>"
    " --datafilename=<FILE_NAME>"
    " --scriptfilename=<FILE_NAME> \n\n"
	" Example: "+string(modulename)+" -v -p --output=o.txt --stokesparameter=1 --method=2 --numberofexposures=4 --ordernumber=34 --plotfilename=plot.eps --datafilename=data.dat --scriptfilename=script.gnu -v -p \n\n"
    "  -o, --output=<FILE_NAME>,  Output file name \n"
    "  -s, --stokesparameter=<UNS_VALUE>, Which Stokes parameter the module is calculating \n"
    "                              Available options are = 0, 1, 2 or 3, where: \n"
    "                              0. Stokes I (default)\n"
    "                              1. Stokes Q \n"
    "                              2. Stokes U \n"
    "                              3. Stokes V \n"
    "  -m, --method=<UNS_VALUE>, Method for calculation of polarisation \n"
    "                              Available options are = 1, 2, where: \n"
    "                              1. Difference \n"
    "                              2. Ratio (default) \n"
    "  -c, --numberofexposures=<UNS_VALUE>, Number of input file to use \n"
    "                              Available options are = 2, 4 \n"
    "                              2. Use the first 2 input files \n"
    "                              4. Use all 4 input files (default) \n"
    "  -r, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all) \n"
    "  -N, --numberofamplifiers=<UNS_VALUE>, ... \n"
    "  -P, --plotfilename=<FILE_NAME>, Output plot eps file name \n"
    "  -F, --datafilename=<FILE_NAME>, Output data file name \n"
    "  -S, --scriptfilename=<FILE_NAME>, Output gnuplot script file name \n\n"
    "  -p, --plot,  Turn on plotting \n"
    "  -D, --display,  Turn on display of plotting \n"
	"  -v, --verbose,  Turn on message sending \n"
	"  -d, --debug,  Turn on debug messages \n"
	"  -t, --trace,  Turn on trace messages \n"
    "  -h, --help,  display help message \n";
}

void GenerateExtractionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName,"w");
    
    fprintf(fgnu,"reset\n");
    fprintf(fgnu,"set xlabel 'Distance (pixels)'\n");
    fprintf(fgnu,"set ylabel 'Degree of polarization'\n");
    
    fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
    fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
    
    fprintf(fgnu,"plot \"%s\" u 2:3 t \"Stokes Q\"  w p\n",dataFileName);
    
    if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
		fprintf(fgnu,"replot\n");
		fclose(fgnu);
		systemf("gnuplot -persist %s",gnuScriptFileName);
    } else {
		fclose(fgnu);
	}
}
