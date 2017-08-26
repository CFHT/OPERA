/********************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaPolar
 Version: 1.0
 Description: This module calculates the statistics.
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaFluxVector.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaPolarimetry.h"
#include "libraries/operaStats.h"

#define NOTPROVIDED -999

/*!
 * \file operaStatistics.cpp
 * \brief Calculate statistics.
 * \details This file holds the implementation for the calculation of statistics.
 */

using namespace std;

/*!
 * \brief Definition of the value of each PolarizationType.
 */
typedef enum { Polarization=0, DegreeOfPolarization, FirstNullPolarization, SecondNullPolarization} polarization_type_t;

static void printUsageSyntax(char *modulename);

void GenerateExtractionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display);

/*!
 * \author Andre Venne
 * \brief Calculate statistics.
 * \param argc
 * \param argv
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \sa class operaFluxVector, class operaStokesVector, class operaPolarimetry
 * \sa class operaSpectralElements, class operaSpectralOrder, class operaSpectralOrderVector
 * \ingroup tools
 */
int main(int argc, char *argv[])
{
	int opt;
	
	int plot=0, display=0, verbose=0, debug=0, trace=0;
	
	string inputfilename;
	string outputfilename;
    
    double UpperLimit = 1.0;
    double LowerLimit = -1.0;
    
    unsigned NumberOfSegments = 3;
    double IncrementOfDivision = 0.0001;
    
    operaSpectralOrder_t format = Polarimetry;
    polarization_type_t PolarizationType = SecondNullPolarization;
    
    stokes_parameter_t StokesParameter = StokesI;
    
	int ordernumber = NOTPROVIDED;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
    	
	struct option longopts[] = {
		{"input",1, NULL, 'i'},
		{"output",1, NULL, 'o'},
		{"ordernumber",1, NULL, 'r'},
		{"upperlimit",1, NULL, 'u'},
		{"lowerlimit",1, NULL, 'l'},
		{"numberofsegments",1, NULL, 's'},
		{"incrementofdivision",1, NULL, 'n'},
		{"format",1, NULL, 'f'},
		{"polarizationtype",1, NULL, 'm'},
        {"stokesparameter",1, NULL, 'q'},
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
	
	while((opt = getopt_long(argc, argv, "i:o:r:u:l:s:n:f:m:q:P:F:S:p::D::v::d::t::h", longopts, NULL))  != -1)
	{
		switch(opt)
		{
			case 'i':
				inputfilename = optarg;
				break;
			case 'o':
				outputfilename = optarg;
				break;
			case 'r':
				ordernumber = atoi(optarg);
				break;
			case 'u':
				UpperLimit = atof(optarg);
				break;
			case 'l':
				LowerLimit = atof(optarg);
				break;
			case 's':
				NumberOfSegments = atoi(optarg);
				break;
			case 'n':
				IncrementOfDivision = atof(optarg);
				break;
			case 'f':
				format = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'm':
				PolarizationType = (polarization_type_t)atoi(optarg);
				break;
            case 'q':
				StokesParameter = (stokes_parameter_t)atoi(optarg);
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
				display = 1;
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
        // we need an inputs
        if (inputfilename.empty()) {
            throw operaException("operaStatistics: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
        }
        // we need an output
		if (outputfilename.empty()) {
			throw operaException("operaStatistics: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
        ofstream *outputfile = new ofstream();
        outputfile->open(outputfilename.c_str());
        outputfile->precision(6);
        *outputfile << fixed;
        *outputfile << "# operaStatistics: <order number> <Stokes parameter> <polarization type> <index> <MedSigma> <Sigma>\n";
        
        ofstream *fdata = NULL ;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
            fdata->precision(6);
            *fdata << fixed;
            *fdata << "# operaStatistics: <center segment> <number of points segment 1> <number of points segment 2> ... <number of points segment N>\n";
        }
        
        switch (format) {
            case Polarimetry: {
                
                /* Create the spectral order vector based on inputs */
                operaSpectralOrderVector spectralOrders;
                operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputfilename);
                
                unsigned minorder = spectralOrders.getMinorder();
                unsigned maxorder = spectralOrders.getMaxorder();
                
                if(ordernumber != NOTPROVIDED) {
                    minorder = ordernumber;
                    maxorder = ordernumber;
                }
                
                for (unsigned order=minorder; order<=maxorder; order++) {
                    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                    operaPolarimetry *polarimetry = spectralOrder->getPolarimetry();
                    
                    unsigned length = polarimetry->getLength();
                    
                    if (verbose)
                        cout << "operaStatistics: Processing order number: " << order << endl;
                    
                    unsigned UnsignedStokesParameterLowerLimit, UnsignedStokesParameterUpperLimit;
                    
                    if (StokesParameter != StokesI) {
                        UnsignedStokesParameterLowerLimit = (unsigned)StokesParameter;
                        UnsignedStokesParameterUpperLimit = (unsigned)StokesParameter;
                    }
                    else {
                        UnsignedStokesParameterLowerLimit = (unsigned)StokesQ;
                        UnsignedStokesParameterUpperLimit = (unsigned)StokesV;
                    }
                    
                    for (unsigned UnsignedStokesParameter = UnsignedStokesParameterLowerLimit ; UnsignedStokesParameter <= UnsignedStokesParameterUpperLimit ; UnsignedStokesParameter++) {
                        StokesParameter = (stokes_parameter_t)UnsignedStokesParameter;
                        if (StokesParameter == StokesQ)
                            if(!polarimetry->getHasStokesQ())
                                continue;
                        if (StokesParameter == StokesU)
                            if(!polarimetry->getHasStokesU())
                                continue;
                        if (StokesParameter == StokesV)
                            if(!polarimetry->getHasStokesV())
                                continue;
                        
                        if (verbose)
                            cout << "operaStatistics: Processing Stokes parameter: " << StokesParameter << endl;
                        
                        /*
                         * length and NumberOfSegments are two integer, which means that the division will drop any remaining fraction.
                         * This way, NumberOfSegments times LengthOfSegments will always be smaller or equal to length.
                         *      NumberOfSegments * LengthOfSegments <= length
                         * If it's smaller, the remaining segment will be shorter than the others and will have no statistical meaning.
                         * Thus it's neglected.
                         */
                        unsigned LengthOfSegments = length / NumberOfSegments;
                        
                        unsigned LengthOfDivision = (unsigned)ceil((UpperLimit - LowerLimit)/IncrementOfDivision);
                        
                        operaFluxVector MedSigma(NumberOfSegments);
                        
                        operaFluxVector **NumberOfPointsInEachDivision = new operaFluxVector*[NumberOfSegments];
                        for (unsigned index = 0 ; index < NumberOfSegments; index++) {
                            NumberOfPointsInEachDivision[index] = new operaFluxVector(LengthOfDivision);
                            *NumberOfPointsInEachDivision[index] = 0.0;
                        }
                        
                        operaFluxVector Sigma(NumberOfSegments);
                        
                        operaFluxVector Segment(LengthOfSegments);
                        double median;
                        
                        operaFluxVector FluxVector(length);
                        
                        switch (PolarizationType) {
                            case Polarization:
                                FluxVector = polarimetry->getStokesParameter(StokesParameter);        // should not be used for now
                                break;
                            case DegreeOfPolarization:
                                FluxVector = polarimetry->getDegreeOfPolarization(StokesParameter);
                                break;
                            case FirstNullPolarization:
                                FluxVector = polarimetry->getFirstNullPolarization(StokesParameter);
                                break;
                            case SecondNullPolarization:
                                FluxVector = polarimetry->getSecondNullPolarization(StokesParameter);
                                break;
                            default:
                                break;
                        }
                        
                        if (verbose)
                            cout << "operaStatistics: Polarization Type: " << PolarizationType << endl;
                        
                        unsigned indexLength = 0;
                        for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                            
                            if (verbose)
                                cout << "operaStatistics: Processing segment: " << indexSegment << endl;
                            
                            Segment = 0.0;
                            
                            for (unsigned index = 0 ; index < LengthOfSegments ; index++) {
                                Segment.setflux(FluxVector.getflux(indexLength), index);
                                indexLength++;
                            }
                            
                            median = (double)operaArrayMedian(LengthOfSegments, (float*)Segment.getfluxpointer());     // is that a good idea ?
                            MedSigma.setflux((double)operaArrayMedianSigma(LengthOfSegments, (float*)Segment.getfluxpointer(), median), indexSegment);
                            
                            for (unsigned index = 0 ; index < LengthOfSegments ; index++) {
                                double LowerValueOfDivision = LowerLimit;
                                for (unsigned indexOfDivision = 0 ; indexOfDivision < LengthOfDivision ; indexOfDivision++) {
                                    double SegmentValue = Segment.getflux(index);
                                    if (SegmentValue >= LowerValueOfDivision && SegmentValue < (LowerValueOfDivision + IncrementOfDivision)) {
                                        double temp = NumberOfPointsInEachDivision[indexSegment]->getflux(indexOfDivision);
                                        temp += 1.0;
                                        NumberOfPointsInEachDivision[indexSegment]->setflux(temp, indexOfDivision);
                                    }
                                    LowerValueOfDivision += IncrementOfDivision;
                                }
                            }
                            
                            // doesn't work
                            median = (double)operaArrayMedian(LengthOfDivision, (float*)NumberOfPointsInEachDivision[indexSegment]->getfluxpointer());     // is that a good idea ?
                            Sigma.setflux((double)operaArrayMedianSigma(LengthOfDivision, (float*)NumberOfPointsInEachDivision[indexSegment]->getfluxpointer(), median), indexSegment);
                        }
                        
                        for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                            *outputfile << "operaStatistics:\t"
                            << order << "\t"
                            << StokesParameter << "\t"
                            << PolarizationType << "\t"
                            << indexSegment << "\t"
                            << MedSigma.getflux(indexSegment) << "\t"
                            << Sigma.getflux(indexSegment) << "\t"
                            << endl;
                        }
                        
                        if (fdata != NULL) {
                            double LowerValueOfDivision = LowerLimit;
                            for (unsigned indexOfDivision = 0 ; indexOfDivision < LengthOfDivision ; indexOfDivision++) {
                                *fdata << "operaStatistics:\t"
                                << LengthOfSegments << "\t"
                                << LowerValueOfDivision + IncrementOfDivision/2 << "\t";
                                for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                                    *fdata << NumberOfPointsInEachDivision[indexSegment]->getflux(indexOfDivision) << "\t";
                                }
                                *fdata << endl;
                                LowerValueOfDivision += IncrementOfDivision;
                            }
                        }
                        
                        if (NumberOfPointsInEachDivision) {
                            for (unsigned index = 0 ; index < NumberOfSegments; index++) {
                                if (NumberOfPointsInEachDivision[index])
                                    delete NumberOfPointsInEachDivision[index];
                                NumberOfPointsInEachDivision[index] = NULL;
                            }
                            delete NumberOfPointsInEachDivision;
                        }
                        NumberOfPointsInEachDivision = NULL;
                    }
                }
            }
                break;
            /*case LibreEspritSpectrum: {
                
                unsigned minorder = 23;
                unsigned maxorder = 60;
                
                unsigned count;
                
                // Create the spectral order vector based on inputs
                operaSpectralOrderVector orderVector;
                orderVector.readOrdersFromLibreEspritPolarimetry(inputfilename, StokesParameter, count, maxorder);
                
                if(ordernumber != NOTPROVIDED) {
                    minorder = ordernumber;
                    maxorder = ordernumber;
                }
                
                if (StokesParameter == StokesI) {
                    cout << "operaStatistics: Please specify a Stokes parameter other than Stokes I (value = 0)." << endl
                    << "Options are: Stokes Q (value = 1), Stokes U (value = 2), Stokes V (value = 3)" << endl;
                }
                
                for (unsigned order=minorder; order<=maxorder; order++) {
                    if (StokesParameter == StokesI)
                        continue;
                    
                    operaSpectralOrder *spectralOrder = orderVector.GetSpectralOrder(order);
                    operaPolarimetry *polarimetry = spectralOrder->getPolarimetry();
                    
                    unsigned length = polarimetry->getLength();
                    
                    if (verbose)
                        cout << "operaStatistics: Processing order number: " << order << endl;
                    
                    if (StokesParameter == StokesQ)
                        if(!polarimetry->getHasStokesQ())
                            continue;
                    if (StokesParameter == StokesU)
                        if(!polarimetry->getHasStokesU())
                            continue;
                    if (StokesParameter == StokesV)
                        if(!polarimetry->getHasStokesV())
                            continue;
                    
                    if (verbose)
                        cout << "operaStatistics: Processing Stokes parameter: " << StokesParameter << endl;
                    
                    // length and NumberOfSegments are two integer, which means that the division will drop any remaining fraction.
                    // This way, NumberOfSegments times LengthOfSegments will always be smaller or equal to length.
                    //      NumberOfSegments * LengthOfSegments <= length
                    // If it's smaller, the remaining segment will be shorter than the others and will have no statistical meaning.
                    // Thus it's neglected.
                    unsigned LengthOfSegments = length / NumberOfSegments;
                    
                    unsigned LengthOfDivision = (unsigned)ceil((UpperLimit - LowerLimit)/IncrementOfDivision);
                    
                    operaFluxVector MedSigma(NumberOfSegments);
                    
                    operaFluxVector **NumberOfPointsInEachDivision = new operaFluxVector*[NumberOfSegments];
                    for (unsigned index = 0 ; index < NumberOfSegments; index++) {
                        NumberOfPointsInEachDivision[index] = new operaFluxVector(LengthOfDivision);
                        *NumberOfPointsInEachDivision[index] = 0.0;
                    }
                    
                    operaFluxVector Sigma(NumberOfSegments);
                    
                    operaFluxVector Segment(LengthOfSegments);
                    double median;
                    
                    operaFluxVector FluxVector(length);
                    
                    switch (PolarizationType) {
                        case Polarization:
                            FluxVector.setVector(*polarimetry->getStokesParameter(StokesParameter));        // should not be used for now
                            break;
                        case DegreeOfPolarization:
                            FluxVector.setVector(*polarimetry->getDegreeOfPolarization(StokesParameter));
                            break;
                        case FirstNullPolarization:
                            FluxVector.setVector(*polarimetry->getFirstNullPolarization(StokesParameter));
                            break;
                        case SecondNullPolarization:
                            FluxVector.setVector(*polarimetry->getSecondNullPolarization(StokesParameter));
                            break;
                        default:
                            break;
                    }
                    
                    if (verbose)
                        cout << "operaStatistics: Polarization Type: " << PolarizationType << endl;
                    
                    unsigned indexLength = 0;
                    for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                        
                        if (verbose)
                            cout << "operaStatistics: Processing segment: " << indexSegment << endl;
                        
                        Segment = 0.0;
                        
                        for (unsigned index = 0 ; index < LengthOfSegments ; index++) {
                            Segment.setflux(FluxVector.getflux(indexLength), index);
                            indexLength++;
                        }
                        
                        median = (double)operaArrayMedian(LengthOfSegments, (float*)Segment.getfluxpointer());     // is that a good idea ?
                        MedSigma.setflux((double)operaArrayMedianSigma(LengthOfSegments, (float*)Segment.getfluxpointer(), median), indexSegment);
                        
                        for (unsigned index = 0 ; index < LengthOfSegments ; index++) {
                            double LowerValueOfDivision = LowerLimit;
                            for (unsigned indexOfDivision = 0 ; indexOfDivision < LengthOfDivision ; indexOfDivision++) {
                                double SegmentValue = Segment.getflux(index);
                                if (SegmentValue >= LowerValueOfDivision && SegmentValue < (LowerValueOfDivision + IncrementOfDivision)) {
                                    double temp = NumberOfPointsInEachDivision[indexSegment]->getflux(indexOfDivision);
                                    temp += 1.0;
                                    NumberOfPointsInEachDivision[indexSegment]->setflux(temp, indexOfDivision);
                                }
                                LowerValueOfDivision += IncrementOfDivision;
                            }
                        }
                        
                        // doesn't work
                        median = (double)operaArrayMedian(LengthOfDivision, (float*)NumberOfPointsInEachDivision[indexSegment]->getfluxpointer());     // is that a good idea ?
                        Sigma.setflux((double)operaArrayMedianSigma(LengthOfDivision, (float*)NumberOfPointsInEachDivision[indexSegment]->getfluxpointer(), median), indexSegment);
                    }
                    
                    for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                        *outputfile << "operaStatistics:\t"
                        << order << "\t"
                        << StokesParameter << "\t"
                        << PolarizationType << "\t"
                        << indexSegment << "\t"
                        << MedSigma.getflux(indexSegment) << "\t"
                        << Sigma.getflux(indexSegment) << "\t"
                        << endl;
                    }
                    
                    if (fdata != NULL) {
                        double LowerValueOfDivision = LowerLimit;
                        for (unsigned indexOfDivision = 0 ; indexOfDivision < LengthOfDivision ; indexOfDivision++) {
                            *fdata << "operaStatistics:\t"
                            << LengthOfSegments << "\t"
                            << LowerValueOfDivision + IncrementOfDivision/2 << "\t";
                            for (unsigned indexSegment = 0 ; indexSegment < NumberOfSegments ; indexSegment++) {
                                *fdata << NumberOfPointsInEachDivision[indexSegment]->getflux(indexOfDivision) << "\t";
                            }
                            *fdata << endl;
                            LowerValueOfDivision += IncrementOfDivision;
                        }
                    }
                    
                    if (NumberOfPointsInEachDivision) {
                        for (unsigned index = 0 ; index < NumberOfSegments; index++) {
                            if (NumberOfPointsInEachDivision[index])
                                delete NumberOfPointsInEachDivision[index];
                            NumberOfPointsInEachDivision[index] = NULL;
                        }
                        delete NumberOfPointsInEachDivision;
                    }
                    NumberOfPointsInEachDivision = NULL;
                }
            }
                break;*/
            default:
                break;
        }
        
        outputfile->close();
        
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateExtractionPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(),(bool)display);
            }
        }
	}
	catch (operaException e) {
		cerr << "operaStatistics: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaStatistics: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
} 

/*! 
 * \brief Print out the proper program usage syntax.
 * \param modulename Name of the module, which is operaPolar
 */
static void printUsageSyntax(char *modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-pDvdth]"
    " --input=<FILE_NAME>"
    " --output=<FILE_NAME>"
    " --ordernumber=<INT_VALUE>"
    " --plotfilename=<FILE_NAME>"
    " --datafilename=<FILE_NAME>"
    " --scriptfilename=<FILE_NAME> \n\n"
	" Example: "+string(modulename)+" --output=o.txt --input=input --ordernumber=34 --plotfilename=plot.eps --datafilename=data.dat --scriptfilename=script.gnu -v -p \n\n"
    "  -i, --input=<FILE_NAME>,  Polar input file name \n"
    "  -o, --output=<FILE_NAME>,  Output file name \n"
    "  -r, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all) \n"
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
