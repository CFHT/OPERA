/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaPolarimetryTest
 Version: 1.0
 Description: Perform various tests on the operaGeometricShapes class.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
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

#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaPolarimetry.h"
#include "libraries/operaMuellerMatrix.h"
#include "libraries/operaStokesVector.h"

static void printUsageSyntax();

/*! \file operaPolarimetryTest.cpp */
/*! \package operaPolarimetryTest - Perform various tests on the operaPolarimetry class.*/

using namespace std;

/*! 
 * operaPolarimetryTest
 * \author Andre Venne
 * \brief Perform various tests on the operaPolarimetry, operaMuellerMatrix and operaStokesVector class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
    int opt;
	int verbose = 0;
	
	struct option longopts[] = {       
		{"verbose",		optional_argument, NULL, 'v'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "v::h", longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'v':
				verbose = 1;
				break;
            case 'h':
				printUsageSyntax();
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax();
				exit(EXIT_SUCCESS);
				break;
		}
	}
    try {
        
        unsigned length = 1;
        
        double PRhomb1 = 0;
        double PRhomb3 = 0;
        
        double Determinant;
        
        double initializationMatrix[4][4] = { {2, 3, 4, 5}, {0, -1, 2, 1}, {8, 0, 2, 4}, {0, 3, -6, 0} };
        operaMuellerMatrix testMatrix;
        testMatrix.setMuellerMatrix(initializationMatrix);

        operaMuellerMatrix FresnelRhomb1;   //Half-wave
        operaMuellerMatrix FresnelRhomb2;   //Quarter-wave
        operaMuellerMatrix FresnelRhomb3;   //Half-wave
        operaMuellerMatrix FresnelRhombs;
        
        operaMuellerMatrix WollastonO;      //Ordinary axis | Perpendicular to the optical axis
        operaMuellerMatrix WollastonE;      //Extraordinary axis | Parallel to the optical axis
        
        operaMuellerMatrix PolarimetryModuleO;
        operaMuellerMatrix PolarimetryModuleE;
        
        operaMuellerMatrix CofactorMatrix;
        operaMuellerMatrix AdjointMatrix;
        operaMuellerMatrix InverseMatrix;

        FresnelRhomb1.createRotatedRetarder(M_PI,0,PRhomb1,0);
        FresnelRhomb2.createRotatedRetarder(M_PI/2,0,0,0);
        FresnelRhomb3.createRotatedRetarder(M_PI,0,PRhomb3,0);
        
        WollastonO.createRotatedPolarizer(1,0,0,0,0,0);
        WollastonE.createRotatedPolarizer(1,0,M_PI/2,0,0,0);
        
        FresnelRhombs = FresnelRhomb3 * FresnelRhomb2 * FresnelRhomb1;
        
        PolarimetryModuleO = WollastonO * FresnelRhombs;
        PolarimetryModuleE = WollastonE * FresnelRhombs;
        
        Determinant = testMatrix.matrixDeterminant().value;
        CofactorMatrix = testMatrix.matrixCofactor();
        AdjointMatrix = testMatrix.matrixAdjoint();
        InverseMatrix = testMatrix.matrixInverse();
        
        if (verbose) {
            cout << "FresnelRhomb1" << endl;
            FresnelRhomb1.printMuellerMatrix();
            
            cout << "FresnelRhomb2" << endl;
            FresnelRhomb2.printMuellerMatrix();
            
            cout << "FresnelRhomb3" << endl;
            FresnelRhomb3.printMuellerMatrix();
            
            cout << "FresnelRhombs" << endl;
            FresnelRhombs.printMuellerMatrix();
            
            cout << "PolarimetryModuleO" << endl;
            PolarimetryModuleO.printMuellerMatrix();
            
            cout << "PolarimetryModuleE" << endl;
            PolarimetryModuleE.printMuellerMatrix();
            
            cout << "Test Matrix" << endl;
            testMatrix.printMuellerMatrix();
            
            cout << "Matrix Determinant = " << Determinant << endl;
            
            cout << "Cofactor Matrix" << endl;
            CofactorMatrix.printMuellerMatrix();
            
            cout << "Adjoint Matrix" << endl;
            AdjointMatrix.printMuellerMatrix();
            
            cout << "Inverse Matrix" << endl;
            InverseMatrix.printMuellerMatrix();
        }
        
        operaStokesVector *testStokes1 = new operaStokesVector(length);
        operaStokesVector *testStokes2 = new operaStokesVector(length);
        
        testStokes1->getStokesParameter(StokesI)->setflux(1.0,0);

        *testStokes2 = FresnelRhomb1 * *testStokes1;
        
        if (verbose) {
            cout << "__________________" << endl;
            
            for (unsigned StokesIndex=0; StokesIndex<4; StokesIndex++) {
                cout << testStokes1->getStokesParameterFlux((stokes_parameter_t)StokesIndex, 0) << " " << testStokes1->getStokesParameterVariance((stokes_parameter_t)StokesIndex, 0) << endl;
            }
            
            cout << "__________________" << endl;
            
            for (unsigned StokesIndex=0; StokesIndex<4; StokesIndex++) {
                cout << testStokes2->getStokesParameterFlux((stokes_parameter_t)StokesIndex, 0) << " " << testStokes2->getStokesParameterVariance((stokes_parameter_t)StokesIndex, 0) << endl;
            }
        }
        
        operaPolarimetry *testPolar = new operaPolarimetry(length);
        
        if (verbose)
            cout << "number of data points = " << testPolar->getLength() << endl;
        
        testPolar->setStokesParameter(StokesI,1.234, 0.2, 0);
        testPolar->setStokesParameter(StokesI,5.678, 0.4, 1);
        
        if (verbose) {
            for (unsigned i=0; i<length; i++) {
                cout << "value of Stokes I index " << i << " = " << testPolar->getStokesParameter(StokesI)->getflux(i) << endl;
                cout << "variance of Stokes I index " << i << " = " << testPolar->getStokesParameter(StokesI)->getvariance(i) << endl;
            }
        }
        
        if (testStokes1)
            delete testStokes1;
        if (testStokes2)
            delete testStokes2;
        if (testPolar)
            delete testPolar;
	}
	catch (operaException e) {
		cerr << "operaPolarimetryTest: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaPolarimetryTest: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cerr << " Usage: operaPolarimetryTest -[vh]\n";
}	
