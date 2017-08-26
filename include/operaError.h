#ifndef OPERAERROR_H
#define OPERAERROR_H

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaError
 Version: 1.0
 Description: Error definitions  
 to start up with an OPERA module.
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
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

/*! \brief Common error definitions. */
/*! \file operaError.h */
/*! \ingroup libraries */

#define MAXERRSTRINGSIZE 4096
typedef int operaErrorCode;

/*
 * errno uses 1-131
 * do not use those ranges
 */

/*
 * Standard opera error codes.
 */

#define		operaErrorCodeOK  0

#define		operaErrorZeroLength 32750
#define		operaErrorDivideByZeroError 32751
#define		operaErrorCodeBadInstrumentModeError 32752
#define		operaErrorCodeFileDoesNotExistError 32753
#define		operaErrorExtensionHasNotBeenRead 32754
#define		operaErrorInvalidInput 32755
#define		operaErrorNoMemory 32756
#define		operaErrorCodeNoMemory 32757
#define		operaErrorInvalidParameter 32758
#define		operaErrorIsInf 32759
#define		operaErrorIsNaN 32760
#define		operaErrorHeaderProblem 32761
#define		operaErrorNoBadPixelMask 32762
#define		operaErrorPickOutofRange 32763
#define		operaErrorNoInput 32764
#define		operaErrorNoOutput 32765
#define		operaErrorCodeNOTIMPLEMENTED 32766
#define		operaErrorCodeNULL 32767
#define		operaErrorCodeNULLString 32768
#define		operaErrorInputValuesDisagree 32770
#define		operaErrorNegativeValues 32771
#define		operaErrorIndexOutOfRange 32772

/* 
 * opera libraries.
 */

#define operaErrorCodeIncorrectFileType 700
#define	operaErrorCodeEnvironmentnotset 701
#define	operaErrorCodeDatatypeNotSupported 702
#define	operaErrorCodeNoFilename 703
#define operaErrorCodeTileError 704
#define operaErrorCodeBracketingError 705
#define operaErrorCodeChangeREADONLYError 706
#define operaErrorCodeNULLfptr 707
#define operaErrorInvalidGuideWindowIndex 708
#define operaErrorLengthMismatch 709
#define operaErrorNotLazy 710
#define operaErrorDifferingLaziness 711
#define operaErrorExtensionOutOfRange 712
#define operaErrorSliceOutOfRange 713

/*
 * matrix
 */
#define MatrixNotSquare 730
#define MatrixInvalidDimensions 731
#define MatrixZeroDeterminant 732

#define operaSplineInterpolationAxesEqual 740

/*
 * operaReductionSet 800-820
 */
#define		operaErrorReductionSetEtypeNotDefined 800
#define		operaErrorReductionSetQualiKeyNotDefined 801
#define		operaErrorReductionSetQualiValNotDefined 802
#define		operaErrorReductionSetObstypeKeyNotDefined 803
#define		operaErrorReductionSetEtypeFailed 804
#define		operaErrorReductionSetInputNotFound 805

/*
 * modules 900 - 1000
 */
#define		operaErrorGeometryNoOrdersFound 900
#define		operaErrorGeometryBadFit 901
#define		operaErrorNoAnchor 902

/*
 * Spectroscopic Libraries 1001 - 1120
 */

#define operaErrorInstrumentProfileImproperCoordinateRequest 1001
#define operaErrorInstrumentProfileImproperDimensions 1002
#define operaErrorNoSpectralLineDetected 1003
#define operaErrorHasNoSpectralElements 1004
#define operaErrorHasNoInstrumentProfile 1005
#define operaErrorHasNoXCorrelation 1006
#define operaErrorHasNoGeometry 1007
#define operaErrorHasNoExtractionAperture 1008
/*
 * cfitsio uses 101-159, 201-264, 301-264, 401-436, 501-599
 * do not use those ranges
 */

// defined by libraries that just want the error codes
// and not these routines
#include <fitsio.h>

#ifndef OPERAERRRORCODESONLY

#ifdef __cplusplus
#include <iostream>
#include <string>

using namespace std;

/*! 
 * void operaStrError(const operaErrorCode errcode)
 * \brief {This function generates and reports a C-style string, 
 * containing an error message derived from the error code passed in with errcode.}
 * \param errcode is an operaErrorCode
 * \note non-reentrant
 * \return char *
 */
static string operaErrorString = "";
static string operaStrError(const operaErrorCode errcode) {
	if (errcode <= 100) {
		operaErrorString = strerror(errcode);
	} else if (errcode < 600) {
		char fitserrbuff[80];
		fits_get_errstatus(errcode, fitserrbuff);
		operaErrorString = string(fitserrbuff);
	} else switch (errcode) {
		case 0:
			operaErrorString = string("");
			break;

		/* Generic */
		case operaErrorCodeNULL:
			operaErrorString = string("NULL value");
			break;
		case operaErrorCodeNOTIMPLEMENTED:
			operaErrorString = string("not implemented");
			break;
		case operaErrorCodeNULLString:
			operaErrorString = string("NULL string");
			break;
		case operaErrorInputValuesDisagree:
			operaErrorString = string("Two or more input values not in agreement");
			break;
		case operaErrorNegativeValues:
			operaErrorString = string("Values must be positive");
			break;
		case operaErrorIndexOutOfRange:
			operaErrorString = string("Index out of range");
			break;
		case operaErrorPickOutofRange:
			operaErrorString = string("pick out of range");
			break;
		case operaErrorNoInput:
			operaErrorString = string("no inputs specified");
			break;
		case operaErrorNoOutput:
			operaErrorString = string("no output specified");
			break;
		case operaErrorHeaderProblem:
			operaErrorString = string("header problem");
			break;
		case operaErrorIsNaN:
			operaErrorString = string("Invalid floating point number (NaN)");
			break;
		case operaErrorIsInf:
			operaErrorString = string("Invalid floating point number (Inf)");
			break;
		case operaErrorNoMemory:
			operaErrorString = string("No memory");
			break;
		case operaErrorInvalidInput:
			operaErrorString = string("invalid input");
			break;
		case operaErrorExtensionHasNotBeenRead:
			operaErrorString = string("Extension Has Not Been Read");
			break;
		case operaErrorCodeFileDoesNotExistError:
			operaErrorString = string("File does not exist");
			break;
		case operaErrorCodeBadInstrumentModeError:
			operaErrorString = string("Bad instrument mode");
			break;
		case operaErrorDivideByZeroError:
			operaErrorString = string("Divide By Zero");
			break;
		case operaErrorZeroLength:
			operaErrorString = string("Zero length object");
			break;
			

			/* Libraries */
        case operaErrorCodeIncorrectFileType:
            operaErrorString = string("incorrect file type submitted");
            break;
		case operaErrorCodeEnvironmentnotset:
			operaErrorString = string("environment variable \"opera\" not set");
			break;
		case operaErrorCodeDatatypeNotSupported:
			operaErrorString = string("datatype not supported");
			break;
		case operaErrorCodeNoFilename:
			operaErrorString = string("no filename");
			break;
		case operaErrorCodeTileError:
			operaErrorString = string("tiling error");
			break;
		case operaErrorCodeBracketingError:
			operaErrorString = string("improper bracketing of operands, got a bool, need an operaFITSImage*");
			break;
		case operaErrorCodeNoMemory:
			operaErrorString = string("not enough memory available");
			break;
		case operaErrorCodeChangeREADONLYError:
			operaErrorString = string("attempt to modify READONLY image");
			break;
		case operaErrorCodeNULLfptr:
			operaErrorString = string("filepointer may not be NULL");
			break;
		case operaErrorInvalidGuideWindowIndex:
			operaErrorString = string("is an invalid guide window index");
			break;
		case operaErrorLengthMismatch:
			operaErrorString = string("length mismatch");
			break;
		case operaErrorNotLazy:
			operaErrorString = string("virtual image must be a lazy read");
			break;
		case operaErrorDifferingLaziness:
			operaErrorString = string("both images must have equivalent lazy status");
			break;
		case operaErrorExtensionOutOfRange:
			operaErrorString = string("extension out of range");
			break;
		case operaErrorSliceOutOfRange:
			operaErrorString = string("slice out of range");
			break;
																							
		/* Modules */
			
		case operaErrorGeometryNoOrdersFound:
			operaErrorString = string(" No orders found");
			break;			
		case operaErrorGeometryBadFit:
			operaErrorString = string("Bad polynomialfit.");
			break;			
		case operaErrorNoAnchor:
			operaErrorString = string("No anchor found.");
			break;			
									
		/* reductionset */
		case operaErrorReductionSetEtypeNotDefined:
			operaErrorString = string("Error: ETYPE option is not defined in config file");
			break;			
		case operaErrorReductionSetQualiKeyNotDefined:
			operaErrorString = string("Error: QUALIFIER_HEADERKEY is not defined in config file");
			break;	
		case operaErrorReductionSetQualiValNotDefined:
			operaErrorString = string("Error: QUALIFIER_HEADERVALUE is not defined in config file");
			break;
		case operaErrorReductionSetObstypeKeyNotDefined:
			operaErrorString = string("Error: OBSTYPE_HEADERKEY is not defined in config file");
			break;			
		case operaErrorReductionSetEtypeFailed:
			operaErrorString = string("Error: Failed to read ETYPE_HEADERVALUE");
			break;	
		case operaErrorReductionSetInputNotFound:
			operaErrorString = string("Error: Failed to access input file");
			break;

			/* Spectroscopic Libraries */
			
		case operaErrorInstrumentProfileImproperCoordinateRequest:
			operaErrorString = string("Error: improper coordinate request");
			break;				
		case operaErrorInstrumentProfileImproperDimensions:
			operaErrorString = string("Error: improper dimensions");
			break;				
		case operaErrorNoSpectralLineDetected:
			operaErrorString = string("Error: no spectral lines have been detected");
			break;						
		case operaErrorHasNoSpectralElements:
			operaErrorString = string("Error: no spectral elements have been created");
			break;   
		case operaErrorHasNoInstrumentProfile:
			operaErrorString = string("Error: no instrument profile has been created");
			break;       
		case operaErrorHasNoXCorrelation:
			operaErrorString = string("Error: no x-correlation has been calculated");
            break;
		case operaErrorHasNoGeometry:
			operaErrorString = string("Error: no geometry has been created");
            break;
		case operaErrorHasNoExtractionAperture:
			operaErrorString = string("Error: no extraction aperture has been created");
            break;
		default:
			break;
	}
	return operaErrorString;
}
/*! 
 * void operaPError(const char* prefix, const operaErrorCode errcode)
 * \brief {operaPErrorwill first print prefix followed by a colon and a space to standard error. 
 * Then, it will print the result of strerror to standard error, followed by a newline character.}
 * \param prefix is a char pointer to the prefix to use
 * \param errcode is an operaErrorCode
 * \return void
 */
static inline void operaPError(string prefix, const operaErrorCode errcode) {
	cerr << prefix << ": " << operaStrError(errcode) << '\n';
}

#else
/*! 
 * void operaStrError(const operaErrorCode errcode)
 * \brief {This function generates and reports a C-style string, 
 * containing an error message derived from the error code passed in with errcode.}
 * \param errcode is an operaErrorCode
 * \note non-reentrant
 * \return char *
 */
static char operaErrorString[MAXERRSTRINGSIZE];
static char *operaStrError(const operaErrorCode errcode) {
	if (errcode <= 100) {
		strncpy(operaErrorString, strerror(errcode), sizeof(operaErrorString));
	} else if (errcode < 600) {
		fits_get_errstatus(errcode, operaErrorString);
	} else switch (errcode) {
		case 0:
			operaErrorString[0] = '\0';
			break;
			
		/* generic */
		case operaErrorCodeNULL:
			strncpy(operaErrorString, "NULL value", sizeof(operaErrorString));
			break;
		case operaErrorCodeNULLString:
			strncpy(operaErrorString, "NULL string", sizeof(operaErrorString));
			break;
		case operaErrorInputValuesDisagree:
			strncpy(operaErrorString, "Two or more input values not in agreement", sizeof(operaErrorString));
			break;
		case operaErrorNegativeValues:
			strncpy(operaErrorString, "Values must be positive", sizeof(operaErrorString));
			break;
		case operaErrorIndexOutOfRange:
			strncpy(operaErrorString, "Index out of range", sizeof(operaErrorString));
			break;
		case operaErrorPickOutofRange:
			strncpy(operaErrorString, "pick out of range", sizeof(operaErrorString));
			break;
		case operaErrorNoInput:
			strncpy(operaErrorString, "no inputs specified", sizeof(operaErrorString));
			break;
		case operaErrorNoOutput:
			strncpy(operaErrorString, "no output specified", sizeof(operaErrorString));
			break;
		case operaErrorHeaderProblem:
			strncpy(operaErrorString, "header problem", sizeof(operaErrorString));
			break;
		case operaErrorIsNaN:
			strncpy(operaErrorString, "Invalid floating point number (NaN)", sizeof(operaErrorString));
			break;
		case operaErrorIsInf:
			strncpy(operaErrorString, "Invalid floating point number (Inf)", sizeof(operaErrorString));
			break;
		case operaErrorNoMemory:
			strncpy(operaErrorString, "No memory", sizeof(operaErrorString));
			break;
		case operaErrorInvalidInput:
			strncpy(operaErrorString, "Invalid input", sizeof(operaErrorString));
			break;
		case operaErrorExtensionHasNotBeenRead:
			strncpy(operaErrorString, "Extension Has Not Been Read", sizeof(operaErrorString));
			break;
		case operaErrorCodeFileDoesNotExistError:
			strncpy(operaErrorString, "File does not exist", sizeof(operaErrorString));
			break;
		case operaErrorCodeBadInstrumentModeError:
			strncpy(operaErrorString, "Bad instrument mode", sizeof(operaErrorString));
			break;
		case operaErrorDivideByZeroError:
			strncpy(operaErrorString, "Divide By Zero", sizeof(operaErrorString));
			break;
		case operaErrorZeroLength:
			strncpy(operaErrorString, "Zero length object", sizeof(operaErrorString));
			break;
			
			/* Modules */
			
		case operaErrorGeometryNoOrdersFound:
			strncpy(operaErrorString, "No orders found", sizeof(operaErrorString));
			break;	
		case operaErrorGeometryBadFit:
			strncpy(operaErrorString, " Bad polynomialfit.", sizeof(operaErrorString));
			break;			
		case operaErrorNoAnchor:
			strncpy(operaErrorString, "No anchor found.", sizeof(operaErrorString));
			break;			
			
			/* Libraries */

		case operaErrorCodeIncorrectFileType:
            strncpy(operaErrorString, "incorrect file type submitted", sizeof(operaErrorString));
            break;
		case operaErrorCodeEnvironmentnotset:
			strncpy(operaErrorString, "environment variable \"opera\" not set", sizeof(operaErrorString));
			break;
		case operaErrorCodeDatatypeNotSupported:
			strncpy(operaErrorString, "datatype not supported", sizeof(operaErrorString));
			break;
		case operaErrorCodeNoFilename:
			strncpy(operaErrorString, "no filename", sizeof(operaErrorString));
			break;
		case operaErrorCodeTileError:
			strncpy(operaErrorString, "tiling error", sizeof(operaErrorString));
			break;
		case operaErrorCodeBracketingError:
			strncpy(operaErrorString, "improper bracketing of operands, got a bool, need an operaFITSImage*", sizeof(operaErrorString));
			break;
		case operaErrorCodeNoMemory:
			strncpy(operaErrorString, "not enough memory available", sizeof(operaErrorString));
			break;
		case operaErrorCodeChangeREADONLYError:
			strncpy(operaErrorString, "attempt to modify READONLY image", sizeof(operaErrorString));
			break;
		case operaErrorCodeNULLfptr:
			strncpy(operaErrorString, "filepointer may not be NULL", sizeof(operaErrorString));
			break;
		case operaErrorInvalidGuideWindowIndex:
			strncpy(operaErrorString, "is an invalid guide window index", sizeof(operaErrorString));
			break;
		case MatrixNotSquare:
			strncpy(operaErrorString, "MatrixDeterminant: error: input matrix is not square", sizeof(operaErrorString));
			break;
		case MatrixInvalidDimensions:
			strncpy(operaErrorString, "MatrixDeterminant: error: dimension of input matrices are not valid", sizeof(operaErrorString));
			break;
		case MatrixZeroDeterminant:
			strncpy(operaErrorString, "MatrixDeterminant: error: determinant of input matrices is zero", sizeof(operaErrorString));
			break;
		case operaErrorLengthMismatch:
			strncpy(operaErrorString, "length mismatch", sizeof(operaErrorString));
			break;
		case operaErrorNotLazy:
			strncpy(operaErrorString, "virtual image must be a lazy read", sizeof(operaErrorString));
			break;
		case operaErrorDifferingLaziness:
			strncpy(operaErrorString, "both images must have equivalent lazy status", sizeof(operaErrorString));
			break;
		case operaErrorExtensionOutOfRange:
			strncpy(operaErrorString, "extension out of range", sizeof(operaErrorString));
			break;
		case operaErrorSliceOutOfRange:
			strncpy(operaErrorString, "slice out of range", sizeof(operaErrorString));
			break;
			
		/* reductionset */
		case operaErrorReductionSetEtypeNotDefined:
			strncpy(operaErrorString, "Error: ETYPE option is not defined in config file", sizeof(operaErrorString));
			break;
		case operaErrorReductionSetQualiKeyNotDefined:
			strncpy(operaErrorString, "Error: QUALIFIER_HEADERKEY is not defined in config file", sizeof(operaErrorString));
			break;			
		case operaErrorReductionSetQualiValNotDefined:
			strncpy(operaErrorString, "Error: QUALIFIER_HEADERVALUE is not defined in config file", sizeof(operaErrorString));
			break;				
		case operaErrorReductionSetObstypeKeyNotDefined:
			strncpy(operaErrorString, "Error: OBSTYPE_HEADERKEY is not defined in config file", sizeof(operaErrorString));
			break;				
		case operaErrorReductionSetEtypeFailed:
			strncpy(operaErrorString, "Error: Failed to read ETYPE_HEADERVALUE", sizeof(operaErrorString));
			break;	
		case operaErrorReductionSetInputNotFound:
			break;

			/* Spectroscopic Libraries */
			
		case operaErrorInstrumentProfileImproperCoordinateRequest:
			strncpy(operaErrorString, "Error: improper coordinate request", sizeof(operaErrorString));
			break;				
		case operaErrorInstrumentProfileImproperDimensions:
			strncpy(operaErrorString, "Error: improper dimensions", sizeof(operaErrorString));
			break;	
		case operaErrorNoSpectralLineDetected:
			strncpy(operaErrorString, "Error: no spectral lines have been detected", sizeof(operaErrorString));
			break;            
		case operaErrorHasNoSpectralElements:
			strncpy(operaErrorString, "Error: no spectral elements have been created", sizeof(operaErrorString));
			break;	
		case operaErrorHasNoInstrumentProfile:
			strncpy(operaErrorString, "Error: no instrument profile has been created", sizeof(operaErrorString));
			break;    
		case operaErrorHasNoXCorrelation:
			strncpy(operaErrorString, "Error: no x-correlation has been calculated", sizeof(operaErrorString));
			break;                          
		case operaErrorHasNoGeometry:
			strncpy(operaErrorString, "Error: no geometry has been created", sizeof(operaErrorString));
			break;
		case operaErrorHasNoExtractionAperture:
			strncpy(operaErrorString, "Error: no extraction aperture has been created", sizeof(operaErrorString));
			break;
		default:
			break;
	}
	return operaErrorString;
}
/*! 
 * void operaPError(const char* prefix, const operaErrorCode errcode)
 * \brief {operaPErrorwill first print prefix followed by a colon and a space to standard error. 
 * Then, it will print the result of strerror to standard error, followed by a newline character.}
 * \param prefix is a char pointer to the prefix to use
 * \param errcode is an operaErrorCode
 * \return void
 */
static inline void operaPError(const char *prefix, const operaErrorCode errcode) {
	if (prefix == NULL) {
		fprintf(stderr, "%s\n", operaStrError(errcode));
	} else {
		fprintf(stderr, "%s: %s\n", prefix, operaStrError(errcode));
	}
}

#endif // __cplusplus
#endif // OPERAERRRORCODESONLY
#endif // OPERAERROR_H
