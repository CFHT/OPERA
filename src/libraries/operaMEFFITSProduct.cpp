/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMEFFITSProduct
 Version: 1.0
 Description: class encapsulates a FITS image.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
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


/*!
 * operaMEFFITSProduct
 * \author Doug Teeple
 * \brief This class extends the MultiExtensionFITS image FITSProduct image.
 * \file operaMEFFITSProduct.cpp
 * \ingroup libraries
 */

#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMEFFITSProduct.h"

using namespace std;

/* 
 * \sa operaMEFFITSProduct()
 * \brief Basic operaMEFFITSProduct constructor.
 * \extends operaFITSImage
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(void) : operaMultiExtensionFITSImage(),
instrumentmode(MODE_UNKNOWN)
{
}

/* 
 * \sa operaMEFFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector)
 * \brief operaMEFFITSProduct constructconsturcted from a spectral order vector.
 * \extends operaFITSImage
 * \oaram Filename
 * \oaram spectralOrerVector
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector) {
}
/* 
 * \sa operaMEFFITSProduct()
 * \brief Basic unnames operaMEFFITSProduct constructor with a size.
 * \extends operaFITSImage
 * \oaram Columns
 * \oaram Rows
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(unsigned Columns, unsigned Rows, unsigned Extensions) : 
operaMultiExtensionFITSImage(Columns, Rows, Extensions, tfloat),
instrumentmode(MODE_UNKNOWN)
{
}

/* 
 * \sa operaMEFFITSProduct(string Filename, int mode=READWRITE|READONLY)
 * \brief Constructor for readng a FITS file and creating the corresponding object.
 * \param Filename containing a product (i.fits or p.fits)
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(string Filename, int mode, unsigned Compression) : 
operaMultiExtensionFITSImage(Filename, tfloat, mode, Compression, true),
instrumentmode(MODE_UNKNOWN)
{
}

/* 
 * \sa operaMEFFITSProduct(string Filename, string baseOnFilename, unsigned Rows=0)
 * \brief Constructor for readng a FITS object file and creating the corresponding product.
 * \param Filename the product file to create
 * \param baseOnFilename the object file from which to get the headers
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(string Filename, string baseOnFilename, unsigned Columns, unsigned Rows, unsigned Extensions, unsigned Compression) : 
operaMultiExtensionFITSImage(Filename, Columns, Rows, Extensions, tfloat, Compression, true),
instrumentmode(MODE_UNKNOWN)
{
	int status = 0;
	int nkeys = 0;
	char card[81];
	
	operaFITSImage basedon(baseOnFilename, tfloat, READONLY, cNone, true);
	/* copy all the user keywords (not the structural keywords) */
	int hdutype = ANY_HDU;
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_get_hdrspace(basedon.getfitsfileptr(), &nkeys, NULL, &status)) {
		throw operaException("operaMEFFITSProduct: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(basedon.getfitsfileptr(), ii, card, &status)) {
			throw operaException("operaMEFFITSProduct: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
			if (fits_write_record(fptr, card, &status)) {
				throw operaException("operaMEFFITSProduct: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	operaFITSDeleteHeaderKey("DETSIZE");
	operaFITSDeleteHeaderKey("DATASEC");
}
/* 
 * \sa operaMEFFITSProduct(string Filename, , unsigned Columns=0, unsigned Rows=0)
 * \brief Constructor for readng a FITS object file and creating the corresponding product.
 * \param Filename the product file to create
 * \param baseOnFilename the object file from which to get the headers
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(string Filename, unsigned Columns, unsigned Rows, unsigned Extensions, unsigned Compression) : 
operaMultiExtensionFITSImage(Filename, Columns, Rows, Extensions, tfloat, Compression, true),
instrumentmode(MODE_UNKNOWN)
{	
}
/* 
 * \sa operaMEFFITSProduct(string Filename, instrumentmode_t, Instrumentmode, int mode=READWRITE|READONLY);
 * \brief Constructor for creating a FITS product file with the correct column width.
 * \param Filename
 * \param mode
 * \param Instrumentmode
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaMEFFITSProduct::operaMEFFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Compression) : 
operaMultiExtensionFITSImage(Filename, tfloat, READWRITE, Compression, true),
instrumentmode(MODE_UNKNOWN)
{
	instrumentmode = Instrumentmode;
}
/* 
 * \sa operaMEFFITSProductSave()
 * \brief Save Image
 * \return void
 */
void operaMEFFITSProduct::operaMEFFITSProductSave(void) {
	operaMultiExtensionFITSImageSave();
}

/* 
 * \sa operaMEFFITSProduct()
 * \brief destructor
 * \return void
 */
operaMEFFITSProduct::~operaMEFFITSProduct(void) {
}

/*
 * Helper functions
 */

/* 
 * instrumentmode_t getmode() 
 * \brief returns intrument mode.
 * \return string
 */
instrumentmode_t operaMEFFITSProduct::getmode() {
	return instrumentmode;
}

