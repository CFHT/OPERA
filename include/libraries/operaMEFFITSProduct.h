#ifndef OPERAMEFFITSPRODUCT_H
#define OPERAMEFFITSPRODUCT_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaMEFFITSProduct
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jul/2011
 
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

#include <string>
#include <fitsio.h>
#include "libraries/operaFITSImage.h"
#include "libraries/operaEspadonsImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralOrderVector.h"

/*! 
 * operaMEFFITSProduct
 * \author Doug Teeple
 * \brief This class extands the FITS image to an opera FITS Product image.
 * \file operaMEFFITSProduct.h
 * \package operaMEFFITSProduct
 * \ingroup libraries
 */

class operaMEFFITSProduct : public operaMultiExtensionFITSImage {
	
private:
	typedef operaMultiExtensionFITSImage& super;	// a way of referring to the super class
	instrumentmode_t instrumentmode;
	
public:
	/*! 
	 * \sa class operaMEFFITSProduct()
	 * \brief Basic operaMEFFITSProduct constructor.
	 * \extends operaFITSImage
	 * \return none
	 */
	operaMEFFITSProduct(void);
	
	/*! 
	 * \sa class operaMEFFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector)
	 * \brief operaMEFFITSProduct constructconsturcted from a spectral order vector.
	 */
	operaMEFFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector);
	
	/*! 
	 * \sa class operaMEFFITSProduct(unsigned Columns, unsigned Rows, unsigned Extensions)
	 * \brief Basic unnames operaMEFFITSProduct constructor with a size.
	 */
	operaMEFFITSProduct(unsigned Columns, unsigned Rows, unsigned Extensions);
	
	/*! 
	 * \sa class operaMEFFITSProduct(string Filename, int mode=READWRITE|READONLY)
	 * \brief Constructor for readng a FITS file and creating the corresponding object.
	 */
	operaMEFFITSProduct(string Filename, int mode=READWRITE/*READONLY*/, unsigned Compression = 0);
	/*! 
	 * \sa class operaMEFFITSProduct(string Filename, instrumentmode_t, Instrumentmode, int mode=READWRITE|READONLY);
	 * \brief Constructor for creating a FITS product file with the correct column width.
	 */
	operaMEFFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Compression = 0);
	/*! 
	 * \sa class operaMEFFITSProduct(string Filename, string baseOnFilename, unsigned Rows, unsigned Compression = 0)
	 * \brief Constructor for readng a FITS object file and creating the corresponding product.
	 */
	operaMEFFITSProduct(string Filename, string baseOnFilename, unsigned Columns, unsigned Rows, unsigned Extensions, unsigned Compression = 0);
	/*! 
	 * \sa class operaMEFFITSProduct(string Filename, int mode=READWRITE|READONLY)
	 * \brief Constructor for readng a FITS object file and creating the corresponding product.
	 */
	operaMEFFITSProduct(string Filename, unsigned Columns, unsigned Rows, unsigned Extensions, unsigned Compression = 0);
	/* 
	 * \sa operaMEFFITSProductSave()
	 * \brief Save Image
	 * \return void
	 */
	void operaMEFFITSProductSave(void);
		/*! 
	 * \sa class ~operaMEFFITSProduct()
	 * \brief destructor
	 */
	~operaMEFFITSProduct(void);
	
	/*
	 * Helper functions
	 */
	/*! 
	 * unsigned getNeextRow() 
	 * \brief returns current row.
	 */
	instrumentmode_t getmode();
};

#endif
