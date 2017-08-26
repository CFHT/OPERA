#ifndef OPERAFITSPRODUCT_H
#define OPERAFITSPRODUCT_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaFITSProduct
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

#include "libraries/operaFITSImage.h"

/*! 
 * operaFITSProduct
 * \author Doug Teeple
 * \brief This class extands the FITS image to an opera FITS Product image.
 * \file operaFITSProduct.h
 * \ingroup libraries
 */

class operaFITSProduct : public operaFITSImage {
private:
	unsigned headercolumns;
public:
	/*! 
	 * \brief Creates a FITS product in memory.
	 * \note extends operaFITSImage
	 */
	operaFITSProduct(void);
	
	/*! 
	 * \brief Creates a FITS product in memory with the specified dimensions.
	 * \note extends operaFITSImage
	 * \param Columns
	 * \param Rows
	 */
	operaFITSProduct(unsigned Columns, unsigned Rows);
	
	/*! 
	 * \brief Creates a FITS product file with the specified dimensions.
	 * \note extends operaFITSImage
	 * \param Filename
	 * \param Columns
	 * \param Rows
	 * \param Compression
	 * \throws operaException operaErrorHeaderProblem
	 */
	operaFITSProduct(string Filename, unsigned Columns, unsigned Rows, unsigned Compression = 0);
	
	/*! 
	 * \brief Creates a FITS product from an existing file.
	 * \param Filename
	 * \param mode READONLY or READWRITE
	 * \param Compression
	 * \throws operaException operaErrorHeaderProblem
	 */
	operaFITSProduct(string Filename, int mode=READWRITE, unsigned Compression = 0);
	
	/*! 
	 * \brief Copies the FITS header into the FITS product from an existing file.
	 * \param filename The existing FITS images
	 * \throws operaException operaErrorHeaderProblem
	 */
	void CopyHeader(string filename);
	
	/*! 
	 * \brief Adds a header keyword describing a column of data.
	 * \details The keyword name will be COLn where n is the number of columns that have been added, including the current one.
	 * \param name Name of the column, inserted as the value of the header keyword.
	 * \param desc Description of the column, inserted as the comment of the header keyword.
	 * \throws operaException operaErrorHeaderProblem
	 */
	void AddColumnToHeader(string name, string desc);
	
	/*! 
	 * \brief Gets the number of column keywords that have been added the header.
	 * \return The number of times AddColumnToHeader has been called.
	 */
	unsigned HeaderColumnCount();
};

#endif
