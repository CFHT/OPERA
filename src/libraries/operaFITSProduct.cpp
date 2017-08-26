/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFITSProduct
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
 * operaFITSProduct
 * \author Doug Teeple
 * \brief This class encapsulates the OPERA products.
 * \file operaFITSProduct.cpp
 * \ingroup libraries
 */

#include "globaldefines.h"
#include "libraries/operaFITSProduct.h"

using namespace std;

operaFITSProduct::operaFITSProduct(void) : operaFITSImage(), headercolumns(0) { }

operaFITSProduct::operaFITSProduct(unsigned Columns, unsigned Rows) : operaFITSImage(Columns, Rows, tfloat), headercolumns(0) { }

operaFITSProduct::operaFITSProduct(string Filename, unsigned Columns, unsigned Rows, unsigned Compression) : operaFITSImage(Filename, Columns, Rows, tfloat, Compression, false), headercolumns(0) { }

operaFITSProduct::operaFITSProduct(string Filename, int mode, unsigned Compression) : operaFITSImage(Filename, tfloat, mode, Compression, false), headercolumns(0)
{
	if (getnaxis1() > getnaxis2()) rotate90();
}

void operaFITSProduct::CopyHeader(string filename) {
	operaFITSImage basedon(filename, tfloat, READONLY, cNone, true);
	operaFITSImageCopyHeader(&basedon);
}

void operaFITSProduct::AddColumnToHeader(string name, string desc) {
	ostringstream ss;
	ss << "COL" << ++headercolumns;
	operaFITSSetHeaderValue(ss.str(), name, desc);
}

unsigned operaFITSProduct::HeaderColumnCount() {
	return headercolumns;
}
