#ifndef OPERAESPADONSSUBIMAGE_H
#define OPERAESPADONSSUBIMAGE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: liboperaespadonsSubImage
 Class: operaEspadonsSubImage
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

#include <libraries/operaEspadonsSubImage.h>

using namespace std;

/*!
 * operaEspadonsSubImage class
 * \author Doug Teeple
 * \brief This class encapsulates an espadons subimage and extends operaFITSSubImage.
 * \file operaEspadonsSubImage.h
 * \ingroup libraries
 */
class operaEspadonsSubImage {
		
public:
	/*!
	 * \class operaEspadonsSubImage()
	 * \brief Basic operaEspadonsSubImage constructor.
	 * \return void
	 */
	operaEspadonsSubImage();														// simply construct a default image
};

#endif
