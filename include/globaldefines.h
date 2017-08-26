#ifndef OPERAGLOBALDEFINES_H
#define OPERAGLOBALDEFINES_H

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: globaldefines
 Version: 1.0
 Description: Global definitions that apply to all libraries and modules.
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

/*! \brief Global definitions that apply to all libraries and modules. */
/*! \file globaldefines.h */
/*! \ingroup libraries */

/*
 * Print debug info.
 */
#undef PRINT_DEBUG

/*
 * Print put of bounds conditions to stderr.
 */
#undef PRINT_OUTOFBOUNDS

/*
 * Range check indices or run fast and loose.
 */
#undef RANGE_CHECK
#define FITS_RANGE_CHECK

/*
 * When calculating indices of extensions and slices we
 * have a choice of fast but less general access
 * or full-blown (but slow) access.
 * Fast means you can't index two different slices or 
 * two different extensions within a single expression
 * using the [] operator.
 */
#undef FAST_EXTENSIONS 
#undef FAST_SLICES

/*
 * Compress intermediate products.
 */
#define GZSTREAM

#ifdef GZSTREAM
#define operaostream ogzstream
#define operaistream igzstream
#else
#define operaostream ofstream
#define operaistream ifstream
#endif

#endif // OPERAGLOBALDEFINES_H
