/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaPolarimetryCorrection
 Version: 1.0
 Description: Create polar corrected spectrum
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

#include "core-espadons/operaExtendedSpectrumCreation.h"

/*! \file operaPolarimetryCorrection.cpp */

/*! 
 * operaPolarimetryCorrection
 * \author Eder Martioli
 * \brief Create polar corrected spectrum.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --polar=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	return ExtendedSpectrumCreation(argc, argv, "operaPolarimetryCorrection", false, true);
}

/*
 Column #3 gives the unnormalizationd Stokes parameters, in the pn.s and the
 pu.s files.
 
 To get a normalizationd Stokes parameter that goes between -1.0 and +1.0 (or
 -100% and +100%), from a pn.s file or a pu.s file, divide column #3 by
 column #2.
 
 EXAMPLE:
 
 1609679pn.s
 369.0888  2.0876e-01  6.8397e-04 -8.1249e-03  6.7757e-03  4.4399e-03
 369.0911  2.0598e-01  5.0268e-03  4.3682e-03 -1.9970e-03  3.9532e-03
 369.0935  2.0679e-01 -3.9196e-05 -7.5297e-04 -1.5541e-03  4.2366e-03
 
 1609679pu.s
 369.0888  1.5895e+00  5.2079e-03 -6.1864e-02  5.1592e-02  3.3806e-02
 369.0911  1.5683e+00  3.8274e-02  3.3259e-02 -1.5205e-02  3.0099e-02
 369.0935  1.5745e+00 -2.9844e-04 -5.7330e-03 -1.1833e-02  3.2257e-02
 
 If you try the numbers, you get the same results when you divide col #3 by
 col #2.
 
 The 'n' and the 'u' refer to the intensity (if it's normalizationd to 1.0 or
 not). But in all cases, the Stokes given is not normalizationd to anything.
 The user has to do the normalization in both cases.
 
 If the data are reduced WITHOUT the continuum polarization subtracted
 (not recommended, rarely used), Upena gives the normalizationd Stokes
 parameter. So column #3 goes between -1.0 and +1.0 (or -100% and +100%).
 */
