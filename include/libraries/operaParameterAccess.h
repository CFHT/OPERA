#ifndef OPERAPARAMETERACCESSLIB_H
#define OPERAPARAMETERACCESSLIB_H

/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaImage
 Version: 1.0
 Description: ThisC library implements low level image routines..
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
 * operaParameterAccess
 * \author Doug Teeple
 * \brief This class interfaces to the parameters store.
 * \file operaParameterAccess.h
 * \ingroup libraries
 */
#ifdef __cplusplus
extern "C" {
#endif
	// the maximum value length that a parameter can have
#define MAXPARAMETERVALUELENGTH 4096

operaErrorCode operaParameterAccessSetParamaterFilepath(const char *filepath);
operaErrorCode operaParameterAccessGet(const char *name, char **value);
operaErrorCode operaParameterAccessSet(const char *name, const char *value);
operaErrorCode operaParameterAccessAdd(const char *name, const char *value);
operaErrorCode operaParameterAccessDelete(const char *name, const char *value);
operaErrorCode operaParameterAccessRemove(const char *name);

#ifdef __cplusplus
}
#endif

#endif

