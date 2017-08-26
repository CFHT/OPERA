#ifndef OPERACONFIGURATIONACCESSLIB_H
#define OPERACONFIGURATIONACCESSLIB_H

/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaConfigurationAccess
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
 * operaConfigurationAccess
 * \author Doug Teeple
 * \brief This class interfaces to the configuration store.
 * \file operaConfigurationAccess.h
 * \ingroup libraries
 */
#ifdef __cplusplus
extern "C" {
#endif
	// the maximum value length that a configuration can have
#define MAXCONFIGURATIONVALUELENGTH 4096

operaErrorCode operaConfigurationAccessSetConfigurationFilepath(const char *filepath);
operaErrorCode operaConfigurationAccessGet(const char *name, char **value);
operaErrorCode operaConfigurationAccessSet(const char *name, const char *value);
operaErrorCode operaConfigurationAccessAdd(const char *name, const char *value);
operaErrorCode operaConfigurationAccessDelete(const char *name, const char *value);
operaErrorCode operaConfigurationAccessRemove(const char *name);

#ifdef __cplusplus
}
#endif

#endif
