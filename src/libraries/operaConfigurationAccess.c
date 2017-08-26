/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Module name: opera
 Version: 1.0
 Description: configuration access library 
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

/*!
 * \brief access the configuration file.
 * \file operaConfigurationAccess.c
 * \ingroup libraries
*/

#include "globaldefines.h"
#include "operaError.h"
#include <stdlib.h>
#include <regex.h>

/*!
 * operaConfigurationAccess
 * \author Doug Teeple
 * \brief access the configuration files.
 * \details {This library accesses the Configuration file, which stores instrument Configurations.
 * The entries consist of name := value [list] entries,
 * The Configuration file may contain comments, and make-style variables and the \  
 * character at end of line signifying continuation.
 *
 * The routines take a const char*name as the first Configuration. The second
 * paramter is a char *&value list in the case of a "get", or a const char*
 * in the case of a "set", "add" or "delete". "remove" removes the entire
 * name/value entry.}
*/

/*!
 * \defgroup libraries Libraries
 */

// do not link in perror and strerror...
#define OPERAERRRORCODESONLY
#include "operaError.h"
#include "libraries/operaConfigurationAccess.h"
#include "libraries/operaLibCommon.h"	// for startsWith

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif
static char configurationfilename[PATH_MAX];	// known to be all NULLs
static FILE *stream = NULL;					// the stream pointer to the file
// the default location of the configuration file. Now how did we find this?
const char *configurationfilebasename = "/harness/espadons/Makefile.configuration";

/* 
 * operaErrorCode operaConfigurationAccessSetParamaterFilepath(const char *filepath)
 * \brief This function sets the static filepath for a module calling this library.
 * \note Modules should really ony use the default except under extenuatin circumstances.
 * \param filepath is a char pointer to the filepath or NULL to set to default
 * \return operaErrorCode
 */
operaErrorCode operaConfigurationAccessSetConfigurationFilepath(const char *filepath) {
	char *prefix = getenv("opera");
	if (filepath == NULL) {
		filepath = configurationfilebasename;
	}
	if (prefix == NULL) {	// means set the default
		strncpy(configurationfilename, "..", sizeof(configurationfilename)-1);
		strncat(configurationfilename, filepath, sizeof(configurationfilename)-1);	
		return operaErrorCodeEnvironmentnotset;
	} else {
		strncpy(configurationfilename, prefix, sizeof(configurationfilename)-1);
		strncat(configurationfilename, filepath, sizeof(configurationfilename)-1);
	}		
	return operaErrorCodeOK;
}

/* 
 * operaErrorCode operaConfigurationAccessGet(const char *name, char *&value)
 * \brief This function gets a value [list] for a given name.
 * \note Note that this function allocates storage, which must be freed by the caller.
 * \param name is a char pointer to the name
 * \param value is a char pointer address to the value [list]
 * \return value = value of name or NULL if not known
 */
operaErrorCode operaConfigurationAccessGet(const char *name, char **value) {
	regex_t regex;
	char scanbuff[MAXCONFIGURATIONVALUELENGTH];
	char namebuff[MAXCONFIGURATIONVALUELENGTH];
	
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (configurationfilename[0] == '\0') 
		operaConfigurationAccessSetConfigurationFilepath(NULL);	// set in the default
	
	if (stream == NULL) {
		if ((stream = fopen(configurationfilename, "r")) == NULL) 
			return errno;
	}
	
	*value = NULL;
	snprintf(namebuff, sizeof(namebuff), "%s[[:space:]]*[?:]=[[:space:]]*[[:print:]]*", name);
	regcomp(&regex, namebuff, 0);
	while (fgets(scanbuff, MAXCONFIGURATIONVALUELENGTH, stream) != NULL) {
		if (!regexec(&regex, scanbuff, 0, NULL, 0)) 
			if (startsWith(scanbuff, name) > 0) {						// we got the name
				*value = (char *)malloc(MAXCONFIGURATIONVALUELENGTH);
				char *start = strstr(scanbuff, "=")+1;							// after the "="
				while (*start == ' ' || *start == '\t') {						// get rid of the ":="
					start++;
				}
				char *end = start;
				while (*end != '\n' && *end != '#') {							// get rid of the "\n"
					end++;
				}
				*end-- = '\0';
				while (*end == '\t' || *end == ' ') {							// get rid of trailing ws
					*end-- = '\0';
				}
				strncpy(*value, start, sizeof(scanbuff)-1); 
				break;
			}
		if (feof(stream))
			break;
	}
	if (stream)
		fclose(stream);
	stream = NULL;
	regfree(&regex);
	return operaErrorCodeOK;
}
/* 
 * operaConfigurationAccessSet(const char *name, const char *value)
 * \brief This function sets a value for a given name.
 * \param name is a char pointer to the name
 * \param value is a char pointer address to the value
 * \return operaErrorCode
 */
operaErrorCode operaConfigurationAccessSet(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (configurationfilename[0] == '\0') 
		operaConfigurationAccessSetConfigurationFilepath(NULL);	// set in the default
	
	return operaErrorCodeNOTIMPLEMENTED;
}
/* 
 * operaConfigurationAccessAdd(const char *name, const char *value)
 * \brief This function adds a name value listentry .
 * \param name is a char pointer to the name
 * \param value is a char pointer to the value [list]
 * \return operaErrorCode
 */
operaErrorCode operaConfigurationAccessAdd(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (configurationfilename[0] == '\0') 
		operaConfigurationAccessSetConfigurationFilepath(NULL);	// set in the default
	
	return operaErrorCodeNOTIMPLEMENTED;
}
/*
 * operaConfigurationAccessDelete(const char *name, const char *value)
 * \brief This function deletes a value from a given name.
 * \param name is a char pointer to the name
 * \param value is a char pointer to the value
 * \return operaErrorCode
 */
operaErrorCode operaConfigurationAccessDelete(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (configurationfilename[0] == '\0') 
		operaConfigurationAccessSetConfigurationFilepath(NULL);	// set in the default
	
	return operaErrorCodeNOTIMPLEMENTED;
}
/*
 * operaConfigurationAccessRemove(const char *name)
 * \brief This function removes a name value entry.
 * \param name is a char pointer to the name
 * \return operaErrorCode
 */
operaErrorCode operaConfigurationAccessRemove(const char *name) {
	if (name == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (configurationfilename[0] == '\0') 
		operaConfigurationAccessSetConfigurationFilepath(NULL);	// set in the default
	
	return operaErrorCodeNOTIMPLEMENTED;
}

