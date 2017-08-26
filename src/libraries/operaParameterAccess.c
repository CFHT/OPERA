/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Module name: opera
 Version: 1.0
 Description: This is a template for helping developers 
 to start up with an OPERA module.
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

/*! operaParameterAccess
 * \brief access the parameter file
 * \file operaParameterAccess.c
 * \ingroup libraries
 */

/*
 * operaParameterAccess
 * \author Doug Teeple
 * \brief access the parameter file.
 * \details {This library access the parameter file, swhich stores instrument parameters.
 * The entries consist of name := value [list] entries,
 * The parameter file is a restricted form of configuration access file,
 * in that 1. variables are not supported and 2. line line is restricted to 80
 * characters.
 *
 * The routines take a const char*name as the first parameter. The second
 * paramter is a char *&value list in the case of a "get", or a const char*
 * in the case of a "set", "add" or "delete". "remove" removes the entire
 * name/value entry.}
 */

#include <stdlib.h>
#include <regex.h>

// do not link in perror and strerror...
#define OPERAERRRORCODESONLY
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaParameterAccess.h"
#include "libraries/operaLibCommon.h"	// for startsWith

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif
static char parameterfilename[PATH_MAX];	// known to be all NULLs
static FILE *stream = NULL;					// the stream pointer to the file
// the default location of the parameter file. Now how did we find this?
static const char *parameterfilebasename = "/harness/espadons/Makefile.parameters";

/* 
 * operaErrorCode operaParameterAccessSetParamaterFilepath(const char *filepath)
 * \brief This function sets the static filepath for a module calling this library.
 * \note Modules should really ony use the default except under extenuatin circumstances.
 * \param filepath is a char pointer to the filepath or NULL to set to default
 * \return operaErrorCode
 */
operaErrorCode operaParameterAccessSetParamaterFilepath(const char *filepath) {
	if (filepath == NULL || *filepath == '\0') {
		char *prefix = getenv("opera");
		if (prefix == NULL) {	// means set the default
			strncpy(parameterfilename, "..", sizeof(parameterfilename)-1);
			strncat(parameterfilename, parameterfilebasename, sizeof(parameterfilename)-1);	
			return operaErrorCodeEnvironmentnotset;
		} else {
			strncpy(parameterfilename, prefix, sizeof(parameterfilename)-1);
			strncat(parameterfilename, parameterfilebasename, sizeof(parameterfilename)-1);
		}
	} else {
		strncpy(parameterfilename, filepath, sizeof(parameterfilename)-1);
	}
	return operaErrorCodeOK;
}

/* 
 * operaErrorCode operaParameterAccessGet(const char *name, char *&value)
 * \brief This function gets a value [list] for a given name.
 * \note Note that this function allocates storage, which must be freed by the caller.
 * \param name is a char pointer to the name
 * \param value is a char pointer address to the value [list]
 * \return operaErrorCode or errno
 * \return value = value of name or NULL if not known
 */
operaErrorCode operaParameterAccessGet(const char *name, char **value) {
	regex_t regex;
	char scanbuff[MAXPARAMETERVALUELENGTH];
	char namebuff[MAXPARAMETERVALUELENGTH];

	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;
	
	// check to see that we have a filepath set
	if (parameterfilename[0] == '\0') 
		operaParameterAccessSetParamaterFilepath(NULL);	// set in the default
	
	if (stream == NULL) {
		if ((stream = fopen(parameterfilename, "r")) == NULL) 
			return errno;
	}
	*value = NULL;
	snprintf(namebuff, sizeof(namebuff), "%s[[:space:]]*[?:]=[[:space:]]*[[:print:]]*", name);
	regcomp(&regex, namebuff, 0);
	while (fgets(scanbuff, MAXPARAMETERVALUELENGTH, stream) != NULL) {
		if (!regexec(&regex, scanbuff, 0, NULL, 0)) 
			if (startsWith(scanbuff, name) > 0) {						// we got the name
				*value = (char *)malloc(MAXPARAMETERVALUELENGTH);
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
 * operaParameterAccessSet(const char *name, const char *value)
 * \brief This function sets a value for a given name.
 * \param name is a char pointer to the name
 * \param value is a char pointer address to the value
 * \return operaErrorCode
 */
operaErrorCode operaParameterAccessSet(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;

	// check to see that we have a filepath set
	if (parameterfilename[0] == '\0') 
		operaParameterAccessSetParamaterFilepath(NULL);	// set in the default

	return operaErrorCodeNOTIMPLEMENTED;
}
/* 
 * operaParameterAccessAdd(const char *name, const char *value)
 * \brief This function adds a name value listentry .
 * \param name is a char pointer to the name
 * \param value is a char pointer to the value [list]
 * \return operaErrorCode
 */
operaErrorCode operaParameterAccessAdd(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;

	// check to see that we have a filepath set
	if (parameterfilename[0] == '\0') 
		operaParameterAccessSetParamaterFilepath(NULL);	// set in the default

	return operaErrorCodeNOTIMPLEMENTED;
}
/*
 * operaParameterAccessDelete(const char *name, const char *value)
 * \brief This function deletes a value from a given name.
 * \param name is a char pointer to the name
 * \param value is a char pointer to the value
 * \return operaErrorCode
 */
operaErrorCode operaParameterAccessDelete(const char *name, const char *value) {
	if (name == NULL || value == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;

	// check to see that we have a filepath set
	if (parameterfilename[0] == '\0') 
		operaParameterAccessSetParamaterFilepath(NULL);	// set in the default

	return operaErrorCodeNOTIMPLEMENTED;
}
/*
 * operaParameterAccessRemove(const char *name)
 * \brief This function removes a name value entry.
 * \param name is a char pointer to the name
 * \return operaErrorCode
 */
operaErrorCode operaParameterAccessRemove(const char *name) {
	if (name == NULL)
		return operaErrorCodeNULL;
	if (*name == '\0')
		return operaErrorCodeNULLString;

	// check to see that we have a filepath set
	if (parameterfilename[0] == '\0') 
		operaParameterAccessSetParamaterFilepath(NULL);	// set in the default

	return operaErrorCodeNOTIMPLEMENTED;
}

