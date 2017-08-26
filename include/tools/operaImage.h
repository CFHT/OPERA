#ifndef OPERAIMAGE_H
#define OPERAIMAGE_H
/******************************************************************
 ****                  MODULE FOR OPERA v1.0                   ****
 ******************************************************************
 Module name: operaImage
 Version: 1.0
 Description: Return Image Statitistics.
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

/*! \brief Return Image Statitistics. */
/*! \file operaImage.h */
/*! \ingroup tools */

/* prototypes */
static void printUsageSyntax(char *prgname);
static void print_operaImage_err(int operaImage_status);
static void operaImage_report_error(int status, char *errtext);
#endif
