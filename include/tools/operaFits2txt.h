#ifndef OPERAFITS2TXT_H
#define OPERAFITS2TXT_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFits2txt
 Version: 1.0
 Description: operaFits2txt unpacks data stored in upena fits 
 files to recover the original .s data in text format.
 The default operation is to process all of the fits files
 on the command line producing the associated .s files.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: eder@cfht.hawaii.edu
 
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
 http://www.gnu.or/licenses/gpl-3.0.html
 ********************************************************************/
// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

/*! \brief operaFits2txt unpacks data stored in upena products. */
/*! \file operaFits2txt.h */
/*! \ingroup tools */

#define MAXCOLS 10
/* prototypes */
static void printUsageSyntax(char *prgname);
static void print_operaFits2txt_err(int operaFits2txt_status);
static void operaFits2txt_report_error(int status, char *errtext);
static void printerror( int status);

#endif
