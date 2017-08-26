/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaImage
 Version: 1.0
 Description: Various image manipulation functions.
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
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "globaldefines.h"
#include "operaError.h"

#include "tools/operaImage.h"

/*! \brief Various image manipulation functions. */
/*! \file operaImage.cpp */
/*! \author Eder Martioli */

/*! 
 * operaImage
 * \brief Various image manipulation functions.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup tools
 */

int main(int argc, char *argv[])
{
	int operaImage_status = 0;
	int opt,i;
	int np=0, ni=0;
	char *output = NULL; // Output product 
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"param",		1, NULL, 'P'},
		{"input",		1, NULL, 'i'},
		{"output",		1, NULL, 'o'},
		
		{"verbose",		0, NULL, 'v'},
		{"plot",		0, NULL, 'p'},
		{"debug",		0, NULL, 'd'},
		{"trace",		0, NULL, 't'},
		{"help",		0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "P:i:o:vpdth", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'P':
				np++;
				break;
			case 'i':
				ni++;
				break;    
			case 'o':
				output = optarg;
				break;
				
			case 'v':
				verbose = 1;
				break;
			case 'p':
				plot = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 't':
				trace = 1;
				break;         
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	if (trace) {
		for(i=0;i<argc;i++) {
			if(i==0 || argv[i][0] == '-') {
				fprintf(stdout,"%s ",argv[i]);
			} else {
				fprintf(stdout,"\"%s\" ",argv[i]);
			}
		} 
		fprintf(stdout,"\n"); 
	}   
	
	if (argc < 4 || np == 0 || ni == 0 || output == NULL) {      
		operaImage_status = 2;
		fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__); 
		printUsageSyntax(argv[0]);
		print_operaImage_err(operaImage_status);
		exit(EXIT_FAILURE);
	} 
	
	if(verbose) {    
		fprintf(stdout,"   ....................................................... \n");	
		fprintf(stdout,"   Running module: %s\n", argv[0]);    
		fprintf(stdout,"\n   Total number of parameters: %d\n",np);
		fprintf(stdout,"\n   Total number of inputs: %d\n",ni);
		fprintf(stdout,"\n   output: %s\n", output);    
		fprintf(stdout,"   ....................................................... \n\n"); 
		
	}
	/*Start the module here*/
	
	if (debug) { 
		fprintf(stdout, "\n(%s:%s:%d)\n",__FILE__, __func__, __LINE__); 
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char *modulename) {
	
	fprintf(stderr,
			"\n"
			" Usage: %s  [-vdth] --param=<PARAMETER_VALUE_1> --param=<PARAMETER_VALUE_2> ... --output=<PRODUCT_FILE_NAME> --input=<INPUT_FILE_1> --input=<INPUT_FILE_2> ... \n\n"
			" Example: %s  -v -p 10 -p 2 --output=o.fits -i 001.fits -i 002.fits -i bad_pix.dat \n\n"
			"  -h, --help  display help message\n"
			"  -v, --verbose,  Turn on message sending\n"
			"  -d, --debug,  Turn on debug messages\n"
			"  -t, --trace,  Turn on trace messages\n"
			"  -p, --param=<PARAMETER_VALUE>, Input parameters  \n"
			"  -o, --output=<PRODUCT_FILE_NAME>, Output product file  \n"
			"  -i, --input=<INPUT_FILE_NAME>, Input files  \n\n"
			,modulename,modulename);
}

static void print_operaImage_err(int operaImage_status)
{
    /*****************************************************/
    /* Print out operaImage error messages and exit program  */
    /*****************************************************/
    char reperr[31];
    
    if (operaImage_status)
    {
		operaImage_report_error(operaImage_status,reperr); /* print error report */
		fprintf(stderr,"\n Error:operaImage: %s\n",reperr);
		
		exit(operaImage_status);    /* terminate the program, returning error status */
    }
    return;
}

static void operaImage_report_error(int status,     /* integer error status value */
							 char *errtext)  /* O - error message (max 30 char long + null) */
/*
 Return a short descriptive error message that corresponds to the input
 error status value.  The message may be up to 30 characters long, plus
 the terminating null character.
 */
{
	errtext[0] = '\0';
	
	if (status >= 0)
	{
		switch (status) {
			case 0:
				strcpy(errtext, "OK - no error\n");
				break;
			case 1:
				strcpy(errtext, "non-MODULE error\n");
				break;
			case 2:   
				strcpy(errtext, "Invalid set of arguments\n");
				break;
			case 3:   
				strcpy(errtext, "Invalid parameter value\n");
				break;       
			case 4:  
				strcpy(errtext, "Input file not found\n");
				break; 
			default:
				strcpy(errtext, "Unknown error status\n");
				break;
		}
	} else {
		strcpy(errtext, "Unknown error status\n");
	}
	return;
}

