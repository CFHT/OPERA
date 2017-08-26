/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFits2txt
 Version: 1.0
 Description: operaFits2txt unpacks data stored in upena fits 
 files to recover the original .s and .sn data in text format.
 The default operation is to process all of the fits files
 on the command line producing the associated text files.
 
 Has dependencies on funpack and gunzip being in the user path.
 
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2011
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
#include <string.h>
#include <unistd.h>			// R_OK, access
#include <getopt.h>

#include "fitsio.h"
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLib.h"
#include "tools/operaFits2txt.h"

/*! \file operaFits2txt.cpp */
/*! \author Eder Maritoli */

/*! 
 * operaFits2txt
 * \brief command line interface to create Libre-Esprit-compatible text files from packed FITS.
 * \arg argc
 * \arg argv inputfile
 * \return EXIT_STATUS
 * \ingroup tools
 */
int main(int argc, char *argv[])
{
	int opt, i;
	char *input = NULL;			// input FITS file name 
	char *input_uncomp = NULL;	// the uncompressed filename
	char *outputu = NULL;		// output u.s file name
	char *outputn = NULL;		// output n.s file name
	char *outputuw = NULL;		// output uw.s file name
	char *outputnw = NULL;		// output nw.s file name  
	char *outputsn = NULL;		// output nw.s file name   
	unsigned int polar=0;		// polarization (=1) or intensity (=0) spectra?
	unsigned int col=0;			// print colunm names
	unsigned int debug=0, verbose=0, trace=0;
	unsigned int isCompressed = 0; // we may received either compressed or uncompressed files as input
	
	int status = 0;
	int naxis = 0, ii = 0, ncols = 0;
	long naxes[2] = {1,1};
	long npixels = 1, firstpix[2] = {1,1};
	fitsfile *fptr = NULL;
	char keyobjname[FLEN_VALUE], comment[FLEN_COMMENT], keyobsmode[FLEN_VALUE];
	char pol_str[FLEN_VALUE], sp1_str[FLEN_VALUE], sp2_str[FLEN_VALUE];
	char keycolwl[FLEN_VALUE], keycol[MAXCOLS][FLEN_VALUE];
	char colnamewl[FLEN_VALUE], colname[MAXCOLS][FLEN_VALUE];
	FILE *fout = NULL;
	double *wlpix = NULL; 
	double *pix[MAXCOLS];    
	int option=0, nopt=0;  
	int curKey=0, numKey=0, ord=0, nord=0, wl=0, snr1=0, snr2=0;
	char card[FLEN_CARD], skip[FLEN_CARD];
	
	struct option longopts[] = {     
		{"col",		0, NULL, 'c'},
		{"verbose",	0, NULL, 'v'},
		{"debug",	0, NULL, 'd'},
		{"trace",	0, NULL, 't'},
		{"help",	0, NULL, 'h'},
		{0,0,0,0}};
	
	while ((opt = getopt_long(argc, argv, "cvdt", 
							  longopts, NULL))  != -1) { 
		switch (opt) {
			case 'c':
				col = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 't':
				trace = 1;
				break;         
			default:
				printUsageSyntax(argv[0]);
				exit(EXIT_FAILURE);
				break;
		}
	}	
    
	if (trace) {
		for(i=0;i<argc;i++) {
			if (i==0 || argv[i][0] == '-') {
				fprintf(stdout,"%s ",argv[i]);
			} else {
				fprintf(stdout,"\"%s\" ",argv[i]);
			}
		} 
		fprintf(stdout,"\n"); 
	}   
	
	if (argc < 2) {      
		printUsageSyntax(argv[0]);
		exit(EXIT_FAILURE);
	} 
	
	if (verbose) {    
		fprintf(stdout,"   .................................... \n");	
		fprintf(stdout,"   Running tool: %s\n",argv[0]);    		
		fprintf(stdout,"   .................................... \n\n");			
	}	
	
	for ( ; optind < argc; optind++) {	

		unsigned int lastinputchar = 0, 
			lastuncompressedchar = 0, 
			paddedinputlength = 0,
			paddeduncompressedlength = 0;
		
		input = argv[optind];
        
		lastinputchar = strlen(input);
		paddedinputlength = lastinputchar + 1;
		
        input_uncomp = new char[paddedinputlength]; 	
        strncpy(input_uncomp, input, paddedinputlength); 
		lastuncompressedchar = strlen(input_uncomp);
		paddeduncompressedlength = lastuncompressedchar + 1;

		status = 0;
		
		if (verbose) {    
			fprintf(stdout,"   ....................................................... \n");		
			fprintf(stdout,"   Input file: %s\n",input);		
		}	
		
		if (input[lastinputchar-3] == '.' && 
			input[lastinputchar-2] == 'f' && 
			input[lastinputchar-1] == 'z') {
			
			isCompressed = 1;		     
			input_uncomp[lastinputchar-3] = '\0'; 
			lastuncompressedchar = strlen(input_uncomp);
			paddeduncompressedlength = lastuncompressedchar + 1;
			
			if (verbose) { 			
				fprintf(stdout,"   Input uncompressed file: %s\n",input_uncomp);				
			}			
			if (access(input_uncomp, R_OK)) {
				systemf("funpack %s", input);
				if (verbose)      
					fprintf(stdout,"funpack: created new uncompressed file %s\n", input_uncomp);         
			} else {
				if (verbose)      
					fprintf(stdout,"Using existing uncompressed file %s\n", input_uncomp);       
			}      
			
		} else if (input[lastinputchar-3] == '.' && 
				   input[lastinputchar-2] == 'g' && 
				   input[lastinputchar-1] == 'z') {	
			
			isCompressed = 1;
			input_uncomp[lastinputchar-3] = '\0';       
			lastuncompressedchar = strlen(input_uncomp);
			paddeduncompressedlength = lastuncompressedchar + 1;
			
			if (access(input_uncomp, R_OK)) {     
				systemf("gunzip %s -cf >%s", input, input_uncomp);
				if (verbose)      
					fprintf(stdout,"gunzip: uncompressed file %s\n", input_uncomp);         
			} else {
				if (verbose)      
					fprintf(stdout,"Using existing uncompressed file %s\n", input_uncomp);       
			}                
		} else if (input[lastinputchar-5] == '.' && 
				   input[lastinputchar-4] == 'f' && 
				   input[lastinputchar-3] == 'i' && 
				   input[lastinputchar-2] == 't' && 
				   input[lastinputchar-1] == 's') {	
			
	        isCompressed = 0;		  
			
		}  else if (input[lastinputchar-4] == '.' && 
					input[lastinputchar-3] == 'f' && 
					input[lastinputchar-2] == 'i' && 
					input[lastinputchar-1] == 't') {
			
	        isCompressed = 0;
			
		} else {
			if (debug)  
				fprintf(stdout, "\n(%s:%s:%d)\n",__FILE__, __func__, __LINE__); 
			fprintf(stdout, "Input file %s format is not supported, file has been ignored.\n",input); 
			free(input_uncomp);
			continue;		  
		}
		
		if (verbose) {    
			fprintf(stdout,"   Uncompressed file: %s\n",input_uncomp);
			fprintf(stdout,"   ....................................................... \n");			
		}			
		
		if (access(input_uncomp,R_OK)) {
			if (debug)  
				fprintf(stdout, "\n(%s:%s:%d)\n",__FILE__, __func__, __LINE__); 
			fprintf(stdout, "Input file %s has not been found.\n",input_uncomp); 
			free(input_uncomp);
			continue;
		} 
		
		if (strlen(input_uncomp) > FLEN_FILENAME) {
			if (debug)  
				fprintf(stdout, "\n(%s:%s:%d)\n",__FILE__, __func__, __LINE__);  
			fprintf(stdout, "Filename too long, file has been ignored.\n");  
			free(input_uncomp);
			continue;
		}
		
		if (fits_open_file(&fptr, input_uncomp, READONLY, &status)) {  // Open input image 
			if (debug)  
				fprintf(stdout, "\n(%s:%s:%d)\n",__FILE__, __func__, __LINE__); 
			fprintf(stdout, "Input file %s is not FITS, file has been ignored.\n",input_uncomp); 		
			free(input_uncomp);			
			continue;
		}		
		
		
		if (input_uncomp[lastuncompressedchar-6] == 'p') {
			polar = 1;
		} else {
			polar = 0;
		}
		
		outputu = new char[paddeduncompressedlength];
		outputn = new char[paddeduncompressedlength];
		outputsn = new char[paddeduncompressedlength];
		
		strncpy(outputu, input_uncomp, paddeduncompressedlength);
		strncpy(outputn, input_uncomp, paddeduncompressedlength);
		strncpy(outputsn, input_uncomp, paddeduncompressedlength);

		outputu[lastuncompressedchar-5] = 'u';    
		outputu[lastuncompressedchar-4] = '.';
		outputu[lastuncompressedchar-3] = 's';
		outputu[lastuncompressedchar-2] = '\0'; 
		
		outputn[lastuncompressedchar-5] = 'n';
		outputn[lastuncompressedchar-4] = '.';
		outputn[lastuncompressedchar-3] = 's';
		outputn[lastuncompressedchar-2] = '\0'; 
		
		outputsn[lastuncompressedchar-4] = 's';
		outputsn[lastuncompressedchar-3] = 'n';
		outputsn[lastuncompressedchar-2] = '\0'; 
		
		if (!polar) {
			outputuw = new char[paddeduncompressedlength];
			outputnw = new char[paddeduncompressedlength];
			
			strncpy(outputuw, input_uncomp, paddeduncompressedlength);
			strncpy(outputnw, input_uncomp, paddeduncompressedlength);
			
			outputuw[lastuncompressedchar-5] = 'u';     
			outputuw[lastuncompressedchar-4] = 'w';
			outputuw[lastuncompressedchar-3] = '.';
			outputuw[lastuncompressedchar-2] = 's';
			outputuw[lastuncompressedchar-1] = '\0';
			
			outputnw[lastuncompressedchar-5] = 'n';    
			outputnw[lastuncompressedchar-4] = 'w';
			outputnw[lastuncompressedchar-3] = '.';
			outputnw[lastuncompressedchar-2] = 's';
			outputnw[lastuncompressedchar-1] = '\0';     
		}
		
		if (verbose) {    
			fprintf(stdout,"   Output unormalized spectra: %s\n", outputu); 
			fprintf(stdout,"   Output normalized spectra: %s\n", outputn); 
			if (polar == 0) {    
				fprintf(stdout,"   Output unormalized spectra without autowave: %s\n", outputuw); 
				fprintf(stdout,"   Output normalized spectra without autowave: %s\n", outputnw); 
				fprintf(stdout,"   Polarization spectra?: NO\n");       
			} else if (polar == 1) {
				fprintf(stdout,"   Polarization spectra?: YES\n"); 
			}    
			fprintf(stdout,"   Output S/N file: %s\n", outputsn);     
			fprintf(stdout,"\n"); 
		}
		
		/*The module starts here*/	
		
		if (fits_get_img_dim(fptr, &naxis, &status)) // read dimensions
			printerror( status ); 
		if (fits_get_img_size(fptr, 2, naxes, &status))
			printerror( status );
		
		npixels = naxes[0];
		
		if ( fits_read_keyword(fptr, "OBJECT", keyobjname, comment, &status) )
			printerror( status );    
		if (fits_read_keyword(fptr, "INSTMODE", keyobsmode, comment, &status) )
			printerror( status );
		
		strncpy(pol_str,"\'Polarimetry, R=65,000\'\0",FLEN_VALUE);
		strncpy(sp1_str,"\'Spectroscopy, star+sky, R=65,000\'\0",FLEN_VALUE);  
		strncpy(sp2_str,"\'Spectroscopy, star only, R=80,000\'\0",FLEN_VALUE);
		
		if (polar == 0) { // here we can deal with polar(i), sp1, and sp2
			if (strcmp(keyobsmode,pol_str) == 0 || strcmp(keyobsmode,sp2_str) == 0) {
				ncols = 3;
			} else if (strcmp(keyobsmode,sp1_str) == 0) {
				ncols = 7;  
			} else {
				if (debug)
					fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
				print_operaFits2txt_err(8);  
				return(EXIT_FAILURE);           
			}
		} else if (polar == 1) { // and here we deal with polar(p)
			if (strcmp(keyobsmode,pol_str) == 0) {
				ncols = 6;
			} else if (strcmp(keyobsmode,sp1_str) == 0 || strcmp(keyobsmode,sp2_str) == 0) {
				if (debug)
					fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
				print_operaFits2txt_err(5);  
				return(EXIT_FAILURE);    
			} else {
				if (debug)
					fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
				print_operaFits2txt_err(8);  
				return(EXIT_FAILURE);       
			}
		} 
		if (verbose) {
			fprintf(stdout,"   Observing mode detected: %s\n",keyobsmode);
			fprintf(stdout,"   ....................................................... \n\n"); 			
		}			
		
		if (ncols*4 != naxes[1]) {
			if (debug)
				fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
			print_operaFits2txt_err(6);  
			return(EXIT_FAILURE);       
		}
		
		wlpix = new double[npixels]; /* mem for 1 row */
		if (wlpix == NULL) {
			if (debug)
				fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
			print_operaFits2txt_err(7);  
			return(EXIT_FAILURE); 
		}
		
		for (i=1; i<ncols; i++) {
			pix[i] = new double[npixels]; /* mem for 1 row */
			if (pix[i] == NULL) {
				if (debug)
					fprintf(stderr, "\nError: (%s:%s:%d)\n", __FILE__, __func__, __LINE__);    
				print_operaFits2txt_err(7);  
				return(EXIT_FAILURE); 
			}    
		}
		
		if (polar == 0) {
			nopt = 4;
		} else if (polar == 1) {
			nopt = 2;
		}   
		
		if (fits_get_hdrspace(fptr, &numKey, NULL, &status))
			printerror( status );
		
		nord = 0;
		fout = fopen(outputsn,"w");
		fprintf(fout,"***\n");  
		fprintf(fout,"%d %d\n",40,1);    
		
		for(curKey = 1; curKey <= numKey; curKey++)  { 
			// read the current keyword 
			if (fits_read_record(fptr, curKey, card, &status))
				printerror( status );   
			
			if (polar == 0) {  
				if (!strncmp("COMMENT SNR", card, 11)) {       
					sscanf(card,"COMMENT SNR per spec/ccd pxl in order # %d (%d nm): I> %d / %d",&ord,&wl,&snr1,&snr2);
					fprintf(fout,"%.6f %.6f\n",(float)wl,(float)snr2);          
					nord++;          
					if (nord == 40)
						break;          
				}         
			} else if (polar == 1) {
				if (!strncmp("COMMENT SNR", card, 11)) {       
					sscanf(card,"COMMENT SNR per spec/ccd pxl in order # %d (%d nm)",&ord,&wl);
				}      
				if (!strncmp("COMMENT  I>", card, 11)) {       
					sscanf(card,"COMMENT  I> %d / %d %[^\n]",&snr1,&snr2,skip);
					fprintf(fout,"%.6f %.6f\n",(float)wl,(float)snr2);           
					nord++;          
					if (nord == 40)
						break;
				}           
			}      
		}
		fclose(fout);
		
		for(option=0;option<nopt;option++)
		{
			/********* Read data *************/  
			firstpix[1] = option*ncols + 1;
			if (fits_read_pix(fptr, TDOUBLE, firstpix, npixels, NULL, wlpix,NULL, &status)) 
				printerror( status );  
			
			if (col == 1) {
				sprintf(colnamewl,"COL%ld",firstpix[1]);
				if (fits_read_keyword(fptr,colnamewl, keycolwl, comment, &status) )
					printerror( status );
			}
			
			for(i=1; i<ncols; i++) {  
				firstpix[1] = option*ncols + i + 1;
				if (fits_read_pix(fptr, TDOUBLE, firstpix, npixels, NULL, pix[i],NULL, &status)) 
					printerror( status );
				
				if (col == 1) {
					sprintf(colname[i],"COL%ld",firstpix[1]);
					if (fits_read_keyword(fptr,colname[i], keycol[i], comment, &status) )
						printerror( status ); 
				}     
			}       
			
			/*********************************/   
			//--- Write output file ---------              
			if (option == 0) { //norm, wave
				fout = fopen(outputn,"w");        
			} else if (option == 1) { //no-norm, wave
				fout = fopen(outputu,"w");
			} else if (option == 2) { //norm; no-wave
				fout = fopen(outputnw,"w");
			} else if (option == 3) { //no-norm; no-wave
				fout = fopen(outputuw,"w");
			}      
			
			fprintf(fout,"***Reduced spectrum of %s\n",keyobjname);  
			fprintf(fout," %ld %d\n",npixels,ncols-1);
			if (col == 1) {      
				fprintf(fout,"%s",keycolwl);
				for(i=1;i<ncols;i++)
					fprintf(fout," %s",keycol[i]);
				fprintf(fout,"\n");      
			}  
			//-------------------------------
			for(ii=0; ii< npixels; ii++) {
				fprintf(fout,"  %8.4f",wlpix[ii]);
				for(i=1;i<ncols;i++)
					fprintf(fout," %11.4e",pix[i][ii]);
				fprintf(fout,"\n");      
			}
			fclose(fout);
			//--------------------------------  
		}
		
		free(wlpix);
		for(i=1;i<ncols;i++)
			free(pix[i]);
		
		if ( fits_close_file(fptr, &status) )
			printerror( status );
		
		if (isCompressed) {
			unlink(input_uncomp);
		}
		
		delete[] input_uncomp;
		delete[] outputu;
		delete[] outputn;
		
		if (!polar) {    
			delete[] outputuw;
			delete[] outputnw;
		}   
		
		delete[] outputsn;
		
	} //end of for ( ; optind < argc; optind++)
	
	return EXIT_SUCCESS;
}  


/*************** Module functions start here *****************/

/* Print out the proper program usage syntax */
void printUsageSyntax(char *modulename) {
	
	fprintf(stderr,
			"\n"
			" Usage: %s  [-cvdth] <INPUT_FILE_NAME(s)> \n\n"
			" Example: %s -t *i.fits \n\n"
			"  -c, --col,      Print column names\n"
			"  -v, --verbose,  Turn on message sending\n"
			"  -d, --debug,    Turn on debug messages\n"
			"  -t, --trace,    Turn on trace message\n"
			"  -h, --help      Display help message\n\n"      
			,modulename,modulename);
}

/* Print out operaFits2txt error messages and exit program  */
void print_operaFits2txt_err(int operaFits2txt_status)
{
	char reperr[1024];
	
	if (operaFits2txt_status)
	{
		operaFits2txt_report_error(operaFits2txt_status,reperr); /* print error report */
		fprintf(stderr,"\n Error:operaFits2txt: %s\n",reperr);
		
		exit(operaFits2txt_status);    /* terminate the program, returning error status */
	}
	return;
}

/*
 Return a short descriptive error message that corresponds to the
 input error status value.  The message may be up to 30 characters long. 
 */
void operaFits2txt_report_error(int status,     /* integer error status value */
								char *errtext)  /* O - error message (max 30 char long + null) */
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
				strcpy(errtext, "Input file not found\n");
				break;       
			case 4:  
				strcpy(errtext, "Invalid obsmode. Use only \"pol\", \"sp1\", or \"sp2\".\n");
				break; 
			case 5:  
				strcpy(errtext, "Option \"--polar\" (-p) is only valid for obsmode=\"pol\"\n");
				break;     
			case 6:  
				strcpy(errtext, "FITS file does not match dimensions for this obsmode\n");
				break; 
			case 7:  
				strcpy(errtext, "Memory allocation error\n");
				break;   
			case 8:  
				strcpy(errtext, "Instrument mode in header not recognized.\n");
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

void printerror( int status)
{
	/*****************************************************/
	/* Print out cfitsio error messages and exit program */
	/*****************************************************/
	if (status)
	{
		fits_report_error(stderr, status); /* print error report */
		
		exit( status );    /* terminate the program, returning error status */
	}
	return;
}

