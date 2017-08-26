/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaQueryImageInfo
 Version: 1.0
 Description: This module produces a table with information from images that match 
 a set of pre-defined qualifiers
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
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
#include <getopt.h>
#include <dirent.h>
#include <unistd.h>			// R_OK, access
#include <regex.h>
#include <iostream>

#include "fitsio.h"

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaQueryImageInfo.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaException.h"
#include "libraries/operaLibCommon.h"	// for startsWith
#include "libraries/operaConfigurationAccess.h"

/*! \file operaQueryImageInfo.cpp */
/*! \ingroup tools */

using namespace std;

/* Print out the proper program usage syntax */
static void printUsageSyntax(){	
	
	cout <<
					"\n"
					" Usage: ./operaQueryImageInfo [-hvdt] [-p] [--input=<FILENAME>] [--output=<FILENAME>] [--directory=<DATA_DIR>] --extractkeys=<FITS_KEYWORDS> --qualifierkeys=<FITS_KEYWORDS> --splitkey=<FITS_KEYWORD>\n\n"
					" Example: ./operaQueryImageInfo --directory=/Users/espadons/data1  -r /Users/espadons/data2 -q \"ORIGIN OBSTYPE\" -e \"OBSTYPE EXPTIME\" ORIGIN=CFHT OBSTYPE=FLAT\n"
					"\n"
					"  -h, --help  display help message\n"
					"  -v, --verbose,  Turn on message sending\n"
					"  -d, --debug,  Turn on debug messages\n"
					"  -t, --trace,  Turn on trace messages\n"
					"  -p, --printheader,  Print header for output extracted values\n"	
					"\n" 
					"  -i, --input, Input file with list of paths (optional: override --directory)\n"
					"  -o, --ouput, Output file \n" 
					"  -r, --directory, Input data directory(ies) from command line\n"
					"  -q, --qualifierkeys, list of FITS header keywords to select dataset\n" 
					"  -e, --extractkeys, list of FITS header keywords to extract values\n"
					"  -s, --splitkey,  Define header keyword for which to insert break point upon change\n"	
					"\n";
}

/* 
 * operaQueryImageInfo
 * \author Eder Martioli
 * \brief Creates a data table with info from a list of input FITS file names.
 * \arg argc
 * \arg argv
 * \note --directory=$(DATADIR)
 * \note --qualifierkeys="$(DETECTOR) $(MODE) $(SPEED) ..."
 * \note --exractkeys="$(DETECTOR) $(MODE) $(SPEED) ..." 
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
int main(int argc, char *argv[])
{
	int opt,i,j;
	bool debug=false, verbose=false, trace=false;
	bool printheader=false;
	operaErrorCode errorcode = operaErrorCodeOK;;
	
	string splitkey = ""; // define keywords to set a splitkey in the output table upon every change in they keyword value 
	
	char *listofheaderkeystoextract = NULL; // list of FITS header keyword names to extract values onto output table
	char *listofheaderkeyqualifierstomatch = NULL; // this qualifier matches the keyword value

	FILE *fout = NULL;
	
	char *input = NULL; // input list of file paths
	char *output = NULL; // output file - default is to stdout
	
	char *dirs[MAXNDIRS];
	int ni=0;
	int ndirs = 0;
	
	struct option longopts[] = {       
		{"input",1, NULL, 'i'}, 
		{"output",1, NULL, 'o'}, 		
		{"directory",1, NULL, 'r'}, 			
		{"qualifierkeys",1, NULL, 'q'}, 
		{"extractkeys",1, NULL, 'e'},
		{"splitkey",1,NULL,'s'},
		{"printheader",0,NULL,'p'},		
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	i = 1;
	
	while((opt = getopt_long(argc, argv, "i:o:r:q:e:s:pvdth", 
													 longopts, NULL))  != -1)
	{
		switch(opt) 
		{      
			case 'i':
				input = optarg;
				break;
			case 'o':
				output = optarg;
				break;				
			case 'r':
				dirs[ni++] = optarg;
				break;
			case 'q':
				listofheaderkeyqualifierstomatch = optarg;
				break;
			case 's':
				splitkey = optarg;
				break;            
			case 'e':
				listofheaderkeystoextract = optarg;
				break;   
			case 'p':
				printheader = true;
				break;				
			case 'v':
				verbose = true;
				break;
			case 'd':
				debug = true;
				break;
			case 't':
				trace = true;
				break;         
			case 'h':
				printUsageSyntax();
				exit(EXIT_FAILURE);
				break;
			case '?':
				printUsageSyntax();
				exit(EXIT_FAILURE);
				break;
		}
	}	
	
	try {
		if (ni != 0 && input == NULL) {
			ndirs = ni;
		}
		if (input != NULL) {
			ndirs = 0;
		}
		
		if (argc < 3) {      
			printUsageSyntax();
			return(EXIT_FAILURE);
		} 
		
		if (debug) {
			for (i=0;i<argc;i++) {
				if (i==0 || argv[i][0] == '-') {
					cout << argv[i] << " ";
				} else {
					cout << "\"" << argv[i] << "\" ";
				}
			} 
			cerr<<"\n"; 
		} 
		
		
		if (verbose) {  
			cout << 
			"   ....................................................... \n"	
			"   Running module: " << argv[0] << "\n" 
			"   Data directories:\n";
			for (i=0;i<ndirs;i++) {
				cout << "      " << dirs[i] << "\n";
			} 	
			if (listofheaderkeystoextract != NULL)
				cout << " List of keys to extract: " << listofheaderkeystoextract << "\n"; 		
			if (listofheaderkeyqualifierstomatch != NULL)
				cout << " List of qualifier keys: " << listofheaderkeyqualifierstomatch << "\n";      
			cout << "   ....................................................... \n";	    
		}
		
		//----- read keys for which the values will be extracted -----------  		
		char *extractheaderkey[MAXCONFIGVALUES];
		int nextractkeys=0; 		
	  // search for command line options
		if(listofheaderkeystoextract != NULL) {
			nextractkeys = SplitValues(listofheaderkeystoextract, extractheaderkey, MAXCONFIGVALUES);  		
		}
		//----- read qualifier keys -----------  
		char *qualiheaderkey[MAXCONFIGVALUES];
		char *qualiheaderkeyvalues[MAXCONFIGVALUES];
		int nqualikeys=0; 
		
	  // search for command line qualifiers
		if (listofheaderkeyqualifierstomatch != NULL){
			nqualikeys = SplitValues(listofheaderkeyqualifierstomatch, qualiheaderkey, MAXCONFIGVALUES);  
			
			// parse the command line arguments to search for values of each qualifier to be matched 
			regex_t regex;
			
			char namebuff[MAXCONFIGURATIONVALUELENGTH];	
			
			for (j=0;j<nqualikeys;j++) {
				qualiheaderkeyvalues[j] = NULL;
				snprintf(namebuff, sizeof(namebuff), "%s=[[:print:]]", qualiheaderkey[j]);
				regcomp(&regex, namebuff, REG_ICASE);		
				for (i=1;i<argc;i++) {		
					if (!regexec(&regex, argv[i], 0, NULL, 0)) {		// we got the name
						qualiheaderkeyvalues[j] = new char[MAXCONFIGURATIONVALUELENGTH];
						char *start = strstr(argv[i], "=")+1;							// after the "="
						strncpy(qualiheaderkeyvalues[j], start, MAXCONFIGURATIONVALUELENGTH); 
						break;
					}
				}
			}
			
			if (debug) {
				for (i=0;i<nqualikeys;i++) {
					cout << "operaQueryImageInfo: \n"
					<< "QUALIFIER # " << i << endl
					<< "KEYWORD:         " << qualiheaderkey[i] << endl
					<< "VALUE TO MATCH:    " << qualiheaderkeyvalues[i] << endl;
				}
			}
		}
		//------------------------------------------
		// Open input directories and get file names
		//------------------------------------------
		if (debug && ndirs)
			cout << "operaQueryImageInfo: Opening input directory and reading file names...\n";
		if (debug && input)
			cout << "operaQueryImageInfo: Opening input file list and reading file names... " << input << "\n";
		
		struct dirent *entry;
		DIR *dp[MAXNDIRS];
		
		for (i=0;i<ndirs;i++) {
			dp[i] = opendir(dirs[i]);
			if (dp[i] == NULL) {
				throw operaException("operaQueryImageInfo: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
		
		int nfiles = 0;
		unsigned maxsize = 0;
		
		for (i=0;i<ndirs;i++) {
			// get number of files and maximum size 
			while((entry = readdir(dp[i])))
			{
				if (maxsize <= (strlen(entry->d_name)+strlen(dirs[i])+3))
					maxsize = (strlen(entry->d_name)+strlen(dirs[i])+3);
				nfiles++;     
			}
			closedir(dp[i]);
		}
		
		char fullpathname[MAXNFILES][MAXDIRNAMESIZE];
		char filename[MAXNFILES][MAXDIRNAMESIZE];
		
		FILE *flist;
		char line[MAXDIRNAMESIZE];		
		j=0;
		//---------------------------------------------------------------
		// Try to open input list and get file paths and ignore directory
		//---------------------------------------------------------------
		if (input == NULL) {
			for (int i=0;i<ndirs;i++) {
				// reopen directory    
				dp[i] = opendir(dirs[i]);
				// get file names
				while((entry = readdir(dp[i])))
				{
					snprintf(filename[j],strlen(entry->d_name)+1,"%s",entry->d_name);
					snprintf(fullpathname[j],strlen(entry->d_name)+strlen(dirs[i])+2,"%s%s",dirs[i],entry->d_name); 
					j++;
				}
				closedir(dp[i]);
			}		
		} else {  //get file paths from input file	
			if (access(input,R_OK)) {
				throw operaException("operaQueryImageInfo: ", operaErrorReductionSetInputNotFound, __FILE__, __FUNCTION__, __LINE__);	
			} else {		
				flist = fopen(input,"r");
				nfiles = 0;
				while(fgets(line,sizeof(line),flist) != NULL){nfiles++;}
				fseek(flist,0L,SEEK_SET);	
				for (j=0;j<nfiles;j++) {
					fscanf(flist,"%s\n",fullpathname[j]);
				}
				fclose(flist);
			}
		}
		//----------------------------------------------		
		
		if (debug) {
			for (i=0;i<nfiles;i++)
				cout << "operaQueryImageInfo: Trying to open file "<<fullpathname[i]<<"\n";  
		}
		
		char qualival_from_FITS[MAXCONFIGVALUES][FLEN_VALUE],extractval_from_FITS[MAXCONFIGVALUES][FLEN_VALUE];
		int FitStatus = 0, accept_quali;
		char comment[FLEN_COMMENT];  
		fitsfile *fptr;
		
		if (output != NULL) {
			fout = fopen(output,"w");
		}
		
		if (debug)
			cout << "operaQueryImageInfo: Parsing files...\n";
		
		unsigned char split = 0;			// set to 1 when we find this split key value in the stack
		unsigned modechangecount = 0;
		char currentsplitkeyvalue[FLEN_VALUE] = {'\0'};
		char splitkeyvaluestack[MAX_MODE_CHANGES][FLEN_VALUE] = 
		{{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},
			{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},
			{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'}};
		
		
		if(printheader){
			if (output == NULL) {
				cout << "FILENAME";	
				for(i=0;i<nextractkeys;i++) {
					cout << "\t" << extractheaderkey[i];
				}
				cout << endl;			
			} else {
				fprintf(fout,"FILENAME");
				for(i=0;i<nextractkeys;i++) {
					fprintf(fout,"\t%s",extractheaderkey[i]);
				}
				fprintf(fout,"\n");				
			}
		}
		
		for (int i=0;i<nfiles;i++) {

			FitStatus = 0;
			char asplitkeyvalue[FLEN_VALUE];
			
			if (debug)
				cout << "operaQueryImageInfo: Trying to open file: " << fullpathname[i] << "\n";
			
			if (fits_open_file(&fptr, fullpathname[i], READONLY, &FitStatus)) {
				if (debug && FitStatus)	
					operaPError("operaQueryImageInfo", FitStatus); 
				continue;        
			}
			
			if (debug)
				cout << "operaQueryImageInfo: Open succeeded, now reading FITS file: " << fullpathname[i] << "\n";
			
			//
			// if we have a split key then watch for triggers
			//
			if (!splitkey.empty()) {
				if (fits_read_keyword(fptr, splitkey.c_str(), asplitkeyvalue, comment, &FitStatus) ) {
					if (FitStatus)
						operaPError("operaQueryImageInfo: Could not find split key "+splitkey+" in header ", FitStatus); 	 
				} else {
					if (strlen(currentsplitkeyvalue) == 0) {
						strcpy(currentsplitkeyvalue, asplitkeyvalue); // first one, grab it
						strcpy(splitkeyvaluestack[modechangecount++], asplitkeyvalue);
					} else {
						if (!strcmp(asplitkeyvalue, currentsplitkeyvalue)) {	// same = no change
						} else {	// whoopsie -- the split key  value changed! Note the split
							for (unsigned k=0; k<modechangecount; k++) {
								if (!strcmp(splitkeyvaluestack[k], currentsplitkeyvalue)) {	// we have seen this mode before
									split = 1;
									break;
								}
							}
							strcpy(currentsplitkeyvalue, asplitkeyvalue); // get the new trigger value
							strcpy(splitkeyvaluestack[modechangecount++], asplitkeyvalue);
						}
					}
				}
			}		
			for(unsigned jj=0;jj<(unsigned)nextractkeys;jj++) {
				if (fits_read_keyword(fptr, extractheaderkey[jj], extractval_from_FITS[jj], comment, &FitStatus) ) {
					if (FitStatus){
						sprintf(extractval_from_FITS[jj],"NONEXISTENT_KEYWORD");
						if(debug)		
							cout << "operaQueryImageInfo: WARNING: keyword " << extractheaderkey[jj] << " could not be found. \n";	 
					}	
				}
			}			
	
			accept_quali = 1;
			for (int ii=0;ii<nqualikeys;ii++) {
				if (qualiheaderkeyvalues[ii] != NULL) {
					qualival_from_FITS[ii][0] = '\0';
					FitStatus = 0;
					
					if (fits_read_keyword(fptr, qualiheaderkey[ii], qualival_from_FITS[ii], comment, &FitStatus) ) {
						if (FitStatus){
							if(debug)
								cout << "operaQueryImageInfo: WARNING: keyword " << qualiheaderkey[ii] << " could not be found and has been ignored. \n";	 
							continue;
						}
					}
					// condition based on qualifiers: quali1 && quali2 && quali3 &&..
					if (debug)
						cout << "operaQueryImageInfo: QUALIFIER: " << qualiheaderkey[ii] << " testing " << qualiheaderkeyvalues[ii] << " to match " << qualival_from_FITS[ii] << "\n";
					
					if (qualival_from_FITS[ii][0] != '\0' && qualiheaderkeyvalues[ii][0] != '\0') {
						// handle wildcards - a qualifier ending in %, also get rid of the '
						if (qualiheaderkeyvalues[ii][strlen(qualiheaderkeyvalues[ii])-1] == '%') {
							char temp[FLEN_VALUE];
							strncpy(temp,qualiheaderkeyvalues[ii],FLEN_VALUE);
							temp[strlen(temp)-1] = '\0';	// remove the %
							if (startsWith(qualival_from_FITS[ii], temp) > 0) {
								if (debug)
									cout << "operaQueryImageInfo: qualifier: " << qualiheaderkeyvalues[ii] << " accepted == " << qualival_from_FITS[ii] << " temp=" << temp << "\n";	  
							} else {
								if (debug)
									cout << "operaQueryImageInfo: qualifier: " << qualiheaderkeyvalues[ii] << " NOT accepted  != " << qualival_from_FITS[ii] << "\n";	  
								accept_quali = 0;
								break;
							}
						} else {	// exact match case							
							char tempstr[FLEN_VALUE];								
							cleanFITSHeaderValue(qualival_from_FITS[ii],tempstr);

							if (!strcmp(qualiheaderkeyvalues[ii], tempstr)) {
								if (debug)
									cout << "operaQueryImageInfo: qualifier: " << qualiheaderkeyvalues[ii] << " accepted == " << tempstr << "\n";	  
							} else {
								if (debug)
									cout << "operaQueryImageInfo: qualifier: " << qualiheaderkeyvalues[ii] << " NOT accepted  != " << tempstr << "\n";	  
								accept_quali = 0;
								break;
							}
						}
						
					} else {
						if (debug)
							cout << "operaQueryImageInfo: Warning: could not find header value for qualifier "<< qualiheaderkey[ii] << " in file "<<filename[i]<<"\n";
					}
				}
			}
			
			if (accept_quali==1) {
				if (output == NULL) {
					if (split) {
						split = 0;
						cout << "##########\n";				
					}
					if(listofheaderkeystoextract != NULL) {
						cout << filename[i];	
						for(j=0;j<nextractkeys;j++) {
							if(extractval_from_FITS[j]!=NULL){
								
								char cleanstr[FLEN_VALUE];								
								cleanFITSHeaderValue(extractval_from_FITS[j],cleanstr);
								
								cout << "\t" << cleanstr;
							} else {
								cout << "\t" << "BLANK_VALUE";
							}
						}
						cout << endl;							
					} else {
						cout << fullpathname[i] << endl;				
					}
				} else {
					if (split) {
						split = 0;
						fprintf(fout,"##########\n");
					}
					if(listofheaderkeystoextract != NULL) {
						fprintf(fout,"%s",filename[i]);
						for(j=0;j<nextractkeys;j++) {
							if(extractval_from_FITS[j]!=NULL){
								char cleanstr[FLEN_VALUE];								
								cleanFITSHeaderValue(extractval_from_FITS[j],cleanstr);								
								fprintf(fout,"\t%s",cleanstr);
							} else {
								fprintf(fout,"\tBLANK_VALUE");
							}
						}
						fprintf(fout,"\n");	
					} else {
						fprintf(fout,"%s\n",fullpathname[i]);				
					}					
				}
			}	
			
			if (fits_close_file(fptr, &FitStatus)) {
				if (debug)
					operaPError("operaQueryImageInfo", FitStatus );
			}	
		} // end of for "nfiles" loop
		if (fout) {
			fclose(fout);
		}
	} catch (operaErrorCode errorcode) {
		if (fout) {
			fclose(fout);
		}
		operaPError("operaQueryImageInfo", errorcode );
	}
	return errorcode;
}

unsigned SplitValues(char *parvalues, char *values[], unsigned maxvalues) {
	unsigned i, fullstr;
	unsigned npars;
	unsigned spc=0,j;
	char limchar = '\0';
	
	fullstr = 0;
	j=0; 	
	for (i=0;i<=strlen(parvalues);i++) {
		if ((parvalues[i] == '\'' || parvalues[i] == '\"') && fullstr==0) {
			limchar = parvalues[i];
			fullstr = 1;					
		} else if (fullstr == 1 && parvalues[i] == limchar) {
			fullstr = 0;
		} else {	
			if (j==0){
				values[spc] = new char[strlen(parvalues)];
			}	
			if (fullstr == 0) {
				if (parvalues[i] == ' ' || parvalues[i] == '\0') {	    				
					values[spc][j] = '\0';
					j=0;
					spc++;
					if (spc > maxvalues)
						break;
					continue;
				}
			}
			values[spc][j] = parvalues[i];
			j++;
		}
	}
	npars = spc;
	return npars;	
}

unsigned operaQueryImageInfoConfigurationAccess(char *par, char *values[], unsigned maxvalues, operaErrorCode *errorcode) {
	
	unsigned npars = 0;
	char *parvalues = NULL;
	
	*errorcode = operaConfigurationAccessGet(par, &parvalues);
	
	if (*errorcode == operaErrorCodeOK) {	
		npars = SplitValues(parvalues, values, maxvalues);
	}
	free(parvalues);
	return npars;
}

void cleanFITSHeaderValue(char *value_from_FITS,char *output_clean_value) {
	unsigned inew=0;
	unsigned inilen = strlen(value_from_FITS);
	// removing single quote							
	for(unsigned ichar=0;ichar<inilen;ichar++) {
		if(value_from_FITS[ichar] != '\''){
			output_clean_value[inew++] = value_from_FITS[ichar];
		}
	}
	output_clean_value[inew] = '\0';
	
	// removing blank space
	unsigned len = strlen(output_clean_value) - 1;
	while (isspace(output_clean_value[len])) {len--;}
	output_clean_value[len + 1] = '\0';
}
