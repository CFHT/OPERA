/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaReductionSet
 Version: 1.0
 Description: This module produces a list of FITS files that matches 
 a set of pre-defined qualifiers and also that matches a 
 given observation type (like FLAT, OBJECT, etc.)
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

#include <stdio.h>
#include <getopt.h>
#include <unistd.h>			// R_OK, access
#include <dirent.h>
#include <regex.h>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaReductionSet.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLibCommon.h"	// for startsWith
#include "libraries/operaException.h"
#include "libraries/operaConfigurationAccess.h"

/*! \file operaReductionSet.cpp */

using namespace std;

/* Print out the proper program usage syntax */
static void printUsageSyntax(){	
	cout <<
					"\n"
					" Usage: ./operaReductionSet [-hvdt] [--input=<LIST_OF_PATHS>] [--output=<OUTPUT_FILE>] [--directory=<DATA_DIR>] --etype=<OBSTYPE> --qualifiers=<QUALIFIERS>\n\n"
					" Example: ./operaReductionSet --directory=/Users/espadons/data --etype=ALIGN --qualifiers=\"EEV1 pol Fast\"\n"
					"\n"
					"  -h, --help  display help message\n"
					"  -v, --verbose,  Turn on message sending\n"
					"  -d, --debug,  Turn on debug messages\n"
					"  -t, --trace,  Turn on trace messages\n"
					"\n" 
					"  -i, --input, Input file with list of paths (optional - overrides --directory)\n"
					"  -o, --ouput, Output with list of selected files (print on stdout if ouptut not provided)\n" 
					"  -r, --directory, Input data directories from command line (overrides config file)\n"
					"  -q, --qualifier, qualifiers to select dataset\n" 
					"  -e, --etype, obstype to select dataset\n"
					"\n";
}

/*! 
 * operaReductionSet
 * \author Eder Martioli
 * \brief Creates reduction sets from a list of input FITS file names.
 * \arg argc
 * \arg argv
 * \note --directory=$(DATADIR)
 * \note --qualifiers="$(DETECTOR) $(MODE) $(SPEED)"
 * \note --etype=FLAT
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */

char *filepath = NULL;

int main(int argc, char *argv[])
{
	int opt,i,j;
	int debug=0, verbose=0, trace=0, plot=0;
	operaErrorCode errorcode = operaErrorCodeOK;;
	string commandbuffer = "";
	string modifiers = "";
	string splitkey = "";
	FILE *fout = NULL;
	
	char *input = NULL; // input list of file paths
	char *output = NULL; // output list of files
	
	char *clquals = NULL; // specify qualifiers
	
	char *etype[MAXCONFIGVALUES]; // obstype of images to be listed
	int ne=0;
	
	char *dirs[MAXNDIRS];
	int ni=0;
	try {
		int ndirs = 0;
		struct option longopts[] = {       
			{"input",1, NULL, 'i'}, 
			{"output",1, NULL, 'o'}, 		
			{"directory",1, NULL, 'r'}, 
			{"qualifiers",1, NULL, 'q'}, 
			{"etype",1, NULL, 'e'},
			{"splitkey",1,NULL,'s'},
			{"path",1,NULL,'P'},
			
			{"plot",		optional_argument, NULL, 'p'},       
			{"verbose",		optional_argument, NULL, 'v'},
			{"debug",		optional_argument, NULL, 'd'},
			{"trace",		optional_argument, NULL, 't'},
			{"help",		no_argument, NULL, 'h'},
			{0,0,0,0}};
		
		i = 1;
		
		while((opt = getopt_long(argc, argv, "i:o:r:q:e:s:P:v::d::t::p::h", 
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
					clquals = optarg;
					break;
				case 's':
					splitkey = optarg;
					break;            
				case 'e':
					etype[ne++] = optarg;
					break;   
				case 'P':
					filepath = optarg;
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
					printUsageSyntax();
					exit(EXIT_FAILURE);
					break;
				case '?':
					printUsageSyntax();
					exit(EXIT_FAILURE);
					break;
			}
		}	

		if (ni != 0 && input == NULL) {
			ndirs = ni;
		}
		if (input != NULL) {
			ndirs = 0;
		}
		
		if (argc < 3 || ne == 0 || clquals == NULL) {      
			printUsageSyntax();
			return(EXIT_FAILURE);
		} 
#ifdef PRINT_DEBUG    
		
		if (debug) {
			for (i=0;i<argc;i++) {
				if (i==0 || argv[i][0] == '-') {
					cout << argv[i] << " ";
				} else {
					cout << "\"" << argv[i] << "\" ";
				}
			} 
			cout <<"\n"; 
		} 

		if (debug) {  
			cout << 
					"   ....................................................... \n"	
					"   Running module: " << argv[0] << "\n" 
					"   Data directories:\n";
			for (i=0;i<ndirs;i++) {
				cout << "      " << dirs[i] << "\n";
			} 
			cout << "   Observation types:\n      ";		
			for (i=0;i<ne;i++) {
				cout << etype[i] << " ";
			}  		
			cout << "\n";		
			if (clquals != NULL)
				cout << "   Qualifiers: " << clquals << "\n";      
			cout << "   ....................................................... \n";	    
		}
#endif
		//----- handle etype options -----------  
		char *etypeopts_tmp[MAXCONFIGVALUES],*etypeopts[MAXCONFIGVALUES];
		char *etypevalues[MAXCONFIGVALUES];
		int netypeopts=0,etype_ok[MAXCONFIGVALUES];
		
		// Read ETYPE options listed in the config file (load etypeopts_tmp[]) 
		char etypes_config[MAX_ETYPES_CONFIG];
		snprintf(etypes_config,sizeof(etypes_config),"ETYPES");
		netypeopts = operaReductionSetConfigurationAccess(etypes_config, etypeopts_tmp, MAXCONFIGVALUES, &errorcode);
	
		if (errorcode || netypeopts == 0) {
			throw operaException("operaReductionSet: ", operaErrorReductionSetEtypeNotDefined, __FILE__, __FUNCTION__, __LINE__);	
		}	
		
		for (int ii=0;ii<ne;ii++)	{
			etype_ok[ii] = 1;
			
			for (i=0;i<netypeopts;i++){
				etypeopts[i] = (char *)malloc(strlen(etypeopts_tmp[i])+13);
				if (!etypeopts[i]) {
					throw operaException("operaReductionSet: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
				snprintf(etypeopts[i],strlen(etypeopts_tmp[i])+13,"%s_HEADERVALUE",etypeopts_tmp[i]);
				if (strncmp(etype[ii], etypeopts_tmp[i],sizeof(etype[ii])) == 0) {  	
					errorcode = operaConfigurationAccessGet(etypeopts[i], &etypevalues[ii]);
					if (errorcode) {
						throw operaException("operaReductionSet: ", operaErrorReductionSetEtypeFailed, __FILE__, __FUNCTION__, __LINE__);	
					}	
					etype_ok[ii] = 0;
					break;
				}
			}
			if (etype_ok[ii]) {	
				throw operaException("operaReductionSet: ", operaErrorReductionSetEtypeNotDefined, __FILE__, __FUNCTION__, __LINE__);	
			}		
		}
		
		// get obstypekey
		char *obstypekey;
		errorcode = operaConfigurationAccessGet("OBSTYPE_HEADERKEY", &obstypekey);
		if (errorcode) {
			throw operaException("operaReductionSet: ", operaErrorReductionSetObstypeKeyNotDefined, __FILE__, __FUNCTION__, __LINE__);	
		}		
		//-------------------------------------  
		
		
		//----- handle qualifier option -----------  
		char *qualiopts[MAXCONFIGVALUES],*qualiopts_default[MAXCONFIGVALUES];
		int nqualiopts=0; 
		
		// read default qualifier options from configuration file.
		char qualifiernamelist_config[18];
		snprintf(qualifiernamelist_config,sizeof(qualifiernamelist_config),"QUALIFIERNAMELIST");
		nqualiopts = operaReductionSetConfigurationAccess(qualifiernamelist_config, qualiopts_default, MAXCONFIGVALUES, &errorcode);	
		if (errorcode) {
			throw operaException("operaReductionSet: ", errorcode, __FILE__, __FUNCTION__, __LINE__);	
		}	
		
	  // search for command line options; keep default if no input
		if (clquals == NULL) {
			for (i=0;i<nqualiopts;i++)
				qualiopts[i] = qualiopts_default[i];
		} else {
			nqualiopts = SplitValues(clquals, qualiopts, MAXCONFIGVALUES);  
		}

		// search keyword definitions for each qualifier option
		char *qualioptkeys[MAXCONFIGVALUES],*qualiopt_headerkey[MAXCONFIGVALUES];
		for (i=0;i<nqualiopts;i++) {
			qualiopt_headerkey[i] = (char *)malloc(strlen(qualiopts[i])+11);
			if (!qualiopt_headerkey[i]) {
				throw operaException("operaReductionSet: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			snprintf(qualiopt_headerkey[i],strlen(qualiopts[i])+11,"%s_HEADERKEY",qualiopts[i]);
			errorcode = operaConfigurationAccessGet(qualiopt_headerkey[i], &qualioptkeys[i]);
			if (errorcode) {
				throw operaException("operaReductionSet: ", operaErrorReductionSetQualiKeyNotDefined, __FILE__, __FUNCTION__, __LINE__);	
			}	
		}
		
		// search default selected choice for each qualifier option	
		char *qualioptchoice_default[MAXCONFIGVALUES];
		for (i=0;i<nqualiopts;i++) {		
			errorcode = operaConfigurationAccessGet(qualiopts[i], &qualioptchoice_default[i]);	
			if (qualioptchoice_default[i] == NULL) {
				if (verbose)
					cout << "operaReductionSet: Warning: could not find default value for qualifier: "<< qualiopts[i] << "\n";
			}	
		}
		
		// parse the command line arguments to search for choices for each qualifier option 
		regex_t regex;
		char *qualioptchoice[MAXCONFIGVALUES];
		char namebuff[MAXCONFIGURATIONVALUELENGTH];	
		
		for (j=0;j<nqualiopts;j++) {
			qualioptchoice[j] = NULL;
			snprintf(namebuff, sizeof(namebuff), "%s=[[:print:]]", qualiopts[j]);
			regcomp(&regex, namebuff, REG_ICASE);		
			for (i=1;i<argc;i++) {		
				if (!regexec(&regex, argv[i], 0, NULL, 0)) {		// we got the name
					qualioptchoice[j] = (char*)malloc(MAXCONFIGURATIONVALUELENGTH);
					if (!qualioptchoice[j]) {
						throw operaException("operaReductionSet: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
					}
					char *start = strstr(argv[i], "=")+1;							// after the "="
					strncpy(qualioptchoice[j], start, MAXCONFIGURATIONVALUELENGTH); 
					break;
				}
			}
			if (qualioptchoice[j] == NULL && qualioptchoice_default[j] != NULL) {
				qualioptchoice[j] = qualioptchoice_default[j];
			}
		}

		// read header values for each qualifier option choice
		char *qualioptchoice_headervalue[MAXCONFIGVALUES],*qualioptchoice_value[MAXCONFIGVALUES];	

		for (i=0;i<nqualiopts;i++) {
			qualioptchoice_headervalue[i] = (char *)malloc(strlen(qualioptchoice[i])+13);
			if (!qualioptchoice_headervalue[i]) {
				throw operaException("operaReductionSet: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			snprintf(qualioptchoice_headervalue[i],strlen(qualioptchoice[i])+13,"%s_HEADERVALUE",qualioptchoice[i]);		
			errorcode = operaConfigurationAccessGet(qualioptchoice_headervalue[i],&qualioptchoice_value[i]);
			if (errorcode) {
				throw operaException("operaReductionSet: ", operaErrorReductionSetQualiValNotDefined, __FILE__, __FUNCTION__, __LINE__);	
			}		
		}	 
#ifdef PRINT_DEBUG    
		if (debug) {
			for (i=0;i<nqualiopts;i++) {
				cout << "operaReductionSet: \n"
					<< "QUALIFIER # " << i << ": " << qualiopts[i] << "\n"
					<< "KEYWORD:         " << qualioptkeys[i] << "\n"
					<< "DEFAULT OPTION:  " << qualioptchoice_default[i] << "\n"
					<< "OPTION SELECTED: " << qualioptchoice[i] << "\n"
					<< "HEADER VALUE:    " << (qualioptchoice_value[i]==NULL?"NULL":qualioptchoice_value[i]) << "\n";
			}
		}
		//------------------------------------------
		// Open input directories and get file names
		//------------------------------------------
		if (debug && ndirs)
			cout << "operaReductionSet: Opening input directory and reading file names...\n";
		if (debug && input)
			cout << "operaReductionSet: Opening input file list and reading file names... " << input << "\n";
#endif		
		struct dirent *entry;
		DIR *dp[MAXNDIRS];
		
		for (i=0;i<ndirs;i++) {
			dp[i] = opendir(dirs[i]);
			if (dp[i] == NULL) {
				throw operaException("operaReductionSet: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
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
					snprintf(fullpathname[j],strlen(entry->d_name)+strlen(dirs[i])+2,"%s/%s",dirs[i],entry->d_name); 
					j++;
				}
				closedir(dp[i]);
			}		
		} else {  //get file paths from input file	
			if (access(input,R_OK)) {
				throw operaException("operaReductionSet: ", operaErrorReductionSetInputNotFound, __FILE__, __FUNCTION__, __LINE__);	
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
#ifdef PRINT_DEBUG    		
		if (debug) {
			for (i=0;i<nfiles;i++)
				cout << "operaReductionSet: Trying to open file "<<fullpathname[i]<<"\n";  
		}
#endif
		char qualival_from_FITS[MAXCONFIGVALUES][FLEN_VALUE],etypeval_from_FITS[MAXCONFIGVALUES][FLEN_VALUE];
		int FitStatus = 0, accept_quali, accept_etype;
		int notfits[MAXNFILES];
		char comment[FLEN_COMMENT];  
		fitsfile *fptr;

		if (output != NULL) {
			fout = fopen(output,"w");
		}
#ifdef PRINT_DEBUG    
		if (debug)
			cout << "operaReductionSet: Parsing files...\n";
#endif
		unsigned char split = 0;			// set to 1 when we find this split key value in the stack
		unsigned modechangecount = 0;
		char currentsplitkeyvalue[FLEN_VALUE] = {'\0'};
		char splitkeyvaluestack[MAX_MODE_CHANGES][FLEN_VALUE] = 
			{{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},
				{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},
				{'\0'},{'\0'},{'\0'},{'\0'},{'\0'},{'\0'}};

		for (int i=0;i<nfiles;i++) {
			FitStatus = 0;
			char asplitkeyvalue[FLEN_VALUE];
#ifdef PRINT_DEBUG    
			if (debug)
				cout << "operaReductionSet: Trying to open file: " << fullpathname[i] << "\n";
#endif			
			notfits[i] = 0; 
			if (fits_open_file(&fptr, fullpathname[i], READONLY, &FitStatus)) {
				if (debug && FitStatus)	
					operaPError("operaReductionSet", FitStatus); 
				notfits[i] = 1;
				continue;        
			}
#ifdef PRINT_DEBUG    
			if (debug)
				cout << "operaReductionSet: Open succeeded, now reading FITS file: " << fullpathname[i] << "\n";
#endif			
			accept_etype = 0;
			for (int ii=0;ii<ne;ii++) {
				FitStatus = 0;
				if (fits_read_keyword(fptr, obstypekey, etypeval_from_FITS[ii], comment, &FitStatus) ) {
					if (FitStatus)
						operaPError("operaReductionSet "+string(obstypekey), FitStatus); 	 
				}	
				// condition based on obstype: obstype1 || obstype2 || ..
#ifdef PRINT_DEBUG    
				if (debug)
					cout << "operaReductionSet: ETYPE: testing " << etypevalues[ii] << " to match " << etypeval_from_FITS[ii] << "\n";	  
#endif
				if (!strcmp(etypevalues[ii], etypeval_from_FITS[ii])) {
					accept_etype = 1;
#ifdef PRINT_DEBUG    
					if (debug)
						cout << "operaReductionSet: ETYPE: accepted " << etypevalues[ii] << " to match " << etypeval_from_FITS[ii] << "\n";	  
#endif
					break;
				}    
			}
			//
			// if we have a split key then watch for triggers
			//
			if (!splitkey.empty()) {
				if (fits_read_keyword(fptr, splitkey.c_str(), asplitkeyvalue, comment, &FitStatus) ) {
					if (FitStatus)
						operaPError("operaReductionSet: Could not find split key "+splitkey+" in header ", FitStatus); 	 
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
			
			accept_quali = 1;
			for (int ii=0;ii<nqualiopts;ii++) 
				if (qualioptchoice_value[ii] != NULL)
				{
					qualival_from_FITS[ii][0] = '\0';
					FitStatus = 0;

					if (fits_read_keyword(fptr, qualioptkeys[ii], qualival_from_FITS[ii], comment, &FitStatus) ) {
						if (FitStatus)
							operaPError("operaReductionSet "+string(qualioptkeys[ii])+" ", FitStatus ); 	 
					}

					// condition based on qualifiers: quali1 && quali2 && quali3 &&..
#ifdef PRINT_DEBUG    
					if (debug)
						cout << "operaReductionSet: QUALIFIER: " << qualioptkeys[ii] << " testing " << qualioptchoice_value[ii] << " to match " << qualival_from_FITS[ii] << "\n";
#endif
					if (qualival_from_FITS[ii][0] != '\0' && qualioptchoice_value[ii][0] != '\0') {
						//cout << ":::" << qualioptchoice_value[ii] << ": :" << qualival_from_FITS[ii] << ":\n";
						
						// handle wildcards - a qualifier ending in %, also get rid of the '
						if (qualioptchoice_value[ii][strlen(qualioptchoice_value[ii])-2] == '%') {
							char temp[FLEN_VALUE];
							strncpy(temp,qualioptchoice_value[ii],FLEN_VALUE);
							temp[strlen(temp)-2] = '\0';	// remove the %
							if (startsWith(qualival_from_FITS[ii], temp) > 0) {
#ifdef PRINT_DEBUG    
								if (debug)
									cout << "operaReductionSet: qualifier: " << qualioptchoice_value[ii] << " accepted == " << qualival_from_FITS[ii] << " temp=" << temp << "\n";	  
#endif
							} else {
#ifdef PRINT_DEBUG    
								if (debug)
									cout << "operaReductionSet: qualifier: " << qualioptchoice_value[ii] << " NOT accepted  != " << qualival_from_FITS[ii] << "\n";	  
#endif
								accept_quali = 0;
								break;
							}
						} else {	// exact match case
							if (!strcmp(qualioptchoice_value[ii], qualival_from_FITS[ii])) {
#ifdef PRINT_DEBUG    
								if (debug)
									cout << "operaReductionSet: qualifier: " << qualioptchoice_value[ii] << " accepted == " << qualival_from_FITS[ii] << "\n";	  
#endif
							} else {
#ifdef PRINT_DEBUG    
								if (debug)
									cout << "operaReductionSet: qualifier: " << qualioptchoice_value[ii] << " NOT accepted  != " << qualival_from_FITS[ii] << "\n";	  
#endif
								accept_quali = 0;
								break;
							}
						}
						
					} else {
#ifdef PRINT_DEBUG    
						if (debug)
							cout << "operaReductionSet: Warning: could not find header value for qualifier "<< qualiopts[ii] << " in file "<<filename[i]<<"\n";
#endif
					}
				}
			if (accept_quali==1 && accept_etype==1) {
				if (output == NULL) {
					if (split) {
						split = 0;
						cout << "##########\n";				
					}
					cout << fullpathname[i] << "\n";				
				} else {
					if (split) {
						split = 0;
						fprintf(fout,"##########\n");
					}
					fprintf(fout,"%s\n",fullpathname[i]);
				}
			}	
			
			if (fits_close_file(fptr, &FitStatus)) {
#ifdef PRINT_DEBUG    
				if (debug)
					operaPError("operaReductionSet", FitStatus );
#endif
			}	
		} // end of for "nfiles" loop
		if (fout) {
			fclose(fout);
		}
	} catch (operaErrorCode errorcode) {
		if (fout) {
			fclose(fout);
		}
		operaPError("operaReductionSet", errorcode );
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
				values[spc] = (char *)malloc(strlen(parvalues));
				if (!values[spc]) {
					throw operaException("operaReductionSet: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
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

unsigned operaReductionSetConfigurationAccess(char *par, char *values[], unsigned maxvalues, operaErrorCode *errorcode) {
	
	unsigned npars = 0;
	char *parvalues = NULL;
	
	operaConfigurationAccessSetConfigurationFilepath(filepath);
	*errorcode = operaConfigurationAccessGet(par, &parvalues);
	
	if (*errorcode == operaErrorCodeOK) {	
		npars = SplitValues(parvalues, values, maxvalues);
	}
	free(parvalues);
	return npars;
}

