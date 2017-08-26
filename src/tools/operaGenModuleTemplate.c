/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaGenModuleTemplate
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
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/*! \brief Generate a default module template.  */
/*! \file operaGenModuleTemplate.c */
/*! \package operaGenModuleTemplate */
/*! \author Eder Martioli */

/*! 
 * operaGenModuleTemplate
 * \brief Generate a default module template.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup tools
 */
void Usage();
void alluppercase(char *p);

int gen_cfile (char *modname);
int gen_hfile (char *modname);

int main(int argc, char *argv[])
{
  int status = 0;

  fprintf(stdout, "\n"
                  "***************************************\n"
                  "*** OPERA Module Template Generator ***\n"
                  "***************************************\n");

  if(argc != 2)
  {
    fprintf(stderr, "\n  Error: invalid number of arguments! \n");
    
    Usage();
	exit(1);
  }

  status = gen_cfile(argv[1]);
  status = gen_hfile(argv[1]);  
//  status = gen_corefile(argv[1]);
//  status = gen_makefile(argv[1]);

  fprintf(stdout, "\nThe following files have been created:\n"
                  "\t1. %s.c\n"
                  "\t2. %s.h\n"                
                  ,argv[1],argv[1]); 

  fprintf(stdout, "\nTo compile do the following: \n"
                    "\t gcc -o %s %s.c -lm -lcfitsio ... [-libraries]\n",argv[1],argv[1]);  
                  
  return status;
}
/*-----------------------------------------------------------------------------*/
void Usage()
{
  fprintf(stderr, "\nUsage: ./gen_template MODULE_NAME\n\n");
} 
/*-----------------------------------------------------------------------------*/
int gen_cfile (char *modname)
{
  FILE *fp;
  char fname[300];

  sprintf(fname,"%s.c",modname);  
  fp = fopen(fname,"w");


  fprintf(fp,"/*******************************************************************\n");
  fprintf(fp,"****                  MODULE FOR OPERA v1.0                     ****\n");
  fprintf(fp,"********************************************************************\n");
  fprintf(fp,"Module name: %s\n",modname);
  fprintf(fp,"Version: 1.0\n");
  fprintf(fp,"Description: This is a template for helping developers \n");
  fprintf(fp,"             to start up with an OPERA module.\n");
  fprintf(fp,"Author(s): CFHT OPERA team\n");
  fprintf(fp,"Affiliation: Canada France Hawaii Telescope \n");
  fprintf(fp,"Location: Hawaii USA\n");
  fprintf(fp,"Date: Jan/2011\n");
  fprintf(fp,"Contact: eder@cfht.hawaii.edu\n");
  fprintf(fp,"\n");
  fprintf(fp,"Copyright (Unpublished--all rights reserved under the copyright laws \n");
  fprintf(fp,"of the United States).\n");
  fprintf(fp,"\n");
  fprintf(fp,"Permission to freely use, copy, modify, and distribute this software\n");
  fprintf(fp,"and its documentation without fee is hereby granted, provided that \n");
  fprintf(fp,"this copyright notice and disclaimer of warranty appears in all \n");
  fprintf(fp,"copies.\n");
  fprintf(fp,"\n");
  fprintf(fp,"DISCLAIMER:\n");
  fprintf(fp,"\n");
  fprintf(fp,"THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,\n");
  fprintf(fp,"EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED \n");
  fprintf(fp,"TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, \n");
  fprintf(fp,"ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR\n");
  fprintf(fp,"PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE\n");
  fprintf(fp,"DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE\n");
  fprintf(fp,"SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL ANYBODY BE LIABLE\n");
  fprintf(fp,"FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, \n");
  fprintf(fp,"SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR\n");
  fprintf(fp,"IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON \n");
  fprintf(fp,"WARRANTY, CONTRACT, TORT , OR OTHERWISE, WHETHER OR NOT INJURY WAS \n");
  fprintf(fp,"SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT \n");
  fprintf(fp,"LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, \n");
  fprintf(fp,"THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.\"\n");
  fprintf(fp,"********************************************************************/\n");
  fprintf(fp,"// $Date$\n");
  fprintf(fp,"// $Id$\n");
  fprintf(fp,"// $Revision$\n");
  fprintf(fp,"// $Locker$\n");
  fprintf(fp,"// $Log$\n");
  fprintf(fp,"\n");
  fprintf(fp,"#include \"%s.h\"\n",modname);
  fprintf(fp,"\n");
  fprintf(fp,"int main(int argc, char *argv[])\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  int %s_status = 0;\n",modname);
  fprintf(fp,"  int opt,i;\n");
  fprintf(fp,"  int np=0,ni=0;\n"); 
  fprintf(fp,"  char *param[argc]; // parameter values\n");
  fprintf(fp,"  char *input[argc]; // Input files \n"); 
  fprintf(fp,"  char *product = NULL; // Output product \n");
  fprintf(fp,"  int debug=0, verbose=0, trace=0;\n");
  fprintf(fp,"\n");  
  fprintf(fp,"  struct option longopts[] = {\n");
  fprintf(fp,"      {\"param\",1, NULL, \'p\'},\n");
  fprintf(fp,"      {\"input\",1, NULL, \'i\'},\n");
  fprintf(fp,"      {\"product\",1, NULL, \'o\'},\n");   
  fprintf(fp,"      {\"verbose\",0, NULL, \'v\'},\n");
  fprintf(fp,"      {\"debug\",0, NULL, \'d\'},\n");
  fprintf(fp,"      {\"trace\",0, NULL, \'t\'},\n");    
  fprintf(fp,"      {\"help\",0, NULL, \'h\'},\n");
  fprintf(fp,"      {0,0,0,0}};\n");
  fprintf(fp,"\n");            
  fprintf(fp,"  while((opt = getopt_long(argc, argv, \"p:i:o:vdth\", \n");
  fprintf(fp,"  		longopts, NULL))  != -1)\n");
  fprintf(fp,"  {\n");
  fprintf(fp,"    switch(opt) \n");
  fprintf(fp,"    {\n");
  fprintf(fp,"      case \'p\':\n");
  fprintf(fp,"            param[np] = optarg;\n");
  fprintf(fp,"            np++;\n");  
  fprintf(fp,"            break;\n");
  fprintf(fp,"      case \'i\':\n");
  fprintf(fp,"            input[ni] = optarg;\n");
  fprintf(fp,"            ni++;\n");
  fprintf(fp,"            break;    \n");
  fprintf(fp,"      case \'o\':\n");
  fprintf(fp,"            product = optarg;\n");
  fprintf(fp,"            break;            \n");
  fprintf(fp,"      case \'v\':\n");
  fprintf(fp,"            verbose = 1;\n");
  fprintf(fp,"            break;\n");
  fprintf(fp,"      case \'d\':\n");
  fprintf(fp,"            debug = 1;\n");
  fprintf(fp,"            break;\n");
  fprintf(fp,"      case \'t\':\n");
  fprintf(fp,"            trace = 1;\n");
  fprintf(fp,"            break;         \n");   
  fprintf(fp,"      case \'h\':\n");
  fprintf(fp,"            printUsageSyntax(argv[0]);\n");
  fprintf(fp,"            exit(EXIT_FAILURE);\n");
  fprintf(fp,"            break;\n");
  fprintf(fp,"      case \'?\':\n");
  fprintf(fp,"            printUsageSyntax(argv[0]);\n");
  fprintf(fp,"            exit(EXIT_FAILURE);\n");
  fprintf(fp,"            break;\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"  }	\n");
  fprintf(fp,"\n");  
  fprintf(fp,"  if(trace) {\n");
  fprintf(fp,"    for(i=0;i<argc;i++) {\n");
  fprintf(fp,"      if(i==0 || argv[i][0] == '-') {\n");
  fprintf(fp,"        fprintf(stdout,\"%%s \",argv[i]);\n");
  fprintf(fp,"      } else {\n");
  fprintf(fp,"        fprintf(stdout,\"\\\"%%s\\\" \",argv[i]);\n");
  fprintf(fp,"      }\n");
  fprintf(fp,"    } \n");
  fprintf(fp,"    fprintf(stdout,\"\\n\"); \n");
  fprintf(fp,"  }   \n");
  fprintf(fp,"\n");
  fprintf(fp,"  if (argc < 4 || np == 0 || ni == 0 || product == NULL) {      \n");
  fprintf(fp,"    %s_status = 2;\n",modname);
  fprintf(fp,"    fprintf(stderr, \"\\nError: (%%s:%%s:%%d)\\n\", __FILE__, __func__, __LINE__); \n");
  fprintf(fp,"    printUsageSyntax(argv[0]);\n");
  fprintf(fp,"    print_%s_err(%s_status);\n",modname,modname);   
  fprintf(fp,"    exit(EXIT_FAILURE);\n");
  fprintf(fp,"  } \n");
  fprintf(fp,"\n");
  fprintf(fp,"  float par[np];\n");  
  fprintf(fp,"\n");  
  fprintf(fp,"  if(verbose) {    \n");
  fprintf(fp,"    fprintf(stdout,\"   ....................................................... \\n\");	\n");
  fprintf(fp,"    fprintf(stdout,\"   Running module: %%s\\n\", argv[0]);    \n");
  fprintf(fp,"    fprintf(stdout,\"\\n   Total number of parameters: %%d\\n\",np);\n");  
  fprintf(fp,"    for(i=0;i<np;i++) { \n");  
  fprintf(fp,"      par[i] = atof(param[i]);\n");  
  fprintf(fp,"      fprintf(stdout,\"   Value of param #%%d: %%f\\n\",i,par[i]);\n");
  fprintf(fp,"    }\n");   
  fprintf(fp,"    fprintf(stdout,\"\\n   Total number of inputs: %%d\\n\",ni);\n");
  fprintf(fp,"    for(i=0;i<ni;i++)\n");
  fprintf(fp,"      fprintf(stdout,\"   input #%%d: %%s\\n\",i,input[i]);  \n");
  fprintf(fp,"    fprintf(stdout,\"\\n   product: %%s\\n\", product);    \n");
  fprintf(fp,"    fprintf(stdout,\"   ....................................................... \\n\\n\"); \n");
  fprintf(fp,"\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"/*The module starts here*/\n\n");	  
  fprintf(fp,"  if (debug) { \n");  
  fprintf(fp,"    fprintf(stdout, \"\\n(%%s:%%s:%%d)\\n\",__FILE__, __func__, __LINE__); \n"); 
  fprintf(fp,"  }\n");  
  fprintf(fp,"  return 0;\n");
  fprintf(fp,"}  \n");
  fprintf(fp,"\n\n");
  fprintf(fp,"/*************** Module functions start here *****************/\n");  
  fprintf(fp,"\n");
  fprintf(fp,"/* Print out the proper program usage syntax */\n");
  fprintf(fp,"void printUsageSyntax(char *modulename) {\n");
  fprintf(fp,"\n");
  fprintf(fp,"   fprintf(stderr,\n");
  fprintf(fp,"      \"\\n\"\n");
  fprintf(fp,"      \" Usage: %%s  [-vdth] --param=<PARAMETER_VALUE_1> --param=<PARAMETER_VALUE_2> ... --product=<PRODUCT_FILE_NAME> --input=<INPUT_FILE_1> --input=<INPUT_FILE_2> ... \\n\\n\"\n");
  fprintf(fp,"      \" Example: %%s  -v -p 10 -p 2 --product=o.fits -i 001.fits -i 002.fits -i bad_pix.dat \\n\\n\"\n");
  fprintf(fp,"      \"  -h, --help  display help message\\n\"\n");
  fprintf(fp,"      \"  -v, --verbose,  Turn on message sending\\n\"\n");
  fprintf(fp,"      \"  -d, --debug,  Turn on debug messages\\n\"\n");
  fprintf(fp,"      \"  -t, --trace,  Turn on trace messages\\n\"\n");     
  fprintf(fp,"      \"  -p, --param=<PARAMETER_VALUE>, Input parameters  \\n\"\n");
  fprintf(fp,"      \"  -o, --product=<PRODUCT_FILE_NAME>, Output product file  \\n\"\n");
  fprintf(fp,"      \"  -i, --input=<INPUT_FILE_NAME>, Input files  \\n\\n\"\n");  
  fprintf(fp,"      ,modulename,modulename);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");
  fprintf(fp,"/* Print out %s error messages and exit program  */\n",modname);  
  fprintf(fp,"void print_%s_err(int %s_status)\n",modname,modname);
  fprintf(fp,"{\n");
  fprintf(fp,"    char reperr[31];\n");
  fprintf(fp,"    \n");
  fprintf(fp,"    if (%s_status)\n",modname);
  fprintf(fp,"    {\n");
  fprintf(fp,"       %s_report_error(%s_status,reperr); /* print error report */\n",modname,modname);
  fprintf(fp,"       fprintf(stderr,\"\\n Error:%s: %%s\\n\",reperr);\n",modname);
  fprintf(fp,"\n");
  fprintf(fp,"       exit(%s_status);    /* terminate the program, returning error status */\n",modname);
  fprintf(fp,"    }\n");
  fprintf(fp,"    return;\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");
  fprintf(fp,"/*\n");
  fprintf(fp,"  %s_report_error returns a short descriptive error message that corresponds\n",modname);
  fprintf(fp,"  to the input error status value.  The message may be up to 30 characters long. \n");
  fprintf(fp,"*/\n");
  fprintf(fp,"void %s_report_error(int status,     /* integer error status value */\n",modname);
  fprintf(fp,"            char *errtext)  /* O - error message (max 30 char long + null) */\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  errtext[0] = \'\\0\';\n");
  fprintf(fp,"\n");
  fprintf(fp,"  if (status >= 0)\n");
  fprintf(fp,"  {\n");
  fprintf(fp,"    switch (status) {\n");
  fprintf(fp,"    case 0:\n");
  fprintf(fp,"       strcpy(errtext, \"OK - no error\\n\");\n");
  fprintf(fp,"       break;\n");
  fprintf(fp,"    case 1:\n");
  fprintf(fp,"       strcpy(errtext, \"non-MODULE error\\n\");\n");
  fprintf(fp,"       break;\n");
  fprintf(fp,"    case 2:   \n");
  fprintf(fp,"       strcpy(errtext, \"Invalid set of arguments\\n\");\n");
  fprintf(fp,"       break;\n");
  fprintf(fp,"    case 3:   \n");
  fprintf(fp,"       strcpy(errtext, \"Invalid parameter value\\n\");\n");
  fprintf(fp,"       break;       \n");
  fprintf(fp,"    case 4:  \n"); 
  fprintf(fp,"       strcpy(errtext, \"Input file not found\\n\");\n");
  fprintf(fp,"       break; \n");        
  fprintf(fp,"    default:\n");
  fprintf(fp,"       strcpy(errtext, \"Unknown error status\\n\");\n");
  fprintf(fp,"       break;\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"  } else {\n");
  fprintf(fp,"     strcpy(errtext, \"Unknown error status\\n\");\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  return;\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");
  fprintf(fp,"/* Print out cfitsio error messages and exit program */\n");  
  fprintf(fp,"void PrintFitsError(int FitStatus)\n");
  fprintf(fp,"{\n");
  fprintf(fp,"    if (FitStatus)\n");
  fprintf(fp,"    {\n");
  fprintf(fp,"       fits_report_error(stderr, FitStatus); /* print error report */\n");
  fprintf(fp,"\n");
  fprintf(fp,"       exit(FitStatus); /* terminate the program, returning CFITSIO error status */\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"    return;\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");  
  fclose(fp);
  return 0;
}
/*-----------------------------------------------------------------------------*/
int gen_hfile (char *modname)
{
  FILE *fp;
  char fname[300], madname_uc[300];
  
  sprintf(fname,"%s.h",modname); 
  fp = fopen(fname,"w");
  
  sprintf(madname_uc,"%s",modname); 
  alluppercase(madname_uc);
 
  fprintf(fp,"#ifndef %s_H\n",madname_uc);
  fprintf(fp,"#define %s_H\n",madname_uc);
  fprintf(fp,"#endif\n");  
  fprintf(fp,"/* libraries */\n");  
  fprintf(fp,"#include <stdio.h>\n");
  fprintf(fp,"#include <stdlib.h>\n");
  fprintf(fp,"#include <getopt.h>\n");    
  fprintf(fp,"#include <ctype.h>\n");
  fprintf(fp,"#include <string.h>\n");
  fprintf(fp,"#include <math.h>\n");  
  fprintf(fp,"#include <fitsio.h>\n");  
  fprintf(fp,"\n");    
  fprintf(fp,"/* prototypes */\n");
  fprintf(fp,"void printUsageSyntax(char *prgname);\n");
  fprintf(fp,"void print_%s_err(int %s_status);\n",modname,modname);
  fprintf(fp,"void %s_report_error(int status, char *errtext);\n",modname);
  fprintf(fp,"void PrintFitsError(int FitStatus);\n");

  fclose(fp);
  return 0;
}

void alluppercase(char *p)
{
  while( *p ) {
    *p = toupper(*p); // to upper case
    p++;    
  }
}
  