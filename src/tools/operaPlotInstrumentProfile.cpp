/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaPlotInstrumentProfile
 Version: 1.0
 Description: Tool to plot instrument profile from calibration
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
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

#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile

#include "libraries/operaLib.h"						// systemf
#include "libraries/operaMath.h"                    // for LengthofPolynomial

#include "libraries/operaLibCommon.h"

void GenerateInstrumentProfile3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned minorderWithIP, unsigned maxorderWithIP, unsigned IPxsize, unsigned IPysize, bool display);

static void printUsageSyntax(char * modulename);

#define NOTPROVIDED -999

/*! \file operaPlotInstrumentProfile.cpp */

using namespace std;

/*!
 * operaPlotInstrumentProfile
 * \author Eder Martioli
 * \brief Tool to plot the instrument profile.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup tools
 */

int main(int argc, char *argv[])
{
	int opt;

    unsigned minorder = 22;
    bool minorderprovided = false;
    unsigned maxorder = 62;    
    bool maxorderprovided = false;
    
	string inputgeom; 	
	string inputprof;  	
    
	int ordernumber = NOTPROVIDED;		
    unsigned pickImageRow = 0;
	int debug=0, verbose=0, trace=0, plot=0;

	bool interactive = false;
    string plotfilename;	
	string datafilename;	
	string scriptfilename;	
    
	struct option longopts[] = {      
		{"inputgeom",1, NULL, 'g'},
		{"inputprof",1, NULL, 'p'}, 
		{"pickImageRow",1, NULL, 'R'},        
		{"ordernumber",1, NULL, 'O'},
		{"minorder",1, NULL, 'M'},
		{"maxorder",1, NULL, 'X'},        
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},  
		{"interactive",0, NULL, 'I'},        
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "g:p:R:O:M:X:P:F:S:I:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{                            
			case 'g':
				inputgeom = optarg;
				break;				
			case 'p':
				inputprof = optarg;
				break;		
			case 'R':
				pickImageRow = atoi(optarg);
				break;                    
			case 'O':
				ordernumber = atoi(optarg);
				break;   
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;                 
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				datafilename = optarg;
				break; 	
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
				interactive = true;
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
	
	/*Start the module here*/
	
	try {
		// we need a profile file...
		if (inputprof.empty()) {
			throw operaException("operaPlotInstrumentProfile: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
	
		if (verbose) {
			cout << "operaPlotInstrumentProfile: inputgeom = " << inputgeom << endl; 
			cout << "operaPlotInstrumentProfile: inputprof = " << inputprof << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaPlotInstrumentProfile: ordernumber = " << ordernumber << endl;            
            }
            if(plot) {
                cout << "operaPlotInstrumentProfile: plotfilename = " << plotfilename << endl;
                cout << "operaPlotInstrumentProfile: datafilename = " << datafilename << endl;
                cout << "operaPlotInstrumentProfile: scriptfilename = " << scriptfilename << endl;                
            }            
		}
  	
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }        
        
        operaSpectralOrderVector spectralOrders;
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputprof);
        
		if (!inputgeom.empty()) {
            operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgeom); 
            if(!minorderprovided) {
                minorder = spectralOrders.getMinorder();
            }
            if(!maxorderprovided) {
                maxorder = spectralOrders.getMaxorder();            
            }
		}
		
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}

        int minorderWithIP = NOTPROVIDED;
        int maxorderWithIP = NOTPROVIDED;
        
        unsigned IPxsize = 0;
        unsigned IPysize = 0;
        
 		for (unsigned order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if(spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                
                operaGeometry *Geometry = spectralOrder->getGeometry();
                operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                
                if(IPxsize==0 && IPysize==0) {
                    IPxsize = instrumentProfile->getxsize();
                    IPysize = instrumentProfile->getysize();
                }
                
                float ysample;
                
                if(pickImageRow) {
                    ysample = (float)pickImageRow;
                } else {
                    ysample = (float)(Geometry->getYmax() - Geometry->getYmin())/2;
                }
                
                float distd = (float)Geometry->CalculateDistance(Geometry->getYmin(), ysample);
                
                instrumentProfile->printModel(distd,order,fdata);
                
                if(minorderWithIP == NOTPROVIDED){
                    minorderWithIP = (int)order;
                }
                maxorderWithIP = (int)order;
            }
        }

        if (fdata != NULL) {
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateInstrumentProfile3DPlot(scriptfilename,plotfilename,datafilename,minorderWithIP,maxorderWithIP, IPxsize, IPysize, interactive);
            }
            
        }
	}
	catch (operaException e) {
		cerr << e.getFormattedMessage() << '\n';
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaPlotInstrumentProfile: " << s << '\n';
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPlotInstrumentProfile: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --inputgeom=<GEOM_FILE>"
	" --inputprof=<PROF_FILE>"
	" --pickImageRow=<UNS_VALUE>"
	" --ordernumber=<INT_VALUE>"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"       
	" --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>" 
	" --interactive=<BOOL>\n\n"       
	" Example: "+string(modulename)+" --inputprof=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.prof --inputgeom=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.geom -R 1000 --minorder=22 --maxorder=61  -P testipplot.eps -F testipplot.dat -S testipplot.gnu \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -g, --inputgeom=<GEOM_FILE>, Input geometry file\n"
	"  -p, --inputprof=<PROF_FILE>, Input instrument profile file\n"
	"  -R, --pickImageRow=<UNS_VALUE>, Pick row number to plot IP model (default = naxis2/2)\n"
	"  -O, --ordernumber=<INT_VALUE>, Pick order number to plot IP model (default = all)\n"
	"  -M, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"       
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateInstrumentProfile3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned minorderWithIP, unsigned maxorderWithIP, unsigned IPxsize, unsigned IPysize, bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "set view 0,0" << endl;
    
    *fgnu << "set palette color" << endl;
    *fgnu << "set palette gamma 2.5" << endl;
    *fgnu << "set pm3d map" << endl;
    *fgnu << "unset ztics" << endl;
    *fgnu << "set cblabel \"flux fraction\"" << endl;
    *fgnu << endl;
    *fgnu << "set xrange[-5:(5*" << IPxsize << "+5)]" << endl;
    *fgnu << "set yrange[-5:(("<< (maxorderWithIP - minorderWithIP) <<"*" << IPysize << "/5)+5)]" << endl;
    
    
    for(unsigned order=minorderWithIP;order<=maxorderWithIP;order++) {
        double xlabelpos = ((double)order-((double)minorderWithIP+(double)(5*floor((order-minorderWithIP)/5))))*(double)IPxsize+1;
        double ylabelpos = (double)floor((order-minorderWithIP)/5)*(double)IPysize+1;
        *fgnu << "set label \"" << order << "\" at " << xlabelpos << "," << ylabelpos << " front font \"Helvetica,7\"" << endl;
    }
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}
