/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirDetrend
 Version: 1.0
 Description: Module to detrend a wircam image
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

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "operaError.h"
#include "core-wircam/wirDetrend.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"	
#include "libraries/operaWIRCamImage.h"	
#include "libraries/operaFITSSubImage.h"	

/* \file wirDetrend.cpp */
/* \package core_wircam */

using namespace std;

/* 
 * wirDetrend
 * \author Doug Teeple, Megan Tannock
 * \brief Module to detrend a wircam image. Detrending consist basically of
 * subtracting the dark and dividing by the flat. There are other nuances
 * such as optionally subtracting the reference pixels from the overscan
 * area of the chips. 
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */


int main(int argc, char *argv[])
{
	int opt;
	string name_raw; 
	string name_dark;
	string name_flat;
	string name_badpix;
	string name_detrended;
	string name_refpixels;
	string param_iiwiversion = "3.0";
	string param_procdate;
	
	bool linearize= false, 
	subtractreferencepixels = false, 
	guidewindowcrosstalk = false;
	
	float maskingthreshold = 0.0;
	
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	
	string nlcorrection_date;
	Polynomial chip1(3);
	Polynomial chip2(3);
	Polynomial chip3(3);
	Polynomial chip4(3);
	
	chip1.Set(1.0, 0);
	chip2.Set(1.0, 0);
	chip3.Set(1.0, 0);
	chip4.Set(1.0, 0);
    
	struct option longopts[] = {
		{"name_raw",1, NULL, 'R'},
		{"name_dark",1, NULL, 'D'},		
		{"name_flat",1, NULL, 'T'},	
		{"name_badpix",1, NULL, 'B'},	
		{"name_detrended",1, NULL, 'O'},	
		{"name_refpixels",1, NULL, 'i'},	
		{"param_procdate",1, NULL, 'A'},	
		{"param_iiwiversion",1, NULL, 'V'},	
		{"param_linearize",1, NULL, 'L'},	
		{"nlcorrection_date",1, NULL, 'c'},	
		{"nlcorrection_chip1_parg0",1, NULL, '1'},	
		{"nlcorrection_chip1_parg1",1, NULL, '2'},	
		{"nlcorrection_chip1_parg2",1, NULL, '3'},	
		{"nlcorrection_chip2_parg0",1, NULL, '4'},	
		{"nlcorrection_chip2_parg1",1, NULL, '5'},	
		{"nlcorrection_chip2_parg2",1, NULL, '6'},	
		{"nlcorrection_chip3_parg0",1, NULL, '7'},	
		{"nlcorrection_chip3_parg1",1, NULL, '8'},	
		{"nlcorrection_chip3_parg2",1, NULL, '9'},	
		{"nlcorrection_chip4_parg0",1, NULL, '0'},	
		{"nlcorrection_chip4_parg1",1, NULL, 'a'},	
		{"nlcorrection_chip4_parg2",1, NULL, 'b'},	
		{"param_subrefpix",1, NULL, 'E'},	
		{"param_gwinxtalk",1, NULL, 'X'},	
		{"param_maskingthreshold",1, NULL, 'M'},	
		
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},  
		{"interactive",0, NULL, 'I'},
		
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "R:D:T:B:O:i:A:V:L:c:E:X:M:P:F:S:I:1:2:3:4:5:6:7:8:9:0:a:b:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'R':
				name_raw = optarg;
				break;   
			case 'D':
				name_dark = optarg;
                break; 				
			case 'T':
				name_flat = optarg;
				break;
			case 'B':
				name_badpix = optarg;
				break;
			case 'i':
				name_refpixels = optarg;
				break;
			case 'O':
				name_detrended = optarg;
				break;
			case 'A':
				param_procdate = optarg;
				break;
			case 'V':
				param_iiwiversion = optarg;
				break;
			case 'L':
				linearize = atoi(optarg) == 1;
				break;
			case 'E':
				subtractreferencepixels = atoi(optarg) == 1;
				break;
			case 'X':
				guidewindowcrosstalk = atoi(optarg) == 1;
				break;
			case 'M':
				maskingthreshold = atof(optarg);
				break;
				
			case 'c':
				nlcorrection_date = optarg;
				break;
			case '1':
				chip1.Set(atof(optarg), 0);
				break;
			case '2':
				chip1.Set(atof(optarg), 1);
				break;
			case '3':
				chip1.Set(atof(optarg), 2);
				break;
				
			case '4':
				chip2.Set(atof(optarg), 0);
				break;
			case '5':
				chip2.Set(atof(optarg), 1);
				break;
			case '6':
				chip2.Set(atof(optarg), 2);
				break;
				
			case '7':
				chip3.Set(atof(optarg), 0);
				break;
			case '8':
				chip3.Set(atof(optarg), 1);
				break;
			case '9':
				chip3.Set(atof(optarg), 2);
				break;
				
			case '0':
				chip4.Set(atof(optarg), 0);
				break;
			case 'a':
				chip4.Set(atof(optarg), 1);
				break;
			case 'b':
				chip4.Set(atof(optarg), 2);
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
    
    // DETREND STEPS:
    // 1. reference pixels correction
    // 2. Dark subtraction
    // 3. Apply non-linearity correction for each chip
    // 4. Flat field correction
    // 5. Bad pixel mask
	
	try {
		if (name_raw.empty()) {
			throw operaException("wirDetrend: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_dark.empty()) {
			throw operaException("wirDetrend: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_flat.empty()) {
			throw operaException("wirDetrend: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_detrended.empty()) {
			throw operaException("wirDetrend: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cerr << "wirDetrend: input image = " << name_raw << endl; 
			cerr << "wirDetrend: detrended = " << name_detrended << endl;
			cerr << "wirDetrend: dark = " << name_dark << endl;							
			cerr << "wirDetrend: flat = " << name_flat << endl;
			cerr << "wirDetrend: refpixels = " << name_refpixels << endl;
			cerr << "wirDetrend: badpix = " << name_badpix << endl;
			cerr << "wirDetrend: version = " << param_iiwiversion << endl;
			cerr << "wirDetrend: proccessing date = " << param_procdate << endl;
			cerr << "wirDetrend: linearize = " << linearize << endl;
			cerr << "wirDetrend: subtractreferencepixels = " << subtractreferencepixels << endl;
			cerr << "wirDetrend: guidewindowcrosstalk = " << guidewindowcrosstalk << endl;
			cerr << "wirDetrend: maskingthreshold = " << maskingthreshold << endl;
			
            if (plot) {
                cerr << "wirDetrend: plotfilename = " << plotfilename << endl;
                cerr << "wirDetrend: datafilename = " << datafilename << endl;
                cerr << "wirDetrend: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cerr << "wirDetrend: interactive = YES" << endl; 
                } else {
                    cerr << "wirDetrend: interactive = NO" << endl; 
                }
            }            
		}
		
        ofstream *fdata = NULL;
        
        cout << "Reading files in to memory" << endl;
 		operaWIRCamImage in(name_raw, tfloat, READONLY);
		operaWIRCamImage out(name_detrended, in.getXDimension(), in.getYDimension(), in.getZDimension(), in.getNExtensions(), tfloat, READWRITE, cNone);
		operaWIRCamImage flat(name_flat, tfloat, READONLY);
		operaWIRCamImage dark(name_dark, tfloat, READONLY);
		operaWIRCamImage mask(name_badpix, tfloat, READONLY);
		operaWIRCamImage refpixels(name_refpixels, tfloat, READONLY);
		
#if 0         
        cout << "Masking GuideWindows" << endl;
		out.maskGuideWindows(out);
#endif
		
        // Copy headers to output files
		out.operaMultiExtensionFITSCubeCopyHeader(&in);
		refpixels.operaMultiExtensionFITSCubeCopyHeader(&in);
		
		// Copy raw pixels to out
        out = in;
		
        // 1. Reference pixel correction
		if (subtractreferencepixels) {
			cout << "Creating reference pixel image from overscan pixels " << endl;
			refpixels.createReferencePixelImage(in);
			cout << "Subtracting reference pixels from input image " << endl;
			out -= refpixels;
			refpixels.operaFITSSetHeaderValue("REFPXCOR", "yes", "Ref. pixel correction done?");
		}
        
        // 2. Dark subtraction
        cout << "Applying dark subtraction." << endl;
        out -= dark;
        out.operaFITSSetHeaderValue("DARKNAME", dark.operaFITSGetFilename(),	"Dark Name");
        out.operaFITSSetHeaderValue("DARKSUB", "yes", "Dark Subtraction done?");
		
        // 3. Apply non-linearity correction for each chip, polynomials declared above
        cout << "Applying non-linear correction " << endl;
		out.applyNonLinearCorrection(1, chip1);
		out.applyNonLinearCorrection(2, chip2);
		out.applyNonLinearCorrection(3, chip3);
		out.applyNonLinearCorrection(4, chip4);
        // Set header values
		out.operaFITSSetHeaderValue("NLCORR", "yes", "Non-linearity correction applied?");
		out.operaFITSSetHeaderValue("NLC_FUNC","xc/xm=a0+a1*xm+a2*xm^2","Non-linearity function");
		out.operaFITSSetHeaderValue("NLC_NAME", nlcorrection_date, "Non-linearity solution name");
		out.operaFITSSetHeaderValue("NLC_A0",chip1.Get(0),"Non-linearity function parameters", 1);
		out.operaFITSSetHeaderValue("NLC_A1",chip1.Get(1),"Non-linearity function parameters", 1);
		out.operaFITSSetHeaderValue("NLC_A2",chip1.Get(2),"Non-linearity function parameters", 1);
		out.operaFITSSetHeaderValue("NLC_A0",chip2.Get(0),"Non-linearity function parameters", 2);
		out.operaFITSSetHeaderValue("NLC_A1",chip2.Get(1),"Non-linearity function parameters", 2);
		out.operaFITSSetHeaderValue("NLC_A2",chip2.Get(2),"Non-linearity function parameters", 2);
		out.operaFITSSetHeaderValue("NLC_A0",chip3.Get(0),"Non-linearity function parameters", 3);
		out.operaFITSSetHeaderValue("NLC_A1",chip3.Get(1),"Non-linearity function parameters", 3);
		out.operaFITSSetHeaderValue("NLC_A2",chip3.Get(2),"Non-linearity function parameters", 3);
		out.operaFITSSetHeaderValue("NLC_A0",chip4.Get(0),"Non-linearity function parameters", 4);
		out.operaFITSSetHeaderValue("NLC_A1",chip4.Get(1),"Non-linearity function parameters", 4);
		out.operaFITSSetHeaderValue("NLC_A2",chip4.Get(2),"Non-linearity function parameters", 4);
		
		// 4. Flat field correction
        cout << "Applying flat field correction." << endl;
        out /= flat;
        out.operaFITSSetHeaderValue("FLATDIV", "yes", "Flat field division done?");
        
        // 5. Bad pixel mask
        cout << "Setting bad pixels in outImage to 0" << endl;
		out *= mask; // Set bad pixel values in outImage to zero
        // Update header
        out.operaFITSSetHeaderValue("BPIXNAME", mask.operaFITSGetFilename(), "Badpixelmask Name");
        out.operaFITSSetHeaderValue("BDPIXVAL", 0.0, "Bad pixel value");
		
        // Save detrended image
        cout << "Saving out" << name_detrended << endl;
 		out.operaWIRCamImageSave();
		
        // Close files
		in.operaFITSImageClose();
		dark.operaFITSImageClose();
		flat.operaFITSImageClose();
		mask.operaFITSImageClose();
		out.operaFITSImageClose();
		refpixels.operaFITSImageClose();
		
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }          
		
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateDetrendPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), interactive);
            }
        }                
	}
	catch (operaException e) {
		cerr << "wirDetrend: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirDetrend: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";		
}

void GenerateDetrendPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName,"w");
    
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"set view 0,0\n");
    fprintf(fgnu,"set iso 100\n");
    fprintf(fgnu,"set samples 100\n");
    fprintf(fgnu,"set pm3d at s\n");
    fprintf(fgnu,"set ticslevel 0\n");   
    
    fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
    fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
	
    if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
		fprintf(fgnu,"replot\n");       
		fclose(fgnu);   
		systemf("gnuplot -persist %s",gnuScriptFileName);
    } else {
		fclose(fgnu);  
	}
}
