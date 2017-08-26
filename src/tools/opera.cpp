/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: opera
 Version: 1.0
 Description: Main driver for opera.
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaConfigurationAccess.h"
#include "libraries/operaLibCommon.h" // for startsWith
#include "libraries/operaLib.h"     // for itos
#include "tools/opera.h"

/*! \brief Main driver for opera. */
/*! \file opera.cpp */

using namespace std;

/*! 
 *  opera
 * \author Doug Teeple
 * \brief Main driver for opera.
 * \arg argc
 * \arg argv
 * \note -v --verbose
 * \note -d --debug
 * \note -t --trace
 * \note -p --plot
 * \note -n --noexecute
 * \note -a --valgrind
 * \note -h --help
 * \note For now it uses getenv, but will eventually use getconfiguration
 * \return EXIT_STATUS
 * \ingroup tools
 ********************************************************************/
int debug=0, verbose=0, plot=0, trace=0, execute=1, valgrind=0;

int main(int argc, char *argv[])
{
	char *prefix = NULL;
	string commandbuffer = "";
	string optargs = "";
	string sandbox = "opera-1.0/";
	int debugvalue=0, verbosevalue=0, plotvalue=0, tracevalue=0;
	
	try {
		prefix = getenv("opera");
		//operaErrorCode operaError = operaConfigurationAccessGet("prefix-dir", &prefix);
		//if (operaError) {
		//	operaPError("opera", operaError);
		//}
		if (prefix == NULL || (prefix != NULL && strlen(prefix) == 0)) {
			char *home = getenv("HOME");
			sandbox = string(home)+"/"+sandbox;
			if (directoryexists(sandbox)) {
				prefix = (char *)sandbox.c_str();
			}
		}
		if (prefix != NULL && strlen(prefix) > 0) {
			commandbuffer = "export opera="+string(prefix)+" ; make -f "+string(prefix)+"/harness/Makefile prefix="+string(prefix)+" --no-print-directory --jobs ";
		} else {
			cerr << "Please specify the environment variable \"opera\" pointing to the location of your OPERA installation.\n";
			cerr << "e.g. export opera=$HOME/opera-1.0\n";
			return EXIT_FAILURE;
		}
		for (int i=1; i<argc; i++) {
			if (startsWith(argv[i], "-v") || startsWith(argv[i], "--verbose")) {
				verbose = 1;
				if (sscanf(argv[i], "-v%d", &verbosevalue) == 1 && verbosevalue != 0)
					verbose = verbosevalue;
				if (sscanf(argv[i], "--verbose=%d", &verbosevalue) == 1 && verbosevalue != 0)
					verbose = verbosevalue;
				optargs += "--verbose"+(verbose>1?itos(verbose):"")+" ";
			} else if (startsWith(argv[i], "-d") || startsWith(argv[i], "--debug")) {
				debug = 1;
				if (sscanf(argv[i], "-d%d", &debugvalue) && debugvalue != 0)
					debug = debugvalue;
				if (sscanf(argv[i], "--debug=%d", &debugvalue) && debugvalue != 0)
					debug = debugvalue;
				optargs += "--debug"+(debug>1?itos(debug):"")+" ";
			} else if (startsWith(argv[i], "-t") || startsWith(argv[i], "--trace")) {
				trace = 1;
				if (sscanf(argv[i], "-t%d", &tracevalue) && tracevalue != 0)
					trace = tracevalue;
				if (sscanf(argv[i], "--trace=%d", &tracevalue) && tracevalue != 0)
					trace = tracevalue;
				optargs += "--trace"+(trace>1?itos(trace):"")+" ";
			} else if (startsWith(argv[i], "-p") || startsWith(argv[i], "--plot")) {
				plot = 1;
				if (sscanf(argv[i], "-p%d", &plotvalue) && plotvalue != 0)
					plot = plotvalue;
				if (sscanf(argv[i], "--plot=%d", &plotvalue) && plotvalue != 0)
					plot = plotvalue;
				optargs += "--plot"+(plot>1?itos(plot):"")+" ";
			} else if (startsWith(argv[i], "-a") || startsWith(argv[i], "--valgrind")) {
				valgrind = 1;
			} else if (startsWith(argv[i], "-n") || startsWith(argv[i], "--noexecute")) {
				trace = 1;
				execute = 0;
			} else if (startsWith(argv[i], "-h") || startsWith(argv[i], "--help")) {
				printUsageSyntax();
				exit(EXIT_SUCCESS);
			} else {
				int skip = 0;
				if (argv[i][skip] == '-') skip++;
				if (argv[i][skip] == '-') skip++;
				string arg = (skip ? string(&argv[i][skip]) : string(argv[i]));
				if (arg.find(' ') != string::npos)
					commandbuffer += arg.substr(0,arg.find('=')+1) + '"' + arg.substr(arg.find('=')+1) + '"'  + ' ';
				else commandbuffer += arg + ' ';
			}

		}
		
		if (verbose) {
			commandbuffer += "VERBOSE="+itos(verbose)+" ";
		} else {
			commandbuffer += "VERBOSE=0 ";
		}
		if (debug) {
			commandbuffer += "DEBUG="+itos(debug)+" ";
		} else {
			commandbuffer += "DEBUG=0 ";
		}
		if (plot) {
			commandbuffer += "PLOT="+itos(plot)+" ";
		} else {
			commandbuffer += "PLOT=0 ";
		}
		if (trace) {
			if (execute) {
				if (valgrind) {
					commandbuffer += "TRACE=4 ";
				} else {
					if (debug) {
						commandbuffer += "TRACE=3 ";
					} else {
						commandbuffer += "TRACE=2 ";
					}
				}
			} else {
				commandbuffer += "TRACE=1 ";
			}
		} else {
			commandbuffer += "TRACE=0 ";
		}
		commandbuffer += "optargs=\""+optargs+"\"";
		if (trace) {
			cout << commandbuffer << endl;
		}
		systemf(commandbuffer.c_str());
		
	}
	catch (operaException e) {
		cerr << "opera: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "opera: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n" <<
	" Usage: opera <options> <commands>\n" <<
	"\n" <<
	" Basic Commands:\n" <<
	"     opera --help\n" <<
	"     opera --version\n" <<
	"\n" <<
	" Options:\n" <<
	"     opera -v[n] | --verbose[=n]	-- verbose output\n" <<
	"     opera -t[n] | --trace[=n]		-- trace execution steps\n" <<
	"     opera -d[n] | --debug[=n]		-- debug output and a debug execution trace if -v is also on\n" <<
	"     opera -p[n] | --plot[=n]		-- generate plots (some modules also take --interactive -I)\n" <<
	"     opera -n    | --noexecute		-- trace only, do not execute instructions (generates a script)\n" <<
	"     in all cases the optional number \"n\" is a modifier relevant to a particular module\n" <<
	"     in the case of -p, for example, the modifier would naturally be the order to plot\n" <<
	"     in the case of -v - or -d, or -t, the modifier may indicate the level of verbosity, debug output, and trace resp.\n" <<
	"\n" <<
	"*********************************************\n" <<
	" Calibration Commands:\n" <<
	"     opera <specifier> calibrations		-- parallel calibrations\n" <<
	(verbose?"     opera <specifier> calibrate calibration mastercalibration		-- serial calibrations\n":"") <<
	(verbose?"     opera <specifier> masterimages			-- image calibration step only (masterbias, masterflat, ...)\n":"") <<
	(verbose?"     opera <specifier> reductionset			-- create a reduction set\n":"") <<
	(verbose?"     opera <specifier> masterflats [--pick=0|1|2]		-- create master flat images\n":"") <<
	(verbose?"     opera <specifier> masterbiases [--pick=0|1|2]	-- create master bias images\n":"") <<
	(verbose?"     opera <specifier> mastercomparisons [--pick=0|1|2]	-- create master comparison images\n":"") <<
	(verbose?"     opera <specifier> masterfabryperots [--pick=0|1|2]	-- create master Fabry-Perot images\n":"") <<
	(verbose?"     opera <specifier> normalizedflats		-- create master flat images\n":"") <<
	(verbose?"     opera <specifier> gains				-- gain calculations\n":"") <<
	(verbose?"     opera <specifier> biases				-- bias level calculations\n":"") <<
	(verbose?"     opera <specifier> geoms | geometries			-- geometry calibrations\n":"") <<
	(verbose?"     opera <specifier> wcals | wavelengthcalibrations	-- wavelength calibrations\n":"") <<
	(verbose?"     opera <specifier> profs | profiles			-- instrument profiles\n":"") <<
	(verbose?"     opera <specifier> lines | spectrallines		-- spectral line calibrations\n":"") <<
	(verbose?"     opera <specifier> apers | apertures			-- aperture calibrations\n":"") <<
	(verbose?"     opera <specifier> fcals | fluxcalibrations		-- flux calibrations\n":"") <<
	(verbose?"     opera <specifier> disps | dispersioncalibrations	-- dispersion polynomial calibrations\n":"") <<
	(verbose?"     opera <specifier> badpixelmask | badpix			-- create a badpixel FITS mask from a map\n":"") <<
	(verbose?"     opera <specifier> eveningcalibrations | afternooncalibrations	-- perform limited calibrations for quicklook (up to geom)\n":"") <<
	"*********************************************\n" <<
	" Reduction Commands:\n" <<
	"     opera <specifier> reduce			-- parallel intensity and polarimetry reduction\n" <<
	(verbose?"     opera <specifier> intensity		-- parallel intensity reduction\n" :"") <<
	(verbose?"     opera <specifier> polarimetry	-- parallel polarimetry reduction\n":"")  <<
	(verbose?"*********************************************\n":"") <<
	(verbose?" Distribution Commands:\n":"") <<
	(verbose?"     opera <specifier> distribute		-- separate images by RUNID for distribution\n":"") <<
	(verbose?"     opera <specifier> approve			-- approve distribution for PIs download\n":"") <<
	"*********************************************\n" <<
	" Analysis Commands:\n" <<
	"     opera <specifier> polarbinning			-- perform binning on polar data\n" <<
	"     opera <specifier> requestedsnr			-- get the SNR value requested by the PI (CFHT-only)\n" <<
	"     opera DATADIRS=... [TEMPLATE=...] [OUT=...] OBJECT=... radialvelocitysequence -- Perform a series of radial velocity calculations over several datasets\n" <<
	"*********************************************\n" <<
	" Utility Commands:\n" <<
	"     opera <specifier> unlock					-- unlock a NIGHT/DATADIR directory\n" <<
	"     opera <specifier> clean					-- cleans all\n" <<
	(verbose?"     opera <specifier> cleanreduction				-- cleans intensity and polarimetry\n":"") <<
	(verbose?"     opera <specifier> cleanpol[arimetry]			-- cleans polarimetry\n":"") <<
	(verbose?"     opera <specifier> cleanint[ensity]			-- cleans intensity\n":"") <<
	(verbose?"     opera <specifier> cleanallintensity			-- also cleans the .e files\n":"") <<
	(verbose?"     opera <specifier> cleancal[ibrations]		-- clean all calibrations\n":"") <<
	(verbose?"     opera <specifier> clean[prof|gain|geom|disp|wcal|aper|prof|rvel|fcal|rlist|tell]\n":"") <<
	(verbose?"     opera <specifier> cleanvisuals\n":"") <<
	(verbose?"     opera <specifier> cleanproducts				-- clean m.fits, p.fits, i.fits\n":"") <<
	(verbose?"     opera <specifier> cleanbyproducts			-- clean reduction lists, visualization data files\n":"") <<
	(verbose?"     opera <specifier> cleanlogs\n":"") <<
	(verbose?"     opera <specifier> cleansnr\n":"") <<
	(verbose?"     opera <specifier> cleantmp\n":"") <<
	(verbose?"     opera <specifier> cleansemester SEMESTER=<>	-- clean an entire semester\n":"") <<
	(verbose?"     opera <specifier> cleanodometer ODO=<>		-- clean all products and byproducts associated with an odometer\n":"") <<
	(verbose?"     opera <specifier> cleanit WHAT=<extension>		-- i.e. geom gain prof eps s sn ...etc\n":"") <<
	"     opera <specifier> status				-- returns reduction status\n" <<
	"     opera <specifier> targets				-- lists targets and PI/RUNIDS\n" <<
	"     opera <specifier> configuration			-- returns the current configuration\n"
	"     opera getorder WL=<wavelength in nm>		-- returns the order(s) covering the given wavelength\n" <<
	"     opera version					-- returns the current opera version and subversion\n" <<
	(verbose?"*********************************************\n":"") <<
	(verbose?" Setup / Quicklook Commands:\n":"") <<
	(verbose?"     opera DATADIR=... setup					-- perform setup calibrations and reduction\n":"") <<
	(verbose?"     opera DATADIR=... cleansetup				-- cleanup setup calibrations and reduction\n":"") <<
	(verbose?"     opera <specifier> spectrum | quicklook FILE=<filename>	-- reduction of a single quicklook image\n":"") <<
	(verbose?"     opera <specifier> cleanquicklook			-- cleanup reduction of a single quicklook image\n":"") <<
	(verbose?"     opera telescopesetup						-- setup for a new telescope\n":"") <<
	(verbose?"     opera spectrographsetup					-- setup for a new spectrograph, for all modes. speeds\n":"") <<
	"*********************************************\n" <<
	"\n"
	"Where <specifier> is of the form: "
	"   NIGHT=QRUNID-MMMDD		-- e.g. 08AQ04-Mar16\n" <<
	"or DATADIR=<directory>		-- full directory path containing the data to reduce.\n" <<
	"or DATADIRS=\"<directory> <directory>...\"	-- full directory paths containing the data to reduce.\n" <<
	"\n" <<
	"WHERE=\"OBJECT=... PI_NAME=... RUNID=...\"	-- selects a subset of the objects.\n" <<
	"\n" <<
	(verbose?"Note: When issuing individual clean targets during setup, specify TYPE=setup.\n":"") <<
	(verbose?"e.g.: opera DATADIR=... cleanlogs TYPE=setup\n":"") <<
	//	"Setup commands:\n"
	//	"     opera DIR=<directory> -pcmwf[g] <command>\n"
	//	"Where DIR is that part of the directory not containing /data/niele/espadons/\n"
	//	"e.g. opera DIR=setup/Jul30 -pcmwf calibrate\n"
	//	"The data directory will be set to /data/$sessionhost/espadons/opera\n"
	"Specify -v option for more commands." << 
	endl;
}

