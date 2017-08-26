/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaConfigurationAccess
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <getopt.h>
#include <stdio.h>

#include <fstream>
#include <map>

/*! \brief command line interface to the configuration file.*/
/*! \file operaConfigurationAccess.cpp */
/*! \author Doug Teeple */

using namespace std;

/*! 
 * operaConfigurationAccess
 * \brief command line interface to the configuration file.
 * \arg argc
 * \arg argv
 * \note [get|set|add] keyword[=value]
 * \return EXIT_STATUS
 * \ingroup tools
 */

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaConfigurationAccess.h"

#include "libraries/operaParameterAccess.h"
#include "libraries/operaConfigurationAccess.h"
#include "libraries/operaLib.h"

map<string, string> table;

static string strip(string stringIn, char whitespace)
{
	string::size_type pos = 0;
	bool spacesLeft = true;
	
	while( spacesLeft )
	{
		pos = stringIn.find(whitespace);
		if( pos != string::npos )
			stringIn.erase( pos, 1 );
		else
			spacesLeft = false;
	}
	
	return stringIn;
}

static string stripleft(string stringIn, char whitespace)
{
	string::size_type pos = 0;
	bool spacesLeft = true;
	
	while( spacesLeft )
	{
		if (stringIn[pos] == whitespace)
			stringIn.erase( pos, 1 );
		else
			spacesLeft = false;
	}
	
	return stringIn;
}

//
// add a line of the form foo=bar, foo+=bar foo:=bar
// replace if the var (foor) exists
//
static void tableAdd(string line) {
	size_t pos;
	string var, value;
	map<string,string>::iterator it;
	
	pos=line.find(":=");
	if (pos != string::npos) {
		var = line.substr(0,pos);
		var = strip(var, ' ');
		var = strip(var, '\t');
		value = line.substr(pos+2);
		value = stripleft(value, ' ');
		value = stripleft(value, '\t');
		it = table.find(var);
		if (it != table.end()) {
			table.erase(it);
		}
		table.insert(table.end(), pair<string, string>(var, value));
	} else {
		pos=line.find("+=");
		if (pos != string::npos) {
			var = line.substr(0,pos);
			var = strip(var, ' ');
			var = strip(var, '\t');
			value = line.substr(pos+2);
			value = stripleft(value, ' ');
			value = stripleft(value, '\t');
			it = table.find(var);
			if (it == table.end()) {
				// couldn't find it, just add it
				table.insert(table.end(), pair<string, string>(var, value));
			} else {
				value += " " + it->second;
				table.erase(it);
				table.insert(it, pair<string, string>(var, value));
			}
		} else {
			pos=line.find("-=");
			if (pos != string::npos) {
				var = line.substr(0,pos);
				var = strip(var, ' ');
				var = strip(var, '\t');
				value = line.substr(pos+2);
				value = stripleft(value, ' ');
				value = stripleft(value, '\t');
				it = table.find(var);
				if (it != table.end()) {
					pos = it->second.find(value);
					value = it->second.substr(pos, pos+value.size()+1);
					table.erase(it);
					table.insert(it, pair<string, string>(var, value));	// bug, leaves only var, not everything else
				}
			} else {
				pos=line.find("=");
				if (pos != string::npos) {
					var = line.substr(0,pos);
					var = strip(var, ' ');
					var = strip(var, '\t');
					value = line.substr(pos+1);
					value = stripleft(value, ' ');
					value = stripleft(value, '\t');
					it = table.find(var);
					if (it != table.end()) {
						table.erase(it);
					}
					table.insert(table.end(), pair<string, string>(var, value));
				}
			}
		}
	}
}

int main(int argc, char *argv[])
{
	int opt;
	int debug = 0, verbose = 0, trace = 0, plot = 0;
	char *prefix = NULL;
	char *value = NULL;
	string command = "get";
	string filename;
	string line;
	const string configurationfilebasename = "/harness/espadons/Makefile.configuration";
	
	struct option longopts[] = {
		{"add",			0, NULL, 'a'},
		{"get",			0, NULL, 'g'},
		{"rem",			0, NULL, 'r'},
		{"prefix",		1, NULL, 'P'},
		{"list",		0, NULL, 'l'},
		{"path",		0, NULL, 'A'},
		
		{"verbose",		0, NULL, 'v'},
		{"plot",		0, NULL, 'p'},
		{"debug",		0, NULL, 'd'},
		{"trace",		0, NULL, 't'},
		{"help",		0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "P:vpdthlagri", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'a':
				command = "add";
				break;
			case 'g':
				command = "get";
				break;
			case 'r':
				command = "rem";
				break;
			case 'P':
				prefix = optarg;
				break;
			case 'A':
				operaConfigurationAccessSetConfigurationFilepath(optarg);
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
				exit(EXIT_SUCCESS);
				break;
			default:
				break;
		}
	}
	table.clear();
	prefix = getenv("opera");
	//operaErrorCode operaError = operaConfigurationAccessGet("prefix-dir", &prefix);
	//if (operaError) {
	//	operaPError("operaConfigurationAccess", operaError);
	//}
	if (prefix != NULL) {
		filename = prefix + configurationfilebasename;
	} else {
		filename = ".." + configurationfilebasename;
	}
	// read in the configuration table
	fstream configstream;
	configstream.open(filename.c_str());
	if (configstream.is_open()) {
		while ( configstream.good() )
		{
			getline(configstream, line);
			if (verbose) {
				cout << line << endl;
			}
			if (line.length() > 0) {
				tableAdd(line);
			}
		}
		configstream.close();
	} else {
		cerr << "Unable to open configuration file " << filename << endl;
		exit(EXIT_FAILURE);
	}
	
    for ( ; optind < argc; optind++) {
		string arg = argv[optind];
		if (!arg.compare("get") || !arg.compare("add") || !arg.compare("rem"))
			command = argv[optind];
		else {
			line = argv[optind];
			size_t pos = line.find("=");
			if (pos != string::npos) {
				command = "add";
				tableAdd(line);
			}
		}
	}	
	// if get then print the value
	if (!command.compare("get") && line.length() > 0) {
		if (table.count(line) > 0) {
			operaErrorCode operaError = operaConfigurationAccessGet(table.find(line)->first.c_str(), &value);
			if (operaError) {
				operaPError("operaConfigurationAccess", operaError);
			} else {
				cout << value << endl;
			}
			//cout << table.find(line)->second << endl;
		} 
	}
	// if rem then print the value
	if (!command.compare("rem")) {
		map<string,string>::iterator it = table.find(line);
		if (it != table.end()) {
			table.erase(it);
		}
		else cerr << line << "not found." << endl;
	}
	// if add then write out the new table
	if (!command.compare("add") || !command.compare("rem")) {
		// make a backup copy...
		systemf("cp -f %s %s.bak", filename.c_str(), filename.c_str());
		map<string,string>::iterator it = table.begin();
		ofstream configstream;
		configstream.open(filename.c_str());
		if (configstream.is_open() ) {
			for ( ; it != table.end(); it++ ) {
				configstream << it->first << "\t:=\t" << it->second << endl;
			}
			configstream.close();
		}
			else cerr << "Unable to save configuration file " << filename << endl;
	}
	
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaConfigurationAccesss  [-vdth] [--prefix=...] [--instantiate] [--<var>[=|:=|+=|-=]<value>] [add|get|rem]\n"
	"  -h, --help\n"
	"  -v, --verbose\n"
	"  -d, --debug\n"
	" Example: operaConfigurationAccesss src=foo.c add			# adds new variable value pairvalue\n"
	" Example: operaConfigurationAccesss src:=foo.c add			# adds new variable value pair\n"
	" Example: operaConfigurationAccesss src=foo.c				# adds new value entry (the default action)\n"
	" Example: operaConfigurationAccesss src+=foo.c				# adds new value to the variable \"src\"\n"
	" Example: operaConfigurationAccesss src-=foo.c				# removes the value foo.c from the variable \"src\"\n"
	" Example: operaConfigurationAccesss src get				# returns value of variable \"src\"\n"
	" Example: operaConfigurationAccesss src rem				# removes variable\n"
	" Example: operaConfigurationAccesss src					# returns value of \"src\" (the default action)\n"
	" Example: operaConfigurationAccesss --instantiate				# instantiates the parameters for use by opera\n"
	" Example: operaConfigurationAccesss prefix=./opera-1.0 --instantiate # instantiates the parameters for use by opera\n"
	;
}

