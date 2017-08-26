#ifndef OPERAINSTRUMENTENVIRONMENTSETUP_H
#define OPERAINSTRUMENTENVIRONMENTSETUP_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaInstrumentEnvironmentSetup
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Feb/2013
 
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

#include "operaError.h"
#include "operaLibCommon.h"
#include "libraries/operaObservingConditions.h"		// for operaObservingConditions
#include "libraries/operaObjectInTheSky.h"          // for operaObjectInTheSky
#include "libraries/operaSpectrograph.h"            // for operaSpectrograph
#include "libraries/operaTelescope.h"               // for operaTelescope


using namespace std;

/*! 
 * \sa class operaInstrumentEnvironmentSetup
 * \brief Stores all instrument and environment setup. 
 * \return none
 * \file operaInstrumentEnvironmentSetup.h
 * \ingroup libraries
 */
class operaInstrumentEnvironmentSetup {
	
private:
    operaObservingConditions *ObservingConditions;
    operaObjectInTheSky *ObjectInTheSky;
    operaSpectrograph *Spectrograph;
    operaTelescope *Telescope;
    
public:
    /*
	 * Constructors
	 */
	
	/*! 
	 * \sa class operaInstrumentEnvironmentSetup();
	 * \brief Base constructor.
	 * \return void
	 */
	operaInstrumentEnvironmentSetup(void);
    
	/*! 
	 * \sa class operaInstrumentEnvironmentSetup(string Filename);
	 * \details Base constructor, read a setup object from a filename
	 * \details Filename can contain either a instrument, target, or environment file as given by the format
	 * \param Filename - string Filename to save to
	 * \return void
	 */
	operaInstrumentEnvironmentSetup(string Filename);
    
	/*
	 * Destructor
	 */
	~operaInstrumentEnvironmentSetup(void);

	/*
	 * Methods
	 */
	
	/*!
	 * \sa method operaObservingConditions *getObservingConditions();
	 * \brief returns a pointer to the ObservingConditions class instance.
	 * \return operaObservingConditions pointer.
	 */
	operaObservingConditions *getObservingConditions(void);
    
	/*!
	 * \sa method operaObjectInTheSky *getObjectInTheSky();
	 * \brief returns a pointer to the ObjectInTheSky class instance.
	 * \return operaObjectInTheSky pointer.
	 */
	operaObjectInTheSky *getObjectInTheSky(void);
    
	/*!
	 * \sa method operaSpectrograph *getSpectrograph();
	 * \brief returns a pointer to the Spectrograph class instance.
	 * \return operaSpectrograph pointer.
	 */
	operaSpectrograph *getSpectrograph(void);
    
	/*!
	 * \sa method operaTelescope *getTelescope();
	 * \brief returns a pointer to the Telescope class instance.
	 * \return operaTelescope pointer.
	 */
	operaTelescope *getTelescope(void);

	/*!
	 * \sa method void ReadInstrumentEnvironmentSetup(string Filename);
	 * \brief augment an existing vector with information from a file
	 * \param Filename - string.
	 * \return none.
	 */
	void ReadInstrumentEnvironmentSetup(string Filename);
    
	/*
	 * void operaInstrumentEnvironmentSetup::ReadInstrumentEnvironmentSetup(string Filename, operaSpectralOrder_t Format)
	 * \brief Reads from an m.fits product to create Instrument/Environment Setup
	 */
	bool ReadInstrumentEnvironmentSetup(string Filename, InstrumentEnvironment_t Format);
    
	/*!
	 * \sa method void WriteInstrumentEnvironmentSetup(string Filename, operaSpectralOrder_t Format);
	 * \details Writes an Instrument/Environment setup to a File
	 * \return operaInstrumentEnvironmentSetup* - pointer to the updated vector.
	 */
	void WriteInstrumentEnvironmentSetup(string Filename, InstrumentEnvironment_t Format);
	/*
	 * Read support routines...
	 */
	void readInstrumentEnvironmentFromObsCond(string filename);
    
	void readInstrumentEnvironmentFromSkyObj(string filename);
	
	void readInstrumentEnvironmentFromSpectrograph(string filename);
	
    void readInstrumentEnvironmentFromTelescope(string filename);
        
};

#endif
