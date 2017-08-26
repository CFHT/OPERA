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

#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaException.h"

#include "libraries/operaObservingConditions.h"		// for operaObservingConditions
#include "libraries/operaObjectInTheSky.h"          // for operaObjectInTheSky
#include "libraries/operaSpectrograph.h"            // for operaSpectrograph
#include "libraries/operaTelescope.h"               // for operaTelescope

#include "libraries/operaInstrumentEnvironmentSetup.h"

#include "libraries/operaLibCommon.h"
#include "libraries/gzstream.h"

/*!
 * operaInstrumentEnvironmentSetup
 * \author Doug Teeple and Eder Martioli
 * \brief instrument and environment setup.
 * \details {This library contains serialization and deserialization of instrument and environment setup.}
 * \file operaInstrumentEnvironmentSetup.cpp
 * \ingroup libraries
 */

using namespace std;

/*
 * Constructors
 */

/*!
 * \sa class operaInstrumentEnvironmentSetup();
 * \brief Base constructor.
 * \return void
 */
operaInstrumentEnvironmentSetup::operaInstrumentEnvironmentSetup(void) :
ObservingConditions(NULL),
ObjectInTheSky(NULL),
Spectrograph(NULL),
Telescope(NULL)
{
    ObservingConditions = new operaObservingConditions();
    ObjectInTheSky = new operaObjectInTheSky();
    Spectrograph = new operaSpectrograph();
    Telescope = new operaTelescope();
}

/*!
 * \sa class operaInstrumentEnvironmentSetup(string Filename);
 * \details Base constructor, read a setup object from a filename
 * \details Filename can contain either a instrument, target, or environment file as given by the format
 * \param Filename - string Filename to save to
 * \return void
 */
operaInstrumentEnvironmentSetup::operaInstrumentEnvironmentSetup(string Filename) :
ObservingConditions(NULL),
ObjectInTheSky(NULL),
Spectrograph(NULL),
Telescope(NULL)
{
    ObservingConditions = new operaObservingConditions();
    ObjectInTheSky = new operaObjectInTheSky();
    Spectrograph = new operaSpectrograph();
    Telescope = new operaTelescope();
    
    ReadInstrumentEnvironmentSetup(Filename);
}

/*
 * Destructor
 */
operaInstrumentEnvironmentSetup::~operaInstrumentEnvironmentSetup(void) {
    if(ObservingConditions)
        delete ObservingConditions;
    if(ObjectInTheSky)
        delete ObjectInTheSky;
    if(Spectrograph)
        delete Spectrograph;
    if(Telescope)
        delete Telescope;
}

/*
 * Methods
 */

/*!
 * \sa method operaObservingConditions *getObservingConditions();
 * \brief returns a pointer to the ObservingConditions class instance.
 * \return operaObservingConditions pointer.
 */
operaObservingConditions *operaInstrumentEnvironmentSetup::getObservingConditions(void) {
    return ObservingConditions;
}

/*!
 * \sa method operaObjectInTheSky *getObjectInTheSky();
 * \brief returns a pointer to the ObjectInTheSky class instance.
 * \return operaObjectInTheSky pointer.
 */
operaObjectInTheSky *operaInstrumentEnvironmentSetup::getObjectInTheSky(void) {
    return ObjectInTheSky;
}

/*!
 * \sa method operaSpectrograph *getSpectrograph();
 * \brief returns a pointer to the Spectrograph class instance.
 * \return operaSpectrograph pointer.
 */
operaSpectrograph *operaInstrumentEnvironmentSetup::getSpectrograph(void) {
    return Spectrograph;
}

/*!
 * \sa method operaTelescope *getTelescope();
 * \brief returns a pointer to the Telescope class instance.
 * \return operaTelescope pointer.
 */
operaTelescope *operaInstrumentEnvironmentSetup::getTelescope(void) {
    return Telescope;
}

/*!
 * \sa method void ReadInstrumentEnvironmentSetup(string Filename);
 * \brief augment an existing vector with information from a file
 * \param Filename - string.
 * \return none.
 */
void operaInstrumentEnvironmentSetup::ReadInstrumentEnvironmentSetup(string Filename) {
	InstrumentEnvironment_t format = InstrumentEnvironmentUnkownFormat;
    
	operaistream fin(Filename.c_str());
	if (fin.is_open()) {
		string dataline;
		if (fin.good()) {
			getline(fin, dataline);
			if (!dataline.compare("#!obscond")) {
				format = ObsCond;
			}
			if (!dataline.compare("#!skyobj")) {
				format = SkyObj;
			}
			if (!dataline.compare("#!spectrograph")) {
				format = spectrograph;
			}
			if (!dataline.compare("#!telescope")) {
				format = telescope;
			}
		}
	}
	fin.close();
    
	switch (format) {
		case ObsCond:
			readInstrumentEnvironmentFromObsCond(Filename);
			break;
		case SkyObj:
			readInstrumentEnvironmentFromSkyObj(Filename);
			break;
		case spectrograph:
			readInstrumentEnvironmentFromSpectrograph(Filename);
			break;
		case telescope:
			readInstrumentEnvironmentFromTelescope(Filename);
			break;
		default:
			throw operaException("operaInstrumentEnvironmentSetup: unkown content type in "+Filename+' ', operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
			break;
	}
}

/*
 * void operaInstrumentEnvironmentSetup::ReadInstrumentEnvironmentSetup(string Filename, operaSpectralOrder_t Format)
 * \brief Reads from an m.fits product to create Instrument/Environment Setup
 */
bool operaInstrumentEnvironmentSetup::ReadInstrumentEnvironmentSetup(string Filename, InstrumentEnvironment_t Format) {
    // Up to you Doug ...
	return false;
}

/*!
 * \sa method void WriteInstrumentEnvironmentSetup(string Filename, operaSpectralOrder_t Format);
 * \details Writes an Instrument/Environment setup to a File
 * \return operaInstrumentEnvironmentSetup* - pointer to the updated vector.
 */
void operaInstrumentEnvironmentSetup::WriteInstrumentEnvironmentSetup(string Filename, InstrumentEnvironment_t Format) {
	operaostream fout;
	if (!Filename.empty()) {
		fout.open(Filename.c_str());
		switch (Format) {
			case ObsCond: {
				fout << "#!obscond\n";
				fout << "######################################################################\n";
				fout << "# Observing Conditions format is:\n";
				fout << "# <JDTime (day)>\n";
				fout << "# <exposureTime (s)>\n";
				fout << "# <imageQuality (arcsec)>\n";
				fout << "# <airmass> \n";
				fout << "# <photometric (boolean 0 or 1)>\n";
				fout << "# <moonphase (supported values: {newmoon=0, crescentmoon=1, quartermoon=2, gibbousmoon=3, fullmoon=4})>\n";
				fout << "# <zenithalDistofMoon (deg)>\n";
				fout << "# <angularDistFromMoon (deg)>\n";
				fout << "#\n";
				fout << "# observercomments: " << getObservingConditions()->getobservercomments() << endl;
				fout << "# qccomments: " << getObservingConditions()->getqccomments() << endl;
				fout << "#\n";
				fout << "######################################################################\n";
				fout << getObservingConditions()->getJDTime() << endl;
				fout << getObservingConditions()->getexposureTime() << endl;
				fout << getObservingConditions()->getimageQuality() << endl;
				fout << getObservingConditions()->getairmass() << endl;
				fout << getObservingConditions()->getphotometric() << endl;
				fout << getObservingConditions()->getmoonphase() << endl;
				fout << getObservingConditions()->getzenithalDistofMoon() << endl;
				fout << getObservingConditions()->getangularDistFromMoon() << endl;
			}
				break;
			case SkyObj: {
				fout << "#!skyobj\n";
				fout << "######################################################################\n";
				fout << "# Object in the sky format is:\n";
				fout << "# <sourceID>\n";
				fout << "# <RA (deg)> <Dec (deg)>\n";
				fout << "# <ProperMotionRA (mas/yr)> <ProperMotionDec (mas/yr)> <Parallax (mas)>\n";
				fout << "# <V magnitude in Vega system (V-Johnson)>\n";
				fout << "# <EffectiveTemperature (K)>\n";
				fout << "# <RadialVelocity (km/s)>\n";
                fout << "# <SpectralType (Supported types O=1,B=2,A=3,F=4,G=5,K=6,M=7,L=8,T=9,Y=10) >\n";
				fout << "#\n";
				fout << "######################################################################\n";
				fout << getObjectInTheSky()->getsourceID() << endl;
				fout << getObjectInTheSky()->getRA() << ' ' << getObjectInTheSky()->getDec() << endl;
				fout << getObjectInTheSky()->getProperMotionRA() << ' ' << getObjectInTheSky()->getProperMotionDec() << ' ' << getObjectInTheSky()->getParallax() << endl;
				fout << getObjectInTheSky()->getV_magnitude() << endl;
				fout << getObjectInTheSky()->getEffectiveTemperature() << endl;
				fout << getObjectInTheSky()->getRadialVelocity() << endl;
				fout << getObjectInTheSky()->getSpectralType() << endl;
			}
				break;
			case spectrograph: {
				fout << "#!spectrograph\n";
				fout << "######################################################################\n";
				fout << "# Spectrograph format is:\n";
				fout << "# <InjectionHoleDiameter (arcsec)> \n";
				fout << "# <OpticalFiber (Supported types: FBPPolymicro=1, STUPolymicro=2)> \n";
				fout << "# <fiberLength (m)> \n";
				fout << "# <fiberCoreDiameter (microns)> \n";
				fout << "# <numberOfInputFibers> \n";
				fout << "# <numberOfSlices> \n";
				fout << "# <spectral Resolution (R=lambda/dlambda)> \n";
				fout << "# <Spectrograph CCD (Supported types: Olapa=1, EEV1=2)> \n";
				fout << "# <Espadons CCD Readout Speed (Supported modes: slowmode=1, normalmode=2, fastmode=3)> \n";
				fout << "# <Espadons Instrument Mode (Supported modes: polarimetric=1, staronly=2, starplussky=3, GRACES_staronly=4, GRACES_starplussky=5)> \n";
				fout << "# <pixel x-size (microns)> <pixel y-size (microns)>\n";
                fout << "#\n";
				fout << "######################################################################\n";
				fout << getSpectrograph()->getInjectionHoleDiameter() << endl;
				fout << getSpectrograph()->getOpticalFiber() << endl;
				fout << getSpectrograph()->getfiberLength() << endl;
				fout << getSpectrograph()->getfiberCoreDiameter() << endl;
				fout << getSpectrograph()->getnumberOfInputFibers() << endl;
				fout << getSpectrograph()->getnumberOfSlices() << endl;
				fout << getSpectrograph()->getspectralResolution() << endl;
				fout << getSpectrograph()->getSpectrographCCD() << endl;
				fout << getSpectrograph()->getEspadonsCCDReadoutSpeed() << endl;
				fout << getSpectrograph()->getEspadonsInstrumentMode() << endl;
				fout << getSpectrograph()->getx_pixelsize() << ' ' << getSpectrograph()->gety_pixelsize() << endl;
			}
				break;
            case telescope: {
				fout << "#!telescope\n";
				fout << "######################################################################\n";
				fout << "# Telescope format is:\n";
				fout << "# <latitude (deg, (+) North, (-) South)>\n";
				fout << "# <longitude (deg, +East)>\n";
				fout << "# <elevation (m)>\n";
				fout << "# <CollectingArea (m^2)>\n";
				fout << "# <Aperture (m)>\n";
				fout << "# <FocalRatio (f-number:  f/focalratio)>\n";
				fout << "# <Telescope Mount (supported types: equatorial=1, altazimuth=2)>\n";
				fout << "# <Optical Coating (supported types: aluminium(Al)=1, silver(Ag)=2)>\n";
				fout << "#\n";
				fout << "######################################################################\n";
				fout << getTelescope()->getlatitude() << endl;
				fout << getTelescope()->getlongitude() << endl;
				fout << getTelescope()->getelevation() << endl;
				fout << getTelescope()->getCollectingArea() << endl;
				fout << getTelescope()->getAperture() << endl;
				fout << getTelescope()->getFocalRatio() << endl;
				fout << getTelescope()->getTelescopeMount() << endl;
				fout << getTelescope()->getOpticalCoating() << endl;
			}
				break;
            default:
				break;
		}
		fout.close();
	}
}

/*
 * Read support routines...
 */
void operaInstrumentEnvironmentSetup::readInstrumentEnvironmentFromObsCond(string filename) {
	operaistream fobscond(filename.c_str());
	if (fobscond.is_open()) {
		string dataline;
        unsigned line = 0;
        
        double JDTime = 0.0;
        double exposureTime = 0.0;
        double imageQuality = 0.0;
        double airmass = 1.0;
        //bool photometric = true;
        //moonphase_t moonphase = quartermoon;
        double zenithalDistofMoon = 0.0;
        double angularDistFromMoon = 0.0;
		int scanin = 0;
        
		while (fobscond.good()) {
			getline(fobscond, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments, except observer/qc comments
					string otag = "observercomments:";
					string qtag = "qccomments:";
					if (dataline.find(otag) != string::npos) {
						string comments = dataline.substr(dataline.find(otag)+otag.length());
						getObservingConditions()->setobservercomments(comments);
					}
					if (dataline.find(qtag) != string::npos) {
						string comments = dataline.substr(dataline.find(qtag)+qtag.length());
						getObservingConditions()->setqccomments(comments);
					}
				} else if (line==0) {
					sscanf(dataline.c_str(), "%lf", &JDTime);
					getObservingConditions()->setJDTime(JDTime);
					line++;
				} else if (line==1) {
					sscanf(dataline.c_str(), "%lf", &exposureTime);
					getObservingConditions()->setexposureTime(exposureTime);
					line++;
				} else if (line==2) {
					sscanf(dataline.c_str(), "%lf", &imageQuality);
					getObservingConditions()->setimageQuality(imageQuality);
					line++;
				} else if (line==3) {
					sscanf(dataline.c_str(), "%lf", &airmass);
					getObservingConditions()->setairmass(airmass);
					line++;
				} else if (line==4) {
					int photom = 1;
					sscanf(dataline.c_str(), "%d", &photom);
					getObservingConditions()->setphotometric(photom==1);
					line++;
				} else if (line==5) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getObservingConditions()->setmoonphase((moonphase_t)scanin);
					line++;
				} else if (line==6) {
					sscanf(dataline.c_str(), "%lf", &zenithalDistofMoon);
					getObservingConditions()->setzenithalDistofMoon(zenithalDistofMoon);
					line++;
				} else if (line==7) {
					sscanf(dataline.c_str(), "%lf", &angularDistFromMoon);
					getObservingConditions()->setangularDistFromMoon(angularDistFromMoon);
					line++;
				}
			}
		}
		fobscond.close();
	}
}

void operaInstrumentEnvironmentSetup::readInstrumentEnvironmentFromSkyObj(string filename) {
	operaistream fskyobj(filename.c_str());
	if (fskyobj.is_open()) {
		string dataline;
        unsigned line = 0;
        
        string sourceID = NULL;
        double RA = 0.0, Dec = 0.0;
        double ProperMotionRA = 0.0, ProperMotionDec = 0.0, Parallax = 0.0;
        double V_magnitude = 0.0;
        double EffectiveTemperature = 0.0;
        double RadialVelocity = 0.0;
        //operaSpectralType_t SpectralType = A_type;
        
		while (fskyobj.good()) {
			getline(fskyobj, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line==0) {
					getObjectInTheSky()->setsourceID(string(dataline.c_str()));
					line++;
				} else if (line==1) {
					sscanf(dataline.c_str(), "%lf %lf", &RA, &Dec);
					getObjectInTheSky()->setRA(RA);
					getObjectInTheSky()->setDec(Dec);
					line++;
				} else if (line==2) {
					sscanf(dataline.c_str(), "%lf %lf %lf", &ProperMotionRA, &ProperMotionDec, &Parallax);
					getObjectInTheSky()->setProperMotionRA(ProperMotionRA);
					getObjectInTheSky()->setProperMotionDec(ProperMotionDec);
					getObjectInTheSky()->setParallax(Parallax);
					line++;
				} else if (line==3) {
					sscanf(dataline.c_str(), "%lf", &V_magnitude);
					getObjectInTheSky()->setV_magnitude(V_magnitude);
					line++;
				} else if (line==4) {
					sscanf(dataline.c_str(), "%lf", &EffectiveTemperature);
					getObjectInTheSky()->setEffectiveTemperature(EffectiveTemperature);
					line++;
				} else if (line==5) {
					sscanf(dataline.c_str(), "%lf", &RadialVelocity);
					getObjectInTheSky()->setRadialVelocity(RadialVelocity);
					line++;
				} else if (line==6) {
					int scanin;
					sscanf(dataline.c_str(), "%d", &scanin);
					getObjectInTheSky()->setSpectralType((operaSpectralType_t)scanin);
					line++;
				}
			}
		}
		fskyobj.close();
	}
}

void operaInstrumentEnvironmentSetup::readInstrumentEnvironmentFromSpectrograph(string filename) {
	operaistream fspectrograph(filename.c_str());
	if (fspectrograph.is_open()) {
		string dataline;
        unsigned line = 0;
        
        double InjectionHoleDiameter = 0.0;
        //opticalFiber_t OpticalFiber = undefinedFiber;
        double fiberLength=0.0;
        double fiberCoreDiameter=0.0;
        unsigned numberOfInputFibers=0;
        unsigned numberOfSlices=0;
        double spectralResolution=0.0;
        //spectrographCCD_t SpectrographCCD = undefinedCCD;
        //EspadonsCCDReadoutSpeed_t EspadonsCCDReadoutSpeed = undefinedReadoutMode;
        //EspadonsInstrumentMode_t EspadonsInstrumentMode = undefinedInstrumentMode;
        double x_pixelsize=0.0, y_pixelsize=0.0;
		int scanin = 0;
        
		while (fspectrograph.good()) {
			getline(fspectrograph, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line==0) {
					sscanf(dataline.c_str(), "%lf", &InjectionHoleDiameter);
					getSpectrograph()->setInjectionHoleDiameter(InjectionHoleDiameter);
					line++;
				} else if (line==1) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getSpectrograph()->setOpticalFiber((opticalFiber_t)scanin);
					line++;
				} else if (line==2) {
					sscanf(dataline.c_str(), "%lf", &fiberLength);
					getSpectrograph()->setfiberLength(fiberLength);
					line++;
				} else if (line==3) {
					sscanf(dataline.c_str(), "%lf", &fiberCoreDiameter);
					getSpectrograph()->setfiberCoreDiameter(fiberCoreDiameter);
					line++;
				} else if (line==4) {
					sscanf(dataline.c_str(), "%u", &numberOfInputFibers);
					getSpectrograph()->setnumberOfInputFibers(numberOfInputFibers);
					line++;
				} else if (line==5) {
					sscanf(dataline.c_str(), "%u", &numberOfSlices);
					getSpectrograph()->setnumberOfSlices(numberOfSlices);
					line++;
				} else if (line==6) {
					sscanf(dataline.c_str(), "%lf", &spectralResolution);
					getSpectrograph()->setspectralResolution(spectralResolution);
					line++;
				} else if (line==7) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getSpectrograph()->setSpectrographCCD((spectrographCCD_t)scanin);
					line++;
				} else if (line==8) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getSpectrograph()->setEspadonsCCDReadoutSpeed((EspadonsCCDReadoutSpeed_t)scanin);
					line++;
				} else if (line==9) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getSpectrograph()->setEspadonsInstrumentMode((EspadonsInstrumentMode_t)scanin);
					line++;
				} else if (line==10) {
					sscanf(dataline.c_str(), "%lf %lf", &x_pixelsize, &y_pixelsize);
					getSpectrograph()->setx_pixelsize(x_pixelsize);
					getSpectrograph()->sety_pixelsize(y_pixelsize);
					line++;
				}
			}
		}
		fspectrograph.close();
	}
}

void operaInstrumentEnvironmentSetup::readInstrumentEnvironmentFromTelescope(string filename) {
	operaistream ftelescope(filename.c_str());
	if (ftelescope.is_open()) {
		string dataline;
        unsigned line = 0;
        
        double latitude=0.0;
        double longitude=0.0;
        double elevation=0.0;
        double CollectingArea=0.0;
        double Aperture=0.0;
        double FocalRatio=0.0;
        //telescopeMount_t TelescopeMount = undefinedMount;
        //opticalCoating_t OpticalCoating = undefinedCoating;
        int scanin = 0;
        
		while (ftelescope.good()) {
			getline(ftelescope, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line==0) {
					sscanf(dataline.c_str(), "%lf", &latitude);
					getTelescope()->setlatitude(latitude);
					line++;
				} else if (line==1) {
					sscanf(dataline.c_str(), "%lf", &longitude);
					getTelescope()->setlongitude(longitude);
					line++;
				} else if (line==2) {
					sscanf(dataline.c_str(), "%lf", &elevation);
					getTelescope()->setelevation(elevation);
					line++;
				} else if (line==3) {
					sscanf(dataline.c_str(), "%lf", &CollectingArea);
					getTelescope()->setCollectingArea(CollectingArea);
					line++;
				} else if (line==4) {
					sscanf(dataline.c_str(), "%lf", &Aperture);
					getTelescope()->setAperture(Aperture);
					line++;
				} else if (line==5) {
					sscanf(dataline.c_str(), "%lf", &FocalRatio);
					getTelescope()->setFocalRatio(FocalRatio);
					line++;
				} else if (line==6) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getTelescope()->setTelescopeMount((telescopeMount_t)scanin);
					line++;
				} else if (line==7) {
					sscanf(dataline.c_str(), "%d", &scanin);
					getTelescope()->setOpticalCoating((opticalCoating_t)scanin);
					line++;
				}
			}
            
		}
		ftelescope.close();
	}
}



