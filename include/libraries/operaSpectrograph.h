#ifndef OPERASPECTROGRAPH_H
#define OPERASPECTROGRAPH_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectrograph
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

enum spectrographCCD_t {undefinedCCD=0, Olapa=1, EEV1=2}; // could add more ...

enum EspadonsCCDReadoutSpeed_t {undefinedReadoutMode=0, xslowmode=1, slowmode=2, normalmode=3, fastmode=4}; // specific to Espadons

enum EspadonsInstrumentMode_t {undefinedInstrumentMode=0, polarimetric=1, staronly=2, starplussky=3, GRACES_staronly=4, GRACES_starplussky=5};   // specific to Espadons

enum opticalFiber_t {undefinedFiber=0, FBPPolymicro=1, STUPolymicro=2}; // could add more ...

/*!
 * \sa class operaSpectrograph
 * \brief Encapsulation of spectrograph instrumental information.
 * \ingroup libraries
 */

class operaSpectrograph {
    
private:    

    double InjectionHoleDiameter;               // in arcsec
    
    opticalFiber_t  OpticalFiber;               // Supported types: FBPPolymicro=1, STUPolymicro=2
    double fiberLength;                         // in meters
    double fiberCoreDiameter;                   // in microns
    unsigned numberOfInputFibers;               // 1 for star-only, 2 for polar or star+sky
    unsigned numberOfSlices;                    // Espadons currently only uses 3 or 6 slices
    
    double spectralResolution;                  // R=lambda/dlambda

    spectrographCCD_t SpectrographCCD;                  // Supported types: Olapa, EEV1
    EspadonsCCDReadoutSpeed_t EspadonsCCDReadoutSpeed;  // Supported modes: slowmode, normalmode, fastmode
    EspadonsInstrumentMode_t EspadonsInstrumentMode;    // polarimetric, staronly, starplussky, GRACES_staronly, GRACES_starplussky
    
    double x_pixelsize;                         // in microns
    double y_pixelsize;                         // in microns
    
public:
		
	/*
	 * Constructors
	 */
    
	operaSpectrograph();
    
	/*
	 * Destructor
	 */   
    
	~operaSpectrograph();
	
	/*
	 * Methods for managing data
	 */    
    
    void setInjectionHoleDiameter(double injectionHoleDiameter);
    double getInjectionHoleDiameter(void);

    void setOpticalFiber(opticalFiber_t opticalFiberType);
    opticalFiber_t getOpticalFiber(void);
    
    void setfiberLength(double FiberLength);
    double getfiberLength(void);
    
    void setfiberCoreDiameter(double FiberCoreDiameter);
    double getfiberCoreDiameter(void);

    void setnumberOfInputFibers(unsigned NumberOfInputFibers);
    unsigned getnumberOfInputFibers(void);

    void setnumberOfSlices(unsigned NumberOfSlices);
    unsigned getnumberOfSlices(void);

    void setspectralResolution(double SpectralResolution);
    double getspectralResolution(void);
    
    void setSpectrographCCD(spectrographCCD_t spectrographCCD);
    spectrographCCD_t getSpectrographCCD(void);
    
    void setEspadonsCCDReadoutSpeed(EspadonsCCDReadoutSpeed_t espadonsCCDReadoutSpeed);
    EspadonsCCDReadoutSpeed_t getEspadonsCCDReadoutSpeed(void);
    
    void setEspadonsInstrumentMode(EspadonsInstrumentMode_t espadonsInstrumentMode);
    EspadonsInstrumentMode_t getEspadonsInstrumentMode(void);
    
    void setx_pixelsize(double x);
    double getx_pixelsize(void);

    void sety_pixelsize(double y);
    double gety_pixelsize(void);
    
    /* 
     * Other Methods     
     */
    
    double getSaturationExptime(double Vmag);
    
    double CalculateIPIE(double seeing);
    
    double OpticsThroughput(double wavelength_nm);

    double CCDQuantumEfficiency(double wavelength_nm);

    string getCCDName(void);
    
    double getNominalResolution (void);
    
    double getNominalAperture (void);
    
    unsigned getNumberOfExposures (void);

    string getFullModeName(void);
    
    double getNominalReadoutTime(void);
        
    double getNominalGain(void);
    
    double getNominalNoise(void);
    
    double FiberThroughput(double wavelength_nm);

    double STUDPolymicroFiber100m_throughput(double wavelength_nm);
    
    double FBPPolymicroFiber_throughput(double wavelength_nm);
    
    double FBPPolymicroFiber_attenuation_db_per_km(double wavelength_nm);
 
    double FiberTransmissionFromAttenuation(double attenuation_db_per_km);

    void setOpticalFiberFromInstrumentMode(void);
    
    double EEV1QuantumEfficiency(double wavelength_nm);
    
    double OlapaQuantumEfficiency(double wavelength_nm);

};

#endif
