#ifndef OPERAOBSERVINGCONDITIONS_H
#define OPERAOBSERVINGCONDITIONS_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaObservingConditions
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

#ifndef MAXNUMBEROFPOINTSINATMTRANSMISSION
#define MAXNUMBEROFPOINTSINATMTRANSMISSION 10000
#endif

enum moonphase_t {newmoon=0, crescentmoon=1, quartermoon=2, gibbousmoon=3, fullmoon=4};

/*! 
 * \sa class operaObservingConditions
 * \brief Encapsulation of observing conditions information.
 * \ingroup libraries
 */

class operaObservingConditions {
    
private:    

    double JDTime;                  // Time in Julian days at start of exposure
    double exposureTime;            // Exposure time in seconds
    double imageQuality;            // seeing in arcsec
    double airmass;                 // airmass@zenith = 1
    bool photometric;               // based on extinction (true or false)
    moonphase_t moonphase;          // 5 phases: new, crescent, quarter, gibbous, full
    double zenithalDistofMoon;      // in degrees
    double angularDistFromMoon;     // in degrees
    string observercomments;		// ...
    string qccomments;				// ...
        
public:
		
	/*
	 * Constructors
	 */
    
	operaObservingConditions();
    
	/*
	 * Destructor
	 */   
    
	~operaObservingConditions();
	
	/*
	 * Methods for managing data
	 */    
    void setJDTime(double jd);
    double getJDTime(void);
    
    void setexposureTime(double exptime);
    double getexposureTime(void);
    
    void setimageQuality(double IQ);
    double getimageQuality(void);

    void setairmass(double Airmass);
    double getairmass(void);
    
    void setphotometric(bool Photometric);
    bool getphotometric(void);
    
    void setmoonphase(moonphase_t MoonPhase);
    moonphase_t getmoonphase(void);

    void setzenithalDistofMoon(double ZenithalDistofMoon);
    double getzenithalDistofMoon(void);
    
    void setangularDistFromMoon(double AngularDistFromMoon);
    double getangularDistFromMoon(void);
    
    void setobservercomments(string Comments);
    string getobservercomments(void);
    
    void setqccomments(string Comments);
    string getqccomments(void);
    
    /*
     * Other Methods     
     */
    double getMoonDay(void);
    
    double getMoonAlpha(void);

    double getZenithalDistanceOfSource(void);
    
    double getMoonDeltaV(void);
    
    double getSkyZeroVMagnitude(double wavelength_nm);
    
    double getAtmosphericTransmissionAtMaunaKea(string filename, double wavelength_nm);
    
    double calculateSkyFlux(double wavelength_nm, double dwl, double InjectionHoleDiameter);

};

#endif
