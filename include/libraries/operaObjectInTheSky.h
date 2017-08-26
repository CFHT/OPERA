#ifndef OPERAOBJECTINTHESKY_H
#define OPERAOBJECTINTHESKY_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaObjectInTheSky
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

enum operaSpectralType_t {
    O_type=1,
    B_type=2,
    A_type=3,
    F_type=4,
    G_type=5,
    K_type=6,
    M_type=7,
    L_type=8,
    T_type=9,
    Y_type=10
};

/*!
 * \sa class operaObjectInTheSky
 * \brief This class encapsulates the information of an object in the sky.
 * \ingroup libraries
 */

class operaObjectInTheSky {
	
private:    
    string sourceID;                        // source identifier: name, catalog number, etc
    double RA;                              // in degrees
    double Dec;                             // in degrees
    double ProperMotionRA;                  // in mas/yr
    double ProperMotionDec;                 // in mas/yr
    double Parallax;                        // in mas
    double V_magnitude;                     // in Vega magnitude in V-band (Johnson)
    double EffectiveTemperature;            // in Kelvin
    double RadialVelocity;                  // in km/s
    operaSpectralType_t SpectralType;       // main sequence (V): O,B,A,F,G,K,M,L,T,Y

public:
		
	/*
	 * Constructors
	 */
    
	operaObjectInTheSky();

	/*
	 * Destructor
	 */   
    
	~operaObjectInTheSky();
	
	/*
	 * Methods for managing data
	 */    
    
    void setsourceID(string SourceID);
    string getsourceID(void);

    void setRA(double alpha);
    double getRA(void);
    
    void setDec(double delta);
    double getDec(void);

    void setProperMotionRA(double muRA);
    double getProperMotionRA(void);

    void setProperMotionDec(double muDec);
    double getProperMotionDec(void);

    void setParallax(double pi);
    double getParallax(void);
    
    void setV_magnitude(double Vmag);
    double getV_magnitude(void);
    
    void setEffectiveTemperature(double Teff);
    double getEffectiveTemperature(void);
    
    void setRadialVelocity(double RV);
    double getRadialVelocity(void);

    void setSpectralType(operaSpectralType_t SpecTy);
    operaSpectralType_t getSpectralType(void);
    
    /*
     * Other Methods
     */
    
    double getBlackBodyVFlux();
    double calculateVFlux(double Temperature);
    double calculatePhiFactor(void);
    double getSpectralBinFlux(double wl0, double wlf);
    
};

#endif
