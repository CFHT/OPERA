#ifndef OPERATELESCOPE_H
#define OPERATELESCOPE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaTelescope
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

enum opticalCoating_t {undefinedCoating=0, aluminium=1, silver=2};
enum telescopeMount_t {undefinedMount=0, equatorial=1, altazimuth=2};

/*! 
 * \sa class operaTelescope
 * \brief Encapsulation of telescope information.
 * \ingroup libraries
 */

class operaTelescope {
    
private:    
    double latitude;                    // in deg, (+) North, (-) South
    double longitude;                   // in deg, +East
    double elevation;                   // in meters
    
    double CollectingArea;              // in m^2
    double Aperture;                    // in meters
    double FocalRatio;                  // f-number:  f/focalratio

    telescopeMount_t TelescopeMount;    // supported types: equatorial, altazimuth
    opticalCoating_t OpticalCoating;    // supported types: aluminium (Al), silver (Ag), 
    
    
public:
		
	/*
	 * Constructors
	 */
    
	operaTelescope();
    
	/*
	 * Destructor
	 */   
    
	~operaTelescope();
	
	/*
	 * Methods for managing data
	 */    
    void setlatitude(double Latitude);
    double getlatitude(void);
    
    void setlongitude(double Longitude);
    double getlongitude(void);
    
    void setelevation(double Elevation);
    double getelevation(void);
    
    void setCollectingArea(double area);
    double getCollectingArea(void);
   
    void setAperture(double diameter);
    double getAperture(void);
    
    void setFocalRatio(double fnumber);
    double getFocalRatio(void);
    
    void setTelescopeMount(telescopeMount_t telescopeMountType);
    telescopeMount_t getTelescopeMount(void);
    
    void setOpticalCoating(opticalCoating_t coatingType);
    opticalCoating_t getOpticalCoating(void);
    
    /* 
     * Other Methods     
     */
    
};

#endif
