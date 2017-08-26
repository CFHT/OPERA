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
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "operaError.h"
#include "globaldefines.h"

#include "libraries/operaFit.h"               // for cubicsplineDouble and splineinterpolateDouble

#include "libraries/operaObjectInTheSky.h"  // for operaObjectInTheSky

#include "libraries/operaException.h"
#include "libraries/operaObservingConditions.h"  // for operaObservingConditions

#include "libraries/operaLibCommon.h"
#include "libraries/gzstream.h"

/*!
 * operaObservingConditions
 * \author Eder Martioli
 * \brief This class encapsulates the spectral energy distribution object.
 * \file operaObservingConditions.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaObservingConditions
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

/*
 * Constructors
 */

operaObservingConditions::operaObservingConditions() :
JDTime(0.0),
exposureTime(0.0),
imageQuality(0.0),
airmass(1.0),
photometric(true),
moonphase(quartermoon),
zenithalDistofMoon(0.0),
angularDistFromMoon(0.0)
{
}

/*
 * Destructor
 */

operaObservingConditions::~operaObservingConditions() {
  
}

/*
 * Methods for managing data
 */
void operaObservingConditions::setJDTime(double jd) {
    JDTime = jd;
}

double operaObservingConditions::getJDTime(void) {
    return JDTime;
}

void operaObservingConditions::setexposureTime(double exptime) {
    exposureTime = exptime;
}

double operaObservingConditions::getexposureTime(void) {
    return exposureTime;
}

void operaObservingConditions::setimageQuality(double IQ) {
    imageQuality = IQ;
}

double operaObservingConditions::getimageQuality(void) {
    return imageQuality;
}

void operaObservingConditions::setairmass(double Airmass) {
    airmass = Airmass;
}

double operaObservingConditions::getairmass(void) {
    return airmass;
}

void operaObservingConditions::setphotometric(bool Photometric) {
    photometric = Photometric;
}

bool operaObservingConditions::getphotometric(void) {
    return photometric;
}

void operaObservingConditions::setmoonphase(moonphase_t MoonPhase) {
    moonphase = MoonPhase;
}

moonphase_t operaObservingConditions::getmoonphase(void) {
    return moonphase;
}

void operaObservingConditions::setzenithalDistofMoon(double ZenithalDistofMoon) {
    zenithalDistofMoon = ZenithalDistofMoon;
}

double operaObservingConditions::getzenithalDistofMoon(void) {
    return zenithalDistofMoon;
}

void operaObservingConditions::setangularDistFromMoon(double AngularDistFromMoon) {
    angularDistFromMoon = AngularDistFromMoon;
}

double operaObservingConditions::getangularDistFromMoon(void) {
    return angularDistFromMoon;
}

void operaObservingConditions::setobservercomments(string Comments) {
	observercomments = Comments;
}

string operaObservingConditions::getobservercomments(void) {
	return observercomments;
}

void operaObservingConditions::setqccomments(string Comments) {
	qccomments = Comments;
}

string operaObservingConditions::getqccomments(void) {
	return qccomments;
}


double operaObservingConditions::getMoonDay(void) {
    switch (moonphase) {
        case newmoon:
            return 0;
            break;
        case crescentmoon:
            return 3.5;
            break;
        case quartermoon:
            return 7.0;
            break;
        case gibbousmoon:
            return 10.5;
            break;
        case fullmoon:
            return 14.0;
            break;
       default:
            break;
    }
    return 7.0;
}

double operaObservingConditions::getMoonAlpha(void) {
    switch (moonphase) {
        case newmoon:
            return 180;
            break;
        case crescentmoon:
            return 135;
            break;
        case quartermoon:
            return 90;
            break;
        case gibbousmoon:
            return 45;
            break;
        case fullmoon:
            return 0;
            break;
        default:
            break;
    }
    return 90;
}

double operaObservingConditions::getZenithalDistanceOfSource(void) {
    return ((180./M_PI)*acos(1./airmass));
}

double operaObservingConditions::getMoonDeltaV(void) {
    
    // reference:  Krisciunas K. & Schaefer B. E., PASP, v. 103, p. 1033-1039, Sep 1991.
    double BMoon, BSky, B0;
    double Iflux, frho;
    double XZ, XZMoon, kextc;
    kextc = 0.172;
    B0 = 79.0;
    
    Iflux = pow(10,(-0.4*(3.84 + 0.026*getMoonAlpha() + (4e-9)*pow(getMoonAlpha(),4))));
    
    frho = pow(10,5.36)*(1.06 + cos(angularDistFromMoon*M_PI/180)*cos(angularDistFromMoon*M_PI/180)) + pow(10, 6.15 - (angularDistFromMoon/40));
    
    XZ = 1/sqrt(1 - 0.96*sin(getZenithalDistanceOfSource()*M_PI/180)*sin(getZenithalDistanceOfSource()*M_PI/180));
    XZMoon = 1/sqrt(1 - 0.96*sin(zenithalDistofMoon*M_PI/180)*sin(zenithalDistofMoon*M_PI/180));
    
    BSky = B0*pow(10,-0.4*kextc*(XZ - 1))*XZ;
    
    BMoon = frho*Iflux*pow(10,-0.4*kextc*XZMoon)*(1 - pow(10,-0.4*kextc*XZ));
    
    return (-2.5*log10((BMoon + BSky)/BSky));
}


double operaObservingConditions::getSkyZeroVMagnitude(double wavelength_nm) {
    
    double wl[40],skyzeromag[40];
    
    wl[0]=372;
    wl[1]=378;
    wl[2]=384;
    wl[3]=391;
    wl[4]=398;
    wl[5]=405;
    wl[6]=412;
    wl[7]=420;
    wl[8]=428;
    wl[9]=436;
    wl[10]=444;
    wl[11]=453;
    wl[12]=463;
    wl[13]=472;
    wl[14]=482;
    wl[15]=493;
    wl[16]=504;
    wl[17]=515;
    wl[18]=527;
    wl[19]=540;
    wl[20]=553;
    wl[21]=567;
    wl[22]=581;
    wl[23]=597;
    wl[24]=613;
    wl[25]=630;
    wl[26]=648;
    wl[27]=667;
    wl[28]=687;
    wl[29]=709;
    wl[30]=732;
    wl[31]=756;
    wl[32]=782;
    wl[33]=810;
    wl[34]=840;
    wl[35]=872;
    wl[36]=907;
    wl[37]=945;
    wl[38]=986;
    wl[39]=1031;
    
    skyzeromag[0]=23.20;
    skyzeromag[1]=23.30;
    skyzeromag[2]=23.20;
    skyzeromag[3]=23.30;
    skyzeromag[4]=23.25;
    skyzeromag[5]=23.00;
    skyzeromag[6]=23.10;
    skyzeromag[7]=23.15;
    skyzeromag[8]=23.15;
    skyzeromag[9]=23.10;
    skyzeromag[10]=23.03;
    skyzeromag[11]=22.91;
    skyzeromag[12]=23.05;
    skyzeromag[13]=23.03;
    skyzeromag[14]=22.79;
    skyzeromag[15]=22.90;
    skyzeromag[16]=22.90;
    skyzeromag[17]=22.84;
    skyzeromag[18]=22.84;
    skyzeromag[19]=22.84;
    skyzeromag[20]=22.60;
    skyzeromag[21]=22.70;
    skyzeromag[22]=22.50;
    skyzeromag[23]=22.50;
    skyzeromag[24]=22.40;
    skyzeromag[25]=22.25;
    skyzeromag[26]=22.20;
    skyzeromag[27]=22.30;
    skyzeromag[28]=22.30;
    skyzeromag[29]=22.30;
    skyzeromag[30]=22.20;
    skyzeromag[31]=22.20;
    skyzeromag[32]=22.25;
    skyzeromag[33]=22.00;
    skyzeromag[34]=21.90;
    skyzeromag[35]=21.60;
    skyzeromag[36]=21.10;
    skyzeromag[37]=20.70;
    skyzeromag[38]=20.20;
    skyzeromag[39]=20.40;
    
    unsigned nin = 40;
    
    double yp1 = (skyzeromag[1] - skyzeromag[0])/(wl[1] - wl[0]);
	double ypn = (skyzeromag[nin-1] - skyzeromag[nin-2])/(wl[nin-1] - wl[nin-2]);
	double *y2 = (double *)malloc(nin*sizeof(double));
	
	// Call cubicspline to get second derivatives
	cubicsplineDouble(wl, skyzeromag, nin, yp1, ypn, y2);
	
    double outputskyzeromag;
    
	// Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, skyzeromag, y2, nin, wavelength_nm , &outputskyzeromag);

    free(y2);
    
    return outputskyzeromag;
}

double operaObservingConditions::getAtmosphericTransmissionAtMaunaKea(string filename, double wavelength_nm) {
    
    double *wavelength = new double[MAXNUMBEROFPOINTSINATMTRANSMISSION];
    double *transmission = new double[MAXNUMBEROFPOINTSINATMTRANSMISSION];
	
    igzstream astream;
	string dataline;
    
	double tmpwl = -1.0;
	double tmpi = -1.0;
	unsigned np = 0;
   
    /*
     * Read the the full atmospheric transmission spectrum
     */
	astream.open(filename.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    
                    wavelength[np] = (tmpwl/10.0);
                    transmission[np] = tmpi;
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open()
    
    unsigned nin = np;
    
    double yp1 = (transmission[1] - transmission[0])/(wavelength[1] - wavelength[0]);
	double ypn = (transmission[np-1] - transmission[np-2])/(wavelength[np-1] -wavelength[np-2]);
	double *y2 = new double[np];
	
	// Call cubicspline to get second derivatives
	cubicsplineDouble(wavelength, transmission, np, yp1, ypn, y2);
	
    double TransmissionAtZenith;
    
	// Call splineinterpolate for interpolations
    splineinterpolateDouble(wavelength, transmission, y2, nin, wavelength_nm , &TransmissionAtZenith);
    
    delete[] y2;
    delete[] wavelength;
    delete[] transmission;
    
    double outputTransmission = exp(log(TransmissionAtZenith)*airmass);

    return outputTransmission;
}

double operaObservingConditions::calculateSkyFlux(double wavelength_nm, double dwl_m, double InjectionHoleDiameter) {
    double wavelength_m = wavelength_nm*1e-9;

    double wl0 = wavelength_m - dwl_m/2;
    double wlf = wavelength_m + dwl_m/2;
    /*
     * Figure out sky brightness due to the Moon
     */
    operaObjectInTheSky Moon;
    Moon.setEffectiveTemperature(SUN_EFF_TEMPERATURE_K);
    Moon.setSpectralType(G_type);
    double SkyVMag = 21.1 + getMoonDeltaV();
    double MoonSkyVMag = -2.5*log10(pow(10,-0.4*(SkyVMag - VEGA_VBAND_MAGNITUDE)) - pow(10,-0.4*(21.1 - VEGA_VBAND_MAGNITUDE)));
    Moon.setV_magnitude(MoonSkyVMag);
    
    double SkyMoonEmit = Moon.getSpectralBinFlux(wl0,wlf);
    double SkyMag = getSkyZeroVMagnitude(wavelength_nm);
    double SkyFlux = pow(10,-(SkyMag - VEGA_VBAND_MAGNITUDE)/2.5)*(VEGA_PHOTONFLUX_AT_548NM)*(M_PI*(InjectionHoleDiameter/2)*(InjectionHoleDiameter/2))*dwl_m + SkyMoonEmit;
    
    return SkyFlux;
}
