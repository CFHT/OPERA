#ifndef OPERASPECTRALENERGYDISTRIBUTION_H
#define OPERASPECTRALENERGYDISTRIBUTION_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralEnergyDistribution
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

#include "libraries/operaFluxVector.h"

/*! 
 * \sa class operaSpectralEnergyDistribution
 * \brief Encapsulation of Spectral Energy Distribution information.
 * \ingroup libraries
 */

class operaSpectralOrder;

class operaSpectralEnergyDistribution {
	
private:    
    
    unsigned nDataPoints;
    
    operaVector distanceData;
    operaVector wavelengthData;
    operaVector fluxData; // flux measurements in instrumental units (ADU/pixel^2/spectralElement) - no point in using a flux vector since we never keep track of the variance
    double wavelengthForNormalization;
    
    operaVector calibrationDist;
    operaVector calibrationWavelength;
    operaFluxVector uncalibratedFlux;
    operaFluxVector calibratedFlux;
    operaFluxVector fluxCalibration;
    operaFluxVector instrumentThroughput;
    
    bool hasFluxData;
	
	bool hasUncalibratedFlux;
	bool hasCalibratedFlux;
	bool hasFluxCalibration;
   	bool hasInstrumentThroughput;      

public:
		
	operaSpectralEnergyDistribution();
    
	operaSpectralEnergyDistribution(unsigned NDataPoints);   

    operaSpectralEnergyDistribution(unsigned NDataPoints, unsigned nElements);  
    
    /*
	 * Methods for managing data
	 */    
    
    unsigned getnDataPoints() const;
    
    void resizeDataVectors(unsigned NDataPoints);
    
    void resizeCalibrationVectors(unsigned nElements);
    
    void setwavelengthForNormalization(double WavelengthForNormalization);
    
    double getwavelengthForNormalization() const;
    
	void setdistanceData(double Distance, unsigned index);
    
    void setwavelengthData(double Wavelength, unsigned index);
    
    void setfluxData(double Flux, unsigned index);
    
    double getdistanceData(unsigned index) const;
    
    double getwavelengthData(unsigned index) const;
    
    double getfluxData(unsigned index) const;
    
    bool getHasFluxData() const { return hasFluxData; }
    
	void setHasFluxData(bool HasFluxData) { hasFluxData = HasFluxData; }
    
	/*
	 * Methods for managing flux calibration elements
	 */
	  
	void setCalibrationDist(const operaVector& CalibrationDist);
	
	void setCalibrationWavelength(const operaVector& CalibrationWavelength);
	
	void setUncalibratedFlux(const operaFluxVector& UncalibratedFluxElements);

    void setCalibratedFlux(const operaFluxVector& CalibratedFluxElements);
    
    void setFluxCalibration(const operaFluxVector& FluxCalibrationElements);
    
    void setThroughput(const operaFluxVector& ThroughputElements);
    
    operaVector& getCalibrationWavelength(){ return calibrationWavelength; }
    
    const operaVector& getCalibrationWavelength() const { return calibrationWavelength; }
    
    /*operaFluxVector& getUncalibratedFlux(){ return uncalibratedFlux; }
    
    const operaFluxVector& getUncalibratedFlux() const { return uncalibratedFlux; }*/
    
    operaFluxVector& getUncalibratedFlux(){ return uncalibratedFlux; }
    
    const operaFluxVector& getUncalibratedFlux() const { return uncalibratedFlux; }
    
    operaFluxVector& getCalibratedFlux(){ return calibratedFlux; }
    
    const operaFluxVector& getCalibratedFlux() const { return calibratedFlux; }
    
    operaFluxVector& getFluxCalibration(){ return fluxCalibration; }
    
    const operaFluxVector& getFluxCalibration() const { return fluxCalibration; }
    
    operaFluxVector& getThroughput(){ return instrumentThroughput; }
    
    const operaFluxVector& getThroughput() const { return instrumentThroughput; }
        
    operaFluxVector& getCalibration(bool useThroughput) { if(useThroughput) return instrumentThroughput; return fluxCalibration; }
    
    const operaFluxVector& getCalibration(bool useThroughput) const { if(useThroughput) return instrumentThroughput; return fluxCalibration; }
    
    bool getHasUncalibratedFlux() const { return hasUncalibratedFlux; };

	void setHasUncalibratedFlux(bool HasUncalibratedFlux) { hasUncalibratedFlux = HasUncalibratedFlux; }

	bool getHasCalibratedFlux() const { return hasCalibratedFlux; }

	void setHasCalibratedFlux(bool HasCalibratedFlux) { hasCalibratedFlux = HasCalibratedFlux; }

	bool getHasFluxCalibration() const { return hasFluxCalibration; }

	void setHasFluxCalibration(bool HasFluxCalibration) { hasFluxCalibration = HasFluxCalibration; }
    
	bool getHasInstrumentThroughput() const { return hasInstrumentThroughput; }

	void setHasInstrumentThroughput(bool HasInstrumentThroughput) { hasInstrumentThroughput = HasInstrumentThroughput; }
    
    /*
     * Other Methods     
     */
    
    void measureUncalibratedContinuum(unsigned binsize, unsigned nsigcut);
    
    void populateUncalibratedFluxFromContinuumData();
    
    void calculateUncalibratedFlux(unsigned binsize, unsigned nsigcut);
};
#endif
