#ifndef OPERASPECTRUMSIMULATION_H
#define OPERASPECTRUMSIMULATION_H

/*******************************************************************
 ****                  LIBRARY FOR OPERA v1.0                   ****
 *******************************************************************
 Module name: operaSpectrumSimulation
 Version: 1.0
 Description: Generate various simulated spectra for polarimetry.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jun/2012
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

/*! 
 * \brief Spectrum Simulation.
 * \file operaSpectrumSimulation.h
 * \ingroup libraries
 */

using namespace std;

typedef enum { BlackBody=0, GalileanMoon, ZeemanLine} simulation_t;

typedef enum { NoError=0, UniformError, GaussianError} error_type_t;

typedef struct SpectrumVariables
{
    stokes_parameter_t StokesParameter;
    
    unsigned Beam;
    double Rhomb1Angle;
    double Rhomb2Angle;
    
    unsigned length;
    double MinWavelength;
    double MaxWavelength;
    double WavelengthIncrement;
    
    error_type_t ErrorType;
    double ErrorPercentage;
    
    /*
     * BlackBody
     */
    
    double Temperature;
    double StokesParameterFractionPolarizationOfBlackBody;
    
    /*
     * GalileanMoon
     */
    
    operaFluxVector *SunSpectrum;
    double ReflectedJupiterFractionOfSunSpectrum;
    
    unsigned NumberOfPolarizationPeaks;
    double* AmplitudeOfPolarizationVector;
    double* SigmaOfPolarizationVector;
    double* CenterOfPolarizationVector;
  
    /*
     * ZeemanLine
     */
    double ZeemanSplit;
    double linedepth;
    double linewidth; 
    double linecenter;
    
} SpectrumVariables_t;

/*
 * BlackBody
 */

operaFluxVector& SimulatedSpectrum(simulation_t SimulationType, SpectrumVariables_t SpectrumVariables, operaFluxVector& ESPaDOnSBeamIntensity);

operaFluxVector& BlackBodySimulation(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, operaFluxVector& SpectralRadianceVector);
operaFluxVector& BlackBodySimulationUniformError(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, double ErrorPercentage, operaFluxVector& SpectralRadianceVector);
operaFluxVector& BlackBodySimulationGaussianError(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, double ErrorPercentage, operaFluxVector& SpectralRadianceVector);

operaStokesVector& GenerateBlackBodyStokesVector(unsigned length, error_type_t ErrorType, double ErrorPercentage, double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, stokes_parameter_t StokesParameter, double StokesParameterFractionPolarizationOfBlackBody, operaStokesVector& BlackBodyStokesVector);

/*
 * GalileanMoons
 */

operaFluxVector& GaussianPolarizationSimulation(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, operaFluxVector& SimulatedPolarizationVector);
operaFluxVector& GaussianPolarizationSimulationUniformError(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, double ErrorPercentage, operaFluxVector& SimulatedPolarizationVector);
operaFluxVector& GaussianPolarizationSimulationGaussianError(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, double ErrorPercentage, operaFluxVector& SimulatedPolarizationVector);

operaStokesVector& GenerateGalileanMoonStokesVector(unsigned length, error_type_t ErrorType, double ErrorPercentage, double MinWavelength, double MaxWavelength, double WavelengthIncrement, operaFluxVector& SunSpectrum, double ReflectedJupiterFractionOfSunSpectrum, stokes_parameter_t StokesParameter, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, operaStokesVector& GalileanMoonStokesVector);

/*
 * ZeemanLine
 */
operaStokesVector& GenerateZeemanLineStokesVector(SpectrumVariables_t SpectrumVariables, operaStokesVector& ZeemanLineStokesVector);
operaFluxVector& IntroduceGaussianError(operaFluxVector& SimulatedPolarizationVector, double ErrorPercentage);
operaFluxVector& IntroduceUniformError(operaFluxVector& SimulatedPolarizationVector, double ErrorPercentage);
/*
 * BeamCreator
 */

operaStokesVector& ESPaDOnSBeamCreator(operaStokesVector &InputStokesVector, operaStokesVector& OutputStokesVector, double Rhomb1Angle, double Rhomb2Angle, unsigned Beam);

#endif
