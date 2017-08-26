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

#include <time.h>       // for random errors
#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaMath.h"    // for PlanckLaw
#include "libraries/operaFluxVector.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaMuellerMatrix.h"
#include "libraries/operaStats.h"   // for Uniform and Gaussian error
#include "libraries/Gaussian.h"     // for Gaussian class
#include "libraries/operaSpectrumSimulation.h" 

/*! 
 * SimulatedSpectrum
 * \author Eder Martioli
 * \brief This class simulates a spectrum
 * \file SimulatedSpectrum.cpp
 * \ingroup libraries
 */

using namespace std;

operaFluxVector& SimulatedSpectrum(simulation_t SimulationType, SpectrumVariables_t SpectrumVariables, operaFluxVector& ESPaDOnSBeamIntensity)
{
    switch (SimulationType) {
        case BlackBody: {
            operaStokesVector BlackBodyStokesVector = GenerateBlackBodyStokesVector(SpectrumVariables.length, SpectrumVariables.ErrorType, SpectrumVariables.ErrorPercentage, SpectrumVariables.Temperature, SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, SpectrumVariables.StokesParameter, SpectrumVariables.StokesParameterFractionPolarizationOfBlackBody, BlackBodyStokesVector);
            operaStokesVector ESPaDOnSBeam = ESPaDOnSBeamCreator(BlackBodyStokesVector, ESPaDOnSBeam, SpectrumVariables.Rhomb1Angle, SpectrumVariables.Rhomb2Angle, SpectrumVariables.Beam);
            
            ESPaDOnSBeamIntensity = *ESPaDOnSBeam.getStokesParameter(StokesI);
        }
            break;
        case GalileanMoon: {
            operaStokesVector GalileanMoonStokesVector = GenerateGalileanMoonStokesVector(SpectrumVariables.length, SpectrumVariables.ErrorType, SpectrumVariables.ErrorPercentage, SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, *SpectrumVariables.SunSpectrum, SpectrumVariables.ReflectedJupiterFractionOfSunSpectrum, SpectrumVariables.StokesParameter, SpectrumVariables.NumberOfPolarizationPeaks, SpectrumVariables.AmplitudeOfPolarizationVector, SpectrumVariables.SigmaOfPolarizationVector, SpectrumVariables.CenterOfPolarizationVector, GalileanMoonStokesVector);
            operaStokesVector ESPaDOnSBeam = ESPaDOnSBeamCreator(GalileanMoonStokesVector, ESPaDOnSBeam, SpectrumVariables.Rhomb1Angle, SpectrumVariables.Rhomb2Angle, SpectrumVariables.Beam);
            
			ESPaDOnSBeamIntensity = *ESPaDOnSBeam.getStokesParameter(StokesI);
        }
            break;
        case ZeemanLine: {
            SpectrumVariables.StokesParameter = StokesV;
            
            operaStokesVector ZeemanLineStokesVector = GenerateZeemanLineStokesVector(SpectrumVariables,ZeemanLineStokesVector);
            operaStokesVector ESPaDOnSBeam = ESPaDOnSBeamCreator(ZeemanLineStokesVector, ESPaDOnSBeam, SpectrumVariables.Rhomb1Angle, SpectrumVariables.Rhomb2Angle, SpectrumVariables.Beam);

			ESPaDOnSBeamIntensity = *ESPaDOnSBeam.getStokesParameter(StokesI);
        }            
            break;
        default:
            break;
    }
	return ESPaDOnSBeamIntensity;
}

/*
 * BlackBody
 */

operaFluxVector& BlackBodySimulation(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, operaFluxVector& SpectralRadianceVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    double SpectralRadiance;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SpectralRadiance = (double)PlanckLaw((float)Temperature, (float)Wavelength);
        SpectralRadianceVector.setflux(SpectralRadiance,index);
        Wavelength += WavelengthIncrement;
    }
    
    return SpectralRadianceVector;
}

operaFluxVector& BlackBodySimulationUniformError(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, double ErrorPercentage, operaFluxVector& SpectralRadianceVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    double SpectralRadiance;
    double SpectralRadianceUniformError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SpectralRadiance = (double)PlanckLaw((float)Temperature, (float)Wavelength);
        SpectralRadianceUniformError = (double)operaUniformRand((float)(SpectralRadiance - SpectralRadiance * ErrorPercentage), (float)(SpectralRadiance + SpectralRadiance * ErrorPercentage));
        SpectralRadianceVector.setflux(SpectralRadianceUniformError,index);
        SpectralRadianceVector.setvariance(pow(SpectralRadiance * ErrorPercentage,2),index);
        Wavelength += WavelengthIncrement;
    }
    
    return SpectralRadianceVector;
}

operaFluxVector& BlackBodySimulationGaussianError(double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, double ErrorPercentage, operaFluxVector& SpectralRadianceVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    double SpectralRadiance;
    double SpectralRadianceGaussianError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SpectralRadiance = (double)PlanckLaw((float)Temperature, (float)Wavelength);
        SpectralRadianceGaussianError = (double)operaGaussRand((float)SpectralRadiance, (float)(SpectralRadiance * ErrorPercentage));
        SpectralRadianceVector.setflux(SpectralRadianceGaussianError,index);
        SpectralRadianceVector.setvariance(pow(SpectralRadiance * ErrorPercentage,2),index);
        Wavelength += WavelengthIncrement;
    }
    
    return SpectralRadianceVector;
}

operaStokesVector& GenerateBlackBodyStokesVector(unsigned length, error_type_t ErrorType, double ErrorPercentage, double Temperature, double MinWavelength, double MaxWavelength, double WavelengthIncrement, stokes_parameter_t StokesParameter, double StokesParameterFractionPolarizationOfBlackBody, operaStokesVector& BlackBodyStokesVector)
{
#ifdef RANGE_CHECK
    if (StokesParameter >= 4) {
        throw operaException("operaSpectrumSimulation: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    operaFluxVector SpectralRadianceVector(length);
    
    if (ErrorType == NoError)
        SpectralRadianceVector = BlackBodySimulation(Temperature, MinWavelength, MaxWavelength, WavelengthIncrement, SpectralRadianceVector);
    else if (ErrorType == UniformError)
        SpectralRadianceVector = BlackBodySimulationUniformError(Temperature, MinWavelength, MaxWavelength, WavelengthIncrement, ErrorPercentage, SpectralRadianceVector);
    else if (ErrorType == GaussianError)
        SpectralRadianceVector = BlackBodySimulationGaussianError(Temperature, MinWavelength, MaxWavelength, WavelengthIncrement, ErrorPercentage, SpectralRadianceVector);
    
    operaFluxVector StokesIVector(length);
    operaFluxVector StokesQVector(length);
    operaFluxVector StokesUVector(length);
    operaFluxVector StokesVVector(length);
    
    StokesIVector = SpectralRadianceVector;
    if (StokesParameter == StokesQ) {
        StokesQVector = SpectralRadianceVector * StokesParameterFractionPolarizationOfBlackBody;
        StokesUVector = 0.0;
        StokesVVector = 0.0;
    }
    else if (StokesParameter == StokesU) {
        StokesQVector = 0.0;
        StokesUVector = SpectralRadianceVector * StokesParameterFractionPolarizationOfBlackBody;
        StokesVVector = 0.0;
    }
    else if (StokesParameter == StokesV) {
        StokesQVector = 0.0;
        StokesUVector = 0.0;
        StokesVVector = SpectralRadianceVector * StokesParameterFractionPolarizationOfBlackBody;
    }
    
    BlackBodyStokesVector.setStokesParameters(StokesIVector, StokesQVector, StokesUVector, StokesVVector);
    
    return BlackBodyStokesVector;
}

/*
 * GalileanMoons
 */

operaFluxVector& GaussianPolarizationSimulation(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, operaFluxVector& SimulatedPolarizationVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    Gaussian SimulatedPolarization(NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector);
    
    double SimulatedPolarizationValue;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SimulatedPolarizationValue = SimulatedPolarization.EvaluateGaussian(Wavelength);
        SimulatedPolarizationVector.setflux(SimulatedPolarizationValue,index);
        Wavelength += WavelengthIncrement;
    }
    
    return SimulatedPolarizationVector;
}

operaFluxVector& GaussianPolarizationSimulationUniformError(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, double ErrorPercentage, operaFluxVector& SimulatedPolarizationVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    Gaussian SimulatedPolarization(NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector);
    
    double SimulatedPolarizationValue;
    double SimulatedPolarizationValueUniformError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SimulatedPolarizationValue = SimulatedPolarization.EvaluateGaussian(Wavelength);
        SimulatedPolarizationValueUniformError = (double)operaUniformRand((float)(SimulatedPolarizationValue - SimulatedPolarizationValue * ErrorPercentage), (float)(SimulatedPolarizationValue + SimulatedPolarizationValue * ErrorPercentage));
        SimulatedPolarizationVector.setflux(SimulatedPolarizationValueUniformError,index);
        SimulatedPolarizationVector.setvariance(pow(SimulatedPolarizationValue * ErrorPercentage,2),index);

        Wavelength += WavelengthIncrement;
    }
    
    return SimulatedPolarizationVector;
}

operaFluxVector& GaussianPolarizationSimulationGaussianError(double MinWavelength, double MaxWavelength, double WavelengthIncrement, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, double ErrorPercentage, operaFluxVector& SimulatedPolarizationVector)
{
    unsigned length = (unsigned)( (MaxWavelength - MinWavelength) / WavelengthIncrement );
    double Wavelength = MinWavelength;
    
    Gaussian SimulatedPolarization(NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector);
    
    double SimulatedPolarizationValue;
    double SimulatedPolarizationValueGaussianError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SimulatedPolarizationValue = SimulatedPolarization.EvaluateGaussian(Wavelength);
        SimulatedPolarizationValueGaussianError = (double)operaGaussRand((float)SimulatedPolarizationValue, (float)(SimulatedPolarizationValue * ErrorPercentage));
        SimulatedPolarizationVector.setflux(SimulatedPolarizationValueGaussianError,index);
        SimulatedPolarizationVector.setvariance(pow(SimulatedPolarizationValue * ErrorPercentage,2),index);
        Wavelength += WavelengthIncrement;
    }
    
    return SimulatedPolarizationVector;
}
operaFluxVector& IntroduceUniformError(operaFluxVector& SimulatedPolarizationVector, double ErrorPercentage)
{
    unsigned length = SimulatedPolarizationVector.getlength();
    
    double SimulatedPolarizationValue;
    double SimulatedPolarizationValueGaussianError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SimulatedPolarizationValue = SimulatedPolarizationVector.getflux(index);
        SimulatedPolarizationValueGaussianError = (double)operaUniformRand((float)(SimulatedPolarizationValue - SimulatedPolarizationValue * ErrorPercentage), (float)(SimulatedPolarizationValue + SimulatedPolarizationValue * ErrorPercentage));
        SimulatedPolarizationVector.setflux(SimulatedPolarizationValueGaussianError,index);
        SimulatedPolarizationVector.setvariance(pow(SimulatedPolarizationValue * ErrorPercentage,2),index);
    }
    
    return SimulatedPolarizationVector;
}

operaFluxVector& IntroduceGaussianError(operaFluxVector& SimulatedPolarizationVector, double ErrorPercentage)
{
    unsigned length = SimulatedPolarizationVector.getlength();
    
    double SimulatedPolarizationValue;
    double SimulatedPolarizationValueGaussianError;
    
    for (unsigned index = 0 ; index < length ; index++) {
        SimulatedPolarizationValue = SimulatedPolarizationVector.getflux(index);
        SimulatedPolarizationValueGaussianError = (double)operaGaussRand((float)SimulatedPolarizationValue, (float)(SimulatedPolarizationValue * ErrorPercentage));
        SimulatedPolarizationVector.setflux(SimulatedPolarizationValueGaussianError,index);
        SimulatedPolarizationVector.setvariance(pow(SimulatedPolarizationValue * ErrorPercentage,2),index);
    }
    
    return SimulatedPolarizationVector;
}



operaStokesVector& GenerateGalileanMoonStokesVector(unsigned length, error_type_t ErrorType, double ErrorPercentage, double MinWavelength, double MaxWavelength, double WavelengthIncrement, operaFluxVector& SunSpectrum, double ReflectedJupiterFractionOfSunSpectrum, stokes_parameter_t StokesParameter, unsigned NumberOfPolarizationPeaks, double* AmplitudeOfPolarizationVector, double* SigmaOfPolarizationVector, double* CenterOfPolarizationVector, operaStokesVector& GalileanMoonStokesVector)
{
#ifdef RANGE_CHECK
    if (StokesParameter >= 4) {
        throw operaException("operaSpectrumSimulation: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    operaFluxVector SimulatedJupiterSpectrum(length);
    SimulatedJupiterSpectrum = SunSpectrum * ReflectedJupiterFractionOfSunSpectrum;
    
    operaFluxVector SimulatedPolarizationVector(length);

    if (ErrorType == NoError) {
        SimulatedPolarizationVector = GaussianPolarizationSimulation(MinWavelength, MaxWavelength, WavelengthIncrement, NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector, SimulatedPolarizationVector);
    } else if (ErrorType == UniformError) {
        IntroduceUniformError(SimulatedJupiterSpectrum, ErrorPercentage);
        IntroduceUniformError(SunSpectrum, ErrorPercentage);
        SimulatedPolarizationVector = GaussianPolarizationSimulationUniformError(MinWavelength, MaxWavelength, WavelengthIncrement, NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector, ErrorPercentage, SimulatedPolarizationVector);
    } else if (ErrorType == GaussianError) {
        IntroduceGaussianError(SimulatedJupiterSpectrum, ErrorPercentage); 
        IntroduceGaussianError(SunSpectrum, ErrorPercentage); 
        SimulatedPolarizationVector = GaussianPolarizationSimulationGaussianError(MinWavelength, MaxWavelength, WavelengthIncrement, NumberOfPolarizationPeaks, AmplitudeOfPolarizationVector, SigmaOfPolarizationVector, CenterOfPolarizationVector, ErrorPercentage, SimulatedPolarizationVector);
    }
    operaFluxVector SimulatedGalileanMoonSpectrum(length);
    
    // line below needs to be fixed
    SimulatedGalileanMoonSpectrum = SunSpectrum + (SimulatedJupiterSpectrum - SimulatedPolarizationVector);                   
    
    operaFluxVector StokesIVector(length);
    operaFluxVector StokesQVector(length);
    operaFluxVector StokesUVector(length);
    operaFluxVector StokesVVector(length);
    
    StokesIVector = SimulatedGalileanMoonSpectrum;
    
    //double polarization = 0.10;
    
    if (StokesParameter == StokesQ) {
        // the operator below doesn't work
//        StokesQVector = polarization*(SimulatedJupiterSpectrum - SimulatedPolarizationVector);        
        StokesQVector = (SimulatedJupiterSpectrum - SimulatedPolarizationVector);
        StokesUVector = 0.0;
        StokesVVector = 0.0;
    }
    else if (StokesParameter == StokesU) {
        StokesQVector = 0.0;
        StokesUVector = (SimulatedJupiterSpectrum - SimulatedPolarizationVector);
        StokesVVector = 0.0;
    }
    else if (StokesParameter == StokesV) {
        StokesQVector = 0.0;
        StokesUVector = 0.0;
        StokesVVector = (SimulatedJupiterSpectrum - SimulatedPolarizationVector);
    }
    
    GalileanMoonStokesVector.setStokesParameters(StokesIVector, StokesQVector, StokesUVector, StokesVVector);
    
    return GalileanMoonStokesVector;
}

operaStokesVector& GenerateZeemanLineStokesVector(SpectrumVariables_t SpectrumVariables, operaStokesVector& ZeemanLineStokesVector)
{
#ifdef RANGE_CHECK
    if (SpectrumVariables.StokesParameter >= 4) {
        throw operaException("operaSpectrumSimulation: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    operaFluxVector continuumSpectrum(SpectrumVariables.length);
    continuumSpectrum = 0.5;
    
    unsigned NumberOfPolarizationPeaks1 = 1;
    double* AmplitudeOfPolarizationVector1 = new double[NumberOfPolarizationPeaks1];
    double* SigmaOfPolarizationVector1 = new double[NumberOfPolarizationPeaks1];
    double* CenterOfPolarizationVector1 = new double[NumberOfPolarizationPeaks1];    

    SigmaOfPolarizationVector1[0] =  SpectrumVariables.linewidth/2.0; 
    AmplitudeOfPolarizationVector1[0] = SpectrumVariables.linedepth;
    CenterOfPolarizationVector1[0] = SpectrumVariables.linecenter + 2.0*SigmaOfPolarizationVector1[0]*(SpectrumVariables.ZeemanSplit/2.0);
    
    unsigned NumberOfPolarizationPeaks2 = 1;
    double* AmplitudeOfPolarizationVector2 = new double[NumberOfPolarizationPeaks2];
    double* SigmaOfPolarizationVector2 = new double[NumberOfPolarizationPeaks2];
    double* CenterOfPolarizationVector2 = new double[NumberOfPolarizationPeaks2];      
    
    SigmaOfPolarizationVector2[0] =  SpectrumVariables.linewidth/2.0; 
    AmplitudeOfPolarizationVector2[0] = SpectrumVariables.linedepth;
    CenterOfPolarizationVector2[0] = SpectrumVariables.linecenter - 2.0*SigmaOfPolarizationVector2[0]*(SpectrumVariables.ZeemanSplit/2.0);
        
    operaFluxVector SimulatedNegativePolarizationVector(SpectrumVariables.length);
    operaFluxVector SimulatedPositivePolarizationVector(SpectrumVariables.length);
    
    if (SpectrumVariables.ErrorType == NoError){
        SimulatedNegativePolarizationVector = GaussianPolarizationSimulation(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks1, AmplitudeOfPolarizationVector1, SigmaOfPolarizationVector1, CenterOfPolarizationVector1, SimulatedNegativePolarizationVector);
        SimulatedPositivePolarizationVector = GaussianPolarizationSimulation(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks2, AmplitudeOfPolarizationVector2, SigmaOfPolarizationVector2, CenterOfPolarizationVector2, SimulatedPositivePolarizationVector);        
    } else if (SpectrumVariables.ErrorType == UniformError) {
        IntroduceUniformError(continuumSpectrum, SpectrumVariables.ErrorPercentage);
        SimulatedNegativePolarizationVector = GaussianPolarizationSimulationUniformError(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks1, AmplitudeOfPolarizationVector1, SigmaOfPolarizationVector1, CenterOfPolarizationVector1, SpectrumVariables.ErrorPercentage, SimulatedNegativePolarizationVector);
        SimulatedPositivePolarizationVector = GaussianPolarizationSimulationUniformError(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks2, AmplitudeOfPolarizationVector2, SigmaOfPolarizationVector2, CenterOfPolarizationVector2, SpectrumVariables.ErrorPercentage, SimulatedPositivePolarizationVector);
    } else if (SpectrumVariables.ErrorType == GaussianError) {
        IntroduceGaussianError(continuumSpectrum, SpectrumVariables.ErrorPercentage);        
        SimulatedNegativePolarizationVector = GaussianPolarizationSimulationGaussianError(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks1, AmplitudeOfPolarizationVector1, SigmaOfPolarizationVector1, CenterOfPolarizationVector1, SpectrumVariables.ErrorPercentage, SimulatedNegativePolarizationVector);
        SimulatedPositivePolarizationVector = GaussianPolarizationSimulationGaussianError(SpectrumVariables.MinWavelength, SpectrumVariables.MaxWavelength, SpectrumVariables.WavelengthIncrement, NumberOfPolarizationPeaks2, AmplitudeOfPolarizationVector2, SigmaOfPolarizationVector2, CenterOfPolarizationVector2, SpectrumVariables.ErrorPercentage, SimulatedPositivePolarizationVector);        
    }
    operaFluxVector SimulatedZeemanNegativeLineSpectrum(SpectrumVariables.length);
    operaFluxVector SimulatedZeemanPositiveLineSpectrum(SpectrumVariables.length);

    SimulatedZeemanNegativeLineSpectrum = continuumSpectrum - SimulatedNegativePolarizationVector*continuumSpectrum;    
    SimulatedZeemanPositiveLineSpectrum = continuumSpectrum - SimulatedPositivePolarizationVector*continuumSpectrum;            

    operaFluxVector StokesIVector(SpectrumVariables.length);
    operaFluxVector StokesQVector(SpectrumVariables.length);
    operaFluxVector StokesUVector(SpectrumVariables.length);
    operaFluxVector StokesVVector(SpectrumVariables.length);
    
    StokesIVector = SimulatedZeemanNegativeLineSpectrum + SimulatedZeemanPositiveLineSpectrum;
    
    StokesQVector = 0.0;
    
    StokesUVector = 0.0;
    
    StokesVVector = SimulatedZeemanPositiveLineSpectrum - SimulatedZeemanNegativeLineSpectrum;
    
    ZeemanLineStokesVector.setStokesParameters(StokesIVector, StokesQVector, StokesUVector, StokesVVector);
    
    delete[] AmplitudeOfPolarizationVector1;
    delete[] SigmaOfPolarizationVector1;
    delete[] CenterOfPolarizationVector1;
    
    delete[] AmplitudeOfPolarizationVector2;
    delete[] SigmaOfPolarizationVector2;
    delete[] CenterOfPolarizationVector2;    
    
    return ZeemanLineStokesVector;
}


/*
 * BeamCreator
 */

operaStokesVector& ESPaDOnSBeamCreator(operaStokesVector &InputStokesVector, operaStokesVector& OutputStokesVector, double Rhomb1Angle, double Rhomb2Angle, unsigned Beam)
{    
    operaMuellerMatrix WollastonO;      //Ordinary axis | Perpendicular to the optical axis
    operaMuellerMatrix WollastonE;      //Extraordinary axis | Parallel to the optical axis
    
    WollastonO.createRotatedPolarizer(1,0,0,0,0,0);
    WollastonE.createRotatedPolarizer(1,0,M_PI/2,0,0,0);
    
    operaMuellerMatrix WollastonBeam[2] = {WollastonO, WollastonE};
    
    operaMuellerMatrix FresnelHalfWaveRhomb1;   //Half-wave
    operaMuellerMatrix FresnelQuaterWaveRhomb;  //Quarter-wave
    operaMuellerMatrix FresnelHalfWaveRhomb2;   //Half-wave
    
    FresnelHalfWaveRhomb1.createRotatedRetarder(M_PI,0,Rhomb1Angle,0);
    FresnelQuaterWaveRhomb.createRotatedRetarder(M_PI/2,0,0,0);
    FresnelHalfWaveRhomb2.createRotatedRetarder(M_PI,0,Rhomb2Angle,0);
    
    OutputStokesVector = WollastonBeam[Beam] * FresnelHalfWaveRhomb2 * FresnelQuaterWaveRhomb * FresnelHalfWaveRhomb1 * InputStokesVector;
    
    return OutputStokesVector;
}

