/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaPolarimetry
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

#ifndef OPERAPOLARIMETRY_H
#define OPERAPOLARIMETRY_H

#include "libraries/operaStokesVector.h"

/*!
 * \file operaPolarimetry.h
 * \brief This file holds the declaration of the class operaPolarimetry.
 * \ingroup libraries
 */

using namespace std;

/*!
 * \brief Definition of the method used to calculate polarization.
 */
typedef enum { Difference=1, Ratio, DifferenceWithBeamSwapped, NewMethod} method_t;

/*!
 * \author Andre Venne
 * \brief This class encapsulates the polarimetry results.
 * \sa class operaStokesVector
 * 
 * This class holds in operaStokesVector classes the 4 Stokes parameters, the 4 associated degrees of polarization
 * and the 2 null polarization spectra for each Stokes parameter.
 */
class operaPolarimetry {
	
private:
	unsigned length;
    method_t method;
    
	operaStokesVector stokesParameter;			// The 4 Stokes parameters I, Q, U, V
    operaStokesVector degreeOfPolarization;		// The 4 degrees of polarization for each Stokes parameters I, Q, U, V
    operaStokesVector continuumRemoved;			// The 4 degrees of polarization after removal of continuum polarization
    operaStokesVector firstNullPolarization;	// The 4 null polarization for the first null spectrum
    operaStokesVector secondNullPolarization;	// The 4 null polarization for the second null spectrum
	operaVector wavelength;						// wavelength in nm 

    bool hasStokesI;
    bool hasStokesQ;
    bool hasStokesU;
    bool hasStokesV;
    bool hasDegreeOfStokesI;
    bool hasDegreeOfStokesQ;
    bool hasDegreeOfStokesU;
    bool hasDegreeOfStokesV;
    bool hasContinuumRemoved;
    bool hasFirstNullPolarization;
    bool hasSecondNullPolarization;
	bool hasWavelength;

public:
	/*!
     * \brief Creates empty set of polarimetry elements.
     * \return void
     */
	operaPolarimetry();
    
    /*!
     * \brief Creates empty set of operaSpectralElements.
     * \param Length The number of polarimetry elements.
     */
	operaPolarimetry(unsigned Length);
    
    /*!
     * \brief Resizes to the specified number of polarimetry elements.
     * \param Length The new number of polarimetry elements.
     * \details Existing elements are preserved up to Length.
     */
	void resize(unsigned Length);
	
	/*!
     * \brief Removes all polarimetry elements outside of the specified range.
     * \param range The index range of elements to keep.
     */
	void trim(operaIndexRange range);
	
    /*!
     * \brief Gets the number of polarimetry elements.
     * \return The number of polarimetry elements.
     */
	unsigned getLength(void) const;

    /*!
     * \brief Get the method used to calculate polarization.
     * \return The method used.
     * \details Supported methods are Difference=1, Ratio=2, DifferenceWithBeamSwapped=3, NewMethod=4.
     */
	method_t getmethod(void) const {return method;}
    
    /*!
     * \brief Set the method used to calculate polarization.
     * \param Method The method used.
     * \details Supported methods are Difference=1, Ratio=2, DifferenceWithBeamSwapped=3, NewMethod=4.
     */
	void setmethod(method_t Method) {method = Method;};
	
	/*!
	 * \brief Gets the wavelength of a polarimetry element.
	 * \param indexElem The index of the element.
	 * \return Wavelength The wavelength at that index.
	 */
	double getwavelength(unsigned indexElem) const;
	
	/*!
	 * \brief Sets a wavelength of a polarimetry element.
	 * \param Wavelength The new wavelength value.
	 * \param indexElem The index of the element.
	 */
	void setwavelength(double Wavelength, unsigned indexElem);
	
    /*!
     * \brief Gets a Stokes parameter.
     * \param StokesIndex Which Stokes paramter.
     * \return An immutable reference to the operaFluxVector of that Stokes parameter.
     */
    const operaFluxVector &getStokesParameter(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Sets a Stokes parameter.
     * \param StokesIndex Which Stokes parameter.
     * \param FluxVector An operaFluxVector to set as that Stokes parameter.
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The value of the element.
     */
    double getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The variance of the element.
     */
    double getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Sets a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter.
     * \param SecondNullPolarizationValue The value to set for the element.
     * \param Variance The variance to set for the element.
     * \param index The index of the element.
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter of the degree of polarization vector.
     * \param StokesIndex Which Stokes paramter.
     * \return An immutable reference to the operaFluxVector of the degree of polarization for that Stokes.
     */
    const operaFluxVector &getDegreeOfPolarization(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Sets a Stokes parameter of the degree of polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param FluxVector An operaFluxVector to set as the degree of polarization for that Stokes.
     */
    void setDegreeOfPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element of the degree of polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The value of the element.
     */
    double getDegreeOfPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element of the degree of polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The variance of the element.
     */
    double getDegreeOfPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Sets a Stokes parameter element of the degree of polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param DegreeOfPolarizationValue The value to set for the element.
     * \param Variance The variance to set for the element.
     * \param index The index of the element.
     */
    void setDegreeOfPolarization(stokes_parameter_t StokesIndex, double DegreeOfPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter of the degree of polarization with continuum polarization removed.
     * \param StokesIndex Which Stokes paramter.
     * \return An immutable reference to the operaFluxVector of the degree of polarization with continuum polarization removed for that Stokes.
     */
    const operaFluxVector &getContinuumRemoved(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Sets a Stokes parameter of the degree of polarization with continuum polarization removed.
     * \param StokesIndex Which Stokes parameter.
     * \param FluxVector An operaFluxVector to set as the degree of polarization with continuum polarization removed for that Stokes.
     */
    void setContinuumRemoved(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element of the degree of polarization with continuum polarization removed.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The value of the element.
     */
    double getContinuumRemovedFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element of the degree of polarization with continuum polarization removed.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The variance of the element.
     */
    double getContinuumRemovedVariance(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Sets a Stokes parameter element of the degree of polarization with continuum polarization removed.
     * \param StokesIndex Which Stokes parameter.
     * \param ContinuumRemovedValue The value to set for the element.
     * \param Variance The variance to set for the element.
     * \param index The index of the element.
     */
    void setContinuumRemoved(stokes_parameter_t StokesIndex, double ContinuumRemovedValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter of the first null polarization vector.
     * \param StokesIndex Which Stokes paramter.
     * \return An immutable reference to the operaFluxVector of the first null polarization for that Stokes.
     */
    const operaFluxVector &getFirstNullPolarization(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Sets a Stokes parameter of the first null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param FluxVector An operaFluxVector to set as the first null polarization for that Stokes.
     */
    void setFirstNullPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element of the first null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The value of the element.
     */
    double getFirstNullPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element of the first null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The variance of the element.
     */
    double getFirstNullPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Sets a Stokes parameter element of the first null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param SecondNullPolarizationValue The value to set for the element.
     * \param Variance The variance to set for the element.
     * \param index The index of the element.
     */
    void setFirstNullPolarization(stokes_parameter_t StokesIndex, double FirstNullPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter of the second null polarization vector.
     * \param StokesIndex Which Stokes paramter.
     * \return An immutable reference to the operaFluxVector of the second null polarization for that Stokes.
     */
    const operaFluxVector &getSecondNullPolarization(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Sets a Stokes parameter of the second null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param FluxVector An operaFluxVector to set as the second null polarization for that Stokes.
     */
    void setSecondNullPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element of the second null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The value of the element.
     */
    double getSecondNullPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element of the second null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param index The index of the element.
     * \return The variance of the element.
     */
    double getSecondNullPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Sets a Stokes parameter element of the second null polarization vector.
     * \param StokesIndex Which Stokes parameter.
     * \param SecondNullPolarizationValue The value to set for the element.
     * \param Variance The variance to set for the element.
     * \param index The index of the element.
     */
    void setSecondNullPolarization(stokes_parameter_t StokesIndex, double SecondNullPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Get the boolean flag indicating if there are values for a given Stokes parameter.
     * \param StokesIndex Which Stokes parameter.
     * \return Value of the flag.
     */
    bool getHasStokes(stokes_parameter_t StokesIndex) const;
    
    /*!
     * \brief Set the boolean flag indicating if there are values for a given Stokes parameter.
     * \param StokesIndex Which Stokes parameter.
     * \param HasDegreeOfStokes Value to set the flag to.
     */
    void setHasStokes(stokes_parameter_t StokesIndex, bool HasStokes);
	
    /*!
     * \brief Get the boolean flag indicating if there are values for the degree of polarization of a given Stokes parameter.
     * \param StokesIndex Which Stokes parameter.
     * \return Value of the flag.
     */
    bool getHasDegreeOfStokes(stokes_parameter_t StokesIndex) const;

    /*!
     * \brief Set the boolean flag indicating if there are values for the degree of polarization of a given Stokes parameter.
     * \param StokesIndex Which Stokes parameter.
     * \param HasDegreeOfStokes Value to set the flag to.
     */
    void setHasDegreeOfStokes(stokes_parameter_t StokesIndex, bool HasDegreeOfStokes);
    
    void copyFROMcontinuumremoved() { degreeOfPolarization = continuumRemoved; }
    
    bool getHasStokesI(void) const { return hasStokesI; }
    bool getHasStokesQ(void) const { return hasStokesQ; }
    bool getHasStokesU(void) const { return hasStokesU; }
    bool getHasStokesV(void) const { return hasStokesV; }
    
    void setHasStokesI(bool HasStokesI) { hasStokesI = HasStokesI; }
    void setHasStokesQ(bool HasStokesQ) { hasStokesQ = HasStokesQ; }
    void setHasStokesU(bool HasStokesU) { hasStokesU = HasStokesU; }
    void setHasStokesV(bool HasStokesV) { hasStokesV = HasStokesV; }
    
    bool getHasDegreeOfStokesI(void) const { return hasDegreeOfStokesI; }
    bool getHasDegreeOfStokesQ(void) const { return hasDegreeOfStokesQ; }
    bool getHasDegreeOfStokesU(void) const { return hasDegreeOfStokesU; }
    bool getHasDegreeOfStokesV(void) const { return hasDegreeOfStokesV; }
    
    void setHasDegreeOfStokesI(bool HasDegreeOfStokesI) { hasDegreeOfStokesI = HasDegreeOfStokesI; }
    void setHasDegreeOfStokesQ(bool HasDegreeOfStokesQ) { hasDegreeOfStokesQ = HasDegreeOfStokesQ; }
    void setHasDegreeOfStokesU(bool HasDegreeOfStokesU) { hasDegreeOfStokesU = HasDegreeOfStokesU; }
    void setHasDegreeOfStokesV(bool HasDegreeOfStokesV) { hasDegreeOfStokesV = HasDegreeOfStokesV; }
    
    bool getHasContinuumRemoved(void) const { return hasContinuumRemoved; }
    void setHasContinuumRemoved(bool HasContinuumRemoved) { hasContinuumRemoved = HasContinuumRemoved; }
    
    bool getHasFirstNullPolarization(void) const { return hasFirstNullPolarization; }
    void setHasFirstNullPolarization(bool HasFirstNullPolarization) { hasFirstNullPolarization = HasFirstNullPolarization; }
    
    bool getHasSecondNullPolarization(void) const { return hasSecondNullPolarization; }
    void setHasSecondNullPolarization(bool HasSecondNullPolarization) { hasSecondNullPolarization = HasSecondNullPolarization; }
    
    bool getHasWavelength() const { return hasWavelength; }
	void setHasWavelength(bool HasWavelength) { hasWavelength = HasWavelength; }
	
    /*!
     * \brief Calculates the polarization.
     * \details A function that calculates the polarization from the degree of polarization and stores it in the Stokes vector.
     * \return void
     */
    void calculatePolarization(void);
    
    /*!
     * \brief Calculates the polarization.
     * \details A function that calculates the polarization from the degree of polarization and stores it in the Stokes vector.
     * \param StokesIndex is the selected Stokes parameter
     * \return void
     */
    void calculatePolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Calculates the degree of polarization.
     * \details A function that calculates the degree of polarization from the polarization and stores it in the degree of polarization vector.
     * \return void
     */
    void calculateDegreeOfPolarization(void);
    
    /*!
     * \brief Calculates the degree of polarization.
     * \details A function that calculates the degree of polarization from the polarization and stores it in the degree of polarization vector.
     * \param StokesIndex is the selected Stokes parameter
     * \return void
     */
    void calculateDegreeOfPolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Calculates the degree of polarization.
     * \details A function that calculates the degree of polarization from the observed ordinary and extraordinary beam fluxes.
     * \details This function accepts 2 or 4 input pairs of fluxes (polarimetric exposures).
     * \param StokesIndex is the selected Stokes parameter
     * \param *iE[4] is a set of flux vectors for the input beams with a given state of polarization (ordinary beams)
     * \param *iA[4] is a set of flux vectors for the input beams with a given orthogonal state of polarization (extra-ordinary beams)
     * \param NumberOfExposures is the number of input exposures (accepts only 2 or 4)
     * \return void
     */
    void calculateDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures);

    /*!
     * \brief Calculates Stokes I and another given Stokes parameter.
     * \details A function that calculates the polarized flux for a given Stokes and the 
     * \details total flux for Stokes I given the observed ordinary and extra-ordinary beam fluxes.
     * \details This function accepts either 2 or 4 input pairs of fluxes (polarimetric exposures).
     * \param StokesIndex is the selected Stokes parameter
     * \param *iE[4] is a set of flux vectors for the input beams with a given state of polarization (ordinary beams)
     * \param *iA[4] is a set of flux vectors for the input beams with a given orthogonal state of polarization (extra-ordinary beams)
     * \param NumberOfExposures is the number of input exposures (accepts only 2 or 4)
     * \return void
     */
    void calculateStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures);
};
#endif
