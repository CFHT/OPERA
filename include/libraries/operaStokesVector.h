/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaStokesVector
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

#ifndef OPERASTOKESVECTOR_H
#define OPERASTOKESVECTOR_H

#include "libraries/operaFluxVector.h"

/*!
 * \file operaStokesVector.h
 * \brief This file holds the declaration of the class operaStokesVector.
 * \ingroup libraries
 */

/*!
 * \brief Definition of the value of each Stokes parameter.
 */
typedef enum { StokesI=0, StokesQ, StokesU, StokesV } stokes_parameter_t;

/*!
 * \brief Convert a Stokes parameter index value to the corresponding name.
 * \param StokesParameter Which Stokes parameter.
 * \return The name of the Stokes parameter.
 */
string StokesName(stokes_parameter_t StokesParameter);

/*!
 * \author Andre Venne
 * \brief This class encapsulates the Stokes vector.
 * \sa class operaFluxVector, class operaMuellerMatrix, class operaPolarimetry
 * 
 * This class contains 4 operaFluxVectors to store each Stokes parameter or
 * any results with the same format, such as the degree of polarization or
 * the 2 null polarization spectra.
 */
class operaStokesVector {

friend class operaMuellerMatrix;

private:

protected:
    unsigned length;
    std::vector<operaFluxVector> stokesVector; // The 4 Stokes parameters I, Q, U, V
	
public:
	/*!
     * \brief Creates an empty operaStokesVector.
     */
    operaStokesVector();
    
    /*!
     * \brief Creates an operaStokesVector with the specified size.
     * \param Length The number of elements in each operaFluxVector.
     */
	operaStokesVector(unsigned Length);
    
    /*!
     * \brief Creates an operaStokesVector from 4 operaFluxVectors.
     * \param StokesVectorI Stokes I.
     * \param StokesVectorQ Stokes Q.
     * \param StokesVectorU Stokes U.
     * \param StokesVectorV Stokes V.
     */
    operaStokesVector(const operaFluxVector &StokesVectorI, const operaFluxVector &StokesVectorQ, const operaFluxVector &StokesVectorU, const operaFluxVector &StokesVectorV);
    
    /*!
     * \brief Resizes the operaStokesVector.
     * \details Resizes each of the 4 operaFluxVectors.
	 * \param Length The length to resize to.
     */
    void resize(unsigned Length);
    
    /*!
     * \brief Removes all elments outside of the specified range from the operaStokesVector.
     * \param range The index range of elements to keep.
     */
	void trim(operaIndexRange range);
    
    /*!
     * \brief Gets the length of the operaStokesVector.
     * \return The length of the 4 operaFluxVectors.
     */
	unsigned getLength(void) const;
    
    /*!
     * \brief Sets the operaFluxVector for each Stokes parameter.
     * \param StokesIFluxVector Stokes I.
     * \param StokesQFluxVector Stokes Q.
     * \param StokesUFluxVector Stokes U.
     * \param StokesVFluxVector Stokes V.
     */
    void setStokesParameters(const operaFluxVector &StokesIFluxVector, const operaFluxVector &StokesQFluxVector, const operaFluxVector &StokesUFluxVector, const operaFluxVector &StokesVFluxVector);
    
    /*!
     * \brief Gets the operaFluxVector for a Stokes parameter.
     * \param StokesIndex Which Stokes parameter to get.
     * \return An immutable reference to that operaFluxVector.
     */
    const operaFluxVector &getStokesParameter(stokes_parameter_t StokesIndex) const;
	
    /*!
     * \brief Sets the operaFluxVector for a single Stokes parameter.
     * \param StokesIndex Which Stokes parameter to set.
     * \param FluxVector The new operaFluxVector.
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector);
    
    /*!
     * \brief Gets the value of a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter to get.
     * \param index The index of the element.
     * \return The value at that index.
     */
    double getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index) const;
    
    /*!
     * \brief Gets the variance of a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter to get.
     * \param index The index of the element.
     * \return The variance at that index.
     */
    double getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index) const;

    /*!
     * \brief Sets the value and variance of a Stokes parameter element.
     * \param StokesIndex Which Stokes parameter to set.
     * \param StokesValue The value of the element.
     * \param Variance The variance of the element.
     * \param index The index of the element.
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index);
};

#endif
