#ifndef GAINBIASNOISE_H
#define GAINBIASNOISE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: GainBiasNoise
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
 
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

#define MAXAMPS 2

/*!
 * \file GainBiasNoise.h
 */

/*!
 * \author Doug Teeple
 * \brief This class encapsulates espadons gain, bias and noise values.
 * \ingroup libraries
 * \sa class GainBiasNoise
 */

#include "libraries/operaFITSImage.h"	// for DATASEC_t


class GainBiasNoise {
	
private:
	unsigned namps;
	double gains[MAXAMPS];
	double noises[MAXAMPS];
	double gainerrors[MAXAMPS];
	double biases[MAXAMPS];
    
    DATASEC_t ampsDataSec[MAXAMPS];
    
public:
	/*
	 * Constructors / Destructors
	 */
	GainBiasNoise();	
	GainBiasNoise(unsigned Namps);	

	~GainBiasNoise();
	
	/*! 
	 * \sa method double getBias();
	 * \brief returns the Bias of amp
	 */
	double getBias(unsigned amp) const { return biases[amp]; };
	/*! 
	 * \sa method void setBias(unsigned amp, double noise);
	 * \brief sets the Bias of amp
	 */
	void setBias(unsigned amp, double bias) { biases[amp] = bias; };
	/*! 
	 * \sa method double getNoise();
	 * \brief returns the Noise of amp
	 */
	double getNoise(unsigned amp) const { return noises[amp]; };
    
	/*!
	 * \sa method double getNoise();
	 * \brief returns the Noise of amp with respect to pixel x,y
	 */
    double getNoise(unsigned x, unsigned y) const;

	/*!
	 * \sa method void setNoise(unsigned amp, double noise);
	 * \brief sets the Noise of amp
	 */
	void setNoise(unsigned amp, double noise) { noises[amp] = noise; };
	
    /*!
	 * \sa method double getGain();
	 * \brief returns the Gain of amp
	 */
	double getGain(unsigned amp) const { return gains[amp]; };
    
    /*!
     * \sa method double getGain();
     * \brief returns the Gain of amp with respect to pixel x,y
     */
    double getGain(unsigned x, unsigned y) const;
        
	/*! 
	 * \sa method void setGain(unsigned amp, double gain);
	 * \brief sets the Gain of amp
	 */
	void setGain(unsigned amp, double gain){ gains[amp] = gain; };
	/*! 
	 * \sa method double getDatasec();
	 * \brief returns the datasec ot amp
	 */
	void getDatasec(unsigned amp, DATASEC_t &datasec) const;
	/*! 
	 * \sa method void setDatasec(unsigned amp, DATASEC_t &datasec);
	 * \brief sets the datsec of amp
	 */
	void setDatasec(unsigned amp, DATASEC_t &datasec);
	/*! 
	 * \sa method double getGainError();
	 * \brief returns the Gain Error of amp
	 */
	double getGainError(unsigned amp) const { return gainerrors[amp]; };
	/*! 
	 * \sa method void setGain(unsigned amp, double gain);
	 * \brief sets the Gain Error of amp
	 */
	void setGainError(unsigned amp, double gainerror){ gainerrors[amp] = gainerror; };
	/*! 
	 * \sa method unsigned getAmps(void);
	 * \brief sets the number of amps
	 */
	unsigned getAmps(void) const { return namps; };
	/*! 
	 * \sa method void setAmps(unsigned Amps);
	 * \brief sets the number of amps
	 */
	void setAmps(unsigned Amps) { namps = Amps; };

};
#endif
