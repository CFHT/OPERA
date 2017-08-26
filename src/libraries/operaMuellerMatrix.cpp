/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaMuellerMatrix
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

#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaMuellerMatrix.h"

/*!
 * \file operaMuellerMatrix.cpp
 * \brief This file holds the implementation of the class operaMuellerMatrix.
 */

/*
 * \author Andre Venne
 * \brief This class encapsulates the Mueller matrix.
 * \sa class operaStokesVector, class operaPolarimetry
 * 
 * This class holds all the parameters and their variance to create a Mueller matrix. It also holds the Mueller matrix itself and its variance matrix.
 * The variances are propagated as followed :
 * 
 * F = F(a,b)
 * DF = Pow(dF/da,2) * Da + Pow(dF/db,2) *Db
 * 
 * where DF is the resulting variance, Da and Db are the variance of the variables a and b, dF/da and dF/db are the partial derivatives of F.
 * The variables are supposed uncorrelated.
 */

/*
 * Constructors / Destructors
 */

/*
 * \brief Basic operaMuellerMatrix constructor.
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaMuellerMatrix::operaMuellerMatrix(bool Istemp) :
p(0.0),
alpha(0.0),
phi(0.0),
theta(0.0),
pVariance(0.0),
alphaVariance(0.0),
phiVariance(0.0),
thetaVariance(0.0)
{
    istemp = Istemp;
    
    for (unsigned row = 0 ; row < 4 ; row++) {
        for (unsigned column = 0 ; column < 4 ; column++) {
            muellerMatrix[row][column] = 0;
            varianceMatrix[row][column] = 0;
        }
    }
}

/*
 * \brief Basic operaMuellerMatrix destructor.
 * \return void
 */
operaMuellerMatrix::~operaMuellerMatrix()
{
    
}

/*
 * Creators
 */

/*
 * \brief Creates a rotated polarizer.
 * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotated polarizer.
 * \param P A double value, 0 <= P <= 1
 * \param PVariance A double value
 * \param Alpha A double value, 0 <= Alpha <= 90 deg
 * \param AlphaVariance A double value
 * \param Theta A double value
 * \param ThetaVariance A double value
 * \return void
 */
void operaMuellerMatrix::createRotatedPolarizer(double P, double PVariance, double Alpha, double AlphaVariance, double Theta, double ThetaVariance)     // for a Rotated Polarizer
{
    p = P;
    double pSquare = pow(P,2);
    alpha = Alpha;
    phi = 0.0;
    theta = Theta;
    
    pVariance = PVariance;
    double pSquareVariance = 2.0 * P * PVariance;
    alphaVariance = AlphaVariance;
    phiVariance = 0.0;
    thetaVariance = ThetaVariance;
    
    /*
     * muellerMatrix =    [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    muellerMatrix[0][0] = pSquare / 2.0 * 1.0;                                                                          // 00
    muellerMatrix[0][1] = pSquare / 2.0 * cos(2.0 * alpha) * cos(2.0 * theta);                                          // 01
    muellerMatrix[0][2] = pSquare / 2.0 * cos(2.0 * alpha) * sin(2.0 * theta);                                          // 02
    muellerMatrix[0][3] = 0.0;                                                                                          // 03
    
    muellerMatrix[1][0] = pSquare / 2.0 * cos(2.0 * alpha) * cos(2.0 * theta);                                          // 10
    muellerMatrix[1][1] = pSquare / 2.0 * (pow(cos(2.0 * theta),2) + sin(2.0 * alpha) * pow(sin(2.0 * theta),2));       // 11
    muellerMatrix[1][2] = pSquare / 2.0 * (1.0 - sin(2.0 * alpha)) * sin(2.0 * theta) * cos(2.0 * theta);               // 12
    muellerMatrix[1][3] = 0.0;                                                                                          // 13
    
    muellerMatrix[2][0] = pSquare / 2.0 * cos(2.0 * alpha) * sin(2.0 * theta);                                          // 20
    muellerMatrix[2][1] = pSquare / 2.0 * (1.0 - sin(2.0 * alpha)) * sin(2.0 * theta) * cos(2.0 * theta);               // 21
    muellerMatrix[2][2] = pSquare / 2.0 * (pow(sin(2.0 * theta),2) + sin(2.0 * alpha) * pow(cos(2.0 * theta),2));       // 22
    muellerMatrix[2][3] = 0.0;                                                                                          // 23
    
    muellerMatrix[3][0] = 0.0;                                                                                          // 30
    muellerMatrix[3][1] = 0.0;                                                                                          // 31
    muellerMatrix[3][2] = 0.0;                                                                                          // 32
    muellerMatrix[3][3] = pSquare / 2.0 * sin(2.0 * alpha);                                                             // 33
    
    
    /*
     * varianceMatrix =   [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    varianceMatrix[0][0] = pow(1.0 / 2.0 * 1.0,2) * pSquareVariance;                                                    // 00
    varianceMatrix[0][1] = pow(1.0 / 2.0 * cos(2.0 * alpha) * cos(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * sin(2.0 * alpha) * cos(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * cos(2.0 * alpha) * sin(2.0 * theta),2) * thetaVariance;                        // 01
    varianceMatrix[0][2] = pow(1.0 / 2.0 * cos(2.0 * alpha) * sin(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * sin(2.0 * alpha) * sin(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * cos(2.0 * alpha) * cos(2.0 * theta),2) * thetaVariance;                        // 02
    varianceMatrix[0][3] = 0.0;                                                                                         // 03
    
    varianceMatrix[1][0] = pow(1.0 / 2.0 * cos(2.0 * alpha) * cos(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * sin(2.0 * alpha) * cos(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * cos(2.0 * alpha) * sin(2.0 * theta),2) * thetaVariance;                        // 10
    varianceMatrix[1][1] = pow(1.0 / 2.0 * (pow(cos(2.0 * theta),2) + sin(2.0 * alpha) * pow(sin(2.0 * theta),2)),2) * pSquareVariance
                         + pow(pSquare * cos(2.0 * alpha) * pow(sin(2.0 * theta),2),2) * alphaVariance
                         + pow(pSquare * (sin(2.0 * alpha) - 1.0) * sin(4.0 * theta),2) * thetaVariance;                // 11
    varianceMatrix[1][2] = pow(1.0 / 2.0 * (1.0 - sin(2.0 * alpha)) * sin(2.0 * theta) * cos(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * cos(2.0 * alpha) * sin(2.0 * theta) * cos(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * (1.0 - sin(2.0 * alpha)) * cos(4.0 * theta),2) * thetaVariance;                // 12
    varianceMatrix[1][3] = 0.0;                                                                                         // 13
    
    varianceMatrix[2][0] = pow(1.0 / 2.0 * cos(2.0 * alpha) * sin(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * sin(2.0 * alpha) * sin(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * cos(2.0 * alpha) * cos(2.0 * theta),2) * thetaVariance;                        // 20
    varianceMatrix[2][1] = pow(1.0 / 2.0 * (1.0 - sin(2.0 * alpha)) * sin(2.0 * theta) * cos(2.0 * theta),2) * pSquareVariance
                         + pow(pSquare * cos(2.0 * alpha) * sin(2.0 * theta) * cos(2.0 * theta),2) * alphaVariance
                         + pow(pSquare * (1.0 - sin(2.0 * alpha)) * cos(4.0 * theta),2) * thetaVariance;                // 21
    varianceMatrix[2][2] = pow(1.0 / 2.0 * (pow(cos(2.0 * theta),2) + sin(2.0 * alpha) * pow(sin(2.0 * theta),2)),2) * pSquareVariance
                         + pow(pSquare * cos(2.0 * alpha) * pow(sin(2.0 * theta),2),2) * alphaVariance
                         + pow(pSquare * (1.0 - sin(2.0 * alpha)) * sin(4.0 * theta),2) * thetaVariance;                // 22
    varianceMatrix[2][3] = 0.0;                                                                                         // 23
    
    varianceMatrix[3][0] = 0.0;                                                                                         // 30
    varianceMatrix[3][1] = 0.0;                                                                                         // 31
    varianceMatrix[3][2] = 0.0;                                                                                         // 32
    varianceMatrix[3][3] = pow(1.0 / 2.0 * sin(2.0 * alpha),2) * pSquareVariance
                         + pow(pSquare * cos(2.0 * alpha),2) * alphaVariance;                                           // 33
}

/*
 * \brief Creates a rotated retarder.
 * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotated retarder.
 * \param Phi A double value
 * \param PhiVariance A double value
 * \param Theta A double value
 * \param ThetaVariance A double value
 * \return void
 */
void operaMuellerMatrix::createRotatedRetarder(double Phi, double PhiVariance, double Theta, double ThetaVariance)      // for a Rotated Retarder
{
    p = 0.0;
    alpha = 0.0;
    phi = Phi;
    theta = Theta;
    
    pVariance = 0.0;
    alphaVariance = 0.0;
    phiVariance = PhiVariance;
    thetaVariance = ThetaVariance;
    
    /*
     * muellerMatrix =    [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    muellerMatrix[0][0] = 1.0;                                                                                          // 00
    muellerMatrix[0][1] = 0.0;                                                                                          // 01
    muellerMatrix[0][2] = 0.0;                                                                                          // 02
    muellerMatrix[0][3] = 0.0;                                                                                          // 03
    
    muellerMatrix[1][0] = 0.0;                                                                                          // 10
    muellerMatrix[1][1] = pow(cos(2.0 * theta),2) + cos(phi) * pow(sin(2.0 * theta),2);                                 // 11
    muellerMatrix[1][2] = (1.0 - cos(phi)) * sin(2.0 * theta) * cos(2.0 * theta);                                       // 12
    muellerMatrix[1][3] = sin(phi) * sin(2.0 * theta);                                                                  // 13
                           
    muellerMatrix[2][0] = 0.0;                                                                                          // 20
    muellerMatrix[2][1] = (1.0 - cos(phi)) * sin(2.0 * theta) * cos(2.0 * theta);                                       // 21
    muellerMatrix[2][2] = pow(sin(2.0 * theta),2) + cos(phi) * pow(cos(2.0 * theta),2);                                 // 22
    muellerMatrix[2][3] = -sin(phi) * cos(2.0 * theta);                                                                 // 23
                                                  
    muellerMatrix[3][0] = 0.0;                                                                                          // 30
    muellerMatrix[3][1] = -sin(phi) * sin(2.0 * theta);                                                                 // 31
    muellerMatrix[3][2] = sin(phi) * cos(2.0 * theta);                                                                  // 32
    muellerMatrix[3][3] = cos(phi);                                                                                     // 33
    
    
    /*
     * varianceMatrix =   [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    varianceMatrix[0][0] = 0.0;                                                                                         // 00
    varianceMatrix[0][1] = 0.0;                                                                                         // 01
    varianceMatrix[0][2] = 0.0;                                                                                         // 02
    varianceMatrix[0][3] = 0.0;                                                                                         // 03
    
    varianceMatrix[1][0] = 0.0;                                                                                         // 10
    varianceMatrix[1][1] = pow(sin(phi) * pow(sin(2.0 * theta),2),2) * phiVariance
                         + pow(2.0 * (cos(phi) - 1.0) * sin(4.0 * theta),2) * thetaVariance;                            // 11
    varianceMatrix[1][2] = pow(sin(phi) * sin(2.0 * theta) * cos(2.0 * theta),2) * phiVariance
                         + pow(2.0 * (1.0 - cos(phi)) * cos(4.0 * theta),2) * thetaVariance;                            // 12
    varianceMatrix[1][3] = pow(cos(phi) * sin(2.0 * theta),2) * phiVariance
                         + pow(2.0 * sin(phi) * cos(2.0 * theta),2) * thetaVariance;                                    // 13
    
    varianceMatrix[2][0] = 0.0;                                                                                         // 20
    varianceMatrix[2][1] = pow(sin(phi) * sin(2.0 * theta) * cos(2.0 * theta),2) * phiVariance
                         + pow(2.0 * (1.0 - cos(phi)) * cos(4.0 * theta),2) * thetaVariance;                            // 21
    varianceMatrix[2][2] = pow(sin(phi) * pow(cos(2.0 * theta),2),2) * phiVariance
                         + pow(2.0 * (1.0 - cos(phi)) * sin(4.0 * theta),2) * thetaVariance;                            // 22
    varianceMatrix[2][3] = pow(cos(phi) * cos(2.0 * theta),2) * phiVariance
                         + pow(2.0 * sin(phi) * sin(2.0 * theta),2) * thetaVariance;                                    // 23
    
    varianceMatrix[3][0] = 0.0;                                                                                         // 30
    varianceMatrix[3][1] = pow(cos(phi) * sin(2.0 * theta),2) * phiVariance
                         + pow(2.0 * sin(phi) * cos(2.0 * theta),2) * thetaVariance;                                    // 31
    varianceMatrix[3][2] = pow(cos(phi) * cos(2.0 * theta),2) * phiVariance
                         + pow(2.0 * sin(phi) * sin(2.0 * theta),2) * thetaVariance;                                    // 32
    varianceMatrix[3][3] = pow(sin(phi),2) * phiVariance;                                                               // 33
}

/*
 * \brief Creates a rotator.
 * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents a rotator.
 * \param Theta A double value
 * \param ThetaVariance A double value
 * \return void
 */
void operaMuellerMatrix::createRotator(double Theta, double ThetaVariance)                                              // for a Rotator
{
    p = 0.0;
    alpha = 0.0;
    phi = 0.0;
    theta = Theta;
    
    pVariance = 0.0;
    alphaVariance = 0.0;
    phiVariance = 0.0;
    thetaVariance = ThetaVariance;
    
    /*
     * muellerMatrix =    [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    muellerMatrix[0][0] = 1.0;                                                                                          // 00
    muellerMatrix[0][1] = 0.0;                                                                                          // 01
    muellerMatrix[0][2] = 0.0;                                                                                          // 02
    muellerMatrix[0][3] = 0.0;                                                                                          // 03
    
    muellerMatrix[1][0] = 0.0;                                                                                          // 10
    muellerMatrix[1][1] = cos(2.0 * theta);                                                                             // 11
    muellerMatrix[1][2] = sin(2.0 * theta);                                                                             // 12
    muellerMatrix[1][3] = 0.0;                                                                                          // 13
    
    muellerMatrix[2][0] = 0.0;                                                                                          // 20
    muellerMatrix[2][1] = -sin(2.0 * theta);                                                                            // 21
    muellerMatrix[2][2] = cos(2.0 * theta);                                                                             // 22
    muellerMatrix[2][3] = 0.0;                                                                                          // 23
    
    muellerMatrix[3][0] = 0.0;                                                                                          // 30
    muellerMatrix[3][1] = 0.0;                                                                                          // 31
    muellerMatrix[3][2] = 0.0;                                                                                          // 32
    muellerMatrix[3][3] = 1.0;                                                                                          // 33
    
    
    /*
     * varianceMatrix =   [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    varianceMatrix[0][0] = 0.0;                                                                                         // 00
    varianceMatrix[0][1] = 0.0;                                                                                         // 01
    varianceMatrix[0][2] = 0.0;                                                                                         // 02
    varianceMatrix[0][3] = 0.0;                                                                                         // 03
    
    varianceMatrix[1][0] = 0.0;                                                                                         // 10
    varianceMatrix[1][1] = pow(2.0 * sin(2.0 * theta),2) * thetaVariance;                                               // 11
    varianceMatrix[1][2] = pow(2.0 * cos(2.0 * theta),2) * thetaVariance;                                               // 12
    varianceMatrix[1][3] = 0.0;                                                                                         // 13
    
    varianceMatrix[2][0] = 0.0;                                                                                         // 20
    varianceMatrix[2][1] = pow(2.0 * cos(2.0 * theta),2) * thetaVariance;                                               // 21
    varianceMatrix[2][2] = pow(2.0 * sin(2.0 * theta),2) * thetaVariance;                                               // 22
    varianceMatrix[2][3] = 0.0;                                                                                         // 23
    
    varianceMatrix[3][0] = 0.0;                                                                                         // 30
    varianceMatrix[3][1] = 0.0;                                                                                         // 31
    varianceMatrix[3][2] = 0.0;                                                                                         // 32
    varianceMatrix[3][3] = 0.0;                                                                                         // 33
}

/*
 * \brief Creates an attenuator.
 * \details A function that sets the 16 elements of the Mueller matrix and its variance matrix so that it represents an attenuator.
 * \param P A double value, 0 <= P <= 1
 * \param PVariance A double value
 * \return void
 */
void operaMuellerMatrix::createAttenuator(double P, double PVariance)                                                   // for an Attenuator
{
    p = P;
    alpha = 0.0;
    phi = 0.0;
    theta = 0.0;
    
    pVariance = PVariance;
    alphaVariance = 0.0;
    phiVariance = 0.0;
    thetaVariance = 0.0;
    
    /*
     * muellerMatrix =    [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    muellerMatrix[0][0] = p;                                                                                            // 00
    muellerMatrix[0][1] = 0.0;                                                                                          // 01
    muellerMatrix[0][2] = 0.0;                                                                                          // 02
    muellerMatrix[0][3] = 0.0;                                                                                          // 03
    
    muellerMatrix[1][0] = 0.0;                                                                                          // 10
    muellerMatrix[1][1] = p;                                                                                            // 11
    muellerMatrix[1][2] = 0.0;                                                                                          // 12
    muellerMatrix[1][3] = 0.0;                                                                                          // 13
    
    muellerMatrix[2][0] = 0.0;                                                                                          // 20
    muellerMatrix[2][1] = 0.0;                                                                                          // 21
    muellerMatrix[2][2] = p;                                                                                            // 22
    muellerMatrix[2][3] = 0.0;                                                                                          // 23
    
    muellerMatrix[3][0] = 0.0;                                                                                          // 30
    muellerMatrix[3][1] = 0.0;                                                                                          // 31
    muellerMatrix[3][2] = 0.0;                                                                                          // 32
    muellerMatrix[3][3] = p;                                                                                            // 33
    
    
    /*
     * varianceMatrix =   [0] [1] [2] [3]
     *                [0] 00  01  02  03
     *                [1] 10  11  12  13
     *                [2] 20  21  22  23
     *                [3] 30  31  32  33
     */
    
    varianceMatrix[0][0] = pVariance;                                                                                   // 00
    varianceMatrix[0][1] = 0.0;                                                                                         // 01
    varianceMatrix[0][2] = 0.0;                                                                                         // 02
    varianceMatrix[0][3] = 0.0;                                                                                         // 03
    
    varianceMatrix[1][0] = 0.0;                                                                                         // 10
    varianceMatrix[1][1] = pVariance;                                                                                   // 11
    varianceMatrix[1][2] = 0.0;                                                                                         // 12
    varianceMatrix[1][3] = 0.0;                                                                                         // 13
    
    varianceMatrix[2][0] = 0.0;                                                                                         // 20
    varianceMatrix[2][1] = 0.0;                                                                                         // 21
    varianceMatrix[2][2] = pVariance;                                                                                   // 22
    varianceMatrix[2][3] = 0.0;                                                                                         // 23
    
    varianceMatrix[3][0] = 0.0;                                                                                         // 30
    varianceMatrix[3][1] = 0.0;                                                                                         // 31
    varianceMatrix[3][2] = 0.0;                                                                                         // 32
    varianceMatrix[3][3] = pVariance;                                                                                   // 33
}

/*
 * Getters/Setters
 */

/*
 * \brief Sets the amplitude of attenuation.
 * \details A function that sets the amplitude of attenuation, which is defined as: px^2 + py^2 = p^2.
 * \param P A double value, 0 <= P <= 1
 * \return void
 */
void operaMuellerMatrix::setP(double P)
{
    p = P;
}

/*
 * \brief Gets the amplitude of attenuation.
 * \details A function that gets the amplitude of attenuation, which is defined as: px^2 + py^2 = p^2.
 * \return A double value
 */
double operaMuellerMatrix::getP(void)
{
    return p;
}

/*
 * \brief Sets the variance of the amplitude of attenuation.
 * \details A function that sets the variance of the amplitude of attenuation.
 * \param PVariance A double value
 * \return void
 */
void operaMuellerMatrix::setPVariance(double PVariance)
{
    pVariance = PVariance;
}

/*
 * \brief Gets the variance of the amplitude of attenuation.
 * \details A function that gets the variance of the amplitude of attenuation.
 * \return A double value
 */
double operaMuellerMatrix::getPVariance(void)
{
    return pVariance;
}

/*
 * \brief Sets the trigonometric angle.
 * \details A function that sets the trigonometric angle, which is defined as: px = p cos(alpha), py = p sin(alpha).
 * \param Alpha A double value, 0 <= Alpha <= 90 deg
 * \return void
 */
void operaMuellerMatrix::setAlpha(double Alpha)
{
    alpha = Alpha;
}

/*
 * \brief Gets the trigonometric angle.
 * \details A function that gets the trigonometric angle, which is defined as: px = p cos(alpha), py = p sin(alpha).
 * \return A double value
 */
double operaMuellerMatrix::getAlpha(void)
{
    return alpha;
}

/*
 * \brief Sets the variance of the trigonometric angle.
 * \details A function that sets the variance of the trigonometric angle.
 * \param AlphaVariance A double value
 * \return void
 */
void operaMuellerMatrix::setAlphaVariance(double AlphaVariance)
{
    alphaVariance = AlphaVariance;
}

/*
 * \brief Gets the variance of the trigonometric angle.
 * \details A function that gets the variance of the trigonometric angle.
 * \return A double value
 */
double operaMuellerMatrix::getAlphaVariance(void)
{
    return alphaVariance;
}

/*
 * \brief Sets the phase shift between orthogonal components.
 * \details A function that sets the phase shift between orthogonal components px and py.
 * \param Phi A double value
 * \return void
 */
void operaMuellerMatrix::setPhi(double Phi)
{
    phi = Phi;
}

/*
 * \brief Gets the phase shift between orthogonal components.
 * \details A function that gets the phase shift between orthogonal components px and py.
 * \return A double value
 */
double operaMuellerMatrix::getPhi(void)
{
    return phi;
}

/*
 * \brief Sets the variance of the phase shift between orthogonal components.
 * \details A function that sets the variance of the phase shift between orthogonal components px and py.
 * \param PhiVariance A double value
 * \return void
 */
void operaMuellerMatrix::setPhiVariance(double PhiVariance)
{
    phiVariance = PhiVariance;
}

/*
 * \brief Gets the variance of the phase shift between orthogonal components.
 * \details A function that gets the variance of the phase shift between orthogonal components px and py.
 * \return A double value
 */
double operaMuellerMatrix::getPhiVariance(void)
{
    return phiVariance;
}

/*
 * \brief Sets the rotation angle of components.
 * \details A function that sets the rotation angle of components px and py.
 * \param Theta A double value
 * \return void
 */
void operaMuellerMatrix::setTheta(double Theta)
{
    theta = Theta;
}

/*
 * \brief Gets the rotation angle of components.
 * \details A function that gets the rotation angle of components px and py.
 * \return A double value
 */
double operaMuellerMatrix::getTheta(void)
{
    return theta;
}

/*
 * \brief Sets the variance of the rotation angle of components.
 * \details A function that sets the variance of the rotation angle of components px and py.
 * \param ThetaVariance A double value
 * \return void
 */
void operaMuellerMatrix::setThetaVariance(double ThetaVariance)
{
    thetaVariance = ThetaVariance;
}

/*
 * \brief Gets the variance of the rotation angle of components.
 * \details A function that gets the variance of the rotation angle of components px and py.
 * \return A double value
 */
double operaMuellerMatrix::getThetaVariance(void)
{
    return thetaVariance;
}

/*
 * \brief Sets the Mueller matrix.
 * \details A function that sets the Mueller matrix from a 4x4 matrix.
 * \param MuellerMatrix A 4x4 matrix of double values
 * \return void
 */
void operaMuellerMatrix::setMuellerMatrix(double MuellerMatrix[4][4])
{
    for (unsigned row=0; row<4; row++) {
        for (unsigned column = 0; column < 4; column++) {
            muellerMatrix[row][column] = MuellerMatrix[row][column];
        }
    }
}

/*
 * \brief Gets a Mueller matrix element.
 * \details A function that gets a Mueller matrix element at a specified row and column.
 * \param indexX An unsigned index to the row position
 * \param indexY An unsigned index to the column position
 * \return A double value
 */
double operaMuellerMatrix::getMuellerMatrixElement(unsigned indexX, unsigned indexY)
{
#ifdef RANGE_CHECK
    if (indexX >= 4) {
        throw operaException("operaMuellerMatrix: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
    if (indexY >= 4) {
        throw operaException("operaMuellerMatrix: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return muellerMatrix[indexX][indexY];
}

/*
 * \brief Sets the variance matrix.
 * \details A function that sets the variance matrix from a 4x4 matrix.
 * \param VarianceMatrix A 4x4 matrix of double values
 * \return void
 */
void operaMuellerMatrix::setVarianceMatrix(double VarianceMatrix[4][4])
{
    for (unsigned row=0; row<4; row++) {
        for (unsigned column = 0; column < 4; column++) {
            varianceMatrix[row][column] = VarianceMatrix[row][column];
        }
    }
}

/*
 * \brief Gets a variance matrix element.
 * \details A function that gets a variance matrix element at a specified row and column.
 * \param indexX An unsigned index to the row position
 * \param indexY An unsigned index to the column position
 * \return A double value
 */
double operaMuellerMatrix::getVarianceMatrixElement(unsigned indexX, unsigned indexY)
{
#ifdef RANGE_CHECK
    if (indexX >= 4) {
        throw operaException("operaMuellerMatrix: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
    if (indexY >= 4) {
        throw operaException("operaMuellerMatrix: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return varianceMatrix[indexX][indexY];
}

/*
 * Matrix Operations
 */

/*
 * \brief Prints the Mueller matrix.
 * \details A function that prints the Mueller matrix to the terminal in a 4x4 matrix representation.
 * \return void
 */
void operaMuellerMatrix::printMuellerMatrix(void)
{
	for (unsigned row=0; row<4; row++) {
        for (unsigned column=0; column<4; column++) {
            printf("%6.3f\t",muellerMatrix[row][column]);
        }
        printf("\n");
    }
}

/*
 * \brief Prints the variance matrix.
 * \details A function that prints the variance matrix to the terminal in a 4x4 matrix representation.
 * \return void
 */
void operaMuellerMatrix::printVarianceMatrix(void)
{
	for (unsigned row=0; row<4; row++) {
        for (unsigned column=0; column<4; column++) {
            printf("%6.3f\t",varianceMatrix[row][column]);
        }
        printf("\n");
    }
}

/*
 * \brief Calculates the determinant of the Mueller matrix.
 * \details A function that calculates the determinant of the Mueller matrix. It also calculates the variance of the determinant using the variance matrix.
 * \return A doubleValue structure.
 */
doubleValue operaMuellerMatrix::matrixDeterminant(void)
{
	doubleValue det = {0.0, 0.0};
    operaMuellerMatrix minorMatrix;     // for a 3x3 matrix
    operaMuellerMatrix subMinorMatrix;  // for a 2x2 matrix
    
    /* To find the 3x3 minor matrix of the 4x4 matrix */
	for (unsigned columnFirstRow = 0 ; columnFirstRow < 4 ; columnFirstRow++) {
        unsigned columnMinorMatrix = 0;
		for (unsigned columnSecondRow = 0 ; columnSecondRow < 4 ; columnSecondRow++) {
            if (columnFirstRow == columnSecondRow)
                continue;
            for (unsigned row = 1 ; row < 4 ; row++) {
                minorMatrix.muellerMatrix[row-1][columnMinorMatrix] = this->muellerMatrix[row][columnSecondRow];
                minorMatrix.varianceMatrix[row-1][columnMinorMatrix] = this->varianceMatrix[row][columnSecondRow];
            }
            columnMinorMatrix++;
        }
        
        /* To find the 2x2 minor matrix of the 3x3 matrix */
        for (unsigned columnMinorFirstRow = 0 ; columnMinorFirstRow < 3 ; columnMinorFirstRow++) {
            unsigned columnSubMinorMatrix = 0;
            for (unsigned columnMinorSecondRow = 0 ; columnMinorSecondRow < 3 ; columnMinorSecondRow++) {
                if (columnMinorFirstRow == columnMinorSecondRow)
                    continue;
                for (unsigned rowMinor = 1 ; rowMinor < 3 ; rowMinor++) {
                    subMinorMatrix.muellerMatrix[rowMinor-1][columnSubMinorMatrix] = minorMatrix.muellerMatrix[rowMinor][columnMinorSecondRow];
                    subMinorMatrix.varianceMatrix[rowMinor-1][columnSubMinorMatrix] = minorMatrix.varianceMatrix[rowMinor][columnMinorSecondRow];
                }
                columnSubMinorMatrix++;
            }
            /* Calculate the determinant of the Mueller matrix */
            det.value += pow(-1.0,1.0 + (columnFirstRow + 1.0)) * this->muellerMatrix[0][columnFirstRow] * pow(-1.0,1.0 + (columnMinorFirstRow + 1.0)) * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * (subMinorMatrix.muellerMatrix[0][0] * subMinorMatrix.muellerMatrix[1][1] - subMinorMatrix.muellerMatrix[0][1] * subMinorMatrix.muellerMatrix[1][0]);
            det.error += pow(minorMatrix.muellerMatrix[0][columnMinorFirstRow] * (subMinorMatrix.muellerMatrix[0][0] * subMinorMatrix.muellerMatrix[1][1] - subMinorMatrix.muellerMatrix[0][1] * subMinorMatrix.muellerMatrix[1][0]),2) * this->varianceMatrix[0][columnFirstRow]
            + pow(this->muellerMatrix[0][columnFirstRow] * (subMinorMatrix.muellerMatrix[0][0] * subMinorMatrix.muellerMatrix[1][1] - subMinorMatrix.muellerMatrix[0][1] * subMinorMatrix.muellerMatrix[1][0]),2) * minorMatrix.varianceMatrix[0][columnMinorFirstRow] 
            + pow(this->muellerMatrix[0][columnFirstRow] * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[1][1],2) * subMinorMatrix.varianceMatrix[0][0] 
            + pow(this->muellerMatrix[0][columnFirstRow] * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[0][0],2) * subMinorMatrix.varianceMatrix[1][1] 
            + pow(this->muellerMatrix[0][columnFirstRow] * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[1][0],2)  * subMinorMatrix.varianceMatrix[0][1] 
            + pow(this->muellerMatrix[0][columnFirstRow] * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[0][1],2) * subMinorMatrix.varianceMatrix[1][0];
        }
    }
    
    return det;
}

/*
 * \brief Calculates the cofactor matrix of the Mueller matrix.
 * \details A function that calculates the cofactor matrix of the Mueller matrix. It also calculates the variance of the cofactor matrix elements using the variance matrix.
 * \param Istemp Optional bool defaults to false
 * \return An operaMuellerMatrix address.
 */
operaMuellerMatrix& operaMuellerMatrix::matrixCofactor(bool Istemp)
{
    istemp = Istemp;
    
    operaMuellerMatrix *cofactorMatrix = new operaMuellerMatrix(istemp);
    operaMuellerMatrix minorMatrix;     // for a 3x3 matrix
    operaMuellerMatrix subMinorMatrix;  // for a 2x2 matrix
    
    /* To find the 3x3 minor matrix of the 4x4 matrix */
    for (unsigned row = 0 ; row < 4 ; row++) {
        for (unsigned column = 0 ; column < 4 ; column++) {
            unsigned columnMinorMatrix = 0;
            for (unsigned columnMuellerMatrix = 0 ; columnMuellerMatrix < 4 ; columnMuellerMatrix++) {
                if (column == columnMuellerMatrix)
                    continue;
                unsigned rowMinorMatrix = 0;
                for (unsigned rowMuellerMatrix = 0 ; rowMuellerMatrix < 4 ; rowMuellerMatrix++) {
                    if (row == rowMuellerMatrix)
                        continue;
                    minorMatrix.muellerMatrix[rowMinorMatrix][columnMinorMatrix] = this->muellerMatrix[rowMuellerMatrix][columnMuellerMatrix];
                    minorMatrix.varianceMatrix[rowMinorMatrix][columnMinorMatrix] = this->varianceMatrix[rowMuellerMatrix][columnMuellerMatrix];
                    rowMinorMatrix++;
                }
                columnMinorMatrix++;
            }
            
            /* To find the 2x2 minor matrix of the 3x3 matrix */
            for (unsigned columnMinorFirstRow = 0 ; columnMinorFirstRow < 3 ; columnMinorFirstRow++) {
                unsigned columnSubMinorMatrix = 0;
                for (unsigned columnMinorSecondRow = 0 ; columnMinorSecondRow < 3 ; columnMinorSecondRow++) {
                    if (columnMinorFirstRow == columnMinorSecondRow)
                        continue;
                    for (unsigned rowSubMinor = 1 ; rowSubMinor < 3 ; rowSubMinor++) {
                        subMinorMatrix.muellerMatrix[rowSubMinor-1][columnSubMinorMatrix] = minorMatrix.muellerMatrix[rowSubMinor][columnMinorSecondRow];
                        subMinorMatrix.varianceMatrix[rowSubMinor-1][columnSubMinorMatrix] = minorMatrix.varianceMatrix[rowSubMinor][columnMinorSecondRow];
                    }
                    columnSubMinorMatrix++;
                }
                /* Calculate the cofactor matrix elements */
                cofactorMatrix->muellerMatrix[row][column] += pow(-1.0,(row + 1.0) + (column + 1.0)) * pow(-1.0,1.0 + (columnMinorFirstRow + 1.0)) * minorMatrix.muellerMatrix[0][columnMinorFirstRow] * (subMinorMatrix.muellerMatrix[0][0] * subMinorMatrix.muellerMatrix[1][1] - subMinorMatrix.muellerMatrix[0][1] * subMinorMatrix.muellerMatrix[1][0]);
                cofactorMatrix->varianceMatrix[row][column] += pow(subMinorMatrix.muellerMatrix[0][0] * subMinorMatrix.muellerMatrix[1][1] - subMinorMatrix.muellerMatrix[0][1] * subMinorMatrix.muellerMatrix[1][0],2) * minorMatrix.varianceMatrix[0][columnMinorFirstRow] + pow(minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[1][1],2) * subMinorMatrix.varianceMatrix[0][0] + pow(minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[0][0],2) * subMinorMatrix.varianceMatrix[1][1] + pow(minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[1][0],2)  * subMinorMatrix.varianceMatrix[0][1] + pow(minorMatrix.muellerMatrix[0][columnMinorFirstRow] * subMinorMatrix.muellerMatrix[0][1],2) * subMinorMatrix.varianceMatrix[1][0];
            }
        }
    }
    
    return *cofactorMatrix;
}

/*
 * \brief Calculates the transpose matrix of the Mueller matrix.
 * \details A function that calculates the transpose matrix of the Mueller matrix. It also calculates the variance of the transpose matrix elements using the variance matrix.
 * \param Istemp Optional bool defaults to false
 * \return An operaMuellerMatrix address.
 */
operaMuellerMatrix& operaMuellerMatrix::matrixTranspose(bool Istemp)
{
    istemp = Istemp;
    
	operaMuellerMatrix *transposeMatrix = new operaMuellerMatrix(istemp);
	
    for (unsigned row = 0 ; row < 4 ; row++) {
        for (unsigned column = 0 ; column < 4 ; column++) {
            transposeMatrix->muellerMatrix[row][column] = this->muellerMatrix[column][row];
            transposeMatrix->varianceMatrix[row][column] = this->varianceMatrix[column][row];
        }
    }
	
	return *transposeMatrix;
}

/*
 * \brief Calculates the adjoint matrix of the Mueller matrix.
 * \details A function that calculates the adjoint matrix of the Mueller matrix. It also calculates the variance of the adjoint matrix elements using the variance matrix.
 * \details It calls matrixCofactor and matrixTranspose.
 * \param Istemp Optional bool defaults to false
 * \return An operaMuellerMatrix address.
 */
operaMuellerMatrix& operaMuellerMatrix::matrixAdjoint(bool Istemp)
{
    istemp = Istemp;
    
    operaMuellerMatrix temporaryMatrix;
    operaMuellerMatrix *adjointMatrix = new operaMuellerMatrix(istemp);
    
    temporaryMatrix = this->matrixCofactor(true);
    *adjointMatrix = temporaryMatrix.matrixTranspose(true);
    
	return *adjointMatrix;
}

/*
 * \brief Calculates the inverse matrix of the Mueller matrix.
 * \details A function that calculates the inverse matrix of the Mueller matrix. It also calculates the variance of the inverse matrix elements using the variance matrix.
 * \details It calls matrixDeterminant and matrixAdjoint.
 * \param Istemp Optional bool defaults to false
 * \return An operaMuellerMatrix address.
 */
operaMuellerMatrix& operaMuellerMatrix::matrixInverse(bool Istemp)
{
    istemp = Istemp;
    
	doubleValue det = this->matrixDeterminant();
    operaMuellerMatrix *inverseMatrix = new operaMuellerMatrix(istemp);
    operaMuellerMatrix *EMPTY = new operaMuellerMatrix(istemp);     // for an empty operaMuellerMatrix in case of a determinant equal to 0
	
	if(det.value == 0) {
		operaPError("operaMuellerMatrix::matrixInverse ", MatrixZeroDeterminant);
		return *EMPTY;
	}
    
    *inverseMatrix = this->matrixAdjoint(true) / det;
    
    if (EMPTY)
        delete EMPTY;
	
	return *inverseMatrix;
}
