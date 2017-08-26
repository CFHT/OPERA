/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectrograph
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

#include "operaError.h"
#include "globaldefines.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectrograph.h"    // for operaSpectrograph
#include "libraries/operaFit.h"             // for operaFitSplineDouble

/*!
 * operaSpectrograph
 * \author Eder Martioli
 * \brief Encapsulation of spectrograph instrumental information. 
 * \file operaSpectrograph.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectrograph
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

/*
 * Constructors
 */

operaSpectrograph::operaSpectrograph() :
InjectionHoleDiameter(0.0),
OpticalFiber(undefinedFiber),
fiberLength(0.0),
fiberCoreDiameter(0.0),
numberOfInputFibers(1),
numberOfSlices(1),
spectralResolution(0.0),
SpectrographCCD(undefinedCCD),
EspadonsCCDReadoutSpeed(undefinedReadoutMode),
EspadonsInstrumentMode(undefinedInstrumentMode),
x_pixelsize(0.0),
y_pixelsize(0.0)
{
}

/*
 * Destructor
 */

operaSpectrograph::~operaSpectrograph() {
  
}

/*
 * Methods for managing data
 */

void operaSpectrograph::setInjectionHoleDiameter(double injectionHoleDiameter) {
    InjectionHoleDiameter = injectionHoleDiameter;
}

double operaSpectrograph::getInjectionHoleDiameter(void) {
    return InjectionHoleDiameter;
}

void operaSpectrograph::setOpticalFiber(opticalFiber_t opticalFiberType) {
    OpticalFiber = opticalFiberType;
}

opticalFiber_t operaSpectrograph::getOpticalFiber(void) {
    return OpticalFiber;
}

void operaSpectrograph::setfiberLength(double FiberLength) {
    fiberLength = FiberLength;
}

double operaSpectrograph::getfiberLength(void) {
    return fiberLength;
}

void operaSpectrograph::setfiberCoreDiameter(double FiberCoreDiameter) {
    fiberCoreDiameter = FiberCoreDiameter;
}

double operaSpectrograph::getfiberCoreDiameter(void) {
    return fiberCoreDiameter;
}

void operaSpectrograph::setnumberOfInputFibers(unsigned NumberOfInputFibers) {
    numberOfInputFibers = NumberOfInputFibers;
}

unsigned operaSpectrograph::getnumberOfInputFibers(void) {
    return numberOfInputFibers;
}

void operaSpectrograph::setnumberOfSlices(unsigned NumberOfSlices) {
    numberOfSlices = NumberOfSlices;
}

unsigned operaSpectrograph::getnumberOfSlices(void) {
    return numberOfSlices;
}

void operaSpectrograph::setspectralResolution(double SpectralResolution) {
    spectralResolution = SpectralResolution;
}

double operaSpectrograph::getspectralResolution(void) {
    return spectralResolution;
}

void operaSpectrograph::setSpectrographCCD(spectrographCCD_t spectrographCCD) {
    SpectrographCCD = spectrographCCD;
}

spectrographCCD_t operaSpectrograph::getSpectrographCCD(void) {
    return SpectrographCCD;
}

void operaSpectrograph::setEspadonsCCDReadoutSpeed(EspadonsCCDReadoutSpeed_t espadonsCCDReadoutSpeed) {
    EspadonsCCDReadoutSpeed = espadonsCCDReadoutSpeed;
}

EspadonsCCDReadoutSpeed_t operaSpectrograph::getEspadonsCCDReadoutSpeed(void) {
    return EspadonsCCDReadoutSpeed;
}

void operaSpectrograph::setEspadonsInstrumentMode(EspadonsInstrumentMode_t espadonsInstrumentMode) {
    EspadonsInstrumentMode = espadonsInstrumentMode;
    setOpticalFiberFromInstrumentMode();
}

EspadonsInstrumentMode_t operaSpectrograph::getEspadonsInstrumentMode(void) {
    return EspadonsInstrumentMode;
}

void operaSpectrograph::setx_pixelsize(double x) {
    x_pixelsize = x;
}

double operaSpectrograph::getx_pixelsize(void) {
    return x_pixelsize;
}

void operaSpectrograph::sety_pixelsize(double y) {
    y_pixelsize = y;
}

double operaSpectrograph::gety_pixelsize(void) {
    return y_pixelsize;
}


double operaSpectrograph::getSaturationExptime(double Vmag) {
    
    double texp[9];
    double mag[9];
    double tsat;
    int i;
    double xout[1], yout[1];
    
    for(i=0;i<9;i++) {
        mag[i] = (double)i;
    }

    if(SpectrographCCD == Olapa) {
        texp[0] = 2;
        texp[1] = 5;
        texp[2] = 12;
        texp[3] = 31;
        texp[4] = 78;
        texp[5] = 195;
        texp[6] = 488;
        texp[7] = 1220;
        texp[8] = 3050;
    } else if (SpectrographCCD == EEV1) {
        texp[0] = 4;
        texp[1] = 10;
        texp[2] = 25;
        texp[3] = 63;
        texp[4] = 156;
        texp[5] = 390;
        texp[6] = 976;
        texp[7] = 2440;
        texp[8] = 6100;
    }
    
    xout[0] = Vmag;
    yout[0] = 0;
    
    operaFitSplineDouble(9, mag, texp, 1, xout, yout);
    
    tsat = yout[0];
    
    return tsat;
}

double operaSpectrograph::CalculateIPIE(double seeing) {
    
    double sig = seeing/(2.*sqrt(2.*log(2.)));
    
    double ipie = 1. - exp(- ((InjectionHoleDiameter/2.0)*(InjectionHoleDiameter/2.0))/(2.*sig*sig));
    
    return ipie;
}

double operaSpectrograph::OpticsThroughput(double wavelength_nm) {
    double throughput = 1;
    return throughput;
}

double operaSpectrograph::FiberThroughput(double wavelength_nm) {
    double throughput = 1.0;  // initialize with clear throughput
    
    switch (OpticalFiber) {
        case STUPolymicro:   // Transmission data available fro 100m only - should be good enough for Espadons (l=33m)
            throughput = STUDPolymicroFiber100m_throughput(wavelength_nm);
            break;
        case FBPPolymicro:  // calculate throughput using fiberLength - for GRACES
            throughput = FBPPolymicroFiber_throughput(wavelength_nm);
            break;
        default:
            break;
    }
    
    return throughput;
}

double operaSpectrograph::STUDPolymicroFiber100m_throughput(double wavelength_nm) {
    double wl[56],throughput[56];
    
    wl[0]=304.478; throughput[0]=0.375940;
    wl[1]=349.254; throughput[1]=0.375940;
    wl[2]=367.164; throughput[2]=2.63158;
    wl[3]=380.597; throughput[3]=10.5263;
    wl[4]=394.030; throughput[4]=18.0451;
    wl[5]=411.940; throughput[5]=28.9474;
    wl[6]=447.761; throughput[6]=48.4962;
    wl[7]=465.672; throughput[7]=54.8872;
    wl[8]=505.970; throughput[8]=66.1654;
    wl[9]=550.746; throughput[9]=75.1880;
    wl[10]=573.134; throughput[10]=76.3158;
    wl[11]=595.522; throughput[11]=72.9323;
    wl[12]=608.955; throughput[12]=65.7895;
    wl[13]=631.343; throughput[13]=60.5263;
    wl[14]=649.254; throughput[14]=57.8947;
    wl[15]=662.687; throughput[15]=63.9098;
    wl[16]=685.075; throughput[16]=74.8120;
    wl[17]=694.030; throughput[17]=83.8346;
    wl[18]=711.940; throughput[18]=88.7218;
    wl[19]=738.806; throughput[19]=90.9774;
    wl[20]=770.149; throughput[20]=92.8571;
    wl[21]=846.269; throughput[21]=94.3609;
    wl[22]=904.478; throughput[22]=95.8647;
    wl[23]=935.821; throughput[23]=95.8647;
    wl[24]=953.731; throughput[24]=94.7368;
    wl[25]=976.119; throughput[25]=95.8647;
    wl[26]=1052.24; throughput[26]=96.9925;
    wl[27]=1101.49; throughput[27]=96.6165;
    wl[28]=1155.22; throughput[28]=96.2406;
    wl[29]=1208.96; throughput[29]=97.3684;
    wl[30]=1231.34; throughput[30]=95.8647;
    wl[31]=1244.78; throughput[31]=94.7368;
    wl[32]=1280.60; throughput[32]=95.4887;
    wl[33]=1307.46; throughput[33]=95.8647;
    wl[34]=1334.33; throughput[34]=93.9850;
    wl[35]=1356.72; throughput[35]=90.9774;
    wl[36]=1365.67; throughput[36]=86.4662;
    wl[37]=1365.67; throughput[37]=80.0752;
    wl[38]=1374.63; throughput[38]=72.9323;
    wl[39]=1388.06; throughput[39]=66.5414;
    wl[40]=1401.49; throughput[40]=75.5639;
    wl[41]=1405.97; throughput[41]=82.3308;
    wl[42]=1414.93; throughput[42]=88.3459;
    wl[43]=1437.31; throughput[43]=93.6090;
    wl[44]=1468.66; throughput[44]=95.8647;
    wl[45]=1517.91; throughput[45]=96.2406;
    wl[46]=1576.12; throughput[46]=95.4887;
    wl[47]=1625.37; throughput[47]=94.7368;
    wl[48]=1679.10; throughput[48]=94.7368;
    wl[49]=1741.79; throughput[49]=92.4812;
    wl[50]=1795.52; throughput[50]=89.4737;
    wl[51]=1840.30; throughput[51]=84.5865;
    wl[52]=1876.12; throughput[52]=79.3233;
    wl[53]=1920.90; throughput[53]=69.9248;
    wl[54]=1952.24; throughput[54]=59.7744;
    wl[55]=1983.58; throughput[55]=51.5038;
    
    unsigned nin = 56;
    
    double yp1 = (throughput[1] - throughput[0])/(wl[1] - wl[0]);
	double ypn = (throughput[nin-1] - throughput[nin-2])/(wl[nin-1] - wl[nin-2]);
	double *y2 = new double[nin];
	
	// Call cubicspline to get second derivatives
	cubicsplineDouble(wl, throughput, nin, yp1, ypn, y2);
	
    double outputthroughput;
    
	// Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, throughput, y2, nin, wavelength_nm , &outputthroughput);
    
    delete[] y2;
    
    return outputthroughput/100;
}

double operaSpectrograph::FBPPolymicroFiber_throughput(double wavelength_nm) {
    return FiberTransmissionFromAttenuation(FBPPolymicroFiber_attenuation_db_per_km(wavelength_nm));
}

double operaSpectrograph::FiberTransmissionFromAttenuation(double attenuation_db_per_km) {
// based on formula: Attenuation (dB) = 10 x log10 ( Input Intensity (W) /Output Intensity (W) )
    double fiberLength_km = fiberLength/1000;
    return pow(10,-attenuation_db_per_km*fiberLength_km/10);
}

double operaSpectrograph::FBPPolymicroFiber_attenuation_db_per_km(double wavelength_nm) {
    double wl[61],attenuation[61];
    
    wl[0]=404.645; attenuation[0]=357.925;
    wl[1]=414.412; attenuation[1]=322.781;
    wl[2]=434.840; attenuation[2]=264.654;
    wl[3]=455.277; attenuation[3]=230.852;
    wl[4]=485.488; attenuation[4]=179.473;
    wl[5]=536.149; attenuation[5]=125.374;
    wl[6]=596.599; attenuation[6]=92.8868;
    wl[7]=657.943; attenuation[7]=75.2639;
    wl[8]=711.289; attenuation[8]=67.1078;
    wl[9]=793.084; attenuation[9]=48.1152;
    wl[10]=809.088; attenuation[10]=46.7494;
    wl[11]=839.317; attenuation[11]=39.9654;
    wl[12]=877.550; attenuation[12]=39.9310;
    wl[13]=926.457; attenuation[13]=49.3464;
    wl[14]=943.365; attenuation[14]=84.4663;
    wl[15]=955.806; attenuation[15]=66.8875;
    wl[16]=977.138; attenuation[16]=46.5980;
    wl[17]=996.697; attenuation[17]=39.8236;
    wl[18]=1032.26; attenuation[18]=37.0889;
    wl[19]=1071.38; attenuation[19]=35.7023;
    wl[20]=1107.84; attenuation[20]=37.0208;
    wl[21]=1141.63; attenuation[21]=39.6930;
    wl[22]=1173.64; attenuation[22]=39.6642;
    wl[23]=1203.87; attenuation[23]=43.6910;
    wl[24]=1227.00; attenuation[24]=67.9945;
    wl[25]=1234.12; attenuation[25]=90.9611;
    wl[26]=1239.47; attenuation[26]=120.686;
    wl[27]=1245.69; attenuation[27]=130.140;
    wl[28]=1252.80; attenuation[28]=116.620;
    wl[29]=1263.46; attenuation[29]=97.6914;
    wl[30]=1284.80; attenuation[30]=76.0506;
    wl[31]=1302.58; attenuation[31]=80.0886;
    wl[32]=1325.71; attenuation[32]=115.203;
    wl[33]=1333.73; attenuation[33]=148.979;
    wl[34]=1345.30; attenuation[34]=192.212;
    wl[35]=1355.99; attenuation[35]=231.392;
    wl[36]=1361.35; attenuation[36]=301.657;
    wl[37]=1366.72; attenuation[37]=393.544;
    wl[38]=1371.22; attenuation[38]=535.432;
    wl[39]=1373.94; attenuation[39]=659.754;
    wl[40]=1376.66; attenuation[40]=800.292;
    wl[41]=1377.58; attenuation[41]=869.211;
    wl[42]=1386.48; attenuation[42]=907.040;
    wl[43]=1390.93; attenuation[43]=913.793;
    wl[44]=1396.26; attenuation[44]=897.572;
    wl[45]=1398.91; attenuation[45]=859.732;
    wl[46]=1402.43; attenuation[46]=767.837;
    wl[47]=1408.62; attenuation[47]=667.831;
    wl[48]=1415.67; attenuation[48]=502.960;
    wl[49]=1424.52; attenuation[49]=396.195;
    wl[50]=1433.38; attenuation[50]=315.106;
    wl[51]=1449.35; attenuation[51]=229.956;
    wl[52]=1466.22; attenuation[52]=179.941;
    wl[53]=1482.22; attenuation[53]=151.548;
    wl[54]=1518.66; attenuation[54]=123.137;
    wl[55]=1554.22; attenuation[55]=101.484;
    wl[56]=1581.78; attenuation[56]=96.0533;
    wl[57]=1618.24; attenuation[57]=97.3718;
    wl[58]=1653.81; attenuation[58]=105.448;
    wl[59]=1684.93; attenuation[59]=116.231;
    wl[60]=1699.16; attenuation[60]=116.218;
    
    unsigned nin = 61;
    
    double yp1 = (attenuation[1] - attenuation[0])/(wl[1] - wl[0]);
    double ypn = (attenuation[nin-1] - attenuation[nin-2])/(wl[nin-1] - wl[nin-2]);
    double *y2 = new double[nin];
    
    // Call cubicspline to get second derivatives
    cubicsplineDouble(wl, attenuation, nin, yp1, ypn, y2);
    
    double outputattenuation;
    
    // Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, attenuation, y2, nin, wavelength_nm , &outputattenuation);
    
    delete[] y2;
    
    return outputattenuation;
}

double operaSpectrograph::CCDQuantumEfficiency(double wavelength_nm) {
    double qe = 1.0;
        
    switch (SpectrographCCD) {
        case Olapa:
            qe = OlapaQuantumEfficiency(wavelength_nm);
            break;
        case EEV1:
            qe = EEV1QuantumEfficiency(wavelength_nm);
            break;
        default:
            break;
    }
    
    return qe;
}
double operaSpectrograph::OlapaQuantumEfficiency(double wavelength_nm) {
    
    double wl[43],Olapa_QE[43];
    
    wl[	0	]=	300.000	; Olapa_QE[	0	]=	0.361757	;
    wl[	1	]=	309.977	; Olapa_QE[	1	]=	0.413437	;
    wl[	2	]=	319.048	; Olapa_QE[	2	]=	0.459948	;
    wl[	3	]=	329.025	; Olapa_QE[	3	]=	0.498708	;
    wl[	4	]=	339.002	; Olapa_QE[	4	]=	0.534884	;
    wl[	5	]=	349.887	; Olapa_QE[	5	]=	0.563307	;
    wl[	6	]=	359.864	; Olapa_QE[	6	]=	0.614987	;
    wl[	7	]=	369.841	; Olapa_QE[	7	]=	0.674419	;
    wl[	8	]=	379.819	; Olapa_QE[	8	]=	0.739018	;
    wl[	9	]=	389.796	; Olapa_QE[	9	]=	0.780362	;
    wl[	10	]=	398.866	; Olapa_QE[	10	]=	0.811370	;
    wl[	11	]=	409.751	; Olapa_QE[	11	]=	0.832041	;
    wl[	12	]=	420.635	; Olapa_QE[	12	]=	0.837209	;
    wl[	13	]=	430.612	; Olapa_QE[	13	]=	0.842377	;
    wl[	14	]=	440.590	; Olapa_QE[	14	]=	0.844961	;
    wl[	15	]=	449.660	; Olapa_QE[	15	]=	0.847545	;
    wl[	16	]=	461.451	; Olapa_QE[	16	]=	0.847545	;
    wl[	17	]=	470.522	; Olapa_QE[	17	]=	0.844961	;
    wl[	18	]=	490.476	; Olapa_QE[	18	]=	0.839793	;
    wl[	19	]=	520.408	; Olapa_QE[	19	]=	0.826873	;
    wl[	20	]=	541.270	; Olapa_QE[	20	]=	0.819121	;
    wl[	21	]=	577.551	; Olapa_QE[	21	]=	0.801034	;
    wl[	22	]=	616.553	; Olapa_QE[	22	]=	0.785530	;
    wl[	23	]=	643.764	; Olapa_QE[	23	]=	0.775194	;
    wl[	24	]=	688.209	; Olapa_QE[	24	]=	0.762274	;
    wl[	25	]=	725.397	; Olapa_QE[	25	]=	0.749354	;
    wl[	26	]=	758.957	; Olapa_QE[	26	]=	0.736434	;
    wl[	27	]=	789.796	; Olapa_QE[	27	]=	0.713178	;
    wl[	28	]=	821.542	; Olapa_QE[	28	]=	0.669251	;
    wl[	29	]=	843.311	; Olapa_QE[	29	]=	0.625323	;
    wl[	30	]=	860.544	; Olapa_QE[	30	]=	0.583979	;
    wl[	31	]=	871.429	; Olapa_QE[	31	]=	0.550388	;
    wl[	32	]=	885.941	; Olapa_QE[	32	]=	0.501292	;
    wl[	33	]=	901.361	; Olapa_QE[	33	]=	0.449612	;
    wl[	34	]=	925.850	; Olapa_QE[	34	]=	0.361757	;
    wl[	35	]=	948.526	; Olapa_QE[	35	]=	0.276486	;
    wl[	36	]=	976.644	; Olapa_QE[	36	]=	0.183463	;
    wl[	37	]=	1006.58	; Olapa_QE[	37	]=	0.103359	;
    wl[	38	]=	1024.72	; Olapa_QE[	38	]=	0.0620155	;
    wl[	39	]=	1041.04	; Olapa_QE[	39	]=	0.0387597	;
    wl[	40	]=	1052.83	; Olapa_QE[	40	]=	0.0232558	;
    wl[	41	]=	1070.07	; Olapa_QE[	41	]=	0.0129199	;
    wl[	42	]=	1094.56	; Olapa_QE[	42	]=	0.00775194	;
    
    unsigned nin = 43;
    
    double yp1 = (Olapa_QE[1] - Olapa_QE[0])/(wl[1] - wl[0]);
    double ypn = (Olapa_QE[nin-1] - Olapa_QE[nin-2])/(wl[nin-1] - wl[nin-2]);
    double *y2 = new double[nin];
    
    // Call cubicspline to get second derivatives
    cubicsplineDouble(wl, Olapa_QE, nin, yp1, ypn, y2);
    
    double outputOlapa_QE;
    
    // Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, Olapa_QE, y2, nin, wavelength_nm , &outputOlapa_QE);
    
    delete[] y2;
    
    return outputOlapa_QE;
}

double operaSpectrograph::EEV1QuantumEfficiency(double wavelength_nm) {
    double wl[38],eev1_QE[38];
    
    wl[0]=321.185; eev1_QE[0]=0.515642;
    wl[1]=330.668; eev1_QE[1]=0.540279;
    wl[2]=341.639; eev1_QE[2]=0.502132;
    wl[3]=352.743; eev1_QE[3]=0.560397;
    wl[4]=362.257; eev1_QE[4]=0.607456;
    wl[5]=371.796; eev1_QE[5]=0.672452;
    wl[6]=381.276; eev1_QE[6]=0.694847;
    wl[7]=390.808; eev1_QE[7]=0.755358;
    wl[8]=401.860; eev1_QE[8]=0.775507;
    wl[9]=431.790; eev1_QE[9]=0.782149;
    wl[10]=452.359; eev1_QE[10]=0.851599;
    wl[11]=471.222; eev1_QE[11]=0.826882;
    wl[12]=502.684; eev1_QE[12]=0.802130;
    wl[13]=532.571; eev1_QE[13]=0.777382;
    wl[14]=551.496; eev1_QE[14]=0.797509;
    wl[15]=571.956; eev1_QE[15]=0.788483;
    wl[16]=600.284; eev1_QE[16]=0.774951;
    wl[17]=631.764; eev1_QE[17]=0.763652;
    wl[18]=650.631; eev1_QE[18]=0.741177;
    wl[19]=671.088; eev1_QE[19]=0.729909;
    wl[20]=700.963; eev1_QE[20]=0.696193;
    wl[21]=730.871; eev1_QE[21]=0.687140;
    wl[22]=733.981; eev1_QE[22]=0.657983;
    wl[23]=751.217; eev1_QE[23]=0.595155;
    wl[24]=771.680; eev1_QE[24]=0.588371;
    wl[25]=801.536; eev1_QE[25]=0.541202;
    wl[26]=831.370; eev1_QE[26]=0.478338;
    wl[27]=851.749; eev1_QE[27]=0.411016;
    wl[28]=870.573; eev1_QE[28]=0.357151;
    wl[29]=900.373; eev1_QE[29]=0.269623;
    wl[30]=930.192; eev1_QE[30]=0.195548;
    wl[31]=950.599; eev1_QE[31]=0.148406;
    wl[32]=971.015; eev1_QE[32]=0.107990;
    wl[33]=999.312; eev1_QE[33]=0.0720358;
    wl[34]=1010.29; eev1_QE[34]=0.0383726;
    wl[35]=1030.73; eev1_QE[35]=0.0181358;
    wl[36]=1060.64; eev1_QE[36]=0.00684121;
    wl[37]=1100.00; eev1_QE[37]=4.41E-06;
    
    unsigned nin = 38;
    
    double yp1 = (eev1_QE[1] - eev1_QE[0])/(wl[1] - wl[0]);
    double ypn = (eev1_QE[nin-1] - eev1_QE[nin-2])/(wl[nin-1] - wl[nin-2]);
    double *y2 = new double[nin];
    
    // Call cubicspline to get second derivatives
    cubicsplineDouble(wl, eev1_QE, nin, yp1, ypn, y2);
    
    double outputeev1_QE;
    
    // Call splineinterpolate for interpolations
    splineinterpolateDouble(wl, eev1_QE, y2, nin, wavelength_nm , &outputeev1_QE);
    
    delete[] y2;
    
    return outputeev1_QE;
}

string operaSpectrograph::getCCDName(void) {
    string outputname;
    switch (SpectrographCCD) {
        case Olapa:
            outputname.assign("Olapa");
            break;
        case EEV1:
            outputname.assign("EEV1");
            break;
        default:
            break;
    }
    return outputname;
}

double operaSpectrograph::getNominalResolution (void) {
    switch (EspadonsInstrumentMode) {
        case polarimetric:
            return 65000;
            break;
        case staronly:
            return 80000;
            break;
        case starplussky:
            return 65000;
            break;
        case GRACES_staronly:
            return 55000;
            break;
        case GRACES_starplussky:
            return 33000;
            break;
        default:
            break;
    }
    return 0;
}

double operaSpectrograph::getNominalAperture (void) {
    switch (EspadonsInstrumentMode) {
        case polarimetric:
            return 20;
            break;
        case staronly:
            return 10;
            break;
        case starplussky:
            return 20;
            break;
        case GRACES_staronly:
            return 10;
            break;
        case GRACES_starplussky:
            return 20;
            break;
        default:
            break;
    }
    return 0;
}

unsigned operaSpectrograph::getNumberOfExposures (void) {
    switch (EspadonsInstrumentMode) {
        case polarimetric:
            return 4;
            break;
        case staronly:
            return 1;
            break;
        case starplussky:
            return 1;
            break;
        case GRACES_staronly:
            return 1;
            break;
        case GRACES_starplussky:
            return 1;
            break;
        default:
            break;
    }
    return 0;
}

string operaSpectrograph::getFullModeName(void) {
    string outputname;
    switch (EspadonsInstrumentMode) {
        case polarimetric:
            outputname.assign("polarimetric (res=65k)");
            break;
        case staronly:
            outputname.assign("spectroscopic / star+sky (res=65k)");
            break;
        case starplussky:
            outputname.assign("spectroscopic / star only (res=80k)");
            break;
        case GRACES_staronly:
            outputname.assign("GRACES spectroscopic / star only (res=55k)");
            break;
        case GRACES_starplussky:
            outputname.assign("GRACES spectroscopic / star + sky (res=38k)");
            break;
        default:
            break;
    }
    return outputname;
}

double operaSpectrograph::getNominalNoise(void) {
    double nominalNoise = 0;
    switch (SpectrographCCD) {
        case Olapa:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalNoise = 2.99;
                    break;
                case normalmode:
                    nominalNoise = 3.8;
                    break;
                case fastmode:
                    nominalNoise = 4.12;
                    break;
                default:
                    break;
            }
            break;
        case EEV1:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalNoise = 3.19;
                    break;
                case normalmode:
                    nominalNoise = 4.36;
                    break;
                case fastmode:
                    nominalNoise = 7.64;
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    return nominalNoise;
}

double operaSpectrograph::getNominalGain(void) {
    double nominalGain = 0;
    switch (SpectrographCCD) {
        case Olapa:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalGain = 1.09;
                    break;
                case normalmode:
                    nominalGain = 1.22;
                    break;
                case fastmode:
                    nominalGain = 1.51;
                    break;
                default:
                    break;
            }
            break;
        case EEV1:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalGain = 1.15;
                    break;
                case normalmode:
                    nominalGain = 1.27;
                    break;
                case fastmode:
                    nominalGain = 1.69;
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    return nominalGain;
}


double operaSpectrograph::getNominalReadoutTime(void) {
    double nominalReadout = 0;
    switch (SpectrographCCD) {
        case Olapa:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalReadout = 60;
                    break;
                case normalmode:
                    nominalReadout = 35.5;
                    break;
                case fastmode:
                    nominalReadout = 31.5;
                    break;
                default:
                    break;
            }
            break;
        case EEV1:
            switch (EspadonsCCDReadoutSpeed) {
                case slowmode:
                    nominalReadout = 67;
                    break;
                case normalmode:
                    nominalReadout = 40;
                    break;
                case fastmode:
                    nominalReadout = 25;
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    return nominalReadout;
}


void operaSpectrograph::setOpticalFiberFromInstrumentMode(void) {
    switch (EspadonsInstrumentMode) {
        case polarimetric:
            OpticalFiber= STUPolymicro;
            fiberLength = 33;
            break;
        case staronly:
            OpticalFiber= STUPolymicro;
            fiberLength = 33;
            break;
        case starplussky:
            OpticalFiber= STUPolymicro;
            fiberLength = 33;
            break;
        case GRACES_staronly:
            OpticalFiber= FBPPolymicro;
            fiberLength = 270;
            break;
        case GRACES_starplussky:
            OpticalFiber= FBPPolymicro;
            fiberLength = 270;
            break;
        default:
            break;
    }
}
