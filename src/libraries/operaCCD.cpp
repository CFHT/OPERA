/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCCD
 Version: 1.0
 Description: CCD-related library routines.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "globaldefines.h"
#include "operaError.h"

/*!
 * \brief CCD library.
 * \file operaCCD.c
 */

#include "libraries/operaLibCommon.h"
#include "libraries/operaCCD.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/ladfit.h"	
#include "libraries/operaFit.h"
#include "libraries/operaStats.h"

/*!
 * operaCCD
 * \author Doug Teeple & Eder Martioli
 * \brief CCD library.
 * \details {This library supports operaGain and operaIMShift}
 * \ingroup libraries
 */

/*
 *  Routines used by operaGain
 */

/* 
 * void operaMaskPixbyCountRange(unsigned npixels, float *array, float *badpixmask, float *newmask, float minvalue, float maxvalue, bool invert)
 * \brief This function creates a new mask based on a previous badpixel mask and on a defined range of values
 * \param npixels is an unsigned for the number of elements in the array
 * \param array is a float pointer with the data
 * \param badpixmask is a float pointer with the input badpixel mask data
 * \param newmask is a float pointer for the output mask
 * \param minvalue is a float that defines the minimum data value to keep mask on (value=1) 
 * \param maxvalue is a float that defines the maximum data value to keep mask on (value=1)
 * \param invert is a bool that allows one to invert the mask values for the selected elements
 * \return void
 */
void operaMaskPixbyCountRange(unsigned npixels, float *array, float *badpixmask, float *newmask, float minvalue, float maxvalue, unsigned char invert) {
	while (npixels--) {
		if(*array > minvalue && *array < maxvalue) {
			if(invert) {
				*newmask = 0.0;
			} else {	
				*newmask = *badpixmask;				
			}	
		} else {
			if(invert) {
				*newmask = *badpixmask;
			} else {	
				*newmask = 0.0;
			}	
		}	
		newmask++;
		array++;
		badpixmask++;
	}			
}


/* 
 * void operaMarkPixbyCountRange(unsigned npixels, float *array, float *previousmark, float *newmark, float minvalue, float maxvalue, unsigned index2mark)
 * \brief This function creates a mask where it marks pixels with the same number based on a defined range of values. 
 * \note It keeps the previous mark mask values.
 * \param npixels is an unsigned for the number of elements in the array
 * \param array is a float pointer with the data
 * \param previousmarks is a float pointer with the input badpixel mask data
 * \param newmarks is a float pointer for the output mask
 * \param minvalue is a float that defines the minimum data value to keep mask on (value=1) 
 * \param maxvalue is a float that defines the maximum data value to keep mask on (value=1)
 * \param index2mark is an unsigned to mark values on the selected elements
 * \return void
 */
void operaMarkPixbyCountRange(unsigned npixels, float *array, float *previousmarks, float *newmarks, float minvalue, float maxvalue, float index2mark) {
	while (npixels--) {
		if(*array > minvalue && *array < maxvalue) {
			*newmarks = index2mark;	
		} else {
			*newmarks = *previousmarks;
		}	
		newmarks++;
		array++;
		previousmarks++;
	}			
}


/* 
 * void operaCCDGainNoise(unsigned npixels, unsigned nbias, float *biasdata[], unsigned nflat, float *flatdata[], float *badpixdata, float lowcount, float highcount, float maxbins, float minnpixelsperbin, float *gain, float *bias, float *noise)
 * \brief This function calculates the CCD gain and noise given a set of bias and a set of flat arrays. 
 * \param npixels is an unsigned for the number of pixels in the arrays,
 * \param nbias is an unsigned for the number of bias arrays,
 * \param biasdata is an array of float pointers that contain the bias data,
 * \param nflat is an unsigned for the number of flat arrays,
 * \param badpixdata is a float pointer that contains the badpixel mask data,
 * \param lowcount is a float for the lowest pixel count value to be considered
 * \param highcount is a float for the highest pixel count value to be considered
 * \param maxbins is a float for the number of bins to divide the range considered
 * \param minnpixelsperbin is a float for the minimum number of pixels per bin allowed in order to consider the bin
 * \param gain is a float pointer that returns the calculated gain value
 * \param bias is a float pointer that returns the calculated bias value 
 * \param noise is a float pointer that returns the calculated noise value 
 * \return void
 */

void operaCCDGainNoise(unsigned npixels, unsigned nbias, float *biasdata[], unsigned nflat, float *flatdata[], float *badpixdata, float lowcount, float highcount, float maxbins, float minnpixelsperbin, float *gain, float *gainerror, float *bias, float *noise) 
{	
	float *masterflat = NULL;
	float *flatdif = NULL; 
	float *biasdif = NULL;	
	float *masktmp = NULL;
	float *indexmask = NULL;
	float bin = (highcount - lowcount)/maxbins;
	double *dxbin = NULL;
	double *dybin = NULL;		
	
	masterflat = (float *) malloc (sizeof(float)*npixels);
	operaImMeanQuick(nflat, npixels, masterflat, flatdata); 					
	
	flatdif = (float *) malloc (sizeof(float)*npixels); 
	biasdif = (float *) malloc (sizeof(float)*npixels);		
	masktmp = (float *) malloc (sizeof(float)*npixels);
	indexmask = (float *) malloc (sizeof(float)*npixels);
	
	memset(indexmask, 0, sizeof(float)*npixels);
	
	unsigned npb = 0;
	
	for (unsigned k = 0; k < maxbins ; k++) { 
		if (operaCountPixels(npixels, masterflat, bin*(float)k, bin*(float)(k+1)) > minnpixelsperbin) {
			memcpy(masktmp, indexmask, sizeof(float)*npixels); 
			operaMarkPixbyCountRange(npixels, masterflat, masktmp, indexmask, bin*(float)k, bin*(float)(k+1), (float)npb+1);										
			npb++;
		}
	}	
	
	float *sigbiasdifarr = (float *) malloc (sizeof(float)*npb*nbias);
	
	unsigned npsigbias = 0;
	
	float *xbin = (float *) malloc (sizeof(float)*npb);
	float *ybin = (float *) malloc (sizeof(float)*npb);	
	
	float *meanbias1_tmp = (float *) malloc (sizeof(float)*npb);
	float *meanbias2_tmp = (float *) malloc (sizeof(float)*npb);
	float *sigbiasdif_tmp = (float *) malloc (sizeof(float)*npb);
	float *meanflat1_tmp = (float *) malloc (sizeof(float)*npb);
	float *meanflat2_tmp = (float *) malloc (sizeof(float)*npb);
	float *sigflatdif_tmp = (float *) malloc (sizeof(float)*npb);
	
	float *meanbias1 = (float *) malloc (sizeof(float)*npb);
	float *meanbias2 = (float *) malloc (sizeof(float)*npb);		
	float *sigbiasdif = (float *) malloc (sizeof(float)*npb);
	float *meanflat1 = (float *) malloc (sizeof(float)*npb);
	float *meanflat2 = (float *) malloc (sizeof(float)*npb);
	float *sigflatdif = (float *) malloc (sizeof(float)*npb);	
	
	memset(meanbias2, 0, sizeof(float)*npb);		
	memset(meanbias1, 0, sizeof(float)*npb);
	memset(sigbiasdif, 0, sizeof(float)*npb);
	memset(meanflat1, 0, sizeof(float)*npb);
	memset(meanflat2, 0, sizeof(float)*npb);
	memset(sigflatdif, 0, sizeof(float)*npb);
	
    float biasinADU = 0;
    
	if(nbias > 1) {
		for (unsigned ibias = 0; ibias < (nbias-1); ibias++) { 
			
			memcpy(biasdif, biasdata[ibias],sizeof(float)*npixels);
			operaImSubtractIm(npixels, biasdif, biasdata[ibias+1]);	// Note: this generates a lot of negative numbers...	
			
			operaArrayIndexedMeanQuick(npixels, biasdata[ibias], indexmask, npb, meanbias1_tmp);
			operaArrayIndexedMeanQuick(npixels, biasdata[ibias+1], indexmask, npb, meanbias2_tmp);
			operaArrayIndexedSigmaQuick(npixels, biasdif, indexmask, npb, sigbiasdif_tmp);	
			
            biasinADU += operaArrayMedian(npixels, biasdata[ibias])/(float)nbias;
            
			for (unsigned k = 0; k < npb ; k++) {
				meanbias1[k] += meanbias1_tmp[k]/(float)(nbias-1);
				meanbias2[k] += meanbias2_tmp[k]/(float)(nbias-1);
				sigbiasdif[k] += sigbiasdif_tmp[k]/(float)(nbias-1);
				sigbiasdifarr[npsigbias++] = sigbiasdif_tmp[k];
			}	
            
            if(ibias == nbias-2) {
                biasinADU += operaArrayMedian(npixels, biasdata[ibias+1])/(float)nbias;
            }
		} 
	} else if (nbias == 1) {
		operaArrayIndexedMeanQuick(npixels, biasdata[0], indexmask, npb, meanbias1);
		operaArrayIndexedSigmaQuick(npixels, biasdata[0], indexmask, npb, sigbiasdif);		
		
        biasinADU = operaArrayMedian(npixels, biasdata[0]);
        
		for (unsigned k = 0; k < npb ; k++) {         
			sigbiasdifarr[npsigbias++] = sigbiasdif[k];
		}	
	}	
	
	for (unsigned iflat = 0; iflat < (nflat-1); iflat++) {
		
		memcpy(flatdif, flatdata[iflat], sizeof(float)*npixels);
		operaImSubtractIm(npixels, flatdif, flatdata[iflat+1]);	// note that we have to index - 1 in the loop	
		
		operaArrayIndexedMeanQuick(npixels, flatdata[iflat], indexmask, npb, meanflat1_tmp);
		operaArrayIndexedMeanQuick(npixels, flatdata[iflat+1], indexmask, npb, meanflat2_tmp);
		operaArrayIndexedSigmaQuick(npixels, flatdif, indexmask, npb, sigflatdif_tmp);
		
		for (unsigned k = 0; k < npb ; k++) {
			meanflat1[k] += meanflat1_tmp[k]/(float)(nflat-1);
			meanflat2[k] += meanflat2_tmp[k]/(float)(nflat-1);
			sigflatdif[k] += sigflatdif_tmp[k]/(float)(nflat-1);	
		}			
	}
	
	if(nbias > 1) {		
		for (unsigned k = 0; k < npb ; k++) {
			xbin[k] = - sigflatdif[k]*sigflatdif[k]/(sqrt(2.0)*sigbiasdif[k]);
			ybin[k] = - (meanflat1[k] + meanflat2[k] - meanbias1[k] - meanbias2[k])/(sqrt(2.0)*sigbiasdif[k]);
		}
	} if (nbias == 1) {
		for (unsigned k = 0; k < npb ; k++) {
			xbin[k] = - sigflatdif[k]*sigflatdif[k]/(sqrt(2.0)*sigbiasdif[k]);
			ybin[k] = - (meanflat1[k] + meanflat2[k] - 2*meanbias1[k])/(sqrt(2.0)*sigbiasdif[k]);
		}		
	}
	
	float sigbiasdifcomb = (npsigbias>0?operaArrayMedianQuick(npsigbias, sigbiasdifarr):0.0); // destructive
	
	*bias = biasinADU;
    
#ifndef LMFIT			
	float chisqr = 0.0;
	float noiseerror = 0.0;
	if (npb > 0) {
		ladfitWithError(xbin, ybin, npb, noise, &noiseerror, gain, gainerror, &chisqr); /* robust linear fit: f(x) =  a + b*x */
	} else {
		*gain = 1.0;
	}
	
#else
	// below in case one wants to use LM fit instead of medfit. it also works! 
	dxbin = (double *) malloc (sizeof(double)*npb);
	dybin = (double *) malloc (sizeof(double)*npb);		
	for(unsigned i=0; i<npb; i++) {
		dxbin[i] = (double)xbin[i];
		dybin[i] = (double)ybin[i];
	}			
	double dchisqr;
	double par[2] = {0,1};
	int npar = 2;
	operaLMFitPolynomial(npb, dxbin, dybin, npar, par, &dchisqr,1);
	*gain = par[1];
#endif
	if(nbias > 1) {			
		*noise = (*gain) * sigbiasdifcomb / sqrt(2); // the fit above cannot provide a good estimate for noise, so we recalculate it
	} else if ( nbias == 1) {
		*noise = (*gain) * sigbiasdifcomb;
	}
	
	free(meanbias1_tmp);
	free(meanbias2_tmp);
	free(sigbiasdif_tmp);	
	free(meanflat1_tmp);
	free(meanflat2_tmp);
	free(sigflatdif_tmp);
	
	free(meanbias1);
	free(meanbias2);	
	free(sigbiasdif);
	free(meanflat1);
	free(meanflat2);
	free(sigflatdif);	
	
	free(masterflat);
	free(flatdif); 
	free(biasdif);		
	free(masktmp);
	free(xbin);
	free(ybin);
	free(sigbiasdifarr);			
	free(dxbin);
	free(dybin);	
}

/*
 *
 * This section is for operaGeometry/operaFindOrders
 *
 */
void operaCCDFitIP(unsigned np, float *x,float *y, unsigned nords, float *xmean, float *ymean, float *ipfunc, float *ipx,float *iperr, unsigned slit)
{
	unsigned i,j,k;
	
	float *my = (float *) malloc (nords * sizeof(float));
	memset(my, 0, sizeof(float)*nords);
	
	float *mx = (float *) malloc (nords * sizeof(float));
	memset(mx, 0, sizeof(float)*nords);
	
	float *bkg, *xbkg;
	bkg = (float *) malloc (2*slit * sizeof(float));
	xbkg = (float *) malloc (2*slit * sizeof(float));
	float a,b,abdev;
	
	float *yin, *xin;
	yin = (float *) malloc (2*slit * sizeof(float));
	xin = (float *) malloc (2*slit * sizeof(float));
	
	unsigned nbkg=0, nin=0;
	
	float ipf[MAXIPPOINTS][MAXORDERS];
	float ipxout[MAXIPPOINTS][MAXORDERS];
	
	j=0;
	
	// loop over given cross-section
	for(i=0;i<np;i++){
		
		if (x[i] >= x[np-1] - 2*(float)slit || j>nords-1) {
			break;
		}
		if((x[i] >= xmean[j] - (float)slit && x[i] <= xmean[j] - (float)slit/2) ||
		   (x[i] >= xmean[j] + (float)slit/2 && x[i] <= xmean[j] + (float)slit)) {
			xbkg[nbkg] = x[i];
			bkg[nbkg] = y[i];
			nbkg++;
		} else if (x[i] > xmean[j] - (float)(slit/2) && x[i] < xmean[j] + (float)(slit/2)) {
			if(nin < (float)slit) {
				xin[nin] = x[i];
				yin[nin] = y[i];
				nin++;
			}
		} else if (x[i] > xmean[j] + (float)slit && x[i] < x[np-1] - 2*(float)slit) {
			
			ladfit(xbkg,bkg,nbkg,&a,&b,&abdev);
			
			float totredflux=0;
			for(k=0;k<nin;k++) {
				totredflux += (yin[k] - (a + b*xin[k]));
			}
			
			for(k=0;k<nin;k++) {
				mx[j] += xin[k]*(yin[k] - (a + b*xin[k]))/totredflux;
				ipf[k][j] = (yin[k] - (a + b*xin[k]))/totredflux;
				//				printf("%u\t%u\t%f\t%f\n",j,k,xin[k],ipf[j][k]);
			}
			for(k=0;k<nin;k++) {
				ipxout[k][j] = xin[k] - mx[j];
			}
			
			nbkg = 0;
			nin = 0;
			j++;
			i -= slit/2;
		}
	}
	
	float totflux = 0;
	for(k=0;k<slit;k++) {
		ipfunc[k] = operaArrayMedian(nords,ipf[k]);
		ipx[k] = operaArrayMean(nords,ipxout[k]);
		iperr[k] = operaArrayMedianSigma(nords,ipf[k],ipfunc[k]);
		if(ipfunc[k] < 0) {
			ipfunc[k] = 0;
		}
		totflux += ipfunc[k];
		//		printf("%u\t%f\t%f\n",k,ipfunc[k],iperr[k]); 
	}
	
	for(k=0;k<slit;k++) {
		ipfunc[k] /= totflux;
		//		printf("%d\t%f\n",(int)k-(int)slit/2,ipfunc[k]); 		
	}
}

/*
 * This function constructs an order map based on the crosscorrelation between a reference oder
 * map and the current one. This has proven to be unreliable for espadons, since the spaces between
 * orders vary very slowly and sometimes the orders matched can be misleading. -- E. Martioli Mar 11, 2015
 */

int operaCCDDetectMissingOrdersUsingRefMap(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,unsigned nrefs,float *xref,float *yref,int *AbsRefOrdNumber,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,float xrange,float xstep) {
    
    double dx = operaXCorrelationWithRefOrders(np,fx,fy,nrefs,xref,yref,nords,xmean,ymean,slit,xrange,xstep);

    unsigned NumberOfOrdersInRow = 0;
    unsigned j0 = 0;

    for (unsigned i=0; i<nrefs; i++) {
        bool refMatched = false;
        
        // Avoid exceeding the image bound.
        if(xref[i] + dx > fx[np-1] - slit) {
            break;
        }
        
        // Try to match input reference map
        for (unsigned j=j0; j<nords; j++) {
            if(fabs(xmean[j] - (xref[i] + dx)) < slit/2) {
                /*
                if(i==0 && j>0) {
                    for (jj==j-1; j>=0; j--) {
                        float predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);

                        if(fabs(xmean[j] -)
                    }
                }
                */
                float xtmp = xmean[j];
                float ytmp = ymean[j];
                float xerr = slit/4;
                int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
                if(isItAboveNoise) {
                    xord[NumberOfOrdersInRow] = xtmp;
                    yord[NumberOfOrdersInRow] = ytmp;
                    xerrord[NumberOfOrdersInRow] = xerr;
                } else {
                    xord[NumberOfOrdersInRow] = xref[i] + dx;
                    yord[NumberOfOrdersInRow] = yref[i];
                    xerrord[NumberOfOrdersInRow] = slit/4;
                }
                AbsOrdNumber[NumberOfOrdersInRow] = AbsRefOrdNumber[i];
                NumberOfOrdersInRow++;
                j0 = j+1;
                refMatched = true;
                break;
            }
        }
        
        // If order couldn't be matched then try to detect again using the ref map position or
        // force point to exist by inputing map-corrected coordinates
        if(refMatched == false) {
            float xtmp = xref[i] + dx;
            float ytmp = yref[i];
            float xerr = slit/4;
            int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
            if(isItAboveNoise) {
                xord[NumberOfOrdersInRow] = xtmp;
                yord[NumberOfOrdersInRow] = ytmp;
                xerrord[NumberOfOrdersInRow] = xerr;
            } else {
                xord[NumberOfOrdersInRow] = xref[i] + dx;
                yord[NumberOfOrdersInRow] = yref[i];
                xerrord[NumberOfOrdersInRow] = slit/4;
            }
            AbsOrdNumber[NumberOfOrdersInRow] = AbsRefOrdNumber[i];
            NumberOfOrdersInRow++;
            refMatched = true;
        }
        
        // If we are out of detected orders then continue to input orders in the map
        if(j0==nords) {
            for (unsigned ii=i+1; ii<nrefs; ii++) {
                float xtmp = xref[ii] + dx;
                float ytmp = yref[ii];
                float xerr = slit/4;
                int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
                if(isItAboveNoise) {
                    xord[NumberOfOrdersInRow] = xtmp;
                    yord[NumberOfOrdersInRow] = ytmp;
                    xerrord[NumberOfOrdersInRow] = xerr;
                } else {
                    xord[NumberOfOrdersInRow] = xref[ii] + dx;
                    yord[NumberOfOrdersInRow] = yref[ii];
                    xerrord[NumberOfOrdersInRow] = slit/4;
                }
                AbsOrdNumber[NumberOfOrdersInRow] = AbsRefOrdNumber[i];
                NumberOfOrdersInRow++;
            }
            break;
        
        } else if (j0 < nords-1 && i == nrefs-1) {
            
            // Here the map is out of orders but there is still detcted orders, so we continue
            // inputting them as long as they agree within an error of slit/2 with the predicted
            // order position using the spacing polynomial.
            
            for (unsigned j=j0; j<nords; j++) {
                float predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);
                if(predX < fx[np-1] - slit) {
                    if(fabs(xmean[j] - predX) < slit/2) {
                        float xtmp = xmean[j];
                        float ytmp = ymean[j];
                        float xerr = slit/4;
                        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
                        if(isItAboveNoise) {
                            xord[NumberOfOrdersInRow] = xtmp;
                            yord[NumberOfOrdersInRow] = ytmp;
                            xerrord[NumberOfOrdersInRow] = xerr;
                        } else {
                            xord[NumberOfOrdersInRow] = predX;
                            yord[NumberOfOrdersInRow] = yord[NumberOfOrdersInRow-1];
                            xerrord[NumberOfOrdersInRow] = slit/4;
                        }
                        NumberOfOrdersInRow++;
                    } else {
                        float xtmp = predX;
                        float ytmp = yord[NumberOfOrdersInRow-1];
                        float xerr = slit/4;
                        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
                        if(isItAboveNoise) {
                            xord[NumberOfOrdersInRow] = xtmp;
                            yord[NumberOfOrdersInRow] = ytmp;
                            xerrord[NumberOfOrdersInRow] = xerr;
                        } else {
                            xord[NumberOfOrdersInRow] = predX;
                            yord[NumberOfOrdersInRow] = yord[NumberOfOrdersInRow-1];
                            xerrord[NumberOfOrdersInRow] = slit/4;
                        }
                        NumberOfOrdersInRow++;
                    }
                } else {
                    break;
                }
            }
            
            float predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);
            while (predX < fx[np-1] - slit) {
                float xtmp = predX;
                float ytmp = yord[NumberOfOrdersInRow-1];
                float xerr = slit/4;
                int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
                if(isItAboveNoise) {
                    xord[NumberOfOrdersInRow] = xtmp;
                    yord[NumberOfOrdersInRow] = ytmp;
                    xerrord[NumberOfOrdersInRow] = xerr;
                } else {
                    xord[NumberOfOrdersInRow] = predX;
                    yord[NumberOfOrdersInRow] = yord[NumberOfOrdersInRow-1];
                    xerrord[NumberOfOrdersInRow] = slit/4;
                }
                NumberOfOrdersInRow++;

                predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);
            }
        }
    }
    
    //cout << dx << " " << xord[0] << " "  << NumberOfOrdersInRow << endl;
    
    return NumberOfOrdersInRow;
}

/*
 * This is now the main function to construct an order map, where the order numbers must match 
 * the numbers of orders in a reference map. Note that the reference map must be displaced by 
 * at most slit/2 +/- xerror from the current map, otherwise it won't work. This function should
 * be used interactively where the current map should be used as the reference map to detect 
 * orders in an adjacent row or sample. This is the most reliable approach, since any other 
 * attempt to match order numbers eventually fails. -- E. Martioli Mar 11, 2015
 */
int operaCCDDetectMissingOrdersUsingNearMap(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,unsigned nrefs,float *xref,float *yref,int *AbsRefOrdNumber,float *xord,float *yord, float *xerrord, int *AbsOrdNumber) {
    
    // Find the FIRST order to match a near map and save the reference order
    float mindist = BIG;
    bool neverin = true;
    unsigned jref = 0;
    unsigned iref = 0;
    
    for (unsigned i=0; i<nrefs; i++) {
        for (unsigned j=0; j<nords; j++) {
            float dist = fabs(xmean[j] - (xref[i]));
            if(dist < slit/2) {
                if(dist < mindist) {
                    mindist = dist;
                    jref = j;
                    iref = i;
                }
                neverin = false;
            } else if (dist > slit/2 && neverin == false) {
                break;
            }
        }
    }
    
    int refOrderNumber = AbsRefOrdNumber[iref];

    //cout << jref << " " << iref << " " << refOrderNumber << " " << xmean[jref] << " " << xref[iref] << endl;

    unsigned NumberOfOrdersInRow = 0;
    
    /*
     * First figure out orders lying before reference order:
     */
    float *xord_b = new float[MAXORDERS];
    float *yord_b = new float[MAXORDERS];
    float *xerr_b = new float[MAXORDERS];
    float *AbsOrdNumber_b = new float[MAXORDERS];
    
    float predX = xmean[jref] - (float)PolynomialFunction((double)refOrderNumber,par,npar);
    float yRef = ymean[jref];

    unsigned nOrdsBefore = 0;
    
    while(predX > fx[0] + slit) {
        float xtmp = predX;
        float ytmp = yRef;
        float xerr = slit/4;
        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
        if(isItAboveNoise) {
            xord_b[nOrdsBefore] = xtmp;
            yord_b[nOrdsBefore] = ytmp;
            xerr_b[nOrdsBefore] = xerr;
        } else {
            xord_b[nOrdsBefore] = predX;
            yord_b[nOrdsBefore] = yRef;
            xerr_b[nOrdsBefore] = slit/4;
        }
        AbsOrdNumber_b[nOrdsBefore] = refOrderNumber - 1 - nOrdsBefore;
        predX = xord_b[nOrdsBefore] - (float)PolynomialFunction((double)AbsOrdNumber_b[nOrdsBefore],par,npar);
        nOrdsBefore++;
    }
    
    for (unsigned i=0; i<nOrdsBefore; i++) {
        unsigned backwrd_i = nOrdsBefore-i-1;
        
        xord[NumberOfOrdersInRow] = xord_b[backwrd_i];
        yord[NumberOfOrdersInRow] = yord_b[backwrd_i];
        xerrord[NumberOfOrdersInRow] = xerr_b[backwrd_i];
        AbsOrdNumber[NumberOfOrdersInRow] = AbsOrdNumber_b[backwrd_i];
        
        //cout << refOrderNumber << " " << AbsOrdNumber_b[backwrd_i] << " " << xord_b[backwrd_i] << " " << yord_b[backwrd_i] << endl;
        NumberOfOrdersInRow++;
    }
    
    /*
     * Then assign reference order:
     */
    float xtmp = xmean[jref];
    float ytmp = ymean[jref];
    float xerr = xmeanerr[jref];
    int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
    if(isItAboveNoise) {
        xord[NumberOfOrdersInRow] = xtmp;
        yord[NumberOfOrdersInRow] = ytmp;
        xerrord[NumberOfOrdersInRow] = xerr;
    } else {
        xord[NumberOfOrdersInRow] = xmean[jref];
        yord[NumberOfOrdersInRow] = ymean[jref];
        xerrord[NumberOfOrdersInRow] = xmeanerr[jref];
    }
    AbsOrdNumber[NumberOfOrdersInRow] = refOrderNumber;
    NumberOfOrdersInRow++;

    /*
     * Finally figure out orders lying after reference order:
     */
    predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);
    float predY = yord[NumberOfOrdersInRow-1];
    float predXerr = xerrord[NumberOfOrdersInRow-1];

    unsigned jcurr = jref+1;
    if(fabs(xmean[jcurr] - predX) < slit/2) {
        predX = xmean[jcurr];
        predY = ymean[jcurr];
        predXerr = xmeanerr[jcurr];
    }
    
    while(predX < fx[np-1] - slit) {
        float xtmp = predX;
        float ytmp = predY;
        float xerr = predXerr;
        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);
        if(isItAboveNoise) {
            xord[NumberOfOrdersInRow] = xtmp;
            yord[NumberOfOrdersInRow] = ytmp;
            xerrord[NumberOfOrdersInRow] = xerr;
        } else {
            xord[NumberOfOrdersInRow] = predX;
            yord[NumberOfOrdersInRow] = predY;
            xerrord[NumberOfOrdersInRow] = predXerr;
        }
        AbsOrdNumber[NumberOfOrdersInRow] = AbsOrdNumber[NumberOfOrdersInRow-1] + 1;

        // -- increment general order count
        NumberOfOrdersInRow++;
        // --------------

        predX = xord[NumberOfOrdersInRow-1] + (float)PolynomialFunction((double)AbsOrdNumber[NumberOfOrdersInRow-1],par,npar);
        predY = yord[NumberOfOrdersInRow-1];
        predXerr = xerrord[NumberOfOrdersInRow-1];
        
        // -- use existing detected orders as long as they last
        jcurr++;
        if(jcurr < nords) {
            for(unsigned j=jcurr;j<nords;j++) {
                if(fabs(xmean[j] - predX) < slit/2) {
                    predX = xmean[j];
                    predY = ymean[j];
                    predXerr = xmeanerr[j];
                    jcurr = j;
                    break;
                }
            }
        }
    }
    
    return NumberOfOrdersInRow;
}

/*
 *  This constructs an order map based on an input set of detected orders and 
 *  on the order spacing polynomial. -- E. Martioli Mar 11, 2015
 */
int operaCCDDetectOrderMapBasedOnSpacingPolynomial(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,unsigned minordertouse, unsigned maxorders)
{
    float *orderMap = new float[MAXORDERS];
    float *ordSepPred = new float[MAXORDERS];
    float *ordXPosPred = new float[MAXORDERS];
    float *ordYValue = new float[MAXORDERS];
    unsigned *ordIndex = new unsigned[MAXORDERS];

    double leastSqr = BIG;
    unsigned jmin = 0;
    unsigned nAccHops = 5;
    float xiSqr = 0;
    
    for (unsigned j=0;j<nords;j++) {
        
        xiSqr = 0;
        
        matchMeasuredOrdersWithMap(np,fx,fy,slit,npar,par,nords,xmean,ymean,orderMap,ordSepPred,ordXPosPred,ordYValue,ordIndex,j,nAccHops,minordertouse,maxorders,&xiSqr);
        
        // cout << j << " " << xiSqr << endl;
        
        if(xiSqr < leastSqr){
            leastSqr = xiSqr;
            jmin = j;
        }
    }
    
    unsigned NumberOfOrdersInRow = matchMeasuredOrdersWithMap(np,fx,fy,slit,npar,par,nords,xmean,ymean,orderMap,ordSepPred,ordXPosPred,ordYValue,ordIndex,jmin,nAccHops,minordertouse,maxorders,&xiSqr);
    
    for (unsigned o=0;o<NumberOfOrdersInRow;o++) {
        
        float xtmp = ordXPosPred[o];
        float ytmp = ordYValue[o];
        float xerr = slit/4;

        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xtmp,&ytmp,&xerr);

        if(isItAboveNoise) {
            xord[o] = xtmp;
            yord[o] = ytmp;
            xerrord[o] = xerr;
        } else {
            xord[o] = ordXPosPred[o];
            yord[o] = ordYValue[o];
            xerrord[o] = slit/4;
        }

        AbsOrdNumber[o] = (int)orderMap[o];
        // cout << orderMap[o] << " " << ordSepPred[o] << " " << ordXPosPred[o] << " " << ordIndex[o] << endl;
    }
    
    //cout << orderMap[0] << " " << ordXPosPred[0] << " " << ordIndex[0] << endl;
    return NumberOfOrdersInRow;
}

/*
 * This function will match detected orders with a reference map, which is obtained 
 * from a first pass using the spacing polynomial and a reference order, where it has better snr.
 */
unsigned matchMeasuredOrdersWithMap(unsigned np,float *fx,float *fy,float slit, unsigned npar,double *par,unsigned nords,float *xmean,float *ymean, float *orderMap, float *ordSepPred, float *ordXPosPred,float *ordYValue, unsigned *ordIndex, unsigned j, unsigned nAccHops,unsigned minordertouse, unsigned maxorders, float *xiSqr) {
    
    float sumOfsqrs = 0;
    unsigned Nsqrs = 0;
    unsigned NumberOfOrdersInRow = 0;
    unsigned nhops = 0;
    
    for (unsigned o=0;o<=maxorders;o++) {
        
        orderMap[o] = minordertouse + o + j;
        ordSepPred[o] = (float)PolynomialFunction((double)orderMap[o],par,npar);
        
        if (o==0) {
            ordXPosPred[o] = xmean[0];
            ordYValue[o] = ymean[0];
            ordIndex[o] = o;
            NumberOfOrdersInRow++;
        } else if (o > nhops && o-nhops < nords) {
            
            unsigned jj = (unsigned)(o-nhops);
            // The condition below tests if next point is a valid order
            if (fabs(xmean[jj] - (ordXPosPred[o-1] + ordSepPred[o-1])) < slit/2.0) {
                // If order is valid then save order position
                ordXPosPred[o] = xmean[jj];
                ordYValue[o] = ymean[jj];
                sumOfsqrs += (xmean[jj] - (ordXPosPred[o-1] + ordSepPred[o-1]))*(xmean[jj] - (ordXPosPred[o-1] + ordSepPred[o-1]));
                ordIndex[o] = jj;
                NumberOfOrdersInRow++;
                Nsqrs++;
            } else {
                // Here we don't trust the measurements so we use the model and consider
                // that there was a hop
                ordXPosPred[o] = ordXPosPred[o-1] + ordSepPred[o-1];
                ordYValue[o] = ymean[jj];
                ordIndex[o] = MAXORDERS;
                nhops++; // This will make it to reuse the same point
                
                if(ordXPosPred[o] > fx[np-1] - slit) {
                    break; // break loop and don't increment NumberOfOrdersInRow
                } else {
                    NumberOfOrdersInRow++;
                    
                    // Test whether the current measured xmean matches any further predicted order,
                    // which is to make sure the hop was an actual hop
                    bool match = false;
                    float predPos = ordXPosPred[o];
                    
                    for (unsigned hops=0; hops<nAccHops; hops++) {
                        predPos += (float)PolynomialFunction((double)(orderMap[o]+hops),par,npar);
                        if(predPos > fx[np-1] - slit) {
                            break;
                        }
                        if (fabs(xmean[o-nhops] - predPos) < slit/2) {
                            match = true;
                            break;
                        }
                    }
                    
                    // If current measurent doesn't match any model position
                    // then it was not a hop, so we throw away the current point and try the next one
                    if (match==false) {
                        if((int)(o+1)-(nhops-2)>=0){
                            nhops -= 2;
                        } else {
                            nhops--;
                        }
                        
                        if (o-nhops >= nords) {
                            // There is no more points left, so accepts model position
                            // and leave loop.
                            break;
                        } else {
                            // Otherwise go back to previous point (o--) and try again,
                            // since it's possible to find another
                            NumberOfOrdersInRow--;
                            o--;
                        }
                    }
                }
            }
        } else {
            break;
        }
    }

    *xiSqr = sumOfsqrs/(float)Nsqrs;
    return NumberOfOrdersInRow;
}

/*
 * The function below calculates the shift between two order maps in a given data for the
 * maximum cross-correlation between these two maps. -- E. Martioli Mar 11, 2015
 */
float operaXCorrelationWithRefOrders(unsigned np,float *fx,float *fy,unsigned nrefs, float *xref,float *yref,unsigned nords, float *xmean,float *ymean, float slit, float xrange, float xstep) {

    // First create a set of Gaussian parameters for reference orders and
    // input measured orders
    int n_par = 3;
    
    int n_totalpar_ref = (int)nrefs*n_par;
    double *par_ref = new double[n_totalpar_ref];
    for (unsigned i=0; i<nrefs; i++) {
// the line below is commented because without a good flat-field correction
// the order will have different flux levels at different detector positions
// so, in order to keep only the information on the spacing between orders the
// gaussian amplitudes were set to 1.0. But one could revert this by using the
// the actual order fluxes> yref[i] and ymean[i]
//        par_ref[0+i*n_par] = (double)yref[i];
        par_ref[0+i*n_par] = 1.0;
        par_ref[1+i*n_par] = (double)xref[i];
        par_ref[2+i*n_par] = (double)slit/8.0;
    }
    
    int n_totalpar_ords = (int)nords*n_par;
    double *par_ords = new double[n_totalpar_ords];
    for (unsigned i=0; i<nords; i++) {
//        par_ords[0+i*n_par] = (double)ymean[i];
        par_ords[0+i*n_par] = 1.0;
        par_ords[1+i*n_par] = (double)xmean[i];
        par_ords[2+i*n_par] = (double)slit/8.0;
    }
    
    // Then create a synthetic data set for each order map
    double *reffy = new double[np];
    double *ordsfy = new double[np];

    for (unsigned i=0; i<np; i++) {
        reffy[i] = GaussianFunction((double)fx[i],par_ref,n_totalpar_ref);
        ordsfy[i] = GaussianFunction((double)fx[i],par_ords,n_totalpar_ords);
        //cout << fx[i] << " " << fy[i] << " " << reffy[i] << " " << ordsfy[i] << endl;
    }
    double maxcorrelation = -BIG;
    float maxdx = 0;
    
    for(float dx = -xrange/2.0; dx < +xrange/2.0; dx+=xstep) {
        
        for (unsigned i=0; i<np; i++) {
            reffy[i] = GaussianFunction((double)(fx[i] + dx),par_ref,n_totalpar_ref);
        }
    
        double corr = operaCrossCorrelation(np,reffy,ordsfy);
        
        if (corr > maxcorrelation) {
            maxcorrelation = corr;
            maxdx = dx;
        }
        // cout << dx << " " << corr << endl;
    }
    
    /*
    for (unsigned i=0; i<np; i++) {
        reffy[i] = GaussianFunction((double)(fx[i] + maxdx),par_ref,n_totalpar_ref);
        ordsfy[i] = GaussianFunction((double)fx[i],par_ords,n_totalpar_ords);
        cout << fx[i] << " " << fy[i] << " " << reffy[i] << " " << ordsfy[i] << endl;
    }
    */
    return maxdx;
}

/*
 * This function is the first reliable version of opera geometry algorithms to construct
 * a numbered order map. This is still a realiable tool, altough there is an issue that
 * turned into a problem for Espadons. It relies on the model for the spacing between orders
 * versus detector position. For Espadons this function cannot be modeled by a simple poynomial
 * and on previous versions of opera we approximated it to a linear model. This works well
 * for **most** of the case, but in the very few cases where it doesn't work it was assigning 
 * misleading order numbers and consequently causing errors in the geometry calibration. 
 * -- E. Martioli Mar 11, 2015
 */
int operaCCDDetectMissingOrders(unsigned np,float *fx,float *fy,unsigned npip,float *ipiny,float *ipinx,float slit,float noise,float gain, unsigned npar,double *par,unsigned nords, float *xmean,float *ymean,float *xmeanerr,float *xord,float *yord, float *xerrord, int *AbsOrdNumber,int AbsPrevOrdNumber,float x0prev)
{
    /*
     * The portion below is a fix to convert the dependency of spacing polynomial from detector x-position 
     * to ordernumber. the definition of this polynomial used to be spacing versus detector postion and now 
     * order spacing versus order number.
     */
    float *OrderPosition = new float[MAXORDERS];
    float *OrderSeparation = new float[MAXORDERS];
    unsigned npts=0;
    for(unsigned i=0; i<nords-1; i++) {
        OrderPosition[npts] = float(xmean[i] + (xmean[i+1] - xmean[i]) / 2.0);
        OrderSeparation[npts] = float(xmean[i+1] - xmean[i]);
        npts++;
    }
    float am,bm,abdevm;
    //--- Robust linear fit
    ladfit(OrderPosition,OrderSeparation,npts,&am,&bm,&abdevm); // robust linear fit: f(x) =  a + b*x
    double *par_tmp = new double[npar];
    unsigned npars_tmp = npar;
    for (unsigned i=0;i<npars_tmp; i++) {
        par_tmp[i] = par[i];
    }
    par[0] = (double)am;
    par[1] = (double)bm;
    npar = 2;
    //--- end of fix to spacing polynomial
    
	unsigned i;
	float xloword[MAXORDERS],yloword[MAXORDERS],xerrloword[MAXORDERS];
	float xhiord[MAXORDERS],yhiord[MAXORDERS],xerrhiord[MAXORDERS];
	int ordNumberinRow=0;
	unsigned nloword, nhiord;		
	float xmmin, xmmintmp, xerrmin;
	float ymmin;	
	int acceptorder = 0;			
	unsigned initord;
	float xmiss,xmisstmp, xerrmiss;
	float ymiss;	
	int firstorderinrow=0,nd;
	
	nloword=0;
	nhiord=0;
	
	xmmin = xmean[0] - (float)PolynomialFunction((double)xmean[0],par,npar);
	xmmintmp = xmmin;
	
    // Measure orders lying before the first order detected and after fx[0]+slit 
	nloword = 0;
	while(xmmin > fx[0] + slit && nloword < nords && nloword < MAXORDERS) {	// D.T. May 6 2014 check for out of bounds
        int isItAboveNoise = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xmmintmp,&ymmin,&xerrmin);
        
		if(isItAboveNoise) {
			xmmin = xmmintmp;
			xloword[nloword] = xmmintmp;
			yloword[nloword] = ymmin;
			xerrloword[nloword] = xerrmin;
			nloword++;
		} else {
			xloword[nloword] = xmmin;
			yloword[nloword] = noise/gain;
			xerrloword[nloword] = slit/4;
			nloword++;
		}
		//recalculate predicted order position				
		xmmin -= (float)PolynomialFunction((double)xmmin,par,npar);		
		xmmintmp = xmmin;				
	}		
	
	initord=0;
	xmiss = xmean[0];
	xmisstmp = xmiss;
	acceptorder = 0;
	nhiord = 0;
	
	while(xmiss < fx[np-1] - slit && nhiord < nords && nhiord < MAXORDERS) {	// D.T. May 6 2014 check for out of bounds	
		/*
		 * test if there is already a detection for the xmiss predicted order position:
		 */
		// first consider only predicted orders which falls before the last detected order
		if(xmiss < (xmean[nords-1] + (float)PolynomialFunction((double)xmean[nords-1],par,npar) - slit/2)) {
			// loop over detected orders
			for(i=initord;i<nords;i++) {
				// if predicted order falls within the same range as a detected order then pick detected one						
				if(xmiss > (xmean[i] - slit/2) &&	xmiss < (xmean[i] + slit/2)) {
					xhiord[nhiord] = xmean[i];
					yhiord[nhiord] = ymean[i];
                    xerrhiord[nhiord] = xmeanerr[i];                    
					xmiss = xmean[i];
					ymiss = ymean[i];
					nhiord++;							
					initord = i+1; //skip this order next time
					break; //get out of the loop
					// if predicted order falls before next detected order then accept it as a missing order							
				} else if (xmiss < xmean[i] - (float)PolynomialFunction((double)xmean[i],par,npar) + slit/2) {
					// recalculate order position and mean flux. 
					// if level is higher than noise then just take predicted value.
					acceptorder = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xmisstmp,&ymiss,&xerrmiss);
					if(acceptorder && fabs(xmisstmp - xmiss) < slit/2) {
						xmiss = xmisstmp;
						xhiord[nhiord] = xmisstmp;
						yhiord[nhiord] = ymiss;
						xerrhiord[nhiord] = xerrmiss;
						nhiord++;
					} else {
						xhiord[nhiord] = xmiss;
						yhiord[nhiord] = noise/gain;
						xerrhiord[nhiord] = slit/4;
						nhiord++;
					}
					initord = i; // try order i again on next pass of "while" loop
					break;
				} 
			}
		}	else {
			acceptorder = operaCCDRecenterOrderUsingXCorrWithIP(np,fx,fy,npip,ipiny,ipinx,noise,gain,&xmisstmp,&ymiss,&xerrmiss);
			if(acceptorder && fabs(xmisstmp - xmiss) < slit/2) {	
				xmiss = xmisstmp;
				xhiord[nhiord] = xmisstmp;
				yhiord[nhiord] = ymiss;
				xerrhiord[nhiord] = xerrmiss;
				nhiord++;
			} else {
				xhiord[nhiord] = xmiss;
				yhiord[nhiord] = noise/gain;
				xerrhiord[nhiord] = slit/4;
				nhiord++;
			}
			
		}
		//recalculate predicted order position
		xmiss += (float)PolynomialFunction((double)xmiss,par,npar); // set new predicted position		
		xmisstmp = xmiss;				
	}	
	
	ordNumberinRow = 0;
	for(i=0;i<nloword;i++) {
		xord[ordNumberinRow] = xloword[nloword-i-1];
		yord[ordNumberinRow] = yloword[nloword-i-1];
		xerrord[ordNumberinRow] = xerrloword[nloword-i-1];		
		ordNumberinRow++;
	}
	for(i=0;i<nhiord;i++) {
		xord[ordNumberinRow] = xhiord[i];
		yord[ordNumberinRow] = yhiord[i];
		xerrord[ordNumberinRow] = xerrhiord[i];		
		ordNumberinRow++;			
	}	
	// figure out absolute order number for first order in row
	if(x0prev == 0) {
		firstorderinrow = AbsPrevOrdNumber;	
	}	else {
		nd = round((xord[0] - x0prev)/(float)PolynomialFunction((double)((x0prev+xord[0])/2),par,npar));
		firstorderinrow = AbsPrevOrdNumber + nd;						
	}
	
	for(i=0;i<(unsigned)ordNumberinRow;i++) {				
		AbsOrdNumber[i] = firstorderinrow + (int)i;				
		//		cout << AbsOrdNumber[i] << " " << xord[i] << " " << y << " " << yord[i] << endl;		
	}
	
    npar = npars_tmp;
    for (unsigned i=0;i<npar; i++) {
        par[i] = par_tmp[i];
    }
	
	return ordNumberinRow;
}

int operaCCDRecenterOrderUsingXCorrWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float *ipx, float noise, float gain, float *xmean, float *ymean, float *xmeanerr)
{
    float scanrange = (fabs(ipx[0])+fabs(ipx[nip-1]))/FRACTIONOFIPTOSCANFORRECENTERORDER;    
    
    float xmin = (float)floor(*xmean - (0.5*(fabs(ipx[0])+fabs(ipx[nip-1])) + scanrange/2)); // ipx[0] should be a negative number
    float xmax = (float)ceil(*xmean + (0.5*(fabs(ipx[0])+fabs(ipx[nip-1])) + scanrange/2));
	
    unsigned npin = 0;
    double *yin = (double *) malloc (np * sizeof(double));
    double *xin = (double *) malloc (np * sizeof(double));
    unsigned *index = (unsigned *) malloc (np * sizeof(unsigned));
    
    for(unsigned i=0;i<np;i++) {
        if(x[i] >= xmin && x[i] <= xmax) {
            yin[npin] = (double)y[i];
            index[npin] = i;
            npin++;
			//            printf("%d\t%f\t%f\n",i,x[i],y[i]);
        }
	}
	
    if(npin==0) {
        return 0;
    }
    
    double *y_withipsamp = (double *) malloc (nip * sizeof(double));
    double *ipx_double = (double *) malloc (nip * sizeof(double));
    double *ipfunc_double = (double *) malloc (nip * sizeof(double));
    
    for(unsigned i=0;i<nip;i++) {
        ipx_double[i] = (double)ipx[i];
        ipfunc_double[i] = (double)ipfunc[i];
    }
	
    float trial_x_center = *xmean - scanrange/2;
    float trial_x_step = scanrange/(NUMBEROFSAMPLESTORECENTERORDER);
    
    double maxcorrelation = -BIG;
    float newxcenter = 0;
    
    for(unsigned k=0; k<NUMBEROFSAMPLESTORECENTERORDER; k++) {
        for(unsigned i=0;i<npin;i++) {
            xin[i] = (double)(x[index[i]] - trial_x_center);
        }
        operaFitSplineDouble(npin,xin,yin,nip,ipx_double,y_withipsamp);
        
        double crosscorrelation = operaCrossCorrelation(nip,ipfunc_double,y_withipsamp);
        
        if (crosscorrelation > maxcorrelation) {
            maxcorrelation = crosscorrelation;
            newxcenter = trial_x_center;
        }
        
		//        printf("%f\t%f\t%lf\n",trial_x_center,trial_x_step,crosscorrelation);
        
        trial_x_center += trial_x_step;
    }
    
    for(unsigned i=0;i<npin;i++) {
        xin[i] = (double)(x[index[i]] - newxcenter);
    }
    operaFitSplineDouble(npin,xin,yin,nip,ipx_double,y_withipsamp);
    
	float intflux = 0;    
    for(unsigned i=0;i<nip;i++) {
        intflux += ipfunc[i]*(float)y_withipsamp[i];
    }
    
    *xmean = newxcenter;
    *ymean = intflux;
    *xmeanerr = trial_x_step;
    
    free(yin);
    free(xin);
    
    free(y_withipsamp);
    free(ipx_double);
    free(ipfunc_double);
	
	int detectprob = 0;
	
	if(maxcorrelation > MINIMUMCROSSCORRELATIONTORECENTERORDER) {
		detectprob = 1;
	} else {
		detectprob = 0;
	}
	return detectprob; 
}

int operaCCDRecenterOrderWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float *ipx, float noise, float gain, float *xmean, float *ymean, float *xmeanerr)
{  
	unsigned i,j;
	unsigned ix=0;
	unsigned slit=0;
	
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slit++;
		}	
	}
	float *slitfunc = (float *) malloc (slit * sizeof(float));
    
	slit=0;
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slitfunc[slit] = ipfunc[i];
            //printf("%u\t%f\n",slit,slitfunc[slit]);
			slit++;
		}
	}
	
	if(*xmean < x[slit/2] || *xmean > x[np - slit/2]){
		return 0;
	}
	
	float dxmin = fabs(x[np-1] - x[0]);
	
	for(i=0;i<np;i++) {		
		if(fabs(x[i]	- *xmean) < dxmin) {
			dxmin = fabs(x[i]	- *xmean);
			ix = i;
		}
        //printf("%u\t%f\t%f\n",i,x[i],y[i]);
	}
	
	//printf("%u\t%f\t%f\t%f\n",ix,dxmin,x[ix],*xmean);
	
	float my = 0;
	float mx = 0;
	float mxerr = 0;
	
	float intflux = 0;
	
    float *yvar = (float *) malloc (slit * sizeof(float));
    float yvarsum = 0;
	
	for(j=0;j<slit;j++) {
		my += y[ix-slit/2+j]*slitfunc[j];
		intflux += y[ix-slit/2+j];
        yvar[j] = (noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain;
        yvarsum += yvar[j];        
	}
	
	float avgx = 0;
	
    for(j=0;j<slit;j++) {
		mx += x[ix-slit/2+j]*y[ix-slit/2+j]/intflux;		
		avgx += x[ix-slit/2+j]/(float)slit;
        mxerr += ((noise/gain)*(noise/gain) + fabs(y[ix-slit/2+j])/gain); // method 1
        //	mxerr += ((intflux-y[ix-slit/2+j])/(intflux*intflux))*((intflux-y[ix-slit/2+j])/(intflux*intflux))*yvar[j] + (y[ix-slit/2+j]/(intflux*intflux))*(y[ix-slit/2+j]/(intflux*intflux))*yvarsum;  // method 2	
    }		
	
    float avgy = intflux/(float)slit;
    mxerr = sqrt(mxerr*(fabs(avgx - mx)/avgy));  // method 1
    //		mxerr = sqrt(mxerr); // method 2
    
	float depth;
	int detectprob = 0;
	
	depth = my - (y[ix+slit/2] + y[ix-slit/2])/2;
	
	if(depth > noise/gain) {
		detectprob = 1;
	} else {
		detectprob = 0;
	}		
	
	*xmean = mx;
	*ymean = my;
	*xmeanerr = mxerr;
	//	printf("%f\t%f\t%f\n",mx,my,mxerr);
    
	return detectprob;
}


/*** Algorithm to detect orders using an IP function ***/
unsigned operaCCDDetectPeaksWithErrorsUsingIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	float *mxerr = (float *) malloc (np * sizeof(float));	
	memset(mxerr, 0, sizeof(float)*np);		
	
	unsigned slit=0;
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slit++;
		}	
	}
	float *slitfunc = (float *) malloc (slit * sizeof(float));	
	slit=0;
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slitfunc[slit] = ipfunc[i];
			//			printf("%d\t%f\n",slit,slitfunc[slit]);			
			slit++;
		}	
	}
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
        
        float *yvar = (float *) malloc (slit * sizeof(float));
        float yvarsum = 0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]*slitfunc[j];
			intflux += y[i-slit/2+j];
            yvar[j] = (noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain;
            yvarsum += yvar[j];
		}
		
		float avgx = 0;
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux;
			avgx += x[i-slit/2+j]/(float)slit;
			mxerr[i] += ((noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain); // method 1
			//			mxerr[i] += ((intflux-y[i-slit/2+j])/(intflux*intflux))*((intflux-y[i-slit/2+j])/(intflux*intflux))*yvar[j] + (y[i-slit/2+j]/(intflux*intflux))*(y[i-slit/2+j]/(intflux*intflux))*yvarsum; // method 2	
		}		
		
		float avgy = intflux/(float)slit;
		
        mxerr[i] = sqrt(mxerr[i]*(fabs(avgx - mx[i])/avgy));  // method 1
		//		mxerr[i] = sqrt(mxerr[i]); // method 2
		
		//		printf("%d\t%f\t%f\t%f\t%f\t%f\n",i,x[i],y[i],mx[i],my[i],mxerr[i]);
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center and it decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			xmeanerr[nords] = sqrt(mxerr[i]*mxerr[i] + (1/(float)slit)*(1/(float)slit));			
			ymean[nords] = my[i];
			nords++;
		}
		//		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	free(mxerr);
	return nords;
}

/*** Algorithm to detect orders using a Gaussian function ***/
unsigned operaCCDDetectPeaksWithErrorsUsingGaussian(unsigned np, float *x,float *y,float sigma, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	float *mxerr = (float *) malloc (np * sizeof(float));	
	memset(mxerr, 0, sizeof(float)*np);		
	
	unsigned slit;
	
	slit = 4.0*(unsigned)sigma;
	
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
        
        float *yvar = (float *) malloc (slit * sizeof(float));
        float yvarsum = 0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]*exp(-((x[i-slit/2+j] - x[i])*(x[i-slit/2+j] - x[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma); 
			intflux += y[i-slit/2+j];
            yvar[j] = (noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain;
            yvarsum += yvar[j];
		}
		
		float avgx = 0;
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux;
			avgx += x[i-slit/2+j]/(float)slit;
			mxerr[i] += ((noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain); // method 1
			//			mxerr[i] += ((intflux-y[i-slit/2+j])/(intflux*intflux))*((intflux-y[i-slit/2+j])/(intflux*intflux))*yvar[j] + (y[i-slit/2+j]/(intflux*intflux))*(y[i-slit/2+j]/(intflux*intflux))*yvarsum; // method 2	
		}		
        
		float avgy = intflux/(float)slit;
		
        mxerr[i] = sqrt(mxerr[i]*(fabs(avgx - mx[i])/avgy));  // method 1
		//		mxerr[i] = sqrt(mxerr[i]); // method 2
        
		//		printf("%d\t%f\t%f\t%f\t%f\t%f\n",i,x[i],y[i],mx[i],my[i],mxerr[i]);
		
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center or decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
			
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			xmeanerr[nords] = sqrt(mxerr[i]*mxerr[i] + (1/(float)slit)*(1/(float)slit));			
			ymean[nords] = my[i];
			nords++;
		}
		//		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	free(mxerr);
	return nords;
}

unsigned operaCCDDetectPeaksWithErrorsUsingGaussianDouble(unsigned np, double *x,double *y,double sigma, double noise, double gain, double threshold,double *xmean, double *ymean, double *xmeanerr)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	double *my = (double *) malloc (np * sizeof(double));
	memset(my, 0, sizeof(double)*np);
	
	double *mx = (double *) malloc (np * sizeof(double));	
	memset(mx, 0, sizeof(double)*np);	
	
	double *mxerr = (double *) malloc (np * sizeof(double));	
	memset(mxerr, 0, sizeof(double)*np);		
	
	unsigned slit;
	
	slit = 4.0*(unsigned)sigma;
	
	double intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
        
        double *yvar = (double *) malloc (slit * sizeof(double));
        double yvarsum = 0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]*exp(-((x[i-slit/2+j] - x[i])*(x[i-slit/2+j] - x[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma); 
			intflux += y[i-slit/2+j];
            yvar[j] = (noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain;
            yvarsum += yvar[j];
		}
		
		double avgx = 0;
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux;
			avgx += x[i-slit/2+j]/(double)slit;
			mxerr[i] += ((noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain); // method 1
			//			mxerr[i] += ((intflux-y[i-slit/2+j])/(intflux*intflux))*((intflux-y[i-slit/2+j])/(intflux*intflux))*yvar[j] + (y[i-slit/2+j]/(intflux*intflux))*(y[i-slit/2+j]/(intflux*intflux))*yvarsum; // method 2	
		}		
        
		double avgy = intflux/(double)slit;
		
        mxerr[i] = sqrt(mxerr[i]*(fabs(avgx - mx[i])/avgy));  // method 1
		//		mxerr[i] = sqrt(mxerr[i]); // method 2
        
		//		printf("%d\t%f\t%f\t%f\t%f\t%f\n",i,x[i],y[i],mx[i],my[i],mxerr[i]);
		
	}
	
	unsigned npts=0, isitmax=1;
	double rate, depth;
	double detectprob,peakprob;
	
	double x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center or decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
			
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((double)npts/(double)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			xmeanerr[nords] = sqrt(mxerr[i]*mxerr[i] + (1/(double)slit)*(1/(double)slit));			
			ymean[nords] = my[i];
			nords++;
		}
		//		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	free(mxerr);
	return nords;
}

/*** Algorithm to detect orders using a flat top hat function ***/
unsigned operaCCDDetectPeaksWithErrorsUsingTopHat(unsigned np, float *x,float *y,unsigned width, float noise, float gain, float threshold,float *xmean, float *ymean, float *xmeanerr)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	float *mxerr = (float *) malloc (np * sizeof(float));	
	memset(mxerr, 0, sizeof(float)*np);		
	
	unsigned slit;
	
	slit = width;
	
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
        
        float *yvar = (float *) malloc (slit * sizeof(float));
        float yvarsum = 0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]/(float)slit; 
			intflux += y[i-slit/2+j];
            yvar[j] = (noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain;
            yvarsum += yvar[j];
		}
		
		float avgx = 0;
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux;
			avgx += x[i-slit/2+j]/(float)slit;
			mxerr[i] += ((noise/gain)*(noise/gain) + fabs(y[i-slit/2+j])/gain); // method 1
            //			mxerr[i] += ((intflux-y[i-slit/2+j])/(intflux*intflux))*((intflux-y[i-slit/2+j])/(intflux*intflux))*yvar[j] + (y[i-slit/2+j]/(intflux*intflux))*(y[i-slit/2+j]/(intflux*intflux))*yvarsum; // method 2	
		}		
        
		float avgy = intflux/(float)slit;
		
        mxerr[i] = sqrt(mxerr[i]*(fabs(avgx - mx[i])/avgy));  // method 1
        //		mxerr[i] = sqrt(mxerr[i]); // method 2
		
		//		printf("%d\t%f\t%f\t%f\t%f\t%f\n",i,x[i],y[i],mx[i],my[i],mxerr[i]);		
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center or decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			xmeanerr[nords] = sqrt(mxerr[i]*mxerr[i] + (1/(float)slit)*(1/(float)slit));			
			ymean[nords] = my[i];
			nords++;
		}
		//		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	free(mxerr);
	return nords;
}

/*** Algorithm to detect orders using an IP function ***/
unsigned operaCCDDetectPeaksWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	unsigned slit=0;
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slit++;
		}	
	}
	float *slitfunc = (float *) malloc (slit * sizeof(float));	
	slit=0;
	for(i=0;i<nip;i++) {
		if(ipfunc[i]) {
			slitfunc[slit] = ipfunc[i];
			//			printf("%d\t%f\n",slit,slitfunc[slit]);			
			slit++;
		}	
	}
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]*slitfunc[j];
			intflux += y[i-slit/2+j];
		}
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux; 
		}		
		
//		printf("%d\t%f\t%f\t%f\t%f\n",i,x[i],y[i],mx[i],my[i]);
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center and it decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);

		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			ymean[nords] = my[i];
			nords++;
		}
        //printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,x[i],y[i],mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}
	free(my);
	free(mx);
	return nords;
}

/*** Algorithm to detect orders using an IP function ***/
unsigned operaCCDDetectPeaksByXCorrWithIP(unsigned np, float *x,float *y,unsigned nip, float *ipfunc, float noise, float gain, float threshold,float *xmean, float *ymean)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
    float intflux=0;
	
    float mx[MAXESPADONSX];
    float my[MAXESPADONSX];
    float xcorr[MAXESPADONSX];
    float stdev[MAXESPADONSX];
    
    float datatmp[MAXIPPOINTS];
	
    unsigned firstj, lastj;
    
    unsigned halfnip = (unsigned)floor(nip/2);
    
    for(i=0;i<np;i++) {
        
        if( i>=halfnip || i < np-halfnip) {

            firstj = i-halfnip;
            lastj = firstj+nip;
            
            intflux=0;
            
            unsigned k=0;
            my[i] = 0;
            for(j=firstj;j<lastj;j++) {
                datatmp[k] = y[j];
                intflux += y[j];
                my[i] += y[j]*ipfunc[k];
                k++;
            }
            
            stdev[i] = operaArraySigma(nip,datatmp);
            
            mx[i]=0;
            for(j=firstj;j<lastj;j++) {
                mx[i] += x[j]*y[j]/intflux;
            }
            
            xcorr[i] = operaCrossCorrelation_f(nip,ipfunc,datatmp);
            
            // printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,x[i],y[i],mx[i],my[i],stdev[i],xcorr[i]);
        } else {
            xcorr[i] = 0.0;
        }
    }
	
    for(i=0;i<np;i++) {
        if( i>=nip || i < np-nip) {
            
            firstj = i-2*nip/3;
            lastj = firstj+2*2*nip/3;
            
            float locmax = -BIG;
            unsigned jmax = np+1;
            for(j=firstj;j<lastj;j++) {
                if(xcorr[j] > locmax) {
                    locmax = xcorr[j];
                    jmax = j;
                }
            }
            
            if(jmax==i && xcorr[i]>threshold && stdev[i] > noise/gain) {
                if(nords > MAXORDERS) { // run out of memory
                    return nords;
                }
                xmean[nords] = mx[i];
                ymean[nords] = my[i];
                //printf("%d\t%lf\t%lf\n",i,xmean[nords],ymean[nords]);
                nords++;
            }
        }
        //printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,x[i],y[i],mx[i],my[i],stdev[i]);
    }
    
	return nords;
}

/*** Algorithm to detect orders using a Gaussian function ***/
unsigned operaCCDDetectPeaksWithGaussian(unsigned np, float *x,float *y,float sigma, float noise, float gain, float threshold,float *xmean, float *ymean)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	unsigned slit;
	
	slit = 4.0*(unsigned)sigma;
	
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]*exp(-((x[i-slit/2+j] - x[i])*(x[i-slit/2+j] - x[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma); 
			intflux += y[i-slit/2+j];
		}
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux; 
		}		
		
		//	printf("%d\t%lf\t%lf\n",i,mx[i],my[i]);
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center or decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
			
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			ymean[nords] = my[i];
			nords++;
		}
		//	printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	return nords;
}

/*** Algorithm to detect orders using a flat top hat function ***/
unsigned operaCCDDetectPeaksWithTopHat(unsigned np, float *x,float *y,unsigned width, float noise, float gain, float threshold,float *xmean, float *ymean)
{
	/* the algorithm below detects the orders, their photocenter and mean values */
	unsigned i,j,nords=0;
	
	float *my = (float *) malloc (np * sizeof(float));
	memset(my, 0, sizeof(float)*np);
	
	float *mx = (float *) malloc (np * sizeof(float));	
	memset(mx, 0, sizeof(float)*np);	
	
	unsigned slit;
	
	slit = width;
	
	float intflux;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		intflux=0;
		for(j=0;j<slit;j++) {
			my[i] += y[i-slit/2+j]/(float)slit; 
			intflux += y[i-slit/2+j];
		}
		
		for(j=0;j<slit;j++) {
			mx[i] += x[i-slit/2+j]*y[i-slit/2+j]/intflux; 
		}		
		
		//	printf("%d\t%lf\t%lf\n",i,mx[i],my[i]);
	}
	
	unsigned npts=0, isitmax=1;
	float rate, depth;
	float detectprob,peakprob;
	
	float x0,y0,xf,yf,slope0,slope1;
	
	for(i=slit/2;i<(np-slit/2);i++) {
		
		x0 = mx[i-slit/2];
		y0 = my[i-slit/2];
		xf = mx[i+slit/2];
		yf = my[i+slit/2];
		
		npts=0;
		isitmax=1;
		
		for(j=0;j<slit;j++) {
			
			slope0 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j] - x0); 
			slope1 = y0 + ((yf-y0)/(xf-x0))*(mx[i-slit/2+j+1] - x0);
			
			//check whether it is a local maximum
			if(my[i-slit/2+j] > my[i])
				isitmax = 0;
			
			//calculate the growth rate
			rate = (my[i-slit/2+j+1]-slope0) - (my[i-slit/2+j]-slope1);
			
			//count points whenever it grows before center or decreases after center
			if(j<slit/2) {
				if(rate > 0) {
					npts++;
				}
			} else if (j>=slit/2) {
				if(rate < 0) {
					npts++;
				} 
			}
		}	
		
		depth = my[i] - ((my[i+slit/2] + my[i-slit/2])/2);
		
		if(depth > noise/gain) {
			detectprob = 1;
		} else {
			detectprob = 0;
		}	
		
		peakprob = isitmax*detectprob*((float)npts/(float)(slit));
		
		if(peakprob > threshold) {
			xmean[nords] = mx[i];
			ymean[nords] = my[i];
			nords++;
		}
		//	printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%u\t%u\n",i,mx[i],my[i],peakprob, depth, noise/gain, npts, slit);
	}	
	free(my);
	free(mx);
	return nords;
}


void operaMedianWidthFromSetOfLines(unsigned np, float *mx, float *my, float *myerr, unsigned nlines, float *xlines, float *ylines, float *medianWidth){
    
    unsigned slit = 6.0*(unsigned)(*medianWidth) + 2;    
    // The above slit is define according to the following criteria:
    // 4*sigma points contain >99% of the integrated flux of a gaussian. 
    // The remaining 2*sigmas+2 points is used to estimate the background
    // The 2 extra points are meant to make sure one would have enough points for background estimation
    // which uses a linear fit with 2 parameters (at least 2 points univocally define a slope). 
    
    unsigned *ii = (unsigned *) malloc (nlines * sizeof(unsigned));
    
    unsigned lasti = 0;
    float xminseparation = (unsigned)(*medianWidth);
    
    for(unsigned j=0;j<nlines;j++) {
        
        for(unsigned i=lasti;i<np;i++) {
            
            float xseparation = fabs(mx[i] - xlines[j]);
            
            if(xseparation < xminseparation) {
                xminseparation = xseparation;
                ii[j] = i;
            } else if (xseparation > xminseparation && xminseparation < (unsigned)(*medianWidth)) {
                lasti = i;
                xminseparation = (unsigned)(*medianWidth); 
                break;
            }
        }
    }
	/*
	 for(unsigned i=0;i<np;i++){
	 printf("%d\t%f\t%f\t%f\n",i,mx[i],my[i],myerr[i]);
	 }
	 for(unsigned j=0;j<nlines;j++) {
	 printf("%u\t%f\t%f\t%f\n",ii[j],xlines[j],ylines[j]);
	 }
     */
	
	/*  
	 float *bkg, *xbkg;
	 bkg = (float *) malloc (2*slit * sizeof(float));
	 xbkg = (float *) malloc (2*slit * sizeof(float));
	 float a,b,abdev;
	 
	 unsigned nbkg=0, nin=0;    
	 
	 
	 if((x[i] >= xmean[j] - (float)slit && x[i] <= xmean[j] - (float)slit/2) ||
	 (x[i] >= xmean[j] + (float)slit/2 && x[i] <= xmean[j] + (float)slit)) {
	 xbkg[nbkg] = x[i];
	 bkg[nbkg] = y[i];			
	 nbkg++;    
	 
	 ladfit(xbkg,bkg,nbkg,&a,&b,&abdev);
	 
	 */
	for(unsigned j=0;j<nlines;j++) {
        double a=ylines[j],x0=xlines[j],sig=(*medianWidth);
        double ea,ex0,esig;
        double chisqr;
        
        unsigned imin = ii[j] - slit/2;
        unsigned imax = ii[j] + slit/2;
        
        ////DT May 16 2014 NOTE: imin is unsigned so this always fails -- an unsigned cn't be less than zero: if(imin < 0) {imin = 0;}
        if(slit/2 > ii[j]) {imin = 0;}
        if(imax >= np) {imax = np;}
        
        if(fabs(imax-imin) < slit) {
            continue;
        }
        
        double *x = (double *) malloc (slit * sizeof(double));
        double *y = (double *) malloc (slit * sizeof(double));
        double *ey = (double *) malloc (slit * sizeof(double));
        
        unsigned npts = 0;
		for(unsigned i=imin;i<imax;i++) {
			x[npts] = (double)mx[i];
            y[npts] = (double)my[i];
            ey[npts] = (double)myerr[i];            
            if(npts == slit) {
                break;
            }
            npts++;
		}
		
        operaMPFitGaussian(npts, x, y, ey,&a,&ea,&x0,&ex0,&sig,&esig,&chisqr);
        
		//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",j,a,ea,x0,ex0,sig,esig,chisqr);
        
        npts = 0;
		for(unsigned i=imin;i<imax;i++) {
            double ymodel = a*exp(-((double)mx[i]-x0)*((double)mx[i]-x0)/(2*sig*sig));
            printf("%u\t%f\t%f\t%f\t%f\n",i,(double)mx[i],(double)my[i],ymodel,(double)myerr[i]);
            if(npts == slit) {
                break;
            }
            npts++;
		}        
	}    
	
}
