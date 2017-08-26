/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFitLibTest
 Version: 1.0
 Description: This module is for testing and exemplify the operaFit 
 library.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: eder@cfht.hawaii.edu
 
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFit.h"
#include "libraries/operaStats.h"

/*! \file operaFFitLibTest.c */

/*! 
 * operaFFitLibTest
 * \author Doug Teeple
 * \brief Perform various tests on the fitting routines.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */
int main(int argc, char** argv)
{
    unsigned i,np=30;	
    double *x, *y;
	x = (double*) malloc(np * sizeof(double));
	y = (double*) malloc(np * sizeof(double));	
	
/*** Test polynomial fit ***/
//-- generate simulated data
	srand(time(NULL));
	for(i=0;i<np;i++) {
		x[i] = (double)i;
		y[i] = (5 - 3.8*x[i] + 2.5*x[i]*x[i] - 0.5*x[i]*x[i]*x[i]) + operaUniformRand(-0.25, 0.25);
	}
	printf("\nInput polynomial:\n");	
	printf("f(x) = 5 - 3.8*x + 2.5*x**(2) - 0.5*x**(3)\n");
//--------------------------
	int npar = 4;
	double par[4] = {1,1,1,1};
	double chisqr;
	
	operaLMFitPolynomial(np, x, y, npar, par, &chisqr);

	printf("\nFit polynomial:\n");	
	printf("f(x) =");
	for(i=0;i<npar;i++) {
		if(i==0)
			printf(" %lf",par[i]);
		else if (i==1)
			printf(" + %lf*x",par[i]);
		else
			printf(" + %lf*x**(%d)",par[i],i);			
	}
	printf("\nChi-square = %lf\n",chisqr);

/*** Test gaussian fit ***/
//-- generate simulated data
	srand(time(NULL));
	for(i=0;i<np;i++) {
		x[i] = (double)i;
		y[i] = (25*exp(-(x[i] - 13)*(x[i] - 13)/(2*7*7))) + operaUniformRand(-1, 1);
	}
	printf("\nInput Gaussian:\n");	
	printf("f(x) = 25*exp(-(x[i] - 13)*(x[i] - 13)/(2*7*7))\n");	
	//--------------------------
	double a=1,x0=1,sig=1;
	
	operaLMFitGaussian(np, x, y, &a, &x0, &sig, &chisqr);
	
	printf("\nFit Gaussian:\n");
	printf("f(x) =%.2lf*exp(-(x-%.2lf)*(x-%.2lf)/(2*%.2lf*%.2lf))\n",a,x0,x0,sig,sig); 	
	printf("Chi-square = %lf\n",chisqr);

	
/*** Test 2D functions ***/
	
	unsigned nx = 30, ny = 30;	
	unsigned mp;
	mp = nx*ny;
	double *xx, *yy, *z;
	xx = (double*) malloc(mp * sizeof(double));
	yy = (double*) malloc(mp * sizeof(double));	
	z = (double*) malloc(mp * sizeof(double));		

	/*** Test 2D Gaussian ***/
	//-- generate simulated data
	double aa=20,xx0=13,yy0=11,sigx=5,sigy=4;	
	srandom(time(NULL));
	int k,j;
	k=0;
	for(i=0;i<nx;i++) {
		for(j=0;j<ny;j++) {
			xx[k] = (double)i;			
			yy[k] = (double)j;
			z[k] = aa*exp(-((xx[k]-xx0)*(xx[k]-xx0)/(2*sigx*sigx)+(yy[k]-yy0)*(yy[k]-yy0)/(2*sigy*sigy))) + operaUniformRand(-3.0, 3.0);
			k++;
//			printf("%lf\t%lf\t%lf\n",xx[k-1],yy[k-1],z[k-1]);
		}	
	}
	printf("\nInput 2D Gaussian:\n");	
	printf("f(x,y) = %.2f*exp(-((x - %.2f)*(x - %.2f)/(2*%.2f*%.2f)" 
				 " + (y - %.2f)*(y - %.2f)/(2*%.2f*%.2f)))\n",
				 aa,xx0,xx0,sigx,sigx,yy0,yy0,sigy,sigy);	
	
	aa=1;xx0=1;yy0=1;sigx=10;sigy=10;
	
	//--------------------------
	operaLMFit2DGaussian(mp,xx,yy,z,&aa,&xx0,&yy0,&sigx,&sigy,&chisqr);

	printf("\nFit 2D Gaussian:\n");	
	printf("f(x,y) = %.2f*exp(-((x - %.2f)*(x - %.2f)/(2*%.2f*%.2f)" 
				 " + (y - %.2f)*(y - %.2f)/(2*%.2f*%.2f)))\n",
				 aa,xx0,xx0,sigx,sigx,yy0,yy0,sigy,sigy);	
	printf("Chi-square = %lf\n",chisqr);	
	
/*** Test polynomial surface fit ***/
	//-- generate simulated data
	npar = 5;
	double pars[5] = {1,2,0.5,-0.03,0.2};
	unsigned kk,jj,ii;
	int n = npar/2;
	int m = npar - n;
	double fxy;

	k=0;
	
	srandom(time(NULL));
	
	for(i=0;i<nx;i++) {
		for(j=0;j<ny;j++) {
			xx[k] = (double)i;			
			yy[k] = (double)j;
			
			fxy=0;
			kk=0;
			for(jj=0;jj<m;jj++){
				for(ii=0;ii<n;ii++){	
					fxy += pars[kk++]*pow(xx[k],(double)ii)*pow(yy[k],(double)jj);
				}
			}		
			z[k++] = fxy + operaUniformRand(-2.0, 2.0);
		//	printf("%lf\t%lf\t%lf\n",xx[k-1],yy[k-1],z[k-1]);			
		}	
	}
	
	kk=0;
	printf("\nInput 2D polynomial:\n");	
	printf("f(x,y) = ");
	for(jj=0;jj<m;jj++){
		for(ii=0;ii<n;ii++){	
			printf("+ %lf*x**(%d)*y**(%d) ",pars[kk],ii,jj);
			pars[kk] = pars[kk] + operaUniformRand(-2.0, 2.0);
			kk++;
		}
	}		
	printf("\n");
				 
	//--------------------------
	operaLMFit2DPolynomial(mp,xx,yy,z,npar,pars,&chisqr);
	
	kk=0;
	printf("\nFit 2D polynomial:\n");	
	printf("f(x,y) = ");
	for(jj=0;jj<m;jj++){
		for(ii=0;ii<n;ii++){	
			printf("+ %lf*x**(%u)*y**(%u) ",pars[kk],ii,jj);
			kk++;
		}
	}		
	printf("\n");	
	printf("Chi-square = %lf\n",chisqr);

	
	//--------------------------
	/*** Test lad fit ***/	
	printf("\nTest ladfit:\n");
	float *xf, *yf;
	xf = (float*) malloc(np * sizeof(float));
	yf = (float*) malloc(np * sizeof(float));	
	
	//-- generate simulated data
	srandom(time(NULL));
	for(i=0;i<np;i++) {
		xf[i] = (float)i;
		x[i] = (double)i;
		yf[i] = 10 + 0.55*xf[i] + operaGaussRand(0, 1);
		y[i] = (double)yf[i];
//		printf("%lf\t%lf\n",x[i],y[i]);
	}
	printf("Input Line:\n");	
	printf("f(x) = 10 + 0.55 x\n");	
	//--------------------------

	float meda=0,medb=1,abdev;
	ladfit(xf, yf, np, &meda, &medb, &abdev);	

	int mednpar = 2;
	double medpar[2] = {1,1};	// used to be 1, fixed! DT Nov 1 2011
	
	operaLMFitPolynomial(np, x, y, mednpar, medpar, &chisqr);	
	
	printf("\nLinear robust ladfit :\n");
	printf("f(x) = %.2f + %.2f*x \n",meda,medb); 	
	printf("Chi-square = %f\n",abdev);	
	printf("\nLinear LM polynomial fit :\n");	
	printf("f(x) = %.2lf + %.2lf*x \n",medpar[0],medpar[1]); 	
	printf("Chi-square = %lf\n",chisqr);	
	//--------------------------
	
	/*** Test spline interpolation ***/		
	printf("\nTest spline interpolation:\n");	
	unsigned npsp = 75;
	float *xsp, *ysp;
	xsp = (float*) malloc(npsp * sizeof(float));
	ysp = (float*) malloc(npsp * sizeof(float));		

	for(i=0;i<npsp;i++) {
		xsp[i] = (float)i*(float)(np-1)/(float)npsp ;
	}	
	operaFitSpline(np,xf,yf,npsp,xsp,ysp);
	
/*	for(i=0;i<npsp;i++) {
		printf("%f\t%f\n",xsp[i],ysp[i]);
	}		
	*/
	//operaFit2DSpline(unsigned nxin, float *xin, unsigned nyin, float *yin,double *fxyin, unsigned nxout, float *xout, unsigned nyout, float *yout, float *fxyout);
	
	return EXIT_SUCCESS;
}
