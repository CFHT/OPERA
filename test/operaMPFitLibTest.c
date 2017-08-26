/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
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

/*! \file operaMPFitLibTest.cpp */

/*! 
 * operaMPFitLibTest
 * \author Eder Maritoli
 * \brief Perform various tests on the MP Fit functions.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char** argv)
{
	unsigned i,np=60;	
	double *x, *y, *ey;
	x = (double*) malloc(np * sizeof(double));
	y = (double*) malloc(np * sizeof(double));	
	ey = (double*) malloc(np * sizeof(double));	 
	
/*** Test polynomial fit ***/
//-- generate simulated data
	srand(time(NULL));
	for(i=0;i<np;i++) {
		x[i] = (double)i;
		y[i] = (5 - 3.8*x[i] + 2.5*x[i]*x[i] - 0.5*x[i]*x[i]*x[i]) + operaUniformRand(-0.25, 0.25);
		ey[i] = fabs((5 - 3.8*x[i] + 2.5*x[i]*x[i] - 0.5*x[i]*x[i]*x[i]) - y[i]);
	}
	printf("---------------------------------------------------------------\n");	
	printf("\nInput polynomial:\n");	
	printf("f(x) = 5 - 3.8*x + 2.5*x**(2) - 0.5*x**(3)\n");
//--------------------------
	int npar = 4;
	double par[4] = {1,1,1,1};
	double epar[4] = {0,0,0,0};	
	double chisqr = 1e30;
	
	operaMPFitPolynomial(np, x, y, ey, npar, par, epar, &chisqr);
	
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

	printf("\nResult parameters:\n");	
	for(i=0;i<npar;i++) {
		printf("p[%d] = %.5lf +/- %.5lf (%.1lf %%)\n",i,par[i],epar[i],100*epar[i]/(double)fabs((float)par[i]));
	}
	printf("---------------------------------------------------------------\n");		
	
/*** Test gaussian fit ***/
//-- generate simulated data
	srand(time(NULL));
	for(i=0;i<np;i++) {
		x[i] = (double)i;
		y[i] = (25*exp(-(x[i] - 13)*(x[i] - 13)/(2*7*7))) + operaUniformRand(-1, 1);
		ey[i] = 1.0;
	}
	printf("---------------------------------------------------------------\n");	
	printf("\nInput Gaussian:\n");	
	printf("f(x) = 25*exp(-(x[i] - 13)*(x[i] - 13)/(2*7*7))\n");	
	//--------------------------
	double a=1,x0=1,sig=1;
	double ea=0,ex0=0,esig=0;
    chisqr = 1e30;
    
	operaMPFitGaussian(np, x, y, ey,&a,&ea, &x0,&ex0, &sig,&esig, &chisqr);
	
	printf("\nFit Gaussian:\n");
	printf("f(x) =%.2lf*exp(-(x-%.2lf)*(x-%.2lf)/(2*%.2lf*%.2lf))\n",a,x0,x0,sig,sig); 	
	printf("Chi-square = %lf\n",chisqr);

	printf("\nResult parameters:\n");	
	printf("Amplitude = %.5lf +/- %.5lf (%.1lf %%)\n",a,ea,100*ea/(double)fabs((float)a));
	printf("Center = %.5lf +/- %.5lf (%.1lf %%)\n",x0,ex0,100*ex0/(double)fabs((float)x0));
	printf("Sigma = %.5lf +/- %.5lf (%.1lf %%)\n",sig,esig,100*esig/(double)fabs((float)sig));		
	printf("---------------------------------------------------------------\n");	
	
/*** Test multiple gaussian fit ***/
#define NGAUSSIANS 2 
    unsigned ngauss = NGAUSSIANS; // number of gaussians to fit
//-- generate simulated data
	for(i=0;i<np;i++) {
		x[i] = (double)i;
		y[i] = (10*exp(-(x[i] - 35)*(x[i] - 35)/(2*5*5))) + (15*exp(-(x[i] - 13)*(x[i] - 13)/(2*5*5))) + operaUniformRand(-1, 1);
//		y[i] = (10*exp(-(x[i] - 35)*(x[i] - 35)/(2*5*5))) + operaUniformRand(-1, 1); // uncomment this line to test ngauss=1
		ey[i] = 1.0;
//        printf("%lf\t%lf\t%lf\n",x[i],y[i],ey[i]); // uncomment this line to print simulated data
	}
    
	printf("---------------------------------------------------------------\n");	
	printf("\nInput Gaussian:\n");	
	printf("f(x) = (10*exp(-(x - 35)*(x - 35)/(2*5*5))) + (15*exp(-(x - 13)*(x - 13)/(2*5*5)))\n");	
	//--------------------------
	double *am = (double*) malloc(ngauss * sizeof(double));
    double *x0m = (double*) malloc(ngauss * sizeof(double));
    double *sigm = (double*) malloc(ngauss * sizeof(double));
    
	double *eam = (double*) malloc(ngauss * sizeof(double));
    double *ex0m = (double*) malloc(ngauss * sizeof(double));
    double *esigm = (double*) malloc(ngauss * sizeof(double));  
    
    for(i=0;i<ngauss;i++){
        am[i] = 10;
        sigm[i] = 10;
        x0m[i] = 20;
        eam[i] = 0;
        ex0m[i] = 0;
        esigm[i] = 0;
    }
    chisqr = 1e30;
    
    //operaLMFitMultipleGaussian(np, x, y,ngauss,am,x0m,sigm,&chisqr); // this doesn't work. I couldn't figure out why.
    
    operaMPFitMultipleGaussian(np, x, y, ey,ngauss,am,eam,x0m,ex0m,sigm,esigm,&chisqr);
    
	printf("\nFit Gaussian:\n");
    printf("f(x) = ");
    for(i=0;i<ngauss;i++){    
        printf("+ %.2lf*exp(-(x-%.2lf)*(x-%.2lf)/(2*%.2lf*%.2lf))\n",am[i],x0m[i],x0m[i],sigm[i],sigm[i]); 	
    }
	printf("Chi-square = %lf\n",chisqr);
    
	printf("\nResult parameters:\n");	
    for(i=0;i<ngauss;i++){      
        printf("Amplitude[%d] = %.5lf +/- %.5lf (%.1lf %%)\n",i,am[i],eam[i],100*eam[i]/(double)fabs((float)am[i]));
        printf("Center[%d] = %.5lf +/- %.5lf (%.1lf %%)\n",i,x0m[i],ex0m[i],100*ex0m[i]/(double)fabs((float)x0m[i]));
        printf("Sigma[%d] = %.5lf +/- %.5lf (%.1lf %%)\n\n",i,sigm[i],esigm[i],100*esigm[i]/(double)fabs((float)sigm[i]));	
    }        
	printf("---------------------------------------------------------------\n");	    
    
    free(am);
    free(x0m);
    free(sigm);
    free(eam);
    free(ex0m);
    free(esigm); 
	
    return EXIT_SUCCESS;
}
