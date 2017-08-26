/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaMathLibTest
 Version: 1.0
 Description: Test the operaMath library functions.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaMath.h"
#include "libraries/operaFit.h"

/*! \file operaMathLibTest.c */

/*! 
 * operaMathLibTest
 * \author Eder Martioli
 * \brief Test the operaMath library functions.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
	
	//--- Derivative of polynomial
	int npar = 4;
	double par[4] = {1,2,3,1};
	printf("\n********************\n");		
	printf("f(x) =");
	for(unsigned i=0;i<npar;i++) {
		if(i==0)
			printf(" %.0lf",par[i]);
		else if (i==1)
			printf(" + %.0lf*x",par[i]);
		else
			printf(" + %.0lf*x**(%d)",par[i],i);			
	}
	printf("\n");
	
	printf("f(2) = %lf\n",PolynomialFunction(2,par,npar));
	printf("\n");		
	printf("f'(x) =");
	for(unsigned i=1;i<npar;i++) {
		if(i==1)
			printf(" %.0lf",(double)i*par[i]);
		else if (i==2)
			printf(" + %.0lf*x",(double)i*par[i]);
		else
			printf(" + %.0lf*x**(%d)",(double)i*par[i],i-1);			
	}
	printf("\n");	
	printf("f'(2) = %lf\n",DiffPolynomialFunction(2,par,npar));  
	
	printf("********************\n\n");	
	//---------------------------------
	
	float a=0.5,b=2.3;
	float polylength;
	polylength = LengthofPolynomial(a,b,par,npar); // this method uses Simpson's rule
												   //	polylength = PolyLengthTrapzd(par,npar,a,b,1); // this method uses the trapezoidal rule
	printf("********************\n");
	printf("Length = int_a^b{sqrt[1+f'(x)^2]dx}\n");	
	printf("Length of polynomial from x=%.1f to x=%.1f is: %f\n",a,b,polylength); 	
	printf("********************\n");
	
	
	
	
	return EXIT_SUCCESS;
}  






