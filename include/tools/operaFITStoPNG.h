#ifndef OPERAFITSTOPNG_H
#define OPERAFITSTOPNG_H
/******************************************************************
 ****                  MODULE FOR OPERA v1.0                   ****
 ******************************************************************
 Module name: operaFITStoPNG
 Version: 1.0
 Description:  Convert an FITS image to PNG.
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

/*! \brief  Convert an FITS image to PNG. */
/*! \file operaFITStoPNG.h */
/*! \ingroup tools */

/* prototypes */

#define histDepth 256
#define index(a,x,y,w) a[(x)+((y)*(w))]
#define ELEM_SWAP(a,b) {register float t=(a);(a)=(b);(b)=t; }

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} colour;

enum outputkinds {
	E_JPG, E_PNG, E_MPG
} outputkinds;

static colour lingray(float val, float z1, float z2, char neg);
static colour lingraycutcol(float val, float z1, float z2, char neg);
static colour irafloggray(float val, float z1, float z2, char neg);
static colour irafloggraycutcol(float val, float z1, float z2, char neg);
static colour plinlog(float val, float z1, float z2, char neg);
static colour logcontours(float val, float z1, float z2, char neg);
static colour loghsv(float val, float z1, float z2, char neg);
static colour lincol1(float val, float z1, float z2, char neg);
static colour lingrayupena(float val, float z1, float z2, char neg);
static colour lingrayupenalog(float val, float z1, float z2, char neg);

static colour HSVtoRGB(float h, float s, float v);
static colour rainbow(float x);

static void GetMedian(const float *inarr, unsigned int n, float *a1, float *a2);
static void GetMinAndMax(const float *inarr, unsigned int n, float *min, float *max);
static float *shrink(const float *image, const int dimx, const int dimy, const int ratio, const float bias, float *smallimage);
static void histogramEqualize(float *image, unsigned int width, unsigned int height, float min, float max);
	
static void printUsageSyntax();

#endif
