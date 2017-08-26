/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirpreview
 Version: 1.0
 Description: Tool to create a wircam thumbnail 
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
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

/* \file wirpreview.cpp */

using namespace std;

/*!
 * wirpreview
 * \author Doug Teeple
 * \brief Tool to create a wircam thumbnail.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOutput
 * \return EXIT_STATUS
 * \ingroup tools
 */

#include <getopt.h>

#include "fitsio.h"
#include "png.h"
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaWIRCamImage.h"			// WIRCAM_BIAS
#include "libraries/operaImage.h"
#include "libraries/operaLib.h"					// systemf
#include "libraries/operaStats.h"				// median

#define MEANSHRINK
#define MAX_FILEPATH_SIZE 1024
#define histDepth 256

#define index(a,x,y,w) a[(x)+((y)*(w))]

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} colour;

enum outputkinds {
	E_JPG, E_PNG, E_MPG
} outputkinds;

#ifdef MEANSHRINK
inline float *shrink(const float *image, const int dim, const int ratio, float *smallimage) {
	int y, x, locx, xx, yy;
	const int smallsize = dim / ratio;
	int locy = 0;
	float mean;
	for (y=0; y<smallsize; y++) {			// iterate line-by-line through the thumbnail
		locx = 0;
		for (x=0; x<smallsize; x++) {
			mean = 0.0;							// find the mean of the subwindow in image
			for (yy=locy; yy<(locy+ratio); yy++) {
				for (xx=locx; xx<(locx+ratio); xx++) {
					mean = (index(image,xx,yy,dim)-WIRCAM_CHIPBIAS+mean)/2.0;
				}
			}
			index(smallimage,x,y,smallsize) = mean;
			locx += ratio;
		}
		locy += ratio;
	}
	return smallimage;
}
#else
inline float *shrink(const float *image, const int dim, const int ratio, float *smallimage) {
	int y, x, locx;
	int locy = 0;
	const int smallsize = dim / ratio;
	const int subindex = ratio / 2;
	for (y=0; y<smallsize; y++) {					// iterate line-by-line through the thumbnail
		locx = 0;
		for (x=0; x<smallsize; x++) {
			index(smallimage,x,y,smallsize) = index(image,locx+subindex/2,locy+subindex/2,dim)-7000.0;
			locx += ratio;
		}
		locy += ratio;
	}
	return smallimage;
}
#endif
// ztrans functions

/*********************************************/
/*  HISTOGRAM EQUALIZATION of A GRAY IMAGE   */
/*********************************************/
static void histogramEqualize(float *image, unsigned int width, unsigned int height, float z1, float z2)
{	
	unsigned int hist[histDepth];
	unsigned short *img_data = NULL;
	unsigned short *img_base;
	float *image_base = image;
	float s_hist_eq[histDepth]={0.0}, sum_of_hist[histDepth]={0.0};
	unsigned int i, npixels;
	const float scale = (float)(histDepth-1)/(z2-z1);
	
	npixels = width * height;
	img_base = img_data = (unsigned short *)malloc(npixels*sizeof(unsigned short));
	for(i=0; i<histDepth; i++)
	{
		hist[i] = 0;
	}
	for(i=0; i<npixels; i++)
	{
		*img_data = (unsigned short)((*image-z1)*scale);
		if (*img_data < histDepth) {
			hist[*img_data++]++;
		} else {
			fprintf(stderr, "Invalid access to hist: img_data=%d image=%4.2f z1=%4.2f z2=%4.2f scale=%4.2f\n", *img_data++, *image, z1, z2, scale);
		}
		image++;
	}
	for (i=0; i<histDepth; i++)  // pdf of image
	{
		s_hist_eq[i] = (double)hist[i]/(double)npixels;
	}
	sum_of_hist[0] = s_hist_eq[0];
	for (i=1; i<histDepth; i++)	 // cdf of image
	{
		sum_of_hist[i] = sum_of_hist[i-1] + s_hist_eq[i];
	}
	img_data = img_base;
	image = image_base;
	for(i=0; i<npixels; i++)
	{
		*image++ = (float)round( sum_of_hist[*img_data++] * (float)(histDepth-1) );
	}
	free(img_base);
}

static colour lingray(float val, float z1, float z2, char neg) {
	colour col;
	if (val < z1) {
		col.r = 0;
	} else if (val > z2) {
		col.r = 255;
	} else {
		col.r = (unsigned char) floor(255.9 * (val-z1)/(z2-z1));
	}
	if (neg == 'n'){
		col.r = 255 - col.r;
	}
	col.g = col.r;
	col.b = col.r;
	return(col);
}

static colour lingraycutcol(float val, float z1, float z2, char neg) {
	colour col;
	if (val < z1) {
		col.r = 0; col.g = 0; col.b = 255;
	} else if (val > z2) {
		col.r = 255; col.g = 0; col.b = 0;
	} else {
		col.r = (unsigned char) floor(255.9 * (val-z1)/(z2-z1));
		if (neg == 'n'){
			col.r = 255 - col.r;
		}
		col.g = col.r;
		col.b = col.r;
	}
	return(col);
}

static colour irafloggray(float val, float z1, float z2, char neg) {
	colour col;
	float x = 1.0;
	if (val < z1) {
		col.r = 0;
	} else if (val > z2) {
		col.r = 255;
	} else {
		x = 1.0 + 1000.0 * ((val-z1)/(z2-z1));
		x = log10(x);
		col.r = (unsigned char) floor(255.9 * x/3.0);
	}
	if (neg == 'n'){
		col.r = 255 - col.r;
	}
	col.g = col.r;
	col.b = col.r;
	return(col);
}

static colour irafloggraycutcol(float val, float z1, float z2, char neg) {
	colour col;
	float x = 1.0;
	if (val < z1) {
		col.r = 0; col.g = 0; col.b = 255;
	} else if (val > z2) {
		col.r = 255; col.g = 0; col.b = 0;
	} else {
		x = 1.0 + 1000.0 * ((val-z1)/(z2-z1));
		x = log10(x);
		col.r = (unsigned char) floor(255.9 * x/3.0);
		if (neg == 'n'){
			col.r = 255 - col.r;
		}
		col.g = col.r;
		col.b = col.r;	
	}
	return(col);
}

static colour plinlog(float val, float z1, float z2, char neg) {
	colour col;
	float x = 1.0;
	unsigned char sep = 80;
	if (val < z1) {
		col.r = 0; col.g = 0; col.b = 0;
	} else if (val > z2) {
		col.r = 255; col.g = 255; col.b = 255;
	} else {
		if (val < -z1) {
			col.r = (unsigned char)floor((sep + 0.9) * (val-z1)/(-z1*2.0));
		} else {
			x = 1.0 + 1000.0 * ((val+z1)/(z2+z1));
			x = log10(x);
			col.r = (unsigned char)(sep + 1 + (unsigned char)floor((255-sep-0.1) * x/3.0));
		}
		if (neg == 'n'){
			col.r = 255 - col.r;
		}
		col.g = col.r;
		col.b = col.r;	
	}
	return(col);
}

static colour logcontours(float val, float z1, float z2, char neg) {
	colour col;
	float x;
	if (val < z1) {
		col.r = 255; col.g = 218; col.b = 15;
	} else if (val > z2) {
		col.r = 255; col.g = 218; col.b = 15;
	} else {
		x = 1.0 + 1000.0 * ((val-z1)/(z2-z1));
		x = log10(x);
		col.r = (unsigned char) floor(255.9 * x/3.0);
		col.g = col.r;
		col.b = col.r;
		if (fabs(val-4.0) < 0.15) {			// courbe rouge a 4.0
			if (col.r < 215) col.r = col.r + 40;
			if (col.g > 20) col.g = col.g - 20;
			if (col.b > 20) col.b = col.b - 20;
		}
		if (fabs(val-1.0) < 0.06) {			// courbe verte a 1.0
			if (col.r > 60) col.r = col.r - 60;
			if (col.g < 235) col.g = col.g + 20;
			if (col.b > 60) col.b = col.b - 60;
		}
		if (fabs(val-0.5) < 0.03) {			// courbe bleue a 0.5
			if (col.r > 10) col.r = col.r - 10;
			if (col.g > 10) col.g = col.g - 10;
			if (col.b < 225) col.b = col.b + 30;
		}
		
	}	
	return(col);
}

// This function is originally from http://www.cs.rit.edu/~ncs/color/t_convert.html, by Eugene Vishnevsky
// I have adapted it to use the colour structure, and changed a few details.
//	h is from 0 to 360 (hue)
//	s from 0 to 1 (saturation)
//	v from 0 to 1 (brightness)

static colour HSVtoRGB(float h, float s, float v ) {
	int i;
	float f, p, q, t;
	colour col;
	if( s == 0 ) {		// grey
		col.r = (unsigned char) floor(255.9 * v);
		col.g = col.r;
		col.b = col.r;
		return(col);
	}
	h /= 60;			// sector 0 to 5
	i = (int)floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	switch( i ) {
		case 0:
			col.r = (unsigned char) floor(255.9 * v);
			col.g = (unsigned char) floor(255.9 * t);
			col.b = (unsigned char) floor(255.9 * p);
			break;
		case 1:
			col.r = (unsigned char) floor(255.9 * q);
			col.g = (unsigned char) floor(255.9 * v);
			col.b = (unsigned char) floor(255.9 * p);
			break;
		case 2:
			col.r = (unsigned char) floor(255.9 * p);
			col.g = (unsigned char) floor(255.9 * v);
			col.b = (unsigned char) floor(255.9 * t);
			break;
		case 3:
			col.r = (unsigned char) floor(255.9 * p);
			col.g = (unsigned char) floor(255.9 * q);
			col.b = (unsigned char) floor(255.9 * v);
			break;
		case 4:
			col.r = (unsigned char) floor(255.9 * t);
			col.g = (unsigned char) floor(255.9 * p);
			col.b = (unsigned char) floor(255.9 * v);
			break;
		default:		// case 5:
			col.r = (unsigned char) floor(255.9 * v);
			col.g = (unsigned char) floor(255.9 * p);
			col.b = (unsigned char) floor(255.9 * q);
			break;
	}
	return(col);
}

static colour loghsv(float val, float z1, float z2, char neg) {
	colour col;
	float x;
	float h, s, v;
	if (val < z1) {
		col.r = 120; col.g = 120; col.b = 120;
	} else if (val > z2) {
		col.r = 120; col.g = 120; col.b = 120;
	} else {
		x = 1.0 + 1000.0 * ((val-z1)/(z2-z1));
		x = log10(x);
		h = 359.9 * x/3.0;
		s = 0.9;
		v = 0.8;
		col = HSVtoRGB(h,s,v);
		
		if (neg == 'n'){
			col.r = 255 - col.r;
			col.g = 255 - col.g;
			col.b = 255 - col.b;
		}
		
	}	
	return(col);
}

static colour rainbow(float x) {		// gives a rainbow. x between 0 and 1
										// rainbow = 6 segments in 3D rgb space
										// I start in blue corner, going towards cyan and green
	x = x * 6.0;
	colour col;
	
	if ((x > 0.0) && (x <= 1.0)) {
		col.r = (unsigned char) floor(255.9 * (1.0 - x));
		col.g = 0;
		col.b = 255;
	} else if ((x > 1.0) && (x <= 2.0)) {
		col.r = 0;
		col.g = (unsigned char) floor(255.9 * (x - 1.0));
		col.b = 255;
	} else if ((x > 2.0) && (x <= 3.0)) {
		col.r = 0;
		col.g = 255;
		col.b = (unsigned char) floor(255.9 * (3.0 - x));
	} else if ((x > 3.0) && (x <= 4.0)) {
		col.r = (unsigned char) floor(255.9 * (x - 3.0));
		col.g = 255;
		col.b = 0;
	} else if ((x > 4.0) && (x <= 5.0)) {
		col.r = 255;
		col.g = (unsigned char) floor(255.9 * (5.0 - x));
		col.b = 0;
	} else {
		col.r = 255;
		col.g = 0;
		col.b = (unsigned char) floor(255.9 * (x - 5.0));
	}
	
	return(col);
}

static colour lincol1(float val, float z1, float z2, char neg) {
	colour col;
	float x;
	if (val < z1) {
		col.r = 0; col.g = 0; col.b = 0;
	} else if (val > z2) {
		col.r = 255; col.g = 255; col.b = 255;
	} else {
		//x = 1.0 + 1000.0 * ((val-z1)/(z2-z1));
		//x = log10(x);	// x/3.0 est donc entre 0 et 1
		
		x = ((val-z1)/(z2-z1));
		
		x = 0.1666 + x * 0.6666;  // ira du bleu au rouge (voir figure)
		
		col = rainbow(x);
		
		if (neg == 'n'){
			col.r = 255 - col.r;
			col.g = 255 - col.g;
			col.b = 255 - col.b;
		}
		
	}	
	return(col);
}


static void GetMinAndMax(const float *inarr, unsigned int n, float *min, float *max)
{
	*min = *inarr;
	*max = *inarr;
	for (unsigned int i=0; i<n; i++) {
		if (*inarr < *min) {
			*min = *inarr;
		}
		if (*inarr > *max) {
			*max = *inarr;
		}
		inarr++;
	}
}
/* This function is the quick_select routine *
 * arr[] is the image
 * n is the number of pixels in the image
 * min is the min pixel value
 * max is the max pixel value
 * a1 is the median
 * a2 is the 95%ile
 */

static void GetMedianA1A2(const float *inarr, unsigned int n, float *a1, float *a2) 
{
	int low, high; 
	int median; 
	int middle, ll, hh; 
	
	short odd = (n & 1);
	
	float *arr = (float *)malloc(n*sizeof(float));
	memcpy(arr, inarr, n*sizeof(float));
	
	low = 0; high = n-1; median = (low + high) / 2; 
	while (1) { 
		
		if (high <= low) {/* One element only */ 
			*a1 = arr[median];
			*a2 = arr[(int)floor((n-1) * 0.95)];
			free(arr);
			return; 
		}
		if (high == low + 1) { /* Two elements only */ 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]) ; 
			if (odd) {
				*a1 = arr[median];
				*a2 = arr[(int)floor((n-1) * 0.95)];
				free(arr);
				return; 
			} else {
				*a1 = (arr[low] + arr[high]) / 2.0;
				*a2 = arr[(int)floor((n-1) * 0.95)];
				free(arr);
				return; 
			}
		} 
		
		/* Find median of low, middle and high items; swap into position low */ 
		middle = (low + high) / 2; 
		if (arr[middle] > arr[high]) 
			ELEM_SWAP(arr[middle], arr[high]) 
			if (arr[low] > arr[high]) 
				ELEM_SWAP(arr[low], arr[high]) 
				if (arr[middle] > arr[low]) 
					ELEM_SWAP(arr[middle], arr[low]) 
					
				/* Swap low item (now in position middle) into position (low+1) */ 
					ELEM_SWAP(arr[middle], arr[low+1]) ; 
		
		/* Nibble from each end towards middle, swapping items when stuck */ 
		ll = low + 1; 
		hh = high; 
		for (;;) { 
			do ll++; while (arr[low] > arr[ll]) ; 
			do hh--; while (arr[hh] > arr[low]) ; 
			
			if (hh < ll)
				break;
			ELEM_SWAP(arr[ll], arr[hh])
		} 
		
		/* Swap middle item (in position low) back into correct position */ 
		ELEM_SWAP(arr[low], arr[hh]) 
		
		/* Re-set active partition */ 
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1; 
	} 
} 

/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
	fprintf(stderr,
			"usage: wirpreview  [<optional arguments>] <fitsfilename list>\n"
			"  Optional arguments:\n"
			"  -o, --output=<output filename>\n"
			"  -D, --dir=<output directory>\n"
			"  -s, --src=<source directory>\n"
			"  -r, --ratio=<integer from 1 to 64> (default 26) -- smount by which to reduce image size\n"
			"  -a, --args=\"<list of convert jpg/png conversion args>\"\n"
			"  -m, --margs=\"<list of ffmpeg movie conversion args>\"\n"
			"  -k, --kind= jpg | png | mpg\n"
			"  -f, --filename, Annotate the filename to the image\n"
			"  -p, --pi,       Annotate the PI name to the image\n"
			"  -c, --coords,   Annotate the coordinates to the image\n"
			"  -n, --name,     Annotate the object name to the image\n"
			"  -h, --help,     Display this help message\n"
			"  -v, --verbose,  Turn on verbose message output.\n"
			"  -l, --list,     Turn on listing of filenames processed in to movie.\n"
			"  -d, --debug,    Turn on debug message output and retain intermediate temp files.\n"
			"  \n"
			"  wirpreview converts a list WIRCam FITS images to jpg or png image thumbnails.\n"
			"  The ratio arguments determines the amount by which the image is reduced (default 16).\n"
			"  If the output filename is not given, then it is constructed from the input filename list.\n"
			"  The --args= option sends all arguments to \"convert\". The most useful option is -equalize which\n"
			"  stretches the image by histogram equalization. The default output kind is determined by the\n"
			"  output filename file extension if known, or as given by the --kind argument. (default is jpg)\n"
			"  Examples:\n"
			"  wirpreview 1244561o.fits\n"
			"  wirpreview --kind=png 1244561o.fits\n"
			"  wirpreview --ouput=/data/wircam/iiwi2/09BQ04-Sep22.mpg /data/lokahi/wircam/09BQ04-Sep22/*o.fits\n"
			"  wirpreview --filename --pi --name --ouput=night.avi *o.fits\n"
			"  wirpreview --output=night.m4v *o.fits\n"
			"  wirpreview --output=night.mpg margs=\"-r 4\" *o.fits\n"
			"  wirpreview --dir=/data/uhane4/1244/ 1244560o.fits 1244561o.fits 1244562o.fits\n"
			"  wirpreview --args=\"-equalize\" 1244561o.fits\n"
			"  wirpreview --output=./out.png 1244561o.fits\n"
			"  wirpreview 1244561o.fits[2]\n"
			"  wirpreview 1244561o.fits[*]\n");
}

//-------------------- main code --------------------------
static char tmpfilename[MAX_FILEPATH_SIZE];
static char mpgfilename[MAX_FILEPATH_SIZE];
static char fitsfilename[MAX_FILEPATH_SIZE];
static char outfilename[MAX_FILEPATH_SIZE];
static char dirname[MAX_FILEPATH_SIZE];
static char srcdirname[MAX_FILEPATH_SIZE];
static char odometer[512];
static char userargs[2048];
static char usermargs[2048];
static char args[4096];
static char piname[256];
static char objname[256];
static char RA[256];
static char DEC[256];

int main(int argc, char *const argv[])
{
	unsigned short ext;
	char *filebase = NULL;
	char *tmpbase = NULL;
	
	int verbose = 0, debug = 0, plot=0;
	
	int list = 0;
	int slice = 1;
	int slices = 1;
	int ratio = 16;
	short addpiname = FALSE;
	short addcoords = FALSE;
	short addobjectname = FALSE;
	short addfilename = FALSE;
	
	enum outputkinds outputkind = E_JPG;
	
	float *image = NULL;
	float *constructed = NULL;
	float *smallimage = NULL;
	
	// PNG related
	FILE *outfptr = NULL;		// to write the png
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep row_ptr;	    
	float *fitsrow;				// a line of pixels
	unsigned char *rgb;			// a png is a table of bytes
	char *ztrans = (char *)"l";			// linear by default
	const int borderwidth = 2;
	
	colour col;
	
	float min = 0.0;
	float max = 65535.0;
	
	float a1, a2;				// automatic cutoffs
	
	float z=0, z1=0, z2=0, logz1=0, logz2=0, logval=0;
	char neg = 'p';				// 'n' means black on white	
	char *lower = (char *)"0", *upper = (char *)"65535";
	
	fitsfile *fptr = NULL, *tptr = NULL;
	
	int status = 0;
	int framecount = 0;
	long opixel[2];
	int hdutype = ANY_HDU;
	int naxis = 0;
	long naxes[3];  
	long nconstructedaxes[2];
	long npixels;
	long nsmallpixels;
	long nconstructedpixels;
	
	int gap;
	int smallsize; 
	int constructedsize;
	int pixelindex = 0;
	char hasextension = FALSE;
	
	args[0] = 0;
	userargs[0] = 0;
	usermargs[0] = 0;
	RA[0] = 0;
	DEC[0] = 0;
	objname[0] = 0;
	piname[0] = 0;
	outfilename[0] = 0;
	mpgfilename[0] = 0;
	srcdirname[0] = 0;
	odometer[0] = 0;
	
	int opt;
	
	/* Go through and handle the filename and path information */
	struct option longopts[] = {
		{"output",  required_argument, NULL, 'o'},
		{"args",    optional_argument, NULL, 'a'},
		{"margs",   optional_argument, NULL, 'm'},
		{"dir",     optional_argument, NULL, 'D'},
		{"kind",    optional_argument, NULL, 'k'},
		{"ratio",   optional_argument, NULL, 'r'},
		{"src",     optional_argument, NULL, 's'},
		{"pi",      no_argument, NULL, 'i'},
		{"coords",  no_argument, NULL, 'c'},
		{"name",    no_argument, NULL, 'n'},
		{"filename",no_argument, NULL, 'f'},
		{"ztrans",	optional_argument, NULL, 'z'},
		{"lower",	optional_argument, NULL, 'l'},
		{"upper",	optional_argument, NULL, 'u'},
		{"list",    optional_argument, NULL, 't'},
		
		{"plot",    optional_argument, NULL, 'p'},
		{"verbose", optional_argument, NULL, 'v'},
		{"debug",   optional_argument, NULL, 'd'},
		{"help",    optional_argument, NULL, 'h'},
		{0,0,0,0}
	};
	
	strcpy(dirname, ".");
	
	while ((opt = getopt_long(argc, argv, "o:a:m:D:k:r:s:z:l:u:pcfntv::d::p::h::", longopts, NULL))  != -1) {
		switch(opt) {
			case 'o':
				strncpy(outfilename, optarg, sizeof(outfilename));
				break;
			case 'D':
				strncpy(dirname, optarg, sizeof(dirname));
				if (strlen(outfilename) > 0 && strchr(outfilename, '/') != NULL) {	// --dir overrides dirpath in outfilename
					strncpy(tmpfilename, strrchr(outfilename, '/')+1, sizeof(tmpfilename));
					strncpy(outfilename, dirname, sizeof(outfilename));
					strncat(outfilename, tmpfilename, sizeof(outfilename));
				}
				break;
			case 'a':
				strncpy(userargs, optarg, sizeof(userargs));
				break;
			case 'm':
				strncpy(usermargs, optarg, sizeof(usermargs));
				break;
			case 'k':
				if (strstr(optarg, "jpg"))
					outputkind = E_JPG;
				else if (strstr(optarg, "png"))
					outputkind = E_PNG;
				else if (strstr(optarg, "mpg"))
					outputkind = E_MPG;
				else
					fprintf(stderr, "Invalid output kind %s, ignored, jpg used.\n", optarg);
				break;
			case 'r':
				ratio = (abs(atoi(optarg))>>1)<<1; //  (x & (y − 1))    ceil(ilog2(value))
				if (ratio <= 0)
					ratio = 1;
				break;
			case 's':
				strncpy(srcdirname, optarg, sizeof(srcdirname));
				break;
			case 'i':
				addpiname = TRUE;
				addfilename = TRUE;
				break;
			case 'c':
				addcoords = TRUE;
				addfilename = TRUE;
				break;
			case 'n':
				addobjectname = TRUE;
				addfilename = TRUE;
				break;
			case 'f':
				addfilename = TRUE;
				break;
			case 't':
				list = 1;
				break;
			case 'z':
				ztrans = optarg;
				if (!strcmp(optarg, "l")) {				// LINEAR SCALE
					if (verbose) printf("wirpreview: using linear scale\n");
				} else if (!strcmp(optarg, "lc")) {		// LINEAR SCALE WITH COLOURED CUTOFFS
					if (verbose) printf("wirpreview: using linear scale with coloured cutoffs\n");
				} else if (!strcmp(optarg, "lm")) {		// LINEAR SCALE WITH CENTRAL MARK
					if (verbose) printf("wirpreview: using linear scale with crosshair\n");	
				} else if (!strcmp(optarg, "e")) {		// LOGARITHMIC SCALE "CENTERED ON 0.0"
					if (verbose) printf("wirpreview: using log scale centered on 0.0\n");
				} else if (!strcmp(optarg, "f")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY
					if (verbose) printf("wirpreview: using IRAF log scale\n");
				} else if (!strcmp(optarg, "fc")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY WITH COLOURED CUTOFFS
					if (verbose) printf("wirpreview: using IRAF log scale with coloured cutoffs\n");
				} else if (!strcmp(optarg, "p")) {		// LIN FROM Z1 TO -Z1, THEN LOG
					if (verbose) printf("wirpreview: using lin in the [z1,-z1] range, then log\n");
				} else if (!strcmp(optarg, "c")) {		// LOG WITH COLOUR CONTOURS
					if (verbose) printf("wirpreview: using log with colour contours\n");
					if (neg != 'p') { 
						if (verbose) fprintf(stderr, "wirpreview: please use positive output for ztrans=c\n");
						return(EXIT_FAILURE);
					}
				} else if (!strcmp(optarg, "r")) {		// PSEUDOCOLOUR LOG HSV
					if (verbose) printf("wirpreview: using pseudocolour log hsv\n");
				} else if (!strcmp(optarg, "s")) {		// PSEUDOCOLOUR LINCOL1
					if (verbose) printf("wirpreview: using rainbow\n");
				} else if (!strcmp(optarg, "h")) {		// HISTOGRAM EQUALIZATION
					if (verbose) printf("wirpreview: using histogram equalization\n");
				} else if (!strcmp(optarg, "u")) {		// HISTOGRAM EQUALIZATION
					if (verbose) printf("wirpreview: using upena equalization\n");
				} else {
					fprintf(stderr, "wirpreview: unknown scale\n");
					return(EXIT_FAILURE);
				}
				break;
			case 'l':
				z1 = 1;
				lower = optarg;
				break;
			case 'u':
				upper = optarg;
				z2 = 1;
				break;
				
			case 'v':
				verbose = 1;
				break;
			case 'p':
				plot = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	if (optind == argc) {
		fprintf(stderr, "error: wirpreview: (%s:%s:%d) Please specify a fits input filename\n", __FILE__, __func__, __LINE__);
		printUsageSyntax(argv[0]);
		exit(EXIT_FAILURE);
	}
	if (argc < 1 ) {
		printUsageSyntax(argv[0]);    
		exit(EXIT_FAILURE);
	}
	if (strlen(outfilename) > 0) {
		if (strstr(outfilename, ".jpg"))
			outputkind = E_JPG;
		else if (strstr(outfilename, ".png"))
			outputkind = E_PNG;
		else if (strstr(outfilename, ".mpg"))
			outputkind = E_MPG;
		else if (strstr(outfilename, ".mov"))
			outputkind = E_MPG;
		else if (strstr(outfilename, ".m4v"))
			outputkind = E_MPG;
		else if (strstr(outfilename, ".avi"))
			outputkind = E_MPG;
	}
	
	if (outputkind == E_MPG) {
		if (strlen(outfilename) > 0) {
			strncpy(mpgfilename, outfilename, sizeof(mpgfilename));
			outfilename[0] = 0;
			if (!strcmp(dirname, ".")) {
				strcpy(dirname, "/tmp/");
			}
		}
		addfilename = TRUE;
	}
	/* Process a filename list. */
	
	try {
		while (optind < argc) {
			
			strncpy(fitsfilename, argv[optind++], sizeof(fitsfilename));
			
			if (outputkind == E_MPG) {
				if (strlen(mpgfilename) == 0) {
					strncpy(mpgfilename, fitsfilename, sizeof(mpgfilename));
					if ( ( tmpbase = strchr(mpgfilename, '[') ) ) {
						*tmpbase = 0;
					}
					if (strchr(mpgfilename, '/') != NULL) {	// --dir overrides dirpath in outfilename
						strncpy(tmpfilename, strrchr(mpgfilename, '/')+1, sizeof(tmpfilename));
						strncpy(mpgfilename, dirname, sizeof(mpgfilename));
						strncat(mpgfilename, tmpfilename, sizeof(mpgfilename));
					}
					mpgfilename[strlen(mpgfilename)-5] = 0;
					strncat(mpgfilename, ".mpg", sizeof(mpgfilename));
				}
			}
			
			//---------------------------------------------------------------
			filebase = strrchr(fitsfilename, '/');
			if ( filebase == NULL ) {
				filebase = fitsfilename;
			} else {
				filebase++;
			}
			
			if ( ( tmpbase = strchr(filebase, '[') ) ) {
				*tmpbase = 0;
				hasextension = TRUE;
			}
			
			if (strlen(odometer) == 0 ) {
				strncpy(odometer, filebase, sizeof(odometer));
				odometer[strlen(odometer)-5] = 0;
			}
			
			//-- Open the fits source file      
			if ( fits_open_image(&fptr, fitsfilename, READONLY, &status) ) {
				fits_report_error(stderr, status); /* print error report and continue */
			} else {			
				
				//-- Get the dimensions      
				hdutype = ANY_HDU;
				if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
					throw operaException("wirpreview: cfitsio error Moving to HDU ", status, __FILE__, __FUNCTION__, __LINE__);	
				if (fits_get_img_dim(fptr, &naxis, &status))	// read current ext img dimensions
					throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
				if (fits_get_img_size(fptr, naxis, naxes, &status)) // read size of each dimension
					throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
				
				if (addpiname) {
					fits_read_key(fptr, TSTRING, "PI_NAME", piname, NULL, &status);
				}
				if (addcoords) {
					fits_read_key(fptr, TSTRING, "RA", RA, NULL, &status);
					fits_read_key(fptr, TSTRING, "DEC", DEC, NULL, &status);
				}
				if (addobjectname) {
					fits_read_key(fptr, TSTRING, "OBJNAME", objname, NULL, &status);
				}
				
				if (hasextension) {
					*tmpbase = '[';
				}
				if ( strstr(filebase, "[*]") ) {
					if (naxis == 3)
						slices = naxes[2];
					if ( ( tmpbase = strchr(filebase, '[') ) ) {
						*tmpbase = 0;
					}
				} else {
					if ( ( tmpbase = strchr(filebase, '[') ) ) {
						slices = atoi(tmpbase+1);
						if (naxis == 3) {
							if (slices > naxes[2]) {
								slices = naxes[2];
							}
						} else {
							slices = 1;
						}
						*tmpbase = 0;
					}
				}
				
				while (slice <= slices) {
					
					if (!strlen(outfilename)) {
						strncpy(outfilename, dirname, sizeof(outfilename));
						strncat(outfilename, "/", sizeof(outfilename));
						strncat(outfilename, filebase, sizeof(outfilename));
						outfilename[strlen(outfilename)-5] = 0;
						if (slice > 1) {
							char tmp[24];
							sprintf(tmp, "-slice%.2d", slice);
							strncat(outfilename, tmp, sizeof(outfilename));
						}
						if (outputkind == E_JPG)
							strncat(outfilename, ".jpg", sizeof(outfilename));
						else if (outputkind == E_PNG)
							strncat(outfilename, ".png", sizeof(outfilename));
						else if (outputkind == E_MPG)
							strncat(outfilename, ".jpg", sizeof(outfilename));
					}
					
					if (verbose) fprintf(stderr, "wirpreview: processing %s->%s verbose=%d slice=%d ratio=%d\n", fitsfilename, outfilename, verbose, slice, ratio);
					
					strcpy(tmpfilename, "/tmp/");
					strncat(tmpfilename, filebase, sizeof(tmpfilename));
					remove(tmpfilename);
					
					if (fits_create_file(&tptr, tmpfilename, &status)) /* create new FITS file */         
						throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
					
					npixels = naxes[0] * naxes[1];      /* number of pixels in the image */	
					
					gap = 150 / ratio; 
					smallsize = naxes[0] / ratio;
					constructedsize = (2 * smallsize) + gap;
					nconstructedaxes[0] = nconstructedaxes[1] = constructedsize;
					nsmallpixels = smallsize * smallsize; 
					nconstructedpixels = nconstructedaxes[0] * nconstructedaxes[1];
					
					//-----------------------------------------------------------------------  
					
					if (fits_create_img(tptr, FLOAT_IMG, 2, nconstructedaxes, &status)) 
						throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
					//-- Allocate storage  
					image = (float *)malloc(npixels*sizeof(float));      
					constructed = (float *)malloc(nconstructedpixels*sizeof(float));
					bzero(constructed, nconstructedpixels*sizeof(float));
					smallimage = (float *)malloc(nsmallpixels*sizeof(float));
					
					if (verbose) fprintf(stderr, "wirpreview: small: %dx%d constructed: %dx%d\n", smallsize, smallsize, smallsize+gap, smallsize+gap);
					
					//-- Loop over extensions inside each input file  
					for (ext=0; ext<WIRCAM_EXTENSIONS; ext++) { 
						int xx, yy;
						float offmed;
						float offdev;
						long fpixel[3], lpixel[3], inc[3];
						int anynull;
						float nullval = 0.0;   
						inc[0] = 1; 
						inc[1] = 1;
						inc[2] = 1;
						fpixel[0] = 1; 
						fpixel[1] = 1;
						fpixel[2] = slice;
						lpixel[0] = naxes[0]; 
						lpixel[1] = naxes[1];
						lpixel[2] = slice;
						opixel[0] = 1; 
						opixel[1] = 1;
						
						/* select extension */
						if ( fits_movabs_hdu(fptr, ext+2, &hdutype, &status) )
							throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
						
						/* load image slice */
						if (fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &nullval, image, &anynull, &status))
							throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
						
						smallimage = shrink(image, naxes[0], ratio, smallimage);	
						
						offmed = operaArrayMedian(nsmallpixels, (float *)smallimage);
						offdev = operaArrayMedianSigma(nsmallpixels, (float *)smallimage, offmed);
						
						if (verbose) fprintf(stdout, "wirpreview: %s[%d]: image median level: %4.2f deviation: %4.2f\n",fitsfilename,ext+1,offmed,offdev);
						//- Rescale the image and compute stats of the image background
						for (yy=0; yy<smallsize; yy++) {
							for (xx=0; xx<smallsize; xx++) {
								if ( (index(smallimage,xx,yy,smallsize) > (offmed+15.0*offdev)) || (index(smallimage,xx,yy,smallsize) < (offmed-15.0*offdev)) ) {
									index(smallimage,xx,yy,smallsize) = 0.0;
								}
								if (offdev != 0.0) {
									float z1 = -15.0;
									float z2 = +15.0;
									float m = 255.0/((z2-z1)*offdev);
									float devd = m*(index(smallimage,xx,yy,smallsize)-(offmed+z1*offdev));
									if (index(smallimage,xx,yy,smallsize) != 0.0) {
										index(smallimage,xx,yy,smallsize) = devd; // this highlights dead amps
									}
								}
							}
						}
						if (ext == 0) {
							for (yy=0; yy<smallsize; yy++) {
								for (xx=0; xx<smallsize; xx++) {
									index(constructed,xx+smallsize+gap,yy+smallsize+gap,constructedsize) = index(smallimage,xx,yy,smallsize);
								}
							}
						} else if (ext == 1) {
							for (yy=0; yy<smallsize; yy++) {
								for (xx=0; xx<smallsize; xx++) {
									index(constructed,xx+smallsize+gap,yy,constructedsize) = index(smallimage,xx,yy,smallsize);
								}
							}
						} else if (ext == 2) {
							for (yy=0; yy<smallsize; yy++) {
								for (xx=0; xx<smallsize; xx++) {
									index(constructed,xx,yy,constructedsize) = index(smallimage,xx,yy,smallsize);
								}
							}
						} else if (ext == 3) {
							for (yy=0; yy<smallsize; yy++) {
								for (xx=0; xx<smallsize; xx++) {
									index(constructed,xx,yy+smallsize+gap,constructedsize) = index(smallimage,xx,yy,smallsize);
								}
							}
						}
					} // end of loop over extensions   
					
					if (outputkind == E_PNG) {
						
						GetMinAndMax(constructed, nconstructedpixels, &min, &max);
						
						if (z2 == 1) { 
							GetMedianA1A2(constructed, nconstructedpixels, &a1, &a2);
							if (!strcmp(lower, "m")) {
								z1 = min;
								if (verbose) printf("wirpreview: z1 = min = %f\n", z1);
							} else if (!strcmp(lower, "a")) {
								z1 = a1;
								if (verbose) printf("wirpreview: z1 = auto = %f\n", z1);
								if (verbose) printf("wirpreview: (min value is %f)\n", min);
							} else {
								z1 = atof(lower);
								if (verbose) printf("z1 = provided = %f\n", z1);
							}
							if (!strcmp(upper, "m")) {
								z2 = max;
								if (verbose) printf("wirpreview: z2 = max = %f\n", z2);
							} else if (!strcmp(upper, "a")) {
								z2 = a2;
								if (verbose) printf("wirpreview: z2 = auto = %f\n", z2);
								if (verbose) printf("wirpreview: (max value is %f)\n", max);
							} else {
								z2 = atof(upper);
								if (verbose) printf("wirpreview: z2 = provided = %f\n", z2);
							}
						} else {
							GetMedianA1A2(constructed, nconstructedpixels, &a1, &a2);
							z1 = min;
							z2 = max;
						}	
						if (z1 > z2) {
							if (verbose) fprintf(stderr, "wirpreview: z1 > z2 : switching...\n");
							z = z1;
							z1 = z2;
							z2 = z;
						}
						if (z1 == 0.0 && z2 == 0.0) {
							throw operaException("wirpreview: both cutoffs = 0, aborting... ", 0, __FILE__, __FUNCTION__, __LINE__);	
						}
						
						if (!strcmp(ztrans, "e")) {		// LOGARITHMIC SCALE "CENTERED ON 0.0"
							logz1=log10(fabs(z1)+1);
							logz2=log10(fabs(z2)+1);
							if (z1 < 0) {
								logz1 = -logz1;
							}
							if (z2 < 0) {
								logz2 = -logz2;
							}
							if (verbose) printf("wirpreview: logz1=%f logz2=%f \n", logz1, logz2);
						} else if (!strcmp(ztrans, "p")) {		// LIN FROM Z1 TO -Z1, THEN LOG
							if (verbose) printf("wirpreview: using lin in the [z1,-z1] range, then log\n");
							if (z1 >= 0.0) { 
								throw operaException("wirpreview: z1 must be negative...", 0, __FILE__, __FUNCTION__, __LINE__);	
							}
							if (z2 <= 0.0) { 
								throw operaException("wirpreview: z1 must be positive...", 0, __FILE__, __FUNCTION__, __LINE__);	
							}
						}  else if (!strcmp(ztrans, "h"))  {
							histogramEqualize(constructed, constructedsize, constructedsize, z1, z2);
							ztrans = (char *)"l";
						}
						if (verbose)  printf("wirpreview: z1=%4.2f, z2=%4.2f, min=%4.2f, max=%4.2f\n", z1, z2, min, max);
						
						png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL);
						if (!png_ptr) {
							throw operaException("wirpreview: Could not create png write structure...", 0, __FILE__, __FUNCTION__, __LINE__);	
						}
						info_ptr = png_create_info_struct(png_ptr);
						if (!info_ptr)
						{
							png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
							throw operaException("wirpreview: Could not create png info structure...", 0, __FILE__, __FUNCTION__, __LINE__);	
						}
						
						outfptr = fopen(outfilename, "w");
						png_init_io(png_ptr, outfptr);
						
						if (verbose) fprintf(stderr, "wirpreview: Constructed size = %d x %d\n", constructedsize, constructedsize);
						png_set_IHDR(png_ptr, info_ptr, constructedsize, constructedsize, 8, PNG_COLOR_TYPE_RGB,
									 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
									 PNG_FILTER_TYPE_DEFAULT);
						png_write_info(png_ptr, info_ptr);
						
						// allocate vector for the lines of the fits image
						fitsrow = (float *) malloc(constructedsize * sizeof(float)); 	// a single line
						if (fitsrow == NULL) {
							throw operaException("wirpreview: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
						}
						rgb = (unsigned char *) malloc(constructedsize * 3 * sizeof(unsigned char));
						if (rgb == NULL) {
							throw operaException("wirpreview: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
						}
						pixelindex = nconstructedpixels;
						for (int j = 0; j < constructedsize; j++) {
							
							pixelindex -= constructedsize;
							memcpy((char *)fitsrow, (char *)(&constructed[pixelindex]), constructedsize*sizeof(float));
							
							for (int i = 0; i < constructedsize; i++) { 
								
								col.r = 255; col.g = 255; col.b = 255;
								
								if (i < borderwidth || i > (constructedsize-borderwidth) || j < borderwidth || j > (constructedsize-borderwidth)) {
								} else if (!strcmp(ztrans, "l")) {				// LINEAR SCALE
									col = lingray(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "lc")) {				// LINEAR SCALE WITH COLOURED CUTOFFS
									col = lingraycutcol(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "lm")) {				// LINEAR SCALE WITH CROSSHAIR
									int centi = (int)floor(constructedsize/2);
									int centj = (int)floor(constructedsize/2);
									col = lingray(fitsrow[i], z1, z2, neg);
									if (((abs(j-centj) == 10) && (i == centi)) || ((abs(i-centi) == 10) && (j == centj))) {
										col.r = 0; col.g = 255; col.b = 0;
									}
									if (((abs(j-centj) == 15) && (i == centi)) || ((abs(i-centi) == 15) && (j == centj))) {
										col.r = 0; col.g = 255; col.b = 0;
									}
								} else if (!strcmp(ztrans, "e")) {				// LOGARITHMIC SCALE "CENTERED ON 0.0"
									logval = log10(fabs(fitsrow[i])+1);
									if (fitsrow[i] < 0) {
										logval = -logval;
									}
									col = lingray(logval, logz1, logz2, neg);
									
								} else if (!strcmp(ztrans, "f")) {				// LOGARITHMIC SCALE LIKE IRAF DISPLAY
									col = irafloggray(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "fc")) {				// LOGARITHMIC SCALE LIKE IRAF DISPLAY WITH COLOURED CUTOFFS
									col = irafloggraycutcol(fitsrow[i], z1, z2, neg);
									
								} else if (!strcmp(ztrans, "p")) {				// LIN FROM Z1 TO -Z1, THEN LOG
									col = plinlog(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "c")) {				// LOG WITH COLOUR CONTOURS
									col = logcontours(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "r")) {				// PSEUDOCOLOUR LOG HSV
									col = loghsv(fitsrow[i], z1, z2, neg);
								} else if (!strcmp(ztrans, "s")) {				// PSEUDOCOLOUR LINCOL1
									col = lincol1(fitsrow[i], z1, z2, neg);
								}
								
								rgb[0 + 3*i] = col.r;
								rgb[1 + 3*i] = col.g;
								rgb[2 + 3*i] = col.b;			
							}
							row_ptr = rgb;
							png_write_rows(png_ptr, &row_ptr, 1);
						}
						
						png_write_end(png_ptr, info_ptr);
						png_destroy_write_struct(&png_ptr, &info_ptr);
						
						free(rgb);
						free(fitsrow);
						
						fclose (outfptr);
					} else {	//- Do conversion to jpg
						/* Write data into the tmp file*/
						if (verbose) fprintf(stderr, "wirpreview: Writing data to %s\n", tmpfilename);	
						if (fits_write_pix(tptr, TFLOAT, opixel, nconstructedpixels, constructed, &status))
							throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
						
						free(constructed); 
						free(image);
						free(smallimage);
						
						if ( fits_close_file(tptr, &status) )                /* close the temp fits file */
							throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
						
						//- Annotate image
						strncpy(args, userargs, sizeof(userargs));
						if (addfilename) {
							int pointsize = 10;
							if (ratio <= 8) pointsize = 16;
							if (ratio <= 4) pointsize = 20;
							if (ratio <= 2) pointsize = 32;
							snprintf(&args[strlen(args)], sizeof(args), " -gravity east -pointsize %d -fill white -annotate 0 '%s'", pointsize, filebase);
							if (addpiname || addcoords || addobjectname) {
								snprintf(&args[strlen(args)], sizeof(args), " -gravity west -fill white -annotate 0 '%s %s %s %s'", piname, objname, RA, DEC);
							}
						}
						//- Do conversion to jpg
						if (verbose) fprintf(stderr, "wirpreview: convert %s %s %s\n", args, tmpfilename, outfilename);	
						systemf("convert %s %s %s", args, tmpfilename, outfilename);
					}
					//- Rename to sequential filenames for ffmpeg
					if (outputkind == E_MPG) {
						systemf("mkdir -p /tmp/%s/", odometer);
						systemf("mv %s /tmp/%s/%04d.jpg", outfilename, odometer, framecount);
					}
					if (!debug) remove(tmpfilename);   /* remove the temp */
					framecount++;
					slice++;
					outfilename[0] = 0;
					if (outputkind == E_MPG && verbose == 0) {
						if (list) {
							printf("%s\n", filebase);
						} else {
							printf(".");
							fflush(stdout);
						}
					}
				} // while <slice < slices
				if ( fits_close_file(fptr, &status) )                /* close the source fits file */
					throw operaException("wirpreview: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
				slice = 1;
			} // if fits_open_image
		} // while optind < argc
		
		//- now make the movie...
		if (outputkind == E_MPG && framecount > 0) {
			printf("\n");
			if (verbose) printf("wirpreview: ffmpeg -y -an %s -i /tmp/%s/%%04d.jpg %s\n", usermargs, odometer, mpgfilename);
			systemf("ffmpeg -y -an %s -i /tmp/%s/%%04d.jpg %s &>/dev/null", usermargs, odometer, mpgfilename);
			if (!debug) {
				systemf("rm -rf /tmp/%s/ 2>/dev/null", odometer);
				systemf("rmdir /tmp/%s/ 2>/dev/null", odometer);
			}
		}
		if (verbose && framecount == 0) fprintf(stderr, "wirpreview: No images were found.\n");
	}
	catch (operaException e) {
		cerr << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirpreview: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
