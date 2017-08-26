/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFITStoPNG
 Version: 1.0
 Description: Convert a FITS image to PNG.
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

/*
 Based on:
 f2n	-	fits 2 net
 A tiny program to "convert" fits images into png files
 with options about the grayscale (or colour) ztrans
 Malte Tewes, last update : December 2008 
 */

#include <stdio.h>
#include <math.h>
#include <getopt.h>

#include "fitsio.h"
#include "png.h"

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaFITStoPNG.h"

#define MEANSHRINK

/*! \brief Create a PNG thumbnail of a FITS Image.  */
/*! \file operaFITStoPNG.c */
/*! \author Doug Teeple */

/*! 
 * operaFITStoPNG - FITS to PNG
 * \brief Create a PNG thumbnail of a FITS Image.
 * \brief Based on f2n: version 0.9.1, December 2008 by malte.tewes@epfl.ch
 * \verbatim
 * Usage :   operaFITStoPNG --negate ztrans= in.fits --output=out.png [lower= [upper=]]
 * --negate: negative (inverse video) output
 * --ztrans =
 *    l: lin ztrans  (lc : same but with coloured cutoffs)
 *    e: log centered on 0.0
 *    f: log ztrans as in IRAF's display (fc : coloured cutoffs)
 *    p: lin from z1 to -z1, then log up to z2 (my fav for dec.fits)
 *    c: log with color countours around 1.0 (my fav for sm residuals)
 *    r: pseudocolour log hsv
 *    s: somewhere over the rainbow (exprimental)
 *    h: histogram equalization
 * --lower= --upper=
 *    the lower and upper cutoffs to use. Give either none or both.
 *    (not given) : min and max cutoffs (same as m m, see below)
 *    -1.2e-5 : hard cutoff (for example)
 *    m : extremal (min or max) cutoff
 *    a : auto cutoff
 *    You can also mix them.
 *    Examples : \"--lower=-10.0 --upper=a\"
 *    Examples : \"--lower=a --upper=m\"
 *    Examples : \"--lower=a --upper=a\"
 *    Examples : \"--lower=0 --upper=1\"
 * \endverbatim
 * \note   You can use file.fits[a:b,c:d] to select a region of the input file.
 * \note   (In case of trouble use quotes : \"file.fits[a:b,c:d]\")
 * \return EXIT_STATUS
 * \ingroup tools
 */
int main(int argc, char *argv[])
{
	
	char *infilename = NULL, *outfilename = NULL;
	char *lower="0", *upper="65535";
	
	fitsfile *infptr = NULL;	// cfitsio input
	FILE *outfptr = NULL;		// to write the png
	float *image = NULL;		// the whole image
	float *thumbnail = NULL;	// the thumbnail
	int anynull;
	float nullval = 0.0;   
	int status = 0;				// CFITSIO status value MUST be initialized to zero..
	unsigned int imagesizex = 0;
	unsigned int imagesizey = 0;
	unsigned int npixels = 0;
	unsigned int thumbsizex = 0;
	unsigned int thumbsizey = 0;
	unsigned int nthumbpixels = 0;
	unsigned int pixelindex = 0;
	float chipbias = 0.0;
	
	// FITS related
	long size[2] = {1,1}; 		// cfitsio size vector
	long firstpix[2] = {1,1}; 	// coordinates of first pixels, 1 based... !
	float *fitsrow;				// a line of pixels
								//int hdutype = ANY_HDU;
	
	// PNG related
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep row_ptr;	    
	unsigned char *rgb;			// a png is a table of bytes
	char *ztrans = "l";			// linear by default
	const int borderwidth = 2;
	
	int i, j;
	colour col;
	enum outputkinds outputkind = E_JPG;
	
	float min = 0.0;
	float max = 65535.0;
	
	float a1, a2;				// automatic cutoffs
	
	float z, z1=0, z2=0, logz1=0, logz2=0, logval=0;
	char neg = 'p';				// 'n' means black on white	
	unsigned int shrinksize = 100;	// % = do not reduce size, 50 would mean 50%, should be a multiple of 2
	unsigned int ratio = 1;		// calculated from size
	
	int opt;
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"output",	1, NULL, 'o'},
		{"negate",	0, NULL, 'n'},
		{"ztrans",	1, NULL, 'z'},
		{"lower",	1, NULL, 'l'},
		{"upper",	1, NULL, 'u'},
		{"kind",	1, NULL, 'k'},		// not used
		{"size",	1, NULL, 's'},
		{"ratio",	1, NULL, 'r'},
		
		{"plot",	0, NULL, 'p'},
		{"verbose",	0, NULL, 'v'},
		{"debug",	0, NULL, 'd'},
		{"trace",	0, NULL, 't'},
		{"help",	0, NULL, 'h'},
		{0,0,0,0}};
	
	if (argc < 1 ) {
		printUsageSyntax();    
		exit(EXIT_FAILURE);
	}
	
	while((opt = getopt_long(argc, argv, "o:z:l:u:s:k:r:pvndth", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'o':
				outfilename = optarg;													
				break;
			case 'n':
				neg = opt;
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
				// now set ratio to be a multiple of 2
				ratio = (abs(atoi(optarg)>>1))<<1; //  (x & (y − 1))    ceil(ilog2(value))
				if (ratio <= 0)
					ratio = 1;
				if (verbose) printf("ratio=1/%d\n", ratio);
				break;
			case 's':
				shrinksize = atoi(optarg);
				if (shrinksize < 2) {
					fprintf(stderr, "shrinksize may not be less that 2%%\n");
					return(EXIT_FAILURE);
				}
				if (shrinksize > 100) {
					fprintf(stderr, "shrinksize may not be more that 100%%, setting to 100%%...\n");
					shrinksize = 100;
				}
				// now set shrinksize to be a multiple of 2
				ratio = 100/shrinksize;
				ratio = (abs(ratio>>1))<<1; //  (x & (y − 1))    ceil(ilog2(value))
				if (ratio <= 0)
					ratio = 1;
				if (verbose) printf("shrinksize=%d%% ratio=1/%d\n", shrinksize, ratio);
				break;
			case 'z':
				ztrans = optarg;
				if (!strcmp(optarg, "l")) {				// LINEAR SCALE
					if (verbose) printf("using linear scale\n");
				} else if (!strcmp(optarg, "lc")) {		// LINEAR SCALE WITH COLOURED CUTOFFS
					if (verbose) printf("using linear scale with coloured cutoffs\n");
				} else if (!strcmp(optarg, "lm")) {		// LINEAR SCALE WITH CENTRAL MARK
					if (verbose) printf("using linear scale with crosshair\n");	
				} else if (!strcmp(optarg, "e")) {		// LOGARITHMIC SCALE "CENTERED ON 0.0"
					if (verbose) printf("using log scale centered on 0.0\n");
				} else if (!strcmp(optarg, "f")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY
					if (verbose) printf("using IRAF log scale\n");
				} else if (!strcmp(optarg, "fc")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY WITH COLOURED CUTOFFS
					if (verbose) printf("using IRAF log scale with coloured cutoffs\n");
				} else if (!strcmp(optarg, "p")) {		// LIN FROM Z1 TO -Z1, THEN LOG
					if (verbose) printf("using lin in the [z1,-z1] range, then log\n");
				} else if (!strcmp(optarg, "c")) {		// LOG WITH COLOUR CONTOURS
					if (verbose) printf("using log with colour contours\n");
					if (neg != 'p') { 
						if (verbose) fprintf(stderr, "please use positive output for ztrans=c\n");
						return(EXIT_FAILURE);
					}
				} else if (!strcmp(optarg, "r")) {		// PSEUDOCOLOUR LOG HSV
					if (verbose) printf("using pseudocolour log hsv\n");
				} else if (!strcmp(optarg, "s")) {		// PSEUDOCOLOUR LINCOL1
					if (verbose) printf("using rainbow\n");
				} else if (strstr(optarg, "h")) {			// HISTOGRAM EQUALIZATION
					if (verbose) printf("using histogram equalization\n");
				} else if (!strcmp(optarg, "u")) {		// upena colours
					if (verbose) printf("using upena equalization\n");
				} else {
					fprintf(stderr, "unknown scale\n");
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
			case 't':
				trace = 1;
				break;         
			case 'h':
				printUsageSyntax();
				exit(EXIT_SUCCESS);
				break;
			default:
				break;
		}
	}	
	
	if ((argc-optind) > 1) {
		fprintf(stderr, "Only one input FITS file may be specified, using the last one.\n");
	}
	for ( ; optind < argc; optind++) {
		infilename = argv[optind];
	}
	if (outfilename == NULL) {
		outfilename = (char *)malloc(strlen(infilename));
		strncpy(outfilename, infilename, strlen(infilename));
		outfilename[strlen(outfilename)-4] = 'p';
		outfilename[strlen(outfilename)-3] = 'n';
		outfilename[strlen(outfilename)-2] = 'g';
		outfilename[strlen(outfilename)-1] = '\0';
		if (verbose) fprintf(stdout, "outfilename=%s\n", outfilename);
	}
	// note: doesn't handle MEF
	fits_open_file(&infptr, infilename, READONLY, &status);	
	fits_get_img_size(infptr, 2, size, &status);
	
	if (status) {
		fprintf(stderr, "Could not open/read the fits file %s!\n", infilename);
		fits_report_error(stderr, status);
		return(status);
	}
	
	imagesizex = size[0];				// width
	imagesizey = size[1];				// height
	npixels = imagesizex * imagesizey;	// number of pixels in the image	
	thumbsizex = size[0] / ratio;
	thumbsizey = size[1] / ratio;
	nthumbpixels = thumbsizex * thumbsizey; 
	
	if (verbose) printf("%s : %i x %i pixels\n", infilename, imagesizex, imagesizey);
	
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL);
	if (!png_ptr) {
		fprintf(stderr, "Could not create png write structure...\n");
		return(EXIT_FAILURE);
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		fprintf(stderr, "Could not create png info structure...\n");
		png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
		return(EXIT_FAILURE);
	}
	
	outfptr = fopen(outfilename, "w");
	png_init_io(png_ptr, outfptr);
	
	png_set_IHDR(png_ptr, info_ptr, thumbsizex, thumbsizey, 8, PNG_COLOR_TYPE_RGB,
				 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
				 PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
	
	image = malloc(npixels*sizeof(float));      
	if (ratio == 1)
		thumbnail = image;
	else
		thumbnail = malloc(nthumbpixels*sizeof(float));
	
	if (verbose) fprintf(stderr, "thumbnail: %dx%d \n", thumbsizex, thumbsizey);
	
	rgb = (unsigned char *) malloc(thumbsizex * 3 * sizeof(unsigned char));
	if (rgb == NULL) {
		fprintf(stderr, "rgb : memory allocation error\n");
		return(EXIT_FAILURE);
	}
	
	/* load image  */
	if (fits_read_pix(infptr, TFLOAT, firstpix, npixels, &nullval, image, &anynull, &status)) {
		fprintf(stderr, "Can't read fits file %s\n", infilename);
		return(EXIT_FAILURE);
	}
	
	if (ratio != 1) {
		thumbnail = shrink(image, imagesizex, imagesizey, ratio, chipbias, thumbnail);	
		free(image);	// no longer needed, if ratio == 1 the thumbnail simply points to the image
	}
	// allocate vector for the lines of the fits image
	fitsrow = (float *) malloc(thumbsizex * sizeof(float)); 	// a single line
	if (fitsrow == NULL) {
		fprintf(stderr, "fitsrow : memory allocation error\n");
		return(EXIT_FAILURE);
	}
	
	GetMinAndMax(thumbnail, nthumbpixels, &min, &max);
	
	if (z2 == 1) { 
		GetMedian(thumbnail, nthumbpixels, &a1, &a2);
		if (!strcmp(lower, "m")) {
			z1 = min;
			if (verbose) printf("z1 = min = %f\n", z1);
		} else if (!strcmp(lower, "a")) {
			z1 = a1;
			if (verbose) printf("z1 = auto = %f\n", z1);
			if (verbose) printf("(min value is %f)\n", min);
		} else {
			z1 = atof(lower);
			if (verbose) printf("z1 = provided = %f\n", z1);
		}
		if (!strcmp(upper, "m")) {
			z2 = max;
			if (verbose) printf("z2 = max = %f\n", z2);
		} else if (!strcmp(upper, "a")) {
			z2 = a2;
			if (verbose) printf("z2 = auto = %f\n", z2);
			if (verbose) printf("(max value is %f)\n", max);
		} else {
			z2 = atof(upper);
			if (verbose) printf("z2 = provided = %f\n", z2);
		}
	} else {
		GetMedian(thumbnail, nthumbpixels, &a1, &a2);
		z1 = min;
		z2 = max;
	}	
	if (z1 > z2) {
		if (verbose) fprintf(stderr, "z1 > z2 : switching...\n");
		z = z1;
		z1 = z2;
		z2 = z;
	}
	if (z1 == 0.0 && z2 == 0.0) {
		if (verbose) fprintf(stderr, "both cutoffs = 0, aborting...\n");
		return(EXIT_FAILURE);
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
		if (verbose) printf("logz1=%f logz2=%f \n", logz1, logz2);
	} else if (!strcmp(ztrans, "p")) {		// LIN FROM Z1 TO -Z1, THEN LOG
		if (verbose) printf("using lin in the [z1,-z1] range, then log\n");
		if (z1 >= 0.0) { 
			if (verbose) fprintf(stderr, "z1 must be negative !\n");
			return(EXIT_FAILURE);
		}
		if (z2 <= 0.0) { 
			if (verbose) fprintf(stderr, "z2 must be positive !\n");
			return(EXIT_FAILURE);
		}
	}  else if (strstr(ztrans, "h"))  {
		histogramEqualize(thumbnail, thumbsizex, thumbsizey, z1, z2);
		GetMinAndMax(thumbnail, nthumbpixels, &z1, &z2);
		GetMedian(thumbnail, nthumbpixels, &a1, &a2);
		if (strlen(ztrans) == 1) ztrans = "l";
	}
	if (verbose)  printf("z1=%4.2f, z2=%4.2f, min=%4.2f, max=%4.2f\n", z1, z2, min, max);
	
	for (j = 0; j < thumbsizey; j++) {
		
		memcpy((char *)fitsrow, (char *)(&thumbnail[pixelindex]), thumbsizex*sizeof(float));
		pixelindex += thumbsizex;
		
		for (i = 0; i < thumbsizex; i++) {
			
			col.r = 255; col.g = 255; col.b = 255;
			
			if (i < borderwidth || i > (thumbsizex-borderwidth) || j < borderwidth || j > (thumbsizey-borderwidth)) {
			} else if (!strcmp(ztrans, "l")) {				// LINEAR SCALE
				col = lingray(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "lc")) {		// LINEAR SCALE WITH COLOURED CUTOFFS
				col = lingraycutcol(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "lm")) {		// LINEAR SCALE WITH CROSSHAIR
				unsigned int centi = floor(thumbsizex/2);
				unsigned int centj = floor(thumbsizey/2);
				col = lingray(fitsrow[i], z1, z2, neg);
				if (((abs(j-centj) == 10) && (i == centi)) || ((abs(i-centi) == 10) && (j == centj))) {
					col.r = 0; col.g = 255; col.b = 0;
				}
				if (((abs(j-centj) == 15) && (i == centi)) || ((abs(i-centi) == 15) && (j == centj))) {
					col.r = 0; col.g = 255; col.b = 0;
				}
			} else if (!strcmp(ztrans, "e")) {		// LOGARITHMIC SCALE "CENTERED ON 0.0"
				logval = log10(fabs(fitsrow[i])+1);
				if (fitsrow[i] < 0) {
					logval = -logval;
				}
				col = lingray(logval, logz1, logz2, neg);
				
			} else if (!strcmp(ztrans, "f")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY
				col = irafloggray(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "fc")) {		// LOGARITHMIC SCALE LIKE IRAF DISPLAY WITH COLOURED CUTOFFS
				col = irafloggraycutcol(fitsrow[i], z1, z2, neg);
				
			} else if (!strcmp(ztrans, "p")) {		// LIN FROM Z1 TO -Z1, THEN LOG
				col = plinlog(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "c")) {		// LOG WITH COLOUR CONTOURS
				col = logcontours(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "r")) {		// PSEUDOCOLOUR LOG HSV
				col = loghsv(fitsrow[i], z1, z2, neg);
			} else if (!strcmp(ztrans, "s")) {		// PSEUDOCOLOUR LINCOL1
				col = lincol1(fitsrow[i], z1, z2, neg);
			} else if (strstr(ztrans, "u")) {		// PSEUDOCOLOUR UPENA
				if (strstr(ztrans, "h"))
					col = lingrayupena(fitsrow[i], z1, z2, neg);
				else
					col = lingrayupenalog(fitsrow[i], z1, z2, neg);
			}
			
			rgb[0 + 3*i] = col.r;
			rgb[1 + 3*i] = col.g;
			rgb[2 + 3*i] = col.b;			
		}
		row_ptr = rgb;
		png_write_rows(png_ptr, &row_ptr, 1);
	}
	
	free(rgb);
	free(fitsrow);
	free(thumbnail);
	
	fits_close_file(infptr, &status);
	if (status) {
		fprintf(stderr, "I can't close the fits file !");
		fits_report_error(stderr, status);
	}
	
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	
	fclose (outfptr);
	
	return(EXIT_SUCCESS);
	
}

// ztrans functions

/*********************************************/
/*  HISTOGRAM EQUALIZATION of A GRAY IMAGE   */
/*********************************************/
static void histogramEqualize(float *image, unsigned int width, unsigned int height, float z1, float z2)
{	
	unsigned int pdf[histDepth];
	float cdf[histDepth];
	unsigned short *img_data = NULL;
	unsigned short *img_base;
	float *image_base = image;
	unsigned int i, npixels;
	const float scale = (float)(histDepth-1)/(z2-z1);
	
	const unsigned char verbose = 0;
	
	npixels = width * height;
	img_base = img_data = (unsigned short *)malloc(npixels*sizeof(unsigned short));
 	memset(pdf, 0, histDepth*sizeof(float));
 	memset(cdf, 0, histDepth*sizeof(float));
	if (verbose) fprintf(stdout, "\nScaled Image table\n");
	for(i=0; i<npixels; i++)
	{
		*img_data = (unsigned short)((*image-z1)*scale);	// *img_data is now scaled between 0 and 255
		if (verbose) {
			fprintf(stdout, "%d ", *img_data);
			if (i == width)
				fprintf(stdout, "\n");
		}
		if (*img_data < histDepth) {
			pdf[*img_data++]++;
		} else {
			fprintf(stderr, "Invalid access to hist: img_data=%d image=%4.2f z1=%4.2f z2=%4.2f scale=%4.2f\n", *img_data++, *image, z1, z2, scale);
		}
		image++;
	}
	if (verbose) {
		fprintf(stdout, "\nPDF table\n");
		for(i=0; i<(histDepth-1); i++)
		{
			fprintf(stdout, "%d ", pdf[i]);
				if (i == 8)
					fprintf(stdout, "\n");
		}
	}
	cdf[0] = (float)pdf[0];
	if (verbose) {
		fprintf(stdout, "\nCDF table\n");
		fprintf(stdout, "%4.2f ", cdf[0]);
	}
	for (i=1; i<histDepth; i++)	 // cdf of image
	{
		if (verbose) {
			fprintf(stdout, "%4.2f ", cdf[i]);
			if (i ==width)
				fprintf(stdout, "\n");
		}
		cdf[i] = cdf[i-1] + pdf[i];
		if (verbose) {
			fprintf(stdout, "%4.2f ", cdf[i]);
			if (i ==width)
				fprintf(stdout, "\n");
		}
	}
	img_data = img_base;
	image = image_base;
	if (verbose) fprintf(stdout, "\nHistogram Image\n");
	for(i=0; i<npixels; i++)
	{
		*image = ((float)(histDepth-1)*(cdf[*img_data] - cdf[0])/(npixels - cdf[0]))/scale;
		if (verbose) {
			fprintf(stdout, "%4.2f ", *image);
			if (i == width)
				fprintf(stdout, "\n");
		}
		img_data++;
		image++;
	}
	if (verbose) fprintf(stdout, "\n");
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

static colour lingrayupena(float val, float z1, float z2, char neg) {
	colour col;
	if (val < z1) {
		col.g = col.r = 0;
		col.b = 255;
	} else if (val > z2) {
		col.g = col.r = 255;
		col.b = 0;
	} else {
		col.g = col.r = (unsigned char) floor(255.9 * (val-z1)/(z2-z1));
		col.b = (unsigned char)(255 - col.r);
	}
	if (neg == 'n'){
		col.r = 255 - col.r;
		col.g = 255 - col.g;
		col.b = 255 - col.b;
	}
	return(col);
}

static colour lingrayupenalog(float val, float z1, float z2, char neg) {
	colour col;
	if (val < z1) {
		col.g = col.r = 0;
		col.b = 255;
	} else if (val > z2) {
		col.g = col.r = 255;
		col.b = 0;
	} else {
		col.g = col.r = (unsigned char) floor(255.9 * log10(1.0 + 1000.0 * ((val-z1)/(z2-z1)))/3.0);
		col.b = (unsigned char)(255 - col.r);
	}
	if (neg == 'n'){
		col.r = 255 - col.r;
		col.g = 255 - col.g;
		col.b = 255 - col.b;
	}
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
			col.r = (unsigned char) floor((sep + 0.9) * (val-z1)/(-z1*2.0));
		} else {
			x = 1.0 + 1000.0 * ((val+z1)/(z2+z1));
			x = log10(x);
			col.r = (unsigned char) sep + 1 + floor((255-sep-0.1) * x/3.0);
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
	i = floor( h );
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

/* This function is the quick_select routine based on the algorithm found
 * in Numerical Recipes in C.
 *
 * arr[] is the image
 * n is the number of pixels in the image
 * min is the min pixel value
 * max is the max pixel value
 * a1 is the median
 * a2 is the 95%ile
 */

static void GetMedian(const float *inarr, unsigned int n, float *a1, float *a2) 
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


#ifdef MEANSHRINK
static float *shrink(const float *image, const int dimx, const int dimy, const int ratio, const float bias, float *smallimage) {
	int y, x, locx, xx, yy;
	const int smallsizex = dimx / ratio;
	const int smallsizey = dimy / ratio;
	int locy = 0;
	float mean;
	for (y=0; y<smallsizey; y++) {			// iterate line-by-line through the thumbnail
		locx = 0;
		for (x=0; x<smallsizex; x++) {
			mean = 0.0;							// find the mean of the subwindow in image
			for (yy=locy; yy<(locy+ratio); yy++) {
				for (xx=locx; xx<(locx+ratio); xx++) {
					mean = (index(image,xx,yy,dimx)-bias+mean)/2.0;
				}
			}
			index(smallimage,x,y,smallsizex) = mean;
			locx += ratio;
		}
		locy += ratio;
	}
	return smallimage;
}
#else
static float *shrink(const float *image, const int dimx, const int dimy, const int ratio, const float bias, float *smallimage) {
	int y, x, locx;
	int locy = 0;
	const int smallsizex = dimx / ratio;
	const int smallsizey = dimy / ratio;
	const int subindex = ratio / 2;
	for (y=0; y<smallsizey; y++) {					// iterate line-by-line through the thumbnail
		locx = 0;
		for (x=0; x<smallsizex; x++) {
			index(smallimage,x,y,smallsizex) = index(image,locx+subindex/2,locy+subindex/2,dim)-bias;
			locx += ratio;
		}
		locy += ratio;
	}
	return smallimage;
}
#endif

static void printUsageSyntax() {
	printf("\n");
	printf("       operaFITStoPNG - FITS to PNG\n");
	printf("	Based on f2n: version 0.9.1, December 2008 by malte.tewes@epfl.ch\n");
	printf("\n");
	printf("Usage :   operaFITStoPNG --negate ztrans= in.fits --output=out.png [lower= [upper=]]\n");
	printf("	--negate: negative (inverse video) output\n");
	printf("	--ztrans =\n");
	printf("		l: lin ztrans  (lc : same but with coloured cutoffs)\n");
	printf("		e: log centered on 0.0\n");
	printf("		f: log ztrans as in IRAF's display (fc : coloured cutoffs)\n");
	printf("		p: lin from z1 to -z1, then log up to z2 (my fav for dec.fits)\n");
	printf("		c: log with color countours around 1.0 (my fav for sm residuals)\n");
	printf("		r: pseudocolour log hsv\n");
	printf("		s: somewhere over the rainbow (exprimental)\n");
	printf("		h: histogram equalization\n");
	printf("	--lower= --upper=\n");
	printf("		the lower and upper cutoffs to use. Give either none or both.\n");
	printf("		(not given) : min and max cutoffs (same as m m, see below)\n");
	printf("		-1.2e-5 : hard cutoff (for example)\n");
	printf("		m : extremal (min or max) cutoff\n");
	printf("		a : auto cutoff\n");
	printf("		You can also mix them.\n");
	printf("		Examples : \"--lower=-10.0 --upper=a\"\n");
	printf("		Examples : \"--lower=a --upper=m\"\n");
	printf("		Examples : \"--lower=a --upper=a\"\n");
	printf("		Examples : \"--lower=0 --upper=1\"\n");
	printf("You can use file.fits[a:b,c:d] to select a region of the input file.\n");
	printf("(In case of trouble use quotes : \"file.fits[a:b,c:d]\")\n");
	printf("\n");
	printf("-------------------------------------------------------------------------\n");
}

