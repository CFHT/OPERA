/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaPNG
 Version: 1.0
 Description: This library implements basic png plotting.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: eder@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 Based on pngwriter by Paul Blackburn
 
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
#include "libraries/operaException.h"
#include "libraries/operaPNG.h"

#include <math.h>
#include <png.h>

/* 
 * From libpngwriter 
 * by Paul Blackburn
 */

/*!
 * operaPNG
 * \author Doug Teeple
 * \brief png plotting and text library.
 * \file operaPNG.cpp
 * \ingroup libraries
 */

/*
 * Default Constructor
 */
operaPNG::operaPNG(void) :
width(250),
height(250),
backgroundcolour(65535),
compressionlevel(DEFAULT_COMPRESSION),
bitdepth(16),
colortype(2),
rowbytes(0),
filegamma(0.5),
screengamma(2.2) {
	
	filename = "out.png";
	init(width, height);	
}

/*
 *
 * Constructor for int, const char * filename
 *
 */
operaPNG::operaPNG(const char *Filename, int x, int y, int Backgroundcolour = 65535) :
compressionlevel(DEFAULT_COMPRESSION),
bitdepth(16),
colortype(2),
rowbytes(0),
filegamma(0.5),
screengamma(2.2) {
	
	width = x;
	height = y;
	filename = Filename;
	backgroundcolour = Backgroundcolour;
	init(width, height);	
}

/*
 *
 * Constructor for double, const char * filename
 *
 */
operaPNG::operaPNG(const char *Filename, int x, int y, double Backgroundcolour = 65535) :
compressionlevel(DEFAULT_COMPRESSION),
bitdepth(16),
colortype(2),
rowbytes(0),
filegamma(0.5),
screengamma(2.2) {
	
	backgroundcolour = int(Backgroundcolour*65535);
	
	if (backgroundcolour > 65535) {
		backgroundcolour = 65535;
	}
	
	if (backgroundcolour < 0) {
		backgroundcolour = 0;
	}
	width = x;
	height = y;
	filename = Filename;
	init(width, height);	
}

/*
 *
 * Common initialization code
 *
 */
void operaPNG::init(int width, int height) {
	
	graph = (png_bytepp)malloc(height * sizeof(png_bytep));
	if (graph == NULL) {
		throw operaException("operaPNG: ", operaErrorCodeNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	for (int k = 0; k < height; k++)
	{
        graph[k] = (png_bytep)malloc(6*width * sizeof(png_byte));
		if (graph[k] == NULL) {
			throw operaException("operaPNG: ", operaErrorCodeNoMemory, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	
	for (int w = 0; w<width; w++) {
		for (int h = 0; h<height; h++) {
			graph[h][6*w]   = (char)floor(((double)backgroundcolour)/256);
			graph[h][6*w+1] = (char)backgroundcolour%256;
			graph[h][6*w+2] = (char)floor(((double)backgroundcolour)/256);
			graph[h][6*w+3] = (char)backgroundcolour%256;
			graph[h][6*w+4] = (char)floor(((double)backgroundcolour)/256);
			graph[h][6*w+5] = (char)backgroundcolour%256;
		}
	}
}

/*
 * Destructor
 */
operaPNG::~operaPNG(void)
{
	for (int j = 0; j < height; j++) 
		free(graph[j]);
	free(graph);
}

/*
 * open file for writing
 */
void open(const char *Filename, int x, int y, int Backgroundcolour = 65535) {
	operaPNG(Filename, x, y, Backgroundcolour);
}

/*
 * save and close file
 */
void operaPNG::close(void)
{
	
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		throw operaException("operaPNG: Error creating file "+filename+" ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	
	png_init_io(png_ptr, fp);
	
	png_set_compression_level(png_ptr, compressionlevel);
	
	png_set_IHDR(png_ptr, info_ptr, width, height,
				 bitdepth, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
				 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	
	if (filegamma < 1.0e-1) {
		filegamma = 0.5;
	}
	
	png_set_gAMA(png_ptr, info_ptr, filegamma);
	
	time_t          gmt;
	png_time        mod_time;
	time(&gmt);
	png_convert_from_time_t(&mod_time, gmt);
	png_set_tIME(png_ptr, info_ptr, &mod_time);
	
	png_write_info(png_ptr, info_ptr);
	png_write_image(png_ptr, graph);
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fp);
}

/*
 * plot a point
 */
void operaPNG::plot(int x, int y, int red, int green, int blue)
{	
	if (red > 65535) {
		red = 65535;
	}
	if (green > 65535) {
		green = 65535;
	}
	if (blue > 65535) {
		blue = 65535;
	}
	
	if ((bitdepth == 16)) {
		if ((y<=height) && (y>0) && (x>0) && (x<=width) ) {
			graph[height-y][6*x-6] = (char) floor(((double)red)/256);
			graph[height-y][6*x-6+1] = (char)(red%256);
			graph[height-y][6*x-6+2] = (char) floor(((double)green)/256);
			graph[height-y][6*x-6+3] = (char)(green%256);
			graph[height-y][6*x-6+4] = (char) floor(((double)blue)/256);
			graph[height-y][6*x-6+5] = (char)(blue%256);
		}
	} else if ((bitdepth == 8)) {
		if ((y<height+1) && (y>0) && (x>0) && (x<width+1) ) {
			graph[height-y][3*x-3] = (char)(floor(((double)red)/257.0));
			graph[height-y][3*x-3+1] = (char)(floor(((double)green)/257.0));
			graph[height-y][3*x-3+2] = (char)(floor(((double)blue)/257.0));
			
		}
	}
}

void operaPNG::plot(int x, int y, double red, double green, double blue)
{
	plot(x,y,int(red*65535),int(green*65535),int(blue*65535));
}

/*
 * plot a line
 */
void operaPNG::line(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue)
{
	//  Bresenham Algorithm.
	//
	int dy = yto - yfrom;
	int dx = xto - xfrom;
	int stepx, stepy;
	
	if (dy < 0) {
		dy = -dy;  stepy = -1;
	} else {
		stepy = 1;
	}
	
	if (dx < 0) {
		dx = -dx;  stepx = -1;
	} else {
		stepx = 1;
	}
	dy <<= 1;     // dy is now 2*dy
	dx <<= 1;     // dx is now 2*dx
	
	plot(xfrom, yfrom, red, green, blue);
	
	if (dx > dy) {
		int fraction = dy - (dx >> 1);
		
		while (xfrom != xto) {
			if (fraction >= 0) {
				yfrom += stepy;
				fraction -= dx;
			}
			xfrom += stepx;
			fraction += dy;
			plot(xfrom,yfrom,red,green,blue);
		}
	} else {
		int fraction = dx - (dy >> 1);
		while (yfrom != yto) {
			if (fraction >= 0) {
				xfrom += stepx;
				fraction -= dy;
			}
			yfrom += stepy;
			fraction += dx;
			plot(xfrom,yfrom,red,green,blue);
		}
	}
	
}

void operaPNG::line(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue)
{
	line(xfrom, yfrom, xto, yto, (int)(red*65535), (int)(green*65535), (int)(blue*65535));
}

/*
 * plot a square
 */
void operaPNG::square(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue)
{
	line(xfrom, yfrom, xfrom, yto, red, green, blue);
	line(xto, yfrom, xto, yto, red, green, blue);
	line(xfrom, yfrom, xto, yfrom, red, green, blue);
	line(xfrom, yto, xto, yto, red, green, blue);
}

void operaPNG::square(int xfrom, int yfrom, int xto, int yto, double red, double green, double blue)
{
	square(xfrom,  yfrom,  xto,  yto, int(red*65535), int(green*65535), int(blue*65535));
}

void operaPNG::filledsquare(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue)
{
	for (int x = xfrom; x <xto+1; x++) {
		line(x, yfrom, x, yto, red, green, blue);
	}
}

void operaPNG::filledsquare(int xfrom, int yfrom, int xto, int yto, double red, double green, double blue)
{
	filledsquare(xfrom,  yfrom,  xto,  yto, int(red*65535), int(green*65535), int(blue*65535));
}

/*
 * plot a circle
 */
void operaPNG::circle(int xcenter, int ycenter, int radius, int red, int green, int blue)
{
	int x = 0;
	int y = radius;
	int p = (5 - radius*4)/4;
	
	circleaux(xcenter, ycenter, x, y, red, green, blue);
	while (x < y) {
		x++;
		if (p < 0) {
			p += 2*x+1;
		} else {
			y--;
			p += 2*(x-y)+1;
		}
		circleaux(xcenter, ycenter, x, y, red, green, blue);
	}
}

void operaPNG::circle(int xcenter, int ycenter, int radius, double red, double green, double blue)
{
	circle(xcenter,ycenter,radius, int(red*65535), int(green*65535), int(blue*65535));
}

void operaPNG::filledcircle(int xcenter, int ycenter, int radius, int red, int green, int blue)
{
	for (int j = ycenter-radius; j< ycenter+radius+1; j++) {
		line(xcenter - int(sqrt((double)(radius*radius) - (-ycenter + j)*(-ycenter + j ))), j,
			 xcenter + int(sqrt((double)(radius*radius) - (-ycenter + j)*(-ycenter + j ))),j,red,green,blue);
	}
}

void operaPNG::filledcircle(int xcenter, int ycenter, int radius, double red, double green, double blue)
{
	filledcircle(xcenter, ycenter,  radius, int(red*65535), int(green*65535), int(blue*65535));
}

/*
 * plot a triangle
 */
void operaPNG::triangle(int x1, int y1, int x2, int y2, int x3, int y3, int red, int green, int blue)
{
	line(x1, y1, x2, y2, red, green, blue);
	line(x2, y2, x3, y3, red, green, blue);
	line(x3, y3, x1, y1, red, green, blue);
}

void operaPNG::triangle(int x1, int y1, int x2, int y2, int x3, int y3, double red, double green, double blue)
{
	
	line(x1, y1, x2, y2, ((int)65535*red), ((int)65535*green), ((int)65535*blue));
	line(x2, y2, x3, y3, ((int)65535*red), ((int)65535*green), ((int)65535*blue));
	line(x3, y3, x1, y1, ((int)65535*red), ((int)65535*green), ((int)65535*blue));
	
}

void operaPNG::filledtriangle(int x1,int y1,int x2,int y2,int x3,int y3, int red, int green, int blue)
{
	if ((x1==x2 && x2==x3) || (y1==y2 && y2==y3)) 
		return;
	
	if (y2<y1) {
		// x2^=x1^=x2^=x1;
		x2^=x1;
		x1^=x2;
		x2^=x1;
		// y2^=y1^=y2^=y1;
		y2^=y1;
		y1^=x2;
		y2^=y1;
	}
	
	if (y3<y1) {
		//x3^=x1^=x3^=x1;
		x3^=x1;
		x1^=x3;
		x3^=x1;
		//y3^=y1^=y3^=y1;
		y3^=y1;
		y1^=y3;
		y3^=y1;
	}
	
	if (y3<y2) {
		//x2^=x3^=x2^=x3;
		x2^=x3;
		x3^=x2;
		x2^=x3;
		//y2^=y3^=y2^=y3;
		y2^=y3;
		y3^=y2;
		y2^=y3;
	}
	
	if (y2==y3) {
		drawtop(x1, y1, x2, y2, x3, red, green, blue);
	} else {
		if (y1==y3 || y1==y2) {
			drawbottom(x1, y1, x2, x3, y3, red, green, blue);
		} else {
			int newx = x1 + (int)((double)(y2-y1)*(double)(x3-x1)/(double)(y3-y1));
			drawtop(x1, y1, newx, y2, x2, red, green, blue);
			drawbottom(x2, y2, newx, x3, y3, red, green, blue);
		}
	}
	
}

void operaPNG::filledtriangle(int x1,int y1,int x2,int y2,int x3,int y3, double red, double green, double blue)
{
	filledtriangle(x1, y1, x2, y2, x3, y3, (int) (red*65535), (int) (green*65535),  (int) (blue*65535)); 
}

/*
 * public interface
 */
int operaPNG::getheight(void)
{
	return height;
}

int operaPNG::getwidth(void)
{
	return width;
}


int operaPNG::getbitdepth(void)
{
	return bitdepth;
}

int operaPNG::getcolortype(void)
{
	return colortype;
}

void operaPNG::setcompressionlevel(int level)
{
	if ((level < -1) || (level > 9) ) {
		level = 0;
	}
	compressionlevel = level;
}

#ifndef NOFREETYPE

/*
 * Freetype-based text rendering functions.
 */

void operaPNG::plottext(char * facepath, int fontsize, int xstart, int ystart, double angle, char * text, double red, double green, double blue)
{
	FT_Library  library;
	FT_Face     face;
	FT_Matrix   matrix;      // transformation matrix
	FT_Vector   pen;
	
	FT_UInt glyphindex;
	FT_Error error;
	
	FT_Bool usekerning;
	FT_UInt previous = 0;
	
	/* Set up transformation Matrix */
	matrix.xx = (FT_Fixed)(cos(angle)*0x10000);   /* It would make more sense to do this (below), but, bizzarely, */
	matrix.xy = (FT_Fixed)(-sin(angle)*0x10000);   /* if one does, FT_LoadGlyph fails consistently.               */
	matrix.yx = (FT_Fixed)(sin(angle)*0x10000);  //   matrix.yx = - matrix.xy;
	matrix.yy = (FT_Fixed)(cos(angle)*0x10000);  //   matrix.yy = matrix.xx;
	
	/* Place starting coordinates in adequate form. */
	pen.x = xstart*64 ;
	pen.y =   (int)(ystart/64.0);
	
	/*Count the length of the string */
	int numchars = strlen(text);
	
	/* Initialize FT Library object */
	error = FT_Init_FreeType(&library);
	
	/* Initialize FT face object */
	error = FT_New_Face(library,facepath,0,&face);
	
	/* Set the Char size */
	error = FT_Set_Char_Size(face,          /* handle to face object           */
							 0,             /* charwidth in 1/64th of points  */
							 fontsize*64,   /* charheight in 1/64th of points */
							 100,           /* horizontal device resolution    */
							 100);         /* vertical device resolution      */
	
	/* A way of accesing the glyph directly */
	FT_GlyphSlot  slot = face->glyph;  // a small shortcut
	
	/* Does the font file support kerning? */
	usekerning = FT_HAS_KERNING(face);
	
	for (int n = 0; n < numchars; n++ ) {
		/* Convert character code to glyph index */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Retrieve kerning distance and move pen position */
		if (usekerning && previous && glyphindex ) {
			FT_Vector  delta;
			FT_Get_Kerning(face,
						   previous,
						   glyphindex,
						   ft_kerning_default, //FT_KERNING_DEFAULT,
						   &delta);
			
			/* Transform this kerning distance into rotated space */
			pen.x += (int) (((double) delta.x)*cos(angle));
			pen.y +=  (int) (((double) delta.x)*(sin(angle)));
		}
		
		/* Set transform */
		FT_Set_Transform(face, &matrix, &pen);
		
		/*set char size*/
		
		/* Retrieve glyph index from character code */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Load glyph image into the slot (erase previous one) */
		error = FT_Load_Glyph(face, glyphindex, FT_LOAD_DEFAULT);
		
		/* Convert to an anti-aliased bitmap */
		//	error = FTRenderGlyph(face->glyph, FTRENDERMODENORMAL);
		error = FT_Render_Glyph(face->glyph, ft_render_mode_normal);
		
		/* Now, draw to our target surface */
		mydrawbitmap(&slot->bitmap,
					 slot->bitmap_left,
					 ystart + slot->bitmap_top,
					 red,
					 green,
					 blue);
		
		/* Advance to the next position */
		pen.x += slot->advance.x;
		pen.y += slot->advance.y;
		
		/* record current glyph index */
		previous = glyphindex;
	}
	
	/* Free the face and the library objects */
	FT_Done_Face    (face);
	FT_Done_FreeType(library);
}


void operaPNG::plottext(char * facepath, int fontsize, int xstart, int ystart, double angle, char * text, int red, int green, int blue)
{
	plottext(facepath, fontsize, xstart, ystart,  angle,  text,  ((double) red)/65535.0,  ((double) green)/65535.0,  ((double) blue)/65535.0  );
}

void operaPNG::plotcenteredtext(char * facepath, int fontsize, int xcenter, int ystart, double angle, char * text, double red, double green, double blue)
{
	FT_Library  library;
	FT_Face     face;
	FT_Matrix   matrix;      // transformation matrix
	FT_Vector   pen;
	
	FT_UInt glyphindex;
	FT_Error error;
	
	FT_Bool usekerning;
	FT_UInt previous = 0;
	
	/* Set up transformation Matrix */
	matrix.xx = (FT_Fixed)(cos(angle)*0x10000);   /* It would make more sense to do this (below), but, bizzarely, */
	matrix.xy = (FT_Fixed)(-sin(angle)*0x10000);   /* if one does, FT_LoadGlyph fails consistently.               */
	matrix.yx = (FT_Fixed)(sin(angle)*0x10000);  //   matrix.yx = - matrix.xy;
	matrix.yy = (FT_Fixed)(cos(angle)*0x10000);  //   matrix.yy = matrix.xx;
	
	/* Place starting coordinates in adequate form. */
	pen.x = xcenter*64 ;
	pen.y =   (int)(ystart/64.0);
	
	/*Count the length of the string */
	int numchars = strlen(text);
	
	/* Initialize FT Library object */
	error = FT_Init_FreeType(&library);
	
	/* Initialize FT face object */
	error = FT_New_Face(library,facepath,0,&face);
	
	/* Set the Char size */
	error = FT_Set_Char_Size(face,          /* handle to face object           */
							 0,             /* charwidth in 1/64th of points  */
							 fontsize*64,   /* charheight in 1/64th of points */
							 100,           /* horizontal device resolution    */
							 100);         /* vertical device resolution      */
	
	/* A way of accesing the glyph directly */
	FT_GlyphSlot  slot = face->glyph;  // a small shortcut
	
	/* Does the font file support kerning? */
	usekerning = FT_HAS_KERNING(face);
	// Pass 1 calc distane
	for (int n = 0; n < numchars; n++ ) {
		/* Convert character code to glyph index */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Retrieve kerning distance and move pen position */
		if (usekerning && previous&& glyphindex ) {
			FT_Vector  delta;
			FT_Get_Kerning(face,
						   previous,
						   glyphindex,
						   ft_kerning_default, //FT_KERNING_DEFAULT,
						   &delta);
			
			/* Transform this kerning distance into rotated space */
			pen.x += (int) (delta.x);
			pen.y +=  0;
		}
		
		/* Set transform */
		FT_Set_Transform(face, &matrix, &pen);
		
		/*set char size*/
		
		/* Retrieve glyph index from character code */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Load glyph image into the slot (erase previous one) */
		error = FT_Load_Glyph(face, glyphindex, FT_LOAD_DEFAULT);
		
		/* Convert to an anti-aliased bitmap */
		//	error = FTRenderGlyph(face->glyph, FTRENDERMODENORMAL);
		error = FT_Render_Glyph(face->glyph, ft_render_mode_normal);
		
		/* Advance to the next position */
		pen.x += slot->advance.x;
		pen.y += slot->advance.y;
		
		/* record current glyph index */
		previous = glyphindex;
	}
	pen.x -= (pen.x - xcenter*64) / 2;
	// Pass 2 render
	for (int n = 0; n < numchars; n++ ) {
		/* Convert character code to glyph index */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Retrieve kerning distance and move pen position */
		if (usekerning && previous && glyphindex ) {
			FT_Vector  delta;
			FT_Get_Kerning(face,
						   previous,
						   glyphindex,
						   ft_kerning_default, //FT_KERNING_DEFAULT,
						   &delta);
			
			/* Transform this kerning distance into rotated space */
			pen.x += (int) (((double) delta.x)*cos(angle));
			pen.y +=  (int) (((double) delta.x)*(sin(angle)));
		}
		
		/* Set transform */
		FT_Set_Transform(face, &matrix, &pen);
		
		/*set char size*/
		
		/* Retrieve glyph index from character code */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Load glyph image into the slot (erase previous one) */
		error = FT_Load_Glyph(face, glyphindex, FT_LOAD_DEFAULT);
		
		/* Convert to an anti-aliased bitmap */
		//	error = FTRenderGlyph(face->glyph, FTRENDERMODENORMAL);
		error = FT_Render_Glyph(face->glyph, ft_render_mode_normal);
		
		/* Now, draw to our target surface */
		mydrawbitmap(&slot->bitmap,
					 slot->bitmap_left,
					 ystart + slot->bitmap_top,
					 red,
					 green,
					 blue);
		
		/* Advance to the next position */
		pen.x += slot->advance.x;
		pen.y += slot->advance.y;
		
		/* record current glyph index */
		previous = glyphindex;
	}
	
	/* Free the face and the library objects */
	FT_Done_Face    (face);
	FT_Done_FreeType(library);
}


void operaPNG::plotcenteredtext(char * facepath, int fontsize, int xcenter, int ystart, double angle, char * text, int red, int green, int blue)
{
	plotcenteredtext(facepath, fontsize, xcenter, ystart,  angle,  text,  ((double) red)/65535.0,  ((double) green)/65535.0,  ((double) blue)/65535.0  );
}

/*
 * Get text width
 */
int operaPNG::gettextwidth(char * facepath, int fontsize, char * text)
{
	
	FT_Library  library;
	FT_Face     face;
	FT_Matrix   matrix;      // transformation matrix
	FT_Vector   pen;
	
	FT_UInt glyphindex;
	FT_Error error;
	
	FT_Bool usekerning;
	FT_UInt previous = 0;
	
	/* Set up transformation Matrix */
	matrix.xx = (FT_Fixed)(1.0*0x10000);   /* It would make more sense to do this (below), but, bizzarely, */
	matrix.xy = (FT_Fixed)(0.0*0x10000);   /* if one does, FT_LoadGlyph fails consistently.               */
	matrix.yx = (FT_Fixed)(0.0*0x10000);  //   matrix.yx = - matrix.xy;
	matrix.yy = (FT_Fixed)(1.0*0x10000);  //   matrix.yy = matrix.xx;
	
	/* Place starting coordinates in adequate form. */
	pen.x = 0;
	pen.y = 0;
	
	/*Count the length of the string */
	int numchars = strlen(text);
	
	/* Initialize FT Library object */
	error = FT_Init_FreeType(&library);
	
	/* Initialize FT face object */
	error = FT_New_Face(library,facepath,0,&face);
	
	/* Set the Char size */
	error = FT_Set_Char_Size(face,       /* handle to face object           */
							 0,             /* charwidth in 1/64th of points  */
							 fontsize*64,   /* charheight in 1/64th of points */
							 100,           /* horizontal device resolution    */
							 100);         /* vertical device resolution      */
	
	/* A way of accesing the glyph directly */
	FT_GlyphSlot  slot = face->glyph;  // a small shortcut
	
	/* Does the font file support kerning? */
	usekerning = FT_HAS_KERNING(face);
	
	for (int n = 0; n < numchars; n++ ) {
		/* Convert character code to glyph index */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Retrieve kerning distance and move pen position */
		if (usekerning && previous&& glyphindex ) {
			FT_Vector  delta;
			FT_Get_Kerning(face,
						   previous,
						   glyphindex,
						   ft_kerning_default, //FT_KERNING_DEFAULT,
						   &delta);
			
			/* Transform this kerning distance into rotated space */
			pen.x += (int) (delta.x);
			pen.y +=  0;
		}
		
		/* Set transform */
		FT_Set_Transform(face, &matrix, &pen);
		
		/*set char size*/
		
		/* Retrieve glyph index from character code */
		glyphindex = FT_Get_Char_Index(face, text[n]);
		
		/* Load glyph image into the slot (erase previous one) */
		error = FT_Load_Glyph(face, glyphindex, FT_LOAD_DEFAULT);
		
		/* Convert to an anti-aliased bitmap */
		//	error = FTRenderGlyph(face->glyph, FTRENDERMODENORMAL);
		error = FT_Render_Glyph(face->glyph, ft_render_mode_normal);
		
		/* Advance to the next position */
		pen.x += slot->advance.x;
		pen.y += slot->advance.y;
		
		/* record current glyph index */
		previous = glyphindex;
	}
	/* Free the face and the library objects */
	FT_Done_Face(face);
	FT_Done_FreeType(library);
	
	return (int)(((double)pen.x)/64.0);
}
#endif

/*
 * private helper routines
 */

#ifndef NOFREETYPE

void operaPNG::mydrawbitmap(FT_Bitmap *bitmap, int x, int y, double red, double green, double blue)
{
	double temp;
	for (int j=1; j<bitmap->rows+1; j++) {
		for (int i=1; i< bitmap->width + 1; i++) {
			temp = (double)(bitmap->buffer[(j-1)*bitmap->width + (i-1)] )/255.0;
			if (temp) {
				plot(x + i, y  - j,
					 temp*red + (1-temp)*(doubleread(x+i,y-j,1)),
					 temp*green + (1-temp)*(doubleread(x+i,y-j,2)),
					 temp*blue + (1-temp)*(doubleread(x+i,y-j,3)));
			}
		}
	}
}
#endif

int operaPNG::read(int x, int y, int colour) {
	int temp1,temp2;
	
	if ((colour !=1 ) && (colour !=2 ) && (colour !=3 )) {
		return 0;
	}
	
	if ((x>0 ) && (x <= (width) ) && (y>0 ) && (y <= (height) ) ) {
		
		if (bitdepth == 16) {
			temp2=6*(x-1);
			if (colour == 1) {
				temp1 = (graph[(height-y)][temp2])*256 + graph[height-y][temp2+1];
				return temp1;
			}
			
			if (colour == 2) {
				temp1 = (graph[height-y][temp2+2])*256 + graph[height-y][temp2+3];
				return temp1;
			}
			
			if (colour == 3) {
				temp1 = (graph[height-y][temp2+4])*256 + graph[height-y][temp2+5];
				return temp1;
			}
		}
		
		if (bitdepth == 8) {
			temp2=3*(x-1);
			if (colour == 1) {
				temp1 = graph[height-y][temp2];
				return temp1*256;
			}
			
			if (colour == 2) {
				temp1 =  graph[height-y][temp2+1];
				return temp1*256;
			}
			
			if (colour == 3) {
				temp1 =  graph[height-y][temp2+2];
				return temp1*256;
			}
		}
	} 
	return 0;
}

int operaPNG::read(int x, int y)
{
	int temp1,temp2,temp3,temp4,temp5;
	
	if ((x>0 ) && (x <= (width) ) && (y>0 ) && (y <= (height) )) {
		if (bitdepth == 16) {
			temp5=6*x;
			temp1 = (graph[(height-y)][temp5-6])*256 + graph[height-y][temp5-5];
			temp2 = (graph[height-y][temp5-4])*256 + graph[height-y][temp5-3];
			temp3 = (graph[height-y][temp5-2])*256 + graph[height-y][temp5-1];
			temp4 =  int((temp1+temp2+temp3)/3.0);
		} else if (bitdepth == 8) {
			temp5 = 3*x;
			temp1 = graph[height-y][temp5-3];
			temp2 = graph[height-y][temp5-2];
			temp3 = graph[height-y][temp5-1];
			temp4 = int((temp1+temp2+temp3)/3.0);
		} else {
			temp4 = 0;
		}
		return temp4;
		
	}
	return 0;
}

double operaPNG::doubleread(int x, int y, int colour)
{
	return double(read(x,y,colour))/65535.0;
}

double operaPNG::doubleread(int x, int y)
{
	return double(read(x,y))/65535.0;
}

/*
 * drawtop(), drawbottom() and filledtriangle() were contributed by Gurkan Sengun
 * (<gurkan@linuks.mine.nu>, http://www.linuks.mine.nu/ )
 */
void operaPNG::drawtop(long x1,long y1,long x2,long y2,long x3, int red, int green, int blue)
{
	if (x2>x3) {
		x2^=x3;
		x3^=x2;
		x2^=x3;
	}
	
	long posl = x1*256;
	long posr = posl;
	
	long cl=((x2-x1)*256)/(y2-y1);
	long cr=((x3-x1)*256)/(y2-y1);
	
	for (int y=y1; y<y2; y++) {
		line(posl/256, y, posr/256, y, red, green, blue);
		posl+=cl;
		posr+=cr;
	}
}

void operaPNG::drawbottom(long x1,long y1,long x2,long x3,long y3, int red, int green, int blue)
{
	if (x1>x2) {
		x2^=x1;
		x1^=x2;
		x2^=x1;
    }
	
	long posl=x1*256;
	long posr=x2*256;
	
	long cl=((x3-x1)*256)/(y3-y1);
	long cr=((x3-x2)*256)/(y3-y1);
	
	for (int y=y1; y<y3; y++) {
		line(posl/256, y, posr/256, y, red, green, blue);
		
		posl+=cl;
		posr+=cr;
	}
}

void operaPNG::circleaux(int xcenter, int ycenter, int x, int y, int red, int green, int blue)
{
	if (x == 0) {
		plot(xcenter, ycenter + y, red, green, blue);
		plot(xcenter, ycenter - y, red, green, blue);
		plot(xcenter + y, ycenter, red, green, blue);
		plot(xcenter - y, ycenter, red, green, blue);
	} else if (x == y) {
		plot(xcenter + x, ycenter + y, red, green, blue);
		plot(xcenter - x, ycenter + y, red, green, blue);
		plot(xcenter + x, ycenter - y, red, green, blue);
		plot(xcenter - x, ycenter - y, red, green, blue);
	} else if (x < y) {
		plot(xcenter + x, ycenter + y, red, green, blue);
		plot(xcenter - x, ycenter + y, red, green, blue);
		plot(xcenter + x, ycenter - y, red, green, blue);
		plot(xcenter - x, ycenter - y, red, green, blue);
		plot(xcenter + y, ycenter + x, red, green, blue);
		plot(xcenter - y, ycenter + x, red, green, blue);
		plot(xcenter + y, ycenter - x, red, green, blue);
		plot(xcenter - y, ycenter - x, red, green, blue);
	}
	
}


