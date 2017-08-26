#ifndef LIBOPERAPNGH
#define LIBOPERAPNGH
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
 Contact: opera@cfht.hawaii.edu
 
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

#include <png.h>

#ifndef NOFREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#endif

#define PNG_BYTES_TO_CHECK (4)
#define DEFAULT_COMPRESSION (6)

/*! 
 * \sa class operaPNG
 * \brief Encapsulates a PNG image. 
 * \file operaPNG.h
 * \ingroup libraries
 */

class operaPNG {

private:
	
	string filename;   
	
	int width;
	int height;
	int backgroundcolour;
	int compressionlevel;
	int bitdepth;
	int colortype;
	int rowbytes;
	
	double filegamma;
	double screengamma;
	
	unsigned char **graph;

	void init(int width, int height);
	void floodfillinternal(int xstart, int ystart,  double startred, double startgreen, double startblue, double fillred, double fillgreen, double fillblue);
	void floodfillinternalblend(int xstart, int ystart, double opacity,  double startred, double startgreen, double startblue, double fillred, double fillgreen, double fillblue);
	
#ifndef NOFREETYPE
	void mydrawbitmap(FT_Bitmap * bitmap, int x, int y, double red, double green, double blue);
#endif
	
	/* drwatop(), drawbottom() and filledtriangle() were contributed by Gurkan Sengun
	 * ( <gurkan@linuks.mine.nu>, http://www.linuks.mine.nu/ )
	 *
	 */
	void drawtop(long x1,long y1,long x2,long y2,long x3, int red, int green, int blue);
	void drawbottom(long x1,long y1,long x2,long x3,long y3, int red, int green, int blue);
	void drawbottomblend(long x1,long y1,long x2,long x3,long y3, double opacity, int red, int green, int blue);
	void drawtopblend(long x1,long y1,long x2,long y2,long x3, double opacity, int red, int green, int blue);
	void circleaux(int xcenter, int ycenter, int x, int y, int red, int green, int blue);
	int  read(int x, int y, int colour);
	int  read(int x, int y);
	double doubleread(int x, int y, int colour);	
	double doubleread(int x, int y);
	
public:
	
	operaPNG();   
	operaPNG(const char *filename, int width, int height, int backgroundcolour);   
	operaPNG(const char *filename, int x, int y, double backgroundcolour);
	~operaPNG();  
	
	void open(const char *filename, int x, int y, int backgroundcolour); 
	void close(void); 
	
	void  plot(int x, int y, int red, int green, int blue); 
	void  plot(int x, int y, double red, double green, double blue); 
	
	void line(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
	void line(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);
	
	void triangle(int x1, int y1, int x2, int y2, int x3, int y3, int red, int green, int blue);
	void triangle(int x1, int y1, int x2, int y2, int x3, int y3, double red, double green, double blue);
	
	void square(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
	void square(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);
	
	void filledsquare(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
	void filledsquare(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);
	
	void circle(int xcentre, int ycentre, int radius, int red, int green, int blue);
	void circle(int xcentre, int ycentre, int radius, double red, double green, double blue);
	
	void filledcircle(int xcentre, int ycentre, int radius, int red, int green, int blue);
	void filledcircle(int xcentre, int ycentre, int radius, double red, double green, double blue);
	
	void filledtriangle(int x1,int y1,int x2,int y2,int x3,int y3, int red, int green, int blue);
	void filledtriangle(int x1,int y1,int x2,int y2,int x3,int y3, double red, double green, double blue);
	
	void plottext(char * facepath, int fontsize, int xstart, int ystart, double angle, char * text, double red, double green, double blue);
	void plottext(char * facepath, int fontsize, int xstart, int ystart, double angle, char * text, int red, int green, int blue);
	void plotcenteredtext(char * facepath, int fontsize, int xcenter, int ystart, double angle, char * text, double red, double green, double blue);
	void plotcenteredtext(char * facepath, int fontsize, int xcenter, int ystart, double angle, char * text, int red, int green, int blue);
	int gettextwidth(char * facepath, int fontsize,  char * text);
	
	int getheight(void);
	int getwidth(void);
	
    void setcompressionlevel(int level);
	
	int getbitdepth(void);
	int getcolortype(void);
	
};

#endif

