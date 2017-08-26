/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaGeometricShapesTest
 Version: 1.0
 Description: Perform various tests on the operaGeometricShapes class.
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

#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaGeometricShapes.h"

#define MAXIMAGES 100

/*! \file operaGeometricShapesTest.cpp */

using namespace std;

/*! 
 * operaGeometricShapesTest
 * \author Eder Martioli
 * \brief Perform various tests on the operaGeometricShapes class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
    try {
        cout << "Test Module for operaGeometricShapes" << endl << endl;
      
        /*
         * Test Point
         */        
        operaPoint TestPoint(-3.9999,0.0);
    
        cout << "Test Point: (" << TestPoint.getXcoord() << "," << TestPoint.getYcoord() << ")" << endl << endl;

        /*
         * Test Circle
         */
        operaPoint circleCenter(-3,4);
        float radius = 4;
        
        Circle *mycircle = new Circle(radius, circleCenter);
        
        cout << "Circle radius = " << mycircle->getRadius() << endl;
        cout << "Circle center = (" << mycircle->getCenter()->getXcoord()  << "," << mycircle->getCenter()->getYcoord()  << ")" << endl;        
        
        if(mycircle->pointInCircle(TestPoint)) {
            cout << "Test Point is INside of circle" << endl;  
        } else if(!mycircle->pointInCircle(TestPoint)) {
            cout << "Test Point is OUTside of circle" << endl;        
        }
        
        cout << endl;
        
        delete mycircle;
        
        /*
         * Test Rectangle
         */  
        float width = 4;
        float height = 3;
        float rectangleAngle = 0;
        operaPoint rectangleCenter(-4,0);
        
        Rectangle *myrectangle = new Rectangle(width,height,rectangleAngle,rectangleCenter);
        
        cout << "Rectangle width = " << myrectangle->getWidth() << endl;
        cout << "Rectangle height = " << myrectangle->getHeight() << endl;
        cout << "Rectangle angle = " << myrectangle->getAngle() << endl;         
        cout << "Rectangle center = (" << myrectangle->getCenter()->getXcoord()  << "," << myrectangle->getCenter()->getYcoord()  << ")" << endl;        
        
        operaPoint *corners[FOURSIDES];
        for(unsigned i=0;i<FOURSIDES;i++) {
            corners[i] = myrectangle->getCorner(i);
        }
        cout << "Corners of Rectangle are: " << endl;
        myrectangle->printCorners();
        if(myrectangle->pointInRectangle(TestPoint)) {
            cout << "Test Point is INside of rectangle" << endl;  
        } else if(!myrectangle->pointInRectangle(TestPoint)) {
            cout << "Test Point is OUTside of rectangle" << endl;        
        }
        
        cout << endl;    
        delete myrectangle;
        
        /*
         * Test Polygon
         */        
        unsigned nsides = 4;
        operaPoint polygoncorners[4];
 
        polygoncorners[0].setPoint(-4,0);
        polygoncorners[1].setPoint(-1,6);
        polygoncorners[2].setPoint(5,2);
        polygoncorners[3].setPoint(1,-2);
        
        Polygon *myPolygon  = new Polygon(nsides,polygoncorners); 

        cout << "Corners of Polygon are: " << endl;
        myPolygon->printVertexCoordinates();
        
        if(myPolygon->pointInPolygon(TestPoint)) {
            cout << "Test Point is INside of polygon" << endl;  
        } else if(!myPolygon->pointInPolygon(TestPoint)) {
            cout << "Test Point is OUTside of polygon" << endl;        
        }
        cout << endl;
  
        delete myPolygon;
 
        /*
         * Test Line
         */  
        float lineSlope = -0.54;
        float lineWidth = 2;
        float lineLength = 6;
        operaPoint lineMidPoint(-1,0);
        
        Line *myline = new Line(lineSlope,lineWidth,lineLength,lineMidPoint);  
        
        cout << "Line slope = " << myline->getSlope() << endl;
        cout << "Line width = " << myline->getWidth() << endl;
        cout << "Line length = " << myline->getLength() << endl;         
        cout << "Line Middle Point = (" << myline->getMidPoint()->getXcoord()  << "," << myline->getMidPoint()->getYcoord()  << ")" << endl;        

        cout << "Line equation = ";
        myline->printLineEquation();
        
        if(myline->pointOnLine(TestPoint)) {
            cout << "Test Point is ON the line" << endl;  
        } else if(!myline->pointOnLine(TestPoint)) {
            cout << "Test Point is NOT ON the line" << endl;        
        }        
        
        if(myline->pointInLineWidth(TestPoint, *myline)) {
            cout << "Test Point is INside the line width" << endl;  
        } else if(!myline->pointInLineWidth(TestPoint, *myline)) {
            cout << "Test Point is OUTside the line width" << endl;        
        }           
        
		Line topline(*myline, line_top);
		Line bottomline(*myline, line_bottom);
		Line leftline(*myline, line_left);
		Line rightline(*myline, line_right);
		Line perpendicularline(*myline, line_perpendicular);

        cout << "topline slope = " << topline.getSlope() << endl;
        cout << "topline width = " << topline.getWidth() << endl;
        cout << "topline length = " << topline.getLength() << endl;         
        cout << "topline Middle Point = (" << topline.getMidPoint()->getXcoord()  << "," << topline.getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "bottomline slope = " << bottomline.getSlope() << endl;
        cout << "bottomline width = " << bottomline.getWidth() << endl;
        cout << "bottomline length = " << bottomline.getLength() << endl;         
        cout << "bottomline Middle Point = (" << bottomline.getMidPoint()->getXcoord()  << "," << bottomline.getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "leftline slope = " << leftline.getSlope() << endl;
        cout << "leftline width = " << leftline.getWidth() << endl;
        cout << "leftline length = " << leftline.getLength() << endl;         
        cout << "leftline Middle Point = (" << leftline.getMidPoint()->getXcoord()  << "," << leftline.getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "rightline slope = " << rightline.getSlope() << endl;
        cout << "rightline width = " << rightline.getWidth() << endl;
        cout << "rightline length = " << rightline.getLength() << endl;         
        cout << "rightline Middle Point = (" << rightline.getMidPoint()->getXcoord()  << "," << rightline.getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "perpendicularline slope = " << perpendicularline.getSlope() << endl;
        cout << "perpendicularline width = " << perpendicularline.getWidth() << endl;
        cout << "perpendicularline length = " << perpendicularline.getLength() << endl;         
        cout << "perpendicularline Middle Point = (" << perpendicularline.getMidPoint()->getXcoord()  << "," << perpendicularline.getMidPoint()->getYcoord()  << ")" << endl;        
		
		Line *topline_p = myline->newTopLine(*myline);
		Line *bottomline_p = myline->newBottomLine(*myline);
		Line *leftline_p = myline->newLeftLine(*myline);
		Line *rightline_p = myline->newRightLine(*myline);
		Line *perpendicularline_p = myline->newPerpendicularLine(*myline);
		
        cout << "topline_p slope = " << topline_p->getSlope() << endl;
        cout << "topline_p width = " << topline_p->getWidth() << endl;
        cout << "topline_p length = " << topline_p->getLength() << endl;         
        cout << "topline_p Middle Point = (" << topline_p->getMidPoint()->getXcoord()  << "," << topline_p->getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "bottomline_p slope = " << bottomline_p->getSlope() << endl;
        cout << "bottomline_p width = " << bottomline_p->getWidth() << endl;
        cout << "bottomline_p length = " << bottomline_p->getLength() << endl;         
        cout << "bottomline Middle Point = (" << bottomline_p->getMidPoint()->getXcoord()  << "," << bottomline_p->getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "leftline_p slope = " << leftline_p->getSlope() << endl;
        cout << "leftline_p width = " << leftline_p->getWidth() << endl;
        cout << "leftline_p length = " << leftline_p->getLength() << endl;         
        cout << "leftline_p Middle Point = (" << leftline_p->getMidPoint()->getXcoord()  << "," << leftline_p->getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "rightline_p slope = " << rightline_p->getSlope() << endl;
        cout << "rightline_p width = " << rightline_p->getWidth() << endl;
        cout << "rightline_p length = " << rightline_p->getLength() << endl;         
        cout << "rightline_p Middle Point = (" << rightline_p->getMidPoint()->getXcoord()  << "," << rightline_p->getMidPoint()->getYcoord()  << ")" << endl;        
		
        cout << "perpendicularline_p slope = " << perpendicularline_p->getSlope() << endl;
        cout << "perpendicularline_p width = " << perpendicularline_p->getWidth() << endl;
        cout << "perpendicularline_p length = " << perpendicularline_p->getLength() << endl;         
        cout << "perpendicularline_p Middle Point = (" << perpendicularline_p->getMidPoint()->getXcoord()  << "," << perpendicularline_p->getMidPoint()->getYcoord()  << ")" << endl;        
		
		cout << endl;        
	}
	catch (operaException e) {
		cerr << "operaGeometricShapes: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaGeometricShapes: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
//static void printUsageSyntax() {
//
//	cerr << " Usage: operaGeometricShapesTest [--image=<filename>]+ --output=<file name> -[dvth]\n";
//}	
