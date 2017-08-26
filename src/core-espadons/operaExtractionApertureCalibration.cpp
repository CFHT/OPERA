/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaExtractionApertureCalibration
 Version: 1.0
 Description: Calibrate aperture for extraction
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

#include <fstream>
#include "libraries/operaIOFormats.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaExtractionApertureCalibration.cpp */

using namespace std;

/*!
 * operaExtractionApertureCalibration
 * \author Eder Martioli
 * \brief Tool to calibrate the aperture for extraction.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

void GenerateExtractionAperturePlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display);

void GenerateExtractionApertureTiltPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName1, string dataFileName2, float tiltAngle, float tiltAngleError, bool display);

bool CalculateTilt(operaSpectralOrder *spectralOrder, double width, double height, unsigned nRowSamples, unsigned pickImageRow, unsigned xbin, double& tiltAngle, double& tiltAngleError, double& yintercept, double& fluxFraction, ofstream& ftiltdata1, bool debug);

operaExtractionAperture<Line> CreateLineAperture(operaPoint apertureMidPoint, double apertureSlope, double apertureHeight, double apertureWidth, operaInstrumentProfile* instrumentProfile, double yCenter, double normalizationFactor);

operaExtractionAperture<Line> ProcessAperture(string apertureName, const Line& extractionLine, double apertureXoffset, double apertureWidth, operaInstrumentProfile* instrumentProfile, double yCenter, double normalizationFactor, DMatrix& ipImage, unsigned ordernumber, bool verbose);

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
    string outputApertureFile;
	string inputgeom;
	string inputprof;
	string inputorderspacing;
    unsigned nRowSamples = 1;
	unsigned pickImageRow = 0;
    unsigned xbin = 10;
    unsigned numberOfBeams = 1; // 1 for star-only;  2 for polar/s+s
    double gapBetweenBeams = 0;
    double apertureWidth = 26.0;
    double apertureHeight = 0.6;
    double backgroundAperture = 2.0;
    bool constantTilt = false;
	double tiltValue = 0;
	int ordernumber = NOTPROVIDED;
	int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    string plotfilename;
	string datafilename;
	string scriptfilename;
	bool interactive = false;
    string tiltplotfilename;
    string tiltdata1filename;
    string tiltdata2filename;
    string tiltscriptfilename;
	bool applyoffset;
    
    args.AddRequiredArgument("outputApertureFile", outputApertureFile, "Output aperture file name");
    args.AddRequiredArgument("inputgeom", inputgeom, "Input geometry file");
    args.AddRequiredArgument("inputprof", inputprof, "Input instrument profile file");
    args.AddOptionalArgument("inputorderspacing", inputorderspacing, "", "Input order spacing file");
    args.AddRequiredArgument("nRowSamples", nRowSamples, "Number of equally spaced rows to sample IP");
    args.AddOptionalArgument("pickImageRow", pickImageRow, 0, "Specific single row to use for IP model (overrides nRowSamples)");
    args.AddRequiredArgument("xbin", xbin, "Number of IP points to bin in x-direction");
    args.AddRequiredArgument("numberOfBeams", numberOfBeams, "Number of beams to split aperture");
    args.AddRequiredArgument("gapBetweenBeams", gapBetweenBeams, "Gap between beams in pixel units");
    args.AddRequiredArgument("apertureWidth", apertureWidth, "Aperture width in pixel units");
    args.AddRequiredArgument("apertureHeight", apertureHeight, "Aperture height in pixel units");
    args.AddRequiredArgument("backgroundAperture", backgroundAperture, "Aperture width for background in pixel units");
    args.AddOptionalArgument("constantTilt", constantTilt, false, "Use the provided tilt angle for all orders");
	args.AddOptionalArgument("tilt", tiltValue, NAN, "Tilt angle to set for all orders");
	args.AddOptionalArgument("applyoffset", applyoffset, false, "Apply a subpixel shift to aperture to match IP line fit");
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
    args.AddOptionalArgument("tiltplotfilename", tiltplotfilename, "", "Tilt plot eps file name");
    args.AddOptionalArgument("tiltdata1filename", tiltdata1filename, "", "Tilt first data file name");
    args.AddOptionalArgument("tiltdata2filename", tiltdata2filename, "", "Tilt second data file name");
    args.AddOptionalArgument("tiltscriptfilename", tiltscriptfilename, "", "Tilt gnuplot script file name");
	
	try {
		args.Parse(argc, argv);
		
		// we need a geometry file...
		if (inputgeom.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an instrument profile file...
		if (inputprof.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		if (args.verbose) {
			cout << "operaExtractionApertureCalibration: outputApertureFile = " << outputApertureFile << endl;
			cout << "operaExtractionApertureCalibration: inputorderspacing = " << inputorderspacing << endl;
			cout << "operaExtractionApertureCalibration: inputgeom = " << inputgeom << endl;
			cout << "operaExtractionApertureCalibration: inputprof = " << inputprof << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaExtractionApertureCalibration: ordernumber = " << ordernumber << endl;
            if(pickImageRow) cout << "operaExtractionApertureCalibration: pickImageRow = " << pickImageRow << endl;
            else cout << "operaExtractionApertureCalibration: nRowSamples = " << nRowSamples << endl;
            cout << "operaExtractionApertureCalibration: xbin = " << xbin << endl;
			cout << "operaExtractionApertureCalibration: numberOfBeams = " << numberOfBeams << endl;
			cout << "operaExtractionApertureCalibration: apertureWidth = " << apertureWidth << endl;
			cout << "operaExtractionApertureCalibration: apertureHeight = " << apertureHeight << endl;
			cout << "operaExtractionApertureCalibration: backgroundAperture = " << backgroundAperture << endl;
			cout << "operaExtractionApertureCalibration: constantTilt = " << constantTilt << endl;
			if (constantTilt) cout << "operaExtractionApertureCalibration: tilt = " << tiltValue << endl;
            if(args.plot) {
                cout << "operaExtractionApertureCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaExtractionApertureCalibration: datafilename = " << datafilename << endl;
                cout << "operaExtractionApertureCalibration: scriptfilename = " << scriptfilename << endl;
                cout << "operaExtractionApertureCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
                cout << "operaExtractionApertureCalibration: tiltplotfilename = " << tiltplotfilename << endl;
                cout << "operaExtractionApertureCalibration: tiltdata1filename = " << tiltdata1filename << endl;
                cout << "operaExtractionApertureCalibration: tiltdata2filename = " << tiltdata2filename << endl;
                cout << "operaExtractionApertureCalibration: tiltscriptfilename = " << tiltscriptfilename << endl;
            }
		}
        
        ofstream fdata;
        ofstream ftiltdata1;
        ofstream ftiltdata2;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        if (!tiltdata1filename.empty()) ftiltdata1.open(tiltdata1filename.c_str());
        if (!tiltdata2filename.empty()) ftiltdata2.open(tiltdata2filename.c_str());
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgeom);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputprof);
        if(!inputorderspacing.empty()) operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputorderspacing);
        
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		
        // Calculate the tilt for each order
        operaVector tiltAngleGood;
		operaVector offsetGood;
        std::map <int, double> offsets;
        for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			double tilt, tilterror, yintercept, fluxfraction;
			if(CalculateTilt(spectralOrder, apertureWidth, apertureHeight, nRowSamples, pickImageRow, xbin, tilt, tilterror, yintercept, fluxfraction, ftiltdata1, args.debug)) {
				if(args.verbose) cout << "operaExtractionApertureCalibration: # order = " << order << " tilt = " << tilt << " +/- " << 180*tilterror/M_PI << " FluxFraction=" << fluxfraction << endl;
				tiltAngleGood.insert(tilt);
				offsetGood.insert(yintercept);
				spectralOrder->setTiltInDegrees(tilt, tilterror);
				offsets[order] = yintercept;
				if(!tiltdata2filename.empty()) {
					ftiltdata2 << order << " " << tilt << " " << tilterror  << " " << fluxfraction << endl;
				}
			} else {
				cout << "operaExtractionApertureCalibration: WARNING can't calculate tilt for order " << order << endl;
				spectralOrder->setTiltInDegrees(NAN, NAN);
			}
		}
		
		double mediantilt = Median(tiltAngleGood);
		double mediantilterror = MedianStdDev(tiltAngleGood, mediantilt);
		double medianoffset = Median(offsetGood);
		if(args.verbose) cout << "operaExtractionApertureCalibration: # Final tilt = " << mediantilt << " +/- " << mediantilterror << endl;
		
		// Assign the median final tilt to orders that have been skipped, or set a constant tilt for all orders
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if(constantTilt && !isnan(tiltValue)) {
				spectralOrder->setTiltInDegrees(tiltValue, 0);
				offsets[order] = 0;
			}
			else if(constantTilt || isnan(spectralOrder->getTiltInDegreesValue())) {
				spectralOrder->setTiltInDegrees(mediantilt, mediantilterror);
				offsets[order] = medianoffset;
			}
		}
		
		// Calculate apertures and save the information
        unsigned ystackIndex = 0;
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                if(args.verbose) cout << "operaExtractionApertureCalibration: Calculating the aperture of order " << order << endl;
                spectralOrder->setnumberOfBeams(numberOfBeams);
                operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                
                // Set up FITS image to store IP data - consider IP at center of order
                unsigned NXPoints = instrumentProfile->getNXPoints();
                unsigned NYPoints = instrumentProfile->getNYPoints();
                DMatrix ipImage(NYPoints, NXPoints);
                
                operaGeometry *Geometry = spectralOrder->getGeometry();
                float ycenter = (Geometry->getYmax() + Geometry->getYmin())/2 + 0.5;
                if(pickImageRow) ycenter = (float)pickImageRow + 0.5;
                
                double normalizationFactor = 0;
                for (unsigned j=0; j<NYPoints; j++) {	 
                    for (unsigned i=0; i<NXPoints; i++) {
                        ipImage[j][i] = instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                        normalizationFactor += ipImage[j][i];
                    }
                }
                
                // Set up line for full aperture
                double tiltRadians = spectralOrder->getTiltInDegreesValue()*M_PI/180;
                operaPoint midpoint;
				if(applyoffset) midpoint.setPoint(0, offsets[order]);
				Line extractionLine(tan(tiltRadians), apertureHeight, apertureWidth, midpoint);
                double widthX = apertureWidth*cos(tiltRadians);
                
                // Set up beam apertures
                for(unsigned k=0;k<numberOfBeams;k++) {
                    double beamWidth = apertureWidth/numberOfBeams - gapBetweenBeams;
                    double aperXsize = widthX/numberOfBeams;
					ostringstream ss;
					ss << "beam = " << k;
					operaExtractionAperture<Line> beamAperture = ProcessAperture(ss.str(), extractionLine, -widthX/2 + aperXsize*(k+0.5), beamWidth, instrumentProfile, ycenter, normalizationFactor, ipImage, order, args.verbose);
					// DT May 14 2014, if we have no beam flux, we have no aperture...
					if (beamAperture.getFluxFraction() > 0) {
						spectralOrder->setExtractionApertures(k, new operaExtractionAperture<Line>(beamAperture));
						spectralOrder->sethasExtractionApertures(true);
					} else {
						spectralOrder->sethasExtractionApertures(false);
						cout << "operaExtractionApertureCalibration: WARNING unable to calculate aperture for order " << order << endl;
						break;
					}
                }
                
                // Set up background apertures
                operaExtractionAperture<Line> leftBackgroundAperture = ProcessAperture("Left Background", extractionLine, -(widthX + backgroundAperture)/2, backgroundAperture, instrumentProfile, ycenter, normalizationFactor, ipImage, order, args.verbose);
				spectralOrder->setBackgroundApertures(0, new operaExtractionAperture<Line>(leftBackgroundAperture));
                
                operaExtractionAperture<Line> rightBackgroundAperture = ProcessAperture("Right Background", extractionLine, (widthX + backgroundAperture)/2, backgroundAperture, instrumentProfile, ycenter, normalizationFactor, ipImage, order, args.verbose);
				spectralOrder->setBackgroundApertures(1, new operaExtractionAperture<Line>(rightBackgroundAperture));
                if(args.verbose) cout << endl;
                
                // For plotting
                if(fdata.is_open()) {
					if(order - minorder >= 5) {
						ystackIndex++;
						minorder = order;
					}
					float shiftX = 0, shiftY = 0;
					if(ordernumber == NOTPROVIDED) {
						shiftX = (float)((order - minorder)*instrumentProfile->getxsize());
						shiftY = (float)(ystackIndex*instrumentProfile->getysize());
					}
					for (unsigned j=0; j<NYPoints; j++) {	 
						for (unsigned i=0; i<NXPoints; i++) {
							fdata << order << " " 
							<< instrumentProfile->getIPixXCoordinate(i) + shiftX << " " 
							<< instrumentProfile->getIPixYCoordinate(j) + shiftY << " " 
							<< ipImage[j][i]/normalizationFactor << endl;
						}
						fdata << endl;
					}
					fdata << endl;
				}
            } else if(args.verbose) {
				cout << "operaExtractionApertureCalibration:aperture of order " << order << " skipped." << endl;
            }
        }        
        
        // output aperture to file...
        operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputApertureFile,Aperture);
        
        if (ftiltdata1.is_open()) ftiltdata1.close();
        if (ftiltdata2.is_open()) ftiltdata2.close();

        if (!tiltdata1filename.empty() && !tiltdata2filename.empty() && !scriptfilename.empty()) {
            GenerateExtractionApertureTiltPlot(tiltscriptfilename, tiltplotfilename, tiltdata1filename, tiltdata2filename, mediantilt, mediantilterror, interactive);
        }
 
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateExtractionAperturePlot(scriptfilename,plotfilename,datafilename,interactive);
            }
        }
	}
    
	catch (operaException e) {
		cerr << "operaExtractionApertureCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtractionApertureCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

operaExtractionAperture<Line> ProcessAperture(string apertureName, const Line& extractionLine, double apertureXoffset, double apertureWidth, operaInstrumentProfile* instrumentProfile, double yCenter, double normalizationFactor, DMatrix& ipImage, unsigned ordernumber, bool verbose) {
    if(verbose) cout << "operaExtractionApertureCalibration: # order = " << ordernumber << " " << apertureName << " ";
    double apertureXcenter = extractionLine.getMidPoint().getXcoord() + apertureXoffset;
	double apertureYcenter = extractionLine.getYcoord(apertureXcenter);
    operaExtractionAperture<Line> lineAperture = CreateLineAperture(operaPoint(apertureXcenter, apertureYcenter), extractionLine.getSlope(), extractionLine.getWidth(), apertureWidth, instrumentProfile, yCenter, normalizationFactor);
    if (verbose) cout << "Flux Fraction = " << lineAperture.getFluxFraction() << " # out-of-bound points = " << 0 << endl;
	const PixelSet *aperturePixels = lineAperture.getSubpixels();
	for(unsigned i=0; i<aperturePixels->getNPixels(); i++) ipImage[(unsigned)aperturePixels->getjIndex(i)][(unsigned)aperturePixels->getiIndex(i)] = 0;
	return lineAperture;
}

operaExtractionAperture<Line> CreateLineAperture(operaPoint apertureMidPoint, double apertureSlope, double apertureHeight, double apertureWidth, operaInstrumentProfile* instrumentProfile, double yCenter, double normalizationFactor) {
	Line extractionLine(apertureSlope, apertureHeight, apertureWidth, apertureMidPoint);
	operaExtractionAperture<Line> aperture(&extractionLine, instrumentProfile, yCenter);
	const PixelSet *aperturePixels = aperture.getSubpixels();
	double fluxFraction = 0;
	for(unsigned i=0; i<aperturePixels->getNPixels(); i++) {
		int pixX = aperturePixels->getiIndex(i);
		int pixY = aperturePixels->getjIndex(i);
		if(pixX >= 0 && pixX < (int)instrumentProfile->getNXPoints() && pixY >= 0 && pixY < (int)instrumentProfile->getNYPoints()) {
			fluxFraction += instrumentProfile->getipDataFromPolyModel(yCenter,(unsigned)pixX,(unsigned)pixY);
		} else {
			//This should never happen now -- see operaExtractionAperture::setSubpixels(operaInstrumentProfile*) for details.
			throw operaException("operaExtractionApertureCalibration: out of bound coordinates.",operaErrorInstrumentProfileImproperCoordinateRequest, __FILE__, __FUNCTION__, __LINE__);
		}
	}
	aperture.setFluxFraction(fluxFraction/normalizationFactor);
	return aperture;
}

bool CalculateTilt(operaSpectralOrder *spectralOrder, double width, double height, unsigned nRowSamples, unsigned pickImageRow, unsigned xbin, double& tiltAngle, double& tiltAngleError, double& yintercept, double& fluxFraction, ofstream& ftiltdata1, bool debug) {
	if (!spectralOrder->gethasGeometry() || !spectralOrder->gethasInstrumentProfile()) return false;
	
	operaGeometry *Geometry = spectralOrder->getGeometry();
	operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
	
	unsigned NXPoints = instrumentProfile->getNXPoints();
	unsigned NYPoints = instrumentProfile->getNYPoints();
	
	if (NXPoints < xbin * 4) {
		throw operaException("operaExtractionApertureCalibration: the IP needs at least 4*xbin X points", operaErrorInvalidInput , __FILE__, __FUNCTION__, __LINE__);
	}
	
	double ymiddle = (Geometry->getYmax() + Geometry->getYmin())/2;
	double ystep = (Geometry->getYmax() - Geometry->getYmin())/(double)nRowSamples;
	
	if(pickImageRow) {
		nRowSamples = 1;
	}

	double normalizationFactor = 0;
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			normalizationFactor += instrumentProfile->getipDataFromPolyModel(ymiddle,i,j);
			// for 3D plot
			if(debug) cout << instrumentProfile->getIPixXCoordinate(i) << " " << instrumentProfile->getIPixYCoordinate(j) << " " << instrumentProfile->getipDataFromPolyModel(ymiddle,i,j) << endl;
		}
	}
	
	operaVector xphotocenter, yphotocenter, photosum;
	
	for(unsigned k=0; k<nRowSamples; k++) {
		
		double ycenter = ymiddle - (float)(nRowSamples-1)*ystep/2 + ystep*(float)k + 0.5;
		
		if(pickImageRow) ycenter = pickImageRow + 0.5;
		
		double photosum_tmp = 0;
		double yphotocenter_tmp = 0;
		double xphotocenter_tmp = 0;
		unsigned npoints = 0;
		
		operaVector xcoords;
		operaVector ycoords;
		
		//For tilt plot data file
		operaVector xphotocenterBin, yphotocenterBin;
		
		for (unsigned i=xbin; i<NXPoints-xbin; i++) {
			
			double xcoord = instrumentProfile->getIPixXCoordinate(i);
			xcoords.insert(xcoord);
			
			double localphotosum = 0;
			double ysum = 0;
			for (unsigned j=0; j<NYPoints; j++) {
				double ip = instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
				if(!isnan(ip) && ip>0) {
					localphotosum += ip;
					double ycoord = instrumentProfile->getIPixYCoordinate(j);
					ysum += ycoord*ip;
					yphotocenter_tmp += ycoord*ip;
					npoints++;
				}
			}
			if(localphotosum) {
				ycoords.insert(ysum/localphotosum);
				xphotocenter_tmp += localphotosum*xcoord;
				photosum_tmp += localphotosum;
			}
			
			if(xcoords.size() == xbin) {
				xphotocenter.insert(xphotocenter_tmp/photosum_tmp);
				yphotocenter.insert(yphotocenter_tmp/photosum_tmp);
				photosum.insert(photosum_tmp);
				
				// For plot tilt data file:
				xphotocenterBin.insert(xphotocenter_tmp/photosum_tmp);
				yphotocenterBin.insert(yphotocenter_tmp/photosum_tmp);
				
				if(debug) cout << spectralOrder->getorder() << " " << xphotocenter.last() << " " << yphotocenter.last() << " " <<  StdDev(xcoords) << " " << StdDev(ycoords) << " " << photosum_tmp/npoints << endl;

				xcoords.clear();
				ycoords.clear();
				xphotocenter_tmp = 0;
				yphotocenter_tmp = 0;
				photosum_tmp = 0;
				npoints = 0;
			}
		}
		
		//This debug section is useful to produce a plot that indicates any variation of the tilt angle along the detector rows or along orders.
		if (ftiltdata1.is_open()) {
			double amBin,bmBin,abdevmBin;
			LinearFit(xphotocenterBin, yphotocenterBin, amBin, bmBin, abdevmBin);
			double tiltAngleBin = 180*atan(bmBin)/(M_PI);
			double tiltAngleBinError = abdevmBin/0.674433;
			ftiltdata1 << spectralOrder->getorder() << " " << ycenter << " " << tiltAngleBin << " " << tiltAngleBinError << endl;
		}
	}
	
	double am,bm,abdevm;
	LinearFit(xphotocenter, yphotocenter, am, bm, abdevm); //robust linear fit: f(x) = a + b*x
	
	// This debug section is useful to produce plots containing individual x and y-centroids in the IP as well as the model obtained from measurements above
	if(debug) {
		for (unsigned i=0; i<NXPoints; i++) {
			float ipXcoord = instrumentProfile->getIPixXCoordinate(i);
			float ipsum = 0;
			float ycentroid = 0;
			unsigned nnp = 0;
			for (unsigned j=0; j<NYPoints; j++) {
				float ipYcoord = instrumentProfile->getIPixYCoordinate(j);
				float ip = (float)instrumentProfile->getipDataFromPolyModel(ymiddle,i,j);
				if(ip>0) {
					ipsum += ip;
					ycentroid += ipYcoord*ip;
					nnp++;
				}
			}
			float ipAvg = ipsum/(float)nnp;
			ycentroid /= ipsum;
			float ipYmodel = am + bm*ipXcoord;
			cout << ipXcoord << " " << ipYmodel << " " << ycentroid << " " << ipAvg << endl;
		}
	}
	
	// Calculate the tilt and tilt angle using the values obtained from the fit
	if (!isnan(bm) && bm) {
		tiltAngle = 180*atan(bm)/M_PI;
		tiltAngleError = abdevm/0.674433;
		yintercept = am;
		operaExtractionAperture<Line> lineAperture = CreateLineAperture(operaPoint(0, am), bm, height, width, instrumentProfile, ymiddle, normalizationFactor);
		fluxFraction = lineAperture.getFluxFraction();
		return true;
	}
	return false;
}

void GenerateExtractionAperturePlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "set view 0,0" << endl;
    
    fgnu << "set palette color" << endl;
    fgnu << "set palette gamma 2.5" << endl;
    fgnu << "set pm3d map" << endl;
    fgnu << "unset ztics" << endl;
    fgnu << "set cblabel \"flux fraction\"" << endl;
    fgnu << endl;
	 
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
        if (display) {
            fgnu << "\nset terminal x11" << endl;
            fgnu << "set output" << endl;
            fgnu << "replot" << endl;
        } else {
            fgnu << "\n#set terminal x11" << endl;
            fgnu << "#set output" << endl;
            fgnu << "#replot" << endl;
        }
    } else {
        fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

void GenerateExtractionApertureTiltPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName1, string dataFileName2, float tiltAngle, float tiltAngleError, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());

    fgnu << "reset" << endl;

    fgnu << "set xlabel \"order number\"" << endl;
    fgnu << "set ylabel \"tilt angle (deg)\"" << endl;
    
    fgnu << "set yrange[" << tiltAngle-tiltAngleError*6 << ":" << tiltAngle+tiltAngleError*6 << "]" << endl;
    
    fgnu << "tilt(x) = " << tiltAngle << endl;
    fgnu << "maxtilt(x) = " << tiltAngle+tiltAngleError << endl;
    fgnu << "mintilt(x) = " << tiltAngle-tiltAngleError << endl;
    
    fgnu << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << dataFileName1 << "\" u 1:3:4 t \"individual samples\" w yerr pt 6 ps 1 lw 1 "
        << ",\"" << dataFileName2 << "\" u 1:2 notitle w l lw 2 lt 3"
        << ",\"" << dataFileName2 << "\" u 1:2:3 t \"order median tilt\" w yerr pt 7 lt 3 ps 2 lw 2"
        << ",tilt(x) t \"final median tilt\" w l lt -1 lw 2.5, maxtilt(x) notitle w l lt -1, mintilt(x) notitle w l lt -1" << endl;

        if (display) {
            fgnu << "\nset terminal x11" << endl;
            fgnu << "set output" << endl;
            fgnu << "replot" << endl;
        } else {
            fgnu << "\n#set terminal x11" << endl;
            fgnu << "#set output" << endl;
            fgnu << "#replot" << endl;
        }
    } else {
        fgnu << "\nplot \"" << dataFileName1 << "\" u 1:3:4 t \"individual samples\" w yerr pt 6 ps 1 lw 1 "
        << ",\"" << dataFileName2 << "\" u 1:2 notitle w l lw 2 lt 3"
        << ",\"" << dataFileName2 << "\" u 1:2:3 t \"order median tilt\" w yerr pt 7 lt 3 ps 2 lw 2"
        << ",tilt(x) t \"final median tilt\" w l lt -1 lw 2.5, maxtilt(x) notitle w l lt -1, mintilt(x) notitle w l lt -1" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}
