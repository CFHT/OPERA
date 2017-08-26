/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaGeometryCalibration
 Version: 1.0
 Description: Finds the location of orders.
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

#include <fstream>
#include "libraries/operaIOFormats.h"
#include "libraries/operaCCD.h"
#include "libraries/operaFFT.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

#define MINIMUMORDERTOCONSIDER 15

/*! \file operaGeometryCalibration.cpp */

using namespace std;

/*! 
 * operaGeometryCalibration
 * \author Doug Teeple
 * \brief CCD geometry calculations.
 * \arg argc
 * \arg argv
 * \note --outputGeomFile=...
 * \note --masterbias=...
 * \note --masterflat=...
 * \note --badpixelmask=...
 * \note --subformat=...
 * \note --aperture=...
 * \note --detectionMethod=...
 * \note --FFTfilter=...
 * \note --referenceOrderSamplePosition=...
 * \note --minordertouse=...
 * \note --orderOfTracingPolynomial=...
 * \note --colDispersion=...
 * \note --invertOrders=...
 * \note --binsize=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

void GenerateGeometryPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display, unsigned col0, unsigned colf, unsigned row0, unsigned rowf, int AbsOrdNumber[]);

float calculateCenterOfSymmetry(unsigned np, float *ipx, float *ipfunc, float *iperr, unsigned totalNumberOfSlices, bool applyXcenterCorrection);

float calculatePhotoCenter(unsigned np, float *ipx, float *ipfunc, float *iperr, bool applyXcenterCorrection);

unsigned getRowBinnedData(operaFITSImage& flat,unsigned x1,unsigned x2,unsigned nx,unsigned y1,unsigned y2,unsigned ny,float *fx,float *fy,float *yout, bool FFTfilter);

unsigned geometryDetectOrders(unsigned np,float *fx,float *fy,unsigned uslit,float *ipfunc, unsigned binsize, float noise,float gain,float *xmean,float *ymean,float *xmeanerr,int detectionMethod, bool witherrors, bool graces);

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
    
	string outputGeomFile;
	string masterbias;
	string masterflat;
	string badpixelmask;
	string inputGainFile;
	string inputOrderSpacing;
    string subformat_str;
    
    unsigned referenceOrderSamplePosition = 1;  // This position is with respect to the dispersion direction (rows for Espadons)
	unsigned minordertouse = 0;
	unsigned maxorders = MAXORDERS;
	unsigned orderOfTracingPolynomial = 4;
    bool recenterIPUsingSliceSymmetry = TRUE; // This flag will allow the spatial profile to be recentered such as to divide half number of slices on each side
    unsigned totalNumberOfSlices = 6; // Number of slices used in routine to recenter spatial profile
    int detectionMethod = 1; // 1. Gaussian, 2. IP, 3. Top-hat
    bool FFTfilter = false;
    bool graces = false;
    unsigned colDispersion_val;
	int aperture = 20;
	bool invertOrders = true;
	unsigned binsize = 1;
    unsigned nsamples = 3;
    bool witherrors = false;

    string plotfilename;
	string datafilename;
	string scriptfilename;
    bool interactive = false;
    
    args.AddRequiredArgument("outputGeomFile", outputGeomFile, "Geometry output calibration file name");
    args.AddOptionalArgument("masterbias", masterbias, "", "FITS image with masterbias");
    args.AddRequiredArgument("masterflat", masterflat, "FITS image with masterflat");
    args.AddOptionalArgument("badpixelmask", badpixelmask, "", "FITS image with badpixel mask");
    args.AddOptionalArgument("inputGainFile", inputGainFile, "", "Input noise/gain/bias file");
    args.AddRequiredArgument("inputOrderSpacing", inputOrderSpacing, "Order spacing input file name");
    
    args.AddRequiredArgument("subformat", subformat_str, "Image subformat to be inspected \"x1 x2 y1 y2\"");
    args.AddRequiredArgument("referenceOrderSamplePosition", referenceOrderSamplePosition, "Detector position to pick samples, along the dispersion direction (rows for Espadons)");
    args.AddRequiredArgument("minordertouse", minordertouse, "Number of first useful order");
    args.AddRequiredArgument("maxorders", maxorders, "Maximum number of orders to use");
    args.AddRequiredArgument("orderOfTracingPolynomial", orderOfTracingPolynomial, "Degree of polynomial to trace order positions");
    args.AddRequiredArgument("recenterIPUsingSliceSymmetry", recenterIPUsingSliceSymmetry, "Whether to allow the IP to be recentered placing half number of slices on each side");
    args.AddRequiredArgument("totalNumberOfSlices", totalNumberOfSlices, "Total number of slices, used in routine to recenter spatial profile");
    args.AddOptionalArgument("detectionMethod", detectionMethod, 1, "Method for detecting orders: 1 = Gaussian, 2 = IP, 3 = Top-hat");
    args.AddRequiredArgument("FFTfilter", FFTfilter, "Activate Fourier smoothing filter");
    args.AddOptionalArgument("graces", graces, false, "Use alternative method that gives improved results with GRACES");
    args.AddRequiredArgument("colDispersion", colDispersion_val, "Define dispersion direction: 1 = up; 2 = down");
    args.AddRequiredArgument("aperture", aperture, "Aperture size in pixel units");
    args.AddRequiredArgument("invertOrders", invertOrders, "Select this option to invert the counting of order numbers");
    args.AddRequiredArgument("binsize", binsize, "Number of rows to bin in order to detect order positions");
    args.AddRequiredArgument("nsamples", nsamples, "Number of row samples for detecting orders");
    args.AddSwitch("witherrors", witherrors, "Use error bars for polynomial fit");
    
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
	
	try {
        args.Parse(argc, argv);
        
        dispersiondirection_t colDispersion = dispersiondirection_t(colDispersion_val);
        
		struct subformat {
			unsigned x1, x2;
			unsigned y1, y2;
		} subformat = {8, 2040, 3, 4600 };
        SplitStringIntoVals(subformat_str, subformat.x1, subformat.x2, subformat.y1, subformat.y2);
        
		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaGeometryCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need some place to put data...
		if (outputGeomFile.empty()) {
			throw operaException("operaGeometryCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input order spacing calibration file
		if (inputOrderSpacing.empty()) {
			throw operaException("operaGeometryCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (args.verbose) {
			cout << "operaGeometryCalibration: outputGeomFile = " << outputGeomFile << endl;
			cout << "operaGeometryCalibration: masterbias = " << masterbias << endl;
			cout << "operaGeometryCalibration: masterflat = " << masterflat << endl;
			cout << "operaGeometryCalibration: badpixelmask = " << badpixelmask << endl;
			cout << "operaGeometryCalibration: inputGainFile = " << inputGainFile << endl;
			cout << "operaGeometryCalibration: inputOrderSpacing = " << inputOrderSpacing << endl;
			cout << "operaGeometryCalibration: subformat = " << subformat.x1 << " " << subformat.x2 << " "<< subformat.y1  << " "<< subformat.y2 << "\n";
			cout << "operaGeometryCalibration: referenceOrderSamplePosition = " << referenceOrderSamplePosition << endl;
			cout << "operaGeometryCalibration: minordertouse = " << minordertouse << endl;
			cout << "operaGeometryCalibration: maxorders = " << maxorders << endl;
			cout << "operaGeometryCalibration: orderOfTracingPolynomial = " << orderOfTracingPolynomial << endl;
			cout << "operaGeometryCalibration: recenterIPUsingSliceSymmetry = " << recenterIPUsingSliceSymmetry << endl;
			cout << "operaGeometryCalibration: totalNumberOfSlices = " << totalNumberOfSlices << endl;
			cout << "operaGeometryCalibration: detectionMethod = " << detectionMethod << endl;
			cout << "operaGeometryCalibration: FFTfilter = " << FFTfilter << endl;
			cout << "operaGeometryCalibration: graces = " << graces << endl;
			cout << "operaGeometryCalibration: colDispersion = " << colDispersion << endl;
			cout << "operaGeometryCalibration: aperture = " << aperture << endl;
			cout << "operaGeometryCalibration: invertOrders = " << invertOrders << endl;
			cout << "operaGeometryCalibration: binsize = " << binsize << endl;
			cout << "operaGeometryCalibration: witherrors = " << witherrors << endl;
			cout << "operaGeometryCalibration: nsamples = " << nsamples << endl;
            if(args.plot) {
                cout << "operaGeometryCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaGeometryCalibration: datafilename = " << datafilename << endl;
                cout << "operaGeometryCalibration: scriptfilename = " << scriptfilename << endl;
                cout << "operaGeometryCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
            }
		}
		
        ofstream fdata;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        
		// Open input images and load data into an in-memory image
		unsigned x1 = subformat.x1;
		unsigned y1 = subformat.y1;
		unsigned nx = subformat.x2 - subformat.x1;
		unsigned ny = subformat.y2 - subformat.y1;
		if (args.verbose) {
			cout << "operaGeometryCalibration: x1,y1,nx,ny = " << x1 << ' ' << y1 << ' ' << nx  << ' ' << ny << '\n';
		}
        
        operaFITSImage flat(masterflat, tfloat, READONLY);
		if (!masterbias.empty()){
			operaFITSImage bias(masterbias, tfloat, READONLY);
			flat -= bias; // remove bias from masterflat
        }
        
        // Open badpixel mask
        operaFITSImage *badpix = NULL;
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *badpix = 1.0;
        }
        		
        // Setup slit aperture, and create spectral order vector object.
        double slit = (double)aperture;
        unsigned uslit = (unsigned)aperture;
        
		operaSpectralOrderVector spectralOrders(MAXORDERS, ny, ny, 0);
        
        // Read order spacing polynomial
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputOrderSpacing);
        unsigned npars = spectralOrders.getOrderSpacingPolynomial()->getOrderOfPolynomial();
        double *par = spectralOrders.getOrderSpacingPolynomial()->getVector();
        
		// Read gain and noise
		double gain = 1;
		double noise = 5;
        if (!inputGainFile.empty()) {
			unsigned amp = 0;
            operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputGainFile);
            gain = spectralOrders.getGainBiasNoise()->getGain(amp);
            noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		}
		if (args.verbose) cout << "operaGeometryCalibration: gain="<< gain << " noise=" << noise << endl;
        
        // Set up Geometry variables
        for (unsigned order=0; order<MAXORDERS; order++) {
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setapertureWidth(aperture);
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setdispersionDirection(colDispersion);
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setNumberofPointsToBinInYDirection(binsize);
            spectralOrders.GetSpectralOrder(order)->getGeometry()->setYmin(subformat.y1);
            spectralOrders.GetSpectralOrder(order)->getGeometry()->setYmax(subformat.y2);
		}
        
        // Measure instrument spatial profile from several samples of rows
        float *ipfunc = new float[uslit];
        float *ipx = new float[uslit];
        float *iperr = new float[uslit];
        spectralOrders.measureIPAlongRowsFromSamples(flat, *badpix, slit, nsamples, FFTfilter,(float)gain, (float)noise, subformat.x1, subformat.x2, subformat.y1, subformat.y2, ipfunc, ipx, iperr);
        
        // Apply centering corrections to the measured spatial profile
        bool applyXcenterCorrection = true;
        calculatePhotoCenter(uslit, ipx, ipfunc, iperr, applyXcenterCorrection);
        if(recenterIPUsingSliceSymmetry) {
            calculateCenterOfSymmetry(uslit, ipx, ipfunc, iperr, totalNumberOfSlices,applyXcenterCorrection);
#ifdef PRINT_DEBUG
            for(unsigned i=0;i<uslit;i++) {
                cout << ((float)i - slit/2) <<  " " << ipx[i] <<  " " << ipfunc[i] << " " << iperr[i] << endl;
            }
#endif
        }
#ifdef PRINT_DEBUG
        for(unsigned i=0;i<uslit;i++) {
            cout << ((float)i - slit/2) <<  " " << ipx[i] <<  " " << ipfunc[i] << " " << iperr[i] << endl;
        }
#endif
        
        // Set up binsize, number of samples, and allocate memory
		unsigned NumberOfySamples;
		unsigned NumberofPointsToBinInYDirection = binsize;
        unsigned NumberofPointsToBinInYForReference = 2*NumberofPointsToBinInYDirection;
        float *fxref = new float[nx];
        float *fyref = new float[ny];
        float xmean[MAXORDERS],xmeanerr[MAXORDERS],ymean[MAXORDERS];
		NumberOfySamples = (unsigned)ceil((float)(ny-y1)/(float)NumberofPointsToBinInYDirection);
        
        // First pass only using rows around reference row to figure out a reference order map.
        // Note: for the reference sample we use twice as many points to obain higher SNR
        unsigned firstY = y1;
        unsigned lastY = ny;
        if(referenceOrderSamplePosition >= NumberofPointsToBinInYForReference/2)
            firstY = referenceOrderSamplePosition - NumberofPointsToBinInYForReference/2;
        if(referenceOrderSamplePosition + NumberofPointsToBinInYForReference/2 < ny)
            lastY = referenceOrderSamplePosition + NumberofPointsToBinInYForReference/2;
        
        float y = 0;
        unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fxref,fyref,&y,FFTfilter);
        
        if(args.debug) for (unsigned i=0;i<np;i++) cout << fxref[i] << " " << fyref[i] << endl;
        
        unsigned nords = geometryDetectOrders(np,fxref,fyref,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);

        float *xref = new float[MAXORDERS];
        float *yref = new float[MAXORDERS];
        float *xreferr = new float[MAXORDERS];
        int *AbsRefOrdNumber = new int[MAXORDERS];
        
        unsigned nrefs = operaCCDDetectOrderMapBasedOnSpacingPolynomial(np,fxref,fyref,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,xref,yref,xreferr,AbsRefOrdNumber,minordertouse,maxorders);
        
        if(args.debug) for (unsigned i=0;i<nrefs;i++) cout << AbsRefOrdNumber[i] << " " << xref[i] << " " << yref[i] << " " << xreferr[i] << endl;
        
        // Allocate memory to save data obtained from samples
        float *fx = new float[nx];
        float *fy = new float[ny];
        
        float *xord_tmp = new float[MAXORDERS];
        float *yord_tmp = new float[MAXORDERS];
        float *xerrord_tmp = new float[MAXORDERS];
        int *AbsOrdNumber_tmp = new int[MAXORDERS];
        
        float *ypos = new float[NumberOfySamples];
        unsigned *newnords = new unsigned[NumberOfySamples];

        float **xords = new float*[NumberOfySamples];
        float **yords = new float*[NumberOfySamples];
        float **xerrords = new float*[NumberOfySamples];
        int **AbsOrdNumbers = new int*[NumberOfySamples];
        
        for(unsigned k=0;k<NumberOfySamples;k++){
            xords[k] = new float[MAXORDERS];
            yords[k] = new float[MAXORDERS];
            xerrords[k] = new float[MAXORDERS];
            AbsOrdNumbers[k] = new int[MAXORDERS];
            newnords[k] = 0;
        }
        
        /*
         * Once first order is figured out, then we save the data for each order. Note that
         * the absolute order number is important to make sure the data are always associated
         * to the correct order number. 
         */
        
        // Figure out which bin contains the reference order
        unsigned kref = (unsigned)round(float(referenceOrderSamplePosition - y1)/(float)NumberofPointsToBinInYDirection);
        
        // Start detecting orders in samples ABOVE reference row
        for (unsigned i=0;i<nrefs;i++) {
            xord_tmp[i] = xref[i];
            yord_tmp[i] = yref[i];
            xerrord_tmp[i] = xreferr[i];
            AbsOrdNumber_tmp[i] = AbsRefOrdNumber[i];
        }
        for(unsigned k=kref;k<NumberOfySamples;k++){
            unsigned firstY = y1 + NumberofPointsToBinInYDirection*(k);
            unsigned lastY =  y1 + NumberofPointsToBinInYDirection*(k+1);
            if(lastY >= ny) break;
            
            unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fx,fy,&ypos[k],FFTfilter);
            unsigned nords = geometryDetectOrders(np,fx,fy,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);
            
            newnords[k] = operaCCDDetectMissingOrdersUsingNearMap(np,fx,fy,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,nrefs,xord_tmp,yord_tmp,AbsOrdNumber_tmp,xords[k],yords[k],xerrords[k],AbsOrdNumbers[k]);
            
            if(args.debug) {
                for (unsigned i=0;i<newnords[k];i++) {
                    cout << xords[k][i] << " " << yords[k][i] << " " << xerrords[k][i]<< endl;
                }
            }

            // Save current sample in the tmp to be used in the next loop around
            for (unsigned i=0;i<nrefs;i++) {
                xord_tmp[i] = xords[k][i];
                yord_tmp[i] = yords[k][i];
                xerrord_tmp[i] = xerrords[k][i];
                AbsOrdNumber_tmp[i] = AbsOrdNumbers[k][i];
            }
        }
        
        // Then detect orders in samples BELOW reference row
        for (unsigned i=0;i<nrefs;i++) {
            xord_tmp[i] = xref[i];
            yord_tmp[i] = yref[i];
            xerrord_tmp[i] = xreferr[i];
            AbsOrdNumber_tmp[i] = AbsRefOrdNumber[i];
        }
        for(unsigned k=kref-1;k>0;k--){
            unsigned firstY = y1 + NumberofPointsToBinInYDirection*(k);
            unsigned lastY =  y1 + NumberofPointsToBinInYDirection*(k+1);
            if(lastY >= ny) break;
            
            unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fx,fy,&ypos[k],FFTfilter);
            unsigned nords = geometryDetectOrders(np,fx,fy,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);
        
            newnords[k] = operaCCDDetectMissingOrdersUsingNearMap(np,fx,fy,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,nrefs,xord_tmp,yord_tmp,AbsOrdNumber_tmp,xords[k],yords[k],xerrords[k],AbsOrdNumbers[k]);
            
            // Save current sample in the tmp to be used in the next loop around
            for (unsigned i=0;i<nrefs;i++) {
                xord_tmp[i] = xords[k][i];
                yord_tmp[i] = yords[k][i];
                xerrord_tmp[i] = xerrords[k][i];
                AbsOrdNumber_tmp[i] = AbsOrdNumbers[k][i];
            }
        }
        //---
        
        
        /*
         * Now put data into plot file and geometry class, so it can fit the tracing polynomials
         */
		for(unsigned k=0;k<NumberOfySamples;k++){
            if (fdata.is_open()) {
                for(unsigned i=0;i<newnords[k];i++) {
                    fdata << k << " " <<  AbsOrdNumbers[k][i] << " " << xords[k][i] << " " << ypos[k] << " " << xerrords[k][i] << " " << yords[k][i] << endl;
                }
            }

			unsigned lastorder = 0;
			unsigned order = 0;
			for (unsigned i=0;i<newnords[k];i++) {
				operaSpectralOrder *SpectralOrder = NULL;
				for (order=lastorder; order<MAXORDERS; order++) {
					if (spectralOrders.GetSpectralOrder(order)->getorder() == (unsigned)AbsOrdNumbers[k][i]) {
						SpectralOrder = spectralOrders.GetSpectralOrder(order);
						lastorder = order+1;
						break;
					}
				}
				// SpectralOrder should always be found, but...
				if (SpectralOrder) {
					SpectralOrder->getGeometry()->addOrderCenterValue(xords[k][i], ypos[k], yords[k][i], xerrords[k][i]);
				}
			}
        }
        //---
        
        /*
         * Now calculate the polynomials for each order
         */
        Polynomial *polynomials[MAXORDERS]; // for plot

        int *LastAbsOrdNumbers = new int[MAXORDERS];
        
		unsigned ordercount = 0;
		for (unsigned i=0;i<maxorders;i++) {
			double chisqr;	// NOTE: The chisqr is not used here, but it is saved in the polynomial itself
            
            LastAbsOrdNumbers[i] = minordertouse + (int)i;
            
			if (spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->getNdatapoints() > orderOfTracingPolynomial) {	// i.e. we did add data to this order
                spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->traceOrder(orderOfTracingPolynomial, chisqr, witherrors);		// i.e. fit polynomial to the data of a single order
				spectralOrders.setMaxorder(LastAbsOrdNumbers[i]);
				spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->sethasGeometry(true);
				
                if (ordercount == 0) {
					spectralOrders.setMinorder(LastAbsOrdNumbers[i]);
				}
                polynomials[ordercount] = spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->getCenterPolynomial();
                ordercount++;
			}
        }
                
		spectralOrders.setCount(ordercount);
        //---
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateGeometryPlot(scriptfilename, plotfilename, datafilename, ordercount, polynomials, interactive,subformat.x1,subformat.x2,subformat.y1,subformat.y2,LastAbsOrdNumbers);
            }
        }
		
        // Write out geometry calibration file
		if (!outputGeomFile.empty()) operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputGeomFile, Geom);
        
        if(badpix) delete badpix;
		flat.operaFITSImageClose();
		
		for(unsigned k=0;k<NumberOfySamples;k++){
            delete[] xords[k];
            delete[] yords[k];
            delete[] xerrords[k];
            delete[] AbsOrdNumbers[k];
        }
        delete[] xords;
        delete[] yords;
        delete[] xerrords;
        delete[] AbsOrdNumbers;
		
        delete[] ipfunc;
        delete[] ipx;
        delete[] iperr;  
        delete[] fxref;
        delete[] fyref;
        delete[] xref;
        delete[] yref;
        delete[] xreferr;
        delete[] AbsRefOrdNumber;
        delete[] fx;
        delete[] fy;
		delete[] xord_tmp;
        delete[] yord_tmp;
        delete[] xerrord_tmp;
        delete[] AbsOrdNumber_tmp;
        delete[] ypos;
        delete[] newnords;       
	}
	catch (operaException e) {
		cerr << "operaGeometryCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaGeometryCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaGeometryCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

void GenerateGeometryPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display, unsigned col0, unsigned colf, unsigned row0, unsigned rowf, int AbsOrdNumber[])
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "\nset xlabel \"image rows (pixels)\"" << endl;
    fgnu << "set ylabel \"image cols (pixels)\"" << endl;
    
    fgnu << "set pointsize 0.5" << endl;
    
    fgnu << "set xrange[" << row0 << ":" << rowf << "]" << endl;
    fgnu << "set yrange[" << col0 << ":" << colf << "]" << endl;
    
    
    for(unsigned k=0;k<npolynomials;k++) {
        fgnu << "set label \"" << AbsOrdNumber[k] << "\" at " << rowf+50 << "," << polynomials[k]->Evaluate((double)rowf) << " rotate by 90 font \"Helvetica,4\"" << endl;
        fgnu << "poly" << k;
        polynomials[k]->printEquation(&fgnu);
    }
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";
        
        for(unsigned k=0;k<npolynomials;k++) fgnu << ", poly" << k << "f(x)";
        fgnu << endl;
        
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
        fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";
        
        for(unsigned k=0;k<npolynomials;k++) fgnu << ", poly" << k << "f(x)";
        fgnu << endl;
        
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
    
float calculateCenterOfSymmetry(unsigned np, float *ipx, float *ipfunc, float *iperr, unsigned totalNumberOfSlices, bool applyXcenterCorrection) {
    float xphotocenter = 0;
    unsigned sliceSize = (unsigned)ceil((float)(np/totalNumberOfSlices));
    float minflux = BIG;
    float minx = 0;
    for(unsigned i=(np-sliceSize)/2; i<(np+sliceSize)/2; i++) {
        if(ipfunc[i] < minflux){
            minflux = ipfunc[i];
            minx = ipx[i];
        }
    }
    if(minx) xphotocenter = minx;
    if(applyXcenterCorrection) for(unsigned i=0; i<np; i++) ipx[i] -= xphotocenter;
    return xphotocenter;
}

float calculatePhotoCenter(unsigned np, float *ipx, float *ipfunc, float *iperr, bool applyXcenterCorrection) {
    float xphotocenter = 0;
    float totalflux = 0;
    for(unsigned i=0; i<np; i++) totalflux += ipfunc[i];
    for(unsigned i=0; i<np; i++) xphotocenter += ipx[i]*ipfunc[i]/totalflux;
    if(applyXcenterCorrection) {
        for(unsigned i=0; i<np; i++) ipx[i] -= xphotocenter;
    }
    return xphotocenter;
}


unsigned getRowBinnedData(operaFITSImage& flat,unsigned x1,unsigned x2,unsigned nx,unsigned y1,unsigned y2,unsigned ny,float *fx,float *fy,float *yout, bool FFTfilter) {
    
    unsigned ybinsize = (unsigned)fabs(y2 - y1);
    
    float *fytmp = new float[ny];
    float *fysample = new float[ybinsize];

    unsigned np = 0;
    
    for (unsigned xx=x1; xx<x2; xx++) {
        unsigned ns=0;
        float ysample=0.0;
        
        for (unsigned yy=y1; yy<y2 ; yy++) {
            fysample[ns++] = flat[yy][xx];
            ysample += (float)yy + 0.5;
        }
        ysample /= (float)ns;
        fx[np] = (float)xx + 0.5;
        
        if(FFTfilter){
            fytmp[np] = operaArrayMedian(ns,fysample);
        } else {
            fy[np] = operaArrayMedian(ns,fysample);
        }
        
#ifdef PRINT_DEBUG
        if(FFTfilter){
            cout << fx[np]  << " " << fytmp[np] << endl;
        } else {
            cout << fx[np]  << " " << fy[np] << endl;
        }
#endif
        *yout = ysample;
        np++;
    }
    
    if(FFTfilter){
        operaFFTLowPass(np,fytmp,fy,0.1);
    }
    
    delete[] fysample;
    delete[] fytmp;
    
    return np;
}

unsigned geometryDetectOrders(unsigned np,float *fx,float *fy,unsigned uslit,float *ipfunc, unsigned binsize, float noise,float gain,float *xmean,float *ymean,float *xmeanerr,int detectionMethod, bool witherrors, bool graces) {
    
    unsigned nords = 0;
    
    double slit = (double)uslit;
    double sigma = slit/4.0;
    double threshold = DETECTTHRESHOLD;
    
    if(detectionMethod == 1) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,(float)noise,(float)gain,threshold,xmean,ymean,xmeanerr);
        } else {
            nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,(float)noise,(float)gain,threshold,xmean,ymean);
        }
    } else if (detectionMethod == 2) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingIP(np,fx,fy,uslit,ipfunc,(float)noise,(float)gain,threshold/2,xmean,ymean,xmeanerr);
        } else {
            if(graces) {
                nords = operaCCDDetectPeaksByXCorrWithIP(np,fx,fy,uslit,ipfunc,(float)noise/sqrt((float)binsize),(float)gain,threshold,xmean,ymean);
            } else {
                /*
                 * The function below does not work on GRACES data. The one above should be better but
                 * it hasn't been tested yet for ESPaDOnS@CFHT. E. Martioli May 14 2014.
                 */
                nords = operaCCDDetectPeaksWithIP(np,fx,fy,uslit,ipfunc,(float)noise,(float)gain,threshold/2,xmean,ymean);
            }
        }
    } else if (detectionMethod == 3) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingTopHat(np,fx,fy,uslit,(float)noise,(float)gain,threshold,xmean,ymean,xmeanerr);
        } else {
            nords = operaCCDDetectPeaksWithTopHat(np,fx,fy,uslit,(float)noise,(float)gain,threshold,xmean,ymean);
        }
    }

#ifdef PRINT_DEBUG
    for(unsigned i=0;i<nords;i++) {
        cout << xmean[i] << " " << y << " " << xmeanerr[i] << " " << ymean[i] << endl;
    }
#endif
    return nords;
}
