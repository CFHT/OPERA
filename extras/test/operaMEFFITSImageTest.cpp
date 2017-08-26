
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaMETFITSImageTest 
 Version: 1.0
 Description: Perform various tests on the operaFITSImage class.
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

#include <getopt.h>
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaWIRCamImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"

/*! \file operaMEFFITSImageTest.cpp */

using namespace std;

/*! 
 * operaMEFFITSImageTest
 * \author Megan Tannock
 * \brief Perform various tests on the operaMEFFITSImageTest class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaMEFFITSImageTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
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
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
			}	// switch
		}	// while
		
#ifdef __APPLE__
		string basepath = "/data/";
#else
		string basepath = "/data/uhane5/opera/";
#endif
		string mefimage = basepath + "/wircam/12AQ03-Mar01/1529925o.fits";

		unsigned XDimension, YDimension, ZDimension, Extensions;
		edatatype Datatype;
		long Npixels;
		
		getFITSImageInformation(mefimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        cout << "operaMEFFITSImageTest: " << mefimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;

		if (Extensions == 0) {
			throw operaException("operaMEFFITSImageTest: "+mefimage+" is not a MEF image. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (ZDimension > 1) {
			throw operaException("operaMEFFITSImageTest: "+mefimage+" is a Cube. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
#if 0		
		string geminiimage = basepath + "/data/gemini/NIFSTestImage.fits.gz";
		getFITSImageInformation(geminiimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        cout << "operaMEFFITSImageTest: " << geminiimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
		operaFITSImage gemini(geminiimage, tfloat, READONLY, cNone);
		cout << "gemini[1024][1024] = " << gemini[1024][1024] << endl;
#endif
		// make a default sky
		if (1) {
			operaWIRCamImage outImage("defaultskyallones.fits", 2048, 2048, 1, 4, tfloat, cNone, false);
			outImage = 1.0;	// all ones
			((operaWIRCamImage &)outImage).operaWIRCamImageSave();
			outImage.operaFITSImageClose();
        }
		/*
		 * Step 1 Basic save/saveas AGRESSIVE READ!!!!
		 */
		
		if (0) {
			operaMultiExtensionFITSImage inImage(mefimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSImage outImage("mefimage.fits", inImage.getXDimension(), inImage.getYDimension(), inImage.getNExtensions(), tfloat, cNone, false);
			outImage.operaMultiExtensionFITSImageCopyHeader(&inImage);    // This line is copying pixels in addition to the header.
	 
			cout << "operaMEFFITSImageTest: inImage npixels= " << inImage.getnpixels() << endl;
			cout << "operaMEFFITSImageTest: outImage npixels= " << outImage.getnpixels() << endl;
			
			cout << "operaMEFFITSImageTest: Copying inImage pixels to outImage..." << endl;
			outImage = inImage;	// copy entire image pixels
			cout << "operaMEFFITSImageTest: Copied pixels to outImage..." << endl;
            
			// print pixel values of centre pixel for each extension
            unsigned lastextension = inImage.getNExtensions();
			for (unsigned extension=1; extension<=lastextension; extension++) {
				cout << "operaMEFFITSImageTest: outImage["<<extension<<"][1024][1024]=" << outImage[extension][1024][1024] << " inImage["<<extension<<"][1024][1024]=" << inImage[extension][1024][1024] << endl;
			}
			
            // Save the images using operaMultiExtensionFITSImageSaveAs
			cout << "operaMEFFITSImageTest: Saving inImage as inImageMEFImageSaveAs.fits..." << endl;
			inImage.operaMultiExtensionFITSImageSaveAs("inImageMEFImageSaveAs.fits"); 
            cout << "operaMEFFITSImageTest: Saving outImage as outImageMEFImageSaveAs.fits..." << endl;
			outImage.operaMultiExtensionFITSImageSaveAs("outImageMEFImageSaveAs.fits");

			cout << "operaMEFFITSImageTest: Subtracting bias from outImage..." << endl;
			outImage -= 7000;
			cout << "operaMEFFITSImageTest: bias subtracted outImage[1][1024][1024]=" << outImage[1][1024][1024] << " inImage[1][1024][1024]=" << inImage[1][1024][1024] << endl;
			
			cout << "operaMEFFITSImageTest: Copying extension 4 from inImage to 3 of outImage..." << endl;
			outImage[3] = inImage[4];	// copy one extension
			cout << "operaMEFFITSImageTest: bias back in outImage[3][1024][1024]=" << outImage[3][1024][1024] << " inImage[4][1024][1024]=" << inImage[4][1024][1024] << endl;
			 
			// Sum the pixels from 4 extensions and take mean (this is the "sum of pixels" algorithm in the wircam guider...)
			outImage[1] = (inImage[1] + inImage[2] + inImage[3] + inImage[4]) / 4.0;
			cout << "operaMEFFITSImageTest: mean outImage[1][1024][1024]=" << outImage[1][1024][1024] << endl;
            // Print centre pixel values of each extension for in and out
			cout << "operaMEFFITSImageTest: outImage[1][1024][1024]=" << outImage[1][1024][1024] << " inImage[1][1024][1024]=" << inImage[1][1024][1024] << " (outImage[1] = mean of inImage)" << endl;
			cout << "operaMEFFITSImageTest: outImage[2][1024][1024]=" << outImage[2][1024][1024] << " inImage[2][1024][1024]=" << inImage[2][1024][1024] << " (outImage[2] = inImage[2] - bias)" << endl;
			cout << "operaMEFFITSImageTest: outImage[3][1024][1024]=" << outImage[3][1024][1024] << " inImage[3][1024][1024]=" << inImage[3][1024][1024] << " (outImage[3] = inImage[4])" << endl;
			cout << "operaMEFFITSImageTest: outImage[4][1024][1024]=" << outImage[4][1024][1024] << " inImage[4][1024][1024]=" << inImage[4][1024][1024] << " (outImage[4] = inImage[4] - bias)" << endl;

			// Now, create a FITS Image from a single extension (extension 3)
			cout << "operaMEFFITSImageTest: Downsizing to a FITS Image..." << endl;
			operaFITSImage inImageDownSize(inImage[3], false, true);  // bitpix incorrect in header, likely a cfitsio issue
			cout << "operaMEFFITSImageTest: inImageDownSize[1024][1024]=" << inImageDownSize[1024][1024] << endl;
			operaFITSImage outImageDownSize(outImage[3], false, true);
			cout << "operaMEFFITSImageTest: outImageDownSize[1024][1024]=" << outImageDownSize[1024][1024] << endl;
			
			cout << "operaMEFFITSImageTest: Saving extension 3 (inImageDownSize) as inImageDownSize.fits...  (operaFITSImageSaveAs)" << endl;
			inImageDownSize.operaFITSImageSaveAs("inImageDownSize.fits");  // bitpix incorrect in header
			cout << "operaMEFFITSImageTest: Saving extension 3 (outImageDownSize) as outImageDownSize.fits...  (operaFITSImageSaveAs)" << endl;
			outImageDownSize.operaFITSImageSaveAs("outImageDownSize.fits");
			
            cout << "operaMEFFITSImageTest: Saving outImage as outImageSaveAs.fits..." << endl;
            outImage.operaMultiExtensionFITSImageSaveAs("outImageSaveAs.fits");

			cout << "operaMEFFITSImageTest: Closing file(s)..." << endl;
			inImageDownSize.operaFITSImageClose();
			outImageDownSize.operaFITSImageClose();
			inImage.operaFITSImageClose();
			outImage.operaFITSImageSave();
			outImage.operaFITSImageClose();
			cout << "operaMEFFITSImageTest: basic save/saveas test complete." << endl;
        }
		
		/*  
		 * Step 2 extension save/saveas -- LAZY READ
		 */
		// REQUIRES TESTING 
		if (0) {
			cout << "operaMEFFITSImageTest: Testing lazy extension read routines." << endl;
			operaMultiExtensionFITSImage in(mefimage, tfloat, READONLY, cNone);
			operaMultiExtensionFITSImage out("mefextensionimage.fits", in.getXDimension(), in.getYDimension(), in.getNExtensions(), tfloat, cNone);
			operaMultiExtensionFITSImage out2("mefextensionimageall.fits", in.getXDimension(), in.getYDimension(), in.getNExtensions(), tfloat, cNone);
			
			cout << "operaMEFFITSImageTest: " << mefimage << " npixels= " << in.getnpixels() << endl;
			cout << "operaMEFFITSImageTest: mefextensionimage.fits npixels= " << out.getnpixels() << endl;

            // Copy headers 
			out.operaMultiExtensionFITSImageCopyHeader(&in);
			out2.operaMultiExtensionFITSImageCopyHeader(&in);
			out = 0.0;
			out2 = in;
			unsigned lastextension = in.getNExtensions();
			for (unsigned extension=1; extension<=lastextension; extension++) {
				cout << "operaMEFFITSImageTest: copying extension "<< extension << endl;					
				out[extension] = in[extension]; // Exercises readExtension, saveExtension inside of []
				cout << "operaMEFFITSImageTest: copied extension " << extension << endl;					
				cout << "operaMEFFITSImageTest: out["<<extension<<"][1024][1024]=" << out[extension][1024][1024] << " in["<<extension<<"][1024][1024]=" << in[extension][1024][1024] << endl;
			}
		// Save file
			cout << "operaMEFFITSImageTest: Closing file(s)..." << endl;
			out.operaMultiExtensionFITSImageSave();
			out2.operaMultiExtensionFITSImageSave();
			in.operaFITSImageClose();
			out.operaFITSImageClose();
			out2.operaFITSImageClose();
			cout << "operaMEFFITSImageTest: lazy extension read routines test complete." << endl;
		}
	} 
	catch (operaException e) {
		cerr << "operaMEFFITSImageTest: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMEFFITSImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

