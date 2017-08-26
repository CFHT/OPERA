/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaCreateProduct
 Version: 1.0
 Description: Bundle .s and .sn files into a product
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

#include "libraries/gzstream.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operastringstream.h"		// for Double, Float
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaIOFormats.h"

/*! \file operaCreateProduct.cpp */

using namespace std;

/*!
 * operaCreateProduct
 * \author Doug Teeple
 * \brief Bundle files into an i.fits, p.fits, m.fits Product.
 * \arg argc
 * \arg argv
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

template <typename T>
class DataRow {
private:
	std::vector<T> data;
public:
	unsigned Size() const { return data.size(); }
	void Insert(const T& element) { data.push_back(element); }
	const T& operator[](unsigned index) const { return data[index]; }
};

template <typename T>
class DataMatrix {
private:
	std::vector<DataRow<T> > data;
public:
	unsigned Rows() const { return data.size(); }
	unsigned Cols() const { return data.empty() ? 0 : data[0].Size(); }
	void Insert(const DataRow<T>& datarow) { data.push_back(datarow); }
	const DataRow<T>& operator[](unsigned index) const { return data[index]; }
};

// Reads in a table from a file into matrix, skipping the first skiplines of the file.
void GetMatrixFromDataFile(string filename, DataMatrix<float>& matrix, unsigned skiplines);

// Updates the FITS product to contain the values in matrix. Starts at coloffset.
// Product must have at least as many rows as matrix and at least coloffset more columns than matrix.
void UpdateProductFromMatrix(operaFITSProduct& Product, const DataMatrix<float>& matrix, unsigned coloffset = 0);

void SetHeaderColumnsLE(operaFITSProduct& Product, operaSpectralOrder_t spectralOrderType);

void SetHeaderColumnsExtended(operaFITSProduct& Product, operaSpectralOrder_t spectralOrderType, unsigned cols);

/*!
 * \verbatim
 *	   upena-compatible headers for column names in extension 1
 *	 
 *     (1) Spectroscopy Star only mode
 *         First column = wavelength in nanometres
 *         Second column = intensity
 *         Third column = error bar
 *     
 *     (2) Polarimetry
 *         1st col = wavelength in nanometres
 *         2d  col = intensity
 *         3rd col = polarisation (Q or U or V or W)
 *         4th col = Check Spectra #1
 *         5th col = Check Spectra #2
 *         6th col = error bar
 *     
 *     (3) Spectroscopy Star + Sky
 *         1st col = wavelength
 *         2d  col = star spectra (sky subtracted)
 *         3rd col = star + sky spectra
 *         4th col = sky spectra
 *         5, 6, 7 = error bars for each column 2, 3, 4
 *
 * The headers spectoscopy star-only look like this:
 *
 * REDUCTIO= 'Intensity'          / Type of reduction                              
 * NORMAL  =                    2 / Normalized and Un-normalized Data              
 * COMMENT  File contains automatic wavelength correction and uncorrected data.    
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Intensity'          / Normalized                                     
 * COL3    = 'ErrorBar'           / Normalized                                     
 * COL4    = 'Wavelength'         / UnNormalized                                   
 * COL5    = 'Intensity'          / UnNormalized                                   
 * COL6    = 'ErrorBar'           / UnNormalized                                   
 * COL7    = 'Wavelength'         / Normalized, no autowave correction             
 * COL8    = 'Intensity'          / Normalized, no autowave correction             
 * COL9    = 'ErrorBar'           / Normalized, no autowave correction             
 * COL10   = 'Wavelength'         / UnNormalized, no autowave correction           
 * COL11   = 'Intensity'          / UnNormalized, no autowave correction           
 * COL12   = 'ErrorBar'           / UnNormalized, no autowave correction  
 *
 * or, for polar:
 *
 * REDUCTIO= 'Polar   '           / Type of reduction                              
 * NORMAL  =                    2 / Normalized and Un-normalized Data              
 * COMMENT  File contains automatic wavelength correction and uncorrected data.    
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Intensity'          / Normalized                                     
 * COL3    = 'Stokes  '           / Normalized                                     
 * COL4    = 'CheckN1 '           / Normalized                                     
 * COL5    = 'CheckN2 '           / Normalized                                     
 * COL6    = 'ErrorBar'           / Normalized                                     
 * COL7    = 'Wavelength'         / UnNormalized                                   
 * COL8    = 'Intensity'          / UnNormalized                                   
 * COL9    = 'Stokes  '           / UnNormalized                                   
 * COL10   = 'CheckN1 '           / UnNormalized                                   
 * COL11   = 'CheckN2 '           / UnNormalized                                   
 * COL12   = 'ErrorBar'           / UnNormalized                                   
 * COL13   = 'Distance'           / Normalized, no autowave correction             
 * COL14   = 'Intensity'          / Normalized, no autowave correction             
 * COL15   = 'Stokes  '           / Normalized, no autowave correction             
 * COL16   = 'CheckN1 '           / Normalized, no autowave correction             
 * COL17   = 'CheckN2 '           / Normalized, no autowave correction             
 * COL18   = 'ErrorBar'           / Normalized, no autowave correction             
 * COL19   = 'Distance'           / UnNormalized, no autowave correction           
 * COL20   = 'Intensity'          / UnNormalized, no autowave correction           
 * COL21   = 'Stokes  '           / UnNormalized, no autowave correction           
 * COL22   = 'CheckN1 '           / UnNormalized, no autowave correction           
 * COL23   = 'CheckN2 '           / UnNormalized, no autowave correction           
 * COL24   = 'ErrorBar'           / UnNormalized, no autowave correction           
 * COMMENT For Stokes Q, V, and W, keep the Stokes parameter sign as is            
 * COMMENT For Stokes U, invert the sign of the Stokes parameter 
 *
 * or for spectoscopy, star+sky"
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Star    '           / Normalized                                     
 * COL3    = 'Star+sky'           / Normalized                                     
 * COL4    = 'Sky     '           / Normalized                                     
 * COL5    = 'ErrorBar1'          / Normalized                                     
 * COL6    = 'ErrorBar2'          / Normalized                                     
 * COL7    = 'ErrorBar3'          / Normalized                                     
 * COL8    = 'Wavelength'         / UnNormalized                                   
 * COL9    = 'Star    '           / UnNormalized                                   
 * COL10   = 'Star+sky'           / UnNormalized                                   
 * COL11   = 'Sky     '           / UnNormalized                                   
 * COL12   = 'ErrorBar1'          / UnNormalized                                   
 * COL13   = 'ErrorBar2'          / UnNormalized                                   
 * COL14   = 'ErrorBar3'          / UnNormalized                                   
 * COL15   = 'Distance'           / Normalized, no autowave correction             
 * COL16   = 'Star    '           / Normalized, no autowave correction             
 * COL17   = 'Star+sky'           / Normalized, no autowave correction             
 * COL18   = 'Sky     '           / Normalized, no autowave correction             
 * COL19   = 'ErrorBar1'          / Normalized, no autowave correction             
 * COL20   = 'ErrorBar2'          / Normalized, no autowave correction             
 * COL21   = 'ErrorBar3'          / Normalized, no autowave correction             
 * COL22   = 'Distance'           / UnNormalized, no autowave correction           
 * COL23   = 'Star    '           / UnNormalized, no autowave correction           
 * COL24   = 'Star+sky'           / UnNormalized, no autowave correction           
 * COL25   = 'Sky     '           / UnNormalized, no autowave correction           
 * COL26   = 'ErrorBar1'          / UnNormalized, no autowave correction           
 * COL27   = 'ErrorBar2'          / UnNormalized, no autowave correction           
 * COL28   = 'ErrorBar3'          / UnNormalized, no autowave correction           
 * 
 * 
 * \endverbatim
 */

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string version, date;
	string inputfilename;
	string outputfilename;
	string ufile, nfile, uwfile, nwfile;
	string spectrumfile;
	string snrfilename;
	string rvelfilename;
	string tellfilename;
	string sresfilename;
	string object;
	unsigned spectralOrderType_val = LibreEspritsp2Spectrum;
	int compressionVal;
	args.AddOptionalArgument("version", version, "", "");
	args.AddOptionalArgument("date", date, "", "");
	args.AddRequiredArgument("input", inputfilename, "input file (o.fits)");
	args.AddRequiredArgument("output", outputfilename, "output file (i.fits/m.fits)");
	args.AddOptionalArgument("ufile", ufile, "", "unnormalized spectrum (iu.s/pu.s)");
	args.AddOptionalArgument("nfile", nfile, "", "normalized spectrum (in.s/pn.s)");
	args.AddOptionalArgument("uwfile", uwfile, "", "unnormalized wavelength corrected spectrum (iuw.s/puw.s)");
	args.AddOptionalArgument("nwfile", nwfile, "", "normalized wavelength corrected spectrum (inw.s/pnw.s)");
	args.AddOptionalArgument("spectrumfile", spectrumfile, "", "extended spectrum (.spc)");
	args.AddRequiredArgument("spectrumtype", spectralOrderType_val, "spectral order type");
	args.AddOptionalArgument("snr", snrfilename, "", ".sn");
	args.AddOptionalArgument("rvel", rvelfilename, "", "i.rvel");
	args.AddOptionalArgument("tell", tellfilename, "", "i.tell");
	args.AddOptionalArgument("sres", sresfilename, "", ".sres");
	args.AddOptionalArgument("object", object, "", "object name, needed for Libre-Esprit output");
	args.AddOptionalArgument("compressiontype", compressionVal, cNone, "compression type");
		
	try {
		args.Parse(argc, argv);
		
		operaSpectralOrder_t spectralOrderType = operaSpectralOrder_t(spectralOrderType_val);
		
		if (inputfilename.empty()) // we need an input...
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		if (outputfilename.empty()) // we need an output...
			throw operaException("operaCreateProduct: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		if ((ufile.empty() || nfile.empty() || uwfile.empty() || nwfile.empty()) && spectrumfile.empty()) // we need either a spectrum or all four LE files...
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		if ((!ufile.empty() || !nfile.empty() || !uwfile.empty() || !nwfile.empty()) && !spectrumfile.empty()) // ...but we don't wont both
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		
		eCompression compression = (eCompression)compressionVal;
		if (args.verbose) {
			cout << "operaCreateProduct: input= " << inputfilename << endl;
			cout << "operaCreateProduct: output= " << outputfilename << endl;
			cout << "operaCreateProduct: ufile= " << ufile << endl;
			cout << "operaCreateProduct: nfile= " << nfile << endl;
			cout << "operaCreateProduct: uwfile= " << uwfile << endl;
			cout << "operaCreateProduct: nwfile= " << nwfile << endl;
			cout << "operaCreateProduct: snrfilename= " << snrfilename << endl;
			cout << "operaCreateProduct: rvelfilename= " << rvelfilename << endl;
			cout << "operaCreateProduct: tellfilename= " << tellfilename << endl;
			cout << "operaCreateProduct: OPERA version= " << version << endl;
			cout << "operaCreateProduct: Reduction date= " << date << endl;
			cout << "operaCreateProduct: compression= " << compression << endl;
			cout << "operaCreateProduct: spectrumtype= " << spectralOrderType << endl;
		}
		
		operaFITSProduct Product(outputfilename, 0, 0, compression);
		Product.CopyHeader(inputfilename);
		
		Product.operaFITSDeleteHeaderKey("DATASEC");
		Product.operaFITSDeleteHeaderKey("DETSEC");
		Product.operaFITSAddComment("----------------------------------------------------");
		Product.operaFITSAddComment("| Processed by the CFHT OPERA Open Source Pipeline |");
		Product.operaFITSAddComment("----------------------------------------------------");
		Product.operaFITSAddComment(version);
		Product.operaFITSAddComment("Processing Date");
		Product.operaFITSAddComment("---------------");
		Product.operaFITSAddComment(date);
		Product.operaFITSAddComment("------------------------------------------------------------------------");
		Product.operaFITSAddComment(itos(spectralOrderType));
		
		if (!snrfilename.empty()) {
			FormatData snrdata;
			operaIOFormats::ReadCustomFormat("SNR", snrdata, snrfilename);
			unsigned ordercount = snrdata.extract<unsigned>();
			for(unsigned i=0; i < ordercount; i++) {
				string name = "SNR" + snrdata.extract<string>();
				snrdata.skip(1);
				double snrspec, snrccd;
				snrdata >> snrspec >> snrccd;
				ostringstream val;
				val << snrspec << " / " <<  snrccd;
				Product.operaFITSSetHeaderValue(name, val.str(), "snr per spectral / ccd bin");
			}
		}
		if (!sresfilename.empty()) {
			FormatData sresdata;
			operaIOFormats::ReadCustomFormat("spectralres", sresdata, sresfilename);
			unsigned ordercount = sresdata.extract<unsigned>();
			for(unsigned i=0; i < ordercount; i++) {
				string orderstr, sresval, sresvar;
				sresdata >> orderstr >> sresval >> sresvar;
				Product.operaFITSSetHeaderValue("SPCRES"+orderstr, sresval + " +/- " + sresvar, "spectral resolution  +/- dispersion");
				Product.operaFITSSetHeaderValue("RVPREC"+orderstr, sresdata.extract<double>(), "radial velocity precision (m/s)");
				Product.operaFITSSetHeaderValue("NLINES"+orderstr, sresdata.extract<unsigned short>(), "spectral lines used");
			}
		}
		if (!rvelfilename.empty()) {
			FormatData rveldata;
			operaIOFormats::ReadCustomFormat("rvel", rveldata, rvelfilename);
			Product.operaFITSSetHeaderValue("HRV", rveldata.extract<double>(), "Heliocentric RV correction (km/s)");
			Product.operaFITSSetHeaderValue("HRVLUNAR", rveldata.extract<double>(), "lunar component of HRV correction (km/s)");
			Product.operaFITSSetHeaderValue("HRVORBIT", rveldata.extract<double>(), "orbital component of HRV correction (km/s)");
			Product.operaFITSSetHeaderValue("HRVDIURN", rveldata.extract<double>(), "diurnal component of HRV correction (km/s)");
			Product.operaFITSSetHeaderValue("HJDUTC", rveldata.extract<double>(), "Heliocentric Julian date (UTC) mid-exposure");
			Product.operaFITSSetHeaderValue("HJDTT", rveldata.extract<double>(), "Heliocentric Julian date (TT) mid-exposure");
		}
		if (!tellfilename.empty()) {
			FormatData telldata;
			operaIOFormats::ReadCustomFormat("tell", telldata, tellfilename);
			Product.operaFITSSetHeaderValue("TELLRV", telldata.extract<double>(), "telluric RV correction (km/s)");
			Product.operaFITSSetHeaderValue("TELLERR", telldata.extract<double>(), "telluric RV correction error (km/s)");
		}
		
		if (!ufile.empty() && !nfile.empty() && ! uwfile.empty() && !nwfile.empty()) {
			string inputfiles[4] = {nfile, ufile, nwfile, uwfile};
			for(unsigned i = 0; i < 4; i++) {
				DataMatrix<float> readdata;
				GetMatrixFromDataFile(inputfiles[i], readdata, 2);
				if(i == 0) {
					Product.resize(readdata.Rows(), readdata.Cols()*4);
					SetHeaderColumnsLE(Product, spectralOrderType);
				}
				UpdateProductFromMatrix(Product, readdata, readdata.Cols()*i);
			}
        }
        else if (!spectrumfile.empty()) {
			DataMatrix<float> readdata;
			GetMatrixFromDataFile(spectrumfile, readdata, 1);
			Product.resize(readdata.Rows(), readdata.Cols());
			SetHeaderColumnsExtended(Product, spectralOrderType, readdata.Cols());
			UpdateProductFromMatrix(Product, readdata);
		}
		
		Product.operaFITSImageSave();
		Product.operaFITSImageClose();
		
		if (args.verbose && spectralOrderType == LibreEspritpolarimetry) cout << "operaCreateProduct: done polarimetry " << endl;
		else if (args.verbose) cout << "operaCreateProduct: done intensity " << endl;
    }
    catch (operaException e) {
        cerr << "operaCreateProduct: " << e.getFormattedMessage() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "operaCreateProduct: " << operaStrError(errno) << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

void GetMatrixFromDataFile(const string filename, DataMatrix<float>& matrix, const unsigned skiplines) {
	operaistream fin(filename.c_str());
	if (fin.is_open()) {
		string dataline;
		unsigned line = 0;
		while (getline(fin, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') {
				if (line >= skiplines) {
					istringstream ss (dataline);
					DataRow<float> datarow;
					for (Float NanTolerantFloat = 0.0; ss >> NanTolerantFloat; datarow.Insert(NanTolerantFloat.f));
					matrix.Insert(datarow);
				}
				line++;
			}
		}
		fin.close();
	}
}

void UpdateProductFromMatrix(operaFITSProduct& Product, const DataMatrix<float>& matrix, const unsigned coloffset) {
	for (unsigned row = 0; row < matrix.Rows(); row++) {
		for (unsigned col = 0; col < matrix.Cols(); col++) {
			Product[col+coloffset][row] = matrix[row][col];
		}
	}
}

const string Wavelength = "Wavelength";
const string Intensity = "Intensity";
const string Polar = "Polar";
const string Star = "Star";
const string Starplussky = "Star+sky";
const string Sky = "Sky";
const string ErrorBar = "ErrorBar";
const string ErrorBar1 = "ErrorBar1";
const string ErrorBar2 = "ErrorBar2";
const string ErrorBar3 = "ErrorBar3";
const string Stokes = "Stokes";
const string CheckN1 = "CheckN1";
const string CheckN2 = "CheckN2";

void SetHeaderColumnsLE(operaFITSProduct& Product, operaSpectralOrder_t spectralOrderType) {
	vector<string> colnames;
	switch (spectralOrderType) {
		case LibreEspritsp2Spectrum:
		case LibreEspritpolSpectrum:
			Product.operaFITSSetHeaderValue("REDUCTIO", Intensity, "Type of reduction");                               
			Product.operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");  
			Product.operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			
			colnames.push_back(Wavelength);
			colnames.push_back(Intensity);
			colnames.push_back(ErrorBar);
			break;
		case LibreEspritsp1Spectrum:
			Product.operaFITSSetHeaderValue("REDUCTIO", Intensity, "Type of reduction");                               
			Product.operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");                               
			Product.operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			
			colnames.push_back(Wavelength);
			colnames.push_back(Star);
			colnames.push_back(Starplussky);
			colnames.push_back(Sky);
			colnames.push_back(ErrorBar1);
			colnames.push_back(ErrorBar2);
			colnames.push_back(ErrorBar3);
			break;
		case LibreEspritpolarimetry:
			Product.operaFITSSetHeaderValue("REDUCTIO", Polar, "Type of reduction");                               
			Product.operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");                               
			Product.operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			Product.operaFITSAddComment("For Stokes Q, V, and W, keep the Stokes parameter sign as is");
			Product.operaFITSAddComment("For Stokes U, invert the sign of the Stokes parameter");
			
			colnames.push_back(Wavelength);
			colnames.push_back(Intensity);
			colnames.push_back(Stokes);
			colnames.push_back(CheckN1);
			colnames.push_back(CheckN2);
			colnames.push_back(ErrorBar);
			break;
		default:
			throw operaException("operaCreateProduct: ", operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
			break;
	}
	
	vector<string> colgroups;
	colgroups.push_back("Normalized");
	colgroups.push_back("UnNormalized");
	colgroups.push_back("Normalized, no autowave correction");
	colgroups.push_back("UnNormalized, no autowave correction");
	
	for(unsigned i = 0; i < colgroups.size(); i++) {
		for(unsigned j = 0; j < colnames.size(); j++) {
			Product.AddColumnToHeader(colnames[j], colgroups[i]);
		}
	}
}

void SetHeaderColumnsExtended(operaFITSProduct& Product, operaSpectralOrder_t spectralOrderType, unsigned cols) {
	switch (spectralOrderType) {
		case LibreEspritsp2Spectrum:
		case LibreEspritpolSpectrum:
		case LibreEspritsp1Spectrum:
			Product.operaFITSSetHeaderValue("REDUCTIO", Intensity, "Type of reduction");                               
			Product.operaFITSAddComment("File contains data used to create Libre-Esprit compatible spectra.");
			
			Product.AddColumnToHeader("Order", "Order number");
			Product.AddColumnToHeader("NElements", "Number of elements in order");
			Product.AddColumnToHeader("NBeams", "Number of beams");
			Product.AddColumnToHeader("ElementIndex", "Index of this element");
			Product.AddColumnToHeader("Wavelength", "Uncorrected wavelength (nm)");
			Product.AddColumnToHeader("Tell", "Wavelength w/ telluric correction applied (nm)");
			Product.AddColumnToHeader("RVel", "Heliocentric wavelength correction (nm)");
			Product.AddColumnToHeader("XCorr", "Cross-correlation");
			Product.AddColumnToHeader("RawFlux", "Raw flux (electron)");
			Product.AddColumnToHeader("RawFluxVar", "Raw flux variance (electron)");
			Product.AddColumnToHeader("NormalizedFlux", "Normalized flux");
			Product.AddColumnToHeader("NormalizedFluxVar", "Normalized flux variance");
			Product.AddColumnToHeader("FcalFlux", "UnNormalized flux (electron)");
			Product.AddColumnToHeader("FcalFluxVar", "UnNormalized flux variance (electron)");
			while (Product.HeaderColumnCount() < cols) {
				Product.AddColumnToHeader("Beam", "Beam number (starting from 0)");
				Product.AddColumnToHeader("BeamRawFlux", "Raw beam flux for the beam (electron)");
				Product.AddColumnToHeader("BeamRawFluxVar", "Raw beam flux variance (electron)");
				Product.AddColumnToHeader("BeamNormalizedFlux", "Normalized beam flux");
				Product.AddColumnToHeader("BeamNormalizedFluxVar", "Normalized beam flux variance");
				Product.AddColumnToHeader("BeamFcalFlux", "UnNormalized beam flux (electron)");
				Product.AddColumnToHeader("BeamFcalFluxVar", "UnNormalized beam flux variance (electron)");
			}
			break;
		case LibreEspritpolarimetry:
			Product.operaFITSSetHeaderValue("REDUCTIO", Polar, "Type of reduction");                               
			Product.operaFITSAddComment("File contains data used to create Libre-Esprit compatible polarimetry.");
			Product.operaFITSAddComment("For Stokes Q, V, and W, keep the Stokes parameter sign as is");
			Product.operaFITSAddComment("For Stokes U, invert the sign of the Stokes parameter");
			
			Product.AddColumnToHeader("Order", "Order number");
			Product.AddColumnToHeader("NElements", "Number of elements in order");
			Product.AddColumnToHeader("ElementIndex", "Index of this element");
			Product.AddColumnToHeader("Wavelength", "Uncorrected wavelength (nm)");
			Product.AddColumnToHeader("Tell", "Wavelength w/ telluric correction applied (nm)");
			Product.AddColumnToHeader("RVel", "Heliocentric wavelength correction (nm)");
			Product.AddColumnToHeader("XCorr", "Cross-correlation");
			Product.AddColumnToHeader("StokesI", "Raw Stokes I");
			Product.AddColumnToHeader("StokesIVar", "Raw Stokes I variance");
			Product.AddColumnToHeader("NormalizedStokesI", "Normalized Stokes I");
			Product.AddColumnToHeader("NormalizedStokesIVar", "Normalized Stokes I variance");
			Product.AddColumnToHeader("FcalStokesI", "Calibrated Stokes I");
			Product.AddColumnToHeader("FcalStokesIVar", "Calibrated Stokes I variance");
			Product.AddColumnToHeader("DegreeOfPol", "Degree of polarization = StokesParam/Stokes I");
			Product.AddColumnToHeader("DegreeOfPolVar", "Degree of polarization variance");
			Product.AddColumnToHeader("ContPolRemoved", "Degree of polarization w/out continuum");
			Product.AddColumnToHeader("ContPolRemovedVar", "Degree of polarization variance w/out continuum");
			Product.AddColumnToHeader("FirstNullPol", "First null polarization");
			Product.AddColumnToHeader("SecondNullPol", "Second null polarization");
			break;
		default:
			throw operaException("operaCreateProduct: ", operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
			break;
	}
}
