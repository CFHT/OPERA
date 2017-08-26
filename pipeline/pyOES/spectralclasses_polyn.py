# -*- coding: utf-8 -*-

"""
Spectral Classes
---------------------------
Created on Nov 18 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import os
import gzip
from numpy import *
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pylab
from scipy.interpolate import interp1d
import pandas as pd
from astropy.io import fits

import matplotlib.pyplot as plt
#from numpy.polynomial import polynomial as P
from scipy.optimize import curve_fit
from scipy import stats
########## SPECTRUM CLASS ############
class Spectrum :
    'Common base class for a spectrum'
    
    def __init__(self, Filename):

        """
        Create a Spectrum object.
        
        Parameters
        ----------
        filename : string
            File to read the spectrum from.

        Examples
        --------

        >>> sp = Spectrum("spectrumfile.spc.gz")
        """
        self.filename = Filename
        self.load_opera_spc(self.filename)

    def load_opera_spc(self,Filename):
       self.order,self.wl,self.rvcorr,self.xcorr,self.rawflux,self.rawvar,self.normflux,self.calflux = np.loadtxt(Filename, unpack=True, comments='#', skiprows=11, usecols=(0,4,6,7,8,9,10,12), delimiter=' ')

    def printCleanCalibratedSpectrum (self) :
        for i in range(0,len(self.cleanwl)) :
            print self.cleanwl[i], self.cleannormflux[i], self.cleanrvcorrnormflux[i]

    def getorderspectrum(self, ordernumber, spectrumsize) :
        mask_order = np.where(self.order == ordernumber)
        spectrum = self.rawflux[mask_order]
        
        outspectrum = np.zeros(spectrumsize)

        for i in range(len(outspectrum)) :
            if i < len(spectrum) :
                outspectrum[i] = spectrum[i]      
        return outspectrum
            
    def getorders(self) : 
        orders = []        
        for o in self.order :
            if o not in orders:
                orders.append(o)
        return orders
    
    def makeWavePixfile(self, wavePixfile):
        f = open(wavePixfile, "w")
        last_order = self.order[0]
        pix = 0
        for i in range(len(self.order)):
            if last_order != self.order[i]:
                pix = 0
                last_order = self.order[i]    
            pix += 1
            f.write("%d %f %d\n" % (self.order[i], self.wl[i] * 10,pix))
        f.close()

    def exportToFits(self, originalfitsfile, outputfitsfilename, spectrumsize, wavePixfile) :
        orders = self.getorders()
        
        data = np.zeros((len(orders),spectrumsize))

        for i in range(len(orders)) :
            spectrum = self.getorderspectrum(orders[i],spectrumsize)
            for j in range(len(spectrum)) :
                data[i][j] = spectrum[j]
        
        # Create file containing orders, wavelengths, pixels.
        #if not os.path.exists(wavePixfile):
        """        
        f = open(wavePixfile, "w")
        last_order = self.order[0]
        pix = 0
        for i in range(len(self.order)):
            if last_order != self.order[i]:
                pix = 0
                last_order = self.order[i]    
            pix += 1
            f.write("%d %f %d\n" % (self.order[i], self.wl[i] * 10,pix))
        f.close()
        """
        # Open wavePixfile and add to each variable column from this file 
        self.makeWavePixfile(wavePixfile)
        open_file = open(wavePixfile)
        data1 = genfromtxt(open_file)        
        o = data1[:,0]   
        y = data1[:,1]
        z = data1[:,2]
        o_list = list(set(data1[:,0]))

        # Variable "ap" stands for aperture number and variable "N" stands for numbering for "specN" and "beam" for beam number.
        wat = "wtype=multispec "
        N = 0
        ap = 34
        beam = 0
        for current_o in o_list:
            k = 0
            w1 = [y[j] for j in range(len(y)) if current_o == o[j]   ]
            z1 = [z[j] for j in range(len(z)) if current_o == o[j]   ]
            len_w1 = len(w1) -1
            interval_dw = [(w1[k+1]-w1[k]) for k in range(0,len_w1)]
            dw = np.average(interval_dw)    
            ## Define polynomial fit from numpy ###
            ###Polynomial coeficients             
            #polyn_coefs = poly1d(polyfit(z1,w1,6))
            ###Legendre coeficients
            legendre_coefs = np.polynomial.legendre.legfit(z1, w1, 6)  
            N = N + 1       
            ap = ap + 1
            beam = beam + 1

            ##Produce each "specN =" parameters for 6th order polynomial fit##
            #wat = wat + "spec%d = \"%d %d 0  %f %f 2048 0. 0 0 1. 0. 2 6 0. 2047. %f %f %f %f %f %f\" " % (N, ap, beam, w1[0], dw, polyn_coefs[0], polyn_coefs[1],polyn_coefs[2],polyn_coefs[3],polyn_coefs[4],polyn_coefs[5])

            ##Produce each "specN =" parameters for 6th order Legendre fit##
            wat = wat + "spec%d = \"%d %d 0  %f %f 2048 0. 0 0 1. 0. 2 6 0. 2047. %f %f %f %f %f %f\" " % (N, ap, beam, w1[0], dw, legendre_coefs[0], legendre_coefs[1],legendre_coefs[2], legendre_coefs[3], legendre_coefs[4], legendre_coefs[5])     

        ## Add WAT_0** strings to the header of fits file.
        wat = wat.strip()
        hdu = fits.open(originalfitsfile,ignore_missing_end=True)
        prihdr = hdu[0].header
        prihdr = fits.getheader(originalfitsfile, ignore_missing_end=True)      
        wat_range = [wat[i:i+68] for i in range(0, len(wat), 68)]
        for i in range(0, len(wat_range)):
            prihdr['WAT2_%03d' % (i + 1)] = wat_range[i]

        #Insert to header necessary parameters for fits file.
        prihdr.insert(75, ('WCSDIM', 2))
        prihdr.insert(76, ('CD1_1', 1.))
        prihdr.insert(78, ('LTM1_1', 1.))
        prihdr.insert(79, ('LTM2_2', 1.))
        prihdr.insert(77, ('CD2_2', 1.))
        prihdr.insert(80, ('WAT0_001', 'system=multispec'))
        prihdr.insert(81, ('WAT1_001', 'wtype=multispec label=Wavelength units=angstroms'))
        prihdr.insert(87, ('CCDPROC', 'Feb 22 17:47 CCD processing done')) 
        prihdr.insert(88, ('APSCATTE', 'Scattered ligth substracted'))
        prihdr.insert(89, ('BANDID1', 'spectrum - background non, weights none, clean no'))
        prihdr.append(('DCLOG1', 'REFSPEC1 = dumtd'))
        prihdr.append(('CTYPE1', 'MULTISPE', "Pixel coordinate system"))
        prihdr.append(('CTYPE2', 'MULTISPE', "Pixel coordinate system"))
        prihdr.append(('CDELT1', 1.0, "Binning factor along x"))
        prihdr.append(('CDELT2', 1.0, "Binning factor alog y"))
        del prihdr["GAINM"]
        del prihdr["DATE"] 
        del prihdr["CHECKSUM"]
        del prihdr["DATASUM"]       
        fits.writeto(outputfitsfilename, data, prihdr, clobber=True)
        #print prihdr
