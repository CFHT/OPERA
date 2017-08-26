#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
"""
Created on Oct 14 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Oct 14, 2014
"""

import Tkinter
import sys,os
import espadonspipeline
import espadons

class pipelineFrontEnd_tk(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.grid()

########### CONFIG FILE ##################
        # Below it initializes the string for espadons pipeline (opera) system directory
        self.operaDirVariable = Tkinter.StringVar()
        self.operaDirVariable.set("/Users/edermartioli/opera-1.0")
        # Below is the button to load CONFIG FILE
        loadConfigButton = Tkinter.Button(self,text="Load Config", command=self.LoadConfigButtonClick)
        loadConfigButton.grid(column=0,row=0)
        # Below is the text box for config file path entry
        self.filepathVariable = Tkinter.StringVar()
        self.filepathVariable.set(self.operaDirVariable.get() + "/pipeline/pyScripts/config.espadonspipeline")
        self.filepathentry = Tkinter.Entry(self,textvariable=self.filepathVariable)
        self.filepathentry.grid(column=1,row=0,columnspan=6,sticky='EW')
        self.filepathentry.bind("<Return>", self.OnPressEnterConfigPathfield)
###########################################

########### DATA DIR ##################
        # Below is the sidetext for DATADIR
        self.sidetextdd = Tkinter.StringVar()
        sidetextlabeldd = Tkinter.Label(self,textvariable=self.sidetextdd, anchor="w",fg="blue",bg="white")
        sidetextlabeldd.grid(column=0,row=1,columnspan=1)
        self.sidetextdd.set("DATA ROOT DIR:")
        # Below is the text box for string entry
        self.ddentryVariable = Tkinter.StringVar()
        self.ddentry = Tkinter.Entry(self,textvariable=self.ddentryVariable)
        self.ddentry.grid(column=1,row=1,columnspan=6,sticky='EW')
        self.ddentry.bind("<Return>", self.OnPressEnterDDfield)
###########################################

########### PRODUCT DIR ##################
        # Below is the sidetext for PRODUCTDIR
        self.sidetextpd = Tkinter.StringVar()
        sidetextlabelpd = Tkinter.Label(self,textvariable=self.sidetextpd, anchor="w",fg="blue",bg="white")
        sidetextlabelpd.grid(column=0,row=2,columnspan=1)
        self.sidetextpd.set("PRODUCT ROOT DIR:")
        # Below is the text box for string entry
        self.pdentryVariable = Tkinter.StringVar()
        self.pdentry = Tkinter.Entry(self,textvariable=self.pdentryVariable)
        self.pdentry.grid(column=1,row=2,columnspan=6,sticky='EW')
        self.pdentry.bind("<Return>", self.OnPressEnterPDfield)
###########################################
        ########### INSTRUMENT MODE ##############
        # Below is the sidetext for INSTRUMENT MODE
        self.sidetextpd2 = Tkinter.StringVar()
        sidetextlabelpd2 = Tkinter.Label(self,textvariable=self.sidetextpd2, anchor="w",fg="black",bg="lightgray")
        sidetextlabelpd2.grid(column=0,row=3,columnspan=1, sticky='EW')
        self.sidetextpd2.set("INSTRUMENT MODE:")
        # Below are the radio buttons to select instrument mode
        self.instmode = Tkinter.IntVar()
        
        self.R4 = Tkinter.Radiobutton(self, text="ALL", variable=self.instmode, value=0, anchor="e")
        self.R4.grid(column=1,row=3,columnspan=1)
        
        self.R1 = Tkinter.Radiobutton(self, text="Star Only", variable=self.instmode,value=1, anchor="e")
        self.R1.grid(column=2,row=3,columnspan=1)
        
        self.R2 = Tkinter.Radiobutton(self, text="Star+Sky", variable=self.instmode, value=2, anchor="e")
        self.R2.grid(column=3,row=3,columnspan=1)
        
        self.R3 = Tkinter.Radiobutton(self, text="Polar", variable=self.instmode, value=3, anchor="e")
        self.R3.grid(column=4,row=3,columnspan=1)

        # Default is Star-Only FOURSLICE mode
        self.R4.select();
        
        # Below it sets the header keyword value for chosen mode
        self.instmodekeyword = Tkinter.StringVar()
        if (self.instmode.get() == 1) :
            self.instmodekeyword.set("Spectroscopy, star only, R=80,000")
        elif (self.instmode.get() == 2) :
            self.instmodekeyword.set("Spectroscopy, star+sky, R=65,000")
        elif (self.instmode.get() == 3) :
            self.instmodekeyword.set("Polarimetry, R=65,000")
        
        ###########################################
        
        ############# READOUT SPEED ###############
        # Below is the sidetext for READOUT SPEED
        self.sidetextpd3 = Tkinter.StringVar()
        sidetextlabelpd3 = Tkinter.Label(self,textvariable=self.sidetextpd3, anchor="w",fg="black",bg="lightgray")
        sidetextlabelpd3.grid(column=0,row=4,columnspan=1, sticky='EW')
        self.sidetextpd3.set("READOUT SPEED:")
        # Below are the radio buttons to select readout speed
        self.readoutmode = Tkinter.IntVar()
        self.rb1 = Tkinter.Radiobutton(self, text="Fast", variable=self.readoutmode, value=1)
        self.rb1.grid(column=2,row=4,columnspan=1)
        
        self.rb2 = Tkinter.Radiobutton(self, text="Normal", variable=self.readoutmode, value=2)
        self.rb2.grid(column=3,row=4,columnspan=1)
        
        self.rb3 = Tkinter.Radiobutton(self, text="Slow", variable=self.readoutmode, value=3)
        self.rb3.grid(column=4,row=4,columnspan=1)

        self.rb4 = Tkinter.Radiobutton(self, text="ALL", variable=self.readoutmode, value=0)
        self.rb4.grid(column=1,row=4,columnspan=1)
    
        # Default is Normal readout speed
        self.rb4.select();
 
        # Below it sets the header keyword value for chosen mode
        self.readspeedkeyword = Tkinter.StringVar()
        if (self.readoutmode.get() == 1) :
            self.readspeedkeyword.set("Fast: 4.70e noise, 1.60e/ADU, 32s")
        elif (self.readoutmode.get() == 2) :
            self.readspeedkeyword.set("Normal: 4.20e noise, 1.30e/ADU, 38s")
        elif (self.readoutmode.get() == 3) :
            self.readspeedkeyword.set("Slow: 2.90e noise, 1.20e/ADU, 60s")
        
        ###########################################

################# NIGHT ###################
       # Below is the sidetext for NIGHT
        self.sidetextpd5 = Tkinter.StringVar()
        sidetextlabelpd5 = Tkinter.Label(self,textvariable=self.sidetextpd5, anchor="w",fg="white",bg="darkblue")
        sidetextlabelpd5.grid(column=0,row=5,columnspan=1)
        self.sidetextpd5.set("NIGHT:")
        # Below is the text box for NIGHT string entry
        self.nightentryVariable = Tkinter.StringVar()
        self.nightentry = Tkinter.Entry(self,textvariable=self.nightentryVariable)
        self.nightentry.grid(column=1,row=5,columnspan=2, sticky='EW')
        self.nightentry.bind("<Return>", self.OnPressEnterNightfield)
        # Below is the button to CREATE REDUCTION NIGHT DIR
        self.createNightButton = Tkinter.Button(self,text=u"Create Night DIR", command=self.OnCreateNightButtonClick, state="normal")
        self.createNightButton.grid(column=3,row=5,columnspan=1)
###########################################

############# MODE SELECTION OPTIONS ###############
        # Below is the sidetext for MODE SELECTION OPTIONS
        self.sidetextpd4 = Tkinter.StringVar()
        sidetextlabelpd4 = Tkinter.Label(self,textvariable=self.sidetextpd4, anchor="w",fg="black",bg="lightgray")
        sidetextlabelpd4.grid(column=0,row=6,columnspan=1, sticky='EW')
        self.sidetextpd4.set("MODE SELECTION:")
        
        # Below is the check button to select OBJECTS only
        self.withobjectscheck = Tkinter.IntVar()
        self.cb1 = Tkinter.Checkbutton(self, text="Object only", variable=self.withobjectscheck, onvalue = 1, offvalue = 0)
        self.cb1.grid(column=1,row=6,columnspan=1)
        self.cb1.select();
        
        # Below it initializes the string for DEFAULT CALIBRATION directory
        self.defaultCalDirVariable = Tkinter.StringVar()
        self.defaultCalDirVariable.set("DEFAULT CALIBRATION dir: /Users/edermartioli/espadonspipeline-1.0/DefaultCalibration/")

        # Below is the check button for DEFAULT CALIBRATION  option
        self.defaultcalibrationcheck = Tkinter.IntVar()
        self.cb4 = Tkinter.Checkbutton(self, text="Default Calib", variable=self.defaultcalibrationcheck, onvalue = 1, offvalue = 0)
        self.cb4.grid(column=2,row=6,columnspan=1)
        #self.cb4.select();

        # Below is the check button for ANY READOUT option
        self.anyreadoutcheck = Tkinter.IntVar()
        self.cb5 = Tkinter.Checkbutton(self, text="Any Readout", variable=self.anyreadoutcheck, onvalue = 1, offvalue = 0)
        self.cb5.grid(column=3,row=6,columnspan=1)
        #self.cb4.select();

###########################################

############# REDUCTION OPTIONS ###############
        # Below is the sidetext for REDUCTION OPTIONS
        self.sidetextredop = Tkinter.StringVar()
        sidetextlabelredop = Tkinter.Label(self,textvariable=self.sidetextredop, anchor="w",fg="black",bg="lightgray")
        sidetextlabelredop.grid(column=0,row=7,columnspan=1, sticky='EW')
        self.sidetextredop.set("REDUCTION OPTIONS:")
        
        # Below are the check buttons to select REDUCTION OPTIONS
        self.plotcheck = Tkinter.IntVar()
        self.ro1 = Tkinter.Checkbutton(self, text="plot", variable=self.plotcheck, onvalue = 1, offvalue = 0)
        self.ro1.grid(column=1,row=7,columnspan=1)
        self.ro1.select();
        
        self.verbosecheck = Tkinter.IntVar()
        self.ro2 = Tkinter.Checkbutton(self, text="verbose", variable=self.verbosecheck, onvalue = 1, offvalue = 0)
        self.ro2.grid(column=2,row=7,columnspan=1)
        self.ro2.select();
        
        self.tracecheck = Tkinter.IntVar()
        self.ro3 = Tkinter.Checkbutton(self, text="trace", variable=self.tracecheck, onvalue = 1, offvalue = 0)
        self.ro3.grid(column=3,row=7,columnspan=1)
        self.ro3.select();
        
        self.cleancheck = Tkinter.IntVar()
        self.ro4 = Tkinter.Checkbutton(self, text="clean", variable=self.cleancheck, onvalue = 1, offvalue = 0)
        self.ro4.grid(column=4,row=7,columnspan=1)
        #self.ro4.select();
###########################################

        ############# TARGET PRODUCTS ###############
        # Below is the sidetext for TARGET PRODUCTS
        self.sidetexttp = Tkinter.StringVar()
        sidetextlabeltp = Tkinter.Label(self,textvariable=self.sidetexttp, anchor="w",fg="white",bg="darkblue")
        sidetextlabeltp.grid(column=0,row=8,columnspan=1)
        self.sidetexttp.set("CALIB. PROUCT:")

        # Below are the radio buttons to select final products

        self.targetproduct = Tkinter.StringVar()

        # Button for Geometry
        self.tp1 = Tkinter.Radiobutton(self, text="Geometry", variable=self.targetproduct, value="GEOMETRYPRODUCT", anchor="w")
        self.tp1.grid(column=1,row=8,columnspan=1)

        # Button for Inst. Profile
        self.tp2 = Tkinter.Radiobutton(self, text="Inst. Profile", variable=self.targetproduct, value="INSTRUMENTPROFILEPRODUCT", anchor="w")
        self.tp2.grid(column=2,row=8,columnspan=1)

        # Button for Aperture
        self.tp3 = Tkinter.Radiobutton(self, text="Aperture", variable=self.targetproduct, value="APERTUREPRODUCT", anchor="w")
        self.tp3.grid(column=3,row=8,columnspan=1)

        # Button for Wavelength
        self.tp4 = Tkinter.Radiobutton(self, text="Wavelength", variable=self.targetproduct, value="WAVELENGTHPRODUCT", anchor="w")
        self.tp4.grid(column=4,row=8,columnspan=1)

        self.sidetexttrp = Tkinter.StringVar()
        sidetextlabeltrp = Tkinter.Label(self,textvariable=self.sidetexttrp, anchor="c",fg="white",bg="darkblue")
        sidetextlabeltrp.grid(column=0,row=9,columnspan=1)
        self.sidetexttrp.set("REDUCE PRODUCT:")

        # Button for Reduce Objects (OPERA products)
        self.tp5 = Tkinter.Radiobutton(self, text="Object Spectra", variable=self.targetproduct, value="OBJECTS", anchor="w")
        self.tp5.grid(column=1,row=9,columnspan=1)

        # Button for Libre-Esprit products
        self.tp6 = Tkinter.Radiobutton(self, text="Libre-Esprit", variable=self.targetproduct, value="LIBRE-ESPRIT", anchor="w")
        self.tp6.grid(column=2,row=9,columnspan=1)
    
        # Button for Libre-Esprit products
        # self.tp7 = Tkinter.Radiobutton(self, text="FITS", variable=self.targetproduct, value="FITS", anchor="w")
        # self.tp7.grid(column=3,row=9,columnspan=1)
                
        # Default is ...
        self.tp5.select();
        ###########################################

############### NIGHT LOG  ################
        # Below is the button to generate NIGHT LOG
        nightLogButton = Tkinter.Button(self,text="Gen Night Log", command=self.OnGenNightLogButtonClick)
        nightLogButton.grid(column=0,row=14,columnspan=1, sticky='EW')
###########################################

############### CHECK DATA ################
        # Below is the button to CHECK DATA
        checkDataButton = Tkinter.Button(self,text="Check Data", command=self.OnCheckDataButtonClick)
        checkDataButton.grid(column=1,row=14,columnspan=2, sticky='EW')
###########################################

############### CLEAN Products #################
        # Below is the button to clean products
        self.cleanallFlag = False
        self.cleanProductsButton = Tkinter.Button(self,text="Clean Products", state="active", command=self.OnCleanProductsButtonClick)
        self.cleanProductsButton.grid(column=3,row=14,columnspan=1, sticky='EW')
###########################################
   
############### RUN SIMULATION #################
        # Below is the button to run a processing simulation
        self.simulationcheck = False
        self.runSimulationButton = Tkinter.Button(self,text="Run Simulation", state="active", command=self.OnRunSimulationButtonClick)
        self.runSimulationButton.grid(column=0,row=15,columnspan=1, sticky='EW')
###########################################
                
############### RUN ESPaDOnS PIPELINE #################
        # Below is the button to run ESPaDOnS PIPELINE
        self.runOperaButton = Tkinter.Button(self,text="Run Pipeline", state="active", command=self.OnRunOperaButtonClick)
        self.runOperaButton.grid(column=1,row=15,columnspan=2, sticky='EW')
###########################################

########### Exit button ##################
        # Below is supposed to be a button to abort execution, but so far it just closes python
        self.abortButton = Tkinter.Button(self,text="Exit", state="active", command=self.OnAbortButtonClick)
        self.abortButton.grid(column=3,row=15,columnspan=1, sticky='EW')
###########################################

################# CONSOLE #################
        # Strings for console output
        self.headbottombarlabelVariable = Tkinter.StringVar()
        self.bottombarlabelVariable1 = Tkinter.StringVar()
        self.bottombarlabelVariable2 = Tkinter.StringVar()
        self.bottombarlabelVariable3 = Tkinter.StringVar()

        # Construct labels for console
        headlabel = Tkinter.Label(self,textvariable=self.headbottombarlabelVariable, anchor="c",fg="darkgrey",bg="lightgray")
        self.headbottombarlabelVariable.set("ESPaDOnS Pipeline 1.0 by E. Martioli")
        headlabel.grid(column=0,row=16,columnspan=5, sticky='EW')

        label1 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable1, anchor="c",fg="black",bg="lightgray")
        label2 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable2, anchor="c",fg="black",bg="lightgray")
        #label3 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable3, anchor="c",fg="black",bg="lightgray")

        self.bottombarlabelVariable2.set(self.operaDirVariable.get())

        label1.grid(column=0,row=17,columnspan=5, sticky='EW')
        label2.grid(column=0,row=18,columnspan=5, sticky='EW')
        #label3.grid(column=0,row=19,columnspan=5, sticky='EW')
###########################################

        self.grid_columnconfigure(4,weight=1)
        self.resizable(True,False)
        self.update()
        self.geometry(self.geometry())

###########################################
#### End of initialize function ###########
###########################################

########### ABORT BUTTON ##################
    def OnAbortButtonClick(self):
        self.quit()
###########################################

######## CLEAN PRODUCTS BUTTON ############
    def OnCleanProductsButtonClick(self):
        self.cleanallFlag = True
        self.OnRunOperaButtonClick()
        self.cleanallFlag = False
###########################################
            
######## RUN SIMULATION BUTTON ############
    def OnRunSimulationButtonClick(self):
        self.simulationcheck = True
        self.OnRunOperaButtonClick()
        self.simulationcheck = False
###########################################
            
########### OPERA BUTTON ##################
    def OnRunOperaButtonClick(self):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("Running ESPaDOnS Pipeline for " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")
            self.update()
            print " "
            print "-----------------------------------------------------------------"
            print "Running ESPaDOnS Pipeline for NIGHT: " + self.nightentryVariable.get()
            print "-----------------------------------------------------------------\n"
            pipelinehomedir = self.operaDirVariable.get()
            datarootdir =self.ddentryVariable.get()
            productrootdir = self.pdentryVariable.get()
            
            specificProduct = self.targetproduct.get()
            
            if(len(specificProduct) > 1) :
                print "TARGET PRODUCT is: ", specificProduct, "\n"
            else :
                specificProduct = ""
                print "TARGET PRODUCT is: REDUCE ALL \n"
            
            """
            Set up directories:
            """
            Dirs = espadons.Directories(pipelinehomedir,datarootdir,productrootdir,self.nightentryVariable.get())
            """
            Set up config files:
            """
            config = espadons.ConfigFiles(Dirs)
            """
            Set up Espadons keywords:
            """
            keywords = espadons.Keywords()
        
            forcecalibration = True
            if self.withobjectscheck.get() :
                forcecalibration = False # This may be set to "True" when one wants to run calibrations even when there is no object files

            modes = espadons.ReductionModes(Dirs, keywords, self.anyreadoutcheck.get(), forcecalibration)
                                
            for mode in modes.getInstReadModes() :
                instrumentmode = mode[0]
                readoutspeed = mode[1]

                if (self.readoutmode.get() == readoutspeed or self.readoutmode.get() == 0) and \
                   (self.instmode.get() == instrumentmode or self.instmode.get() == 0) :
                       
                    modes.displayModeStats(instrumentmode,readoutspeed)
                    input = [self.nightentryVariable.get(),instrumentmode,readoutspeed,self.cleancheck.get(),self.simulationcheck,self.plotcheck.get(),self.verbosecheck.get(),self.tracecheck.get(),self.anyreadoutcheck.get(),specificProduct,self.cleanallFlag]
                    espadonspipeline.executePipeline(input, Dirs, config, keywords)
            
            print "-----------------------------------------------------------------"
            print " "
                        
            modes.cleanModes()
                        
            self.bottombarlabelVariable2.set("DONE!")

        elif not(os.path.exists(self.ddentryVariable.get())) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("Error: non-existent DATA DIR : " + self.ddentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
        elif os.path.exists(self.ddentryVariable.get()) and not(os.path.exists(self.pdentryVariable.get())) :
            self.bottombarlabelVariable1.set("Error: non-existent PRODUCT DIR : " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)
        else :
            self.bottombarlabelVariable1.set("Error: non-existent directories")
            self.bottombarlabelVariable2.set("")
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
###########################################


########### NIGHT LOG BUTTON #############
    def OnGenNightLogButtonClick(self):
        # test if all files exist
        nightdir = self.ddentryVariable.get()+self.nightentryVariable.get()+"/"
        if os.path.exists(nightdir) :
            self.bottombarlabelVariable1.set("Generating Night Log for " + nightdir + " ...")
            self.bottombarlabelVariable2.set("** See data list output on terminal **")

            print " "
            print "-----------------------------------------------------------------"
            print "NIGHT LOG for DATA in DIR: " + nightdir
            print "-----------------------------------------------------------------"
           # command line below prints out information on data available in nightdir
            commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE INSTMODE EREADSPD OBSTYPE EXPTIME OBJECT" -p'
            os.system(commandline)
            print "-----------------------------------------------------------------"
            print " "

        else :
            self.bottombarlabelVariable1.set("Error: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###########################################

########### CHECK DATA BUTTON #############
    def OnCheckDataButtonClick(self):
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("Checking data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")
            self.update()
            
            print " "
            print "-----------------------------------------------------------------"
            print "Running Check DATA for NIGHT: " + self.nightentryVariable.get()
            print "-----------------------------------------------------------------"
            
            pipelinehomedir = self.operaDirVariable.get()
            datarootdir =self.ddentryVariable.get()
            productrootdir = self.pdentryVariable.get()
            
            """
                Set up directories:
                """
            Dirs = espadons.Directories(pipelinehomedir,datarootdir,productrootdir,self.nightentryVariable.get())
            """
                Set up config files:
                """
            config = espadons.ConfigFiles(Dirs)
            """
                Set up Espadons keywords:
                """
            keywords = espadons.Keywords()
            
            forcecalibration = True
            if self.withobjectscheck.get() :
                forcecalibration = False # This may be set to "True" when one wants to run calibrations even when there is no object files
            modes = espadons.ReductionModes(Dirs, keywords, self.anyreadoutcheck.get(), forcecalibration)
        
            modes.displayOverallStats()
            #modes.displayObjectData()
            #modes.displayCalibrationData()

            modes.cleanModes()
            self.bottombarlabelVariable2.set("DONE!")
        
        elif not(os.path.exists(self.ddentryVariable.get())) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("Error: non-existent DATA DIR : " + self.ddentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
        elif os.path.exists(self.ddentryVariable.get()) and not(os.path.exists(self.pdentryVariable.get())) :
            self.bottombarlabelVariable1.set("Error: non-existent PRODUCT DIR : " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)
        else :
            self.bottombarlabelVariable1.set("Error: non-existent directories")
            self.bottombarlabelVariable2.set("")
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
###########################################

############## DATA DIR ###################
    def OnPressEnterDDfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) :
            self.bottombarlabelVariable1.set("DATA DIR is OK")
            self.bottombarlabelVariable2.set("")
        else:
            self.bottombarlabelVariable1.set("Error: non-existent directory: " + self.ddentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
###########################################

############ PRODUCT DIR ##################
    def OnPressEnterPDfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("PRODUCT DIR is OK")
            self.bottombarlabelVariable2.set("")

        else:
            self.bottombarlabelVariable1.set("Error: non-existent directory: " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)
###########################################

########### CONFIG FILE ##################
    def OnPressEnterConfigPathfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.isfile(self.filepathVariable.get()) :

            configFile = open(self.filepathVariable.get(), 'r')

            for configdatum in configFile:
                col0,col1 = configdatum.split(' = ')
                if(col0=="ESPADONSPIPELINEDIR") :
                    self.operaDirVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTDATADIR") :
                    self.ddentryVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTPRODUCTDIR") :
                    self.pdentryVariable.set(col1.rstrip('\n'))
                if(col0=="NIGHT") :
                    self.nightentryVariable.set(col1.rstrip('\n'))
                if(col0=="DEFAULTCALIBRATIONDIR") :
                    self.defaultCalDirVariable.set(col1.rstrip('\n'))

            configFile.close()
            self.bottombarlabelVariable1.set("Config File loaded: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("ESPaDOnS Pipeline home dir:" + self.operaDirVariable.get())

        else:
            self.bottombarlabelVariable1.set("Error: non-existent Config File: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("")
            self.filepathentry.focus_set()
            self.filepathentry.selection_range(0, Tkinter.END)

    def LoadConfigButtonClick(self):
        # test if directory exists if not it will ask for it.
        if os.path.isfile(self.filepathVariable.get()) :

            configFile = open(self.filepathVariable.get(), 'r')

            for configdatum in configFile:
                col0,col1 = configdatum.split(' = ')
                if(col0=="ESPADONSPIPELINEDIR") :
                    self.operaDirVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTDATADIR") :
                    self.ddentryVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTPRODUCTDIR") :
                    self.pdentryVariable.set(col1.rstrip('\n'))
                if(col0=="NIGHT") :
                    self.nightentryVariable.set(col1.rstrip('\n'))
                if(col0=="DEFAULTCALIBRATIONDIR") :
                    self.defaultCalDirVariable.set(col1.rstrip('\n'))

            configFile.close()
            self.bottombarlabelVariable1.set("Config File laoded: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("ESPaDOnS Pipeline home dir:" + self.operaDirVariable.get())
        else:
            self.bottombarlabelVariable1.set("Error: non-existent file: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("")
            self.filepathentry.focus_set()
            self.filepathentry.selection_range(0, Tkinter.END)
###########################################

############### NIGHT #####################
    def OnPressEnterNightfield(self,event):
        # test if NIGHT directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get() + "/" + self.nightentryVariable.get()) :
            self.bottombarlabelVariable1.set("NIGHT DATA DIR is OK: " + self.ddentryVariable.get() + "/" + self.nightentryVariable.get())
            if os.path.exists(self.pdentryVariable.get() + "/" + self.nightentryVariable.get()) :
                self.bottombarlabelVariable2.set("NIGHT PRODUCT DIR already exists:" + self.pdentryVariable.get() + "/" + self.nightentryVariable.get())
            else :
                self.bottombarlabelVariable2.set("Click on <Night DIR> to create NIGHT PRODUCT dir:" + self.pdentryVariable.get() + "/" + self.nightentryVariable.get())

        else:
            self.bottombarlabelVariable1.set("Error: non-existent directory: " + self.ddentryVariable.get() + "/" + self.nightentryVariable.get())
            self.bottombarlabelVariable2.set("")

            self.nightentry.focus_set()
            self.nightentry.selection_range(0, Tkinter.END)

    def OnCreateNightButtonClick(self):
        nightdir = self.pdentryVariable.get() + "/" + self.nightentryVariable.get()
        if os.path.exists(nightdir) :
            self.bottombarlabelVariable1.set("Error: night directory already exists -> " + nightdir)
            self.bottombarlabelVariable2.set("")
        else :
            # create new night directory
            self.bottombarlabelVariable1.set("New directory created: " + nightdir)
            self.bottombarlabelVariable2.set("")

            createNightDirCommand = "mkdir " + nightdir
            os.system(createNightDirCommand)

###########################################

if __name__ == "__main__":
    app = pipelineFrontEnd_tk(None)
    app.title('ESPaDOnS Pipeline v. 1.0')
    app.mainloop()
