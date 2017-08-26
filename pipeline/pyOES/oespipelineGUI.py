#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
"""
Created on Wed May 28 21:36:03 2014
@author: Eder Martioli, Romana Grossová
Laboratorio Nacional de Astrofisica, Brazil; Masaryk University, Brno, Czech Republic
Last update on Aug 24, 2016
"""
import Tkinter
import Tkinter, tkFileDialog
import os
import glob

###########################################
############ MESSAGE COMMANDS #############
###########################################

######### Print message command ###########
def print_message(message):
    print " "                       
    print "-----------------------------------------------------------------"
    print message
    print "-----------------------------------------------------------------"
###########################################

####### Print check message command #######
def print_infomessage(infomessage):
    print " "                       
    print "****************************************************************************"
    print infomessage
    print "****************************************************************************"
###########################################

###### Print ERROR message command ########
def print_ERRORmessage(ERRORmessage):
    print " "
    print "----------------------------------------------------------------------------"                  
    print "!!! " + ERRORmessage + " !!!"
    print "----------------------------------------------------------------------------"
###########################################

###########################################
######### SET PATH TO CONFIG FILE #########
###########################################
configfilepath = os.path.realpath(__file__).strip('oespipelineGUI.py') + 'config.oespipeline'
if os.path.isfile(configfilepath) :
    print_infomessage("RUNNING OES PIPELINE GUI FRAMEWORK.")   
    print("CONFIG FILE: " + configfilepath)
else :
    print_ERRORmessage( "ERROR: Non-existent CONFIG file \'" + configfilepath + "\'")
###########################################

###########################################
############ MAIN FRONT END ###############
###########################################

class pipelineFrontEnd_tk(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

############ CONFIG FILE ##################   
    def OnPressEnterConfigPathfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.isfile(self.filepathVariable.get()) :    
            
            configFile = open(self.filepathVariable.get(), 'r')

            for configdatum in configFile:
                col0,col1 = configdatum.split(' = ')
                if(col0=="OESPIPELINEDIR") :
                    self.operaDirVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTDATADIR") :
                    self.ddentryVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTPRODUCTDIR") :
                    self.pdentryVariable.set(col1.rstrip('\n'))
                if(col0=="READOUT") :
                    self.readoutVariable.set(col1.rstrip('\n'))
                if(col0=="PLOTVIEWER") :
                    self.plotViewerVariable.set(col1.rstrip('\n'))
                if(col0=="NIGHT") :
                    self.nightentryVariable.set(col1.rstrip('\n'))
                if(col0=="FINALPRODDIR") :
                    self.finalProdDirVariable.set(col1.rstrip('\n'))
                

            configFile.close()
                                     
            self.bottombarlabelVariable1.set("Config File loaded: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("OES Pipeline home dir:" + self.operaDirVariable.get())
            
        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent Config File: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("")
            self.filepathentry.focus_set()
            self.filepathentry.selection_range(0, Tkinter.END)          

########### LOAD CONFIG BUTTON ###########
    def LoadConfigButtonClick(self):
        # test if directory exists if not it will ask for it.
        if os.path.isfile(self.filepathVariable.get()) :  

            configFile = open(self.filepathVariable.get(), 'r')

            for configdatum in configFile:
                col0,col1 = configdatum.split(' = ')
                if(col0=="OESPIPELINEDIR") :
                    self.operaDirVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTDATADIR") :
                    self.ddentryVariable.set(col1.rstrip('\n'))
                if(col0=="ROOTPRODUCTDIR") :
                    self.pdentryVariable.set(col1.rstrip('\n'))
                if(col0=="READOUT") :
                    self.readoutVariable.set(col1.rstrip('\n'))
                if(col0=="PLOTVIEWER") :
                    self.plotViewerVariable.set(col1.rstrip('\n'))
                if(col0=="NIGHT") :
                    self.nightentryVariable.set(col1.rstrip('\n'))
                if(col0=="FINALPRODDIR") :
                    self.finalProdDirVariable.set(col1.rstrip('\n'))
                

            configFile.close() 
                                   
            self.bottombarlabelVariable1.set("Config File laoded: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("OES Pipeline home dir:" + self.operaDirVariable.get())
            self.fileloaded = True

        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent file: " + self.filepathVariable.get())
            self.bottombarlabelVariable2.set("LOAD CONFIG")            
            self.filepathentry.focus_set()
            self.filepathentry.selection_range(0, Tkinter.END)
        
###########################################

    def initialize(self):
        self.grid()
        # primary set up of loaded files to False
        self.fileloaded = False        
        
########### CONFIG FILE ##################   
        # Below is the button to load CONFIG FILE
        loadConfigButton = Tkinter.Button(self,text="Load Config", command=self.LoadConfigButtonClick)        
        loadConfigButton.grid(column=0,row=0)
        
        # Below it initializes the string for opera system directory
        self.operaDirVariable = Tkinter.StringVar()
        self.operaDirVariable.set(self.operaDirVariable.get())

        # Below is the text box for config file path entry
        self.filepathVariable = Tkinter.StringVar()
        self.filepathVariable.set(configfilepath)
        self.filepathentry = Tkinter.Entry(self,textvariable=self.filepathVariable)
        self.filepathentry.grid(column=1,row=0,columnspan=3,sticky='EW')
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
        self.ddentry.grid(column=1,row=1,columnspan=3,sticky='EW')
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
        self.pdentry.grid(column=1,row=2,columnspan=3,sticky='EW')
        self.pdentry.bind("<Return>", self.OnPressEnterPDfield)
###########################################

########### IMAGE VIEWER ##################   
        # Below is the sidetext for image viewer
        self.sidetextdd = Tkinter.StringVar()      
        # Below is the text box for string entry
        self.plotViewerVariable = Tkinter.StringVar()
        self.vventry = Tkinter.Entry(self,textvariable=self.plotViewerVariable)
        self.vventry.bind("<Return>", self.OnPressViewerVariable)
###########################################

########### READOUT SPEED ##################   
        # Below is the sidetext for readout speed
        self.sidetextdd = Tkinter.StringVar()      
        # Below is the text box for string entry
        self.readoutVariable = Tkinter.StringVar()
        self.vventry = Tkinter.Entry(self,textvariable=self.readoutVariable)
        self.vventry.bind("<Return>", self.OnPressReadoutVariable)
###########################################

################# NIGHT ###################           
       # Below is the sidetext for NIGHT
        self.sidetextpd3 = Tkinter.StringVar()        
        sidetextlabelpd3 = Tkinter.Label(self,textvariable=self.sidetextpd3, anchor="w",fg="white",bg="darkblue")
        sidetextlabelpd3.grid(column=0,row=3,columnspan=1)        
        self.sidetextpd3.set("NIGHT:")
        # Below is the text box for NIGHT string entry
        self.nightentryVariable = Tkinter.StringVar()
        self.nightentry = Tkinter.Entry(self,textvariable=self.nightentryVariable)
        self.nightentry.grid(column=1,row=3,columnspan=1, sticky='EW')
        self.nightentry.bind("<Return>", self.OnPressEnterNightfield)
        # Below is the button to CREATE REDUCTION NIGHT DIR
        self.createNightButton = Tkinter.Button(self,text=u"Create Night DIR", command=self.OnCreateNightButtonClick, state="normal")        
        self.createNightButton.grid(column=2,row=3,columnspan=1)
########################################### 

##### SELECT REDUCTION PRODUCT ############ 
        # Below is the sidetext for REDUCTION PRODUCT       
        self.sidetextpd = Tkinter.StringVar()        
        sidetextlabelpd = Tkinter.Label(self,textvariable=self.sidetextpd, anchor="w",fg="black",bg="white")
        sidetextlabelpd.grid(column=0,row=5,columnspan=1, rowspan=2)        
        self.sidetextpd.set("SELECT PRODUCT:")
        # Below are the check buttons to select REDUCTION PRODUCT
        self.products = Tkinter.StringVar()
        self.p1 = Tkinter.Radiobutton(self, text="OBJECT", variable=self.products, value="'OBJECT'")
        self.p1.grid(column=2,row=5,columnspan=1)        
                
        self.p2 = Tkinter.Radiobutton(self, text="CALIBRATION", variable=self.products, value="'CALIBRATION'")
        self.p2.grid(column=1,row=5,columnspan=1)        
        
        self.p3 = Tkinter.Radiobutton(self, text="1D SPECTRUM", variable=self.products, value="'ONEDSPEC'")
        self.p3.grid(column=1,row=6,columnspan=1)  

        self.p4 = Tkinter.Radiobutton(self, text="STACK", variable=self.products, value = "'STACK'")
        self.p4.grid(column=2,row=6,columnspan=1)

        # Default is "OBJECT"
        self.p1.select();

#############################################
     
############### CHECK PLOTS ###############
        # Below is the button to check plots of the products
        self.sidetextpd = Tkinter.StringVar()        
        sidetextlabelpd = Tkinter.Label(self,textvariable=self.sidetextpd, anchor="w",fg="white",bg="black")
        sidetextlabelpd.grid(column=0,row=7,columnspan=1, rowspan=2)        
        self.sidetextpd.set("CHECK PLOTS:")
        self.GeomPlotsButton = Tkinter.Button(self,text=u"Geometry Plots", state="active", command=self.OnGeomPlotsButtonClick)        
        self.GeomPlotsButton.grid(column=1,row=7,columnspan=1,sticky='EW')

        self.AperPlotsButton = Tkinter.Button(self,text=u"Aperture Plots", state="active", command=self.OnAperPlotsButtonClick, width=20)        
        self.AperPlotsButton.grid(column=2,row=7,columnspan=1,sticky='EW')

        self.CalibPlotsButton = Tkinter.Button(self,text=u"Calibration Plots", state="active", command=self.OnCalibPlotsButtonClick)        
        self.CalibPlotsButton.grid(column=1,row=8,columnspan=1,sticky="EW")
    
        self.FcalPlotsButton = Tkinter.Button(self,text=u"Flat Calib Plots", state="active", command=self.OnFcalPlotsButtonClick)        
        self.FcalPlotsButton.grid(column=2,row=8,columnspan=1,sticky="EW")
########################################### 

############## FINAL PLOTS ################
        self.sidetextpd = Tkinter.StringVar()        
        sidetextlabelpd = Tkinter.Label(self,textvariable=self.sidetextpd, anchor="w",fg="white",bg="black")
        sidetextlabelpd.grid(column=0,row=10,columnspan=1,rowspan=2)        
        self.sidetextpd.set("FINAL SPECTRUM PLOTS:")
        self.RawPlotsButton = Tkinter.Button(self,text=u"Raw spectrum plot", state="active", command=self.OnRawPlotButtonClick,height=1, width=22)
        self.RawPlotsButton.grid(column=1,row=10,columnspan=1,sticky='EW')

        self.FlatPlotsButton = Tkinter.Button(self,text=u"Calibrated spectrum plot", state="active", command=self.OnFlatPlotButtonClick,height=1, width=20)
        self.FlatPlotsButton.grid(column=2,row=10,columnspan=1,sticky='EW')

        self.NormButton = Tkinter.Button(self,text=u"Normalized spectrum plot", state="active", command=self.OnNormButtonClick,height=1, width=23)
        self.NormButton.grid(column=1,row=11,columnspan=1,sticky='EW')

        self.ONEDButton = Tkinter.Button(self,text=u"1D spectrum plot", state="active", command=self.OnONEDButtonClick,height=1, width=22)
        self.ONEDButton.grid(column=2,row=11,columnspan=1,sticky='EW')
#############################################

 ############### NIGHT LOG  ################        
        # Below is the button to generate NIGHT LOG
        nightLogButton = Tkinter.Button(self,text=u"Generate Night Log", command=self.OnGenNightLogButtonClick)        
        nightLogButton.grid(column=0,row=15,columnspan=1,sticky='EW')
########################################### 

############### RUN OES PIPELINE ######       
        # Below is the button to run OES PIPELINE
        self.runOperaButton = Tkinter.Button(self,text=u"Run Pipeline", state="active", command=self.OnRunOperaButtonClick)        
        self.runOperaButton.grid(column=1,row=15,columnspan=1,sticky='EW') 
########################################### 

######### CHECK RAW IMAGES IN DS9 #########     
        # Below is the button to check raw images in ds9
        self.Onds9Button = Tkinter.Button(self,text=u"Check IMAGE (ds9)", state="active", command=self.Onds9ButtonClick)        
        self.Onds9Button.grid(column=0,row=16,columnspan=1,sticky='EW')
########################################### 

##### RUN OES PIPELINE SIMULATION #####       
        # Below is the button to run OES PIPELINE
        self.runOperaSButton = Tkinter.Button(self,text=u"Run Simulation", state="active", command=self.OnRunOperaSimButtonClick)        
        self.runOperaSButton.grid(column=1,row=16,columnspan=1,sticky='EW') 
###########################################
        
########### Abort button ##################   
        # Below is button to close GUI 
        self.abortButton = Tkinter.Button(self,text="Exit", state="active", command=self.OnAbortButtonClick,height=1, width=20)        
        self.abortButton.grid(column=2,row=15,columnspan=1,rowspan=2,sticky='EW')
###########################################

### FINAL PRODUCT DIR - ONLY *.SPC ########         
        # Below it initializes the string for DEFAULT CALIBRATION directory       
        self.finalProdDirVariable = Tkinter.StringVar()        
        self.finalProdDirVariable.set(self.pdentryVariable.get() + "products" + self.nightentryVariable.get() + "/")
###########################################

### DEFINE MODE ########         
        # Below it initializes the string for DEFAULT CALIBRATION directory       
        self.modeName = Tkinter.StringVar()        
        self.modeName.set("OESBLADE")
###########################################

################# CONSOLE #################        
        # Strings for console output
        self.headbottombarlabelVariable = Tkinter.StringVar()        
        self.bottombarlabelVariable1 = Tkinter.StringVar()  
        self.bottombarlabelVariable2 = Tkinter.StringVar()  
        self.bottombarlabelVariable3 = Tkinter.StringVar()  

        # Construct labels for console       
        label1 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable1, anchor="c",fg="darkgray",bg="lightgray")  
        self.bottombarlabelVariable1.set("OES Pipeline 1.0 Console by E. Martioli, R.Grossova")      
        label2 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable2, anchor="c",fg="black",bg="lightgray")
        self.bottombarlabelVariable2.set(self.operaDirVariable.get())
             
        label1.grid(column=0,row=13,columnspan=3,sticky='EW')        
        label2.grid(column=0,row=12,columnspan=3,sticky='EW')                      
###########################################
        self.grid_columnconfigure(4,weight=1)
        self.resizable(True,False)
        self.update()
        self.geometry(self.geometry())     
###########################################

###########################################
#### End of initialize function ###########
###########################################

########### ABORT BUTTON ##################   
    def OnAbortButtonClick(self):
        self.quit()
###########################################
               
###########################################
############## RUN PIPELINE ###############
###########################################

########### OPERA BUTTON ##################   
    def OnRunOperaButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return

        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :  
            self.bottombarlabelVariable1.set("Running OES Pipeline for " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()     
            print_infomessage("Running OES Pipeline for NIGHT: " + self.nightentryVariable.get())
            
            # command line to run pipeline script             
            commandline = self.operaDirVariable.get() + '/pipeline/pyOES/operaoes.py' \
            ' ' + "--datarootdir=" + self.ddentryVariable.get() + \
            ' ' + "--pipelinehomedir=" + self.operaDirVariable.get() + \
            ' ' + "--productrootdir=" + self.pdentryVariable.get() + \
            ' ' + "--night=" + self.nightentryVariable.get() +\
            ' ' + "--product=" + str(self.products.get()) + \
            ' ' + "-pvt"
            print "Issuing command: " + commandline
            os.system(commandline)
            print_message("DONE!")
            self.bottombarlabelVariable2.set("DONE!")         

        elif not(os.path.exists(self.ddentryVariable.get())) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("ERROR: non-existent DATA DIR : " + self.ddentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
        elif os.path.exists(self.ddentryVariable.get()) and not(os.path.exists(self.pdentryVariable.get())) :
            self.bottombarlabelVariable1.set("ERROR: non-existent PRODUCT DIR : " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent directories")
            self.bottombarlabelVariable2.set("")            
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)      
###########################################
            
########### OPERA SIMULATION BUTTON #######
    def OnRunOperaSimButtonClick(self):
        
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :  
            self.bottombarlabelVariable1.set("Running OES Pipeline Simulation for " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()    
            print_infomessage("Running OES Pipeline Simulation for NIGHT: " + self.nightentryVariable.get() + " " + self.operaDirVariable.get() )
            # command line to run pipeline simulation script
            commandline = self.operaDirVariable.get() + '/pipeline/pyOES/operaoes.py' \
            ' ' + "--datarootdir=" + self.ddentryVariable.get() + \
            ' ' + "--pipelinehomedir=" + self.operaDirVariable.get() + \
            ' ' + "--productrootdir=" + self.pdentryVariable.get() + \
            ' ' + "--night=" + self.nightentryVariable.get() +\
            ' ' + "--product=" + str(self.products.get()) + \
            ' ' + "-pvts"
            print "Issuing command: " + commandline
            os.system(commandline)
            print_message("SIMULATION DONE!")
            self.bottombarlabelVariable2.set("DONE!")       

        elif not(os.path.exists(self.ddentryVariable.get())) and os.path.exists(self.pdentryVariable.get()) :
            self.bottombarlabelVariable1.set("ERROR: non-existent DATA DIR : " + self.ddentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)
        elif os.path.exists(self.ddentryVariable.get()) and not(os.path.exists(self.pdentryVariable.get())) :
            self.bottombarlabelVariable1.set("ERROR: non-existent PRODUCT DIR : " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent directories")
            self.bottombarlabelVariable2.set("")            
            self.ddentry.focus_set()
            self.ddentry.selection_range(0, Tkinter.END)      
##########################################

##########################################
###### CHECK CALIBRATION PRODUCTS ########          
##########################################

############# GEOMETRY PLOTS ############# 
    def OnGeomPlotsButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the products in " + self.pdentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the geometry and spacing products for night: " + self.nightentryVariable.get())

            # command line below vizualizes calibration result 
            spcprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_spcplot.eps" 
            geomprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_geomplot.eps" 
            # check existens of geometry plots    
            if os.path.exists(spcprod):
                commandline = self.plotViewerVariable.get() + ' ' + spcprod

                print "Issuing command to plot order spacing plot: " + commandline
                os.system(commandline)
                print "-------------------------------------"
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + spcprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"GEOMETRY\" to produce geometry plots.")
            pass

            if  os.path.exists(geomprod):
                    commandline_1 = self.plotViewerVariable.get() + ' ' + geomprod
                    print "Issuing command to plot geometry calibration plot: " + commandline_1
                    os.system(commandline_1)
                    print_message("Printing plots finished.")
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + geomprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"GEOMETRY\" to produce geometry plots.")
            pass

        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
##########################################
            
############# APERTURE PLOTS #############      
    def OnAperPlotsButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the products in " + self.pdentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the aperture and profile products for night: " + self.nightentryVariable.get())
            
            profprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_profplot.eps"
            aperprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_aperplot.eps"
            tiltprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_tiltplot.eps"
            #check existens of calibration plots
            if os.path.exists(profprod):
                commandline = self.plotViewerVariable.get() + ' ' + profprod
                print "Issuing command to plot instrumental profile calibration plot: " + commandline
                os.system(commandline)
                print "-------------------------------------"
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + profprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce calibration plots.")
            pass

            if os.path.exists(aperprod):
                commandline_1 = self.plotViewerVariable.get() + ' ' + aperprod
                print "Issuing command to plot aperture calibration plot: " + commandline_1
                os.system(commandline_1)
                print "-------------------------------------"
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + aperprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce calibration plots.")
            pass

            if os.path.exists(tiltprod):
                commandline_2 = self.plotViewerVariable.get() + ' ' + tiltprod
                print "Issuing command to plot aperture titl plot: " + commandline_2
                os.system(commandline_2)
                print_message("Printing plots finished.")
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + tiltprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce calibration plots.")
            pass

        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################

########## WAVELENGHT CALIB PLOTS ############# 
    def OnCalibPlotsButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the products in " + self.pdentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the wavelength solution products for night: " + self.nightentryVariable.get())

            # command line below vizualizes calibration result 
            wavespecprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_wavespecplot.eps" 
            waveordsprod = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_waveordsplot.eps"
            # check existens of wavelength calibration plots    
            if os.path.exists(wavespecprod):
                commandline = self.plotViewerVariable.get() + ' ' + wavespecprod
                print "Issuing command to of wavelength solution as echellogram: " + commandline
                os.system(commandline)
                print "-------------------------------------"
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + wavespecprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce wavelength calibration plots.")
            pass

            if  os.path.exists(waveordsprod):
                commandline_1 = self.plotViewerVariable.get() + ' ' + waveordsprod
                print "Issuing command to plot wavelength solution statistics: " + commandline
                os.system(commandline_1)
                print_message("Printing plots finished.")
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + waveordsprod + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent files.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce wavelength calibration plots.")
            pass

        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################

############ FLAT FIELD CALIB PLOTS ###########
    def OnFcalPlotsButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the products in " + self.pdentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the wavelength solution products for night: " + self.nightentryVariable.get())
            # command line below vizualizes calibration result 
            gzffcal = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_flat.fcal.gz"            
            ffile = self.pdentryVariable.get() + self.nightentryVariable.get() + '/' + self.nightentryVariable.get() + "_" + self.modeName.get() + "_" + self.readoutVariable.get() + "_flat.fcal"
            name = ffile.split('/')[-1]
            # check existens of flat flux calibration plot
            if os.path.exists(gzffcal):
                os.system('gunzip -c %s  > %s'  % (gzffcal, ffile))
                print_message("Uncompressed file \'%s\'" % gzffcal)                
                os.system("echo 'set title \"Flat flux calibration\" \n set xlabel \"Wavelength (nm)\" \n set ylabel \"order number + normalized flux\" \n  plot \"%s\" u 6:7 lt -1 t \"%s\"' |  gnuplot -p"  % (ffile, name)) 
                       
                print ("Printing plot \'" + name +  "\'...")
           
            else:
                print_ERRORmessage( "ERROR: Non-existent file \'" + gzffcal + "\'")
                self.bottombarlabelVariable1.set("ERROR: Non-existent file.")
                self.bottombarlabelVariable2.set("Select product: \"CALIBRATION\" to produce flat flux calibration plot.")
            pass

        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
##########################################

##########################################
########### CHECK FINAL PLOTS ############          
##########################################

############## RAW PLOT ##################            
    def OnRawPlotButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            # command line below creates "/products" dir and copies the final products with extention *.spc.gz to this dir                      
            FinalProdDir = self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() 
            if not os.path.exists(FinalProdDir):
                os.system("mkdir -p "  + FinalProdDir) 
                print_message("NEW directory with only final products " + FinalProdDir + " has been produced.")  
            spc_dir = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() + "/*.spc.gz")
            gz_spc = FinalProdDir + "/*spc.gz"
            final_dir = glob.glob(gz_spc)
            add_files = list(set(spc_dir) - set(final_dir))
            for gzfile in add_files:  
                os.system("cp " + gzfile + " " + self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() + "/")
            # check if dir containing just final products exists and uncompress files
            if not os.listdir(FinalProdDir) == []:
                os.system("gunzip -f %s " % gz_spc)   
                # command line below prints out information on data available in nightdir  
                nightdir = self.ddentryVariable.get() + self.nightentryVariable.get()+"/"
                commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE-OBS SGH-OIC IMAGETYP OBJECT"' 
                print "    FILENAME\t\t   DATE\t     IOD.CELL*  TYPE  OBJECT_NAME"   
                os.system(commandline)
                print_message( "*1 = with, 2 = without")
                # choose file from product directory with OPERA extension *.spc
                chooseFile = tkFileDialog.askopenfilename(title='Choose a file', defaultextension='.spc', filetypes =[("SPC extension",".spc"),("all files",".*")], initialdir=FinalProdDir + "/") 
                if chooseFile:
                    FileName = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() +"/*.spc")
                    os.system("cd %s" % FinalProdDir)  
                    openFile = open(chooseFile)
                    name = chooseFile.split('/')[-1]
                    print "Printing '%s'... " % chooseFile               
                    if openFile != None:
                        data = openFile.read()                
                        openFile.close()
                    os.system("echo 'set title \"Raw spectrum reduced by OPERA\" \n set xlabel \"Wavelength (nm)\" \n set ylabel \"Raw flux\" \n  plot \"%s\" u 5:9 lt -1 t \"%s\"' |  gnuplot -p"  % (chooseFile, name)) 
                    print ("Raw spectrum of \"%s\" has been plotted." % name)    
                else:
                    print_ERRORmessage( "ERROR: No file has been chosen."  )
                    self.bottombarlabelVariable1.set("ERROR: No file has been chosen.")
                    self.bottombarlabelVariable2.set("")    
                pass           
            else:        
                print_ERRORmessage( "ERROR: No files in directory. Run pipeline - select product: \"OBJECT\"." )
                self.bottombarlabelVariable1.set("ERROR: No files in directory " + FinalProdDir)
                self.bottombarlabelVariable2.set("Select product: \"OBJECT\" to produces final spectrum files.")
            pass
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################     

######### FLAT FIELD CALIB SPECTRUM ###########            

    def OnFlatPlotButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")            
            
            FinalProdDir = self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() 
            if not os.path.exists(FinalProdDir):
                os.system("mkdir -p "  + FinalProdDir) 
                print_message("NEW directory with only final products " + FinalProdDir + " has been produced.")  
            spc_dir = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() + "/*.spc.gz")
            gz_spc = FinalProdDir + "/*spc.gz"
            final_dir = glob.glob(gz_spc)
            add_files = list(set(spc_dir) - set(final_dir))
            for gzfile in add_files:  
                os.system("cp " + gzfile + " " + self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() + "/")
            if not os.listdir(FinalProdDir) == []:
                os.system("gunzip -f %s " % gz_spc)   
                
                # command line below prints out information on data available in nightdir  
                nightdir = self.ddentryVariable.get() + self.nightentryVariable.get()+"/"
                commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE-OBS SGH-OIC IMAGETYP OBJECT"' 
                print "    FILENAME\\t\\t   DATE\\t     IOD.CELL*  TYPE  OBJECT_NAME"   
                os.system(commandline)
                print_message( "*1 = with, 2 = without")
                # choose file from product directory with OPERA extension *.spc
                chooseFile = tkFileDialog.askopenfilename(title='Choose a file', defaultextension='.spc', filetypes =[("SPC extension",".spc"),("all files",".*")], initialdir=FinalProdDir + "/") 
                if chooseFile:
                    FileName = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() +"/*.spc")
                    os.system("cd %s" % FinalProdDir)  
                    openFile = open(chooseFile)
                    name = chooseFile.split('/')[-1]
                    print "Printing '%s'... " % chooseFile               
                    if openFile != None:
                        data = openFile.read()                
                        openFile.close()
                    os.system("cd %s" % FinalProdDir)
                    os.system("echo 'set title \"Calibrated spectrum reduced by OPERA\" \n set xlabel \"Wavelength (nm)\" \n set ylabel \"Raw flux / Flat-field (arbitrary units)f\" \n plot \"%s\" u 5:13 lt -1 t \"%s\"' |  gnuplot -p"  % (chooseFile, name)) 
                    print ("Calibrated spectrum of \"%s\" has been plotted." % name)
                else:
                    print_ERRORmessage( "ERROR: No file has been chosen."  )
                    self.bottombarlabelVariable1.set("ERROR: No file has been chosen.")
                    self.bottombarlabelVariable2.set("")    
                pass           

            else:        
                print_ERRORmessage( "ERROR: No files in directory. Run pipeline - select product: \"OBJECT\"." )
                self.bottombarlabelVariable1.set("ERROR: No files in directory " + FinalProdDir)
                self.bottombarlabelVariable2.set("Select product: \"OBJECT\" to produces final spectrum files.")
            pass
                           
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################

########### NORMALIZED SPECTRUM PLOT ##########            

    def OnNormButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")

            FinalProdDir = self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get()
            if not os.path.exists(FinalProdDir):
                os.system("mkdir -p "  + FinalProdDir) 
                print_message("NEW directory with only final products " + FinalProdDir + " has been produced.")  
            spc_dir = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() + "/*.spc.gz")
            gz_spc = FinalProdDir + "/*spc.gz"
            final_dir = glob.glob(gz_spc)
            add_files = list(set(spc_dir) - set(final_dir))
            for gzfile in add_files:  
                os.system("cp " + gzfile + " " + self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() + "/")
            if not os.listdir(FinalProdDir) == []:
                os.system("gunzip -f %s " % gz_spc)   
                
                # command line below prints out information on data available in nightdir  
                nightdir = self.ddentryVariable.get() + self.nightentryVariable.get()+"/"
                commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE-OBS SGH-OIC IMAGETYP OBJECT"' 
                print "    FILENAME\t\t   DATE\t     IOD.CELL*  TYPE  OBJECT_NAME"   
                os.system(commandline)
                print_message( "*1 = with, 2 = without")
                # choose file from product directory with OPERA extension *.spc
                chooseFile = tkFileDialog.askopenfilename(title='Choose a file', defaultextension='.spc', filetypes =[("SPC extension",".spc"),("all files",".*")], initialdir=FinalProdDir + "/") 
                if chooseFile:
                    FileName = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() +"/*.spc")
                    os.system("cd %s" % FinalProdDir)  
                    openFile = open(chooseFile)
                    name = chooseFile.split('/')[-1]
                    print "Printing '%s'... " % chooseFile               
                    if openFile != None:
                        data = openFile.read()                
                        openFile.close()
                    os.system("echo 'set title \"Normalized spectrum reduced by OPERA\" \n set xlabel \"Wavelength (nm)\" \n set ylabel \"Normalized flux\" \n plot \"%s\" u 5:11 lt -1 t \"%s\"' |  gnuplot -p"  % (chooseFile, name)) 
                    print ("Normalized spectrum of \"%s\" has been plotted." % name)  
                else:
                    print_ERRORmessage( "ERROR: No file has been chosen." )
                    self.bottombarlabelVariable1.set("ERROR: No file has been chosen.")
                    self.bottombarlabelVariable2.set("")    
                pass           

            else:        
                print_ERRORmessage( "ERROR: No files in directory. Run pipeline - select product: \"OBJECT\"." )
                self.bottombarlabelVariable1.set("ERROR: No files in directory " + FinalProdDir)
                self.bottombarlabelVariable2.set("Select product: \"OBJECT\" to produces final spectrum files.")
            pass   

        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################

############## 1D PLOT ##################            
    def OnONEDButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            # command line below creates "/products" dir and copies the final products with extention *.spc.gz to this dir                      
            ONEdProdDir = self.pdentryVariable.get() + self.finalProdDirVariable.get() + self.nightentryVariable.get() + "_ONED/"
            if not os.path.exists(ONEdProdDir):
                os.system("mkdir -p "  + ONEdProdDir ) 
                print_message("NEW directory with only final products " + ONEdProdDir + " has been produced.")  
            spc_dir = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() + "/*1d.spc")
            gz_spc = ONEdProdDir + "/*1d.spc"
            add_files = list(set(spc_dir) - set(gz_spc))
            for oneDfile in add_files:  
                os.system("cp " + oneDfile + " " + ONEdProdDir+ "/")
            # check if dir containing just final products exists and uncompress files
            if not os.listdir(ONEdProdDir) == []:
  
                # command line below prints out information on data available in nightdir  
                nightdir = self.ddentryVariable.get() + self.nightentryVariable.get()+"/"
                commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE-OBS SGH-OIC IMAGETYP OBJECT"' 
                print "    FILENAME\t\t   DATE\t     IOD.CELL*  TYPE  OBJECT_NAME"   
                os.system(commandline)
                print_message( "*1 = with, 2 = without")
                # choose file from product directory with OPERA extension *.spc
                chooseFile = tkFileDialog.askopenfilename(title='Choose a file', defaultextension='.spc', filetypes =[("SPC extension",".1d.spc"),("all files",".*")], initialdir=ONEdProdDir + "/") 
                if chooseFile:
                    FileName = glob.glob(self.pdentryVariable.get() + self.nightentryVariable.get() +"/*.spc")
                    os.system("cd %s" % ONEdProdDir)  
                    openFile = open(chooseFile)
                    name = chooseFile.split('/')[-1]
                    print "Printing '%s'... " % chooseFile               
                    if openFile != None:
                        data = openFile.read()                
                        openFile.close()
                    os.system("echo 'set title \"1D spectrum reduced by OPERA\" \n set xlabel \"Wavelength (nm)\" \n set ylabel \"Normalized flux\" \n  plot \"%s\" u 1:2 lt -1 t \"%s\"' |  gnuplot -p"  % (chooseFile, name)) 
                    print ("1D spectrum of \"%s\" has been plotted." % name)    
                else:
                    print_ERRORmessage( "ERROR: No file has been chosen."  )
                    self.bottombarlabelVariable1.set("ERROR: No file has been chosen.")
                    self.bottombarlabelVariable2.set("")    
                pass           
            else:        
                print_ERRORmessage( "ERROR: No files in directory. Run pipeline - select product: \"1D SPEC\"." )
                self.bottombarlabelVariable1.set("ERROR: No files in directory " + ONEdProdDir)
                self.bottombarlabelVariable2.set("Select product: \"1D SPEC\" to produces final spectrum files.")
            pass
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###############################################     


###############################################
############ OPEN FILE WITH DS9 ###############          
###############################################
    def Onds9ButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()
            print_infomessage("Checking the data in " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            
	    self.dataNightVariable = Tkinter.StringVar() 	
	    self.dataNightVariable.set(self.ddentryVariable.get() + self.nightentryVariable.get())
            print_message ("Choose a file from data directory: " + self.dataNightVariable.get())
               
            for file in glob.glob(self.ddentryVariable.get() + ".gz"):
		os.system('gunzip %s'  % file)
		print file   
              
            chooseFile = tkFileDialog.askopenfilename(title='Choose a file', defaultextension='.fit', filetypes =[("FITS extension",".fit"),("all files",".*")],initialdir=self.dataNightVariable.get()) 
            if chooseFile: 
                os.system("cd " +  self.dataNightVariable.get())
                openFile = open(chooseFile)
                print_infomessage ("Opening ds9 with file '%s' ... " % chooseFile )
                name = chooseFile.split('/')[-1]                     
                if openFile != None:
                     data = openFile.read()                
                     openFile.close()
                os.system("ds9 \"%s\"& "  % chooseFile)
                print_message("Image \'" + name + "\' was opened by Ds9.") 
            else:
                print_ERRORmessage( "ERROR: No data in directory.")
                self.bottombarlabelVariable1.set("ERROR: No data in directory " +  self.dataNightVariable.get())
                self.bottombarlabelVariable2.set("")
            pass
###########################################

###########################################
############# NIGHT LOG BUTTON ############        
###########################################
    def OnGenNightLogButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
        # test if all files exist
        nightdir = self.ddentryVariable.get()+self.nightentryVariable.get()+"/"
        if os.path.exists(nightdir) :    
            self.bottombarlabelVariable1.set("Generating Night Log for " + nightdir + " ...")
            self.bottombarlabelVariable2.set("** See data list output on terminal **")          
           
            print_infomessage("NIGHT LOG for DATA in DIR: " + nightdir)
            
           # command line below prints out information on data available in nightdir            
            commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE-OBS SGH-OIC IMAGETYP OBJECT"'
            print "    FILENAME\t\t   DATE\t     IOD.CELL*  TYPE  OBJECT_NAME"   
            os.system(commandline)
            print_message( "*1 = with, 2 = without")
        else :
            self.bottombarlabelVariable1.set("ERROR: non-existent data directory " + nightdir)
            self.bottombarlabelVariable2.set("")
        pass
###########################################            

###########################################
########### CONFIG SET UP #################
###########################################

############## DATA DIR ###################
    def OnPressEnterDDfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) :    
            self.bottombarlabelVariable1.set("DATA DIR is OK")
            self.bottombarlabelVariable2.set("")            
        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.ddentryVariable.get())
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
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.pdentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.pdentry.focus_set()
            self.pdentry.selection_range(0, Tkinter.END)                    
###########################################

############ Select PRODUCT ###############           
    def OnPressEnterCPfield(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.spentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Selected PRODUCT is OK")
            self.bottombarlabelVariable2.set("")
            
        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.spentryVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.spentry.focus_set()
            self.spentry.selection_range(0, Tkinter.END)                    
###########################################

############## Viewer ###################
    def OnPressViewerVariable(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.plotViewerVariable.get()) :    
            self.bottombarlabelVariable1.set("DATA DIR is OK")
            self.bottombarlabelVariable2.set("")            
        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.plotViewerVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.vventry.focus_set()
            self.vventry.selection_range(0, Tkinter.END)        
###########################################

############## Readout speed ##############
    def OnPressReadoutVariable(self,event):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.readoutVariable.get()) :    
            self.bottombarlabelVariable1.set("DATA DIR is OK")
            self.bottombarlabelVariable2.set("")            
        else:
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.readoutVariable.get())
            self.bottombarlabelVariable2.set("")            
            self.vventry.focus_set()
            self.vventry.selection_range(0, Tkinter.END)        
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
            self.bottombarlabelVariable1.set("ERROR: non-existent directory: " + self.ddentryVariable.get() + "/" + self.nightentryVariable.get())
            self.bottombarlabelVariable2.set("")

            self.nightentry.focus_set()
            self.nightentry.selection_range(0, Tkinter.END)     

###### CREATE NIGHT DIRECTORY #############          
    def OnCreateNightButtonClick(self):
        # check if config is loaded
        if self.fileloaded != True:
            print_ERRORmessage("ERROR: Load config.")
            self.bottombarlabelVariable1.set("ERROR: No config loaded.")
            self.bottombarlabelVariable2.set("") 
            return
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
    app.title('OES Pipeline v. 1.0')
    app.mainloop()

