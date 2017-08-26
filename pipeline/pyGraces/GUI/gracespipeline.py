#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
"""
Created on Wed May 28 21:36:03 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Jul 30, 2014
"""

import Tkinter
import os 

class pipelineFrontEnd_tk(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.grid()
        
########### Abort button ##################   
        # Below is supposed to be a button to abort execution, but so far it just closes python
        self.abortButton = Tkinter.Button(self,text="Exit", state="active", command=self.OnAbortButtonClick)        
        self.abortButton.grid(column=3,row=8,columnspan=1, sticky='EW')
###########################################        
        
########### CONFIG FILE ##################   
        # Below it initializes the string for opera system directory       
        self.operaDirVariable = Tkinter.StringVar()
        self.operaDirVariable.set("/home/user/opera-1.0")
        # Below is the button to load CONFIG FILE
        loadConfigButton = Tkinter.Button(self,text="Load Config", command=self.LoadConfigButtonClick)        
        loadConfigButton.grid(column=0,row=0)  
        # Below is the text box for config file path entry
        self.filepathVariable = Tkinter.StringVar()
        self.filepathVariable.set(self.operaDirVariable.get() + "/GUI/config.gracespipeline")
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
        self.R1 = Tkinter.Radiobutton(self, text="Star Only", variable=self.instmode,value=1, anchor="e")
        self.R1.grid(column=1,row=3,columnspan=1)
        
        self.R2 = Tkinter.Radiobutton(self, text="Star+Sky", variable=self.instmode, value=2, anchor="e")
        self.R2.grid(column=2,row=3,columnspan=1)

        #self.R3 = Tkinter.Radiobutton(self, text="Polar", variable=self.instmode, value=3, anchor="e", state="disabled")
        #self.R3.grid(column=3,row=3,columnspan=1)
        #self.R4 = Tkinter.Radiobutton(self, text="Any", variable=self.instmode, value=0, anchor="e")
        #self.R4.grid(column=4,row=3,columnspan=1)
        
        # Default is Star-Only FOURSLICE mode
        self.R1.select();
                
        # Below it sets the header keyword value for chosen mode
        self.instmodekeyword = Tkinter.StringVar() 
        if (self.instmode.get() == 1) :
            self.instmodekeyword.set("FOURSLICE")
        elif (self.instmode.get() == 2) :
            self.instmodekeyword.set("TWOSLICE")
        
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
        self.rb1.grid(column=1,row=4,columnspan=1)

        self.rb2 = Tkinter.Radiobutton(self, text="Normal", variable=self.readoutmode, value=2)
        self.rb2.grid(column=2,row=4,columnspan=1)

        self.rb3 = Tkinter.Radiobutton(self, text="Slow", variable=self.readoutmode, value=3)
        self.rb3.grid(column=3,row=4,columnspan=1)

        # Default is Normal readout speed
        self.rb2.select();
        
        self.anyreadoutcheck = Tkinter.IntVar()
        self.cb_anyreadout = Tkinter.Checkbutton(self, text="Allow any", variable=self.anyreadoutcheck, onvalue = 1, offvalue = 0)
        self.cb_anyreadout.grid(column=4,row=4,columnspan=1)
        self.cb_anyreadout.select();
        
        # Below it sets the header keyword value for chosen mode
        self.readspeedkeyword = Tkinter.StringVar()      
        if (self.readoutmode.get() == 1) :
            self.readspeedkeyword.set("Fast: 4.70e noise, 1.60e/ADU, 32s")
        elif (self.readoutmode.get() == 2) :
            self.readspeedkeyword.set("Normal: 4.20e noise, 1.30e/ADU, 38s")         
        elif (self.readoutmode.get() == 3) :
            self.readspeedkeyword.set("Slow: 2.90e noise, 1.20e/ADU, 60s") 
            
###########################################

############# REDUCTION STEPS ###############           
        # Below is the sidetext for REDUCTION STEPS
        self.sidetextpd4 = Tkinter.StringVar()        
        sidetextlabelpd4 = Tkinter.Label(self,textvariable=self.sidetextpd4, anchor="w",fg="black",bg="lightgray")
        sidetextlabelpd4.grid(column=0,row=5,columnspan=1, sticky='EW')        
        self.sidetextpd4.set("REDUCTION STEPS:")
        
        # Below are the check buttons to select REDUCTION STEPS
        self.calibrationcheck = Tkinter.IntVar()
        self.cb1 = Tkinter.Checkbutton(self, text="Calibration", variable=self.calibrationcheck, onvalue = 1, offvalue = 0)
        self.cb1.grid(column=1,row=5,columnspan=1)
        self.cb1.select();

        self.reductioncheck = Tkinter.IntVar()
        self.cb2 = Tkinter.Checkbutton(self, text="Reduction", variable=self.reductioncheck, onvalue = 1, offvalue = 0)
        self.cb2.grid(column=2,row=5,columnspan=1)
###########################################
     
################# NIGHT ###################           
       # Below is the sidetext for NIGHT
        self.sidetextpd5 = Tkinter.StringVar()        
        sidetextlabelpd5 = Tkinter.Label(self,textvariable=self.sidetextpd5, anchor="w",fg="white",bg="darkblue")
        sidetextlabelpd5.grid(column=0,row=6,columnspan=1)        
        self.sidetextpd5.set("NIGHT:")
        # Below is the text box for NIGHT string entry
        self.nightentryVariable = Tkinter.StringVar()
        self.nightentry = Tkinter.Entry(self,textvariable=self.nightentryVariable)
        self.nightentry.grid(column=1,row=6,columnspan=2, sticky='EW')
        self.nightentry.bind("<Return>", self.OnPressEnterNightfield)
        # Below is the button to CREATE REDUCTION NIGHT DIR
        self.createNightButton = Tkinter.Button(self,text=u"Create Night DIR", command=self.OnCreateNightButtonClick, state="normal")        
        self.createNightButton.grid(column=3,row=6,columnspan=1)
###########################################            
      
############# DEFAULT CALIBRATION ###############         
        # Below it initializes the string for DEFAULT CALIBRATION directory       
        self.defaultCalDirVariable = Tkinter.StringVar()        
        self.defaultCalDirVariable.set("DEFAULT CALIBRATION dir: /home/usr/gracespipeline-1.0/DefaultCalibration/")

        # Below is the check button to allow use of DEFAULT CALIBRATION
        self.defaultcalibrationcheck = Tkinter.IntVar()
        self.cb4 = Tkinter.Checkbutton(self, text="Allow use of DEFAULT Calibration", variable=self.defaultcalibrationcheck, onvalue = 1, offvalue = 0)
        self.cb4.grid(column=1,row=7,columnspan=2)
        #self.cb4.select();
###########################################
    
    
############### NIGHT LOG  ################        
        # Below is the button to generate NIGHT LOG
        nightLogButton = Tkinter.Button(self,text=u"Gen Night Log", command=self.OnGenNightLogButtonClick)        
        nightLogButton.grid(column=0,row=8,columnspan=1)
###########################################      
   
############### CHECK DATA ################        
        # Below is the button to CHECK DATA
        checkDataButton = Tkinter.Button(self,text=u"Check Data", command=self.OnCheckDataButtonClick)        
        checkDataButton.grid(column=1,row=8,columnspan=1)
###########################################      
        
############### RUN GRACES PIPELINE #################        
        # Below is the button to run GRACES PIPELINE
        self.runOperaButton = Tkinter.Button(self,text=u"Run Pipeline", state="active", command=self.OnRunOperaButtonClick)        
        self.runOperaButton.grid(column=2,row=8,columnspan=1, sticky='EW')
###########################################               

################# CONSOLE #################        
        # Strings for console output
        self.headbottombarlabelVariable = Tkinter.StringVar()        
        self.bottombarlabelVariable1 = Tkinter.StringVar()  
        self.bottombarlabelVariable2 = Tkinter.StringVar()  
        self.bottombarlabelVariable3 = Tkinter.StringVar()  

        # Construct labels for console 
        headlabel = Tkinter.Label(self,textvariable=self.headbottombarlabelVariable, anchor="c",fg="darkgrey",bg="lightgray")
        self.headbottombarlabelVariable.set("GRACES Pipeline 1.0 Console by E. Martioli")
        headlabel.grid(column=0,row=9,columnspan=5, sticky='EW')   
        
        label1 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable1, anchor="c",fg="black",bg="lightgray")        
        label2 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable2, anchor="c",fg="black",bg="lightgray")
        #label3 = Tkinter.Label(self,textvariable=self.bottombarlabelVariable3, anchor="c",fg="black",bg="lightgray")

        self.bottombarlabelVariable2.set(self.operaDirVariable.get())
             
        label1.grid(column=0,row=10,columnspan=5, sticky='EW')        
        label2.grid(column=0,row=11,columnspan=5, sticky='EW')        
        #label3.grid(column=0,row=12,columnspan=5, sticky='EW')               
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
               
########### OPERA BUTTON ##################   
    def OnRunOperaButtonClick(self):
        # test if directory exists if not it will ask for it.
        if os.path.exists(self.ddentryVariable.get()) and os.path.exists(self.pdentryVariable.get()) :    
            self.bottombarlabelVariable1.set("Running GRACES Pipeline for " + self.ddentryVariable.get() + self.nightentryVariable.get() + "/ ...")
            self.bottombarlabelVariable2.set("")                        
            self.update()

            print " "                       
            print "-----------------------------------------------------------------"
            print "Running GRACES Pipeline for NIGHT: " + self.nightentryVariable.get()
            print "-----------------------------------------------------------------"
           # command line below prints out information on data available in nightdir            

            commandline = self.operaDirVariable.get() + '/scripts/graces_reduction.sh' \
            ' ' + self.operaDirVariable.get() + \
            ' ' + self.ddentryVariable.get() + \
            ' ' + self.pdentryVariable.get() + \
            ' ' + self.nightentryVariable.get() + \
            ' ' + str(self.instmode.get()) + \
            ' ' + str(self.readoutmode.get()) + \
            ' ' + str(self.calibrationcheck.get())  + \
            ' ' + str(self.reductioncheck.get()) + \
            ' ' + str(self.anyreadoutcheck.get()) + \
            ' ' + str(self.defaultcalibrationcheck.get())
            
            print "Issuing command: " + commandline
            os.system(commandline)
            print "-----------------------------------------------------------------"
            print " "  
            
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
            commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "DATE GSLICER EREADSPD OBSTYPE EXPTIME OBJECT" -p'   
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
        
        defaultmasterbias = ''
        defaultmasterflat = ''
        defaultmastercomp = ''
        
        # Below it sets the header keyword value for chosen mode
        if (self.instmode.get() == 1) :
            self.instmodekeyword.set("FOURSLICE")
            defaultmasterflat = self.operaDirVariable.get() + 'DefaultCalibration/StarOnly_masterflat.fits.gz'
            defaultmastercomp = self.operaDirVariable.get() + 'DefaultCalibration/StarOnly_mastercomp.fits.gz'
        elif (self.instmode.get() == 2) :
            self.instmodekeyword.set("TWOSLICE")
            defaultmasterflat = self.operaDirVariable.get() + 'DefaultCalibration/StarPlusSky_masterflat.fits.gz'
            defaultmastercomp = self.operaDirVariable.get() + 'DefaultCalibration/StarPlusSky_mastercomp.fits.gz'            
         
        # Below it sets the header keyword value for chosen mode
        if (self.readoutmode.get() == 1) :
            self.readspeedkeyword.set("Fast: 4.70e noise, 1.60e/ADU, 32s")
            defaultmasterbias = self.operaDirVariable.get() + 'DefaultCalibration/fast_masterbias.fits.gz'
        elif (self.readoutmode.get() == 2) :
            self.readspeedkeyword.set("Normal: 4.20e noise, 1.30e/ADU, 38s")    
            defaultmasterbias = self.operaDirVariable.get() + 'DefaultCalibration/normal_masterbias.fits.gz'            
        elif (self.readoutmode.get() == 3) :
            self.readspeedkeyword.set("Slow: 2.90e noise, 1.20e/ADU, 60s") 
            defaultmasterbias = self.operaDirVariable.get() + 'DefaultCalibration/slow_masterbias.fits.gz'            
            
        # test if all files exist
        nightdir = self.ddentryVariable.get()+self.nightentryVariable.get()+"/"
        if os.path.exists(nightdir) :    
            self.bottombarlabelVariable1.set("Checking available DATA for " + nightdir + " ...")
            self.bottombarlabelVariable2.set("** See data list output on terminal **")                      
            
           # start of command line            
            commandline = self.operaDirVariable.get() + '/bin/operaQueryImageInfo --directory=' + nightdir + ' -e "OBSTYPE GSLICER EREADSPD OBJECT"'   
            
            print " "                       
            print "-----------------------------------------------------------------"
            print "Checking DATA in DIR: " + nightdir
            print "Instrument Mode: " + self.instmodekeyword.get()
            print "Read out Speed: " + self.readspeedkeyword.get()
            print "-----------------------------------------------------------------"            
            
            if self.calibrationcheck.get() :
                # The three command lines below print out information only for  calibration data available in selected mode
                #biascheck = commandline + ' -q "INSTRUME EREADSPD" INSTRUME=GRACES' + ' EREADSPD="' + self.readspeedkeyword.get() + '" | grep "N"'+ self.nightentryVariable.get() +'"B"'
                #flatcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD" INSTRUME=GRACES GSLICER=' + self.instmodekeyword.get() +' EREADSPD="' + self.readspeedkeyword.get() + '" | grep "N"'+ self.nightentryVariable.get() +'"F"'   
                #compcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD" INSTRUME=GRACES GSLICER=' + self.instmodekeyword.get() +' EREADSPD="' + self.readspeedkeyword.get() + '" | grep "N"'+ self.nightentryVariable.get() +'"C"' 
                
                biascheck = commandline + ' -q "INSTRUME EREADSPD OBSTYPE" INSTRUME=GRACES EREADSPD="' + self.readspeedkeyword.get() + '" OBSTYPE=BIAS'

                if self.anyreadoutcheck.get() :                
                    flatcheck = commandline + ' -q "INSTRUME GSLICER OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" OBSTYPE=FLAT'   
                    compcheck = commandline + ' -q "INSTRUME GSLICER OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" OBSTYPE=COMPARISON'   
                else:
                    flatcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" EREADSPD="' + self.readspeedkeyword.get() + '" OBSTYPE=FLAT'   
                    compcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" EREADSPD="' + self.readspeedkeyword.get() + '" OBSTYPE=COMPARISON'   
                    
                print " "                       
                print "BIAS files available for selected mode:"   
                print "---------------------------------------"
                #print biascheck                     
                if self.defaultcalibrationcheck.get() :
                   print defaultmasterbias
                                            
                os.system(biascheck)              
                print " "                       
                print "FLAT files available for selected mode:"            
                print "---------------------------------------" 
                #print flatcheck                     
                if self.defaultcalibrationcheck.get() :
                   print defaultmasterflat                
                
                os.system(flatcheck)
                print " "                       
                print "COMPARISON files available for selected mode:"            
                print "---------------------------------------------"     
                #print compcheck                                 
                if self.defaultcalibrationcheck.get() :
                   print defaultmastercomp
                   
                os.system(compcheck) 
                            
            if self.reductioncheck.get() :   
                # The command lines below prints out information only for object data available in selected mode                
                if self.anyreadoutcheck.get() :                
                    #objcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=GRACES GSLICER=' + self.instmodekeyword.get() +' EREADSPD="' + self.readspeedkeyword.get() + '" | grep "N"'+ self.nightentryVariable.get() +'"O"'                       
                    objcheck = commandline + ' -q "INSTRUME GSLICER OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" OBSTYPE=OBJECT'    
                else:
                    objcheck = commandline + ' -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=GRACES GSLICER="' + self.instmodekeyword.get() +'" EREADSPD="' + self.readspeedkeyword.get() + '" OBSTYPE=OBJECT'       
            
                print " "                       
                print "OBJECT files available for selected mode:"            
                print "---------------------------------------------"     
                #print objcheck
                os.system(objcheck)      
            print " "                       
                            
        else :
            self.bottombarlabelVariable1.set("Error: non-existent data directory " + nightdir + "...")
            self.bottombarlabelVariable2.set("")
        pass
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
                if(col0=="GRACESPIPELINEDIR") :
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
            self.bottombarlabelVariable2.set("GRACES Pipeline home dir:" + self.operaDirVariable.get())
            
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
                if(col0=="GRACESPIPELINEDIR") :
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
            self.bottombarlabelVariable2.set("GRACES Pipeline home dir:" + self.operaDirVariable.get())
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
    app.title('GRACES Pipeline v. 1.0')
    app.mainloop()
