#!/bin/csh
# Description: script to fix header keywords of GRACES commissioning data
# Author: Eder Martioli
# Laboratorio Nacional de Astrofisica, Brazil
# July 2014
#

set DATAROOTDIR="/data/GRACES/"
set SETKEYEXE="/Users/edermartioli/opera-1.0/pipeline/pyGraces/setGracesKeywords.py"

####### NIGHT dir ########
set NIGHT="20140505" 
##########################

ls $DATAROOTDIR/$NIGHT/*b.fits > "bias.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*f.fits > "flat.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*c.fits > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*o.fits > "object.$NIGHT.list"

$SETKEYEXE -i "bias.$NIGHT.list" --instmode=starsky --readspeed=fast --obstype=bias
$SETKEYEXE -i "flat.$NIGHT.list" --instmode=starsky --readspeed=fast --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --instmode=starsky --readspeed=fast --obstype=arc
$SETKEYEXE -i "object.$NIGHT.list" --instmode=starsky --readspeed=fast --obstype=object

####### NIGHT dir ########
set NIGHT="20140506" 
##########################

ls $DATAROOTDIR/$NIGHT/*b.fits > "bias.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*f.fits > "flat.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*c.fits > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*o.fits > "object.$NIGHT.list"

$SETKEYEXE -i "bias.$NIGHT.list" --obstype=bias
$SETKEYEXE -i "flat.$NIGHT.list" --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --obstype=arc
$SETKEYEXE -i "object.$NIGHT.list" --obstype=object

####### NIGHT dir ########
set NIGHT="20140513" 
##########################

ls $DATAROOTDIR/$NIGHT/*C*.fits > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*F*.fits > "flat.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*O*.fits > "object.$NIGHT.list"

# Flats and arcs were taken in fast mode, therefore the --readspeed option artificially 
# changes the readout in the header to allow the pipeline to identify flats and arcs for calibration

$SETKEYEXE -i "flat.$NIGHT.list" --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --obstype=arc
$SETKEYEXE -i "object.$NIGHT.list" --obstype=object

# Below we could borrow a few bias in slow mode from the night of 20140515 
#cp $DATAROOTDIR/20140515/*B001[2-8].fits $DATAROOTDIR/$NIGHT/

####### NIGHT dir ########
set NIGHT="20140514" 
##########################

ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C000[1-5].fits" > "flat.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C001[1-6].fits" >> "flat.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C002[2-6].fits" >> "flat.$NIGHT.list"

ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C000[6-9].fits" > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C0010.fits" >> "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C001[7-9].fits" >> "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*C*.fits | grep "N20140514C002[0-1].fits" >> "arc.$NIGHT.list"

$SETKEYEXE -i "flat.$NIGHT.list" --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --obstype=arc

####### NIGHT dir ########
set NIGHT="20140515" 
##########################

ls $DATAROOTDIR/$NIGHT/*C*.fits > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*F*.fits > "flat.$NIGHT.list"

$SETKEYEXE -i "flat.$NIGHT.list" --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --obstype=arc

####### NIGHT dir ########
set NIGHT="20140516" 
##########################

ls $DATAROOTDIR/$NIGHT/*C*.fits > "arc.$NIGHT.list"
ls $DATAROOTDIR/$NIGHT/*F*.fits > "flat.$NIGHT.list"

$SETKEYEXE -i "flat.$NIGHT.list" --obstype=flat
$SETKEYEXE -i "arc.$NIGHT.list" --obstype=arc

exit


