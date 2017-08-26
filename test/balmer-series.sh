#!/bin/csh
# script to reduce ESPaDOnS data
# Eder Martioli
# CFHT Corporation, Jul 2013
#
echo "Running OPERA pipeline..."
echo " "

set INPUTEXTRACTEDSPECTRUM="$argv[1]"
set INPUTWAVEFILE="$argv[2]"
set OUTPUTROOTNAME="$argv[3]"

set GNUFILENORM="balmer-series.gnu"
set GNUFILE="balmer-series_continuum.gnu"

echo "reset" > $GNUFILENORM
echo "set xlabel '{/Symbol l} - {/Symbol l}_c (nm)'" >> $GNUFILENORM
echo "set ylabel 'relative flux + cte'" >> $GNUFILENORM
echo "set xrange[-5:5]" >> $GNUFILENORM
echo "set yrange[0.5:3.6]" >> $GNUFILENORM

echo "set title 'Balmer Series'" >> $GNUFILENORM

echo "Reducing the Balmer Series: 656 nm, 486 nm, 434 nm, 410 nm, 397 nm, 389 nm, and 383 nm"

####### H-alpha ######
set WAVELENGTH="656.3"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"  

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND

echo "plot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" > $GNUFILE

echo "set label 'H-alpha ($WAVELENGTH nm)' at 2,1.1" >> $GNUFILENORM
echo "plot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10) notitle w d" >> $GNUFILENORM
#################################

####### H-Beta ######
set WAVELENGTH="486.1"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"  

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND

echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-beta ($WAVELENGTH nm)' at 2,1.6" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+0.5) notitle w d" >> $GNUFILENORM
#################################

####### H-gamma ######
set WAVELENGTH="434.1"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"  

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND

echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-gamma ($WAVELENGTH nm)' at 2,2.0" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+0.9) notitle w d" >> $GNUFILENORM
#################################

####### H-delta ######
set WAVELENGTH="410.2"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"  

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND

echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-delta ($WAVELENGTH nm)' at 2,2.4" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+1.3) notitle w d" >> $GNUFILENORM
#################################

####### H-epsilon ######
set WAVELENGTH="397.0"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"   

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND

echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-epsilon ($WAVELENGTH nm)' at 2,2.8" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+1.7) notitle w d" >> $GNUFILENORM
#################################

####### H-zeta ######
set WAVELENGTH="388.9"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"   

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND
echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-zeta ($WAVELENGTH nm)' at 2,3.15" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+2.1) notitle w d" >> $GNUFILENORM
#################################

####### H-eta ######
set WAVELENGTH="383.5"
set BINSIZE="250"
set ORDERBIN="2"

set OUTPUTSPECTRUMFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".s.gz"  
set SPECTRUMDATAFILENAME=$OUTPUTROOTNAME"_"$WAVELENGTH".spec"   

set COMMAND = "operaNormalizeAcrossOrders --inputUncalibratedSpectrum=$INPUTEXTRACTEDSPECTRUM --inputWaveFile=$INPUTWAVEFILE --outputNormalizedSpectrum=$OUTPUTSPECTRUMFILENAME --binsize=$BINSIZE --orderBin=$ORDERBIN --nsigcut=3 --wavelength=$WAVELENGTH --spectrumDataFilename=$SPECTRUMDATAFILENAME"

echo "$COMMAND"
$COMMAND
echo "replot '$SPECTRUMDATAFILENAME' u 7:8 w d, ''  u 7:9 w l" >> $GNUFILE

echo "set label 'H-eta ($WAVELENGTH nm)' at 2,3.4" >> $GNUFILENORM
echo "replot '$SPECTRUMDATAFILENAME' u ("'$'"7-$WAVELENGTH):("'$'"10+2.4) notitle w d" >> $GNUFILENORM

set PLOTFILENAME=$OUTPUTROOTNAME"_BalmerSeries.eps"   

echo 'set terminal postscript enhanced color solid lw 1.5 "Helvetica" 14' >> $GNUFILENORM
echo 'set output "'"$PLOTFILENAME" >> $GNUFILENORM
echo "replot" >> $GNUFILENORM
echo "set terminal x11" >> $GNUFILENORM
echo "set output" >> $GNUFILENORM

#################################

gnuplot -persist $GNUFILENORM

exit