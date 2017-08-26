#! /bin/bash
#
# test 1 amp to 2 amp conversion
#
bindir=$HOME/opera-1.0/bin/
upenabindir=$HOME/upena-1.1/
byproductsdir=/data/uhane5/opera/byproducts
rootdir=/data/niele/espadons/
files="
1310030c.fits
1310031f.fits
1310032f.fits
1310033f.fits
1310034f.fits
1310035f.fits
1310036f.fits
1310037f.fits
1310038f.fits
1310039f.fits
1310040f.fits
1310041f.fits
1310042f.fits
1310043f.fits
1310044f.fits
1310045f.fits
1310046f.fits
1310047f.fits
1310048f.fits
1310049f.fits
1310050f.fits
1310051a.fits
1310052b.fits
1310053b.fits
1310054b.fits
1310058o.fits
1310059o.fits
1310060o.fits
1310061o.fits
1310062o.fits
1310063o.fits
1310064o.fits
1310065o.fits
1310066o.fits
1310067o.fits
1310068o.fits
1310069o.fits
"
#############################################################################
# Algorithm:
# calculate the gains from the b.fits and f.fits images
# apply the gain ratio to all images (calibrations and objects)
# DO NOT measure the bias for each gain-modified b.fits in each mode
# subtract the bias difference from each object file
#############################################################################

# where the data resides
NIGHT=11AQ12-Jun16
# where the gain-corrected data resides
TMPNIGHT=01AQ12-Jun16
# here the gain-corrected, bias corrected data resides
TONIGHT=00AQ12-Jun16
# cleanup
#${bindir}/opera NIGHT=$NIGHT cleanit WHAT=gain
${bindir}/opera NIGHT=$NIGHT cleanit WHAT=bias
rm -f $TMPNIGHT/*
rm -f $TONIGHT/*
# calculate the gains for each mode
#echo opera NIGHT=$NIGHT gains -v
#${bindir}/opera NIGHT=$NIGHT gains -v
echo opera NIGHT=$NIGHT biases -v
${bindir}/opera NIGHT=$NIGHT biases -v
# apply the gain mods to all images
for f in $files 
do
	echo operaConvert2ampTo1amp --outputdir=$rootdir/$TONIGHT/ --images=$rootdir/$NIGHT/$f -v
	${bindir}/operaConvert2ampTo1amp --outputdir=$rootdir/$TONIGHT/ --images=$rootdir/$NIGHT/$f -v
done
exit
gainA_Normal=1.1548
gainB_Normal=1.1404
for f in $files 
do
	speed=`${upenabindir}/espgetspeed $rootdir/$NIGHT/$f`
	if [[ "$speed" == "Normal" ]]
	then
		echo operaConvert2ampTo1amp --gaina=$gainA_Normal --gainb=$gainB_Normal --outputdir=$rootdir/$TMPNIGHT/ --images=$rootdir/$NIGHT/$f -v
		${bindir}/operaConvert2ampTo1amp --gaina=$gainA_Normal --gainb=$gainB_Normal --outputdir=$rootdir/$TMPNIGHT/ --images=$rootdir/$NIGHT/$f -v
	else
		echo "$f not reduced as we only have data for speed Normal"
	fi
done
# calculate the bias difference
echo opera NIGHT=$NIGHT biases -v
${bindir}/opera NIGHT=$NIGHT biases -v
exit
# apply the bias difference to each image
for f in $files 
do
	mode=`${upenabindir}/espgetmode $rootdir/$NIGHT/$f`
	speed=`${upenabindir}/espgetspeed $rootdir/$NIGHT/$f`
	if [[ "$speed" == "Normal" ]]
	then
		biases=`cat <${byproductsdir}/$NIGHT/OLAPAab_${mode}_${speed}.bias`
		biasA=`${bindir}/operagetword 1 $biases`
		biasB=`${bindir}/operagetword 2 $biases`
		echo operaConvert2ampTo1amp --noshift --biasa=$biasA --biasb=$biasB --outputdir=$rootdir/$TONIGHT/ --images=$rootdir/$TMPNIGHT/$f -v
		${bindir}/operaConvert2ampTo1amp --noshift --biasa=$biasA --biasb=$biasB --outputdir=$rootdir/$TONIGHT/ --images=$rootdir/$TMPNIGHT/$f -v
	else
		echo "$f not reduced as we only have data for speed Normal"
	fi
done
exit


