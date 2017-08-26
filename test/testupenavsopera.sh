#!/bin/bash
#
# Compare upena spectra with opera equivalaents
#
if (( $# == 0 ))
then
	echo" usage: $(basename $0) <NIGHT> [ --keep ]"
	exit
fi
operadir=/data/uhane5/opera/spectra/
upenadir=/data/uhane5/upena/spectra/
night=$1
keep=$2

#
# first make sure everything is there....
#
uspectra=`ls $upenadir/$night/*.s`
for spectrum in $uspectra 
do
	ospectrum=$(basename $spectrum)
	if [ ! -e $operadir/$night/$ospectrum.gz ]
	then
		echo "$ospectrum is missing from opera."
	fi
done
ospectra=`ls $operadir/$night/*.s.gz`
for spectrum in $ospectra 
do
	uspectrum=$(basename $spectrum .gz)
	if [ ! -e $upenadir/$night/$uspectrum ]
	then
		echo "$uspectrum is missing from upena."
	fi
done

#
# next check one by one
#
for spectrum in $ospectra 
do
	uspectrum=$(basename $spectrum .gz)
	ospectrum=$(basename $spectrum)
	if [ -e $upenadir/$night/$uspectrum ]
	then
		cd $operadir/$night/ >/dev/null
		if [ ! -e $uspectrum ]
		then
			cp -p $ospectrum $ospectrum.bak
			gunzip -f -d $ospectrum
			mv $ospectrum.bak $ospectrum
		fi
		cd - >/dev/null
		if [[ -e $operadir/$night/$uspectrum && -e $upenadir/$night/$uspectrum ]]
		then
			if [[ "`diff -w -b -q $operadir/$night/$uspectrum $upenadir/$night/`" != "" ]]
			then
				echo "$uspectrum differs"
				diff -w -b -s $operadir/$night/$uspectrum $upenadir/$night/ >$uspectrum.difflog
				if [[ "$keep" == "" ]]
				then
					rm -f $operadir/$night/$uspectrum
				fi
			else
				echo "$uspectrum OK"
				rm -f $operadir/$night/$uspectrum
			fi
		else
			echo "skipping $uspectrum"
		fi
	fi
done
exit

