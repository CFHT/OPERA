#! /bin/bash
###############################################################################
#
# test espadons quicklook
#
###############################################################################
bindir=$HOME/opera-1.0/bin/
night=$1
sessiondir=/data/niele/espadons/
indir=$sessiondir/$night/
outdir=/data/uhane5/opera/$night/
o_only=0
if [ "$night" == "" ]
then
	echo "usage: $(basename $0) <night directory> [-o] (o.fits only)"
	exit
fi
shift
if [ "$1" == "-o" ]
then
	o_only=1
fi
if [ ! -d $indir ]
then
	echo "usage: $(basename $0) $indir does not exist."
	exit
fi

echo "$(basename $0): working in directory $indir output to $outdir"
echo "$(basename $0): start quicklook like this: espqlh --night=$night -v"
echo "$(basename $0): in another terminal..."
cd $indir
if (( o_only == 1 ))
then
	files=`ls *o.fits`
else
	files=`ls *.fits *f.fits`
fi

for file in $files
do
	echo "$(basename $0): $file"
	rm -f ${sessiondir}current.fits
	ln -s ${indir}$file ${sessiondir}current.fits
	rm -f ${indir}current.fits
	ln -s ${indir}$file ${indir}current.fits
	sleep 20
done
echo "$(basename $0): done $night"
exit


