
		 operaFits2txt Version 1.0
	  Canada France Hawaii Telescope
		      Eder Martioli 
			  
			   April 2011

operaFits2txt unpacks data stored in upena fits 
files to recover the original tabular data in text format.
The default operation is to process all of the fits files
on the command line.

Examples:

./operaFits2txt 1287849p.fits
Produces:
1287849pn.s 1287849pu.s 1287849p.sn

./operaFits2txt 1287849i.fits
Produces:
1287849inw.s  1287849iuw.s  1287849i.sn
1287849in.s   1287849iu.s 

Wildcards also work:
./operaFits2txt *p.fits
./operaFits2txt *i.fits

You may also add the column names with the -c option:

./operaFits2txt 1287849p.fits -c

operaFits2txt also works with Rice or gzip compressed 
images, as stored in an archive:

./operaFits2txt 1287849p.fits.fz
./operaFits2txt 1287849p.fits.gz

-------------------------------------
Warning: in order to use operaFits2txt with 
compressed files one needs to have either 
fpack/funpack or gzip/gunzip installed. 
--------------------------------------

operaFits2txt builds on MacOSX and Linux:

make operaFits2txt

-------------------------------------
Warning: operaFits2txt requires cfitsio
installed in the computer and optionally
fpack and/or gzip. 
--------------------------------------

Enjoy!
The CFHT OPERA Team
