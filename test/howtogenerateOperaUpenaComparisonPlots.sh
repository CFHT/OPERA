#!/bin/bash
#
# Graphically compare upena spectra with opera equivalaents
#
# In order to compare LE with opera, a parallel directory structure is assumed.
# i,e, the directory: /data/uhane5/upena/spectra/10AQ00-Jul31/
# If the upena directory ds not exists then only the opera plot is done.
#
# Arguments:
#
#	-I=1	immediate plot
#	linestyle=impulses | dots etc
#	xrange='wl0:wlf' or 'd0:dn'
#	yrange='y0:yf'
#	order=n (only for e.eps, p.eps, pn.eps)
#	comment="..."
#	withtelluriclines is defined in Makefile.parameters
#	withtelluricspectrum is defined in Makefile.parameters
basedir="/data/uhane5/opera/"
basedir="/data/espadons/"
# opera extended spectrum, order 35 only, distance 0 - 2500, interactive plot
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221.e.eps order=35 -I=1 xrange='0:2500'

# opera polarimetry, wl range, (interactive plot not available for a multiplot such as .p)
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221.p.eps xrange='655.5:657'

# opera normalized polarimetry, wl range, (interactive plot not available for a multiplot)
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221.pol.eps type=normalized wl=1 xrange='655.5:657'

# LE-compatible normalized intensity spectra, wl range, interactive plot
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221in.s.eps xrange='655.5:657' -I=1

# LE-compatible normalized, wl-calibrated polarimetry, wl range
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221pnw.s.eps xrange='655.5:657'

# LE-compatible SNR, interactive plot
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221i.sn.eps -I=1

# opera extended spectrum, intensity, wl range, interactive plot, order 35 only
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221.spc.eps type=unnormalized wl=1 xrange='655.5:657' -I=1 order=35

# opera extended spectrum, normalized intensity, wl range, interactive plot
opera DATADIR=${basedir}/spectra/10AQ00-Jul31/ 1219221.spc.eps type=normalized wl=1 xrange='655.5:657' -I=1

# Svetlana Berdyugina
opera DATADIR=${basedir}/spectra/13AQ02-Feb28  1610460.pol.eps -I=1 order=39 type=normalized wl=1 xrange='588.8:589.6' comment="Svetlana Berdyugina"

# Normalized, beam spectrum
opera DATADIR=${basedir}/spectra/10AQ00-Jul31 1219221.spc.eps type=unnormalized wl=1 xrange='655.5:657' -I=1 order=35

# plot an image
opera /data/espadons/11AQ14-Jul08/1315170o.fits.eps -I=1
# green
opera /data/espadons/11AQ14-Jul08/1315170o.fits.eps -I=1 color=22,13,5
# color default palette
opera /data/espadons/11AQ14-Jul08/1315170o.fits.eps -I=1 color=1 
# overplot geometry
opera /data/espadons/11AQ14-Jul08/1315170o.fits.eps -I=1 geom=... order=35,44,50

exit

