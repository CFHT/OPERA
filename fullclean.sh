make uninstall
rm lib/*.dylib
make clean
make distclean
rm -rf autom4te.cache
rm -f m4/*
rm -f aclocal.m4 config.guess config.sub configure depcomp install-sh ltmain.sh missing
find -name Makefile.in -type f -delete
find -name .DS_Store -type f -delete
