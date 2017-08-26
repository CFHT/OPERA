#! /bin/sh
# This script runs Autotools, for use before running ./configure.
# Only use if unable to obtain a version of autotools with autoreconf.
# Otherwise, "autoreconf -vfi" is a better alternative to this script.
(libtoolize || glibtoolize) \
&& aclocal \
&& autoconf \
&& automake --add-missing
