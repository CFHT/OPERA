#########################################################################################
#
# Process this file with autoconf to produce a configure script.
# If you add directories, list them in AC_CONFIG_FILES
#
#########################################################################################

#                                               -*- Autoconf -*-

#AC_PREREQ([2.63])
AC_INIT([opera-pipeline], [1.0], [opera@cfht.hawaii.edu])
AM_INIT_AUTOMAKE([1.10 no-define foreign])
AC_REVISION([$Revision: 1.0 $])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_HEADER_STDC
AM_DISABLE_SHARED
#AC_PROG_RANLIB
AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AC_PREFIX_DEFAULT($HOME/opera-1.0)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h fitsio.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_ARG_ENABLE([LIBPNG],
	[ --enable-LIBPNG	Enable use of libpng in the plot tools ],
	[case "${enableval}" in
		true)	libpng=true ;;
		false)  libpng=false ;;
		yes)	libpng=true ;;
		no)		libpng=false ;;
		1)		libpng=true ;;
		0)		libpng=false ;;
		*)		AC_MSG_ERROR([bad value ${enableval} for --enable-LIBPNG]) ;;
	esac], [libpng=false])
	
AM_CONDITIONAL([LIBPNG], [test x$libpng = xtrue])

AC_ARG_ENABLE([LIBWCS],
	[ --enable-LIBWCS	Enable use of libwcs world coordinate system ],
	[case "${enableval}" in
		true)	libwcs=true ;;
		false)  libwcs=false ;;
		yes)	libwcs=true ;;
		no)		libwcs=false ;;
		1)		libwcs=true ;;
		0)		libwcs=false ;;
		*)		AC_MSG_ERROR([bad value ${enableval} for --enable-LIBWCS]) ;;
	esac], [libwcs=false])
	
AM_CONDITIONAL([LIBWCS], [test x$libwcs = xtrue])

#
# we need cfitsio
#
AC_ARG_WITH([cfitsio], [AS_HELP_STRING([--without-cfitsio],
            [eschew CFITSIO])], [])
if test "x$with_cfitsio" = xno ; then
  AC_MSG_WARN([CFITSIO disabled])
else
  AC_ARG_WITH([cfitsiolib], [AS_HELP_STRING([--with-cfitsiolib=DIR],
              [directory containing cfitsio library])], [])
  if test "x$with_cfitsiolib" != x ; then
    LIBDIRS="$LIBDIRS $with_cfitsiolib"
  fi

  AC_ARG_WITH([cfitsioinc], [AS_HELP_STRING([--with-cfitsioinc=DIR],
              [directory containing cfitsio header files])], [])
  if test "x$with_cfitsioinc" != x ; then
    CFITSIO_INCDIRS="$with_cfitsioinc"
  fi

  CFITSIO_INCDIRS="$CFITSIO_INCDIRS   \
           /usr/local/cfitsio/include \
           /local/cfitsio/include"

  CFITSIO_LIBDIRS="$LIBDIRS           \
           /usr/local/cfitsio/lib \
           /local/cfitsio/lib"
		   
  LIBDIRS="$LIBDIRS           \
           /usr/local/lib     \
           /lib               \
           /local/lib         \
           /opt/local/lib"

  # Generic include directories.
  INCDIRS="/usr/local/include \
           /local/include     \
           /include           \
           /opt/local/include \
           /usr/include"

  for LIBDIR in $CFITSIO_LIBDIRS $LIBDIRS ; do
    AC_CHECK_FILE([$LIBDIR], [CFITSIOLIB="$CFITSIOLIB -L$LIBDIR"], [continue])
  done

    # Search for CFITSIO.
    for INCDIR in $CFITSIO_INCDIRS $INCDIRS ; do
      AC_CHECK_FILE([$INCDIR/cfitsio/fitsio.h], [CFITSIOINC="-I$INCDIR/cfitsio"; break])
      AC_CHECK_FILE([$INCDIR/fitsio.h], [CFITSIOINC="-I$INCDIR"; break])
    done

    if test "x$CFITSIOINC" = x -o "x$CFITSIOLIB" = x; then
      AC_MSG_WARN([CFITSIO NOT AVAILABLE!!!!!! PLEASE READ THE INSTALLATION GUIDE. DO NOT ATTEMPT A BUILD.])
    else
      AC_MSG_NOTICE([CFITSIO appears to be available.])
      AC_DEFINE([HAVE_CFITSIO], [1], [Define to 1 if CFITSIO is available.])
   fi
fi

AC_SUBST([CFITSIOINC])
AC_SUBST([CFITSIOLIB])

AC_MSG_NOTICE([Configuring files...])

AC_CONFIG_FILES([Makefile
                 src/libraries/Makefile
				 scripts/Makefile
				 scripts/espadons/Makefile
				 scripts/guest/Makefile
                 src/core-espadons/Makefile
                 src/core-guest/Makefile
                 src/analysis-espadons/Makefile
                 src/analysis-guest/Makefile
                 src/tools/Makefile])

AC_OUTPUT




