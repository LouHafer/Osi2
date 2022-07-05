# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: configure.ac 1752 2011-06-19 14:49:45Z stefan $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Osi2],[devel],[osi2@list.coin-or.org],[coin-or-osi2],
	[https://projects.coin-or.org/Osi2])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Osi2 which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there.

AC_CONFIG_SRCDIR(src/Osi2Plugin/Osi2Plugin.hpp)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.

AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.

AC_LANG_PUSH(C++)

# Initialize libtool.

AC_COIN_PROG_LIBTOOL

# Assume we have native support for nullptr. We can remove this from the code
# later, if this turns out to be a valid assumption  --lh, 190529 --

AC_DEFINE(COIN_HAS_NULLPTR,[1],
          [Define to 1 if the C++ compiler supports nullptr])

#############################################################################
#                              COIN-OR components                           #
#############################################################################

# First the plugin handling library.

AC_COIN_CHK_PKG(CoinUtils,[Osi2Plugin])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_FINALIZE_FLAGS(Osi2Plugin)

# Next, do up flags for the main Osi2 code base, then tack on info from
# Osi2Plugin.

AC_COIN_FINALIZE_FLAGS(Osi2Lib)
AC_COIN_CHK_HERE(Osi2Plugin,Osi2Lib)

# All we need from Data projects are the COIN_HAS_XXX and XXX_DATA variables.

AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],,dataonly)

# Can we find CxxTest? We can use it in the unit test. If there's a ThirdParty
# project, use it. Otherwise, look for a system installation.

AC_COIN_CHK_PKG(CxxTest,[CxxUnitTestDeps],[coincxxtest])
if test $coin_has_cxxtest = yes ; then
  AC_SUBST(CXXTESTGEN)
  if test -x "${prefix}/bin/cxxtestgen" ; then
    CXXTESTGEN="\$(prefix)/bin/cxxtestgen"
  else
    AC_MSG_WARN([cxxtestgen not in ${prefix}/bin as expected.])
    AC_MSG_WARN([Check that ThirdParty/CxxTest is properly installed?])
    coin_has_cxxtest=no
    COIN_HAS_CXXTEST_TRUE='#'
    COIN_HAS_CXXTEST_FALSE=''
  fi
else
  AC_PATH_PROG(CXXTESTGEN,cxxtestgen)
  if test -n "$CXXTESTGEN" ; then
    coin_has_cxxtest=yes
    COIN_HAS_CXXTEST_TRUE=''
    COIN_HAS_CXXTEST_FALSE='#'
  fi
fi
if test x${coin_has_cxxtest} = xyes ; then
  AC_MSG_NOTICE([Using $CXXTESTGEN])
else
  AC_MSG_NOTICE([Not using CxxTest.])
fi

#############################################################################
#                            Shims                                          #
#############################################################################

# One can build ClpShim without libClp, but ClpHeavyShim must be linked with
# libOsiClp and dependencies. But see below for separate additions to CFLAGS.

AC_COIN_CHK_PKG(OsiClp,[Osi2ClpHeavyShim],[osi-clp])
if test $coin_has_osiclp != yes ; then
  AC_MSG_ERROR([Required package OsiClp not available.])
fi

# Similarly, GlpkHeavyShim needs libOsiGlpk and dependencies.

AC_COIN_CHK_PKG(OsiGlpk,[Osi2GlpkHeavyShim],[osi-glpk])

# Finalize flags, then tack on the local stuff

AC_COIN_FINALIZE_FLAGS([Osi2ClpShim Osi2ClpHeavyShim Osi2GlpkHeavyShim])
AC_COIN_CHK_HERE(Osi2Lib,[Osi2ClpShim Osi2ClpHeavyShim Osi2GlpkHeavyShim],osi2)

# Create some artificial variables that'll be useful in the makefile when
# building unit tests. Finalize, then tack on Osi2.

AC_COIN_CHK_PKG(Osi,[UnitTestDeps CxxUnitTestDeps])
AC_COIN_FINALIZE_FLAGS([UnitTestDeps CxxUnitTestDeps])
AC_COIN_CHK_HERE(Osi2Lib,[UnitTestDeps CxxUnitTestDeps],osi2)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 test/Makefile
                 cxxtest/Makefile
		 src/Osi2Plugin/Makefile
		 src/Osi2/Makefile
		 src/Osi2Shims/Makefile
                 osi2.pc
                 osi2plugin.pc
                 osi2clpshim.pc
                 osi2clpheavyshim.pc
                 osi2-unittests.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/Osi2/config.h src/Osi2/config_osi2.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
