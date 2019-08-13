# Copyright (C) 2010 Lou Hafer
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Lou Hafer        SFU        2011-06-30

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                               libOsi2                                 #
########################################################################

# Name of the library compiled in this directory.

lib_LTLIBRARIES = libOsi2.la

# List all source files for this library, including headers.

libOsi2_la_SOURCES = \
	Osi2API.hpp \
	Osi2ControlAPI.hpp \
	Osi2ControlAPI_Imp.hpp Osi2ControlAPI_Imp.cpp \
	Osi2CtrlAPIMessages.hpp Osi2CtrlAPIMessages.cpp \
	Osi2ParamMgmtAPI.hpp \
	Osi2ParamMgmtAPI_Imp.hpp Osi2ParamMgmtAPI_Imp.cpp \
	Osi2ParamMgmtAPIMessages.hpp Osi2ParamMgmtAPIMessages.cpp \
	Osi2ParamBEAPI_Imp.hpp \
	Osi2ParamBEAPIMessages.hpp Osi2ParamBEAPIMessages.cpp \
	Osi2RunParamsAPI.hpp \
	Osi2RunParamsAPI_Imp.hpp Osi2RunParamsAPI_Imp.cpp

# List all additionally required libraries.

libOsi2_la_LIBADD = $(OSI2LIB_LFLAGS)

# This is for libtool.

libOsi2_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags.

AM_CPPFLAGS = -DOSI2_BUILD \
	      -I$(srcdir)/../Osi2Plugin \
	      $(OSI2LIB_CFLAGS)
	      
########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir.

includecoindir = $(pkgincludedir)
includecoin_HEADERS = \
	Osi2API.hpp \
	Osi2ControlAPI.hpp \
	Osi2ProbMgmtAPI.hpp

