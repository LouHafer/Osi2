/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  Include file for the configuration of Osi2.
  
  On systems where the code is configured with the configure script
  (i.e., compilation is always done with HAVE_CONFIG_H defined), this
  header file includes the automatically generated header file, and
  undefines macros that might configure with other config.h files.
  
  On systems that are compiled in other ways (e.g., with MS Visual
  Studio), a header file is included to define those macros that depend
  on the operating system and the compiler.  The macros that define
  the configuration of the particular build environment (e.g., presence
  of other COIN-OR packages or third party code) are set in the files
  config_*default.h. As a user, you can hope that these will be correct for
  your environment, but you really should check.
*/
/*! \file Osi2Config.h

  Top level for Osi2 configuration headers. Include this file to access Osi2
  configuration information and allow it to find the right configuration
  header files.
*/


#ifndef Osi2Config_H
# define Osi2Config_H

#ifdef HAVE_CONFIG_H

#  ifdef OSI2_BUILD
#    include "config.h"
#    ifdef OSI2LIB_EXPORT
#      undef OSI2LIB_EXPORT
#    endif
#    ifdef DLL_EXPORT
#      define OSI2LIB_EXPORT __declspec(dllexport)
#    elif defined(__GNUC__) && __GNUC__ >= 4
#      define OSI2LIB_EXPORT __attribute__((__visibility__("default")))
#    else
#      define OSI2LIB_EXPORT
#    endif
#  else /* OSI2_BUILD */
#    include "config_osi2.h"
#  endif

#else /* HAVE_CONFIG_H */

#  ifndef OSI2LIB_EXPORT
#    if defined(_WIN32) && defined(DLL_EXPORT)
#      ifdef OSI2_BUILD
#        define OSI2LIB_EXPORT __declspec(dllexport)
#      else
#        define OSI2LIB_EXPORT __declspec(dllimport)
#      endif
#    elif defined(__GNUC__) && __GNUC__ >= 4
#      ifdef OSI2_BUILD
#        define OSI2LIB_EXPORT __attribute__((__visibility__("default")))
#      else
#        define OSI2LIB_EXPORT
#      endif
#    else
#      define OSI2LIB_EXPORT
#    endif
#  endif
#  ifdef OSI2LIB_BUILD
#    include "config_default.h"
#  else
#    include "config_osi2_default.h"
#  endif

#endif /* HAVE_CONFIG_H */

#endif  /* Osi2Config_H */
