/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  Based on original design and code by Gigi Sayfan published in five parts
  in Dr. Dobbs, starting November 2007.
*/
/*! \file Osi2DynamicLibrary.cpp

  This file contains the methods that handle loading and unloading of
  libraries.

  Currently works only for Linux & Solaris. There are hooks for Windows but
  they are untested and/or disabled. The eventual plan is to make this work
  with libtldl.  -- lh, 110825 --
*/

#include "Osi2Config.h"
#include "Osi2nullptr.hpp"
#include "Osi2DynamicLibrary.hpp"
#include <sstream>
#include <iostream>

/*
  Default for OSI2PLUGINDIR is $(libdir) (library installation directory).
  Specified with -D in the Makefile.
*/
#ifndef OSI2PLUGINDIR
# define OSI2DFLTPLUGINDIR "/usr/local/lib"
#else
# define OSI2DFLTPLUGINDIR OSI2PLUGINDIR
#endif

namespace Osi2 {

DynamicLibrary::DynamicLibrary (void *handle)
    : handle_(handle)
{  }

DynamicLibrary::~DynamicLibrary ()
{
  if (handle_) {
      if (::dlclose(handle_)) {
	  std::string errorString ;
	  errorString += "Failed to unload library \"" + fullPath_ + '"' ;
	  const char *zErrorString = ::dlerror() ;
	  if (zErrorString)
	      errorString = errorString + ": " + zErrorString ;
	  std::cerr << errorString << std::endl ;
      }
  }
}

DynamicLibrary *DynamicLibrary::load (const std::string &name,
                                      std::string &errorString)
{
  if (name.empty()) {
      errorString = "Empty path." ;
      return (nullptr) ;
  }

  void *handle = nullptr ;

  handle = ::dlopen(name.c_str(),RTLD_NOW) ;
  if (handle == nullptr) {
      if (!errorString.empty()) {
        errorString += "\n" ;
      }
      errorString += "Failed to load library \"" + name + '"' ;
      const char *zErrorString = ::dlerror() ;
      if (zErrorString)
	  errorString = errorString + ": " + zErrorString ;
      return (nullptr) ;
  }

  DynamicLibrary *dynLib = new DynamicLibrary(handle) ;
  dynLib->fullPath_ = name ;
  return (dynLib) ;
}

std::ostream &operator<< (std::ostream &ostr, const DynamicLibrary *lib)
{
  ostr
    << "  DynamicLibrary (" << (static_cast<const void *>(lib)) << "):"
    << std::endl ;
  ostr << "    handle: " << lib->handle_ << std::endl ;
  ostr << "    path: " << lib->fullPath_ ;
  return (ostr) ;
}

}  // end namespace Osi2
