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

#ifdef WIN32
#include <Windows.h>
#else
#include <dlfcn.h>
#endif

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
{
    dfltPluginDir_ = std::string(OSI2DFLTPLUGINDIR) ;
}

DynamicLibrary::~DynamicLibrary ()
{
    if (handle_) {
#   ifdef WIN32
        ::FreeLibrary((HMODULE)handle_) ;
#   else
        if (::dlclose(handle_)) {
            std::string errorString ;
            errorString += "Failed to unload library \"" + fullPath_ + '"' ;
            const char *zErrorString = ::dlerror() ;
            if (zErrorString)
                errorString = errorString + ": " + zErrorString ;
            std::cerr << errorString << std::endl ;
        }
#   endif
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

# ifdef WIN32
    handle = ::LoadLibraryA(name.c_str()) ;
    if (handle == nullptr) {
        DWORD errorCode = ::GetLastError() ;
        std::stringstream ss ;
        ss << "LoadLibrary(" << name << ") Failed. errorCode: " << errorCode ;
        errorString = ss.str() ;
    }
# else
    handle = ::dlopen(name.c_str(),RTLD_NOW) ;
    if (handle == nullptr) {
        errorString += "Failed to load library \"" + name + '"' ;
        const char *zErrorString = ::dlerror() ;
        if (zErrorString)
            errorString = errorString + ": " + zErrorString ;
        return (nullptr) ;
    }
# endif

    DynamicLibrary *dynLib = new DynamicLibrary(handle) ;
    dynLib->fullPath_ = name ;
    return (dynLib) ;
}

void *DynamicLibrary::getSymbol (const std::string &symbol,
                                 std::string &errorString)
{
    if (handle_ == nullptr) return (nullptr) ;

#ifdef WIN32
    return (::GetProcAddress((HMODULE)handle_, symbol.c_str())) ;
#else
    void *sym = ::dlsym(handle_,symbol.c_str()) ;
    if (sym == nullptr) {
        errorString += "Failed to load symbol \"" + symbol + '"' ;
        const char *zErrorString = ::dlerror() ;
        if (zErrorString)
            errorString = errorString + ": " + zErrorString ;
    }
    return (sym) ;
#endif
}

}  // end namespace Osi2
