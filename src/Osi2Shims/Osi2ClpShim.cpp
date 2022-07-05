/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ClpShim.cpp
    \brief Method definitions for ClpShim.

  This shim is written to dynamically load libClp. As such, it doesn't
  need to be statically linked with libClp. The tradeoff is that it must
  work through clp's C interface. The shim, and any objects it constructs,
  must dynamically load any methods they wish to use.
*/

#include <iostream>
#include <vector>

#include "ClpConfig.h"
#include "Osi2ClpShim.hpp"
#include "Clp_C_Interface.h"

#include "Osi2Config.h"
#include "Osi2nullptr.hpp"

#include "Osi2ClpLite_Wrap.hpp"

typedef Clp_Simplex *(*ClpFactory)() ;

using namespace Osi2 ;

/*
  Default constructor.
*/
ClpShim::ClpShim ()
    : services_(0),
      libClp_(0),
      ourID_(0),
      verbosity_(1)
{ }

/*! \brief Object factory

  Create clp-specific objects to satisfy the OSI2 API specified as the
  \p objectType member of of \p params.

  The name "WildClpSimplex" is used for testing wildcard object creation.
*/
void *ClpShim::create (const ObjectParams *params)
{

  std::string what = reinterpret_cast<const char *>(params->apiStr_) ;
  void *retval = nullptr ;

  std::cout << "Clp create: " << what << " API requested." << std::endl ;
/*
  For an object that supports the ClpSimplex API, create a ClpLite_Wrap object
  and add a ClpSimplexAPI_ClpLite object to it.
*/
  if (what == "ClpSimplex" || what == "WildClpSimplex") {
    std::cout
      << "  " << what << " API is  supported." << std::endl ;
    ClpShim *shim = static_cast<ClpShim*>(params->ctrlObj_) ;
    DynamicLibrary *libClp = shim->libClp_ ;
    ClpLite_Wrap *clpliteWrapper = new ClpLite_Wrap() ;
    clpliteWrapper->addClpSimplex(libClp) ;
    retval = clpliteWrapper ;
  } else {
    std::cout
      << "    " << what << " API is not supported." << std::endl ;
  }

  return (retval) ;
}

/*! \brief Object destructor

  Given that ClpShim only hands out C++ objects that are derived from
  Osi2::API, we can simply invoke the destructor with delete.
*/
int32_t ClpShim::destroy (void *victim, const ObjectParams *objParms)
{
    std::string what = reinterpret_cast<const char *>(objParms->apiStr_) ;
    std::cout
            << "Request to destroy " << what << " recognised." << std::endl ;
    API *api = static_cast<API *>(victim) ;
    delete api ;

    return (0) ;
}


/*
  Plugin cleanup method. Does whatever is needed to clean up after the plugin
  prior to unloading the library.
*/
extern "C"
int32_t cleanupPlugin (const PlatformServices *services)
{
  std::cout << "Executing ClpShim cleanupPlugin." << std::endl ;
    return (0) ;
}

/*
  Plugin initialisation method. The first task is to find and load libClp.
  Once that's done, register the APIs we provide.  The return value is the
  exit method to be called before the plugin is unloaded.
*/
extern "C" OSI2LIB_EXPORT
ExitFunc initPlugin (PlatformServices *services)
{
  std::string version = CLP_VERSION ;
  std::cout
      << "Executing ClpShim::initPlugin, clp version "
      << version << "." << std::endl ;
/*
  Attempt to load clp. The plugin search path is composed of directories
  separated by ':' characters. Stop with the first successful load.
*/
  std::string libClpName = "libClp.so.0" ;
  const char *tmp = reinterpret_cast<const char*>(services->plugSrchPath_) ;
  std::string searchDirs(tmp) ;
  std::string libPath ;
  std::string::size_type sepPos = 0 ;
  std::string::size_type startPos = 0 ;
  std::vector<std::pair<std::string,std::string>> errMsgs ;
  DynamicLibrary *libClp = nullptr ;
  while (sepPos < std::string::npos && libClp == nullptr) {
    sepPos = searchDirs.find_first_of(':',startPos) ;
    libPath = searchDirs.substr(startPos,sepPos-startPos) ;

    std::string errMsg ;
    std::string fullPath = libPath + "/" + libClpName ;
    libClp = DynamicLibrary::load(fullPath, errMsg) ;
    if (libClp != nullptr) break ;
    errMsgs.push_back(std::pair<std::string,std::string>(fullPath,errMsg)) ;
    startPos = sepPos+1 ;
  }
  if (libClp == nullptr) {
    std::cout
	<< "Apparent failure opening " << libClpName << "." << std::endl ;
    std::cout << "Error(s):" << std::endl ;
    std::vector<std::pair<std::string,std::string>>::const_iterator iter ;
    for (iter = errMsgs.begin() ; iter != errMsgs.end() ; iter++) {
      std::pair<std::string,std::string> errPair = *iter ;
      std::cout
        << "  Path '" << errPair.first << "', error '" << errPair.second
	<< "." << std::endl ;
    }
    return (nullptr) ;
  }
/*
  Create the plugin library state object, ClpShim.  Arrange to remember the
  handle to libClp, and our unique ID from the plugin manager.  Then stash a
  pointer to the shim in PlatformServices to return it to the plugin manager.
  This allows the plugin manager to hand back the shim object with each call,
  which in turn allows us to remember what we're doing.
*/
  ClpShim *shim = new ClpShim() ;
  shim->setLibClp(libClp) ;
  shim->setPluginID(services->pluginID_) ;
  services->ctrlObj_ = static_cast<PluginState *>(shim) ;
/*
  Register our APIs.

  For each API, fill in the the registration parameters and invoke the
  registration function. We could specify a separate state object for each
  API, but so far that doesn't seem necessary --- just use the library's
  state object.
*/
  APIRegInfo reginfo ;
  reginfo.version_.major_ = 1 ;
  reginfo.version_.minor_ = CLP_VERSION_MINOR ;
  reginfo.pluginID_ = shim->getPluginID() ;
  reginfo.lang_ = Plugin_CPP ;
  reginfo.ctrlObj_ = static_cast<APIState *>(shim) ;
  reginfo.createFunc_ = ClpShim::create ;
  reginfo.destroyFunc_ = ClpShim::destroy ;
  int retval = services->registerAPI_(
	  reinterpret_cast<const CharString*>("ClpSimplex"),&reginfo) ;
  if (retval < 0) {
    std::cout
	<< "Apparent failure to register ClpSimplex API." << std::endl ;
    return (nullptr) ;
  }
/*
  Register a wildcard object (just for testing)
*/
  retval =
      services->registerAPI_(
	  reinterpret_cast<const CharString*>("*"), &reginfo) ;
  if (retval < 0) {
    std::cout
	<< "Apparent failure to register wildcard plugin." << std::endl ;
    return (nullptr) ;
  }
/*
  std::cout
    << "   in ClpShim initPlugin, cleanupPlugin is "
    << std::hex << (void *) cleanupPlugin << std::dec
    << std::endl ;
*/
  return (cleanupPlugin) ;
}
