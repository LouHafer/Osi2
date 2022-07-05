/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2GlpkHeavyShim.cpp
    \brief Method definitions for GlpkHeavyShim.

  This shim is written to work directly with glpk objects (in particular,
  OsiGlpkSolverInterface). As such, it must be linked with the appropriate
  glpk libraries.
*/

#include <iostream>

#include "glpk.h"

#include "Osi2GlpkHeavyShim.hpp"

#include "Osi2Config.h"
#include "Osi2nullptr.hpp"
#include "Osi2DynamicLibrary.hpp"

#include "Osi2Osi1API_GlpkHeavy.hpp"

using namespace Osi2 ;

/*
  Default constructor.
*/
GlpkHeavyShim::GlpkHeavyShim ()
    : services_(nullptr),
      verbosity_(1)
{ }

/*! \brief Object factory

  Create glpk-specific objects to satisfy the Osi2 API specified as the
  \p objectType member of of \p params.
*/
void *GlpkHeavyShim::create (const ObjectParams *params)
{
    std::string what = reinterpret_cast<const char *>(params->apiStr_) ;
    void *retval = nullptr ;

    std::cout << "GlpkHeavy create: type " << what << "." << std::endl ;

    if (what == "Osi1") {
      std::cout
	  << "Request to create " << what << " recognised." << std::endl ;
      retval = new Osi1API_GlpkHeavy() ;
    } else {
      std::cout
	  << "GlpkHeavy create: unrecognised type " << what << "."
	  << std::endl ;
    }

    return (retval) ;
}

/*! \brief Object destructor

  Given that GlpkHeavyShim only hands out C++ objects that are derived from
  Osi2::API, we can simply invoke the destructor with delete.
*/
int32_t GlpkHeavyShim::destroy (void *victim, const ObjectParams *objParms)
{
    std::string what = reinterpret_cast<const char *>(objParms->apiStr_) ;
    std::cout
            << "Request to destroy " << what << " recognised." << std::endl ;
    API *api = static_cast<API *>(victim) ;
    delete api ;

    return (0) ;
}


/*
  Plugin initialisation method. The job here is to construct a parameter
  structure which can be passed to the plugin manager's registration method
  (supplied in the services parameter). The return value is the exit method to
  be called before the plugin is unloaded.
*/
extern "C" OSI2LIB_EXPORT
ExitFunc initPlugin (PlatformServices *services)
{
  std::string version = glp_version() ;
  std::cout
      << "Executing GlpkHeavyShim::initPlugin, glpk version "
      << version << "." << std::endl ;
/*
  Create the plugin library state object, GlpkHeavyShim.  Arrange to remember
  our unique ID from the plugin manager.  Then stash a pointer to the shim
  in PlatformServices to return it to the plugin manager.  This allows the
  plugin manager to hand back the shim object with each call, which in turn
  allows us to remember what we're doing.
*/
  GlpkHeavyShim *shim = new GlpkHeavyShim() ;
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
  reginfo.version_.minor_ = GLP_MINOR_VERSION ;
  reginfo.pluginID_ = shim->getPluginID() ;
  reginfo.lang_ = Plugin_CPP ;
  reginfo.ctrlObj_ = static_cast<PluginState *>(shim) ;
  reginfo.createFunc_ = GlpkHeavyShim::create ;
  reginfo.destroyFunc_ = GlpkHeavyShim::destroy ;
  int retval = services->registerAPI_(
	  reinterpret_cast<const unsigned char*>("Osi1"), &reginfo) ;
  if (retval < 0) {
    std::cout
	<< "Apparent failure to register Osi1 plugin." << std::endl ;
    return (nullptr) ;
  }

  return (cleanupPlugin) ;
}

/*
  Plugin cleanup method. Does whatever is needed to clean up after the plugin
  prior to unloading the library.
*/
extern "C" int32_t cleanupPlugin (const PlatformServices *services)
{
    return (0) ;
}

