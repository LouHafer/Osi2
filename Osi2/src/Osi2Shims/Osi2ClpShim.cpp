
/*
  Method definitions for ClpShim.

  $Id$
*/

#include <iostream>

#include "Osi2ClpShim.hpp"
#include "Clp_C_Interface.h"
#include "ClpConfig.h"

#include "Osi2DynamicLibrary.hpp"

#include "ClpSimplex.hpp"
#include "Osi2ProbMgmtAPI_Clp.hpp"

using namespace Osi2 ;

/*
  Default constructor.
*/
ClpShim::ClpShim ()
  : services_(0),
    libClp_(0),
    verbosity_(1)
{ }

typedef Clp_Simplex *(*ClpFactory)() ;
typedef ClpSimplex *(*ClpSimplexFactory)(Clp_Simplex *clp) ;

/*! \brief Object factory

  Create clp-specific objects to satisfy the Osi2 API specified as the
  \p objectType member of of \p params.

  Object "WildProbMgmt" is used for testing wildcard object creation.
*/
void *ClpShim::create (const ObjectParams *params)
{
  std::string what = reinterpret_cast<const char *>(params->apiStr_) ;
  void *retval = 0 ;

  std::cout << "Clp create: type " << what << "." << std::endl ;

  if (what == "ClpSimplex" || what == "ProbMgmt" || what == "WildProbMgmt") {
    std::cout
      << "Request to create " << what << " recognised." << std::endl ;
    ClpShim *shim = static_cast<ClpShim*>(params->ctrlObj_) ;
    DynamicLibrary *libClp = shim->libClp_ ;
    std::string errStr ;
    ClpFactory factory =
      reinterpret_cast<ClpFactory>(libClp->getSymbol("Clp_newModel",errStr)) ;
    if (factory == 0) {
      std::cout << "Apparent failure to find Clp_newModel." << std::endl ;
      std::cout << errStr << std::endl ;
      return (0) ;
    }
    Clp_Simplex *wrapper = factory() ;
    ClpSimplexFactory underlyingModel =
	reinterpret_cast<ClpSimplexFactory>
	    (libClp->getSymbol("Clp_model",errStr)) ;
    if (underlyingModel == 0) {
      std::cout << "Apparent failure to find Clp_model." << std::endl ;
      std::cout << errStr << std::endl ;
      return (0) ;
    }
    ClpSimplex *retval = underlyingModel(wrapper) ;
    if (what == "ProbMgmt" || what == "WildProbMgmt") {
      ProbMgmtAPI *probMgmt = new ProbMgmtAPI_Clp(libClp,wrapper) ;
      return (probMgmt) ;
    } else {
      return (retval) ;
    }
  } else {
      std::cout
	<< "Clp create: unrecognised type " << what << "." << std::endl ;
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
  Plugin initialisation method. The job here is to construct a parameter
  structure which can be passed to the plugin manager's registration method
  (supplied in the services parameter). The return value is the exit method to
  be called before the plugin is unloaded.
*/
extern "C"
ExitFunc initPlugin (PlatformServices *services)
{
  std::string version = CLP_VERSION ;
  std::cout
    << "Executing ClpShim::initPlugin, clp version "
    << version << "." << std::endl ;
/*
  Attempt to load clp.
*/
  std::string libClpName = "libClp.so.0" ;
  const char *tmp = reinterpret_cast<const char*>(services->dfltPluginDir_) ;
  std::string libPath(tmp) ;
  std::string errMsg ;
  std::string fullPath = libPath+"/"+libClpName ;
  DynamicLibrary *libClp = DynamicLibrary::load(fullPath,errMsg) ;
  if (libClp == 0) {
    std::cout
      << "Apparent failure opening " << fullPath << "." << std::endl ;
    std::cout
      << "Error is " << errMsg << "." << std::endl ;
    return (0) ;
  }
/*
  Create the plugin library state object, ClpShim.  Arrange to remember the
  handle to libClp, and our unique ID from the plugin manager.  Then stash a
  pointer to the shim in PlatformServices to return it to the plugin
  manager.  This allows the plugin manager to hand back the shim object with
  each call, which in turn allows us to remember what we're doing.
*/
  ClpShim *shim = new ClpShim() ;
  shim->setLibClp(libClp) ;
  shim->setPluginID(services->pluginID_) ;
  services->ctrlObj_ = static_cast<PluginState *>(shim) ;
/*
  RegisterParams.
*/
/*
  Fill in the rest of the registration parameters and invoke the registration
  method. We could specify a separate state object for each API, but so far
  that doesn't seem necessary --- just use the library's state object.
*/
  RegisterParams reginfo ;
  reginfo.ctrlObj_ = static_cast<PluginState *>(shim) ;
  reginfo.version_.major_ = 1 ;
  reginfo.version_.minor_ = CLP_VERSION_MINOR ;
  reginfo.lang_ = Plugin_CPP ;
  reginfo.pluginID_ = shim->getPluginID() ;
  reginfo.createFunc_ = ClpShim::create ;
  reginfo.destroyFunc_ = ClpShim::destroy ;
  int retval =
      services->registerObject_(
	  reinterpret_cast<const unsigned char*>("ClpSimplex"),&reginfo) ;
  if (retval < 0) {
    std::cout
      << "Apparent failure to register ClpSimplex plugin." << std::endl ;
    return (0) ;
  }
  retval =
      services->registerObject_(
	  reinterpret_cast<const unsigned char*>("ProbMgmt"),&reginfo) ;
  if (retval < 0) {
    std::cout
      << "Apparent failure to register ProgMgmt plugin." << std::endl ;
    return (0) ;
  }
/*
  Register a wildcard object (for testing)
*/
  retval =
      services->registerObject_(
	  reinterpret_cast<const unsigned char*>("*"),&reginfo) ;
  if (retval < 0) {
    std::cout
      << "Apparent failure to register wildcard plugin." << std::endl ;
    return (0) ;
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

