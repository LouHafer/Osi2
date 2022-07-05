/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  Based on and extended from original design and code by Gigi Sayfan published
  in five parts in Dr. Dobbs, starting November 2007.
*/
/*! \file Osi2PluginManager.cpp
    \brief Definition of methods for Osi2::PluginManager.
*/

#include <cstring>
#include <cassert>
#include <string>
#include <iostream>

#include "CoinHelperFunctions.hpp"

#include "Osi2PluginManager.hpp"
#include "Osi2DynamicLibrary.hpp"
#include "Osi2ObjectAdapter.hpp"

#ifndef OSI2PLUGINDIR
# define OSI2DFLTPLUGINDIR "/usr/local/lib"
#else
# define OSI2DFLTPLUGINDIR OSI2PLUGINDIR
#endif

using namespace Osi2 ;

namespace {
/*
  Utilities for the PlatformServices struct. This can't be an honest class
  if we want to support C language plugins because a pointer to the struct is
  passed to the plugin for its use.
*/

/*
  A utility method to convert the std::vector<std::string> form for plugin
  search directories to a persistent char string allocated on the heap. This
  can be safely passed to plugin libraries.
*/

char *safePlugSrchPath (const std::vector<std::string> &plugSrchDirs)
{ std::string plugSrchPath ;
  plugSrchPath.clear() ;
  std::vector<std::string>::const_iterator iter ;
  for (iter = plugSrchDirs.cbegin() ; iter != plugSrchDirs.cend() ; iter++) {
    if (!plugSrchPath.empty()) plugSrchPath += ':' ;
    plugSrchPath += *iter ;
  }
  std::string::size_type pathLen = plugSrchPath.length()+1 ;
  char *tmp = new char[pathLen] ;
  std::strncpy(tmp,plugSrchPath.c_str(),pathLen) ;
  return (tmp) ;
}

/*
  PlatformServices `constructor'

  Initialises the PlatformServices struct. The one thing that takes work is
  creating a null-terminated character string that's persistent, hence safe to
  pass a pointer to a plugin library.
*/
void initPlatformServices (PlatformServices &platServ,
			   const std::vector<std::string> &plugSrchDirs,
			   APIRegFunc regFunc)
{ platServ.version_ = { 1, 0 } ;
  platServ.plugSrchPath_ =
    reinterpret_cast<const CharString*>(safePlugSrchPath(plugSrchDirs)) ;
  platServ.pluginID_ = 0 ;
  platServ.ctrlObj_ = nullptr ;
  platServ.registerAPI_ = regFunc ;
  platServ.invokeService_ = nullptr ;
}

/*
  PlatformServices `destructor'

  We own the string in plugSrchPath_, but that's all. In particular, we're not
  responsible for the control object ctrlObj_.
*/
void clearPlatformServices (PlatformServices &platServ)
{ platServ.version_ = { 0, 0 } ;
  delete[] platServ.plugSrchPath_ ;
  platServ.plugSrchPath_ = nullptr ;
  platServ.pluginID_ = 0 ;
  platServ.ctrlObj_ = nullptr ;
  platServ.registerAPI_ = nullptr ;
  platServ.invokeService_ = nullptr ;
}

/*
  PlatformServices update plugSrchPath_
*/
void updatePlatSrvSrchPath (PlatformServices &platServ,
			    const std::vector<std::string> &plugSrchDirs)
{ delete[] platServ.plugSrchPath_ ;
  platServ.plugSrchPath_ =
    reinterpret_cast<const CharString*>(safePlugSrchPath(plugSrchDirs)) ;
}

}  // end anonymous namespace

/*
  Given an initialisation function, do the necessary bookkeeping to register
  the library and its APIs with the plugin manager.

  The method works for innate (compiled in) and dynamic libraries. The sole
  purpose of the dynLib parameter is to correctly initialise the LibraryInfo
  structure; it should be null for an innate library.

  The expected flow of events is that the library's initFunc will
  register the APIs provided by the library via callbacks using
  PluginManager::registerAPI. Once this is complete, the initFunc will
  return and the plugin manager bookkeeping can be completed.

  The method returns the library's unique ID, or 0 for failure.
*/

PluginUniqueID PluginManager::initOneLib (std::string fullPath,
	InitFunc initFunc, DynamicLibrary *dynLib)
{

/*
  Invoke the initialisation function.

  Set initialisingPlugin_ to true so that the APIs will be registered into
  the temporary wildcard and exact match vectors. If initialisation is
  successful, we'll transfer them to the permanent vectors.
*/
  initialisingPlugin_ = true ;
  libInInit_ = genUniqueID() ;
  pathInInit_ = fullPath ;
  tmpExactMatchMap_.clear() ;
  tmpWildCardVec_.clear() ;
  platformServices_.pluginID_ = libInInit_ ;
  ExitFunc exitFunc = initFunc(&platformServices_) ;
  if (exitFunc == nullptr) {
    msgHandler_->message(PLUGMGR_LIBINITFAIL, msgs_)
	<< fullPath << CoinMessageEol ;
    initialisingPlugin_ = false ;
    return (0) ;
  }
  msgHandler_->message(PLUGMGR_LIBINITOK, msgs_)
      << ((dynLib == nullptr)?"innate":"plugin")
      << fullPath << CoinMessageEol ;
/*
  We have happiness: the library is initialised. Do the bookkeeping.
  Enter the library in the library map. Add the exit function to the vector
  of exit functions, and copy information from the temporary wildcard and
  exact match vectors to the permanent vectors.

  The library may well be dynamic, but it's the responsibility of the calling
  method to deal with that question.

  \todo
  How paranoid do we want to be? Given that we check entries for uniqueness
  as they're going into tmpExactMatchMap_, and that you can't load the same
  lib twice, should we check tmpExactMatchMap_ against exactMatchMap_?
  -- lh, 111013 --
*/
  libPathToIDMap_[fullPath] = libInInit_ ;
  LibraryInfo &info = libraryMap_[libInInit_] ;
  info.id_ = libInInit_ ;
  info.isDynamic_ = ((dynLib == nullptr)?false:true) ;
  info.dynLib_ = dynLib ;
  info.ctrlObj_ = platformServices_.ctrlObj_ ;
  info.exitFunc_ = exitFunc ;
  exactMatchMap_.insert(tmpExactMatchMap_.begin(), tmpExactMatchMap_.end()) ;
  tmpExactMatchMap_.clear() ;
  wildCardVec_.insert(wildCardVec_.end(),
		      tmpWildCardVec_.begin(), tmpWildCardVec_.end()) ;
  tmpWildCardVec_.clear() ;
  initialisingPlugin_ = false ;
  libInInit_ = 0 ;

  msgHandler_->message(PLUGMGR_LIBLDOK, msgs_)
      << ((info.isDynamic_ == true)?"plugin":"innate")
      << fullPath << CoinMessageEol ;

  return (info.id_) ;
}

/*
  Plugin manager constructor.

  There can only be one! Correct operation of the PluginManager depends on
  the fact that there is exactly one static instance, created at the first
  invocation of PluginManager::getInstance.
*/
PluginManager::PluginManager()
  : currentID_(0),
    initialisingPlugin_(false),
    dfltInnateDir_("innate"),
    dfltHandler_(true),
    logLvl_(7)
{
  msgHandler_ = new CoinMessageHandler() ;
  msgs_ = PlugMgrMessages() ;
  msgHandler_->setLogLevel(logLvl_) ;
  msgHandler_->message(PLUGMGR_INIT, msgs_) << CoinMessageEol ;
  plugSrchDirs_ = std::vector<std::string>() ;
  plugSrchDirs_.push_back(std::string(OSI2DFLTPLUGINDIR)) ;
  initPlatformServices(platformServices_,plugSrchDirs_,registerAPI) ;
}

/*
  Plugin manager destructor

  Call shutdown to close out the plugins, then destroy the pieces of ourself.
*/
PluginManager::~PluginManager ()
{
/*
  Just in case it wasn't called earlier.
*/
  shutdown() ;
/*
  Clean out the PlatformServices struct.
*/
  clearPlatformServices(platformServices_) ;
/*
  If this is our message handler, delete it. Otherwise it's the client's
  responsibility.
*/
  if (dfltHandler_) {
      delete msgHandler_ ;
      msgHandler_ = nullptr ;
  }
}

/*
  Set the plugin search directories using a vector of strings.
*/
void PluginManager::setPluginDirs(const std::vector<std::string> searchDirs)
{
  plugSrchDirs_ = searchDirs ;
  updatePlatSrvSrchPath(platformServices_,plugSrchDirs_) ;
}
/*
  Methods to set and get the plugin search path using a single string.
*/
void PluginManager::setPluginDirsStr (const std::string searchDirs)
{
  std::string::size_type sepPos = 0 ;
  std::string::size_type startPos = 0 ;

  plugSrchDirs_.clear() ;
  while (sepPos < std::string::npos) {
    sepPos = searchDirs.find_first_of(':',startPos) ;
    plugSrchDirs_.push_back(searchDirs.substr(startPos,sepPos-startPos)) ;
    startPos = sepPos+1 ;
  }
  updatePlatSrvSrchPath(platformServices_,plugSrchDirs_) ;
}

std::string PluginManager::getPluginDirsStr () const
{
  std::string dirs ;
  dirs.clear() ;
  for (std::vector<std::string>::const_iterator iter = plugSrchDirs_.cbegin() ;
       iter != plugSrchDirs_.cend() ;
       iter++)
  { if (!dirs.empty()) dirs += ':' ;
    dirs += *iter ; }

  return (dirs) ;
}
     

/*
  There's no particular reason to trust a plugin, so we need to do some
  validation of the parameters it gives to the plugin manager when registering
  an API.

  \parm Registration string supplied by plugin
  \parm Registration registration information block supplied by plugin

  \returns true if the parameter block validates, false otherwise
*/

bool PluginManager::validateAPIRegInfo (const CharString *apiStr,
		const APIRegInfo *params) const
{
    bool retval = true ;
    std::string errStr = "" ;
    retval = true ;

    // The API to be registered must have a name
    if (apiStr == nullptr || !(*apiStr)) {
        errStr += "; bad API ID string" ;
        retval = false ;
    }

    // If we're missing the parameter block, nothing more to say.
    if (params == nullptr) {
        errStr += "; no registration parameter block" ;
        retval = false ;
    } else {
        /*
          The pluginID must be valid. If we're initialising the library,
	  it will not appear in libraryMap. The library is allowed to register
	  new APIs at any time, so allow for that.
        */
        if (initialisingPlugin_) {
            if (libInInit_ != params->pluginID_) {
                errStr += "; bad plugin library ID" ;
                retval = false ;
            } 
        } else {
            LibraryMap::const_iterator lmIter ;
	    lmIter = libraryMap_.find(params->pluginID_) ;
            if (lmIter == libraryMap_.end()) {
                errStr += "; bad plugin library ID" ;
                retval = false ;
            }
        }
        // There must be a constructor
        if (params->createFunc_ == nullptr) {
            errStr += "; missing constructor" ;
            retval = false ;
        }
        // There must be a destructor
        if (params->destroyFunc_ == nullptr) {
            errStr += "; missing destructor" ;
            retval = false ;
        }
        // The major version must match
        PluginAPIVersion ver = platformServices_.version_ ;
        if (ver.major_ != params->version_.major_) {
            msgHandler_->message(PLUGMGR_BADVER, msgs_)
                    << params->version_.major_ << ver.major_ << CoinMessageEol ;
            errStr += "; version mismatch" ;
            retval = false ;
        }
    }

    if (!retval) {
        errStr = errStr.substr(2) ;
        msgHandler_->message(PLUGMGR_APIBADPARM, msgs_)
	    << errStr << CoinMessageEol ;
    }

    return (retval) ;
}


/*
  Check for an entry that matches key and libID. Returns an iterator for the
  entry, if it exists, else regMap.end(). A libID of 0 is a wildcard.
*/
PluginManager::APIRegMap::const_iterator
PluginManager::apiEntryExists (const APIRegMap &regMap,
                               const std::string &key,
                               PluginUniqueID libID)
{
    typedef APIRegMap::const_iterator RMCI ;
    /*
      See if we find anything. If we have a hit, iterPair.first->first will be
      the key we're looking for; otherwise, we've struck out.
    */
    std::pair<RMCI, RMCI> iterPair = regMap.equal_range(key) ;
    if (iterPair.first->first != key) return (regMap.end()) ;
    /*
      For the case where the libID is a wildCard (0), we're perfectly happy
      to use the first entry found.
    */
    if (libID == 0) return (iterPair.first) ;
    /*
      There's at least one entry matching the API (key). Check to see if any of
      them are from the same plugin library. If not, return failure.
    */
    for (RMCI iter = iterPair.first ; iter != iterPair.second ; iter++) {
        PluginUniqueID id = iter->second.id_ ;
        if (id == libID) return (iter) ;
    }

    return (regMap.end()) ;
}



/*
  This method implements a callback from the plugin to register the APIs
  it supports.  A pointer to this method is passed to the plugin in a
  PlatformServices object.

  Returns 0 for success, -1 for failure.
*/
int32_t PluginManager::registerAPI (const CharString *apiStr,
				    const APIRegInfo *params)
{
    PluginManager &pm = getInstance() ;

    // Validate the parameter block
    bool val = pm.validateAPIRegInfo(apiStr, params) ;
    if (!val) return (-1) ;

    const std::string api(reinterpret_cast<const char *>(apiStr)) ;
    PluginUniqueID id = params->pluginID_ ;
    int retval = 0 ;
/*
  If this call comes while we're initialising a plugin, add the APIs
  to temporary vectors for merge into permanent vectors if and when
  initialisation is successful.

  Pointers to members. Always a treat.
*/
    APIRegVec PluginManager:: *wcVec = nullptr ;
    APIRegMap PluginManager:: *regMap = nullptr ;
    if (pm.initialisingPlugin_) {
      wcVec = &PluginManager::tmpWildCardVec_ ;
      regMap = &PluginManager::tmpExactMatchMap_ ;
    } else {
      wcVec = &PluginManager::wildCardVec_ ;
      regMap = &PluginManager::exactMatchMap_ ;
    }
/*
  If the registration is a wild card, add it to the wild card vector,
  otherwise add to the exact match vector that tracks specific APIs.

  Duplicates are not allowed --- each <api,lib> pair must be unique.
*/
    if (api == std::string("*")) {
      APIRegVec::const_iterator rvIter ;
      for (APIRegVec::const_iterator rvIter = (pm.*wcVec).begin() ;
           rvIter != (pm.*wcVec).end() ;
	   rvIter++) {
        if (rvIter->api_ == api) {
	  retval = -1 ;
	  break ;
	}
      }
      if (retval != -1) {
        APIInfo tmp ;
	tmp.api_ = "*" ;
	tmp.id_ = id ;
	tmp.ctrlObj_ = params->ctrlObj_ ;
	tmp.lang_ = params->lang_ ;
	tmp.createFunc_ = params->createFunc_ ;
	tmp.destroyFunc_ = params->destroyFunc_ ;
        (pm.*wcVec).push_back(tmp) ;
      }
    } else {
      APIRegMap::const_iterator rmIter ;
      rmIter = pm.apiEntryExists(pm.*regMap,api,id) ;
      if (rmIter != (pm.*regMap).end()) {
        retval = -1 ;
      } else {
        APIInfo tmp ;
	tmp.api_ = api ;
	tmp.id_ = id ;
	tmp.ctrlObj_ = params->ctrlObj_ ;
	tmp.lang_ = params->lang_ ;
	tmp.createFunc_ = params->createFunc_ ;
	tmp.destroyFunc_ = params->destroyFunc_ ;
	(pm.*regMap).insert(std::pair<std::string,APIInfo>(tmp.api_,tmp)) ;
      }
    }
    if (!retval) {
      pm.msgHandler_->message(PLUGMGR_APIREGOK, pm.msgs_)
	  << api << pm.getLibPath(id) << CoinMessageEol ;
    } else {
      pm.msgHandler_->message(PLUGMGR_APIREGDUP,pm.msgs_)
	  << api << CoinMessageEol ; }

    return (retval) ;
}

/*
  Instantiate a single instance of the manager.
*/
PluginManager &PluginManager::getInstance()
{
    static PluginManager instance ;

    return (instance) ;
}


/*
  Load a plugin.

  If uniqueID is supplied, it will be loaded with the PluginUniqueID assigned
  to the library.

  Returns:
    -3: library failed to initialise
    -2: failed to find initFunc
    -1: library failed to load
     0: library loaded and initialised without error
     1: library is already loaded

  The ability to `load' an innate plugin is limited. It's expected that innate
  plugins will preregister --- otherwise there's no way to know the address of
  the initFunc. The only time this method would get involved is if the user
  actually `unloaded' the innate plugin and now wants it back.

  \todo Implement platform-independent path handling.
  \todo Implement resolution of symbolic links.
  \todo Implement conversion to absolute path.
*/
int PluginManager::loadOneLib (const std::string &libName,
			       const std::string *dir,
                               PluginUniqueID *uniqueID)
{
  if (uniqueID != 0) (*uniqueID) = 0 ;
/*
  If no directory is specified, consider both the plugin search path and
  the innate `directory'. If a directory is specified, use that exclusively.
*/
  std::vector<std::string> plugDirs ;
  std::string pluginPath ;
  std::string innatePath ;
  std::string fullPath = "" ;
  char dirSep = CoinFindDirSeparator() ;

  if (dir == nullptr || (dir->compare("") == 0)) {
    plugDirs = getPluginDirs() ;
    innatePath = getInnateDir() ;
  } else {
    plugDirs.push_back(*dir) ;
    innatePath = "" ;
  }
/*
  Is this library already loaded? If so, don't do it again.
*/
  bool foundIt = false ;
  LibPathToIDMap::iterator thisLib ;
  for (std::vector<std::string>::const_iterator iter = plugDirs.cbegin() ;
       iter != plugDirs.cend() && !foundIt ;
       iter++) {
    pluginPath = *iter+dirSep+libName ;
    if (libPathToIDMap_.find(pluginPath) != libPathToIDMap_.end()) {
      foundIt = true ;
      thisLib = libPathToIDMap_.find(pluginPath) ;
      fullPath = pluginPath ;
    }
  }
  if (!foundIt && innatePath != "") {
    pluginPath = innatePath+dirSep+libName ;
    if (libPathToIDMap_.find(pluginPath) != libPathToIDMap_.end()) {
      foundIt = true ;
      thisLib = libPathToIDMap_.find(pluginPath) ;
      fullPath = pluginPath ;
    }
  }
  if (foundIt) {
    msgHandler_->message(PLUGMGR_LIBLDDUP, msgs_)
	<< fullPath << CoinMessageEol ;
    if (uniqueID != 0) { *uniqueID = thisLib->second ; }
    return (1) ; }
/*
  The library isn't loaded. Look for an honest plugin first (so that the user
  can override an innate plugin if they choose). If we don't find one, try for
  an innate plugin if we can.
*/
  DynamicLibrary *dynLib = nullptr ;
  InitFunc initFunc = nullptr ;
  std::string errStr ;

  for (std::vector<std::string>::const_iterator iter = plugDirs.cbegin() ;
       iter != plugDirs.cend() && dynLib == nullptr ;
       iter++) { 
    pluginPath = *iter+dirSep+libName ;
    dynLib = DynamicLibrary::load(pluginPath, errStr) ;
  }
  if (dynLib != nullptr) {
    fullPath = pluginPath ;
    initFunc = dynLib->getFunc<InitFunc>("initPlugin",errStr) ;
  } else if (innatePath != "") {
    pluginPath = innatePath+dirSep+libName ;
    PreloadMap::const_iterator plmIter = preloadLibs_.find(pluginPath) ;
    if (plmIter != preloadLibs_.end()) {
      fullPath = pluginPath ;
      initFunc = plmIter->second ;
    }
  }
/*
  Check for errors: no library or no init function.
*/
  if (fullPath == "") {
    msgHandler_->message(PLUGMGR_LIBLDFAIL, msgs_)
	<< libName << errStr << CoinMessageEol ;
    return (-1) ;
  }
  if (initFunc == nullptr) {
    msgHandler_->message(PLUGMGR_SYMLDFAIL, msgs_)
	<< "function" << "initPlugin" << fullPath << errStr
	<< CoinMessageEol ;
    return (-2) ;
  }
/*
  Invoke the initialisation function to complete the initialisation and do
  PluginManager bookkeeping.
*/
  PluginUniqueID id = initOneLib(fullPath,initFunc,dynLib) ;
  if (!id) {
    delete dynLib ;
    return (-3) ;
  }

  if (uniqueID != 0) (*uniqueID) = id ;

  return (0) ;
}


/*
  We're not going to implement loadAll functionality in the first round.
  Make this method into a guaranteed failure, so we'll catch it when it's
  called. Unclear we'll ever want this.
*/

int PluginManager::loadAllLibs (const std::string &libDir,
                                InvokeServiceFunc func)
{
    assert(false) ;

    return (-1) ;
}

#if defined(OSI2_IMPLEMENT_LOADALL)

int32_t PluginManager::loadAll(const std::string &libDir,
                               InvokeServiceFunc func)
{

#if defined(OSI2PLATFORM_MAC)
  std::string dynamicLibraryExtension("dylib") ;
#elif defined(OSI2PLATFORM_LINUX)
  std::string dynamicLibraryExtension("so") ;
#elif defined(OSI2PLATFORM_WINDOWS)
  std::string dynamicLibraryExtension("dll") ;
#endif

  if (libDir.empty()) // Check that the path is non-empty.
      return -1 ;

  platformServices_.invokeService = func ;

  Path dir_path(libDir) ;
  if (!dir_path.exists() || !dir_path.isDirectory())
      return -1 ;

  Directory::Entry e ;
  Directory::Iterator di(dir_path) ;
  while (di.next(e)) {
      Path full_path(dir_path + Path(e.path)) ;

      // Skip directories
      if (full_path.isDirectory())
	  continue ;

      // Skip files with the wrong extension
      std::string ext = std::string(full_path.getExtension()) ;
      if (ext != dynamicLibraryExtension)
	  continue ;

      // Ignore return value
      /*int32_t res = */
      loadOneLib(std::string(full_path)) ;
  }

  return 0 ;
}
#endif		// OSI2_IMPLEMENT_LOADALL


/*
  Unload a single library specified by name. The name must exactly match the
  name used to load the library. Issue a warning if the library isn't loaded.

  Returns:  1 if the library isn't loaded
	    0 if the library unloads successfully
	   -1 exit function failed
*/
int PluginManager::unloadOneLib (const std::string &libName,
                                 const std::string *dir)
{
  int result = 0 ;
/*
  If no directory is specified, consider both the default plugin directory and
  the innate `directory'. If a directory is specified, use that exclusively.
*/
  std::vector<std::string> plugDirs ;
  std::string pluginPath ;
  std::string innatePath ;
  std::string fullPath ;
  char dirSep = CoinFindDirSeparator() ;

  if (dir == nullptr || (dir->compare("") == 0)) {
    plugDirs = getPluginDirs() ;
    innatePath = getInnateDir() + dirSep + libName ;
  } else {
    plugDirs.push_back(*dir) ;
    innatePath = "" ;
  }
/*
  Look up the library's unique ID. It's possible (though unlikely) that we've
  been asked to unload an innate library. If we don't find anything at all,
  complain and return.
*/
  bool foundIt = false ;
  LibPathToIDMap::iterator lpiIter ;
  for (std::vector<std::string>::const_iterator iter = plugDirs.cbegin() ;
       iter != plugDirs.cend() && !foundIt ;
       iter++) {
    pluginPath = *iter+dirSep+libName ;
    lpiIter = libPathToIDMap_.find(pluginPath) ;
    if (lpiIter != libPathToIDMap_.end()) {
      foundIt = true ;
      fullPath = pluginPath ;
    }
  }
  if (!foundIt && innatePath != "") {
    lpiIter = libPathToIDMap_.find(innatePath) ;
    if (lpiIter != libPathToIDMap_.end()) {
      foundIt = true ;
      fullPath = innatePath ;
    }
  }
  if (!foundIt) {
    msgHandler_->message(PLUGMGR_LIBNOTFOUND, msgs_)
	<< libName << pluginPath << innatePath << CoinMessageEol ;
    return (1) ; }
  PluginUniqueID id = lpiIter->second ;
  LibraryInfo &lib = libraryMap_[id] ;
/*
  Step through the exact match map and remove any entries that are registered
  to this library. Regrettably, erase for a map is defined to invalidate
  the iterator pointing to the element. But not iterators pointing to other
  elements, so we need to advance prior to erasing.
*/
  APIRegMap::iterator rmIter = exactMatchMap_.begin() ;
  while (rmIter != exactMatchMap_.end()) {
    if (rmIter->second.id_ == id) {
      APIRegMap::iterator tmpIter = rmIter ;
      rmIter++ ;
      msgHandler_->message(PLUGMGR_APIUNREG, msgs_)
	  << tmpIter->first << fullPath << CoinMessageEol ;
      exactMatchMap_.erase(tmpIter) ;
    } else {
      rmIter++ ;
    }
  }
/*
  See if there's an entry in the wildcard vector. Vectors don't have the
  same problem as maps (erase returns a valid iterator), but there's only
  one entry so it's irrelevant.
*/
  for (APIRegVec::iterator rvIter = wildCardVec_.begin() ;
       rvIter != wildCardVec_.end() ;
       rvIter++) {
    if (rvIter->id_ == id) {
      msgHandler_->message(PLUGMGR_APIUNREG, msgs_)
	<< "wildcard"
	<< fullPath << CoinMessageEol ;
      wildCardVec_.erase(rvIter) ;
      break ;
    }
  }
  /*
    Execute the exit function for the library.
  */
  bool threwError = false ;
  ExitFunc func = lib.exitFunc_ ;
  platformServices_.pluginID_ = id ;
  platformServices_.ctrlObj_ = lib.ctrlObj_ ;
  try {
    result = (*func)(&platformServices_) ;
  } catch (...) {
    threwError = true ;
  }
  if (threwError || result != 0) {
    msgHandler_->message(PLUGMGR_LIBEXITFAIL, msgs_)
      << fullPath << CoinMessageEol ;
    result = -1 ;
  } else {
    msgHandler_->message(PLUGMGR_LIBEXITOK, msgs_)
      << fullPath << CoinMessageEol ;
  }
/*
  If this is an actual dynamic library, unload the library via the
  destructor. Then remove the map entry.
*/
  if (lib.isDynamic_) {
    msgHandler_->message(PLUGMGR_LIBCLOSE, msgs_)
      << fullPath << CoinMessageEol ;
    delete lib.dynLib_ ;
  }
  LibraryMap::iterator lmIter = libraryMap_.find(id) ;
  libraryMap_.erase(lmIter) ;
  libPathToIDMap_.erase(lpiIter) ;

  return (result) ;
}

/*
  This is a hook for compiled-in libraries to register their APIs with the
  PluginManager. The trick is to get a one-time registration. The strategy 
  is that a compiled-in API can create a single static instance of itself,
  with minimal initialisation, whose sole purpose in life is to invoke
  addPreloadLib and provide a pointer to the plugin's InitFunc. The
  PluginManager remembers the initFunc so that it can provide the illusion of
  unloading an innate plugin and reloading it.
*/

void PluginManager::addPreloadLib (std::string lib, InitFunc initFunc)
{
  char dirSep = CoinFindDirSeparator() ;
  std::string fullPath = dfltInnateDir_ + dirSep + lib ;
  preloadLibs_[fullPath] = initFunc ;
  PluginUniqueID id = initOneLib(fullPath,initFunc) ;
}

/*
  Execute the exit function for each library, then clear out the maps in the
  manager. Executing the destructor for the DynamicLibrary object will unload
  the library.

  \todo: The exit functions can throw? Why isn't there a catch block for all
	 the others (init function, etc.)
*/
int PluginManager::shutdown()
{
  int overallResult = 0 ;

  for (LibraryMap::iterator lmIter = libraryMap_.begin() ;
       lmIter != libraryMap_.end() ;
       lmIter++) {

    LibraryInfo &libInfo = lmIter->second ;
    ExitFunc func = libInfo.exitFunc_ ;
    platformServices_.pluginID_ = libInfo.id_ ;
    platformServices_.ctrlObj_ = libInfo.ctrlObj_ ;

    int result = 0 ;
    bool threwError = false ;
/*
  Figure out the full path for the library. For honest plugins we can get this
  from the dynamic library entry of libInfo. For innate plugins, we have to
  search through the libPathToIDMap_
*/
    std::string fullPath ;
    if (libInfo.isDynamic_) {
      fullPath = libInfo.dynLib_->getLibPath() ;
    } else {
      for (LibPathToIDMap::const_iterator lptiIter = libPathToIDMap_.begin() ;
           lptiIter != libPathToIDMap_.end() ;
	   lptiIter++) {
	if (lptiIter->second == libInfo.id_) {
	  fullPath = lptiIter->first ;
	  break ;
	}
      }
    }
    try {
      result = (*func)(&platformServices_) ;
    } catch (...) {
      threwError = true ;
    }
    if (threwError || result != 0) {
      msgHandler_->message(PLUGMGR_LIBEXITFAIL, msgs_)
	  <<  fullPath << CoinMessageEol ;
      overallResult-- ;
    } else {
      msgHandler_->message(PLUGMGR_LIBEXITOK, msgs_)
	  << fullPath << CoinMessageEol ;
    }
  }
/*
  Clear out the maps. Before we clear libraryMap_, go through and
  delete the DynamicLibrary objects; the destructor will unload the library.
  Otherwise, the maps contain actual objects and can simply be cleared.
*/
  for (LibraryMap::iterator lmIter = libraryMap_.begin() ;
       lmIter != libraryMap_.end() ;
       lmIter++) {
    LibraryInfo &libInfo = lmIter->second ;
    if (libInfo.isDynamic_) {
      DynamicLibrary *dynLib = libInfo.dynLib_ ;
      msgHandler_->message(PLUGMGR_LIBCLOSE, msgs_)
	  << dynLib->getLibPath() << CoinMessageEol ;
      delete dynLib ;
    }
  }
  libPathToIDMap_.clear() ;
  libraryMap_.clear() ;
  exactMatchMap_.clear() ;
  wildCardVec_.clear() ;

  return (overallResult) ;
}

/*
  Replace the current handler with a new handler. The current handler may or
  may not be our responsibility. If newHandler is null, create a default
  handler. (We can't be without one; the code isn't prepared for that.)
*/
void PluginManager::setMsgHandler (CoinMessageHandler *newHandler)
{
    if (dfltHandler_) {
        delete msgHandler_ ;
        msgHandler_ = nullptr ;
    }
    if (newHandler) {
        dfltHandler_ = false ;
        msgHandler_ = newHandler ;
    } else {
        dfltHandler_ = true ;
        msgHandler_ = new CoinMessageHandler() ;
        msgHandler_->setLogLevel(logLvl_) ;
    }
}


// ---------------------------------------------------------------

/*
  Construct an ObjectParams block to pass to the plugin.

  This method hides the details of constructing the object creation
  parameters.  In particular, it hides the business of finding the correct
  library state objects. Recall that the state object in PlatFormServices
  is state for the plugin as a whole, while the state object in ObjectParams
  is state for the API.
*/
ObjectParams *PluginManager::buildObjectParams (const std::string apiStr,
        const APIInfo &apiInfo)
{
  ObjectParams *objParms = new ObjectParams() ;
  objParms->apiStr_ = reinterpret_cast<const CharString *>(apiStr.c_str()) ;

  LibraryMap::iterator lmIter = libraryMap_.find(apiInfo.id_) ;
  LibraryInfo &libInfo = lmIter->second ;
  platformServices_.pluginID_ = apiInfo.id_ ;
  // Library control object
  platformServices_.ctrlObj_ = libInfo.ctrlObj_ ;
  objParms->platformServices_ = &platformServices_ ;
  // API control object
  objParms->ctrlObj_ = apiInfo.ctrlObj_ ;

  return (objParms) ;
}


void *PluginManager::createObject (const std::string &apiStr,
                                   PluginUniqueID &libID,
                                   IObjectAdapter &adapter)
{
/*
  "*" is not a valid type for createObject --- a specific type is needed.
*/
  if (apiStr == std::string("*")) {
    msgHandler_->message(PLUGMGR_APICREATEFAIL, msgs_)
	<< apiStr << "wildcard is invalid for createObject"
	<< CoinMessageEol ;
    return (nullptr) ;
  }
/*
  Check for an exact (string) match. If so, add the plugin's management object
  to the parameter block and ask for an object. If we're successful, we need
  one last step for a C plugin --- wrap it for C++ use.
*/
  APIRegMap::const_iterator apiIter =
      apiEntryExists(exactMatchMap_,apiStr,libID) ;
  if (apiIter != exactMatchMap_.end()) {
    const APIInfo &apiInfo = apiIter->second ;
    ObjectParams *objParms = buildObjectParams(apiStr, apiInfo) ;
    void *object = apiInfo.createFunc_(objParms) ;
    delete objParms ;
    if (object) {
      msgHandler_->message(PLUGMGR_APICREATEOK, msgs_)
	  << apiStr << "exact" << CoinMessageEol ;
      if (libID == 0) libID = apiInfo.id_ ;
      if (apiInfo.lang_ == Plugin_C)
	object = adapter.adapt(object, apiInfo.destroyFunc_) ;
      return (object) ;
    }
  }
/*
  No exact match. Try for a wildcard match. If some plugin volunteers an
  object, check to see if the plugin registered the API `on demand', so to
  speak. If so, we're good. If not, construct an exact match entry based on
  the wildcard apiInfo.
*/
  for (size_t i = 0 ; i < wildCardVec_.size() ; ++i) {
    APIInfo &wcInfo = wildCardVec_[i] ;
    if (libID && wcInfo.id_ != libID) continue ;
    ObjectParams *objParms = buildObjectParams(apiStr,wcInfo) ;
    void *object = wcInfo.createFunc_(objParms) ;
    delete objParms ;
    if (object) {
      msgHandler_->message(PLUGMGR_APICREATEOK, msgs_)
	  << apiStr << "wildcard" << CoinMessageEol ;
      const APIInfo *apiInfo = nullptr ;
      APIRegMap::const_iterator apiIter =
          apiEntryExists(exactMatchMap_,apiStr,libID) ;
      if (apiIter != exactMatchMap_.end()) {
        apiInfo = &apiIter->second ;
      } else {
	APIInfo tmp = wcInfo ;
	tmp.api_ = apiStr ;
	exactMatchMap_.insert(std::pair<std::string,APIInfo>(tmp.api_,tmp)) ;
	apiIter = apiEntryExists(exactMatchMap_,apiStr,libID) ;
	apiInfo = &apiIter->second ;
      }
      if (libID == 0) libID = apiInfo->id_ ;
      if (apiInfo->lang_ == Plugin_C)
	object = adapter.adapt(object, apiInfo->destroyFunc_) ;
      return (object) ;
    }
  }
  /*
    No plugin volunteered. We can't create this object.
  */
  msgHandler_->message(PLUGMGR_APICREATEFAIL, msgs_)
	  << apiStr << "no capable plugin" << CoinMessageEol ;
  return (nullptr) ;
}

#ifdef UNDEFINED

/*
  This needs to be updated to current code if we ever start to work with a
  straight C plugin.

  \todo Create a trivial C plugin for testing.  -- lh, 190711 --
*/
template <typename T, typename U>
T *PluginManager::createObject(const std::string & apiStr, IObjectAdapter<T, U> & adapter)
{
    // "*" is not a valid object type
    if (apiStr == std::string("*"))
        return (nullptr) ;

    // Prepare object params
    ObjectParams np ;
    np.apiStr = apiStr.c_str() ;
    np.platformServices = &ServiceProvider::getInstance() ;

    // Exact match found
    if (exactMatchMap_.find(apiStr) != exactMatchMap_.end()) {
        APIRegInfo & rp = exactMatchMap_[apiStr] ;
        IObject * object = createObject(rp, np, adapter) ;
        if (object) // great, it worked
            return object ;
    }

    for (Size i = 0; i < wildCardVec_.size(); ++i) {
        APIRegInfo & rp = wildCardVec_[i] ;
        IObject * object = createObject(rp, np, adapter) ;
        if (object) { // great, it worked
            // promote registration to exactMatc_
            // (but keep also the wild card registration for other object types)
            int32_t res = registerAPI(np.apiStr, &rp) ;
            if (res < 0) {
                THROW << "PluginManager::createObject(" << np.apiStr << "), registration failed" ;
                delete object ;
                return (nullptr) ;
            }
            return object ;
        }
    }

    // Too bad no one can create this objectType
    return (nullptr) ;
}
#endif


/*
  Destroy an object.

  Given that anything handed out by this manager is a C++ object, the client
  can invoke standard C++ delete and the destructor will do its work. But
  suppose that the plugin wants to do more than just invoke the destructor?
  This provides a way to invoke an arbitrary function. The only requirement is
  that the specified object be destroyed as part of the execution.
*/
int PluginManager::destroyObject (const std::string &apiStr,
                                  PluginUniqueID libID, void *victim)
{
  int result = 0 ;
  APIRegMap::const_iterator apiIter =
      apiEntryExists(exactMatchMap_,apiStr,libID) ;
  if (apiIter == exactMatchMap_.end()) {
    msgHandler_->message(PLUGMGR_APIDELFAIL, msgs_)
	<< apiStr << "no such API" << CoinMessageEol ;
    result = -1 ;
  } else {
    const APIInfo &apiInfo = apiIter->second ;
    ObjectParams *objParms = buildObjectParams(apiStr,apiInfo) ;
    result = apiInfo.destroyFunc_(victim,objParms) ;
    if (result < 0) {
      msgHandler_->message(PLUGMGR_APIDELFAIL, msgs_)
	  << apiStr << "DestroyFunc failed" << CoinMessageEol ;
    }
  }
  if (result >= 0)
    msgHandler_->message(PLUGMGR_APIDELOK, msgs_)
        << apiStr << CoinMessageEol ;

  return (result) ;
}


PlatformServices &PluginManager::getPlatformServices ()
{
    return (platformServices_) ;
}

/*
  Find the full path for the library specified by libID. Keep in mind that
  the library could be in initialisation.
*/
std::string PluginManager::getLibPath (PluginUniqueID libID)
{
  typedef LibPathToIDMap::const_iterator LPTIMI ;

  if (libID == libInInit_)
  { return (pathInInit_) ; }
  else
  { for (LPTIMI iter = libPathToIDMap_.begin() ;
	 iter != libPathToIDMap_.end() ;
	 iter++) {
    if (iter->second == libID) return (iter->first) ;
    }
  }
  return ("<library not loaded>") ;
}

