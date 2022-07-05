/*
  Copyright 2019 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  This file contains the implementation of the methods for the Osi2 Parameter
  Managememt API.
*/

#include <cstring>

#include "Osi2Config.h"

#include "Osi2ParamMgmtAPI_Imp.hpp"


namespace {

/*
  The `registration' instance. The sole purpose of this instance of
  ParamMgmtAPI is to register the initFunc with the PluginManager. Since we
  won't have to look these up with dlsym, then can be plain static (hence
  mangled) names. No need for "C" linkage.
*/

Osi2::ParamMgmtAPI_Imp regInstance("ParamMgmt") ;

/*
  Plugin cleanup function.
*/
static int32_t cleanupPlugin (const Osi2::PlatformServices *services)
{
  return (0) ;
}

/*
  Plugin initialisation function.
*/

static Osi2::ExitFunc initPlugin (Osi2::PlatformServices *services)
{
  return (cleanupPlugin) ;
}


/*
  Utility method to check that the object to be enrolled actually supports
  parameter management.
*/
bool supportsMgmt (Osi2::API *obj)
{
  const char **apiVec = 0 ;
  int apiCnt = obj->getAPIs(apiVec) ;
  if (apiCnt == 0 || apiVec == nullptr) { return (false) ; }
  std::string paramBEID = Osi2::ParamBEAPI::getAPIIDString() ;

  for (int ndx = 0 ; ndx < apiCnt ; ndx++) {
    std::string apiName = apiVec[ndx] ;
    if (apiName == paramBEID) { return (true) ; }
  }

  return (false) ;
}
    



}  // end file-local namespace


namespace Osi2 {

/*
  Registration function to enroll an object for parameter management.
*/

bool ParamMgmtAPI_Imp::enroll (std::string objIdent, Osi2::API *enrollee)
{
/*
  Check that the object supports parameter management; if not, we're done.
*/
  if (!supportsMgmt(enrollee)) {
    msgHandler_->message(PMMGAPI_NOSUPPORT, msgs_)
        << objIdent << CoinMessageEol ;
    return (false) ; }
/*
  Check that the prefix is available (no clashes allowed).
*/
  if (indexMap_.find(objIdent) != indexMap_.end()) {
    msgHandler_->message(PMMGAPI_DUPIDENT, msgs_)
        << objIdent << CoinMessageEol ;
    return (false) ; }
/*
  Get the enrollee's back end parameter management object.
*/
  const char *paramBEID = Osi2::ParamBEAPI::getAPIIDString() ;
  void *tmp = enrollee->getAPIPtr(paramBEID) ;
  Osi2::ParamBEAPI *handler = reinterpret_cast<Osi2::ParamBEAPI *>(tmp) ;
/*
  Ask for the set of exported parameters and set up the index map entry.
*/
  ObjData entry ;
  entry.enrolledObject_ = enrollee ;
  entry.paramHandler_ = handler ;
  const char **exportedParams = nullptr ;
  int exportCnt = handler->reportParams(exportedParams) ;

  for (int ndx = 0 ; ndx < exportCnt ; ndx++) {
    const char *param = exportedParams[ndx] ;
    entry.paramNames_.push_back(param) ;
  }
  indexMap_[objIdent] = entry ;

  return (true) ;
}

/*
  Removal function, to remove an object from parameter management.
*/
bool ParamMgmtAPI_Imp::remove (std::string objIdent)
{
/*
  Check that the prefix is registered and complain if it isn't.
*/
  IndexMap::iterator iter = indexMap_.find(objIdent) ;
  if (iter == indexMap_.end()) {
    msgHandler_->message(PMMGAPI_UNREG, msgs_)
        << objIdent << CoinMessageEol ;
    return (false) ;
  }
  indexMap_.erase(iter) ;

  return (true) ;
}

/*
  Write the basic void* method, to get us started.
*/
bool ParamMgmtAPI_Imp::get (std::string objIdent, std::string param,
			     void *blob)
{
/*
  Check that the prefix is registered and complain if it isn't.
*/
  IndexMap::iterator iter = indexMap_.find(objIdent) ;
  if (iter == indexMap_.end()) {
    msgHandler_->message(PMMGAPI_UNREG, msgs_)
        << objIdent << CoinMessageEol ;
    return (false) ;
  }
  const ObjData &objData = iter->second ;
  ParamBEAPI *hdlr = objData.paramHandler_ ;
/*
  Check that the parameter is registered and complain if it isn't.
*/
  const char *pStr = param.c_str() ;
  bool foundIt = false ;
  for (int ndx = 0 ; ndx < objData.paramNames_.size() ; ndx++) {
    if (std::strcmp(pStr,objData.paramNames_[ndx]) == 0) {
      foundIt = true ;
      break ;
    }
  }
  if (!foundIt) {
    msgHandler_->message(PMMGAPI_PARMUNREG,msgs_)
      << objIdent << param << CoinMessageEol ;
    return (false) ;
  }
/*
  Invoke the parameter handler's get method.
*/
  if (!hdlr->get(param.c_str(),blob)) {
    msgHandler_->message(PMMGAPI_OPFAIL,msgs_)
      << objIdent << "get" << param << CoinMessageEol ;
    return (false) ;
  }
/*
  Success!
*/
  return (true) ;
}


/*
  We can get here directly or via convenience wrappers that handle standard
  types. The job here is to make sure the parameter is valid, find the
  parameter handler for the specified object, and invoke the ParamBE set
  method.
*/
bool ParamMgmtAPI_Imp::set (std::string objIdent, std::string param,
			    const void *blob)
{
/*
  Check that the object's identifying prefix is registered and complain if
  it isn't.
*/
  IndexMap::iterator iter = indexMap_.find(objIdent) ;
  if (iter == indexMap_.end()) {
    msgHandler_->message(PMMGAPI_UNREG, msgs_)
        << objIdent << CoinMessageEol ;
    return (false) ;
  }
  const ObjData &objData = iter->second ;
  ParamBEAPI *hdlr = objData.paramHandler_ ;
/*
  Check that the parameter is registered and complain if it isn't.
*/
  const char *pStr = param.c_str() ;
  bool foundIt = false ;
  for (int ndx = 0 ; ndx < objData.paramNames_.size() ; ndx++) {
    if (std::strcmp(pStr,objData.paramNames_[ndx]) == 0) {
      foundIt = true ;
      break ;
    }
  }
  if (!foundIt) {
    msgHandler_->message(PMMGAPI_PARMUNREG,msgs_)
      << objIdent << param << CoinMessageEol ;
    return (false) ;
  }
/*
  Invoke the parameter handler's set method.
*/
  if (!hdlr->set(param.c_str(),blob)) {
    msgHandler_->message(PMMGAPI_OPFAIL,msgs_)
      << objIdent << "set" << param << CoinMessageEol ;
    return (false) ;
  }
/*
  Success!
*/
  return (true) ;
}

/*
  Boilerplate: Constructors, destructors, & such like
*/

/*
  Default constructor
*/
ParamMgmtAPI_Imp::ParamMgmtAPI_Imp ()
  : logLvl_(7)
{
  indexMap_.clear() ;
  msgHandler_ = new CoinMessageHandler() ;
  msgs_ = ParamMgmtAPIMessages() ;
  msgHandler_->setLogLevel(logLvl_) ;
  msgHandler_->message(PMMGAPI_INIT, msgs_)
      << "default" << CoinMessageEol ;
}

/*
  Registration constructor

  The sole purpose of this constructor is to register the ParamMgmt API
  with the plugin manager during program startup. The parameter exists to
  allow the compiler to choose the proper overload. Use this constructor
  to construct a static, file local instance that has absolutely no other use.
*/
ParamMgmtAPI_Imp::ParamMgmtAPI_Imp (std::string name)
  : dfltHandler_(true),
    msgHandler_(0),
    msgs_(0),
    logLvl_(0)
{
#ifndef NDEBUG
/*
  Useful for debugging.
*/
  msgHandler_ = new CoinMessageHandler() ;
  msgs_ = ParamMgmtAPIMessages() ;
  msgHandler_->setLogLevel(7) ;
  msgHandler_->message(PMMGAPI_INIT, msgs_)
      << "registration" << CoinMessageEol ;
#endif

  PluginManager &pmgr = PluginManager::getInstance() ;
  pmgr.addPreloadLib(name,initPlugin) ;
}

/*
  Copy constructor
*/
ParamMgmtAPI_Imp::ParamMgmtAPI_Imp (const ParamMgmtAPI_Imp &rhs)
    : dfltHandler_(rhs.dfltHandler_),
      logLvl_(rhs.logLvl_)
{
  indexMap_ = rhs.indexMap_ ;
/*
  If this is our handler, make an independent copy. If it's the client's
  handler, we can't make an independent copy because the client won't know
  about it and won't delete it.
*/
  if (dfltHandler_) {
    msgHandler_ = new CoinMessageHandler(*rhs.msgHandler_) ;
  } else {
    msgHandler_ = rhs.msgHandler_ ;
  }
  msgs_ = rhs.msgs_ ;
  msgHandler_->setLogLevel(logLvl_) ;
  msgHandler_->message(PMMGAPI_INIT, msgs_)
    << "copy" << CoinMessageEol ;
}

/*
  Assignment
*/
ParamMgmtAPI_Imp &ParamMgmtAPI_Imp::operator= (const ParamMgmtAPI_Imp &rhs)
{
/*
  Self-assignment requires no work.
*/
  if (this == &rhs) return (*this) ;
/*
  Otherwise, get to it.
*/
  indexMap_ = rhs.indexMap_ ;
/*
  If it's our handler, we need to delete the old and replace with the new.
  If it's the user's handler, it's the user's problem. We just assign the
  pointer.
*/
  if (dfltHandler_) {
    delete msgHandler_ ;
    msgHandler_ = nullptr ;
  }
  dfltHandler_ = rhs.dfltHandler_ ;
  if (dfltHandler_) {
    msgHandler_ = new CoinMessageHandler(*rhs.msgHandler_) ;
  } else {
    msgHandler_ = rhs.msgHandler_ ;
  }
  msgs_ = rhs.msgs_ ;
  msgHandler_->setLogLevel(logLvl_) ;

  return (*this) ;
}


/*
  Destructor
*/
ParamMgmtAPI_Imp::~ParamMgmtAPI_Imp ()
{
/*
  If this is our handler, delete it. Otherwise it's the client's
  responsibility.
*/
  if (dfltHandler_) {
    delete msgHandler_ ;
    msgHandler_ = nullptr ;
  }
}

/*
  Factory constructor
*/
ParamMgmtAPI *ParamMgmtAPI_Imp::create ()
{
  ParamMgmtAPI *api = new ParamMgmtAPI_Imp() ;
  return (api) ;
}

/*
  Clone (factory copy)
*/
ParamMgmtAPI *ParamMgmtAPI_Imp::clone ()
{
  ParamMgmtAPI *api = new ParamMgmtAPI_Imp(*this) ;
  return (api) ;
}



} // end namespace Osi2

