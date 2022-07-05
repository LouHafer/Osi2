
/*! \file Osi2RunParamsAPI_Imp.cpp

  Method definitions for an implementation of Osi2::RunParamsAPI.
*/

#include <string>
#include <map>


#ifndef Provide_RunParamsAPI_TemplateMethodDefs

#include "Osi2PluginManager.hpp"
#include "Osi2RunParamsAPI_Imp.hpp"


namespace {

/*
  Registration instance for RunParamsAPI.
*/

static Osi2::RunParamsAPI_Imp regObj(Osi2::RunParamsAPI::getAPIIDString()) ;

/*
  The plugin management methods. These need to be static so we can pass
  pointers to PluginManager. They need to be in the file-local namespace so
  that they don't pollute the global namespace. (Particularly if using the
  conventional names, as below.)
*/

/*
  The initial implementations of create and destroy simply invoke the
  constructor and destructor, respectively.
*/
static void *create (const Osi2::ObjectParams *params)
{
  std::string what = reinterpret_cast<const char *>(params->apiStr_) ;
  void *retval = nullptr ;
  if (what == Osi2::RunParamsAPI::getAPIIDString()) {
    retval = new Osi2::RunParamsAPI_Imp() ;
  }
  return (retval) ;
}

/*
  Everything RunParamsAPI knows about is descended from API, so we can simply
  invoke the virtual destructor.
*/
static int destroy (void *victim, const Osi2::ObjectParams *objParms)
{
  std::string what = reinterpret_cast<const char *>(objParms->apiStr_) ;
  int retval = 1 ;
  if (what == Osi2::RunParamsAPI::getAPIIDString()) {
    Osi2::API *api = static_cast<Osi2::API *>(victim) ;
    delete api ;
    retval = 0 ;
  }

  return (retval) ;
}

/*
  Nothing to be done to close out an innate plugin.
*/
static int cleanup (const Osi2::PlatformServices *services)
{
  return (0) ;
}

/*
  Plugin initialisation method. In theory, we should construct a plugin state
  object (ctrlObj) here, but so far I don't see the need.
*/
static Osi2::ExitFunc initPlugin (Osi2::PlatformServices *services)
{
  services->ctrlObj_ = nullptr ;
/*
  Register our API.
*/
  Osi2::APIRegInfo reginfo ;
  reginfo.version_.major_ = 1 ;
  reginfo.version_.minor_ = 0 ;
  reginfo.pluginID_ = services->pluginID_ ;
  reginfo.lang_ = Osi2::Plugin_CPP ;
  reginfo.ctrlObj_ = nullptr ;
  reginfo.createFunc_ = create ;
  reginfo.destroyFunc_ = destroy ;
  const char *apiName = Osi2::RunParamsAPI::getAPIIDString() ;
  int retval =
    services->registerAPI_(reinterpret_cast<const Osi2::CharString *>(apiName),
    			   &reginfo) ;

  if (retval < 0) { return (nullptr) ; }

  return (cleanup) ;
}

}  // end file-local namespace


namespace Osi2 {

/*
  Default constructor. All we need to do is initialise empty maps.
*/
RunParamsAPI_Imp::RunParamsAPI_Imp ()
  : intParams_(),
    dblParams_(),
    strParams_(),
    voidParams_()
{ /* nothing more to do */ }

/*
  Registration constructor

  The sole purpose of this constructor is to register the RunParams API
  with the plugin manager during program startup. The parameter exists to
  allow the compiler to choose the proper overload. Use this constructor
  to construct a static, file local instance that has absolutely no other use.
*/
RunParamsAPI_Imp::RunParamsAPI_Imp (std::string name)
{
  PluginManager *pluginMgr = &PluginManager::getInstance() ;

  pluginMgr->addPreloadLib(name,initPlugin) ;
}

RunParamsAPI_Imp::~RunParamsAPI_Imp ()
{ /* nothing more to do */ }

}  // end namespace Osi2

#endif  // Provide_RunParamsAPI_TemplateMethodDefs not defined


#ifdef Provide_RunParamsAPI_TemplateMethodDefs

/*
  Template methods. These must be present when instantiated. They'll be
  included within namespace Osi2, so don't wrap them again here.
*/

template <class ParamType>
void RunParamsAPI_Imp::addParam (ParamMap<ParamType> &paramMap,
			    std::string id, ParamType val)
{
  paramMap[id] = val ;
  return ;
}

template <class ParamType>
void RunParamsAPI_Imp::setParam (ParamMap<ParamType> &paramMap,
			    std::string id, ParamType val)
{
  typename ParamMap<ParamType>::iterator paramIter ;
  paramIter = paramMap.find(id) ;
  if (paramIter != paramMap.end()) {
    paramIter->second = val ;
  }
  return ;
}

template <class ParamType>
ParamType RunParamsAPI_Imp::getParam (const ParamMap<ParamType> &paramMap,
				      std::string id) const
{
  typename ParamMap<ParamType>::const_iterator paramIter ;
  paramIter = paramMap.find(id) ;
  if (paramIter != paramMap.end()) {
    return (paramIter->second) ;
  } else {
    return (ParamType()) ;
  }
}

template <class ParamType>
std::vector<std::string> RunParamsAPI_Imp::getParamIds
  (const ParamMap<ParamType> &paramMap) const
{
  std::vector<std::string> ids ;
  for (typename ParamMap<ParamType>::const_iterator paramIter =
  							paramMap.begin() ;
       paramIter != paramMap.end() ;
       paramIter++)
  { ids.push_back(paramIter->first) ; }
  return (ids) ;
}

#endif    // Provide_RunParamsAPI_TemplateMethodDefs defined

