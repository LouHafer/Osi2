/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$

  This file contains the unit test for the OSI2 PluginManager and OSI2 APIs.
*/

#include "CoinHelperFunctions.hpp"


#include "Osi2Config.h"
#include "Osi2nullptr.hpp"
#include "Osi2PluginManager.hpp"
#include "Osi2DynamicLibrary.hpp"
#include "Osi2ObjectAdapter.hpp"

#include "Osi2ParamMgmtAPI_Imp.hpp"
#include "Osi2ControlAPI_Imp.hpp"
#include "Osi2ProbMgmtAPI.hpp"
#include "Osi2Osi1API.hpp"
#include "Osi2ClpSimplexAPI.hpp"
#include "Osi2ClpLite_Wrap.hpp"
#include "Osi2ClpSolveParamsAPI.hpp"

#include "Osi2RunParamsAPI.hpp"

using namespace Osi2 ;

namespace {

/*
  Test the bare PluginManager API:
    * Initialise the PluginManager.
    * Load a plugin library.
    * Create ClpSimplex objects: exact match and wild card. Also check that
      we fail correctly for a nonexistent object.

  The test is (sort of) clp-specific, but only in the sense that the test clp
  plugin will return a ClpSimpex object via the wildcard mechanism when asked
  for a WildClpSimpex object.

  The test also checks that innate plugins work properly, using the ParamMgmt
  plugin as the test case. By the time the unit test executes, ParamMgmt will
  already be loaded.

  Over time, this test should be expanded and broken out into a separate file.
  Arguably, it should be a separate executable.
*/
int testPluginManager (const std::string libName)
{
  int errcnt = 0 ;

  std::cout << "Instantiating PluginManager." << std::endl ;

  PluginManager &plugMgr = PluginManager::getInstance() ;
/*
  Try to load ParamMgmt. This should provoke an error as it's an innate plugin
  and should have registered itself during program startup.
*/
  std::string innateName = "ParamMgmt" ;
  int retval = plugMgr.loadOneLib(innateName) ;
  if (retval == 0) {
    std::cout
      << "Unexpected success loading" << innateName << "." << std::endl ;
    std::cout
      << "Should report " << innateName << " already loaded." << std::endl ;
    std::cout
      << "Error code is " << retval << "." << std::endl ;
  }
/*
  Now unload ParamMgmt and reload it, to test that this works for an innate
  plugin.
*/
  retval = plugMgr.unloadOneLib(innateName) ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
      << "Apparent failure to unload " << innateName << "." << std::endl ;
    std::cout
      << "Error code is " << retval << "." << std::endl ;
    return (errcnt) ;
  }
  retval = plugMgr.loadOneLib(innateName) ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
      << "Apparent failure to load " << innateName << "." << std::endl ;
    std::cout
      << "Error code is " << retval << "." << std::endl ;
    return (errcnt) ;
  }
/*
  Now let's try to load the shim.
*/
  std::string dfltDir = plugMgr.getDfltPluginDir() ;
  char dirSep = CoinFindDirSeparator() ;
  std::string shimPath = dfltDir + dirSep + libName ;

  retval = plugMgr.loadOneLib(libName) ;

  if (retval != 0) {
    errcnt++ ;
    std::cout
      << "Apparent failure to load " << shimPath << "." << std::endl ;
    std::cout
      << "Error code is " << retval << "." << std::endl ;
    return (errcnt) ;
  }
/*
  Invoke createObject. If it works, try to invoke a nontrivial method. Which
  will fail, because exmip1 is not available in the current directory,
  but that's not the point. Then destroy the object.
*/
  DummyAdapter dummy ;
  PluginUniqueID libID = 0 ;
  ClpLite_Wrap *clpWrap =
    static_cast<ClpLite_Wrap *>(plugMgr.createObject("ClpSimplex",libID,dummy)) ;
  if (clpWrap == nullptr) {
      errcnt++ ;
      std::cout
	<< "Apparent failure to create a ClpLite_Wrap object." << std::endl ;
  } else {
      ClpSimplexAPI *clp =
          static_cast<ClpSimplexAPI *>(clpWrap->getAPIPtr("ClpSimplex")) ;
      clp->readMps("exmip1.mps", true) ;
      int retval = plugMgr.destroyObject("ClpSimplex", 0, clp) ;
      if (retval < 0) {
	  errcnt++ ;
	  std::cout
	    << "Apparent failure to destroy a ClpSimplex object."
	    << std::endl ;
      }
      clp = nullptr ;
  }
/*
  Ask for a nonexistent API and check that we (correctly) fail to provide
  one. We will fail twice: once on exact match, once on an attempt at a
  wildcard match.
*/
  libID = 0 ;
  ProbMgmtAPI *bogus =
    static_cast<ProbMgmtAPI *>(plugMgr.createObject("BogusAPI",libID,dummy)) ;
  if (bogus == nullptr) {
      std::cout
	<< "Apparent failure to create a BogusAPI object (expected)."
	<< std::endl ;
  } else {
      errcnt++ ;
      std::cout
	<< "Eh? We shouldn't be able to create a BogusAPI object!"
	<< std::endl ;
  }
  /*
    Check that we can create an object through the wildcard mechanism. The Clp
    shim doesn't register a WildClpSimplex object, so the only way we can get
    it is by the wild card mechanism (and a special hack in Osi2ClpShim,
    placed there for testing).
  */
  libID = 0 ;
  clpWrap = static_cast<ClpLite_Wrap *>(plugMgr.createObject("WildClpSimplex",
				   libID, dummy)) ;
  if (clpWrap == nullptr) {
      errcnt++ ;
      std::cout
	<< "Apparent failure to create a ClpLite_Wrap object."
	<< std::endl ;
  } else {
      ClpSimplexAPI *clp =
          static_cast<ClpSimplexAPI *>(clpWrap->getAPIPtr("ClpSimplex")) ;
      clp->readMps("exmip1.mps", true) ;
      int retval = plugMgr.destroyObject("WildClpSimplex", 0, clp) ;
      if (retval < 0) {
	std::cout
	  << "Apparent failure to destroy a WildProbMgmt object."
	  << std::endl ;
      }
      clp = nullptr ;
  }
  /*
    Unload the plugin library.
  */
  retval = plugMgr.unloadOneLib(libName) ;
  if (retval != 0) {
      errcnt++ ;
      std::cout
	      << "Apparent failure to unload " << shimPath << "." << std::endl ;
      std::cout
	      << "Error code is " << retval << "." << std::endl ;
  }
  /*
    Shut down the plugin manager. This will call the plugin library exit
    functions and unload the libraries.
  */
  std::cout
	  << "Shutting down PluginManager (executing exit functions)." << std::endl ;
  retval = plugMgr.shutdown() ;
  if (retval < 0) {
      errcnt++ ;
      std::cout
	      << "Apparent failure of PluginManager shutdown." << std::endl ;
  }

  return (errcnt) ;
}

/*
  Test the control API and various aspects of object creation and destruction.
  The method creates various objects, works with them, and finally destroys
  all of them. Part of the test is that there are multiple objects in
  existence, and multiple objects of the same API.
*/
int testControlAPI (const std::string &shortName,
		    const std::string &dfltSampleDir)
{
  int retval = 0 ;
  int errcnt = 0 ;
  std::vector<API *> apiObjs ;
  char dirSep = CoinFindDirSeparator() ;

  std::cout
      << "Attempting to instantiate a ControlAPI_Imp object." << std::endl ;
  ControlAPI_Imp ctrlAPI ;
  std::cout << "Log level is " << ctrlAPI.getLogLvl() << std::endl ;
/*
  Load the shim.
*/
  retval = ctrlAPI.load(shortName) ;
  if (retval != 0) {
    std::cout
	<< "Apparent failure to load " << shortName << "." << std::endl ;
    std::cout
	<< "Error code is " << retval << "." << std::endl ;
  }
/*
  Create a ProbMgmt API object and invoke a nontrivial method.
*/
  API *apiObj = nullptr ;
  retval = ctrlAPI.createObject(apiObj, "ProbMgmt") ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
	<< "Apparent failure to create a ProbMgmt object." << std::endl ;
  } else {
    apiObjs.push_back(apiObj) ;
    ProbMgmtAPI *clp = dynamic_cast<ProbMgmtAPI *>(apiObj) ;
    std::string exmip1Path = dfltSampleDir+dirSep+"brandy.mps" ;
    clp->readMps(exmip1Path.c_str(), true) ;
    clp->initialSolve() ;
  }
/*
  Create an Osi1 object and invoke a nontrivial method.  The clone test
  here is checking several things: First, that the cloned object is viable
  (i.e., the original clone method still works) and, second, that the
  control information block has been stripped. The second destroy should
  fail because there's no control information block.  (The ControlAPI
  didn't create this object and takes no responsibility for it!)  Later,
  the destruction of the original object should succeed, because it still
  possesses the control information block.
*/
  apiObj = nullptr ;
  retval = ctrlAPI.createObject(apiObj, "Osi1") ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
	<< "Apparent failure to create an Osi1 object." << std::endl ;
  } else {
    apiObjs.push_back(apiObj) ;
    Osi1API *osi = dynamic_cast<Osi1API *>(apiObj) ;
    std::string exmip1Path = dfltSampleDir+dirSep+"brandy.mps" ;
    osi->readMps(exmip1Path.c_str()) ;
    std::cout << "    cloning ... " << std::endl ;
    Osi1API *o2 = osi->clone() ;
    o2->initialSolve() ;
    if (o2->isProvenOptimal())
      std::cout << "    solved to optimality." << std::endl ;
    apiObj = o2 ;
    retval = ctrlAPI.destroyObject(apiObj) ;
    if (retval < 0) {
	std::cout
	    << "Apparent failure to destroy a cloned Osi1 object"
	    << " (expected)." << std::endl ;
	delete o2 ;
    } else {
      errcnt++ ;
      std::cout
	<< "Invoking ControlAPI::destroyObject on a cloned object"
	<< " should fail!" << std::endl ;
    }
  }
  /*
    Create a restricted ProbMgmt object and invoke a nontrivial method.
  */
  apiObj = nullptr ;
  retval = ctrlAPI.createObject(apiObj, "ProbMgmt", &shortName) ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
	<< "Apparent failure to create a restricted ProbMgmt object."
	<< std::endl ;
  } else {
    apiObjs.push_back(apiObj) ;
    ProbMgmtAPI *clp = dynamic_cast<ProbMgmtAPI *>(apiObj) ;
    clp->readMps("exmip1.mps", true) ;
  }
/*
  Now destroy all the objects.
*/
  for (std::vector<API *>::iterator iter = apiObjs.begin() ;
       iter != apiObjs.end() ;
       iter++) {
    API *apiObj = *iter ;
    int retval = ctrlAPI.destroyObject(apiObj) ;
    if (retval < 0) {
      errcnt++ ;
      std::cout
	<< "Apparent failure to destroy a ProbMgmt object." << std::endl ;
    }
  }
  apiObjs.clear() ;
/*
  Unload the shim.
*/
  retval = ctrlAPI.unload(shortName) ;
  if (retval != 0) {
    errcnt++ ;
    std::cout
	<< "Apparent failure to unload " << shortName << "." << std::endl ;
    std::cout
	<< "Error code is " << retval << "." << std::endl ;
  }
  return (errcnt) ;
}

int testParamMgmtAPI ()
{
  int errCnt = 0 ;

  ParamMgmtAPI_Imp mgmtAPI ;
  std::cout
      << "Instantiated a ParamMgmtAPI_Imp object "
      << std::hex << &mgmtAPI << std::dec << "." << std::endl ;

  ControlAPI_Imp ctrlAPI1 ;
  std::cout
      << "Instantiated a ControlAPI_Imp object "
      << std::hex << &ctrlAPI1 << std::dec << "." << std::endl ;
/*
  Check that the ControlAPI_Imp object thinks that it can implement
  ControlAPI.
*/
  std::vector<const char *>apiStrings ;
  apiStrings.push_back(Osi2::ControlAPI::getAPIIDString()) ;
  apiStrings.push_back(Osi2::ParamBEAPI::getAPIIDString()) ;
  apiStrings.push_back("BogusAPI") ;
  for (int ndx = 0 ; ndx < apiStrings.size() ; ndx++) {
    const char *apiStr = apiStrings[ndx] ;
    void *obj = ctrlAPI1.getAPIPtr(apiStr) ;
    if (obj == nullptr) {
      std::cout
	  << "Control object " << std::hex << &ctrlAPI1 << std::dec
	  << " believes it does not implement " << apiStr << "."
	  << std::endl ;
    } else {
      std::cout
	  << "Control object " << std::hex << &ctrlAPI1 << std::dec
	  << " believes it implements " << apiStr << " with "
	  << std::hex << obj << std::dec << "."
	  << std::endl ;
    }
  }

  std::string ctrlAPI1ID = "CtrlAPI1" ;
  std::cout
      << "Enrolling CtrlAPI1 as \"" << ctrlAPI1ID << "\"" << std::endl ;
  if (!mgmtAPI.enroll(ctrlAPI1ID,&ctrlAPI1)) {
    std::cout
        << "Failed to enroll CtrlAPI1 (" << ctrlAPI1ID << ")" << std::endl ;
    errCnt++ ;
  }
/*
  Attempt to retrieve a parameter value.
*/
  {
    int blob ;
    const char *paramID = "log level" ;
    if (mgmtAPI.get(ctrlAPI1ID,paramID,&blob))
    { std::cout
	<< "The value of " << ctrlAPI1ID << ":" << paramID
	<< " is " << blob << std::endl ;
    } else {
      std::cout
	<< "Failed to get the value of " << ctrlAPI1ID << ":" << paramID
	<< std::endl ;
    }
  /*
    Attempt to set the parameter and retrieve it again.
  */
    blob = blob-2 ;
    if (!mgmtAPI.set(ctrlAPI1ID,paramID,&blob)) {
      std::cout
	<< "Failed to set the value of " << ctrlAPI1ID << ":" << paramID
	<< std::endl ;
    }
    if (mgmtAPI.get(ctrlAPI1ID,paramID,&blob))
    { std::cout
	<< "The value of " << ctrlAPI1ID << ":" << paramID
	<< " is " << blob << std::endl ;
    } else {
      std::cout
	<< "Failed to get the value of " << ctrlAPI1ID << ":" << paramID
	<< std::endl ;
    }
  }
/*
  Attempt to retrieve a parameter value.
*/
  {
    std::string blob ;
    const char *paramID = "DfltPlugDir" ;
    ctrlAPI1.setDfltPluginDir("/home/Coin.Mix/Split-FedGCC/lib") ;
    if (mgmtAPI.get(ctrlAPI1ID,paramID,&blob))
    { std::cout
	<< "The value of " << ctrlAPI1ID << ":" << paramID
	<< " is " << blob << std::endl ;
    } else {
      std::cout
	<< "Failed to get the value of " << ctrlAPI1ID << ":" << paramID
	<< std::endl ;
    }
  }
/*
  Use the ControlAPI object to load ClpLite. Find the ClpSimplex object and
  register it with parameter management.
*/
  std::string shortName = "ClpLite" ;
  int retval = ctrlAPI1.load(shortName,"libOsi2ClpShim.so") ;
  if (retval != 0) {
    errCnt++ ;
    std::cout
      << "Apparent failure to load " << shortName << "." << std::endl ;
    std::cout
      << "Error code is " << retval << "." << std::endl ;
    return (errCnt) ;
  }
  API *apiObj = nullptr ;
  retval = ctrlAPI1.createObject(apiObj,"ClpSimplex") ;
  if (retval != 0) {
    errCnt++ ;
    std::cout
      << "Apparent failure to create a ClpSimplex wrapper object."
      << std::endl ;
    return (errCnt) ;
  }
  ClpLite_Wrap *wrap = static_cast<ClpLite_Wrap *>(apiObj) ;
  ClpSimplexAPI *clp =
      static_cast<ClpSimplexAPI *>(wrap->getAPIPtr("ClpSimplex")) ;
  clp->readMps("../../share/coin-or-sample/exmip1.mps") ;
  mgmtAPI.enroll("ClpLite",clp) ;
  double dblblob = 0.0 ;
  mgmtAPI.get("ClpLite","primal tolerance",&dblblob) ;
  std::cout
    << "The primal zero tolerance is " << dblblob << "." << std::endl ;
/*
  Work with the ClpSolveParams object.
*/
  ClpSolveParamsAPI *clpsolve =
      static_cast<ClpSolveParamsAPI *>(wrap->getAPIPtr("ClpSolveParams")) ;
  const char **idents = nullptr ;
  int apiCnt = clpsolve->getAPIs(idents) ;
  std::cout << "  clpsolve implements " << apiCnt << " APIs:" ;
  for (int ndx = 0 ; ndx < apiCnt ; ndx++)
    std::cout << " " << idents[ndx] ;
  std::cout << std::endl ;
  mgmtAPI.enroll("ClpP",clpsolve) ;
  std::string inStrBlob = "Hi, mom!" ;
  std::string outStrBlob = std::string() ;
  mgmtAPI.set("ClpP","problem name",&inStrBlob) ;
  mgmtAPI.get("ClpP","problem name",&outStrBlob) ;
  std::cout << " and the name is " << outStrBlob << std::endl ;

  return (errCnt) ;
}

int testRunParamsAPI ()

{ int errCnt = 0 ;
  std::map<std::string,API *> extantObjs ;

  ControlAPI_Imp ctrlAPI ;
  std::cout
      << "Instantiated a ControlAPI_Imp object "
      << std::hex << &ctrlAPI << std::dec << "." << std::endl ;
/*
  Reload the RunParams plugin (unloaded as part of the plugin manager test).
*/
  std::string rpDir = "" ;
  std::string rpShortName = "RunParams" ;
  int retval = ctrlAPI.load(rpShortName,rpShortName,&rpDir) ;
  if (retval != 0) {
    std::cout
      << "Error " << retval << " loading RunParams plugin."
      << std::endl ;
    errCnt++ ;
    return (errCnt) ;
  } else {
    std::cout << "Loaded RunParams plugin." << std::endl ;
  }
/*
  Instantiate a RunParams API object.
*/
  API *apiObj = nullptr ;
  std::string rp1Name = "RunParams #1" ;
  retval = ctrlAPI.createObject(apiObj,RunParamsAPI::getAPIIDString()) ;
  if (retval != 0) {
    std::cout
      << "Error " << retval << " creating RunParamsAPI object."
      << std::endl ;
    errCnt++ ;
    return (errCnt) ;
  } else {
    std::cout
      << "Created RunParamsAPI object " << std::hex << apiObj << std::dec
      << "." << std::endl ;
    extantObjs[rp1Name] = apiObj ;
  }
  RunParamsAPI *rpObj = dynamic_cast<RunParamsAPI *>(apiObj) ;
/*
  Instantiate a ClpSimplex API object and load in a model.
*/
  apiObj = nullptr ;
  std::string clpName = "ClpC" ;
  retval = ctrlAPI.createObject(apiObj,ClpSimplexAPI::getAPIIDString()) ;
  if (retval != 0) {
    std::cout
      << "Error " << retval << " creating ClpSimplexAPI object."
      << std::endl ;
    errCnt++ ;
    return (errCnt) ;
  } else {
    std::cout
      << "Created ClpSimplexAPI object " << std::hex << apiObj << std::dec
      << "." << std::endl ;
    extantObjs[clpName] = apiObj ;
  }
  ClpLite_Wrap *wrap = static_cast<ClpLite_Wrap *>(apiObj) ;
  ClpSimplexAPI *clpObj =
    static_cast<ClpSimplexAPI *>(wrap->getAPIPtr("ClpSimplex")) ;
/*
  Try to load a RunParams object. See what it holds.
*/
  clpObj->exposeParams(*rpObj) ;
  rpObj->setStrParam("problem name","pilot") ;
  rpObj->setDblParam("obj sense",1.0) ;
  std::vector<std::string> paramNames = rpObj->getDblParamIds() ;
  for (std::vector<std::string>::const_iterator iter = paramNames.begin() ;
       iter != paramNames.end() ;
       iter++) {
    std::cout
      << "  parameter " << *iter << " has value "
      << rpObj->getDblParam(*iter) << "." << std::endl ;
  }
  std::cout
    << "The problem name is " << rpObj->getStrParam("problem name")
    << "." << std::endl ;
  clpObj->readMps("../../share/coin-or-netlib/pilot.mps") ;
  clpObj->loadParams(*rpObj) ;
  clpObj->initialSolve() ;
  std::cout
    << "Solve took " << clpObj->numberIterations() << " iterations."
    << std::endl ;
  rpObj->setDblParam("obj sense",-1.0) ;
  clpObj->readMps("../../share/coin-or-netlib/pilot.mps") ;
  clpObj->loadParams(*rpObj) ;
  clpObj->initialSolve() ;
  std::cout
    << "Solve took " << clpObj->numberIterations() << " iterations."
    << std::endl ;
/*
  Destroy the objects we've created.
*/
  for (std::map<std::string,API *>::iterator iter = extantObjs.begin() ;
       iter != extantObjs.end() ;
       iter++) {
    std::cout
      << "Attempting to destroy " << iter->first << "." << std::endl ;
    API *tmp = iter->second ;
    retval = ctrlAPI.destroyObject(iter->second) ;
    if (retval != 0) {
      std::cout
	<< "Error " << retval << " destroying RunParamsAPI object "
	<< std::hex << tmp << std::dec << "."
	<< std::endl ;
      errCnt++ ;
      return (errCnt) ;
    } else {
      std::cout
	<< "Destroyed RunParamsAPI object " << std::hex << tmp << std::dec
	<< "." << std::endl ;
    }
  }

  return errCnt ;
}

} // end unnamed file-local namespace


int main(int argC, char* argV[])
{

  std::cout << "START UNIT TEST" << std::endl ;

  std::string dfltSampleDir = "../../Data/Sample" ;
/*
  Test the bare PluginManager. There's no sense proceeding to the API
  tests if the PluginManager isn't working.
*/
  std::cout << "Testing bare PluginManager." << std::endl ;
  int retval = testPluginManager("libOsi2ClpShim.so") ;
  std::cout
    << "End test of bare PluginManager, " << retval << " errors."
    << std::endl << std::endl ;
  if (retval != 0) {
      std::cout
	<< "Aborting unitTest; errors in PluginManager." << std::endl ;
      return (retval) ;
  }
/*
  Now let's try the Osi2 control API.
*/
  typedef std::pair<std::string,int> TestVec ;
  std::vector<TestVec> solvers ;
#if 0
  Temporarily disable while ClpLite is under renovation -- lh, 190802 --
  solvers.push_back(TestVec("clp",1)) ;
#endif
  solvers.push_back(TestVec("clpHeavy",0)) ;
# ifdef COIN_HAS_OSIGLPK
  solvers.push_back(TestVec("glpkHeavy",2)) ;
# endif
  std::vector<TestVec>::const_iterator iter ;
  int totalErrs = 0 ;
  for (iter = solvers.begin() ; iter != solvers.end() ; iter++) {
    std::string solverName = iter->first ;
    std::cout
      << "Testing ControlAPI (" << solverName << ")." << std::endl ;
    retval = testControlAPI(solverName,dfltSampleDir) ;
    std::cout
      << "End test of ControlAPI (" << solverName << "), "
      << retval << " errors, expected " << iter->second << "."
      << std::endl << std::endl ;
    totalErrs += retval-(iter->second) ;
  }
/*
  Test the parameter management API.
*/
  int errCnt = 0 ;
  int expectedErrs = 0 ;
  std::cout << "Testing ParamMgmtAPI." << std::endl ;
  errCnt = testParamMgmtAPI() ;
  std::cout
    << "End test of ParamMgmtAPI, " << errCnt << " errors, expected "
    << expectedErrs << "." << std::endl ;
/*
  Test the RunParams API.
*/
  std::cout << std::endl << std::endl ;
  std::cout << "Testing RunParamsAPI." << std::endl ;
  errCnt = testRunParamsAPI() ;
  expectedErrs = 0 ;
  std::cout
    << "End test of RunParamsAPI, " << errCnt << " errors, expected "
    << expectedErrs << "." << std::endl ;
/*
  Shut down the plugin manager. This will call the plugin library exit
  functions and unload the libraries.
*/
  PluginManager &plugMgr = PluginManager::getInstance() ;
  std::cout
    << "Shutting down plugins (executing exit functions)." << std::endl ;
  plugMgr.shutdown() ;

  std::cout << "END UNIT TEST" << std::endl ;

  return (totalErrs) ;

}

